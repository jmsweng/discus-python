"""Regression tests for native DISCUS 3D output wiring and crystallographic validity.

Real C3/C4 HDF5 fixtures were produced by DISCUS_SUITE 6.23.01.
A tiny fake executable is used to test runner plumbing without claiming a
fresh Fortran calculation occurred in CI.
"""
from __future__ import annotations

import os
from pathlib import Path
import shutil
import subprocess
import sys

import h5py
import numpy as np
import pytest

from discuspdf import Diffuse3DRunner, FFTValidationError
from discuspdf.discus_output import DiscusVolume, delta_pdf_fractional, read_discus_volume


ROOT = Path(__file__).resolve().parents[1]
FIX = ROOT / "tests" / "fixtures" / "discus_6.23.01"


def make_native(path, *, values=None, corner=None, inc=None,
                lengths=(5.0, 6.0, 7.0), angles=(90.0, 105.0, 90.0)):
    n = 9
    if values is None:
        values = np.ones((n, n, n))
    if corner is None:
        corner = [-1., -1., -1.]
    if inc is None:
        inc = np.eye(3) * .25
    with h5py.File(path, "w") as f:
        g = f.create_group("entry/data")
        for name, text in {
            "audit_conform_dict_name": "Disorder unified data",
            "coordinate_unit": "basecell_fractional",
            "data_type_axes": "hkl",
            "data_type_reciprocal": "reciprocal",
            "data_type_content": "intensity",
            "data_type_with_bragg": "bragg_subtracted",
        }.items():
            g.create_dataset(name, data=np.bytes_(text))
        for name, value in {
            "data_values": values, "data_dimension": np.array(values.shape),
            "data_axes": np.array([1, 2, 3]),
            "data_corner": corner, "data_increment_vector": inc,
            "unit_cell_lengths": lengths, "unit_cell_angles": angles,
        }.items():
            g.create_dataset(name, data=value)


@pytest.mark.parametrize("lengths,angles,error", [
    ((0, 6, 7), (90, 90, 90), "lengths"),
    ((-5, 6, 7), (90, 90, 90), "lengths"),
    ((5, 6, 7), (0, 90, 90), "angles"),
    ((5, 6, 7), (180, 90, 90), "angles"),
    ((5, 6, 7), (10, 10, 179), "positive definite"),
    ((5, 6, 7), (float("nan"), 90, 90), "angles"),
])
def test_reject_unphysical_cell(tmp_path, lengths, angles, error):
    p = tmp_path / "bad.h5"
    make_native(p, lengths=lengths, angles=angles)
    with pytest.raises(ValueError, match=error):
        read_discus_volume(p)


@pytest.mark.parametrize("case,corner,inc,error", [
    ("badcorner", [0, 0], np.eye(3)*.25, "data_corner"),
    ("nan", [0, np.nan, 0], np.eye(3)*.25, "data_corner"),
    ("zero", [-1]*3, np.zeros((3, 3)), "singular"),
    ("flat", [-1]*3, np.diag([.25, .25, 0]), "singular"),
    ("shape", [-1]*3, np.eye(2), "data_increment_vector"),
])
def test_reject_malformed_reciprocal_grid(tmp_path, case, corner, inc, error):
    p = tmp_path / (case + ".h5")
    make_native(p, corner=corner, inc=inc)
    with pytest.raises(ValueError, match=error):
        read_discus_volume(p)


def test_directly_constructed_volume_is_validated():
    with pytest.raises(ValueError, match="lengths"):
        DiscusVolume(
            np.ones((3, 3, 3)), np.zeros(3), np.eye(3)*.5,
            np.array([4., -1., 6.]), np.array([90., 90., 90.]), "bragg_subtracted")


def test_fractional_fft_error_preserves_data(tmp_path):
    n = 9
    values = np.zeros((n, n, n))
    values[n//2 + 1, n//2, n//2] = 1.0  # no Friedel partner
    p = tmp_path / "asym.h5"
    make_native(p, values=values)
    vol = read_discus_volume(p)
    with pytest.raises(FFTValidationError) as caught:
        delta_pdf_fractional(vol, imaginary_tolerance=1e-8)
    result = caught.value.result
    assert result.complex_volume.shape == (n, n, n)
    assert result.imaginary_ratio > 1e-8
    assert result.friedel_relative_error > 0
    assert len(result.u_axes) == 3
    assert result.provenance["sha256"] == vol.provenance["sha256"]
    with pytest.warns(RuntimeWarning):
        alternative = delta_pdf_fractional(vol, imaginary_tolerance=1e-8,
                                           raise_on_imaginary=False)
    assert np.allclose(alternative.complex_volume, result.complex_volume)
    with pytest.raises(ValueError, match="imaginary_tolerance"):
        delta_pdf_fractional(vol, imaginary_tolerance=-1)


def make_fake_runner(tmp_path, output_name, python_body):
    script = tmp_path / "fake_engine.py"
    script.write_text(python_body, encoding="utf-8")
    macro = tmp_path / "input.mac"
    macro.write_text("# staged simulation macro\n", encoding="utf-8")
    runner = Diffuse3DRunner(runs_directory=tmp_path / "runs",
                             output_filename=output_name,
                             launcher=(sys.executable, str(script)))
    return runner, macro


def test_native_runner_loads_real_fixture_and_preserves_metadata(tmp_path):
    output = "C3_total.h5"
    body = "import shutil\nshutil.copyfile(" + repr(str(FIX / output)) + ", " + repr(output) + ")\n"
    runner, macro = make_fake_runner(tmp_path, output, body)
    job = runner.prepare(macro)
    done = runner.execute(job, native_nexus=True)
    assert done.output.name == output
    volume = done.load_discus_volume()
    assert volume.values.shape == (33, 41, 49)
    assert volume.with_bragg == "bragg_present"
    status = __import__("json").loads((done.directory / "status.json").read_text())
    assert status["status"] == "complete"
    assert status["output_filename"] == output
    assert status["output_validation"]["native_nexus"] is True
    assert status["output_validation"]["physical_validation"] is False
    with pytest.raises(ValueError, match="bragg_present"):
        delta_pdf_fractional(volume)


def test_native_runner_reports_wrong_schema_failed(tmp_path):
    body = ("import h5py\n"
            "with h5py.File('native.h5', 'w') as f:\n"
            "    f.create_dataset('intensity', data=[1,2,3])\n")
    runner, macro = make_fake_runner(tmp_path, "native.h5", body)
    job = runner.prepare(macro)
    with pytest.raises(RuntimeError, match="Invalid 3D output"):
        runner.execute(job, native_nexus=True)
    import json
    status = json.loads((job.directory / "status.json").read_text())
    assert status["status"] == "failed"
    assert "common/nexus" in status["reason"]


@pytest.mark.parametrize("name", [
    "../escape.h5", "nested/foo.h5", r"C:\\escape.h5", "run.mac", ".", "",
])
def test_output_filename_is_restricted(name):
    with pytest.raises(ValueError, match="output_filename"):
        Diffuse3DRunner(output_filename=name)


def test_conflicting_native_validation_modes(tmp_path):
    from discuspdf import OutputSpec
    runner, macro = make_fake_runner(tmp_path, "native.h5", "print('unused')\n")
    job = runner.prepare(macro)
    with pytest.raises(ValueError, match="OR"):
        runner.execute(job, native_nexus=True,
                       output_spec=OutputSpec("i", ("qx", "qy", "qz")))


def test_probe_case_arguments_with_stub_discus(tmp_path):
    """Run the real shell entrypoint with an HDF5-producing stub engine.

    This is the original test, updated to use the default stub that writes
    valid HDF5 files so the strict output-integrity checks pass.
    """
    copied, env = _probe_env(tmp_path, None)
    selected = tmp_path / "probe_selected"
    subprocess.run(["bash", str(copied / "run_cases.sh"), str(selected), "C3", "C4"],
                   env=env, check=True, capture_output=True, text=True)
    assert (selected / "C3.mac").exists()
    assert (selected / "C4.mac").exists()
    assert not (selected / "C0.mac").exists()
    default = tmp_path / "probe_default"
    subprocess.run(["bash", str(copied / "run_cases.sh"), str(default)],
                   env=env, check=True, capture_output=True, text=True)
    assert sorted(x.stem for x in default.glob("C?.mac")) == ["C0", "C1", "C2", "C3", "C4"]
    assert not (default / "probe_default.mac").exists()
    unknown = subprocess.run(["bash", str(copied / "run_cases.sh"), str(tmp_path / "invalid"), "BAD"],
                             env=env, capture_output=True, text=True)
    assert unknown.returncode == 2
    assert "Unknown case" in unknown.stderr


# --- Gate 2: run_cases.sh safety regression tests --------------------------
#
# These exercise the real shell entrypoint with stub executables (a no-op
# engine, a failing engine, and a macro-producing engine) so that the
# fail-fast, overwrite-safety, and per-macro output-integrity behavior is
# pinned without requiring a live DISCUS installation.

def _probe_env(tmp_path, engine_body):
    """Copy the spike dir + examples and install a stub discus_suite on PATH.

    The default engine_body writes a minimal valid HDF5 file for each ``outf``
    name declared in the macro, so the inspector and HDF5 checks pass. Pass a
    custom body to simulate failure modes (nonzero exit, missing output, etc.).
    """
    source = ROOT / "spike_discus_3d"
    copied = tmp_path / "spike_discus_3d"
    shutil.copytree(source, copied)
    (tmp_path / "examples").mkdir(exist_ok=True)
    shutil.copy2(ROOT / "examples" / "gaas.cll", tmp_path / "examples" / "gaas.cll")
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir(exist_ok=True)
    engine = fake_bin / "discus_suite"
    if engine_body is None:
        # Default: write a minimal valid HDF5 file for each outf in the macro.
        # Use the test's Python interpreter (sys.executable) which has h5py.
        engine_body = (
            "#!/usr/bin/env python3\n"
            "import re, sys\n"
            "import h5py\n"
            "macro = sys.argv[-1]\n"
            "text = open(macro).read()\n"
            "for name in re.findall(r'outf\\s+(\\S+)', text):\n"
            "    with h5py.File(name, 'w') as f:\n"
            "        f.create_dataset('data', data=[1.0, 2.0, 3.0])\n"
            "    print('stub engine wrote', name)\n"
        )
        # The stub must run with the same Python that has h5py installed.
        # Replace the shebang with the test's interpreter.
        engine_body = engine_body.replace("#!/usr/bin/env python3", f"#!{sys.executable}")
    engine.write_text(engine_body, encoding="utf-8")
    engine.chmod(0o755)
    env = dict(os.environ, PATH=str(fake_bin) + os.pathsep + os.environ.get("PATH", ""),
               PYTHON=sys.executable)
    return copied, env


def _run_probe(copied, env, *args):
    return subprocess.run(["bash", str(copied / "run_cases.sh"), *args],
                          env=env, capture_output=True, text=True)


def test_probe_explicit_c3_c4_selection(tmp_path):
    copied, env = _probe_env(tmp_path, None)
    out = tmp_path / "sel_c3c4"
    r = _run_probe(copied, env, str(out), "C3", "C4")
    assert r.returncode == 0, r.stderr
    assert (out / "C3.mac").exists()
    assert (out / "C4.mac").exists()
    # Only the selected cases are generated.
    assert not (out / "C0.mac").exists()
    assert not (out / "C1.mac").exists()


def test_probe_default_case_selection(tmp_path):
    copied, env = _probe_env(tmp_path, None)
    out = tmp_path / "default"
    r = _run_probe(copied, env, str(out))
    assert r.returncode == 0, r.stderr
    assert sorted(x.stem for x in out.glob("C?.mac")) == ["C0", "C1", "C2", "C3", "C4"]
    assert not (out / "C5.mac").exists()


def test_probe_unknown_case_fails_before_writing(tmp_path):
    copied, env = _probe_env(tmp_path, "#!/bin/sh\nexit 0\n")
    out = tmp_path / "invalid"
    r = _run_probe(copied, env, str(out), "BAD")
    assert r.returncode == 2
    assert "Unknown case" in r.stderr
    # No files were written for the bad case.
    assert not (out / "BAD.mac").exists()


def test_probe_missing_executable_fails(tmp_path):
    copied, env = _probe_env(tmp_path, "#!/bin/sh\nexit 0\n")
    # Point the engine at a binary that does not exist.
    env = dict(env, DISCUS_ENGINE="definitely_not_a_binary_xyz")
    out = tmp_path / "noengine"
    r = _run_probe(copied, env, str(out), "C0")
    assert r.returncode != 0
    assert "engine" in r.stderr.lower()


def test_probe_engine_nonzero_exit_fails(tmp_path):
    copied, env = _probe_env(tmp_path, "#!/bin/sh\nexit 7\n")
    out = tmp_path / "badexit"
    r = _run_probe(copied, env, str(out), "C0")
    assert r.returncode != 0
    assert "C0" in r.stderr
    assert "7" in r.stderr


def test_probe_engine_exit0_missing_output_fails(tmp_path):
    # Engine exits 0 but writes no HDF5 at all.
    copied, env = _probe_env(tmp_path, "#!/bin/sh\nexit 0\n")
    out = tmp_path / "nooutput"
    r = _run_probe(copied, env, str(out), "C0")
    assert r.returncode != 0
    assert "missing" in r.stderr.lower() or "empty" in r.stderr.lower()


def test_probe_malformed_hdf5_fails(tmp_path):
    # Engine writes a non-HDF5 file with the expected name.
    body = ("#!/bin/sh\n"
            "echo 'not-hdf5-bytes' > C0_total.h5\n"
            "echo 'not-hdf5-bytes' > C0_diffuse.h5\n"
            "echo 'not-hdf5-bytes' > C0_diffuse_3dpdf.h5\n"
            "exit 0\n")
    copied, env = _probe_env(tmp_path, body)
    out = tmp_path / "malformed"
    r = _run_probe(copied, env, str(out), "C0")
    assert r.returncode != 0
    assert "HDF5" in r.stderr or "hdf5" in r.stderr


def test_probe_refuses_occupied_directory(tmp_path):
    copied, env = _probe_env(tmp_path, "#!/bin/sh\nexit 0\n")
    out = tmp_path / "occupied"
    out.mkdir()
    (out / "prior_evidence.txt").write_text("do not delete")
    r = _run_probe(copied, env, str(out), "C0")
    assert r.returncode == 3
    assert "not empty" in r.stderr
    # Prior evidence is preserved.
    assert (out / "prior_evidence.txt").exists()


def test_probe_overwrite_flag_clobbers(tmp_path):
    copied, env = _probe_env(tmp_path, None)
    out = tmp_path / "overwrite_me"
    out.mkdir()
    (out / "prior_evidence.txt").write_text("will be clobbered")
    r = _run_probe(copied, env, str(out), "--overwrite", "C0")
    assert r.returncode == 0, r.stderr
    # Prior evidence was removed by the explicit overwrite.
    assert not (out / "prior_evidence.txt").exists()
    assert (out / "C0.mac").exists()


# --- Gate 2b: run_cases.sh final engineering fixes --------------------------
#
# These tests cover the three defects found during independent review:
#   1. Multi-case macro hashing (bash array expansion)
#   2. Archive creation (outside the archived tree, no self-inclusion)
#   3. --overwrite guard (protected paths, symlink aliases)

def test_probe_hash_manifest_contains_all_selected_macros(tmp_path):
    """The SHA-256 manifest must contain every selected macro exactly once."""
    copied, env = _probe_env(tmp_path, None)
    out = tmp_path / "hash_multi"
    r = _run_probe(copied, env, str(out), "C3", "C4")
    assert r.returncode == 0, r.stderr
    hashes = (out / "hashes.txt").read_text()
    # Both selected macros must appear.
    assert "C3.mac" in hashes
    assert "C4.mac" in hashes
    # Unselected macros must NOT appear.
    assert "C0.mac" not in hashes
    assert "C1.mac" not in hashes
    assert "C2.mac" not in hashes
    assert "C5.mac" not in hashes
    # Shared inputs must appear.
    assert "gaas.cll" in hashes
    assert "mono2.cll" in hashes


def test_probe_hash_manifest_single_case(tmp_path):
    """A single-case selection produces a manifest with exactly that macro."""
    copied, env = _probe_env(tmp_path, None)
    out = tmp_path / "hash_single"
    r = _run_probe(copied, env, str(out), "C0")
    assert r.returncode == 0, r.stderr
    hashes = (out / "hashes.txt").read_text()
    assert "C0.mac" in hashes
    assert "C1.mac" not in hashes
    assert "C3.mac" not in hashes


def test_probe_hash_manifest_default_selection(tmp_path):
    """The default 5-case selection produces a manifest with all five macros."""
    copied, env = _probe_env(tmp_path, None)
    out = tmp_path / "hash_default"
    r = _run_probe(copied, env, str(out))
    assert r.returncode == 0, r.stderr
    hashes = (out / "hashes.txt").read_text()
    for c in ("C0", "C1", "C2", "C3", "C4"):
        assert f"{c}.mac" in hashes
    assert "C5.mac" not in hashes


def test_probe_archive_created_outside_tree(tmp_path):
    """The archive must be created in the parent directory, not inside the run dir."""
    copied, env = _probe_env(tmp_path, None)
    out = tmp_path / "arch_out"
    r = _run_probe(copied, env, str(out), "C3")
    assert r.returncode == 0, r.stderr
    # Find the archive in the run directory (moved there after creation).
    archives = list(out.glob("probe_3d_*.tgz"))
    assert len(archives) == 1, f"Expected exactly one archive, found {len(archives)}"
    archive = archives[0]
    # The archive must be readable.
    import tarfile
    with tarfile.open(archive, "r:gz") as tf:
        names = tf.getnames()
    # It must contain the run directory contents.
    assert any("C3.mac" in n for n in names)
    assert any("hashes.txt" in n for n in names)
    # It must NOT contain itself (no self-reference).
    assert not any(archive.name in n for n in names)
    # It must NOT contain .h5 files (excluded).
    assert not any(n.endswith(".h5") for n in names)


def test_probe_archive_relative_path(tmp_path):
    """Archive creation must work with a relative output directory path."""
    copied, env = _probe_env(tmp_path, None)
    # Use a relative path: the run dir is relative to CWD.
    rel_out = "rel_arch_test"
    r = subprocess.run(
        ["bash", str(copied / "run_cases.sh"), rel_out, "C3"],
        env=env, capture_output=True, text=True,
        cwd=str(tmp_path),
    )
    assert r.returncode == 0, r.stderr
    out = tmp_path / rel_out
    assert (out / "C3.mac").exists()
    # Archive should exist in the run directory.
    archives = list(out.glob("probe_3d_*.tgz"))
    assert len(archives) == 1


def test_probe_overwrite_refuses_home_directory(tmp_path):
    """--overwrite must refuse to delete the home directory."""
    copied, env = _probe_env(tmp_path, "#!/bin/sh\nexit 0\n")
    # Create a fake "home" directory and point HOME at it.
    fake_home = tmp_path / "fake_home"
    fake_home.mkdir()
    (fake_home / "important.txt").write_text("do not delete")
    env = dict(env, HOME=str(fake_home))
    r = _run_probe(copied, env, str(fake_home), "--overwrite", "C0")
    assert r.returncode == 3
    assert "Refusing" in r.stderr
    # The data must survive.
    assert (fake_home / "important.txt").exists()


def test_probe_overwrite_refuses_repo_root(tmp_path):
    """--overwrite must refuse to delete the repository root."""
    copied, env = _probe_env(tmp_path, "#!/bin/sh\nexit 0\n")
    # The repo root is the parent of the copied spike dir.
    repo_root = copied.parent
    r = _run_probe(copied, env, str(repo_root), "--overwrite", "C0")
    assert r.returncode == 3
    assert "Refusing" in r.stderr


def test_probe_overwrite_refuses_script_dir(tmp_path):
    """--overwrite must refuse to delete the script/source directory."""
    copied, env = _probe_env(tmp_path, "#!/bin/sh\nexit 0\n")
    r = _run_probe(copied, env, str(copied), "--overwrite", "C0")
    assert r.returncode == 3
    assert "Refusing" in r.stderr


def test_probe_overwrite_refuses_symlink_to_protected(tmp_path):
    """--overwrite must refuse a symlink that resolves into a protected dir."""
    copied, env = _probe_env(tmp_path, "#!/bin/sh\nexit 0\n")
    # Create a symlink pointing to the copied spike dir (protected).
    symlink = tmp_path / "link_to_spike"
    symlink.symlink_to(copied)
    r = _run_probe(copied, env, str(symlink), "--overwrite", "C0")
    assert r.returncode == 3
    assert "Refusing" in r.stderr
    # The script must still exist.
    assert (copied / "run_cases.sh").exists()


def test_probe_overwrite_allows_normal_directory(tmp_path):
    """--overwrite must still work for a normal, non-protected directory."""
    copied, env = _probe_env(tmp_path, None)
    out = tmp_path / "ok_overwrite"
    out.mkdir()
    (out / "stale.txt").write_text("old data")
    r = _run_probe(copied, env, str(out), "--overwrite", "C0")
    assert r.returncode == 0, r.stderr
    assert not (out / "stale.txt").exists()
    assert (out / "C0.mac").exists()
