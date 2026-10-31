"""DISCUS output files (volumes, planes, lines) straight into ScatteringVolume, checked against the
FourierConfig that wrote them. All files are real DISCUS 6.23.01 output from the B probes."""
from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pytest

h5py = pytest.importorskip("h5py")

from discuspdf import Diffuse3DRunner, ScatteringVolume, read_discus_volume
from discuspdf.fourier_macro import FourierConfig, prepare_fourier_run
from tests.test_fourier_macro import EVIDENCE, G32, GAAS, probe_configs

ROOT = Path(__file__).resolve().parents[1]
C3 = ROOT / "tests" / "fixtures" / "discus_6.23.01" / "C3_total.h5"


def ev(run_id, value="int"):
    return EVIDENCE / run_id / f"{run_id}_{value}.h5"


# --------------------------------------------------------------------------- planes and lines
@pytest.mark.parametrize("run_id, shape, ndim, corner", [
    ("p1a_hk0", (41, 41, 1), 2, (-2, -2, 0)),
    ("p1b_h00", (41, 1, 1), 1, (-2, 0, 0)),
    ("p1e_h0l", (41, 1, 41), 2, (-2, 0, -2)),
])
def test_plane_and_line_files(run_id, shape, ndim, corner):
    sv = ScatteringVolume.from_discus_file(ev(run_id))
    assert sv.intensity.shape == shape and sv.grid.ndim == ndim
    assert np.allclose(sv.grid.corner, corner)
    hkl = sv.grid.hkl()
    last = tuple(n - 1 for n in shape)
    expected_last = np.array(corner, float) + 4.0 * (np.array(shape) > 1)
    assert np.allclose(hkl[last], expected_last)            # far corner is +4 r.l.u. along active axes
    assert sv.bragg == "bragg_present" and sv.provenance["content"] == "intensity"
    with h5py.File(ev(run_id), "r") as f:
        assert np.array_equal(sv.intensity, f["/entry/data/data_values"][()])
    with pytest.raises(ValueError, match="singular"):
        read_discus_volume(ev(run_id))                      # the older reader cannot do this


def test_tl_offset_plane_reads_as_the_same_grid():
    a = ScatteringVolume.from_discus_file(ev("p1a_hk0"))
    d = ScatteringVolume.from_discus_file(ev("p1d_tloff"))  # DISCUS wrote (0,0,1) for the 1-point axis
    with h5py.File(ev("p1d_tloff"), "r") as f:
        assert np.allclose(f["/entry/data/data_increment_vector"][()][:, 2], (0, 0, 1))
    assert d.grid == a.grid and np.array_equal(d.intensity, a.intensity)


@pytest.mark.parametrize("name", ["C3_total.h5", "C4_total.h5"])   # C4: skewed increments
def test_volume_matches_the_older_reader(name):
    path = C3.parent / name
    new = ScatteringVolume.from_discus_file(path)
    old = ScatteringVolume.from_discus_volume(read_discus_volume(path))
    assert new.grid.isclose(old.grid, rtol=1e-12, atol=1e-12)
    assert np.array_equal(new.intensity, old.intensity) and new.bragg == old.bragg
    assert new.provenance["sha256"] == old.provenance["sha256"]


def test_real_and_imaginary_parts():
    re_ = ScatteringVolume.from_discus_file(ev("p2_aver100", "real"))
    im_ = ScatteringVolume.from_discus_file(ev("p2_aver100", "imag"))
    assert re_.provenance["content"] == "real part" and im_.provenance["content"] == "imaginary part"
    assert re_.bragg == "bragg_subtracted" and re_.intensity.min() < 0


def test_refuses_real_space_files(tmp_path):
    q = tmp_path / "uvw.h5"
    shutil.copyfile(ev("p2_aver0"), q)
    with h5py.File(q, "r+") as f:
        del f["/entry/data/data_type_axes"]
        f["/entry/data/data_type_axes"] = np.bytes_(b"uvw")
    with pytest.raises(ValueError, match="data_type_axes"):
        ScatteringVolume.from_discus_file(q)


# --------------------------------------------------------------------------- config cross-check
def test_config_is_checked_and_kept():
    structure, cfg = probe_configs()["p2_aver50_s8a"]
    sv = ScatteringVolume.from_discus_file(ev("p2_aver50_s8a"), config=cfg)
    fc = sv.provenance["fourier_config"]
    assert fc["seed"] == [11, 22, 33, 44, 55, 66, 77, 88] and fc["aver_percent"] == 50.0
    real = ScatteringVolume.from_discus_file(ev("p2_aver50_s8a", "real"), config=cfg)
    assert real.provenance["content"] == "real part"


@pytest.mark.parametrize("cfg, match", [
    (FourierConfig(G32, (("int", "p2_aver100_int.h5"),), aver_percent=0), "implies 'bragg_present'"),
    (FourierConfig(G32, (("real", "p2_aver100_int.h5"),), aver_percent=100), "content"),
    (FourierConfig(G32, (("int", "other.h5"),), aver_percent=100), "not one of"),
])
def test_config_mismatches_are_refused(cfg, match):
    with pytest.raises(ValueError, match=match):
        ScatteringVolume.from_discus_file(ev("p2_aver100"), config=cfg)


def test_grid_mismatch_is_refused():
    from discuspdf import ReciprocalGrid, Lattice
    shifted = ReciprocalGrid.from_ranges(Lattice.from_parameters(4, 4, 4, 90, 90, 90),
                                         *[(-1.875, 2.0, 32)] * 3)
    cfg = FourierConfig(shifted, (("int", "p2_aver100_int.h5"),), aver_percent=100)
    with pytest.raises(ValueError, match="corner"):
        ScatteringVolume.from_discus_file(ev("p2_aver100"), config=cfg)


# --------------------------------------------------------------------------- runner
def fake_engine(tmp_path, source: Path, target: str):
    script = tmp_path / "engine.py"
    script.write_text(f"import shutil\nshutil.copyfile({str(source)!r}, {target!r})\n", encoding="utf-8")
    return (sys.executable, str(script))


def test_runner_completes_a_plane_and_checks_its_config(tmp_path):
    _, cfg = probe_configs()["p1a_hk0"]
    runner = Diffuse3DRunner(runs_directory=tmp_path / "runs", output_filename="p1a_hk0_int.h5",
                             launcher=fake_engine(tmp_path, ev("p1a_hk0"), "p1a_hk0_int.h5"))
    job = prepare_fourier_run(runner, cfg, GAAS, ROOT / "tests/fixtures/discus_6_23_01/gaas.cll")
    runner.execute(job, native_nexus=True)
    st = json.loads((job.directory / "status.json").read_text())
    assert st["status"] == "complete"
    rep = st["output_validation"]
    assert rep["grid_dimensionality"] == 2 and rep["fourier_config_checked"] and rep["content"] == "intensity"
    sv = job.load_scattering_volume()
    assert sv.grid.shape == (41, 41, 1) and sv.provenance["fourier_config"]["labels"] == ["h", "k", "l"]


def test_runner_fails_when_output_contradicts_config(tmp_path):
    _, cfg = probe_configs()["p2_aver100"]                  # expects bragg_subtracted
    runner = Diffuse3DRunner(runs_directory=tmp_path / "runs", output_filename="p2_aver100_int.h5",
                             launcher=fake_engine(tmp_path, ev("p2_aver0"), "p2_aver100_int.h5"))
    job = prepare_fourier_run(runner, cfg, probe_configs()["p2_aver100"][0],
                              ROOT / "benchmarks/discus_6_23_01/gate2_prereg/inputs/gate2_random.stru")
    with pytest.raises(RuntimeError, match="Invalid 3D output"):
        runner.execute(job, native_nexus=True)
    st = json.loads((job.directory / "status.json").read_text())
    assert st["status"] == "failed" and "bragg_present" in st["reason"]
