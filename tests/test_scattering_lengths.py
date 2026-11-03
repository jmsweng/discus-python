"""Own neutron scattering lengths (``scat``) for Fourier and powder-PDF runs, and the ``.stru`` source
of :class:`PDFConfig`.

The expected macros are written out by hand here, and the engine-parity tests rebuild the text that the
disorder builder's ``engine.py`` (discus-disorder-builder @ 92fe94e) produced by splicing strings into
the rendered macros, so the move into discuspdf is checked byte for byte.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest

from discuspdf import FourierConfig, Lattice, PDFConfig, PDFRunner, ReciprocalGrid, Seed, StructureSource
from discuspdf.fourier_macro import normalize_scattering_lengths, scat_commands
from discuspdf.pdf import render_pdf_macro

SEED = Seed((11, 22, 33, 44, 55, 66, 77, 88))


def _grid():
    lat = Lattice.from_parameters(3.52, 3.52, 3.52, 90, 90, 90)
    return ReciprocalGrid.centered(lat, (9, 9, 9), (0.25, 0.25, 0.25))


def _fourier(**kw):
    kw.setdefault("radiation", "neutron")
    return FourierConfig(_grid(), outputs=(("int", "diffuse.h5"),), aver_percent=100.0, seed=SEED, **kw)


# --------------------------------------------------------------------------- normalisation
def test_normalize_mapping_and_pairs():
    assert normalize_scattering_lengths({"ni62": 10.3, "C13": 6.19}) == (("C13", 6.19), ("NI62", 10.3))
    assert normalize_scattering_lengths([("Ni62", -8.7)]) == (("NI62", -8.7),)
    assert normalize_scattering_lengths(None) == ()
    assert normalize_scattering_lengths({}) == ()


@pytest.mark.parametrize("bad", [
    {"Ni62": float("nan")}, {"Ni62": "10"}, {"Ni62": True}, {"": 1.0}, {"NI625": 1.0}, {"6LI": 2.0},
    {"VOID": 1.0}, {"Ni62": 1.0, "NI62": 2.0}, [("Ni62",)], {3: 1.0},
])
def test_normalize_rejects(bad):
    with pytest.raises(ValueError):
        normalize_scattering_lengths(bad)


def test_scat_command_format():
    assert scat_commands((("NI62", 10.3),)) == ["  scat NI62, 0,0,0,0,0,0,0,0, 10.3"]
    assert scat_commands((("B11", 6.65e-1), ("LI7", -2.22))) == [
        "  scat B11, 0,0,0,0,0,0,0,0, 0.665", "  scat LI7, 0,0,0,0,0,0,0,0, -2.22"]
    assert scat_commands((("X", 1 / 3),)) == ["  scat X, 0,0,0,0,0,0,0,0, 0.33333333"]


# --------------------------------------------------------------------------- FourierConfig
def test_fourier_default_macro_unchanged_by_the_new_field():
    src = StructureSource.stru("supercell.stru", crystal_cells=(4, 4, 4))
    plain = _fourier().render(src, title="t")
    assert "scat" not in plain
    assert _fourier(scattering_lengths={}).render(src, title="t") == plain
    assert "scattering_lengths" not in _fourier().to_dict()


def test_fourier_scat_lines_before_show():
    cfg = _fourier(scattering_lengths={"Ni62": 10.3, "C13": 6.19})
    lines = cfg.fourier_lines()
    i = lines.index("  show")
    assert lines[i - 2:i] == ["  scat C13, 0,0,0,0,0,0,0,0, 6.19", "  scat NI62, 0,0,0,0,0,0,0,0, 10.3"]
    assert lines[i - 3] == "  set aver,100.0"
    assert lines[i:] == ["  show", "  run", "exit"]


def test_fourier_scat_after_symmetry_line():
    lines = _fourier(scattering_lengths={"Ni62": 10.3}, symmetry="ignore").fourier_lines()
    i = lines.index("  show")
    assert lines[i - 2:i] == ["  set symm:ignore", "  scat NI62, 0,0,0,0,0,0,0,0, 10.3"]


def test_fourier_xray_refuses_scattering_lengths():
    with pytest.raises(ValueError, match="neutron"):
        _fourier(radiation="xray", scattering_lengths={"Ni62": 10.3})


def test_fourier_dict_round_trip():
    cfg = _fourier(scattering_lengths={"Ni62": 10.3})
    d = cfg.to_dict()
    assert d["scattering_lengths"] == {"NI62": 10.3}
    back = FourierConfig.from_dict(json.loads(json.dumps(d)))
    assert back == cfg
    src = StructureSource.stru("supercell.stru", crystal_cells=(4, 4, 4))
    assert back.render(src) == cfg.render(src)


def _engine_fourier_text(cfg_plain: FourierConfig, src, title, extra):
    """engine.py stage_fourier_run (builder 92fe94e): extra lines spliced in before show/run."""
    text = cfg_plain.render(src, title=title)
    marker = "\n  show\n  run\nexit\n"
    assert text.count(marker) == 1
    return text.replace(marker, "\n" + "\n".join(extra) + marker)


@pytest.mark.parametrize("symmetry", [None, "friedel"])
def test_fourier_matches_builder_engine(symmetry):
    custom = {"NI62": 10.3, "C13": 6.19}
    src = StructureSource.stru("supercell.stru", crystal_cells=(4, 4, 4))
    extra = [f"  scat {nm}, 0,0,0,0,0,0,0,0, {custom[nm]:.8g}" for nm in sorted(custom)]   # engine scat_lines
    old = _engine_fourier_text(_fourier(symmetry=symmetry), src, "model diffuse", extra)
    new = _fourier(symmetry=symmetry, scattering_lengths=custom).render(src, title="model diffuse")
    assert new == old


# --------------------------------------------------------------------------- PDFConfig
def test_pdf_cll_macro_unchanged(tmp_path):
    cfg = PDFConfig("model.cll", supercell=(3, 4, 5))
    text = render_pdf_macro(cfg)
    assert "\ncell model.cll,3,4,5\npowder\n" in text
    assert "fourier" not in text and "scat" not in text
    assert cfg.staged_name == "model.cll" and cfg.structure_kind == "cell"


def test_pdf_stru_macro():
    cfg = PDFConfig("supercell.stru", supercell=(1, 1, 1), rmax=10, dr=0.02)
    text = render_pdf_macro(cfg)
    assert text.splitlines()[1:5] == ["discus", "read", "stru model.stru", "powder"]
    assert cfg.staged_name == "model.stru" and cfg.structure_kind == "stru"


def test_pdf_stru_needs_supercell_one():
    with pytest.raises(ValueError, match=r"supercell=\(1, 1, 1\)"):
        PDFConfig("supercell.stru")


@pytest.mark.parametrize("name", ["model.xyz", "model.cif", "model"])
def test_pdf_other_suffixes_refused(name):
    with pytest.raises(ValueError, match=r"\.cll unit cell or a \.stru"):
        PDFConfig(name, supercell=(1, 1, 1))


def test_pdf_scat_block_before_powder():
    cfg = PDFConfig("supercell.stru", supercell=(1, 1, 1), radiation="neutron", qmax=7.0,
                    scattering_lengths={"Ni62": 10.3})
    lines = render_pdf_macro(cfg).splitlines()
    i = lines.index("powder")
    assert lines[i - 4:i + 2] == ["stru model.stru", "fourier", "  scat NI62, 0,0,0,0,0,0,0,0, 10.3",
                                  "exit", "powder", "reset"]


def test_pdf_xray_refuses_scattering_lengths():
    with pytest.raises(ValueError, match="neutron"):
        PDFConfig("supercell.stru", supercell=(1, 1, 1), scattering_lengths={"Ni62": 10.3})


@pytest.mark.parametrize("radiation,custom", [("xray", {}), ("neutron", {}), ("neutron", {"NI62": 10.3, "C13": 6.19})])
def test_pdf_matches_builder_engine(radiation, custom):
    """engine.py run_powder_pdf (builder 92fe94e) with the staged file name model.stru."""
    common = dict(rmax=20.0, dr=0.01, qmax=7.0, qmin=0.5, dq=0.01, radiation=radiation)
    macro = render_pdf_macro(PDFConfig("model.cll", supercell=(1, 1, 1), **common))
    macro = macro.replace("cell model.cll,1,1,1", "stru model.stru")
    extra = [f"  scat {nm}, 0,0,0,0,0,0,0,0, {custom[nm]:.8g}" for nm in sorted(custom)]
    if extra:
        assert macro.count("\npowder\n") == 1
        macro = macro.replace("\npowder\n", "\nfourier\n" + "\n".join(extra) + "\nexit\npowder\n")
    new = render_pdf_macro(PDFConfig("x.stru", supercell=(1, 1, 1), scattering_lengths=custom, **common))
    assert new == macro


# --------------------------------------------------------------------------- PDFRunner with a .stru
def _fake(tmp_path, output="0.0 1.0\n0.1 1.2\n0.2 0.8\n"):
    script = tmp_path / "fake_discus.py"
    script.write_text("from pathlib import Path\nassert Path('model.stru').is_file()\n"
                      "Path('result.pdf').write_text(" + repr(output) + ")\n")
    return PDFRunner(runs_directory=tmp_path / "runs", launcher=[sys.executable, str(script)])


def test_runner_stages_stru_with_lf_and_runs(tmp_path):
    stru = tmp_path / "supercell.stru"
    stru.write_bytes(b"title t\r\ncell 1,1,1,90,90,90\r\natoms\r\nNI62 0,0,0,1\r\n")
    cfg = PDFConfig(stru, supercell=(1, 1, 1), rmax=0.2, dr=0.1, qmax=7.0, radiation="neutron",
                    scattering_lengths={"Ni62": 10.3})
    runner = _fake(tmp_path)
    run = runner.prepare(cfg)
    staged = run.directory / "model.stru"
    assert staged.read_bytes() == stru.read_bytes().replace(b"\r\n", b"\n")
    assert not (run.directory / "model.cll").exists()
    params = json.loads((run.directory / "config.json").read_text())
    assert params["scattering_lengths"] == [["NI62", 10.3]]
    runner.execute(run)
    r, g = run.load()
    assert np.allclose(r, [0, 0.1, 0.2]) and np.allclose(g, [1, 1.2, 0.8])
    assert json.loads((run.directory / "status.json").read_text())["status"] == "complete"


def test_runner_stru_hash_checked(tmp_path):
    stru = tmp_path / "supercell.stru"
    stru.write_text("title t\n")
    runner = _fake(tmp_path)
    run = runner.prepare(PDFConfig(stru, supercell=(1, 1, 1), rmax=0.2, dr=0.1, qmax=7.0))
    (run.directory / "model.stru").write_text("tampered\n")
    with pytest.raises(ValueError, match="hash"):
        runner.execute(run)


def test_runner_cll_config_json_has_no_new_key(tmp_path):
    cell = tmp_path / "cell.cll"
    cell.write_bytes(b"dummy\r\n")
    run = PDFRunner(runs_directory=tmp_path / "runs").prepare(PDFConfig(cell, rmax=0.2, dr=0.1, qmax=7.0))
    params = json.loads((run.directory / "config.json").read_text())
    assert "scattering_lengths" not in params
    assert (run.directory / "model.cll").read_bytes() == b"dummy\r\n"   # .cll staging unchanged
