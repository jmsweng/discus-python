"""discuspdf.builder: the disorder-builder back end.

The fixtures under tests/fixtures/builder/ are models exported by the discus-disorder-builder page (@ 92fe94e)
and the SHA-256 of the files its self-contained engine.py 1.0 wrote for them (``--dry-run``). The moved code must
reproduce those files byte for byte.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pytest

from discuspdf.builder import EXPORT_DEFAULTS, ModelError, export_formats, load_model, run
from discuspdf.builder.workflow import cli, main

FIX = Path(__file__).parent / "fixtures" / "builder"
EXPECTED = json.loads((FIX / "engine_1_0_expected.json").read_text())["expected"]


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


@pytest.mark.parametrize("name", sorted(EXPECTED))
def test_structure_files_identical_to_engine_1_0(name, tmp_path):
    res = run(FIX / f"{name}.json", out=tmp_path / "out", dry_run=True, log=lambda m: None)
    out = res["out"]
    for f in ("supercell.stru", "supercell.xyz"):
        assert sha(out / f) == EXPECTED[name][f], f
    summary = json.loads((out / "structure_summary.json").read_text())
    assert summary.pop("engine") == "2.0"
    assert summary.pop("outputs") == EXPORT_DEFAULTS
    assert summary == EXPECTED[name]["summary"]


def test_dry_run_stages_runs_without_executing(tmp_path):
    res = run(FIX / "test_co2_asis.json", out=tmp_path / "o", dry_run=True, log=lambda m: None)
    out = res["out"]
    macros = sorted((out / "runs_3d").glob("*/run.mac"))
    assert macros and all("set technique:turbo" in m.read_text() for m in macros)
    assert "diffuse" not in res and "pdf" not in res
    assert (out / "run_log.txt").read_text().rstrip().endswith(f"Results in {out}")


def test_neutron_model_writes_scat_only_for_present_names(tmp_path):
    cfg = load_model(FIX / "nh4cl_neutron.json")
    cfg["scattering"] = dict(cfg["scattering"], XX9=1.0)           # a name not in the structure
    res = run(cfg, out=tmp_path / "o", dry_run=True, log=lambda m: None)
    texts = [p.read_text() for p in res["out"].rglob("run.mac")]
    assert texts and all("scat" in t for t in texts)
    assert not any("XX9" in t for t in texts)
    assert any("\nfourier\n  scat " in t and "\npowder\n" in t for t in texts)   # 1D PDF run


@pytest.mark.parametrize("bad,match", [
    ({"schema": "something-else"}, "not a disorder-builder model"),
    ({"version": 2}, "newer than this discuspdf"),
])
def test_load_model_refuses(bad, match):
    cfg = load_model(FIX / "cuau_asis.json")
    with pytest.raises(ModelError, match=match):
        load_model({**cfg, **bad})


def test_load_model_needs_keys():
    with pytest.raises(ModelError, match="'sites'"):
        load_model({"cell": {}, "supercell": [1, 1, 1]})


def test_load_model_from_text_and_path():
    p = FIX / "cuau_asis.json"
    assert load_model(p) == load_model(str(p)) == load_model(p.read_text())


def test_export_formats():
    assert export_formats({}) == EXPORT_DEFAULTS
    assert export_formats({"outputs": {"vtk": True, "png": False}}) == {**EXPORT_DEFAULTS, "vtk": True, "png": False}
    with pytest.raises(ModelError, match="unknown export"):
        export_formats({"outputs": {"pdf": True}})


def test_existing_output_refused(tmp_path):
    (tmp_path / "o").mkdir()
    (tmp_path / "o" / "x").write_text("x")
    with pytest.raises(FileExistsError):
        run(FIX / "cuau_asis.json", out=tmp_path / "o", dry_run=True, log=lambda m: None)


def test_cli_usage_and_dry_run(tmp_path, capsys):
    assert cli([]) == 2
    assert cli(["--help"]) == 0
    assert cli(["build", "x.json"]) == 2
    assert cli(["run", str(FIX / "cuau_asis.json"), "--dry-run", "--structure-only", "--out", str(tmp_path / "c")]) == 0
    assert (tmp_path / "c" / "supercell.stru").is_file() and not (tmp_path / "c" / "runs_3d").exists()


def test_main_with_embedded_model(tmp_path):
    cfg = load_model(FIX / "cuau_asis.json")
    assert main(["--dry-run", "--structure-only", "--out", str(tmp_path / "e")], model=cfg) == 0
    assert sha(tmp_path / "e" / "supercell.stru") == EXPECTED["cuau_asis"]["supercell.stru"]


def _pdf_only(cfg: dict, formats: dict) -> dict:
    cfg = json.loads(json.dumps(cfg))
    cfg["calc"]["fourier"]["enabled"] = False
    cfg["mmc"] = {"enabled": False}
    cfg["calc"]["pdf"] = {"enabled": True, "rmax": 0.2, "dr": 0.1, "qmax": 7.0, "qmin": 0.5, "dq": 0.01,
                          "radiation": "xray"}
    cfg["outputs"] = formats
    return cfg


@pytest.mark.parametrize("formats,files,absent", [
    ({}, {"pdf_1d.csv", "pdf_1d.png"}, {"pdf_1d.txt"}),
    ({"txt": True, "csv": False, "png": False}, {"pdf_1d.txt"}, {"pdf_1d.csv", "pdf_1d.png"}),
])
def test_pdf_outputs_follow_formats(tmp_path, formats, files, absent):
    fake = tmp_path / "fake_discus.py"
    fake.write_text("from pathlib import Path\nassert 'stru model.stru' in Path('run.mac').read_text()\n"
                    "Path('result.pdf').write_text('0.0 1.0\\n0.1 1.2\\n0.2 0.8\\n')\n")
    cfg = _pdf_only(load_model(FIX / "cuau_asis.json"), formats)
    res = run(cfg, out=tmp_path / "o", launcher=[sys.executable, str(fake)], log=lambda m: None)
    out = res["out"]
    names = {p.name for p in out.iterdir()}
    assert files <= names and not (absent & names)
    r, g = res["pdf"]
    assert np.allclose(r, [0, 0.1, 0.2]) and np.allclose(g, [1, 1.2, 0.8])
    if "pdf_1d.csv" in names:
        assert (out / "pdf_1d.csv").read_text().splitlines()[0] == "r (Angstrom),G(r)"


@pytest.mark.parametrize("radiation", ["xray", "neutron"])
def test_numpy_preview(tmp_path, radiation):
    cfg = load_model(FIX / "test_co2_asis.json")
    cfg["calc"]["fourier"]["radiation"] = radiation
    cfg["calc"]["preview"] = {"u": [1, 0, 0], "v": [0, 1, 0], "range": 1}
    res = run(cfg, out=tmp_path / "o", structure_only=True, numpy_preview=True, log=lambda m: None)
    d = np.load(res["out"] / "numpy_preview_plane.npz")
    I = d["intensity"]
    assert I.shape == (len(d["t"]), len(d["t"])) and np.isfinite(I).all()
    assert ("neutron" in str(d["factors"])) == (radiation == "neutron")
    # Bragg-subtracted: the diffuse intensity vanishes at integer hkl only if the cells were identical;
    # here it must at least be non-negative within rounding and not identically zero
    assert I.min() > -1e-6 * I.max() and I.max() > 0


def test_anchor_mass_rules():
    from discuspdf.builder.model import anchor_mass
    assert anchor_mass("H") == 1.008 and anchor_mass("D") == 2.014
    assert anchor_mass("C") == 12.0 and anchor_mass("Br") == 70.5      # 2Z, +0.5 above Z = 20
    assert anchor_mass("O18") == 18.0 and anchor_mass("Ni62") == 62.0


def test_site_rotation_cartesian():
    from discuspdf.builder.model import ModelError, site_rotation
    A = np.diag([6.0, 6.0, 6.0])
    R = site_rotation({"label": "S", "sym": "-y,x,z+1/2"}, A)
    assert np.allclose(R, [[0, -1, 0], [1, 0, 0], [0, 0, 1]])
    assert site_rotation({"label": "S"}, A) is None
    hexa = np.array([[4.0, 0, 0], [-2.0, 2.0 * np.sqrt(3), 0], [0, 0, 6.0]])
    R6 = site_rotation({"label": "S", "sym": "x-y,x,z"}, hexa)             # 6-fold about c
    assert np.allclose(R6 @ R6.T, np.eye(3)) and np.isclose(np.trace(R6), 2.0)
    with pytest.raises(ModelError, match="does not fit"):
        site_rotation({"label": "S", "sym": "x-y,x,z"}, A)                   # not a symmetry of a cube
