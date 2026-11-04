"""Tests for FourierConfig (Package B, Features 8 and 9) and the hardened Diffuse3DRunner.

Independent references:
* the 35 probe macros that DISCUS 6.23.01 actually ran (benchmarks/discus_6_23_01/b_probes/macros,
  written by make_macros.py, a separate generator): rendered configs must equal them byte for byte;
* the committed engine evidence (…/b_probes/evidence/20261031T182936Z_LAPTOP-NGU3BT0N): the engine
  behaviour this API relies on is asserted on the real output files.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np
import pytest

from discuspdf import Diffuse3DRunner, Lattice, ReciprocalGrid
from discuspdf.fourier_macro import (FourierConfig, Lots, Seed, StructureSource, axis_labels,
                                     prepare_fourier_run)

ROOT = Path(__file__).resolve().parents[1]
PROBES = ROOT / "benchmarks" / "discus_6_23_01" / "b_probes"
EVIDENCE = PROBES / "evidence" / "20261031T182936Z_LAPTOP-NGU3BT0N"
LAT = Lattice.from_parameters(4, 4, 4, 90, 90, 90)
SEED_A = Seed((11, 22, 33, 44, 55, 66, 77, 88))
SEED_B = Seed((81, 72, 63, 54, 45, 36, 27, 18))


def cubic(lo, step, n):
    return ReciprocalGrid.from_ranges(LAT, *[(lo, lo + (n - 1) * step, n)] * 3)


G32, G33, GNC = cubic(-2.0, 0.125, 32), cubic(-2.0, 0.125, 33), cubic(-2.0, 0.2, 22)
RANDOM = StructureSource.stru("gate2_random.stru", (16, 16, 16))
BLOCK = StructureSource.stru("gate2_block.stru", (16, 16, 16))
GAAS = StructureSource.cell("gaas.cll", (4, 4, 4))


def outs(run_id, values=("int",)):
    return tuple((v, f"{run_id}_{v}.h5") for v in values)


def probe_configs():
    """Configs that should reproduce the probe macros DISCUS ran (excluding the refused cases)."""
    three = ("int", "real", "imag")
    c = {
        "p1a_hk0": (GAAS, FourierConfig(ReciprocalGrid.from_ranges(LAT, (-2, 2, 41), (-2, 2, 41), (0, 0, 1)),
                                        outs("p1a_hk0"))),
        "p1b_h00": (GAAS, FourierConfig(ReciprocalGrid.from_ranges(LAT, (-2, 2, 41), (0, 0, 1), (0, 0, 1)),
                                        outs("p1b_h00"))),
        "p1e_h0l": (GAAS, FourierConfig(ReciprocalGrid.from_ranges(LAT, (-2, 2, 41), (0, 0, 1), (-2, 2, 41)),
                                        outs("p1e_h0l"))),
        "p2_aver0": (RANDOM, FourierConfig(G32, outs("p2_aver0", three))),
        "p2_aver100": (RANDOM, FourierConfig(G32, outs("p2_aver100", three), aver_percent=100)),
        "p2_aver50_s8a": (RANDOM, FourierConfig(G32, outs("p2_aver50_s8a", three), aver_percent=50, seed=SEED_A)),
        "p2_aver50_s8b": (RANDOM, FourierConfig(G32, outs("p2_aver50_s8b", three), aver_percent=50, seed=SEED_A)),
        "p2_aver50_s8c": (RANDOM, FourierConfig(G32, outs("p2_aver50_s8c", three), aver_percent=50, seed=SEED_B)),
        "p3_off": (BLOCK, FourierConfig(G32, outs("p3_off"))),
        "p3_box4_all_y": (BLOCK, FourierConfig(G32, outs("p3_box4_all_y"), lots=Lots("box", (4, 4, 4), "all", True))),
        "p3_box4_all_n": (BLOCK, FourierConfig(G32, outs("p3_box4_all_n"), lots=Lots("box", (4, 4, 4), "all", False))),
        "p3_ell4_all_y": (BLOCK, FourierConfig(G32, outs("p3_ell4_all_y"),
                                               lots=Lots("ellipsoid", (4, 4, 4), "all", True))),
        "p3_box4_all_y_aver100": (BLOCK, FourierConfig(G32, outs("p3_box4_all_y_aver100"), aver_percent=100,
                                                       lots=Lots("box", (4, 4, 4), "all", True))),
        "p3_box16_n2_y": (BLOCK, FourierConfig(G32, outs("p3_box16_n2_y"), seed=SEED_A,
                                               lots=Lots("box", (16, 16, 16), 2, True))),
        "p5_default": (BLOCK, FourierConfig(G33, outs("p5_default"), aver_percent=100)),
    }
    for tag, seed in (("a", SEED_A), ("b", SEED_A), ("c", SEED_B)):
        rid = f"p3_box4_n64_s8{tag}"
        c[rid] = (BLOCK, FourierConfig(G32, outs(rid), seed=seed, lots=Lots("box", (4, 4, 4), 64, True)))
    for opt in ("ignore", "friedel", "apply"):
        rid = f"p5_symm_{opt}"
        c[rid] = (BLOCK, FourierConfig(G33, outs(rid), aver_percent=100, symmetry=opt))
    for gname, g in (("g32", G32), ("gnc", GNC)):
        for a in (0, 100):
            rid = f"p4_turbo_{gname}_a{a}"
            c[rid] = (BLOCK, FourierConfig(g, outs(rid), aver_percent=a))
    return c


# --------------------------------------------------------------------------- byte equality with probes
@pytest.mark.parametrize("run_id", sorted(probe_configs()))
def test_rendered_macro_equals_the_macro_discus_ran(run_id):
    structure, cfg = probe_configs()[run_id]
    committed = (PROBES / "macros" / f"{run_id}.mac").read_bytes().replace(b"\r\n", b"\n").decode()
    body = committed.split("\n", 1)[1]                       # drop the generator's comment line
    assert cfg.render(structure) == body
    assert cfg.labels_verified


def test_every_successful_turbo_probe_is_covered():
    status = json.loads((EVIDENCE / "runs_status.json").read_text())
    ok = {r["id"] for r in status["runs"] if r["outputs_present"] and "nufft" not in r["id"]
          and r["id"] not in ("p1d_tloff", "p2_aver50_s1a", "p2_aver50_s1b")}
    # p1d (tl != ll, refused by ReciprocalGrid on purpose) and the 1-value seed runs are not offered.
    assert ok == set(probe_configs())


# --------------------------------------------------------------------------- refusals
def test_seed_rules():
    with pytest.raises(ValueError, match="8 values"):
        Seed((1234,))
    with pytest.raises(ValueError, match="1..99999"):
        Seed((0, 1, 2, 3, 4, 5, 6, 7))
    with pytest.raises(ValueError, match="1..99999"):
        Seed((12345678, 1, 2, 3, 4, 5, 6, 7))
    with pytest.raises(ValueError, match="Seed"):
        FourierConfig(G32, outs("x"), aver_percent=50)
    with pytest.raises(ValueError, match="Seed"):
        FourierConfig(G32, outs("x"), lots=Lots("box", (4, 4, 4), 10, True))
    FourierConfig(G32, outs("x"), aver_percent=100)            # deterministic: no seed needed
    FourierConfig(G32, outs("x"), lots=Lots("box", (4, 4, 4), "all", True))


@pytest.mark.parametrize("kwargs, match", [
    (dict(outputs=(("amp", "a.h5"),)), "amp"),
    (dict(outputs=(("int", "a.h5"), ("int", "a.h5"))), "twice"),
    (dict(outputs=(("int", "has space.h5"),)), "plain file name"),
    (dict(outputs=(("int", "a.txt"),)), "end with"),
    (dict(outputs=(("real", "a.h5"),), lots=Lots("box", (4, 4, 4), "all", True)), "only 'int'"),
    (dict(aver_percent=101), r"\[0, 100\]"),
    (dict(temperature="on"), "temperature"),
    (dict(symmetry="mirror"), "symmetry"),
    (dict(wavelength="cuka"), "wavelength"),
    (dict(wavelength=0.0), "wavelength"),
    (dict(wavelength=True), "wavelength"),
    (dict(radiation="gamma"), "radiation"),
    (dict(labels=("h", "h", "l")), "permutation"),
    (dict(outputs=()), "at least one"),
])
def test_invalid_configs(kwargs, match):
    args = dict(grid=G32, outputs=outs("x"))
    args.update(kwargs)
    with pytest.raises(ValueError, match=match):
        FourierConfig(**args)


def test_no_nufft_option():
    import inspect
    assert "technique" not in inspect.signature(FourierConfig).parameters
    assert "set technique:turbo" in FourierConfig(G32, outs("x")).render(BLOCK)


def test_lots_size_against_crystal():
    cfg = FourierConfig(G32, outs("x"), lots=Lots("box", (17, 17, 17), "all", True))
    with pytest.raises(ValueError, match="-101"):
        cfg.render(BLOCK)
    with pytest.raises(ValueError, match="crystal size"):
        cfg.render(StructureSource.stru("gate2_block.stru"))
    assert Lots("box", (4, 4, 4), "all", True).n_lots((16, 16, 16)) == 4096
    assert Lots("box", (4, 4, 4), "all", False).n_lots((16, 16, 16)) == 2197
    with pytest.raises(ValueError):
        Lots("sphere", (4, 4, 4), "all", True)
    with pytest.raises(ValueError):
        Lots("box", (4, 4, 0), "all", True)


def test_lots_counts_match_engine_log():
    for rid, lots in (("p3_box4_all_y", Lots("box", (4, 4, 4), "all", True)),
                      ("p3_box4_all_n", Lots("box", (4, 4, 4), "all", False)),
                      ("p3_box4_n64_s8a", Lots("box", (4, 4, 4), 64, True))):
        log = (EVIDENCE / rid / "run.log").read_text(errors="replace")
        n = int(re.search(r"Fourier volume\s*:\s*(\d+)\s+box shaped lots", log).group(1))
        assert lots.n_lots((16, 16, 16)) == n


def test_structure_source_rules():
    with pytest.raises(ValueError):
        StructureSource.stru("my file.stru")
    with pytest.raises(ValueError):
        StructureSource.stru("a,b.stru")
    with pytest.raises(ValueError):
        StructureSource("cell", "gaas.cll")                 # read cell needs n1,n2,n3
    with pytest.raises(ValueError):
        StructureSource.stru("gaas.cll")


# --------------------------------------------------------------------------- axis labels
def test_oblique_plane_gets_labels_discus_accepts():
    hhl = ReciprocalGrid.from_corners(LAT, (-2, -2, -2), (2, 2, -2), (-2, -2, 2), (-2, -2, -2), (41, 41, 1))
    assert axis_labels(hhl) == ("h", "l", "k")
    cfg = FourierConfig(hhl, outs("hhl"))
    text = cfg.render(GAAS)
    assert "  abs h\n  ord l\n  top k\n" in text
    assert not cfg.labels_verified                          # follow-up live probe
    with pytest.raises(ValueError, match="Component of increment vector is zero"):
        FourierConfig(hhl, outs("hhl"), labels=("h", "k", "l"))  # what p1c_hhl ran, and DISCUS refused


# --------------------------------------------------------------------------- engine facts on the evidence
h5py = pytest.importorskip("h5py")


def values(run_id, v="int"):
    with h5py.File(EVIDENCE / run_id / f"{run_id}_{v}.h5", "r") as f:
        return np.asarray(f["/entry/data/data_values"][()], float)


def bragg_mask():
    h = -2.0 + 0.125 * np.arange(32)
    on = np.abs(h - np.round(h)) < 1e-12
    return on[:, None, None] & on[None, :, None] & on[None, None, :]


def test_engine_subtracts_amplitudes():
    # DISCUS writes value real/imag rounded to single precision (in float64 datasets), so |F|^2 from
    # them matches value int only to float32 rounding: |I - |F|^2| <= 4 eps32 (|Re| + |Im|)^2.
    eps = np.finfo(np.float32).eps
    for rid in ("p2_aver0", "p2_aver100", "p2_aver50_s8a"):
        R, M = values(rid, "real"), values(rid, "imag")
        assert np.all(np.abs(values(rid) - (R * R + M * M)) <= 4 * eps * (np.abs(R) + np.abs(M)) ** 2 + 1e-6)
    b = bragg_mask()
    I0, I100 = values("p2_aver0"), values("p2_aver100")
    assert b.sum() == 64 and np.all(I100[b] == 0) and I0[b].min() > 1e9
    assert np.allclose(I100[~b], I0[~b], rtol=1e-12, atol=1e-6)


def test_engine_seed_reproducibility():
    assert np.array_equal(values("p2_aver50_s8a"), values("p2_aver50_s8b"))
    assert not np.array_equal(values("p2_aver50_s8a"), values("p2_aver50_s8c"))
    assert not np.array_equal(values("p2_aver50_s1a"), values("p2_aver50_s1b"))   # 1-value seed
    assert np.array_equal(values("p3_box4_n64_s8a"), values("p3_box4_n64_s8b"))


def test_engine_lots_scaled_to_whole_crystal():
    off, whole = values("p3_off"), values("p3_box16_n2_y")
    assert np.max(np.abs(whole - off)) <= 1e-15 * off.max()
    b = bragg_mask()
    ratio = values("p3_box4_all_y")[b].sum() / off[b].sum()
    assert ratio == pytest.approx((4 / 16) ** 3, rel=0.01)


# --------------------------------------------------------------------------- runner hardening
def fake(tmp_path, body, output="out.h5"):
    script = tmp_path / "engine.py"
    script.write_text(body, encoding="utf-8")
    macro = tmp_path / "in.mac"
    macro.write_text("# m\n")
    runner = Diffuse3DRunner(runs_directory=tmp_path / "runs", output_filename=output,
                             launcher=(sys.executable, str(script)))
    return runner, runner.prepare(macro)


C3 = ROOT / "tests" / "fixtures" / "discus_6.23.01" / "C3_total.h5"


@pytest.mark.parametrize("line", [
    " ***APPL*** Output value NOT allowed using lots           *** -60 ***",
    " ***MAC *** Error in macro line:    31 Level   1          ***",
    "FINUFFT invokeGuru: setpts error (ier=4)!",
])
def test_runner_fails_on_discus_error_lines(tmp_path, line):
    body = f"import shutil\nshutil.copyfile({str(C3)!r}, 'out.h5')\nprint({line!r})\n"
    runner, job = fake(tmp_path, body)
    with pytest.raises(RuntimeError, match="failed"):
        runner.execute(job, native_nexus=True)
    st = json.loads((job.directory / "status.json").read_text())
    assert st["status"] == "failed" and st["engine_error"]


def test_runner_banner_is_not_an_error(tmp_path):
    banner = ("          ***********************************************************\n"
              "          *         D I S C U S - S U I T E  Version 6.23.01        *\n")
    body = f"import shutil\nshutil.copyfile({str(C3)!r}, 'out.h5')\nprint({banner!r})\n"
    runner, job = fake(tmp_path, body)
    runner.execute(job, native_nexus=True)
    assert json.loads((job.directory / "status.json").read_text())["status"] == "complete"


def test_runner_refuses_all_zero_output(tmp_path):
    body = (f"import shutil, h5py\nshutil.copyfile({str(C3)!r}, 'out.h5')\n"
            "with h5py.File('out.h5', 'r+') as f:\n    f['/entry/data/data_values'][...] = 0\n")
    runner, job = fake(tmp_path, body)
    with pytest.raises(RuntimeError, match="Invalid 3D output"):
        runner.execute(job, native_nexus=True)
    st = json.loads((job.directory / "status.json").read_text())
    assert st["status"] == "failed" and "zero" in st["reason"]


# --------------------------------------------------------------------------- staging and persistence
def test_prepare_fourier_run_stages_macro_config_and_structure(tmp_path):
    cfg = FourierConfig(G32, outs("job", ("int",)) + (("3DPDF", "job_3dpdf.h5"),), aver_percent=100)
    src = ROOT / "benchmarks/discus_6_23_01/gate2_prereg/inputs/gate2_block.stru"
    runner = Diffuse3DRunner(runs_directory=tmp_path / "runs", output_filename="job_int.h5")
    job = prepare_fourier_run(runner, cfg, BLOCK, src, title="block, aver 100")
    staged = (job.directory / "run.mac").read_text()
    assert staged == cfg.render(BLOCK, title="block, aver 100")
    hashes = json.loads((job.directory / "inputs.json").read_text())
    assert set(hashes) == {"run.mac", "gate2_block.stru", "fourier_config.json"}
    saved = json.loads((job.directory / "fourier_config.json").read_text())
    assert saved["bragg"] == "bragg_subtracted" and saved["technique"] == "turbo"
    assert FourierConfig.from_dict(saved).render(BLOCK) == cfg.render(BLOCK)
    with pytest.raises(ValueError, match="output_filename"):
        prepare_fourier_run(Diffuse3DRunner(runs_directory=tmp_path / "r2", output_filename="other.h5"),
                            cfg, BLOCK, src)
    with pytest.raises(ValueError, match="does not match"):
        prepare_fourier_run(runner, cfg, RANDOM, src)
    # CRLF input (Windows checkout) is staged with LF; the bytes DISCUS reads equal the git blob.
    crlf = tmp_path / "crlf" / "gate2_block.stru"
    crlf.parent.mkdir()
    crlf.write_bytes(src.read_bytes().replace(b"\n", b"\r\n"))
    job2 = prepare_fourier_run(runner, cfg, BLOCK, crlf)
    assert (job2.directory / "gate2_block.stru").read_bytes() == src.read_bytes()
    assert json.loads((job2.directory / "fourier_config.json").read_text())["structure"]["staged_with_lf"]


def test_round_trip_and_bragg_label():
    cfg = FourierConfig(G32, outs("x"), aver_percent=50, seed=SEED_A, symmetry="ignore",
                        lots=Lots("ellipsoid", (4, 5, 6), 9, False))
    again = FourierConfig.from_dict(json.loads(json.dumps(cfg.to_dict())))
    assert again == cfg
    assert cfg.bragg == "bragg_subtracted" and FourierConfig(G32, outs("x")).bragg == "bragg_present"
    with pytest.raises(ValueError):
        FourierConfig.from_dict({**cfg.to_dict(), "version": 3})


# --------------------------------------------------------------------------- phase 2c-A options
def _lines(cfg):
    text = cfg.render(BLOCK)
    from discuspdf.suite_catalogue import load_catalogue
    assert load_catalogue().check(text) == []
    return [l.strip() for l in text.splitlines()]


def test_electron_energy_dispersion_delf_and_3dpdf_limit_lines():
    from discuspdf.engine import fork_reasons
    e = _lines(FourierConfig(G32, outs("e", ("int",)), radiation="electron", energy=200.0))
    assert e[e.index("fourier") + 1:e.index("fourier") + 3] == ["electron", "energy 200"]
    assert not any(l.startswith("wvle") for l in e)
    x = _lines(FourierConfig(G32, outs("x", ("int",)), dispersion=True))
    assert "disp anom" in x and "wvle moa1" in x
    d = _lines(FourierConfig(G32, outs("d", ("int",)), wavelength=0.7093, dispersion=True,
                             delf={"Cu": (0.26, 1.6), "au": (-2.0, 8.5)}))
    assert d.index("delf CU,0.26,1.6") < d.index("show") and "delf AU,-2,8.5" in d
    assert fork_reasons("\n".join(d)) == [] and fork_reasons("\n".join(e)) == []


@pytest.mark.parametrize("kwargs, match", [
    (dict(energy=0.0), "energy"),
    (dict(energy=17.0, wavelength="cua1"), "either"),
    (dict(dispersion=True, radiation="neutron"), "X-rays only"),
    (dict(dispersion=True, wavelength=0.71), "K-alpha1"),
    (dict(dispersion=True, wavelength="moa2"), "K-alpha1"),
    (dict(dispersion=True, energy=17.4), "K-alpha1"),
    (dict(delf={"CU": (0.1, 0.2)}), "dispersion=True"),
    (dict(radiation="electron"), "electrons need"),
    (dict(wavelength=None), "needs energy"),
    (dict(dispersion=True, delf={"CU": (0.1,)}), "f'"),
    (dict(dispersion=True, delf={"CU": (0, 1), "cu": (0, 1)}), "twice"),
])
def test_phase2c_refusals(kwargs, match):
    with pytest.raises((ValueError, TypeError), match=match):
        FourierConfig(G32, outs("r", ("int",)), **kwargs)


def test_3dpdf_must_be_the_single_last_output():
    with pytest.raises(ValueError, match="last"):
        FourierConfig(G32, (("3DPDF", "a.h5"), ("int", "b.h5")), aver_percent=100)
    with pytest.raises(ValueError, match="last"):
        FourierConfig(G32, (("3DPDF", "a.h5"), ("3DPDF", "b.h5")), aver_percent=100)


def test_phase2c_round_trip():
    cfg = FourierConfig(G32, outs("j", ("int",)) + (("3DPDF", "j_3dpdf.h5"),), aver_percent=100, wavelength=0.7093,
                        dispersion=True, delf={"CU": (0.26, 1.6)})
    back = FourierConfig.from_dict(json.loads(json.dumps(cfg.to_dict())))
    assert back.render(BLOCK) == cfg.render(BLOCK)
    e = FourierConfig(G32, outs("e", ("int",)), radiation="electron", energy=300.0)
    assert FourierConfig.from_dict(json.loads(json.dumps(e.to_dict()))).render(BLOCK) == e.render(BLOCK)
    plain = FourierConfig(G32, outs("q", ("int",)))
    assert not {"energy", "dispersion", "delf"} & set(plain.to_dict())


def test_energy_runs_record_no_wavelength():
    d = FourierConfig(G32, outs("e", ("int",)), radiation="electron", energy=200.0).to_dict()
    assert d["wavelength"] is None and d["energy"] == 200.0
    assert FourierConfig(G32, outs("n", ("int",)), radiation="electron", wavelength=0.0251).render(BLOCK)
