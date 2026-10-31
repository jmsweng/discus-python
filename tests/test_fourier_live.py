"""Live DISCUS 6.23.01 tests for FourierConfig + Diffuse3DRunner (opt-in: DISCUS_LIVE=1).

Each test generates a macro with FourierConfig, stages it with prepare_fourier_run, runs the real
engine through Diffuse3DRunner and compares the output with the committed probe evidence
(benchmarks/discus_6_23_01/b_probes/evidence/20261031T182936Z_LAPTOP-NGU3BT0N). The engine was
deterministic in Gates 2–3, so equal macros and inputs must give byte-identical data values.

Environment:
    DISCUS_LIVE=1             enable
    DISCUS_LIVE_WORK          scratch directory without spaces (default C:\\tmp\\b3live on Windows,
                              /tmp/b3live elsewhere); results.json is written there
    DISCUS_WSL_DISTRO         WSL distribution on Windows (default Ubuntu-22.04)
"""
from __future__ import annotations

import datetime as dt
import json
import os
import platform
import shutil
from pathlib import Path

import numpy as np
import pytest

from discuspdf import Diffuse3DRunner, Lattice, ReciprocalGrid
from discuspdf.fourier_macro import FourierConfig, StructureSource, prepare_fourier_run
from tests.test_fourier_macro import EVIDENCE, probe_configs

pytestmark = pytest.mark.discus_live
if os.environ.get("DISCUS_LIVE") != "1":
    pytest.skip("live DISCUS test; set DISCUS_LIVE=1 to enable", allow_module_level=True)

h5py = pytest.importorskip("h5py")
ROOT = Path(__file__).resolve().parents[1]
INPUTS = {"gate2_random.stru": ROOT / "benchmarks/discus_6_23_01/gate2_prereg/inputs/gate2_random.stru",
          "gate2_block.stru": ROOT / "benchmarks/discus_6_23_01/gate2_prereg/inputs/gate2_block.stru",
          "gaas.cll": ROOT / "tests/fixtures/discus_6_23_01/gaas.cll"}
WORK = Path(os.environ.get("DISCUS_LIVE_WORK", r"C:\tmp\b3live" if platform.system() == "Windows" else "/tmp/b3live"))
RESULTS: dict = {}


def launcher():
    if platform.system() == "Windows":
        # wsl.exe translates the Windows working directory into /mnt/<drive>/... for the command.
        return ("wsl.exe", "-d", os.environ.get("DISCUS_WSL_DISTRO", "Ubuntu-22.04"), "--",
                "discus_suite", "-noautorun", "-macro")
    return ("discus_suite", "-noautorun", "-macro")


@pytest.fixture(scope="module", autouse=True)
def work_dir():
    if " " in str(WORK):
        pytest.fail(f"DISCUS_LIVE_WORK must not contain spaces: {WORK}")
    WORK.mkdir(parents=True, exist_ok=True)
    yield
    RESULTS["finished_utc"] = dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")
    (WORK / "results.json").write_text(json.dumps(RESULTS, indent=1) + "\n", encoding="utf-8", newline="\n")


def run(cfg: FourierConfig, structure: StructureSource, name: str, *, native: bool = True):
    runner = Diffuse3DRunner(runs_directory=WORK / name, output_filename=cfg.outputs[0][1],
                             launcher=launcher(), timeout_seconds=1800)
    job = prepare_fourier_run(runner, cfg, structure, INPUTS[structure.file], title=f"b3 live {name}")
    runner.execute(job, native_nexus=native)
    status = json.loads((job.directory / "status.json").read_text())
    RESULTS.setdefault(name, {})["status"] = status["status"]
    RESULTS[name]["directory"] = str(job.directory)
    return job


def data(path: Path) -> np.ndarray:
    with h5py.File(path, "r") as f:
        return np.asarray(f["/entry/data/data_values"][()], float)


@pytest.mark.parametrize("run_id", ["p2_aver100", "p2_aver50_s8a", "p3_box4_n64_s8a",
                                    "p3_box4_all_y_aver100", "p5_symm_ignore"])
def test_generated_macro_reproduces_probe_evidence(run_id):
    structure, cfg = probe_configs()[run_id]
    job = run(cfg, structure, run_id)
    assert json.loads((job.directory / "status.json").read_text())["status"] == "complete"
    for _, name in cfg.outputs:
        new, old = data(job.directory / name), data(EVIDENCE / run_id / name)
        same = bool(np.array_equal(new, old))
        RESULTS[run_id][name] = {"identical_to_evidence": same,
                                 "max_abs_diff": float(np.max(np.abs(new - old))), "max": float(np.max(np.abs(old)))}
        assert same, f"{name}: differs from evidence by up to {np.max(np.abs(new - old)):.3g}"


def test_oblique_plane_with_generated_labels():
    """Follow-up to p1c_hhl (refused with 'abs h, ord k'): the generated labels are 'abs h, ord l'.
    The l = 0 row of the (h h l) plane must equal the h = k diagonal of the probe's hk0 plane."""
    lat = Lattice.from_parameters(5.6537, 5.6537, 5.6537, 90, 90, 90)
    hhl = ReciprocalGrid.from_corners(lat, (-2, -2, -2), (2, 2, -2), (-2, -2, 2), (-2, -2, -2), (41, 41, 1))
    cfg = FourierConfig(hhl, (("int", "hhl_int.h5"),))
    assert cfg.labels == ("h", "l", "k")
    gaas = StructureSource.cell("gaas.cll", (4, 4, 4))
    # Planes cannot be validated by read_discus_volume (singular increments): run without native check.
    job = run(cfg, gaas, "hhl_plane", native=False)
    with h5py.File(job.directory / "hhl_int.h5", "r") as f:
        g = f["/entry/data"]
        plane = np.asarray(g["data_values"][()], float)
        inc = np.asarray(g["data_increment_vector"][()], float)
        corner = np.asarray(g["data_corner"][()], float)
    assert plane.shape == (41, 41, 1)
    assert np.allclose(inc[:, 0], (0.1, 0.1, 0)) and np.allclose(inc[:, 1], (0, 0, 0.1))
    assert np.allclose(corner, (-2, -2, -2))
    hk0 = data(EVIDENCE / "p1a_hk0" / "p1a_hk0_int.h5")
    diag = hk0[np.arange(41), np.arange(41), 0]
    row = plane[:, 20, 0]                                     # l = -2 + 20*0.1 = 0
    RESULTS["hhl_plane"].update(rel_diff_vs_hk0_diagonal=float(np.linalg.norm(row - diag) / np.linalg.norm(diag)),
                                increment_columns=inc.T.tolist())
    assert np.allclose(row, diag, rtol=1e-10, atol=1e-6 * diag.max())


def test_runner_rejects_engine_error():
    """A lot larger than the crystal: DISCUS reports error -101; the run must be 'failed'."""
    from discuspdf.fourier_macro import Lots
    _, cfg = probe_configs()["p3_off"]
    bad = FourierConfig(cfg.grid, (("int", "big_int.h5"),), lots=Lots("box", (8, 8, 8), "all", True))
    # Lie about the crystal size so the Python check passes and DISCUS has to refuse it.
    structure = StructureSource.stru("gate2_block.stru", (16, 16, 16))
    text = bad.render(structure).replace("lots box,8,8,8", "lots box,17,17,17")
    WORK.joinpath("err").mkdir(exist_ok=True)
    macro = WORK / "err" / "run.mac"
    macro.write_text(text, encoding="ascii", newline="\n")
    staged = WORK / "err" / "gate2_block.stru"
    staged.write_bytes(INPUTS["gate2_block.stru"].read_bytes().replace(b"\r\n", b"\n"))
    runner = Diffuse3DRunner(runs_directory=WORK / "err_runs", output_filename="big_int.h5",
                             launcher=launcher(), timeout_seconds=600)
    job = runner.prepare(macro, assets=(staged,))
    with pytest.raises(RuntimeError):
        runner.execute(job, native_nexus=True)
    status = json.loads((job.directory / "status.json").read_text())
    RESULTS["engine_error"] = {"status": status["status"], "engine_error": status.get("engine_error")}
    assert status["status"] == "failed" and status.get("engine_error")
