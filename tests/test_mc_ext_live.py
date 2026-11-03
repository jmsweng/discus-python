"""Live checks of the phase 2b MMC extensions (opt-in: DISCUS_LIVE=1). Engine-ran level, small models.

* E1: ``Convergence(diff=0.02)`` on one occupational feedback target stops before the full number of
  moves, and the final achieved correlation is within 0.02 of the target.
* E2: ``Shift(..., directions={"Cu": (1, 0, 0)})`` moves atoms along x only.
* E3: repulsive and Buckingham targets are accepted by DISCUS and their rows are read from the log.

Structures: the committed Package D probe crystals (``benchmarks/discus_6_23_01/d_mc_probes``).
Environment as in test_mc_live.py (DISCUS_LIVE_WORK default C:\\tmp\\mcextlive or /tmp/mcextlive).
"""
from __future__ import annotations

import os
import platform
from pathlib import Path

import numpy as np
import pytest

from discuspdf.discus_io import read_stru
from discuspdf.fourier_macro import Seed
from discuspdf.mc import (BuckinghamTarget, Convergence, MCConfig, MCRunner, Neighborhood, OccupationalTarget,
                          RepulsiveTarget, Shift, SpringTarget, Swap, Vector)

pytestmark = pytest.mark.discus_live
if os.environ.get("DISCUS_LIVE") != "1":
    pytest.skip("live DISCUS test; set DISCUS_LIVE=1 to enable", allow_module_level=True)

ROOT = Path(__file__).resolve().parents[1]
EV = ROOT / "benchmarks/discus_6_23_01/d_mc_probes/evidence/run1"
WORK = Path(os.environ.get("DISCUS_LIVE_WORK",
                           r"C:\tmp\mcextlive" if platform.system() == "Windows" else "/tmp/mcextlive"))
SEED = Seed((11, 22, 33, 44, 55, 66, 77, 88))
SIX = Neighborhood(tuple(Vector(1, 1, t) for t in [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]))


def launcher():
    if platform.system() == "Windows":
        return ("wsl.exe", "-d", os.environ.get("DISCUS_WSL_DISTRO", "Ubuntu-22.04"), "--",
                "discus_suite", "-noautorun", "-macro")
    return ("discus_suite", "-noautorun", "-macro")


def run(name, cfg, stru, title):
    if " " in str(WORK):
        pytest.fail(f"DISCUS_LIVE_WORK must not contain spaces: {WORK}")
    structure = read_stru(EV / "work" / stru).structure
    runner = MCRunner(WORK / name, launcher=launcher(), timeout_seconds=900)
    return structure, runner.execute(runner.prepare(cfg, structure, title=title))


def test_e1_convergence_stops_early():
    moves = 30
    cfg = MCConfig((SIX,), (OccupationalTarget(0, "Cu", "Au", correlation=-0.3),), (Swap(),), SEED,
                   moves_per_atom=moves, feedback_per_atom=1, finish=Convergence(diff=0.02))
    _, res = run("e1_converge", cfg, "cuau.stru", "D MMC probe Cu0.5Au0.5 sc 16^3")
    n_feed = len(res.log["feedback"])
    assert n_feed < moves, n_feed
    rep = res.target_report(tol=0.02)[0]
    assert rep["reached"], rep


def test_e2_directional_shift_moves_along_x_only():
    cfg = MCConfig((SIX,), (SpringTarget(0, "Cu", "Cu", 3.0, 1.0),),
                   (Shift({"Cu": 0.01}, directions={"Cu": (1, 0, 0)}),), SEED, moves_per_atom=5, feedback_per_atom=1)
    before, res = run("e2_direction", cfg, "cu_disp.stru", "D MMC probe Cu sc 12^3 displaced")
    d = np.asarray(res.structure.frac) - np.asarray(before.frac)
    d -= np.round(d)
    assert np.max(np.abs(d[:, 1:])) < 1e-6, np.max(np.abs(d[:, 1:]))
    assert np.max(np.abs(d[:, 0])) > 1e-4


def test_e3_repulsive_and_buckingham_rows():
    cfg = MCConfig((SIX,), (RepulsiveTarget(0, "Cu", "Cu", 1.0, b_scale=2.5, m=12.0),
                           BuckinghamTarget(0, "Cu", "Cu", 1000.0, 0.3, 10.0)),
                   (Shift({"Cu": (0.005, 0.005, 0.005)}),), SEED, moves_per_atom=3, feedback_per_atom=1)
    _, res = run("e3_rep_buck", cfg, "cu_disp.stru", "D MMC probe Cu sc 12^3 displaced")
    kinds = {r["kind"] for r in res.log["final"]["rows"]}
    assert {"repulsive", "buckingham"} <= kinds, kinds
