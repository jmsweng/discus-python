"""Probe: what does stock DISCUS 6.23.01 RMC do without 'set const, p, p'?

rmc_constrain is allocated as 0 for every plane, and rmc_calcskal copies the scale and background of
plane rmc_constrain(p) when it differs from p (rmc_sup_mod.f90:2581-2584), i.e. element 0 of a 1-based
array. This script runs the same short fit (fixed scale 1, fixed background 0, sigma 0) twice, with
and without the 'set const,1,1' line that discuspdf writes, and prints the scale/background lines.
With the line, every block must print scal 1 / back 0. Usage (from a discus-python checkout):

    python tools/probes/rmc_constrain_probe.py WORKDIR
"""
from __future__ import annotations

import os
import platform
import sys
from pathlib import Path

from discuspdf.discus_ops import ReplaceRandom, StructureJob
from discuspdf.discus_rmc import FourierPlane, PlaneGrid, RmcFit, RmcPlane, RmcSwap, parse_rmc_log
from discuspdf.engine import Engine
from discuspdf.fourier_macro import Seed, StructureSource
from discuspdf.suite_runner import SuiteRunner

ROOT = Path(__file__).resolve().parents[2]
CU = ROOT / "tests" / "fixtures" / "discus_ops" / "cu_b05.cll"
SEED = Seed((11, 22, 33, 44, 55, 66, 77, 88))
GRID = PlaneGrid((0, -1, 0), (3, -1, 0), (0, 1, 0), 31, 21)
SRC = StructureSource.cell("cu_b05.cll", (4, 4, 4))


def main(work: str) -> None:
    if platform.system() == "Windows":
        engine = Engine(("wsl.exe", "-d", os.environ.get("DISCUS_WSL_DISTRO", "Ubuntu-22.04"), "--",
                         "discus_suite", "-noautorun", "-macro"), timeout_seconds=900)
    else:
        engine = Engine(timeout_seconds=900)
    runner = SuiteRunner(work, engine)
    target = StructureJob(SRC, [ReplaceRandom("CU", "AU", 1.0, site=1), FourierPlane(GRID, "target.nipl")],
                          seed=SEED, readback="none").run(runner, {"cu_b05.cll": CU}).run.output("target.nipl")
    fit = RmcFit(planes=(RmcPlane("target.nipl", GRID),), moves=(RmcSwap(),), select=("CU", "AU"),
                 tries=500, sigma=0.0, display=100, fit_scale=False, scale=1.0)
    text = StructureJob(SRC, [ReplaceRandom("CU", "AU", 0.25), fit], seed=SEED).macro().render()
    without = "\n".join(line for line in text.splitlines() if line != "set const,1,1") + "\n"
    for label, macro in (("with set const", text), ("without set const", without)):
        run = runner.prepare(macro, {"cu_b05.cll": CU, "target.nipl": target}, ["result.stru"])
        try:
            runner.execute(run)
            status = "complete"
        except RuntimeError as exc:
            status = f"failed: {exc}"
        log = parse_rmc_log(run.log.read_text(errors="replace") if run.log.is_file() else "")
        print(f"== {label}: {status}; run directory {run.directory}")
        for p in log.progress:
            print(f"   gen {p.generated} tries {p.tries} good {p.accepted_downhill} chi2 {p.chi2:.8g} "
                  f"planes {p.planes} final={p.final}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "rmc_constrain_probe")
