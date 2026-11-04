"""Probes for two suspected DISCUS 6.23.01 output-menu bugs (for issue #73). Predictions written before running.

P1 (qval.f90:188-253): the <f>^2-type values (f2aver, faver2, faver, ..., and the 'type:sharp' 3D-PDF, which
    divides by faver2) look up the form factor at istl(i,j,1), i.e. at |Q| of the FIRST layer of a 3D grid,
    for every layer. Prediction on a 3D grid h, k, l in [-2, 2]: `value faver2` is the same in every layer l
    (max relative spread along l < 1e-12), while it changes along h by more than 10 %.
P2 (graphic.f90:2861-2889, output_limit): `value 3DPDF, qmax:` tests |h| inside the loop over the three
    components of h and zeroes dsi(i,j,k) with the component index k, so only layers 1-3 can be zeroed.
    Prediction: after `value 3DPDF, qmax:4` the intensity written by a following `value int` has zeros only
    in layers 1-3 of the third axis and is unchanged (equal to the first `value int`) in layers 4 and up.

Usage, from a discus-python checkout with DISCUS installed (Windows: through WSL as in the live tests):

    python tools/probes/fourier_output_probes.py WORKDIR
"""
from __future__ import annotations

import os
import platform
import sys
from pathlib import Path

import numpy as np

from discuspdf import Lattice, ReciprocalGrid
from discuspdf.engine import Engine
from discuspdf.fourier_macro import FourierConfig, StructureSource
from discuspdf.suite_runner import SuiteRunner

ROOT = Path(__file__).resolve().parents[2]
CU = ROOT / "tests" / "fixtures" / "discus_ops" / "cu_b05.cll"
A = 3.615
SRC = StructureSource.cell("cu_b05.cll", (4, 4, 4))
GRID = ReciprocalGrid.from_ranges(Lattice.from_parameters(A, A, A, 90, 90, 90), (-2, 2, 17), (-2, 2, 17), (-2, 2, 17))


def engine() -> Engine:
    if platform.system() == "Windows":
        return Engine(("wsl.exe", "-d", os.environ.get("DISCUS_WSL_DISTRO", "Ubuntu-22.04"), "--",
                       "discus_suite", "-noautorun", "-macro"), timeout_seconds=900)
    return Engine(timeout_seconds=900)


def macro(output_block: list[str]) -> str:
    """The FourierConfig macro with its output block replaced (values FourierConfig does not offer)."""
    # aver 0: with Bragg subtraction a perfect crystal scatters ~0 everywhere on this grid, and P2 could not
    # tell the bug from correct behaviour. four_weight's scalef = 64 * 0.25^3 = 1 here, so dsi3d = dsi.
    text = FourierConfig(GRID, (("int", "x.h5"),)).render(SRC)
    lines = text.splitlines()
    start = lines.index("output")
    end = start + lines[start:].index("exit")
    return "\n".join(lines[:start] + ["output", *output_block, "exit"] + lines[end + 1:]) + "\n"


def main(work: str) -> None:
    runner = SuiteRunner(work, engine())
    # P1
    m1 = macro(["  value int", "  form nexus", "  outf p1_int.h5", "  run",
                "  value faver2", "  form nexus", "  outf p1_faver2.h5", "  run"])
    r1 = runner.run(m1, {"cu_b05.cll": CU}, ["p1_int.h5", "p1_faver2.h5"])
    f = _values(r1.output("p1_faver2.h5"))
    spread_l = float(np.max(np.ptp(f, axis=2) / np.maximum(np.abs(f).max(axis=2), 1e-300)))
    along_h = float(np.ptp(f[:, 8, 8]) / np.abs(f[:, 8, 8]).max())
    print(f"P1: max relative spread of faver2 along the third axis = {spread_l:.3e} (predicted < 1e-12); "
          f"relative change along the first axis = {along_h:.3f} (predicted > 0.1)")
    # P2
    m2 = macro(["  value int", "  form nexus", "  outf p2_before.h5", "  run",
                "  value 3DPDF,qmax:4", "  form nexus", "  outf p2_3dpdf.h5", "  run",
                "  value int", "  form nexus", "  outf p2_after.h5", "  run"])
    r2 = runner.run(m2, {"cu_b05.cll": CU}, ["p2_before.h5", "p2_3dpdf.h5", "p2_after.h5"])
    before, after = _values(r2.output("p2_before.h5")), _values(r2.output("p2_after.h5"))
    zero_new = (after == 0) & (before != 0)
    per_layer = zero_new.sum(axis=(0, 1))
    same_rest = bool(np.array_equal(before[:, :, 3:], after[:, :, 3:]))
    print(f"P2: newly zero points per layer of the third axis = {per_layer.tolist()} "
          f"(predicted: nonzero only in layers 1-3); layers 4+ unchanged = {same_rest} (predicted True)")
    print(f"run directories: {r1.directory} {r2.directory}")


def _values(path):
    import h5py
    with h5py.File(path, "r") as h:
        return np.asarray(h["entry/data/data_values"], dtype=float)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "fourier_output_probes")
