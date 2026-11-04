"""Probes for two suspected DISCUS 6.23.01 `pdf`-menu defects (for issue #73). Predictions written before running.

Q1 (pdf_mod.f90:86/90, fourier_sup.f90:2947-2971): without `set radiation`, pdf_radiation is X-ray (the X-ray
    table is loaded) but pdf_lxray is false, so `form` returns only the constant term c of the table.
    Prediction for L1_2 Cu3Au (3x3x3, periodic, qmax 0): the default G(r) differs from the one with
    `set radiation,xray,0`, and equals within 1e-6 the G(r) computed in Python with weights c_i c_j / <c>^2,
    c_Cu = -3.254477, c_Au = -6.279078 (Waasmaier per_wc); the xray run equals the one with f(0) weights.
Q2 (pdf.f90:552-562, 2723-2733): `set qsig` multiplies G by pdf_exp, which is filled only when `set therm gauss`
    or qalpha > 0. Prediction: with `set qsig,0.05` and `set therm,off`, every G value written is 0 for r > dr.

Usage (from a discus-python checkout with DISCUS; Windows through WSL as in the live tests):

    python tools/probes/pdf_menu_probes.py WORKDIR
"""
from __future__ import annotations

import os
import platform
import sys
from pathlib import Path

import numpy as np

from discuspdf.discus_ops import ReplaceRandom, StructureJob
from discuspdf.discus_pdfmenu import PdfCalc, read_pdf_menu_file
from discuspdf.engine import Engine
from discuspdf.fourier_macro import Seed, StructureSource
from discuspdf.suite_runner import SuiteRunner

ROOT = Path(__file__).resolve().parents[2]
CU = ROOT / "tests" / "fixtures" / "discus_ops" / "cu_b05.cll"
SRC = StructureSource.cell("cu_b05.cll", (3, 3, 3))
SEED = Seed((11, 22, 33, 44, 55, 66, 77, 88))
F0 = {"CU": 28.992065, "AU": 78.967456}
C = {"CU": -3.254477, "AU": -6.279078}


def engine() -> Engine:
    if platform.system() == "Windows":
        return Engine(("wsl.exe", "-d", os.environ.get("DISCUS_WSL_DISTRO", "Ubuntu-22.04"), "--",
                       "discus_suite", "-noautorun", "-macro"), timeout_seconds=600)
    return Engine(timeout_seconds=600)


def macro(drop=(), add_before_calc=()):
    job = StructureJob(SRC, [ReplaceRandom("CU", "AU", 1.0, site=1), PdfCalc("g.dat", rmax=6.0, dr=0.01)],
                       seed=SEED, readback="none")
    lines = [l for l in job.macro().render().splitlines() if not any(l.startswith(d) for d in drop)]
    i = lines.index("calc")
    return "\n".join(lines[:i] + list(add_before_calc) + lines[i:]) + "\n"


def predicted_g(f, rmax=6.0, dr=0.01, a=3.615, nc=3):
    """Periodic weighted pair histogram -> G, as in tests/test_powder_pdfmenu_live.py (Au on site 0)."""
    import itertools
    basis = np.array([[0, 0, 0], [0, .5, .5], [.5, 0, .5], [.5, .5, 0]])
    cells = np.array(list(itertools.product(range(nc), repeat=3)))
    xyz = ((cells[:, None, :] + basis[None]).reshape(-1, 3)) * a
    site = np.tile(np.arange(4), len(cells))
    mean = (f["AU"] + 3 * f["CU"]) / 4
    w = np.where(site == 0, f["AU"], f["CU"]) / mean
    L, n = nc * a, len(xyz)
    H = np.zeros(int(rmax / dr) + 2)
    for sh in itertools.product((-1, 0, 1), repeat=3):
        d = np.linalg.norm(xyz[:, None] - (xyz[None] + np.array(sh) * L), axis=2)
        keep = (d > 1e-6) & (d <= rmax)
        np.add.at(H, np.floor((d[keep] + dr / 2) / dr).astype(int), (w[:, None] * w[None] * np.ones_like(d))[keep])
    r = np.arange(len(H)) * dr
    with np.errstate(divide="ignore", invalid="ignore"):
        return r, H / (n * dr * r) - 4 * np.pi * r * n / L ** 3


def main(work: str) -> None:
    runner = SuiteRunner(work, engine())
    out = {}
    for label, text in (("xray", macro()),
                        ("default", macro(drop=("set radiation",))),
                        ("qsig", macro(add_before_calc=("set qsig,0.05",)))):
        run = runner.run(text, {"cu_b05.cll": CU}, ["g.dat"])
        out[label] = read_pdf_menu_file(run.output("g.dat"))
        print(f"{label}: run directory {run.directory}")
    r, gx = out["xray"]
    _, gd = out["default"]
    sel = (r > 0.5) & (r < 5.9)
    idx = np.rint(r[sel] / 0.01).astype(int)
    rp, g_c = predicted_g(C)
    _, g_f = predicted_g(F0)
    print(f"Q1: max |G(default) - G(c-only weights)| = {np.max(np.abs(gd[sel] - g_c[idx])):.3e} (predicted < 1e-6); "
          f"max |G(xray) - G(f(0) weights)| = {np.max(np.abs(gx[sel] - g_f[idx])):.3e} (predicted < 1e-6); "
          f"max |G(default) - G(xray)| = {np.max(np.abs(gd[sel] - gx[sel])):.3e} (predicted > 0.01)")
    _, gq = out["qsig"]
    print(f"Q2: max |G| with qsig and therm off = {np.max(np.abs(gq[r > 0.011])):.3e} (predicted 0)")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "pdf_menu_probes")
