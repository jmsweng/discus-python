"""Probes of the DISCUS 6.23.01 ``exp2pdf`` menu before wrapping it. Predictions written before running (2026-11-05).

The input is a computed powder pattern: DISCUS ``powder`` I(Q) of Cu fcc 5x5x5 cells (500 atoms, X-ray at 0.5 A, Debye,
Q 0.5..17.9, dq 0.01). Every probe runs ``exp2pdf`` on that pattern (``discus_exp2pdf_*.f90``; XS = _supp,
XL = _load, XR = _run).

X1 limits (XS:445-647): ``limits inst:[1.0,17.0]`` alone fails with "Upper fourier limit is larger than upper
   instrumental limit" (the fourier default [0, 1e9] is still checked). With ``fourier:[1.2,16.5]`` added the run
   completes and writes G(r).
X2 composition (XS:9-232): the log's "ATOM NAME" lines. Predicted: ``Cu`` -> Cu 1; ``Cu3 Au1`` -> Cu 3, Au 1;
   ``Cu:3 Au:1`` -> Cu 1, Au 1 (the value after ':' is overwritten by 1.0; code-read, uncertain);
   ``Cu3Au`` -> an error (one token, '3AU' is not a number; uncertain).
X3 output format (XS:408-415, XD:74): without ``format:`` the G(r) file is HDF5 (starts with \\x89HDF); with
   ``format:stan`` it is text with '#' header lines.
X4 transform (XR:843, PF:7-112): G(r) written equals (2/pi) * sum_k dq F_k sin(Q_k r) over the F(Q) file
   written by the same run, plus the linear ramp from Q = 0 to its first point, within 1 % of max |G|.
X5 electron energy (XR:375, pitfall 22): ``radiation electron, energy:300`` and ``energy:200`` give the same
   G(r) (data are on a Q scale), within 1e-9 of max |G|.
X6 Q range: the F(Q) file runs from the Fourier Qmin to Qmax given in X1 (1.2 .. 16.5) within 0.01.
X7 background (XL:52, XR:166): data I + 2 B with ``back xy, b.xy, scale:2`` gives the same G(r) as data I with
   no background, within 1e-6 of max |G|.
X8 loader: a data file with a '#' first line and three columns (Q, I, dI) loads; the G(r) equals the two-column
   run's within 1e-6 of max |G| (sigma does not enter the unit-weight poly fit, XR:515).

Usage (from a discus-python checkout with DISCUS; Windows through WSL as in the live tests):

    python tools/probes/exp2pdf_probes.py WORKDIR
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np

from discuspdf.discus_ops import StructureJob
from discuspdf.discus_powder import PowderPattern, read_powder_file
from discuspdf.fourier_macro import StructureSource
from discuspdf.suite_runner import SuiteRunner
from tests.live_engine import live_engine

ROOT = Path(__file__).resolve().parents[2]
CU = ROOT / "tests" / "fixtures" / "discus_ops" / "cu_b05.cll"
SRC = StructureSource.cell("cu_b05.cll", (5, 5, 5))


def exp2pdf_macro(lines) -> str:
    return "discus\nexp2pdf\n" + "\n".join(lines) + "\nexit\nexit\n"


def run(runner, label, lines, inputs, outputs):
    try:
        r = runner.run(exp2pdf_macro(lines), inputs, outputs)
        print(f"{label}: complete, {r.directory}")
        return r, r.log.read_text(errors="replace")
    except RuntimeError as exc:
        d = Path(str(exc).rsplit("(see ", 1)[-1].rstrip(")"))
        log = (d / "run.log").read_text(errors="replace") if (d / "run.log").is_file() else ""
        print(f"{label}: FAILED: {exc}")
        for line in log.splitlines():
            if "rror" in line or "limit" in line:
                print(f"    log: {line.strip()}")
        return None, log


def read_text_xy(path):
    rows = [l.split() for l in Path(path).read_text().splitlines() if l.strip() and not l.lstrip().startswith("#")]
    a = np.array([[float(v) for v in row[:2]] for row in rows])
    return a[:, 0], a[:, 1]


def base(data="i.xy", back=None, comp="Cu", rad="radiation xray, lambda:0.1", gr="g.dat", extra=()):
    lines = ["reset", f"data xy, {data}"]
    if back:
        lines.append(back)
    lines += [f"composition {comp}", rad, *extra,
              f"outfile gr:{gr}, rmin:0.5, rmax:20.0, rstep:0.01, format:stan, fq:f.dat", "run mode:silent"]
    return lines


def main(work: str) -> None:
    work = Path(work)
    runner = SuiteRunner(work, live_engine(600))
    step = PowderPattern((("I(Q)", "i.dat"),), qmin=0.5, qmax=17.9, dq=0.01, radiation="xray",
                         wavelength=0.5, stepr=0.001)  # Mo Ka reaches only 17.7
    pr = StructureJob(SRC, [step], readback="none").run(runner, {"cu_b05.cll": CU})
    q, i = read_powder_file(pr.run.output("i.dat"))
    files = work / "inputs"
    files.mkdir(parents=True, exist_ok=True)
    b = 0.05 * i.mean() * (1 + 0.1 * q)
    np.savetxt(files / "i.xy", np.c_[q, i])
    np.savetxt(files / "ib.xy", np.c_[q, i + 2 * b])
    np.savetxt(files / "b.xy", np.c_[q, b])
    np.savetxt(files / "i3.xy", np.c_[q, i, np.sqrt(np.abs(i))], header="Q I dI")
    inp = {n: files / n for n in ("i.xy", "ib.xy", "b.xy", "i3.xy")}
    lim = "limits inst:[1.0,17.0], fourier:[1.2,16.5]"

    print("== X1")
    run(runner, "X1 inst only", base(extra=("limits inst:[1.0,17.0]",)), inp, ["g.dat"])
    r1, _ = run(runner, "X1 inst+fourier", base(extra=(lim,)), inp, ["g.dat", "f.dat"])

    print("== X2")
    for comp in ("Cu", "Cu3 Au1", "Cu:3 Au:1", "Cu3Au"):
        _, log = run(runner, f"X2 {comp!r}", base(comp=comp), inp, ["g.dat"])
        print("    " + " | ".join(l.strip() for l in log.splitlines() if "ATOM NAME" in l or "FOUND" in l))

    print("== X3")
    lines = base()
    lines[-2] = "outfile gr:g.h5, rmin:0.5, rmax:20.0, rstep:0.01"
    r3, _ = run(runner, "X3 default format", lines, inp, ["g.h5"])
    if r3:
        print(f"    first bytes {r3.output('g.h5').read_bytes()[:8]!r} (predicted b'\\x89HDF...')")

    print("== X4 / X6")
    if r1:
        qf, f = read_text_xy(r1.output("f.dat"))
        r, g = read_text_xy(r1.output("g.dat"))
        dq = np.median(np.diff(qf))
        qq = np.concatenate([np.arange(1, int(round(qf[0] / dq))) * dq, qf])
        ff = np.concatenate([np.arange(1, int(round(qf[0] / dq))) * dq * f[0] / qf[0], f])
        pred = 2 / np.pi * dq * (ff[None] * np.sin(qq[None] * r[:, None])).sum(axis=1)
        print(f"X4: max|G - pred| / max|G| = {np.max(np.abs(g - pred)) / np.max(np.abs(g)):.3e} (predicted < 1e-2); "
              f"dq of F file {dq:.4g}")
        print(f"X6: F(Q) file Q range {qf[0]:.4f} .. {qf[-1]:.4f} (predicted 1.2 .. 16.5 within 0.01)")

    print("== X5")
    gs = {}
    for e in (200, 300):
        rr, _ = run(runner, f"X5 electron {e}", base(rad=f"radiation electron, energy:{e}"), inp, ["g.dat"])
        if rr:
            gs[e] = read_text_xy(rr.output("g.dat"))[1]
    if len(gs) == 2:
        print(f"X5: max|G300 - G200| / max|G| = {np.max(np.abs(gs[300] - gs[200])) / np.max(np.abs(gs[200])):.3e} "
              f"(predicted < 1e-9)")

    print("== X7 / X8")
    g0 = run(runner, "X7 no background", base(), inp, ["g.dat"])[0]
    g7 = run(runner, "X7 background x2", base(data="ib.xy", back="back xy, b.xy, scale:2"), inp, ["g.dat"])[0]
    g8 = run(runner, "X8 three columns", base(data="i3.xy"), inp, ["g.dat"])[0]
    if g0:
        ref = read_text_xy(g0.output("g.dat"))[1]
        for label, rr in (("X7", g7), ("X8", g8)):
            if rr:
                other = read_text_xy(rr.output("g.dat"))[1]
                n = min(len(ref), len(other))
                print(f"{label}: max|dG| / max|G| = {np.max(np.abs(other[:n] - ref[:n])) / np.max(np.abs(ref)):.3e} "
                      f"(predicted < 1e-6); lengths {len(ref)}, {len(other)}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "exp2pdf_probes")
