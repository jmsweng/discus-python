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

Round 1 result (2026-11-05, 2cee59e): X1 inst-only failed as predicted (message: "Lower fourier limit is less
than lower instrumental limit"); X2 composition parse: ``Cu`` -> Cu 1 and ``Cu3 Au1`` -> Cu 3, Au 1 as predicted,
``Cu:3 Au:1`` -> error -12 (not as predicted: "Cu:3" is read as an optional parameter name), ``Cu3Au`` -> error as
predicted. Every ``run`` then failed with -169 "Fractions add up to 0.000" (phases_average), so X3-X8 gave nothing.

Round 2 (``--round 2``, written before running, 2026-11-05): which set-up lets ``run`` finish. Code-read: the
fraction is set to 1 in powder_run only when powder_run reaches its calculation (powder.f90:376-377); exp2pdf does
not check the error code after powder_run (XR:429-437), so an earlier powder_run error shows up as -169.
Y1 a ``powder`` calculation of a Cu crystal earlier in the same session, then exp2pdf (``composition Cu``);
Y2 a Cu crystal read, then exp2pdf with ``composition comp:current``; Y3 both; Y4 exp2pdf alone with
``radiation neutron``; Y5 exp2pdf alone with ``radiation xray`` (no lambda); Y6 as Y5 without ``reset``.
No prediction of which finish; each line reports complete or the error.

Round 2 result (2026-11-05, 88b0bed): only Y3 finished (G(r) first peak at 2.550 A; Cu a/sqrt2 = 2.556). Y1:
"Unknown element, no Fourier calculated" (-20; the synthetic crystal of ``composition <string>`` is broken, XR:315-355);
Y2: SIGABRT; Y4-Y6: -169 as in round 1. So a run needs a crystal read, a powder calculation in the same session,
then ``composition comp:current``.

Round 3 (``--round 3``, predictions written before running, 2026-11-05), all with the Y3 set-up:
Z1 on the Cu fcc 2x2x2 crystal: X3, X4, X6 (with ``limits inst:[1.0,17.0], fourier:[1.2,16.5]``), X7, X8 as in
   round 1, and X5 (electron 200 vs 300 keV) and a neutron run (completes; first G peak within 0.02 of 2.556).
Z2 composition from a small P1 crystal (cell 20 A): 3 Cu + 1 Au atoms -> the G(r) header's composition line
   gives Cu 0.75, Au 0.25 (or 3 and 1); 1 Cu (occupancy 0.75) + 1 Au (0.25) -> the same G(r) within 1e-6 of
   max |G| (code-read: <f>^2 comes from powder_run on the current crystal, which weights by occupancy; uncertain).
Z3 the powder calculation on the Cu crystal, then ``read`` of the Z2 crystal, then exp2pdf: the same G(r) as Z2
   (3 Cu + 1 Au) within 1e-6 of max |G| (uncertain: whether a new read keeps the phases set-up).

Usage (from a discus-python checkout with DISCUS; Windows through WSL as in the live tests):

    python tools/probes/exp2pdf_probes.py WORKDIR [--round 2]
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


POWDER_BLOCK = ["read", "cell cu_b05.cll,2,2,2", "powder", "reset", "xray", "set wvle,0.5", "set calc,debye",
                "set qmin,0.5", "set qmax,17.9", "set dq,0.01", "set stepr,0.001", "set profile,off",
                "set lpcor,none", "set temp,ignore", "set disp,off", "set period,off", "run", "exit"]
READ_ONLY = ["read", "cell cu_b05.cll,2,2,2"]


def round2(work: str) -> None:
    work = Path(work)
    runner = SuiteRunner(work, live_engine(600))
    step = PowderPattern((("I(Q)", "i.dat"),), qmin=0.5, qmax=17.9, dq=0.01, radiation="xray",
                         wavelength=0.5, stepr=0.001)
    pr = StructureJob(SRC, [step], readback="none").run(runner, {"cu_b05.cll": CU})
    q, i = read_powder_file(pr.run.output("i.dat"))
    files = work / "inputs2"
    files.mkdir(parents=True, exist_ok=True)
    np.savetxt(files / "i.xy", np.c_[q, i])
    inp = {"i.xy": files / "i.xy", "cu_b05.cll": CU}

    def e2p(comp="Cu", rad="radiation xray, lambda:0.1", reset=True):
        return ((["reset"] if reset else []) + ["data xy, i.xy", f"composition {comp}", rad,
                "outfile gr:g.dat, rmin:0.5, rmax:20.0, rstep:0.01, format:stan, fq:f.dat", "run mode:silent"])

    cases = {
        "Y1 powder first": POWDER_BLOCK + ["exp2pdf"] + e2p(),
        "Y2 comp:current": READ_ONLY + ["exp2pdf"] + e2p(comp="comp:current"),
        "Y3 powder first + comp:current": POWDER_BLOCK + ["exp2pdf"] + e2p(comp="comp:current"),
        "Y4 neutron": ["exp2pdf"] + e2p(rad="radiation neutron"),
        "Y5 xray no lambda": ["exp2pdf"] + e2p(rad="radiation xray"),
        "Y6 no reset": ["exp2pdf"] + e2p(rad="radiation xray", reset=False),
    }
    for label, lines in cases.items():
        text = "discus\n" + "\n".join(lines) + "\nexit\nexit\n"
        try:
            r = runner.run(text, inp, ["g.dat", "f.dat"])
            rr, g = read_text_xy(r.output("g.dat"))
            print(f"{label}: complete, {r.directory}; G(r) {len(g)} points, r {rr[0]:.2f}..{rr[-1]:.2f}, "
                  f"max|G| {np.max(np.abs(g)):.4g}, first peak at r = {rr[np.argmax(g)]:.3f}")
        except RuntimeError as exc:
            print(f"{label}: FAILED: {exc}")


def comp_stru(atoms) -> str:
    """A P1 crystal, cell 20 A, atoms (name, occupancy) far apart along x."""
    lines = ["title exp2pdf composition crystal", "spcgr P1",
             "cell    20.000000,  20.000000,  20.000000, 90.000000, 90.000000, 90.000000",
             f"ncell        1,        1,        1,       {len(atoms)}", "format numbers,XYZBPMMOS",
             "atoms      x,              y,              z,             Biso,    Property,  MoleNo,  MoleAt,   Occ,"
             "     St,  Sh,  Sk,  Sl"]
    for k, (name, occ) in enumerate(atoms):
        lines.append(f"{name:<10} {0.2 * k:.6f},       0.000000,       0.000000,      0.000000,       1,       0,"
                     f"       0,   {occ:.6f}, _,   0,   0,   0")
    return "\n".join(lines) + "\n"


def round3(work: str) -> None:
    work = Path(work)
    runner = SuiteRunner(work, live_engine(600))
    step = PowderPattern((("I(Q)", "i.dat"),), qmin=0.5, qmax=17.9, dq=0.01, radiation="xray",
                         wavelength=0.5, stepr=0.001)
    pr = StructureJob(SRC, [step], readback="none").run(runner, {"cu_b05.cll": CU})
    q, i = read_powder_file(pr.run.output("i.dat"))
    files = work / "inputs3"
    files.mkdir(parents=True, exist_ok=True)
    b = 0.05 * i.mean() * (1 + 0.1 * q)
    np.savetxt(files / "i.xy", np.c_[q, i])
    np.savetxt(files / "ib.xy", np.c_[q, i + 2 * b])
    np.savetxt(files / "b.xy", np.c_[q, b])
    np.savetxt(files / "i3.xy", np.c_[q, i, np.sqrt(np.abs(i))], header="Q I dI")
    (files / "c31.stru").write_text(comp_stru([("Cu", 1), ("Cu", 1), ("Cu", 1), ("Au", 1)]))
    (files / "cocc.stru").write_text(comp_stru([("Cu", 0.75), ("Au", 0.25)]))
    inp = {n: files / n for n in ("i.xy", "ib.xy", "b.xy", "i3.xy", "c31.stru", "cocc.stru")}
    inp["cu_b05.cll"] = CU
    powder_cu = POWDER_BLOCK
    powder_on = lambda stru: ["read", f"stru {stru}"] + POWDER_BLOCK[2:]  # noqa: E731

    def go(label, prefix, data="i.xy", back=None, rad="radiation xray, lambda:0.1", extra=(), gr="g.dat",
           fmt=", format:stan", outs=("g.dat", "f.dat")):
        lines = list(prefix) + ["exp2pdf", "reset", f"data xy, {data}"] + ([back] if back else []) + [
            "composition comp:current", rad, *extra,
            f"outfile gr:{gr}, rmin:0.5, rmax:20.0, rstep:0.01{fmt}, fq:f.dat", "run mode:silent"]
        text = "discus\n" + "\n".join(lines) + "\nexit\nexit\n"
        try:
            r = runner.run(text, inp, list(outs))
            print(f"{label}: complete, {r.directory}")
            return r
        except RuntimeError as exc:
            print(f"{label}: FAILED: {exc}")
            return None

    def g_of(r):
        return read_text_xy(r.output("g.dat"))

    def rel(a, b):
        n = min(len(a), len(b))
        return np.max(np.abs(a[:n] - b[:n])) / np.max(np.abs(b))

    print("== Z1")
    base = go("Z1 base", powder_cu)
    lim = go("Z1 X4/X6 limits", powder_cu, extra=("limits inst:[1.0,17.0], fourier:[1.2,16.5]",))
    h5 = go("Z1 X3 default format", powder_cu, gr="g.h5", fmt="", outs=("g.h5", "f.dat"))
    if h5:
        print(f"    X3 first bytes {h5.output('g.h5').read_bytes()[:8]!r} (predicted b'\\x89HDF...')")
    if lim:
        qf, f = read_text_xy(lim.output("f.dat"))
        r, g = g_of(lim)
        dq = np.median(np.diff(qf))
        n0 = int(round(qf[0] / dq))
        qq = np.concatenate([np.arange(1, n0) * dq, qf])
        ff = np.concatenate([np.arange(1, n0) * dq * f[0] / qf[0], f])
        pred = 2 / np.pi * dq * (ff[None] * np.sin(qq[None] * r[:, None])).sum(axis=1)
        print(f"    X4: max|G - pred| / max|G| = {np.max(np.abs(g - pred)) / np.max(np.abs(g)):.3e} (predicted < 1e-2)"
              f"; F file dq {dq:.4g}, {len(qf)} points")
        print(f"    X6: F(Q) file Q range {qf[0]:.4f} .. {qf[-1]:.4f} (predicted 1.2 .. 16.5 within 0.01)")
        head = [l for l in lim.output("g.dat").read_text().splitlines() if l.startswith("#")]
        print("    G(r) header:\n      " + "\n      ".join(head[:40]))
    bg = go("Z1 X7 background x2", powder_cu, data="ib.xy", back="back xy, b.xy, scale:2")
    c3 = go("Z1 X8 three columns", powder_cu, data="i3.xy")
    e2 = go("Z1 X5 electron 200", powder_cu, rad="radiation electron, energy:200")
    e3 = go("Z1 X5 electron 300", powder_cu, rad="radiation electron, energy:300")
    ne = go("Z1 neutron", powder_cu, rad="radiation neutron")
    if base:
        r0, g0 = g_of(base)
        print(f"    base: first peak r = {r0[np.argmax(g0)]:.3f}, max|G| {np.max(np.abs(g0)):.4g}")
        for label, rr in (("X7", bg), ("X8", c3)):
            if rr:
                print(f"    {label}: max|dG|/max|G| = {rel(g_of(rr)[1], g0):.3e} (predicted < 1e-6)")
    if e2 and e3:
        print(f"    X5: max|dG|/max|G| = {rel(g_of(e3)[1], g_of(e2)[1]):.3e} (predicted < 1e-9)")
    if ne:
        rn, gn = g_of(ne)
        print(f"    neutron: first peak r = {rn[np.argmax(gn)]:.3f} (predicted 2.556 within 0.02)")

    print("== Z2 / Z3")
    z31 = go("Z2 3 Cu + 1 Au", powder_on("c31.stru"))
    zoc = go("Z2 occupancies", powder_on("cocc.stru"))
    z3 = go("Z3 powder on Cu, then read 3Cu+1Au", powder_cu + ["read", "stru c31.stru"])
    for label, rr in (("Z2 3+1", z31), ("Z2 occ", zoc), ("Z3", z3)):
        if rr:
            comp = [l for l in rr.output("g.dat").read_text().splitlines() if "composition" in l.lower()]
            print(f"    {label} header: {comp}")
    if z31 and zoc:
        print(f"    Z2: occupancy vs atoms max|dG|/max|G| = {rel(g_of(zoc)[1], g_of(z31)[1]):.3e} (predicted < 1e-6)")
    if z31 and z3:
        print(f"    Z3: max|dG|/max|G| = {rel(g_of(z3)[1], g_of(z31)[1]):.3e} (predicted < 1e-6)")


if __name__ == "__main__":
    if "--round" in sys.argv and sys.argv[-1] == "3":
        round3(sys.argv[1])
        sys.exit()
    args = [a for a in sys.argv[1:] if a != "--round"]
    if "--round" in sys.argv and args[-1] == "2":
        round2(args[0])
    else:
        main(args[0] if args else "exp2pdf_probes")
