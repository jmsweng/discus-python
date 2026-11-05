"""Probes of the DISCUS 6.23.01 ``pdf`` menu fit (``data`` + ``run``) and partial PDFs, before wrapping them.
Predictions written before running (2026-11-05); sections refer to docs/DISCUS_PDFMENU_EXP2PDF_FACTS.md.

Crystal: Cu fcc 3x3x3 cells (108 atoms, B 0.5, ``tests/fixtures/discus_ops/cu_b05.cll``); ordered L1_2 Cu3Au = Au on
site 1; random Cu0.75Au0.25 from the seed.

F1 target: ``calc`` of the ordered crystal (X-ray, periodic, Gaussian peaks from B, qmax 0, ``set range, 8, 0.01``),
   ``save pdf`` -> the first r is 0.01 and the step 0.01, so the file is valid ``data`` input as it stands (3.1, 3.4).
F2 fit: the random crystal, ``data target.dat``, swaps only (``set mode, 1, swchem, all`` + ``set mode, 0, shift,
   all``), sigma 0, 2000 tries, display 200: the log has ``Gen: .. try: .. acc: (good/bad): .. s2x2:`` lines that
   :func:`discuspdf.discus_rmc.parse_rmc_log` reads, with a final block; no uphill acceptances; chi2 never rises.
F3 the two accumulation paths (pitfall 18): ``save pdf`` right after ``run`` equals ``calc`` + ``save pdf`` on the
   fitted crystal within 1e-5 of max |G| (uncertain).
F4 the fit orders the crystal: first-shell Warren-Cowley alpha(Cu-Au) of the result below the start's (printed).
F5 dropped before running: the catalogue already refuses the full word ``select`` in ``discus/pdf`` (the macro
   builder raises), so the wrapper can only write ``sele``.
F6 partials (neutron, periodic, no Gaussian, qmax 0, range 8, 0.001) on the ordered crystal: ``idese all; isele AU;
   jdese all; jsele CU`` and the Au-Au pair. The script fits G = a H_p/(N dr r) + b r to each (H_p the ordered
   pair count of the selected pairs, both orders) and prints a, b and the residual; no prediction of a, b.
F7 selection state (pitfall 15): Au-Au computed after Au-Cu in the same session equals Au-Au in a fresh session.

Round 1 result (2026-11-05, 2d92218): F1 as predicted. F2 failed at ``data``: "Too many points in PDF" (-1). Code
(pdf.f90:944-981): after ``reset all`` PDF_MAXDAT is 1 (:3568), ``data`` allocates exactly n_dat points, and the read loop
tests ``ip > PDF_MAXDAT`` before the read that meets the end of file, so a file of n lines needs room for n + 1.
F6 (neutron weights w_ab = b_a b_b / <b>^2: w_AuAu 0.9829, w_AuCu 0.9943): Au-Au a = 0.98292 = w_AuAu, b/(-4 pi rho0)
= 0.06143 = c_Au^2 w_AuAu; Au-Cu a = 0.49713 = w_AuCu / 2 against the pair count of both orders, b/(-4 pi rho0)
= 0.37285 = 2 c_Au c_Cu w_AuCu. So the cross partial carries half the pairs its baseline assumes (residuals 1e-13).
F7: no state leaks between selections.

Round 2 (``--round 2``, predictions written before running, 2026-11-05):
F2b-F4b as F2-F4, with ``set range, 9, 0.01`` and one ``calc`` before ``data`` (allocates more than 800 + 1 points;
   uncertain). F2b: the log parses, no uphill, chi2 never rises; F3b: run vs calc within 1e-5 of max |G| (uncertain);
   F4b: alpha1 falls.
F8 sum rule (neutron, as F6): G_total - (G_AuAu + G_CuCu + G_AuCu) equals w_AuCu H_AuCu / (2 N dr r) (the missing half
   of the cross pairs, H_AuCu over both orders) within 1e-6 of max |G_total|, and G_CuCu has a = w_CuCu and
   b/(-4 pi rho0) = c_Cu^2 w_CuCu.

Usage: python tools/probes/pdf_fit_probes.py WORKDIR [--round 2]
"""
from __future__ import annotations

import itertools
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from discuspdf.discus_ops import ReplaceRandom, Step, StructureJob
from discuspdf.discus_rmc import parse_rmc_log
from discuspdf.fourier_macro import Seed, StructureSource
from discuspdf.macro import Macro
from discuspdf.suite_runner import SuiteRunner
from tests.live_engine import live_engine

ROOT = Path(__file__).resolve().parents[2]
CU = ROOT / "tests" / "fixtures" / "discus_ops" / "cu_b05.cll"
SRC = StructureSource.cell("cu_b05.cll", (3, 3, 3))
SEED = Seed((11, 22, 33, 44, 55, 66, 77, 88))
A, NC = 3.615, 3
ORDER = ReplaceRandom("CU", "AU", 1.0, site=1)
RANDOM = ReplaceRandom("CU", "AU", 0.25)


@dataclass(frozen=True)
class Lines(Step):
    lines: tuple
    outs: tuple = ()
    ins: tuple = ()
    random: bool = False

    @property
    def uses_random_numbers(self) -> bool:
        return self.random

    def input_names(self):
        return self.ins

    def output_names(self):
        return self.outs

    def apply(self, m: Macro) -> None:
        for line in self.lines:
            m.raw(line)


COMMON_X = ["reset all", "set radiation, xray, 0", "set boundary, periodic, 3D", "set qmax, 0", "set finite, periodic",
            "set density, auto", "set therm, gauss"]
COMMON_N = ["reset all", "set radiation, neutron", "set boundary, periodic, 3D", "set range, 8, 0.001", "set qmax, 0",
            "set finite, periodic", "set density, auto", "set therm, off"]


def read_cols(path):
    rows = [l.split() for l in Path(path).read_text().splitlines() if l.strip() and not l.startswith("#")]
    a = np.array([[float(x) for x in r[:2]] for r in rows])
    return a[:, 0], a[:, 1]


def alpha1(s) -> float:
    xyz = np.asarray(s.cartesian)
    box = NC * A
    d = xyz[:, None, :] - xyz[None, :, :]
    d -= box * np.round(d / box)
    nn = np.abs(np.linalg.norm(d, axis=2) - A / np.sqrt(2)) < 0.2
    sp = np.array([str(x).upper() for x in s.species])
    p = nn[sp == "CU"][:, sp == "AU"].sum() / nn[sp == "CU"].sum()
    return float(1 - p / np.mean(sp == "AU"))


def pair_hist(pair, dr=0.001, rmax=8.0):
    basis = np.array([[0, 0, 0], [0, .5, .5], [.5, 0, .5], [.5, .5, 0]])
    cells = np.array(list(itertools.product(range(NC), repeat=3)))
    xyz = ((cells[:, None, :] + basis[None]).reshape(-1, 3)) * A
    sp = np.where(np.tile(np.arange(4), len(cells)) == 0, "AU", "CU")
    L = NC * A
    H = np.zeros(int(rmax / dr) + 2)
    sel_i = sp == pair[0]
    sel_j = sp == pair[1]
    for sh in itertools.product((-1, 0, 1), repeat=3):
        d = np.linalg.norm(xyz[:, None] - (xyz[None] + np.array(sh) * L), axis=2)
        m = (sel_i[:, None] & sel_j[None]) | (sel_j[:, None] & sel_i[None])
        keep = m & (d > 1e-6) & (d <= rmax)
        np.add.at(H, np.floor((d[keep] + dr / 2) / dr).astype(int), 1)
    return np.arange(len(H)) * dr, H, len(xyz)


def run(runner, label, job, inputs):
    try:
        r = job.run(runner, inputs)
        print(f"{label}: complete, {r.run.directory}")
        return r
    except Exception as exc:  # noqa: BLE001 - a probe reports and goes on
        print(f"{label}: FAILED: {exc}")
        return None


def main(work: str) -> None:
    runner = SuiteRunner(Path(work), live_engine(900))
    cu = {"cu_b05.cll": CU}

    print("== F1")
    t = run(runner, "F1 target", StructureJob(SRC, [ORDER, Lines(("pdf", *COMMON_X, "set range, 8, 0.01", "calc",
                                                                   "save pdf, target.dat", "exit"), ("target.dat",))],
                                               seed=SEED, readback="none"), cu)
    if t is None:
        return
    target = t.run.output("target.dat")
    r0, g0 = read_cols(target)
    print(f"    first r {r0[0]:.4f}, step {np.median(np.diff(r0)):.4f}, {len(r0)} points, max|G| {np.abs(g0).max():.4g}"
          f" (predicted 0.0100, 0.0100)")
    print("    first lines: " + " | ".join(target.read_text().splitlines()[:2]))

    print("== F2 / F3 / F4")
    start = run(runner, "F4 start", StructureJob(SRC, [RANDOM], seed=SEED), cu)
    fit_lines = ("pdf", *COMMON_X, "data target.dat", "dese all", "sele CU, AU", "set mode, 0, shift, all",
                 "set mode, 1, swchem, all", "set scal, off, 1", "set sigm, 0", "set cycl, 2000", "set disp, 200",
                 "run", "save pdf, fit.dat", "calc", "save pdf, recalc.dat", "exit")
    fit = run(runner, "F2 fit", StructureJob(SRC, [RANDOM, Lines(fit_lines, ("fit.dat", "recalc.dat"), ("target.dat",),
                                                                 random=True)], seed=SEED),
              {**cu, "target.dat": target})
    if fit is not None:
        log = parse_rmc_log(fit.run.log.read_text(errors="replace"))
        chi = [p.chi2 for p in log.progress]
        print(f"    F2: {len(log.progress)} blocks, final {log.final}, uphill {[p.accepted_uphill for p in log.progress]}")
        print(f"    F2: chi2 {chi}; never rises: {all(b <= a for a, b in zip(chi, chi[1:]))}")
        rf, gf = read_cols(fit.run.output("fit.dat"))
        rc, gc = read_cols(fit.run.output("recalc.dat"))
        n = min(len(gf), len(gc))
        print(f"    F3: fit.dat {len(gf)} pts r {rf[0]:.3f}..{rf[-1]:.3f}; recalc {len(gc)} pts r {rc[0]:.3f}..{rc[-1]:.3f}; "
              f"max|dG|/max|G| = {np.max(np.abs(gf[:n] - gc[:n])) / np.max(np.abs(gc)):.3e} (predicted < 1e-5)")
        n2 = min(len(gc), len(g0))
        print(f"    fit vs target: max|dG|/max|G| = {np.max(np.abs(gc[:n2] - g0[:n2])) / np.max(np.abs(g0)):.3e}")
        if start is not None:
            print(f"    F4: alpha1 start {alpha1(start.structure):.4f}, fitted {alpha1(fit.structure):.4f}; "
                  f"Au start {sum(str(x).upper() == 'AU' for x in start.structure.species)}, "
                  f"fitted {sum(str(x).upper() == 'AU' for x in fit.structure.species)}")

    print("== F6 / F7")
    sel_aucu = ("idese all", "isele AU", "jdese all", "jsele CU")
    sel_auau = ("idese all", "isele AU", "jdese all", "jsele AU")
    p1 = run(runner, "F6 Au-Cu", StructureJob(SRC, [ORDER, Lines(("pdf", *COMMON_N, *sel_aucu, "calc", "save pdf, p.dat",
                                                                  "exit"), ("p.dat",))], seed=SEED, readback="none"), cu)
    p2 = run(runner, "F6 Au-Au", StructureJob(SRC, [ORDER, Lines(("pdf", *COMMON_N, *sel_auau, "calc", "save pdf, p.dat",
                                                                  "exit"), ("p.dat",))], seed=SEED, readback="none"), cu)
    p3 = run(runner, "F7 Au-Cu then Au-Au", StructureJob(SRC, [ORDER, Lines(
        ("pdf", *COMMON_N, *sel_aucu, "calc", "save pdf, p1.dat", *sel_auau, "calc", "save pdf, p2.dat", "exit"),
        ("p1.dat", "p2.dat"))], seed=SEED, readback="none"), cu)
    for label, res, pair in (("Au-Cu", p1, ("AU", "CU")), ("Au-Au", p2, ("AU", "AU"))):
        if res is None:
            continue
        r, g = read_cols(res.run.output("p.dat"))
        rh, H, n = pair_hist(pair)
        idx = np.rint(r / 0.001).astype(int)
        sel = (r > 0.5) & (r <= 7.9)
        X = np.c_[H[idx[sel]] / (n * 0.001 * r[sel]), r[sel]]
        (a, b), *_ = np.linalg.lstsq(X, g[sel], rcond=None)
        resid = np.max(np.abs(X @ [a, b] - g[sel])) / np.max(np.abs(g[sel]))
        print(f"    F6 {label}: G = a H/(N dr r) + b r with a = {a:.6g}, b = {b:.6g} (b / (-4 pi rho0) = "
              f"{b / (-4 * np.pi * n / (NC * A) ** 3):.6g}); max residual / max|G| = {resid:.2e}")
    if p2 is not None and p3 is not None:
        _, ga = read_cols(p2.run.output("p.dat"))
        _, gb = read_cols(p3.run.output("p2.dat"))
        _, g1 = read_cols(p3.run.output("p1.dat"))
        _, g1f = read_cols(p1.run.output("p.dat")) if p1 is not None else (None, g1)
        print(f"    F7: Au-Au after Au-Cu vs fresh: max|dG|/max|G| = {np.max(np.abs(gb - ga)) / np.max(np.abs(ga)):.3e} "
              f"(predicted 0); Au-Cu in that session vs fresh: {np.max(np.abs(g1 - g1f)) / np.max(np.abs(g1f)):.3e}")


def round2(work: str) -> None:
    runner = SuiteRunner(Path(work), live_engine(900))
    cu = {"cu_b05.cll": CU}
    t = run(runner, "target", StructureJob(SRC, [ORDER, Lines(("pdf", *COMMON_X, "set range, 8, 0.01", "calc",
                                                                "save pdf, target.dat", "exit"), ("target.dat",))],
                                            seed=SEED, readback="none"), cu)
    if t is None:
        return
    target = t.run.output("target.dat")
    _, g0 = read_cols(target)
    start = run(runner, "F4b start", StructureJob(SRC, [RANDOM], seed=SEED), cu)
    fit_lines = ("pdf", *COMMON_X, "set range, 9, 0.01", "calc", "data target.dat", "dese all", "sele CU, AU",
                 "set mode, 0, shift, all", "set mode, 1, swchem, all", "set scal, off, 1", "set sigm, 0",
                 "set cycl, 2000", "set disp, 200", "run", "save pdf, fit.dat", "calc", "save pdf, recalc.dat", "exit")
    fit = run(runner, "F2b fit", StructureJob(SRC, [RANDOM, Lines(fit_lines, ("fit.dat", "recalc.dat"), ("target.dat",),
                                                                  random=True)], seed=SEED),
              {**cu, "target.dat": target})
    if fit is not None:
        text = fit.run.log.read_text(errors="replace")
        log = parse_rmc_log(text)
        chi = [p.chi2 for p in log.progress]
        print(f"    F2b: {len(log.progress)} blocks, final {log.final}, uphill {[p.accepted_uphill for p in log.progress]}")
        print(f"    F2b: chi2 {chi}; never rises: {all(b <= a for a, b in zip(chi, chi[1:]))}; delta {log.delta_chi2}")
        keep = [l.strip() for l in text.splitlines() if any(k in l for k in ("Gen:", "Final", "Delta", "Data", "data"))]
        print("    F2b log lines: " + " || ".join(keep[:30]))
        rf, gf = read_cols(fit.run.output("fit.dat"))
        rc, gc = read_cols(fit.run.output("recalc.dat"))
        n = min(len(gf), len(gc))
        print(f"    F3b: fit.dat {len(gf)} pts r {rf[0]:.3f}..{rf[-1]:.3f}; recalc {len(gc)} pts r {rc[0]:.3f}..{rc[-1]:.3f}; "
              f"max|dG|/max|G| = {np.max(np.abs(gf[:n] - gc[:n])) / np.max(np.abs(gc)):.3e} (predicted < 1e-5)")
        n2 = min(len(gf), len(g0))
        print(f"    fit vs target: max|dG|/max|G| = {np.max(np.abs(gf[:n2] - g0[:n2])) / np.max(np.abs(g0)):.3e}")
        if start is not None:
            print(f"    F4b: alpha1 start {alpha1(start.structure):.4f}, fitted {alpha1(fit.structure):.4f}; "
                  f"Au start {sum(str(x).upper() == 'AU' for x in start.structure.species)}, "
                  f"fitted {sum(str(x).upper() == 'AU' for x in fit.structure.species)}")
    print("== F8")
    sels = {"AuAu": ("idese all", "isele AU", "jdese all", "jsele AU"),
            "CuCu": ("idese all", "isele CU", "jdese all", "jsele CU"),
            "AuCu": ("idese all", "isele AU", "jdese all", "jsele CU"),
            "total": ()}
    lines = ["pdf", *COMMON_N, "calc", "save pdf, total.dat"]  # total first, then the partials
    for k in ("AuAu", "CuCu", "AuCu"):
        lines += [*sels[k], "calc", f"save pdf, {k}.dat"]
    outs = [f"{k}.dat" for k in sels]
    res = run(runner, "F8 partials", StructureJob(SRC, [ORDER, Lines((*lines, "exit"), tuple(outs))], seed=SEED,
                                                  readback="none"), cu)
    if res is None:
        return
    g = {k: read_cols(res.run.output(f"{k}.dat")) for k in sels}
    r = g["total"][0]
    rh, H, n = pair_hist(("AU", "CU"))
    idx = np.rint(r / 0.001).astype(int)
    b = {"CU": 7.718, "AU": 7.63}
    mean = (3 * b["CU"] + b["AU"]) / 4
    w = b["AU"] * b["CU"] / mean ** 2
    sel = (r > 0.5) & (r <= 7.9)
    with np.errstate(divide="ignore", invalid="ignore"):
        missing = w * H[idx] / (2 * n * 0.001 * r)
    diff = g["total"][1] - (g["AuAu"][1] + g["CuCu"][1] + g["AuCu"][1])
    top = np.max(np.abs(g["total"][1][sel]))
    print(f"    F8: max|G_total - sum| / max|G_total| = {np.max(np.abs(diff[sel])) / top:.3e}; after adding the "
          f"missing half: {np.max(np.abs(diff[sel] - missing[sel])) / top:.3e} (predicted < 1e-6); b_Cu, b_Au assumed "
          f"{b['CU']}, {b['AU']}")
    rh, Hc, _ = pair_hist(("CU", "CU"))
    X = np.c_[Hc[idx[sel]] / (n * 0.001 * r[sel]), r[sel]]
    (a, bb), *_ = np.linalg.lstsq(X, g["CuCu"][1][sel], rcond=None)
    print(f"    F8 Cu-Cu: a = {a:.6g} (w_CuCu {b['CU'] ** 2 / mean ** 2:.6g}), b/(-4 pi rho0) = "
          f"{bb / (-4 * np.pi * n / (NC * A) ** 3):.6g} (c_Cu^2 w_CuCu {0.5625 * b['CU'] ** 2 / mean ** 2:.6g})")


if __name__ == "__main__":
    if "--round" in sys.argv and sys.argv[-1] == "2":
        round2(sys.argv[1])
    else:
        main(sys.argv[1] if len(sys.argv) > 1 else "pdf_fit_probes")
