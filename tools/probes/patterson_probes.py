"""Probes of the DISCUS 6.23.01 ``patterson``, ``inverse`` and ``diff-four`` menus (all run ``patters.f90``) before
wrapping them. Predictions written before running (2026-11-05), from the source (facts gathered by a subagent,
checked against patters.f90 by the overseer for Q2 and Q4).

Input: a random Cu0.75Au0.25 crystal, 2x2x2 fcc cells (32 atoms, so F is complex), X-ray. The hk0 plane h, k = -3..3 in
steps 0.5 (13 x 13) is written by the ``fourier`` + ``output`` menus as NIPL (``form stan``) files: intensity, real
and imaginary part. The real-space plane is x, y = 0..2 cells in 21 x 21 points, z = 0.

Q1 Patterson (``format stan``, ``file a`` = intensity, ``type a, intensity``): rho(r) = vr sum_h I(h) exp(+2 pi i h.r)
   with vr = 1/V_cell (cr_vr, calc_patters patters.f90:1813-1890, the A values times patt_scale * vr). The script
   compares DISCUS's real part with that sum over the file values, within 1e-4 of max |rho|, everywhere except
   Q2's points.
Q2 row-phase bug (CODE-READ, patters.f90:1872-1885): after the inner loops ``iarg = iarg0 + iincu (j-1) + iincv i``
   is not reset when j advances, so for j >= 2 the first point i = 1 uses the phase of the previous row's end. The
   prediction: the first ordinate point of every abscissa column j >= 2 is wrong (more than 1e-2 of max |rho|),
   and every other point is right.
Q3 inverse (``file a`` real, ``file b`` imaginary, ``type a, real``, ``type b, imag``): DISCUS's real part equals
   Re[vr sum F(h) exp(s 2 pi i h.r)] for s = +1 (the code's sign) within 1e-4 of max, outside Q2's points; the script
   also prints the error for s = -1 (the sign that gives the density at +x if F = sum f exp(+2 pi i h.x)).
Q4 diff-four with ``format stan`` (CODE-READ: set_patt_value has no INV_DIFF branch, patters.f90:1905-1992): the map
   is all zero (uncertain: the value may be left from an earlier reflection instead).
Q5 ``sign`` has no effect (CODE-READ, patters.f90:1840 sets patt_sign = 1): Q1 with ``sign -1`` gives the same file.

Usage: python tools/probes/patterson_probes.py WORKDIR
"""
from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from discuspdf.discus_ops import ReplaceRandom, Step, StructureJob
from discuspdf.discus_rmc import FourierPlane, PlaneGrid, read_nipl
from discuspdf.fourier_macro import Seed, StructureSource
from discuspdf.macro import Macro
from discuspdf.suite_runner import SuiteRunner
from tests.live_engine import live_engine

ROOT = Path(__file__).resolve().parents[2]
CU = ROOT / "tests" / "fixtures" / "discus_ops" / "cu_b05.cll"
SRC = StructureSource.cell("cu_b05.cll", (2, 2, 2))
SEED = Seed((11, 22, 33, 44, 55, 66, 77, 88))
A = 3.615
GRID = PlaneGrid(ll=(-3, -3, 0), lr=(3, -3, 0), ul=(-3, 3, 0), nx=13, ny=13)
RHO = dict(ll=(0, 0, 0), lr=(2, 0, 0), ul=(0, 2, 0), na=21, no=21)


@dataclass(frozen=True)
class Lines(Step):
    lines: tuple
    outs: tuple = ()
    ins: tuple = ()

    def input_names(self):
        return self.ins

    def output_names(self):
        return self.outs

    def apply(self, m: Macro) -> None:
        for line in self.lines:
            m.raw(line)


def fourier_value(value, out):
    return ("output", f"value {value}", "form stan", f"outf {out}", "run", "exit")


def menu_lines(menu, files, extra=()):
    g = GRID
    lines = [menu, "reset", "format stan"]
    for ab, (name, kind) in zip("ab", files):
        lines += [f"file {ab}, {name}", f"type {ab}, {kind}"]
    lines += [f"ll {g.ll[0]}, {g.ll[1]}, {g.ll[2]}", f"lr {g.lr[0]}, {g.lr[1]}, {g.lr[2]}",
              f"ul {g.ul[0]}, {g.ul[1]}, {g.ul[2]}", f"na {g.nx}", f"no {g.ny}",
              "rholl {}, {}, {}".format(*RHO["ll"]), "rholr {}, {}, {}".format(*RHO["lr"]),
              "rhoul {}, {}, {}".format(*RHO["ul"]), f"rhona {RHO['na']}", f"rhono {RHO['no']}", *extra, "run", "exit"]
    return lines


def hkl_points():
    g = GRID
    ll, lr, ul = (np.array(v, float) for v in (g.ll, g.lr, g.ul))
    i = np.arange(g.nx)
    j = np.arange(g.ny)
    return ll[None, None] + j[:, None, None] * (ul - ll)[None, None] / (g.ny - 1) \
        + i[None, :, None] * (lr - ll)[None, None] / (g.nx - 1)            # (ny, nx, 3), rows = ordinate


def rho_points():
    ll, lr, ul = (np.array(v, float) for v in (RHO["ll"], RHO["lr"], RHO["ul"]))
    a = np.arange(RHO["na"])
    o = np.arange(RHO["no"])
    return ll[None, None] + o[:, None, None] * (ul - ll)[None, None] / (RHO["no"] - 1) \
        + a[None, :, None] * (lr - ll)[None, None] / (RHO["na"] - 1)       # (no, na, 3)


def direct_sum(values, sign=+1):
    h = hkl_points().reshape(-1, 3)
    v = np.asarray(values).reshape(-1)
    r = rho_points().reshape(-1, 3)
    ph = np.exp(sign * 2j * np.pi * (r @ h.T))
    return (ph @ v).reshape(RHO["no"], RHO["na"]) / A ** 3


def compare(label, got, pred):
    top = np.max(np.abs(pred))
    err = np.abs(got - pred) / top
    bad = np.zeros_like(err, dtype=bool)
    bad[0, 1:] = True       # ordinate point 0 of abscissa columns j >= 2 (file row 0 = first ordinate point)
    print(f"    {label}: max err {err.max():.3e}; outside the Q2 points {err[~bad].max():.3e}; "
          f"at the Q2 points min {err[bad].min():.3e} max {err[bad].max():.3e}; worst at (row, col) "
          f"{np.unravel_index(np.argmax(err), err.shape)}")
    return err


def main(work: str) -> None:
    runner = SuiteRunner(Path(work), live_engine(600))
    cu = {"cu_b05.cll": CU}
    data = StructureJob(SRC, [ReplaceRandom("CU", "AU", 0.25), FourierPlane(GRID, "i.nipl"),
                              Lines(fourier_value("real", "re.nipl") + fourier_value("imag", "im.nipl"),
                                    ("re.nipl", "im.nipl"))], seed=SEED, readback="none").run(runner, cu)
    f = {k: data.run.output(f"{k}.nipl") for k in ("i", "re", "im")}
    zi, _ = read_nipl(f["i"])
    zr, _ = read_nipl(f["re"])
    zm, _ = read_nipl(f["im"])
    print(f"input: I {zi.shape}, max {zi.max():.4g}; |I - re^2 - im^2| max {np.max(np.abs(zi - zr**2 - zm**2)):.3e}")
    inputs = {**cu, "i.nipl": f["i"], "re.nipl": f["re"], "im.nipl": f["im"]}

    def run(label, lines, out):
        try:
            r = StructureJob(SRC, [Lines(tuple(lines), (out,), tuple(n for n in ("i.nipl", "re.nipl", "im.nipl")
                                                                         if any(n in l for l in lines)))],
                             readback="none").run(runner, inputs)
            print(f"{label}: complete, {r.run.directory}")
            return read_nipl(r.run.output(out))[0]
        except Exception as exc:  # noqa: BLE001
            print(f"{label}: FAILED: {exc}")
            return None

    print("== Q1 / Q2 Patterson")
    p = run("Q1 patterson", menu_lines("patterson", [("i.nipl", "intensity")]) + list(fourier_value("real", "p.nipl")),
            "p.nipl")
    if p is not None:
        print(f"    output shape {p.shape}")
        compare("Q1 real part vs vr sum I exp(+2 pi i h.r)", p, direct_sum(zi).real)
    print("== Q3 inverse")
    inv = run("Q3 inverse", menu_lines("inverse", [("re.nipl", "real"), ("im.nipl", "imag")])
              + list(fourier_value("real", "rho.nipl")), "rho.nipl")
    if inv is not None:
        F = zr + 1j * zm
        compare("Q3 sign +1", inv, direct_sum(F, +1).real)
        compare("Q3 sign -1", inv, direct_sum(F, -1).real)
    print("== Q4 diff-four")
    d = run("Q4 diff-four stan", menu_lines("diff-four", [("re.nipl", "real")]) + list(fourier_value("real", "d.nipl")),
            "d.nipl")
    if d is not None:
        print(f"    Q4: max |map| = {np.max(np.abs(d)):.3e} (predicted 0)")
    print("== Q5 sign")
    s = run("Q5 patterson with sign -1", menu_lines("patterson", [("i.nipl", "intensity")], ("sign -1",))
            + list(fourier_value("real", "p.nipl")), "p.nipl")
    if s is not None and p is not None:
        print(f"    Q5: max |diff| / max = {np.max(np.abs(s - p)) / np.max(np.abs(p)):.3e} (predicted 0)")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "patterson_probes")
