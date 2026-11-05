"""Live checks of PattersonMap and InverseMap (opt-in: DISCUS_LIVE=1). Pre-registered (verified level).

Data: random Cu0.75Au0.25 crystals (fcc, from ``tests/fixtures/discus_ops/cu_b05.cll`` and the seed), X-ray; the hk0
plane h, k = -3..3 in steps of 0.5 (13 x 13) written by FourierPlane (intensity, real and imaginary part). Map plane
x, y = 0..2 cells in 21 x 21 points, z = 0.

Registered before any run of these steps (2026-11-05; the probe, tools/probes/patterson_probes.py attempt 2, matched to
3e-5 of max away from the wrong row):

* L1 (Patterson, 3x3x3 cells, 108 atoms): read_map equals (4/108) / V_cell * sum_h I(h) exp(+2 pi i h.r), summed over
  the values of the input file, within 1e-4 of max |rho| at every point of the plane (the first row included: the
  wrong row is the extra one the step asks for). The 4/108 (atoms per cell / atoms) is the normalisation the source
  sets at entry (patters.f90:70-71); the probe used 2x2x2 cells (1/8).
* L2 (inverse, 2x2x2 cells): the real and imaginary outputs equal Re and Im of (4/32) / V_cell * sum_h F(h)
  exp(+2 pi i h.r), F = real + i imag from the input files, within 1e-4 of max |rho|; with exp(-2 pi i h.r) instead the
  real part differs by more than 1e-2 of max |rho| (so the sign is tested).

Environment: DISCUS_LIVE=1; DISCUS_LIVE_WORK (default C:\\tmp\\pattlive or /tmp/pattlive); DISCUS_WSL_DISTRO.
"""
from __future__ import annotations

import os
import platform
from pathlib import Path

import numpy as np
import pytest

from discuspdf.discus_ops import ReplaceRandom, StructureJob
from discuspdf.discus_patterson import InverseMap, PattersonMap, RealPlane, read_map, run_map
from discuspdf.discus_rmc import FourierPlane, PlaneGrid, read_nipl
from discuspdf.fourier_macro import Seed, StructureSource
from tests.live_engine import live_runner

pytestmark = pytest.mark.discus_live
if os.environ.get("DISCUS_LIVE") != "1":
    pytest.skip("live DISCUS test; set DISCUS_LIVE=1 to enable", allow_module_level=True)

ROOT = Path(__file__).resolve().parents[1]
CU = ROOT / "tests" / "fixtures" / "discus_ops" / "cu_b05.cll"
WORK = Path(os.environ.get("DISCUS_LIVE_WORK",
                           r"C:\tmp\pattlive" if platform.system() == "Windows" else "/tmp/pattlive"))
SEED = Seed((11, 22, 33, 44, 55, 66, 77, 88))
A = 3.615
GRID = PlaneGrid(ll=(-3, -3, 0), lr=(3, -3, 0), ul=(-3, 3, 0), nx=13, ny=13)
PLANE = RealPlane((0, 0, 0), (2, 0, 0), (0, 2, 0), 21, 21)


def runner():
    return live_runner(WORK, 600)


def hkl():
    ll, lr, ul = (np.array(v, float) for v in (GRID.ll, GRID.lr, GRID.ul))
    i = np.arange(GRID.nx)[None, :, None]
    j = np.arange(GRID.ny)[:, None, None]
    return (ll + i * (lr - ll) / (GRID.nx - 1) + j * (ul - ll) / (GRID.ny - 1)).reshape(-1, 3)


def direct(values, scale, sign=+1):
    r = PLANE.points().reshape(-1, 3)
    return (np.exp(sign * 2j * np.pi * r @ hkl().T) @ np.asarray(values).reshape(-1)).reshape(PLANE.no, PLANE.na) \
        * scale / A ** 3


def test_l1_patterson_normalisation_and_every_point():
    src = StructureSource.cell("cu_b05.cll", (3, 3, 3))
    step = PattersonMap("i.nipl", GRID, PLANE, "p.nipl")
    data = StructureJob(src, [ReplaceRandom("CU", "AU", 0.25), FourierPlane(GRID, "i.nipl")], seed=SEED,
                        readback="none").run(runner(), {"cu_b05.cll": CU})
    res = run_map(StructureJob(src, [ReplaceRandom("CU", "AU", 0.25), step], seed=SEED, readback="none"), runner(),
                  {"cu_b05.cll": CU, "i.nipl": data.run.output("i.nipl")}, WORK / "l1_inputs")
    z, _ = read_nipl(data.run.output("i.nipl"))
    got = read_map(res.run.output("p.nipl"), PLANE)
    pred = direct(z, 4 / 108).real
    err = np.max(np.abs(got - pred)) / np.max(np.abs(pred))
    assert err < 1e-4, err


def test_l2_inverse_sign_and_imaginary_part():
    src = StructureSource.cell("cu_b05.cll", (2, 2, 2))
    data = StructureJob(src, [ReplaceRandom("CU", "AU", 0.25), FourierPlane(GRID, "re.nipl", value="real"),
                              FourierPlane(GRID, "im.nipl", value="imag")], seed=SEED, readback="none").run(
        runner(), {"cu_b05.cll": CU})
    step = InverseMap("re.nipl", "im.nipl", GRID, PLANE, "rho.nipl", imag_output="rhoi.nipl")
    res = run_map(StructureJob(src, [ReplaceRandom("CU", "AU", 0.25), step], seed=SEED, readback="none"), runner(),
                  {"cu_b05.cll": CU, "re.nipl": data.run.output("re.nipl"), "im.nipl": data.run.output("im.nipl")},
                  WORK / "l2_inputs")
    zr, _ = read_nipl(data.run.output("re.nipl"))
    zi, _ = read_nipl(data.run.output("im.nipl"))
    pred = direct(zr + 1j * zi, 4 / 32)
    top = np.max(np.abs(pred.real))
    got_r = read_map(res.run.output("rho.nipl"), PLANE)
    got_i = read_map(res.run.output("rhoi.nipl"), PLANE)
    assert np.max(np.abs(got_r - pred.real)) / top < 1e-4
    assert np.max(np.abs(got_i - pred.imag)) / top < 1e-4
    other = direct(zr + 1j * zi, 4 / 32, sign=-1).real
    assert np.max(np.abs(got_r - other)) / top > 1e-2
