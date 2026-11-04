"""Live checks of PowderPattern and PdfCalc (opt-in: DISCUS_LIVE=1). Pre-registered (verified level).

Model: Cu fcc 3x3x3 cells (108 atoms) from ``tests/fixtures/discus_ops/cu_b05.cll``; for the two-species check
the same crystal with Au on site 1 (L1_2, 27 Au + 81 Cu).

Registered before any run (2026-11-03):

* D1 (powder, neutron, Debye): the S(Q) DISCUS writes at label L equals, for every label but the first,
  S(L - dq) = 1 + (2/N) sum_{i<j} sinetab(r~ (L - dq)) / ((L - dq) r~) with r~ = round(r/stepr) stepr and
  DISCUS's 65536-entry sine table (floor of the argument), within 1e-5; and the textbook sinc version within 2e-3.
* D2: F(Q) at label L equals L (S(L - dq) - 1) from the same emulation, within 1e-4.
* G1 (pdf menu, neutron, periodic, qmax 0): G(r_k) = H_k / (N dr r_k) - 4 pi r_k rho0 with H_k the ordered pair
  count over periodic images in bin k = int((d + dr/2)/dr), rho0 = N/V, within 1e-6 for r_k <= rmax - 0.1.
* G2 (pdf menu, X-ray, Cu3Au): the same with weights f_i f_j / <f>^2, f = Waasmaier f(0) = c + sum a_i of
  lib_f90/element_data_mod.f90, within 1e-6.

Environment: DISCUS_LIVE=1; DISCUS_LIVE_WORK (default C:\\tmp\\powpdflive or /tmp/powpdflive);
DISCUS_WSL_DISTRO (default Ubuntu-22.04).
"""
from __future__ import annotations

import itertools
import math
import os
import platform
from pathlib import Path

import numpy as np
import pytest

from discuspdf.discus_ops import ReplaceRandom, StructureJob
from discuspdf.discus_pdfmenu import PdfCalc, read_pdf_menu_file
from discuspdf.discus_powder import PowderPattern, read_powder_file
from discuspdf.engine import Engine
from discuspdf.fourier_macro import Seed, StructureSource
from discuspdf.suite_runner import SuiteRunner

pytestmark = pytest.mark.discus_live
if os.environ.get("DISCUS_LIVE") != "1":
    pytest.skip("live DISCUS test; set DISCUS_LIVE=1 to enable", allow_module_level=True)

ROOT = Path(__file__).resolve().parents[1]
CU = ROOT / "tests" / "fixtures" / "discus_ops" / "cu_b05.cll"
WORK = Path(os.environ.get("DISCUS_LIVE_WORK",
                           r"C:\tmp\powpdflive" if platform.system() == "Windows" else "/tmp/powpdflive"))
A, NC = 3.615, 3
SRC = StructureSource.cell("cu_b05.cll", (NC, NC, NC))
SEED = Seed((11, 22, 33, 44, 55, 66, 77, 88))
# Waasmaier f(0) = c + sum a_i (lib_f90/element_data_mod.f90, per_wc + per_wa1..5)
F0 = {"CU": -3.254477 + 14.014192 + 4.784577 + 5.056806 + 1.457971 + 6.932996,    # 28.992065
      "AU": -6.279078 + 16.777389 + 19.317156 + 32.979682 + 5.595453 + 10.576854}  # 78.967456


def runner() -> SuiteRunner:
    if " " in str(WORK):
        pytest.skip("DISCUS_LIVE_WORK must not contain spaces")
    if platform.system() == "Windows":
        engine = Engine(("wsl.exe", "-d", os.environ.get("DISCUS_WSL_DISTRO", "Ubuntu-22.04"), "--",
                         "discus_suite", "-noautorun", "-macro"), timeout_seconds=600)
    else:
        engine = Engine(timeout_seconds=600)
    return SuiteRunner(WORK, engine)


def fcc_positions():
    basis = np.array([[0, 0, 0], [0, .5, .5], [.5, 0, .5], [.5, .5, 0]])
    cells = np.array(list(itertools.product(range(NC), repeat=3)))
    return ((cells[:, None, :] + basis[None]) .reshape(-1, 3)) * A, np.tile(np.arange(4), len(cells))


# ----------------------------------------------------------------------------------------------- powder
DQ, QMIN, QMAX, STEPR = 0.01, 0.5, 15.0, 0.01


def emulated_s(q, exact=False):
    xyz, _ = fcc_positions()
    d = np.linalg.norm(xyz[:, None] - xyz[None], axis=2)[np.triu_indices(len(xyz), 1)]
    j = np.floor((d + 0.5 * STEPR) / STEPR).astype(int)
    bins, counts = np.unique(j, return_counts=True)
    rt = bins * STEPR
    x = q[:, None] / (2 * np.pi)
    arg = rt[None] * x
    if exact:
        s = np.sin(2 * np.pi * arg)
    else:
        idx = np.floor(arg * 65536).astype(np.int64) & 65535
        s = np.sin(2 * np.pi * idx / 65536)
    return 1 + 2.0 / len(xyz) * (counts[None] * s / (2 * np.pi * arg)).sum(axis=1)


@pytest.fixture(scope="module")
def powder():
    step = PowderPattern((("S(Q)", "s.dat"), ("F(Q)", "f.dat")), qmin=QMIN, qmax=QMAX, dq=DQ,
                         radiation="neutron", wavelength=0.5, stepr=STEPR)
    r = StructureJob(SRC, [step], readback="none").run(runner(), {"cu_b05.cll": CU})
    return read_powder_file(r.run.output("s.dat")), read_powder_file(r.run.output("f.dat"))


def test_d1_debye_s_of_q(powder):
    (q, s), _ = powder
    assert np.allclose(np.diff(q), DQ, atol=1e-6) and abs(q[0] - QMIN) < 1e-6
    q, s = q[1:], s[1:]
    emul = emulated_s(q - DQ)
    assert np.max(np.abs(s - emul)) < 1e-5, np.max(np.abs(s - emul))
    assert np.max(np.abs(s - emulated_s(q - DQ, exact=True))) < 2e-3


def test_d2_f_of_q(powder):
    _, (q, f) = powder
    q, f = q[1:], f[1:]
    pred = q * (emulated_s(q - DQ) - 1)
    assert np.max(np.abs(f - pred)) < 1e-4, np.max(np.abs(f - pred))


# ----------------------------------------------------------------------------------------------- pdf menu
RMAX, DR = 8.0, 0.01


def periodic_g(weights):
    xyz, site = fcc_positions()
    n = len(xyz)
    L = NC * A
    w = weights(site)
    nb = int(RMAX / DR) + 2
    H = np.zeros(nb)
    shifts = np.array(list(itertools.product((-1, 0, 1), repeat=3))) * L
    for sh in shifts:
        d = np.linalg.norm(xyz[:, None] - (xyz[None] + sh), axis=2)
        ww = w[:, None] * w[None] * np.ones_like(d)
        keep = (d > 1e-6) & (d <= RMAX)
        k = np.floor((d[keep] + DR / 2) / DR).astype(int)
        np.add.at(H, k, ww[keep])
    r = np.arange(nb) * DR
    rho0 = n / L ** 3
    with np.errstate(divide="ignore", invalid="ignore"):
        g = H / (n * DR * r) - 4 * np.pi * r * rho0
    return r, g


def compare(path, weights):
    r, g = read_pdf_menu_file(path)
    rp, gp = periodic_g(weights)
    sel = (r > 0.5) & (r <= RMAX - 0.1)
    idx = np.rint(r[sel] / DR).astype(int)
    assert np.allclose(r[sel], rp[idx], atol=1e-4)
    err = np.max(np.abs(g[sel] - gp[idx]))
    assert err < 1e-6, err


def test_g1_neutron_single_species():
    step = PdfCalc("g.dat", rmax=RMAX, dr=DR, radiation="neutron")
    r = StructureJob(SRC, [step], readback="none").run(runner(), {"cu_b05.cll": CU})
    compare(r.run.output("g.dat"), lambda site: np.ones(len(site)))


def test_g2_xray_two_species_weights():
    step = PdfCalc("g.dat", rmax=RMAX, dr=DR, radiation="xray")
    job = StructureJob(SRC, [ReplaceRandom("CU", "AU", 1.0, site=1), step], seed=SEED, readback="none")
    r = job.run(runner(), {"cu_b05.cll": CU})
    fau, fcu = F0["AU"], F0["CU"]
    mean = (3 * fcu + fau) / 4
    compare(r.run.output("g.dat"), lambda site: np.where(site == 0, fau, fcu) / mean)

