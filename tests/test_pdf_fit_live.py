"""Live checks of PdfFit and PdfPartials (opt-in: DISCUS_LIVE=1). Pre-registered (verified level).

Crystal: Cu fcc 3x3x3 cells (108 atoms, B 0.5, ``tests/fixtures/discus_ops/cu_b05.cll``).

Registered before any run of these steps (2026-11-05; numbers from probe rounds 2-3, ``tools/probes/pdf_fit_probes.py``):

* P1 (fit): target = G(r) of ordered L1_2 Cu3Au (X-ray, Gaussian peaks, qmax 0, r 0.01..8 step 0.01); start =
  random Cu0.75Au0.25 from the seed; PdfFit with swaps only, sigma 0, 2000 tries, fixed scale 1, fit range
  0.01..7.5: no uphill acceptance, chi2 never rises and ends below its first value, the Au count is kept, and the
  first-shell Warren-Cowley alpha(Cu-Au) of the result is below -0.25 (perfect order -1/3; the probe without a fit
  range reached -0.318).
* P2 (two code paths): the G(r) saved after the fit and PdfCalc of the fitted crystal with the same settings agree
  within 2e-3 of max |G| for r <= 7.5 (probe: 8e-4).
* P3 (partials): ternary Cu2AuAg (Au site 1, Ag site 2), neutron, no Gaussian, qmax 0, r 0..8 step 0.001: each of
  G_AuAu, G_CuCu, G_AuCu, G_AuAg from read_pdf_partials equals w_ab P_ab / (N dr r) - (2 - delta_ab) c_a c_b w_ab
  4 pi r rho0 (P_ab the pair count over both orders, w_ab = b_a b_b / <b>^2) within 1e-6 of max |G_ab| for
  0.5 < r <= 7.9; and a total PdfCalc after the partials in the same session equals a total in a fresh one within
  1e-12 of max |G| (the selections are restored).

Environment: DISCUS_LIVE=1; DISCUS_LIVE_WORK (default C:\\tmp\\pdffitlive or /tmp/pdffitlive); DISCUS_WSL_DISTRO.
"""
from __future__ import annotations

import itertools
import os
import platform
from pathlib import Path

import numpy as np
import pytest

from discuspdf.discus_ops import ReplaceRandom, StructureJob
from discuspdf.discus_pdfmenu import (PdfCalc, PdfData, PdfFit, PdfPartials, read_pdf_menu_file, read_pdf_partials,
                                      run_pdf_fit)
from discuspdf.discus_rmc import RmcSwap
from discuspdf.fourier_macro import Seed, StructureSource
from tests.live_engine import live_runner

pytestmark = pytest.mark.discus_live
if os.environ.get("DISCUS_LIVE") != "1":
    pytest.skip("live DISCUS test; set DISCUS_LIVE=1 to enable", allow_module_level=True)

ROOT = Path(__file__).resolve().parents[1]
CU = ROOT / "tests" / "fixtures" / "discus_ops" / "cu_b05.cll"
WORK = Path(os.environ.get("DISCUS_LIVE_WORK",
                           r"C:\tmp\pdffitlive" if platform.system() == "Windows" else "/tmp/pdffitlive"))
SRC = StructureSource.cell("cu_b05.cll", (3, 3, 3))
SEED = Seed((11, 22, 33, 44, 55, 66, 77, 88))
A, NC = 3.615, 3
XRAY = dict(radiation="xray", thermal=True, rmax=8.0, dr=0.01)


def runner():
    return live_runner(WORK, 900)


def alpha1(s) -> float:
    xyz = np.asarray(s.cartesian)
    box = NC * A
    d = xyz[:, None, :] - xyz[None, :, :]
    d -= box * np.round(d / box)
    nn = np.abs(np.linalg.norm(d, axis=2) - A / np.sqrt(2)) < 0.2
    sp = np.array([str(x).upper() for x in s.species])
    p = nn[sp == "CU"][:, sp == "AU"].sum() / nn[sp == "CU"].sum()
    return float(1 - p / np.mean(sp == "AU"))


@pytest.fixture(scope="module")
def fitted():
    t = StructureJob(SRC, [ReplaceRandom("CU", "AU", 1.0, site=1), PdfCalc("target.dat", **XRAY)], seed=SEED,
                     readback="none").run(runner(), {"cu_b05.cll": CU})
    data = PdfData.from_file(t.run.output("target.dat"))
    fit = PdfFit(PdfCalc("fit.dat", **XRAY), data, (RmcSwap(),), ("CU", "AU"), 2000, sigma=0.0, display=200,
                 fit_scale=False, scale=1.0, fit_range=(0.01, 7.5))
    start = StructureJob(SRC, [ReplaceRandom("CU", "AU", 0.25)], seed=SEED).run(runner(), {"cu_b05.cll": CU}).structure
    out = run_pdf_fit(StructureJob(SRC, [ReplaceRandom("CU", "AU", 0.25), fit], seed=SEED), runner(),
                      {"cu_b05.cll": CU, "target.dat": t.run.output("target.dat")})
    return start, out


def test_p1_fit_orders_the_crystal(fitted):
    start, out = fitted
    log = out.log
    assert log.final.tries == 2000
    assert all(p.accepted_uphill == 0 for p in log.progress)
    chi = [p.chi2 for p in log.progress]
    assert all(b <= a for a, b in zip(chi, chi[1:])) and chi[-1] < chi[0], chi
    sp = lambda s: sum(str(x).upper() == "AU" for x in s.species)  # noqa: E731
    assert sp(out.result.structure) == sp(start)
    a0, a1 = alpha1(start), alpha1(out.result.structure)
    assert a1 < -0.25, (a0, a1, chi[0], chi[-1])


def test_p2_fit_path_equals_calc_path(fitted):
    _, out = fitted
    r, g_fit = read_pdf_menu_file(out.result.run.output("fit.dat"))
    calc = StructureJob(out.result.structure, [PdfCalc("recalc.dat", **XRAY)], readback="none").run(runner())
    rc, g_calc = read_pdf_menu_file(calc.run.output("recalc.dat"))
    n = min(len(r), len(rc))
    assert np.allclose(r[:n], rc[:n], atol=1e-6)
    sel = r[:n] <= 7.5
    err = np.max(np.abs(g_fit[:n][sel] - g_calc[:n][sel])) / np.max(np.abs(g_calc[:n][sel]))
    assert err <= 2e-3, err


# ----------------------------------------------------------------------------------------------- partials
SITES = ("AU", "AG", "CU", "CU")
B = {"CU": 7.718, "AU": 7.63, "AG": 5.922}  # neutron b (fm), DISCUS 6.23.01 table
CONC = {"AU": 0.25, "AG": 0.25, "CU": 0.5}
NEUTRON = dict(radiation="neutron", rmax=8.0, dr=0.001)


def pair_count(a, b, dr=0.001, rmax=8.0):
    basis = np.array([[0, 0, 0], [0, .5, .5], [.5, 0, .5], [.5, .5, 0]])
    cells = np.array(list(itertools.product(range(NC), repeat=3)))
    xyz = ((cells[:, None, :] + basis[None]).reshape(-1, 3)) * A
    sp = np.array(SITES * len(cells))
    L = NC * A
    H = np.zeros(int(rmax / dr) + 2)
    si, sj = sp == a, sp == b
    for sh in itertools.product((-1, 0, 1), repeat=3):
        d = np.linalg.norm(xyz[:, None] - (xyz[None] + np.array(sh) * L), axis=2)
        m = (si[:, None] & sj[None]) | (sj[:, None] & si[None])
        keep = m & (d > 1e-6) & (d <= rmax)
        np.add.at(H, np.floor((d[keep] + dr / 2) / dr).astype(int), 1)
    return H


def test_p3_partials_and_restored_selection():
    order = [ReplaceRandom("CU", "AU", 1.0, site=1), ReplaceRandom("CU", "AG", 1.0, site=2)]
    step = PdfPartials(PdfCalc("unused.dat", **NEUTRON), (("AU", "AU"), ("CU", "CU"), ("AU", "CU"), ("AU", "AG")))
    run = StructureJob(SRC, [*order, step, PdfCalc("after.dat", **NEUTRON)], seed=SEED, readback="none").run(
        runner(), {"cu_b05.cll": CU}).run
    fresh = StructureJob(SRC, [*order, PdfCalc("total.dat", **NEUTRON)], seed=SEED, readback="none").run(
        runner(), {"cu_b05.cll": CU}).run
    parts = read_pdf_partials(step, run)
    mean = sum(CONC[x] * B[x] for x in CONC)
    n, rho0 = 108, 108 / (NC * A) ** 3
    for (a, b), (r, g) in parts.items():
        idx = np.rint(r / 0.001).astype(int)
        sel = (r > 0.5) & (r <= 7.9)
        w = B[a] * B[b] / mean ** 2
        model = w * pair_count(a, b)[idx[sel]] / (n * 0.001 * r[sel]) \
            - (1 if a == b else 2) * CONC[a] * CONC[b] * w * 4 * np.pi * r[sel] * rho0
        err = np.max(np.abs(g[sel] - model)) / np.max(np.abs(g[sel]))
        assert err < 1e-6, ((a, b), err)
    _, g_after = read_pdf_menu_file(run.output("after.dat"))
    _, g_total = read_pdf_menu_file(fresh.output("total.dat"))
    assert np.max(np.abs(g_after - g_total)) <= 1e-12 * np.max(np.abs(g_total))
