"""Live checks of the phase 2c-A Fourier options (opt-in: DISCUS_LIVE=1). Pre-registered (verified level).

Model: Cu fcc 4x4x4 cells (256 atoms, ``tests/fixtures/discus_ops/cu_b05.cll``), ``temp ignore``,
``set aver 0``, plane hk0 from (0,0,0) to (4,4,0) in steps of 0.5. At the fcc reflections (h, k even)
the intensity is (256 f)^2, so f = sqrt(I)/256 at each of (200), (220), (400), (420), (440).

Registered before any run (2026-11-03), with the numbers computed from the 6.23.01 tables in
``lib_f90/element_data_mod.f90`` (Cu: Waasmaier per_wa*/per_wb*/per_wc, electron per_ea*/per_eb*,
Mo Ka1 f' = 0.3201, f'' = 1.2651):

* X1: X-ray f at the five reflections equals the Waasmaier 5-Gaussian value within 0.5 % (the
  sin(theta)/lambda table step is 0.001 A^-1).
* E1: electron f equals DISCUS's 5-Gaussian electron table within 0.5 %, and agrees with the
  Mott-Bethe formula f_e = 0.023934 (Z - f_x) / s^2 (A, s in A^-1) from the X-ray run within 5 %
  (computed from the tables: 0.2-1.9 %).
* A1: ``disp anom`` with Mo Ka1 (internal f', f''): I_anom / I = ((f0 + f')^2 + f''^2) / f0^2 within 1e-5.
* A2: ``delf CU, 0.5, 2.0`` with a numeric wavelength: the same formula with (0.5, 2.0), within 1e-5.
* W1: ``energy 17.4793`` keV (X-ray) gives the same intensities as ``wvle moa1`` (no dispersion), within 1e-12
  relative: on an hkl grid the wavelength is metadata only.

Environment as in the other live tests: DISCUS_LIVE=1; DISCUS_LIVE_WORK (default C:\\tmp\\fourxlive or
/tmp/fourxlive); DISCUS_WSL_DISTRO (default Ubuntu-22.04).
"""
from __future__ import annotations

import math
import os
import platform
from pathlib import Path

import numpy as np
import pytest

from discuspdf import Lattice, ReciprocalGrid
from discuspdf.discus_output import read_discus_volume
from discuspdf.engine import Engine
from discuspdf.fourier_macro import FourierConfig, StructureSource
from discuspdf.suite_runner import SuiteRunner

pytestmark = pytest.mark.discus_live
if os.environ.get("DISCUS_LIVE") != "1":
    pytest.skip("live DISCUS test; set DISCUS_LIVE=1 to enable", allow_module_level=True)

ROOT = Path(__file__).resolve().parents[1]
CU = ROOT / "tests" / "fixtures" / "discus_ops" / "cu_b05.cll"
WORK = Path(os.environ.get("DISCUS_LIVE_WORK",
                           r"C:\tmp\fourxlive" if platform.system() == "Windows" else "/tmp/fourxlive"))
A = 3.615
N = 256
SRC = StructureSource.cell("cu_b05.cll", (4, 4, 4))
GRID = ReciprocalGrid.from_ranges(Lattice.from_parameters(A, A, A, 90, 90, 90), (0, 4, 9), (0, 4, 9), (0, 0, 1))
REFLECTIONS = [(2, 0), (2, 2), (4, 0), (4, 2), (4, 4)]

# lib_f90/element_data_mod.f90, Cu (index of 'CU  ' in per_name)
WAAS_A = (14.014192, 4.784577, 5.056806, 1.457971, 6.932996)
WAAS_B = (3.73828, 0.003744, 13.034982, 72.554793, 0.265666)
WAAS_C = -3.254477
ELEC_A = (0.3501, 1.6558, 1.9582, 0.2134, 1.4109)
ELEC_B = (0.1867, 1.9917, 11.3396, 53.2619, 63.252)
MOA1 = (0.3201, 1.2651)


def runner() -> SuiteRunner:
    if " " in str(WORK):
        pytest.skip("DISCUS_LIVE_WORK must not contain spaces")
    if platform.system() == "Windows":
        engine = Engine(("wsl.exe", "-d", os.environ.get("DISCUS_WSL_DISTRO", "Ubuntu-22.04"), "--",
                         "discus_suite", "-noautorun", "-macro"), timeout_seconds=600)
    else:
        engine = Engine(timeout_seconds=600)
    return SuiteRunner(WORK, engine)


def bragg(**options) -> dict:
    """f = sqrt(I)/N at each reflection, from one DISCUS run."""
    cfg = FourierConfig(GRID, (("int", "plane.h5"),), **options)
    run = runner().run(cfg.render(SRC), {"cu_b05.cll": CU}, ["plane.h5"])
    vol = read_discus_volume(run.output("plane.h5"))
    h, k, _ = vol.hkl_axes()
    out = {}
    for hh, kk in REFLECTIONS:
        i, j = int(np.argmin(abs(h - hh))), int(np.argmin(abs(k - kk)))
        out[(hh, kk)] = math.sqrt(vol.values[i, j, 0]) / N
    return out


def s_of(hh, kk):
    return math.hypot(hh, kk) / A / 2


def gauss(a, b, c, s):
    return c + sum(x * math.exp(-y * s * s) for x, y in zip(a, b))


@pytest.fixture(scope="module")
def xray():
    return bragg()


def test_x1_xray_form_factor_is_waasmaier(xray):
    for (hh, kk), f in xray.items():
        ref = gauss(WAAS_A, WAAS_B, WAAS_C, s_of(hh, kk))
        assert abs(f / ref - 1) < 0.005, ((hh, kk), f, ref)


def test_e1_electron_table_and_mott_bethe(xray):
    elec = bragg(radiation="electron", energy=200.0)
    for (hh, kk), f in elec.items():
        s = s_of(hh, kk)
        assert abs(f / gauss(ELEC_A, ELEC_B, 0.0, s) - 1) < 0.005, ((hh, kk), f)
        mott = 0.023934 * (29 - xray[(hh, kk)]) / s**2
        assert abs(f / mott - 1) < 0.05, ((hh, kk), f, mott)


def _anomalous_ratio(f0, fp, fpp):
    return ((f0 + fp) ** 2 + fpp ** 2) / f0 ** 2


def test_a1_internal_mo_ka1_dispersion(xray):
    anom = bragg(dispersion=True)
    for key, f in anom.items():
        assert abs((f / xray[key]) ** 2 / _anomalous_ratio(xray[key], *MOA1) - 1) < 1e-5, (key, f, xray[key])


def test_a2_user_delf_with_a_numeric_wavelength(xray):
    anom = bragg(wavelength=0.7093, dispersion=True, delf={"CU": (0.5, 2.0)})
    for key, f in anom.items():
        assert abs((f / xray[key]) ** 2 / _anomalous_ratio(xray[key], 0.5, 2.0) - 1) < 1e-5, (key, f, xray[key])


def test_w1_energy_equals_wavelength_for_intensities(xray):
    by_energy = bragg(energy=17.4793)
    for key, f in by_energy.items():
        assert abs(f / xray[key] - 1) < 1e-12, (key, f, xray[key])
