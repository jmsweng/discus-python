"""Live checks of ExpPdf / run_exp2pdf (opt-in: DISCUS_LIVE=1). Pre-registered (verified level).

Input: DISCUS ``powder`` (Debye, X-ray at 0.5 A, Q 0.5..17.9 step 0.01) of Cu fcc 5x5x5 cells (500 atoms,
``tests/fixtures/discus_ops/cu_b05.cll``): I(Q) and S(Q) from the same run. The "measured" pattern is
I(Q) + B(Q) with B = 0.05 <I> (1 + 0.1 Q), given to the wrapper **in 2-theta** (lambda 0.5) with B as the
background file, Fourier limits 1.0 .. 16.5.

Registered before any run (2026-11-05; a numpy emulation of the reduction, unit-weight polynomial of order 7,
gave a correlation of 0.998 and the same three peaks within 0.01 A):

* L1 (transform): the G(r) written equals (2/pi) sum dq F(Q) sin(Q r) over the F(Q) written by the same run,
  with a linear ramp from Q = 0 to its first point, within 1e-3 of max |G| (probe round 3: 2.8e-4).
* L2 (shape against the model): over 1.5 <= r <= 15, the Pearson correlation of G(r) with the model
  G_m(r) = (2/pi) sum dq Q (S(Q) - 1) sin(Q r), summed over the DISCUS S(Q) inside the Fourier limits DISCUS used
  (header ``lim_four_app``), is at least 0.98; the first three peaks (near 2.556, 3.615, 4.427 A) are within 0.02 A
  of the model's.
* L3 (composition is used): the same pattern reduced with composition ``Cu3 Au1`` gives a G(r) differing from the
  ``Cu`` one by more than 1 % of max |G|.

Engine-ran only: ``poly_order=5`` is accepted (header npoly 5); neutron and electron runs complete.

Environment: DISCUS_LIVE=1; DISCUS_LIVE_WORK (default C:\\tmp\\e2plive or /tmp/e2plive); DISCUS_WSL_DISTRO.
"""
from __future__ import annotations

import os
import platform
from pathlib import Path

import numpy as np
import pytest

from discuspdf.discus_exp2pdf import ExpPdf, run_exp2pdf
from discuspdf.discus_ops import StructureJob
from discuspdf.discus_powder import PowderPattern, read_powder_file
from discuspdf.fourier_macro import StructureSource
from tests.live_engine import live_runner

pytestmark = pytest.mark.discus_live
if os.environ.get("DISCUS_LIVE") != "1":
    pytest.skip("live DISCUS test; set DISCUS_LIVE=1 to enable", allow_module_level=True)

ROOT = Path(__file__).resolve().parents[1]
CU = ROOT / "tests" / "fixtures" / "discus_ops" / "cu_b05.cll"
WORK = Path(os.environ.get("DISCUS_LIVE_WORK",
                           r"C:\tmp\e2plive" if platform.system() == "Windows" else "/tmp/e2plive"))
LAM, DQ = 0.5, 0.01


def runner():
    return live_runner(WORK, 600)


@pytest.fixture(scope="module")
def pattern():
    step = PowderPattern((("I(Q)", "i.dat"), ("S(Q)", "s.dat")), qmin=0.5, qmax=17.9, dq=DQ, radiation="xray",
                         wavelength=LAM, stepr=0.001)
    r = StructureJob(StructureSource.cell("cu_b05.cll", (5, 5, 5)), [step], readback="none").run(
        runner(), {"cu_b05.cll": CU})
    q, i = read_powder_file(r.run.output("i.dat"))
    _, s = read_powder_file(r.run.output("s.dat"))
    files = WORK / "e2p_case"
    files.mkdir(parents=True, exist_ok=True)
    b = 0.05 * i.mean() * (1 + 0.1 * q)
    tth = 2 * np.degrees(np.arcsin(q * LAM / (4 * np.pi)))
    np.savetxt(files / "meas.xy", np.c_[tth, i + b])
    np.savetxt(files / "back.xy", np.c_[tth, b])
    return q, s, files


def reduce(pattern, **kw):
    _, _, files = pattern
    cfg = ExpPdf(files / "meas.xy", kw.pop("composition", "Cu"), axis="tth", wavelength=LAM,
                 background=files / "back.xy", fourier_qmin=1.0, fourier_qmax=16.5, rmin=0.5, rmax=20.0, **kw)
    return run_exp2pdf(cfg, runner(), files / f"in_{abs(hash(repr(sorted(kw.items())))) % 10**6}")


@pytest.fixture(scope="module")
def cu(pattern):
    return reduce(pattern)


def peak(r, g, near):
    sel = (r > near - 0.2) & (r < near + 0.2)
    return r[sel][np.argmax(g[sel])]


def test_l1_transform(cu):
    qf, f = cu.fq
    dq = np.median(np.diff(qf))
    n0 = int(round(qf[0] / dq))
    qq = np.concatenate([np.arange(1, n0) * dq, qf])
    ff = np.concatenate([np.arange(1, n0) * dq * f[0] / qf[0], f])
    pred = 2 / np.pi * dq * (ff[None] * np.sin(qq[None] * cu.r[:, None])).sum(axis=1)
    err = np.max(np.abs(cu.g - pred)) / np.max(np.abs(cu.g))
    assert err < 1e-3, err


def test_l2_shape_matches_the_model(pattern, cu):
    q, s, _ = pattern
    lo, hi = cu.fourier_limits
    sel = (q >= lo) & (q <= hi)
    r = cu.r[(cu.r >= 1.5) & (cu.r <= 15)]
    g = cu.g[(cu.r >= 1.5) & (cu.r <= 15)]
    gm = 2 / np.pi * DQ * ((q[sel] * (s[sel] - 1))[None] * np.sin(q[sel][None] * r[:, None])).sum(axis=1)
    corr = np.corrcoef(g, gm)[0, 1]
    peaks = [(peak(r, g, a), peak(r, gm, a)) for a in (2.556, 3.615, 4.427)]
    assert corr >= 0.98, (corr, peaks, (lo, hi))
    assert all(abs(a - b) <= 0.02 for a, b in peaks), peaks


def test_l3_composition_is_used(pattern, cu):
    cuau = reduce(pattern, composition="Cu3 Au1")
    n = min(len(cu.g), len(cuau.g))
    diff = np.max(np.abs(cuau.g[:n] - cu.g[:n])) / np.max(np.abs(cu.g))
    assert diff > 0.01, diff


def test_engine_ran_poly_neutron_electron(pattern):
    p5 = reduce(pattern, poly_order=5, outputs=())
    assert p5.header["npoly"] == ["5"], p5.header.get("npoly")
    for rad, extra in (("neutron", {}), ("electron", {"composition": {"Cu": 1}})):
        _, _, files = pattern
        cfg = ExpPdf(files / "meas.xy", extra.get("composition", "Cu"), radiation=rad, axis="tth", wavelength=LAM,
                     background=files / "back.xy", rmin=0.5, rmax=20.0, outputs=())
        res = run_exp2pdf(cfg, runner(), files / f"in_{rad}")
        assert len(res.g) > 100 and np.isfinite(res.g).all()
