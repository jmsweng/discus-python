"""Live checks of FourierPlane and RmcFit on the real engine (opt-in: DISCUS_LIVE=1).

Small on purpose: Cu fcc 4x4x4 cells (256 atoms) and one 31 x 21 plane hk0.

* R1 (engine-ran): the plane DISCUS writes has the size and orientation ``rmc data`` reads
  (superstructure peak of L1_2 Cu3Au at h=1, k=0; fundamental at h=2, k=0).
* R2 (pre-registered, verified level): an RMC fit with swaps only and sigma = 0, started from a
  random Cu0.75Au0.25 crystal and fitted to the L1_2 plane, (a) accepts no uphill move, (b) never
  increases chi2 between log lines and ends below the first logged value, (c) keeps the Au count,
  and (d) gives a first-neighbour Warren-Cowley alpha below -0.15 (perfect L1_2: -1/3; the start:
  |alpha| < 0.1).
* R3 (engine-ran): the same fit with the same seed gives the same file.

Environment as in test_discus_ops_live.py: DISCUS_LIVE=1; DISCUS_LIVE_WORK (no spaces; default
C:\\tmp\\rmclive or /tmp/rmclive); DISCUS_WSL_DISTRO (default Ubuntu-22.04).
"""
from __future__ import annotations

import os
import platform
from pathlib import Path

import numpy as np
import pytest

from discuspdf.discus_ops import ReplaceRandom, StructureJob
from discuspdf.discus_rmc import FourierPlane, PlaneGrid, RmcFit, RmcPlane, RmcSwap, read_nipl, run_rmc
from discuspdf.engine import Engine
from discuspdf.fourier_macro import Seed, StructureSource
from discuspdf.suite_runner import SuiteRunner

pytestmark = pytest.mark.discus_live
if os.environ.get("DISCUS_LIVE") != "1":
    pytest.skip("live DISCUS test; set DISCUS_LIVE=1 to enable", allow_module_level=True)

ROOT = Path(__file__).resolve().parents[1]
CU = ROOT / "tests" / "fixtures" / "discus_ops" / "cu_b05.cll"
WORK = Path(os.environ.get("DISCUS_LIVE_WORK",
                           r"C:\tmp\rmclive" if platform.system() == "Windows" else "/tmp/rmclive"))
SEED = Seed((11, 22, 33, 44, 55, 66, 77, 88))
N = 4
A = 3.615
GRID = PlaneGrid(ll=(0, -1, 0), lr=(3, -1, 0), ul=(0, 1, 0), nx=31, ny=21)
SRC = StructureSource.cell("cu_b05.cll", (N, N, N))


def runner() -> SuiteRunner:
    if " " in str(WORK):
        pytest.skip("DISCUS_LIVE_WORK must not contain spaces")
    if platform.system() == "Windows":
        engine = Engine(("wsl.exe", "-d", os.environ.get("DISCUS_WSL_DISTRO", "Ubuntu-22.04"), "--",
                         "discus_suite", "-noautorun", "-macro"), timeout_seconds=900)
    else:
        engine = Engine(timeout_seconds=900)
    return SuiteRunner(WORK, engine)


def species(s):
    return np.array([str(x).upper() for x in s.species])


def alpha1(s, a="CU", b="AU") -> float:
    """First-neighbour Warren-Cowley alpha = 1 - P(b next to a) / c_b, periodic cubic box."""
    xyz = np.asarray(s.cartesian)
    box = N * A
    d = xyz[:, None, :] - xyz[None, :, :]
    d -= box * np.round(d / box)
    r = np.linalg.norm(d, axis=2)
    nn = np.abs(r - A / np.sqrt(2)) < 0.2
    sp = species(s)
    assert (nn.sum(axis=1) == 12).all()
    p_b_given_a = nn[sp == a][:, sp == b].sum() / nn[sp == a].sum()
    return float(1 - p_b_given_a / np.mean(sp == b))


@pytest.fixture(scope="module")
def target():
    """The L1_2 plane: Au on site 1 of every cell."""
    job = StructureJob(SRC, [ReplaceRandom("CU", "AU", 1.0, site=1), FourierPlane(GRID, "target.nipl")],
                       seed=SEED, readback="none")
    r = job.run(runner(), {"cu_b05.cll": CU})
    return r.run.output("target.nipl")


def fit(tries=3000):
    return RmcFit(planes=(RmcPlane("target.nipl", GRID),), moves=(RmcSwap(),), select=("CU", "AU"),
                  tries=tries, sigma=0.0, display=250)


def test_r1_plane_size_and_orientation(target):
    z, _ = read_nipl(target)
    assert z.shape == (21, 31)
    med = float(np.median(z))
    assert z[10, 10] > 100 * med, (z[10, 10], med)      # h=1, k=0: L1_2 superstructure
    assert z[10, 20] > z[10, 10]                          # h=2, k=0: fundamental
    assert z[10, 5] < 0.01 * z[10, 10]                    # h=0.5: no peak


def test_r2_fit_orders_the_crystal(target):
    start = StructureJob(SRC, [ReplaceRandom("CU", "AU", 0.25)], seed=SEED).run(
        runner(), {"cu_b05.cll": CU}).structure
    a0 = alpha1(start)
    assert abs(a0) < 0.1, a0
    job = StructureJob(SRC, [ReplaceRandom("CU", "AU", 0.25), fit()], seed=SEED)
    out = run_rmc(job, runner(), {"cu_b05.cll": CU, "target.nipl": target})
    log = out.log
    assert log.final.tries == 3000
    assert all(p.accepted_uphill == 0 for p in log.progress)
    chi = [p.chi2 for p in log.progress]
    assert all(b <= a for a, b in zip(chi, chi[1:])), chi
    assert chi[-1] < chi[0], chi
    s = out.result.structure
    assert (species(s) == "AU").sum() == (species(start) == "AU").sum()
    a1 = alpha1(s)
    assert a1 < -0.15, (a0, a1, chi[0], chi[-1])


def test_r3_same_seed_same_fit(target):
    job = StructureJob(SRC, [ReplaceRandom("CU", "AU", 0.25), fit(500)], seed=SEED)
    a = run_rmc(job, runner(), {"cu_b05.cll": CU, "target.nipl": target})
    b = run_rmc(job, runner(), {"cu_b05.cll": CU, "target.nipl": target})
    assert a.result.run.output("result.stru").read_bytes() == b.result.run.output("result.stru").read_bytes()
