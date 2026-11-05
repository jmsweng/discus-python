"""Live checks of discus_ops steps on the real engine (opt-in: DISCUS_LIVE=1). Engine-ran level.

Small crystals on purpose: Cu fcc 6x6x6 cells (864 atoms) from tests/fixtures/discus_ops/cu_b05.cll
(B = 0.5 A^2). Each check compares the saved crystal with what the step should do, with
tolerances from counting statistics where random numbers are involved.

Environment: DISCUS_LIVE=1; DISCUS_LIVE_WORK (no spaces; default C:\\tmp\\opslive on Windows,
/tmp/opslive); DISCUS_WSL_DISTRO (default Ubuntu-22.04).
"""
from __future__ import annotations

import math
import os
import platform
from pathlib import Path

import numpy as np
import pytest

from discuspdf.discus_ops import Purge, Remove, ReplaceAtom, ReplaceRandom, StructureJob, Switch, Thermal
from discuspdf.fourier_macro import Seed, StructureSource
from discuspdf.suite_runner import SuiteRunner
from tests.live_engine import live_runner

pytestmark = pytest.mark.discus_live
if os.environ.get("DISCUS_LIVE") != "1":
    pytest.skip("live DISCUS test; set DISCUS_LIVE=1 to enable", allow_module_level=True)

ROOT = Path(__file__).resolve().parents[1]
CU = ROOT / "tests" / "fixtures" / "discus_ops" / "cu_b05.cll"
WORK = Path(os.environ.get("DISCUS_LIVE_WORK",
                           r"C:\tmp\opslive" if platform.system() == "Windows" else "/tmp/opslive"))
SEED = Seed((11, 22, 33, 44, 55, 66, 77, 88))
N = 6  # cells per edge; 4 atoms per cell -> 864 atoms


def runner() -> SuiteRunner:
    return live_runner(WORK, 600)


def job(steps, seed=SEED):
    return StructureJob(StructureSource.cell("cu_b05.cll", (N, N, N)), steps, seed=seed)


def run(steps, seed=SEED):
    return job(steps, seed).run(runner(), {"cu_b05.cll": CU})


def species(s):
    return np.array([str(x).upper() for x in s.species])


def test_plain_read_and_save_is_the_ideal_crystal():
    r = run([])
    s = r.structure
    assert s.n_atoms == 4 * N**3
    assert set(species(s)) == {"CU"}
    # fractions of the whole crystal times N are integers or halves (fcc)
    x = np.asarray(s.frac) * N
    assert np.allclose(x * 2, np.round(x * 2), atol=1e-6)


def test_replace_random_fraction_and_reproducibility():
    a = run([ReplaceRandom("CU", "AU", 0.25)])
    b = run([ReplaceRandom("CU", "AU", 0.25)])
    n = a.structure.n_atoms
    frac = float(np.mean(species(a.structure) == "AU"))
    sigma = math.sqrt(0.25 * 0.75 / n)
    assert abs(frac - 0.25) < 4 * sigma, (frac, sigma)
    # same seed, same draws
    assert (species(a.structure) == species(b.structure)).all()


def test_thermal_mean_square_displacement():
    ideal = run([]).structure
    hot = run([Thermal()]).structure
    assert hot.n_atoms == ideal.n_atoms
    a = 3.615 * N  # whole-crystal edge in Angstrom (cubic)
    d = (np.asarray(hot.frac) - np.asarray(ideal.frac))
    d -= np.round(d)  # periodic images
    u2 = float(np.mean((d * a) ** 2))  # per Cartesian direction
    expected = 0.5 / (8 * math.pi ** 2)  # B / 8 pi^2
    # 3*864 samples of a Gaussian: relative s.e. of the mean square sqrt(2/2592) = 2.8 %
    assert abs(u2 / expected - 1) < 0.12, (u2, expected)


def test_replace_switch_remove_are_deterministic():
    r = run([ReplaceAtom(1, "AU", b=0.5), Switch(1, 2), Remove(5, 7)]).structure
    sp = species(r)
    assert r.n_atoms == 4 * N**3  # removed atoms stay in the list as vacancies
    assert sp[0] == "CU" and sp[1] == "AU"  # atom 1 became Au, then 1 and 2 switched
    assert (sp == "AU").sum() == 1
    assert list(sp[4:7]) == ["VOID"] * 3, list(sp[:8])


def test_purge_drops_the_vacancies():
    from discuspdf.discus_ops import ResultNotAStructure
    with pytest.raises(ResultNotAStructure) as exc:
        run([Remove(5, 7), Purge()])
    text = exc.value.run.output("result.stru").read_text()
    atoms = text.lower().split("atoms", 1)[1].strip().splitlines()[1:]  # skip the column header line
    assert sum(1 for line in atoms if line.strip()) == 4 * N**3 - 3


def test_longitudinal_wave_displaces_along_the_wave_vector_only():
    from discuspdf.discus_ops import Wave
    ideal = run([]).structure
    amp = 0.05
    w = Wave("longitudinal", (1, 0, 0), 3.615 * N, amplitude=amp, phase=0)
    moved = run([w]).structure
    d = (np.asarray(moved.frac) - np.asarray(ideal.frac))
    d -= np.round(d)
    u = d * 3.615 * N  # Angstrom (cubic)
    assert np.max(np.abs(u[:, 1:])) < 1e-6  # nothing perpendicular to [100]
    assert 0.9 * amp < np.max(np.abs(u[:, 0])) <= amp + 1e-6
