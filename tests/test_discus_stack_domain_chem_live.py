"""Live checks of Stacking, Domains and ChemAnalysis on the real engine (opt-in: DISCUS_LIVE=1).

Engine-ran level, small models (at most 864 atoms). Environment as in test_discus_ops_live.py:
DISCUS_LIVE=1; DISCUS_LIVE_WORK (no spaces; default C:\\tmp\\s2blive or /tmp/s2blive);
DISCUS_WSL_DISTRO (default Ubuntu-22.04).
"""
from __future__ import annotations

import os
import platform
import shutil
from pathlib import Path

import numpy as np
import pytest

from discuspdf.discus_domain import Domains
from discuspdf.discus_ops import Purge, ReplaceRandom, StructureJob, read_box
from discuspdf.discus_stack import Stacking
from discuspdf.engine import Engine
from discuspdf.fourier_macro import Seed, StructureSource
from discuspdf.suite_runner import SuiteRunner

pytestmark = pytest.mark.discus_live
if os.environ.get("DISCUS_LIVE") != "1":
    pytest.skip("live DISCUS test; set DISCUS_LIVE=1 to enable", allow_module_level=True)

ROOT = Path(__file__).resolve().parents[1]
STACK = ROOT / "tests" / "fixtures" / "discus_stack"
CU = ROOT / "tests" / "fixtures" / "discus_ops" / "cu_b05.cll"
WORK = Path(os.environ.get("DISCUS_LIVE_WORK",
                           r"C:\tmp\s2blive" if platform.system() == "Windows" else "/tmp/s2blive"))
SEED = Seed((11, 22, 33, 44, 55, 66, 77, 88))
SEED2 = Seed((12, 23, 34, 45, 56, 67, 78, 89))


def runner() -> SuiteRunner:
    if " " in str(WORK):
        pytest.skip("DISCUS_LIVE_WORK must not contain spaces")
    if platform.system() == "Windows":
        engine = Engine(("wsl.exe", "-d", os.environ.get("DISCUS_WSL_DISTRO", "Ubuntu-22.04"), "--",
                         "discus_suite", "-noautorun", "-macro"), timeout_seconds=600)
    else:
        engine = Engine(timeout_seconds=600)
    return SuiteRunner(WORK, engine)


def species(s):
    return np.array([str(x).upper() for x in s.species])


# --- stacking ------------------------------------------------------------------------------

ABAB = dict(layers=("cu_layer.stru", "ag_layer.stru"), correlation=((0, 1), (1, 0)),
            translations=(((0, 0, 0.5), (0, 0, 0.5)), ((0, 0, 0.5), (0, 0, 0.5))),
            modulus=((1, 0, 0), (0, 1, 0)), number=8, first=1)
LAYERS = {n: STACK / n for n in ("cu_layer.stru", "ag_layer.stru")}


def test_s1_fault_free_abab_stack():
    r = StructureJob(Stacking(**ABAB)).run(runner(), LAYERS)
    s = r.structure
    sp = species(s)
    assert s.n_atoms == 8
    assert list(sp) == ["CU", "AG"] * 4
    z = np.asarray(s.cartesian)[:, 2]
    assert np.allclose(np.diff(z), 2.0, atol=1e-4)  # 0.5 of the 4 A layer cell per layer


def _fcc(alpha, **kw):
    pos = {0: (0.0, 0.0), 1: (1 / 3, 2 / 3), 2: (2 / 3, 1 / 3)}
    trans = tuple(tuple(((pos[j][0] - pos[i][0]) % 1.0, (pos[j][1] - pos[i][1]) % 1.0, 1 / 3)
                        for j in range(3)) for i in range(3))
    args = dict(layers=("fcc_a.stru", "fcc_b.stru", "fcc_c.stru"),
                correlation=((0, 1 - alpha, alpha), (alpha, 0, 1 - alpha), (1 - alpha, alpha, 0)),
                translations=trans, modulus=((1, 0, 0), (0, 1, 0)), average=(0, 0, 1 / 3),
                number=30, first=1, modulus_on=True)
    args.update(kw)
    return Stacking(**args)


FCC_LAYERS = {n: STACK / n for n in ("fcc_a.stru", "fcc_b.stru", "fcc_c.stru")}


def _positions(s):
    """In-plane position class (0 = A, 1 = B, 2 = C) of each layer, from Cartesian x, y."""
    a = 2.863782
    xy = np.asarray(s.cartesian)[:, :2]
    # A at (0,0), B at (1/3,2/3), C at (2/3,1/3) of the hexagonal cell
    basis = np.array([[a, 0.0], [-a / 2, a * np.sqrt(3) / 2]])
    f = np.linalg.solve(basis.T, xy.T).T % 1.0
    ref = np.array([[0, 0], [1 / 3, 2 / 3], [2 / 3, 1 / 3]])
    d = np.abs(((f[:, None, :] - ref[None]) + 0.5) % 1.0 - 0.5).sum(axis=2)
    return d.argmin(axis=1), d.min(axis=1)


def test_s2_fcc_abc_without_faults_and_markov_steps():
    perfect = StructureJob(_fcc(0.0)).run(runner(), FCC_LAYERS).structure
    cls, err = _positions(perfect)
    assert perfect.n_atoms == 30 and np.max(err) < 1e-3
    assert list(cls) == [0, 1, 2] * 10  # ABC... from first = 1
    z = np.sort(np.asarray(perfect.cartesian)[:, 2])
    assert np.allclose(np.diff(z), 7.014805 / 3, atol=1e-3)
    faulted = StructureJob(_fcc(0.3), seed=SEED).run(runner(), FCC_LAYERS).structure
    cls, err = _positions(faulted)
    assert np.max(err) < 1e-3
    steps = (np.diff(cls) % 3)  # 1 = forward (A->B), 2 = reversed (A->C); 0 never
    assert set(steps) <= {1, 2} and (steps == 2).sum() > 0


def test_s3_same_seed_same_stack():
    a = StructureJob(_fcc(0.3), seed=SEED).run(runner(), FCC_LAYERS)
    b = StructureJob(_fcc(0.3), seed=SEED).run(runner(), FCC_LAYERS)
    assert a.run.output("result.stru").read_bytes() == b.run.output("result.stru").read_bytes()


# --- domains -------------------------------------------------------------------------------

# Free-format atom lines on purpose: 'format numbers' makes DISCUS read fixed Fortran columns
# (atom_line.f90 44444/44445), which a hand-written file rarely matches (error -49 live).
PSEUDO = """title pseudo atoms
spcgr F m -3 m
cell    3.615000,  3.615000,  3.615000, 90.000000, 90.000000, 90.000000
atoms
DOM1 0.000000 0.000000 0.000000 0.000000
"""

SPHERE = ((1.5, 0, 0, 0), (0, 1.5, 0, 0), (0, 0, 1.5, 0))


@pytest.fixture(scope="module")
def domain_files(tmp_path_factory):
    """Content: the Cu 6x6x6 crystal saved by DISCUS with every atom renamed AG (DISCUS's own
    centred frame); pseudo atom file: one DOM1 at the origin with the host lattice."""
    d = tmp_path_factory.mktemp("domain")
    job = StructureJob(StructureSource.cell("cu_b05.cll", (6, 6, 6)), [ReplaceRandom("CU", "AG", 1.0)], seed=SEED)
    content = job.run(runner(), {"cu_b05.cll": CU}).run.output("result.stru")
    shutil.copyfile(content, d / "dom_ag.stru")
    (d / "pseudo.stru").write_text(PSEUDO)
    return {"cu_b05.cll": CU, "dom_ag.stru": d / "dom_ag.stru", "pseudo.stru": d / "pseudo.stru"}


def _domain_job(fuzzy, steps=(), seed=SEED, **kw):
    step = Domains(content="dom_ag.stru", inputfile="pseudo.stru", name="DOM1", character="sphere",
                   shape=SPHERE, fuzzy=fuzzy, **kw)
    return StructureJob(StructureSource.cell("cu_b05.cll", (6, 6, 6)), [step, *steps], seed=seed,
                        readback="box")


def test_d1_sphere_domain_replaces_host_atoms(domain_files):
    s = _domain_job(1.0).run(runner(), domain_files).structure
    sp = species(s)
    n_ag, n_void = int((sp == "AG").sum()), int((sp == "VOID").sum())
    assert n_ag > 0
    assert n_void == n_ag  # the same lattice: each inserted Ag sits on one host Cu
    assert s.n_atoms == 864 + n_ag


def test_d2_no_removal_without_fuzzy(domain_files):
    s = _domain_job(0.0).run(runner(), domain_files).structure
    sp = species(s)
    assert (sp == "VOID").sum() == 0 and (sp == "AG").sum() > 0


def test_d3_same_result_with_other_seed_when_no_size_spread(domain_files):
    a = _domain_job(1.0).run(runner(), domain_files)
    b = _domain_job(1.0, seed=SEED2).run(runner(), domain_files)
    assert a.run.output("result.stru").read_bytes() == b.run.output("result.stru").read_bytes()


# --- chem ----------------------------------------------------------------------------------

def test_c1_chem_occupational_correlation_matches_python():
    from discuspdf.discus_chem import (ChemAnalysis, ChemOccupational, CorrelationVector, Neighbourhood,
                                       VectorNeighbourhood, run_chem)
    vec = CorrelationVector(1, 1, 1, (1, 0, 0))
    analysis = ChemAnalysis(vectors=(vec,), neighbourhoods=(Neighbourhood(VectorNeighbourhood(1)),),
                            calculations=(ChemOccupational("CU", "AU"),))
    job = StructureJob(StructureSource.cell("cu_b05.cll", (6, 6, 6)),
                       [ReplaceRandom("CU", "AU", 0.5), analysis], seed=SEED, readback="none")
    out = run_chem(job, runner(), {"cu_b05.cll": CU})
    assert out.results, "no correlation parsed from the log"
    value = out.results[0].value
    assert abs(value) < 0.2  # random alloy: near zero (864 atoms, one vector)
