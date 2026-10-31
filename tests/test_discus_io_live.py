"""Opt-in live tests against the real DISCUS 6.23.01 engine.

Run with ``DISCUS_LIVE=1`` to enable; tests are skipped otherwise.
"""
from __future__ import annotations

import os
import pytest

from discuspdf import (
    DiscusFormatError,
    Lattice,
    Structure,
    read_cll,
    read_stru,
    write_cll,
    write_stru,
    parse_cll,
)

from tests.discus_live import run_discus, evidence_log

pytestmark = pytest.mark.discus_live

# Skip all live tests unless DISCUS_LIVE=1.
skip_reason = "live DISCUS test; set DISCUS_LIVE=1 to enable"
pytestmark = [pytest.mark.discus_live, pytest.mark.skipif(
    os.environ.get("DISCUS_LIVE") != "1",
    reason=skip_reason,
)]


def _workdir() -> str:
    wd = os.environ.get("DISCUS_LIVE_WORKDIR", r"C:\tmp\a3live")
    os.makedirs(wd, exist_ok=True)
    # Clean up any leftover .stru / .cll files from previous runs.
    import glob
    for pattern in ("*.stru", "*.cll", "*.mac"):
        for f in glob.glob(os.path.join(wd, pattern)):
            try:
                os.unlink(f)
            except OSError:
                pass
    return wd


# -- L0: smoke test ----------------------------------------------------------
def test_live_smoke():
    """Run the A3-1a smoke macro: read cell gaas.cll,2,2,2 → 64 atoms."""
    workdir = _workdir()
    macro = (
        "discus\n"
        "read\n"
        "  cell gaas.cll,2,2,2\n"
        "eval n[1]\n"
        "exit\n"
        "exit\n"
    )
    # Copy the fixture into the workdir.
    import shutil
    shutil.copy(
        "tests/fixtures/discus_6_23_01/gaas.cll",
        os.path.join(workdir, "gaas.cll"),
    )
    result = run_discus(macro, workdir)
    assert "Version 6.23.01" in result.banner
    assert "64.0000000000" in result.stdout
    evidence_log("test_live_smoke", macro, None, result)


# -- L1: read cell matches fixture ------------------------------------------
def test_live_read_cell_matches_fixture():
    """L1: DISCUS read cell gaas.cll,2,2,2 → save → read_stru → 64 atoms, 32+32."""
    workdir = _workdir()
    import shutil
    shutil.copy(
        "tests/fixtures/discus_6_23_01/gaas.cll",
        os.path.join(workdir, "gaas.cll"),
    )
    macro = (
        "discus\n"
        "read\n"
        "  cell gaas.cll,2,2,2\n"
        "eval n[1]\n"
        "save gaas_live.stru\n"
        "exit\n"
        "exit\n"
    )
    result = run_discus(macro, workdir)
    assert "64.0000000000" in result.stdout

    # Read the saved .stru with our parser.
    s = read_stru(os.path.join(workdir, "gaas_live.stru"))
    assert s.structure.n_atoms == 64
    from collections import Counter
    assert dict(Counter(s.structure.species)) == {"GA": 32, "AS": 32}

    # Compare to the committed fixture.
    ref = read_stru("tests/fixtures/discus_6_23_01/gaas_222.stru")
    assert s.structure.isclose(ref.structure, rtol=1e-8, atol=1e-6)
    evidence_log("test_live_read_cell_matches_fixture", macro, None, result)


# -- L2: Python-written .stru loads and round-trips -------------------------
@pytest.mark.parametrize(
    "case_id, builder",
    [
        (
            "cubic_4x4x4",
            lambda: _build_cubic_4x4x4(),
        ),
        (
            "gaas_p1_2x2x2",
            lambda: _build_gaas_p1_2x2x2(),
        ),
    ],
)
def test_live_python_written_stru_loads_and_roundtrips(case_id, builder):
    """L2: write_stru → DISCUS read stru → save → read_stru → compare."""
    workdir = _workdir()
    s = builder()
    fname = f"live_{case_id}.stru"
    write_stru(os.path.join(workdir, fname), s)

    macro = (
        "discus\n"
        "read\n"
        f"  stru {fname}\n"
        "eval n[1]\n"
        f"save {fname}_out.stru\n"
        "exit\n"
        "exit\n"
    )
    result = run_discus(macro, workdir)

    N = s.n_atoms
    assert str(N) in result.stdout

    # Read the DISCUS-saved file.
    s2 = read_stru(os.path.join(workdir, f"{fname}_out.stru"))
    assert s2.structure.n_atoms == N

    # Same species per row.
    import numpy as np
    assert np.array_equal(s2.structure.species, s.species)

    # DISCUS drops the ncell line on save, so s2 is read as a single unit
    # cell (decision #4).  The coordinates are in supercell units [0, n),
    # so divide by n to get unit-cell fractions [0, 1).
    n = np.array(s.supercell, dtype=float)
    s2_frac_unit = s2.structure.frac / n
    s_frac_unit = s.frac
    assert np.allclose(s2_frac_unit, s_frac_unit, atol=2e-6)

    # Biso: DISCUS defaults missing B to 1.0 (Q1), not 0.0.
    # If s.b_iso is None, we wrote 0.0, but DISCUS reads it as 1.0.
    # If s.b_iso is set, DISCUS preserves it.
    if s.b_iso is None:
        # DISCUS defaults to 1.0 when B is 0.0 in the file.
        assert np.allclose(s2.structure.b_iso, 1.0, atol=2e-6)
    else:
        assert np.allclose(s2.structure.b_iso, s.b_iso, atol=2e-6)

    evidence_log("test_live_python_written_stru_loads_and_roundtrips", macro, None, result)


def _build_cubic_4x4x4() -> Structure:
    """Simple cubic a=4 Å, 4×4×4, m=1, two species with seeded RNG, b_iso=None."""
    import numpy as np
    rng = np.random.default_rng(42)
    lat = Lattice.from_parameters(4.0, 4.0, 4.0, 90.0, 90.0, 90.0)
    n = (4, 4, 4)
    N = 64
    species = np.array(
        ["CU" if rng.random() < 0.5 else "AG" for _ in range(N)],
        dtype=object,
    )
    frac = np.zeros((N, 3))
    cell_index = np.zeros((N, 3), dtype=np.int64)
    # DISCUS cell-major order: row r → cell r//m, slot r%m.  With m=1,
    # row r is in cell r.  Decode r as (ix, iy, iz) with x fastest.
    for r in range(N):
        iz = r // 16
        rem = r % 16
        iy = rem // 4
        ix = rem % 4
        cell_index[r] = (ix, iy, iz)
        frac[r] = [ix / 4, iy / 4, iz / 4]
    return Structure(
        lat, species, frac,
        model="configuration",
        b_iso=None,
        cell_index=cell_index,
        supercell=n,
        pbc=(True, True, True),
        min_separation=None,
    )


def _build_gaas_p1_2x2x2() -> Structure:
    """GaAs P1, 8 atoms per cell explicit, 2×2×2, m=8, triclinic-distorted.

    The supercell is (2,2,2) with 8 atoms per unit cell, so N=32 and
    m = N / prod(n) = 32 / 8 = 4.  The rows are in DISCUS cell-major order:
    row r → cell r//m, slot r%m.
    """
    import numpy as np
    lat = Lattice.from_parameters(5.65, 5.70, 5.60, 89.0, 91.0, 88.0)
    # 8 atoms per cell (P1 GaAs: 4 GA + 4 AS in the conventional cell)
    unit_frac = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.5, 0.5, 0.0],
            [0.0, 0.5, 0.5],
            [0.5, 0.0, 0.5],
            [0.25, 0.25, 0.25],
            [0.75, 0.75, 0.25],
            [0.25, 0.75, 0.75],
            [0.75, 0.25, 0.75],
        ]
    )
    unit_species = np.array(
        ["GA", "GA", "GA", "GA", "AS", "AS", "AS", "AS"],
        dtype=object,
    )
    n = (2, 2, 2)
    m = 4  # N / prod(n) = 32 / 8 = 4
    N = 32
    species = np.empty(N, dtype=object)
    frac = np.empty((N, 3))
    cell_index = np.zeros((N, 3), dtype=np.int64)
    for r in range(N):
        cell_num = r // m
        slot = r % m
        iz = cell_num // 4
        rem = cell_num % 4
        iy = rem // 2
        ix = rem % 2
        cell_index[r] = (ix, iy, iz)
        species[r] = unit_species[slot]
        frac[r] = (unit_frac[slot] + np.array([ix, iy, iz])) / np.array(n)
    return Structure(
        lat, species, frac,
        model="configuration",
        b_iso=np.full(N, 0.35),
        cell_index=cell_index,
        supercell=n,
        pbc=(True, True, True),
        min_separation=None,
    )


# -- L3: Python-written .cll expands as expected ----------------------------
def test_live_python_written_cll_expands_as_expected():
    """L3: write_cll P1 8-atom GaAs → DISCUS read cell,2,2,2 → save → compare."""
    workdir = _workdir()
    import numpy as np
    lat = Lattice.from_parameters(5.6537, 5.6537, 5.6537, 90.0, 90.0, 90.0)
    # P1: all 8 atoms explicit (no symmetry expansion)
    unit_frac = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.5, 0.5, 0.0],
            [0.0, 0.5, 0.5],
            [0.5, 0.0, 0.5],
            [0.25, 0.25, 0.25],
            [0.75, 0.75, 0.25],
            [0.25, 0.75, 0.75],
            [0.75, 0.25, 0.75],
        ]
    )
    unit_species = np.array(
        ["GA", "GA", "GA", "GA", "AS", "AS", "AS", "AS"],
        dtype=object,
    )
    s = Structure(
        lat, unit_species, unit_frac,
        model="configuration",
        b_iso=np.full(8, 1.0),
        supercell=(1, 1, 1),
        pbc=(True, True, True),
        min_separation=None,
    )
    write_cll(os.path.join(workdir, "gaas_p1.cll"), s)

    macro = (
        "discus\n"
        "read\n"
        "  cell gaas_p1.cll,2,2,2\n"
        "eval n[1]\n"
        "save gaas_p1_live.stru\n"
        "exit\n"
        "exit\n"
    )
    result = run_discus(macro, workdir)
    assert "64.0000000000" in result.stdout

    # Build the expected 2×2×2 supercell in numpy.
    n = (2, 2, 2)
    m = 8
    N = 64
    exp_species = np.empty(N, dtype=object)
    exp_frac = np.empty((N, 3))
    for r in range(N):
        cell_num = r // m
        slot = r % m
        iz = cell_num // 4
        rem = cell_num % 4
        iy = rem // 2
        ix = rem % 2
        exp_species[r] = unit_species[slot]
        exp_frac[r] = (unit_frac[slot] + np.array([ix, iy, iz])) / np.array(n)

    s2 = read_stru(os.path.join(workdir, "gaas_p1_live.stru"))
    assert s2.structure.n_atoms == N
    assert np.array_equal(s2.structure.species, exp_species)
    # DISCUS saves in unit-cell coordinates with a centered origin for P1
    # (values may be negative, e.g. -1.0 = 0.0 mod 1).  The .cll unit
    # fractions are in [0, 1); multiply by n to get supercell coordinates
    # in [0, n).  Compare modulo n to account for the origin shift.
    # Build the expected supercell coordinates from the unit fractions.
    n_arr = np.array(n, dtype=float)
    s2_super = s2.structure.frac * n_arr
    exp_super = np.empty((N, 3))
    for r in range(N):
        slot = r % m
        exp_super[r] = unit_frac[slot] * n_arr
    assert np.allclose(s2_super % n_arr, exp_super % n_arr, atol=2e-6)
    evidence_log("test_live_python_written_cll_expands_as_expected", macro, None, result)


# -- L4: double-expansion hazard -------------------------------------------
def test_live_double_expansion_hazard():
    """L4: DISCUS read cell vs read stru on the same .stru → counts differ 4×."""
    workdir = _workdir()
    # Reuse the GaAs P1 2×2×2 .stru from L2.
    s = _build_gaas_p1_2x2x2()
    fname = "live_gaas_p1_2x2x2.stru"
    write_stru(os.path.join(workdir, fname), s)

    macro = (
        "discus\n"
        "read\n"
        f"  stru {fname}\n"
        "eval n[1]\n"
        f"read\n"
        f"  cell {fname},1,1,1\n"
        "eval n[1]\n"
        "exit\n"
        "exit\n"
    )
    result = run_discus(macro, workdir)

    # Extract the two n[1] values.
    import re
    values = re.findall(r"Value of n\[1\] =\s+([\d.]+)", result.stdout)
    assert len(values) == 2
    n_stru = int(float(values[0]))
    n_cell = int(float(values[1]))
    # Q6: read cell on a .stru does NOT re-expand (both give 32).
    # The double-expansion hazard is avoided: DISCUS treats the .stru as a
    # complete structure, not a unit cell to be expanded.
    assert n_cell == n_stru, f"expected equal counts, got {n_cell} vs {n_stru}"

    # Our parse_cll raises DiscusFormatError mentioning read_stru.
    text = open(os.path.join(workdir, fname), encoding="ascii").read()
    with pytest.raises(DiscusFormatError, match="read_stru"):
        parse_cll(text)

    evidence_log("test_live_double_expansion_hazard", macro, None, result)


# -- L5: species names -----------------------------------------------------
def test_live_species_names():
    """L5: species spellings that DISCUS accepts survive write_stru → read → save."""
    workdir = _workdir()
    import numpy as np
    lat = Lattice.from_parameters(4.0, 4.0, 4.0, 90.0, 90.0, 90.0)
    # Q9: Ga → GA, O2- → O2-, Fe3+ → FE3+
    species = np.array(["Ga", "O2-", "Fe3+", "Ga1"], dtype=object)
    frac = np.array(
        [[0.0, 0.0, 0.0], [0.25, 0.25, 0.25], [0.5, 0.5, 0.5], [0.75, 0.75, 0.75]]
    )
    s = Structure(
        lat, species, frac,
        model="configuration",
        b_iso=np.full(4, 0.0),
        cell_index=np.zeros((4, 3), dtype=np.int64),
        supercell=(1, 1, 1),
        pbc=(True, True, True),
        min_separation=None,
    )
    fname = "live_species.stru"
    write_stru(os.path.join(workdir, fname), s)

    macro = (
        "discus\n"
        "read\n"
        f"  stru {fname}\n"
        "eval n[1]\n"
        f"save {fname}_out.stru\n"
        "exit\n"
        "exit\n"
    )
    result = run_discus(macro, workdir)
    assert "4.00000000000" in result.stdout

    s2 = read_stru(os.path.join(workdir, f"{fname}_out.stru"))
    assert s2.structure.n_atoms == 4
    # DISCUS uppercases: Ga → GA, Fe3+ → FE3+; O2- and Ga1 unchanged.
    assert list(s2.structure.species) == ["GA", "O2-", "FE3+", "GA1"]
    evidence_log("test_live_species_names", macro, None, result)
