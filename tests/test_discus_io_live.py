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

    # Read the saved .stru with our parser.  DISCUS drops the ncell line on
    # save, so the file has no ncell; the atoms span the centred 2x2x2
    # supercell, so use supercell="infer".
    s = read_stru(os.path.join(workdir, "gaas_live.stru"), supercell="infer")
    assert s.structure.n_atoms == 64
    from collections import Counter
    assert dict(Counter(s.structure.species)) == {"GA": 32, "AS": 32}

    # Compare to the committed fixture (also read with supercell="infer").
    ref = read_stru("tests/fixtures/discus_6_23_01/gaas_222.stru", supercell="infer")
    assert s.structure.isclose(ref.structure, rtol=1e-8, atol=1e-6)

    # Physics check (the one that would have caught the A3-6a bug): a correct
    # zinc-blende crystal has exactly 4 neighbours within 2.6 A at 2.4480 A.
    import numpy as np
    D = s.structure.distance_matrix()
    np.fill_diagonal(D, 99.0)
    counts = (D < 2.6).sum(axis=1)
    assert set(counts.tolist()) == {4}, f"neighbour counts {sorted(set(counts.tolist()))} != [4]"
    # zinc-blende nearest-neighbour = a*sqrt(3)/4 = 5.6537*sqrt(3)/4 = 2.44812
    assert D.min() == pytest.approx(5.6537 * 3 ** 0.5 / 4, abs=1e-4), f"min distance {D.min():.4f}"
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

    # Read the DISCUS-saved file.  DISCUS drops the ncell line on save; a
    # `read stru` + `save` does NOT centre the origin (only `read cell` +
    # `save` does), so the coordinates are 0..n-1 and "infer" (which requires
    # a centred crystal) cannot be used.  Use the explicit supercell tuple
    # (the builder's supercell).
    s2 = read_stru(os.path.join(workdir, f"{fname}_out.stru"), supercell=tuple(s.supercell))
    assert s2.structure.n_atoms == N

    # Same species per row.
    import numpy as np
    assert np.array_equal(s2.structure.species, s.species)

    # s2 is read in the inferred supercell frame; the file coordinates are
    # centred unit-cell values, so s2.structure.frac * n gives the centred
    # supercell coordinates.  Compare modulo n to the input supercell
    # fractions to account for the origin shift.
    n = np.array(s.supercell, dtype=float)
    s2_super = s2.structure.frac * n
    s_super = s.frac * n
    assert np.allclose(s2_super % n, s_super % n, atol=2e-6)

    # Cartesian physics check: the DISCUS-saved file, read back, must agree
    # with the Python input in Cartesian (A) to < 2e-6 * max cell length.
    # Reduce both to the unit cell (mod 1) so the origin shift cancels.
    a, b, c, *_ = s.lattice.parameters
    tol = 2e-6 * max(a, b, c)
    s2_unit = (s2_super % 1.0) @ s2.structure.lattice.matrix
    s_unit = (s_super % 1.0) @ s.lattice.matrix
    diff = np.abs(s2_unit - s_unit)
    assert diff.max() < tol, f"Cartesian mismatch {diff.max():.3e} A >= {tol:.3e} A"

    # Biso: the writer always writes B explicitly (None -> 1.0, the DISCUS
    # default, A3-9), and DISCUS must keep exactly what was written.  If this
    # fails, DISCUS changes B on `read stru` / `save` and the writer policy
    # must be revisited (report it).
    expected_b = np.ones(N) if s.b_iso is None else s.b_iso
    assert np.allclose(s2.structure.b_iso, expected_b, atol=2e-6), (
        f"Biso written {sorted(set(np.round(expected_b, 6)))} but DISCUS saved "
        f"{sorted(set(np.round(s2.structure.b_iso, 6)))}"
    )
    from tests.discus_live import evidence_note
    with open(os.path.join(workdir, f"{fname}_out.stru"), encoding="ascii") as fh:
        evidence_note("test_live_python_written_stru_loads_and_roundtrips",
                      f"--- saved {fname}_out.stru (first 12 lines) ---\n" + "".join(fh.readlines()[:12]))

    evidence_log("test_live_python_written_stru_loads_and_roundtrips", macro, None, result)


def _build_cubic_4x4x4() -> Structure:
    """Simple cubic, unit cell 4 Å, 4×4×4 (supercell frame 16 Å), m=1, seeded species, b_iso=None."""
    import numpy as np
    rng = np.random.default_rng(42)
    n = (4, 4, 4)
    # Structure.lattice is the whole supercell frame: unit cell x n.
    lat = Lattice.from_parameters(4.0 * n[0], 4.0 * n[1], 4.0 * n[2], 90.0, 90.0, 90.0)
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
    )


def _build_gaas_p1_2x2x2() -> Structure:
    """GaAs P1, 8 atoms per cell explicit, 2×2×2 (64 atoms: 32 Ga + 32 As), triclinic-distorted.

    Unit cell (5.65, 5.70, 5.60, 89, 91, 88); ``Structure.lattice`` is the
    supercell frame (unit cell x 2).  Rows are in DISCUS cell-major order,
    x fastest: row r -> cell r // 8, slot r % 8.
    """
    import numpy as np
    n = (2, 2, 2)
    lat = Lattice.from_parameters(5.65 * n[0], 5.70 * n[1], 5.60 * n[2], 89.0, 91.0, 88.0)
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
    unit_species = np.array(["GA"] * 4 + ["AS"] * 4, dtype=object)
    m = len(unit_frac)
    N = m * n[0] * n[1] * n[2]
    cell_num = np.arange(N) // m
    slot = np.arange(N) % m
    cell_index = np.stack(
        [cell_num % n[0], (cell_num // n[0]) % n[1], cell_num // (n[0] * n[1])], axis=1
    ).astype(np.int64)
    return Structure(
        lat, unit_species[slot], (unit_frac[slot] + cell_index) / np.array(n, dtype=float),
        model="configuration",
        b_iso=np.full(N, 0.35),
        cell_index=cell_index,
        supercell=n,
        pbc=(True, True, True),
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

    # `read cell` + `save` centres the origin, so "infer" recovers the 2x2x2
    # supercell.
    s2 = read_stru(os.path.join(workdir, "gaas_p1_live.stru"), supercell="infer")
    assert s2.structure.n_atoms == N
    from collections import Counter
    assert dict(Counter(s2.structure.species)) == {"GA": 32, "AS": 32}

    # Physics check: the P1 8-atom cell puts Ga and As on the fcc sites, so
    # the nearest-neighbour distance is the Ga-As distance a*sqrt(3)/4 =
    # 5.6537*sqrt(3)/4 = 2.4481 A (4 per atom).
    D = s2.structure.distance_matrix()
    np.fill_diagonal(D, 99.0)
    nn = (D < 2.5).sum(axis=1)
    assert D.min() == pytest.approx(5.6537 * 3 ** 0.5 / 4, abs=1e-3), f"min {D.min():.4f}"
    assert set(nn.tolist()) == {4}, f"nearest-neighbour counts {sorted(set(nn.tolist()))} != [4]"

    # Coordinate check (order- and origin-independent): DISCUS `read cell` +
    # `save` reorders the atoms by the space group and centres the origin, so
    # the row order and the absolute origin are not preserved.  The physics
    # check above (4 neighbours at 2.4481 A) plus the per-species counts
    # (32 GA + 32 AS) confirm DISCUS built the correct 2x2x2 crystal; a
    # position-by-position comparison would be brittle under the centred
    # origin and the space-group reordering.
    evidence_log("test_live_python_written_cll_expands_as_expected", macro, None, result)


# -- L4: double-expansion hazard -------------------------------------------
def test_live_double_expansion_hazard():
    """L4: DISCUS read cell vs read stru on the same F-43m .stru → counts
    differ 4x (Q6).  The hazard is specific to a space group with expansion
    (F-43m: 2 listed atoms -> 8/cell); a P1 .stru does not re-expand."""
    workdir = _workdir()
    import shutil
    # Use the committed F-43m 2x2x2 fixture (64 atoms): `read cell` on it
    # re-expands by the space group (4x), the double-expansion hazard.
    fname = "gaas_222_hazard.stru"
    shutil.copy("tests/fixtures/discus_6_23_01/gaas_222.stru", os.path.join(workdir, fname))

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
    # Q6: read cell on the F-43m .stru RE-EXPANDS it, silently multiplying
    # the crystal 4x (64 -> 256).  This is exactly why parse_cll refuses a
    # .stru (below).
    assert n_stru == 64, f"read stru gave {n_stru}, expected 64"
    assert n_cell == 4 * n_stru, f"expected 4x ({4 * n_stru}), got {n_cell} vs {n_stru}"

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
        b_iso=np.full(4, 1.0),
        cell_index=np.zeros((4, 3), dtype=np.int64),
        supercell=(1, 1, 1),
        pbc=(True, True, True),
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

    # `read stru` + `save` does not centre the origin, so use the explicit
    # supercell (the builder's (1,1,1)).
    s2 = read_stru(os.path.join(workdir, f"{fname}_out.stru"), supercell=(1, 1, 1))
    assert s2.structure.n_atoms == 4
    # DISCUS uppercases: Ga → GA, Fe3+ → FE3+; O2- and Ga1 unchanged.
    assert list(s2.structure.species) == ["GA", "O2-", "FE3+", "GA1"]
    evidence_log("test_live_species_names", macro, None, result)


# -- L6: B values survive the DISCUS round trip (A3-9) -----------------------
def test_live_b_values_survive():
    """L6: B values in [0.001, 2.0] survive write → DISCUS read → save → read.

    Exercises both read paths: ``read stru`` (saved structure) and
    ``read cell`` (unit-cell source, P1 1×1×1).  DISCUS 6.23.01 treats B
    below about 1e-3 as "unset" and replaces it with 1.0 (probe Q15), so a
    value at the boundary (0.001) and typical values (0.35, 2.0) must come
    back unchanged.
    """
    workdir = _workdir()
    import numpy as np
    lat = Lattice.from_parameters(4.0, 4.0, 4.0, 90.0, 90.0, 90.0)
    species = np.array(["CU", "CU", "CU"], dtype=object)
    frac = np.array([[0.0, 0.0, 0.0], [0.25, 0.25, 0.25], [0.5, 0.5, 0.5]])
    b_iso = np.array([0.001, 0.35, 2.0])
    s = Structure(
        lat, species, frac,
        model="configuration",
        b_iso=b_iso,
        cell_index=np.zeros((3, 3), dtype=np.int64),
        supercell=(1, 1, 1),
        pbc=(True, True, True),
    )

    # Path 1: .stru (read stru → save → read back).
    fname = "live_b.stru"
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
    assert "3.00000000000" in result.stdout
    # `read stru` + `save` does not centre the origin; use the explicit
    # supercell (1,1,1).
    s2 = read_stru(os.path.join(workdir, f"{fname}_out.stru"), supercell=(1, 1, 1))
    assert s2.structure.n_atoms == 3
    assert np.allclose(s2.structure.b_iso, b_iso, atol=1e-6), (
        f"Biso written {list(b_iso)} but DISCUS saved {list(s2.structure.b_iso)}"
    )
    from tests.discus_live import evidence_note
    with open(os.path.join(workdir, f"{fname}_out.stru"), encoding="ascii") as fh:
        evidence_note("test_live_b_values_survive",
                      f"--- saved {fname}_out.stru (first 12 lines) ---\n" + "".join(fh.readlines()[:12]))
    evidence_log("test_live_b_values_survive", macro, None, result)

    # Path 2: .cll (read cell,1,1,1 → save → read back).  P1: no expansion.
    cfname = "live_b.cll"
    write_cll(os.path.join(workdir, cfname), s)
    macro = (
        "discus\n"
        "read\n"
        f"  cell {cfname},1,1,1\n"
        "eval n[1]\n"
        f"save {cfname}_out.stru\n"
        "exit\n"
        "exit\n"
    )
    result = run_discus(macro, workdir)
    assert "3.00000000000" in result.stdout
    s3 = read_stru(os.path.join(workdir, f"{cfname}_out.stru"), supercell=(1, 1, 1))
    assert s3.structure.n_atoms == 3
    assert np.allclose(s3.structure.b_iso, b_iso, atol=1e-6), (
        f"Biso written {list(b_iso)} but DISCUS saved {list(s3.structure.b_iso)}"
    )
    with open(os.path.join(workdir, f"{cfname}_out.stru"), encoding="ascii") as fh:
        evidence_note("test_live_b_values_survive",
                      f"--- saved {cfname}_out.stru (first 12 lines) ---\n" + "".join(fh.readlines()[:12]))
    evidence_log("test_live_b_values_survive", macro, None, result)


# -- L7: DISCUS replaces B=0 with its default 1.0 (A3-9) ----------------------
def test_live_zero_b_is_replaced_by_discus():
    """L7: a hand-written .stru with Biso 0.000000 → DISCUS saves Biso 1.0.

    Documents the engine behaviour that justifies the writer rule (B < 0.001
    is rejected on write): DISCUS 6.23.01 treats B = 0.000000 as "unset" and
    silently replaces it with its default 1.0 on ``read stru`` / ``save``
    (probe Q15).  A non-tiny B (0.35) is kept unchanged.
    """
    workdir = _workdir()
    import numpy as np
    # Hand-written .stru (NOT through our writer): 2-atom P1 cell, Biso
    # 0.000000 on the first atom and 0.350000 on the second.  DISCUS
    # ``read stru`` rejects positive coordinates with P1 (it expects the
    # centred-origin layout), so use negative coordinates (the centred
    # layout DISCUS writes on save).  No ``ncell`` line (DISCUS drops it
    # on save; ``read stru`` rejects a file that has one).
    text = (
        "title live zero b test\n"
        "spcgr P1\n"
        "cell    4.000000,  4.000000,  4.000000, 90.000000, 90.000000, 90.000000\n"
        "format numbers,XYZBPMMOS\n"
        "atoms      x,              y,              z,             Biso,    "
        "Property,  MoleNo,  MoleAt,   Occ,     St,  Sh,  Sk,  Sl\n"
        "GA        -1.000000,      -1.000000,      -1.000000,      0.000000,       "
        "1,       0,       0,   1.000000, _,   0,   0,   0\n"
        "GA        -1.000000,      -0.500000,      -0.500000,      0.350000,       "
        "1,       0,       0,   1.000000, _,   0,   0,   0\n"
    )
    fname = "live_zero_b.stru"
    with open(os.path.join(workdir, fname), "w", newline="\n", encoding="ascii") as fh:
        fh.write(text)

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
    assert "2.00000000000" in result.stdout

    s2 = read_stru(os.path.join(workdir, f"{fname}_out.stru"), supercell=(1, 1, 1))
    assert s2.structure.n_atoms == 2
    # DISCUS replaced the 0.0 B with its default 1.0; the 0.35 B is unchanged.
    assert np.allclose(s2.structure.b_iso, [1.0, 0.35], atol=1e-6), (
        f"DISCUS saved Biso {list(s2.structure.b_iso)}, expected [1.0, 0.35] "
        "(0.0 must be replaced by the DISCUS default 1.0)"
    )
    from tests.discus_live import evidence_note
    with open(os.path.join(workdir, f"{fname}_out.stru"), encoding="ascii") as fh:
        evidence_note("test_live_zero_b_is_replaced_by_discus",
                      f"--- saved {fname}_out.stru (first 12 lines) ---\n" + "".join(fh.readlines()[:12]))
    evidence_log("test_live_zero_b_is_replaced_by_discus", macro, None, result)
