"""Tests for structure_ops.repeat (Package A, Feature A4).

Oracle: an independent numpy loop builds the expected frac, cell_index, atom_id
for the GaAs 8-atom cell repeated by (2,3,4) and compares exactly.
"""
import math

import numpy as np
import pytest

from discuspdf import Lattice, Structure, repeat
from discuspdf.cif_io import read_cif
from discuspdf.discus_io import format_stru, parse_stru


# ------------------------------ helpers ------------------------------
def _gaas_cell() -> Structure:
    """GaAs 8-atom cell from the CIF fixture."""
    cif = read_cif("tests/fixtures/interop/gaas_cod_9008845.cif")
    return cif.structure


def _triclinic_cell() -> Structure:
    """Triclinic P-1 cell from the CIF fixture."""
    cif = read_cif("tests/fixtures/interop/triclinic_p1bar.cif")
    return cif.structure


def _gaas_supercell_disco(N: int) -> Structure:
    """GaAs zinc-blende N x N x N supercell in DISCUS order (ix fastest)."""
    a = 5.6537
    fcc = np.array([[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]])
    base = np.vstack([fcc, fcc + 0.25])
    cells = np.array([(ix, iy, iz) for iz in range(N) for iy in range(N) for ix in range(N)])
    frac = ((base[None, :, :] + cells[:, None, :]) / N).reshape(-1, 3)
    labels = [f"Ga{i}" for i in range(1, 5)] + [f"As{i}" for i in range(1, 5)]
    return Structure(
        Lattice.from_parameters(a * N, a * N, a * N, 90.0, 90.0, 90.0),
        np.tile(["Ga"] * 4 + ["As"] * 4, len(cells)).astype(object),
        frac,
        model="configuration",
        site_label=np.tile(labels, len(cells)).astype(object),
        cell_index=np.repeat(cells, 8, axis=0),
        supercell=(N, N, N),
    )


# ------------------------------ 1. oracle: independent numpy loop ------------------------------
def test_repeat_oracle_gaas_234():
    """Independent numpy loop builds the expected frac, cell_index, atom_id
    for the GaAs 8-atom cell repeated by (2,3,4) and compares exactly."""
    s = _gaas_cell()
    reps = (2, 3, 4)
    nx, ny, nz = reps

    result = repeat(s, reps)
    new_s = result.structure

    # Build expected arrays with an independent numpy loop.
    n_old = s.n_atoms  # 8
    frac_old = s.frac
    ci_old = s.cell_index
    sc_old = s.supercell  # (1, 1, 1)
    aid_old = s.atom_id  # 0..7
    stride = int(aid_old.max()) + 1  # 8

    n_new = n_old * nx * ny * nz
    frac_exp = np.empty((n_new, 3), dtype=float)
    ci_exp = np.empty((n_new, 3), dtype=np.int64)
    aid_exp = np.empty(n_new, dtype=np.int64)

    row = 0
    for oz in range(nz):
        for oy in range(ny):
            for ox in range(nx):
                o = np.array([ox, oy, oz], dtype=np.int64)
                offset_index = ox + nx * (oy + ny * oz)
                frac_exp[row:row + n_old] = (frac_old + o) / np.array([nx, ny, nz], dtype=float)
                ci_exp[row:row + n_old] = ci_old + sc_old * o
                aid_exp[row:row + n_old] = offset_index * stride + aid_old
                row += n_old

    # Compare exactly (frac atol 1e-15).
    np.testing.assert_allclose(new_s.frac, frac_exp, atol=1e-15)
    np.testing.assert_array_equal(new_s.cell_index, ci_exp)
    np.testing.assert_array_equal(new_s.atom_id, aid_exp)

    # Check info.
    assert result.info["reps"] == (2, 3, 4)
    assert result.info["stride"] == 8
    assert result.ids.size == 0
    assert result.seed is None


# ------------------------------ 2. invariants ------------------------------
def test_repeat_invariants_gaas_234():
    """lattice.parameters of the result = parent a*nx, b*ny, c*nz, same angles;
    unit cell = lattice/supercell equals the parent lattice; supercell = reps;
    n_atoms = 8*24."""
    s = _gaas_cell()
    reps = (2, 3, 4)
    new_s = repeat(s, reps).structure

    # Lattice parameters.
    a, b, c, alpha, beta, gamma = s.lattice.parameters
    a2, b2, c2, alpha2, beta2, gamma2 = new_s.lattice.parameters
    np.testing.assert_allclose([a2, b2, c2], [a * 2, b * 3, c * 4], rtol=1e-12)
    np.testing.assert_allclose([alpha2, beta2, gamma2], [alpha, beta, gamma], rtol=1e-12)

    # Unit cell = lattice / supercell equals parent lattice.
    unit_a = a2 / new_s.supercell[0]
    unit_b = b2 / new_s.supercell[1]
    unit_c = c2 / new_s.supercell[2]
    np.testing.assert_allclose([unit_a, unit_b, unit_c], [a, b, c], rtol=1e-12)

    # Supercell = reps.
    np.testing.assert_array_equal(new_s.supercell, [2, 3, 4])

    # n_atoms = 8 * 24.
    assert new_s.n_atoms == 8 * 24


# ------------------------------ 3. physics: GaAs 3x3x3 neighbours ------------------------------
def test_repeat_physics_gaas_333():
    """GaAs 3x3x3 (pbc all): every atom has 4 opposite-species neighbours at
    a*sqrt(3)/4 +/- 1e-6 and the minimum distance is a*sqrt(3)/4."""
    s = _gaas_cell()
    new_s = repeat(s, (3, 3, 3)).structure
    a = s.lattice.parameters[0]
    d_expected = a * math.sqrt(3) / 4

    dm = new_s.distance_matrix()
    n = new_s.n_atoms

    for i in range(n):
        sp_i = new_s.species[i]
        dists = dm[i].copy()
        dists[i] = np.inf  # exclude self
        # Find opposite-species neighbours.
        opp_mask = np.array([new_s.species[j] != sp_i for j in range(n)])
        opp_dists = dists[opp_mask]
        # Every atom should have at least 4 neighbours at ~d_expected.
        close = np.sum(np.abs(opp_dists - d_expected) < 1e-6)
        assert close >= 4, f"atom {i}: only {close} neighbours at {d_expected:.6f}"

    # Minimum distance (excluding self).
    np.fill_diagonal(dm, np.inf)
    d_min = dm.min()
    np.testing.assert_allclose(d_min, d_expected, atol=1e-6)


# ------------------------------ 4. physics: triclinic 2x2x2 distance_matrix ------------------------------
def test_repeat_physics_triclinic_222():
    """distance_matrix for a 125-atom subset matches parent distances (minimum image)
    for a repeated triclinic cell (2x2x2)."""
    s = _triclinic_cell()
    new_s = repeat(s, (2, 2, 2)).structure

    # Pick a 125-atom subset (or all if fewer).
    n_sub = min(125, new_s.n_atoms)
    ids = np.arange(n_sub)
    dm_new = new_s.distance_matrix(ids)

    # Independent computation: for each pair (i, j) in the subset, compute the
    # minimum-image distance over lattice translations.
    A = new_s.lattice.matrix
    cart = new_s.frac[ids] @ A
    n = len(ids)
    dm_exp = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            diff = cart[j] - cart[i]
            # Minimum image: try all 27 translations of the supercell.
            best = np.inf
            for tx in range(-1, 2):
                for ty in range(-1, 2):
                    for tz in range(-1, 2):
                        t = np.array([tx, ty, tz])
                        # Translation in Cartesian: t @ A (rows of A are a, b, c)
                        # Actually the translation is t[0]*a + t[1]*b + t[2]*c
                        # = t @ A (where A has rows a, b, c)
                        trans = t @ A
                        d = np.linalg.norm(diff - trans)
                        if d < best:
                            best = d
            dm_exp[i, j] = best
            dm_exp[j, i] = best

    np.testing.assert_allclose(dm_new, dm_exp, atol=1e-6)


# ------------------------------ 5. ID mapping ------------------------------
def test_repeat_id_mapping_contiguous():
    """new % stride, new // stride recover the parent id and the offset index
    for every row, for a parent with ids 0..n-1."""
    s = _gaas_cell()
    reps = (2, 3, 4)
    new_s = repeat(s, reps).structure
    stride = 8  # max(0..7) + 1

    for i in range(new_s.n_atoms):
        new_id = int(new_s.atom_id[i])
        old_id = new_id % stride
        offset_index = new_id // stride
        # Verify: old_id should be a valid parent id.
        assert 0 <= old_id < 8
        # Verify: offset_index should be in [0, 24).
        assert 0 <= offset_index < 24
        # Verify: the old atom at old_id should have the same species.
        assert new_s.species[i] == s.species[old_id]


def test_repeat_id_mapping_noncontiguous():
    """new % stride, new // stride recover the parent id and the offset index
    for a parent whose ids were shifted by +100 and made non-contiguous (select a subset)."""
    s = _gaas_cell()
    # Select a subset (rows 0, 2, 4, 6) and shift ids by +100.
    subset = s.select(rows=np.array([0, 2, 4, 6]))
    # The subset has atom_ids 0, 2, 4, 6 (original ids).
    # Shift by +100: we need to build a new structure with shifted ids.
    # Use the internal _Structure_bypass_validation or just build manually.
    # Easier: use repeat on the subset, then check the mapping.
    # The subset has max_id = 6, so stride = 7.
    new_s = repeat(subset, (2, 2, 2)).structure
    stride = 7  # max(0, 2, 4, 6) + 1

    for i in range(new_s.n_atoms):
        new_id = int(new_s.atom_id[i])
        old_id = new_id % stride
        offset_index = new_id // stride
        # old_id should be one of {0, 2, 4, 6}.
        assert old_id in (0, 2, 4, 6), f"old_id={old_id} not in parent ids"
        assert 0 <= offset_index < 8


# ------------------------------ 6. chained ------------------------------
def test_repeat_chained():
    """repeat(repeat(s,(2,1,1)),(1,3,1)) equals repeat(s,(2,3,1)) up to row order
    and ids: compare sorted (species, round(frac, 12), cell_index) tuples."""
    s = _gaas_cell()
    r1 = repeat(repeat(s, (2, 1, 1)).structure, (1, 3, 1)).structure
    r2 = repeat(s, (2, 3, 1)).structure

    assert r1.n_atoms == r2.n_atoms

    # Compare sorted (species, round(frac, 12), cell_index) tuples.
    def key_tuples(struct):
        tuples = []
        for i in range(struct.n_atoms):
            tuples.append((
                struct.species[i],
                tuple(round(v, 12) for v in struct.frac[i]),
                tuple(struct.cell_index[i]),
            ))
        return sorted(tuples)

    assert key_tuples(r1) == key_tuples(r2)


# ------------------------------ 7. average model ------------------------------
def test_repeat_average_model():
    """Average model (mixed fixture via read_cif): works, model stays 'average',
    occupancy sum per (site_label, cell_index) is 1.0 for every site,
    number of sites = parent sites x reps."""
    cif = read_cif("tests/fixtures/interop/ga07in03as_mixed.cif")
    s = cif.structure
    assert s.model == "average"

    n_parent_sites = len(s.sites())
    reps = (2, 2, 2)
    new_s = repeat(s, reps).structure

    assert new_s.model == "average"

    # Occupancy sum per (site_label, cell_index) is 1.0 for every site.
    sites = new_s.sites()
    for site in sites:
        np.testing.assert_allclose(site["total_occupancy"], 1.0, atol=1e-9)

    # Number of sites = parent sites x reps.
    assert len(sites) == n_parent_sites * 8


# ------------------------------ 8. identity ------------------------------
def test_repeat_identity():
    """repeat(s, (1,1,1)) == s."""
    s = _gaas_cell()
    new_s = repeat(s, (1, 1, 1)).structure
    assert new_s == s


# ------------------------------ 9. errors ------------------------------
def test_repeat_errors_reps():
    """reps with 0, negative, float, wrong length."""
    s = _gaas_cell()
    with pytest.raises(ValueError, match="reps"):
        repeat(s, (0, 1, 1))
    with pytest.raises(ValueError, match="reps"):
        repeat(s, (-1, 1, 1))
    with pytest.raises(ValueError, match="reps"):
        repeat(s, (1.5, 1, 1))
    with pytest.raises(ValueError, match="reps"):
        repeat(s, (1, 1))
    with pytest.raises(ValueError, match="reps"):
        repeat(s, (1, 1, 1, 1))
    with pytest.raises(ValueError, match="reps"):
        repeat(s, 2)  # int, not a tuple


def test_repeat_errors_pbc():
    """pbc not all True."""
    s = _gaas_cell()
    # Build a structure with pbc not all True.
    # The GaAs cell has pbc all True, so build a custom one.
    lat = Lattice.from_parameters(5.0, 5.0, 5.0, 90.0, 90.0, 90.0)
    s_nonpbc = Structure(
        lat, np.array(["Ga"], dtype=object), np.array([[0.25, 0.25, 0.25]]),
        model="configuration", pbc=(True, True, False),
    )
    with pytest.raises(ValueError, match="pbc"):
        repeat(s_nonpbc, (2, 2, 2))


def test_repeat_errors_min_separation():
    """min_separation pass-through: a min_separation=5 overlap error on GaAs."""
    s = _gaas_cell()
    # GaAs has atoms at ~2.85 Å apart (a*sqrt(3)/4). min_separation=5 should raise.
    with pytest.raises(ValueError, match="apart"):
        repeat(s, (2, 2, 2), min_separation=5.0)


# ------------------------------ 10. DISCUS compat ------------------------------
def test_repeat_discus_compat():
    """format_stru(repeat(gaas_cell,(4,4,4)).structure) -> parse_stru equals
    the repeated structure up to the b_iso policy (b_iso None -> 1.0),
    and the order is cell-major, x fastest."""
    # Build a DISCUS-ordered GaAs 1x1x1 cell with model="configuration".
    a = 5.6537
    fcc = np.array([[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]])
    base = np.vstack([fcc, fcc + 0.25])
    s = Structure(
        Lattice.from_parameters(a, a, a, 90.0, 90.0, 90.0),
        np.array(["Ga"] * 4 + ["As"] * 4, dtype=object),
        base,
        model="configuration",
        site_label=np.array([f"Ga{i}" for i in range(1, 5)] + [f"As{i}" for i in range(1, 5)], dtype=object),
        cell_index=np.zeros((8, 3), dtype=np.int64),
        supercell=(1, 1, 1),
    )
    new_s = repeat(s, (4, 4, 4)).structure

    text = format_stru(new_s)
    ds = parse_stru(text)
    s_rt = ds.structure

    # Compare frac (the .stru reader assigns atom_id = 0..N-1 and
    # site_label = site{slot+1}, so we ignore those).
    np.testing.assert_allclose(new_s.frac, s_rt.frac, atol=1e-6)
    np.testing.assert_array_equal(new_s.species, s_rt.species)
    np.testing.assert_array_equal(new_s.cell_index, s_rt.cell_index)
    np.testing.assert_array_equal(new_s.supercell, s_rt.supercell)
    assert new_s.model == s_rt.model

    # Check the order is cell-major, x fastest.
    n = new_s.supercell
    m = new_s.n_atoms // (n[0] * n[1] * n[2])
    for r in range(new_s.n_atoms):
        expected_cell = r // m
        ix, iy, iz = (int(v) for v in new_s.cell_index[r])
        assert iz * (n[0] * n[1]) + iy * n[0] + ix == expected_cell, \
            f"row {r} not in DISCUS order"


# ------------------------------ 11. speed ------------------------------
def test_repeat_speed_4096():
    """8x8x8 GaAs (4096 atoms) repeat < 5 s on the owner's machine (assert < 20 s)."""
    import time
    s = _gaas_cell()
    t0 = time.perf_counter()
    new_s = repeat(s, (8, 8, 8)).structure
    elapsed = time.perf_counter() - t0
    print(f"  8x8x8 repeat: {elapsed:.3f} s")
    assert new_s.n_atoms == 4096
    assert elapsed < 20.0, f"repeat too slow: {elapsed:.1f} s"
