"""Tests for discuspdf.neighbors (Package C, Feature 12).

Oracles:
* an independent brute-force enumeration (every node pair x every image shift in
  a box that is provably large enough) gives the exact expected pair set for
  cubic, tetragonal, triclinic, small and wrapped/unwrapped cells;
* textbook shell distances and multiplicities for sc, bcc, fcc and a primitive
  tetragonal lattice;
* DISCUS ``set vec`` neighbour arithmetic (cell index + translation, both
  directions) reproduced by index arithmetic on the DISCUS atom order.
"""
from __future__ import annotations

import csv
import itertools
import math
import time

import numpy as np
import pytest

from discuspdf import Lattice, NeighborList, Structure, neighbor_list


# ------------------------------ helpers ------------------------------
def _config(lattice, species, frac, **kw):
    kw.setdefault("model", "configuration")
    return Structure(lattice, np.asarray(species, dtype=object), np.asarray(frac, float), **kw)


def _cubic(kind: str, a: float = 4.0) -> Structure:
    base = {
        "sc": [[0, 0, 0]],
        "bcc": [[0, 0, 0], [0.5, 0.5, 0.5]],
        "fcc": [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]],
    }[kind]
    n = len(base)
    return _config(Lattice.from_parameters(a, a, a, 90, 90, 90), ["Cu"] * n, base,
                   site_label=[f"S{k + 1}" for k in range(n)])


def _brute_force(structure: Structure, cutoff: float) -> dict:
    """{(i, j, s0, s1, s2): distance} by enumerating every shift in a safe box (independent oracle)."""
    frac = structure.frac
    A = structure.lattice.matrix
    n = len(frac)
    out = {}
    if all(structure.pbc):
        B = np.linalg.inv(A).T
        # |Δf_k| <= cutoff |b_k| for any vector of length <= cutoff; positions span < 1 + their range
        span = frac.max(axis=0) - frac.min(axis=0)
        K = [int(math.ceil(cutoff * np.linalg.norm(B[k]) + span[k])) + 1 for k in range(3)]
        shifts = list(itertools.product(*[range(-K[k], K[k] + 1) for k in range(3)]))
    else:
        shifts = [(0, 0, 0)]
    for i in range(n):
        for j in range(n):
            for s in shifts:
                if i == j and s == (0, 0, 0):
                    continue
                d = float(np.linalg.norm((frac[j] + np.array(s) - frac[i]) @ A))
                if d <= cutoff * (1 + 1e-12):
                    out[(i, j) + tuple(s)] = d
    return out


def _as_dict(nl: NeighborList) -> dict:
    return {(int(i), int(j)) + tuple(int(x) for x in s): float(d)
            for i, j, s, d in zip(nl.i, nl.j, nl.shift, nl.distance)}


def _assert_matches_brute_force(structure, cutoff):
    nl = neighbor_list(structure, cutoff)
    got = _as_dict(nl)
    want = _brute_force(structure, cutoff)
    assert set(got) == set(want), (len(got), len(want), sorted(set(want) ^ set(got))[:5])
    for k, d in want.items():
        assert got[k] == pytest.approx(d, rel=1e-12, abs=1e-12)
    # stored vectors are consistent with the definition
    vec = (structure.frac[nl.j] + nl.shift - structure.frac[nl.i]) @ structure.lattice.matrix
    np.testing.assert_allclose(nl.vector, vec, atol=1e-12)
    np.testing.assert_allclose(nl.distance, np.linalg.norm(vec, axis=1), atol=1e-12)
    return nl


def _shell_table(nl, tol=1e-6):
    return [(round(s.distance, 6), s.coordination) for s in nl.shells(tol)]


# ------------------------------ brute-force equality ------------------------------
@pytest.mark.parametrize("cutoff", [2.5, 4.0, 6.1])
def test_fcc_supercell_matches_brute_force(cutoff):
    _assert_matches_brute_force(_cubic("fcc").repeat((2, 2, 2)), cutoff)


def test_tetragonal_random_positions_match_brute_force():
    rng = np.random.default_rng(12)
    lat = Lattice.from_parameters(3.1, 3.1, 4.7, 90, 90, 90)
    frac = rng.random((7, 3))
    s = _config(lat, ["Fe"] * 7, frac, min_separation=None)
    _assert_matches_brute_force(s, 5.5)


@pytest.mark.parametrize("params", [
    (4.1, 5.3, 6.2, 73.0, 81.0, 104.0),
    (3.0, 3.4, 9.0, 60.0, 115.0, 64.0),    # strongly skewed: naive per-axis rounding fails here
])
def test_triclinic_random_positions_match_brute_force(params):
    rng = np.random.default_rng(3)
    lat = Lattice.from_parameters(*params)
    frac = rng.random((6, 3))
    s = _config(lat, ["O"] * 6, frac, min_separation=None)
    _assert_matches_brute_force(s, 7.0)


def test_cutoff_larger_than_cell_lists_all_images_including_self():
    lat = Lattice.from_parameters(2.0, 2.3, 2.9, 88, 95, 101)
    s = _config(lat, ["X", "Y"], [[0.1, 0.2, 0.3], [0.6, 0.55, 0.8]])
    nl = _assert_matches_brute_force(s, 5.0)
    self_pairs = nl.i == nl.j
    assert self_pairs.any() and np.all(np.any(nl.shift[self_pairs] != 0, axis=1))


def test_unwrapped_coordinates_give_same_pairs_up_to_shift():
    """Coordinates outside [0, 1) (DISCUS centred crystals) are handled; shift refers to the stored frac."""
    rng = np.random.default_rng(5)
    lat = Lattice.from_parameters(5.0, 6.0, 7.0, 80, 95, 110)
    frac = rng.random((8, 3))
    offset = rng.integers(-3, 4, size=(8, 3))
    s0 = _config(lat, ["A"] * 8, frac, min_separation=None)
    s1 = _config(lat, ["A"] * 8, frac + offset, min_separation=None)
    nl0 = _assert_matches_brute_force(s0, 6.5)
    nl1 = _assert_matches_brute_force(s1, 6.5)
    np.testing.assert_allclose(np.sort(nl0.distance), np.sort(nl1.distance), atol=1e-12)


def test_open_boundaries_match_brute_force():
    rng = np.random.default_rng(9)
    lat = Lattice.from_parameters(10, 11, 12, 90, 90, 90)
    s = _config(lat, ["C"] * 40, rng.random((40, 3)), pbc=(False, False, False), min_separation=None)
    nl = _assert_matches_brute_force(s, 3.3)
    assert np.all(nl.shift == 0)


# ------------------------------ textbook shells ------------------------------
@pytest.mark.parametrize("kind, expected", [
    ("sc", [(4.0, 6), (4 * math.sqrt(2), 12), (4 * math.sqrt(3), 8), (8.0, 6)]),
    ("bcc", [(2 * math.sqrt(3), 8), (4.0, 6), (4 * math.sqrt(2), 12), (2 * math.sqrt(11), 24), (4 * math.sqrt(3), 8)]),
    ("fcc", [(2 * math.sqrt(2), 12), (4.0, 6), (2 * math.sqrt(6), 24), (4 * math.sqrt(2), 12), (2 * math.sqrt(10), 24)]),
])
@pytest.mark.parametrize("reps", [(1, 1, 1), (3, 3, 3)])
def test_cubic_shell_distances_and_multiplicities(kind, expected, reps):
    s = _cubic(kind).repeat(reps)
    nl = neighbor_list(s, max(d for d, _ in expected) + 0.01)
    got = _shell_table(nl)
    assert len(got) == len(expected)
    for (d, c), (de, ce) in zip(got, expected):
        assert d == pytest.approx(de, abs=1e-6)
        assert c == ce
    assert np.all(nl.coordination() == sum(c for _, c in expected))


def test_tetragonal_shells():
    s = _config(Lattice.from_parameters(3, 3, 4, 90, 90, 90), ["Fe"], [[0, 0, 0]])
    got = _shell_table(neighbor_list(s, 5.01))
    assert got == [(3.0, 4), (4.0, 2), (round(3 * math.sqrt(2), 6), 4), (5.0, 8)]


def test_rotated_lattice_gives_same_distances():
    rng = np.random.default_rng(1)
    lat = Lattice.from_parameters(4.2, 5.1, 6.3, 77, 85, 99)
    q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    q *= np.sign(np.linalg.det(q))
    lat_rot = Lattice.from_matrix(lat.matrix @ q.T)
    frac = rng.random((5, 3))
    a = neighbor_list(_config(lat, ["A"] * 5, frac, min_separation=None), 6.0)
    b = neighbor_list(_config(lat_rot, ["A"] * 5, frac, min_separation=None), 6.0)
    assert _as_dict(a).keys() == _as_dict(b).keys()
    np.testing.assert_allclose(a.distance, b.distance, atol=1e-12)


# ------------------------------ translations and DISCUS vectors ------------------------------
def test_translation_and_cell_vector_are_parent_cell_quantities():
    parent = _cubic("fcc")
    s = parent.repeat((3, 2, 4))
    nl = neighbor_list(s, 4.5)
    parent_A = parent.lattice.matrix
    np.testing.assert_allclose(nl.cell_vector @ parent_A, nl.vector, atol=1e-12)
    # cell_vector = translation + (target site position - origin site position) in the parent cell
    base = {lab: f for lab, f in zip(parent.site_label, parent.frac)}
    site_diff = np.array([base[t] - base[o] for o, t in zip(nl.origin_site, nl.target_site)])
    np.testing.assert_allclose(nl.cell_vector, nl.translation + site_diff, atol=1e-12)


def test_symmetric_vector_selection_reproduces_discus_set_vec():
    """DISCUS: an atom on site is1 in cell c gets site is2 in cell c+t; an atom on is2 gets is1 in c-t."""
    N = (3, 4, 2)
    s = _cubic("fcc").repeat(N)            # DISCUS order: x fastest, 4 sites per cell
    nl = neighbor_list(s, 7.0)             # |(1.5, 0, -0.5)| a = 6.32 Å
    is1, is2, t = "S1", "S3", np.array([1, 0, -1])
    sel = nl.select(origin_site=is1, target_site=is2, translation=t, symmetric=True)

    def index(cell, site):            # 0-based DISCUS celltoindex
        cell = np.mod(cell, N)
        return ((cell[2] * N[1] + cell[1]) * N[0] + cell[0]) * 4 + site
    expected = set()
    for cz, cy, cx in itertools.product(range(N[2]), range(N[1]), range(N[0])):
        c = np.array([cx, cy, cz])
        expected.add((index(c, 0), index(c + t, 2)))
        expected.add((index(c, 2), index(c - t, 0)))
    got = {(int(a), int(b)) for a, b in zip(nl.node_atom_id[sel.i], nl.node_atom_id[sel.j])}
    assert got == expected
    assert sel.n_pairs == 2 * np.prod(N)


def test_select_criteria():
    s = _cubic("fcc").repeat((3, 3, 3))
    nl = neighbor_list(s, 6.0)
    assert np.all(nl.select(shell=2).distance == pytest.approx(4.0))
    assert nl.select(shell=2).coordination().max() == 6
    first = nl.select(shell=1)
    # exactly one first neighbour of every atom along +[110] (fractional direction of the parent cell)
    along = first.select(direction=(1, 1, 0), direction_frame="fractional")
    assert np.all(along.coordination() == 1)
    np.testing.assert_allclose(along.cell_vector, [[0.5, 0.5, 0.0]] * along.n_pairs, atol=1e-12)
    # Cartesian direction with an angular tolerance: the cone of 30° around +x holds the 4 neighbours with x = +a/2
    cone = first.select(direction=(1, 0, 0), angle_tol=46.0)
    assert np.all(cone.coordination() == 4)
    # distance range and value
    assert nl.select(distance=(2.8, 4.1)).n_pairs == nl.select(shell=1).n_pairs + nl.select(shell=2).n_pairs
    assert nl.select(distance=4.0).n_pairs == nl.select(shell=2).n_pairs
    # site selection
    s1s2 = nl.select(origin_site="S1", target_site=["S2"], shell=1)
    assert np.all(s1s2.coordination()[s1s2.i] == 4)
    # cell_vector selection equals the matching translation selection here
    a = nl.select(origin_site="S1", target_site="S2", cell_vector=(0.5, 0.5, 0.0))
    b = nl.select(origin_site="S1", target_site="S2", translation=(0, 0, 0))
    assert a.n_pairs == b.n_pairs == 27


def test_fractional_direction_uses_the_lattice_in_a_skewed_cell():
    lat = Lattice.from_parameters(3.0, 4.0, 5.0, 70, 80, 115)
    s = _config(lat, ["A"], [[0, 0, 0]])
    nl = neighbor_list(s, 12.0)
    sel = nl.select(direction=(1, 1, 0), direction_frame="fractional")
    assert sel.n_pairs >= 1
    cv = sel.cell_vector
    np.testing.assert_allclose(cv[:, 0], cv[:, 1], atol=1e-9)
    np.testing.assert_allclose(cv[:, 2], 0.0, atol=1e-9)
    assert (1, 1, 0) in {tuple(int(x) for x in t) for t in sel.translation}
    # the Cartesian direction (1, 1, 0) is a different line in this cell
    cart = nl.select(direction=(1, 1, 0))
    assert not np.array_equal(np.sort(cart.distance), np.sort(sel.distance))


def test_shell_numbering_follows_the_other_criteria():
    s = _cubic("bcc").repeat((2, 2, 2))
    nl = neighbor_list(s, 4.1)
    same_site = nl.select(origin_site="S1", target_site="S1")
    # on one sublattice of bcc the nearest shell is the cubic edge
    assert same_site.select(shell=1).distance == pytest.approx(4.0)
    with pytest.raises(ValueError, match="shell 2 requested"):
        same_site.select(shell=2)


def test_species_selection():
    lat = Lattice.from_parameters(4, 4, 4, 90, 90, 90)
    s = _config(lat, ["Cs", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]]).repeat((2, 2, 2))
    nl = neighbor_list(s, 3.5)
    cs_cl = nl.select(origin_species="Cs", target_species="Cl")
    assert cs_cl.n_pairs == 8 * 8 and set(cs_cl.target_species) == {"Cl"}
    assert nl.select(origin_species="Cs", target_species="Cs").n_pairs == 0


def test_unique_keeps_one_direction_and_reversal_restores_all():
    s = _config(Lattice.from_parameters(2.0, 2.3, 2.9, 88, 95, 101), ["X", "Y"],
                [[0.1, 0.2, 0.3], [0.6, 0.55, 0.8]])
    nl = neighbor_list(s, 5.0)
    u = nl.unique()
    assert 2 * u.n_pairs == nl.n_pairs
    fwd = set(_as_dict(u))
    rev = {(j, i, -a, -b, -c) for i, j, a, b, c in fwd}
    assert fwd.isdisjoint(rev)
    assert fwd | rev == set(_as_dict(nl))
    # an asymmetric selection loses nothing
    for origin in ("X1", "Y1"):          # Y1 -> X1 pairs are the non-canonical direction (i > j)
        half = nl.select(origin_site=origin)
        n_self = int(np.sum(half.i == half.j))
        assert n_self > 0 and half.unique().n_pairs == half.n_pairs - n_self // 2


def test_average_model_mixed_site_collapses_to_one_node():
    lat = Lattice.from_parameters(3, 3, 3, 90, 90, 90)
    s = Structure(lat, np.array(["Au", "Cu"], object), [[0, 0, 0], [0, 0, 0]], model="average",
                  occupancy=[0.25, 0.75], site_label=["M", "M"])
    nl = neighbor_list(s, 3.0)
    assert nl.n_nodes == 1 and nl.n_pairs == 6
    assert list(nl.node_species) == [None] and list(nl.node_atom_id) == [-1]
    with pytest.raises(ValueError, match="one atom per site"):
        nl.select(origin_species="Au")


def test_split_sites_never_reach_the_neighbor_list():
    """Structure itself refuses members of one site at different positions, so nodes are well defined."""
    lat = Lattice.from_parameters(3, 3, 3, 90, 90, 90)
    with pytest.raises(ValueError, match="differing fractional positions"):
        Structure(lat, np.array(["Cu", "Cu"], object), [[0, 0, 0], [0.05, 0, 0]], model="average",
                  occupancy=[0.5, 0.5], site_label=["M", "M"])


def test_errors():
    s = _cubic("sc")
    with pytest.raises(ValueError, match="cutoff"):
        neighbor_list(s, 0.0)
    with pytest.raises(ValueError, match="cutoff"):
        neighbor_list(s, float("nan"))
    with pytest.raises(TypeError):
        neighbor_list("not a structure", 3.0)
    mixed = _config(s.lattice, ["Cu"], [[0, 0, 0]], pbc=(True, True, False))
    with pytest.raises(NotImplementedError, match="mixed PBC"):
        neighbor_list(mixed, 3.0)
    nl = neighbor_list(s, 4.1)
    with pytest.raises(ValueError, match="translation"):
        nl.select(translation=(0.5, 0, 0))
    with pytest.raises(ValueError, match="direction_frame"):
        nl.select(direction=(1, 0, 0), direction_frame="polar")
    with pytest.raises(ValueError, match="shell"):
        nl.select(shell=0)
    with pytest.raises(AttributeError):
        nl.cutoff = 5.0
    with pytest.raises(ValueError):
        nl.distance[0] = 1.0


def test_structure_neighbors_method_and_summary_and_csv(tmp_path):
    s = _cubic("fcc").repeat((2, 2, 2))
    nl = s.neighbors(3.0, origin_site="S1")
    assert nl.n_pairs == 32 * 12 // 4
    rows = s.neighbors(4.1).summary()
    s1 = [r for r in rows if r["origin_site"] == "S1"]
    assert {(r["target_site"], r["shell"], r["coordination"]) for r in s1} == {
        ("S1", 1, 6.0), ("S2", 1, 4.0), ("S3", 1, 4.0), ("S4", 1, 4.0)}
    path = nl.to_csv(tmp_path / "pairs.csv")
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    assert len(rows) == nl.n_pairs
    assert float(rows[0]["distance"]) == nl.distance[0]
    assert {"tx", "ty", "tz"}.isdisjoint(rows[0]) and {"ta", "tb", "tc", "fa", "shift_a"} <= set(rows[0])


def test_large_supercell_is_fast():
    """62 500 atoms (fcc 25³), first shell: 750 000 directed pairs."""
    s = _cubic("fcc").repeat((25, 25, 25))
    t0 = time.perf_counter()
    nl = neighbor_list(s, 3.0)
    elapsed = time.perf_counter() - t0
    assert nl.n_pairs == 62_500 * 12
    assert elapsed < 30.0, elapsed
