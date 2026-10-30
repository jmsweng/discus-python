"""Tests for the Structure model (Package A, Feature A1).

Expected values come from independent computations (hand-built cases, literal
numbers), never from calling the function under test a second way.
"""
import math

import numpy as np
import pytest

from discuspdf import Lattice, Structure


# ------------------------------ helpers ------------------------------
def _cubic(a=5.0) -> Lattice:
    return Lattice.from_parameters(a, a, a, 90.0, 90.0, 90.0)


def _triclinic() -> Lattice:
    return Lattice.from_parameters(4.1, 5.3, 6.7, 78.0, 95.0, 112.0)


# ------------------------------ 1. construction and validation ------------------------------
def test_construction_shape_mismatch():
    lat = _cubic()
    with pytest.raises(ValueError, match="shape"):
        Structure(lat, np.array(["Ga"]), np.array([[0.25, 0.25]]), model="configuration")


def test_construction_nan_position():
    lat = _cubic()
    with pytest.raises(ValueError, match="NaN|inf|finite"):
        Structure(lat, np.array(["Ga"]), np.array([[0.25, float("nan"), 0.25]]), model="configuration")


def test_construction_inf_position():
    lat = _cubic()
    with pytest.raises(ValueError, match="NaN|inf|finite"):
        Structure(lat, np.array(["Ga"]), np.array([[0.25, float("inf"), 0.25]]), model="configuration")


def test_construction_occupancy_zero():
    lat = _cubic()
    with pytest.raises(ValueError, match="occupancy"):
        Structure(lat, np.array(["Ga"]), np.array([[0.25, 0.25, 0.25]]),
                  model="configuration", occupancy=np.array([0.0]))


def test_construction_occupancy_gt_1():
    lat = _cubic()
    with pytest.raises(ValueError, match="occupancy"):
        Structure(lat, np.array(["Ga"]), np.array([[0.25, 0.25, 0.25]]),
                  model="configuration", occupancy=np.array([1.5]))


def test_construction_site_occupancy_sum_gt_1():
    lat = _cubic()
    with pytest.raises(ValueError, match="occupancy"):
        Structure(lat, np.array(["Cu", "Au"]),
                  np.array([[0.25, 0.25, 0.25], [0.25, 0.25, 0.25]]),
                  model="average",
                  occupancy=np.array([0.7, 0.5]),
                  site_label=np.array(["S1", "S1"]))


def test_construction_site_differing_positions():
    lat = _cubic()
    with pytest.raises(ValueError, match="differing|position"):
        Structure(lat, np.array(["Cu", "Au"]),
                  np.array([[0.25, 0.25, 0.25], [0.30, 0.25, 0.25]]),
                  model="average",
                  occupancy=np.array([0.5, 0.5]),
                  site_label=np.array(["S1", "S1"]))


def test_construction_site_differing_b_iso():
    lat = _cubic()
    with pytest.raises(ValueError, match="b_iso"):
        Structure(lat, np.array(["Cu", "Au"]),
                  np.array([[0.25, 0.25, 0.25], [0.25, 0.25, 0.25]]),
                  model="average",
                  occupancy=np.array([0.5, 0.5]),
                  site_label=np.array(["S1", "S1"]),
                  b_iso=np.array([1.0, 2.0]))


def test_construction_duplicate_atom_id():
    lat = _cubic()
    with pytest.raises(ValueError, match="unique|atom_id"):
        Structure(lat, np.array(["Ga", "As"]),
                  np.array([[0.25, 0.25, 0.25], [0.75, 0.75, 0.75]]),
                  model="configuration",
                  atom_id=np.array([0, 0]))


def test_construction_negative_b_iso():
    lat = _cubic()
    with pytest.raises(ValueError, match="b_iso"):
        Structure(lat, np.array(["Ga"]), np.array([[0.25, 0.25, 0.25]]),
                  model="configuration", b_iso=np.array([-1.0]))


def test_construction_invalid_species_empty():
    lat = _cubic()
    with pytest.raises(ValueError, match="non-empty|species"):
        Structure(lat, np.array([""]), np.array([[0.25, 0.25, 0.25]]), model="configuration")


def test_construction_invalid_species_whitespace():
    lat = _cubic()
    with pytest.raises(ValueError, match="whitespace|species"):
        Structure(lat, np.array(["Ga As"]), np.array([[0.25, 0.25, 0.25]]), model="configuration")


def test_construction_invalid_species_too_long():
    lat = _cubic()
    with pytest.raises(ValueError, match="8 char|species"):
        Structure(lat, np.array(["GaAsO2NaX"]), np.array([[0.25, 0.25, 0.25]]), model="configuration")


def test_construction_missing_model():
    lat = _cubic()
    with pytest.raises((ValueError, TypeError), match="model"):
        Structure(lat, np.array(["Ga"]), np.array([[0.25, 0.25, 0.25]]))


def test_construction_configuration_occupancy_lt_1():
    lat = _cubic()
    with pytest.raises(ValueError, match="configuration|1.0|average"):
        Structure(lat, np.array(["Ga"]), np.array([[0.25, 0.25, 0.25]]),
                  model="configuration", occupancy=np.array([0.5]))


def test_construction_two_atoms_one_configuration_site():
    lat = _cubic()
    with pytest.raises(ValueError, match="configuration|one atom"):
        Structure(lat, np.array(["Ga", "As"]),
                  np.array([[0.25, 0.25, 0.25], [0.25, 0.25, 0.25]]),
                  model="configuration",
                  site_label=np.array(["S1", "S1"]))


def test_construction_cell_index_outside_supercell():
    lat = _cubic()
    with pytest.raises(ValueError, match="cell_index|supercell"):
        Structure(lat, np.array(["Ga"]), np.array([[0.25, 0.25, 0.25]]),
                  model="configuration",
                  cell_index=np.array([[2, 0, 0]]),
                  supercell=(2, 2, 2))


def test_construction_overlap_across_pbc():
    """Two atoms closer than 0.1 Å across a periodic boundary raise."""
    lat = _cubic(5.0)
    # x = 0.001 and x = 0.999 in a 5 Å cell: distance = 0.01 Å < 0.1 Å
    with pytest.raises(ValueError, match="0.1|apart|overlap"):
        Structure(lat, np.array(["Ga", "As"]),
                  np.array([[0.001, 0.5, 0.5], [0.999, 0.5, 0.5]]),
                  model="configuration")


def test_construction_average_close_sites_accepted():
    """Two different average-model sites at 0.05 Å with occupancies 0.5 + 0.5 is accepted.

    The overlap check applies only when two occupancies sum to > 1.
    """
    lat = _cubic(5.0)
    # 0.05 Å apart in a 5 Å cell: 0.01 in fractional units
    # Different sites (S1, S2), so the site-consistency check doesn't apply.
    s = Structure(lat, np.array(["Cu", "Au"]),
                 np.array([[0.25, 0.25, 0.25], [0.26, 0.25, 0.25]]),
                 model="average",
                 occupancy=np.array([0.5, 0.5]),
                 site_label=np.array(["S1", "S2"]))
    assert s.n_atoms == 2


def test_construction_average_close_sites_rejected():
    """Two different average-model sites at 0.05 Å with occupancies 0.6 + 0.6 = 1.2 > 1 raises."""
    lat = _cubic(5.0)
    with pytest.raises(ValueError, match="0.1|apart|overlap"):
        Structure(lat, np.array(["Cu", "Au"]),
                  np.array([[0.25, 0.25, 0.25], [0.26, 0.25, 0.25]]),
                  model="average",
                  occupancy=np.array([0.6, 0.6]),
                  site_label=np.array(["S1", "S2"]))


# ------------------------------ 2. mixed-site preservation ------------------------------
def test_mixed_site_preservation():
    """An average Cu/Au site (0.75/0.25) plus a vacancy site (0.9) keeps rows,
    occupancies, sites and composition() through select, with_frac, replace and concat."""
    lat = _cubic(4.0)
    # S1 at (0.25, 0.25, 0.25) with occ 0.75+0.25=1.0
    # S2 at (0.75, 0.75, 0.75) with occ 0.9
    # Distance S1-S2: 0.5*sqrt(3)*4 = 3.46 Å > 0.1 Å, so overlap check passes
    s = Structure(lat,
                  np.array(["Cu", "Au", "Cu"]),
                  np.array([[0.25, 0.25, 0.25], [0.25, 0.25, 0.25], [0.75, 0.75, 0.75]]),
                  model="average",
                  occupancy=np.array([0.75, 0.25, 0.9]),
                  site_label=np.array(["S1", "S1", "S2"]))

    # select
    s2 = s.select(species="Cu")
    assert s2.n_atoms == 2
    assert s2.composition() == {"Cu": 0.75 + 0.9}

    # with_frac
    s3 = s.with_frac(np.array([[0.3, 0.25, 0.25], [0.3, 0.25, 0.25], [0.75, 0.75, 0.75]]))
    assert s3.n_atoms == 3
    assert s3.composition() == {"Cu": 1.65, "Au": 0.25}

    # replace
    s4 = s.replace(occupancy=np.array([0.7, 0.3, 0.9]))
    assert s4.occupancy[0] == pytest.approx(0.7)
    assert s4.composition() == {"Cu": 1.6, "Au": 0.3}

    # concat
    s5 = s.concat(s, renumber=True)
    assert s5.n_atoms == 6
    assert s5.composition() == {"Cu": 3.3, "Au": 0.5}


# ------------------------------ 3. coordinates ------------------------------
def test_frac_cartesian_round_trip():
    """frac → cartesian → from_cartesian round trip on a triclinic cell."""
    lat = _triclinic()
    rng = np.random.default_rng(42)
    frac = rng.uniform(0, 1, (10, 3))
    s = Structure(lat, np.array(["Ga"] * 10), frac, model="configuration")
    cart = s.cartesian
    s2 = Structure.from_cartesian(lat, np.array(["Ga"] * 10), cart, model="configuration")
    assert np.allclose(s2.frac, frac, rtol=1e-12, atol=1e-14)
    assert np.array_equal(s2.atom_id, s.atom_id)
    assert np.array_equal(s2.species, s.species)


def test_wrapped_literals():
    """wrapped() on positions such as (1.2, -0.3, 0.5) gives expected literals."""
    lat = _cubic(5.0)
    s = Structure(lat, np.array(["Ga", "As", "Cu"]),
                 np.array([[1.2, -0.3, 0.5], [0.7, 0.9, 0.1], [0.0, 0.0, 0.0]]),
                 model="configuration")
    s_w = s.wrapped()
    assert np.allclose(s_w.frac[0], [0.2, 0.7, 0.5], atol=1e-12)
    assert np.allclose(s_w.frac[1], [0.7, 0.9, 0.1], atol=1e-12)
    assert np.allclose(s_w.frac[2], [0.0, 0.0, 0.0], atol=1e-12)
    # cell_index unchanged
    assert np.array_equal(s_w.cell_index, s.cell_index)


def test_wrapped_non_periodic_raises():
    lat = _cubic(5.0)
    s = Structure(lat, np.array(["Ga"]), np.array([[0.5, 0.5, 0.5]]),
                 model="configuration", pbc=(False, False, False))
    with pytest.raises(ValueError, match="periodic"):
        s.wrapped()


# ------------------------------ 4. non-aliasing and immutability ------------------------------
def test_input_array_not_aliased():
    """Editing input arrays after construction does not change the structure."""
    lat = _cubic()
    frac_in = np.array([[0.25, 0.25, 0.25]])
    s = Structure(lat, np.array(["Ga"]), frac_in, model="configuration", min_separation=None)
    frac_in[0, 0] = 0.99
    assert s.frac[0, 0] == pytest.approx(0.25)


def test_returned_array_read_only():
    """Editing returned arrays is blocked."""
    lat = _cubic()
    s = Structure(lat, np.array(["Ga"]), np.array([[0.25, 0.25, 0.25]]), model="configuration")
    with pytest.raises(ValueError):
        s.frac[0, 0] = 0.99


def test_attribute_set_raises():
    lat = _cubic()
    s = Structure(lat, np.array(["Ga"]), np.array([[0.25, 0.25, 0.25]]), model="configuration")
    with pytest.raises(AttributeError):
        s.new_attr = 1
    with pytest.raises(AttributeError):
        del s._frac


def test_functional_edit_leaves_original_untouched():
    lat = _cubic()
    s = Structure(lat, np.array(["Ga", "As"]),
                 np.array([[0.25, 0.25, 0.25], [0.75, 0.75, 0.75]]),
                 model="configuration")
    frac_before = s.frac.copy()
    s2 = s.with_frac(np.array([[0.3, 0.25, 0.25], [0.75, 0.75, 0.75]]))
    assert np.array_equal(s.frac, frac_before)
    assert s2.frac[0, 0] == pytest.approx(0.3)


# ------------------------------ 5. identity ------------------------------
def test_select_keeps_ids_and_order():
    lat = _cubic()
    s = Structure(lat, np.array(["Ga", "As", "Ga"]),
                 np.array([[0.25, 0.25, 0.25], [0.75, 0.75, 0.75], [0.5, 0.5, 0.5]]),
                 model="configuration",
                 atom_id=np.array([5, 3, 7]))
    s2 = s.select(species="Ga")
    assert np.array_equal(s2.atom_id, [5, 7])
    assert np.array_equal(s2.species, ["Ga", "Ga"])


def test_remove_keeps_ids_and_order():
    lat = _cubic()
    s = Structure(lat, np.array(["Ga", "As", "Ga"]),
                 np.array([[0.25, 0.25, 0.25], [0.75, 0.75, 0.75], [0.5, 0.5, 0.5]]),
                 model="configuration",
                 atom_id=np.array([5, 3, 7]))
    s2 = s.remove(np.array([3]))
    assert np.array_equal(s2.atom_id, [5, 7])


def test_concat_id_clash_raises():
    lat = _cubic()
    s1 = Structure(lat, np.array(["Ga"]), np.array([[0.25, 0.25, 0.25]]),
                   model="configuration", atom_id=np.array([0]))
    s2 = Structure(lat, np.array(["As"]), np.array([[0.75, 0.75, 0.75]]),
                   model="configuration", atom_id=np.array([0]))
    with pytest.raises(ValueError, match="clash|renumber"):
        s1.concat(s2)


def test_concat_renumber_assigns_max_plus_1():
    lat = _cubic()
    s1 = Structure(lat, np.array(["Ga"]), np.array([[0.25, 0.25, 0.25]]),
                   model="configuration", atom_id=np.array([5]))
    s2 = Structure(lat, np.array(["As"]), np.array([[0.75, 0.75, 0.75]]),
                   model="configuration", atom_id=np.array([3]))
    s3 = s1.concat(s2, renumber=True)
    # max(5, 3) + 1 = 6, so s2's atom gets id 6
    assert np.array_equal(s3.atom_id, [5, 6])


def test_renumbered_is_only_method_that_renumbers():
    lat = _cubic()
    s = Structure(lat, np.array(["Ga", "As"]),
                 np.array([[0.25, 0.25, 0.25], [0.75, 0.75, 0.75]]),
                 model="configuration",
                 atom_id=np.array([10, 20]))
    s2 = s.renumbered()
    assert np.array_equal(s2.atom_id, [0, 1])
    # Original unchanged
    assert np.array_equal(s.atom_id, [10, 20])


# ------------------------------ 5b. ids mean atom_id ------------------------------
def test_select_ids_are_atom_ids_not_rows():
    """ids= selects by atom_id values, not row numbers."""
    lat = _cubic()
    s = Structure(lat, np.array(["Ga", "As", "Cu"]),
                 np.array([[0.25, 0.25, 0.25], [0.75, 0.75, 0.75], [0.5, 0.5, 0.5]]),
                 model="configuration",
                 atom_id=np.array([10, 20, 30]))
    # id 20 is row 1; id 1 does not exist
    s2 = s.select(ids=[20])
    assert np.array_equal(s2.atom_id, [20])
    with pytest.raises(ValueError, match="unknown atom_id"):
        s.select(ids=[1])


def test_remove_unknown_id_raises():
    lat = _cubic()
    s = Structure(lat, np.array(["Ga", "As"]),
                 np.array([[0.25, 0.25, 0.25], [0.75, 0.75, 0.75]]),
                 model="configuration",
                 atom_id=np.array([1, 2]))
    with pytest.raises(ValueError, match="unknown atom_id"):
        s.remove([99])


def test_remove_selects_by_atom_id_not_row():
    """atom_id = [1, 2, 0]: remove([0]) removes the atom with id 0 (row 2)."""
    lat = _cubic()
    s = Structure(lat, np.array(["Ga", "As", "Cu"]),
                 np.array([[0.25, 0.25, 0.25], [0.75, 0.75, 0.75], [0.5, 0.5, 0.5]]),
                 model="configuration",
                 atom_id=np.array([1, 2, 0]))
    s2 = s.remove([0])
    assert np.array_equal(s2.atom_id, [1, 2])


def test_select_ids_preserve_row_order():
    lat = _cubic()
    s = Structure(lat, np.array(["Ga", "As", "Cu", "Au"]),
                 np.array([[0.1, 0.2, 0.3], [0.4, 0.5, 0.6], [0.2, 0.3, 0.4], [0.5, 0.6, 0.7]]),
                 model="configuration",
                 atom_id=np.array([30, 10, 20, 40]))
    s2 = s.select(ids=[40, 10, 30])
    assert np.array_equal(s2.atom_id, [30, 10, 40])


def test_remove_after_select():
    """remove on a structure that has already been through select."""
    lat = _cubic()
    s = Structure(lat, np.array(["Ga", "As", "Cu"]),
                 np.array([[0.25, 0.25, 0.25], [0.75, 0.75, 0.75], [0.5, 0.5, 0.5]]),
                 model="configuration",
                 atom_id=np.array([10, 20, 30]))
    s2 = s.select(ids=[10, 30])
    assert np.array_equal(s2.atom_id, [10, 30])
    s3 = s2.remove([30])
    assert np.array_equal(s3.atom_id, [10])


def test_select_rows_keyword_selects_by_row_position():
    lat = _cubic()
    s = Structure(lat, np.array(["Ga", "As", "Cu"]),
                 np.array([[0.25, 0.25, 0.25], [0.75, 0.75, 0.75], [0.5, 0.5, 0.5]]),
                 model="configuration",
                 atom_id=np.array([10, 20, 30]))
    s2 = s.select(rows=[2, 0])
    assert np.array_equal(s2.atom_id, [10, 30])


def test_distance_and_displacement_use_atom_ids():
    """distance/displacement take atom_id values, not row numbers."""
    lat = _cubic(5.0)
    s = Structure(lat, np.array(["Ga", "As"]),
                 np.array([[0.05, 0.5, 0.5], [0.95, 0.5, 0.5]]),
                 model="configuration",
                 atom_id=np.array([10, 20]))
    assert s.distance(10, 20) == pytest.approx(0.5, rel=1e-12)
    d = s.displacement(10, 20)
    assert np.allclose(np.abs(d), [0.5, 0.0, 0.0], atol=1e-12)
    with pytest.raises(ValueError, match="unknown atom_id"):
        s.distance(0, 1)


def test_distance_matrix_ids_are_atom_ids():
    lat = _cubic(5.0)
    s = Structure(lat, np.array(["Ga", "As", "Cu"]),
                 np.array([[0.05, 0.5, 0.5], [0.95, 0.5, 0.5], [0.5, 0.5, 0.5]]),
                 model="configuration",
                 atom_id=np.array([10, 20, 30]))
    m = s.distance_matrix(ids=[10, 30])
    assert m.shape == (2, 2)
    assert m[0, 0] == pytest.approx(0.0, abs=1e-12)
    with pytest.raises(ValueError, match="unknown atom_id"):
        s.distance_matrix(ids=[0])


# ------------------------------ 6. periodicity ------------------------------
def test_periodic_distance_across_boundary():
    """0.05 → 0.95 in 5 Å gives 0.5 Å."""
    lat = _cubic(5.0)
    s = Structure(lat, np.array(["Ga", "As"]),
                 np.array([[0.05, 0.5, 0.5], [0.95, 0.5, 0.5]]),
                 model="configuration")
    d = s.distance(0, 1)
    assert d == pytest.approx(0.5, rel=1e-12)


def test_non_periodic_distance():
    """Same atoms with pbc all False give 4.5 Å."""
    lat = _cubic(5.0)
    s = Structure(lat, np.array(["Ga", "As"]),
                 np.array([[0.05, 0.5, 0.5], [0.95, 0.5, 0.5]]),
                 model="configuration",
                 pbc=(False, False, False))
    d = s.distance(0, 1)
    assert d == pytest.approx(4.5, rel=1e-12)


def test_mixed_pbc_raises():
    lat = _cubic(5.0)
    # Use positions far enough apart to pass the overlap check
    s = Structure(lat, np.array(["Ga", "As"]),
                 np.array([[0.1, 0.5, 0.5], [0.9, 0.5, 0.5]]),
                 model="configuration",
                 pbc=(True, False, True),
                 min_separation=None)
    with pytest.raises(NotImplementedError, match="mixed PBC"):
        s.distance(0, 1)


def test_skewed_cell_distance_matches_lattice():
    """A skewed-cell distance matches Lattice.distance_pbc."""
    lat = Lattice.from_parameters(5.0, 5.0, 5.0, 90.0, 90.0, 20.0)
    f1 = np.array([0.1, 0.2, 0.3])
    f2 = np.array([0.9, 0.8, 0.7])
    s = Structure(lat, np.array(["Ga", "As"]),
                 np.array([f1, f2]), model="configuration")
    d_struct = s.distance(0, 1)
    d_lattice = lat.distance_pbc(f1, f2)
    assert d_struct == pytest.approx(d_lattice, rel=1e-12)


# ------------------------------ 7. summaries ------------------------------
def test_counts_composition_formula_gaas():
    lat = _cubic(4.5)
    s = Structure(lat, np.array(["GA", "AS"]),
                 np.array([[0.25, 0.25, 0.25], [0.75, 0.75, 0.75]]),
                 model="configuration")
    assert s.counts() == {"GA": 1, "AS": 1}
    assert s.composition() == {"GA": 1.0, "AS": 1.0}
    assert s.formula() == "GAAS"


def test_atomic_numbers_gaas():
    lat = _cubic(4.5)
    s = Structure(lat, np.array(["GA", "AS"]),
                 np.array([[0.25, 0.25, 0.25], [0.75, 0.75, 0.75]]),
                 model="configuration")
    z = s.atomic_numbers()
    assert z[0] == 31  # Ga
    assert z[1] == 33  # As


def test_atomic_numbers_cupau_average():
    lat = _cubic(4.0)
    s = Structure(lat, np.array(["Cu", "Au"]),
                 np.array([[0.25, 0.25, 0.25], [0.25, 0.25, 0.25]]),
                 model="average",
                 occupancy=np.array([0.75, 0.25]),
                 site_label=np.array(["S1", "S1"]))
    z = s.atomic_numbers()
    assert z[0] == 29  # Cu
    assert z[1] == 79  # Au


def test_atomic_numbers_unknown_raises():
    lat = _cubic()
    s = Structure(lat, np.array(["Xx"]), np.array([[0.25, 0.25, 0.25]]),
                 model="configuration")
    with pytest.raises(ValueError, match="unknown|element"):
        s.atomic_numbers()


def test_composition_cupau_average():
    lat = _cubic(4.0)
    s = Structure(lat, np.array(["Cu", "Au"]),
                 np.array([[0.25, 0.25, 0.25], [0.25, 0.25, 0.25]]),
                 model="average",
                 occupancy=np.array([0.75, 0.25]),
                 site_label=np.array(["S1", "S1"]))
    comp = s.composition()
    assert comp["Cu"] == pytest.approx(0.75)
    assert comp["Au"] == pytest.approx(0.25)


# ------------------------------ 8. scale ------------------------------
def test_scale_1e6_atoms():
    """Construct and validate a 10⁶-atom configuration with min_separation=None."""
    import tracemalloc
    lat = _cubic(5.0)
    n = 1_000_000
    rng = np.random.default_rng(42)
    # Spread atoms on a grid to avoid overlaps
    per_axis = int(round(n ** (1 / 3)))
    coords = np.stack(np.meshgrid(
        np.linspace(0.001, 0.999, per_axis),
        np.linspace(0.001, 0.999, per_axis),
        np.linspace(0.001, 0.999, per_axis),
        indexing="ij",
    ), axis=-1).reshape(-1, 3)
    coords = coords[:n]
    species = np.array(["Ga"] * n, dtype=object)
    tracemalloc.start()
    s = Structure(lat, species, coords, model="configuration", min_separation=None)
    cur, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    assert s.n_atoms == n
    # select and edit
    s2 = s.select(ids=np.arange(100))
    assert s2.n_atoms == 100
    s3 = s.with_frac(coords + 0.001)
    assert s3.n_atoms == n
    print(f"  10^6 atoms: peak {peak / 1e6:.1f} MB")


def test_scale_20k_atoms_skips_separation():
    """Construct 20 000 atoms; the separation check is skipped at n >= 20000."""
    import time
    lat = _cubic(50.0)  # large cell to avoid overlaps
    n = 20_000
    rng = np.random.default_rng(42)
    coords = rng.uniform(0.01, 0.99, (n, 3))
    species = np.array(["Ga"] * n, dtype=object)
    t0 = time.perf_counter()
    s = Structure(lat, species, coords, model="configuration")
    t1 = time.perf_counter()
    assert s.n_atoms == n
    print(f"  20k atoms (separation check skipped): {t1 - t0:.2f} s")


def test_scale_1k_atoms_with_separation():
    """Construct 1 000 atoms with the separation check on (below the 20k limit)."""
    import time
    lat = _cubic(50.0)  # large cell to avoid overlaps
    n = 1_000
    rng = np.random.default_rng(42)
    coords = rng.uniform(0.01, 0.99, (n, 3))
    species = np.array(["Ga"] * n, dtype=object)
    t0 = time.perf_counter()
    s = Structure(lat, species, coords, model="configuration")
    t1 = time.perf_counter()
    assert s.n_atoms == n
    print(f"  1k atoms with separation check: {t1 - t0:.2f} s")


# ------------------------------ 8b. isclose ------------------------------
def _two_atom_base():
    lat = _cubic(5.0)
    return Structure(lat, np.array(["Ga", "As"]),
                    np.array([[0.25, 0.25, 0.25], [0.75, 0.75, 0.75]]),
                    model="configuration",
                    atom_id=np.array([1, 2]),
                    site_label=np.array(["S1", "S2"]))


def test_isclose_requires_matching_species():
    s = _two_atom_base()
    other = s.replace(species=np.array(["Ga", "Xe"]))
    assert not s.isclose(other)


def test_isclose_requires_matching_atom_id():
    s = _two_atom_base()
    other = s.replace(atom_id=np.array([1, 3]))
    assert not s.isclose(other)


def test_isclose_requires_matching_pbc():
    s = _two_atom_base()
    other = s.replace(pbc=(True, True, False))
    assert not s.isclose(other)


def test_isclose_requires_matching_site_label_and_cell_index():
    s = _two_atom_base()
    assert not s.isclose(s.replace(site_label=np.array(["T1", "S2"])))
    assert not s.isclose(s.replace(cell_index=np.array([[0, 0, 0], [0, 0, 1]]), supercell=(1, 1, 2)))


def test_isclose_requires_matching_supercell_and_model():
    s = _two_atom_base()
    assert not s.isclose(s.replace(supercell=(2, 1, 1)))
    other = Structure(_cubic(5.0), np.array(["Ga", "As"]),
                     np.array([[0.25, 0.25, 0.25], [0.75, 0.75, 0.75]]),
                     model="average",
                     occupancy=np.array([0.5, 0.5]),
                     atom_id=np.array([1, 2]),
                     site_label=np.array(["S1", "S2"]))
    assert not s.isclose(other)


def test_isclose_b_iso_presence_must_match():
    s = _two_atom_base()
    with_b = s.replace(b_iso=np.array([1.0, 1.0]))
    assert not s.isclose(with_b)
    assert not with_b.isclose(s)
    # both have b_iso: numeric closeness applies
    assert with_b.isclose(with_b.replace(b_iso=np.array([1.0 + 1e-13, 1.0])))


def test_isclose_tolerance_applies_to_frac():
    s = _two_atom_base()
    other = s.replace(frac=np.array([[0.25 + 1e-13, 0.25, 0.25], [0.75, 0.75, 0.75]]))
    assert s.isclose(other)
    assert not s.isclose(other, rtol=0.0, atol=0.0)


# ------------------------------ 8c. validation ------------------------------
def test_site_label_validation():
    lat = _cubic()
    # non-string
    with pytest.raises(ValueError, match="site_label"):
        Structure(lat, np.array(["Ga"]), np.array([[0.25, 0.25, 0.25]]),
                  model="configuration", site_label=np.array([3.5]))
    # empty
    with pytest.raises(ValueError, match="site_label"):
        Structure(lat, np.array(["Ga"]), np.array([[0.25, 0.25, 0.25]]),
                  model="configuration", site_label=np.array([""]))
    # whitespace
    with pytest.raises(ValueError, match="whitespace|site_label"):
        Structure(lat, np.array(["Ga"]), np.array([[0.25, 0.25, 0.25]]),
                  model="configuration", site_label=np.array(["S 1"]))
    # too long (> 16)
    with pytest.raises(ValueError, match="16|site_label"):
        Structure(lat, np.array(["Ga"]), np.array([[0.25, 0.25, 0.25]]),
                  model="configuration", site_label=np.array(["S" * 17]))
    # internal key separators
    with pytest.raises(ValueError, match="site_label"):
        Structure(lat, np.array(["Ga"]), np.array([[0.25, 0.25, 0.25]]),
                  model="configuration", site_label=np.array(["S|1"]))
    with pytest.raises(ValueError, match="site_label"):
        Structure(lat, np.array(["Ga"]), np.array([[0.25, 0.25, 0.25]]),
                  model="configuration", site_label=np.array(["S,1"]))
    # boundary: 16 characters is accepted
    s = Structure(lat, np.array(["Ga"]), np.array([[0.25, 0.25, 0.25]]),
                 model="configuration", site_label=np.array(["S" * 16]))
    assert s.n_atoms == 1


def test_pbc_validation():
    lat = _cubic()
    for bad in (("yes", 0, None), (True, 1, False), (True, False), (1, 0, 0)):
        with pytest.raises(ValueError, match="pbc"):
            Structure(lat, np.array(["Ga"]), np.array([[0.25, 0.25, 0.25]]),
                      model="configuration", pbc=bad)
    # np.bool_ is accepted
    s = Structure(lat, np.array(["Ga"]), np.array([[0.25, 0.25, 0.25]]),
                 model="configuration", pbc=(np.bool_(True), True, np.bool_(True)))
    assert s.pbc == (True, True, True)
    # mixed PBC is accepted
    s = Structure(lat, np.array(["Ga"]), np.array([[0.25, 0.25, 0.25]]),
                 model="configuration", pbc=(np.bool_(True), False, np.bool_(True)),
                 min_separation=None)
    assert s.pbc == (True, False, True)


def test_site_error_messages_format_site_and_cell():
    """Site-consistency and overlap errors format the site as
    ``site 'S0' in cell (0, 0, 0)``, in both the small-n and large-n paths."""
    lat = _cubic(50.0)

    def _build(n, model, occupancy=None, site_label=None):
        rng = np.random.default_rng(0)
        coords = np.zeros((n, 3))
        # spread atoms so the overlap check passes, then duplicate one site
        for i in range(n):
            coords[i] = [0.01 * (i % 7 + 1), 0.01 * ((i // 7) % 7 + 1), 0.01 * ((i // 49) % 7 + 1)]
        if site_label is None:
            site_label = np.array([f"S{i}" for i in range(n)], dtype=object)
        if occupancy is None:
            occupancy = np.ones(n)
        return Structure(lat, np.array(["Ga"] * n, dtype=object), coords,
                         model=model, occupancy=occupancy, site_label=site_label)

    # small-n path: configuration, two atoms on one site
    n = 10
    sl = np.array([f"S{i}" for i in range(n)], dtype=object)
    sl[1] = sl[0]
    with pytest.raises(ValueError, match=r"site 'S0' in cell \(0, 0, 0\)"):
        _build(n, "configuration", site_label=sl)
    # large-n path (> 1000): same expectation
    n = 1100
    sl = np.array([f"S{i}" for i in range(n)], dtype=object)
    sl[1] = sl[0]
    with pytest.raises(ValueError, match=r"site 'S0' in cell \(0, 0, 0\)"):
        _build(n, "configuration", site_label=sl)


# ------------------------------ 9. example script ------------------------------
def test_example_script_runs():
    """CI keeps examples/structure_basics.py working."""
    import subprocess
    import sys
    import pathlib
    script = pathlib.Path(__file__).parent.parent / "examples" / "structure_basics.py"
    result = subprocess.run([sys.executable, str(script)], capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    assert "gaas" in result.stdout.lower() or "structure" in result.stdout.lower()
