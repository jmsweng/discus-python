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

    # concat with a compatible average model (different sites, no overlap)
    other = Structure(lat,
                      np.array(["Cu", "Au"]),
                      np.array([[0.25, 0.75, 0.25], [0.25, 0.75, 0.25]]),
                      model="average",
                      occupancy=np.array([0.5, 0.5]),
                      site_label=np.array(["S3", "S3"]))
    s5 = s.concat(other, renumber=True)
    assert s5.n_atoms == 5
    assert s5.composition() == pytest.approx({"Cu": 2.15, "Au": 0.75})
    assert list(s5.occupancy) == [0.75, 0.25, 0.9, 0.5, 0.5]  # no row became fully occupied


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
    s3 = s.with_frac(coords + 0.001, min_separation=None)  # grid spacing is far below 0.1 Å
    assert s3.n_atoms == n
    print(f"  10^6 atoms: peak {peak / 1e6:.1f} MB")


def _gaas_supercell(N):
    """GaAs zinc-blende N x N x N supercell as a configuration (8 atoms per cell)."""
    a = 5.6537
    fcc = np.array([[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]])
    base = np.vstack([fcc, fcc + 0.25])
    cells = np.array(np.meshgrid(range(N), range(N), range(N), indexing="ij")).reshape(3, -1).T
    frac = ((base[None, :, :] + cells[:, None, :]) / N).reshape(-1, 3)
    labels = [f"Ga{i}" for i in range(1, 5)] + [f"As{i}" for i in range(1, 5)]
    return dict(
        lattice=Lattice.from_parameters(a * N, a * N, a * N, 90.0, 90.0, 90.0),
        species=np.tile(["Ga"] * 4 + ["As"] * 4, len(cells)).astype(object),
        frac=frac,
        site_label=np.tile(labels, len(cells)).astype(object),
        cell_index=np.repeat(cells, 8, axis=0),
        supercell=(N, N, N),
    )


def test_overlap_check_runs_at_scale():
    """The overlap check is never skipped: a 16x16x16 GaAs supercell (32 768 atoms)
    validates quickly, and a single planted overlap in it is found."""
    import time
    d = _gaas_supercell(16)
    t0 = time.perf_counter()
    s = Structure(d["lattice"], d["species"], d["frac"], model="configuration",
                  site_label=d["site_label"], cell_index=d["cell_index"], supercell=d["supercell"])
    elapsed = time.perf_counter() - t0
    assert s.n_atoms == 32_768
    assert elapsed < 30.0, f"overlap check too slow: {elapsed:.1f} s"
    frac = d["frac"].copy()
    frac[-1] = frac[0] + np.array([0.05, 0.0, 0.0]) / (5.6537 * 16)  # 0.05 Å from atom 0
    with pytest.raises(ValueError, match="apart"):
        Structure(d["lattice"], d["species"], frac, model="configuration",
                  site_label=d["site_label"], cell_index=d["cell_index"], supercell=d["supercell"])


def _brute_min_overlap(model, frac, occ, site, lat, periodic, cut):
    import itertools
    shifts = np.array(list(itertools.product(range(-3, 4), repeat=3)), float)
    sid = np.unique(site, return_inverse=True)[1].reshape(-1)
    site_occ = np.bincount(sid, weights=occ)
    i, j = np.triu_indices(len(frac), 1)
    keep = sid[i] != sid[j]
    if model == "average":
        keep &= (site_occ[sid[i]] + site_occ[sid[j]]) > 1.0 + 1e-9
    i, j = i[keep], j[keep]
    df = frac[j] - frac[i]
    if periodic:
        base = df - np.rint(df)
        dist = np.linalg.norm((base[:, None, :] - shifts[None]) @ lat.matrix, axis=2).min(axis=1)
    else:
        dist = np.linalg.norm(df @ lat.matrix, axis=1)
    dist = dist[dist < cut]
    return float(dist.min()) if len(dist) else None


def test_overlap_cell_list_matches_brute_force():
    """Cell-list overlap check agrees with brute force (all pairs, all images)
    on random structures, incl. skewed and smaller-than-cutoff cells."""
    from discuspdf.structure import _check_overlap
    rng = np.random.default_rng(0)
    cells = [(5, 5, 5, 90, 90, 90), (5, 5, 5, 90, 90, 20), (4.1, 5.3, 6.7, 78, 95, 112),
             (0.25, 0.3, 0.28, 90, 90, 90)]
    for p in cells:
        lat = Lattice.from_parameters(*p)
        inv = np.linalg.inv(lat.matrix)
        for _ in range(10):
            n = int(rng.integers(20, 80))
            frac = rng.uniform(-1.5, 2.5, (n, 3))
            for _ in range(int(rng.integers(0, 4))):  # plant near pairs, some across the boundary
                i, j = rng.integers(n, size=2)
                if i != j:
                    v = rng.normal(size=3)
                    v *= rng.choice([0.03, 0.0999, 0.1001, 0.2]) / np.linalg.norm(v)
                    frac[j] = frac[i] + v @ inv + rng.integers(-2, 3, 3)
            model = str(rng.choice(["configuration", "average"]))
            occ = np.ones(n) if model == "configuration" else rng.choice([0.3, 0.5, 0.6, 1.0], n)
            site = np.array([f"S{k}" for k in range(n)], dtype=object)
            for periodic in (True, False):
                ref = _brute_min_overlap(model, frac, occ, site, lat, periodic, 0.1)
                try:
                    _check_overlap(model, frac, occ, site, np.zeros((n, 3), int), lat, (periodic,) * 3, 0.1)
                    got = None
                except ValueError as err:
                    got = float(str(err).split(" are ")[1].split(" Å")[0])
                assert (ref is None) == (got is None), (p, model, periodic, ref, got)
                if ref is not None:
                    assert abs(ref - got) < 1e-6


def test_distance_matrix_matches_pairwise_distance():
    lat = Lattice.from_parameters(5.0, 5.0, 5.0, 90.0, 90.0, 20.0)
    rng = np.random.default_rng(3)
    n = 60
    frac = rng.uniform(-1, 2, (n, 3))
    s = Structure(lat, np.array(["Cu"] * n, dtype=object), frac, model="configuration",
                  atom_id=np.arange(100, 100 + n), min_separation=None)
    D = s.distance_matrix()
    assert D.shape == (n, n) and np.allclose(D, D.T) and np.all(np.diag(D) == 0)
    for a in range(0, n, 7):
        for b in range(a + 1, n, 5):
            assert D[a, b] == pytest.approx(s.distance(100 + a, 100 + b), abs=1e-12)
    sub = s.distance_matrix(ids=[130, 101, 155])
    assert sub[0, 2] == pytest.approx(s.distance(130, 155), abs=1e-12)


def test_scale_1k_atoms_with_separation():
    """Construct 1 000 atoms with the separation check on."""
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
    other = s.replace(pbc=(True, True, False), min_separation=None)
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


# ------------------------------ review fix-ups: Task 2a ------------------------------
def test_edits_recheck_overlap():
    s = _two_atom_base()
    with pytest.raises(ValueError, match=r"site 'S1' in cell \(0, 0, 0\)"):
        s.with_frac(np.array([[0.25, 0.25, 0.25], [0.2501, 0.25, 0.25]]))
    with pytest.raises(ValueError, match="apart"):
        s.with_cartesian(np.array([[1.25, 1.25, 1.25], [1.25, 1.25, 1.26]]))
    with pytest.raises(ValueError, match="apart"):
        s.replace(frac=np.array([[0.25, 0.25, 0.25], [0.25, 0.25, 0.25]]))
    # explicit opt-out still works
    t = s.with_frac(np.array([[0.25, 0.25, 0.25], [0.2501, 0.25, 0.25]]), min_separation=None)
    assert t.n_atoms == 2
    # the original is untouched
    assert np.array_equal(s.frac, [[0.25, 0.25, 0.25], [0.75, 0.75, 0.75]])


def test_concat_validates_result():
    s = _two_atom_base()
    # configuration concatenated with itself: two atoms per site
    with pytest.raises(ValueError, match="configuration"):
        s.concat(s, renumber=True)
    # average model concatenated with itself: site occupancy sum 2.0
    lat = _cubic(4.0)
    avg = Structure(lat, np.array(["Cu", "Au"]),
                    np.array([[0.25, 0.25, 0.25], [0.25, 0.25, 0.25]]),
                    model="average", occupancy=np.array([0.75, 0.25]),
                    site_label=np.array(["S1", "S1"]))
    with pytest.raises(ValueError, match="occupancy"):
        avg.concat(avg, renumber=True)


def test_concat_rejects_mismatched_metadata():
    lat = _cubic(5.0)
    a = Structure(lat, np.array(["Ga"]), np.array([[0.1, 0.1, 0.1]]),
                  model="configuration", atom_id=np.array([0]), site_label=np.array(["A1"]))
    b_iso = Structure(lat, np.array(["As"]), np.array([[0.6, 0.6, 0.6]]),
                      model="configuration", atom_id=np.array([1]), site_label=np.array(["B1"]),
                      b_iso=np.array([0.5]))
    with pytest.raises(ValueError, match="b_iso"):
        a.concat(b_iso)
    big = Structure(lat, np.array(["As"]), np.array([[0.6, 0.6, 0.6]]),
                    model="configuration", atom_id=np.array([1]), site_label=np.array(["B1"]),
                    supercell=(2, 2, 2))
    with pytest.raises(ValueError, match="supercell"):
        a.concat(big)
    slab = Structure(lat, np.array(["As"]), np.array([[0.6, 0.6, 0.6]]),
                     model="configuration", atom_id=np.array([1]), site_label=np.array(["B1"]),
                     pbc=(True, True, False), min_separation=None)
    with pytest.raises(ValueError, match="pbc"):
        a.concat(slab)
    # compatible structures concatenate, keeping ids and order
    c = Structure(lat, np.array(["As"]), np.array([[0.6, 0.6, 0.6]]),
                  model="configuration", atom_id=np.array([7]), site_label=np.array(["B1"]))
    ab = a.concat(c)
    assert list(ab.atom_id) == [0, 7]
    assert list(ab.species) == ["Ga", "As"]


def test_site_consistency_checks_the_right_rows():
    """Mixed sites are validated using their own rows, whatever order the site keys sort in."""
    lat = _cubic(4.0)
    # S2 (mixed) sorts after S1 (single): its two rows sit at different positions -> must raise
    with pytest.raises(ValueError, match="S2.*differing fractional positions"):
        Structure(lat, np.array(["X", "Cu", "Au"]),
                  np.array([[0.0, 0.0, 0.0], [0.25, 0.25, 0.25], [0.5, 0.5, 0.5]]),
                  model="average", occupancy=np.array([1.0, 0.5, 0.5]),
                  site_label=np.array(["S1", "S2", "S2"]))
    # two mixed sites after a single site: the over-occupied one (S3) is reported
    with pytest.raises(ValueError, match="S3.*total occupancy"):
        Structure(lat, np.array(["X", "Cu", "Au", "Cu", "Au"]),
                  np.array([[0, 0, 0], [0.25, 0.25, 0.25], [0.25, 0.25, 0.25],
                            [0.5, 0.5, 0.5], [0.5, 0.5, 0.5]], dtype=float),
                  model="average", occupancy=np.array([1.0, 0.5, 0.5, 0.7, 0.7]),
                  site_label=np.array(["S1", "S2", "S2", "S3", "S3"]))
    # same structure with valid occupancies is accepted
    ok = Structure(lat, np.array(["X", "Cu", "Au", "Cu", "Au"]),
                   np.array([[0, 0, 0], [0.25, 0.25, 0.25], [0.25, 0.25, 0.25],
                             [0.5, 0.5, 0.5], [0.5, 0.5, 0.5]], dtype=float),
                   model="average", occupancy=np.array([1.0, 0.5, 0.5, 0.6, 0.4]),
                   site_label=np.array(["S1", "S2", "S2", "S3", "S3"]))
    assert ok.n_atoms == 5


@pytest.mark.parametrize("label, z", [("O2-", 8), ("Fe3+", 26), ("Na+", 11), ("Cl-", 17), ("O-2", 8), ("Fe1", 26),
                                      ("D", 1), ("T", 1), ("CU", 29), ("Cu", 29), ("cu", 29)])
def test_atomic_number_strips_charges_and_reads_isotopes(label, z):
    from discuspdf._elements import _strip_suffix, atomic_number
    assert atomic_number(label) == z
    if label in ("O2-", "Na+", "Fe3+"):
        assert _strip_suffix(label) == {"O2-": "O", "Na+": "Na", "Fe3+": "Fe"}[label]
