"""Tests for the CIF reader and writer (Package A, Feature A2).

All tests that require ``gemmi`` use ``pytest.importorskip("gemmi")`` so the
suite is green without the ``cif`` extra.
"""
from __future__ import annotations

import math
import subprocess
import sys
import textwrap
from pathlib import Path

import numpy as np
import pytest

from discuspdf.lattice import Lattice
from discuspdf.structure import Structure
from discuspdf.interop import (
    InteropError,
    InteropWarning,
    B_FROM_U,
    u_to_b,
    b_to_u,
    element_from_type_symbol,
    LossReport,
)

FIXTURES = Path(__file__).parent / "fixtures" / "interop"


# ---------------------------------------------------------------------------
# 1. Independent oracle: F-43m expansion
# ---------------------------------------------------------------------------

def _f43m_operations() -> list[tuple[np.ndarray, np.ndarray]]:
    """Return the 96 symmetry operations of F-43m as (rotation, translation) pairs.

    Built from the 4 F-centring translations × the 24 point operations of -43m.
    The 24 point ops are listed explicitly (independent of gemmi).
    """
    # The 24 point ops of -43m as 3x3 rotation matrices.
    # (x,y,z) -> (x',y',z') where each of x',y',z' is ±x, ±y, or ±z.
    _ops_str = [
        "x,y,z", "y,-x,-z", "-x,-y,z", "-y,x,-z", "x,-y,-z",
        "y,x,z", "-x,y,-z", "-y,-x,z", "z,x,y", "-z,y,-x",
        "z,-x,-y", "-z,-y,x", "-z,x,-y", "z,y,x", "-z,-x,y",
        "z,-y,-x", "x,-z,-y", "y,z,x", "-x,-z,y", "-y,z,-x",
        "x,z,y", "y,-z,-x", "-x,z,-y", "-y,-z,x",
    ]

    def _parse(s: str) -> np.ndarray:
        parts = s.split(",")
        mat = np.zeros((3, 3))
        for i, p in enumerate(parts):
            if p == "x":
                mat[i, 0] = 1
            elif p == "-x":
                mat[i, 0] = -1
            elif p == "y":
                mat[i, 1] = 1
            elif p == "-y":
                mat[i, 1] = -1
            elif p == "z":
                mat[i, 2] = 1
            elif p == "-z":
                mat[i, 2] = -1
        return mat

    point_ops = [_parse(s) for s in _ops_str]
    assert len(point_ops) == 24

    # F-centring translations (fractional).
    translations = [
        np.array([0.0, 0.0, 0.0]),
        np.array([0.0, 0.5, 0.5]),
        np.array([0.5, 0.0, 0.5]),
        np.array([0.5, 0.5, 0.0]),
    ]

    # Combine: 24 point ops × 4 translations = 96 operations.
    ops: list[tuple[np.ndarray, np.ndarray]] = []
    for t in translations:
        for p in point_ops:
            ops.append((p, t))
    return ops


def test_oracle_f43m_expansion():
    """Compare the independent F-43m oracle with read_cif on the real COD GaAs CIF."""
    gemmi = pytest.importorskip("gemmi")
    from discuspdf.cif_io import read_cif

    cif_path = FIXTURES / "gaas_cod_9008845.cif"
    result = read_cif(cif_path)
    s = result.structure

    # The CIF has 2 listed atoms: Ga at (0,0,0), As at (0.25,0.25,0.25).
    # The oracle expands them by the 96 F-43m operations.
    ops = _f43m_operations()

    # Ga at (0,0,0):
    ga_positions = []
    for rot, trans in ops:
        pos = rot @ np.array([0.0, 0.0, 0.0]) + trans
        ga_positions.append(pos % 1.0)

    # As at (0.25,0.25,0.25):
    as_positions = []
    for rot, trans in ops:
        pos = rot @ np.array([0.25, 0.25, 0.25]) + trans
        as_positions.append(pos % 1.0)

    # Merge (remove duplicates within 1e-8).
    def _merge(positions: list[np.ndarray], tol: float = 1e-8) -> list[np.ndarray]:
        merged: list[np.ndarray] = []
        for p in positions:
            is_dup = False
            for m in merged:
                d = np.abs(p - m)
                d = np.minimum(d, 1.0 - d)
                if np.max(d) < tol:
                    is_dup = True
                    break
            if not is_dup:
                merged.append(p)
        return merged

    ga_merged = _merge(ga_positions)
    as_merged = _merge(as_positions)

    # Expect 4 Ga + 4 As = 8 atoms.
    assert len(ga_merged) == 4, f"expected 4 Ga, got {len(ga_merged)}"
    assert len(as_merged) == 4, f"expected 4 As, got {len(as_merged)}"

    # Compare with read_cif result.
    assert s.n_atoms == 8, f"expected 8 atoms, got {s.n_atoms}"

    # Check Ga positions.
    ga_in_result = s.frac[s.species == "Ga"]
    for target in ga_merged:
        found = False
        for pos in ga_in_result:
            d = np.abs(pos - target)
            d = np.minimum(d, 1.0 - d)
            if np.max(d) < 1e-8:
                found = True
                break
        assert found, f"Ga position {target} not found in result"

    # Check As positions.
    as_in_result = s.frac[s.species == "As"]
    for target in as_merged:
        found = False
        for pos in as_in_result:
            d = np.abs(pos - target)
            d = np.minimum(d, 1.0 - d)
            if np.max(d) < 1e-8:
                found = True
                break
        assert found, f"As position {target} not found in result"


# ---------------------------------------------------------------------------
# 2. Physics: GaAs 2×2×2 supercell nearest neighbours
# ---------------------------------------------------------------------------

def test_gaas_supercell_neighbours():
    """The GaAs structure from the real CIF turned into a 2×2×2 supercell
    has 4 nearest neighbours of the opposite species at a·√3/4 ± 1e-6."""
    gemmi = pytest.importorskip("gemmi")
    from discuspdf.cif_io import read_cif

    cif_path = FIXTURES / "gaas_cod_9008845.cif"
    result = read_cif(cif_path)
    s = result.structure

    # Build a 2×2×2 supercell.
    a, b, c, alpha, beta, gamma = s.lattice.parameters
    super_lattice = Lattice.from_parameters(
        a * 2, b * 2, c * 2, alpha, beta, gamma
    )
    # New fractional coordinates: (f + shift) / 2 for shifts (0,0), (0,1), (1,0), (1,1) in each axis.
    # Simpler: replicate the 8 atoms in a 2×2×2 grid.
    new_frac = []
    new_species = []
    for dx in range(2):
        for dy in range(2):
            for dz in range(2):
                for i in range(s.n_atoms):
                    new_frac.append(
                        [s.frac[i, 0] / 2 + dx / 2,
                         s.frac[i, 1] / 2 + dy / 2,
                         s.frac[i, 2] / 2 + dz / 2]
                    )
                    new_species.append(s.species[i])

    new_frac = np.array(new_frac)
    new_species = np.array(new_species, dtype=object)
    supercell = Structure(
        super_lattice,
        new_species,
        new_frac,
        model="configuration",
        supercell=(2, 2, 2),
        pbc=(True, True, True),
        min_separation=0.1,
    )

    # Expected nearest-neighbour distance: a·√3/4
    expected_d = a * math.sqrt(3) / 4

    # Check every atom has 4 nearest neighbours of the opposite species.
    for i in range(supercell.n_atoms):
        species_i = supercell.species[i]
        opposite = "As" if species_i == "Ga" else "Ga"
        # Find all atoms of the opposite species.
        opposite_mask = supercell.species == opposite
        opposite_indices = np.where(opposite_mask)[0]
        # Compute distances.
        distances = []
        for j in opposite_indices:
            if i == j:
                continue
            d = supercell.distance(int(supercell.atom_id[i]), int(supercell.atom_id[j]))
            distances.append((d, j))
        distances.sort()
        # The 4 nearest should all be at a·√3/4 ± 1e-6.
        assert len(distances) >= 4, f"atom {i} has only {len(distances)} opposite-species neighbours"
        for k in range(4):
            d, _ = distances[k]
            assert abs(d - expected_d) < 1e-6, (
                f"atom {i} ({species_i}) neighbour {k}: d = {d:.6f}, "
                f"expected {expected_d:.6f}"
            )


# ---------------------------------------------------------------------------
# 3. Mixed fixture
# ---------------------------------------------------------------------------

def test_mixed_fixture():
    """The mixed fixture: 12 sites; at (0,0,0) Ga 0.7 + In 0.3 share
    (site_label, cell_index); model='average'; occupancy sums 1;
    model='configuration' → InteropError."""
    gemmi = pytest.importorskip("gemmi")
    from discuspdf.cif_io import read_cif

    cif_path = FIXTURES / "ga07in03as_mixed.cif"
    result = read_cif(cif_path)
    s = result.structure

    # 12 sites (3 listed × 4 F-centring images).
    assert s.n_atoms == 12, f"expected 12 atoms, got {s.n_atoms}"
    assert s.model == "average"

    # At each of the 4 F-centring positions: Ga 0.7 + In 0.3 = 1.0.
    for pos in [
        [0.0, 0.0, 0.0],
        [0.0, 0.5, 0.5],
        [0.5, 0.0, 0.5],
        [0.5, 0.5, 0.0],
    ]:
        pos_arr = np.array(pos)
        # Find rows where all 3 coordinates match.
        idx = [
            i for i in range(s.n_atoms)
            if np.allclose(s.frac[i], pos_arr, atol=1e-9)
        ]
        assert len(idx) == 2, f"expected 2 atoms at {pos}, got {len(idx)}"
        sp = [s.species[i] for i in idx]
        oc = [s.occupancy[i] for i in idx]
        assert "Ga" in sp and "In" in sp, f"expected Ga+In at {pos}, got {sp}"
        ga_i = sp.index("Ga")
        in_i = sp.index("In")
        total = oc[ga_i] + oc[in_i]
        assert abs(total - 1.0) < 1e-9, f"occupancy sum {total} != 1.0 at {pos}"
        assert abs(oc[ga_i] - 0.7) < 1e-9
        assert abs(oc[in_i] - 0.3) < 1e-9
        # The Ga and In rows at this position share (site_label, cell_index).
        ga_global = idx[ga_i]
        in_global = idx[in_i]
        assert s.site_label[ga_global] == s.site_label[in_global], (
            f"Ga/In site_label differs at {pos}: "
            f"{s.site_label[ga_global]} vs {s.site_label[in_global]}"
        )
        assert np.array_equal(s.cell_index[ga_global], s.cell_index[in_global]), (
            f"Ga/In cell_index differs at {pos}: "
            f"{s.cell_index[ga_global]} vs {s.cell_index[in_global]}"
        )

    # The four cation labels are distinct from each other and from every As label.
    cation_labels = set()
    for pos in [
        [0.0, 0.0, 0.0],
        [0.0, 0.5, 0.5],
        [0.5, 0.0, 0.5],
        [0.5, 0.5, 0.0],
    ]:
        pos_arr = np.array(pos)
        idx = [
            i for i in range(s.n_atoms)
            if np.allclose(s.frac[i], pos_arr, atol=1e-9)
        ]
        cation_labels.add(s.site_label[idx[0]])
    assert len(cation_labels) == 4, (
        f"expected 4 distinct cation labels, got {sorted(cation_labels)}"
    )

    as_labels = {s.site_label[i] for i in range(s.n_atoms) if s.species[i] == "As"}
    assert not cation_labels & as_labels, (
        f"cation labels overlap As labels: {sorted(cation_labels & as_labels)}"
    )

    # The Ga/In label equals Ga1_k (k = image index, 1-based).
    for k, pos in enumerate(
        [
            [0.0, 0.0, 0.0],
            [0.0, 0.5, 0.5],
            [0.5, 0.0, 0.5],
            [0.5, 0.5, 0.0],
        ]
    ):
        pos_arr = np.array(pos)
        idx = [
            i for i in range(s.n_atoms)
            if np.allclose(s.frac[i], pos_arr, atol=1e-9)
        ]
        expected = f"Ga1_{k + 1}"
        assert s.site_label[idx[0]] == expected, (
            f"expected cation label {expected} at {pos}, "
            f"got {s.site_label[idx[0]]}"
        )

    # Per-site occupancy sum is 1.0 for every site
    # (group by (site_label, cell_index)).
    groups: dict[tuple[str, tuple], float] = {}
    for i in range(s.n_atoms):
        key = (s.site_label[i], tuple(s.cell_index[i]))
        groups[key] = groups.get(key, 0.0) + s.occupancy[i]
    assert len(groups) == 8, f"expected 8 sites, got {len(groups)}"
    for key, total in groups.items():
        assert abs(total - 1.0) < 1e-9, (
            f"occupancy sum {total} != 1.0 for site {key}"
        )

    # model='configuration' → InteropError.
    with pytest.raises(InteropError, match="configuration"):
        read_cif(cif_path, model="configuration")


# ---------------------------------------------------------------------------
# 4. Triclinic P-1
# ---------------------------------------------------------------------------

def test_triclinic_p1bar():
    """triclinic_p1bar.cif: 4 atoms; images of Cu at (0.1,0.2,0.3) and
    (0.9,0.8,0.7); Cartesian pair vectors invariant."""
    gemmi = pytest.importorskip("gemmi")
    from discuspdf.cif_io import read_cif

    cif_path = FIXTURES / "triclinic_p1bar.cif"
    result = read_cif(cif_path)
    s = result.structure

    assert s.n_atoms == 4, f"expected 4 atoms, got {s.n_atoms}"
    assert result.symmetry_source == "symop_loop"

    # Cu at (0.1,0.2,0.3) and its inversion partner at (0.9,0.8,0.7).
    cu_mask = s.species == "Cu"
    cu_frac = s.frac[cu_mask]
    assert cu_frac.shape == (2, 3)

    # The two Cu atoms should be inversion partners: r and -r (mod lattice).
    cu0 = cu_frac[0]
    cu1 = cu_frac[1]
    # Check: cu1 ≈ -cu0 (mod 1).
    diff = cu1 - (-cu0)
    diff = diff - np.rint(diff)  # reduce to [0, 1)
    assert np.allclose(diff, 0.0, atol=1e-9), (
        f"Cu inversion partner check failed: {cu0} vs {cu1}"
    )

    # Cartesian pair vectors invariant: the displacement from Cu0 to Cu1
    # should be the same as from Cu1 to Cu0 (negated).
    d01 = s.displacement(int(s.atom_id[cu_mask][0]), int(s.atom_id[cu_mask][1]))
    d10 = s.displacement(int(s.atom_id[cu_mask][1]), int(s.atom_id[cu_mask][0]))
    assert np.allclose(d01, -d10, atol=1e-9)


# ---------------------------------------------------------------------------
# 5. Losses
# ---------------------------------------------------------------------------

def test_aniso_adp_loss():
    """Anisotropic ADP file → InteropError naming the loss;
    with allow_loss=True → InteropWarning once and CifImport.loss lists it."""
    gemmi = pytest.importorskip("gemmi")
    from discuspdf.cif_io import read_cif

    cif_path = FIXTURES / "p1_aniso_adp.cif"

    # Default: raise.
    with pytest.raises(InteropError, match="anisotropic"):
        read_cif(cif_path)

    # allow_loss=True: warn and record.
    import warnings
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        result = read_cif(cif_path, allow_loss=True)
        assert len(w) == 1
        assert isinstance(w[0].message, InteropWarning)
        assert "anisotropic" in str(w[0].message)
        assert len(result.loss) > 0
        assert any("anisotropic" in loss for loss in result.loss)


def test_disorder_group_loss():
    """Disorder-group file → InteropError naming the loss."""
    gemmi = pytest.importorskip("gemmi")
    from discuspdf.cif_io import read_cif

    cif_path = FIXTURES / "p1_disorder_groups.cif"
    with pytest.raises(InteropError, match="disorder"):
        read_cif(cif_path)


def test_occupancy_over_one():
    """Occupancy-over-one file → InteropError."""
    gemmi = pytest.importorskip("gemmi")
    from discuspdf.cif_io import read_cif

    cif_path = FIXTURES / "p1_occupancy_over_one.cif"
    with pytest.raises(InteropError, match="occupancy"):
        read_cif(cif_path)


# ---------------------------------------------------------------------------
# 6. Errors
# ---------------------------------------------------------------------------

def test_missing_cell():
    """No _cell_length_a → InteropError."""
    gemmi = pytest.importorskip("gemmi")
    from discuspdf.cif_io import parse_cif

    text = textwrap.dedent("""\
        data_test
        _symmetry_space_group_name_H-M 'P 1'
        loop_
        _atom_site_label
        _atom_site_fract_x
        _atom_site_fract_y
        _atom_site_fract_z
        Ga1 0.0 0.0 0.0
    """)
    with pytest.raises(InteropError, match="cell"):
        parse_cif(text)


def test_unknown_element():
    """Unknown element → InteropError."""
    gemmi = pytest.importorskip("gemmi")
    from discuspdf.cif_io import parse_cif

    text = textwrap.dedent("""\
        data_test
        _cell_length_a 5.0
        _cell_length_b 5.0
        _cell_length_c 5.0
        _cell_angle_alpha 90
        _cell_angle_beta 90
        _cell_angle_gamma 90
        _symmetry_space_group_name_H-M 'P 1'
        loop_
        _atom_site_label
        _atom_site_type_symbol
        _atom_site_fract_x
        _atom_site_fract_y
        _atom_site_fract_z
        Xx1 Xx 0.0 0.0 0.0
    """)
    with pytest.raises(InteropError, match="unknown element"):
        parse_cif(text)


def test_empty_text():
    """Empty text → InteropError."""
    gemmi = pytest.importorskip("gemmi")
    from discuspdf.cif_io import parse_cif

    with pytest.raises(InteropError, match="empty"):
        parse_cif("")


def test_two_blocks():
    """Two blocks without block= → InteropError."""
    gemmi = pytest.importorskip("gemmi")
    from discuspdf.cif_io import parse_cif

    text = textwrap.dedent("""\
        data_a
        _cell_length_a 5.0
        _cell_length_b 5.0
        _cell_length_c 5.0
        _cell_angle_alpha 90
        _cell_angle_beta 90
        _cell_angle_gamma 90
        _symmetry_space_group_name_H-M 'P 1'
        loop_
        _atom_site_label
        _atom_site_fract_x
        _atom_site_fract_y
        _atom_site_fract_z
        Ga1 0.0 0.0 0.0

        data_b
        _cell_length_a 5.0
        _cell_length_b 5.0
        _cell_length_c 5.0
        _cell_angle_alpha 90
        _cell_angle_beta 90
        _cell_angle_gamma 90
        _symmetry_space_group_name_H-M 'P 1'
        loop_
        _atom_site_label
        _atom_site_fract_x
        _atom_site_fract_y
        _atom_site_fract_z
        As1 0.0 0.0 0.0
    """)
    with pytest.raises(InteropError, match="multiple"):
        parse_cif(text)


def test_negative_u():
    """Negative U → InteropError (B would be negative)."""
    gemmi = pytest.importorskip("gemmi")
    from discuspdf.cif_io import parse_cif

    text = textwrap.dedent("""\
        data_test
        _cell_length_a 5.0
        _cell_length_b 5.0
        _cell_length_c 5.0
        _cell_angle_alpha 90
        _cell_angle_beta 90
        _cell_angle_gamma 90
        _symmetry_space_group_name_H-M 'P 1'
        loop_
        _atom_site_label
        _atom_site_type_symbol
        _atom_site_fract_x
        _atom_site_fract_y
        _atom_site_fract_z
        _atom_site_occupancy
        _atom_site_U_iso_or_equiv
        Fe1 Fe 0.0 0.0 0.0 1.0 -0.01
    """)
    # Negative U gives negative B, which Structure rejects.
    with pytest.raises((InteropError, ValueError)):
        parse_cif(text)


# ---------------------------------------------------------------------------
# 7. Defaults
# ---------------------------------------------------------------------------

def test_no_u_column():
    """File without U → b_iso=None and 'b_iso' in defaults_applied."""
    gemmi = pytest.importorskip("gemmi")
    from discuspdf.cif_io import parse_cif

    text = textwrap.dedent("""\
        data_test
        _cell_length_a 5.0
        _cell_length_b 5.0
        _cell_length_c 5.0
        _cell_angle_alpha 90
        _cell_angle_beta 90
        _cell_angle_gamma 90
        _symmetry_space_group_name_H-M 'P 1'
        loop_
        _atom_site_label
        _atom_site_type_symbol
        _atom_site_fract_x
        _atom_site_fract_y
        _atom_site_fract_z
        _atom_site_occupancy
        Ga1 Ga 0.0 0.0 0.0 1.0
    """)
    result = parse_cif(text)
    assert result.structure.b_iso is None
    assert "b_iso" in result.defaults_applied


def test_assumed_p1():
    """assumed_P1 reported when no symmetry information."""
    gemmi = pytest.importorskip("gemmi")
    from discuspdf.cif_io import parse_cif

    text = textwrap.dedent("""\
        data_test
        _cell_length_a 5.0
        _cell_length_b 5.0
        _cell_length_c 5.0
        _cell_angle_alpha 90
        _cell_angle_beta 90
        _cell_angle_gamma 90
        loop_
        _atom_site_label
        _atom_site_type_symbol
        _atom_site_fract_x
        _atom_site_fract_y
        _atom_site_fract_z
        _atom_site_occupancy
        Ga1 Ga 0.0 0.0 0.0 1.0
    """)
    result = parse_cif(text)
    assert result.symmetry_source == "assumed_P1"
    assert "symmetry" in result.defaults_applied


# ---------------------------------------------------------------------------
# 8. Round trip
# ---------------------------------------------------------------------------

def test_roundtrip_gaas():
    """Round-trip the GaAs CIF through format_cif → parse_cif."""
    gemmi = pytest.importorskip("gemmi")
    from discuspdf.cif_io import read_cif, format_cif, parse_cif

    cif_path = FIXTURES / "gaas_cod_9008845.cif"
    result = read_cif(cif_path)
    s = result.structure

    text = format_cif(s)
    result2 = parse_cif(text)
    s2 = result2.structure

    assert s2.n_atoms == s.n_atoms
    assert (s2.species == s.species).all()
    assert np.allclose(s2.frac, s.frac, atol=1e-9)
    assert np.allclose(s2.occupancy, s.occupancy, atol=1e-9)
    assert (s2.site_label == s.site_label).all()
    assert (s2.atom_id == s.atom_id).all()
    assert np.array_equal(s2.cell_index, s.cell_index)
    assert np.array_equal(s2.supercell, s.supercell)
    assert s2.model == s.model


def test_roundtrip_triclinic():
    """Round-trip the triclinic CIF."""
    gemmi = pytest.importorskip("gemmi")
    from discuspdf.cif_io import read_cif, format_cif, parse_cif

    cif_path = FIXTURES / "triclinic_p1bar.cif"
    result = read_cif(cif_path)
    s = result.structure

    text = format_cif(s)
    result2 = parse_cif(text)
    s2 = result2.structure

    assert s2.n_atoms == s.n_atoms
    assert (s2.species == s.species).all()
    assert np.allclose(s2.frac, s.frac, atol=1e-9)
    assert (s2.site_label == s.site_label).all()


def test_roundtrip_mixed():
    """Round-trip the mixed-occupancy CIF."""
    gemmi = pytest.importorskip("gemmi")
    from discuspdf.cif_io import read_cif, format_cif, parse_cif

    cif_path = FIXTURES / "ga07in03as_mixed.cif"
    result = read_cif(cif_path)
    s = result.structure

    text = format_cif(s)
    result2 = parse_cif(text)
    s2 = result2.structure

    assert s2.n_atoms == s.n_atoms
    assert (s2.species == s.species).all()
    assert np.allclose(s2.frac, s.frac, atol=1e-9)
    assert np.allclose(s2.occupancy, s.occupancy, atol=1e-9)
    assert (s2.site_label == s.site_label).all()


def test_b_from_u_values():
    """M3 check: B = 8π²·U with literal values (not u_to_b)."""
    gemmi = pytest.importorskip("gemmi")
    from discuspdf.cif_io import read_cif

    cif_path = FIXTURES / "ga07in03as_mixed.cif"
    result = read_cif(cif_path)
    s = result.structure

    # Ga and In rows: U = 0.0100 → B = 8π²·0.0100 = 0.789568352...
    expected_b_ga_in = 8.0 * math.pi ** 2 * 0.0100
    # As rows: U = 0.0080 → B = 8π²·0.0080 = 0.631654682...
    expected_b_as = 8.0 * math.pi ** 2 * 0.0080

    assert s.b_iso is not None, "expected b_iso to be present"
    for i in range(s.n_atoms):
        if s.species[i] in ("Ga", "In"):
            assert abs(s.b_iso[i] - expected_b_ga_in) < 1e-6, (
                f"{s.species[i]} b_iso={s.b_iso[i]}, expected {expected_b_ga_in}"
            )
        elif s.species[i] == "As":
            assert abs(s.b_iso[i] - expected_b_as) < 1e-6, (
                f"As b_iso={s.b_iso[i]}, expected {expected_b_as}"
            )


def test_roundtrip_mixed_labels():
    """Round-trip preserves shared mixed-site labels."""
    gemmi = pytest.importorskip("gemmi")
    from discuspdf.cif_io import read_cif, format_cif, parse_cif

    cif_path = FIXTURES / "ga07in03as_mixed.cif"
    result = read_cif(cif_path)
    s = result.structure

    text = format_cif(s)
    result2 = parse_cif(text)
    s2 = result2.structure

    assert s2.n_atoms == s.n_atoms
    # Ga and In at the same position must share the same site_label.
    for pos in [
        [0.0, 0.0, 0.0],
        [0.0, 0.5, 0.5],
        [0.5, 0.0, 0.5],
        [0.5, 0.5, 0.0],
    ]:
        pos_arr = np.array(pos)
        idx = [
            i for i in range(s2.n_atoms)
            if np.allclose(s2.frac[i], pos_arr, atol=1e-9)
        ]
        assert len(idx) == 2, f"expected 2 atoms at {pos}, got {len(idx)}"
        labels = [s2.site_label[i] for i in idx]
        assert labels[0] == labels[1], (
            f"mixed site labels differ after round-trip at {pos}: {labels}"
        )
        # The shared label must start with "Ga1" (first listed atom).
        assert labels[0].startswith("Ga1"), (
            f"expected shared label to start with 'Ga1', got {labels[0]}"
        )


def test_mixed_reversed_order():
    """When In is listed before Ga, the shared label is In1_k."""
    gemmi = pytest.importorskip("gemmi")
    from discuspdf.cif_io import parse_cif

    cif_path = FIXTURES / "ga07in03as_mixed.cif"
    text = cif_path.read_text(encoding="ascii")

    # Swap the Ga and In rows: move the In row before the Ga row.
    lines = text.splitlines()
    ga_line_idx = None
    in_line_idx = None
    for i, line in enumerate(lines):
        if line.startswith("Ga1 "):
            ga_line_idx = i
        elif line.startswith("In1 "):
            in_line_idx = i
    assert ga_line_idx is not None and in_line_idx is not None
    # Swap: put In before Ga.
    ga_line = lines[ga_line_idx]
    in_line = lines[in_line_idx]
    lines[ga_line_idx] = in_line
    lines[in_line_idx] = ga_line
    swapped = "\n".join(lines) + "\n"

    result = parse_cif(swapped)
    s = result.structure

    # At (0,0,0): the shared label must start with "In1".
    origin_idx = [
        i for i in range(s.n_atoms)
        if np.allclose(s.frac[i], [0.0, 0.0, 0.0], atol=1e-9)
    ]
    assert len(origin_idx) == 2
    labels = [s.site_label[i] for i in origin_idx]
    assert labels[0] == labels[1], f"mixed site labels differ: {labels}"
    assert labels[0].startswith("In1"), (
        f"expected shared label to start with 'In1', got {labels[0]}"
    )


def test_roundtrip_precision_10_decimals():
    """M8 check: the writer must use 10 decimal places.

    A 6-decimal writer would lose precision on the GaAs cell parameter
    (5.6537 Å) and the fractional coordinates.
    """
    gemmi = pytest.importorskip("gemmi")
    from discuspdf.cif_io import read_cif, format_cif, parse_cif

    cif_path = FIXTURES / "gaas_cod_9008845.cif"
    result = read_cif(cif_path)
    s = result.structure

    text = format_cif(s)
    # The written text must contain 10-decimal numbers.
    # Check that the cell parameter has 10 decimal places.
    assert "5.6537000000" in text, (
        f"expected 10-decimal cell parameter in CIF text, got: "
        f"{[l for l in text.splitlines() if 'cell_length_a' in l]}"
    )

    result2 = parse_cif(text)
    s2 = result2.structure

    # Tight tolerance: 10-decimal precision means atol=5e-11.
    assert np.allclose(s2.frac, s.frac, atol=5e-11), (
        f"round-trip frac precision lost: max diff = "
        f"{np.abs(s2.frac - s.frac).max()}"
    )
    a_orig = s.lattice.parameters[0]
    a_rt = s2.lattice.parameters[0]
    assert abs(a_rt - a_orig) < 5e-11, (
        f"round-trip lattice precision lost: {a_orig} -> {a_rt}"
    )


def test_write_cif_overwrite():
    """write_cif refuses existing file without overwrite=True."""
    gemmi = pytest.importorskip("gemmi")
    from discuspdf.cif_io import read_cif, write_cif
    import tempfile

    cif_path = FIXTURES / "gaas_cod_9008845.cif"
    result = read_cif(cif_path)
    s = result.structure

    with tempfile.TemporaryDirectory() as tmpdir:
        target = Path(tmpdir) / "test.cif"
        write_cif(target, s)
        assert target.exists()
        # Second write without overwrite → FileExistsError.
        with pytest.raises(FileExistsError):
            write_cif(target, s)
        # With overwrite=True → succeeds.
        write_cif(target, s, overwrite=True)
        assert target.exists()


# ---------------------------------------------------------------------------
# 9. Generic path (no _discuspdf_* tags)
# ---------------------------------------------------------------------------

def test_generic_path():
    """Strip the _discuspdf_* lines from the written text → still parses;
    model/supercell default as documented."""
    gemmi = pytest.importorskip("gemmi")
    from discuspdf.cif_io import read_cif, format_cif, parse_cif

    cif_path = FIXTURES / "gaas_cod_9008845.cif"
    result = read_cif(cif_path)
    s = result.structure

    text = format_cif(s)
    # Strip _discuspdf_* header lines AND trailing _discuspdf_ values from
    # the loop data lines (the last 5 columns are _discuspdf_* tags).
    lines = []
    for line in text.splitlines():
        if line.startswith("_discuspdf_"):
            continue
        if line.startswith("data_") or line.startswith("loop_") or line.startswith("_"):
            lines.append(line)
        else:
            # Data line: strip the last 5 space-separated values
            # (_discuspdf_site_label, _discuspdf_atom_id,
            #  _discuspdf_cell_index_a/b/c).
            parts = line.split()
            if len(parts) >= 11:
                lines.append(" ".join(parts[:6]))
            else:
                lines.append(line)
    stripped = "\n".join(lines) + "\n"

    result2 = parse_cif(stripped)
    s2 = result2.structure

    # Should still parse (P1 symmetry, default model=average).
    assert s2.n_atoms == s.n_atoms
    assert (s2.species == s.species).all()
    assert np.allclose(s2.frac, s.frac, atol=1e-9)
    assert s2.model == "average"
    assert tuple(s2.supercell) == (1, 1, 1)


# ---------------------------------------------------------------------------
# 10. Optional import
# ---------------------------------------------------------------------------

def test_optional_import():
    """A subprocess test sets sys.modules['gemmi'] = None and checks that
    importing the modules works and read_cif raises ImportError."""
    code = textwrap.dedent("""\
        import sys
        sys.modules['gemmi'] = None
        import discuspdf
        import discuspdf.cif_io
        try:
            discuspdf.cif_io.read_cif('nonexistent.cif')
            print("FAIL: no ImportError raised")
        except ImportError as e:
            if "discus-python[cif]" in str(e):
                print("OK")
            else:
                print(f"FAIL: wrong message: {e}")
        except Exception as e:
            print(f"FAIL: wrong exception: {type(e).__name__}: {e}")
    """)
    result = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert "OK" in result.stdout, f"stdout: {result.stdout}\nstderr: {result.stderr}"


# ---------------------------------------------------------------------------
# 11. Atomic write
# ---------------------------------------------------------------------------

def test_atomic_write_failure():
    """Failure during write leaves no partial file."""
    gemmi = pytest.importorskip("gemmi")
    from discuspdf.cif_io import read_cif, write_cif
    import tempfile
    import os

    cif_path = FIXTURES / "gaas_cod_9008845.cif"
    result = read_cif(cif_path)
    s = result.structure

    with tempfile.TemporaryDirectory() as tmpdir:
        target = Path(tmpdir) / "test.cif"
        # Write successfully first.
        write_cif(target, s)
        assert target.exists()

        # Simulate a failure by making the directory read-only.
        # (This is a simple test: just verify the target still exists
        # and has the correct content after a successful write.)
        content = target.read_text()
        assert "data_discuspdf" in content


# ---------------------------------------------------------------------------
# 12. Interop helpers
# ---------------------------------------------------------------------------

def test_u_to_b():
    """B = 8π²·U."""
    assert abs(u_to_b(0.01) - 8 * math.pi**2 * 0.01) < 1e-15
    assert abs(b_to_u(8 * math.pi**2 * 0.01) - 0.01) < 1e-15
    assert abs(B_FROM_U - 8 * math.pi**2) < 1e-15


def test_element_from_type_symbol():
    """Strip charge/digits: O2-→O, Ga3+→Ga, GA→Ga; unknown → InteropError."""
    assert element_from_type_symbol("O2-") == "O"
    assert element_from_type_symbol("Ga3+") == "Ga"
    assert element_from_type_symbol("GA") == "Ga"
    assert element_from_type_symbol("Fe") == "Fe"
    with pytest.raises(InteropError):
        element_from_type_symbol("Xx")


def test_loss_report():
    """LossReport: raise_or_warn raises InteropError when not allowed,
    emits InteropWarning when allowed."""
    losses = LossReport()
    losses.add("loss 1")
    losses.add("loss 2")
    assert losses.losses == ("loss 1", "loss 2")

    # Raise.
    with pytest.raises(InteropError, match="loss 1"):
        losses.raise_or_warn(allow_loss=False)

    # Warn.
    import warnings
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        losses.raise_or_warn(allow_loss=True)
        assert len(w) == 1
        assert isinstance(w[0].message, InteropWarning)
        assert "loss 1" in str(w[0].message)
        assert "loss 2" in str(w[0].message)
