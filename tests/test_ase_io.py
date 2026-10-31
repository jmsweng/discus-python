"""Tests for the ASE adapter (Package A, Feature A2).

All tests that require ``ase`` use ``pytest.importorskip("ase")`` so the
suite is green without the ``ase`` extra.  ``gemmi`` is additionally
required to build the fixture structures via ``read_cif``.
"""
from __future__ import annotations

import pickle
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
)

FIXTURES = Path(__file__).parent / "fixtures" / "interop"


# ---------------------------------------------------------------------------
# (a) Round trip
# ---------------------------------------------------------------------------

def test_roundtrip_gaas():
    """GaAs CIF → to_ase → from_ase round trip preserves the structure."""
    ase = pytest.importorskip("ase")
    gemmi = pytest.importorskip("gemmi")
    from discuspdf.cif_io import read_cif
    from discuspdf.ase_io import to_ase, from_ase

    s = read_cif(FIXTURES / "gaas_cod_9008845.cif").structure
    atoms = to_ase(s)
    s2 = from_ase(atoms)

    assert s2.isclose(s, rtol=1e-8, atol=1e-12)
    assert s2.model == s.model
    assert s2.n_atoms == s.n_atoms


def test_roundtrip_mixed():
    """Mixed-occupancy CIF → to_ase → from_ase round trip preserves the structure."""
    ase = pytest.importorskip("ase")
    gemmi = pytest.importorskip("gemmi")
    from discuspdf.cif_io import read_cif
    from discuspdf.ase_io import to_ase, from_ase

    s = read_cif(FIXTURES / "ga07in03as_mixed.cif").structure
    atoms = to_ase(s)
    s2 = from_ase(atoms)

    assert s2.isclose(s, rtol=1e-8, atol=1e-12)
    assert s2.model == s.model
    assert s2.n_atoms == s.n_atoms


def test_roundtrip_triclinic():
    """Triclinic P-1 CIF → to_ase → from_ase round trip preserves the structure."""
    ase = pytest.importorskip("ase")
    gemmi = pytest.importorskip("gemmi")
    from discuspdf.cif_io import read_cif
    from discuspdf.ase_io import to_ase, from_ase

    s = read_cif(FIXTURES / "triclinic_p1bar.cif").structure
    atoms = to_ase(s)
    s2 = from_ase(atoms)

    assert s2.isclose(s, rtol=1e-8, atol=1e-12)
    assert s2.model == s.model
    assert s2.n_atoms == s.n_atoms


def test_roundtrip_b_iso_none_and_per_atom():
    """b_iso=None survives as None; per-atom B survives within 1e-6."""
    ase = pytest.importorskip("ase")
    gemmi = pytest.importorskip("gemmi")
    from discuspdf.cif_io import read_cif
    from discuspdf.ase_io import to_ase, from_ase

    # b_iso=None (the real CIF has U=0 → b_iso=None).
    s_none = read_cif(FIXTURES / "gaas_cod_9008845.cif").structure
    assert s_none.b_iso is None
    atoms = to_ase(s_none)
    s2 = from_ase(atoms)
    assert s2.b_iso is None

    # Per-atom B (the mixed fixture has B values from U).
    s_b = read_cif(FIXTURES / "ga07in03as_mixed.cif").structure
    assert s_b.b_iso is not None
    atoms = to_ase(s_b)
    s2 = from_ase(atoms)
    assert s2.b_iso is not None
    assert np.allclose(s2.b_iso, s_b.b_iso, atol=1e-6)


def test_roundtrip_after_copy_and_pickle():
    """The dp_* arrays and info survive atoms.copy() and pickle."""
    ase = pytest.importorskip("ase")
    gemmi = pytest.importorskip("gemmi")
    from discuspdf.cif_io import read_cif
    from discuspdf.ase_io import to_ase, from_ase

    s = read_cif(FIXTURES / "gaas_cod_9008845.cif").structure
    atoms = to_ase(s)

    # copy()
    atoms2 = atoms.copy()
    s2 = from_ase(atoms2)
    assert s2.isclose(s, rtol=1e-8, atol=1e-12)

    # pickle round trip
    atoms3 = pickle.loads(pickle.dumps(atoms))
    s3 = from_ase(atoms3)
    assert s3.isclose(s, rtol=1e-8, atol=1e-12)


# ---------------------------------------------------------------------------
# (b) Distance invariance
# ---------------------------------------------------------------------------

def test_ase_distance_invariance():
    """The multiset of pairwise minimum-image distances is preserved by the
    to_ase conversion (compare Structure.distance_matrix with ASE
    get_all_distances(mic=True))."""
    ase = pytest.importorskip("ase")
    gemmi = pytest.importorskip("gemmi")
    from discuspdf.cif_io import read_cif
    from discuspdf.ase_io import to_ase

    s = read_cif(FIXTURES / "gaas_cod_9008845.cif").structure
    atoms = to_ase(s)

    # (N, N) distance matrix from our Structure.
    d_ours = s.distance_matrix()
    np.fill_diagonal(d_ours, 0.0)

    # (N, N) distance matrix from ASE (minimum image).
    d_ase = atoms.get_all_distances(mic=True)
    np.fill_diagonal(d_ase, 0.0)

    # Both should have the same multiset of pairwise distances.
    flat_ours = np.sort(d_ours.flatten())
    flat_ase = np.sort(d_ase.flatten())
    assert np.allclose(flat_ours, flat_ase, atol=1e-6), (
        f"distance multisets differ: max diff = "
        f"{np.abs(flat_ours - flat_ase).max()}"
    )


# ---------------------------------------------------------------------------
# (c) Foreign object
# ---------------------------------------------------------------------------

def test_foreign_ase():
    """A foreign ase.Atoms (no dp_* metadata) converts with model given,
    occupancy 1.0, and the element symbols as species."""
    ase = pytest.importorskip("ase")
    gemmi = pytest.importorskip("gemmi")
    from discuspdf.ase_io import from_ase

    atoms = ase.Atoms(
        "GaAs",
        cell=[5.6537] * 3,
        pbc=True,
        positions=[[0.0, 0.0, 0.0], [0.25, 0.25, 0.25]],
    )
    s = from_ase(atoms, model="configuration")

    assert s.n_atoms == 2
    assert s.model == "configuration"
    assert np.allclose(s.occupancy, 1.0)
    assert list(s.species) == ["Ga", "As"]


def test_foreign_ase_partial_occupancy_loss():
    """info['occupancy'] with a value < 1 → InteropError by default;
    with allow_loss=True → InteropWarning and it proceeds."""
    ase = pytest.importorskip("ase")
    gemmi = pytest.importorskip("gemmi")
    from discuspdf.ase_io import from_ase

    atoms = ase.Atoms(
        "GaAs",
        cell=[5.6537] * 3,
        pbc=True,
        positions=[[0.0, 0.0, 0.0], [0.25, 0.25, 0.25]],
    )
    atoms.info['occupancy'] = {'0': {'Ga': 0.7}}

    # Default: raise.
    with pytest.raises(InteropError):
        from_ase(atoms, model="configuration")

    # allow_loss=True: warn and proceed.
    with pytest.warns(InteropWarning):
        s = from_ase(atoms, model="configuration", allow_loss=True)
    assert s.n_atoms == 2


# ---------------------------------------------------------------------------
# (d) Losses raise and report with allow_loss=True
# ---------------------------------------------------------------------------

def test_ase_pbc_mixed_loss():
    """to_ase on an average-model structure with partial occupancy always
    succeeds but emits an InteropWarning (ASE atoms are fully occupied)."""
    ase = pytest.importorskip("ase")
    gemmi = pytest.importorskip("gemmi")
    from discuspdf.cif_io import read_cif
    from discuspdf.ase_io import to_ase

    s = read_cif(FIXTURES / "ga07in03as_mixed.cif").structure
    assert s.model == "average"
    assert np.any(s.occupancy < 1.0)

    with pytest.warns(InteropWarning):
        atoms = to_ase(s)
    assert len(atoms) == s.n_atoms


# ---------------------------------------------------------------------------
# (e) Errors
# ---------------------------------------------------------------------------

def test_ase_singular_cell_error():
    """A zero cell (det = 0) → InteropError from from_ase."""
    ase = pytest.importorskip("ase")
    gemmi = pytest.importorskip("gemmi")
    from discuspdf.ase_io import from_ase

    atoms = ase.Atoms(
        "Ga",
        cell=[[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
        pbc=True,
        positions=[[0.0, 0.0, 0.0]],
    )
    with pytest.raises(InteropError):
        from_ase(atoms, model="configuration")


def test_ase_left_handed_cell_error():
    """A left-handed cell (det < 0) → InteropError from from_ase."""
    ase = pytest.importorskip("ase")
    gemmi = pytest.importorskip("gemmi")
    from discuspdf.ase_io import from_ase

    atoms = ase.Atoms(
        "Ga",
        cell=np.array([[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]]),
        pbc=True,
        positions=[[0.0, 0.0, 0.0]],
    )
    with pytest.raises(InteropError):
        from_ase(atoms, model="configuration")


def test_ase_foreign_missing_model_error():
    """A foreign Atoms object without dp_model and no model argument
    → InteropError mentioning 'model'."""
    ase = pytest.importorskip("ase")
    gemmi = pytest.importorskip("gemmi")
    from discuspdf.ase_io import from_ase

    atoms = ase.Atoms(
        "GaAs",
        cell=[5.6537] * 3,
        pbc=True,
        positions=[[0.0, 0.0, 0.0], [0.25, 0.25, 0.25]],
    )
    with pytest.raises(InteropError, match="model"):
        from_ase(atoms)


def test_ase_unknown_element_error():
    """A Structure with an unknown element symbol → InteropError from to_ase
    (element_from_type_symbol)."""
    ase = pytest.importorskip("ase")
    gemmi = pytest.importorskip("gemmi")
    from discuspdf.ase_io import to_ase

    lattice = Lattice.from_parameters(5.0, 5.0, 5.0, 90.0, 90.0, 90.0)
    s = Structure(
        lattice,
        np.array(["Xx"], dtype=object),
        np.array([[0.0, 0.0, 0.0]]),
        model="configuration",
    )
    with pytest.raises(InteropError):
        to_ase(s)


# ---------------------------------------------------------------------------
# (f) Optional import
# ---------------------------------------------------------------------------

def test_optional_import_ase():
    """A subprocess test sets sys.modules['ase'] = None and checks that
    importing the modules works and to_ase raises ImportError."""
    code = textwrap.dedent("""\
        import sys
        sys.modules['ase'] = None
        import discuspdf
        import discuspdf.ase_io
        from discuspdf.lattice import Lattice
        from discuspdf.structure import Structure
        try:
            s = Structure(
                Lattice.from_parameters(5.0, 5.0, 5.0, 90.0, 90.0, 90.0),
                ['Ga'],
                [[0.0, 0.0, 0.0]],
                model='configuration',
            )
            discuspdf.ase_io.to_ase(s)
            print("FAIL: no ImportError raised")
        except ImportError as e:
            if "discus-python[ase]" in str(e):
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
# (g) CIF ↔ ASE cross-check
# ---------------------------------------------------------------------------

def test_cif_ase_cross_check():
    """to_ase output matches the CIF read directly by ASE (cell + fractional
    positions within 1e-8)."""
    ase = pytest.importorskip("ase")
    gemmi = pytest.importorskip("gemmi")
    from discuspdf.cif_io import read_cif
    from discuspdf.ase_io import to_ase

    import ase.io
    s = read_cif(FIXTURES / "gaas_cod_9008845.cif").structure
    atoms = to_ase(s)
    ase_atoms = ase.io.read(str(FIXTURES / "gaas_cod_9008845.cif"))

    assert len(ase_atoms) == 8, f"expected 8 atoms, got {len(ase_atoms)}"
    assert len(atoms) == 8, f"expected 8 atoms, got {len(atoms)}"

    # Cell comparison.
    assert np.allclose(
        atoms.cell.array, ase_atoms.cell.array, rtol=1e-8
    ), (
        f"cell mismatch:\n{atoms.cell.array}\n{ase_atoms.cell.array}"
    )

    # Fractional positions: compare as sorted multisets.
    f_ours = atoms.get_scaled_positions()
    f_ase = ase_atoms.get_scaled_positions()
    # Sort both sets row-wise (by x, then y, then z) for a stable comparison.
    order_ours = np.lexsort((f_ours[:, 2], f_ours[:, 1], f_ours[:, 0]))
    order_ase = np.lexsort((f_ase[:, 2], f_ase[:, 1], f_ase[:, 0]))
    f_ours = f_ours[order_ours]
    f_ase = f_ase[order_ase]
    # Wrap to [0, 1) to normalise any coordinate at 1.0.
    f_ours = f_ours % 1.0
    f_ase = f_ase % 1.0
    assert np.allclose(f_ours, f_ase, atol=1e-8), (
        f"fractional positions mismatch:\n{f_ours}\n{f_ase}"
    )


# ---------------------------------------------------------------------------
# (h) Module-level importorskip
# ---------------------------------------------------------------------------

def test_module_importorskip():
    """The module imports cleanly when ase is present and the functions are
    callable."""
    ase = pytest.importorskip("ase")
    from discuspdf.ase_io import to_ase, from_ase

    assert callable(to_ase)
    assert callable(from_ase)
