"""Tests for the DiffPy adapter (Package A, Feature A2).

All tests that require ``diffpy.structure`` use
``pytest.importorskip("diffpy.structure")`` so the suite is green without
the ``diffpy`` extra.  ``gemmi`` is additionally required to build the
fixture structures via ``read_cif``.
"""
from __future__ import annotations

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
    u_to_b,
)

FIXTURES = Path(__file__).parent / "fixtures" / "interop"


# ---------------------------------------------------------------------------
# (a) Round trip
# ---------------------------------------------------------------------------

def test_roundtrip_gaas():
    """GaAs CIF → to_diffpy → from_diffpy round trip preserves the structure.
    b_iso is None (U=0 in the CIF) and must stay None."""
    diffpy = pytest.importorskip("diffpy.structure")
    gemmi = pytest.importorskip("gemmi")
    from discuspdf.cif_io import read_cif
    from discuspdf.diffpy_io import to_diffpy, from_diffpy

    s = read_cif(FIXTURES / "gaas_cod_9008845.cif").structure
    assert s.b_iso is None

    ds = to_diffpy(s, allow_loss=True)
    s2 = from_diffpy(ds, model=s.model)

    assert s2.isclose(s, rtol=1e-8, atol=1e-12)
    assert s2.n_atoms == s.n_atoms
    assert s2.model == s.model
    # All Uisoequiv are 0 → b_iso must be None after the round trip.
    assert s2.b_iso is None


def test_roundtrip_mixed():
    """Mixed-occupancy CIF → to_diffpy → from_diffpy round trip preserves the
    structure, including per-atom B within 1e-6."""
    diffpy = pytest.importorskip("diffpy.structure")
    gemmi = pytest.importorskip("gemmi")
    from discuspdf.cif_io import read_cif
    from discuspdf.diffpy_io import to_diffpy, from_diffpy

    s = read_cif(FIXTURES / "ga07in03as_mixed.cif").structure
    assert s.b_iso is not None

    ds = to_diffpy(s, allow_loss=True)
    s2 = from_diffpy(ds, model=s.model)

    assert s2.isclose(s, rtol=1e-8, atol=1e-12)
    assert s2.n_atoms == s.n_atoms
    assert s2.model == s.model
    assert s2.b_iso is not None
    assert np.allclose(s2.b_iso, s.b_iso, atol=1e-6)


def test_roundtrip_triclinic():
    """Triclinic P-1 CIF → to_diffpy → from_diffpy round trip preserves the
    structure (lattice parameters and fractional coordinates)."""
    diffpy = pytest.importorskip("diffpy.structure")
    gemmi = pytest.importorskip("gemmi")
    from discuspdf.cif_io import read_cif
    from discuspdf.diffpy_io import to_diffpy, from_diffpy

    s = read_cif(FIXTURES / "triclinic_p1bar.cif").structure
    ds = to_diffpy(s, allow_loss=True)
    s2 = from_diffpy(ds, model=s.model)

    assert s2.isclose(s, rtol=1e-8, atol=1e-12)
    assert s2.n_atoms == s.n_atoms
    assert s2.model == s.model


def test_roundtrip_b_iso_none_and_per_atom():
    """b_iso=None (GaAs) stays None; per-atom B (mixed) survives within 1e-6."""
    diffpy = pytest.importorskip("diffpy.structure")
    gemmi = pytest.importorskip("gemmi")
    from discuspdf.cif_io import read_cif
    from discuspdf.diffpy_io import to_diffpy, from_diffpy

    # b_iso=None (the real CIF has U=0 → b_iso=None).
    s_none = read_cif(FIXTURES / "gaas_cod_9008845.cif").structure
    assert s_none.b_iso is None
    ds = to_diffpy(s_none, allow_loss=True)
    s2 = from_diffpy(ds, model=s_none.model)
    assert s2.b_iso is None

    # Per-atom B (the mixed fixture has B values from U).
    s_b = read_cif(FIXTURES / "ga07in03as_mixed.cif").structure
    assert s_b.b_iso is not None
    ds = to_diffpy(s_b, allow_loss=True)
    s2 = from_diffpy(ds, model=s_b.model)
    assert s2.b_iso is not None
    assert np.allclose(s2.b_iso, s_b.b_iso, atol=1e-6)


# ---------------------------------------------------------------------------
# (b) Distance invariance
# ---------------------------------------------------------------------------

def test_diffpy_distance_invariance():
    """The multiset of pairwise minimum-image distances is preserved by the
    to_diffpy → from_diffpy round trip.  Distances are computed independently
    from the fractional coordinates and the lattice metric (minimum image),
    not via any DiffPy distance API."""
    diffpy = pytest.importorskip("diffpy.structure")
    gemmi = pytest.importorskip("gemmi")
    from discuspdf.cif_io import read_cif
    from discuspdf.diffpy_io import to_diffpy, from_diffpy

    s = read_cif(FIXTURES / "gaas_cod_9008845.cif").structure
    ds = to_diffpy(s, allow_loss=True)
    s2 = from_diffpy(ds, model=s.model)

    # Independent (N, N) distance computation from the round-tripped
    # structure: minimum-image fractional difference mapped through the
    # lattice metric.
    A = s2.lattice.matrix
    frac = s2.frac
    n = s2.n_atoms
    out = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            df_min, _ = s2.lattice.minimum_image(frac[i] - frac[j])
            out[i, j] = float(np.linalg.norm(df_min @ A))

    # (N, N) distance matrix from the original Structure.
    d_ours = s.distance_matrix()

    # Both should have the same multiset of pairwise distances.
    flat_ours = np.sort(d_ours.flatten())
    flat_rt = np.sort(out.flatten())
    assert np.allclose(flat_ours, flat_rt, atol=1e-6), (
        f"distance multisets differ after round trip: max diff = "
        f"{np.abs(flat_ours - flat_rt).max()}"
    )


# ---------------------------------------------------------------------------
# (c) Foreign object
# ---------------------------------------------------------------------------

def test_foreign_diffpy():
    """A hand-built diffpy Structure (no to_diffpy involvement) converts with
    model given, occupancy 1.0, and b_iso=None (all Uisoequiv 0)."""
    diffpy = pytest.importorskip("diffpy.structure")
    gemmi = pytest.importorskip("gemmi")
    from discuspdf.diffpy_io import from_diffpy

    lattice = diffpy.Lattice(5.6537, 5.6537, 5.6537, 90.0, 90.0, 90.0)
    ga = diffpy.Atom(atype="Ga", xyz=(0.0, 0.0, 0.0), label="Ga1")
    as_ = diffpy.Atom(atype="As", xyz=(0.25, 0.25, 0.25), label="As1")
    ds = diffpy.Structure(lattice=lattice, atoms=[ga, as_])

    s = from_diffpy(ds, model="configuration")

    assert s.n_atoms == 2
    assert s.model == "configuration"
    assert np.allclose(s.occupancy, 1.0)
    assert s.b_iso is None


def test_foreign_diffpy_missing_model_error():
    """from_diffpy without the required model keyword → TypeError."""
    diffpy = pytest.importorskip("diffpy.structure")
    gemmi = pytest.importorskip("gemmi")
    from discuspdf.diffpy_io import from_diffpy

    lattice = diffpy.Lattice(5.6537, 5.6537, 5.6537, 90.0, 90.0, 90.0)
    ga = diffpy.Atom(atype="Ga", xyz=(0.0, 0.0, 0.0), label="Ga1")
    ds = diffpy.Structure(lattice=lattice, atoms=[ga])

    with pytest.raises(TypeError):
        from_diffpy(ds)


# ---------------------------------------------------------------------------
# (d) Losses raise and report with allow_loss=True
# ---------------------------------------------------------------------------

def test_diffpy_aniso_loss():
    """An anisotropic atom (distinct U11/U22/U33) → InteropError by default;
    with allow_loss=True → InteropWarning and the isotropic equivalent is kept
    (b_iso = 8π²·Uisoequiv)."""
    diffpy = pytest.importorskip("diffpy.structure")
    gemmi = pytest.importorskip("gemmi")
    from discuspdf.diffpy_io import from_diffpy

    lattice = diffpy.Lattice(5.0, 5.0, 5.0, 90.0, 90.0, 90.0)
    atom = diffpy.Atom(atype="Fe", xyz=(0.0, 0.0, 0.0), label="Fe1")
    # Set anisotropy=True FIRST, then the U values, so that Uisoequiv is the
    # correct lattice-weighted isotropic equivalent (P14).
    atom.anisotropy = True
    atom.U11 = 0.013
    atom.U22 = 0.012
    atom.U33 = 0.011
    assert atom.anisotropy, "expected the atom to be anisotropic"

    ds = diffpy.Structure(lattice=lattice, atoms=[atom])

    # Default: raise.
    with pytest.raises(InteropError, match="anisotropic"):
        from_diffpy(ds, model="configuration")

    # allow_loss=True: warn and keep the isotropic equivalent.
    with pytest.warns(InteropWarning):
        s = from_diffpy(ds, model="configuration", allow_loss=True)
    assert s.n_atoms == 1
    expected_b = u_to_b(0.012)
    assert abs(s.b_iso[0] - expected_b) < 1e-6, (
        f"b_iso={s.b_iso[0]}, expected {expected_b}"
    )


# ---------------------------------------------------------------------------
# (e) Errors
# ---------------------------------------------------------------------------

def test_diffpy_missing_model_error():
    """from_diffpy without the required model keyword → TypeError."""
    diffpy = pytest.importorskip("diffpy.structure")
    gemmi = pytest.importorskip("gemmi")
    from discuspdf.diffpy_io import from_diffpy

    lattice = diffpy.Lattice(5.6537, 5.6537, 5.6537, 90.0, 90.0, 90.0)
    ga = diffpy.Atom(atype="Ga", xyz=(0.0, 0.0, 0.0), label="Ga1")
    ds = diffpy.Structure(lattice=lattice, atoms=[ga])

    with pytest.raises(TypeError):
        from_diffpy(ds)


def test_diffpy_unknown_element():
    """A Structure with an unknown element symbol → InteropError from to_diffpy
    (element_from_type_symbol)."""
    diffpy = pytest.importorskip("diffpy.structure")
    gemmi = pytest.importorskip("gemmi")
    from discuspdf.diffpy_io import to_diffpy

    lattice = Lattice.from_parameters(5.0, 5.0, 5.0, 90.0, 90.0, 90.0)
    s = Structure(
        lattice,
        np.array(["Xx"], dtype=object),
        np.array([[0.0, 0.0, 0.0]]),
        model="configuration",
    )
    with pytest.raises(InteropError):
        to_diffpy(s, allow_loss=True)


# ---------------------------------------------------------------------------
# (f) Optional import
# ---------------------------------------------------------------------------

def test_optional_import_diffpy():
    """A subprocess test sets sys.modules['diffpy'] = None and
    sys.modules['diffpy.structure'] = None and checks that importing the
    modules works and to_diffpy raises ImportError."""
    code = textwrap.dedent("""\
        import sys
        sys.modules['diffpy'] = None
        sys.modules['diffpy.structure'] = None
        import discuspdf
        import discuspdf.diffpy_io
        from discuspdf.lattice import Lattice
        from discuspdf.structure import Structure
        try:
            s = Structure(
                Lattice.from_parameters(5.0, 5.0, 5.0, 90.0, 90.0, 90.0),
                ['Ga'],
                [[0.0, 0.0, 0.0]],
                model='configuration',
            )
            discuspdf.diffpy_io.to_diffpy(s)
            print("FAIL: no ImportError raised")
        except ImportError as e:
            if "discus-python[diffpy]" in str(e):
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
# (g) CIF ↔ DiffPy cross-check
# ---------------------------------------------------------------------------

def test_cif_diffpy_cross_check():
    """to_diffpy preserves the lattice parameters (orientation is NOT
    preserved, but the parameters are)."""
    diffpy = pytest.importorskip("diffpy.structure")
    gemmi = pytest.importorskip("gemmi")
    from discuspdf.cif_io import read_cif
    from discuspdf.diffpy_io import to_diffpy

    s = read_cif(FIXTURES / "gaas_cod_9008845.cif").structure
    ds = to_diffpy(s, allow_loss=True)

    assert np.allclose(
        ds.lattice.cell_parms(), s.lattice.parameters, rtol=1e-8
    ), (
        f"lattice parameters mismatch:\n"
        f"{ds.lattice.cell_parms()}\n{s.lattice.parameters}"
    )


# ---------------------------------------------------------------------------
# (h) Module-level importorskip
# ---------------------------------------------------------------------------

def test_module_importorskip():
    """The module imports cleanly when diffpy is present and the functions are
    callable."""
    diffpy = pytest.importorskip("diffpy.structure")
    from discuspdf.diffpy_io import to_diffpy, from_diffpy

    assert callable(to_diffpy)
    assert callable(from_diffpy)
