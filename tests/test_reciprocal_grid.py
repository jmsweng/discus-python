"""Tests for ReciprocalGrid (Package B, Feature 6).

Independent references used here (not derived from the code under test):
* monoclinic d-spacing textbook formula and the reciprocal angle beta* = 180° - beta;
* the real DISCUS 6.23.01 C3/C4 fixtures (C4 is deliberately skewed) and the corner
  commands in their committed macros, parsed from the macro text;
* a brute-force direct sum for the real-space dual grid of a skewed volume.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pytest

from discuspdf import Lattice, ReciprocalGrid, ReciprocalVolume
from discuspdf.reciprocal_grid import FFTCompatibility

FIX = Path(__file__).parent / "fixtures" / "discus_6.23.01"

MONO = dict(a=5.0, b=6.0, c=7.0, alpha=90.0, beta=105.0, gamma=90.0)


def mono():
    return Lattice.from_parameters(**MONO)


def ortho():
    return Lattice.from_parameters(4.0, 5.0, 6.0, 90, 90, 90)


def monoclinic_d(h, k, l, a, b, c, beta):
    """Textbook monoclinic (unique axis b) d-spacing; independent of Lattice."""
    sb, cb = np.sin(np.deg2rad(beta)), np.cos(np.deg2rad(beta))
    inv_d2 = (h * h / a**2 + k * k * sb**2 / b**2 + l * l / c**2 - 2 * h * l * cb / (a * c)) / sb**2
    return 1.0 / np.sqrt(inv_d2)


def macro_corners(path: Path):
    """Parse ll/lr/ul/tl/na/no/nt from a committed DISCUS macro."""
    out = {}
    for line in path.read_text().splitlines():
        m = re.match(r"\s*(ll|lr|ul|tl)\s+(.*)", line)
        if m:
            out[m.group(1)] = [float(x) for x in m.group(2).split(",")]
        m = re.match(r"\s*(na|no|nt)\s+(\d+)", line)
        if m:
            out[m.group(1)] = int(m.group(2))
    return out


# --------------------------------------------------------------------------- Q calibration
@pytest.mark.parametrize("hkl", [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 0, -1), (2, -1, 3), (-3, 2, 1)])
def test_reflection_q_matches_textbook_d_spacing(hkl):
    lat = mono()
    grid = ReciprocalGrid.from_ranges(lat, (-4, 4, 33), (-5, 5, 41), (-6, 6, 49))
    idx = grid.hkl_to_index(hkl)
    assert np.allclose(idx, np.round(idx))                       # reflection is a grid point
    i = tuple(int(round(x)) for x in idx)
    qmag = grid.q_magnitude()[i]
    d_ref = monoclinic_d(*hkl, MONO["a"], MONO["b"], MONO["c"], MONO["beta"])
    assert qmag == pytest.approx(2 * np.pi / d_ref, rel=1e-12)
    assert np.allclose(grid.hkl()[i], hkl, atol=1e-12)


def test_reciprocal_angle_and_real_space_distance():
    grid = ReciprocalGrid.centered(mono(), (9, 9, 9), (1.0, 1.0, 1.0))
    S = grid.q_step_vectors()
    cos = S[0] @ S[2] / np.linalg.norm(S[0]) / np.linalg.norm(S[2])
    assert np.rad2deg(np.arccos(cos)) == pytest.approx(180.0 - MONO["beta"], abs=1e-10)
    assert abs(S[0] @ S[1]) < 1e-12 and abs(S[1] @ S[2]) < 1e-12
    # A fractional real-space vector u = (0.3, 0.2, 0.1): law-of-cosines length (independent).
    u = np.array([0.3, 0.2, 0.1])
    x, y, z = u * [MONO["a"], MONO["b"], MONO["c"]]
    ref = np.sqrt(x * x + y * y + z * z + 2 * x * z * np.cos(np.deg2rad(MONO["beta"])))
    assert np.linalg.norm(grid.lattice.frac_to_cart(u)) == pytest.approx(ref, rel=1e-12)


def test_orthogonal_cell_q_axes_are_2pi_h_over_a():
    lat = ortho()
    grid = ReciprocalGrid.centered(lat, (5, 6, 7), (0.5, 0.25, 0.2))
    qx, qy, qz = grid.q_axes()
    h, k, l = grid.hkl_axes()
    assert np.allclose(qx, 2 * np.pi * h / 4.0)
    assert np.allclose(qy, 2 * np.pi * k / 5.0)
    assert np.allclose(qz, 2 * np.pi * l / 6.0)
    comp = grid.fft_compatibility()
    assert comp.fractional and comp.cartesian and comp.reasons == ()
    vol = ReciprocalVolume(np.zeros(grid.shape), grid.q_axes(), "cartesian_q", "1/angstrom")
    assert all(np.array_equal(a, b) for a, b in zip(vol.q_axes, grid.q_axes()))


def test_monoclinic_grid_is_fractional_only():
    grid = ReciprocalGrid.centered(mono(), (5, 5, 5), (0.5, 0.5, 0.5))
    comp = grid.fft_compatibility()
    assert comp.fractional and not comp.cartesian
    assert any("cartesian" in r for r in comp.reasons)
    assert not grid.is_q_orthogonal()
    with pytest.raises(ValueError, match="Cartesian"):
        grid.q_axes()


def test_hexagonal_hk_plane_is_not_q_orthogonal():
    hexa = Lattice.from_parameters(3.0, 3.0, 5.0, 90, 90, 120)
    grid = ReciprocalGrid.from_ranges(hexa, (-2, 2, 9), (-2, 2, 9), (0, 0, 1))
    assert grid.ndim == 2 and grid.is_hkl_axis_aligned()
    S = grid.q_step_vectors()
    cos = S[0] @ S[1] / np.linalg.norm(S[0]) / np.linalg.norm(S[1])
    assert np.rad2deg(np.arccos(cos)) == pytest.approx(60.0, abs=1e-10)   # gamma* = 180 - 120
    assert not grid.is_q_orthogonal()
    assert ReciprocalGrid.from_ranges(ortho(), (-1, 1, 3), (-1, 1, 3), (0, 0, 1)).is_q_orthogonal()


# --------------------------------------------------------------------------- real DISCUS fixtures
@pytest.mark.parametrize("case", ["C3", "C4"])
def test_from_corners_matches_real_discus_output(case):
    pytest.importorskip("h5py")
    from discuspdf import read_discus_volume
    vol = read_discus_volume(FIX / f"{case}_total.h5")
    from_file = ReciprocalGrid.from_discus_volume(vol)
    m = macro_corners(FIX / f"{case}.mac")
    from_macro = ReciprocalGrid.from_corners(from_file.lattice, m["ll"], m["lr"], m["ul"], m["tl"],
                                             (m["na"], m["no"], m["nt"]))
    assert from_macro.isclose(from_file, rtol=0, atol=1e-6)
    # Element-wise: the engine's increments are columns, ours rows.
    for k in range(3):
        assert np.allclose(from_file.axis_vectors[k], vol.increments[:, k])
    dc = from_macro.discus_corners()
    for key in ("ll", "lr", "ul", "tl"):
        assert np.allclose(dc[key], m[key], atol=1e-12)
    assert (dc["na"], dc["no"], dc["nt"]) == (m["na"], m["no"], m["nt"])


def test_skewed_c4_grid_geometry():
    pytest.importorskip("h5py")
    from discuspdf import read_discus_volume
    c4 = ReciprocalGrid.from_discus_volume(read_discus_volume(FIX / "C4_total.h5"))
    assert not c4.is_hkl_axis_aligned()
    with pytest.raises(ValueError, match="skewed"):
        c4.hkl_axes()
    comp = c4.fft_compatibility()
    assert not comp.fractional and not comp.cartesian
    # Last point of the ordinate axis is the macro's ul corner (-3, 5, -6).
    assert np.allclose(c4.index_to_hkl([0, 40, 0]), (-3, 5, -6))
    assert c4.handedness() == 1


def test_c3_real_space_steps_match_fractional_transform():
    pytest.importorskip("h5py")
    from discuspdf import read_discus_volume, delta_pdf_fractional
    vol = read_discus_volume(FIX / "C3_total.h5")
    grid = ReciprocalGrid.from_discus_volume(vol)
    assert grid.fft_compatibility().fractional
    R = grid.real_space_step_vectors()
    res = delta_pdf_fractional(vol, require_diffuse=False)
    for k in range(3):
        assert R[k, k] == pytest.approx(np.diff(res.u_axes[k])[0], rel=1e-12)
    assert np.allclose(R - np.diag(np.diag(R)), 0)
    assert grid.lattice.isclose(Lattice.from_matrix(res.direct_basis))


# --------------------------------------------------------------------------- dual grid, brute force
def test_real_space_dual_of_skewed_grid_by_direct_sum():
    rng = np.random.default_rng(7)
    V = np.array([[0.11, 0.02, 0.0], [-0.03, 0.09, 0.01], [0.0, 0.015, 0.13]])
    grid = ReciprocalGrid(mono(), (-0.4, -0.3, -0.5), V, (6, 5, 7))
    I = rng.normal(size=grid.shape)
    F = np.fft.fftn(I)
    R = grid.real_space_step_vectors()
    hkl = grid.hkl().reshape(-1, 3)
    for m in [(0, 0, 0), (1, 0, 0), (2, 3, 1), (5, 4, 6)]:
        u = np.asarray(m, float) @ R
        direct = np.sum(I.ravel() * np.exp(-2j * np.pi * hkl @ u))
        assert direct == pytest.approx(np.exp(-2j * np.pi * grid.corner @ u) * F[m], rel=1e-10, abs=1e-10)
    assert np.allclose(grid.axis_vectors @ R.T, np.diag(1.0 / np.array(grid.shape)))


# --------------------------------------------------------------------------- index <-> hkl
def test_index_hkl_round_trip_and_edges():
    V = np.array([[0.1, 0.0, 0.02], [0.0, 0.2, 0.0], [0.03, 0.0, 0.25]])
    grid = ReciprocalGrid(mono(), (1.0, -2.0, 0.5), V, (4, 3, 5))
    idx = np.array([[0, 0, 0], [3, 2, 4], [1.5, 0.25, 2.0]])
    assert np.allclose(grid.hkl_to_index(grid.index_to_hkl(idx)), idx)
    assert np.allclose(grid.index_to_hkl([-0.5, 0, 0]), grid.corner - 0.5 * V[0])   # voxel edge
    full = grid.hkl()
    assert full.shape == (4, 3, 5, 3)
    assert np.allclose(full[3, 2, 4], grid.corner + 3 * V[0] + 2 * V[1] + 4 * V[2])
    assert np.allclose(grid.q(), grid.lattice.hkl_to_q(full))


def test_even_and_odd_centering():
    for shape in [(4, 5, 6), (7, 7, 7)]:
        grid = ReciprocalGrid.centered(ortho(), shape, (0.1, 0.2, 0.3))
        assert grid.is_centered()
        assert np.allclose(grid.index_to_hkl([n // 2 for n in shape]), 0)
        h = grid.hkl_axes()[0]
        assert np.allclose(h, (np.arange(shape[0]) - shape[0] // 2) * 0.1)
    shifted = ReciprocalGrid.from_ranges(ortho(), (0, 1, 3), (-1, 1, 3), (-1, 1, 3))
    assert not shifted.is_centered()
    comp = shifted.fft_compatibility()
    assert not comp.fractional and not comp.cartesian and any("n//2" in r for r in comp.reasons)


def test_decreasing_axis_is_left_handed_and_not_fft_ready():
    grid = ReciprocalGrid(ortho(), (1, -1, -1), np.diag([-1.0, 1.0, 1.0]), (3, 3, 3))
    assert grid.handedness() == -1
    comp = grid.fft_compatibility()
    assert not comp.fractional and any("positive" in r for r in comp.reasons)
    assert not comp.cartesian


# --------------------------------------------------------------------------- planes and lines
def test_plane_and_line():
    grid = ReciprocalGrid.centered(mono(), (9, 11, 13), (0.25, 0.2, 0.5))
    p = grid.plane(2, 6)
    assert p.ndim == 2 and p.shape == (9, 11, 1) and p.active_axes == (0, 1)
    assert np.allclose(p.corner, grid.index_to_hkl([0, 0, 6]))
    assert np.allclose(p.hkl()[..., 0, :], grid.hkl()[:, :, 6, :])
    with pytest.raises(ValueError, match="3D grid"):
        p.hkl_to_index((0, 0, 0))
    with pytest.raises(ValueError, match="3D"):
        p.real_space_step_vectors()
    with pytest.raises(ValueError):
        p.handedness()
    assert p.fft_compatibility().reasons[0].startswith("grid has 2 active axes")
    line = ReciprocalGrid.from_ranges(mono(), (-2, 2, 5), (1, 1, 1), (0, 0, 1))
    assert line.ndim == 1 and np.allclose(line.hkl()[:, 0, 0, 1], 1.0)


def test_inactive_axis_vector_is_ignored():
    a = ReciprocalGrid(ortho(), (0, 0, 0), [[0.1, 0, 0], [0, 0.1, 0], [0, 0, 0]], (3, 3, 1))
    b = ReciprocalGrid(ortho(), (0, 0, 0), [[0.1, 0, 0], [0, 0.1, 0], [5.0, 7.0, 9.0]], (3, 3, 1))
    assert a == b and hash(a) == hash(b)
    assert np.array_equal(a.hkl(), b.hkl())
    assert np.allclose(b.q_step_vectors()[2], 0)


def test_oblique_plane_from_corners():
    # (h, h, l) plane: abscissa along [1 1 0], ordinate along [0 0 1].
    grid = ReciprocalGrid.from_corners(ortho(), (-2, -2, -3), (2, 2, -3), (-2, -2, 3), (-2, -2, -3), (9, 13, 1))
    assert grid.ndim == 2 and not grid.is_hkl_axis_aligned()
    assert np.allclose(grid.axis_vectors[0], (0.5, 0.5, 0))
    assert grid.is_q_orthogonal()          # [110]* ⟂ [001]* in an orthorhombic cell
    with pytest.raises(ValueError, match="single point"):
        ReciprocalGrid.from_corners(ortho(), (0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (2, 2, 1))


# --------------------------------------------------------------------------- validation
@pytest.mark.parametrize("shape", [(0, 2, 2), (2, 2), (2.0, 2, 2), (True, 2, 2), "abc"])
def test_bad_shape(shape):
    with pytest.raises(ValueError):
        ReciprocalGrid(ortho(), (0, 0, 0), np.eye(3), shape)


def test_bad_vectors():
    with pytest.raises(ValueError, match="zero step"):
        ReciprocalGrid(ortho(), (0, 0, 0), [[0, 0, 0], [0, 1, 0], [0, 0, 1]], (2, 2, 2))
    with pytest.raises(ValueError, match="dependent"):
        ReciprocalGrid(ortho(), (0, 0, 0), [[1, 1, 0], [2, 2, 0], [0, 0, 1]], (2, 2, 2))
    with pytest.raises(ValueError, match="dependent"):      # nearly dependent in a plane
        ReciprocalGrid(ortho(), (0, 0, 0), [[1, 0, 0], [1, 1e-12, 0], [0, 0, 0]], (2, 2, 1))
    with pytest.raises(ValueError, match="finite"):
        ReciprocalGrid(ortho(), (0, np.nan, 0), np.eye(3), (2, 2, 2))
    with pytest.raises(ValueError, match="3x3"):
        ReciprocalGrid(ortho(), (0, 0, 0), np.eye(2), (2, 2, 2))
    with pytest.raises(TypeError):
        ReciprocalGrid(None, (0, 0, 0), np.eye(3), (2, 2, 2))
    with pytest.raises(ValueError):
        ReciprocalGrid.centered(ortho(), (3, 3, 3), (0.1, 0, 0.1))
    with pytest.raises(ValueError, match="start == stop"):
        ReciprocalGrid.from_ranges(ortho(), (0, 1, 1), (0, 1, 2), (0, 1, 2))


def test_immutable():
    grid = ReciprocalGrid.centered(ortho(), (3, 3, 3), (0.1, 0.1, 0.1))
    with pytest.raises(AttributeError):
        grid.shape = (1, 1, 1)
    c = grid.corner
    c[0] = 99
    assert grid.corner[0] != 99
    with pytest.raises(ValueError):
        grid._corner[0] = 1.0


# --------------------------------------------------------------------------- persistence
def test_dict_and_json_round_trip_preserves_orientation():
    rot = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]], float)          # 90° about z
    lat = Lattice.from_matrix(mono().matrix @ rot.T)
    grid = ReciprocalGrid(lat, (0.5, -1, 2), [[0.1, 0.01, 0], [0, 0.2, 0], [0, 0, 0.3]], (3, 4, 5))
    again = ReciprocalGrid.from_json(grid.to_json())
    assert again == grid
    assert np.allclose(again.q(), grid.q())
    d = json.loads(grid.to_json())
    assert d["schema"] == "discuspdf.ReciprocalGrid" and d["version"] == 1


def test_from_dict_refuses_unknown_schema():
    d = ReciprocalGrid.centered(ortho(), (3, 3, 3), (0.1, 0.1, 0.1)).to_dict()
    with pytest.raises(ValueError, match="version"):
        ReciprocalGrid.from_dict({**d, "version": 2})
    with pytest.raises(ValueError, match="schema"):
        ReciprocalGrid.from_dict({**d, "schema": "other"})
    bad = dict(d)
    del bad["corner"]
    with pytest.raises(ValueError, match="missing"):
        ReciprocalGrid.from_dict(bad)


def test_fft_compatibility_is_not_a_truth_value():
    comp = ReciprocalGrid.centered(ortho(), (3, 3, 3), (0.1, 0.1, 0.1)).fft_compatibility()
    assert isinstance(comp, FFTCompatibility)
    with pytest.raises(TypeError):
        bool(comp)
