"""Tests for SampledVolume, slices, isosurfaces and VTK export (Package B, Feature 20).

Independent references: the real DISCUS C3 Patterson peak at u = (0.3, 0.2, 0.1) (from the C3 macro
comment, by construction of mono2.cll), a sphere |r|² = c for isosurfaces, plane geometry from the
lattice, and — when the ``vtk`` package is installed — VTK's own XML readers.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pytest

from discuspdf import (Lattice, ReciprocalGrid, ScatteringVolume, delta_pdf_fractional,
                       delta_pdf_from_diffuse)
from discuspdf.viz import SampledVolume, isosurface, plot_isosurfaces, plot_slice, read_vtk, write_vtk

FIX = Path(__file__).parent / "fixtures" / "discus_6.23.01"
MONO = Lattice.from_parameters(5, 6, 7, 90, 105, 90)
ORTHO = Lattice.from_parameters(4, 5, 6, 90, 90, 90)


def oblique_volume():
    g = ReciprocalGrid(MONO, (-1, -1, -1), [[0.1, 0.02, 0], [0, 0.12, 0], [0.01, 0, 0.15]], (5, 4, 3))
    I = np.arange(60, dtype=float).reshape(5, 4, 3)
    m = np.zeros(I.shape, bool)
    m[1, 2, 0] = True
    sv = ScatteringVolume(g, np.where(m, np.nan, I), "e^2", mask=m, mask_reasons=("gap",),
                          variance=np.full(I.shape, 2.0))
    return sv, SampledVolume.from_scattering_volume(sv)


# --------------------------------------------------------------------------- geometry
def test_reciprocal_points_are_grid_q():
    sv, vol = oblique_volume()
    assert np.allclose(vol.points(), sv.grid.q())
    assert vol.unit == "Å⁻¹" and not vol.is_axis_aligned()
    aligned = SampledVolume.from_scattering_volume(
        ScatteringVolume(ReciprocalGrid.centered(ORTHO, (5, 4, 3), (0.25, 0.25, 0.5)), np.ones((5, 4, 3)), "a"))
    assert aligned.is_axis_aligned()


def test_real_space_of_c3_patterson_peak():
    pytest.importorskip("h5py")
    from discuspdf import read_discus_volume
    vol = read_discus_volume(FIX / "C3_total.h5")
    res = delta_pdf_fractional(vol, require_diffuse=False)
    sv = SampledVolume.from_fractional_fft(res, quantity="patterson")
    assert sv.unit == "Å" and not sv.is_axis_aligned()          # monoclinic: c is not along z
    v = sv.values.copy()
    centre = tuple(n // 2 for n in v.shape)
    v[tuple(slice(c - 2, c + 3) for c in centre)] = -np.inf    # exclude the origin peak
    # Bragg is present, so every lattice translation is a strong Patterson peak; the interatomic
    # vector is searched for inside one unit cell (|u_k| < 0.5).
    U = np.meshgrid(*res.u_axes, indexing="ij")
    v[np.any([np.abs(u) >= 0.5 for u in U], axis=0)] = -np.inf
    peak = np.unravel_index(np.argmax(v), v.shape)
    target = res.direct_basis.T @ np.array([0.3, 0.2, 0.1])
    other = res.direct_basis.T @ np.array([-0.3, -0.2, -0.1])
    pos = sv.points()[peak]
    half = 0.5 * np.linalg.norm(sv.steps, axis=1).sum()        # within one voxel diagonal
    assert min(np.linalg.norm(pos - target), np.linalg.norm(pos - other)) < half
    u = np.array([res.u_axes[k][peak[k]] for k in range(3)])
    assert np.allclose(pos, u @ res.direct_basis)


def test_cartesian_fft_result_is_axis_aligned():
    g = ReciprocalGrid.centered(ORTHO, (8, 9, 10), (0.25, 0.25, 0.25))
    res = delta_pdf_from_diffuse(np.ones(g.shape) * 0, g.q_axes(), raise_on_imaginary=False)
    vol = SampledVolume.from_fft_result(res)
    assert vol.is_axis_aligned()
    assert np.allclose(vol.points()[..., 0][:, 0, 0], res.r_axes[0])


@pytest.mark.parametrize("kwargs, match", [
    (dict(values=np.zeros((2, 2))), "3D"),
    (dict(space="hkl"), "space"),
    (dict(quantity="bad name"), "identifier"),
    (dict(steps=np.zeros((3, 3))), "dependent"),
    (dict(values=np.full((2, 2, 2), np.nan)), "finite"),
    (dict(mask=np.zeros((2, 2, 2), dtype=int)), "boolean"),
    (dict(extra={"v": np.zeros(3)}), "extra"),
])
def test_validation(kwargs, match):
    args = dict(values=np.zeros((2, 2, 2)), origin=np.zeros(3), steps=np.eye(3), space="real", quantity="x")
    args.update(kwargs)
    with pytest.raises(ValueError, match=match):
        SampledVolume(**args)


# --------------------------------------------------------------------------- VTK
def test_vts_round_trip_own_reader(tmp_path):
    sv, vol = oblique_volume()
    p = write_vtk(vol, tmp_path / "a.vts")
    back = read_vtk(p)
    assert back["kind"] == "StructuredGrid" and back["shape"] == (5, 4, 3)
    assert np.allclose(back["points"], sv.grid.q())
    assert np.array_equal(back["arrays"]["intensity"], vol.filled(np.nan), equal_nan=True)
    assert back["arrays"]["mask"][1, 2, 0] == 1 and back["arrays"]["mask"].sum() == 1
    assert np.all(back["arrays"]["variance"] == 2.0)
    side = json.loads((tmp_path / "a.vts.json").read_text())
    assert side["unit"] == "Å⁻¹" and side["metadata"]["bragg"] == "unknown"
    with pytest.raises(FileExistsError):
        write_vtk(vol, p)
    with pytest.raises(ValueError, match="vts"):
        write_vtk(vol, tmp_path / "a.vti")
    with pytest.raises(ValueError, match="end in"):
        write_vtk(vol, tmp_path / "a.vtk")


def test_vti_round_trip_own_reader(tmp_path):
    g = ReciprocalGrid.centered(ORTHO, (5, 4, 3), (0.25, 0.25, 0.5))
    I = np.random.default_rng(0).normal(size=g.shape)
    vol = SampledVolume.from_scattering_volume(ScatteringVolume(g, I, "a", bragg="bragg_subtracted"))
    back = read_vtk(write_vtk(vol, tmp_path / "b.vti"))
    assert back["kind"] == "ImageData"
    assert np.allclose(back["points"], g.q()) and np.array_equal(back["arrays"]["intensity"], I)


def test_vtk_package_reads_our_files(tmp_path):
    vtk = pytest.importorskip("vtk")
    from vtk.util.numpy_support import vtk_to_numpy
    sv, vol = oblique_volume()
    r = vtk.vtkXMLStructuredGridReader()
    r.SetFileName(str(write_vtk(vol, tmp_path / "a.vts")))
    r.Update()
    out = r.GetOutput()
    assert out.GetNumberOfPoints() == 60
    pts = vtk_to_numpy(out.GetPoints().GetData()).reshape(3, 4, 5, 3).transpose(2, 1, 0, 3)
    assert np.allclose(pts, sv.grid.q())
    vals = vtk_to_numpy(out.GetPointData().GetArray("intensity")).reshape(3, 4, 5).transpose(2, 1, 0)
    assert np.array_equal(vals, vol.filled(np.nan), equal_nan=True)
    g = ReciprocalGrid.centered(ORTHO, (5, 4, 3), (0.25, 0.25, 0.5))
    I = np.arange(60, dtype=float).reshape(5, 4, 3)
    v2 = SampledVolume.from_scattering_volume(ScatteringVolume(g, I, "a"))
    r2 = vtk.vtkXMLImageDataReader()
    r2.SetFileName(str(write_vtk(v2, tmp_path / "b.vti")))
    r2.Update()
    img = r2.GetOutput()
    pid = img.ComputePointId([3, 2, 1])
    assert np.allclose(img.GetPoint(pid), g.q()[3, 2, 1], atol=1e-12)
    assert vtk_to_numpy(img.GetPointData().GetArray("intensity"))[pid] == I[3, 2, 1]


# --------------------------------------------------------------------------- slices
def test_slice_keeps_true_geometry():
    sv, vol = oblique_volume()
    ax = plot_slice(vol, axis=2, index=0)
    mesh = ax.collections[0]
    coords = mesh.get_coordinates()                             # (n0+1, n1+1, 2) corner coordinates
    d0 = coords[1, 0] - coords[0, 0]
    d1 = coords[0, 1] - coords[0, 0]
    S = vol.steps
    assert np.linalg.norm(d0) == pytest.approx(np.linalg.norm(S[0]), rel=1e-12)
    assert np.linalg.norm(d1) == pytest.approx(np.linalg.norm(S[1]), rel=1e-12)
    cos_plot = d0 @ d1 / np.linalg.norm(d0) / np.linalg.norm(d1)
    cos_true = S[0] @ S[1] / np.linalg.norm(S[0]) / np.linalg.norm(S[1])
    assert cos_plot == pytest.approx(cos_true, abs=1e-12)
    arr = mesh.get_array()
    assert np.ma.is_masked(arr) and np.ma.getmaskarray(arr).reshape(5, 4)[1, 2]   # masked voxel stays masked
    assert "[1 0.2 0]*" in ax.get_xlabel()
    plt.close("all")


def test_signed_data_symmetric_colour_scale():
    g = ReciprocalGrid.centered(ORTHO, (5, 5, 5), (0.25, 0.25, 0.25))
    I = np.linspace(-3, 1, 125).reshape(5, 5, 5)
    vol = SampledVolume.from_scattering_volume(ScatteringVolume(g, I, "a", bragg="bragg_subtracted"))
    ax = plot_slice(vol, axis=0, index=0)
    norm = ax.collections[0].norm
    assert norm.vmin == -norm.vmax and norm.vmax == pytest.approx(np.abs(I[0]).max())
    with pytest.raises(ValueError, match="signed"):
        plot_slice(vol, log=True)
    with pytest.raises(ValueError, match="index"):
        plot_slice(vol, axis=0, index=5)
    plt.close("all")


# --------------------------------------------------------------------------- isosurfaces
@pytest.mark.parametrize("steps", [np.diag([0.1, 0.1, 0.1]),
                                   np.array([[0.1, 0.02, 0], [0, 0.09, 0.01], [0.015, 0, 0.11]])])
def test_isosurface_of_a_sphere(steps):
    pytest.importorskip("skimage")
    n = 41
    origin = -20 * steps.sum(0)
    idx = np.stack(np.meshgrid(*[np.arange(n)] * 3, indexing="ij"), -1)
    pts = origin + idx @ steps
    vol = SampledVolume(np.sum(pts ** 2, -1), origin, steps, "real", "r2")
    verts, faces = isosurface(vol, 1.0)
    r = np.linalg.norm(verts, axis=1)
    assert len(faces) > 100 and np.abs(r - 1.0).max() < 0.02
    with pytest.raises(ValueError, match="outside"):
        isosurface(vol, 1e9)
    ax = plot_isosurfaces(vol, [1.0, 2.0])
    assert ax.get_xlabel() == "x (Å)"
    plt.close("all")
