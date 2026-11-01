"""Structure pictures and viewer export (Package D, Feature 21)."""
from __future__ import annotations

import subprocess
import sys

import numpy as np
import pytest

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from discuspdf.lattice import Lattice  # noqa: E402
from discuspdf.neighbors import neighbor_list  # noqa: E402
from discuspdf.structure import Structure  # noqa: E402
from discuspdf.structure_ops import repeat  # noqa: E402
from discuspdf.structure_viz import (displacement_vectors, export_structure, plot_bonds,  # noqa: E402
                                     plot_displacements, plot_structure, species_colors)


@pytest.fixture(autouse=True)
def _close():
    yield
    plt.close("all")


def _mono(n=3):
    lat = Lattice.from_parameters(4.1, 5.3, 3.7, 90, 103, 90)
    cell = Structure(lat, ["Ga", "As"], [[0.05, 0.1, 0.2], [0.55, 0.62, 0.71]], model="configuration",
                     site_label=["s1", "s2"])
    return repeat(cell, (n, n, n)).structure


def _displaced(ref, sigma=0.05, seed=0):
    return ref.with_cartesian(ref.cartesian + np.random.default_rng(seed).normal(0, sigma, (ref.n_atoms, 3)))


def _offsets(coll):
    return np.column_stack([np.asarray(v, float) for v in coll._offsets3d])


def test_import_is_headless():
    code = "import sys, discuspdf, discuspdf.structure_viz; print('matplotlib' in sys.modules, 'ase' in sys.modules)"
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True).stdout.split()
    assert out == ["False", "False"]


def test_plot_structure_draws_stored_cartesian_positions():
    st = _mono()
    ax = plot_structure(st)
    got = np.vstack([_offsets(c) for c in ax.collections])
    want = st.cartesian
    # same set of points (drawn per species)
    order = lambda a: a[np.lexsort(a.T)]   # noqa: E731
    np.testing.assert_allclose(order(got), order(want), atol=1e-12)
    assert sorted(t.get_text() for t in ax.get_legend().get_texts()) == ["As", "Ga"]
    assert len(ax.lines) == 12                                     # supercell frame


def test_plot_structure_options():
    st = _mono(2)
    ax = plot_structure(st, color_by="site", select=st.atom_id[:5], cell="parent")
    assert sum(len(_offsets(c)) for c in ax.collections) == 5
    np.testing.assert_allclose(np.ptp(np.vstack([l.get_xydata() for l in ax.lines])[:, 0]) > 0, True)
    with pytest.raises(ValueError, match="colour"):
        plot_structure(st, colors={"Ga": "red"})
    with pytest.raises(ValueError, match="3D"):
        plot_structure(st, ax=plt.figure().add_subplot())
    with pytest.raises(ValueError, match="unknown atom_id"):
        plot_structure(st, select=[10_000])
    with pytest.raises(ValueError, match="color_by"):
        plot_structure(st, color_by="charge")


def test_species_colors_are_deterministic():
    assert species_colors(["Ga", "As", "Ga"]) == species_colors(["As", "Ga"])
    assert species_colors(["As", "Ga"])["As"] != species_colors(["As", "Ga"])["Ga"]


def test_displacement_vectors_minimum_image_across_boundary():
    ref = _mono(3)
    # move one atom by −0.1 Å along x so that it leaves the supercell and is wrapped to the far side
    cart = ref.cartesian.copy()
    k = int(np.argmin(ref.frac[:, 0]))
    cart[k, 0] -= 0.3
    conf = ref.with_cartesian(cart)
    wrapped = conf.wrapped()
    start, d = displacement_vectors(wrapped, ref)
    expect = np.zeros_like(cart)
    expect[k, 0] = -0.3
    order = np.argsort(wrapped.atom_id)
    np.testing.assert_allclose(d[order], expect[np.argsort(ref.atom_id)], atol=1e-9)
    np.testing.assert_allclose(start, ref.cartesian, atol=1e-12)


def test_displacement_vectors_refusals():
    ref = _mono(2)
    other = repeat(Structure(Lattice.from_parameters(4, 4, 4, 90, 90, 90), ["Ga"], [[0, 0, 0]],
                             model="configuration", site_label=["s1"]), (2, 2, 2)).structure
    with pytest.raises(ValueError, match="lattices"):
        displacement_vectors(other, ref)
    with pytest.raises(ValueError, match="two configuration atoms|no reference site"):
        displacement_vectors(_mono(2).remove(np.array([0])), _mono(2).remove(np.array([1])))


def test_plot_displacements_arrows_are_scaled_displacements():
    ref = _mono(2)
    conf = _displaced(ref)
    ax = plot_displacements(conf, ref, scale=10, draw_atoms=False)
    assert ax.get_title() == "displacements × 10"
    start, disp = displacement_vectors(conf, ref)
    for q in ax.collections:                                  # one quiver per species
        sp = q.get_label()
        k = np.asarray(conf.species) == sp
        segs = np.asarray(q._segments3d)
        n = int(k.sum())
        assert segs.shape == (3 * n, 2, 3)                    # shaft + two head lines per arrow
        # matplotlib stores each shaft as (tip, tail)
        np.testing.assert_allclose(segs[:n, 1], start[k], atol=1e-12)
        np.testing.assert_allclose(segs[:n, 0] - segs[:n, 1], 10 * disp[k], atol=1e-12)
    with pytest.raises(ValueError, match="scale"):
        plot_displacements(conf, ref, scale=0)


def test_plot_bonds_segments_match_pair_vectors():
    ref = _mono(2)
    nl = neighbor_list(ref, 4.0).select(shell=1)
    ax = plot_bonds(nl, draw_atoms=False)
    segs = [l.get_data_3d() for l in ax.lines if len(l.get_xdata()) == 2]
    assert len(segs) == nl.n_pairs
    lengths = sorted(float(np.linalg.norm(np.array(s)[:, 1] - np.array(s)[:, 0])) for s in segs)
    np.testing.assert_allclose(lengths, np.sort(nl.distance), atol=1e-12)
    with pytest.raises(TypeError):
        plot_bonds(ref)


@pytest.mark.parametrize("fmt,name", [("extxyz", "a.xyz"), ("vasp", "POSCAR"), ("xsf", "a.xsf"), ("cif", "a.cif")])
def test_export_keeps_cell_and_positions(tmp_path, fmt, name):
    ase_io = pytest.importorskip("ase.io")
    conf = _displaced(_mono(2))
    p = export_structure(conf, tmp_path / name, format=fmt)
    back = ase_io.read(p, format=fmt)
    np.testing.assert_allclose(back.cell.array, conf.lattice.matrix, atol=1e-5)
    np.testing.assert_allclose(back.get_positions(), conf.cartesian, atol=1e-5)
    assert list(back.get_chemical_symbols()) == list(conf.species)


def test_export_refuses_a_format_that_reorients(tmp_path):
    pytest.importorskip("ase.io")
    th = np.deg2rad(30)
    R = np.array([[np.cos(th), -np.sin(th), 0], [np.sin(th), np.cos(th), 0], [0, 0, 1]])
    A = Lattice.from_parameters(4.0, 4.5, 5.0, 90, 95, 90).matrix @ R.T
    st = Structure(Lattice.from_matrix(A), ["Ga", "As"], [[0.1, 0.2, 0.3], [0.6, 0.7, 0.8]],
                   model="configuration")
    export_structure(st, tmp_path / "ok.xyz", format="extxyz")      # extxyz keeps the full cell matrix
    with pytest.raises(ValueError, match="cell changed|positions changed"):
        export_structure(st, tmp_path / "bad.cif", format="cif")     # CIF stores parameters only
