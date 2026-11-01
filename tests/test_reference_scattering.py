"""Analytical benchmarks for the independent structure-factor reference (Package D, Feature 10)."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from discuspdf.lattice import Lattice
from discuspdf.structure import Structure
from discuspdf.structure_ops import repeat
from discuspdf.neighbors import neighbor_list
from discuspdf.correlations import occupational_correlation
from discuspdf.reference_scattering import (
    ScatteringFactors, average_structure_factor, diffuse_intensity, intensity, scattering_vectors,
    structure_factor,
)

ROOT = Path(__file__).resolve().parents[1]
# DABAX Waasmaier-Kirfel rows already transcribed and cross-checked against DISCUS's own table in PR #8.
WK = json.loads((ROOT / "benchmarks/discus_6_23_01/gate2_prereg/wk_coefficients.json").read_text())

NEUTRON = ScatteringFactors.neutron({"Cu": 7.718, "Au": 7.63, "A": 5.0, "B": -3.0}, source="test values")


def _wk(el):
    return WK["elements"][el.upper()]["dabax"]


def _wk_value(el, s):
    c = _wk(el)
    s2 = np.asarray(s, float) ** 2
    return c["c"] + sum(a * np.exp(-b * s2) for a, b in zip(c["a"], c["b"]))


XRAY = ScatteringFactors.xray_wk({"Cu": _wk("CU"), "Au": _wk("AU")}, source="DABAX (PR #8 transcription)")


def _cubic(a=4.0):
    return Lattice.from_parameters(a, a, a, 90, 90, 90)


def _single(species="Cu", frac=(0.1, 0.2, 0.3), b=None, a=4.0):
    return Structure(_cubic(a), [species], [frac], model="configuration",
                     b_iso=None if b is None else [b])


def _random_hkl(n=50, seed=0, scale=3.0):
    return np.random.default_rng(seed).uniform(-scale, scale, (n, 3))


# ------------------------------------------------------------------------------------- conventions
def test_q_convention_and_s():
    st = repeat(_single(frac=(0, 0, 0)), (2, 3, 4)).structure
    Q = scattering_vectors(st, hkl=[[1, 0, 0], [0, 2, 0]])
    np.testing.assert_allclose(Q, [[2 * np.pi / 4, 0, 0], [0, 4 * np.pi / 4, 0]], rtol=1e-14, atol=1e-15)
    # the same physical Q in supercell r.l.u.
    Qs = scattering_vectors(st, hkl=[[2, 0, 0], [0, 6, 0]], cell="supercell")
    np.testing.assert_allclose(Qs, Q, rtol=1e-14, atol=1e-15)
    with pytest.raises(ValueError):
        scattering_vectors(st, hkl=[[1, 0, 0]], q=[[1, 0, 0]])


def test_phase_sign_is_plus_i():
    # one atom at x = 0.25: F(h=1) = f exp(+i pi/2) = i f
    st = _single("A", frac=(0.25, 0, 0))
    F = structure_factor(st, NEUTRON, hkl=[[1, 0, 0]])
    np.testing.assert_allclose(F, [5.0j], atol=1e-12)


# ------------------------------------------------------------------------------------- single atom
@pytest.mark.parametrize("b", [None, 0.7])
def test_single_atom_xray(b):
    st = _single("Cu", b=b)
    hkl = _random_hkl()
    s = np.linalg.norm(scattering_vectors(st, hkl=hkl), axis=1) / (4 * np.pi)
    expected = _wk_value("Cu", s) ** 2 * (np.exp(-2 * b * s ** 2) if b else 1.0)
    got = intensity(st, XRAY, hkl=hkl, debye_waller=b is not None)
    np.testing.assert_allclose(got, expected, rtol=1e-13)
    assert abs(intensity(st, XRAY, hkl=[[0, 0, 0]])[0] - _wk_value("Cu", 0.0) ** 2) < 1e-10


def test_single_atom_neutron_is_flat():
    st = _single("Au")
    np.testing.assert_allclose(intensity(st, NEUTRON, hkl=_random_hkl()), 7.63 ** 2, rtol=1e-14)


# ------------------------------------------------------------------------------------- two atoms
def test_two_atom_interference():
    lat = _cubic(5.0)
    fa, fb = 5.0, -3.0
    st = Structure(lat, ["A", "B"], [[0.1, 0.1, 0.1], [0.37, 0.22, 0.61]], model="configuration")
    hkl = _random_hkl(200, seed=3)
    Q = scattering_vectors(st, hkl=hkl)
    d = st.cartesian[1] - st.cartesian[0]
    expected = fa ** 2 + fb ** 2 + 2 * fa * fb * np.cos(Q @ d)
    np.testing.assert_allclose(intensity(st, NEUTRON, hkl=hkl), expected, rtol=1e-12, atol=1e-11)


def test_occupancy_weights_amplitude():
    st = Structure(_cubic(), ["Cu"], [[0, 0, 0]], model="average", occupancy=[0.25])
    np.testing.assert_allclose(structure_factor(st, NEUTRON, hkl=[[1, 2, 3]]), [0.25 * 7.718], rtol=1e-14)


# ------------------------------------------------------------------------------------- perfect crystal
def _laue(h, n):
    h = np.asarray(h, float)
    num, den = np.sin(np.pi * n * h) ** 2, np.sin(np.pi * h) ** 2
    out = np.full(h.shape, float(n * n))
    ok = den > 1e-24
    out[ok] = num[ok] / den[ok]
    return out


def test_perfect_crystal_laue_function():
    n = (5, 4, 3)
    cell = Structure(_cubic(3.0), ["A"], [[0, 0, 0]], model="configuration")
    st = repeat(cell, n).structure
    hkl = _random_hkl(300, seed=7, scale=2.0)
    expected = 25.0 * _laue(hkl[:, 0], n[0]) * _laue(hkl[:, 1], n[1]) * _laue(hkl[:, 2], n[2])
    np.testing.assert_allclose(intensity(st, NEUTRON, hkl=hkl), expected, rtol=1e-9, atol=1e-9)


def test_fcc_bragg_and_extinctions():
    fcc = Structure(_cubic(3.6), ["Cu"] * 4, [[0, 0, 0], [.5, .5, 0], [.5, 0, .5], [0, .5, .5]],
                    model="configuration", site_label=["s1", "s2", "s3", "s4"])
    st = repeat(fcc, (3, 3, 3)).structure
    F = structure_factor(st, NEUTRON, hkl=[[1, 1, 1], [2, 0, 0], [1, 0, 0], [1, 1, 0]])
    np.testing.assert_allclose(F, [27 * 4 * 7.718, 27 * 4 * 7.718, 0, 0], atol=1e-9)


def test_average_equals_total_for_perfect_crystal():
    cell = Structure(_cubic(3.0), ["A", "B"], [[0, 0, 0], [.5, .5, .5]], model="configuration",
                     site_label=["s1", "s2"])
    st = repeat(cell, (3, 4, 2)).structure
    hkl = _random_hkl(80, seed=11)
    np.testing.assert_allclose(average_structure_factor(st, NEUTRON, hkl=hkl),
                               structure_factor(st, NEUTRON, hkl=hkl), rtol=1e-10, atol=1e-9)
    assert np.max(diffuse_intensity(st, NEUTRON, hkl=hkl)) < 1e-9


# ------------------------------------------------------------------------------------- correlated disorder
def _alloy(n=8, seed=5, c=0.4, order=0.0):
    """Simple-cubic binary A/B alloy on an n^3 periodic crystal, optionally with nearest-neighbour order."""
    rng = np.random.default_rng(seed)
    occ = (rng.random((n, n, n)) < c).astype(float)
    if order:                                     # a few Metropolis sweeps of an Ising swap model
        for _ in range(20 * n ** 3):
            i, j = rng.integers(0, n, (2, 3))
            if occ[tuple(i)] == occ[tuple(j)]:
                continue
            def e(o):
                return sum(o[tuple(p)] * (np.roll(o, s, axis=ax)[tuple(p)])
                           for p in (i, j) for ax in range(3) for s in (1, -1))
            e0 = e(occ)
            occ[tuple(i)], occ[tuple(j)] = occ[tuple(j)], occ[tuple(i)]
            if rng.random() >= np.exp(-order * (e(occ) - e0)):
                occ[tuple(i)], occ[tuple(j)] = occ[tuple(j)], occ[tuple(i)]
    cell = Structure(_cubic(3.0), ["A"], [[0, 0, 0]], model="configuration", site_label=["s1"])
    st = repeat(cell, (n, n, n)).structure
    ci = st.cell_index
    species = np.where(occ[ci[:, 0], ci[:, 1], ci[:, 2]] > 0.5, "A", "B")
    return st.replace(species=species), occ


@pytest.mark.parametrize("order", [0.0, 1.5])
def test_diffuse_equals_warren_cowley_sum(order):
    """|F − ⟨F⟩|² = (f_A − f_B)² Σ_R Σ_i σ_i σ_{i+R} e^{−iQ·R}, with σ = indicator − c.

    Exact on the grid commensurate with the n³ crystal (h·n integer), where the periodic pair sum over
    R ∈ [0, n)³ equals the sum over all real pair separations."""
    n = 6
    st, occ = _alloy(n=n, order=order)
    c = occ.mean()
    sig = occ - c
    hkl = np.random.default_rng(2).integers(-n, n + 1, (60, 3)) / n
    # pair-correlation sums by FFT-free rolling: C[R] = Σ_i σ_i σ_{i+R}
    lag = np.arange(n)
    total = np.zeros(len(hkl), complex)
    for rx in lag:
        for ry in lag:
            for rz in lag:
                C = np.sum(sig * np.roll(sig, (-rx, -ry, -rz), axis=(0, 1, 2)))
                total += C * np.exp(-2j * np.pi * (hkl @ np.array([rx, ry, rz], float)))
    expected = (5.0 - (-3.0)) ** 2 * total.real
    got = diffuse_intensity(st, NEUTRON, hkl=hkl)
    np.testing.assert_allclose(got, expected, rtol=1e-9, atol=1e-8)


def test_nearest_neighbour_correlation_matches_package_c():
    """The roll-based C(100) used above agrees with Package C's Pearson correlation r_AA."""
    n = 6
    st, occ = _alloy(n=n, order=1.5)
    c = occ.mean()
    roll = np.mean((occ - c) * (np.roll(occ, -1, axis=0) - c)) / (c * (1 - c))
    nl = neighbor_list(st, 3.1)
    sro = occupational_correlation(nl.select(translation=(1, 0, 0), symmetric=True), st)
    assert abs(sro.value("pearson", "A", "A") - roll) < 1e-12
    assert roll < -0.05                              # the Ising model above produced real order


def test_bragg_points_of_random_alloy_vanish_after_subtraction():
    st, _ = _alloy(n=5, seed=9)
    hkl = np.array([[1, 0, 0], [2, 1, 0], [0, 0, 0], [1, 1, 1]], float)
    np.testing.assert_allclose(diffuse_intensity(st, NEUTRON, hkl=hkl), 0.0, atol=1e-9)


def test_normalizations():
    st, _ = _alloy(n=4)
    hkl = _random_hkl(10)
    I = intensity(st, NEUTRON, hkl=hkl)
    np.testing.assert_allclose(intensity(st, NEUTRON, hkl=hkl, normalize="per_atom"), I / 64)
    np.testing.assert_allclose(intensity(st, NEUTRON, hkl=hkl, normalize="per_cell"), I / 64)
    with pytest.raises(ValueError):
        intensity(st, NEUTRON, hkl=hkl, normalize="per_volume")


def test_chunking_does_not_change_result():
    st, _ = _alloy(n=4)
    hkl = _random_hkl(37)
    np.testing.assert_array_equal(structure_factor(st, NEUTRON, hkl=hkl, chunk=5),
                                  structure_factor(st, NEUTRON, hkl=hkl, chunk=37))
    np.testing.assert_allclose(structure_factor(st, NEUTRON, hkl=hkl, chunk=1),
                               structure_factor(st, NEUTRON, hkl=hkl), rtol=1e-12)


def test_grid_shape_is_kept():
    st = _single()
    hkl = np.stack(np.meshgrid(*(np.linspace(-1, 1, 4),) * 3, indexing="ij"), axis=-1)
    assert structure_factor(st, NEUTRON, hkl=hkl).shape == (4, 4, 4)


# ------------------------------------------------------------------------------------- refusals
def test_refusals():
    st = _single("Cu")
    with pytest.raises(KeyError):
        structure_factor(_single("Fe"), NEUTRON, hkl=[[1, 0, 0]])
    with pytest.raises(ValueError, match="b_iso"):
        structure_factor(st, NEUTRON, hkl=[[1, 0, 0]], debye_waller=True)
    with pytest.raises(ValueError, match="configuration"):
        average_structure_factor(Structure(_cubic(), ["Cu"], [[0, 0, 0]], model="average"), NEUTRON,
                                 hkl=[[1, 0, 0]])
    with pytest.raises(ValueError, match="radiation"):
        ScatteringFactors("electron", {"Cu": 1.0}, "x")
    with pytest.raises(ValueError, match="source"):
        ScatteringFactors.neutron({"Cu": 1.0}, source="")
    with pytest.raises(ValueError, match="twice"):
        ScatteringFactors.neutron({"Cu": 1.0, "CU": 2.0})
    with pytest.raises(ValueError, match="real"):
        ScatteringFactors.neutron({"Cu": 1 + 1j})
    with pytest.raises(ValueError, match="function"):
        ScatteringFactors("xray", {"Cu": 29.0}, "x")


def test_case_insensitive_species_lookup():
    st = _single("CU")
    np.testing.assert_allclose(intensity(st, NEUTRON, hkl=[[1, 0, 0]]), 7.718 ** 2)


# ------------------------------------------------------------------------------------- optional tables
def test_xraydb_agrees_with_dabax_transcription():
    xraydb = pytest.importorskip("xraydb")
    fac = ScatteringFactors.from_xraydb(["Cu", "AU"])
    s = np.linspace(0, 2.0, 41)
    for label, el in (("Cu", "CU"), ("AU", "AU")):
        np.testing.assert_allclose(fac.factor(label, s), _wk_value(el, s), rtol=2e-6, atol=1e-6)
    assert "xraydb" in fac.source
    with pytest.raises(ValueError, match="ion"):
        ScatteringFactors.from_xraydb(["O2-"])
    del xraydb


def test_periodictable_lengths():
    pytest.importorskip("periodictable")
    fac = ScatteringFactors.from_periodictable(["Cu", "Ni"])
    assert fac.radiation == "neutron"
    assert abs(fac.factor("Cu", 0.0) - 7.718) < 0.01        # Sears (1992) 7.718 fm
    assert fac.factor("Ni", 0.0) > 10                       # Ni 10.3 fm


def test_diffuse_off_grid_equals_brute_force_pair_sum():
    """Off the commensurate grid: |F − ⟨F⟩|² = Δf² Σ_ij σ_i σ_j e^{iQ·(r_i − r_j)} (O(N²) double sum)."""
    st, occ = _alloy(n=4, seed=13, order=1.0)
    c = occ.mean()
    ci = st.cell_index
    sig = occ[ci[:, 0], ci[:, 1], ci[:, 2]] - c
    hkl = _random_hkl(25, seed=4, scale=1.3)
    Q = scattering_vectors(st, hkl=hkl)
    r = st.cartesian
    ph = np.exp(1j * (Q @ r.T))                              # (M, N)
    expected = 8.0 ** 2 * np.abs(ph @ sig) ** 2              # Σ_ij σ_i σ_j e^{iQ(ri−rj)} = |Σ σ e^{iQr}|²
    np.testing.assert_allclose(diffuse_intensity(st, NEUTRON, hkl=hkl), expected, rtol=1e-9, atol=1e-8)


# ------------------------------------------------------------------------------------- general geometry
def _triclinic_crystal(seed=21, vacancies=3):
    """3 sites in a triclinic cell, 3×2×4 cells, B ≠ 0, displacements, a few vacancies."""
    lat = Lattice.from_parameters(4.3, 5.1, 3.9, 77, 103, 95)
    cell = Structure(lat, ["A", "B", "A"], [[0.05, 0.1, 0.2], [0.55, 0.45, 0.6], [0.3, 0.8, 0.1]],
                     model="configuration", site_label=["s1", "s2", "s3"], b_iso=[0.4, 0.9, 0.4])
    st = repeat(cell, (3, 2, 4)).structure
    rng = np.random.default_rng(seed)
    st = st.with_cartesian(st.cartesian + rng.normal(0, 0.08, (st.n_atoms, 3)))
    drop = rng.choice(st.n_atoms, vacancies, replace=False)
    return st.remove(st.atom_id[drop])


def _brute(st, hkl, fac, average=False):
    """Independent loops: parent-fraction phases exp(2πi h·x_parent), s from the reciprocal metric."""
    n = st.supercell
    Ap = st.lattice.matrix / n[:, None]
    Gstar = np.linalg.inv(Ap @ Ap.T)
    xpar = st.frac * n                                    # parent-cell fractional coordinates
    out = np.zeros(len(hkl), complex)
    for m, h in enumerate(hkl):
        s = 0.5 * np.sqrt(h @ Gstar @ h)
        if not average:
            for sp, x, b in zip(st.species, xpar, st.b_iso):
                out[m] += fac.factor(sp, s) * np.exp(-b * s * s) * np.exp(2j * np.pi * h @ x)
            continue
        ncell = int(np.prod(n))
        cellsum = 0j
        for sp, x, b, ci in zip(st.species, xpar, st.b_iso, st.cell_index):
            cellsum += fac.factor(sp, s) * np.exp(-b * s * s) * np.exp(2j * np.pi * h @ (x - ci))
        L = sum(np.exp(2j * np.pi * h @ np.array([i, j, k], float))
                for i in range(n[0]) for j in range(n[1]) for k in range(n[2]))
        out[m] = L * cellsum / ncell
    return out


def test_triclinic_vacancies_against_brute_force():
    st = _triclinic_crystal()
    hkl = _random_hkl(30, seed=8, scale=2.0)
    np.testing.assert_allclose(structure_factor(st, NEUTRON, hkl=hkl, debye_waller=True),
                               _brute(st, hkl, NEUTRON), rtol=1e-10, atol=1e-10)
    np.testing.assert_allclose(average_structure_factor(st, NEUTRON, hkl=hkl, debye_waller=True),
                               _brute(st, hkl, NEUTRON, average=True), rtol=1e-10, atol=1e-10)


def test_supercell_rlu_in_triclinic_cell():
    st = _triclinic_crystal(vacancies=0)
    hkl = _random_hkl(10, seed=9)
    np.testing.assert_allclose(structure_factor(st, NEUTRON, hkl=hkl * st.supercell, cell="supercell"),
                               structure_factor(st, NEUTRON, hkl=hkl), rtol=1e-10, atol=1e-10)


def test_normalizations_multi_site():
    st = _triclinic_crystal(vacancies=3)                  # 72 sites − 3 vacancies, 24 cells
    hkl = _random_hkl(5)
    I = intensity(st, NEUTRON, hkl=hkl)
    np.testing.assert_allclose(intensity(st, NEUTRON, hkl=hkl, normalize="per_atom"), I / 69)
    np.testing.assert_allclose(intensity(st, NEUTRON, hkl=hkl, normalize="per_cell"), I / 24)


def test_single_cell_average_refused():
    st = _alloy(n=1)[0]
    with pytest.raises(ValueError, match="single cell"):
        diffuse_intensity(st, NEUTRON, hkl=[[0.5, 0, 0]])


def test_cell_argument_validated_with_q():
    with pytest.raises(ValueError, match="cell"):
        scattering_vectors(_single(), q=[[1, 0, 0]], cell="bogus")
