"""Analytical checks of the experimental magnetic neutron reference (Package D, Feature 11)."""
from __future__ import annotations

import numpy as np
import pytest

from discuspdf.lattice import Lattice
from discuspdf.magnetic_reference import (P_MAGNETIC_FM, MagneticFormFactor, magnetic_intensity,
                                          magnetic_structure_factor)

P2 = P_MAGNETIC_FM ** 2
ONE = lambda q: np.ones_like(np.asarray(q, float))      # noqa: E731  f = 1 isolates geometry


def test_prefactor_value():
    assert abs(P_MAGNETIC_FM - 2.6954) < 1e-4            # γ r0 / 2 in fm


@pytest.mark.parametrize("theta_deg", [0, 30, 60, 90, 135])
def test_single_moment_sin2(theta_deg):
    th = np.deg2rad(theta_deg)
    mu = 2.5 * np.array([np.cos(th), np.sin(th), 0.0])
    q = np.array([[1.7, 0.0, 0.0]])                       # Q along x
    I, ok = magnetic_intensity([[0.3, -0.2, 0.9]], [mu], q, ONE)
    assert ok.all()
    np.testing.assert_allclose(I, P2 * 2.5 ** 2 * np.sin(th) ** 2, rtol=1e-12, atol=1e-12)


def test_moment_parallel_to_q_in_a_nonorthogonal_cell_gives_zero():
    # catches mixing r.l.u. components with Cartesian |Q| (the Javelin bug): a = 5.3 Å, triclinic
    lat = Lattice.from_parameters(5.3, 6.1, 4.7, 81, 97, 113)
    Q = lat.hkl_to_q([[1, 2, -1], [0.5, 0, 1.5]])
    mus = 1.7 * Q / np.linalg.norm(Q, axis=1)[:, None]
    for q, mu in zip(Q, mus):
        I, _ = magnetic_intensity([[0, 0, 0]], [mu], q[None], ONE)
        assert abs(I[0]) < 1e-24
        perp = np.cross(q, [0.3, 1.0, -0.2])
        perp *= 1.7 / np.linalg.norm(perp)
        I2, _ = magnetic_intensity([[0, 0, 0]], [perp], q[None], ONE)
        np.testing.assert_allclose(I2, P2 * 1.7 ** 2, rtol=1e-12)


def test_ferro_and_antiferro_pairs():
    rng = np.random.default_rng(1)
    Q = rng.normal(size=(40, 3)) * 2
    d = np.array([1.1, 0.4, -0.7])
    mu = np.array([0.0, 0.0, 3.0])
    sin2 = 1 - (Q @ mu / (np.linalg.norm(Q, axis=1) * 3.0)) ** 2
    Ifm, _ = magnetic_intensity([[0, 0, 0], d], [mu, mu], Q, ONE)
    Iaf, _ = magnetic_intensity([[0, 0, 0], d], [mu, -mu], Q, ONE)
    np.testing.assert_allclose(Ifm, P2 * 9 * sin2 * 4 * np.cos(Q @ d / 2) ** 2, rtol=1e-11, atol=1e-11)
    np.testing.assert_allclose(Iaf, P2 * 9 * sin2 * 4 * np.sin(Q @ d / 2) ** 2, rtol=1e-11, atol=1e-11)


def test_antiferromagnetic_chain_peaks_at_half_integers():
    a, N = 3.0, 8
    lat = Lattice.from_parameters(a, 10, 10, 90, 90, 90)
    r = np.column_stack([np.arange(N) * a, np.zeros(N), np.zeros(N)])
    mu = np.zeros((N, 3))
    mu[:, 2] = np.where(np.arange(N) % 2 == 0, 1.0, -1.0)
    h = np.array([[0.5, 0, 0], [1.5, 0, 0], [1.0, 0, 0], [2.0, 0, 0]])
    I, _ = magnetic_intensity(r, mu, lat.hkl_to_q(h), ONE)
    np.testing.assert_allclose(I[:2], P2 * N ** 2, rtol=1e-12)
    np.testing.assert_allclose(I[2:], 0.0, atol=1e-20)
    Ifm, _ = magnetic_intensity(r, np.abs(mu), lat.hkl_to_q(h), ONE)
    np.testing.assert_allclose(Ifm[2:], P2 * N ** 2, rtol=1e-12)
    np.testing.assert_allclose(Ifm[:2], 0.0, atol=1e-20)


def test_powder_average_two_thirds():
    # Q directions on a Fibonacci sphere: mean sin² → 2/3
    n = 20000
    k = np.arange(n) + 0.5
    phi, ct = np.pi * (1 + 5 ** 0.5) * k, 1 - 2 * k / n
    st = np.sqrt(1 - ct ** 2)
    Q = 2.0 * np.column_stack([st * np.cos(phi), st * np.sin(phi), ct])
    I, _ = magnetic_intensity([[0, 0, 0]], [[0.0, 0.0, 1.0]], Q, ONE)
    assert abs(I.mean() / P2 - 2 / 3) < 1e-4


def test_q_zero_is_nan_and_masked():
    Q = np.array([[0.0, 0, 0], [1e-12, 0, 0], [1.0, 0, 0]])
    I, ok = magnetic_intensity([[0, 0, 0]], [[0, 1.0, 0]], Q, ONE)
    assert ok.tolist() == [False, False, True]
    assert np.isnan(I[:2]).all() and np.isfinite(I[2])
    F, ok2 = magnetic_structure_factor([[0, 0, 0]], [[0, 1.0, 0]], Q, ONE)
    assert np.isnan(F[:2]).all() and F.shape == (3, 3)


def test_form_factor_formula_against_periodictable():
    pt = pytest.importorskip("periodictable")
    ff = MagneticFormFactor.from_periodictable("V", 3)
    q = np.linspace(0, 8, 33)
    np.testing.assert_allclose(ff(q), pt.V.magnetic_ff[3].j0_Q(q), rtol=1e-12)
    # Javelin's own test value: j0(V3+, 2π/5 Å) = 0.88829
    assert abs(ff(2 * np.pi / 5) - 0.88829) < 1e-5
    ffg = MagneticFormFactor.from_periodictable("Fe", 3, g=1.9)
    expect = pt.Fe.magnetic_ff[3].j0_Q(q) + (1 - 2 / 1.9) * pt.Fe.magnetic_ff[3].j2_Q(q)
    np.testing.assert_allclose(ffg(q), expect, rtol=1e-12)


def test_form_factor_is_one_at_q_zero_for_normalised_coefficients():
    ff = MagneticFormFactor((0.4, 15.0, 0.6, 5.0, -0.1, 1.0, 0.1), source="test")
    assert abs(ff(0.0) - 1.0) < 1e-12


def test_form_factor_enters_squared():
    ff = MagneticFormFactor((0.4, 15.0, 0.6, 5.0, -0.1, 1.0, 0.1), source="test")
    q = np.array([[0, 2.2, 0]])
    I, _ = magnetic_intensity([[0, 0, 0]], [[1.0, 0, 0]], q, {"Mn": ff}, species=["Mn"])
    np.testing.assert_allclose(I, P2 * ff(2.2) ** 2, rtol=1e-12)


def test_refusals():
    with pytest.raises(ValueError, match="j2"):
        MagneticFormFactor((1, 0, 0, 0, 0, 0, 0), source="x", g=1.9)
    with pytest.raises(ValueError, match="source"):
        MagneticFormFactor((1, 0, 0, 0, 0, 0, 0), source="")
    with pytest.raises(ValueError, match="7"):
        MagneticFormFactor((1, 0, 0), source="x")
    with pytest.raises(KeyError):
        magnetic_intensity([[0, 0, 0]], [[1, 0, 0]], [[1, 0, 0]], {"Fe": ONE}, species=["Mn"])
    with pytest.raises(ValueError, match="species"):
        magnetic_intensity([[0, 0, 0]], [[1, 0, 0]], [[1, 0, 0]], {"Fe": ONE})
    with pytest.raises(ValueError, match="N, 3"):
        magnetic_intensity([[0, 0, 0]], [[1, 0]], [[1, 0, 0]], ONE)
    with pytest.raises(ValueError, match="finite"):
        magnetic_intensity([[0, 0, np.nan]], [[1, 0, 0]], [[1, 0, 0]], ONE)
