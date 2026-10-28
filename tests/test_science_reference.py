"""Analytical tests for the independent scientific reference, not production FFT."""
import os
import sys
from pathlib import Path
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]
                       / "benchmarks" / "discus_6_23_01"))
from reference import (
    hkl_grid, direct_average, factored_average, xray_factor,
    metrics, load_c1_structure, empirical_basis_probabilities,
)


def test_skewed_grid_uses_entire_increment_matrix():
    inc = np.array([[.11, .02, 0], [-.01, .13, .03], [.05, 0, .17]])
    corner = np.array([-.2, .3, .4])
    grid = hkl_grid((2, 3, 4), corner, inc)
    assert grid.shape == (2, 3, 4, 3)
    np.testing.assert_allclose(grid[1, 2, 3],
                               corner + inc @ np.array([1, 2, 3]))
    np.testing.assert_allclose(grid[0, 1, 0], corner + inc[:, 1])
    np.testing.assert_allclose(grid[1, 0, 0], corner + inc[:, 0])
    assert not np.isclose(grid[1, 2, 3, 0], corner[0] + .11)


@pytest.mark.parametrize("increment", [
    np.eye(2), np.zeros((3, 3)),
    [[1, 0, 0], [0, 1, 0], [0, 0, np.nan]],
])
def test_bad_grid_rejected(increment):
    with pytest.raises(ValueError):
        hkl_grid((2, 3, 4), [0, 0, 0], increment)


def test_2x2x2_analytical_oracle_and_origin_shift():
    cells = np.stack(np.meshgrid(*([np.arange(2)] * 3),
                                 indexing="ij"), -1).reshape(-1, 3)
    basis = np.array([.25, .5, .125])
    positions = cells + basis
    hkl = np.array([[0, 0, 0], [1, 0, 0], [.13, .27, .41],
                    [.5, .5, .5], [-.1, .12, 0]])
    direct = direct_average(hkl, positions, np.full(8, .25),
                            a=5.6537, biso=np.ones(8))
    factored = factored_average(hkl, {tuple(basis): .25}, a=5.6537,
                                n_min=0, n_max=1)
    np.testing.assert_allclose(direct, factored, rtol=1e-12, atol=1e-11)
    expected = 8 * (.25 * xray_factor("GA", [0])[0]
                    + .75 * xray_factor("AS", [0])[0])
    assert direct[0] == pytest.approx(expected)
    shift = np.array([.2, -.3, .4])
    shifted = direct_average(hkl, positions + shift, np.full(8, .25),
                             a=5.6537, biso=np.ones(8))
    np.testing.assert_allclose(
        shifted, direct * np.exp(2j * np.pi * (hkl @ shift)),
        rtol=2e-12, atol=2e-11)


def test_prescribed_occupancy_is_not_sample_sum_over_n():
    positions = np.array([[0, 0, 0], [.5, 0, 0.]])
    out = direct_average([[0, 0, 0], [.3, 0, 0]], positions, [1., 0.],
                         a=5.6537, biso=np.ones(2))
    expected = xray_factor("GA", [0])[0] + xray_factor("AS", [0])[0]
    assert out[0] == pytest.approx(expected)
    with pytest.raises(ValueError):
        direct_average([[0, 0, 0]], positions, [1.2, 0],
                       a=5.6537, biso=np.ones(2))


def test_metrics_safely_skip_zero_relative_denominators():
    observed = np.array([[[0., 1.]]])
    predicted = np.array([[[3., 1.]]])
    grid = np.zeros((1, 1, 2, 3))
    m = metrics(predicted, observed, grid)
    assert m["above_floor"] == 1
    assert m["worst_index"] == [0, 0, 0]


def test_real_C1_direct_vs_factored_opt_in():
    path = os.environ.get("PR8_C1_FROZEN_STRU")
    if not path:
        pytest.skip("set PR8_C1_FROZEN_STRU to archived C1_frozen.stru")
    species, positions, biso = load_c1_structure(Path(path))
    basis = empirical_basis_probabilities(species, positions)
    p = np.array([basis[tuple(b)] for b in np.round(positions % 1, 4)])
    hkl = hkl_grid((2, 3, 4), [-.1, .02, .07],
                   np.diag([.1, .06, .11])).reshape(-1, 3)
    hkl = np.vstack([hkl, [[0, 0, 0], [1, 0, 0], [.07, .03, .01]]])
    actual_direct = direct_average(hkl, positions, p, a=5.6537,
                                   biso=biso, chunk=8)
    predicted = factored_average(hkl, basis, a=5.6537)
    np.testing.assert_allclose(actual_direct, predicted,
                               rtol=5e-11, atol=5e-8)
