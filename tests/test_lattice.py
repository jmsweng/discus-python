"""Tests for the Lattice core (Package A, Feature A5).

Expected values come from independent computations (closed-form formulas, literal
numbers, brute-force enumeration), never from calling the function under test a
second way.
"""
import math

import numpy as np
import pytest

from discuspdf import Lattice
from discuspdf.lattice import _canonical_basis
from discuspdf.discus_output import DiscusVolume, _validated_direct_basis


def _discus_volume_with_cell(cell, angles, n=8, d=0.25):
    """A minimal synthetic DiscusVolume (axis-aligned, centered) for parity checks."""
    values = np.ones((n, n, n))
    return DiscusVolume(
        values,
        np.array([-(n // 2) * d] * 3),
        np.eye(3) * d,
        np.asarray(cell, float),
        np.asarray(angles, float),
        "bragg_subtracted",
    )


# ------------------------------ 1. cubic ------------------------------
def test_cubic_q_d_volume():
    lat = Lattice.from_parameters(5.0, 5.0, 5.0, 90.0, 90.0, 90.0)
    assert lat.volume == pytest.approx(125.0, rel=1e-14)
    q100 = lat.hkl_to_q((1, 0, 0))
    assert np.linalg.norm(q100) == pytest.approx(2 * math.pi / 5.0, rel=1e-13)
    assert lat.d_spacing((1, 0, 0)) == pytest.approx(5.0, rel=1e-13)
    assert lat.d_spacing((1, 1, 1)) == pytest.approx(5.0 / math.sqrt(3.0), rel=1e-13)
    # canonical orientation: a ∥ x
    assert np.allclose(lat.matrix, np.diag([5.0, 5.0, 5.0]), atol=1e-14)


# ------------------------------ 2. triclinic ------------------------------
TRIC = (4.1, 5.3, 6.7, 78.0, 95.0, 112.0)


def _triclinic_volume_closed_form(a, b, c, alpha, beta, gamma):
    ca, cb, cg = (math.cos(math.radians(x)) for x in (alpha, beta, gamma))
    return a * b * c * math.sqrt(1 - ca * ca - cb * cb - cg * cg + 2 * ca * cb * cg)


def _triclinic_d2_inv(h, k, l, a, b, c, alpha, beta, gamma):
    """Textbook 1/d² = (h²a*² + k²b*² + l²c*² + 2hka*b* + 2hlc*a*c* + 2klb*c*) / (a b c)²...
    in the standard form: 1/d² = (h² sin²γ + k² sin²α + l² sin²β
                               + 2hk cosγ cosα... ) / V²  — use the metric form:
    1/d² = hkl G* hklᵀ with G* = inv(G), G from the cell."""
    ca, cb, cg = (math.cos(math.radians(x)) for x in (alpha, beta, gamma))
    G = np.array([[1.0, cg, cb], [cg, 1.0, ca], [cb, ca, 1.0]])
    G = (G * np.outer([a, b, c], [a, b, c]))
    hkl = np.array([h, k, l], float)
    return float(hkl @ np.linalg.inv(G) @ hkl)


def test_triclinic_volume_parameters_dspacing():
    a, b, c, alpha, beta, gamma = TRIC
    lat = Lattice.from_parameters(*TRIC)
    v_closed = _triclinic_volume_closed_form(*TRIC)
    assert lat.volume == pytest.approx(v_closed, rel=1e-12)
    # recover all six parameters
    p = lat.parameters
    assert p[0] == pytest.approx(a, rel=1e-12)
    assert p[1] == pytest.approx(b, rel=1e-12)
    assert p[2] == pytest.approx(c, rel=1e-12)
    assert p[3] == pytest.approx(alpha, rel=1e-10)
    assert p[4] == pytest.approx(beta, rel=1e-10)
    assert p[5] == pytest.approx(gamma, rel=1e-10)
    # d_hkl from the textbook 1/d² formula for several hkl
    for hkl in ((1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (1, 1, 1), (2, -1, 3)):
        d = lat.d_spacing(hkl)
        d_closed = 1.0 / math.sqrt(_triclinic_d2_inv(*hkl, *TRIC))
        assert d == pytest.approx(d_closed, rel=1e-11)


def test_triclinic_basis_identities():
    a, b, c, alpha, beta, gamma = TRIC
    A = _validated_direct_basis((a, b, c), (alpha, beta, gamma))
    B = np.linalg.inv(A).T
    assert np.allclose(A @ B.T, np.eye(3), atol=1e-13)
    G = A @ A.T
    Gstar = B @ B.T
    assert np.allclose(G @ Gstar, np.eye(3), atol=1e-10)
    # Lattice agrees with the existing validated helper
    lat = Lattice.from_parameters(*TRIC)
    assert np.allclose(lat.matrix, A, rtol=1e-14, atol=1e-15)
    assert np.allclose(lat.reciprocal_matrix, B, rtol=1e-13, atol=1e-15)
    assert np.allclose(lat.metric, G, rtol=1e-13)
    assert np.allclose(lat.reciprocal_metric, Gstar, rtol=1e-12)


# ------------------------------ 3. 2π identity ------------------------------
def test_twopi_identity():
    rng = np.random.default_rng(42)
    lat = Lattice.from_parameters(*TRIC)
    for _ in range(50):
        hkl = rng.uniform(-5, 5, 3)
        f = rng.uniform(-2, 2, 3)
        lhs = lat.hkl_to_q(hkl) @ lat.frac_to_cart(f)
        rhs = 2 * math.pi * float(hkl @ f)
        assert lhs == pytest.approx(rhs, rel=1e-12)


# ------------------------------ 4. round trips ------------------------------
def test_round_trips():
    lat = Lattice.from_parameters(*TRIC)
    rng = np.random.default_rng(7)
    f = rng.uniform(-3, 3, 3)
    r = lat.frac_to_cart(f)
    assert np.allclose(lat.cart_to_frac(r), f, rtol=1e-12, atol=1e-14)
    hkl = np.array([2.0, -1.0, 3.0])
    q = lat.hkl_to_q(hkl)
    assert np.allclose(lat.q_to_hkl(q), hkl, rtol=1e-12, atol=1e-14)
    # batched shapes
    fb = rng.uniform(-3, 3, (4, 3))
    assert np.allclose(lat.cart_to_frac(lat.frac_to_cart(fb)), fb, rtol=1e-12, atol=1e-14)
    hb = rng.uniform(-5, 5, (2, 5, 3))
    assert np.allclose(lat.q_to_hkl(lat.hkl_to_q(hb)), hb, rtol=1e-12, atol=1e-14)
    # scalar and (N,3) hkl
    q1 = lat.hkl_to_q((1, 0, 0))
    assert q1.shape == (3,)
    qn = lat.hkl_to_q(np.arange(6).reshape(2, 3).astype(float))
    assert qn.shape == (2, 3)
    # (2,N,3) inputs
    fb2 = rng.uniform(-3, 3, (2, 4, 3))
    assert np.allclose(lat.cart_to_frac(lat.frac_to_cart(fb2)), fb2, rtol=1e-12, atol=1e-14)
    hb2 = rng.uniform(-5, 5, (2, 4, 3))
    assert np.allclose(lat.q_to_hkl(lat.hkl_to_q(hb2)), hb2, rtol=1e-12, atol=1e-14)


# ------------------------------ 5. rejections ------------------------------
def test_rejections():
    with pytest.raises(ValueError):
        Lattice.from_parameters(0.0, 5.0, 5.0, 90, 90, 90)
    with pytest.raises(ValueError):
        Lattice.from_parameters(-5.0, 5.0, 5.0, 90, 90, 90)
    with pytest.raises(ValueError):
        Lattice.from_parameters(float("nan"), 5.0, 5.0, 90, 90, 90)
    with pytest.raises(ValueError):
        Lattice.from_parameters(float("inf"), 5.0, 5.0, 90, 90, 90)
    with pytest.raises(ValueError):
        Lattice.from_parameters(5.0, 5.0, 5.0, 0.0, 90, 90)
    with pytest.raises(ValueError):
        Lattice.from_parameters(5.0, 5.0, 5.0, 180.0, 90, 90)
    with pytest.raises(ValueError):
        Lattice.from_parameters(5.0, 5.0, 5.0, -10.0, 90, 90)
    with pytest.raises(ValueError):
        Lattice.from_parameters(5.0, 5.0, 5.0, float("nan"), 90, 90)
    # impossible angle set (10°, 10°, 170°): metric not positive definite
    with pytest.raises(ValueError):
        Lattice.from_parameters(5.0, 5.0, 5.0, 10.0, 10.0, 170.0)
    # near-singular cell
    with pytest.raises(ValueError):
        Lattice.from_parameters(5.0, 5.0, 5.0, 90.0, 90.0, 179.9999)
    # singular and left-handed from_matrix
    with pytest.raises(ValueError):
        Lattice.from_matrix(np.zeros((3, 3)))
    with pytest.raises(ValueError):
        Lattice.from_matrix(-np.eye(3))
    with pytest.raises(ValueError):
        Lattice.from_matrix(np.diag([1.0, 1.0, 1e-15]))
    with pytest.raises(ValueError):
        Lattice.from_matrix(np.eye(3) * float("nan"))
    # wrong shapes
    with pytest.raises(ValueError):
        Lattice.from_matrix(np.ones((2, 2)))
    with pytest.raises((ValueError, TypeError)):
        Lattice.from_parameters(5.0, 5.0, 5.0, 90, 90)
    # d_spacing (0,0,0)
    lat = Lattice.from_parameters(5.0, 5.0, 5.0, 90, 90, 90)
    with pytest.raises(ValueError):
        lat.d_spacing((0, 0, 0))
    # non-finite coordinate inputs
    with pytest.raises(ValueError):
        lat.frac_to_cart([1.0, float("nan"), 0.0])
    with pytest.raises(ValueError):
        lat.hkl_to_q([1.0, float("inf"), 0.0])
    with pytest.raises(ValueError):
        lat.frac_to_cart([1.0, 2.0])


# ------------------------------ 6. parity with existing code ------------------------------
def test_parity_with_discus_output():
    # the cells used in tests/test_discus_output.py
    cells = [
        ((5.0, 5.0, 5.0), (90.0, 90.0, 90.0)),
        ((4.0, 5.0, 8.0), (90.0, 90.0, 90.0)),
        ((5.0, 6.0, 7.0), (90.0, 105.0, 90.0)),
    ]
    rng = np.random.default_rng(123)
    # 20 random valid cells
    for _ in range(20):
        a = rng.uniform(3.0, 10.0)
        b = rng.uniform(3.0, 10.0)
        c = rng.uniform(3.0, 10.0)
        alpha = rng.uniform(70.0, 110.0)
        beta = rng.uniform(70.0, 110.0)
        gamma = rng.uniform(70.0, 110.0)
        cells.append(((a, b, c), (alpha, beta, gamma)))

    worst = 0.0
    for cell, angles in cells:
        A_ref = _validated_direct_basis(cell, angles)
        lat = Lattice.from_parameters(*cell, *angles)
        worst = max(worst, float(np.max(np.abs(lat.matrix - A_ref)) / np.max(np.abs(A_ref))))
        # 2π reciprocal basis equals DiscusVolume.reciprocal_basis(two_pi=True)
        vol = _discus_volume_with_cell(cell, angles)
        assert np.allclose(2 * np.pi * lat.reciprocal_matrix, vol.reciprocal_basis(two_pi=True),
                          rtol=1e-13, atol=1e-15)
    assert worst <= 1e-14

    # rejection set matches: same inputs raise in both
    bad = [
        ((0.0, 5.0, 5.0), (90.0, 90.0, 90.0)),
        ((5.0, -5.0, 5.0), (90.0, 90.0, 90.0)),
        ((5.0, 5.0, 5.0), (0.0, 90.0, 90.0)),
        ((5.0, 5.0, 5.0), (180.0, 90.0, 90.0)),
        ((5.0, 5.0, 5.0), (10.0, 10.0, 170.0)),
        ((5.0, 5.0, 5.0), (90.0, 90.0, 179.9999)),
        ((float("nan"), 5.0, 5.0), (90.0, 90.0, 90.0)),
    ]
    for cell, angles in bad:
        with pytest.raises(ValueError):
            _validated_direct_basis(cell, angles)
        with pytest.raises(ValueError):
            Lattice.from_parameters(*cell, *angles)


# ------------------------------ 7. minimum image ------------------------------
def _brute_force_minimum_image(lat, df, window=4):
    """Independent brute-force over all shifts in a fixed ±window window."""
    df = np.asarray(df, float)
    A = lat.matrix
    best_d, best_n = None, None
    for i in range(-window, window + 1):
        for j in range(-window, window + 1):
            for k in range(-window, window + 1):
                n = np.array([i, j, k], float)
                d = (df - n) @ A
                dist = float(np.linalg.norm(d))
                if best_d is None or dist < best_d - 1e-15:
                    best_d, best_n = dist, n
    return best_d, best_n


def _brute_force_minimum_image_vec(lat, df, window=4):
    """Vectorised independent brute-force over all shifts in a fixed ±window window.

    ``df`` has shape (m, 3); returns (d_ref, n_ref) each of shape (m, ...) .
    """
    df = np.asarray(df, float)
    A = lat.matrix
    shifts = np.stack(np.meshgrid(
        *[range(-window, window + 1) for _ in range(3)], indexing="ij"), axis=-1
    ).reshape(-1, 3).astype(float)  # (K, 3)
    # (m, 1, 3) - (1, K, 3) -> (m, K, 3)
    d = df[:, None, :] - shifts[None, :, :]
    dist = np.linalg.norm(d @ A, axis=2)  # (m, K)
    best = dist.argmin(axis=1)
    return dist.min(axis=1), shifts[best]


def test_minimum_image_skewed_cells():
    # strongly skewed cell where naive per-component rounding is demonstrably wrong
    lat = Lattice.from_parameters(5.0, 5.0, 5.0, 90.0, 90.0, 20.0)
    rng = np.random.default_rng(99)
    n = 10_000
    df = rng.uniform(-2.5, 2.5, (n, 3))
    # counter-example: naive rounding is not the minimum image for at least one point
    naive = np.rint(df)
    naive_dist = np.linalg.norm((df - naive) @ lat.matrix, axis=1)
    df_min, n_shift = lat.minimum_image(df)
    min_dist = np.linalg.norm((df - n_shift) @ lat.matrix, axis=1)
    assert np.any(min_dist < naive_dist - 1e-10)
    # check ALL 10^4 against independent vectorised brute force over a ±4 window.
    # Two-sided: the code may not return a *shorter* vector than the true
    # minimum (impossible) NOR a longer one (the batched tie-resolution bug).
    d_ref, _ = _brute_force_minimum_image_vec(lat, df, window=4)
    assert np.max(np.abs(min_dist - d_ref)) <= 1e-12

    # a cell with all angles ≈ 30°: check all 500 (two-sided, same rule)
    lat30 = Lattice.from_parameters(4.0, 4.0, 4.0, 30.0, 30.0, 30.0)
    df2 = rng.uniform(-1.5, 1.5, (500, 3))
    df_min2, n2 = lat30.minimum_image(df2)
    min_dist2 = np.linalg.norm((df2 - n2) @ lat30.matrix, axis=1)
    d_ref2, _ = _brute_force_minimum_image_vec(lat30, df2, window=4)
    assert np.max(np.abs(min_dist2 - d_ref2)) <= 1e-12


def test_minimum_image_batch_equals_single_and_oracle():
    """Batched minimum_image equals per-vector calls and an independent per-row
    oracle, for the five cells from the PR #11 review.

    Cells: cubic 5 Å; hexagonal (4, 4, 4, 90, 90, 120); γ = 20° (5, 5, 5,
    90, 90, 20); (4, 4, 4, 30, 30, 30); triclinic (4.1, 5.3, 6.7, 78, 95, 112).
    For each: 2000 random df in [−3, 3] plus the 27 vectors v/2 with
    v ∈ {−1, 0, 1}³.  Asserts:
      * batched shifts == per-vector shifts (np.array_equal);
      * both == an independent per-row oracle: all shifts within ±4 of rint(df),
        shortest within 1e-12·R + 1e-15, then lexicographically smallest
        round((df − n)/1e-12).
    Also runs one batch large enough to span several chunks.
    """
    cells = [
        ((5.0, 5.0, 5.0), (90.0, 90.0, 90.0)),      # cubic 5 Å
        ((4.0, 4.0, 4.0), (90.0, 90.0, 120.0)),   # hexagonal
        ((5.0, 5.0, 5.0), (90.0, 90.0, 20.0)),    # γ = 20°
        ((4.0, 4.0, 4.0), (30.0, 30.0, 30.0)),    # all angles 30°
        ((4.1, 5.3, 6.7), (78.0, 95.0, 112.0)),   # triclinic
    ]
    rng = np.random.default_rng(20261030)

    def _oracle_row(lat, df):
        """Independent per-row brute force over shifts within ±4 of rint(df)."""
        df = np.asarray(df, float)
        A = lat.matrix
        n0 = np.rint(df)
        best_d, best_n = None, None
        for i in range(-4, 5):
            for j in range(-4, 5):
                for k in range(-4, 5):
                    n = n0 + np.array([i, j, k], float)
                    d = (df - n) @ A
                    dist = float(np.linalg.norm(d))
                    if best_d is None or dist < best_d - 1e-15:
                        best_d, best_n = dist, n
        # tie rule: among images within atol of the shortest, the reduced
        # vector (df - n) lexicographically smallest wins
        atol = 1e-12 * best_d + 1e-15
        best_key = None
        for i in range(-4, 5):
            for j in range(-4, 5):
                for k in range(-4, 5):
                    n = n0 + np.array([i, j, k], float)
                    d = (df - n) @ A
                    if float(np.linalg.norm(d)) <= best_d + atol:
                        key = tuple(np.round((df - n) / 1e-12))
                        if best_key is None or key < best_key:
                            best_key = key
                            best_n = n
        return best_n

    tie_vecs = np.array(
        [[x, y, z] for x in (-1.0, 0.0, 1.0) for y in (-1.0, 0.0, 1.0)
         for z in (-1.0, 0.0, 1.0)]
    ) / 2.0

    for cell, angles in cells:
        lat = Lattice.from_parameters(*cell, *angles)
        df = np.vstack([rng.uniform(-3.0, 3.0, (2000, 3)), tie_vecs])
        # batched call: 2027 rows
        df_min_b, n_b = lat.minimum_image(df)
        # per-vector calls
        n_single = np.empty_like(n_b)
        for i in range(df.shape[0]):
            _, n_single[i] = lat.minimum_image(df[i])
        # 1. batched == per-vector, exactly
        assert np.array_equal(n_b, n_single), (cell, angles)
        # 2. both equal the independent per-row oracle
        n_oracle = np.empty_like(n_b)
        for i in range(df.shape[0]):
            n_oracle[i] = _oracle_row(lat, df[i])
        assert np.array_equal(n_b, n_oracle), (cell, angles)

    # one batch large enough to span several chunks
    lat = Lattice.from_parameters(5.0, 5.0, 5.0, 90.0, 90.0, 20.0)
    df = rng.uniform(-3.0, 3.0, (50_000, 3))
    n_b = lat.minimum_image(df)[1]
    n_single = np.empty_like(n_b)
    for i in range(0, df.shape[0], 500):  # representative sample of per-vector calls
        _, n_single[i] = lat.minimum_image(df[i])
    assert np.array_equal(n_b[::500], n_single[::500])


def test_minimum_image_cubic_half_cell_tie():
    lat = Lattice.from_parameters(5.0, 5.0, 5.0, 90.0, 90.0, 90.0)
    # Tie rule: lexicographically smallest reduced vector -> half-cell ties map to -1/2,
    # i.e. the half-open interval [-1/2, 1/2).
    # (0.5, 0, 0), (-0.5, 0, 0), (1.5, 0, 0) all give (-0.5, 0, 0)
    for df in (np.array([0.5, 0.0, 0.0]), np.array([-0.5, 0.0, 0.0]), np.array([1.5, 0.0, 0.0])):
        df_min, _ = lat.minimum_image(df)
        assert np.allclose(df_min, [-0.5, 0.0, 0.0], atol=1e-12)
    # (0.5, 0.5, 0) gives (-0.5, -0.5, 0)
    df_min, _ = lat.minimum_image(np.array([0.5, 0.5, 0.0]))
    assert np.allclose(df_min, [-0.5, -0.5, 0.0], atol=1e-12)
    # batched
    df_min, n_shift = lat.minimum_image(np.array([[0.5, 0.0, 0.0], [0.0, 0.5, 0.0]]))
    assert n_shift.shape == (2, 3)
    assert np.allclose(df_min[0], [-0.5, 0.0, 0.0], atol=1e-12)
    assert np.allclose(df_min[1], [0.0, -0.5, 0.0], atol=1e-12)


def test_minimum_image_translation_invariance():
    """minimum_image(df + n)[0] == minimum_image(df)[0] for integer n (|n_i| <= 50),
    including constructed exact ties."""
    rng = np.random.default_rng(2024)
    for cell in ((5.0, 5.0, 5.0, 90.0, 90.0, 20.0), (4.0, 4.0, 4.0, 30.0, 30.0, 30.0)):
        lat = Lattice.from_parameters(*cell)
        for _ in range(200):
            df = rng.uniform(-3, 3, 3)
            n = rng.integers(-50, 51, 3)
            base, _ = lat.minimum_image(df)
            shifted, _ = lat.minimum_image(df + n.astype(float))
            assert np.allclose(base, shifted, atol=1e-12)
        # constructed exact ties: half-cell in a cubic-like axis
        for df in (np.array([0.5, 0.0, 0.0]), np.array([0.5, 0.5, 0.5]), np.array([2.5, -1.5, 0.5])):
            base, _ = lat.minimum_image(df)
            shifted, _ = lat.minimum_image(df + np.array([3, -2, 1], float))
            assert np.allclose(base, shifted, atol=1e-12)


def test_minimum_image_large_df_memory():
    """10^5 displacements with |df| up to 1000 complete and agree with df - rint(df)
    reduction; peak memory is recorded via tracemalloc."""
    import tracemalloc
    lat = Lattice.from_parameters(5.0, 5.0, 5.0, 90.0, 90.0, 20.0)
    rng = np.random.default_rng(7)
    n = 100_000
    df = rng.uniform(-1000.0, 1000.0, (n, 3))
    tracemalloc.start()
    df_min, n_shift = lat.minimum_image(df)
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    # The reduced displacement (df - rint(df)) lies in [-0.5, 0.5]^3; the minimum
    # image of df must equal the minimum image of that reduced displacement.
    reduced = df - np.rint(df)
    ref, _ = lat.minimum_image(reduced)
    assert np.allclose(df_min, ref, atol=1e-9)
    # Peak memory should be far below the naive (n * (2*1000)^3) allocation.
    assert peak < 200 * 1024 * 1024  # < 200 MB
    print(f"  peak memory for 10^5 displacements (|df|<=1000): {peak / 1024 / 1024:.1f} MB")


def test_minimum_image_shift_dtype_and_identity():
    """n_shift is integer (regression: chunking once made it float64) and
    df - n_shift equals the returned reduced vector exactly, for (3,), (N, 3)
    and (2, N, 3) inputs."""
    rng = np.random.default_rng(31337)
    for lat, tag in ((Lattice.from_parameters(5.0, 5.0, 5.0, 90.0, 90.0, 90.0), "cubic"),
                     (Lattice.from_parameters(*TRIC), "triclinic"),
                     (Lattice.from_parameters(5.0, 5.0, 5.0, 90.0, 90.0, 20.0), "skewed")):
        for df in (rng.uniform(-2.5, 2.5, 3),
                   rng.uniform(-2.5, 2.5, (500, 3)),
                   rng.uniform(-2.5, 2.5, (2, 500, 3))):
            df_min, n_shift = lat.minimum_image(df)
            assert np.issubdtype(n_shift.dtype, np.integer), (tag, df.shape)
            assert n_shift.shape == df.shape
            assert np.array_equal(df - n_shift, df_min), (tag, df.shape)


def test_distance_pbc():
    lat = Lattice.from_parameters(5.0, 5.0, 5.0, 90.0, 90.0, 90.0)
    # (0,0,0) to (0.9, 0, 0): PBC distance is 0.1 * 5 = 0.5 Å
    assert lat.distance_pbc((0.0, 0.0, 0.0), (0.9, 0.0, 0.0)) == pytest.approx(0.5, rel=1e-12)
    # (0,0,0) to (0.5, 0, 0): 2.5 Å
    assert lat.distance_pbc((0.0, 0.0, 0.0), (0.5, 0.0, 0.0)) == pytest.approx(2.5, rel=1e-12)
    # triclinic: (0,0,0) to (0.1, 0, 0) -> 0.1 * a = 0.41 Å
    lat_t = Lattice.from_parameters(*TRIC)
    assert lat_t.distance_pbc((0.0, 0.0, 0.0), (0.1, 0.0, 0.0)) == pytest.approx(0.41, rel=1e-12)
    # triclinic: (0,0,0) to (0.9, 0, 0) -> PBC wraps to 0.1 * a = 0.41 Å
    assert lat_t.distance_pbc((0.0, 0.0, 0.0), (0.9, 0.0, 0.0)) == pytest.approx(0.41, rel=1e-12)


# ------------------------------ 8. immutability ------------------------------
def test_immutability():
    lat = Lattice.from_parameters(5.0, 5.0, 5.0, 90.0, 90.0, 90.0)
    m = lat.matrix
    m[0, 0] = 999.0  # mutating the returned array must not change the lattice
    assert lat.matrix[0, 0] == 5.0
    # in-place write to the stored array raises (writeable flag cleared)
    with pytest.raises(ValueError):
        lat._A[0, 0] = 9.0
    # rebinding, deletion, and new attributes all raise AttributeError
    with pytest.raises(AttributeError):
        lat._A = np.eye(3)
    with pytest.raises(AttributeError):
        del lat._A
    with pytest.raises(AttributeError):
        lat.new_attr = 1
    # hash is unchanged afterwards
    assert hash(lat) == hash(Lattice.from_parameters(5.0, 5.0, 5.0, 90.0, 90.0, 90.0))
    # from_matrix stores a copy, not the caller's array
    A = np.diag([5.0, 5.0, 5.0])
    lat2 = Lattice.from_matrix(A)
    A[0, 0] = 77.0
    assert lat2.matrix[0, 0] == 5.0


# ------------------------------ 9. misc ------------------------------
def test_repr_and_comparison():
    lat = Lattice.from_parameters(5.0, 5.0, 5.0, 90.0, 90.0, 90.0)
    assert "Lattice" in repr(lat)
    lat2 = Lattice.from_parameters(5.0, 5.0, 5.0, 90.0, 90.0, 90.0)
    assert lat == lat2
    assert lat.isclose(lat2)
    lat3 = Lattice.from_parameters(5.1, 5.0, 5.0, 90.0, 90.0, 90.0)
    assert lat != lat3
    assert not lat.isclose(lat3, rtol=1e-4, atol=0.0)
    assert lat.isclose(lat3, rtol=1e-1, atol=0.0)
    assert lat != "not a lattice"


def test_example_script_runs():
    """CI keeps examples/lattice_basics.py working."""
    import subprocess
    import sys
    import pathlib
    script = pathlib.Path(__file__).parent.parent / "examples" / "lattice_basics.py"
    result = subprocess.run([sys.executable, str(script)], capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    assert "cubic" in result.stdout.lower()
