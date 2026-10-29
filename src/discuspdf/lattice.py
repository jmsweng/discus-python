"""Crystal-lattice value type: direct/reciprocal basis, Q vectors, d-spacings, PBC distances.

Conventions (shared with the validated ``DiscusVolume`` geometry in ``discus_output.py``):

* Direct basis ``A``: 3x3, **rows** are **a**, **b**, **c** in Å.
* Fractional -> Cartesian: row vectors, ``r = f @ A``.
* Reciprocal basis ``B = inv(A).T``: rows **a\\***, **b\\***, **c\\*** in Å⁻¹, **no 2π**
  (``A @ B.T = I``).
* Cartesian Q: ``Q = 2π · (hkl @ B)`` in Å⁻¹ (2π included); hkl in r.l.u.
* d-spacing: ``d_hkl = 1 / |hkl @ B|`` in Å; ``hkl = (0,0,0)`` raises.
* Metric ``G = A @ A.T``; reciprocal metric ``G* = inv(G) = B @ B.T``.
* Volume ``det(A)`` in Å³, must be > 0.
* Angles in degrees in and out; lengths in Å.

``Lattice.from_parameters`` constructs the canonical orientation (**a** ∥ x, **b** in the xy
plane, right-handed) and uses the same rejection rules as ``discus_output._validated_direct_basis``.
It does **not** reconstruct an experimental orientation matrix: any rotation of the canonical
cell is a valid ``Lattice`` and must be supplied via :meth:`Lattice.from_matrix`.
"""
from __future__ import annotations

import numpy as np


def _validate_lengths(lengths: np.ndarray) -> None:
    if lengths.shape != (3,) or not np.isfinite(lengths).all() or np.any(lengths <= 0):
        raise ValueError("unit cell lengths must be three finite positive values in Å")


def _validate_angles(angles: np.ndarray) -> None:
    if angles.shape != (3,) or not np.isfinite(angles).all() or np.any(angles <= 0) or np.any(angles >= 180):
        raise ValueError("unit cell angles must be finite degrees strictly between 0 and 180")


def _canonical_basis(a: float, b: float, c: float, alpha: float, beta: float, gamma: float) -> np.ndarray:
    """Canonical orientation (a ∥ x, b in xy plane, right-handed); same rejection rules as
    ``discus_output._validated_direct_basis``."""
    alpha_r, beta_r, gamma_r = np.deg2rad((alpha, beta, gamma))
    ca, cb, cg = np.cos((alpha_r, beta_r, gamma_r))
    metric = np.array([[1.0, cg, cb], [cg, 1.0, ca], [cb, ca, 1.0]], dtype=float)
    eigenvalues = np.linalg.eigvalsh(metric)
    if eigenvalues[0] <= 1e-12 * eigenvalues[-1]:
        raise ValueError("Invalid or degenerate unit cell: metric is not positive definite")
    sg = np.sin(gamma_r)
    bx, by = b * cg, b * sg
    cx = c * cb
    cy = c * (ca - cb * cg) / sg
    cz_sq = c * c - cx * cx - cy * cy
    if not np.isfinite(cz_sq) or cz_sq <= 0:
        raise ValueError("Invalid unit cell: nonpositive squared c_z component")
    basis = np.array([[a, 0.0, 0.0], [bx, by, 0.0], [cx, cy, np.sqrt(cz_sq)]], dtype=float)
    if not np.isfinite(basis).all() or np.linalg.det(basis) <= 0:
        raise ValueError("Invalid unit cell basis")
    return basis


def _validate_basis(A: np.ndarray) -> None:
    """Reject singular, left-handed, or non-finite 3x3 basis matrices."""
    if A.shape != (3, 3) or not np.issubdtype(A.dtype, np.number):
        raise ValueError("lattice matrix must be a numeric 3x3 array")
    A = np.asarray(A, dtype=float)
    if not np.isfinite(A).all():
        raise ValueError("lattice matrix must be finite")
    det = np.linalg.det(A)
    if det <= 0:
        raise ValueError(
            f"lattice matrix is singular (det = 0) or left-handed (det < 0); got det = {det:.6g}. "
            "Only right-handed cells (det > 0) are supported."
        )
    if np.linalg.cond(A) > 1e12:
        raise ValueError("lattice matrix is nearly singular")


class Lattice:
    """Immutable crystal lattice.

    The 3x3 direct basis ``A`` (rows **a**, **b**, **c** in Å) is stored read-only.
    Construct with :meth:`from_parameters` (canonical orientation) or :meth:`from_matrix`
    (any right-handed orientation).
    """

    __slots__ = ("_A", "_B", "_Ainv")

    def __init__(self, A: np.ndarray) -> None:
        A = np.array(A, dtype=float, copy=True)
        _validate_basis(A)
        A.flags.writeable = False
        B = np.linalg.solve(A.T, np.eye(3))
        B.flags.writeable = False
        Ainv = np.linalg.inv(A)
        Ainv.flags.writeable = False
        object.__setattr__(self, "_A", A)
        object.__setattr__(self, "_B", B)
        object.__setattr__(self, "_Ainv", Ainv)

    def __setattr__(self, name, value) -> None:
        raise AttributeError(f"Lattice is immutable; cannot set {name!r}")

    def __delattr__(self, name) -> None:
        raise AttributeError(f"Lattice is immutable; cannot delete {name!r}")

    # -- construction -----------------------------------------------------
    @classmethod
    def from_parameters(cls, a: float, b: float, c: float, alpha: float, beta: float, gamma: float) -> "Lattice":
        """Build a lattice from six crystallographic parameters.

        ``a, b, c`` in Å; ``alpha, beta, gamma`` in degrees (alpha = angle between
        b and c; beta = angle between a and c; gamma = angle between a and b).
        The result uses the canonical orientation (**a** ∥ x, **b** in the xy plane,
        right-handed); an experimentally rotated orientation is not reconstructed —
        use :meth:`from_matrix` for that.
        """
        lengths = np.asarray((a, b, c), dtype=float)
        angles = np.asarray((alpha, beta, gamma), dtype=float)
        _validate_lengths(lengths)
        _validate_angles(angles)
        return cls(_canonical_basis(a, b, c, alpha, beta, gamma))

    @classmethod
    def from_matrix(cls, A) -> "Lattice":
        """Build a lattice from a 3x3 direct basis ``A`` (rows **a**, **b**, **c** in Å).

        Any right-handed orientation is preserved. ``det(A) <= 0`` (singular or
        left-handed) raises ``ValueError``; there is no left-handed opt-in.
        """
        A = np.asarray(A, dtype=float)
        _validate_basis(A)
        return cls(A)

    # -- properties -------------------------------------------------------
    @property
    def matrix(self) -> np.ndarray:
        """Direct basis ``A`` (rows a, b, c in Å). Returns a copy; the stored array is read-only."""
        return self._A.copy()

    @property
    def reciprocal_matrix(self) -> np.ndarray:
        """Reciprocal basis ``B = inv(A).T`` (rows a*, b*, c* in Å⁻¹, no 2π). Returns a copy."""
        return self._B.copy()

    @property
    def parameters(self) -> tuple[float, float, float, float, float, float]:
        """(a, b, c, alpha, beta, gamma) in Å and degrees, recovered from the stored basis.

        alpha = angle between b and c; beta = angle between a and c; gamma = angle between a and b.
        """
        A = self._A
        a, b, c = np.linalg.norm(A, axis=1)
        cos_ac = (A[0] @ A[2]) / (a * c)   # beta
        cos_bc = (A[1] @ A[2]) / (b * c)   # alpha
        cos_ab = (A[0] @ A[1]) / (a * b)   # gamma
        # clamp to [-1, 1] to absorb round-off before arccos
        alpha = float(np.rad2deg(np.arccos(np.clip(cos_bc, -1.0, 1.0))))
        beta = float(np.rad2deg(np.arccos(np.clip(cos_ac, -1.0, 1.0))))
        gamma = float(np.rad2deg(np.arccos(np.clip(cos_ab, -1.0, 1.0))))
        return (float(a), float(b), float(c), alpha, beta, gamma)

    @property
    def volume(self) -> float:
        """Cell volume in Å³ (det A, positive by construction)."""
        return float(np.linalg.det(self._A))

    @property
    def metric(self) -> np.ndarray:
        """Direct metric tensor ``G = A @ A.T``. Returns a copy."""
        return self._A @ self._A.T

    @property
    def reciprocal_metric(self) -> np.ndarray:
        """Reciprocal metric tensor ``G* = inv(G) = B @ B.T``. Returns a copy."""
        return self._B @ self._B.T

    # -- transforms -------------------------------------------------------
    def frac_to_cart(self, f) -> np.ndarray:
        """Fractional -> Cartesian (Å): ``r = f @ A``.

        ``f`` has shape ``(3,)`` or ``(..., 3)``; the result keeps the leading shape.
        """
        f = np.asarray(f, dtype=float)
        if f.shape[-1:] != (3,):
            raise ValueError(f"fractional coordinates must have last dimension 3, got shape {f.shape}")
        if not np.isfinite(f).all():
            raise ValueError("fractional coordinates must be finite")
        return f @ self._A

    def cart_to_frac(self, r) -> np.ndarray:
        """Cartesian (Å) -> fractional: ``f = r @ inv(A)``.

        ``r`` has shape ``(3,)`` or ``(..., 3)``; the result keeps the leading shape.
        """
        r = np.asarray(r, dtype=float)
        if r.shape[-1:] != (3,):
            raise ValueError(f"Cartesian coordinates must have last dimension 3, got shape {r.shape}")
        if not np.isfinite(r).all():
            raise ValueError("Cartesian coordinates must be finite")
        return r @ self._Ainv

    def hkl_to_q(self, hkl) -> np.ndarray:
        """hkl (r.l.u.) -> Cartesian Q in Å⁻¹: ``Q = 2π · (hkl @ B)`` (2π included)."""
        hkl = np.asarray(hkl, dtype=float)
        if hkl.shape[-1:] != (3,):
            raise ValueError(f"hkl must have last dimension 3, got shape {hkl.shape}")
        if not np.isfinite(hkl).all():
            raise ValueError("hkl must be finite")
        return 2.0 * np.pi * (hkl @ self._B)

    def q_to_hkl(self, q) -> np.ndarray:
        """Cartesian Q (Å⁻¹, 2π convention) -> hkl in r.l.u.: ``hkl = (Q / 2π) @ B⁻¹ = (Q / 2π) @ Aᵀ``."""
        q = np.asarray(q, dtype=float)
        if q.shape[-1:] != (3,):
            raise ValueError(f"Q must have last dimension 3, got shape {q.shape}")
        if not np.isfinite(q).all():
            raise ValueError("Q must be finite")
        return (q / (2.0 * np.pi)) @ self._A.T

    def d_spacing(self, hkl) -> float:
        """d-spacing in Å: ``d = 1 / |hkl @ B|``.

        ``hkl = (0,0,0)`` raises. Hkl vectors with ``|hkl| < 1e-8`` are treated
        as zero and also raise (``np.allclose`` threshold).
        """
        hkl = np.asarray(hkl, dtype=float)
        if hkl.shape != (3,):
            raise ValueError("d_spacing expects a single hkl triple of shape (3,)")
        if not np.isfinite(hkl).all():
            raise ValueError("hkl must be finite")
        if np.allclose(hkl, 0.0):
            raise ValueError("d_spacing is undefined for hkl = (0, 0, 0)")
        return 1.0 / float(np.linalg.norm(hkl @ self._B))

    # -- periodic boundary conditions -------------------------------------
    def minimum_image(self, df) -> tuple[np.ndarray, np.ndarray]:
        """Minimum-image convention for fractional displacements.

        ``df`` has shape ``(3,)`` or ``(..., 3)`` (fractional). Returns
        ``(df_min, n_shift)`` where ``df_min = df - n_shift`` is the shortest
        periodic image and ``n_shift`` (same leading shape, integers) is the
        shift applied. The result is a function of ``df mod 1`` (translation
        invariant: ``minimum_image(df + n)[0] == minimum_image(df)[0]`` for
        integer ``n``).

        Algorithm (exact bounded search; naive per-component rounding is **not**
        correct in skewed cells):

        1. ``n0 = round(df)``, ``r0 = (df - n0) @ A``, ``R = |r0|``.
        2. Any better image satisfies ``|r| <= R``. Since the fractional
           components of a shift ``n`` are ``f_i = r . b_i`` (rows of ``B``),
           each component of the *offset* ``k = n - n0`` lies in
           ``[-R |b_i|, R |b_i|]``. Enumerate integer offsets in that range per
           axis (vectorised, in chunks so peak memory is independent of the
           number of displacements) and keep the shortest image.
        3. **Ties:** among images within ``atol = 1e-12·R + 1e-15`` Å of the
           shortest, return the one whose **reduced fractional vector
           ``df − n`` is lexicographically smallest** (component 0, then 1,
           then 2; components equal if they differ by < 1e-12). In a cubic
           cell this maps half-cell ties to −½, i.e. the half-open interval
           [−½, ½). The result is deterministic and translation invariant.
        """
        df = np.asarray(df, dtype=float)
        if df.shape[-1:] != (3,):
            raise ValueError(f"fractional displacements must have last dimension 3, got shape {df.shape}")
        if not np.isfinite(df).all():
            raise ValueError("fractional displacements must be finite")
        leading = df.shape[:-1]
        flat = df.reshape(-1, 3)
        A = self._A
        B = self._B

        n0 = np.rint(flat)
        r0 = (flat - n0) @ A
        R = np.linalg.norm(r0, axis=1)
        Rmax = R.max() if R.size else 0.0
        # Bounded search window for the offset k = n - n0: |k_i| <= R |b_i| (rows of B).
        b_norms = np.linalg.norm(B, axis=1)  # (3,)
        bounds = np.ceil(Rmax * b_norms).astype(int) + 1
        # Candidate offsets k in the box [-bounds, +bounds] per axis.
        ranges = [np.arange(-bounds[i], bounds[i] + 1) for i in range(3)]
        grid = np.stack(np.meshgrid(*ranges, indexing="ij"), axis=-1)  # (N1, N2, N3, 3)
        k_all = grid.reshape(-1, 3).astype(float)  # (K, 3)

        # Process in chunks so the temporary (chunk, K, 3) stays bounded.
        chunk = max(1, min(flat.shape[0], 65536 // max(k_all.shape[0], 1)))
        out = np.empty_like(flat)
        n_out = np.empty(flat.shape, dtype=np.int64)
        for start in range(0, flat.shape[0], chunk):
            stop = min(start + chunk, flat.shape[0])
            fc = flat[start:stop]
            nc = n0[start:stop]
            Rc = R[start:stop]
            # Cartesian displacement for each candidate: r = (df - (n0 + k)) @ A
            d = fc[:, None, :] - nc[:, None, :] - k_all[None, :, :]
            dist = np.linalg.norm(d @ A, axis=2)  # (chunk, K)
            dmin = dist.min(axis=1)
            best = dist.argmin(axis=1)
            # Tie rule: among candidates within atol of the shortest, pick the one
            # whose reduced fractional vector (df - n) is lexicographically smallest.
            atol = 1e-12 * Rc + 1e-15
            tie = dist <= dmin[:, None] + atol[:, None]  # (chunk, K)
            if tie.any():
                # Reduced fractional vector for each candidate: (df - (n0 + k))
                red = fc[:, None, :] - nc[:, None, :] - k_all[None, :, :]  # (chunk, K, 3)
                # Lexicographic key with 1e-12 component equality:
                # key = (round(red_0 / 1e-12)) * S^2 + (round(red_1 / 1e-12)) * S + round(red_2 / 1e-12)
                S = 2 * int(bounds.max()) + 3  # > max possible |red_i / 1e-12| in the window
                key = (np.round(red[:, :, 0] / 1e-12) * S**2
                       + np.round(red[:, :, 1] / 1e-12) * S
                       + np.round(red[:, :, 2] / 1e-12))
                key[~tie] = np.inf
                best = key.argmin(axis=1)
            n_best = (nc + k_all[best]).astype(int)
            out[start:stop] = fc - n_best
            n_out[start:stop] = n_best
        return out.reshape(leading + (3,)), n_out.reshape(leading + (3,))

    def distance_pbc(self, f1, f2) -> float:
        """Distance in Å between two fractional points under the minimum-image convention."""
        f1 = np.asarray(f1, dtype=float)
        f2 = np.asarray(f2, dtype=float)
        if f1.shape != (3,) or f2.shape != (3,):
            raise ValueError("distance_pbc expects two fractional points of shape (3,)")
        df_min, _ = self.minimum_image(f1 - f2)
        return float(np.linalg.norm(df_min @ self._A))

    # -- comparison -------------------------------------------------------
    def __eq__(self, other) -> bool:
        if not isinstance(other, Lattice):
            return False
        return np.array_equal(self._A, other._A)

    def isclose(self, other, rtol: float = 1e-12, atol: float = 1e-15) -> bool:
        """Element-wise closeness of the stored basis matrices."""
        if not isinstance(other, Lattice):
            raise TypeError("isclose requires another Lattice")
        return bool(np.allclose(self._A, other._A, rtol=rtol, atol=atol))

    def __hash__(self) -> int:
        return hash(self._A.tobytes())

    def __repr__(self) -> str:
        a, b, c, alpha, beta, gamma = self.parameters
        return (f"Lattice(a={a:.6g} Å, b={b:.6g} Å, c={c:.6g} Å, "
                f"alpha={alpha:.4g}°, beta={beta:.4g}°, gamma={gamma:.4g}°, "
                f"V={self.volume:.6g} Å³)")
