"""Regular (affine) reciprocal-space grids: volumes, oblique planes and lines (Package B, Feature 6).

Conventions (shared with :class:`~discuspdf.lattice.Lattice` and the verified DISCUS 6.23.01 reader):

* A grid point with integer array index ``i = (i0, i1, i2)`` sits at
  ``hkl = corner + i @ axis_vectors`` (r.l.u.). **Row** ``axis_vectors[k]`` is the hkl step along
  array axis ``k``. DISCUS stores the same steps as **columns**:
  ``DiscusVolume.increments[:, k] == axis_vectors[k]`` (verified on the real skewed C4 fixture).
* Grid points are **voxel centres**. Voxel edges lie at half-integer indices; use
  :meth:`ReciprocalGrid.index_to_hkl` with fractional indices (e.g. ``i - 0.5``).
* Cartesian ``Q = 2π · (hkl @ B)`` in Å⁻¹ (2π included) for the grid's :class:`Lattice` and its
  orientation. ``|Q| = 2π/d``.
* An axis of length 1 is *inactive*: its vector never enters any coordinate. A grid with one
  inactive axis is a plane, with two a line. Rank and conditioning are checked over the
  active axes only.

Nothing here resamples data. Whether an existing transform may be applied to data on this grid is
reported by :meth:`ReciprocalGrid.fft_compatibility`; an oblique grid is never relabelled as
Cartesian Q.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any
import json

import numpy as np

from .lattice import Lattice

SCHEMA_NAME = "discuspdf.ReciprocalGrid"
SCHEMA_VERSION = 1

_DEFAULT_TOL = 1e-9


def _readonly(a: np.ndarray) -> np.ndarray:
    a = np.array(a, dtype=float, copy=True)
    a.flags.writeable = False
    return a


def _as_shape(shape) -> tuple[int, int, int]:
    try:
        items = tuple(shape)
    except TypeError as exc:
        raise ValueError("shape must be three positive integers") from exc
    if len(items) != 3:
        raise ValueError(f"shape must have three entries, got {len(items)}")
    out = []
    for n in items:
        if isinstance(n, (bool, np.bool_)) or not isinstance(n, (int, np.integer)):
            raise ValueError(f"shape entries must be integers, got {n!r}")
        if n < 1:
            raise ValueError(f"shape entries must be >= 1, got {n}")
        out.append(int(n))
    return tuple(out)  # type: ignore[return-value]


def _as_vec3(name: str, value) -> np.ndarray:
    a = np.asarray(value, dtype=float)
    if a.shape != (3,) or not np.isfinite(a).all():
        raise ValueError(f"{name} must be three finite numbers")
    return a


@dataclass(frozen=True)
class FFTCompatibility:
    """Which existing 3D transform may be applied to data sampled on a grid, and why not otherwise.

    * ``fractional`` — :func:`discuspdf.delta_pdf_fractional` (hkl axis-aligned, increasing,
      centred on hkl = 0 at index n//2). Valid for any cell metric.
    * ``cartesian`` — :func:`discuspdf.delta_pdf_from_diffuse` / :class:`ReciprocalVolume`
      (Q step k parallel to +Cartesian axis k, centred on Q = 0 at index n//2).
    """
    fractional: bool
    cartesian: bool
    reasons: tuple[str, ...]

    def __bool__(self) -> bool:  # pragma: no cover - guard against `if grid.fft_compatibility():`
        raise TypeError("FFTCompatibility is not a truth value; test .fractional or .cartesian")


class ReciprocalGrid:
    """Immutable regular grid in reciprocal space, tied to a :class:`Lattice`.

    Parameters
    ----------
    lattice:
        Direct lattice (with its Cartesian orientation) that defines hkl -> Q.
    corner:
        hkl (r.l.u.) of array element ``[0, 0, 0]``.
    axis_vectors:
        3x3; row ``k`` is the hkl step between neighbouring points along array axis ``k``.
    shape:
        Number of points ``(n0, n1, n2)``, each >= 1.
    """

    __slots__ = ("_lattice", "_corner", "_axis_vectors", "_shape")

    def __init__(self, lattice: Lattice, corner, axis_vectors, shape) -> None:
        if not isinstance(lattice, Lattice):
            raise TypeError("lattice must be a discuspdf.Lattice")
        shape = _as_shape(shape)
        corner = _as_vec3("corner", corner)
        V = np.asarray(axis_vectors, dtype=float)
        if V.shape != (3, 3) or not np.isfinite(V).all():
            raise ValueError("axis_vectors must be a finite 3x3 array (row k = hkl step of axis k)")
        active = [k for k in range(3) if shape[k] > 1]
        if active:
            Va = V[active]
            norms = np.linalg.norm(Va, axis=1)
            for k, nrm in zip(active, norms):
                if nrm == 0.0:
                    raise ValueError(f"axis {k} has {shape[k]} points but a zero step vector")
            sv = np.linalg.svd(Va, compute_uv=False)
            if sv[-1] <= 1e-10 * sv[0]:
                raise ValueError("active axis vectors are linearly dependent or nearly so "
                                 f"(singular values {sv.tolist()}); the grid would fold onto itself")
        object.__setattr__(self, "_lattice", lattice)
        object.__setattr__(self, "_corner", _readonly(corner))
        object.__setattr__(self, "_axis_vectors", _readonly(V))
        object.__setattr__(self, "_shape", shape)

    def __setattr__(self, name, value) -> None:
        raise AttributeError(f"ReciprocalGrid is immutable; cannot set {name!r}")

    def __delattr__(self, name) -> None:
        raise AttributeError(f"ReciprocalGrid is immutable; cannot delete {name!r}")

    # -- construction -----------------------------------------------------
    @classmethod
    def from_corners(cls, lattice: Lattice, ll, lr, ul, tl, shape) -> "ReciprocalGrid":
        """DISCUS-style definition: corners ``ll`` (element [0,0,0]), ``lr`` (last point of axis 0),
        ``ul`` (last point of axis 1), ``tl`` (last point of axis 2), and point counts
        ``shape = (na, no, nt)``. Mirrors the DISCUS ``fourier`` menu commands of the same names.

        For an axis with one point the matching corner must equal ``ll`` (DISCUS behaviour for
        such corners is not verified, so a contradiction is refused rather than guessed).
        """
        shape = _as_shape(shape)
        ll = _as_vec3("ll", ll)
        ends = [_as_vec3(name, v) for name, v in (("lr", lr), ("ul", ul), ("tl", tl))]
        V = np.zeros((3, 3))
        for k, (name, end) in enumerate(zip(("lr", "ul", "tl"), ends)):
            if shape[k] > 1:
                V[k] = (end - ll) / (shape[k] - 1)
            elif not np.allclose(end, ll, rtol=0, atol=1e-12):
                raise ValueError(f"axis {k} has a single point but {name}={end.tolist()} != ll={ll.tolist()}")
        return cls(lattice, ll, V, shape)

    @classmethod
    def from_ranges(cls, lattice: Lattice, h, k, l) -> "ReciprocalGrid":
        """hkl-axis-aligned grid from inclusive ranges ``(start, stop, n)`` per axis (DISCUS-like:
        both end points are sampled). ``n = 1`` requires ``start == stop``."""
        corner, steps, shape = [], [], []
        for name, spec in (("h", h), ("k", k), ("l", l)):
            try:
                start, stop, n = spec
            except (TypeError, ValueError) as exc:
                raise ValueError(f"{name} range must be (start, stop, n)") from exc
            n = _as_shape((n, 1, 1))[0]
            start, stop = float(start), float(stop)
            if not (np.isfinite(start) and np.isfinite(stop)):
                raise ValueError(f"{name} range must be finite")
            if n == 1:
                if start != stop:
                    raise ValueError(f"{name}: a single point needs start == stop")
                steps.append(0.0)
            else:
                if stop == start:
                    raise ValueError(f"{name}: start == stop with {n} points")
                steps.append((stop - start) / (n - 1))
            corner.append(start)
            shape.append(n)
        return cls(lattice, corner, np.diag(steps), shape)

    @classmethod
    def centered(cls, lattice: Lattice, shape, steps) -> "ReciprocalGrid":
        """hkl-axis-aligned grid with hkl = 0 at index ``n // 2`` on every axis (FFT-ready).

        ``steps`` are positive hkl steps (r.l.u.). Odd ``n`` is symmetric about 0; even ``n``
        has one more point on the negative side (numpy ``fftshift`` convention)."""
        shape = _as_shape(shape)
        steps = _as_vec3("steps", steps)
        if np.any(steps <= 0):
            raise ValueError("steps must be positive")
        corner = -np.array([n // 2 for n in shape], dtype=float) * steps
        return cls(lattice, corner, np.diag(steps), shape)

    @classmethod
    def from_discus_volume(cls, volume) -> "ReciprocalGrid":
        """Grid of a :class:`~discuspdf.discus_output.DiscusVolume` (verified DISCUS 6.23.01 NeXus).

        The lattice is rebuilt from the file's cell parameters in the canonical orientation
        (a ∥ x, b in xy); DISCUS files carry no orientation matrix."""
        lengths = np.asarray(volume.cell_lengths, dtype=float)
        angles = np.asarray(volume.cell_angles, dtype=float)
        lattice = Lattice.from_parameters(*lengths, *angles)
        inc = np.asarray(volume.increments, dtype=float)
        return cls(lattice, volume.corner, inc.T, np.shape(volume.values))

    # -- basic properties --------------------------------------------------
    @property
    def lattice(self) -> Lattice:
        return self._lattice

    @property
    def corner(self) -> np.ndarray:
        """hkl of element [0,0,0]. Returns a copy."""
        return self._corner.copy()

    @property
    def axis_vectors(self) -> np.ndarray:
        """Row k = hkl step along array axis k. Returns a copy."""
        return self._axis_vectors.copy()

    @property
    def shape(self) -> tuple[int, int, int]:
        return self._shape

    @property
    def size(self) -> int:
        return int(np.prod(self._shape))

    @property
    def active_axes(self) -> tuple[int, ...]:
        """Array axes with more than one point."""
        return tuple(k for k in range(3) if self._shape[k] > 1)

    @property
    def ndim(self) -> int:
        """Effective dimensionality: 3 volume, 2 plane, 1 line, 0 single point."""
        return len(self.active_axes)

    def _effective_vectors(self) -> np.ndarray:
        V = self._axis_vectors.copy()
        for k in range(3):
            if self._shape[k] == 1:
                V[k] = 0.0
        return V

    # -- coordinates -------------------------------------------------------
    def index_to_hkl(self, index) -> np.ndarray:
        """hkl (r.l.u.) at array index (shape ``(3,)`` or ``(..., 3)``); fractional indices allowed,
        so ``i ± 0.5`` gives voxel edges. Inactive axes contribute nothing."""
        idx = np.asarray(index, dtype=float)
        if idx.shape[-1:] != (3,):
            raise ValueError(f"index must have last dimension 3, got shape {idx.shape}")
        if not np.isfinite(idx).all():
            raise ValueError("index must be finite")
        return self._corner + idx @ self._effective_vectors()

    def hkl_to_index(self, hkl) -> np.ndarray:
        """Fractional array index of hkl. Only for volumes (3 active axes); a plane or line has no
        inverse for off-grid points."""
        if self.ndim != 3:
            raise ValueError(f"hkl_to_index needs a 3D grid; this grid has {self.ndim} active axes")
        hkl = np.asarray(hkl, dtype=float)
        if hkl.shape[-1:] != (3,):
            raise ValueError(f"hkl must have last dimension 3, got shape {hkl.shape}")
        if not np.isfinite(hkl).all():
            raise ValueError("hkl must be finite")
        return np.linalg.solve(self._axis_vectors.T, (hkl - self._corner).reshape(-1, 3).T).T.reshape(hkl.shape)

    def hkl(self) -> np.ndarray:
        """hkl of every point, shape ``shape + (3,)``."""
        grids = np.meshgrid(*[np.arange(n, dtype=float) for n in self._shape], indexing="ij")
        return self.index_to_hkl(np.stack(grids, axis=-1))

    def q(self) -> np.ndarray:
        """Cartesian Q (Å⁻¹, 2π included) of every point, shape ``shape + (3,)``."""
        return self._lattice.hkl_to_q(self.hkl())

    def q_magnitude(self) -> np.ndarray:
        """|Q| (Å⁻¹) of every point; equals 2π/d at a reflection."""
        return np.linalg.norm(self.q(), axis=-1)

    def q_corner(self) -> np.ndarray:
        """Cartesian Q of element [0,0,0]."""
        return self._lattice.hkl_to_q(self._corner)

    def q_step_vectors(self) -> np.ndarray:
        """Row k = Cartesian Q step (Å⁻¹) along array axis k (zero for inactive axes)."""
        return 2.0 * np.pi * (self._effective_vectors() @ self._lattice.reciprocal_matrix)

    def hkl_axes(self, tol: float = _DEFAULT_TOL) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Independent 1D hkl axes ``(h, k, l)``; only for an hkl-axis-aligned grid."""
        if not self.is_hkl_axis_aligned(tol):
            raise ValueError("grid is skewed in hkl: axis k is not parallel to hkl direction k, so "
                             "independent 1D h, k, l axes do not exist")
        steps = np.diag(self._effective_vectors())
        return tuple(self._corner[k] + steps[k] * np.arange(n) for k, n in enumerate(self._shape))

    def q_axes(self, tol: float = _DEFAULT_TOL) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Independent 1D Cartesian Q axes (Å⁻¹); only when step k is parallel to Cartesian axis k."""
        if not self.is_q_axis_aligned(tol):
            raise ValueError("Q step vectors are not parallel to the Cartesian x, y, z axes; exact 1D "
                             "Cartesian Q axes do not exist for this grid and lattice orientation")
        S = self.q_step_vectors()
        q0 = self.q_corner()
        return tuple(q0[k] + S[k, k] * np.arange(n) for k, n in enumerate(self._shape))

    # -- geometry tests -----------------------------------------------------
    @staticmethod
    def _parallel_to_unit(vectors: np.ndarray, active, tol: float, require_positive: bool) -> bool:
        for k in active:
            v = vectors[k]
            off = np.delete(np.abs(v), k)
            if np.any(off > tol * np.linalg.norm(v)):
                return False
            if require_positive and v[k] <= 0:
                return False
        return True

    def is_hkl_axis_aligned(self, tol: float = _DEFAULT_TOL) -> bool:
        """Active axis k runs along hkl direction k (h, k or l) only, in either sense."""
        return self._parallel_to_unit(self._axis_vectors, self.active_axes, tol, False)

    def is_q_axis_aligned(self, tol: float = _DEFAULT_TOL) -> bool:
        """Active axis k runs along +Cartesian axis k (x, y or z) in this lattice orientation."""
        return self._parallel_to_unit(self.q_step_vectors(), self.active_axes, tol, True)

    def is_q_orthogonal(self, tol: float = _DEFAULT_TOL) -> bool:
        """Active Q step vectors are mutually perpendicular (|cos| <= tol)."""
        S = self.q_step_vectors()
        act = self.active_axes
        for a in range(len(act)):
            for b in range(a + 1, len(act)):
                u, v = S[act[a]], S[act[b]]
                if abs(u @ v) > tol * np.linalg.norm(u) * np.linalg.norm(v):
                    return False
        return True

    def handedness(self) -> int:
        """+1 if the array axes (0, 1, 2) form a right-handed frame in Cartesian Q, -1 if left-handed.
        Only defined for volumes."""
        if self.ndim != 3:
            raise ValueError("handedness is defined only for a 3D grid")
        return 1 if np.linalg.det(self.q_step_vectors()) > 0 else -1

    def is_centered(self, tol: float = _DEFAULT_TOL) -> bool:
        """hkl = 0 at index ``n // 2`` along every axis (inactive axes: the corner component is 0)."""
        centre = self.index_to_hkl([n // 2 for n in self._shape])
        scale = max(float(np.max(np.abs(self._axis_vectors))), 1.0)
        return bool(np.all(np.abs(centre) <= tol * scale))

    def fft_compatibility(self, tol: float = _DEFAULT_TOL) -> FFTCompatibility:
        """Report which existing 3D transform is exact for data on this grid (no resampling)."""
        reasons: list[str] = []
        if self.ndim != 3:
            reasons.append(f"grid has {self.ndim} active axes; the 3D transforms need a volume")
            return FFTCompatibility(False, False, tuple(reasons))
        centred = self.is_centered(tol)
        if not centred:
            reasons.append("hkl = 0 is not at index n//2 on every axis (centred FFT would be mis-phased)")
        frac = True
        if not self.is_hkl_axis_aligned(tol):
            frac = False
            reasons.append("fractional: grid is skewed in hkl (non-diagonal axis_vectors)")
        elif np.any(np.diag(self._axis_vectors) <= 0):
            frac = False
            reasons.append("fractional: hkl steps must be positive (increasing axes)")
        cart = True
        if not self.is_q_axis_aligned(tol):
            cart = False
            reasons.append("cartesian: Q steps are not parallel to +x, +y, +z for this lattice orientation "
                           "(use the fractional transform for non-orthogonal cells)")
        return FFTCompatibility(frac and centred, cart and centred, tuple(reasons))

    def real_space_step_vectors(self) -> np.ndarray:
        """Fractional real-space step vectors (rows) of the grid that a DFT over this volume returns.

        Row k satisfies ``axis_vectors @ row_k = e_k / n_k``; for an hkl-aligned grid this is
        ``1 / (n_k · dh_k)`` along direction k, as in :func:`delta_pdf_fractional`. Lengths in Å
        follow from ``lattice.frac_to_cart``."""
        if self.ndim != 3:
            raise ValueError("real-space dual steps are defined only for a 3D grid")
        R = np.linalg.inv(self._axis_vectors).T
        return R / np.asarray(self._shape, dtype=float)[:, None]

    # -- derived grids -------------------------------------------------------
    def plane(self, axis: int, index: int) -> "ReciprocalGrid":
        """The single-layer grid at ``index`` along ``axis`` (that axis becomes inactive)."""
        if axis not in (0, 1, 2):
            raise ValueError("axis must be 0, 1 or 2")
        if isinstance(index, bool) or not isinstance(index, (int, np.integer)) or not 0 <= index < self._shape[axis]:
            raise ValueError(f"index must be an integer in [0, {self._shape[axis]})")
        corner = self._corner + int(index) * self._effective_vectors()[axis]
        shape = list(self._shape)
        shape[axis] = 1
        return ReciprocalGrid(self._lattice, corner, self._axis_vectors, shape)

    def discus_corners(self) -> dict[str, Any]:
        """``ll, lr, ul, tl`` (hkl lists) and ``na, no, nt`` as for the DISCUS ``fourier`` menu.
        Inverse of :meth:`from_corners`."""
        last = [self.index_to_hkl([n - 1 if j == k else 0 for j, n in enumerate(self._shape)])
                for k in range(3)]
        return {"ll": self._corner.tolist(), "lr": last[0].tolist(), "ul": last[1].tolist(),
                "tl": last[2].tolist(), "na": self._shape[0], "no": self._shape[1], "nt": self._shape[2]}

    # -- persistence -----------------------------------------------------------
    def to_dict(self) -> dict[str, Any]:
        """JSON-safe description, including the full lattice matrix (orientation preserved)."""
        return {"schema": SCHEMA_NAME, "version": SCHEMA_VERSION,
                "lattice_matrix": self._lattice.matrix.tolist(),
                "corner": self._corner.tolist(),
                "axis_vectors": self._axis_vectors.tolist(),
                "shape": list(self._shape),
                "units": {"corner": "r.l.u.", "axis_vectors": "r.l.u. per index",
                          "lattice_matrix": "angstrom, rows a b c"}}

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "ReciprocalGrid":
        """Inverse of :meth:`to_dict`; refuses unknown schemas and versions."""
        if not isinstance(data, dict):
            raise TypeError("from_dict expects a dict")
        if data.get("schema") != SCHEMA_NAME:
            raise ValueError(f"not a {SCHEMA_NAME} description (schema={data.get('schema')!r})")
        if data.get("version") != SCHEMA_VERSION:
            raise ValueError(f"unsupported {SCHEMA_NAME} version {data.get('version')!r}; "
                             f"this code reads version {SCHEMA_VERSION}")
        missing = [k for k in ("lattice_matrix", "corner", "axis_vectors", "shape") if k not in data]
        if missing:
            raise ValueError(f"{SCHEMA_NAME} description is missing {missing}")
        return cls(Lattice.from_matrix(data["lattice_matrix"]), data["corner"],
                   data["axis_vectors"], tuple(data["shape"]))

    def to_json(self) -> str:
        return json.dumps(self.to_dict(), sort_keys=True)

    @classmethod
    def from_json(cls, text: str) -> "ReciprocalGrid":
        return cls.from_dict(json.loads(text))

    # -- comparison ------------------------------------------------------------
    def __eq__(self, other) -> bool:
        if not isinstance(other, ReciprocalGrid):
            return NotImplemented
        return (self._shape == other._shape and self._lattice == other._lattice
                and np.array_equal(self._corner, other._corner)
                and np.array_equal(self._effective_vectors(), other._effective_vectors()))

    def __hash__(self) -> int:
        return hash((self._shape, hash(self._lattice), self._corner.tobytes(),
                     self._effective_vectors().tobytes()))

    def isclose(self, other: "ReciprocalGrid", rtol: float = 1e-12, atol: float = 1e-12) -> bool:
        """Same shape, and lattice, corner and active axis vectors equal within tolerance."""
        if not isinstance(other, ReciprocalGrid):
            raise TypeError("isclose requires another ReciprocalGrid")
        return (self._shape == other._shape
                and self._lattice.isclose(other._lattice, rtol=rtol, atol=atol)
                and np.allclose(self._corner, other._corner, rtol=rtol, atol=atol)
                and np.allclose(self._effective_vectors(), other._effective_vectors(), rtol=rtol, atol=atol))

    def __repr__(self) -> str:
        kind = {3: "volume", 2: "plane", 1: "line", 0: "point"}[self.ndim]
        return (f"ReciprocalGrid({kind}, shape={self._shape}, corner={self._corner.tolist()}, "
                f"axis_vectors={self._axis_vectors.tolist()}, {self._lattice!r})")
