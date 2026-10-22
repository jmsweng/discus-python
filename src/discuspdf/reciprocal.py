"""Typed reciprocal-space volume with explicit physical units and provenance."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any
import numpy as np


@dataclass(frozen=True)
class ReciprocalVolume:
    """Intensity on Cartesian momentum-transfer axes in rad/Å (Å^-1).

    The caller must explicitly identify and validate the source coordinate
    system; arrays in reciprocal-lattice units cannot be treated as Cartesian Q.
    """
    intensity: np.ndarray
    q_axes: tuple[np.ndarray, np.ndarray, np.ndarray]
    coordinate_system: str
    units: str
    preprocessing: tuple[str, ...] = ()
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        if self.coordinate_system != "cartesian_q" or self.units not in ("1/angstrom", "angstrom^-1"):
            raise ValueError("3D transform requires Cartesian Q in Å^-1, not hkl/rlu")
        v = np.asarray(self.intensity, dtype=float)
        if v.ndim != 3 or not np.isfinite(v).all():
            raise ValueError("Intensity must be a finite 3D array")
        if len(self.q_axes) != 3:
            raise ValueError("Exactly 3 Q axes required")
        axes = tuple(np.asarray(q, dtype=float) for q in self.q_axes)
        for i, (q, n) in enumerate(zip(axes, v.shape)):
            if q.ndim != 1 or len(q) != n or n < 2 or not np.isfinite(q).all():
                raise ValueError(f"Invalid Q axis {i}")
            d = np.diff(q)
            if np.any(d <= 0) or not np.allclose(d, d[0], rtol=1e-5, atol=1e-10):
                raise ValueError(f"Q axis {i} must be uniformly increasing")
            centered = (np.arange(n) - n // 2) * d[0]
            if not np.allclose(q, centered, rtol=1e-5, atol=max(1e-8, abs(d[0])*1e-5)):
                raise ValueError(f"Q axis {i} must be centered on zero")
        object.__setattr__(self, "intensity", v)
        object.__setattr__(self, "q_axes", axes)

    def require_bragg_subtracted(self) -> None:
        if "bragg_subtracted" not in self.preprocessing and "bragg_excluded" not in self.preprocessing:
            raise ValueError("3D-Delta PDF requires explicitly Bragg-subtracted/excluded input")
