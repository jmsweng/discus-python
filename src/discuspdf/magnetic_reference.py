"""EXPERIMENTAL magnetic neutron scattering reference (Package D, Feature 11).

DISCUS 6.23.01 computes **no** magnetic scattering: it can hold a moment per atom but never uses it in a
Fourier, powder or PDF calculation, and its .stru reader/writer drop moments (see
``docs/MAGNETIC_SCATTERING.md``). This module is therefore an independent, deliberately small reference
for analytical fixtures, not a wrapper and not a production feature.

Model (elastic, unpolarised neutrons, spin-only dipole approximation, Squires ch. 7):

    dσ/dΩ = p² |F⊥(Q)|²,   F(Q) = Σ_j f_j(|Q|) μ_j exp(i Q·r_j),   F⊥ = F − (F·Q̂) Q̂

* ``μ_j`` — Cartesian moment in Bohr magnetons (μB), the same Cartesian frame as the positions (Å);
* ``f_j`` — magnetic form factor, ``⟨j0⟩ + (1 − 2/g)⟨j2⟩`` with ``f(0) = 1``;
* ``p = γ r0 / 2 = 2.6954 fm`` (γ = 1.91304276, r0 = 2.8179403262 fm); intensities in fm²;
* ``Q`` in Å⁻¹ (2π included). At ``|Q| < q_min`` Q̂ is undefined: the intensity is NaN and the returned
  mask is False there. Nothing is returned for a Q = 0 forward beam as if it were a number.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Mapping, Sequence

import numpy as np

__all__ = ["P_MAGNETIC_FM", "MagneticFormFactor", "magnetic_structure_factor", "magnetic_intensity"]

GAMMA_NEUTRON = 1.91304276          # |γ|, neutron magnetic moment in nuclear magnetons (CODATA 2018)
R0_FM = 2.8179403262                # classical electron radius, fm (CODATA 2018)
P_MAGNETIC_FM = GAMMA_NEUTRON * R0_FM / 2.0     # 2.6954 fm


def _coeffs(name: str, c) -> np.ndarray:
    c = np.asarray(c, float)
    if c.shape != (7,) or not np.isfinite(c).all():
        raise ValueError(f"{name} needs 7 finite coefficients (A, a, B, b, C, c, D)")
    return c


def _radial(c: np.ndarray, s2: np.ndarray) -> np.ndarray:
    A, a, B, b, C, cc, D = c
    return A * np.exp(-a * s2) + B * np.exp(-b * s2) + C * np.exp(-cc * s2) + D


@dataclass(frozen=True)
class MagneticFormFactor:
    """Dipole magnetic form factor ``f(Q) = ⟨j0⟩ + (1 − 2/g)⟨j2⟩``, ``s = |Q|/4π``.

    ``g`` is the spectroscopic g-factor of a spin-dominated (3d) ion, where L = (g − 2)S. Rare-earth ions
    described by a Landé g_J use the opposite sign, ⟨j0⟩ + (2/g_J − 1)⟨j2⟩; that convention is not provided
    (``from_periodictable`` refuses g ≠ 2 for lanthanides and actinides).

    ``⟨j0⟩ = A e^{−a s²} + B e^{−b s²} + C e^{−c s²} + D`` and ``⟨j2⟩ = s² (same form)`` (International
    Tables C, Brown). ``j2`` may be omitted only for ``g = 2`` (spin only)."""
    j0: tuple
    source: str
    j2: tuple | None = None
    g: float = 2.0

    def __post_init__(self) -> None:
        object.__setattr__(self, "j0", tuple(_coeffs("j0", self.j0)))
        if self.j2 is not None:
            object.__setattr__(self, "j2", tuple(_coeffs("j2", self.j2)))
        if not np.isfinite(self.g) or self.g <= 0:
            raise ValueError("g must be positive")
        if self.g != 2.0 and self.j2 is None:
            raise ValueError("g != 2 needs <j2> coefficients")
        if not isinstance(self.source, str) or not self.source:
            raise ValueError("state the source of the form-factor coefficients")

    @classmethod
    def from_periodictable(cls, element: str, charge: int, g: float = 2.0) -> "MagneticFormFactor":
        """Coefficients from ``periodictable`` (Brown, International Tables C 4.4.5)."""
        try:
            import periodictable
        except ImportError as exc:                                     # pragma: no cover - env specific
            raise ImportError("from_periodictable needs the optional package 'periodictable'") from exc
        el = getattr(periodictable, element)
        if g != 2.0 and (57 <= el.number <= 71 or 89 <= el.number <= 103):
            raise ValueError(f"{element}: rare-earth/actinide form factors use the Landé g_J convention "
                             "(2/g_J − 1)⟨j2⟩, which this reference does not provide; use g = 2 or explicit "
                             "coefficients with a stated convention")
        if charge not in el.magnetic_ff:
            raise ValueError(f"periodictable has no magnetic form factor for {element}{charge:+d}")
        ff = el.magnetic_ff[charge]
        return cls(tuple(ff.j0), f"periodictable {periodictable.__version__} {element}{charge:+d}",
                   tuple(ff.j2) if ff.j2 is not None else None, g)

    def __call__(self, q) -> np.ndarray:
        q = np.asarray(q, float)
        s2 = (q / (4.0 * np.pi)) ** 2
        f = _radial(np.asarray(self.j0), s2)
        if self.j2 is not None and self.g != 2.0:
            f = f + (1.0 - 2.0 / self.g) * s2 * _radial(np.asarray(self.j2), s2)
        return f


def _as_form_factors(form_factors, n: int, species: Sequence[str] | None) -> list[Callable]:
    if callable(form_factors):
        return [form_factors] * n
    if isinstance(form_factors, Mapping):
        if species is None:
            raise ValueError("a form-factor mapping needs species= (one label per atom)")
        if len(species) != n:
            raise ValueError("species must have one label per atom")
        missing = sorted({s for s in species if s not in form_factors})
        if missing:
            raise KeyError(f"no magnetic form factor for {missing}")
        return [form_factors[s] for s in species]
    ffs = list(form_factors)
    if len(ffs) != n:
        raise ValueError("give one form factor per atom, a mapping by species, or a single form factor")
    return ffs


def _validate(positions, moments, q):
    r = np.asarray(positions, float)
    mu = np.asarray(moments, float)
    Q = np.asarray(q, float)
    if r.ndim != 2 or r.shape[1] != 3 or mu.shape != r.shape:
        raise ValueError("positions and moments must both be (N, 3) Cartesian arrays (Å, μB)")
    if Q.shape[-1:] != (3,):
        raise ValueError("q must have last dimension 3 (Cartesian Å⁻¹)")
    for name, a in (("positions", r), ("moments", mu), ("q", Q)):
        if not np.isfinite(a).all():
            raise ValueError(f"{name} must be finite")
    return r, mu, Q


def magnetic_structure_factor(positions_cartesian, moments_cartesian, q, form_factors, *,
                              species: Sequence[str] | None = None, q_min: float = 1e-9):
    """``(F⊥, valid)``: the perpendicular magnetic amplitude ``(..., 3)`` complex, in μB, and the mask of
    Q points where it is defined (``|Q| ≥ q_min``; F⊥ is NaN elsewhere). Multiply by ``P_MAGNETIC_FM`` for fm."""
    r, mu, Q = _validate(positions_cartesian, moments_cartesian, q)
    shape = Q.shape[:-1]
    Qf = Q.reshape(-1, 3)
    qn = np.linalg.norm(Qf, axis=1)
    ffs = _as_form_factors(form_factors, len(r), species)
    f = np.column_stack([np.asarray(ff(qn), float) for ff in ffs]) if len(r) else np.zeros((len(Qf), 0))
    if not np.isfinite(f).all():
        raise ValueError("a form factor returned non-finite values")
    phase = np.exp(1j * (Qf @ r.T))                       # (M, N)
    F = (f * phase) @ mu                                  # (M, 3)
    valid = qn >= q_min
    qhat = np.zeros_like(Qf)
    qhat[valid] = Qf[valid] / qn[valid, None]
    Fperp = F - np.sum(F * qhat, axis=1)[:, None] * qhat
    Fperp[~valid] = np.nan
    return Fperp.reshape(shape + (3,)), valid.reshape(shape)


def magnetic_intensity(positions_cartesian, moments_cartesian, q, form_factors, *,
                       species: Sequence[str] | None = None, q_min: float = 1e-9):
    """``(I, valid)``: unpolarised elastic magnetic intensity ``p² |F⊥|²`` in fm² for the whole sample.

    ``form_factors``: one callable for all atoms, one per atom, or a mapping by ``species`` label.
    ``I`` is NaN where ``|Q| < q_min`` (``valid`` False)."""
    Fp, valid = magnetic_structure_factor(positions_cartesian, moments_cartesian, q, form_factors,
                                          species=species, q_min=q_min)
    I = P_MAGNETIC_FM ** 2 * np.sum(np.abs(Fp) ** 2, axis=-1)
    return I, valid
