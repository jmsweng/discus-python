"""Powder patterns from DISCUS's ``powder`` menu (Debye calculation) as a :class:`~discuspdf.discus_ops.Step`.

The 1D PDF route through the powder menu is :mod:`discuspdf.pdf`. This step writes the powder-side
quantities themselves: ``inte`` (intensity with polarisation), ``S(Q)``, ``F(Q)`` = Q (S − 1) and ``I(Q)``.
Facts from the 6.23.01 source (``docs/DISCUS_POWDER_FACTS.md``; the items used here were checked):

* ``set calc, debye``: the crystal as read is one finite cluster (no periodic images); pair distances
  are binned to multiples of ``stepr`` (``set stepr``, start value 0.01 Å, **not** restored by ``reset``, so
  the step always writes it) and summed with a 65536-entry sine table (``powder_pdf_hist.f90``).
* For one species, S(Q) = 1 + (2/N) Σ_bins H(r̃) sin(Q r̃)/(Q r̃) (``phases_average``); S(Q) is divided
  by its mean above 18 Å⁻¹ when more than 50 points lie there, so keep ``qmax`` ≤ 18.
* **Label shift:** the value written at label Q is the one computed at Q − dq (the Debye result is
  stored one index up, ``powder_write_mod.f90:1351``). :func:`read_powder_file` returns the file as written;
  ``shift_labels=True`` moves the labels down by dq.
* ``reset`` restores a pseudo-Voigt profile (w = 0.05); the step writes ``set profile, off``. ``set axis``
  always gives Q in 6.23.01, so no 2θ axis is offered. ``set temp`` and ``set disp`` are not restored by
  ``reset`` either and are always written.
* ``dq`` is rounded by DISCUS to 1e-5 Å⁻¹ and doubled silently if the grid would exceed 2^18 points;
  both are refused here. ``qmax`` · λ / 4π must stay below 1.
"""
from __future__ import annotations

import math
from dataclasses import dataclass


from .discus_ops import Step
from .fourier_macro import ENERGY_UNITS, _RADIATION_KEYWORD, RADIATIONS, WAVELENGTH_ANGSTROM, WAVELENGTHS
from .macro import Macro
from ._validate import PLAIN_FILE_RE as _NAME, real as _num

__all__ = ["POWDER_VALUES", "PowderPattern", "read_powder_file"]

POWDER_VALUES = ("inte", "S(Q)", "F(Q)", "I(Q)")

@dataclass(frozen=True)
class PowderPattern(Step):
    """``powder`` (Debye) then one ``output`` block per entry of ``outputs`` ((value, file), value in
    :data:`POWDER_VALUES`), each written on the Q grid qmin, qmin + dq, …, qmax.

    ``wavelength`` (Å, or a K-alpha1 symbol ``moa1``/``cua1``/``aga1``) or ``energy`` (keV for X-rays and
    electrons, meV for neutrons). ``temperature='use'`` applies exp(−B s²) to the form factors.
    ``polarisation`` (``set lpcor``) affects ``inte`` only: ``none`` (default here) or DISCUS's ``bragg``."""
    outputs: tuple
    qmin: float = 0.5
    qmax: float = 15.0
    dq: float = 0.01
    radiation: str = "xray"
    wavelength: str | float | None = "moa1"
    energy: float | None = None
    stepr: float = 0.01
    temperature: str = "ignore"
    polarisation: str = "none"

    def __post_init__(self) -> None:
        outs = tuple(tuple(o) for o in self.outputs)
        if not outs:
            raise ValueError("at least one output is required")
        names = set()
        for o in outs:
            if len(o) != 2 or o[0] not in POWDER_VALUES:
                raise ValueError(f"each output is (value, file) with value in {POWDER_VALUES}, got {o!r}")
            if not isinstance(o[1], str) or not _NAME.match(o[1]) or o[1] in names:
                raise ValueError(f"output file {o[1]!r} must be a plain, unique file name")
            names.add(o[1])
        object.__setattr__(self, "outputs", outs)
        if self.radiation not in RADIATIONS:
            raise ValueError(f"radiation must be one of {RADIATIONS}")
        qmin, qmax, dq = _num("qmin", self.qmin), _num("qmax", self.qmax), _num("dq", self.dq)
        if not 0 < qmin < qmax:
            raise ValueError("need 0 < qmin < qmax")
        if dq <= 0 or abs(round(dq * 1e5) - dq * 1e5) > 1e-6:
            raise ValueError("dq must be positive and a multiple of 1e-5 (DISCUS rounds it)")
        if (qmax - qmin) / dq + 1 > 2 ** 18:
            raise ValueError("more than 2^18 Q points: DISCUS would double dq silently")
        if abs(round(qmin / dq) - qmin / dq) > 1e-9:
            raise ValueError("qmin must be a multiple of dq (the written grid is then the calculated one)")
        if qmax > 18.0:
            raise ValueError("qmax > 18: DISCUS divides S(Q) by its mean above 18 A^-1")
        if self.energy is not None:
            if self.wavelength not in ("moa1", None):
                raise ValueError("give either wavelength or energy, not both")
            if _num("energy", self.energy) <= 0:
                raise ValueError(f"energy must be positive ({ENERGY_UNITS[self.radiation]})")
            lam = None
        elif isinstance(self.wavelength, str):
            if self.wavelength.lower() not in WAVELENGTHS:
                raise ValueError(f"unknown wavelength symbol {self.wavelength!r}")
            if self.radiation != "xray":
                raise ValueError("give neutrons and electrons a wavelength in A or an energy")
            lam = WAVELENGTH_ANGSTROM[self.wavelength.lower()]
        else:
            lam = _num("wavelength", self.wavelength)
            if lam <= 0:
                raise ValueError("wavelength must be positive")
        if lam is not None and qmax * lam / (4 * math.pi) >= 1:
            raise ValueError(f"qmax {qmax} is beyond 4 pi / lambda for lambda = {lam} A")
        if _num("stepr", self.stepr) <= 0:
            raise ValueError("stepr must be positive")
        if self.temperature not in ("ignore", "use"):
            raise ValueError("temperature must be 'ignore' or 'use'")
        if self.polarisation not in ("none", "bragg"):
            raise ValueError("polarisation must be 'none' or 'bragg'")

    def output_names(self) -> tuple[str, ...]:
        return tuple(name for _, name in self.outputs)

    def apply(self, m: Macro) -> None:
        with m.menu("powder"):
            m.command("reset")
            m.command(_RADIATION_KEYWORD[self.radiation])
            if self.energy is not None:
                m.command("set", "energy", self.energy)
            else:
                m.command("set", "wvle", self.wavelength)
            m.command("set", "calc", "debye")
            m.command("set", "qmin", self.qmin)
            m.command("set", "qmax", self.qmax)
            m.command("set", "dq", self.dq)
            m.command("set", "stepr", self.stepr)
            m.command("set", "profile", "off")
            m.command("set", "lpcor", self.polarisation)
            m.command("set", "temp", self.temperature)
            m.command("set", "disp", "off")
            m.command("set", "period", "off")
            m.command("run")
        for value, name in self.outputs:
            with m.menu("output"):
                m.command("outf", name)
                m.command("form", "powd", "q", self.qmin, self.qmax, self.dq)
                m.command("value", value)
                m.command("run")


def read_powder_file(path, *, dq: float | None = None, shift_labels: bool = False):
    """(Q, y) from a powder output file (two columns, no header). With ``shift_labels=True`` (and ``dq``)
    the labels are moved down by dq, to the Q at which DISCUS computed the values."""
    from .pdf import read_pdf  # the same two-column reader (comments, D exponents, increasing x)

    q, y = read_pdf(path)
    if shift_labels:
        if dq is None:
            raise ValueError("shift_labels needs dq")
        q = q - dq
    return q, y
