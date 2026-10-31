"""Shared interop helpers: B↔U conversion, element-symbol extraction, loss handling.

These are used by the CIF, ASE, and DiffPy adapters.  The B↔U conversion is
``B = 8π²·U`` (gemmi has no ``u_to_b``; see ``docs/INTEROP_FACTS.md`` P12).
"""
from __future__ import annotations

import math
import re

from ._elements import _ELEMENTS


class InteropError(ValueError):
    """An interop operation cannot be completed without data loss or an invalid
    result.  The message lists every loss (one per line)."""


class InteropWarning(UserWarning):
    """An interop operation proceeded with data loss (``allow_loss=True``)."""


#: Conversion factor from isotropic U (Å²) to isotropic B (Å²): B = 8π²·U.
B_FROM_U = 8.0 * math.pi ** 2


def u_to_b(u: float) -> float:
    """Convert isotropic U (Å²) to isotropic B (Å²): ``B = 8π²·U``."""
    return B_FROM_U * u


def b_to_u(b: float) -> float:
    """Convert isotropic B (Å²) to isotropic U (Å²): ``U = B / (8π²)``."""
    return b / B_FROM_U


def element_from_type_symbol(s: str) -> str:
    """Extract the element symbol from a CIF ``_atom_site_type_symbol`` value.

    Strips charge and number suffixes: ``O2-`` → ``O``, ``Ga3+`` → ``Ga``,
    ``GA`` → ``Ga``.  Unknown symbols raise :class:`InteropError`.
    """
    if not isinstance(s, str) or not s:
        raise InteropError(f"empty type_symbol")
    # Strip charge suffixes (e.g. "O2-", "Ga3+") and number suffixes (e.g. "Fe1").
    # The element symbol is the leading run of letters.
    m = re.match(r"([A-Za-z]+)", s)
    if m is None:
        raise InteropError(f"no element symbol found in {s!r}")
    stripped = m.group(1)
    # Try exact match first, then case-insensitive.
    if stripped in _ELEMENTS:
        return stripped
    upper = stripped.upper()
    if upper in _ELEMENTS:
        return upper
    # Try with first letter upper, rest lower (standard element symbol form).
    if len(stripped) >= 2:
        canonical = stripped[0].upper() + stripped[1:].lower()
        if canonical in _ELEMENTS:
            return canonical
    raise InteropError(
        f"unknown element symbol {s!r} (stripped: {stripped!r}); "
        "the element table covers the first 100 elements"
    )


class LossReport:
    """Accumulates loss messages and raises or warns according to policy.

    Usage::

        losses = LossReport()
        if aniso:
            losses.add("anisotropic ADPs not stored")
        losses.raise_or_warn(allow_loss)
    """

    def __init__(self) -> None:
        self._losses: list[str] = []

    def add(self, msg: str) -> None:
        """Record a loss message (one per line)."""
        if msg not in self._losses:
            self._losses.append(msg)

    @property
    def losses(self) -> tuple[str, ...]:
        """The accumulated loss messages as a tuple."""
        return tuple(self._losses)

    def raise_or_warn(self, allow_loss: bool) -> None:
        """Raise :class:`InteropError` listing all losses, or emit one
        :class:`InteropWarning` carrying the lines (``allow_loss=True``)."""
        if not self._losses:
            return
        text = "\n".join(self._losses)
        if allow_loss:
            import warnings
            warnings.warn(text, InteropWarning, stacklevel=3)
        else:
            raise InteropError(text)
