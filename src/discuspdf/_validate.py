"""Shared argument checks for the DISCUS_SUITE wrappers (one copy instead of one per module).

Each function returns the cleaned value or raises ``ValueError`` with the parameter's name in the message.
Booleans are refused wherever a number is expected (``True`` is an ``int`` in Python).
"""
from __future__ import annotations

import math
import numbers
import re

__all__ = ["PLAIN_FILE_RE", "real", "positive", "integer", "probability", "vec3", "int3", "plain_file_name"]

# A file name DISCUS reads as one macro parameter: no blanks, commas, quotes, comment characters or paths.
PLAIN_FILE_RE = re.compile(r"^[A-Za-z0-9_.\-]+$")


def real(name: str, v) -> float:
    """A finite real number."""
    if isinstance(v, bool) or not isinstance(v, numbers.Real) or not math.isfinite(float(v)):
        raise ValueError(f"{name} must be a finite number, got {v!r}")
    return float(v)


def positive(name: str, v, *, allow_zero: bool = False) -> float:
    """A finite number > 0 (or >= 0 with ``allow_zero``)."""
    x = real(name, v)
    if x < 0 or (x == 0 and not allow_zero):
        raise ValueError(f"{name} must be {'>= 0' if allow_zero else '> 0'}, got {v!r}")
    return x


def integer(name: str, v, minimum: int = 1) -> int:
    """An integer >= ``minimum``."""
    if isinstance(v, bool) or not isinstance(v, numbers.Integral) or v < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}, got {v!r}")
    return int(v)


def probability(name: str, v) -> float:
    """A probability in [0, 1] (not a relative move weight: those use :func:`positive`)."""
    if isinstance(v, bool) or not isinstance(v, numbers.Real) or not 0.0 <= float(v) <= 1.0:
        raise ValueError(f"{name} must be a probability in [0, 1], got {v!r}")
    return float(v)


def _three(name: str, v, what: str) -> tuple:
    try:
        items = tuple(v)
    except TypeError:
        raise ValueError(f"{name} must be three {what}, got {v!r}") from None
    if len(items) != 3:
        raise ValueError(f"{name} must be three {what}, got {v!r}")
    return items


def vec3(name: str, v) -> tuple[float, float, float]:
    """Three finite numbers."""
    return tuple(real(f"{name}[{i}]", x) for i, x in enumerate(_three(name, v, "numbers")))


def int3(name: str, v) -> tuple[int, int, int]:
    """Three integers (any sign)."""
    items = _three(name, v, "integers")
    if any(isinstance(x, bool) or not isinstance(x, numbers.Integral) for x in items):
        raise ValueError(f"{name} must be three integers, got {v!r}")
    return tuple(int(x) for x in items)


def plain_file_name(name: str, v) -> str:
    """A plain file name: letters, digits, '_', '-', '.'; not '.' or '..'."""
    if not isinstance(v, str) or not PLAIN_FILE_RE.match(v) or v in (".", ".."):
        raise ValueError(f"{name} must be a plain file name (letters, digits, '_', '-', '.'), got {v!r}")
    return v
