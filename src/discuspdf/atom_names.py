"""Atom names as DISCUS 6.23.01 reads them.

DISCUS identifies an atom by its name (up to 4 characters, case ignored). It has a table of names it knows
(:data:`DISCUS_NAMES`: elements, ions such as ``FE3+`` or ``O2-``, ``D`` for deuterium, ``H.``/``C.``/``SI.``).
Any other name is guessed (``lib_f90`` ``guess_element``): the name up to a ``+``/``-`` sign at its 3rd or 4th
character, then its first two characters, then its first character. Names are compared upper case; DISCUS
keeps only 4 characters, so longer names are refused by the writers (not truncated here). So ``Ni62`` and ``C1`` are labels that scatter as Ni and C unless they are
given a scattering length of their own (:class:`~discuspdf.FourierConfig` and :class:`~discuspdf.PDFConfig`
``scattering_lengths``).

>>> discus_atom_name("Ni62")
AtomName(name='NI62', known=False, read_as='NI', element='Ni', b=10.3)
>>> element_of("D"), mass_number("D"), mass_number("Ni62"), mass_number("C1")
('H', 2, 62, None)
"""
from __future__ import annotations

import re
from dataclasses import dataclass

from ._discus_names import DISCUS_NAME_TABLE, SOURCE_SHA256
from ._elements import _ELEMENTS

__all__ = ["AtomName", "DISCUS_NAMES", "discus_atom_name", "element_of", "atomic_number", "mass_number",
           "is_isotope_label", "discus_name_rows", "SOURCE_SHA256"]

#: name -> (element DISCUS uses for the X-ray form factor, coherent neutron b in fm). Upper case.
DISCUS_NAMES: dict[str, tuple[str, float]] = {n: (e, b) for n, e, b in DISCUS_NAME_TABLE}


@dataclass(frozen=True)
class AtomName:
    """How DISCUS reads one atom name."""
    name: str       #: the name, upper case
    known: bool     #: DISCUS has the name in its table (else it is guessed)
    read_as: str    #: the table entry DISCUS uses
    element: str    #: element of that entry, e.g. 'Ni', 'D' (deuterium), 'Void'
    b: float        #: DISCUS's coherent neutron scattering length of that entry (fm)


def _title(s: str) -> str:
    return s[:1] + s[1:].lower()


def discus_atom_name(name) -> AtomName | None:
    """The table entry DISCUS uses for ``name``, or ``None`` if it cannot read the name at all."""
    up = str(name if name is not None else "").strip().upper()
    if not up:
        return None
    # Names longer than 4 characters are not truncated here (the page does the same); writers and the page
    # refuse them, because DISCUS cuts them to 4 characters.
    if up in DISCUS_NAMES:
        el, b = DISCUS_NAMES[up]
        return AtomName(up, True, up, _title(el), b)
    k = 3 if up[2:3] in ("+", "-") and len(up) > 2 else 4 if up[3:4] in ("+", "-") and len(up) > 3 else 0
    tries = ([up[:k]] if k else []) + [up[:2], up[:1]]
    for t in tries:
        if t in DISCUS_NAMES:
            el, b = DISCUS_NAMES[t]
            return AtomName(up, False, t, _title(el), b)
    return None


def element_of(name) -> str | None:
    """Chemical element of an atom name (``D`` counts as H); ``None`` for VOID and unreadable names."""
    info = discus_atom_name(name)
    if info is None:
        return None
    el = "H" if info.element == "D" else info.element
    return el if el in _ELEMENTS else None


def atomic_number(name) -> int:
    """Atomic number of the element DISCUS reads for ``name`` (D: 1). Raises ``ValueError`` otherwise."""
    el = element_of(name)
    if el is None:
        raise ValueError(f"atom name {name!r} is not an element and does not start with one")
    return _ELEMENTS[el]


_MASS_RE = re.compile(r"[A-Za-z]{1,2}(\d{1,3})")


def mass_number(name) -> int | None:
    """Mass number written after the element (``Ni62``, ``O18``, ``Li7``); ``D`` is 2.

    Numbers outside 2Z−1 … 3Z+20 are read as plain labels (``C1``, ``Fe2``), and hydrogen labels are never
    isotopes (``H2`` is a label; write ``D``)."""
    info = discus_atom_name(name)
    if info is None:
        return None
    if info.element == "D":
        return 2
    m = _MASS_RE.fullmatch(str(name).strip())
    el = element_of(name)
    z = _ELEMENTS.get(el, 0) if el else 0
    if not m or z <= 1:
        return None
    a = int(m.group(1))
    return a if 2 * z - 1 <= a <= 3 * z + 20 else None


def is_isotope_label(name) -> bool:
    """A name DISCUS does not know that carries a mass number, e.g. ``Ni62`` (it is read as plain Ni)."""
    info = discus_atom_name(name)
    return bool(info and not info.known and mass_number(name))


def discus_name_rows(*, scattering_only: bool = True) -> list[tuple[str, str, float]]:
    """The table rows in DISCUS's order; by default without VOID and the non-scattering special names."""
    return [r for r in DISCUS_NAME_TABLE if not scattering_only or r[1] != "VOID"]


# One rule for every macro writer: a letter, then up to 3 letters, digits, '+' or '-' (ions such as O2-).
_DISCUS_NAME_RE = re.compile(r"^[A-Za-z][A-Za-z0-9+\-]{0,3}$")


def check_discus_name(label: str, v, *, allow_number: bool = False, allow_all: bool = False):
    """An atom name as a DISCUS macro command reads it, stripped: a letter, then up to 3 letters, digits, '+'
    or '-' (DISCUS keeps 4 characters). With ``allow_number`` a non-negative integer
    (a scattering-type number) is returned as an ``int``; ``all`` is accepted only with ``allow_all``.
    Shared by every macro writer (``mc``, ``fourier_macro`` and the suite-interface steps)."""
    import numbers

    if allow_number and isinstance(v, numbers.Integral) and not isinstance(v, bool):
        if v < 0:
            raise ValueError(f"{label} must be a type number >= 0, got {v!r}")
        return int(v)
    if not isinstance(v, str):
        raise ValueError(f"{label} must be an atom name{' or type number' if allow_number else ''}, got {v!r}")
    text = v.strip()
    if text.lower() == "all":
        if not allow_all:
            raise ValueError(f"{label}: 'all' is not allowed here")
        return text
    if len(text) > 4:
        raise ValueError(f"{label} {v!r} is not a DISCUS atom name: longer than 4 characters (DISCUS keeps only 4)")
    if not _DISCUS_NAME_RE.match(text):
        raise ValueError(f"{label} must be a DISCUS atom name (a letter, then up to 3 letters, digits, '+' or '-'), "
                         f"got {v!r}")
    return text
