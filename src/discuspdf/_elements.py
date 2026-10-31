"""Small element table for atomic-number lookup.

Maps element symbols (case-insensitive, with optional charge or number suffixes
stripped) to atomic numbers.  Only the first ~100 elements are included;
unknown labels raise ``ValueError`` when ``atomic_numbers`` is requested.
"""
from __future__ import annotations

import re

# Element symbol -> atomic number (first 100 elements + a few common ions).
_ELEMENTS: dict[str, int] = {
    "H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8,
    "F": 9, "Ne": 10, "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15,
    "S": 16, "Cl": 17, "Ar": 18, "K": 19, "Ca": 20, "Sc": 21, "Ti": 22,
    "V": 23, "Cr": 24, "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29,
    "Zn": 30, "Ga": 31, "Ge": 32, "As": 33, "Se": 34, "Br": 35, "Kr": 36,
    "Rb": 37, "Sr": 38, "Y": 39, "Zr": 40, "Nb": 41, "Mo": 42, "Tc": 43,
    "Ru": 44, "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48, "In": 49, "Sn": 50,
    "Sb": 51, "Te": 52, "I": 53, "Xe": 54, "Cs": 55, "Ba": 56, "La": 57,
    "Ce": 58, "Pr": 59, "Nd": 60, "Pm": 61, "Sm": 62, "Eu": 63, "Gd": 64,
    "Tb": 65, "Dy": 66, "Ho": 67, "Er": 68, "Tm": 69, "Yb": 70, "Lu": 71,
    "Hf": 72, "Ta": 73, "W": 74, "Re": 75, "Os": 76, "Ir": 77, "Pt": 78,
    "Au": 79, "Hg": 80, "Tl": 81, "Pb": 82, "Bi": 83, "Po": 84, "At": 85,
    "Rn": 86, "Fr": 87, "Ra": 88, "Ac": 89, "Th": 90, "Pa": 91, "U": 92,
    "Np": 93, "Pu": 94, "Am": 95, "Cm": 96, "Bk": 97, "Cf": 98, "Es": 99,
    "Fm": 100,
}

# Strip charge suffixes like "O2-", "Na+", "Cl-" and number suffixes like "Fe1".
_SUFFIX_RE = re.compile(r"^(.+?)(?:[+-]\d+|\d+)?$")


def _strip_suffix(label: str) -> str:
    """Strip charge or number suffixes from a species label.

    Examples: "O2-" -> "O", "Na+" -> "Na", "Fe1" -> "Fe", "GA" -> "GA".
    """
    m = _SUFFIX_RE.match(label)
    if m:
        return m.group(1)
    return label


def atomic_number(species: str) -> int:
    """Return the atomic number for a species label.

    The label is matched case-insensitively after stripping any charge or
    number suffix.  Unknown labels raise ``ValueError``.
    """
    if not isinstance(species, str) or not species:
        raise ValueError(f"species label must be a non-empty string, got {species!r}")
    stripped = _strip_suffix(species)
    # Try exact match first, then case-insensitive.
    if stripped in _ELEMENTS:
        return _ELEMENTS[stripped]
    upper = stripped.upper()
    if upper in _ELEMENTS:
        return _ELEMENTS[upper]
    # Try with first letter upper, rest lower (standard element symbol form).
    if len(stripped) >= 2:
        canonical = stripped[0].upper() + stripped[1:].lower()
        if canonical in _ELEMENTS:
            return _ELEMENTS[canonical]
    raise ValueError(
        f"unknown element symbol {species!r} (stripped: {stripped!r}); "
        "the element table covers the first 100 elements"
    )
