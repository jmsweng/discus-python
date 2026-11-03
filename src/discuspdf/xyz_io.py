"""Extended XYZ output (readable by ASE, OVITO, VESTA's XYZ import and most viewers), without ASE.

One frame: the atom count, a comment line with the periodic frame (``Lattice="ax ay az bx by bz cx cy cz"``,
Å, rows a, b, c of :attr:`Structure.lattice`), ``Properties=species:S:1:pos:R:3`` and ``pbc``, then one line
per atom with its species and Cartesian position (Å). Coordinates are written as they are (no wrapping).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from .structure import Structure

__all__ = ["format_extxyz", "write_extxyz"]


def format_extxyz(structure: Structure, *, mask=None, decimals: int = 6) -> str:
    """Extended XYZ text of ``structure`` (or of the atoms where boolean ``mask`` is True)."""
    if not isinstance(structure, Structure):
        raise TypeError("structure must be a discuspdf Structure")
    keep = np.ones(structure.n_atoms, bool) if mask is None else np.asarray(mask, bool)
    if keep.shape != (structure.n_atoms,):
        raise ValueError(f"mask must have one entry per atom ({structure.n_atoms})")
    species = [str(s) for s in np.asarray(structure.species, dtype=object)[keep]]
    for s in species:
        if not s or any(c.isspace() for c in s) or '"' in s:
            raise ValueError(f"species {s!r} cannot be written to XYZ (empty, spaces or quotes)")
    cart = np.asarray(structure.cartesian, float)[keep]
    d = int(decimals)
    lat = " ".join(f"{x:.{d}f}" for x in np.asarray(structure.lattice.matrix, float).ravel())
    pbc = " ".join("T" if p else "F" for p in structure.pbc)
    lines = [str(len(species)), f'Lattice="{lat}" Properties=species:S:1:pos:R:3 pbc="{pbc}"']
    lines += [f"{s} {x:.{d}f} {y:.{d}f} {z:.{d}f}" for s, (x, y, z) in zip(species, cart)]
    return "\n".join(lines) + "\n"


def write_extxyz(path, structure: Structure, *, mask=None, decimals: int = 6, overwrite: bool = False) -> Path:
    """Write :func:`format_extxyz` text to ``path`` (ASCII, LF). Refuses an existing file unless ``overwrite``."""
    from ._safe_write import atomic_write_text
    return atomic_write_text(Path(path), format_extxyz(structure, mask=mask, decimals=decimals),
                             overwrite=overwrite)
