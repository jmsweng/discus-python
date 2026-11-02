"""DiffPy adapter (Package A, Feature A2).

Converts between :class:`~discuspdf.structure.Structure` and
``diffpy.structure.Structure``.  ``diffpy.structure`` is imported inside the
functions so that ``import discuspdf`` works without the ``diffpy`` extra.

Orientation is NOT preserved (DiffPy defines its own Cartesian frame);
fractional coordinates, lattice parameters and all distances/angles are.
"""
from __future__ import annotations

import numpy as np

from .lattice import Lattice
from .structure import Structure
from .interop import (
    InteropError,
    element_from_type_symbol,
    b_to_u,
    u_to_b,
    LossReport,
)


def _import_diffpy():
    """Import diffpy.structure, raising a helpful ``ImportError`` if missing."""
    try:
        import diffpy.structure
        return diffpy.structure
    except ImportError:
        raise ImportError(
            "diffpy.structure is required for DiffPy conversion; "
            "install it with: pip install 'discus-python[diffpy]'"
        ) from None


def to_diffpy(structure: Structure, *, allow_loss: bool = False):
    """Convert a :class:`Structure` to a ``diffpy.structure.Structure``.

    Lattice parameters are recovered via ``Lattice.parameters`` and passed to
    ``diffpy.Lattice``.  Orientation is NOT preserved (DiffPy defines its own
    Cartesian frame); fractional coordinates, lattice parameters and all
    distances/angles are.

    The model is not stored; pass it again to :func:`from_diffpy`.

    Losses (design rule 1): atom_id not 0..N-1, non-trivial cell_index/
    supercell, pbc not all True, and species charges.
    """
    diffpy = _import_diffpy()

    a, b, c, alpha, beta, gamma = structure.lattice.parameters
    lattice = diffpy.Lattice(a, b, c, alpha, beta, gamma)

    n = structure.n_atoms
    atoms = []
    for i in range(n):
        element = element_from_type_symbol(structure.species[i])
        xyz = tuple(structure.frac[i])
        label = structure.site_label[i]
        occupancy = float(structure.occupancy[i])
        Uisoequiv = (
            b_to_u(structure.b_iso[i]) if structure.b_iso is not None else 0.0
        )
        atoms.append(diffpy.Atom(
            atype=element,
            xyz=xyz,
            label=label,
            occupancy=occupancy,
            Uisoequiv=Uisoequiv,
        ))

    ds = diffpy.Structure(lattice=lattice, atoms=atoms)

    # Losses (design rule 1).
    losses = LossReport()

    if not np.array_equal(structure.atom_id, np.arange(n)):
        losses.add("atom_id is not 0..N-1")

    if np.any(structure.cell_index != 0) or tuple(structure.supercell) != (1, 1, 1):
        losses.add("cell_index/supercell not trivial")

    if structure.pbc != (True, True, True):
        losses.add("pbc not all True")

    # Species charges: any species whose type symbol differs from the bare
    # element symbol carries a charge (or number) suffix.
    if any(element_from_type_symbol(s) != s for s in structure.species):
        losses.add("species charges not representable")

    losses.raise_or_warn(allow_loss)

    return ds


def from_diffpy(
    dstructure,
    *,
    model: str,
    allow_loss: bool = False,
    min_separation: float | None = 0.1,
) -> Structure:
    """Convert a ``diffpy.structure.Structure`` to a :class:`Structure`.

    ``model`` is required (no default).  The lattice is recovered via
    ``lattice.cell_parms()`` — NEVER ``.base`` (different orientation).
    Anisotropic atoms are a loss (the isotropic equivalent is kept).
    ``min_separation`` (Å, or None to skip) is passed to :class:`Structure`.
    """
    diffpy = _import_diffpy()

    # ``cell_parms`` is the replacement for the deprecated ``abcABG``;
    # fall back to the individual attributes (present in all versions).
    lat = dstructure.lattice
    a, b, c, alpha, beta, gamma = lat.a, lat.b, lat.c, lat.alpha, lat.beta, lat.gamma
    lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)

    atoms = list(dstructure)
    n = len(atoms)

    species = np.empty(n, dtype=object)
    frac = np.empty((n, 3), dtype=float)
    occupancy = np.empty(n, dtype=float)
    labels = np.empty(n, dtype=object)
    u_values = np.zeros(n, dtype=float)

    losses = LossReport()
    for i, atom in enumerate(atoms):
        species[i] = atom.element
        frac[i] = np.array(atom.xyz, dtype=float)
        occupancy[i] = float(atom.occupancy)
        u_values[i] = float(atom.Uisoequiv)
        labels[i] = atom.label
        if atom.anisotropy:
            losses.add(
                "anisotropic atoms not fully representable "
                "(isotropic equivalent kept)"
            )

    # B factor: B = 8π²·U.  If every Uisoequiv is 0 → b_iso = None.
    b_iso = u_to_b(u_values)
    if np.all(b_iso == 0.0):
        b_iso = None

    losses.raise_or_warn(allow_loss)

    return Structure(
        lattice,
        species,
        frac,
        model=model,
        occupancy=occupancy,
        site_label=labels,
        b_iso=b_iso,
        min_separation=min_separation,
    )
