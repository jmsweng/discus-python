"""ASE adapter (Package A, Feature A2).

Converts between :class:`~discuspdf.structure.Structure` and ``ase.Atoms``.
``ase`` is imported inside the functions so that ``import discuspdf`` works
without the ``ase`` extra.
"""
from __future__ import annotations

import warnings

import numpy as np

from .lattice import Lattice
from .structure import Structure
from .interop import (
    InteropError,
    InteropWarning,
    element_from_type_symbol,
    LossReport,
)


def _import_ase():
    """Import ase, raising a helpful ``ImportError`` if missing."""
    try:
        import ase
        return ase
    except ImportError:
        raise ImportError(
            "ase is required for ASE conversion; "
            "install it with: pip install 'discus-python[ase]'"
        ) from None


def _pbc_from_ase(atoms) -> tuple[bool, bool, bool]:
    """Recover the 3-tuple pbc from an ASE Atoms object."""
    pbc = atoms.pbc
    if pbc is None:
        return (True, True, True)
    return tuple(bool(x) for x in pbc)


def to_ase(structure: Structure) -> "ase.Atoms":
    """Convert a :class:`Structure` to an ``ase.Atoms`` object.

    The full DISCUS state is stored in ``atoms.arrays`` (``dp_*`` keys) and
    ``atoms.info`` (``dp_model``, ``dp_supercell``) so that
    :func:`from_ase` can restore it exactly.  Extra arrays are set *after*
    construction because the ``ase.Atoms`` constructor does not accept an
    ``arrays=`` keyword argument.
    """
    ase = _import_ase()

    cell = structure.lattice.matrix
    pbc = structure.pbc
    symbols = [element_from_type_symbol(s) for s in structure.species]
    # Cartesian positions from the (unwrapped) fractional coordinates.
    positions = structure.frac @ structure.lattice.matrix

    atoms = ase.Atoms(symbols, cell=cell, pbc=pbc, positions=positions)

    # Extra arrays (set after construction; no arrays= kwarg).
    atoms.arrays['dp_species'] = np.array(structure.species, dtype=object)
    atoms.arrays['dp_site_label'] = np.array(structure.site_label, dtype=object)
    atoms.arrays['dp_atom_id'] = np.array(structure.atom_id, dtype=np.int64)
    atoms.arrays['dp_occupancy'] = np.array(structure.occupancy, dtype=float)
    atoms.arrays['dp_cell_index'] = np.array(structure.cell_index, dtype=np.int64)
    if structure.b_iso is not None:
        atoms.arrays['dp_b_iso'] = np.array(structure.b_iso, dtype=float)

    atoms.info['dp_model'] = structure.model
    atoms.info['dp_supercell'] = tuple(int(x) for x in structure.supercell)

    # Loss policy (design rule 2): an average-model structure with partial
    # occupancy or a mixed site cannot be represented exactly in ASE, where
    # every atom is fully occupied.  This is a warning, not a loss that
    # raises — to_ase always succeeds.
    if structure.model == "average":
        has_partial = bool(np.any(structure.occupancy < 1.0))
        has_mixed = False
        if not has_partial:
            # Group rows by (site_label, cell_index); any group with > 1 row
            # is a mixed site.
            n = structure.n_atoms
            keys = [
                (structure.site_label[i],
                 int(structure.cell_index[i, 0]),
                 int(structure.cell_index[i, 1]),
                 int(structure.cell_index[i, 2]))
                for i in range(n)
            ]
            counts: dict[tuple, int] = {}
            for k in keys:
                counts[k] = counts.get(k, 0) + 1
            has_mixed = any(c > 1 for c in counts.values())
        if has_partial or has_mixed:
            warnings.warn(
                "ASE consumers will treat these atoms as fully occupied; "
                "use the dp_* arrays",
                InteropWarning,
                stacklevel=2,
            )

    return atoms


def from_ase(
    atoms,
    *,
    model: str | None = None,
    allow_loss: bool = False,
    min_separation: float = 0.1,
) -> Structure:
    """Convert an ``ase.Atoms`` object to a :class:`Structure`.

    Two paths:

    * **Path A** — the Atoms object carries our ``dp_model`` info key
      (i.e. it came from :func:`to_ase`): the state is restored *exactly*
      from the ``dp_*`` arrays.
    * **Path B** — a foreign Atoms object (no ``dp_model``): ``model`` must
      be given, occupancies are assumed 1.0, and unsupported features are
      reported as losses via :class:`LossReport`.
    """
    ase = _import_ase()

    # Cell validation (the ASE constructor does not reject bad cells).
    cell = np.array(atoms.cell.array, dtype=float)
    if np.linalg.det(cell) <= 0:
        raise InteropError(
            "cell is singular or left-handed (det <= 0); no vacuum box is invented"
        )
    lattice = Lattice.from_matrix(cell)

    pbc = _pbc_from_ase(atoms)

    if 'dp_model' in atoms.info:
        # Path A: our own Atoms object — restore exactly.
        if model is not None and model != atoms.info['dp_model']:
            raise InteropError(
                f"model argument {model!r} does not match the stored "
                f"dp_model {atoms.info['dp_model']!r}"
            )
        model_val = atoms.info['dp_model']
        species = np.array(atoms.arrays['dp_species'], dtype=object)
        site_label = np.array(atoms.arrays['dp_site_label'], dtype=object)
        atom_id = np.array(atoms.arrays['dp_atom_id'], dtype=np.int64)
        occupancy = np.array(atoms.arrays['dp_occupancy'], dtype=float)
        cell_index = np.array(atoms.arrays['dp_cell_index'], dtype=np.int64)
        if 'dp_b_iso' in atoms.arrays:
            b_iso = np.array(atoms.arrays['dp_b_iso'], dtype=float)
        else:
            b_iso = None
        supercell = tuple(int(x) for x in atoms.info['dp_supercell'])
        frac = atoms.get_scaled_positions()

        return Structure(
            lattice,
            species,
            frac,
            model=model_val,
            occupancy=occupancy,
            site_label=site_label,
            atom_id=atom_id,
            b_iso=b_iso,
            cell_index=cell_index,
            supercell=supercell,
            pbc=pbc,
            min_separation=min_separation,
        )

    # Path B: foreign Atoms object.
    if model is None:
        raise InteropError(
            "foreign Atoms object: pass model='configuration' or 'average'"
        )

    n = len(atoms)
    species = np.array(atoms.get_chemical_symbols(), dtype=object)
    occupancy = np.ones(n, dtype=float)
    atom_id = np.arange(n, dtype=np.int64)
    cell_index = np.zeros((n, 3), dtype=np.int64)
    frac = atoms.get_scaled_positions()

    # Losses (design rule 4).
    losses = LossReport()

    # info['occupancy']: format is {str(listed_site_index): {symbol: occ}}.
    # Any occupancy value < 1 is a loss (the atom→key mapping is not
    # recoverable from Atoms; see INTEROP_FACTS P10).
    if 'occupancy' in atoms.info:
        occ_info = atoms.info['occupancy']
        for val in occ_info.values():
            if isinstance(val, dict):
                occs = list(val.values())
            else:
                occs = [val]
            if any(occ < 1.0 - 1e-9 for occ in occs):
                losses.add(
                    "partial occupancy in info['occupancy'] cannot be mapped "
                    "to atoms (use read_cif for such CIFs)"
                )
                break

    if atoms.constraints:
        losses.add("constraints not supported")

    # These attributes do not exist on a fresh Atoms object; use getattr.
    momenta = getattr(atoms, "momenta", None)
    if momenta is not None and np.any(np.asarray(momenta) != 0):
        losses.add("non-zero momenta not supported")

    init_mm = getattr(atoms, "init_magnetic_moments", None)
    if init_mm is not None and np.any(np.asarray(init_mm) != 0):
        losses.add("non-zero initial magnetic moments not supported")

    init_ch = getattr(atoms, "init_charges", None)
    if init_ch is not None and np.any(np.asarray(init_ch) != 0):
        losses.add("non-zero initial charges not supported")

    tags = getattr(atoms, "tags", None)
    if tags is not None and np.any(np.asarray(tags) != 0):
        losses.add("non-zero tags not supported")

    if atoms.calc is not None:
        losses.add("a calculator is attached; not supported")

    # info keys other than 'occupancy' (no dp_* keys exist in path B).
    for key in atoms.info:
        if key != 'occupancy':
            losses.add(f"info key {key!r} not supported")

    losses.raise_or_warn(allow_loss)

    return Structure(
        lattice,
        species,
        frac,
        model=model,
        occupancy=occupancy,
        atom_id=atom_id,
        cell_index=cell_index,
        supercell=(1, 1, 1),
        pbc=pbc,
        min_separation=min_separation,
    )
