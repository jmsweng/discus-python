"""Structure editing operations (Package A, Feature A4).

All operations are pure functions that return an ``EditResult``.  The
``Structure`` stays immutable; atom order and ``atom_id``s are never changed
except as stated for ``repeat``.
"""
from __future__ import annotations

import dataclasses

import numpy as np

from .lattice import Lattice
from .structure import Structure


@dataclasses.dataclass(frozen=True)
class EditResult:
    """Result of a structure editing operation.

    Attributes:
        structure: the new ``Structure``.
        ids: sorted ``np.ndarray`` of changed ``atom_id``s (int64).
        seed: the integer seed used, or ``None``.
        info: plain-python dict of operation-specific metadata.
    """
    structure: Structure
    ids: np.ndarray
    seed: int | None
    info: dict


def repeat(structure: Structure, reps, *, min_separation: float | None = 0.1) -> EditResult:
    """Repeat a structure into a supercell.

    Args:
        structure: the parent ``Structure``.
        reps: a 3-tuple of positive ints ``(nx, ny, nz)``.
        min_separation: overlap-check threshold in Å (``None`` to skip).

    Returns:
        An ``EditResult`` with the repeated structure.

    Raises:
        ValueError: if ``reps`` is not a 3-tuple of positive ints, or if the
            parent has ``pbc`` not all True.
    """
    # Validate reps.
    if isinstance(reps, (int, np.integer)):
        raise ValueError(f"reps must be a 3-tuple, got int {reps}")
    if len(reps) != 3:
        raise ValueError(f"reps must have length 3, got {len(reps)}")
    for i, r in enumerate(reps):
        if isinstance(r, bool) or not isinstance(r, (int, np.integer)):
            raise ValueError(f"reps[{i}] must be a positive int, got {r!r}")
        if r < 1:
            raise ValueError(f"reps[{i}] must be >= 1, got {r}")
    nx, ny, nz = (int(r) for r in reps)

    # PBC check.
    if not all(structure.pbc):
        raise ValueError(f"repeat requires all pbc True, got {structure.pbc}")

    # New lattice: scale rows by reps.
    A = structure.lattice.matrix
    A_new = A * np.array([nx, ny, nz], dtype=float)[:, None]
    lat_new = Lattice.from_matrix(A_new)

    # Old arrays.
    n_old = structure.n_atoms
    frac_old = structure.frac
    ci_old = structure.cell_index
    sc_old = structure.supercell
    aid_old = structure.atom_id

    # Stride for atom_id mapping.
    stride = int(aid_old.max()) + 1

    # Build new arrays.
    # Row order: oz outermost, oy, ox innermost (x fastest); old rows in old order inside each offset.
    n_new = n_old * nx * ny * nz
    frac_new = np.empty((n_new, 3), dtype=float)
    ci_new = np.empty((n_new, 3), dtype=np.int64)
    aid_new = np.empty(n_new, dtype=np.int64)
    species_new = np.empty(n_new, dtype=object)
    occ_new = np.empty(n_new, dtype=float)
    sl_new = np.empty(n_new, dtype=object)
    b_iso_new = None
    if structure.b_iso is not None:
        b_iso_new = np.empty(n_new, dtype=float)

    row = 0
    for oz in range(nz):
        for oy in range(ny):
            for ox in range(nx):
                o = np.array([ox, oy, oz], dtype=np.int64)
                offset_index = ox + nx * (oy + ny * oz)
                frac_new[row:row + n_old] = (frac_old + o) / np.array([nx, ny, nz], dtype=float)
                ci_new[row:row + n_old] = ci_old + sc_old * o
                aid_new[row:row + n_old] = offset_index * stride + aid_old
                species_new[row:row + n_old] = structure.species
                occ_new[row:row + n_old] = structure.occupancy
                sl_new[row:row + n_old] = structure.site_label
                if b_iso_new is not None:
                    b_iso_new[row:row + n_old] = structure.b_iso
                row += n_old

    sc_new = tuple(int(s) * r for s, r in zip(sc_old, (nx, ny, nz)))

    new_struct = Structure(
        lat_new,
        species_new,
        frac_new,
        model=structure.model,
        occupancy=occ_new,
        site_label=sl_new,
        atom_id=aid_new,
        b_iso=b_iso_new,
        cell_index=ci_new,
        supercell=sc_new,
        pbc=structure.pbc,
        min_separation=min_separation,
    )

    return EditResult(
        structure=new_struct,
        ids=np.array([], dtype=np.int64),
        seed=None,
        info={"reps": (nx, ny, nz), "stride": stride},
    )
