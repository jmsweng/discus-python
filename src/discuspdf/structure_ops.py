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


def substitute(
    structure: Structure,
    *,
    old: str,
    new: str,
    seed: int,
    fraction: float | None = None,
    count: int | None = None,
    site: str | None = None,
    new_b_iso: float | None = None,
) -> EditResult:
    """Substitute some atoms of species ``old`` with species ``new``.

    Only ``model="configuration"``.  Exactly one of ``fraction`` / ``count``.
    """
    if structure.model != "configuration":
        raise ValueError("substitute requires model='configuration'; use sample_configuration first")

    if (fraction is None) == (count is None):
        raise ValueError("exactly one of fraction or count must be given")
    if fraction is not None:
        if not (0.0 <= fraction <= 1.0):
            raise ValueError(f"fraction must be in [0, 1], got {fraction}")
    if count is not None:
        if count < 0:
            raise ValueError(f"count must be >= 0, got {count}")

    # Validate new species (like Structure does).
    if not isinstance(new, str) or not new:
        raise ValueError(f"new species must be a non-empty string, got {new!r}")
    if len(new) > 8:
        raise ValueError(f"new species exceeds 8 characters: {new!r}")
    if any(c.isspace() for c in new):
        raise ValueError(f"new species contains whitespace: {new!r}")

    # Candidates.
    mask = structure.species == old
    if site is not None:
        mask &= structure.site_label == site
    candidate_rows = np.where(mask)[0]
    if len(candidate_rows) == 0:
        raise ValueError(f"no atoms with species={old!r}" + (f" and site={site!r}" if site else ""))

    # Number to change.
    if count is not None:
        n = count
    else:
        n = int(np.floor(fraction * len(candidate_rows) + 0.5))
    if n > len(candidate_rows):
        raise ValueError(f"n={n} exceeds number of candidates ({len(candidate_rows)})")

    # Choose rows.
    chosen = np.sort(np.random.default_rng(seed).choice(candidate_rows, size=n, replace=False))

    # Build new arrays.
    species_new = structure.species.copy()
    species_new[chosen] = new

    b_iso_new = None
    if new_b_iso is not None:
        if structure.b_iso is None:
            raise ValueError("new_b_iso requires the structure to have b_iso")
        b_iso_new = structure.b_iso.copy()
        b_iso_new[chosen] = new_b_iso

    if b_iso_new is not None:
        new_struct = structure.replace(species=species_new, b_iso=b_iso_new)
    else:
        new_struct = structure.replace(species=species_new)

    ids = np.sort(structure.atom_id[chosen])

    return EditResult(
        structure=new_struct,
        ids=ids,
        seed=seed,
        info={
            "n_candidates": int(len(candidate_rows)),
            "n_changed": int(n),
            "old": old,
            "new": new,
        },
    )


def displace(
    structure: Structure,
    *,
    seed: int | None = None,
    sigma: float | None = None,
    vector=None,
    coords: str = "cartesian",
    ids: np.ndarray | None = None,
    species: str | list[str] | None = None,
    site: str | None = None,
    mask: np.ndarray | None = None,
    min_separation: float | None = 0.1,
) -> EditResult:
    """Displace selected atoms by a Gaussian (``sigma``) or a fixed vector.

    Only ``model="configuration"``.  Exactly one of ``sigma`` / ``vector``.
    """
    if structure.model != "configuration":
        raise ValueError("displace requires model='configuration'")

    if (sigma is None) == (vector is None):
        raise ValueError("exactly one of sigma or vector must be given")
    if sigma is not None:
        if sigma <= 0:
            raise ValueError(f"sigma must be > 0, got {sigma}")
        if seed is None:
            raise ValueError("sigma requires a seed")
    if vector is not None:
        if seed is not None:
            raise ValueError("vector must not use a seed (seed must be None)")
        if coords not in ("cartesian", "fractional"):
            raise ValueError(f"coords must be 'cartesian' or 'fractional', got {coords!r}")

    # Selection: at most one of ids / species / site / mask; none = all rows.
    n_sel_args = sum(x is not None for x in (ids, species, site, mask))
    if n_sel_args > 1:
        raise ValueError("at most one of ids, species, site, mask may be given")

    if ids is not None:
        sel = structure.select(ids=ids)
    elif species is not None:
        sel = structure.select(species=species)
    elif site is not None:
        sel = structure.select(site=site)
    elif mask is not None:
        sel = structure.select(mask=mask)
    else:
        sel = structure

    n_sel = sel.n_atoms
    if n_sel == 0:
        raise ValueError("empty selection")

    # Get the row indices of the selected rows in the original structure.
    # sel.atom_id gives the atom_ids; we need the row positions.
    # Build a map from atom_id -> row index.
    aid_to_row = {int(a): i for i, a in enumerate(structure.atom_id)}
    sel_rows = np.array([aid_to_row[int(a)] for a in sel.atom_id], dtype=np.int64)

    # Displacement.
    frac_new = structure.frac.copy()
    if sigma is not None:
        noise = np.random.default_rng(seed).normal(0, sigma, size=(n_sel, 3))
        # Convert Cartesian noise to fractional: df = dcart @ inv(A)
        A = structure.lattice.matrix
        Ainv = np.linalg.inv(A)
        df = noise @ Ainv.T
        frac_new[sel_rows] += df
    else:
        vec = np.asarray(vector, dtype=float)
        if vec.shape != (3,):
            raise ValueError(f"vector must have 3 components, got shape {vec.shape}")
        if coords == "cartesian":
            A = structure.lattice.matrix
            Ainv = np.linalg.inv(A)
            df = vec @ Ainv.T
        else:  # fractional
            df = vec
        frac_new[sel_rows] += df

    new_struct = structure.replace(frac=frac_new, min_separation=min_separation)

    return EditResult(
        structure=new_struct,
        ids=np.sort(structure.atom_id[sel_rows]),
        seed=seed,
        info={
            "n_displaced": int(n_sel),
            "sigma": sigma,
            "vector": None if vector is None else tuple(vector),
        },
    )


def check_composition(structure: Structure, expected: dict, *, atol: float = 0.0, rtol: float = 0.0) -> None:
    """Check that the structure's composition matches ``expected``.

    Compares ``structure.counts()`` for ``"configuration"`` and
    ``structure.composition()`` for ``"average"``.  Missing/extra species
    count as 0.  Raises ``ValueError`` listing every mismatch.
    """
    if structure.model == "configuration":
        actual = structure.counts()
    else:
        actual = structure.composition()

    all_species = set(expected) | set(actual)
    mismatches = []
    for sp in sorted(all_species):
        exp = float(expected.get(sp, 0.0))
        act = float(actual.get(sp, 0.0))
        if abs(act - exp) > atol + rtol * abs(exp):
            mismatches.append((sp, exp, act))

    if mismatches:
        details = "; ".join(f"{sp}: expected {exp}, got {act}" for sp, exp, act in mismatches)
        raise ValueError(f"composition mismatch: {details}")


def sample_configuration(structure: Structure, *, seed: int) -> EditResult:
    """Draw one realised configuration from an average (site-occupancy) structure.

    Each site's occupancies are treated as a categorical distribution over its
    member rows (in row order); the last member is realised when the total
    occupancy is 1 (within 1e-9), so the site may be a vacancy when the total
    occupancy is < 1.  Draws use ``np.random.default_rng(seed).random(n_sites)``
    once, in site order of first appearance.

    Args:
        structure: a ``model="average"`` ``Structure``.
        seed: integer seed for the random draw.

    Returns:
        An ``EditResult`` with ``model="configuration"`` (occupancies all 1.0),
        the realised rows in their original order, and ``ids`` = the sorted
        ``atom_id``s of removed rows (vacancies and rejected alternatives).

    Raises:
        ValueError: if ``structure.model`` is not ``"average"``, or if the
            result is empty (all sites vacant).
    """
    if structure.model == "configuration":
        raise ValueError(
            "sample_configuration requires model='average'; "
            "a configuration structure is already sampled"
        )
    if structure.model != "average":
        raise ValueError(
            f"sample_configuration requires model='average', got {structure.model!r}"
        )

    # Sites in order of first appearance (smallest member row index).
    sites = sorted(structure.sites(), key=lambda s: min(s["members"]))

    # One draw per site, in that site order.
    u = np.random.default_rng(seed).random(len(sites))

    keep_rows: list[int] = []
    n_vacancies = 0
    for k, site in enumerate(sites):
        members = site["members"]  # in row order
        cum = np.cumsum(structure.occupancy[members])
        if cum[-1] >= 1.0 - 1e-9:
            cum = cum.copy()
            cum[-1] = 1.0
        j = int(np.searchsorted(cum, u[k], side="right"))
        if j < len(members):
            keep_rows.append(members[j])
        else:
            n_vacancies += 1

    if not keep_rows:
        raise ValueError("sample_configuration produced an empty structure (all sites vacant)")

    keep_rows = sorted(keep_rows)

    keep_set = set(keep_rows)
    removed_rows = [i for i in range(structure.n_atoms) if i not in keep_set]

    n = len(keep_rows)
    species_new = structure.species[keep_rows]
    frac_new = structure.frac[keep_rows]
    occ_new = np.ones(n, dtype=float)
    sl_new = structure.site_label[keep_rows]
    aid_new = structure.atom_id[keep_rows]
    ci_new = structure.cell_index[keep_rows]
    b_iso_new = None
    if structure.b_iso is not None:
        b_iso_new = structure.b_iso[keep_rows]

    new_struct = Structure(
        structure.lattice,
        species_new,
        frac_new,
        model="configuration",
        occupancy=occ_new,
        site_label=sl_new,
        atom_id=aid_new,
        b_iso=b_iso_new,
        cell_index=ci_new,
        supercell=structure.supercell,
        pbc=structure.pbc,
    )

    # expected: occupancy-weighted count per species (from the input)
    expected = structure.composition()
    # realised: raw count per species in the output
    realised = new_struct.counts()

    return EditResult(
        structure=new_struct,
        ids=np.sort(structure.atom_id[removed_rows]),
        seed=seed,
        info={
            "n_sites": int(len(sites)),
            "n_vacancies": int(n_vacancies),
            "expected": expected,
            "realised": realised,
        },
    )
