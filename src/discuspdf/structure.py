"""Immutable atomic-structure value type (Package A, Feature A1).

A ``Structure`` holds per-atom NumPy arrays (no per-atom Python objects) plus a
``Lattice``.  It supports two models:

* ``"average"`` — a site-occupancy model (as in a CIF).  Mixed sites are several
  rows with the same ``(site_label, cell_index)``; occupancies sum to ≤ 1.
* ``"configuration"`` — a realised atomistic arrangement (as in a DISCUS
  supercell).  Every occupancy is exactly 1.0; each site holds at most one atom.

All stored arrays are read-only.  ``__setattr__`` and ``__delattr__`` raise.
Functional edits return new ``Structure`` instances; the original is untouched.
"""
from __future__ import annotations

import numpy as np

from .lattice import Lattice
from ._elements import atomic_number


def _validate_species(species: np.ndarray) -> None:
    """Validate a (N,) string array of species labels."""
    if species.ndim != 1:
        raise ValueError(f"species must be a 1-D array, got shape {species.shape}")
    for i, s in enumerate(species):
        if not isinstance(s, str) or not s:
            raise ValueError(f"species[{i}] must be a non-empty string, got {s!r}")
        if len(s) > 8:
            raise ValueError(f"species[{i}] exceeds 8 characters: {s!r}")
        if any(c.isspace() for c in s):
            raise ValueError(f"species[{i}] contains whitespace: {s!r}")


def _validate_frac(frac: np.ndarray) -> None:
    """Validate a (N, 3) float array of fractional coordinates."""
    if frac.ndim != 2 or frac.shape[1] != 3:
        raise ValueError(f"frac must have shape (N, 3), got {frac.shape}")
    if not np.isfinite(frac).all():
        raise ValueError("frac contains NaN or inf values")


def _validate_occupancy(occ: np.ndarray) -> None:
    """Validate a (N,) float array of occupancies in (0, 1]."""
    if occ.ndim != 1:
        raise ValueError(f"occupancy must be a 1-D array, got shape {occ.shape}")
    if np.any(occ <= 0) or np.any(occ > 1):
        raise ValueError("occupancy values must be in (0, 1]")


def _validate_cell_index(ci: np.ndarray, supercell: np.ndarray) -> None:
    """Validate a (N, 3) int array of unit-cell indices inside a supercell."""
    if ci.ndim != 2 or ci.shape[1] != 3:
        raise ValueError(f"cell_index must have shape (N, 3), got {ci.shape}")
    if np.any(ci < 0):
        raise ValueError("cell_index values must be >= 0")
    for axis in range(3):
        if np.any(ci[:, axis] >= supercell[axis]):
            raise ValueError(
                f"cell_index[:, {axis}] exceeds supercell[{axis}] = {supercell[axis]}"
            )


def _validate_b_iso(b_iso: np.ndarray | None) -> None:
    """Validate an optional (N,) float array of isotropic B factors (>= 0)."""
    if b_iso is None:
        return
    if b_iso.ndim != 1:
        raise ValueError(f"b_iso must be a 1-D array, got shape {b_iso.shape}")
    if np.any(b_iso < 0):
        raise ValueError("b_iso values must be >= 0")
    if not np.isfinite(b_iso).all():
        raise ValueError("b_iso contains NaN or inf values")


def _check_site_consistency(
    model: str,
    species: np.ndarray,
    frac: np.ndarray,
    occupancy: np.ndarray,
    site_label: np.ndarray,
    cell_index: np.ndarray,
    b_iso: np.ndarray | None,
) -> None:
    """Check site-level consistency rules for the given model."""
    n = len(species)
    if n == 0:
        return

    # Build a composite key for grouping.
    # For large n, use a tuple-based approach to avoid slow Python string formatting.
    if n > 1000:
        # Use a structured view: (site_label, ci0, ci1, ci2) as a composite key.
        # np.unique works on structured arrays.
        dtype = np.dtype([("sl", object), ("c0", np.int64), ("c1", np.int64), ("c2", np.int64)])
        keys = np.empty(n, dtype=dtype)
        keys["sl"] = site_label
        keys["c0"] = cell_index[:, 0]
        keys["c1"] = cell_index[:, 1]
        keys["c2"] = cell_index[:, 2]
    else:
        keys = np.array(
            [f"{site_label[i]}|{cell_index[i, 0]},{cell_index[i, 1]},{cell_index[i, 2]}"
             for i in range(n)],
            dtype=object,
        )
    unique_keys, inverse, counts = np.unique(keys, return_inverse=True, return_counts=True)

    # Only sites with > 1 member need checking.
    multi_mask = counts > 1
    if not multi_mask.any():
        return

    # Map from original row index to multi-site index (or -1 if not a multi-site).
    multi_site_ids = np.where(multi_mask[inverse], inverse, -1)
    n_multi = int(multi_mask.sum())

    for key_idx in range(n_multi):
        # Find the original unique_keys index for this multi-site.
        orig_idx = np.where(multi_mask)[0][key_idx]
        key = unique_keys[orig_idx]
        members = np.where(multi_site_ids == key_idx)[0]
        if model == "configuration":
            raise ValueError(
                f"configuration model: site {key!r} has "
                f"{len(members)} atoms; a configuration site holds at most one atom"
            )
        # average model: all rows must share the same fractional position (<= 1e-6)
        # and the same b_iso.
        ref_frac = frac[members[0]]
        if np.max(np.abs(frac[members] - ref_frac)) > 1e-6:
            raise ValueError(
                f"average model: site {key!r} has rows with "
                f"differing fractional positions (> 1e-6)"
            )
        if b_iso is not None:
            ref_b = b_iso[members[0]]
            if np.max(np.abs(b_iso[members] - ref_b)) > 1e-12:
                raise ValueError(
                    f"average model: site {key!r} has rows with "
                    f"differing b_iso values"
                )
        # Occupancies must sum to <= 1 + 1e-6.
        occ_sum = float(np.sum(occupancy[members]))
        if occ_sum > 1.0 + 1e-6:
            raise ValueError(
                f"average model: site {key!r} has total occupancy "
                f"{occ_sum:.6f} > 1 + 1e-6"
            )


def _check_overlap(
    model: str,
    frac: np.ndarray,
    occupancy: np.ndarray,
    site_label: np.ndarray,
    cell_index: np.ndarray,
    lattice: Lattice,
    pbc: tuple[bool, bool, bool],
    min_separation: float | None = 0.1,
    n_atoms_limit: int = 20_000,
) -> None:
    """Check for unphysical atom overlaps.

    In a configuration, every pair of atoms in different sites must be >= min_separation.
    In an average model, the check applies only when two occupancies sum to > 1.
    Skipped above n_atoms_limit unless min_separation is explicitly requested.
    """
    if min_separation is None:
        return
    n = len(frac)
    if n >= n_atoms_limit:
        return

    all_pbc = all(pbc)
    no_pbc = not any(pbc)
    if not all_pbc and not no_pbc:
        raise NotImplementedError("mixed PBC not yet supported")

    # Group rows by site.
    if n > 1000:
        dtype = np.dtype([("sl", object), ("c0", np.int64), ("c1", np.int64), ("c2", np.int64)])
        keys = np.empty(n, dtype=dtype)
        keys["sl"] = site_label
        keys["c0"] = cell_index[:, 0]
        keys["c1"] = cell_index[:, 1]
        keys["c2"] = cell_index[:, 2]
    else:
        keys = np.array(
            [f"{site_label[i]}|{cell_index[i, 0]},{cell_index[i, 1]},{cell_index[i, 2]}"
             for i in range(n)],
            dtype=object,
        )
    unique_sites, inverse = np.unique(keys, return_inverse=True)
    n_sites = len(unique_sites)

    # For each pair of distinct sites, check the minimum distance.
    for i in range(n_sites):
        rows_i = np.where(inverse == i)[0]
        for j in range(i + 1, n_sites):
            rows_j = np.where(inverse == j)[0]
            # Determine if this pair should be checked.
            if model == "configuration":
                check = True
            else:
                occ_i = float(np.sum(occupancy[rows_i]))
                occ_j = float(np.sum(occupancy[rows_j]))
                check = (occ_i + occ_j) > 1.0
            if not check:
                continue
            # Vectorised: compute all pairwise distances between rows_i and rows_j.
            df = frac[rows_i][:, None, :] - frac[rows_j][None, :, :]  # (ni, nj, 3)
            if all_pbc:
                # Process in chunks to keep memory bounded.
                chunk = max(1, 1000 // max(len(rows_j), 1))
                found = False
                for start in range(0, len(rows_i), chunk):
                    stop = min(start + chunk, len(rows_i))
                    dc = df[start:stop]
                    df_min, _ = lattice.minimum_image(dc)
                    dists = np.linalg.norm(df_min @ lattice.matrix, axis=2)  # (chunk, nj)
                    if np.any(dists < min_separation):
                        found = True
                        break
                if found:
                    # Find the specific pair for the error message.
                    dmin_idx = np.unravel_index(np.argmin(dists), dists.shape)
                    d = float(dists[dmin_idx])
                    raise ValueError(
                        f"atoms at sites {unique_sites[i]!r} and {unique_sites[j]!r} "
                        f"are {d:.6f} Å apart (< {min_separation} Å)"
                    )
            else:
                dists = np.linalg.norm(df @ lattice.matrix, axis=2)  # (ni, nj)
                if np.any(dists < min_separation):
                    dmin_idx = np.unravel_index(np.argmin(dists), dists.shape)
                    d = float(dists[dmin_idx])
                    raise ValueError(
                        f"atoms at sites {unique_sites[i]!r} and {unique_sites[j]!r} "
                        f"are {d:.6f} Å apart (< {min_separation} Å)"
                    )


def _Structure_bypass_validation(
    lattice: Lattice,
    species: np.ndarray,
    frac: np.ndarray,
    *,
    model: str,
    occupancy: np.ndarray,
    site_label: np.ndarray,
    atom_id: np.ndarray,
    b_iso: np.ndarray | None,
    cell_index: np.ndarray,
    supercell: tuple[int, int, int],
    pbc: tuple[bool, bool, bool],
) -> "Structure":
    """Internal: create a Structure without site-consistency or overlap validation.

    Used by ``concat`` and ``_subset`` where the caller has already validated
    the individual structures.
    """
    species = np.array(species, dtype=object, copy=True)
    frac = np.array(frac, dtype=float, copy=True)
    occupancy = np.array(occupancy, dtype=float, copy=True)
    site_label = np.array(site_label, dtype=object, copy=True)
    atom_id = np.array(atom_id, dtype=np.int64, copy=True)
    cell_index = np.array(cell_index, dtype=np.int64, copy=True)
    supercell_arr = np.array(supercell, dtype=np.int64, copy=True)

    if b_iso is not None:
        b_iso = np.array(b_iso, dtype=float, copy=True)

    species.flags.writeable = False
    frac.flags.writeable = False
    occupancy.flags.writeable = False
    site_label.flags.writeable = False
    atom_id.flags.writeable = False
    if b_iso is not None:
        b_iso.flags.writeable = False
    cell_index.flags.writeable = False
    supercell_arr.flags.writeable = False

    obj = object.__new__(Structure)
    object.__setattr__(obj, "_lattice", lattice)
    object.__setattr__(obj, "_species", species)
    object.__setattr__(obj, "_frac", frac)
    object.__setattr__(obj, "_occupancy", occupancy)
    object.__setattr__(obj, "_site_label", site_label)
    object.__setattr__(obj, "_atom_id", atom_id)
    object.__setattr__(obj, "_b_iso", b_iso)
    object.__setattr__(obj, "_cell_index", cell_index)
    object.__setattr__(obj, "_supercell", supercell_arr)
    object.__setattr__(obj, "_pbc", pbc)
    object.__setattr__(obj, "_model", model)
    return obj


class Structure:
    """Immutable atomic structure.

    Holds per-atom NumPy arrays (read-only) plus a ``Lattice``.  Supports
    ``"average"`` (site-occupancy) and ``"configuration"`` (realised) models.
    """

    __slots__ = (
        "_lattice", "_species", "_frac", "_occupancy", "_site_label",
        "_atom_id", "_b_iso", "_cell_index", "_supercell", "_pbc", "_model",
    )

    def __init__(
        self,
        lattice: Lattice,
        species: np.ndarray,
        frac: np.ndarray,
        *,
        model: str,
        occupancy: np.ndarray | None = None,
        site_label: np.ndarray | None = None,
        atom_id: np.ndarray | None = None,
        b_iso: np.ndarray | None = None,
        cell_index: np.ndarray | None = None,
        supercell: tuple[int, int, int] = (1, 1, 1),
        pbc: tuple[bool, bool, bool] = (True, True, True),
        min_separation: float | None = 0.1,
    ) -> None:
        if model not in ("average", "configuration"):
            raise ValueError(f"model must be 'average' or 'configuration', got {model!r}")

        species = np.array(species, dtype=object, copy=True)
        frac = np.array(frac, dtype=float, copy=True)
        n = len(species)

        _validate_species(species)
        _validate_frac(frac)

        if occupancy is None:
            occupancy = np.ones(n, dtype=float)
        else:
            occupancy = np.array(occupancy, dtype=float, copy=True)
        _validate_occupancy(occupancy)

        if site_label is None:
            # Default: species + 1-based index within that species.
            site_label = np.empty(n, dtype=object)
            counts: dict[str, int] = {}
            for i in range(n):
                s = species[i]
                counts[s] = counts.get(s, 0) + 1
                site_label[i] = f"{s}{counts[s]}"
        else:
            site_label = np.array(site_label, dtype=object, copy=True)
            if site_label.shape != (n,):
                raise ValueError(f"site_label must have shape ({n},), got {site_label.shape}")

        if atom_id is None:
            atom_id = np.arange(n, dtype=np.int64)
        else:
            atom_id = np.array(atom_id, dtype=np.int64, copy=True)
            if atom_id.shape != (n,):
                raise ValueError(f"atom_id must have shape ({n},), got {atom_id.shape}")
            if np.any(atom_id < 0):
                raise ValueError("atom_id values must be >= 0")
            if len(np.unique(atom_id)) != n:
                raise ValueError("atom_id values must be unique")

        if b_iso is not None:
            b_iso = np.array(b_iso, dtype=float, copy=True)
        _validate_b_iso(b_iso)

        if cell_index is None:
            cell_index = np.zeros((n, 3), dtype=np.int64)
        else:
            cell_index = np.array(cell_index, dtype=np.int64, copy=True)
        supercell_arr = np.array(supercell, dtype=np.int64, copy=True)
        _validate_cell_index(cell_index, supercell_arr)

        if supercell_arr.shape != (3,) or np.any(supercell_arr < 1):
            raise ValueError(f"supercell must be (3,) ints >= 1, got {supercell_arr}")
        supercell = supercell_arr

        pbc = tuple(pbc)
        if len(pbc) != 3:
            raise ValueError(f"pbc must be 3 booleans, got {pbc}")

        # Model-specific validation.
        if model == "configuration":
            if np.any(occupancy != 1.0):
                raise ValueError(
                    "configuration model: all occupancies must be exactly 1.0; "
                    "for partial occupancies use model='average'"
                )

        # Site consistency check.
        _check_site_consistency(
            model, species, frac, occupancy, site_label, cell_index, b_iso
        )

        # Overlap check.
        _check_overlap(
            model, frac, occupancy, site_label, cell_index,
            lattice, pbc, min_separation=min_separation,
        )

        # Store read-only.
        species.flags.writeable = False
        frac.flags.writeable = False
        occupancy.flags.writeable = False
        site_label.flags.writeable = False
        atom_id.flags.writeable = False
        if b_iso is not None:
            b_iso.flags.writeable = False
        cell_index.flags.writeable = False
        supercell.flags.writeable = False

        object.__setattr__(self, "_lattice", lattice)
        object.__setattr__(self, "_species", species)
        object.__setattr__(self, "_frac", frac)
        object.__setattr__(self, "_occupancy", occupancy)
        object.__setattr__(self, "_site_label", site_label)
        object.__setattr__(self, "_atom_id", atom_id)
        object.__setattr__(self, "_b_iso", b_iso)
        object.__setattr__(self, "_cell_index", cell_index)
        object.__setattr__(self, "_supercell", supercell)
        object.__setattr__(self, "_pbc", pbc)
        object.__setattr__(self, "_model", model)

    @classmethod
    def from_cartesian(
        cls,
        lattice: Lattice,
        species: np.ndarray,
        cart: np.ndarray,
        *,
        model: str,
        **kwargs,
    ) -> "Structure":
        """Build a Structure from Cartesian coordinates (Å).

        Converts to fractional via ``lattice.cart_to_frac`` and delegates to
        :meth:`__init__`.
        """
        cart = np.asarray(cart, dtype=float)
        if cart.ndim != 2 or cart.shape[1] != 3:
            raise ValueError(f"cart must have shape (N, 3), got {cart.shape}")
        if not np.isfinite(cart).all():
            raise ValueError("cart contains NaN or inf values")
        frac = lattice.cart_to_frac(cart)
        return cls(lattice, species, frac, model=model, **kwargs)

    def __setattr__(self, name, value) -> None:
        raise AttributeError(f"Structure is immutable; cannot set {name!r}")

    def __delattr__(self, name) -> None:
        raise AttributeError(f"Structure is immutable; cannot delete {name!r}")

    # -- properties -------------------------------------------------------
    @property
    def n_atoms(self) -> int:
        """Number of atom rows."""
        return len(self._species)

    @property
    def lattice(self) -> Lattice:
        """The lattice (frame for fractional coordinates)."""
        return self._lattice

    @property
    def species(self) -> np.ndarray:
        """(N,) string array of species labels. Read-only view."""
        return self._species

    @property
    def frac(self) -> np.ndarray:
        """(N, 3) float array of fractional coordinates. Read-only view."""
        return self._frac

    @property
    def cartesian(self) -> np.ndarray:
        """(N, 3) float array of Cartesian coordinates in Å. Returns a copy."""
        return self._frac @ self._lattice.matrix

    @property
    def occupancy(self) -> np.ndarray:
        """(N,) float array of occupancies. Read-only view."""
        return self._occupancy

    @property
    def site_label(self) -> np.ndarray:
        """(N,) string array of site labels. Read-only view."""
        return self._site_label

    @property
    def atom_id(self) -> np.ndarray:
        """(N,) int64 array of atom IDs. Read-only view."""
        return self._atom_id

    @property
    def b_iso(self) -> np.ndarray | None:
        """(N,) float array of isotropic B factors (Å²) or None. Read-only view."""
        return self._b_iso

    @property
    def cell_index(self) -> np.ndarray:
        """(N, 3) int array of unit-cell indices. Read-only view."""
        return self._cell_index

    @property
    def supercell(self) -> np.ndarray:
        """(3,) int array of supercell repetitions. Read-only view."""
        return self._supercell

    @property
    def pbc(self) -> tuple[bool, bool, bool]:
        """3-tuple of periodic-axis booleans."""
        return self._pbc

    @property
    def model(self) -> str:
        """'average' or 'configuration'."""
        return self._model

    # -- summaries --------------------------------------------------------
    def sites(self) -> list[dict]:
        """Table of unique sites with members and total occupancy.

        Returns a list of dicts, one per unique (site_label, cell_index) pair:
        ``{"site_label": str, "cell_index": (i, j, k), "members": [row_idx, ...],
        "total_occupancy": float}``.
        """
        n = self.n_atoms
        if n == 0:
            return []
        keys = [
            f"{self._site_label[i]}|{self._cell_index[i, 0]},{self._cell_index[i, 1]},{self._cell_index[i, 2]}"
            for i in range(n)
        ]
        unique_keys, inverse = np.unique(np.array(keys, dtype=object), return_inverse=True)
        result = []
        for key_idx, key in enumerate(unique_keys):
            parts = key.split("|")
            site_label = parts[0]
            ci = tuple(int(x) for x in parts[1].split(","))
            members = list(np.where(inverse == key_idx)[0])
            total_occ = float(np.sum(self._occupancy[members]))
            result.append({
                "site_label": site_label,
                "cell_index": ci,
                "members": members,
                "total_occupancy": total_occ,
            })
        return result

    def atomic_numbers(self) -> np.ndarray:
        """(N,) int array of atomic numbers.

        Resolves each species label via ``_elements.atomic_number``; unknown
        labels raise ``ValueError`` here (not at construction).
        """
        return np.array([atomic_number(s) for s in self._species], dtype=np.int64)

    def composition(self) -> dict[str, float]:
        """Dict species -> occupancy-weighted count."""
        comp: dict[str, float] = {}
        for i in range(self.n_atoms):
            s = self._species[i]
            comp[s] = comp.get(s, 0.0) + float(self._occupancy[i])
        return comp

    def counts(self) -> dict[str, int]:
        """Dict species -> number of rows."""
        c: dict[str, int] = {}
        for s in self._species:
            c[s] = c.get(s, 0) + 1
        return c

    def formula(self) -> str:
        """Chemical formula string, e.g. 'GaAs' or 'Cu0.75Au0.25'."""
        comp = self.composition()
        parts = []
        for s, count in comp.items():
            if count == 1.0:
                parts.append(s)
            else:
                parts.append(f"{s}{count:.2f}".rstrip("0").rstrip("."))
        return "".join(parts)

    # -- selection --------------------------------------------------------
    def select(
        self,
        mask: np.ndarray | None = None,
        *,
        ids: np.ndarray | None = None,
        species: str | list[str] | None = None,
        site: str | None = None,
        allow_empty: bool = False,
    ) -> "Structure":
        """Return a new Structure with the selected rows (original atom_ids, original order).

        An empty selection raises ``ValueError`` unless ``allow_empty=True``.
        """
        n = self.n_atoms
        if mask is not None:
            mask = np.asarray(mask, dtype=bool)
            if mask.shape != (n,):
                raise ValueError(f"mask must have shape ({n},), got {mask.shape}")
            selected = mask
        elif ids is not None:
            ids = np.asarray(ids, dtype=np.int64)
            if np.any(ids < 0) or np.any(ids >= n):
                raise ValueError(f"ids must be in [0, {n})")
            selected = np.zeros(n, dtype=bool)
            selected[ids] = True
        elif species is not None:
            if isinstance(species, str):
                species = [species]
            selected = np.isin(self._species, species)
        elif site is not None:
            selected = self._site_label == site
        else:
            selected = np.ones(n, dtype=bool)

        if not selected.any() and not allow_empty:
            raise ValueError("empty selection; pass allow_empty=True to permit it")

        idx = np.where(selected)[0]
        return self._subset(idx)

    def _subset(self, idx: np.ndarray) -> "Structure":
        """Internal: create a new Structure from the given row indices."""
        b_iso = self._b_iso[idx].copy() if self._b_iso is not None else None
        return _Structure_bypass_validation(
            lattice=self._lattice,
            species=self._species[idx],
            frac=self._frac[idx],
            model=self._model,
            occupancy=self._occupancy[idx],
            site_label=self._site_label[idx],
            atom_id=self._atom_id[idx],
            b_iso=b_iso,
            cell_index=self._cell_index[idx],
            supercell=tuple(self._supercell),
            pbc=self._pbc,
        )

    # -- functional edits -------------------------------------------------
    def replace(self, **fields) -> "Structure":
        """Return a new Structure with the given fields replaced.

        Supported fields: species, frac, occupancy, site_label, atom_id,
        b_iso, cell_index, supercell, pbc.  The new values are re-validated.
        """
        if not fields:
            return self
        kwargs = {
            "lattice": self._lattice,
            "species": self._species.copy(),
            "frac": self._frac.copy(),
            "model": self._model,
            "occupancy": self._occupancy.copy(),
            "site_label": self._site_label.copy(),
            "atom_id": self._atom_id.copy(),
            "cell_index": self._cell_index.copy(),
            "supercell": tuple(self._supercell),
            "pbc": self._pbc,
            "min_separation": None,
        }
        if self._b_iso is not None:
            kwargs["b_iso"] = self._b_iso.copy()
        for name, value in fields.items():
            if name not in ("species", "frac", "occupancy", "site_label", "atom_id",
                           "b_iso", "cell_index", "supercell", "pbc"):
                raise ValueError(f"unknown field {name!r} for replace()")
            kwargs[name] = value
        return Structure(**kwargs)

    def with_frac(self, frac: np.ndarray) -> "Structure":
        """Return a new Structure with the given fractional coordinates."""
        return self.replace(frac=np.asarray(frac, dtype=float))

    def with_cartesian(self, cart: np.ndarray) -> "Structure":
        """Return a new Structure with the given Cartesian coordinates (Å)."""
        cart = np.asarray(cart, dtype=float)
        frac = self._lattice.cart_to_frac(cart)
        return self.replace(frac=frac)

    def with_occupancy(self, values: np.ndarray) -> "Structure":
        """Return a new Structure with the given occupancies (re-validated)."""
        return self.replace(occupancy=np.asarray(values, dtype=float))

    def remove(self, ids: np.ndarray) -> "Structure":
        """Return a new Structure without the given atom IDs."""
        ids = np.asarray(ids, dtype=np.int64)
        mask = np.ones(self.n_atoms, dtype=bool)
        mask[np.isin(np.arange(self.n_atoms), ids)] = False
        # More efficient: find which rows have the given IDs.
        id_to_row = {int(aid): i for i, aid in enumerate(self._atom_id)}
        for aid in ids:
            if int(aid) in id_to_row:
                mask[id_to_row[int(aid)]] = False
        idx = np.where(mask)[0]
        if len(idx) == 0:
            raise ValueError("removing all atoms; pass allow_empty=True to permit it")
        return self._subset(idx)

    def concat(self, other: "Structure", *, renumber: bool = False) -> "Structure":
        """Return a new Structure that is the concatenation of self and other.

        The lattices must be isclose.  atom_id clashes raise unless
        ``renumber=True``, in which case new IDs start at max(existing) + 1.
        """
        if not self._lattice.isclose(other._lattice):
            raise ValueError("concat: lattices are not isclose")
        if self._model != other._model:
            raise ValueError(
                f"concat: models differ ({self._model!r} vs {other._model!r})"
            )
        n1, n2 = self.n_atoms, other.n_atoms
        species = np.concatenate([self._species, other._species])
        frac = np.concatenate([self._frac, other._frac])
        occ = np.concatenate([self._occupancy, other._occupancy])
        site = np.concatenate([self._site_label, other._site_label])
        ci = np.concatenate([self._cell_index, other._cell_index])
        sc = np.maximum(self._supercell, other._supercell)

        if renumber:
            max_id = max(int(self._atom_id.max()), int(other._atom_id.max()))
            id1 = self._atom_id.copy()
            # Remap other's IDs to start at max_id + 1 (preserving relative order).
            min_other = int(other._atom_id.min())
            id2 = other._atom_id - min_other + (max_id + 1)
            atom_id = np.concatenate([id1, id2])
        else:
            atom_id = np.concatenate([self._atom_id, other._atom_id])
            if len(np.unique(atom_id)) != n1 + n2:
                raise ValueError(
                    "concat: atom_id clash; pass renumber=True to reassign IDs"
                )

        b1 = self._b_iso
        b2 = other._b_iso
        if b1 is None and b2 is None:
            b_iso = None
        elif b1 is None:
            b_iso = np.concatenate([np.full(n1, np.nan), b2])
        elif b2 is None:
            b_iso = np.concatenate([b1, np.full(n2, np.nan)])
        else:
            b_iso = np.concatenate([b1, b2])

        # Note: concat does not re-validate site consistency; the caller is
        # responsible for ensuring the combined structure is valid.
        return _Structure_bypass_validation(
            lattice=self._lattice,
            species=species,
            frac=frac,
            model=self._model,
            occupancy=occ,
            site_label=site,
            atom_id=atom_id,
            b_iso=b_iso,
            cell_index=ci,
            supercell=tuple(sc),
            pbc=self._pbc,
        )

    def renumbered(self) -> "Structure":
        """Return a new Structure with atom_ids renumbered to 0..N-1.

        This is the **only** method that changes atom_ids.
        """
        new_ids = np.arange(self.n_atoms, dtype=np.int64)
        return self.replace(atom_id=new_ids)

    # -- distances --------------------------------------------------------
    def _check_pbc(self) -> None:
        all_pbc = all(self._pbc)
        no_pbc = not any(self._pbc)
        if not all_pbc and not no_pbc:
            raise NotImplementedError("mixed PBC not yet supported")

    def distance(self, id_a: int, id_b: int) -> float:
        """Distance in Å between two atoms (minimum image if periodic)."""
        self._check_pbc()
        i = self._row_for_id(id_a)
        j = self._row_for_id(id_b)
        df = self._frac[i] - self._frac[j]
        if all(self._pbc):
            df_min, _ = self._lattice.minimum_image(df)
            return float(np.linalg.norm(df_min @ self._lattice.matrix))
        return float(np.linalg.norm(df @ self._lattice.matrix))

    def distance_matrix(self, ids: np.ndarray | None = None) -> np.ndarray:
        """(M, M) distance matrix in Å.  Refuses M > 5000."""
        self._check_pbc()
        if ids is None:
            idx = np.arange(self.n_atoms)
        else:
            ids = np.asarray(ids, dtype=np.int64)
            idx = np.array([self._row_for_id(int(i)) for i in ids])
        m = len(idx)
        if m > 5000:
            raise ValueError(
                f"distance_matrix: M = {m} > 5000; use Package C neighbour lists "
                "for large structures"
            )
        cart = self._cartesian_subset(idx)
        if all(self._pbc):
            # Compute all pairwise minimum-image distances.
            out = np.empty((m, m))
            for i in range(m):
                for j in range(i + 1, m):
                    df = self._frac[idx[i]] - self._frac[idx[j]]
                    df_min, _ = self._lattice.minimum_image(df)
                    d = float(np.linalg.norm(df_min @ self._lattice.matrix))
                    out[i, j] = d
                    out[j, i] = d
            np.fill_diagonal(out, 0.0)
            return out
        # Non-periodic: plain Cartesian.
        diff = cart[:, None, :] - cart[None, :, :]
        return np.linalg.norm(diff, axis=2)

    def displacement(self, id_a: int, id_b: int) -> np.ndarray:
        """Cartesian minimum-image displacement vector from atom a to atom b (Å)."""
        self._check_pbc()
        i = self._row_for_id(id_a)
        j = self._row_for_id(id_b)
        df = self._frac[j] - self._frac[i]
        if all(self._pbc):
            df_min, _ = self._lattice.minimum_image(df)
            return df_min @ self._lattice.matrix
        return df @ self._lattice.matrix

    def _row_for_id(self, atom_id: int) -> int:
        """Find the row index for a given atom_id."""
        for i, aid in enumerate(self._atom_id):
            if int(aid) == atom_id:
                return i
        raise ValueError(f"atom_id {atom_id} not found")

    def _cartesian_subset(self, idx: np.ndarray) -> np.ndarray:
        return self._frac[idx] @ self._lattice.matrix

    # -- wrapping ---------------------------------------------------------
    def wrapped(self) -> "Structure":
        """Return a new Structure with fractional coordinates wrapped to [0, 1).

        Only valid on fully periodic structures.  ``cell_index`` is unchanged.
        Wrapping changes coordinates but not physics.
        """
        if not all(self._pbc):
            raise ValueError("wrapped() is only valid on fully periodic structures")
        frac_wrapped = self._frac % 1.0
        return self.replace(frac=frac_wrapped)

    # -- comparison -------------------------------------------------------
    def __eq__(self, other) -> bool:
        if not isinstance(other, Structure):
            return False
        return (
            self._lattice == other._lattice
            and np.array_equal(self._species, other._species)
            and np.array_equal(self._frac, other._frac)
            and np.array_equal(self._occupancy, other._occupancy)
            and np.array_equal(self._site_label, other._site_label)
            and np.array_equal(self._atom_id, other._atom_id)
            and np.array_equal(self._cell_index, other._cell_index)
            and np.array_equal(self._supercell, other._supercell)
            and self._pbc == other._pbc
            and self._model == other._model
            and ((self._b_iso is None and other._b_iso is None)
                 or (self._b_iso is not None and other._b_iso is not None
                     and np.array_equal(self._b_iso, other._b_iso)))
        )

    def isclose(self, other: "Structure", rtol: float = 1e-12, atol: float = 1e-15) -> bool:
        """Element-wise closeness of all numeric fields."""
        if not isinstance(other, Structure):
            raise TypeError("isclose requires another Structure")
        if self._model != other._model:
            return False
        if self.n_atoms != other.n_atoms:
            return False
        if not self._lattice.isclose(other._lattice, rtol=rtol, atol=atol):
            return False
        if not np.allclose(self._frac, other._frac, rtol=rtol, atol=atol):
            return False
        if not np.allclose(self._occupancy, other._occupancy, rtol=rtol, atol=atol):
            return False
        if self._b_iso is not None and other._b_iso is not None:
            if not np.allclose(self._b_iso, other._b_iso, rtol=rtol, atol=atol):
                return False
        return True

    __hash__ = None  # type: ignore

    def __repr__(self) -> str:
        a, b, c, alpha, beta, gamma = self._lattice.parameters
        return (
            f"Structure(n={self.n_atoms}, formula={self.formula()}, "
            f"model={self._model!r}, supercell={tuple(self._supercell)}, "
            f"a={a:.4g} Å, b={b:.4g} Å, c={c:.4g} Å, "
            f"alpha={alpha:.4g}°, beta={beta:.4g}°, gamma={gamma:.4g}°)"
        )
