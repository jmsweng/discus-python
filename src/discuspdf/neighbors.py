"""Periodic neighbor lists on a general lattice (Package C, Feature 12).

A :class:`NeighborList` holds every **directed** pair ``(i -> j, shift)`` of
*nodes* whose distance is at most ``cutoff``.  A node is one crystal site
(``site_label`` + ``cell_index``) of the structure the list was built on; in a
``"configuration"`` structure every site holds exactly one atom, so nodes and
atoms coincide.  ``shift`` is the integer translation, in units of the
structure's own (super)cell, that is added to node *j*::

    vector = (frac[j] + shift - frac[i]) @ A          # Å, A = structure.lattice.matrix

Every unordered pair appears twice (``i -> j`` with ``shift`` and
``j -> i`` with ``-shift``); :meth:`NeighborList.unique` keeps one of each.
All periodic images within the cutoff are listed, so a cutoff larger than half
the cell is handled correctly (a node can even be its own neighbor through a
nonzero shift).  The full lattice metric is used throughout: there is no
per-axis rounding, so skewed (triclinic) cells are exact.

Conventions are written out in ``docs/NEIGHBORS.md``.
"""
from __future__ import annotations

import csv
import dataclasses
from pathlib import Path

import numpy as np

from .structure import Structure

__all__ = ["NeighborList", "neighbor_list", "Shell", "label_sites_by_slot"]

_POSITION_TOL = 1e-3  # Å: members of one site must coincide (Structure already enforces 1e-6 in fractions)


@dataclasses.dataclass(frozen=True)
class Shell:
    """One distance shell of a (selected) neighbor list.

    Attributes:
        index: 1-based shell number (1 = nearest).
        distance: mean pair distance in the shell (Å).
        min_distance, max_distance: extreme pair distances in the shell (Å).
        n_pairs: number of directed pairs in the shell.
        n_origins: number of distinct origin nodes that have a pair in this shell.
        coordination: ``n_pairs / n_origin_nodes`` where ``n_origin_nodes`` is
            the number of origin nodes of the **whole** selection (so origins
            with no neighbor in this shell count as zero).
    """
    index: int
    distance: float
    min_distance: float
    max_distance: float
    n_pairs: int
    n_origins: int
    coordination: float


def _site_nodes(structure: Structure, position_tol: float):
    """Unique sites of ``structure`` in order of their first member row.

    Returns ``(first_row, member_count, site_of_row)``.  Raises if the members of
    one site are more than ``position_tol`` Å apart (split sites are not nodes).
    """
    n = structure.n_atoms
    keys: dict[tuple, int] = {}
    site_of = np.empty(n, dtype=np.int64)
    first_row: list[int] = []
    sl = structure.site_label
    ci = structure.cell_index
    for r in range(n):
        key = (sl[r], int(ci[r, 0]), int(ci[r, 1]), int(ci[r, 2]))
        k = keys.get(key)
        if k is None:
            k = len(first_row)
            keys[key] = k
            first_row.append(r)
        site_of[r] = k
    first_row_arr = np.array(first_row, dtype=np.int64)
    counts = np.bincount(site_of, minlength=len(first_row_arr))
    multi = np.flatnonzero(counts > 1)
    if multi.size:
        frac = structure.frac
        A = structure.lattice.matrix
        rows = np.flatnonzero(np.isin(site_of, multi))
        d = frac[rows] - frac[first_row_arr[site_of[rows]]]
        if all(structure.pbc):
            d, _ = structure.lattice.minimum_image(d)
        dist = np.linalg.norm(d @ A, axis=1)
        bad = np.flatnonzero(dist > position_tol)
        if bad.size:
            r = rows[bad[0]]
            raise ValueError(
                f"site {sl[r]!r} in cell {tuple(int(x) for x in ci[r])} has members "
                f"{dist[bad[0]]:.3g} Å apart (> {position_tol} Å); split sites cannot be "
                "neighbor-list nodes: build the list on a configuration or an undisplaced reference"
            )
    return first_row_arr, counts, site_of


def _node_images(structure: Structure, first_row: np.ndarray) -> np.ndarray:
    """(n_nodes, 3) integer cell offset of each node's stored position, in parent cells.

    ``cell_index`` names the parent cell of a site, but the stored ``frac`` may sit a
    whole number of **supercells** away (``Structure.wrapped()``, DISCUS centred files,
    hand-made offsets).  For every site label a canonical parent-cell position ``r``
    is taken from its first node, moved by whole supercells to lie within half a
    supercell of the structure's first node (along an axis with one cell, the first
    node's position is used as is).  A node whose parent-cell position
    ``frac * supercell - cell_index`` differs from ``r`` by about ``k`` supercells gets
    the image ``k * supercell``.  Displacements smaller than half a supercell are
    therefore never mistaken for a different cell.
    """
    sc = structure.supercell.astype(float)
    par = structure.frac[first_row] * sc - structure.cell_index[first_row]
    labels = structure.site_label[first_row]
    _, first, inverse = np.unique(labels.astype(str), return_index=True, return_inverse=True)
    r = par[first]
    # bring every label's reference within half a supercell of the first node's position
    r = np.where(sc > 1, r - sc * np.rint((r - par[0]) / sc), r)
    k = np.rint((par - r[inverse.reshape(-1)]) / sc)
    return (k * sc).astype(np.int64)


def _offsets(m: np.ndarray) -> np.ndarray:
    ranges = [np.arange(-int(m[k]), int(m[k]) + 1) for k in range(3)]
    return np.stack(np.meshgrid(*ranges, indexing="ij"), axis=-1).reshape(-1, 3)


def _pairs_periodic(frac: np.ndarray, A: np.ndarray, B: np.ndarray, cutoff: float):
    """All directed (i, j, shift) with |(f_j + shift - f_i) A| <= cutoff, excluding i == j with shift 0.

    Cell-list search in fractional coordinates.  Bins are at least ``cutoff``
    thick (perpendicular spacing) where possible; the number of neighbouring
    bins searched per axis, ``m``, is chosen so that every image within the
    cutoff is reached (see docs/NEIGHBORS.md, "Completeness").
    """
    n = len(frac)
    d_planes = 1.0 / np.linalg.norm(B, axis=1)            # interplanar spacings of the cell (Å)
    nb = np.maximum(1, np.floor(d_planes / cutoff)).astype(np.int64)
    max_bins = max(1, 2 * n)
    if np.prod(nb.astype(float)) > max_bins:
        scale = (np.prod(nb.astype(float)) / max_bins) ** (1.0 / 3.0)
        nb = np.maximum(1, np.floor(nb / scale)).astype(np.int64)
    # |Δf_k| <= cutoff * |b_k| = cutoff / d_k for any pair within the cutoff; a fractional
    # extent x spans at most floor(x * nb) + 1 bin boundaries.
    m = (np.floor(cutoff / d_planes * nb * (1.0 + 1e-12) + 1e-12) + 1).astype(np.int64)

    fl = np.floor(frac)
    fw = frac - fl
    bins = np.minimum((fw * nb).astype(np.int64), nb - 1)
    flat = (bins[:, 0] * nb[1] + bins[:, 1]) * nb[2] + bins[:, 2]
    order = np.argsort(flat, kind="stable")
    sorted_flat = flat[order]
    cut2 = cutoff * cutoff * (1.0 + 1e-12) + 1e-24

    out_i, out_j, out_s = [], [], []
    for o in _offsets(m):
        tb = bins + o                                      # target bin, unwrapped
        S = np.floor_divide(tb, nb)                        # image of the wrapped crystal
        tbw = tb - S * nb
        tflat = (tbw[:, 0] * nb[1] + tbw[:, 1]) * nb[2] + tbw[:, 2]
        lo = np.searchsorted(sorted_flat, tflat, side="left")
        hi = np.searchsorted(sorted_flat, tflat, side="right")
        counts = hi - lo
        total = int(counts.sum())
        if total == 0:
            continue
        ii = np.repeat(np.arange(n), counts)
        starts = np.repeat(lo - np.concatenate(([0], np.cumsum(counts)[:-1])), counts)
        jj = order[starts + np.arange(total)]
        Si = S[ii]
        for c0 in range(0, total, 1 << 20):
            i_c, j_c, S_c = ii[c0:c0 + (1 << 20)], jj[c0:c0 + (1 << 20)], Si[c0:c0 + (1 << 20)]
            df = fw[j_c] + S_c - fw[i_c]
            r = df @ A
            keep = np.einsum("ij,ij->i", r, r) <= cut2
            keep &= ~((i_c == j_c) & np.all(S_c == 0, axis=1))
            if np.any(keep):
                i_k, j_k = i_c[keep], j_c[keep]
                # back to the original (unwrapped) coordinates: f_j + shift - f_i == fw_j + S - fw_i
                shift = S_c[keep] - fl[j_k].astype(np.int64) + fl[i_k].astype(np.int64)
                out_i.append(i_k)
                out_j.append(j_k)
                out_s.append(shift)
    if not out_i:
        return np.empty(0, np.int64), np.empty(0, np.int64), np.empty((0, 3), np.int64)
    return np.concatenate(out_i), np.concatenate(out_j), np.concatenate(out_s)


def _pairs_open(cart: np.ndarray, cutoff: float):
    """All directed (i, j) with i != j and |r_j - r_i| <= cutoff (no periodic images)."""
    n = len(cart)
    lo = cart.min(axis=0)
    span = cart.max(axis=0) - lo
    nb = np.maximum(1, np.floor(span / cutoff)).astype(np.int64)
    max_bins = max(1, 2 * n)
    if np.prod(nb.astype(float)) > max_bins:
        scale = (np.prod(nb.astype(float)) / max_bins) ** (1.0 / 3.0)
        nb = np.maximum(1, np.floor(nb / scale)).astype(np.int64)
    width = np.where(span > 0, span / nb, 1.0)
    width = np.maximum(width, cutoff)                      # bins at least one cutoff wide
    bins = np.minimum(((cart - lo) / width).astype(np.int64), nb - 1)
    flat = (bins[:, 0] * nb[1] + bins[:, 1]) * nb[2] + bins[:, 2]
    order = np.argsort(flat, kind="stable")
    sorted_flat = flat[order]
    cut2 = cutoff * cutoff * (1.0 + 1e-12) + 1e-24
    out_i, out_j = [], []
    for o in _offsets(np.ones(3, dtype=np.int64)):
        tb = bins + o
        valid = np.all((tb >= 0) & (tb < nb), axis=1)
        tflat = (tb[:, 0] * nb[1] + tb[:, 1]) * nb[2] + tb[:, 2]
        lo_i = np.searchsorted(sorted_flat, tflat, side="left")
        hi_i = np.searchsorted(sorted_flat, tflat, side="right")
        counts = np.where(valid, hi_i - lo_i, 0)
        total = int(counts.sum())
        if total == 0:
            continue
        ii = np.repeat(np.arange(n), counts)
        starts = np.repeat(lo_i - np.concatenate(([0], np.cumsum(counts)[:-1])), counts)
        jj = order[starts + np.arange(total)]
        r = cart[jj] - cart[ii]
        keep = (np.einsum("ij,ij->i", r, r) <= cut2) & (ii != jj)
        out_i.append(ii[keep])
        out_j.append(jj[keep])
    if not out_i:
        return np.empty(0, np.int64), np.empty(0, np.int64)
    return np.concatenate(out_i), np.concatenate(out_j)


def label_sites_by_slot(structure: Structure, prefix: str = "site") -> Structure:
    """Return ``structure`` with ``site_label`` = ``f"{prefix}{k}"``, ``k`` = 1-based slot in the cell.

    DISCUS numbers the sites of a cell by their order (``set vec, n, is1, is2, …``);
    ``read_stru`` instead labels atoms by species and index, which makes every atom
    its own site.  This relabels a structure stored in DISCUS order (cell-major, the
    same number ``m`` of atoms in every cell, as ``write_stru`` writes and
    ``read_stru`` reads) so that row ``r`` gets slot ``r % m + 1``.  Use it on both
    the reference and the configurations before building neighbor lists or
    computing correlations.
    """
    n_cells = int(np.prod(structure.supercell))
    if structure.n_atoms % n_cells:
        raise ValueError(f"{structure.n_atoms} atoms are not a whole number per cell ({n_cells} cells)")
    m = structure.n_atoms // n_cells
    expected = np.repeat(np.arange(n_cells), m)
    sc = structure.supercell
    ci = structure.cell_index
    cell = (ci[:, 2] * sc[1] + ci[:, 1]) * sc[0] + ci[:, 0]
    if not np.array_equal(cell, expected):
        raise ValueError("rows are not in DISCUS cell-major order (x fastest, m atoms per cell)")
    labels = np.array([f"{prefix}{r % m + 1}" for r in range(structure.n_atoms)], dtype=object)
    return structure.replace(site_label=labels)


def neighbor_list(structure: Structure, cutoff: float, *, position_tol: float = _POSITION_TOL) -> "NeighborList":
    """Build the directed neighbor list of ``structure`` up to ``cutoff`` Å (inclusive).

    Args:
        structure: a :class:`~discuspdf.Structure` with ``pbc`` all True (periodic
            images are listed) or all False (no images).  Mixed PBC raises
            ``NotImplementedError``.
        cutoff: maximum center-to-center distance in Å (> 0).  Pairs at exactly
            ``cutoff`` are included (relative tolerance 1e-12).
        position_tol: members of one site (average model) must coincide within
            this many Å; split sites raise ``ValueError``.

    Returns:
        A :class:`NeighborList`; pairs sorted by origin node, then distance,
        then target node, then shift.
    """
    if not isinstance(structure, Structure):
        raise TypeError(f"structure must be a Structure, got {type(structure).__name__}")
    cutoff = float(cutoff)
    if not np.isfinite(cutoff) or cutoff <= 0:
        raise ValueError(f"cutoff must be a positive number of Å, got {cutoff!r}")
    pbc = tuple(bool(p) for p in structure.pbc)
    if any(pbc) and not all(pbc):
        raise NotImplementedError("mixed PBC not yet supported; use all True or all False")
    if structure.n_atoms == 0:
        raise ValueError("structure has no atoms")

    first_row, member_count, _ = _site_nodes(structure, position_tol)
    frac = structure.frac[first_row]
    A = structure.lattice.matrix
    if all(pbc):
        i, j, shift = _pairs_periodic(frac, A, structure.lattice.reciprocal_matrix, cutoff)
    else:
        i, j = _pairs_open(frac @ A, cutoff)
        shift = np.zeros((len(i), 3), dtype=np.int64)
    vec = (frac[j] + shift - frac[i]) @ A
    dist = np.linalg.norm(vec, axis=1)
    order = np.lexsort((shift[:, 2], shift[:, 1], shift[:, 0], j, np.round(dist, 9), i))
    return NeighborList(
        structure=structure,
        cutoff=cutoff,
        node_row=first_row,
        node_image=_node_images(structure, first_row),
        node_members=member_count,
        i=i[order],
        j=j[order],
        shift=shift[order],
        vector=vec[order],
        distance=dist[order],
    )


def _row_keys(rows: np.ndarray) -> np.ndarray:
    """One scalar key per row of an integer matrix (for set membership of whole rows)."""
    rows = np.ascontiguousarray(rows, dtype=np.int64)
    return rows.view(np.dtype((np.void, rows.dtype.itemsize * rows.shape[1]))).reshape(-1)


def _readonly(a: np.ndarray) -> np.ndarray:
    a = np.array(a, copy=True)
    a.flags.writeable = False
    return a


def _as_set(value, name: str):
    if value is None:
        return None
    if isinstance(value, (str, np.str_)):
        return {str(value)}
    try:
        values = {str(v) for v in value}
    except TypeError:
        raise TypeError(f"{name} must be a string or an iterable of strings, got {value!r}") from None
    if not values:
        raise ValueError(f"{name} must not be empty")
    return values


class NeighborList:
    """Directed periodic neighbor pairs between the sites of a structure.

    Build with :func:`neighbor_list` (or ``Structure.neighbors``).  Instances are
    immutable; :meth:`select` and :meth:`unique` return new lists that keep the
    same node table, so node indices stay comparable between selections.

    Per-node arrays (``n_nodes``): :attr:`node_row`, :attr:`node_site`,
    :attr:`node_cell_index`, :attr:`node_frac`, :attr:`node_atom_id`,
    :attr:`node_species`.  Per-pair arrays (``n_pairs``): :attr:`i`, :attr:`j`,
    :attr:`shift`, :attr:`vector`, :attr:`distance`, :attr:`translation`,
    :attr:`cell_vector`, and the origin/target views of the node arrays.
    """

    __slots__ = ("_structure", "_cutoff", "_node_row", "_node_image", "_node_members", "_i", "_j",
                 "_shift", "_vector", "_distance")

    def __init__(self, *, structure, cutoff, node_row, node_image, node_members, i, j, shift, vector,
                 distance):
        object.__setattr__(self, "_structure", structure)
        object.__setattr__(self, "_cutoff", float(cutoff))
        object.__setattr__(self, "_node_row", _readonly(np.asarray(node_row, dtype=np.int64)))
        object.__setattr__(self, "_node_image", _readonly(np.asarray(node_image, dtype=np.int64).reshape(-1, 3)))
        object.__setattr__(self, "_node_members", _readonly(np.asarray(node_members, dtype=np.int64)))
        object.__setattr__(self, "_i", _readonly(np.asarray(i, dtype=np.int64)))
        object.__setattr__(self, "_j", _readonly(np.asarray(j, dtype=np.int64)))
        object.__setattr__(self, "_shift", _readonly(np.asarray(shift, dtype=np.int64).reshape(-1, 3)))
        object.__setattr__(self, "_vector", _readonly(np.asarray(vector, dtype=float).reshape(-1, 3)))
        object.__setattr__(self, "_distance", _readonly(np.asarray(distance, dtype=float)))

    def __setattr__(self, name, value) -> None:
        raise AttributeError(f"NeighborList is immutable; cannot set {name!r}")

    def _subset(self, mask) -> "NeighborList":
        mask = np.asarray(mask)
        return NeighborList(
            structure=self._structure, cutoff=self._cutoff, node_row=self._node_row,
            node_image=self._node_image, node_members=self._node_members, i=self._i[mask], j=self._j[mask],
            shift=self._shift[mask], vector=self._vector[mask], distance=self._distance[mask],
        )

    # -- basic properties -------------------------------------------------
    @property
    def structure(self) -> Structure:
        """The structure the list was built on (its sites are the nodes)."""
        return self._structure

    @property
    def cutoff(self) -> float:
        """Search radius in Å used to build the list."""
        return self._cutoff

    @property
    def n_pairs(self) -> int:
        return len(self._i)

    def __len__(self) -> int:
        return self.n_pairs

    @property
    def n_nodes(self) -> int:
        return len(self._node_row)

    @property
    def i(self) -> np.ndarray:
        """(P,) origin node index."""
        return self._i

    @property
    def j(self) -> np.ndarray:
        """(P,) target node index."""
        return self._j

    @property
    def shift(self) -> np.ndarray:
        """(P, 3) integer image shift of the target, in units of the structure's (super)cell."""
        return self._shift

    @property
    def vector(self) -> np.ndarray:
        """(P, 3) Cartesian pair vector origin -> target image (Å)."""
        return self._vector

    @property
    def distance(self) -> np.ndarray:
        """(P,) pair distance (Å)."""
        return self._distance

    # -- node table ------------------------------------------------------
    @property
    def node_row(self) -> np.ndarray:
        """(n_nodes,) row of the first member atom of each node in ``structure``."""
        return self._node_row

    @property
    def node_image(self) -> np.ndarray:
        """(n_nodes, 3) whole-supercell offset (in parent cells) of each node's stored position.

        Zero unless the coordinates were wrapped or offset by whole supercells
        relative to ``cell_index``; used to keep :attr:`translation` independent of wrapping.
        """
        return self._node_image

    @property
    def node_site(self) -> np.ndarray:
        """(n_nodes,) site label of each node."""
        return self._structure.site_label[self._node_row]

    @property
    def node_cell_index(self) -> np.ndarray:
        """(n_nodes, 3) parent-cell index of each node."""
        return self._structure.cell_index[self._node_row]

    @property
    def node_frac(self) -> np.ndarray:
        """(n_nodes, 3) fractional position (of the structure's own cell) of each node."""
        return self._structure.frac[self._node_row]

    @property
    def node_atom_id(self) -> np.ndarray:
        """(n_nodes,) ``atom_id`` of the node's atom, or -1 where the site has several members."""
        out = self._structure.atom_id[self._node_row].copy()
        out[self._node_members != 1] = -1
        return out

    @property
    def node_species(self) -> np.ndarray:
        """(n_nodes,) species of the node's atom, or ``None`` where the site has several members."""
        out = self._structure.species[self._node_row].copy()
        out[self._node_members != 1] = None
        return out

    # -- per-pair views --------------------------------------------------
    @property
    def origin_site(self) -> np.ndarray:
        return self.node_site[self._i]

    @property
    def target_site(self) -> np.ndarray:
        return self.node_site[self._j]

    @property
    def origin_species(self) -> np.ndarray:
        return self.node_species[self._i]

    @property
    def target_species(self) -> np.ndarray:
        return self.node_species[self._j]

    @property
    def origin_atom_id(self) -> np.ndarray:
        return self.node_atom_id[self._i]

    @property
    def target_atom_id(self) -> np.ndarray:
        return self.node_atom_id[self._j]

    @property
    def translation(self) -> np.ndarray:
        """(P, 3) integer translation between the two **parent** unit cells.

        ``cell_index[target] + shift * supercell - cell_index[origin]``, where the
        cell index of each node is corrected for stored coordinates that sit whole
        supercells away from it (see :attr:`node_image`): the target site sits in the
        parent cell displaced by this many unit cells from the origin's cell, whatever
        way the coordinates were wrapped.  ``(origin site, target site, translation)`` is the same
        description of a pair that DISCUS uses in ``set vec``.
        """
        ci = self.node_cell_index + self._node_image
        return ci[self._j] + self._shift * self._structure.supercell - ci[self._i]

    @property
    def cell_vector(self) -> np.ndarray:
        """(P, 3) pair vector in fractional coordinates of the **parent** unit cell."""
        frac = self.node_frac
        return (frac[self._j] + self._shift - frac[self._i]) * self._structure.supercell

    # -- selection -------------------------------------------------------
    def _species_needed(self) -> np.ndarray:
        species = self.node_species
        if np.any(species == None):  # noqa: E711 (object array comparison)
            raise ValueError(
                "species selection needs one atom per site; this list was built on an "
                "average-model structure with mixed sites — select species in the "
                "correlation functions, which read species from a configuration"
            )
        return species

    def _direction_unit(self, direction, frame: str) -> np.ndarray:
        d = np.asarray(direction, dtype=float)
        if d.shape != (3,) or not np.all(np.isfinite(d)):
            raise ValueError(f"direction must be 3 finite numbers, got {direction!r}")
        if frame == "cartesian":
            v = d
        elif frame == "fractional":
            # [uvw] of the parent unit cell -> Cartesian
            parent = self._structure.lattice.matrix / self._structure.supercell[:, None]
            v = d @ parent
        else:
            raise ValueError(f"direction_frame must be 'cartesian' or 'fractional', got {frame!r}")
        norm = np.linalg.norm(v)
        if norm == 0:
            raise ValueError("direction must be nonzero")
        return v / norm

    def _mask(self, origin_site, target_site, origin_species, target_species,
              translation, cell_vector, vector_tol, direction, frame, angle_tol,
              distance, distance_tol):
        mask = np.ones(self.n_pairs, dtype=bool)
        osite = _as_set(origin_site, "origin_site")
        tsite = _as_set(target_site, "target_site")
        known = None
        for given, name in ((osite, "origin_site"), (tsite, "target_site")):
            if given is not None:
                known = set(self.node_site.astype(str)) if known is None else known
                unknown = sorted(given - known)
                if unknown:
                    raise ValueError(f"{name}: unknown site label(s) {unknown}; the structure has {sorted(known)[:10]}")
        if osite is not None:
            mask &= np.isin(self.origin_site.astype(str), list(osite))
        if tsite is not None:
            mask &= np.isin(self.target_site.astype(str), list(tsite))
        ospec = _as_set(origin_species, "origin_species")
        tspec = _as_set(target_species, "target_species")
        if ospec is not None or tspec is not None:
            species = self._species_needed().astype(str)
            if ospec is not None:
                mask &= np.isin(species[self._i], list(ospec))
            if tspec is not None:
                mask &= np.isin(species[self._j], list(tspec))
        if translation is not None:
            t = np.asarray(translation)
            if t.shape != (3,) or t.dtype.kind not in "iuf" or not np.all(t == np.round(t)):
                raise ValueError(f"translation must be 3 integers, got {translation!r}")
            mask &= np.all(self.translation == t.astype(np.int64), axis=1)
        if cell_vector is not None:
            cv = np.asarray(cell_vector, dtype=float)
            if cv.shape != (3,):
                raise ValueError(f"cell_vector must be 3 numbers, got {cell_vector!r}")
            mask &= np.all(np.abs(self.cell_vector - cv) <= vector_tol, axis=1)
        if direction is not None:
            if not (0 <= angle_tol <= 180):
                raise ValueError(f"angle_tol must be between 0 and 180 degrees, got {angle_tol!r}")
            u = self._direction_unit(direction, frame)
            with np.errstate(invalid="ignore", divide="ignore"):
                cosang = (self._vector @ u) / self._distance
            mask &= cosang >= np.cos(np.deg2rad(angle_tol)) - 1e-12
        if distance is not None:
            if np.ndim(distance) == 0:
                mask &= np.abs(self._distance - float(distance)) <= distance_tol
            else:
                lo, hi = (float(x) for x in distance)
                if lo > hi:
                    raise ValueError(f"distance range must be (min, max) with min <= max, got {distance!r}")
                mask &= (self._distance >= lo - distance_tol) & (self._distance <= hi + distance_tol)
        return mask

    def select(
        self,
        *,
        origin_site=None,
        target_site=None,
        origin_species=None,
        target_species=None,
        shell: int | None = None,
        shell_tol: float = 1e-3,
        distance=None,
        distance_tol: float = 1e-6,
        translation=None,
        cell_vector=None,
        vector_tol: float = 1e-6,
        direction=None,
        direction_frame: str = "cartesian",
        angle_tol: float = 1e-3,
        symmetric: bool = False,
    ) -> "NeighborList":
        """Return the pairs that satisfy **all** given criteria.

        Args:
            origin_site, target_site: site label or iterable of labels (an
                unknown label raises ``ValueError``, so typos are not silent).
            origin_species, target_species: species or iterable of species
                (needs one atom per site, i.e. a configuration).
            distance: a value (kept if ``|d - value| <= distance_tol``) or a
                ``(min, max)`` range (inclusive, widened by ``distance_tol``), Å.
            translation: exact integer parent-cell translation (see :attr:`translation`).
            cell_vector: pair vector in parent-cell fractions, matched to ``vector_tol``.
            direction: pair direction; kept if the angle to it is <= ``angle_tol``
                degrees.  ``direction_frame="fractional"`` reads it as [uvw] of the
                parent cell, ``"cartesian"`` as a Cartesian vector.
            shell: 1-based distance shell, counted **after** all other criteria
                are applied (see :meth:`shells`; shells are separated by gaps
                larger than ``shell_tol`` Å).
            symmetric: also keep the reversed pairs, i.e. pairs that match the
                criteria with origin and target swapped and translation /
                vector / direction negated.  ``select(origin_site="A",
                target_site="B", translation=t, symmetric=True)`` reproduces the
                pairs DISCUS uses for ``set vec, n, A, B, t``.
        """
        args = (origin_site, target_site, origin_species, target_species, translation,
                cell_vector, vector_tol, direction, direction_frame, angle_tol,
                distance, distance_tol)
        mask = self._mask(*args)
        if symmetric:
            neg = lambda v: None if v is None else -np.asarray(v, dtype=float)  # noqa: E731
            rev = (target_site, origin_site, target_species, origin_species,
                   None if translation is None else -np.asarray(translation),
                   neg(cell_vector), vector_tol, neg(direction), direction_frame, angle_tol,
                   distance, distance_tol)
            mask |= self._mask(*rev)
        out = self._subset(mask)
        if shell is not None:
            if isinstance(shell, bool) or not isinstance(shell, (int, np.integer)) or shell < 1:
                raise ValueError(f"shell must be an integer >= 1, got {shell!r}")
            labels = out._shell_labels(shell_tol)
            if out.n_pairs == 0 or shell > int(labels.max()):
                n = 0 if out.n_pairs == 0 else int(labels.max())
                raise ValueError(f"shell {shell} requested but the selection has {n} shell(s) within the cutoff")
            out = out._subset(labels == shell)
        return out

    def _shell_labels(self, tol: float) -> np.ndarray:
        if self.n_pairs == 0:
            return np.zeros(0, dtype=np.int64)
        order = np.argsort(self._distance, kind="stable")
        d = self._distance[order]
        new = np.concatenate(([True], np.diff(d) > tol))
        labels_sorted = np.cumsum(new)
        labels = np.empty(self.n_pairs, dtype=np.int64)
        labels[order] = labels_sorted
        return labels

    def shells(self, tol: float = 1e-3) -> list[Shell]:
        """Group the pairs into distance shells.

        Pairs sorted by distance start a new shell wherever the gap to the
        previous distance exceeds ``tol`` Å.  The last shell may be incomplete if
        the cutoff cuts through it; build with a larger cutoff to be sure.
        """
        labels = self._shell_labels(tol)
        n_origin_nodes = len(np.unique(self._i))
        out = []
        for k in range(1, int(labels.max(initial=0)) + 1):
            sel = labels == k
            d = self._distance[sel]
            out.append(Shell(
                index=k, distance=float(d.mean()), min_distance=float(d.min()),
                max_distance=float(d.max()), n_pairs=int(sel.sum()),
                n_origins=int(len(np.unique(self._i[sel]))),
                coordination=float(sel.sum()) / n_origin_nodes,
            ))
        return out

    def unique(self) -> "NeighborList":
        """Keep one direction of every unordered pair.

        A pair ``(i, j, s)`` is kept when ``i < j``, or ``i == j`` and ``s`` is
        lexicographically positive; its reverse ``(j, i, -s)`` is dropped.  If a
        pair's reverse is not in this list (an asymmetric selection), the pair is
        kept as it is, so the result never loses an unordered pair.  Costs about as
        much as building the list (a sort over the pair keys).
        """
        i, j, s = self._i, self._j, self._shift
        canonical = (i < j) | ((i == j) & (
            (s[:, 0] > 0) | ((s[:, 0] == 0) & ((s[:, 1] > 0) | ((s[:, 1] == 0) & (s[:, 2] > 0))))))
        # a non-canonical pair whose reverse is absent must be kept
        fwd = _row_keys(np.column_stack([i, j, s]))
        rev = _row_keys(np.column_stack([j, i, -s]))
        reverse_present = np.isin(rev, fwd)
        return self._subset(canonical | ~reverse_present)

    # -- summaries -------------------------------------------------------
    def coordination(self) -> np.ndarray:
        """(n_nodes,) number of pairs in this list that start at each node."""
        return np.bincount(self._i, minlength=self.n_nodes)

    def summary(self, tol: float = 1e-3) -> list[dict]:
        """Coordination table grouped by (origin site, target site, shell).

        Each row: ``origin_site``, ``target_site``, ``shell`` (1-based, counted
        separately for every origin/target site pair), ``distance`` (mean, Å),
        ``n_pairs``, ``n_origins`` (origin nodes of that site in the structure)
        and ``coordination`` (= n_pairs / n_origins).
        """
        if self.n_pairs == 0:
            return []
        names, node_code = np.unique(self.node_site.astype(str), return_inverse=True)
        node_code = node_code.reshape(-1)
        n_orig = np.bincount(node_code, minlength=len(names))
        key = node_code[self._i] * len(names) + node_code[self._j]
        order = np.lexsort((self._distance, key))
        k = key[order]
        d = self._distance[order]
        new_group = np.concatenate(([True], k[1:] != k[:-1]))
        new_shell = new_group | np.concatenate(([False], np.diff(d) > tol))
        starts = np.flatnonzero(new_shell)
        counts = np.diff(np.append(starts, len(d)))
        sums = np.add.reduceat(d, starts)
        group_id = np.cumsum(new_group) - 1
        # 1-based shell number within each (origin, target) group
        first_shell_of_group = np.searchsorted(starts, np.flatnonzero(new_group))
        shell_index = np.arange(len(starts)) - first_shell_of_group[group_id[starts]] + 1
        rows = []
        for s_k, c_k, sum_k, idx in zip(starts, counts, sums, shell_index):
            o, t = divmod(int(k[s_k]), len(names))
            rows.append({
                "origin_site": str(names[o]), "target_site": str(names[t]), "shell": int(idx),
                "distance": float(sum_k / c_k), "n_pairs": int(c_k),
                "n_origins": int(n_orig[o]), "coordination": float(c_k / n_orig[o]),
            })
        return rows

    def to_table(self) -> dict[str, np.ndarray]:
        """Per-pair columns as plain arrays (ready for ``pandas.DataFrame``)."""
        tr = self.translation
        cv = self.cell_vector
        table = {
            "origin_node": self._i.copy(), "target_node": self._j.copy(),
            "origin_atom_id": self.origin_atom_id, "target_atom_id": self.target_atom_id,
            "origin_site": self.origin_site.astype(str), "target_site": self.target_site.astype(str),
            "origin_species": self.origin_species.astype(str),
            "target_species": self.target_species.astype(str),
            "distance": self._distance.copy(),
        }
        for k, ax in enumerate("xyz"):
            table[f"r{ax}"] = self._vector[:, k].copy()
        for k, ax in enumerate("abc"):
            table[f"t{ax}"] = tr[:, k]
        for k, ax in enumerate("abc"):
            table[f"f{ax}"] = cv[:, k]
        for k, ax in enumerate("abc"):
            table[f"shift_{ax}"] = self._shift[:, k].copy()
        return table

    def to_csv(self, path) -> Path:
        """Write :meth:`to_table` as CSV (header + one row per pair); returns the path."""
        path = Path(path)
        table = self.to_table()
        names = list(table)
        with open(path, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh, lineterminator="\n")
            w.writerow(names)
            cols = [table[n] for n in names]
            for r in range(self.n_pairs):
                w.writerow([repr(float(c[r])) if c.dtype.kind == "f" else c[r] for c in cols])
        return path

    def __repr__(self) -> str:
        return (f"NeighborList(n_nodes={self.n_nodes}, n_pairs={self.n_pairs}, "
                f"cutoff={self._cutoff:g} Å, pbc={self._structure.pbc})")
