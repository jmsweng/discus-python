"""Occupational and displacement correlations over neighbor pairs (Package C, Features 13 and 14).

Both kinds of correlation are evaluated on the directed pairs of a
:class:`~discuspdf.neighbors.NeighborList` (usually built on an undisplaced
reference structure, or on the configuration itself) and read the species and
positions from a ``"configuration"`` :class:`~discuspdf.Structure` that has the
same lattice, supercell and site/cell bookkeeping.  Nodes of the list are
matched to atoms of the configuration by (``site_label``, ``cell_index``); a node
with no atom in the configuration is a **vacancy**; an atom that matches no
node is an error (it means the two structures label their sites differently).

Definitions, domains, null-model standard errors and the mapping to DISCUS
``chem corr`` are written out in ``docs/CORRELATIONS.md``.
"""
from __future__ import annotations

import dataclasses
import math
import warnings

import numpy as np

from .neighbors import NeighborList, _row_keys
from .structure import Structure

__all__ = [
    "OccupationalCorrelation",
    "DisplacementCorrelation",
    "node_species",
    "concentrations",
    "occupational_correlation",
    "occupational_correlation_table",
    "displacements",
    "displacement_correlation",
    "displacement_correlation_table",
]

_ZERO_VARIANCE = 1e-20   # Å²: projected displacements with rms below 1e-10 Å count as zero (round-off)
_LARGE_DISPLACEMENT = 0.4  # warn when |u| exceeds this fraction of the smallest supercell spacing


# ------------------------------------------------------------------ node <-> configuration mapping
def _check_compatible(nl: NeighborList, configuration: Structure) -> None:
    if not isinstance(configuration, Structure):
        raise TypeError(f"configuration must be a Structure, got {type(configuration).__name__}")
    ref = nl.structure
    if configuration.model != "configuration":
        raise ValueError(
            f"configuration must have model='configuration' (one realised atom per site), got "
            f"{configuration.model!r}; draw one with sample_configuration"
        )
    if not np.array_equal(configuration.supercell, ref.supercell):
        raise ValueError(
            f"configuration supercell {tuple(configuration.supercell)} differs from the "
            f"neighbor-list structure's {tuple(ref.supercell)}"
        )
    if not configuration.lattice.isclose(ref.lattice, rtol=1e-9, atol=1e-9):
        raise ValueError("configuration and neighbor-list structure have different lattices")
    if tuple(configuration.pbc) != tuple(ref.pbc):
        raise ValueError("configuration and neighbor-list structure have different pbc")


def _node_rows(nl: NeighborList, configuration: Structure | None) -> tuple[Structure, np.ndarray]:
    """(configuration, row of the configuration for each node or -1 for a vacancy).

    Raises if an atom of the configuration matches no node: the structures then
    label their sites differently (e.g. ``read_stru`` labels by species), and a
    silent mismatch would turn real atoms into vacancies.
    """
    if configuration is None:
        configuration = nl.structure
        if configuration.model != "configuration":
            raise ValueError(
                "the neighbor list was built on an average-model structure; pass the "
                "configuration whose species/positions should be analysed"
            )
        return configuration, np.array(nl.node_row, dtype=np.int64)
    _check_compatible(nl, configuration)
    sl = configuration.site_label
    ci = configuration.cell_index
    lookup = {(str(sl[r]), int(ci[r, 0]), int(ci[r, 1]), int(ci[r, 2])): r for r in range(configuration.n_atoms)}
    nsite = nl.node_site
    nci = nl.node_cell_index
    rows = np.array([lookup.get((str(nsite[k]), int(nci[k, 0]), int(nci[k, 1]), int(nci[k, 2])), -1)
                     for k in range(nl.n_nodes)], dtype=np.int64)
    n_unmatched = configuration.n_atoms - int(np.count_nonzero(rows >= 0))
    if n_unmatched:
        extra = sorted({str(x) for x in sl} - {str(x) for x in nsite})[:5]
        raise ValueError(
            f"{n_unmatched} of {configuration.n_atoms} configuration atoms match no node of the neighbor list "
            f"by (site_label, cell_index)"
            + (f"; site labels only in the configuration: {extra}" if extra else "")
            + ". Both structures must label sites the same way — for DISCUS files read with read_stru "
            "(labels by species) relabel both with discuspdf.neighbors.label_sites_by_slot"
        )
    return configuration, rows


def node_species(neighbors: NeighborList, configuration: Structure | None = None, *,
                 vacancy: str | None = "Va") -> np.ndarray:
    """Species at every node of ``neighbors`` in ``configuration`` (``vacancy`` where no atom).

    With ``configuration=None`` the list's own structure is used (it must be a
    configuration).  ``vacancy=None`` makes an empty node an error.
    """
    configuration, rows = _node_rows(neighbors, configuration)
    out = np.empty(neighbors.n_nodes, dtype=object)
    have = rows >= 0
    out[have] = [str(x) for x in configuration.species[rows[have]]]
    if not np.all(have):
        if vacancy is None:
            k = int(np.flatnonzero(~have)[0])
            raise ValueError(
                f"node {k} (site {str(neighbors.node_site[k])!r}, cell "
                f"{tuple(int(x) for x in neighbors.node_cell_index[k])}) has no atom in the configuration "
                "and vacancy=None"
            )
        if vacancy in {str(x) for x in configuration.species}:
            raise ValueError(f"vacancy label {vacancy!r} is also a species of the configuration")
        out[~have] = vacancy
    return out


def concentrations(neighbors: NeighborList, configuration: Structure | None = None, *,
                   vacancy: str | None = "Va") -> dict:
    """Site-resolved and global concentrations over the **nodes** of the list.

    Returns ``{"global": {species: fraction}, "by_site": {site_label: {species: fraction}},
    "n_nodes": int, "n_by_site": {site_label: int}}``.  These are per-site
    occupancies of the configuration, not the pair-weighted marginals that enter
    the correlation formulas (those are reported with each correlation).
    """
    spec = node_species(neighbors, configuration, vacancy=vacancy).astype(str)
    site = neighbors.node_site.astype(str)
    names = sorted(set(spec))
    out = {"global": {s: float(np.mean(spec == s)) for s in names}, "by_site": {}, "n_nodes": len(spec),
           "n_by_site": {}}
    for lab in sorted(set(site)):
        m = site == lab
        out["by_site"][lab] = {s: float(np.mean(spec[m] == s)) for s in names}
        out["n_by_site"][lab] = int(m.sum())
    return out


# ------------------------------------------------------------------ pair multiplicities (standard errors)
@dataclasses.dataclass(frozen=True)
class _PairCounts:
    n: int            # directed pairs
    n_reversed: int   # pairs whose exact reverse (j, i, -shift) is in the list
    n_self: int       # pairs i == j (a node with its own periodic image)
    s_same: int       # Σ m_ot² over directed node pairs o != t (m = multiplicity over shifts)
    s_rev: int        # Σ m_ot · m_to over directed node pairs o != t


def _pair_counts(i: np.ndarray, j: np.ndarray, shift: np.ndarray) -> _PairCounts:
    n = len(i)
    if n == 0:
        return _PairCounts(0, 0, 0, 0, 0)
    fwd = _row_keys(np.column_stack([i, j, shift]))
    rev = _row_keys(np.column_stack([j, i, -shift]))
    n_reversed = int(np.isin(rev, fwd).sum())
    self_pair = i == j
    io, jo = i[~self_pair], j[~self_pair]
    if len(io) == 0:
        return _PairCounts(n, n_reversed, int(self_pair.sum()), 0, 0)
    keys = _row_keys(np.column_stack([io, jo]))
    uniq, mult = np.unique(keys, return_counts=True)
    pairs = uniq.view(np.int64).reshape(-1, 2)
    rkeys = _row_keys(pairs[:, ::-1])
    pos = np.searchsorted(uniq, rkeys)
    pos = np.minimum(pos, len(uniq) - 1)
    has_rev = uniq[pos] == rkeys
    s_rev = int(np.sum(mult[has_rev] * mult[pos[has_rev]]))
    return _PairCounts(n, n_reversed, int(self_pair.sum()), int(np.sum(mult.astype(np.int64) ** 2)), s_rev)


def _null_variance_factor(pc: _PairCounts, rho):
    """Var(r) = (S_same + S_rev ρ²) / n² for independent sites; NaN if a node is paired with itself."""
    if pc.n == 0 or pc.n_self:
        return np.nan * np.asarray(rho, dtype=float)
    return (pc.s_same + pc.s_rev * np.asarray(rho, dtype=float) ** 2) / pc.n ** 2


_SELF_PAIR_NOTE = ("standard errors undefined: the selection pairs nodes with their own periodic images "
                   "(cutoff or translation spans the whole supercell); enlarge the supercell")


# ------------------------------------------------------------------ occupational (Feature 13)
@dataclasses.dataclass(frozen=True)
class OccupationalCorrelation:
    """Occupational pair statistics of one set of directed pairs.

    Matrices are indexed ``[a, b]`` with ``a`` the species at the **origin** and
    ``b`` at the **target**, in the order of :attr:`species`.  Undefined entries
    are NaN and listed with a reason in :attr:`undefined`.

    Attributes:
        species: species (and vacancy label) along both matrix axes.
        n_pairs: directed pairs counted (both ends among ``species``).
        n_unique_pairs: distinct unordered pairs among them.
        n_reversed: pairs whose reverse direction is also counted.
        n_self_pairs: pairs that join a node to its own periodic image.
        counts: ``N[a, b]``.
        joint: ``p[a, b] = N[a, b] / n_pairs``.
        origin_concentration, target_concentration: pair-weighted marginals
            ``c_o[a] = Σ_b p[a, b]``, ``c_t[b] = Σ_a p[a, b]``.
        conditional: ``P(b | a) = N[a, b] / Σ_b N[a, b]``.
        warren_cowley: ``α[a, b] = (P(b|a) − c_t[b]) / (δ_ab − c_t[b])`` (Cowley).
        pearson: ``r[a, b] = (p[a,b] − c_o[a] c_t[b]) / sqrt(c_o[a](1−c_o[a]) c_t[b](1−c_t[b]))``,
            the correlation coefficient of the occupation indicators at the two ends.
        stderr_pearson, stderr_warren_cowley: standard errors of ``r`` and ``α``
            under the null model of independent random occupation (see docs).
        undefined: ``{(a, b): reason}`` for every NaN value of ``warren_cowley`` or ``pearson``.
        stderr_note: why the standard errors are NaN although the values are defined, or ``None``.
    """
    species: tuple
    n_pairs: int
    n_unique_pairs: int
    n_reversed: int
    n_self_pairs: int
    counts: np.ndarray
    joint: np.ndarray
    origin_concentration: np.ndarray
    target_concentration: np.ndarray
    conditional: np.ndarray
    warren_cowley: np.ndarray
    pearson: np.ndarray
    stderr_pearson: np.ndarray
    stderr_warren_cowley: np.ndarray
    undefined: dict
    stderr_note: str | None = None

    def _ix(self, a: str, b: str) -> tuple[int, int]:
        try:
            return self.species.index(a), self.species.index(b)
        except ValueError:
            raise KeyError(f"species must be among {self.species}, got {a!r}, {b!r}") from None

    def value(self, kind: str, a: str, b: str) -> float:
        """One entry: ``kind`` in {"joint", "conditional", "warren_cowley", "pearson", "counts",
        "stderr_pearson", "stderr_warren_cowley"}."""
        ia, ib = self._ix(a, b)
        return float(getattr(self, kind)[ia, ib])

    def discus(self, a: str, b: str) -> float:
        """DISCUS ``chem corr occ, a, b`` value from these pairs (a ≠ b).

        DISCUS restricts the pairs to those with both ends in {a, b}, pools the
        two directions, and uses ``c = (P_aa − θ²) / (θ (1 − θ))`` with
        ``θ = (2 N_aa + N_ab + N_ba) / (2 N)``.  DISCUS reports 0 when θ is 0 or 1;
        here that case is NaN.  For a list that contains both directions of
        every pair this equals ``pearson[a, a]`` (and ``−pearson[a, b]`` for a binary).
        """
        if a == b:
            raise ValueError("DISCUS occupational correlations need two different species")
        ia, ib = self._ix(a, b)
        N = self.counts
        naa, nab, nba, nbb = (int(x) for x in (N[ia, ia], N[ia, ib], N[ib, ia], N[ib, ib]))
        n = naa + nab + nba + nbb
        num = 2 * naa + nab + nba        # θ = num / (2n), exact in integers
        if n == 0 or num == 0 or num == 2 * n:
            return math.nan
        theta = num / (2 * n)
        return float((naa / n - theta ** 2) / (theta * (1 - theta)))


def _occupational_core(o: np.ndarray, t: np.ndarray, names: list, pc: _PairCounts,
                       pooled: np.ndarray) -> OccupationalCorrelation:
    S = len(names)
    n = len(o)
    counts = np.bincount(o * S + t, minlength=S * S).reshape(S, S).astype(np.int64)
    undefined: dict = {}
    if n == 0:
        nan = np.full((S, S), np.nan)
        for a in names:
            for b in names:
                undefined[(a, b)] = "no pairs"
        z = np.full(S, np.nan)
        return OccupationalCorrelation(tuple(names), 0, 0, 0, 0, counts, nan.copy(), z, z.copy(), nan.copy(),
                                       nan.copy(), nan.copy(), nan.copy(), nan.copy(), undefined)
    row_n = counts.sum(axis=1)          # origins of each species (integer)
    col_n = counts.sum(axis=0)          # targets of each species (integer)
    joint = counts / n
    co = row_n / n
    ct = col_n / n
    delta = np.eye(S, dtype=bool)
    # Zero denominators are decided on the integer counts, never on rounded fractions.
    origin_absent = (row_n == 0)[:, None] & np.ones((1, S), bool)
    target_absent = np.ones((S, 1), bool) & (col_n == 0)[None, :]
    origin_pure = (row_n == n)[:, None] & np.ones((1, S), bool)
    target_pure = np.ones((S, 1), bool) & (col_n == n)[None, :]
    wc_zero_den = (~delta & target_absent) | (delta & target_pure)
    wc_undefined = origin_absent | wc_zero_den
    pr_undefined = origin_absent | target_absent | origin_pure | target_pure
    with np.errstate(invalid="ignore", divide="ignore"):
        cond = np.where(row_n[:, None] > 0, counts / np.where(row_n[:, None] > 0, row_n[:, None], 1), np.nan)
        wc_den = delta - ct[None, :]
        wc = np.where(wc_undefined, np.nan, (cond - ct[None, :]) / np.where(wc_undefined, 1.0, wc_den))
        pr_den = np.sqrt(np.outer(co * (1 - co), ct * (1 - ct)))
        pr = np.where(pr_undefined, np.nan,
                      (joint - np.outer(co, ct)) / np.where(pr_undefined, 1.0, pr_den))
        # Null model (independent random occupation, pooled concentrations of the nodes involved):
        # rho_ab = correlation of the a- and b-indicators at one site.
        cp = pooled
        sd_site = np.sqrt(cp * (1 - cp))
        rho_site = (np.diag(cp) - np.outer(cp, cp)) / np.outer(sd_site, sd_site)
        se_r = np.sqrt(_null_variance_factor(pc, rho_site))
        se_r = np.where(pr_undefined, np.nan, se_r)
        factor = pr_den / (co[:, None] * wc_den)
        se_wc = np.where(wc_undefined | pr_undefined, np.nan, np.abs(factor) * se_r)
    for ia, a in enumerate(names):
        for ib, b in enumerate(names):
            if not (wc_undefined[ia, ib] or pr_undefined[ia, ib]):
                continue
            if origin_absent[ia, ib]:
                reason = f"species {a!r} absent at the origins"
            elif wc_zero_den[ia, ib]:
                reason = (f"species {b!r} absent at the targets" if ia != ib
                          else f"pure end member: every target is {b!r} (δ_ab − c_b = 0)")
            elif target_absent[ia, ib]:
                reason = f"species {b!r} absent at the targets (correlation coefficient undefined)"
            else:
                reason = "pure end member (zero variance; correlation coefficient undefined)"
            undefined[(a, b)] = reason
    note = _SELF_PAIR_NOTE if pc.n_self else None
    return OccupationalCorrelation(
        species=tuple(names), n_pairs=n, n_unique_pairs=n - pc.n_reversed // 2, n_reversed=pc.n_reversed,
        n_self_pairs=pc.n_self, counts=counts, joint=joint, origin_concentration=co,
        target_concentration=ct, conditional=cond, warren_cowley=wc, pearson=pr, stderr_pearson=se_r,
        stderr_warren_cowley=se_wc, undefined=undefined, stderr_note=note,
    )


def _species_axis(spec: np.ndarray, species) -> list[str]:
    if species is None:
        return sorted({str(x) for x in spec})
    names = [str(s) for s in ([species] if isinstance(species, str) else species)]
    if len(set(names)) != len(names) or not names:
        raise ValueError(f"species must be distinct and non-empty, got {species!r}")
    return names


def _occupational(neighbors: NeighborList, mask, spec: np.ndarray, names: list) -> OccupationalCorrelation:
    """Core shared by the single call and the tables (``spec``: species per node, already mapped)."""
    names_arr = np.array(names)
    order = np.argsort(names_arr)
    i = neighbors.i if mask is None else neighbors.i[mask]
    j = neighbors.j if mask is None else neighbors.j[mask]
    shift = neighbors.shift if mask is None else neighbors.shift[mask]
    node_code = np.full(len(spec), -1, dtype=np.int64)
    pos = np.searchsorted(names_arr[order], spec)
    pos = np.minimum(pos, len(names) - 1)
    hit = names_arr[order][pos] == spec
    node_code[hit] = order[pos[hit]]
    o, t = node_code[i], node_code[j]
    keep = (o >= 0) & (t >= 0)
    i, j, shift, o, t = i[keep], j[keep], shift[keep], o[keep], t[keep]
    nodes = np.unique(np.concatenate([i, j]))
    pooled = (np.bincount(node_code[nodes], minlength=len(names)) / len(nodes) if nodes.size
              else np.full(len(names), np.nan))
    return _occupational_core(o, t, names, _pair_counts(i, j, shift), pooled)


def occupational_correlation(neighbors: NeighborList, configuration: Structure | None = None, *,
                             species=None, vacancy: str | None = "Va") -> OccupationalCorrelation:
    """Pair probabilities, Warren–Cowley and correlation coefficients for the pairs of ``neighbors``.

    Args:
        neighbors: the (selected) directed pairs to analyse — e.g.
            ``nl.select(shell=1)`` for a shell average or
            ``nl.select(origin_site=..., target_site=..., translation=..., symmetric=True)``
            for one signed neighbor vector.
        configuration: the configuration to read species from (default: the
            list's own structure).  Must match lattice, supercell and pbc, and
            every atom must match a node (see module docstring).
        species: restrict to pairs with **both** ends among these species (as
            DISCUS does); default: all species found at the nodes of the list
            (whether or not they occur in the selected pairs), plus the vacancy
            label if any node is empty.
        vacancy: label of empty nodes; ``None`` makes them an error.
    """
    spec = node_species(neighbors, configuration, vacancy=vacancy).astype(str)
    return _occupational(neighbors, None, spec, _species_axis(spec, species))


def _group_labels(nl: NeighborList, by: str, shell_tol: float):
    """[(keys, boolean mask)] per group, in order of distance."""
    if nl.n_pairs == 0:
        if by not in ("shell", "vector"):
            raise ValueError(f"by must be 'shell' or 'vector', got {by!r}")
        return []
    if by == "shell":
        labels = nl._shell_labels(shell_tol)
        out = []
        for k in range(1, int(labels.max()) + 1):
            m = labels == k
            out.append(({"shell": k, "distance": float(nl.distance[m].mean())}, m))
        return out
    if by == "vector":
        names, code = np.unique(nl.node_site.astype(str), return_inverse=True)
        code = code.reshape(-1)
        rows = np.column_stack([code[nl.i], code[nl.j], nl.translation])
        uniq, first, inverse = np.unique(rows, axis=0, return_index=True, return_inverse=True)
        inverse = inverse.reshape(-1)
        order = np.lexsort((uniq[:, 4], uniq[:, 3], uniq[:, 2], uniq[:, 1], uniq[:, 0],
                            np.round(nl.distance[first], 9)))
        out = []
        for g in order:
            m = inverse == g
            out.append(({"origin_site": str(names[uniq[g, 0]]), "target_site": str(names[uniq[g, 1]]),
                         "translation": tuple(int(x) for x in uniq[g, 2:]),
                         "distance": float(nl.distance[first[g]])}, m))
        return out
    raise ValueError(f"by must be 'shell' or 'vector', got {by!r}")


def occupational_correlation_table(neighbors: NeighborList, configuration: Structure | None = None, *,
                                   by: str = "shell", species=None, vacancy: str | None = "Va",
                                   shell_tol: float = 1e-3) -> list[dict]:
    """One row per (group, origin species, target species).

    ``by="shell"`` groups the pairs into distance shells (shell-averaged values);
    ``by="vector"`` groups them by (origin site, target site, parent-cell
    translation) — one signed neighbor vector per group.  Each row has the group
    keys, ``n_pairs``, ``n_unique_pairs``, ``a``, ``b``, ``count``, ``joint``,
    ``conditional``, ``warren_cowley``, ``pearson`` and both standard errors.
    """
    spec = node_species(neighbors, configuration, vacancy=vacancy).astype(str)
    names = _species_axis(spec, species)
    rows = []
    for key, mask in _group_labels(neighbors, by, shell_tol):
        res = _occupational(neighbors, mask, spec, names)
        for ia, a in enumerate(res.species):
            for ib, b in enumerate(res.species):
                rows.append({**key, "n_pairs": res.n_pairs, "n_unique_pairs": res.n_unique_pairs,
                             "a": a, "b": b, "count": int(res.counts[ia, ib]),
                             "joint": float(res.joint[ia, ib]), "conditional": float(res.conditional[ia, ib]),
                             "warren_cowley": float(res.warren_cowley[ia, ib]),
                             "pearson": float(res.pearson[ia, ib]),
                             "stderr_warren_cowley": float(res.stderr_warren_cowley[ia, ib]),
                             "stderr_pearson": float(res.stderr_pearson[ia, ib])})
    return rows


# ------------------------------------------------------------------ displacement (Feature 14)
def displacements(neighbors: NeighborList, configuration: Structure, *, reference="neighbors") -> np.ndarray:
    """(n_nodes, 3) Cartesian displacement of every node's atom from its reference position (Å).

    Args:
        neighbors: defines the nodes; its structure gives the reference positions.
        configuration: the displaced configuration (same lattice/supercell/pbc).
        reference: ``"neighbors"`` — displacement from the node positions of the
            list's structure (the ideal/undisplaced reference); ``"site_average"``
            — additionally subtract, for every site label, the mean displacement
            of that site over all cells (DISCUS ``chem`` measures displacements
            from the average position of each site in this way); or a
            configuration ``Structure`` with the same sites (positions taken from it).

    The difference is taken as the minimum image in the supercell for periodic
    structures (so wrapping of the configuration does not matter).  A
    displacement longer than half the shortest supercell spacing would be folded
    back; a ``RuntimeWarning`` is issued when any displacement exceeds 0.4 of it.
    Rows of vacant nodes are NaN.
    """
    config, rows = _node_rows(neighbors, configuration)
    if isinstance(reference, Structure):
        if reference.model != "configuration":
            raise ValueError("a reference Structure must be a configuration with the same sites")
        ref_rows = _node_rows(neighbors, reference)[1]
        if np.any(ref_rows < 0):
            raise ValueError("the reference structure lacks some nodes of the neighbor list")
        ref_frac = reference.frac[ref_rows]
        mode = "structure"
    elif isinstance(reference, str) and reference in ("neighbors", "site_average"):
        ref_frac = neighbors.node_frac
        mode = reference
    else:
        raise ValueError(f"reference must be 'neighbors', 'site_average' or a Structure, got {reference!r}")
    u = np.full((neighbors.n_nodes, 3), np.nan)
    have = rows >= 0
    df = config.frac[rows[have]] - ref_frac[have]
    if all(config.pbc):
        df, _ = config.lattice.minimum_image(df)
        spacing = float(np.min(1.0 / np.linalg.norm(config.lattice.reciprocal_matrix, axis=1)))
        longest = float(np.max(np.linalg.norm(df @ config.lattice.matrix, axis=1), initial=0.0))
        if longest > _LARGE_DISPLACEMENT * spacing:
            warnings.warn(
                f"a displacement of {longest:.3g} Å exceeds {_LARGE_DISPLACEMENT} of the smallest supercell "
                f"spacing ({spacing:.3g} Å); displacements beyond half of it are folded back by the "
                "minimum image — check the reference structure", RuntimeWarning, stacklevel=2)
    u[have] = df @ config.lattice.matrix
    if mode == "site_average":
        site = neighbors.node_site.astype(str)
        for lab in set(site):
            m = (site == lab) & have
            if np.any(m):
                u[m] -= u[m].mean(axis=0)
    return u


@dataclasses.dataclass(frozen=True)
class DisplacementCorrelation:
    """Displacement statistics of one set of directed pairs.

    Attributes:
        n_pairs: pairs used (both ends occupied and of the selected species).
        n_skipped: pairs dropped because an end is vacant or of another species.
        n_unique_pairs, n_reversed, n_self_pairs: as for :class:`OccupationalCorrelation`.
        centered: whether means over the pairs were subtracted (Pearson) or not (DISCUS).
        origin_mean, target_mean: (3,) mean displacement of the origins / targets over the pairs (Å).
        covariance: (3, 3) ``<(u_o − m_o)(u_t − m_t)ᵀ>`` over the pairs (Å²; means are 0 if not centered).
        origin_covariance, target_covariance: (3, 3) ``<(u − m)(u − m)ᵀ>`` of the origins / targets.
        origin_direction, target_direction: Cartesian unit vectors, or ``None``.
        correlation: ``<x y> / sqrt(<x²> <y²>)`` of the projections
            ``x = (u_o − m_o)·d_o``, ``y = (u_t − m_t)·d_t``; NaN when undefined.
        stderr_null: standard error of ``correlation`` under the null model of
            independent displacements (see docs); NaN when undefined.
        undefined: reason when ``correlation`` is NaN, else ``None``.
        stderr_note: why ``stderr_null`` is NaN although ``correlation`` is defined, or ``None``.
    """
    n_pairs: int
    n_skipped: int
    n_unique_pairs: int
    n_reversed: int
    n_self_pairs: int
    centered: bool
    origin_mean: np.ndarray
    target_mean: np.ndarray
    covariance: np.ndarray
    origin_covariance: np.ndarray
    target_covariance: np.ndarray
    origin_direction: np.ndarray | None
    target_direction: np.ndarray | None
    correlation: float
    stderr_null: float
    undefined: str | None
    stderr_note: str | None = None

    @property
    def correlation_matrix(self) -> np.ndarray:
        """(3, 3) ``C[a, b] / sqrt(Vo[a, a] Vt[b, b])`` for Cartesian components (NaN where a variance is 0)."""
        vo = np.diag(self.origin_covariance)
        vt = np.diag(self.target_covariance)
        with np.errstate(invalid="ignore", divide="ignore"):
            den = np.sqrt(np.outer(vo, vt))
            ok = den > _ZERO_VARIANCE
            return np.where(ok, self.covariance / np.where(ok, den, 1), np.nan)


def _species_mask(spec: np.ndarray, value, nodes: np.ndarray) -> np.ndarray:
    if value is None:
        return np.ones(len(nodes), dtype=bool)
    vals = [value] if isinstance(value, str) else list(value)
    return np.isin(spec[nodes], [str(v) for v in vals])


def _displacement(neighbors: NeighborList, mask, u: np.ndarray, spec: np.ndarray, do, dt, centered: bool,
                  origin_species, target_species) -> DisplacementCorrelation:
    i = neighbors.i if mask is None else neighbors.i[mask]
    j = neighbors.j if mask is None else neighbors.j[mask]
    shift = neighbors.shift if mask is None else neighbors.shift[mask]
    total = len(i)
    keep = ~np.isnan(u[i, 0]) & ~np.isnan(u[j, 0])
    keep &= _species_mask(spec, origin_species, i) & _species_mask(spec, target_species, j)
    i, j, shift = i[keep], j[keep], shift[keep]
    n = len(i)
    pc = _pair_counts(i, j, shift)
    nan3, nan33 = np.full(3, np.nan), np.full((3, 3), np.nan)
    if n == 0:
        return DisplacementCorrelation(0, total, 0, 0, 0, centered, nan3, nan3, nan33, nan33, nan33,
                                       do, dt, math.nan, math.nan, "no pairs with both ends occupied")
    uo, ut = u[i], u[j]
    mo = uo.mean(axis=0) if centered else np.zeros(3)
    mt = ut.mean(axis=0) if centered else np.zeros(3)
    ao, at = uo - mo, ut - mt
    cov, cov_o, cov_t = ao.T @ at / n, ao.T @ ao / n, at.T @ at / n
    corr, se, reason, note = math.nan, math.nan, None, None
    if do is None:
        reason = "no direction given"
    else:
        x, y = ao @ do, at @ dt
        sxx, syy = float(x @ x / n), float(y @ y / n)
        if sxx <= _ZERO_VARIANCE or syy <= _ZERO_VARIANCE:
            reason = "zero variance of the projected displacements (rms < 1e-10 Å)"
        else:
            corr = float(x @ y / n / math.sqrt(sxx * syy))
            # Null model: independent displacements; a pair counted in both orientations
            # contributes the same-atom correlation of the two projections (rho_self).
            nodes = np.unique(np.concatenate([i, j]))
            un = u[nodes] - (u[nodes].mean(axis=0) if centered else 0.0)
            a1, a2 = un @ do, un @ dt
            den = math.sqrt(float(a1 @ a1) * float(a2 @ a2))
            rho_self = float(a1 @ a2) / den if den > 0 else 0.0
            se = float(np.sqrt(_null_variance_factor(pc, rho_self)))
            if pc.n_self:
                note = _SELF_PAIR_NOTE
    return DisplacementCorrelation(
        n_pairs=n, n_skipped=total - n, n_unique_pairs=n - pc.n_reversed // 2, n_reversed=pc.n_reversed,
        n_self_pairs=pc.n_self, centered=centered,
        origin_mean=mo if centered else uo.mean(axis=0), target_mean=mt if centered else ut.mean(axis=0),
        covariance=cov, origin_covariance=cov_o, target_covariance=cov_t, origin_direction=do,
        target_direction=dt, correlation=corr, stderr_null=se, undefined=reason, stderr_note=note,
    )


def _directions(neighbors, direction, target_direction, frame):
    if direction is None:
        if target_direction is not None:
            raise ValueError("target_direction needs direction")
        return None, None
    do = neighbors._direction_unit(direction, frame)
    dt = do if target_direction is None else neighbors._direction_unit(target_direction, frame)
    return do, dt


def displacement_correlation(neighbors: NeighborList, configuration: Structure, *, direction=None,
                             target_direction=None, direction_frame: str = "cartesian",
                             reference="neighbors", centered: bool = True, origin_species=None,
                             target_species=None) -> DisplacementCorrelation:
    """Displacement covariances and the normalised directional correlation over the pairs.

    Args:
        neighbors: the (selected) directed pairs, ideally built on the undisplaced reference.
        configuration: the displaced configuration.
        direction: projection direction at the origin; ``target_direction``
            defaults to the same.  ``direction_frame="fractional"`` reads them as
            [uvw] of the parent cell.  Without a direction only the covariance
            matrices are computed and ``correlation`` is NaN.  (DISCUS's
            direction-free ``dir all`` mode is not implemented.)
        reference: see :func:`displacements`.
        centered: subtract the pair means (Pearson correlation, default) or not
            (DISCUS ``chem corr disp`` convention ``<xy>/sqrt(<x²><y²>)``).
        origin_species, target_species: keep only pairs whose ends carry these species.
    """
    do, dt = _directions(neighbors, direction, target_direction, direction_frame)
    u = displacements(neighbors, configuration, reference=reference)
    spec = node_species(neighbors, configuration).astype(str)
    return _displacement(neighbors, None, u, spec, do, dt, centered, origin_species, target_species)


def displacement_correlation_table(neighbors: NeighborList, configuration: Structure, *, direction,
                                   target_direction=None, direction_frame: str = "cartesian",
                                   by: str = "shell", reference="neighbors", centered: bool = True,
                                   origin_species=None, target_species=None,
                                   shell_tol: float = 1e-3) -> list[dict]:
    """One row per group (``by="shell"`` or ``by="vector"``, as in
    :func:`occupational_correlation_table`) with ``n_pairs``, ``correlation``,
    ``stderr_null`` and ``covariance`` (the projected ``<x y>`` in Å²)."""
    do, dt = _directions(neighbors, direction, target_direction, direction_frame)
    u = displacements(neighbors, configuration, reference=reference)
    spec = node_species(neighbors, configuration).astype(str)
    rows = []
    for key, mask in _group_labels(neighbors, by, shell_tol):
        res = _displacement(neighbors, mask, u, spec, do, dt, centered, origin_species, target_species)
        proj = (math.nan if do is None or res.n_pairs == 0 else float(do @ res.covariance @ dt))
        rows.append({**key, "n_pairs": res.n_pairs, "n_skipped": res.n_skipped,
                     "correlation": res.correlation, "stderr_null": res.stderr_null, "covariance": proj,
                     "undefined": res.undefined})
    return rows
