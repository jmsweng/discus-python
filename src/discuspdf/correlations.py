"""Occupational and displacement correlations over neighbor pairs (Package C, Features 13 and 14).

Both kinds of correlation are evaluated on the directed pairs of a
:class:`~discuspdf.neighbors.NeighborList` (usually built on an undisplaced
reference structure, or on the configuration itself) and read the species and
positions from a ``"configuration"`` :class:`~discuspdf.Structure` that has the
same lattice, supercell and site/cell bookkeeping.  Nodes of the list are
matched to atoms of the configuration by (``site_label``, ``cell_index``); a node
with no atom in the configuration is a **vacancy**.

Definitions, domains, null-model standard errors and the mapping to DISCUS
``chem corr`` are written out in ``docs/CORRELATIONS.md``.
"""
from __future__ import annotations

import dataclasses
import math

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
    """(configuration, row of the configuration for each node or -1 for a vacancy)."""
    if configuration is None:
        configuration = nl.structure
        if configuration.model != "configuration":
            raise ValueError(
                "the neighbor list was built on an average-model structure; pass the "
                "configuration whose species/positions should be analysed"
            )
        rows = np.full(nl.n_nodes, -1, dtype=np.int64)
        rows[:] = nl.node_row
        return configuration, rows
    _check_compatible(nl, configuration)
    sl = configuration.site_label
    ci = configuration.cell_index
    lookup = {(sl[r], int(ci[r, 0]), int(ci[r, 1]), int(ci[r, 2])): r for r in range(configuration.n_atoms)}
    nsite = nl.node_site
    nci = nl.node_cell_index
    rows = np.array([lookup.get((nsite[k], int(nci[k, 0]), int(nci[k, 1]), int(nci[k, 2])), -1)
                     for k in range(nl.n_nodes)], dtype=np.int64)
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
    out[have] = configuration.species[rows[have]]
    if not np.all(have):
        if vacancy is None:
            k = int(np.flatnonzero(~have)[0])
            raise ValueError(
                f"node {k} (site {neighbors.node_site[k]!r}, cell "
                f"{tuple(int(x) for x in neighbors.node_cell_index[k])}) has no atom in the configuration "
                "and vacancy=None"
            )
        if vacancy in set(configuration.species):
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


def _n_reversed(nl: NeighborList) -> int:
    """Number of directed pairs whose reverse (j, i, -shift) is also in the list."""
    if nl.n_pairs == 0:
        return 0
    fwd = _row_keys(np.column_stack([nl.i, nl.j, nl.shift]))
    rev = _row_keys(np.column_stack([nl.j, nl.i, -nl.shift]))
    return int(np.isin(rev, fwd).sum())


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
        undefined: ``{(a, b): reason}`` for every NaN entry.
    """
    species: tuple
    n_pairs: int
    n_unique_pairs: int
    n_reversed: int
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
        naa, nab, nba, nbb = N[ia, ia], N[ia, ib], N[ib, ia], N[ib, ib]
        n = naa + nab + nba + nbb
        if n == 0:
            return math.nan
        theta = (2 * naa + nab + nba) / (2 * n)
        if theta <= 0 or theta >= 1:
            return math.nan
        return float((naa / n - theta ** 2) / (theta * (1 - theta)))


def _occupational_from_codes(o: np.ndarray, t: np.ndarray, names: list, n_unique: int, n_rev: int,
                             pooled: np.ndarray) -> OccupationalCorrelation:
    S = len(names)
    n = len(o)
    counts = np.zeros((S, S), dtype=np.int64)
    np.add.at(counts, (o, t), 1)
    undefined: dict = {}
    nan = np.full((S, S), np.nan)
    if n == 0:
        for a in names:
            for b in names:
                undefined[(a, b)] = "no pairs"
        z = np.zeros(S)
        return OccupationalCorrelation(tuple(names), 0, 0, 0, counts, nan.copy(), z * np.nan, z * np.nan,
                                       nan.copy(), nan.copy(), nan.copy(), nan.copy(), nan.copy(), undefined)
    joint = counts / n
    co = joint.sum(axis=1)
    ct = joint.sum(axis=0)
    row = counts.sum(axis=1, keepdims=True)
    with np.errstate(invalid="ignore", divide="ignore"):
        cond = np.where(row > 0, counts / np.where(row > 0, row, 1), np.nan)
        delta = np.eye(S)
        wc_den = delta - ct[None, :]
        wc = (cond - ct[None, :]) / wc_den
        var_o = co * (1 - co)
        var_t = ct * (1 - ct)
        pr_den = np.sqrt(np.outer(var_o, var_t))
        pr = (joint - np.outer(co, ct)) / pr_den
        # null-model standard error: Var(r) = (n + n_rev rho_ab^2) / n^2 with rho_ab the
        # correlation of the a- and b-indicators at one site (pooled concentrations)
        cp = pooled
        cov_site = np.diag(cp) - np.outer(cp, cp)
        sd_site = np.sqrt(cp * (1 - cp))
        rho_site = cov_site / np.outer(sd_site, sd_site)
        se_r = np.sqrt((n + n_rev * rho_site ** 2) / n ** 2)
        factor = pr_den / (co[:, None] * wc_den)
        se_wc = np.abs(factor) * se_r
    for ia, a in enumerate(names):
        for ib, b in enumerate(names):
            reason = None
            if co[ia] == 0:
                reason = f"species {a!r} absent at the origins"
            elif ct[ib] == 0:
                reason = f"species {b!r} absent at the targets"
            elif co[ia] == 1 or ct[ib] == 1:
                reason = "pure end member (zero variance)"
            if reason is not None:
                pr[ia, ib] = np.nan
                se_r[ia, ib] = np.nan
                # α stays defined where its own denominators are nonzero
                if co[ia] == 0 or ct[ib] == 0 or (ia == ib and ct[ib] == 1):
                    wc[ia, ib] = np.nan
                    se_wc[ia, ib] = np.nan
                    if ia == ib and ct[ib] == 1 and co[ia] != 0:
                        reason = "pure end member (δ_ab − c_b = 0)"
                else:
                    se_wc[ia, ib] = np.nan
                undefined[(a, b)] = reason
            if np.isnan(cond[ia, ib]) and (a, b) not in undefined:
                undefined[(a, b)] = f"species {a!r} absent at the origins"
    return OccupationalCorrelation(
        species=tuple(names), n_pairs=n, n_unique_pairs=n_unique, n_reversed=n_rev, counts=counts,
        joint=joint, origin_concentration=co, target_concentration=ct, conditional=cond,
        warren_cowley=wc, pearson=pr, stderr_pearson=se_r, stderr_warren_cowley=se_wc, undefined=undefined,
    )


def occupational_correlation(neighbors: NeighborList, configuration: Structure | None = None, *,
                             species=None, vacancy: str | None = "Va") -> OccupationalCorrelation:
    """Pair probabilities, Warren–Cowley and correlation coefficients for the pairs of ``neighbors``.

    Args:
        neighbors: the (selected) directed pairs to analyse — e.g.
            ``nl.select(shell=1)`` for a shell average or
            ``nl.select(origin_site=..., target_site=..., translation=..., symmetric=True)``
            for one signed neighbor vector.
        configuration: the configuration to read species from (default: the
            list's own structure).  Must match lattice, supercell and pbc.
        species: restrict to pairs with **both** ends among these species (as
            DISCUS does); default: all species found at the nodes of the list
            (whether or not they occur in the selected pairs), plus the vacancy
            label if any node is empty.
        vacancy: label of empty nodes; ``None`` makes them an error.
    """
    spec = node_species(neighbors, configuration, vacancy=vacancy).astype(str)
    if species is None:
        names = sorted(set(spec))
    else:
        names = [str(s) for s in ([species] if isinstance(species, str) else species)]
        if len(set(names)) != len(names) or not names:
            raise ValueError(f"species must be distinct and non-empty, got {species!r}")
    code = {s: k for k, s in enumerate(names)}
    so = spec[neighbors.i]
    st = spec[neighbors.j]
    keep = np.isin(so, names) & np.isin(st, names)
    sub = neighbors._subset(keep)
    o = np.array([code[s] for s in so[keep]], dtype=np.int64)
    t = np.array([code[s] for s in st[keep]], dtype=np.int64)
    nodes = np.unique(np.concatenate([sub.i, sub.j]))
    pooled = np.array([np.mean(spec[nodes] == s) if nodes.size else np.nan for s in names])
    n_rev = _n_reversed(sub)
    n_unique = sub.n_pairs - n_rev // 2
    return _occupational_from_codes(o, t, names, n_unique, n_rev, pooled)


def _group_labels(nl: NeighborList, by: str, shell_tol: float):
    if by == "shell":
        labels = nl._shell_labels(shell_tol)
        keys = []
        for k in range(1, int(labels.max(initial=0)) + 1):
            d = nl.distance[labels == k]
            keys.append(({"shell": k, "distance": float(d.mean())}, labels == k))
        return keys
    if by == "vector":
        tr = nl.translation
        osite = nl.origin_site.astype(str)
        tsite = nl.target_site.astype(str)
        groups: dict = {}
        for p in range(nl.n_pairs):
            key = (str(osite[p]), str(tsite[p]), int(tr[p, 0]), int(tr[p, 1]), int(tr[p, 2]))
            groups.setdefault(key, []).append(p)
        out = []
        for key in sorted(groups, key=lambda k: (float(nl.distance[groups[k][0]]), k)):
            m = np.zeros(nl.n_pairs, dtype=bool)
            m[groups[key]] = True
            out.append(({"origin_site": key[0], "target_site": key[1], "translation": key[2:],
                         "distance": float(nl.distance[groups[key][0]])}, m))
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
    rows = []
    for key, mask in _group_labels(neighbors, by, shell_tol):
        res = occupational_correlation(neighbors._subset(mask), configuration, species=species, vacancy=vacancy)
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
_ZERO_VARIANCE = 1e-20   # Å²: projected displacements with rms below 1e-10 Å count as zero (round-off)

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
            ``Structure`` with the same nodes (positions taken from it).

    The difference is taken as the minimum image in the structure's cell for
    periodic structures (so wrapping of the configuration does not matter);
    displacements longer than half the shortest supercell spacing are ambiguous.
    Rows of vacant nodes are NaN.
    """
    config, rows = _node_rows(neighbors, configuration)
    if isinstance(reference, Structure):
        ref_nl_rows = _node_rows(neighbors, reference)[1] if reference.model == "configuration" else None
        if ref_nl_rows is None:
            raise ValueError("a reference Structure must be a configuration with the same sites")
        if np.any(ref_nl_rows < 0):
            raise ValueError("the reference structure lacks some nodes of the neighbor list")
        ref_frac = reference.frac[ref_nl_rows]
        mode = "neighbors"
    elif reference in ("neighbors", "site_average"):
        ref_frac = neighbors.node_frac
        mode = reference
    else:
        raise ValueError(f"reference must be 'neighbors', 'site_average' or a Structure, got {reference!r}")
    u = np.full((neighbors.n_nodes, 3), np.nan)
    have = rows >= 0
    df = config.frac[rows[have]] - ref_frac[have]
    if all(config.pbc):
        df, _ = config.lattice.minimum_image(df)
    u[have] = df @ config.lattice.matrix
    if mode == "site_average":
        site = neighbors.node_site.astype(str)
        for lab in set(site):
            m = (site == lab) & have
            if np.any(m):
                u[m] -= u[m].mean(axis=0)
    return u


def _direction(neighbors: NeighborList, direction, frame: str) -> np.ndarray:
    return neighbors._direction_unit(direction, frame)


@dataclasses.dataclass(frozen=True)
class DisplacementCorrelation:
    """Displacement statistics of one set of directed pairs.

    Attributes:
        n_pairs: pairs used (both ends occupied and of the selected species).
        n_skipped: pairs dropped because an end is vacant or of another species.
        n_unique_pairs, n_reversed: as for :class:`OccupationalCorrelation`.
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
    """
    n_pairs: int
    n_skipped: int
    n_unique_pairs: int
    n_reversed: int
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

    @property
    def correlation_matrix(self) -> np.ndarray:
        """(3, 3) ``C[a, b] / sqrt(Vo[a, a] Vt[b, b])`` for Cartesian components (NaN where a variance is 0)."""
        vo = np.diag(self.origin_covariance)
        vt = np.diag(self.target_covariance)
        with np.errstate(invalid="ignore", divide="ignore"):
            den = np.sqrt(np.outer(vo, vt))
            return np.where(den > 0, self.covariance / np.where(den > 0, den, 1), np.nan)


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
            matrices are computed and ``correlation`` is NaN.
        reference: see :func:`displacements`.
        centered: subtract the pair means (Pearson correlation, default) or not
            (DISCUS ``chem corr disp`` convention ``<xy>/sqrt(<x²><y²>)``).
        origin_species, target_species: keep only pairs whose ends carry these species.
    """
    u = displacements(neighbors, configuration, reference=reference)
    spec = node_species(neighbors, configuration).astype(str)
    keep = ~np.isnan(u[neighbors.i, 0]) & ~np.isnan(u[neighbors.j, 0])

    def _species_mask(value, nodes):
        if value is None:
            return np.ones(len(nodes), dtype=bool)
        vals = [value] if isinstance(value, str) else list(value)
        return np.isin(spec[nodes], [str(v) for v in vals])
    keep &= _species_mask(origin_species, neighbors.i) & _species_mask(target_species, neighbors.j)
    sub = neighbors._subset(keep)
    n = sub.n_pairs
    uo = u[sub.i]
    ut = u[sub.j]
    n_rev = _n_reversed(sub)
    nan3 = np.full(3, np.nan)
    nan33 = np.full((3, 3), np.nan)
    do = None if direction is None else _direction(neighbors, direction, direction_frame)
    dt = None
    if direction is not None:
        dt = do if target_direction is None else _direction(neighbors, target_direction, direction_frame)
    elif target_direction is not None:
        raise ValueError("target_direction needs direction")
    if n == 0:
        return DisplacementCorrelation(0, neighbors.n_pairs, 0, 0, centered, nan3, nan3, nan33, nan33, nan33,
                                       do, dt, math.nan, math.nan, "no pairs with both ends occupied")
    mo = uo.mean(axis=0) if centered else np.zeros(3)
    mt = ut.mean(axis=0) if centered else np.zeros(3)
    ao = uo - mo
    at = ut - mt
    cov = ao.T @ at / n
    cov_o = ao.T @ ao / n
    cov_t = at.T @ at / n
    corr, se, reason = math.nan, math.nan, None
    if do is None:
        reason = "no direction given"
    else:
        x = ao @ do
        y = at @ dt
        sxx, syy = float(x @ x / n), float(y @ y / n)
        if sxx <= _ZERO_VARIANCE or syy <= _ZERO_VARIANCE:
            reason = "zero variance of the projected displacements (rms < 1e-10 Å)"
        else:
            corr = float(x @ y / n / math.sqrt(sxx * syy))
            # Null model: independent displacements; a reversed pair contributes the
            # same-atom correlation of the two projections (rho_self).
            nodes = np.unique(np.concatenate([sub.i, sub.j]))
            un = u[nodes] - (u[nodes].mean(axis=0) if centered else 0.0)
            a1, a2 = un @ do, un @ dt
            den = math.sqrt(float(a1 @ a1) * float(a2 @ a2))
            rho_self = float(a1 @ a2) / den if den > 0 else 0.0
            se = math.sqrt((n + n_rev * rho_self ** 2) / n ** 2)
    return DisplacementCorrelation(
        n_pairs=n, n_skipped=neighbors.n_pairs - n, n_unique_pairs=n - n_rev // 2, n_reversed=n_rev,
        centered=centered, origin_mean=mo if centered else uo.mean(axis=0),
        target_mean=mt if centered else ut.mean(axis=0), covariance=cov, origin_covariance=cov_o,
        target_covariance=cov_t, origin_direction=do, target_direction=dt, correlation=corr,
        stderr_null=se, undefined=reason,
    )


def displacement_correlation_table(neighbors: NeighborList, configuration: Structure, *, direction,
                                   target_direction=None, direction_frame: str = "cartesian",
                                   by: str = "shell", reference="neighbors", centered: bool = True,
                                   origin_species=None, target_species=None,
                                   shell_tol: float = 1e-3) -> list[dict]:
    """One row per group (``by="shell"`` or ``by="vector"``, as in
    :func:`occupational_correlation_table`) with ``n_pairs``, ``correlation``,
    ``stderr_null`` and ``covariance`` (the projected ``<x y>`` in Å²)."""
    rows = []
    for key, mask in _group_labels(neighbors, by, shell_tol):
        res = displacement_correlation(
            neighbors._subset(mask), configuration, direction=direction, target_direction=target_direction,
            direction_frame=direction_frame, reference=reference, centered=centered,
            origin_species=origin_species, target_species=target_species)
        proj = (math.nan if res.origin_direction is None or res.n_pairs == 0
                else float(res.origin_direction @ res.covariance @ res.target_direction))
        rows.append({**key, "n_pairs": res.n_pairs, "n_skipped": res.n_skipped,
                     "correlation": res.correlation, "stderr_null": res.stderr_null, "covariance": proj,
                     "undefined": res.undefined})
    return rows
