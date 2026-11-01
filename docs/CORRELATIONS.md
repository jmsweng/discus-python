# Occupational and displacement correlations (Package C, Features 13 and 14)

`discuspdf.correlations` measures short-range order on the pairs of a
[`NeighborList`](NEIGHBORS.md). It is independent of any scattering calculation.

```python
from discuspdf import neighbor_list
from discuspdf.correlations import (occupational_correlation, occupational_correlation_table,
                                    displacement_correlation, displacement_correlation_table,
                                    concentrations)

nl = neighbor_list(reference, 8.0)                        # reference: undisplaced sites
sro = occupational_correlation(nl.select(shell=1), configuration)
sro.value("warren_cowley", "Au", "Cu"), sro.value("stderr_warren_cowley", "Au", "Cu")
occupational_correlation_table(nl, configuration, by="shell")      # shell-averaged
occupational_correlation_table(nl, configuration, by="vector")     # one row per signed vector

disp = displacement_correlation(nl.select(translation=(1, 0, 0), symmetric=True), configuration,
                                direction=(1, 0, 0), direction_frame="fractional")
disp.correlation, disp.stderr_null, disp.covariance
```

## Inputs and bookkeeping

* **Pairs** come from a `NeighborList`. Build it on the *reference* structure (ideal positions,
  or the average-model structure) so that shells and vectors are sharp, then select the pairs
  you want: a shell (`select(shell=k)`), a signed neighbor vector
  (`select(origin_site=…, target_site=…, translation=…)`), a direction, species, …
* **Species and positions** come from a `"configuration"` `Structure` with the same lattice,
  supercell and pbc (checked). Atoms are matched to the list's nodes by
  (`site_label`, `cell_index`); a node with no atom is a **vacancy** (label `"Va"` by default,
  `vacancy=None` makes it an error). If no configuration is passed, the list's own structure
  is used (it must then be a configuration).
* `concentrations(nl, configuration)` gives the per-site-label (**site-specific**) and the
  overall (**global**) concentrations over the nodes. They are reported separately from the
  **pair-weighted** marginals that enter the formulas below; the two coincide only when every
  node takes part in the same number of pairs.

## Occupational correlations (Feature 13)

For a set of directed pairs (origin species *a*, target species *b*), with `N[a,b]` counts and
`n = Σ N`:

| Quantity | Definition | Domain |
|---|---|---|
| `joint` | `p_ab = N_ab / n` | [0, 1] |
| `origin_concentration`, `target_concentration` | `c_o(a) = Σ_b p_ab`, `c_t(b) = Σ_a p_ab` (pair-weighted) | [0, 1] |
| `conditional` | `P(b|a) = N_ab / Σ_b N_ab` | [0, 1] |
| `warren_cowley` | `α_ab = (P(b|a) − c_t(b)) / (δ_ab − c_t(b))` (Cowley's multicomponent form) | `α_ab ≤ 1`; for a ≠ b `α_ab ≥ 1 − 1/c_t(b)` |
| `pearson` | `r_ab = (p_ab − c_o(a) c_t(b)) / √(c_o(a)(1−c_o(a)) c_t(b)(1−c_t(b)))` — the correlation coefficient of the occupation indicators at the two ends | [−1, 1] |

Signs: α < 0 for *a ≠ b* means more unlike neighbors than at random (ordering); α > 0 for
*a = a* means clustering. For a binary alloy on one sublattice with both directions of every
pair in the list, `α_AB = r_AA = −r_AB`, the classical Warren–Cowley parameter.

Checked values: B2 order on bcc gives α = −1 in the first shell and +1 in the second; L1₂
Cu₃Au gives α(Au,Cu) = −1/3 in the first shell; random occupation gives values within the
standard error (tests in `tests/test_correlations.py`).

**Undefined values** are NaN, with the reason in `result.undefined[(a, b)]`:
no pairs; species *a* absent at the origins; species *b* absent at the targets; pure end
member (zero variance of an indicator — `pearson` is then undefined, while `α_ab` stays
defined where its own denominator is nonzero, e.g. `α = 0` for a sublattice entirely occupied by
*b*).

`species=(…)` restricts the pairs to those with **both** ends among the given species (as
DISCUS does). By default all species present at the list's nodes are used, so tables have the
same axes in every group.

### Standard errors

`stderr_pearson` and `stderr_warren_cowley` are the standard errors under the **null model of
independent random occupation** with the observed concentrations. For centred indicators the
products over different pairs are uncorrelated, so exactly

`Var(r_ab) = (n + n_rev ρ_ab²) / n²`,

where `n_rev` is the number of pairs whose reverse direction is also in the list and
`ρ_ab = (δ_ab c_a − c_a c_b)/√(c_a(1−c_a) c_b(1−c_b))` is the correlation of the *a*- and
*b*-indicators at one site (pooled concentrations of the nodes involved). For a list holding both
directions of every pair and a = b this is `1/√n_unique_pairs`. `stderr_warren_cowley` follows
from `α = r · √(c_o(1−c_o) c_t(1−c_t)) / (c_o (δ − c_t))`. The tests check the formula against the
spread over 300 random crystals (agreement within 15 %). The error describes the deviation of a
random crystal from zero; it is not a confidence interval for an ordered crystal, and it is
approximate when sublattices with different concentrations are pooled.

### DISCUS `chem corr occ`

`result.discus(a, b)` reproduces DISCUS 6.23.01 (`chem.f90`, `chem_corr_occ`):
pairs restricted to {a, b}, `θ = (2 N_aa + N_ab + N_ba) / (2N)`,
`c = (N_aa/N − θ²) / (θ(1 − θ))`. DISCUS reports 0 where θ is 0 or 1; we return NaN. To get
DISCUS's pairs for `set vec, n, is1, is2, t`, select
`nl.select(origin_site=is1, target_site=is2, translation=t, symmetric=(is1 != is2))`
(DISCUS takes both directions for two different sites, only `+t` for the same site). On a list
that holds both directions of every pair, `discus(a, b) = pearson[a, a]`.
The engine comparison is in `benchmarks/discus_6_23_01/c_probes/`.

## Displacement correlations (Feature 14)

`displacements(nl, configuration, reference=…)` gives, for every node, the Cartesian
displacement `u` (Å) of its atom from the reference position (NaN for vacancies):

* `reference="neighbors"` (default): the node positions of the list's own structure — build the
  list on the undisplaced reference;
* `reference="site_average"`: in addition, the mean displacement of each site label over all
  cells is subtracted (the DISCUS `chem` convention: displacement from the average position of
  the site);
* a `Structure` (configuration with the same sites): its positions.

The difference is taken as the **minimum image** in the supercell, so wrapping the
configuration into [0, 1) does not matter; displacements longer than half the shortest
supercell spacing would be ambiguous.

`displacement_correlation(pairs, configuration, direction=d_o, target_direction=d_t, …)`:

| Quantity | Definition |
|---|---|
| `covariance` | `C = ⟨(u_o − m_o)(u_t − m_t)ᵀ⟩` over the pairs (3×3, Å²; Cartesian) |
| `origin_covariance`, `target_covariance` | `⟨(u − m)(u − m)ᵀ⟩` of the origins / targets |
| `correlation_matrix` | `C[α,β] / √(V_o[α,α] V_t[β,β])` |
| `correlation` | `⟨x y⟩ / √(⟨x²⟩⟨y²⟩)` with `x = (u_o − m_o)·d̂_o`, `y = (u_t − m_t)·d̂_t` |

`centered=True` (default) subtracts the pair means `m_o`, `m_t` (Pearson correlation).
`centered=False` uses `m = 0`, which with `reference="site_average"` is the DISCUS
`chem corr disp` definition `⟨x(i)x(j)⟩/√(⟨x(i)²⟩⟨x(j)²⟩)`. Directions are Cartesian or, with
`direction_frame="fractional"`, `[uvw]` of the parent cell normalised with the metric (as DISCUS
does); `target_direction` defaults to `direction`, so cross-correlations (e.g. *x* at the
origin with *y* at the target) are available. `origin_species`/`target_species` select
cross-species pairs; pairs with a vacant end are skipped (`n_skipped`).

Checked values: a longitudinal wave `u_x = A cos(2π q·R)` gives `cos(2π q·t)` for the vector *t*;
a transverse sine adds the cross term `sin(2π q·t)` for one direction of *t*, which cancels when
both directions are pooled; independent random displacements give values within the standard
error; the covariance equals an explicit loop over the pairs.

**Undefined:** `correlation` is NaN when no direction is given, when no pair has two occupied
ends, or when the projected displacements have zero variance (rms below 1e-10 Å, i.e. round-off),
with the reason in `undefined`.

**Standard error** (`stderr_null`): independent displacements, `Var = (n + n_rev ρ_self²)/n²`
with `ρ_self` the correlation of the two projections on the same atom.

## Tables

`occupational_correlation_table` and `displacement_correlation_table` return one row (a plain
dict, ready for `pandas.DataFrame`) per group: `by="shell"` (distance shells, shell-averaged
values) or `by="vector"` (origin site, target site, translation — the signed neighbor
vectors). Occupational rows hold every (a, b) combination with counts, probabilities, α, r and
both standard errors; displacement rows hold the pair count, correlation, standard error and the
projected covariance.

## Not covered

Molecules; magnetic moments; correlations between different neighbor sets (three-body);
ensemble averages over several configurations (average the returned counts/covariances
yourself, weighted by `n_pairs`).
