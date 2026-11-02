# Package C engine probe — DISCUS `chem corr` vs `discuspdf.correlations` (DISCUS_SUITE 6.23.01)

**Result: PASS — 29 of 29 correlation values agree to 4.7 × 10⁻¹³ (DISCUS prints 12 significant
digits), and all 29 pair counts agree exactly.**

| | |
|---|---|
| Registration | inputs, macros and `predictions.json` committed in `9c62d2a` **before** the engine run |
| Engine run | owner's laptop, WSL Ubuntu-22.04, `/usr/local/bin/discus_suite` SHA-256 `d4274f13…45706f3` (same binary as the Package B probes), 2026-11-01 18:36 UTC, local agent task 002 |
| Evidence | `evidence/20261101T183654Z_LAPTOP-NGU3BT0N/` (three logs, `engine.txt`, `exit_codes.txt`, `input_sha256.txt`) |
| Input hashes | 6/6 match `predictions.json` |
| Check | `python benchmarks/discus_6_23_01/c_probes/analyze.py` (also run by `tests/test_c_probes.py`) |

## What was compared

| Crystal | Model | DISCUS command | Neighborhoods |
|---|---|---|---|
| `occ_bcc` | bcc 8³, 1024 atoms, Au/Cu, partial B2 order (site 1: 75 % Au, site 2: 30 % Au) | `corr occ, Au, Cu` and `corr occ, Cu, Au` | 5 `set vec` (two-site and same-site, signed translations), 2 `set neig, dis` shells |
| `disp_sc` | sc 12×6×6, Cu, longitudinal cosine wave + transverse sine wave + noise | `corr disp, Cu, Cu` | 6 vectors, incl. the cross term x(origin)·y(target) |
| `mixed_bcc` | bcc 6³, random Au/Cu, random displacements, Cu shifted along +c | `corr disp, Au, Cu` / `Cu, Cu` / `all, all` | 3 vectors along [001] and [110] |

Our side: `NeighborList.select(origin_site=is1, target_site=is2, translation=t, symmetric=(is1 != is2))`
for `set vec`; `select(distance=(d − fsig/2, d + fsig/2))` for `set neig, dis`;
`OccupationalCorrelation.discus(a, b)`; `displacement_correlation(..., reference="site_average",
centered=False, direction_frame="fractional")`.

## Findings

1. **Definitions confirmed.** DISCUS's occupational correlation is
   `(N_aa/N − θ²)/(θ(1 − θ))` with θ pooled over both ends, over pairs with both ends in {a, b};
   the result does not depend on the order of the two species. The displacement correlation is
   the **uncentred** `⟨x_i x_j⟩/√(⟨x_i²⟩⟨x_j²⟩)` with displacements from the per-site average
   position, directions in lattice coordinates normalised with the metric.
2. **Pair sets confirmed.** `set vec, n, is1, is2, t` uses both directions for two different sites
   and only `+t` per atom for the same site (512 = one pair per atom of the site, not 1024);
   `set neig, dis` uses every atom's neighbors at that distance (8 × 1024 and 6 × 1024 pairs).
3. **Uncentred vs Pearson.** With species filters the DISCUS value can differ strongly from the
   Pearson correlation (e.g. Au→Cu along [001]: −0.237 vs −0.128; Cu→Cu same-site: +0.139 vs
   +0.004), because the per-**site** average does not remove the mean +c shift of the Cu atoms.
   Both are reported by `displacement_correlation` (`centered=True/False`); they agree when the
   pair means vanish (`all, all`: −0.1495 vs −0.1503).
4. **Engine crash at exit.** `occ_bcc` ended with `free(): invalid pointer` / SIGABRT
   (exit 134) **after** `exit` from the suite level; every result had already been printed.
   The two displacement runs exited with 0. Candidate for an upstream DiffuseCode report
   (`chem corr occ` followed by `exit`, DISCUS 6.23.01).

## Scope and limits

One engine version (6.23.01), orthogonal cells, one or two sites per cell, no vacancies (DISCUS
atom order requires every site to be present), atom mode only (`mode mole` not tested),
`set mode, quick, period`. Not compared: `corr` with `lots`, angle neighborhoods, `apply_symmetry`.
