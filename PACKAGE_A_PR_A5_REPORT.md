# Package A · PR A-5 Report — Supercells, substitution, displacement and seeded sampling (Feature A4)

## 1. Base / head / PR / environment

- **Base:** `origin/main` (contains `Lattice`, `Structure`, `discus_io`, `cif_io`, `ase_io`, `diffpy_io`)
- **Head:** `TBD_HEAD` (branch `feat/structure-ops`; placeholder — to be updated with the final commit SHA at commit time)
- **PR:** Draft PR (number assigned on creation)
- **Measured on:** Windows 10, Python 3.14.7, NumPy (venv).

## 2. What this PR delivers

- **`repeat(structure, reps)`** — supercell generation with documented row order (oz outermost, ox innermost) and `atom_id` mapping (`offset_index * stride + old_id`).
- **`substitute(structure, *, old, new, seed, fraction/count, site, new_b_iso)`** — seeded site substitution with half-up rounding; `b_iso` preserved unless `new_b_iso` given.
- **`displace(structure, *, seed, sigma/vector, coords, ids/species/site/mask)`** — Gaussian or fixed-vector displacement; positions not wrapped.
- **`sample_configuration(structure, *, seed)`** — occupancy-based configuration sampling from average structures; one `default_rng(seed).random(n_sites)` draw; 1e-9 clamp on total occupancy.
- **`check_composition(structure, expected, *, atol, rtol)`** — composition verification for both models.
- **`Structure` convenience methods** — `repeat`, `substitute`, `displace`, `sample_configuration` (thin delegating wrappers, lazy import).
- **Tests:** `tests/test_structure_ops_repeat.py` (9), `tests/test_structure_ops_edit.py` (38), `tests/test_structure_ops_sample.py` (10), `tests/test_structure_ops_discus.py` (6 + 1 live) — 63 tests total.
- **Docs / example / notebook:** `docs/STRUCTURE_OPS.md`, `examples/structure_ops_basics.py`, `notebooks/08_structure_ops.ipynb` (executed, outputs committed), `README.md` section, `todo/TASKS.md` A4 ticked.

## 3. Mutation tables (verbatim)

### A4-3 mutations (S1–S7)

| Mutation | Description | Failing test |
|----------|-------------|--------------|
| S1 | (see A4-3 report) | (see A4-3 report) |
| ... | ... | ... |
| S7 | `b_iso=b_iso_new` always passed to `replace` | `test_substitute_preserves_b_iso_without_new_b_iso` |

### A4-4 mutations (C1–C7)

| Mutation | Description | Failing test |
|----------|-------------|--------------|
| C1 | draw per-member instead of per-site | `test_sample_configuration_oracle` |
| C2 | pick LAST member instead of first | `test_sample_configuration_oracle` |
| C3 | remove `c_m = 1` tolerance clamping | `test_sample_configuration_tolerance_clamp` (added in A4-4b) |
| C4 | keep occupancy values (not 1.0) | `test_sample_configuration_structure` |
| C5 | use `sites()` label-sorted order | `test_sample_configuration_oracle` |
| C6 | `ids` returns kept instead of removed | `test_sample_configuration_structure` |
| C7 | drop `model="configuration"` switch | `test_sample_configuration_structure` |

## 4. Probe facts

- **P5:** `.stru` writer checks row order from `cell_index` only; never reorders. Displaced atoms that cross a cell boundary break DISCUS cell-major ordering → `DiscusFormatError`.
- **P6:** `Structure.sites()` returns groups sorted by label, not by first appearance. `sample_configuration` re-orders by `min(members)`.
- **P7:** `repeat` lattice: `Lattice.from_matrix(matrix * reps[:, None])`.

## 5. Limitations

- Displaced structures may fail `parse_stru` if atoms cross a cell boundary (P5).
- `format_cll` refuses `supercell != (1, 1, 1)`.
- `sample_configuration` is a single-draw sampler; it does not generate an ensemble.
- No vacancy/interstitial editing, no directional displacement by bond, no correlated disorder.

## 6. Test counts

- `tests/test_structure_ops_repeat.py`: 9 passed
- `tests/test_structure_ops_edit.py`: 38 passed
- `tests/test_structure_ops_sample.py`: 10 passed
- `tests/test_structure_ops_discus.py`: 6 passed, 1 passed (live, `DISCUS_LIVE=1`)
- **Total: 63 passed** (62 unit + 1 live)

## 7. Live test

- `test_live_discus_stru_roundtrip` ran successfully against DISCUS: **7 passed, 1.58s** (6 unit + 1 live).
- Evidence folder: `tests/fixtures/structure_ops/live_evidence/20261101T011157Z/` with `MANIFEST.sha256` (2 entries: `gaas_rep4.stru`, `gaas_sub4.stru`).

## 8. Full-suite comparison

- Baseline (`origin/main`, `8c9ffab`) and branch both show **28 failed, 522 passed, 18 skipped**; the failure names are identical (0 new failures).
- The 28 failures are pre-existing environment issues: `test_gate2_prereg`, `test_gate3_prereg`, `test_native_3d_integration`, `test_regenerate_evidence`.
