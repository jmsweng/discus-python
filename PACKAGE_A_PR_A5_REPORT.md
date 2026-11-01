# Package A · PR A-5 Report — Supercells, substitution, displacement and seeded sampling (Feature A4)

## 1. Base / head / PR / environment

- **Base:** `origin/main` (contains `Lattice`, `Structure`, `discus_io`, `cif_io`, `ase_io`, `diffpy_io`)
- **Head:** `dd32bf8` (branch `feat/structure-ops`)
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

- `test_live_discus_stru_roundtrip` (marker `@pytest.mark.discus_live`, scoped to the single live test only) ran successfully against DISCUS: **1 passed, 6 deselected in 1.36s**.
- The test writes two `.stru` files (4×4×4 GaAs parent and 25% Ga→In substituted), runs DISCUS `read stru` → `eval n[1]` → `save back_<name>` on each, then reads DISCUS's own saved file with `parse_stru(text, supercell=(4,4,4))` and compares atom count, species counts and centred periodic coordinates.
- DISCUS `n[1]` values: `gaas_rep4.stru` → 512, `gaas_sub4.stru` → 512 (both match Python).
- Species counts (DISCUS-saved, upper-cased): `gaas_rep4.stru` → GA: 256, AS: 256; `gaas_sub4.stru` → GA: 192, AS: 256, IN: 64 (25% of 256 Ga → 64 In).
- Mutation check: substituting `py_counts` with the un-substituted `rep` counts causes the test to **fail** with `DISCUS species counts {'GA': 192, 'AS': 256, 'IN': 64} != Python {'GA': 256, 'AS': 256}`.
- Evidence folder: `tests/fixtures/structure_ops/live_evidence/20261101T0156Z/` with `MANIFEST.sha256` (2 entries: the two `.log` files in the folder; `--check` → all ok).
- Previous evidence folders `20261101T011157Z` and `20261101T0143Z` removed (superseded).

## 8. Full-suite comparison

- Baseline (`origin/main`, `8c9ffab`) and branch both show **28 failed, 522 passed, 18 skipped**; the failure names are identical (0 new failures).
- The 28 failures are pre-existing environment issues: `test_gate2_prereg`, `test_gate3_prereg`, `test_native_3d_integration`, `test_regenerate_evidence`.
