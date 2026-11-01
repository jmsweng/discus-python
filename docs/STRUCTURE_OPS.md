# Structure editing operations (Package A, Feature A4)

## Overview

`discuspdf.structure_ops` provides pure functions for editing `Structure` objects:
repeating into supercells, substituting species, displacing atoms, and sampling
configurations from average (site-occupancy) structures. All functions return an
`EditResult` dataclass; the input `Structure` is never mutated.

## API

### `repeat(structure, reps, *, min_separation=0.1) -> EditResult`

Repeat a structure into a supercell of size `reps = (nx, ny, nz)`.

- **Row order:** offsets in loop order `oz` outermost, `oy`, `ox` innermost (x fastest); old rows in old order inside each offset.
- **atom_id mapping:** `stride = max(old atom_id) + 1`; `atom_id_new = offset_index * stride + atom_id_old` with `offset_index = ox + nx*(oy + ny*oz)`. The inverse is `old_id = new % stride`, `offset_index = new // stride`.
- **Lattice:** `Lattice.from_matrix(old_matrix * reps[:, None])` (rows scaled).
- **PBC:** requires all `pbc` True; otherwise `ValueError`.
- **`ids`:** empty (no atoms changed).
- **`info`:** `{"reps": (nx, ny, nz), "stride": stride}`.

### `substitute(structure, *, old, new, seed, fraction=None, count=None, site=None, new_b_iso=None) -> EditResult`

Substitute some atoms of species `old` with species `new`.

- Only `model="configuration"` (average → `ValueError`).
- Exactly one of `fraction` (0 ≤ f ≤ 1) / `count` (int ≥ 0).
- `n = count` or `floor(fraction * len(candidates) + 0.5)` (half-up rounding).
- Chosen rows: `np.random.default_rng(seed).choice(candidate_rows, size=n, replace=False)`, sorted.
- `b_iso` of chosen rows set to `new_b_iso` if given (requires the structure to have `b_iso`).
- **`ids`:** atom_ids of chosen rows.
- **`info`:** `{"n_candidates", "n_changed", "old", "new"}`.

### `displace(structure, *, seed=None, sigma=None, vector=None, coords="cartesian", ids=None, species=None, site=None, mask=None, min_separation=0.1) -> EditResult`

Displace selected atoms by a Gaussian (`sigma`) or a fixed `vector`.

- Only `model="configuration"`.
- Exactly one of `sigma` (> 0 Å; `seed` required) or `vector` (3 numbers; `seed` must be `None`).
- `coords`: `"cartesian"` (Å) or `"fractional"` (of the structure's lattice).
- Displaced positions are **not wrapped** (no minimum-image applied).
- **`ids`:** selected atom_ids.
- **`info`:** `{"n_displaced", "sigma", "vector"}`.

### `sample_configuration(structure, *, seed) -> EditResult`

Draw one realised configuration from an average (site-occupancy) structure.

- Only `model="average"` (configuration → `ValueError`).
- Sites: `Structure.sites()` groups, re-ordered by their smallest member row index (order of first appearance).
- Draw: `u = np.random.default_rng(seed).random(n_sites)` once, in site order.
- Per site: cumulative occupancies; if `c_m >= 1 - 1e-9`, set `c_m = 1`. First member with `u_k < c_j` is realised; if none, the site is a vacancy.
- Result: `model="configuration"`, `occupancy = 1.0`, original `atom_id`s and row order kept.
- **`ids`:** atom_ids of **removed** rows (vacancies and rejected alternatives).
- **`info`:** `{"n_sites", "n_vacancies", "expected", "realised"}`.

### `check_composition(structure, expected, *, atol=0.0, rtol=0.0) -> None`

Check that the structure's composition matches `expected` (species → count).

## Seed contract

- Every random operation takes an explicit integer `seed` (no default, no global RNG).
- Same `seed` ⇒ identical result (reproducibility).
- Different `seed` ⇒ different result (independence).
- **Same NumPy version ⇒ same result.** NumPy does not guarantee stream stability across versions; the seed contract is version-specific.

## Conventions

- `Structure` is immutable; all operations return a new `Structure`.
- `atom_id`s are never changed except as stated for `repeat`.
- `b_iso` is preserved unless `new_b_iso` is given (substitute) or the model changes.
- Displaced positions are not wrapped; `cell_index`, `site_label`, `atom_id`, order unchanged.
- `sample_configuration` converts `model="average"` → `model="configuration"`, sets `occupancy = 1.0`.

## DISCUS round-trip limitations

- `format_stru` / `parse_stru` round-trip works for repeated and substituted structures.
- **Displaced structures may fail `parse_stru`** if atoms cross a cell boundary: DISCUS `.stru` files require cell-major atom ordering (probe P5), and displaced atoms that move across the cell boundary break this ordering. The `DiscusFormatError` is raised in this case.
- `format_cll` works only for unit-cell parents (`supercell == (1, 1, 1)`); it refuses supercells.

## What is not included

- Vacancies/interstitials editing
- Directional displacement by bond
- Correlated disorder
- Anything on `DiscusVolume`
