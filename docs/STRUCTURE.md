# Structure model (Package A, Feature A1)

`discuspdf.Structure` is an immutable atomic-structure value type. It holds per-atom
NumPy arrays (no per-atom Python objects) plus a `Lattice`, and supports two models:

* **`"average"`** — a site-occupancy model (as in a CIF). Mixed sites are several rows
  with the same `(site_label, cell_index)`; occupancies sum to ≤ 1.
* **`"configuration"`** — a realised atomistic arrangement (as in a DISCUS supercell).
  Every occupancy is exactly 1.0; each site holds at most one atom.

## Fields

| Field | Type / shape | Default | Meaning |
|---|---|---|---|
| `lattice` | `Lattice` | required | frame for fractional coordinates |
| `species` | `(N,)` str | required | label as given (e.g. `"Ga"`, `"O2-"`); non-empty, no whitespace, ≤ 8 chars |
| `frac` | `(N, 3)` float | required | fractional coordinates, **stored as given** (never auto-wrapped) |
| `occupancy` | `(N,)` float | 1.0 | in (0, 1] |
| `site_label` | `(N,)` str | species + 1-based index | crystallographic site the atom belongs to |
| `atom_id` | `(N,)` int64 | 0…N−1 | unique, ≥ 0, stable through edits |
| `b_iso` | `(N,)` float or None | None | isotropic B in Å² (DISCUS uses B, not U); ≥ 0 if given |
| `cell_index` | `(N, 3)` int | zeros | unit-cell index inside a supercell; 0 ≤ index < `supercell` |
| `supercell` | `(3,)` int | (1, 1, 1) | repetitions relative to the parent cell |
| `pbc` | 3 bools | (True, True, True) | periodic axes |
| `model` | `"average"` or `"configuration"` | required, no default | see below |

All stored arrays are read-only (`flags.writeable = False`). `__setattr__` and
`__delattr__` raise. Getters return read-only views or copies (documented per property).

## Average model vs configuration

**`"average"`** — a site-occupancy model, as in a CIF. Mixed sites are several rows with
the same `(site_label, cell_index)`. All rows of one site must share the same fractional
position (≤ 1e-6 per component) and the same `b_iso`. Their occupancies must sum to
≤ 1 + 1e-6. A sum below 1 means vacancies.

**`"configuration"`** — a realised atomistic arrangement, as in a DISCUS supercell.
Every occupancy must be exactly 1.0, and each `(site_label, cell_index)` holds at most
one atom. Displaced atoms keep their site.

**No implicit conversion.** Nothing in this PR turns an average model into a configuration;
that is A4's explicit, seeded sampling step. Attempts such as a configuration with
occupancy 0.5 raise `ValueError` with a message pointing to the average model.

**Overlap check.** Two atoms in different sites closer than `min_separation` (default
0.1 Å) raise `ValueError`, naming both sites and the distance. Periodic structures use the
minimum image. In a configuration, this applies to every such pair. In an average model,
it applies only when the two sites' total occupancies add up to more than 1, so split
sites (e.g. two half-occupied positions 0.05 Å apart) are allowed.

The check runs on construction, on every edit (`replace`, `with_frac`, `with_cartesian`,
`with_occupancy`, `wrapped`, `renumbered`) and on `concat`. It is never skipped
automatically, whatever the size; `min_separation=None` is the only way to turn it off,
and each of those methods accepts it. Structures with mixed `pbc` must pass
`min_separation=None`, because mixed-PBC distances are not implemented yet.

*Method.* A cell list: positions are binned so that every bin is at least `min_separation`
thick in each lattice direction (using the interplanar spacings `1/|b_i|`), so any pair
closer than the cutoff lies in the same or an adjacent bin, also in strongly skewed
cells. When every axis has at least 3 bins, the nearest rounded image is the only one
that can be within the cutoff; otherwise (cells not much larger than the cutoff) candidate
pairs go through `Lattice.minimum_image`. The number of bins is capped at 4 per atom.
A test compares it with a brute-force all-pairs, all-images check on random structures.
Typical cost (cloud runner): a 20 × 20 × 20 GaAs supercell (64 000 atoms) validates in
about 0.4 s, and 10⁵ random atoms in about 0.6 s.

## Periodicity

`pbc` is stored as given.

* **All True:** distances use `Lattice.minimum_image`.
* **All False:** distances are plain Cartesian; this is a molecule in a reference frame.
* **Mixed (slabs):** stored, but distance methods raise `NotImplementedError("mixed PBC not yet supported")`.

`wrapped()` is only valid on fully periodic structures. It returns positions in [0, 1)
and leaves `cell_index` unchanged. Wrapping changes coordinates but not physics.

## Identity rules

* `atom_id` never changes except through `renumbered()`.
* `ids=` arguments (`select`, `remove`, `distance`, `displacement`, `distance_matrix`)
  always mean `atom_id` values, never row numbers; unknown ids raise `ValueError`.
* `select`, `remove` and `with_*` keep row order.
* New ids, needed only by `concat(renumber=True)`, start at `max(existing) + 1`.
* No method sorts or shuffles rows.

## Units

* Lengths in Å.
* `b_iso` in Å² (DISCUS uses B, not U).

## API

* **Construction:** `Structure(lattice, species, frac, *, model, ...)` and
  `Structure.from_cartesian(lattice, species, cart, *, model, ...)`.
* **Views:** `n_atoms`, `cartesian` (Å), `frac`, `sites()`, `atomic_numbers`,
  `composition()`, `counts()`, `formula()`.
* **Selection:** `select(mask=None, *, ids=None, species=None, site=None)`.
* **Functional edits** (each returns a new, fully re-validated `Structure`):
  `replace(*, min_separation=0.1, **fields)`, `with_frac(frac)`, `with_cartesian(cart)`,
  `with_occupancy(values)`, `remove(ids)`, `renumbered()`, and
  `concat(other, *, renumber=False, min_separation=0.1)`. `concat` requires matching
  lattices (`isclose`), model, `supercell` and `pbc`, and `b_iso` on both or neither.
* **Distances:** `distance(id_a, id_b)`, `displacement(id_a, id_b)` and
  `distance_matrix(ids=None)` (vectorised; refuses more than 5 000 atoms — use the
  Package C neighbour lists for larger sets; about 3.5 s for 2 000 atoms on the cloud runner).
* **Wrapping:** `wrapped()`.
* **Comparison:** `__eq__` (exact), `isclose(other, rtol, atol)`, `__hash__ = None`.

## Not yet (later PRs)

* CIF, ASE or DiffPy adapters (A2).
* `.cll`/`.stru` reading or writing (A3).
* `repeat`, substitution, random/directional displacements, occupancy sampling (A4).
* Neighbour lists and correlations (Package C).
* Magnetic moments (Package D11), anisotropic ADPs, mixed-PBC distances.

## Example

`examples/structure_basics.py` (run with `python examples/structure_basics.py`); covered by
`tests/test_structure.py::test_example_script_runs`.

Notebook: [`notebooks/04_structure_basics.ipynb`](../notebooks/04_structure_basics.ipynb) —
worked examples of the API (construction, validation, select/edit, distances, wrapping).
