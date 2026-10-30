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

**Overlap check.** Two atoms in different sites closer than 0.1 Å raise an error. Use the
minimum image on periodic structures. In a configuration, this applies to every such pair.
In an average model, it applies only when the two occupancies sum to more than 1. The
threshold is a documented keyword (`min_separation=0.1`, `None` to skip); the check is
skipped above 20 000 atoms.

## Periodicity

`pbc` is stored as given.

* **All True:** distances use `Lattice.minimum_image`.
* **All False:** distances are plain Cartesian; this is a molecule in a reference frame.
* **Mixed (slabs):** stored, but distance methods raise `NotImplementedError("mixed PBC not yet supported")`.

`wrapped()` is only valid on fully periodic structures. It returns positions in [0, 1)
and leaves `cell_index` unchanged. Wrapping changes coordinates but not physics.

## Identity rules

* `atom_id` never changes except through `renumbered()`.
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
* **Functional edits:** `replace(**fields)`, `with_frac(frac)`, `with_cartesian(cart)`,
  `with_occupancy(values)`, `remove(ids)`, `concat(other, *, renumber=False)`,
  `renumbered()`.
* **Distances:** `distance(id_a, id_b)`, `distance_matrix(ids=None)`, `displacement(id_a, id_b)`.
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
