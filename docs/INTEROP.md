# Interop: CIF, ASE and DiffPy adapters

The `discuspdf` package can read crystallographic CIF files (via `gemmi`) and
convert between its `Structure` type and the `Atoms` (ASE) and
`diffpy.structure.Structure` types.  These adapters are optional extras:

```
pip install 'discus-python[cif]'      # gemmi (CIF reading)
pip install 'discus-python[ase]'     # ASE
pip install 'discus-python[diffpy]'  # DiffPy
pip install 'discus-python[interop]' # all three
```

`import discuspdf` works without any of them; the adapters raise a helpful
`ImportError` naming the extra when the library is missing.

Tested with gemmi 0.7.5, ase 3.29.0, diffpy.structure 3.2.3 (CI) and 3.5.0 (local).

## Quick start

### CIF (gemmi)

```python
from discuspdf.cif_io import read_cif

imp = read_cif("structure.cif")          # CifImport
s = imp.structure                        # Structure (model="average")
print(imp.spacegroup, imp.symmetry_source, imp.loss)
```

`read_cif` expands the asymmetric unit by the space group (symop loop,
space-group symbol, or assumed P1) and wraps coordinates to [0, 1).
A structure written with `write_cif` carries `_discuspdf_*` tags and reads
back exactly (model, supercell, site labels, atom ids, cell indices).

### ASE

```python
from discuspdf.ase_io import to_ase, from_ase

atoms = to_ase(s)                       # ase.Atoms (dp_* arrays attached)
s2 = from_ase(atoms)                   # exact round trip (path A)

# Foreign Atoms object (no dp_* arrays):
s3 = from_ase(foreign_atoms, model="configuration")   # path B
```

`to_ase` stores the full DISCUS state in `atoms.arrays` (`dp_species`,
`dp_site_label`, `dp_atom_id`, `dp_occupancy`, `dp_cell_index`, `dp_b_iso`)
and `atoms.info` (`dp_model`, `dp_supercell`).  `from_ase` restores it
exactly when those keys are present (path A).  For a foreign `Atoms`
object (path B) you must pass `model=`; occupancies are assumed 1.0 and
unsupported features (constraints, momenta, magnetic moments, charges, tags,
calculator, non-occupancy info keys) are reported as losses.  Note that
`from_ase` returns a plain `Structure`: the foreign-Atoms defaults are NOT
reported in the returned object — they are only surfaced as losses/warnings
at conversion time.

### DiffPy

```python
from discuspdf.diffpy_io import to_diffpy, from_diffpy

ds = to_diffpy(s)                      # diffpy.structure.Structure
s2 = from_diffpy(ds, model=s.model)    # model must be passed again
```

The model is not stored in the DiffPy object; pass it again to
`from_diffpy`.  Orientation is NOT preserved (DiffPy defines its own
Cartesian frame); fractional coordinates, lattice parameters and all
distances/angles are.

## What is preserved and what is lost

| Quantity | CIF | ASE | DiffPy |
|---|---|---|---|
| species (element) | exact | exact | exact |
| site labels | exact (via `_discuspdf_site_label`) | exact (`dp_site_label`) | lost (DiffPy labels are regenerated) |
| atom ids | exact (via `_discuspdf_atom_id`) | exact (`dp_atom_id`) | lost |
| cell indices / supercell | exact (via `_discuspdf_cell_index_*`) | exact (`dp_cell_index`) | lost |
| pbc | all True (assumed) | exact | lost |
| fractional coordinates | 10 decimals | exact (in memory) | exact (in memory) |
| lattice parameters | 10 decimals | exact (in memory) | exact (in memory) |
| occupancy | exact | exact (in `dp_occupancy`; ASE consumers see 1.0) | exact |
| B factor (isotropic) | exact | exact (`dp_b_iso`) | via Uisoequiv (B = 8π²U) |
| anisotropic ADPs | **loss** | **loss** | **loss** (isotropic equivalent kept) |
| model | exact (via `_discuspdf_model`) | exact (`dp_model`) | **not stored** — pass to `from_diffpy` |
| species charges | **loss** | **loss** | **loss** (element only) |

## Symmetry rules (CIF reading)

1. `_space_group_symop_operation_xyz` or `_symmetry_equiv_pos_as_xyz` present
   → use the symop loop (`symmetry_source="symop_loop"`).
2. Otherwise, space-group symbol present → use the symbol
   (`symmetry_source="space_group_symbol"`).
3. Neither → assume P1 (`symmetry_source="assumed_P1"`); recorded in
   `CifImport.defaults_applied`.

The asymmetric unit is expanded by all symmetry operations, wrapped to [0, 1),
and periodic images are merged within 1e-4 (fractional, per component,
minimum-image).  Modulated structures are rejected with `InteropError`.

## Occupancy and model mapping

* CIF reading defaults to `model="average"` (site-occupancy model).
* `model="configuration"` requires all occupancies = 1.0 (checked on read and
  on write).
* A mixed site is several rows with the same `(site_label, cell_index)`;
  occupancies at a position sum to ≤ 1 + 1e-6 (checked).
* ASE consumers always see fully occupied atoms; the true occupancies live in
  the `dp_occupancy` array.  `to_ase` emits an `InteropWarning` for an
  average-model structure with partial occupancy or a mixed site.
* `from_ase` path B (foreign Atoms) assumes occupancy 1.0.

## Loss policy (`allow_loss`)

Every adapter takes `allow_loss: bool = False`.  Losses are accumulated in a
`LossReport`:

* `allow_loss=False` (default) → `InteropError` listing every loss (one per line).
* `allow_loss=True` → one `InteropWarning` carrying the lines, and the
  conversion proceeds.

`to_ase` is the one exception: it never raises for a loss (design rule 2);
it always succeeds and at most warns.

## `_discuspdf_*` tags (CIF exact round-trip)

`write_cif` emits, in addition to the standard tags:

* `_discuspdf_model` — `average` or `configuration`
* `_discuspdf_supercell_a/b/c` — the supercell
* `_discuspdf_site_label` (loop) — per-atom site labels
* `_discuspdf_atom_id` (loop) — per-atom ids
* `_discuspdf_cell_index_a/b/c` (loop) — per-atom cell indices

When all of these are present, `read_cif` restores the structure exactly
(model, supercell, site labels, atom ids, cell indices) instead of expanding
by symmetry.  The writer is atomic and verifies the written file by parsing
it back; on mismatch the target is left untouched.

## B ↔ U

`B = 8π²·U` (Å²).  `gemmi` has no `u_to_b`; the conversion is implemented in
`discuspdf.interop` (`u_to_b`, `b_to_u`, `B_FROM_U`).  CIF reading prefers
`_atom_site_B_iso_or_equiv` over `_atom_site_U_iso_or_equiv`; if all U (or B)
values are 0, `b_iso` is `None`.

## Tolerances (predeclared)

| Quantity | Tolerance |
|---|---|
| lattice parameters | rtol 1e-8 |
| fractional coordinates (CIF round trip) | atol 1e-9 (10-decimal file precision) |
| fractional coordinates (ASE, DiffPy, in memory) | atol 1e-12 |
| occupancy | atol 1e-9 |
| B factor | atol 1e-6 |
| distances | atol 1e-6 Å |
| species, site labels, atom ids, cell indices, row order | exact |

## Limitations

* Anisotropic ADPs are a loss everywhere (isotropic equivalent kept, if any).
* Disorder groups, magnetic moments: CIF reading reports them as losses.
* Total occupancy > 1 at a position is rejected (`InteropError`).
* Mixed pbc (not all True) is a loss on CIF writing (the reader assumes all
  periodic).
* ASE `info['occupancy']` is never interpreted (overseer ruling): foreign
  partial occupancies cannot be mapped to atoms.
* DiffPy: orientation is not preserved; site labels, atom ids, cell indices,
  supercell and pbc are not stored; the model must be passed to `from_diffpy`.
* The element table covers the first 100 elements; unknown type symbols raise
  `InteropError`.
