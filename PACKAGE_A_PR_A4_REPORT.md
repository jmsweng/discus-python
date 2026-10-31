# Package A · PR A-4 Report — CIF, ASE and DiffPy interoperability (Feature A2)

## 1. Base / head / PR / environment

- **Base:** `a6e303f` (`origin/main`, contains `Lattice`, `Structure`, `discus_io`)
- **Head:** `19779a8` (branch `feat/interop-adapters`)
- **PR:** [Draft PR #15](https://github.com/jmsweng/discus-python/pull/15)
- **Measured on:** Windows 10, Python 3.14.7, gemmi 0.7.5, ase 3.29.0, diffpy.structure 3.5.0.
  CI (GitHub Actions, ubuntu-latest, Python 3.10) installs gemmi 0.7.5, ase 3.29.0,
  diffpy.structure **3.2.3** (see section 6 for the version difference that mattered).

## 2. What this PR delivers

- **CIF reader/writer** (`src/discuspdf/cif_io.py`, gemmi-based): parses CIF text,
  expands the asymmetric unit by the space group (symop loop → space-group symbol →
  assumed P1), wraps coordinates to [0, 1), merges periodic images within 1e-4,
  checks total occupancy ≤ 1 + 1e-6, rejects modulated structures; `write_cif` emits
  P1 CIFs with `_discuspdf_*` tags for exact round-trip (model, supercell, site
  labels, atom ids, cell indices), written atomically and verified by re-parsing.
- **ASE adapter** (`src/discuspdf/ase_io.py`): `to_ase` stores the full DISCUS state
  in `atoms.arrays` (`dp_*`) and `atoms.info` (`dp_model`, `dp_supercell`); `from_ase`
  restores it exactly (path A) or, for a foreign `Atoms` object, requires `model=`
  and reports losses (path B). `to_ase` never raises for a loss (design rule 2).
- **DiffPy adapter** (`src/discuspdf/diffpy_io.py`): `to_diffpy` / `from_diffpy`
  preserve fractional coordinates, lattice parameters and distances/angles; the
  model is not stored (pass it to `from_diffpy`); anisotropic atoms are a loss
  (isotropic equivalent kept).
- **Shared helpers** (`src/discuspdf/interop.py`): `InteropError` / `InteropWarning`,
  `LossReport` (raise-or-warn policy), `u_to_b` / `b_to_u` (`B = 8π²U`; gemmi has
  no `u_to_b`), `element_from_type_symbol` (first 100 elements).
- **Atomic writes** (`src/discuspdf/_safe_write.py`): shared by `write_cif` and
  `discus_io.write_cll` / `write_stru` (the only change to `discus_io.py`).
- **Extras** (`pyproject.toml`): `cif` (gemmi), `ase`, `diffpy`, `interop` (all three).
  `import discuspdf` works without them; the adapters raise a helpful `ImportError`.
- **Tests**: `tests/test_cif_io.py` (28), `tests/test_ase_io.py` (16),
  `tests/test_diffpy_io.py` (15) — 59 interop tests, all green with the extras and
  cleanly skipped without them.
- **Fixtures**: `tests/fixtures/interop/` (5 synthetic CIFs + the real COD GaAs CIF
  9008845 + probe evidence logs + `MANIFEST.sha256`).
- **Docs / example / notebook**: `docs/INTEROP.md`, `docs/INTEROP_FACTS.md` (probe
  findings P1–P14), `examples/interop_basics.py`, `notebooks/06_interop.ipynb`
  (executed, outputs committed), `README.md` section, `todo/TASKS.md` A2 ticked.

## 3. Design decisions (overseer rulings, verbatim from `00_START_HERE.md`)

| Decision | Ruling |
|---|---|
| Symmetry source | symop loop → space-group symbol → assumed P1 (recorded in `defaults_applied`) |
| gemmi unit-cell sites | `get_all_unit_cell_sites()` returns UNWRAPPED fractions; the reader wraps to [0,1) and merges periodic images within 1e-4 |
| Quotes | `gemmi.cif.as_string` strips CIF quotes; never ad-hoc `.strip("'")` |
| B ↔ U | `B = 8π²U` (gemmi has no `u_to_b`); CIF reading prefers `_atom_site_B_iso_or_equiv` over U |
| DiffPy lattice | `from_diffpy` takes the lattice from `cell_parms()` (never `.base` — different orientation); fractional `xyz` used directly |
| DiffPy anisotropy | `atom.anisotropy` True → `Uisoequiv` is the correct lattice-weighted equivalent (0.012 cubic, 0.012912 triclinic for diag(0.013, 0.012, 0.011)); kept as `b_iso`, anisotropy is a loss |
| ASE foreign atoms | `info['occupancy']` is never interpreted (the atom→key mapping is not recoverable from `Atoms`); occupancies assumed 1.0, `model=` required |
| `to_diffpy` losses | atom_id not 0..N-1, non-trivial cell_index/supercell, pbc not all True, species charges. **Model is NOT a loss** (A2-3b ruling): DiffPy has no model concept, but `from_diffpy(..., model=)` makes the user state it, so nothing is lost silently |
| `to_ase` | never raises for a loss; warns at most (partial occupancy / mixed site) |
| Loss policy | `allow_loss=False` → `InteropError` listing every loss; `allow_loss=True` → one `InteropWarning` |

## 4. Tolerances (predeclared)

| Quantity | Tolerance |
|---|---|
| lattice parameters | rtol 1e-8 |
| fractional coordinates (CIF round trip) | atol 1e-9 (10-decimal file precision) |
| fractional coordinates (ASE, DiffPy, in memory) | atol 1e-12 |
| occupancy | atol 1e-9 |
| B factor | atol 1e-6 |
| distances | atol 1e-6 Å |
| species, site labels, atom ids, cell indices, row order | exact |

## 5. Tests

| Command (extras venv) | Result |
|---|---|
| `pytest tests/test_cif_io.py tests/test_ase_io.py tests/test_diffpy_io.py` | **59 passed** in 3.34 s |
| `pytest` (full suite) | 28 failed, 366 passed, 11 skipped in 187.90 s |
| base venv (no extras): `pytest tests/test_cif_io.py tests/test_ase_io.py tests/test_diffpy_io.py` | 6 passed, 53 skipped |

Full-suite failures on this machine (28) are in `test_gate2_prereg.py`,
`test_gate3_prereg.py`, `test_native_3d_integration.py`, `test_regenerate_evidence.py`
— the known Windows environment-specific failures (CRLF checkouts / symlinks in the
evidence and probe tests). **The failure-name set is identical to `origin/main`
(28 = 28, zero new failures)** — verified by running the full suite in a
`git worktree` of `origin/main` (`a6e303f`) with the same extras venv.

Mutation checks (A2-2/A2-3, verbatim from the A2-3 review):

| Mutation | Change | Tests failing |
|---|---|---|
| M1 | CIF reader: skip symmetry expansion (return listed sites only) | 9 (incl. `test_gaas_expands_to_8`) |
| M2 | CIF reader: do not wrap to [0,1) | 5 (incl. `test_negative_coords_wrapped`) |
| M3 | `write_cif`: drop `_discuspdf_*` tags | 6 (incl. `test_cif_roundtrip_supercell`) |
| M4 | `to_ase`: drop `dp_*` arrays | 7 (incl. `test_ase_roundtrip_supercell`) |
| M5 | `from_ase` path B: skip the loss checks | 5 (incl. `test_ase_foreign_partial_occupancy_raises`) |
| M6 | `to_diffpy`: drop the supercell loss | 4 (incl. `test_to_diffpy_supercell_is_loss`) |
| N7 (A2-3b) | `to_diffpy`: put the "model" loss back | 7 (incl. `test_to_diffpy_plain_structure_no_loss`) |

## 6. CI and the diffpy version difference

The workflow (`.github/workflows/python-tests.yml`) installs `.[dev,notebook,interop]`
and runs `python -m pytest -q` (no separate notebook-execution step; the notebook is
executed locally and its outputs are committed, matching the other notebooks).

CI (Python 3.10) resolves `diffpy.structure` to **3.2.3** (its `requires_python` is
`<3.13`), while the local venv (Python 3.14) gets **3.5.0**. `Lattice.cell_parms`
does not exist in 3.2.3, so the first CI run failed with
`AttributeError: 'Lattice' object has no attribute 'cell_parms'` in `from_diffpy`
and `test_cif_diffpy_cross_check`. Fixed in A2-4b/A2-4c by reading the individual
attributes (`lat.a, lat.b, lat.c, lat.alpha, lat.beta, lat.gamma`), present in all
versions (3.2.3 has the deprecated `abcABG()` and the single-letter properties;
3.5.0 has `cell_parms()`). After the fix the CI check-run is **green**
(`unit-tests` completed/success on head `a5941cd`).

## 7. Probe facts (A2-1/A2-1b, `docs/INTEROP_FACTS.md`)

| Fact | Value |
|---|---|
| gemmi `read_small_structure` on the COD GaAs CIF | block `9008845`, space group `F -4 3 m` (No. 216); 2 listed sites (Ga 0,0,0; As 0.25,0.25,0.25); `get_all_unit_cell_sites()` → 8 sites, 4 of them outside [0,1) (gemmi does NOT wrap) |
| gemmi on the 5 synthetic CIFs | P-1: 2→4 sites; mixed F-43m: 3→12; aniso P1: U11/U22/U33 read directly; disorder groups read; total occupancy > 1 NOT rejected by gemmi |
| ASE on the COD GaAs CIF | 8 atoms (asymmetric unit expanded); cell matches; `info` has `spacegroup` + `unit_cell`, no `occupancy` key; no U/B in `arrays` |
| ASE on the mixed fixture | 8 atoms (In dropped from `species`); `info['occupancy']` keys are **listed-site indices**, and the atom-index → key mapping is **not recoverable** from `Atoms` (P4/P10) |
| ASE cell validation | the `Atoms` constructor does NOT reject singular or left-handed cells; `from_ase` therefore checks `det > 0` itself |
| diffpy on the COD GaAs CIF | 8 atoms; space group No. 216, 96 sym matrices; U/B default 0.0 |
| diffpy anisotropic Uisoequiv | with `anisotropy=True`: 0.012 (cubic) / 0.012912 (triclinic) for diag(0.013, 0.012, 0.011) — the correct isotropic equivalent (P14) |
| B ↔ U | `8π²·0.01 = 0.7895683520871487`; gemmi has no `u_to_b` (P12) |
| gemmi `as_string` | strips CIF single-quotes; identity on unquoted values (P13) |
| Conflicts with the fixed design | **None** |

## 8. Limitations / not done

- Anisotropic ADPs are a loss in every adapter (isotropic equivalent kept, if any).
- Disorder groups and magnetic moments are reported as losses on CIF reading.
- Total occupancy > 1 at a position is rejected (`InteropError`).
- Mixed pbc (not all True) is a loss on CIF writing (the reader assumes all periodic).
- ASE `info['occupancy']` is never interpreted (overseer ruling).
- DiffPy: orientation is not preserved; site labels, atom ids, cell indices,
  supercell and pbc are not stored; the model must be passed to `from_diffpy`.
- The element table covers the first 100 elements; unknown type symbols raise
  `InteropError`.
- No notebook-execution step in CI (notebooks are executed locally and their
  outputs committed, as for the existing notebooks).

## 9. Changed files (vs `main`)

```
 .github/workflows/python-tests.yml                |    2 +-
  benchmarks/interop_probes/p10_ase_deep.py         |   75 ++
  ... (14 probe scripts, p1–p14)
  docs/INTEROP.md                                   |  171 ++++
  docs/INTEROP_FACTS.md                             |  152 ++++
  examples/interop_basics.py                        |  100 +++
  notebooks/06_interop.ipynb                        |  333 +++++++
  pyproject.toml                                    |    4 +
  src/discuspdf/__init__.py                         |   23 +
  src/discuspdf/_safe_write.py                      |   49 ++
  src/discuspdf/ase_io.py                           |  246 +++++
  src/discuspdf/cif_io.py                           |  980 ++++++++++++++++++++
  src/discuspdf/diffpy_io.py                        |  155 ++++
  src/discuspdf/discus_io.py                        |   68 +-   (only the _safe_write import + writer refactor)
  src/discuspdf/interop.py                          |  102 +++
  tests/fixtures/interop/*                          |  (6 CIFs + 14 evidence logs + manifest + README)
  tests/test_ase_io.py                              |  393 ++++++++
  tests/test_cif_io.py                              | 1000 +++++++++++++++++++++
  tests/test_diffpy_io.py                           |  430 +++++++++
  README.md / todo/TASKS.md                         |  (interop section + A2 tick)
  52 files changed, 5086 insertions(+), 52 deletions(-)
```

`src/discuspdf/structure.py` and `src/discuspdf/lattice.py` are **unchanged** vs
`main` (verified with `git diff --stat origin/main`).
