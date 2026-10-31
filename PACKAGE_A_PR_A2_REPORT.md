# Package A · PR A-2 Report — Editable `Structure` model (Feature A1)

## 1. Base / head / PR / environment

- **Base:** `d4b60ad` (`main`, PR #10 merged) · **Head at measurement:** `6872dff`
- **PR:** [Draft PR #11](https://github.com/jmsweng/discus-python/pull/11) · **Branch:** `feat/structure-model`
- **Measured on:** Windows 10, Python 3.14.7, NumPy 2.5.4 (numbers in sections 5–6 come from this machine
  unless marked "cloud").

## 2. What this PR delivers

- `discuspdf.Structure`: an immutable structure value (per-atom NumPy arrays + `Lattice`) with an
  explicit `model`: `"average"` (CIF-like site-occupancy model; mixed sites share
  `(site_label, cell_index)`, occupancy sum ≤ 1) or `"configuration"` (realised arrangement; all
  occupancies 1, at most one atom per site). No implicit conversion between them.
- Validation on construction, on every edit and on `concat`: shapes, finiteness, species and
  site-label strings, `pbc` booleans, occupancies, B factors, site consistency and a 0.1 Å overlap
  check (cell list; never skipped automatically; `min_separation=None` opts out).
- Identity: `ids=` always means `atom_id`; ids change only through `renumbered()`; no method sorts rows.
- Distances (`distance`, `displacement`, vectorised `distance_matrix`) for fully periodic or fully
  non-periodic structures; mixed PBC is stored but distances raise `NotImplementedError`.
- `Lattice.minimum_image` (from PR #10) made faster with identical output (see section 6).
- Docs `docs/STRUCTURE.md`, `examples/structure_basics.py`, `notebooks/04_structure_basics.ipynb`.

## 3. Design decisions

| Decision | Choice |
|---|---|
| Value semantics | Immutable; arrays read-only; edits return new, fully re-validated structures |
| Models | Required `model` keyword; `"average"` vs `"configuration"`; no implicit sampling (A4) |
| Overlap rule | Configuration: all pairs in different sites; average: only if the two sites' occupancy sums exceed 1 |
| Overlap method | Periodic cell list (bins ≥ cutoff thick per lattice direction; ≤ 4 bins per atom); nearest rounded image when every axis has ≥ 3 bins, else `Lattice.minimum_image` |
| `concat` | Validating; requires matching lattice, model, `supercell`, `pbc`, and `b_iso` on both or neither |
| `_subset` (`select`/`remove`) | Skips re-validation: a subset of a valid structure is valid |
| `isclose` | Exact on species, labels, ids, cell_index, supercell, pbc, model; tolerances on lattice, frac, occupancy, b_iso |
| `distance_matrix` | Vectorised, chunked; capped at 5 000 atoms (neighbour lists are Package C) |

## 4. Review history

Overseer review of head `04d1034` found five blockers and four should-fix items; all are resolved:

| Item | Problem | Fixed in |
|---|---|---|
| B1 | Batched `minimum_image` returned wrong (longer) images in non-cubic cells; one-sided test missed it | `ea2088f`, `01584c6` (per-row tie resolution; two-sided and batch-vs-single tests) |
| B2 | Edits skipped the overlap check | `d85ccec` |
| B3 | `concat` bypassed validation | `d85ccec` |
| B4 | `ids=` treated as row numbers | `01584c6` |
| B5 | Overlap check ~O(n^2.5) and silently skipped ≥ 20 000 atoms; per-pair `distance_matrix` | 2b-1 commit (cell list; vectorised `distance_matrix`) |
| S1–S3 | `isclose` too loose; `site_label`/`pbc` unvalidated; unreadable messages | `01584c6` |
| Found in review | `_check_site_consistency` checked the wrong rows when a mixed site did not sort first | `d85ccec` (+ regression test) |

## 5. Tests

| Command | Result |
|---|---|
| `pytest tests/test_lattice.py tests/test_structure.py` | 87 passed in 103.93s (0:01:43) (failures: none) |
| `pytest` (full suite) | 28 failed, 230 passed, 2 skipped in 186.53s (0:03:06) |

Full-suite failures on this machine are in: `tests/test_gate2_prereg.py`, `tests/test_gate3_prereg.py`, `tests/test_native_3d_integration.py`, `tests/test_regenerate_evidence.py`. On native Windows these are the
known environment-specific failures (CRLF checkouts and symlinks in the evidence/probe tests); they
occur identically on `main` and do not involve the files changed here. Standard CI (Linux) is green.

Independent checks run by the overseer (cloud):
- `minimum_image` output identical to `main` (`d4b60ad`) on 6 × 60 000 displacements incl. constructed
  ties; 0 mismatches against an independent brute-force + tie-rule oracle (5 cells × 3 seeds).
- Cell-list overlap check identical to a brute-force all-pairs, all-images check in 960 random cases
  (594 with an overlap) across 6 cells incl. strongly skewed and smaller-than-cutoff cells.

## 6. Performance (this machine)

| Operation | Size | Time |
|---|---|---|
| Construct + validate GaAs 4×4×4 configuration | 512 atoms | 0.02 s |
| Construct + validate GaAs 8×8×8 | 4 096 atoms | 0.03 s |
| Construct + validate GaAs 16×16×16 | 32 768 atoms | 0.26 s |
| Construct + validate GaAs 20×20×20 | 64 000 atoms | 0.58 s |
| Construct + validate random configuration | 100 000 atoms | 0.77 s |
| `distance_matrix` | 2 000 atoms | 3.85 s |
| `Lattice.minimum_image`, γ = 20° cell | 400 000 displacements | 2.94 s |

For comparison (cloud runner): before the fixes, the 512-atom case took 16 s and the 4 096-atom case
was estimated at ~40 min; `minimum_image` on 4×10⁵ displacements went from 18.3 s to 2.8 s.

## 7. Caveats / not done

- Mixed-PBC (slab) distances and overlap checks are not implemented (`NotImplementedError`).
- Magnetic moments (Package D11) and anisotropic ADPs are not represented.
- No CIF/ASE/DiffPy (A2), `.cll`/`.stru` I/O (A3) or supercell/disorder editing (A4) yet.
- `distance_matrix` is O(M²) by design and capped at 5 000 atoms.

## 8. Changed files (vs `main`)

```
PACKAGE_A_PR_A2_REPORT.md           |  107 ++++
 README.md                           |    4 +
 docs/STRUCTURE.md                   |  124 +++++
 examples/structure_basics.py        |   64 +++
 notebooks/04_structure_basics.ipynb |  316 +++++++++++
 src/discuspdf/__init__.py           |    4 +
 src/discuspdf/_elements.py          |   68 +++
 src/discuspdf/lattice.py            |   51 +-
 src/discuspdf/structure.py          | 1005 +++++++++++++++++++++++++++++++++++
 tests/test_lattice.py               |   97 +++-
 tests/test_structure.py             |  942 ++++++++++++++++++++++++++++++++
 todo/TASKS.md                       |    3 +-
 12 files changed, 2765 insertions(+), 20 deletions(-)
```
