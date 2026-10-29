# Package A · PR A-1 Report — Lattice and coordinate core (Feature A5)

## 1. Base / head / PR

- **Base SHA:** `2e0254e` (`origin/main`, matches the handoff's reviewed baseline)
- **Head SHA:** `024773d`
- **PR:** [Draft PR #10](https://github.com/jmsweng/discus-python/pull/10)
- **Branch:** `feat/structure-lattice-core`

## 2. Design decisions

| Decision | Choice | Rationale |
|---|---|---|
| Immutability | `__slots__` class with `object.__setattr__` in `__init__`; all getters return **copies** | A1 (Structure) will depend on this; copies prevent accidental mutation of shared state. |
| Stored arrays | `_A` (direct basis), `_B` (reciprocal basis = `inv(A).T`), `_Ainv` (`inv(A)`) | `cart_to_frac` needs `inv(A)` directly; `hkl_to_q` needs `B`. Storing both avoids per-call `inv()`. |
| Reciprocal precision | `np.linalg.solve(A.T, I)` for B; `np.linalg.inv(A)` for A⁻¹ | ~1e-16 accuracy; round-trips pass at rtol 1e-12. |
| `parameters` recovery | `arccos` of row-dot-products with `clip(cos, -1, 1)` | Absorbs round-off; angles recovered to ~1e-10 relative. |
| `minimum_image` tie rule | Within `atol = 1e-12·R + 1e-15` Å: smallest `|n|₁`, then lexicographically smallest `n` | Deterministic; documented in docstring and `docs/LATTICE.md`. |
| `__eq__` | `np.array_equal` on stored `_A`; returns `False` for non-`Lattice` (not `NotImplemented`) | Simpler; `isclose` is the tolerance-based comparison. |
| `from_matrix` validation | `det ≤ 0` → `ValueError`; `cond(A) > 1e12` → `ValueError` | Matches the handoff's rejection rules; no left-handed opt-in. |

No decisions require owner/overseer ruling.

## 3. Changed files

```
 README.md                  |   4 +
 docs/LATTICE.md            |  67 +
 examples/lattice_basics.py |  29 +
 src/discuspdf/__init__.py  |   4 +
 src/discuspdf/lattice.py   | 321 +
 tests/test_lattice.py      | 355 +
 todo/TASKS.md              |   3 +-
 7 files changed
```

## 4. Test commands and results

| Command | Baseline (before) | After |
|---|---|---|
| `python -m pytest -q` | 28 failed, 143 passed, 2 skipped | 28 failed, 156 passed, 2 skipped |
| `python -m pytest -q tests/test_lattice.py -v` | n/a | 13 passed |

The 28 failures are **pre-existing and environment-specific** (Windows CRLF vs LF + Windows symlink semantics in `test_native_3d_integration.py`, `test_gate2_prereg.py`, `test_gate3_prereg.py`, `test_regenerate_evidence.py`). The failure set is **identical** before and after; only the pass count changed (+13 new lattice tests).

## 5. Numerical table

| Group | Worst error | Tolerance |
|---|---|---|
| 2π identity (`hkl_to_q(hkl) · frac_to_cart(f) == 2π·hkl·f`) | 1.7e-15 (relative) | 1e-12 |
| Round trips (frac↔cart, hkl↔Q, scalar + (N,3) + (2,N,3)) | 5.3e-16 (relative) | 1e-12 |
| Parity with `_validated_direct_basis` (23 cells) | 0.0 (exact) | 1e-14 |
| Minimum image vs brute force (γ=20°, 10⁴ displacements, ±4 window) | 0.0 (Å) | — |
| Minimum image vs brute force (α=β=γ=30°, 500 displacements, ±4 window) | 0.0 (Å) | — |

## 6. Environment

- **OS:** Windows 10 (native, not WSL2)
- **Python:** 3.14.7
- **NumPy:** 2.5.4
- **pytest:** 9.2.0

Note: the handoff specifies WSL2 Ubuntu; this run used native Windows Python 3.14. The 28 pre-existing failures are attributable to the Windows environment (CRLF line-ending conversion via `core.autocrlf=true`, and Windows symlink semantics). CI (Python 3.10, Linux) is expected to pass.

## 7. Caveats and items left undone

## Review fix-ups (overseer review of draft PR #10)

**New head SHA:** `af19c5a` · PR remains **draft**.

### What changed per fix

- **Fix 1 — Immutability enforced:** `__init__` now sets `flags.writeable = False` on `_A`, `_B`, `_Ainv`; `__setattr__`/`__delattr__` overridden to raise `AttributeError` (only `object.__setattr__` used inside `__init__`). `test_immutability` extended: in-place write to `lat._A` raises `ValueError`; `lat._A = ...`, `del lat._A`, `lat.new_attr = 1` raise `AttributeError`; `hash(lat)` unchanged.
- **Fix 2 — Bounded search window + chunking:** `minimum_image` now searches around `n0 = rint(df)` with offsets `k_i ∈ [−ceil(R_max|b_i|), +ceil(R_max|b_i|)]` (the `|df|` term dropped), and processes displacements in chunks (≤ 65536/K rows) so peak memory is independent of `m`.
- **Fix 3 — Translation-invariant tie rule:** ties (within `atol = 1e-12·R + 1e-15` Å) now resolved by the **lexicographically smallest reduced fractional vector `df − n`** (components equal if they differ by < 1e-12). Cubic half-cell ties map to −½ ([−½, ½)). Docstring, `docs/LATTICE.md`, and this report updated. New tests: cubic tie cases and a translation-invariance property test (γ=20° and 30°/30°/30° cells, random `df` + integer `n` with |n_i| ≤ 50, including constructed exact ties).
- **Fix 4 — Full brute-force coverage:** the brute force is now vectorised (729 shifts × all rows) and **all** 10⁴ displacements (γ=20°) and all 500 (30° cell) are checked.

### Test counts before/after

| | Before fix-ups | After fix-ups |
|---|---|---|
| `tests/test_lattice.py` | 13 passed | 15 passed |
| Full suite (Windows) | 28 failed, 156 passed, 2 skipped | 28 failed, 158 passed, 2 skipped |

The 28 failures remain the same pre-existing environment-specific set (no regression).

### Numerical results

- **Worst brute-force error** (`d_test − d_ref`, all 10⁴ γ=20° + all 500 30° displacements, ±4 window): **0.0 Å** (the bounded search matches the independent brute force exactly).
- **Large-|df| performance:** 10⁵ displacements with |df| up to 1000 → **6.4 s**, **peak 18.0 MB** (tracemalloc). Previously the window scaled as `m·(2|df|)³` (≈143 MB for 1 000 displacements at |df|≈30); it is now independent of `|df|` and of `m` (chunked).
- **Translation invariance:** `minimum_image(df + n)[0] == minimum_image(df)[0]` to 1e-12 for random `df` and integer `n` (|n_i| ≤ 50), including constructed exact half-cell ties.

### Notes (not required)

- `d_spacing` treats `|hkl| < 1e-8` as zero via `np.allclose` — now mentioned in the docstring.
- `__eq__` returning `False` for non-`Lattice` left as-is (acceptable).
- `todo/TASKS.md` A5 remains ticked (fixes 1–4 are in).

## Caveats (original)

- **Not run on WSL2:** the handoff specifies WSL2 Ubuntu; this agent ran on native Windows. The 28 pre-existing failures are environment-specific and do not reflect code regressions.
- **`from_parameters` does not reconstruct experimental orientation:** documented in the docstring and `docs/LATTICE.md`. A rotated cell must be supplied via `from_matrix`.
- **No new dependencies:** core stays NumPy only.
- **`todo/TASKS.md` A5 ticked** with a link to the PR and `tests/test_lattice.py`, plus the one-line open item for `DiscusVolume` consolidation.
- **No changes to** `discus_output.py`, `reciprocal.py`, `volume.py`, `pdf.py`, `diffuse3d.py`, `benchmarks/`, `spike_discus_3d/`, `examples/run_sample/`, or evidence files.

## Final fix-ups (overseer review, second pass)

**New head SHA:** `5ea0b78` · PR remains **draft**.

### Changes

- **Fix 1 — `minimum_image` returns integer shifts:** `n_out` is now allocated as `np.empty(flat.shape, dtype=np.int64)` (previously `np.empty_like(flat)` → float64, so `n_shift` came back as `[1. 0. 0.]`). New test `test_minimum_image_shift_dtype_and_identity` asserts `np.issubdtype(n_shift.dtype, np.integer)` for `(3,)`, `(N, 3)` and `(2, N, 3)` inputs (cubic, triclinic, γ=20° cells) and that `df − n_shift` equals the returned reduced vector **exactly** (`np.array_equal`).
- **Fix 2 — example notebook:** added `notebooks/03_lattice_basics.ipynb` (written and executed by the overseer against `071661b`), re-executed with `jupyter nbconvert --to notebook --execute --inplace` after Fix 1 — runs cleanly from a fresh kernel, no errors, and the shift column now prints `[1 0 0]` (integer) instead of `[1. 0. 0.]`. No content edited beyond re-execution output changes. Linked from `docs/LATTICE.md`.

### Test counts

| | Before | After |
|---|---|---|
| `tests/test_lattice.py` | 15 passed | **16 passed** |
| Full suite (Windows) | 28 failed, 158 passed, 2 skipped | **28 failed, 159 passed, 2 skipped** |

The 28 failures are the same pre-existing environment-specific set (Windows CRLF/symlink); no regression.

### Untracked `benchmarks/` files (left unstaged, per instructions)

```
benchmarks/discus_6_23_01/correlated_run/inspect_h5.py
benchmarks/discus_6_23_01/correlated_run/run_test7.sh
benchmarks/discus_6_23_01/correlated_run/run_test8.sh
benchmarks/discus_6_23_01/correlated_run/run_test9.sh
benchmarks/discus_6_23_01/correlated_run/test.cll
benchmarks/discus_6_23_01/correlated_run/test.mac
benchmarks/discus_6_23_01/correlated_run/test.stru
benchmarks/discus_6_23_01/correlated_run/test2.mac
benchmarks/discus_6_23_01/correlated_run/test3.mac
benchmarks/discus_6_23_01/correlated_run/test4.mac
benchmarks/discus_6_23_01/correlated_run/test5.mac
benchmarks/discus_6_23_01/correlated_run/test6.mac
benchmarks/discus_6_23_01/correlated_run/test7.mac
benchmarks/discus_6_23_01/correlated_run/test7_wsL.log
benchmarks/discus_6_23_01/correlated_run/test8.log
benchmarks/discus_6_23_01/correlated_run/test8.mac
benchmarks/discus_6_23_01/correlated_run/test9.log
benchmarks/discus_6_23_01/correlated_run/test9.mac
benchmarks/discus_6_23_01/correlated_run/test_proper.cll
benchmarks/discus_6_23_01/correlated_run/test_stru.cll
benchmarks/discus_6_23_01/correlated_run/test_wsL.cll
benchmarks/discus_6_23_01/run_discus.sh
```
