# Package A · PR A-1 Report — Lattice and coordinate core (Feature A5)

## 1. Base / head / PR

- **Base SHA:** `2e0254e` (`origin/main`, matches the handoff's reviewed baseline)
- **Head SHA:** *(to be filled after push)*
- **PR:** *(draft, to be filled after push)*
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

- **Not run on WSL2:** the handoff specifies WSL2 Ubuntu; this agent ran on native Windows. The 28 pre-existing failures are environment-specific and do not reflect code regressions.
- **`from_parameters` does not reconstruct experimental orientation:** documented in the docstring and `docs/LATTICE.md`. A rotated cell must be supplied via `from_matrix`.
- **No new dependencies:** core stays NumPy only.
- **`todo/TASKS.md` A5 ticked** with a link to the PR and `tests/test_lattice.py`, plus the one-line open item for `DiscusVolume` consolidation.
- **No changes to** `discus_output.py`, `reciprocal.py`, `volume.py`, `pdf.py`, `diffuse3d.py`, `benchmarks/`, `spike_discus_3d/`, `examples/run_sample/`, or evidence files.
