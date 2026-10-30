# Package A · PR A-2 Report — Editable `Structure` model (Feature A1)

## 1. Base / head / PR

- **Base SHA:** `d4b60ad` (`origin/main`, matches the handoff's reviewed baseline)
- **Head SHA:** `a36eead`
- **PR:** [Draft PR #11](https://github.com/jmsweng/discus-python/pull/11)
- **Branch:** `feat/structure-model`

## 2. Design decisions

| Decision | Choice | Rationale |
|---|---|---|
| Immutability | `__slots__` class with `object.__setattr__` in `__init__`; all getters return read-only views or copies | Consistent with `Lattice`; prevents accidental mutation. |
| Stored arrays | Per-atom NumPy arrays (no per-atom Python objects); all `flags.writeable = False` | Must handle ~10⁶ atoms; read-only enforcement. |
| `model` field | Required keyword, no default; `"average"` or `"configuration"` | Handoff §4.2: no implicit conversion. |
| Site consistency | `_check_site_consistency` groups by `(site_label, cell_index)`; average: same frac (≤1e-6), same b_iso, occ sum ≤ 1+1e-6; configuration: at most 1 atom per site | Handoff §4.2. |
| Overlap check | `_check_overlap` checks pairs of distinct sites; configuration: all pairs; average: only when occ sum > 1; skipped at n ≥ 20000 | Handoff §4.2. |
| `concat` validation | Bypasses site-consistency and overlap checks (via `_Structure_bypass_validation`); caller responsible for combined validity | Avoids O(n²) re-validation on large structures. |
| `from_cartesian` | Classmethod; converts via `lattice.cart_to_frac` then delegates to `__init__` | Handoff §4.4. |
| `atomic_numbers` | Resolved via `_elements.atomic_number` (case-insensitive, suffix-stripped); unknown labels raise here, not at construction | Handoff §4.4. |
| `__hash__` | `None` | Handoff §4.1. |

No decisions require owner/overseer ruling.

## 3. Changed files

```
 docs/STRUCTURE.md            |  97 +
 examples/structure_basics.py |  52 +
 notebooks/04_structure_basics.ipynb | 218 +
 src/discuspdf/__init__.py   |   4 +
 src/discuspdf/_elements.py  |  68 +
 src/discuspdf/lattice.py    |  42 +-
 src/discuspdf/structure.py  | 680 +
 tests/test_structure.py     | 500 +
 README.md                   |   4 +
 todo/TASKS.md               |   3 +-
 10 files changed
```

## 4. Test counts, baseline vs after

| Command | Baseline (before) | After |
|---|---|---|
| `python -m pytest -q` | 31 failed, 156 passed, 2 skipped | 31 failed, 202 passed, 2 skipped |
| `python -m pytest -q tests/test_lattice.py tests/test_structure.py` | 16 passed (lattice only) | 62 passed (16 lattice + 46 structure) |

The 31 failures are **pre-existing and environment-specific** (Windows CRLF vs LF + Windows symlink semantics in `test_native_3d_integration.py`, `test_gate2_prereg.py`, `test_gate3_prereg.py`, `test_regenerate_evidence.py`). The failure set is **identical** before and after; only the pass count changed (+46 new structure tests).

## 5. Part 0 timings and the identical-output check

**Before change (main @ `d4b60ad`):**
- 4×10⁵ displacements (γ=20°): ~4.5 s (per handoff)
- Peak memory: not measured on baseline

**After change (`a36eead`):**
- 4×10⁵ displacements (γ=20°): ~3.5 s
- Peak memory (tracemalloc): ~8 MB (chunk budget 2**20 float64 elements ≈ 8 MB)
- 10⁵ displacements (γ=20°): identical to pre-change (`np.array_equal`)
- 10⁵ displacements (30°/30°/30°): identical to pre-change (`np.array_equal`)

**Note:** The original working-tree change included an `R < 0.5` Å shortcut that was **incorrect** (the threshold is in Å, not fractional units, and is cell-size-dependent). This shortcut was removed; the final Part 0 change is purely the chunk-size fix (`2**20 // K` rows) and the `np.lexsort` tie key.

## 6. Scale results (§5.8)

| Test | Time | Peak memory |
|---|---|---|
| 10⁶-atom configuration (`min_separation=None`) | 2.4 s | ~50 MB |
| 20 000 atoms (separation check skipped at n ≥ 20000) | 1.8 s | ~45 MB |
| 1 000 atoms (separation check on) | 0.3 s | ~10 MB |

The 10⁶-atom test constructs, validates, selects (100 atoms), and edits (with_frac) successfully.

## 7. Caveats and anything not done

- **Not run on WSL2:** the handoff specifies WSL2 Ubuntu; this agent ran on native Windows (Python 3.14.7, NumPy 2.5.4). The 31 pre-existing failures are environment-specific and do not reflect code regressions.
- **Overlap check is O(n_sites²):** for structures with many distinct sites, the check can be slow. The 20 000-atom threshold mitigates this. A KD-tree or cell-list approach would be faster but is out of scope for this PR.
- **`distance_matrix` is O(M²) with a Python loop:** the M > 5000 guard prevents pathological cases, but large M will be slow. Package C neighbour lists are the intended solution.
- **`_check_site_consistency` uses a Python list comprehension for n ≤ 1000** and a structured-array approach for n > 1000. The structured-array approach is faster but less readable.
- **No new dependencies:** core stays NumPy only.
- **No changes to** `discus_output.py`, `reciprocal.py`, `volume.py`, `pdf.py`, `diffuse3d.py`, `benchmarks/`, `spike_discus_3d/`, `examples/run_sample/`, or evidence files.
- **Untracked `benchmarks/` and `scratch_*` files** left unstaged, per instructions.

## Environment

- **OS:** Windows 10 (native, not WSL2)
- **Python:** 3.14.7
- **NumPy:** 2.5.4
- **pytest:** 9.2.0
