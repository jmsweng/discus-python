# DiscusVolume / Lattice facts and parity baseline

**Date:** 2026-11-01
**Branch:** `refactor/discusvolume-lattice` (from `origin/main` at `d9ba5cf`)

## Facts

| Fact | Evidence file | Verbatim excerpt |
|---|---|---|
| `DiscusVolume` has its own `_validated_direct_basis` function | `src/discuspdf/discus_output.py:35` | `def _validated_direct_basis(lengths, angles) -> np.ndarray:` |
| `DiscusVolume.direct_basis` calls `_validated_direct_basis` | `src/discuspdf/discus_output.py:108-110` | `def direct_basis(self) -> np.ndarray: return _validated_direct_basis(self.cell_lengths, self.cell_angles)` |
| `DiscusVolume.reciprocal_basis` computes `inv(direct_basis).T` | `src/discuspdf/discus_output.py:112-114` | `def reciprocal_basis(self, two_pi: bool = True) -> np.ndarray: b = np.linalg.inv(self.direct_basis()).T` |
| `DiscusVolume.cell_volume` computes `abs(det(direct_basis))` | `src/discuspdf/discus_output.py:118-119` | `def cell_volume(self) -> float: return float(abs(np.linalg.det(self.direct_basis())))` |
| `Lattice` has a `matrix` property (direct basis) | `src/discuspdf/lattice.py:131-134` | `@property def matrix(self) -> np.ndarray: """Direct basis ``A`` (rows a, b, c in Å). Returns a copy; the stored array is read-only.""" return self._A.copy()` |
| `Lattice` has a `reciprocal_matrix` property (no 2π) | `src/discuspdf/lattice.py:136-139` | `@property def reciprocal_matrix(self) -> np.ndarray: """Reciprocal basis ``B = inv(A).T`` (rows a*, b*, c* in Å⁻¹, no 2π). Returns a copy.""" return self._B.copy()` |
| `Lattice` has a `volume` property | `src/discuspdf/lattice.py:158-159` | `@property def volume(self) -> float:` |
| `Lattice` uses the same rejection rules as `_validated_direct_basis` | `src/discuspdf/lattice.py:16,37` | `plane, right-handed) and uses the same rejection rules as ``discus_output._validated_direct_basis``.` |
| `reciprocal_grid.py` already uses `Lattice` | `src/discuspdf/reciprocal_grid.py:8` | `from .lattice import Lattice` |
| Direct basis comparison: max abs = 0.00e+00 for all 200 random + 3 special cells | `benchmarks/discusvolume_lattice/compare_bases.log` | (see log; all direct basis differences are exactly 0) |
| Cell volume comparison: max abs = 0.00e+00 for all cells | `benchmarks/discusvolume_lattice/compare_bases.log` | (see log; all volume differences are exactly 0) |
| Reciprocal basis (2π) comparison: max abs varies (0.36 to 4.14), max rel = 0.841 | `benchmarks/discusvolume_lattice/compare_bases.log` | (see log; the 2π convention difference is expected) |
| Invalid inputs: both reject the same cases with similar error messages | `benchmarks/discusvolume_lattice/compare_bases.log` | `angles_sum_360: discus_output: ValueError: unit_cell_angles must be finite degrees strictly between 0 and 180; lattice: ValueError: unit cell angles must be finite degrees strictly between 0 and 180` |
| Parity baseline: 79 passed, 1 skipped (test_viz.py excluded due to missing PIL) | `benchmarks/discusvolume_lattice/baseline/pytest.log` | `79 passed, 1 skipped in 104.45s (0:01:44)` |
| Gate 2/3 analysis scripts cannot run offline (HDF5 files not committed) | `benchmarks/discus_6_23_01/gate2_prereg/analyze.py:6-10` | `RUN_DIR must contain the generated inputs (copied unmodified from inputs/), the four engine HDF5 files gate2_{block,random}_{total,diffuse}.h5` |

## Conflicts

None found. The direct basis and cell volume computations are bit-for-bit identical between `DiscusVolume` and `Lattice`. The reciprocal basis differs by the 2π convention (expected).

## Open questions for the overseer

1. Should `DiscusVolume.reciprocal_basis(two_pi=True)` delegate to `Lattice.reciprocal_matrix * 2π`, or should `Lattice` gain a `two_pi` parameter?
2. The `np.linalg.inv` calls in `discus_output.py:114` and `lattice.py:90` are mathematically equivalent but may differ in floating-point rounding for ill-conditioned matrices. Should we use a shared helper?
3. `DiscusVolume.is_orthogonal` (line 94) calls `_validated_direct_basis` and checks the off-diagonal elements. Should this delegate to a `Lattice` method?
4. The `cell_volume` property in `discus_output.py:118-119` computes `abs(det(A))`. `Lattice.volume` (line 158) likely does the same. Should we use a shared implementation?
5. Are there any performance considerations for delegating to `Lattice` (e.g., copying arrays vs. views)?
