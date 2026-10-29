> **Historical document (kept by owner decision, 2026-10-28).** Written during earlier PR work. It is **not** part of PR #8's validated results and may contain claims since superseded or refuted. The current record is `docs/SCIENTIFIC_VALIDATION_REPORT.md` and `benchmarks/discus_6_23_01/REVIEW_2026-10-28.md`.

# Off-Commensurate Structure-Factor Benchmark — DISCUS 6.23.01 (v4.2)

**Date:** 2026-10-28  
**Environment:** WSL2 Ubuntu 22.04, DISCUS_SUITE 6.23.01  
**Structure:** C1_frozen.stru (4096 atoms: 1723 GA + 2373 AS)  
**Structure SHA-256:** `ceefca57e26ab40e03d2adfd443cff559bd7fc3c837402f418f6e8def35406b9`

## 1. Purpose

This benchmark validates the independent Python structure-factor calculation against a **real DISCUS calculation on an off-commensurate grid**, addressing the most important finding from the v4 independent review:

> "No real DISCUS off-commensurate experiment. The new `off_grid_verification` JSON values are generated solely by the Python model."

## 2. Experiment

### 2.1 Grid

- **Dimensions:** 33×33×33
- **Corner:** (−0.16, −0.16, −0.16)
- **Step:** 0.01 rlu
- **Range:** h,k,l from −0.16 to +0.16
- **Key property:** This is **deliberately NOT a multiple-of-1/8 reciprocal-space grid**, so the interference factor is nonzero at off-grid points.

### 2.2 DISCUS Macro

```text
discus
read
  stru C1_frozen.stru
fourier
  xray
  wvle moa1
  ll -0.16,-0.16,-0.16
  lr  0.16,-0.16,-0.16
  ul -0.16, 0.16,-0.16
  tl -0.16,-0.16, 0.16
  na 33
  no 33
  nt 33
  abs h
  ord k
  top l
  temp use
  disp off
  lots off
  set technique:turbo
  set aver,0.0
  show
  run
exit
output
  value int
  form nexus
  outf C1_offcomm_total.h5
  run
exit
fourier
  set aver,100.0
  show
  run
exit
output
  value int
  form nexus
  outf C1_offcomm_diffuse.h5
  run
exit
exit
exit
```

### 2.3 Verification

- **Exit code:** 0
- **Output files:** `C1_offcomm_total.h5` (312 KB), `C1_offcomm_diffuse.h5` (312 KB)
- **Engine log:** Confirms 4096 atoms, 33×33×33 grid, corner (−0.16, −0.16, −0.16), step 0.01
- **Metadata:** `data_type_with_bragg = "bragg"` (total) and `"bragg_subtracted"` (diffuse)
- **All values finite:** Yes
- **No engine errors:** Yes

## 3. Results

### 3.1 Total Intensity

All 8 test points match within 0.5% (max rel L2 = 0.046%):

| hkl | DISCUS total | Independent total | Rel L2 | Match (1%) |
|-----|-------------|------------------|--------|------------|
| (0.05, 0, 0) | 9.932×10⁹ | 9.928×10⁹ | 0.046% | **YES** |
| (0.10, 0, 0) | 9.494×10⁸ | 9.491×10⁸ | 0.030% | **YES** |
| (0.10, 0.15, 0.07) | 6.994×10⁶ | 6.991×10⁶ | 0.040% | **YES** |
| (0, 0, 0) | 1.734×10¹⁰ | 1.733×10¹⁰ | 0.040% | **YES** |
| (0.16, 0, 0) | 6.362×10⁸ | 6.359×10⁸ | 0.043% | **YES** |
| (−0.16, 0, 0) | 6.362×10⁸ | 6.359×10⁸ | 0.043% | **YES** |
| (0.05, 0.05, 0.05) | 3.227×10⁹ | 3.226×10⁹ | 0.028% | **YES** |
| (0.07, 0.03, 0.01) | 4.346×10⁹ | 4.344×10⁹ | 0.033% | **YES** |

**Conclusion:** The independent total intensity calculation is **fully validated** on an off-commensurate grid.

### 3.2 Diffuse Intensity

Only 2 of 8 points match within 1%:

| hkl | DISCUS diffuse | Independent diffuse | Rel L2 | Match (1%) |
|-----|---------------|-------------------|--------|------------|
| (0.05, 0, 0) | 14.08 | 32.36 | 129.9% | **NO** |
| (0.10, 0, 0) | 134.44 | 148.01 | 9.99% | **NO** |
| (0.10, 0.15, 0.07) | 287.11 | 285.52 | 0.55% | **YES** |
| (0, 0, 0) | 0.00 | 0.00 | 0% | **YES** |
| (0.16, 0, 0) | 510.13 | 568.75 | 11.5% | **NO** |
| (−0.16, 0, 0) | 510.13 | 568.75 | 11.5% | **NO** |
| (0.05, 0.05, 0.05) | 310.76 | 342.71 | 10.3% | **NO** |
| (0.07, 0.03, 0.01) | 93.72 | 125.79 | 25.7% | **NO** |

**Conclusion:** The site-resolved average structure model does **NOT** reproduce DISCUS's diffuse intensity at off-commensurate points.

### 3.3 Interpretation

The empirical site occupancies in F_avg were derived **from the same frozen structure** whose F is being calculated. At integer hkl, the species-dependent phase is identical across cell translations, so F and this fitted F_avg are identical *by construction*. This is a useful software-consistency test but does **not independently prove** which averaging/Bragg-subtraction algorithm DISCUS uses away from those Bragg nodes.

The off-commensurate results confirm this: DISCUS's averaging procedure is **not** a simple site-resolved occupancy model. The exact algorithm DISCUS uses for `set aver` at off-commensurate hkl remains **undocumented and unverified**.

## 4. v4.2 Improvements

This v4.2 benchmark addresses all remaining issues from the v4 review:

1. **Full metadata grid construction:** `build_hkl_grid_from_metadata()` uses the complete 3-vector corner and 3×3 increment vector, not a shared `hmin` and `step`.

2. **Direct-site-average oracle:** `calculate_average_structure_factor_direct_site()` independently computes the average amplitude as a direct weighted summation over all 4,096 actual site positions, without calling the factored `interference_factor()` formula.

3. **Explicit off-grid assertions:** The 4 off-grid reference values from the v4 review are now tight, explicit assertions in the regression test suite.

4. **Synthetic supercell test:** A 2×2×2 supercell with analytically known amplitudes validates the core structure-factor calculation.

5. **Origin-shift phase check:** Verifies that shifting the origin by (0.5, 0, 0) introduces the correct phase factor.

6. **Asymmetric grid fixture:** Tests that `build_hkl_grid_from_metadata()` correctly handles asymmetric corner and increment vectors.

7. **Failure for unsupported structures:** The reader now raises `ValueError` for unsupported species and non-4096-atom structures, rather than silently skipping atoms.

8. **Documented restrictions:** All restrictions (8-cell range, 4096 atoms, 8 basis sites, cubic a, uniform Biso) are explicitly documented in the code.

## 5. Regression Tests

All 7 tests pass:

1. ✅ Synthetic supercell (analytically known amplitudes)
2. ✅ Origin-shift phase check
3. ✅ Asymmetric hkl-corner/step fixture
4. ✅ Off-grid reference values (4 explicit assertions)
5. ✅ Direct-site-average oracle (agrees with factored formula at integer hkl)
6. ✅ Unsupported species raises ValueError
7. ✅ Non-4096-atom structure raises ValueError

## 6. Files

- `independent_structure_factor_v4_2.py` — Corrected Python reference implementation (v4.2)
- `test_structure_factor_v4_2.py` — Regression tests (7 tests, all pass)
- `C1_offcomm_total.h5` — DISCUS off-commensurate total intensity (33×33×33)
- `C1_offcomm_diffuse.h5` — DISCUS off-commensurate diffuse intensity (33×33×33)
- `offcommensurate_macro.cl` — DISCUS macro for the off-commensurate calculation
- `offcommensurate.log` — Full engine log
- `results_offcomm.json` — Off-commensurate comparison metrics

## 7. Checksums

| File | SHA-256 |
|------|---------|
| C1_frozen.stru | `ceefca57e26ab40e03d2adfd443cff559bd7fc3c837402f418f6e8def35406b9` |
| C1_offcomm_total.h5 | `b7a36c50f64c0f1f92feea8e0bf1a49c8d7047defb638f6f2ce77c0f3c6fa4fe` |
| C1_offcomm_diffuse.h5 | `d2c2e913ad6b7c123af3d0a8a3b68895bc7112fc180bef8c414bad4b11c25ce4` |
| offcommensurate_macro.cl | `da7cde8605f233323f5ceb89fcbee1fb9cc0817abbee391bb409a82396c6ea83` |
| offcommensurate.log | `7647218b7f0d2156f9d4acdc96f2abd22c787132e7ee973758905ca7f342d3be` |
| results_offcomm.json | `4ea83e947e908f9b5de5c84c7277755e32d7c3a5137c641d3bf3e6c1d8e3c6c9` |

## 8. Conclusion

The v4.2 benchmark **validates the independent total intensity calculation** against a real DISCUS off-commensurate experiment. The total intensity matches within 0.5% at all 8 test points.

The **diffuse intensity does NOT match** at 6 of 8 off-commensurate points, confirming that DISCUS's `set aver` algorithm is not a simple site-resolved occupancy model. The exact averaging procedure DISCUS uses at off-commensurate hkl remains undocumented.

**The full 3D-ΔPDF scientific benchmark is NOT complete.** A separate correlated-disorder model with known off-origin peak amplitudes is still required.
