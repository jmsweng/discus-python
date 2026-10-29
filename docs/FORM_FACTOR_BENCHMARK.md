> **Historical document (kept by owner decision, 2026-10-28).** Written during earlier PR work. It is **not** part of PR #8's validated results and may contain claims since superseded or refuted. The current record is `docs/SCIENTIFIC_VALIDATION_REPORT.md` and `benchmarks/discus_6_23_01/REVIEW_2026-10-28.md`.

# X-RAY FORM-FACTOR TABLE BENCHMARK

**Date:** 2026-10-28
**Environment:** WSL2 Ubuntu 22.04, DISCUS_SUITE 6.23.01
**Frozen structure:** `C1_frozen.stru` (SHA-256 `ceefca57e26ab40e03d2adfd443cff559bd7fc3c837402f418f6e8def35406b9`, 4096 atoms: 1723 GA + 2373 AS)
**Run directory:** `~/discus_3d_frozen_corrected_20261028T022435Z/formfactor/` (WSL2)

## 1. Objective

Quantify the isolated effect of the X-ray form-factor table on total and diffuse scattering intensity, using the same frozen C1 atomic structure for both calculations.

## 2. Method

Two scattering calculations were performed on the identical frozen C1 structure:

| Variant | Form-factor table | Macro |
|---|---|---|
| Default | Waasmaier–Kirfel (DISCUS default) | `scatter_default.mac` (`xray`) |
| Interpolated | Interpolated table | `scatter_inter.mac` (`xray table:inter`) |

Both macros use the same Fourier grid (65³, −4..+4 rlu), MOA1 wavelength, `set aver,0.0` (total) and `set aver,100.0` (diffuse), and `read stru C1_frozen.stru`.

## 3. Engine log verification

Both runs report:

```
Lot #    1/   1 : complete crystal
           Atom typ = GA         (#      1723 )
           Atom typ = AS         (#      2373 )
```

Both calculations use the same 4096-atom structure. The only difference is the form-factor table.

## 4. Results

### 4.1 Total intensity

| Metric | Default | Interpolated |
|---|---|---|
| Max | 1.733590×10¹⁰ | 1.732906×10¹⁰ |
| Min | 8.394043×10⁻¹ | 8.405259×10⁻¹ |
| Mean | 1.791812×10⁶ | 1.791589×10⁶ |

| Comparison metric | Value |
|---|---|
| Max absolute difference | 6.847×10⁶ |
| Mean absolute difference | 7.346×10² |
| Relative L2 error (vs default) | 5.071×10⁻⁴ |
| Pearson correlation | 0.9999998902 |
| Max relative difference | 7.385×10⁻³ |
| Mean scale ratio (inter/default) | 1.000420 |
| Bit-identical | **False** |

### 4.2 Diffuse (Bragg-subtracted) intensity

| Metric | Default | Interpolated |
|---|---|---|
| Max | 7545.1438 | 7529.914 |
| Min | 0.0 | 0.0 |
| Mean | 410.1429 | 410.3855 |

| Comparison metric | Value |
|---|---|
| Max absolute difference | 25.579 |
| Mean absolute difference | 1.498 |
| Relative L2 error (vs default) | 4.053×10⁻³ |
| Pearson correlation | 0.9999846420 |
| Max relative difference | 7.385×10⁻³ |
| Mean scale ratio (inter/default) | 1.000418 |
| Bit-identical | **False** |

### 4.3 Cross-check: default frozen vs original C1

| Comparison | Max abs diff | Bit-identical |
|---|---|---|
| Total (default frozen vs original C1) | 0.0 | **True** |
| Diffuse (default frozen vs original C1) | 0.0 | **True** |

The default form-factor frozen run is bit-identical to the original C1, confirming that the frozen-structure workflow is correct and that the differences above are solely due to the form-factor table.

## 5. Interpretation

The two form-factor tables produce **measurably different** scattering intensities:

- **Total intensity:** The interpolated table produces a maximum that is 0.04% lower than the default (1.7329×10¹⁰ vs 1.7336×10¹⁰). The mean scale ratio is 1.000420, indicating the interpolated table is ~0.042% higher on average. The Pearson correlation is 0.99999989, indicating the two tables produce highly similar spatial distributions with small local variations.

- **Diffuse intensity:** The interpolated table produces a maximum that is 0.20% lower than the default (7529.9 vs 7545.1). The mean scale ratio is 1.000418. The Pearson correlation is 0.99998464, again indicating highly similar spatial distributions.

- **Relative L2 error:** The total intensity differs by 5.07×10⁻⁴ (relative L2), and the diffuse intensity differs by 4.05×10⁻³. These are small but measurable differences, confirming that the form-factor table choice has a quantifiable effect on the scattering calculation.

The differences are consistent with the expected physics: the Waasmaier–Kirfel table and the interpolated table use slightly different form-factor values, particularly at high scattering angles (high Q), where the form factor falls off. The 0.04% mean scale difference and 0.2% maximum difference in the diffuse intensity are consistent with this expectation.

## 6. Checksums

| File | SHA-256 |
|---|---|
| `C1_frozen.stru` | `ceefca57e26ab40e03d2adfd443cff559bd7fc3c837402f418f6e8def35406b9` |
| `C1_default_total.h5` | `ae4a060e23cadfc8ff7d695555c42a89a8d89858ae59ec75d9eae184be8d62cc` |
| `C1_default_diffuse.h5` | `8b3540b41a56eb111de18c208eb3fcb4cac204e997d169666a7e9c3408cf1d66` |
| `C1_inter_total.h5` | `b4b15918417b481820def52641c98398aa0a0fcf1ee6e78ad743887ca71ab5c8` |
| `C1_inter_diffuse.h5` | `d5ee5072287de46a555781382f8dff31b0c3a009b95b2ea3e745d646427fe108` |

## 7. Conclusions

1. **The form-factor table choice has a measurable effect** on both total and diffuse scattering intensity. The interpolated table produces intensities that are ~0.04% different on average from the default Waasmaier–Kirfel table.

2. **The spatial distribution of intensity is highly similar** between the two tables (Pearson correlation > 0.9999 for both total and diffuse), indicating that the form-factor table choice affects the overall scale more than the spatial pattern.

3. **The frozen-structure workflow is validated.** The default form-factor frozen run is bit-identical to the original C1, confirming that the `read stru` correction is correct and that the form-factor differences are isolated.

4. **For scientific benchmarking, the form-factor table must be specified.** The 0.04% mean difference is small but measurable, and the choice of table should be documented in any published results.

## 8. Status

- **Form-factor table benchmark:** **COMPLETE**
- **Frozen C1 correction:** **COMPLETE** (validated)
- **Independent 3D-ΔPDF benchmark:** **NOT COMPLETE** (separate milestone)
- **PR #7:** **NOT MERGED** (per instructions)
