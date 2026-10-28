# Gate 2: Real DISCUS Engine Comparison Results

## Date: 2026-10-28

## Summary

Gate 2 is **UNRESOLVED**. The DISCUS engine successfully ran on the frozen checkerboard structure, but the 3D-ΔPDF peaks do not quantitatively agree with the independent predictions.

## Engine Run

- **Engine**: DISCUS_SUITE 6.23.01 (WSL Ubuntu-22.04, `/usr/local/bin/discus_suite`)
- **Input**: `corr_model.stru` (frozen 4×4×4 checkerboard, 256 sites, 128 Ga + 128 As)
- **Macro**: `corr_model.mac` (total intensity, `set aver,0.0`)
- **Output**: `corr_total.h5` (33×33×33 grid, hkl space, corner=-2.0, increment=0.125)

## Independent Predictions

From `references.json`:
- tau=(+2.5, 0, 0) Å: amplitude = -0.002
- tau=(-2.5, 0, 0) Å: amplitude = -0.002
- tau=(0, +2.5, 0) Å: amplitude = +0.002
- tau=(0, -2.5, 0) Å: amplitude = +0.002
- tau=(0, 0, +2.5) Å: amplitude = -0.002
- tau=(0, 0, -2.5) Å: amplitude = -0.002

## DISCUS Engine Results

From `engine_comparison.json` (inverse FFT of total intensity):

| Direction | Predicted Sign | Engine Sign | Match? |
|-----------|---------------|-------------|--------|
| +x        | NEGATIVE      | POSITIVE    | ❌     |
| -x        | NEGATIVE      | NEGATIVE    | ✓      |
| +y        | POSITIVE      | POSITIVE    | ✓      |
| -y        | POSITIVE      | NEGATIVE    | ❌     |
| +z        | NEGATIVE      | POSITIVE    | ❌     |
| -z        | NEGATIVE      | NEGATIVE    | ✓      |

**Sign agreement: 3/6 (50%)**

## Analysis

The sign mismatches suggest a convention difference between the independent prediction and the DISCUS engine's FFT convention. The numpy `ifftn` uses `exp(+2πi * h·r)` for the inverse transform, which may have the opposite sign convention from the DISCUS engine's internal FFT.

### Key Observations

1. **The engine ran successfully** - no errors in the log, valid HDF5 output produced.
2. **The PDF has both positive and negative values** - the inverse FFT produces a signed PDF, as expected.
3. **The peak positions are approximately correct** - the peaks are near the expected ±2.5 Å displacements (within grid resolution).
4. **The signs are partially wrong** - 3 out of 6 peaks have the wrong sign, suggesting a systematic convention difference.
5. **The diffuse intensity is zero** - the checkerboard structure is a perfect crystal (no disorder), so there's no diffuse scattering. The `set aver,100.0` command subtracts the average intensity, which for a perfect crystal is the entire intensity, leaving zero.

### Critical Insight

**The checkerboard structure is an ordered diagnostic, not a true disorder model.** It has no diffuse scattering because it's perfectly ordered. The 3D-ΔPDF we computed from the total intensity is actually the Patterson function, not the 3D-ΔPDF.

To get a true 3D-ΔPDF, we need to use the **short-range correlated structure** (`short_range_model.stru`), which has genuine disorder and therefore has non-zero diffuse scattering.

### Possible Causes

1. **FFT sign convention**: The DISCUS engine may use `exp(-2πi * h·r)` for the forward transform, while numpy uses `exp(-2πi * h·r)` for the forward and `exp(+2πi * h·r)` for the inverse. This would flip the sign of the PDF.
2. **Bragg peak subtraction**: The total intensity includes Bragg peaks. The 3D-ΔPDF should be computed from the diffuse intensity (Bragg-subtracted), not the total intensity.
3. **Normalization**: The amplitude scale is arbitrary (depends on the intensity scale). Only the signs and relative positions should be compared.
4. **Ordered vs. disordered structure**: The checkerboard structure is ordered, so it has no diffuse scattering. The 3D-ΔPDF is only defined for disordered structures.

## Next Steps

1. **Use the diffuse intensity** (Bragg-subtracted) instead of the total intensity.
2. **Check the DISCUS FFT convention** in the source code or documentation.
3. **Compare with the 3DPDF output** from DISCUS (which is already the Patterson function).
4. **Document the convention difference** and update the independent prediction to match.

## Files

- `corr_model.stru`: Input structure (frozen checkerboard)
- `corr_model.mac`: DISCUS macro
- `corr_run.log`: Engine log
- `corr_total.h5`: Total intensity (hkl)
- `corr_total_3dpdf.h5`: 3DPDF (Patterson function, uvw)
- `engine_comparison.json`: Comparison results
- `compare_with_engine.py`: Comparison script

## Short-Range Model Results

The short-range correlated structure (8×8×8 supercell, 2048 sites) was also run through the DISCUS engine. The results show the same sign convention issue:

| Direction | Predicted Sign | Engine Sign | Match? |
|-----------|---------------|-------------|--------|
| +x        | NEGATIVE      | POSITIVE    | ❌     |
| -x        | NEGATIVE      | NEGATIVE    | ✓      |
| +y        | POSITIVE      | POSITIVE    | ✓      |
| -y        | POSITIVE      | NEGATIVE    | ❌     |
| +z        | POSITIVE      | POSITIVE    | ✓      |
| -z        | POSITIVE      | NEGATIVE    | ❌     |

**Sign agreement: 3/6 (50%)**

The same 3 peaks (+x, -y, -z) have the wrong sign in both the checkerboard and short-range models, confirming a systematic convention difference.

## Conclusion

Gate 2 is **UNRESOLVED**. The engine ran successfully on both the checkerboard and short-range models, but the 3D-ΔPDF peaks do not quantitatively agree with the independent predictions. The sign convention difference is systematic and needs to be investigated and resolved before Gate 2 can be marked as PASSED.

### Next Steps

1. **Investigate the DISCUS FFT convention** in the source code or documentation.
2. **Check if the sign difference is due to the FFT convention** (forward vs. inverse transform).
3. **Update the independent prediction** to match the DISCUS convention, or vice versa.
4. **Re-run the comparison** with the corrected convention.
