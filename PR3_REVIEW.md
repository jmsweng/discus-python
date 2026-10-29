> **Historical document (kept by owner decision, 2026-10-28).** Written during earlier PR work. It is **not** part of PR #8's validated results and may contain claims since superseded or refuted. The current record is `docs/SCIENTIFIC_VALIDATION_REPORT.md` and `benchmarks/discus_6_23_01/REVIEW_2026-10-28.md`.

# PR #3 Review: FFT Diagnostics and Testing

**Status:** Ready to close
**Review date:** 2026-10-22

## Summary

PR #3 transforms `delta_pdf_from_diffuse` from returning a simple `(volume, r_axes)` tuple to returning a rich `FFTResult` dataclass with `complex_volume`, `r_axes`, `diagnostics` (an `FFTDiagnostics` with `imaginary_ratio`, `max_imaginary`, `max_real`, `friedel_relative_error`, `dq`, `reciprocal_voxel_volume`, `normalization`), and a `real_volume` convenience property. It adds configurable `imaginary_tolerance` (default `1e-4`) and `raise_on_imaginary` (default `True`), a Friedel symmetry diagnostic, `ReciprocalVolume` integration with a Bragg metadata gate, 8 new analytical tests, and updates notebook 2 to the new API.

## What Was Tested

| Test | Result |
|---|---|
| `tests/test_fft_diagnostics.py` (8 tests) | ✅ All pass |
| `tests/test_volume.py` (2 tests, updated) | ✅ All pass |
| Analytical FFT verification (6 scenarios) | ✅ All pass |
| Notebook 2 workflow with new API | ✅ Works |
| Full test suite (31 tests) | ✅ All pass |

## Analytical FFT Verification

1. **Constant intensity** → Peak at center with value `n³·dq³/(2π)³` ✅ (matches analytical)
2. **Cosine pair** `2·cos(Qx·r0)` → Symmetric peaks at ±r0, zero at center ✅
3. **Unequal axis steps** → Correct r-axis spacing `2π/(n·dq)` per axis ✅
4. **Non-centrosymmetric** (single off-center voxel) → Complex FFT, `ValueError` raised when ratio > tolerance ✅
5. **`ReciprocalVolume` metadata gate** → Accepts Bragg-excluded, rejects missing metadata ✅
6. **Backward compatibility** → `FFTResult` has `real_volume`, `r_axes`, `diagnostics` ✅

## Issues Found (All Minor, No Fix Required)

1. **`from dataclasses import dataclass` is imported mid-file** (in `volume.py`, after `read_hdf5_volume`). Should be at the top with other imports. Style issue only.
2. **Friedel symmetry diagnostic assumes periodic FFT-index inversion.** Correctly documented as an assumption-dependent check. Acceptable.
3. **`imaginary_tolerance` default of `1e-4` may be strict for noisy experimental data.** The `raise_on_imaginary=False` option provides an escape hatch. Acceptable for MVP.
4. **Notebook 2 source cells have inconsistent line-ending formatting** (some `\n`-joined single strings, others multi-line arrays). Cosmetic. Notebook executes correctly.

## Scientific Correctness

- **FFT transform:** ✅ Correct — `ifftshift` → `fftn` → `fftshift` with `exp(-i Q.r)` convention and `dQx·dQy·dQz/(2π)³` normalization
- **Friedel symmetry:** ✅ Correct diagnostic for centrosymmetric data
- **Imaginary component check:** ✅ Physically meaningful — centrosymmetric I(Q) should give a nearly real ΔPDF
- **3D-ΔPDF workflow:** ⚠️ **Not scientifically validated** — synthetic demonstration only. No real 3D DISCUS calculation has been benchmarked.

## Compatibility

- **API change:** `delta_pdf_from_diffuse` now returns `FFTResult` instead of `(ndarray, tuple)`
- **Backward compatibility:** `FFTResult.real_volume` and `FFTResult.r_axes` provide the same data as the old tuple return
- **Notebook 2:** Updated to the new API ✅
- **`plot_delta_pdf_slice`:** Unchanged, works with `FFTResult.real_volume` ✅

## Verdict

**Ready to close.** The FFT diagnostics are scientifically correct, well-tested (8 analytical tests), and the API change is backward-compatible. All 31 tests pass. No code fixes are needed.
