> **Historical document (kept by owner decision, 2026-10-28).** Written during earlier PR work. It is **not** part of PR #8's validated results and may contain claims since superseded or refuted. The current record is `docs/SCIENTIFIC_VALIDATION_REPORT.md` and `benchmarks/discus_6_23_01/REVIEW_2026-10-28.md`.

# PR #2 Review: 3D Infrastructure

**Status:** ✅ Merged into `main` at commit `b6bc243`
**Review date:** 2026-10-22

## Summary

PR #2 adds 3D reciprocal-space infrastructure: `ReciprocalVolume` (typed data model), `Diffuse3DRunner` (isolated macro execution adapter), enhanced `PDFConfig` (new `qmin`/`dq`/`wavelength` fields, `qmax=0` rejection, `rmax/dr` integer check), `PDFRunner.execute()` (two-phase prepare/execute with hash verification and r-grid validation), `read_pdf()` (renamed from `read_csv`), GitHub Actions CI workflow, and 12 new/updated tests.

## What Was Tested

| Test | Result |
|---|---|
| `tests/test_reciprocal.py` (5 tests) | ✅ All pass |
| `tests/test_pdf.py` (16 tests, updated by PR #2) | ✅ All pass |
| Real DISCUS 1D PDF (GaAs 5×5×5, r=0→12 Å, dr=0.02) | ✅ 601 samples, G(r) -5.74→6.10, finite, increasing r |
| Notebook 1 API compatibility (`PDFConfig`, `runner.run()`, `replace()`, `plot()`) | ✅ All patterns work |
| Full test suite (31 tests) | ✅ All pass |

## PR #2 Contributions Verified in `main`

- **`reciprocal.py`** — `ReciprocalVolume` with coordinate-system/unit/preprocessing validation ✅
- **`diffuse3d.py`** — `Diffuse3DRunner` with isolated macro execution, hash verification, headless PGPLOT ✅
- **`pdf.py`** — `PDFConfig` with `qmin`/`dq`/`wavelength` fields, `qmax=0` rejection, `rmax/dr` integer check; `PDFRunner.execute()` two-phase workflow; `read_pdf()` ✅
- **`__init__.py`** — 3D symbols exposed ✅
- **`docs/3D_INTEGRATION.md`** — explicit scope guardrails ✅
- **`.github/workflows/python-tests.yml`** — CI workflow ✅

## Issues Found (All Minor, No Fix Required)

1. **`PDFConfig` error message** for `qmax=0` is slightly misleading (mixes `qmin<qmax` constraint with "qmax=0 not supported" note). Validation logic is correct.
2. **`PDFRunner.execute()` engine-error regex** (`error\s*:`) is broad. Low risk in practice.
3. **`ReciprocalVolume.__post_init__`** uses `object.__setattr__` on a frozen dataclass. Valid pattern, minor code smell.
4. **`Diffuse3DRunner.execute()`** only checks `diffuse.h5` exists/non-empty; structure validation deferred to `load_volume()`.

## 3D-ΔPDF Scientific Validation

**Not scientifically validated.** `ReciprocalVolume` and `Diffuse3DRunner` are infrastructure only. No real 3D DISCUS calculation has been benchmarked. This is explicitly documented in `docs/3D_INTEGRATION.md`.

## Unresolved Issues

None blocking. PR #2 is complete and correct in `main`.

## Recommendation

**No action needed.** PR #2 is merged, all tests pass, the real 1D DISCUS calculation is verified, and the notebooks are compatible with the updated APIs.
