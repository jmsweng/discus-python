# Repository Cleanup Plan for PR #8

## Summary

Commit `6a390d2` inadvertently included unrelated historical files that should not be part of PR #8. This document proposes a safe cleanup plan.

## Files to Remove from PR #8

### 1. Benchmark Archives (~37 MB total)

These are historical evidence from previous PRs and should remain as external evidence, not in the repository:

- `structure_factor_benchmark.tar.gz` (7.2 MB)
- `structure_factor_benchmark_v2.tar.gz` (7.2 MB)
- `structure_factor_benchmark_v3.tar.gz` (7.2 MB)
- `structure_factor_benchmark_v4.tar.gz` (7.2 MB)
- `structure_factor_benchmark_v4.1.tar.gz` (0.6 MB)
- `structure_factor_benchmark_v4.2.tar.gz` (0.6 MB). **Do not remove (2026-10-28):** this is the
  archived C1 evidence, and `benchmarks/discus_6_23_01/gate1_c1_dyadic/c1_dyadic.py` reads the
  archived volumes from it (SHA-256 checked). Keep it until those files are preserved elsewhere
  with the owner's approval.
- `structure_factor_benchmark_v4.zip` (7.2 MB)

**Action:** Remove from the repository. These should be stored as external evidence (e.g., in a separate archive or cloud storage).

### 2. Empty File

- `checksums.txt` (0 bytes)

**Action:** Remove. An empty file serves no purpose.

### 3. Historical Review Documents

These are review documents from previous PRs and should not be part of PR #8:

- `PR2_REVIEW.md`
- `PR3_REVIEW.md`

**Action:** Remove from PR #8. These can be archived separately if needed.

### 4. Benchmark Documentation

These are detailed benchmark reports from previous work:

- `docs/FORM_FACTOR_BENCHMARK.md`
- `docs/FROZEN_C1_CORRECTION_REPORT.md`
- `docs/INDEPENDENT_REVIEW_VERIFICATION.md`
- `docs/OFFCOMMENSURATE_BENCHMARK.md`
- `docs/STRUCTURE_FACTOR_BENCHMARK.md`

**Action:** These may be relevant as historical context, but they are not part of the current PR #8 scope. Consider:
- Option A: Remove from PR #8 and archive separately
- Option B: Keep in the repository but clearly label as historical (not part of PR #8)

### 5. Spike Tooling Instructions

- `spike_discus_3d/PUSH_INSTRUCTIONS.md`

**Action:** The `spike_discus_3d/` directory is historical tooling from PR #6/7. Its retirement is a separate, owner-approved cleanup step. For now, keep the directory but remove the push instructions from PR #8 (they are not relevant to the current work).

## Files to Keep in PR #8

The following files are directly related to PR #8 (scientific validation gates):

- `benchmarks/discus_6_23_01/correlated_model.py` (the corrected model)
- `benchmarks/discus_6_23_01/correlated_run/corr_model.stru` (the frozen structure)
- `benchmarks/discus_6_23_01/correlated_run/references.json` (the independent references)
- `tests/test_correlated_model.py` (the regression tests)
- `docs/SCIENTIFIC_VALIDATION_REPORT.md` (the validation report)

## Proposed Cleanup Steps

1. **Create a new commit** that removes the unrelated files from the PR #8 branch
2. **Do NOT rewrite Git history** (the user explicitly requested no history rewriting)
3. **Add a `.gitignore` entry** if needed to prevent future accidental inclusion of large archives
4. **Document the cleanup** in a commit message explaining which files were removed and why

## Sensitive Information Check

Before removing files, check for any sensitive information:

- [ ] No API keys or credentials in the benchmark archives
- [ ] No personal information in the review documents
- [ ] No proprietary data in the benchmark reports

If any sensitive information is found, handle it according to the project's security policies.

## Approval Required

This cleanup plan requires explicit user approval before proceeding. Do not delete files or rewrite pushed history without permission.

## Next Steps

1. Get user approval for the cleanup plan
2. Create a commit removing the unrelated files
3. Push the cleanup commit to the PR #8 branch
4. Update the PR description to reflect the cleanup
5. Verify that the PR #8 branch now contains only the relevant files
