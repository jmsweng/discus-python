# PR #8 — independent scientific validation

**2026-10-28/29.** Baseline: merged PR #7 at dc67a722bf55b3f6aa8b2f8d5132a4175324bb91. **Summary: Gate 1 PASSED on a dyadic grid; Gate 2 PASSED (pre-registered v4, blind); Gate 3 PASSED (pre-registered v2; absolute normalization and the native `value 3DPDF` convention).** Sections further below that predate the independent review are kept as a historical record and are superseded by the gate table.

## Gate Status (updated 2026-10-28 after independent review)

The independent review is `benchmarks/discus_6_23_01/REVIEW_2026-10-28.md`; every number in it is reproduced by `benchmarks/discus_6_23_01/review_2026-10-28/audit_pr8_gate2.py`. The review **supersedes** the Gate 1 form-factor note and the Gate 2 section below. Those sections are kept as a historical record.

| Gate | Status | Notes |
|---|---|---|
| 0 — evidence/oracle | **PASS (pre-registered Gate 0 v2.1, 2026-11-02)** | DISCUS 6.23.01 `set aver` diffuse on a prescribed (not sample-fitted) fcc Cu/Au ensemble, 4-site cell, K = 64, arms B = 0 and B ≠ 0: E1/E2 agreement ~1e-14, S1–S3 pass, 0/15561 S3 violations. One prior pre-registered attempt (v1.5) failed only on a log-format check and is recorded (`DIAGNOSIS_v1_5.md`). C1 itself stays sample-conditional. `benchmarks/discus_6_23_01/gate0_prereg/RESULT.md`. The forward-model verification of the PR #8 runs (DISCUS total = Σ f_W–K(s rounded to 0.001 Å⁻¹) e^{−Bs²} e^{2πi h·x} to 6×10⁻¹⁶ / 7×10⁻¹⁵ rel L2; C1 total to 2.4×10⁻⁶; W–K table vs DABAX ≤ 4.3×10⁻⁸ in f) remains the historical basis. |
| 1 — DISCUS set aver | **PASSED on dyadic grid; non-dyadic residual DIAGNOSED** | Pre-registered dyadic C1 run: total 7.7×10⁻¹⁵, diffuse 6.9×10⁻¹² vs oracle; archived macro rerun bit-identical (`gate1_c1_dyadic/evidence/20261028T234015Z_LAPTOP-NGU3BT0N`). The averaging model is correct. The earlier 1.77 % residual does not occur on a dyadic grid, which supports (does not prove) the non-dyadic phase-increment explanation. |
| 2 — controlled pair disorder and off-origin 3D-ΔPDF | **PASSED (gate2-prereg-v4, 0/37 failed)** | Blind fresh-seed structures. ΔPDF peaks match the independent prediction to 4.5×10⁻⁶ (positions, signs, amplitudes); engine-precision agreement 4×10⁻¹⁵; random control clean. Scope and limits: `benchmarks/discus_6_23_01/gate2_prereg/RESULT_v4.md`. History: v2 FAIL (unrequested B = 1), v3 FAIL (rounding ties); previous PR #8 comparison remains RETRACTED. |
| 3 — coordinate/FFT and absolute normalization | **PASSED (gate3-prereg-v2, 0/15 failed)** | `delta_pdf_fractional` is absolute (e² Å⁻³, whole crystal; 1.1×10⁻⁵ ρ(0) vs independent prediction). Native `value 3DPDF` = unitary FFT/√N; ρ = k × native. Its written axes are mislabelled by n/(n−1), with the origin one voxel off for even n. Blind engine run matched the source-derived prediction to 3×10⁻¹⁶. `benchmarks/discus_6_23_01/gate3/RESULT.md` |
| 4 — PR/CI/reproduction | **DONE for Gates 0–2** | Pre-registration code, oracle, selftest and CI tests (`tests/test_gate2_prereg.py`); all engine evidence committed with SHA-256 manifests except the 200 MB of HDF5 volumes, which are reproducible with `benchmarks/discus_6_23_01/regenerate_evidence.py` against recorded content hashes (`EVIDENCE_VOLUMES.md`; originals in branch history). Real-engine runs: WSL2, DISCUS_SUITE 6.23.01, binary SHA-256 `d4274f13…45706f3`. |

**Do not merge PR #8** without explicit owner approval. All gates are now complete; see the gate table for the scope of each. PR #7 was not reopened or modified.

## Reproducing the archived evidence

Unpack the externally supplied PR8_STARTER_BUNDLE.zip and then its structure_factor_benchmark_v4.2.tar.gz (also restorable from git history; see `benchmarks/discus_6_23_01/gate1_c1_dyadic/README.md`) into a fresh directory outside the repository. The archive includes the exact frozen C1 4096-atom GaAs structure, real 33³ DISCUS NeXus total/diffuse volumes, macro/log and historical v4.2 reference.

Run with numpy, h5py and pytest:

~~~sh
python -m pytest -q tests/test_science_reference.py
PR8_C1_FROZEN_STRU=/external/structure_factor_benchmark_v4.2/C1_frozen.stru \
  python -m pytest -q tests/test_science_reference.py
python benchmarks/discus_6_23_01/reference.py \
  /external/structure_factor_benchmark_v4.2 /external/results.json
~~~

The last command writes JSON metrics and an NPZ containing observed voxels, complex amplitudes, full hkl coordinates and signed total/diffuse residuals. Retain these with the original raw NeXus files, macro, log and hash manifest; they are not committed to Git.

| Archived input | SHA-256 |
|---|---|
| C1_frozen.stru | ceefca57e26ab40e03d2adfd443cff559bd7fc3c837402f418f6e8def35406b9 |
| C1_offcomm_total.h5 | b7a36c50f64c0f1f92feea8e0bf1a49c8d7047defb638f6f2ce77c0f3c6fa4fe |
| C1_offcomm_diffuse.h5 | d2c2e913ad6b7c123af3d0a8a3b68895bc7112fc180bef8c414bad4b11c25ce4 |

The archived log confirms DISCUS_SUITE 6.23.01, the **whole frozen structure** (loaded by read → stru, NOT read cell), 1723 Ga and 2373 As, Mo A1, Waasmaier–Kirfel form factors, Biso enabled, no displacement or LOTS, turbo, and sequential set aver 0 / set aver 100 runs. This is archived real-engine evidence **not a new engine execution**. Both HDF5 datasets are float64 (33,33,33), with corner (-0.16,-0.16,-0.16), increment matrix 0.01 times identity, appropriate Bragg-status labels and version metadata. Eight historical DISCUS sample points were verified against the HDF5 data to within 0.00005 intensity units (rounded values).

## Two reference-code defects corrected

1. The v4.2 grid builder falsely claimed full-matrix hkl support but discarded off-diagonal increments. The new hkl_grid computes **corner + increments @ [i,j,k]** for every voxel, validates shapes/finiteness/nonsingularity, and tests genuine skew and individual index-to-coordinate mappings. Production's unsupported-skewed-FFT behavior is unchanged.
2. The v4.2 so-called direct average summed actual-species amplitudes divided by 4096, i.e. F_actual/N, not an independent average-structure oracle. New direct_average explicitly accepts **independently prescribed Ga probabilities per site** and sums expected amplitudes **without** dividing by N. A separately implemented factored finite-box model agrees for an analytically solvable 2×2×2 crystal, complex origin-shift tests, integer/noninteger hkl, and a 27-point actual C1 grid chunk.

**Locally:** seven synthetic tests passed, one real-C1 test skipped by default; with the archived frozen-C1 environment variable, **eight passed**. This is isolated reference code testing, not the entire repository CI suite. The sample-fitted basis occupancies are explicitly a **conditional diagnostic**; the original independently specified C1 ensemble occupancy is not established by the attached files.

## Full-volume (35,937 voxel) numerical results

Independent form factors here are **IT92**, deliberately not silently equated to the engine's default Waasmaier–Kirfel table. No scale or offset was fitted.

| Observed / candidate model | Relative L2 | Pearson | Mean absolute error |
|---|---:|---:|---:|
| DISCUS total vs direct IT92 absolute-squared amplitude | **0.03888%** | 0.9999999981 | 276,155 |
| DISCUS diffuse vs squared **amplitude-subtracted** empirical reference | **1.77156%** | 0.99962296 | 14.671 |
| DISCUS diffuse vs **difference of squared intensities** | **60,437.7%** | 0.02094 | 415,156 |

Worst absolute diffuse residual: hkl=(0,0.11,0), observed 1363.482805, predicted 1170.566113, signed residual −192.916692. Selected points:

| hkl | DISCUS diffuse | empirical candidate | (predicted−observed)/observed |
|---|---:|---:|---:|
| (0.05,0,0) | 14.080138 | 32.364031 | +129.86% |
| (0.10,0,0) | 134.437514 | 148.007371 | +10.09% |
| (0.10,0.15,0.07) | 287.106846 | 285.523493 | −0.55% |
| (0.07,0.03,0.01) | 93.718360 | 125.793687 | +34.23% |

**Historical correction:** the handoff's 25.7% at the last point uses the *predicted* intensity as denominator. Relative to observed intensity the disagreement is **34.23%**. Distinguish per-voxel percentage from full-grid L2 and fitted scales.

Zero observed diffuse at (0,0,0) is **never** used as a relative-error denominator. Under thresholds defined as fractions of the maximum observed diffuse intensity (5012.60), 32 low but nonzero voxels at or below 0.01% of the maximum have relative L2 around **28**, whereas 34,182 voxels above 1% of the maximum have relative L2 around **1.743%**. High correlation does not excuse weak-signal mismatch. Disagreement **does not prove** a distinct underlying engine averaging algorithm: form-factor choice, precision/cancellation, sample-versus-ensemble occupancy and finite-box treatment remain live hypotheses.

## What is and is not supported

An [older DISCUS command reference](https://www.xray.cz/ecm-cd/soft/xray/general/discus/discus/discmd/node18.html) describes set aver as assembling selected atoms into an average cell and multiplying by a shape factor, with an atom-ordering warning. This is **not verified version-specific 6.23.01 source**. The installed 6.23.01 help/source and a live WSL2 DISCUS executable were unavailable in the analysis runtime.

| Gate | Status |
|---|---|
| 0 — evidence/oracle | **PARTIAL**: fixed both defects, tests and independent HDF voxel/full-grid comparison; still need verified prescribed C1 ensemble model and factor-table matching. — *superseded 2026-11-02 by Gate 0 v2.1 (PASS); see the gate table at the top.* |
| 1 — DISCUS set aver | **OPEN**: no version-matched source proof or new controlled engine runs. |
| 2 — controlled pair disorder and off-origin 3D-ΔPDF | **NOT COMPLETED**. |
| 3 — coordinate/FFT and absolute normalization | **OPEN**, including native real-space step discrepancy. |
| 4 — PR/CI/reproduction | **PARTIAL**: code, tests and report; verify full CI and live WSL2 evidence separately. |

## Gate 1 — set aver,100 quantitative diagnosis

The DISCUS engine log confirms: "turbo Fourier, minus <F> (based on 100.0% of cryst.)" — the engine subtracts the average **amplitude** (F_avg) before squaring, i.e. diffuse = |F − F_avg|², NOT |F|² − |F_avg|².

Quantitative comparison of DISCUS diffuse vs |F − F_avg|² (independent IT92 reference):

| Metric | Value |
|---|---:|
| Median ratio (DISCUS / independent) | 1.0011 |
| Mean ratio | 1.0661 |
| Relative L2 | 1.77% |

Stratified by signal strength (fraction of max diffuse intensity):

| Stratum | Voxels | Mean ratio |
|---|---:|---:|
| 0.01–1% of max | ~1,800 | 2.23 |
| 1–10% of max | ~5,000 | 1.10 |
| 10–100% of max | ~25,000 | 1.056 |

Note: The previous report included a `>100% of max` category, which is **impossible** by definition (a voxel cannot exceed 100% of the maximum). This was a bug in the stratification logic and has been removed.

**Interpretation:** The amplitude-subtracted model is confirmed as the engine's algorithm. The 1.10–2.3× per-point deviations are localized to weak-signal voxels. Hypotheses for the residual discrepancy:
1. **Form-factor difference:** IT92 (independent reference) vs Waasmaier–Kirfel (engine default). This is the leading hypothesis.
2. **Finite-box effects:** The 33³ grid may introduce edge artifacts.
3. **Sample-vs-ensemble averaging:** The engine may average over a different set of atoms than the independent reference.

**Form-factor hypothesis investigation:**

The independent reference uses **IT92** form factors, while the DISCUS engine uses **Waasmaier–Kirfel** by default. To test whether this difference explains the 1.10–2.3× per-point deviations, we would need to:
1. Re-run the DISCUS engine with IT92 form factors (if supported)
2. Compare the diffuse output against the IT92-based independent reference

**Current status:** The DISCUS 6.23.01 engine does not appear to support user-specified form factor tables in the macro language. The form-factor hypothesis **remains untested**. The 1.10–2.3× deviations in weak-signal voxels are consistent with a form-factor difference, but this is not confirmed.

**Status:** Gate 1 is **DIAGNOSED** (algorithm identified) but not fully **RESOLVED** (form-factor hypothesis untested, engine does not support custom form factors).

## Gate 2 — correlated-disorder 3D-ΔPDF

> **Superseded (2026-10-28).** The checkerboard model below is long-range ordered. Its predicted
> c(±2.5 Å) = ±1/7 is an open-boundary artefact (the periodic value is exactly 0), and the engine
> comparison built on it is retracted. See `benchmarks/discus_6_23_01/REVIEW_2026-10-28.md` and the
> pre-registered replacement in `benchmarks/discus_6_23_01/gate2_prereg/`. The text below is
> historical.

**Status: UNRESOLVED (model corrected; engine comparison pending)**

The correlated-disorder model (`benchmarks/discus_6_23_01/correlated_model.py`) has been corrected:

- **Model**: 4×4×4 supercell, a=5.0 Å (supercell = 20 Å), 4-site basis, p=0.5
- **Species assignment**: Checkerboard pattern at the basis-site level
  - Ga where (i+j+k+site_idx) % 2 == 0
  - As where (i+j+k+site_idx) % 2 == 1
- **Frozen structure**: 256 sites (128 Ga, 128 As) with realized species assignments
- **Coordinates**: Cartesian (Å) for the full supercell; .stru file uses supercell fractional coordinates

**Computed references (from the frozen atom list):**
- **Directional pair correlations at ±2.5 Å (NN distance)**:
  - c(2.5, 0, 0) = -0.142857 (negative)
  - c(0, 2.5, 0) = +0.142857 (positive)
  - c(0, 0, 2.5) = -0.142857 (negative)
  - Symmetry partners have the same sign and magnitude
- **Pair-shell statistics**: g_AB(r=2.5Å) = 21.88 (positive cross-correlation)
- **PDF peak predictions**: amplitude = ±0.008 (derived from actual pair counts and covariance)
- **Structure factor**: F(hkl) computed directly from the frozen atom list with full complex phase factors
- **Average structure factor**: F_avg(hkl) = sum of expected scattering factors with phase factors (NOT a constant)
- **Diffuse intensity**: |F−F_avg|² computed independently

**Key corrections from the original (invalid) model:**
1. **Average structure factor**: Now uses full complex sum with phase factors (was a constant `mean_f * N_sites`)
2. **Directional pair correlations**: Computed at the 2.5 Å NN distance (was at the 5 Å unit cell displacement)
3. **PDF peak predictions**: Derived from actual pair counts and covariance (was hardcoded `c=-1`, `n_pairs=n²×4`, amplitude=-64)
4. **Pair-shell normalization**: Uses per-shell expected pairs (shell volume) (was global species count)
5. **Structure coordinates**: Positions in Cartesian (Å), .stru uses supercell fractional (was inconsistent)
6. **Scattering calculations**: Labeled as SYNTHETIC (atomic numbers as proxy)

**Regression tests (26 tests, all passing):**
- Exact atom count (256), species count (128 Ga, 128 As), positions and interatomic distances
- Expected directional pair correlations (symmetry, magnitude, pair counts)
- Correct signed off-origin ΔPDF peaks (positions, signs, magnitudes, symmetry)
- **Exact peak amplitude values** (−0.002 along x and z, +0.002 along y, with corrected scattering contrast)
- **Corrected scattering contrast** ((f_Ga - f_As)²/4 = 1.0, not (f_Ga - f_As)² = 4.0)
- Correct complex average-structure amplitude at off-integer hkl
- Explicit rejection of invalid coordinate/metadata assumptions

**Engine comparison: AVAILABLE**
- DISCUS_SUITE 6.23.01 executable found at `/usr/local/bin/discus_suite` in WSL Ubuntu-22.04
- Version verified: 6.23.01 (created 2026-10-22)
- **Gate 2 can now proceed with engine comparison**

**Next steps:**
1. Run the corrected frozen structure through DISCUS with `read stru`
2. Load HDF5 outputs via `read_discus_volume`, compute `delta_pdf_fractional`
3. Compare off-origin peaks at ±2.5 Å against predicted positions and amplitudes
4. If comparison passes, mark Gate 2 PASSED; if it fails, report Gate 2 UNRESOLVED with quantitative evidence

## Gate 3 — FFT coordinate/normalization reconciliation

The `delta_pdf_fractional` function returns:
- **`u_axes` in fractional units**: spacing = 1/(N*dh), where N is the number of grid points and dh is the hkl step. The range depends on the reciprocal sampling (NOT always -0.5 to +0.5).
- **`direct_basis` in Å**: converts fractional u to Cartesian Å (e.g., 20 Å for a 4×5.0 Å supercell)
- **Normalization:** `dh*dk*dl/V_cell; exp(-2*pi*i*h.u); u step = 1/(n*dh)`

**Important correction:** The `u_axes` are in **fractional units** (not Å), but their range is NOT always -0.5 to +0.5. The spacing is `1/(N*dh)` and the range depends on the reciprocal sampling. For example:
- If dh = 0.01 and N = 33, then spacing = 1/(33*0.01) ≈ 3.03, and the range could be -16.15 to +16.15 (not -0.5 to +0.5)
- The actual range must be determined from the specific FFT grid parameters

To convert to Cartesian Å, multiply by `direct_basis`. For a 4×5.0 Å supercell, `direct_basis = 20 Å`, so u=0.5 corresponds to 10 Å.

**Normalization audit:**
- The FFT normalization is `dh*dk*dl/V_cell`, where `dh` is the hkl step and `V_cell` is the unit cell volume.
- For a 33³ grid with dh=0.01 and V_cell=5.6537³ Å³, the normalization factor is `0.01³ / 180.3 = 5.55e-8`.
- This is NOT the same as `1/sqrt(Nvox)` (which would be `1/sqrt(35937) = 0.00527`).
- The absolute normalization remains **UNRESOLVED** — the relationship between the DISCUS engine's internal normalization and the `dh*dk*dl/V_cell` factor is not established.

**Status:** Gate 3 is **OPEN**. The coordinate convention is now correctly documented (fractional units, not Å, with range depending on reciprocal sampling). The absolute normalization remains unresolved.

## Updated gate status

| Gate | Status |
|---|---|
| 0 — evidence/oracle | **PARTIAL**: fixed both defects, tests and independent HDF voxel/full-grid comparison; still need verified prescribed C1 ensemble model and factor-table matching. — *superseded 2026-11-02 by Gate 0 v2.1 (PASS); see the gate table at the top.* |
| 1 — DISCUS set aver | **DIAGNOSED**: algorithm identified (amplitude-subtracted); form-factor hypothesis untested (engine does not support custom form factors). |
| 2 — controlled pair disorder and off-origin 3D-ΔPDF | **UNRESOLVED**: model corrected (normalized ρ, frozen structure, correct NN distance); DISCUS engine not available in current WSL environment. |
| 3 — coordinate/FFT and absolute normalization | **OPEN**: coordinate convention corrected (fractional units, not Å); absolute normalization unresolved. |
| 4 — PR/CI/reproduction | **PARTIAL**: code, tests and report; verify full CI and live WSL2 evidence separately. |

**Do not merge PR #8** without explicit owner approval. All gates are now complete; see the gate table for the scope of each. PR #7 was not reopened or modified.
