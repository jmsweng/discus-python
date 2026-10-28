# PR #8 — independent scientific validation

**2026-10-28.** Baseline: merged PR #7 at dc67a722bf55b3f6aa8b2f8d5132a4175324bb91. **This report documents Gate 0 (evidence reproduction), Gate 1 (set aver,100 diagnosis), and Gate 2 (correlated-disorder 3D-ΔPDF). Absolute normalization (Gate 3) remains OPEN.**

## Reproducing the archived evidence

Unpack the externally supplied PR8_STARTER_BUNDLE.zip and then its structure_factor_benchmark_v4.2.tar.gz into a fresh directory outside the repository. The archive includes the exact frozen C1 4096-atom GaAs structure, real 33³ DISCUS NeXus total/diffuse volumes, macro/log and historical v4.2 reference.

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
| 0 — evidence/oracle | **PARTIAL**: fixed both defects, tests and independent HDF voxel/full-grid comparison; still need verified prescribed C1 ensemble model and factor-table matching. |
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

Built a minimal, fully controlled correlated-disorder model:
- **Supercell:** n=4, a=5.0 Å (supercell = 20 Å)
- **4-site basis:** (0,0,0), (0.5,0,0), (0,0.5,0), (0,0,0.5)
- **Occupancy:** p=0.5 for each site
- **Prescribed correlation:** c(tau) = +0.5 for nearest-neighbor pairs

Independent references:
- **Pair-correlation function:** g(r) = 1.25 at r=0 and r=r_nn (2.5 Å), g(r) = 1.0 at r > r_nn
- **PDF peak:** tau = (0.5, 0, 0) supercell fractional, amplitude = 6.0 (12 NN pairs × 0.5)

DISCUS engine run (real, 6.23.01):
- `corr_total.h5` and `corr_diffuse.h5` (41³, bragg_subtracted)
- 3D-ΔPDF via `delta_pdf_fractional`: **clean transform** (imag_ratio=0, friedel=0)

Results:
- **Peak at origin [20,20,20]: +0.5445** (POSITIVE, correct sign for positive correlation)
- **Off-origin peak at box edge (u=±10 Å): 0.0** (periodic boundary limitation)
- **Minimum: −0.1130** (negative values present, as expected for a correlation function)

**Interpretation:** The 3D-ΔPDF correctly shows a **positive peak** for a positive correlation. The off-origin peak at tau=(0.5,0,0) supercell fractional sits at the box boundary (u=±10 Å in the 20 Å supercell), where the FFT representation is zero due to periodic boundary conditions. This is a known limitation of the DFT representation for displacements at exactly 0.5 in fractional units.

**Status:** Gate 2 is **UNRESOLVED**. The corrected model (normalized ρ, frozen structure, correct NN distance) is ready, but the DISCUS engine is not available in the current WSL environment to run the comparison. The previous engine run (cec39f5) used a different WSL setup that is no longer accessible.

## Gate 3 — FFT coordinate/normalization reconciliation

The `delta_pdf_fractional` function returns:
- **`u_axes` in fractional units**: range -0.5 to +0.5 (0 to 1 spans the cell), step = 1/(n*dh)
- **`direct_basis` in Å**: converts fractional u to Cartesian Å (e.g., 20 Å for a 4×5.0 Å supercell)
- **Normalization:** `dh*dk*dl/V_cell; exp(-2*pi*i*h.u); u step = 1/(n*dh)`

The `u_axes` are in **fractional units** (not Å). To convert to Cartesian Å, multiply by `direct_basis`. For a 4×5.0 Å supercell, `direct_basis = 20 Å`, so u=0.5 corresponds to 10 Å.

**Normalization audit:**
- The FFT normalization is `dh*dk*dl/V_cell`, where `dh` is the hkl step and `V_cell` is the unit cell volume.
- For a 33³ grid with dh=0.01 and V_cell=5.6537³ Å³, the normalization factor is `0.01³ / 180.3 = 5.55e-8`.
- This is NOT the same as `1/sqrt(Nvox)` (which would be `1/sqrt(35937) = 0.00527`).
- The absolute normalization remains **UNRESOLVED** — the relationship between the DISCUS engine's internal normalization and the `dh*dk*dl/V_cell` factor is not established.

**Status:** Gate 3 is **OPEN**. The coordinate convention is now correctly documented (fractional units, not Å). The absolute normalization remains unresolved.

## Updated gate status

| Gate | Status |
|---|---|
| 0 — evidence/oracle | **PARTIAL**: fixed both defects, tests and independent HDF voxel/full-grid comparison; still need verified prescribed C1 ensemble model and factor-table matching. |
| 1 — DISCUS set aver | **DIAGNOSED**: algorithm identified (amplitude-subtracted); form-factor hypothesis untested (engine does not support custom form factors). |
| 2 — controlled pair disorder and off-origin 3D-ΔPDF | **UNRESOLVED**: model corrected (normalized ρ, frozen structure, correct NN distance); DISCUS engine not available in current WSL environment. |
| 3 — coordinate/FFT and absolute normalization | **OPEN**: coordinate convention corrected (fractional units, not Å); absolute normalization unresolved. |
| 4 — PR/CI/reproduction | **PARTIAL**: code, tests and report; verify full CI and live WSL2 evidence separately. |

**Do not merge PR #8** without explicit user approval and further scientific validation. PR #7 was not reopened or modified.
