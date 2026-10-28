# Reproduce and independently validate the real 3D disorder benchmarks

**Status: action required on the local WSL2 DISCUS installation.** These instructions address review item **#3**. The prior [`FINDINGS.md`](../spike_discus_3d/FINDINGS.md) describes C0–C5 experiments and reports ≈0.1% agreement with independent scattering calculations. GitHub does **not** contain the full atom lists, intensity outputs, independent-calculation script or all calculation logs behind those reported comparisons. Do not treat the numbers as reverified by this document.

## Before running (operator)

- [ ] Pull the reviewed branch containing the probe-case argument fix. Use the existing Ubuntu/WSL2 environment with the already installed **DISCUS_SUITE 6.23.01**; do not install or replace binaries without approval.
- [ ] Confirm `discus_suite` version, path, SHA-256, operating system, Python and dependency versions. Run `bash spike_discus_3d/probe_environment.sh <output-dir>` or record equivalent facts.
- [ ] Allocate a dedicated output directory outside the repository. Check available disk and memory: the probe's 65³ datasets are reasonably small but generated records should remain outside Git by default.
- [ ] Write down the physics assumptions and *predeclare* numerical comparison metrics/tolerances before adjusting any scattering parameters. Distinguish shape-normalized relative comparisons from absolute normalization.

## Step 1 — regenerate the actual disordered output files

From the repository root, run:

```bash
bash spike_discus_3d/run_cases.sh "$HOME/discus_benchmark_3d" C0 C1 C2 C3 C4 C5
```

- [ ] Confirm every case ran to completion; preserve **all** generated `.mac`, `.stru`, `.log`, `.time.txt`, `.h5`, inspector JSON and `hashes.txt`. The script's archive deliberately excludes HDF5; **keep the original output directory**.
- [ ] Confirm C0 is a perfect crystal, C1/C2 use substituted GaAs with different seeds, C3/C4 test nonorthogonal and skewed grid geometry, and C5 chooses `xray table:inter`. Verify the actual saved crystals/seed behavior; do not infer distinct structures from seed commands alone.
- [ ] Inspect `/entry/data/data_values`, `data_type_axes`, `data_type_with_bragg`, `data_corner`, `data_increment_vector`, `unit_cell_lengths`, `unit_cell_angles` in each NeXus file with `spike_discus_3d/inspect_output.py`.

## Step 2 — resolve reproducibility before benchmarking

- [ ] Run `bash spike_discus_3d/seed_test.sh "$HOME/discus_seed_test"` and compare structure **atom lists** (not just raw text file hashes) from repeated runs at identical seed/thread settings.
- [ ] Record whether the seed command must occur **before or after** `read`, and whether `OMP_NUM_THREADS=1` changes determinism. If deterministic runs cannot be obtained, keep *all* samples and characterize variability; do not compare nonidentical structures as if they were the same.
- [ ] Record individual species counts, substitution fraction, cell, temperature/displacement settings, form-factor table, and exact input hashes for each run.

## Step 3 — independently reproduce the scattering

- [ ] Implement a separate, reviewable Python script that loads **the actual saved atomic coordinates** (`.stru` from C1/C5), computes `F(hkl)=Σ_j f_j(Q) exp(2π i h·r_j)`, and evaluates `I_total=|F|²` at exactly the stored NeXus grid points. Use an explicitly documented scattering-factor source, Debye-Waller assumptions, form-factor lookup and atom-coordinate convention.
- [ ] Independently calculate the appropriate average-structure amplitude and **subtract amplitudes before squaring** to compare to the actual `set aver,100` diffuse intensity. Explain the finite-crystal box/interference factor and confirm the correct species/occupancy averages.
- [ ] Compare C1 and C5 using the **same saved structure** as the output being benchmarked. Report relative L2 error, correlation coefficient, scale ratio, worst residual voxel with hkl coordinates, and a few representative peaks. Compare both total and Bragg-subtracted diffuse intensity. Report the explained form-factor-table difference and any remaining unexplained residual honestly.
- [ ] Preserve the independent implementation and machine-readable comparison outputs in a new `benchmarks/` directory; do not commit large generated raw outputs without deciding on size/LFS/artifact storage.

## Step 4 — validate ΔPDF physics, not just an FFT implementation

- [ ] Construct a disorder model with a **known nonzero pair correlation** (e.g. a controlled nearest-neighbor occupation relation) and preserve its atomic coordinates. Random uncorrelated antisites alone are not enough to test nonzero off-origin ΔPDF peaks.
- [ ] Compute diffuse scattering with Bragg-subtracted `set aver,100`, load it with `read_discus_volume`, obtain `delta_pdf_fractional`, and compare off-origin peak positions **and signed/relative amplitudes** with the known pair correlations and an independent direct summation/reference result.
- [ ] Compare DISCUS's `value 3DPDF` with the Python FFT only after explicitly reconciling its normalization (prior report: `1/sqrt(Nvox)` relationship), the reported 3DPDF header step versus `1/(N Δh)`, and handling of Bragg-present data.
- [ ] Confirm nonorthogonal distances using the full cell metric, assess finite reciprocal-range truncation, edge/windowing effects and masking/uncertainty. State explicitly if absolute normalization remains unresolved.

## Step 5 — package and report the evidence

- [ ] Commit compact **inputs**, a reproducible `run_benchmark.py`/notebook, small representative fixtures, expected outputs, tests and a report that names the exact DISCUS executable version/SHA, OS, run commands, seeds, file hashes, elapsed time, peak memory, and library versions.
- [ ] Attach large immutable result bundles to an approved artifact location; provide their checksums and a retrieval procedure in the report. Avoid checking binary megabytes into regular Git history by default.
- [ ] Give results in a table: predicted value, DISCUS observation, independently calculated value, tolerance, pass/fail/uncertainty for each quantitative claim. Distinguish **engine-ran**, **schema-validated**, **numerically benchmarked**, and **scientifically interpreted** status.
- [ ] Open a dedicated PR and request review; leave uncertain/failed scientific tests in the report instead of widening tolerances to pass.

### Clear completion criterion

#3 is only complete when another scientist can rerun the saved input/model and independent reference calculation, reproduce the numbers within declared tolerances, and explain remaining normalization uncertainty. A green CI check from synthetic files or a Bragg-present C3 Patterson peak alone does not close #3.
