# Development task tracker

**Review snapshot:** 2026-10-22 (America/New_York) · **Evidence:** repository `main` at [`2e74868`](https://github.com/jmsweng/discus-python/commit/2e74868e3a05ec7b84aff8bceea85a21c309fd72), [PRs](https://github.com/jmsweng/discus-python/pulls?q=is%3Apr) and GitHub Actions. **This is a living checklist, not a promise that an agent is working on a task.**

- `[x]` means the specified deliverable has been **implemented or recorded with evidence**, **not** that all related science is verified.
- `[ ]` means still needed, insufficiently verified, or only partly built. When some groundwork exists, it is explicitly described.
- Keep each item small and independently verifiable. Mark complete only after linking the implementation/tests and, for scientific claims, actual validation evidence.
- Package A–D documents below remain the detailed specifications for the Javelin-inspired additions. This page is the single progress overview; don't duplicate a separate status table in every document.

## Completed foundations — source and integration evidence

### DISCUS / 1D PDF

- [x] Establish working DISCUS_SUITE 6.23.01 command line and the verified 1D powder-to-PDF macro in [`docs/VERIFIED_DISCUS_MACRO.md`](../docs/VERIFIED_DISCUS_MACRO.md). This records prior local verification, **not** an engine test that GitHub Actions can run.
- [x] Implement `PDFConfig`, inspected macro generation and an isolated `PDFRunner` with per-run inputs, hashes, macro, logs and status. [`src/discuspdf/pdf.py`](../src/discuspdf/pdf.py); [PR #1](https://github.com/jmsweng/discus-python/pull/1).
- [x] Expose and validate `qmin`, `qmax`, `dq`, `qsig`, wavelength and real-space grid settings. Ambiguous `qmax=0` is rejected. [`pdf.py`](../src/discuspdf/pdf.py).
- [x] Run the exact prepared 1D macro via `execute(prepared)`; validate hashes, engine failure, expected radial sample count and grid before success. [`pdf.py`](../src/discuspdf/pdf.py).
- [x] Commit a real GaAs 1D sample consisting of a structure, macro, log and numeric `result.pdf`. [`examples/run_sample/`](../examples/run_sample/).
- [x] Add a **locally opt-in** DISCUS-vs-reference 1D regression test. It is skipped in default CI unless `DISCUS_INTEGRATION=1` and the binary is available; fixture consistency is an open item below. [`tests/test_integration_discus.py`](../tests/test_integration_discus.py).

### 3D infrastructure and FFT prototype

- [x] Add `Diffuse3DRunner` to stage prevalidated input macros and assets and execute them in isolated directories. **This does not generate or verify a real DISCUS 3D macro.** [`src/discuspdf/diffuse3d.py`](../src/discuspdf/diffuse3d.py); [PR #2 integration history](../docs/INTEGRATION_HISTORY.md).
- [x] Add structural HDF5 validation with optional `OutputSpec` (expected datasets, 3D shape, 1D axis lengths, finite numbers), log-error detection, input-preflight checks and `produced_unverified` status where no spec is provided. **Physical units and Bragg subtraction are not verified.** [PR #4](https://github.com/jmsweng/discus-python/pull/4); [`docs/3D_INTEGRATION.md`](../docs/3D_INTEGRATION.md).
- [x] Add a `ReciprocalVolume` container that requires declared Cartesian-Q units, centered/uniform axes and explicit Bragg-preprocessing labels. **These declarations are not verified against source-file metadata.** [`src/discuspdf/reciprocal.py`](../src/discuspdf/reciprocal.py).
- [x] Implement the synthetic, uniform-Cartesian-Q FFT prototype returning a complex `FFTResult`, real-space axes and numerical diagnostics, plus a plotting helper. [`src/discuspdf/volume.py`](../src/discuspdf/volume.py).
- [x] Add configurable imaginary-component tolerance, diagnostic reporting and `FFTValidationError.result` to retain the computed complex data on failure. [PR #4](https://github.com/jmsweng/discus-python/pull/4).
- [x] Add analytical and negative-case Python tests for FFT behavior, basic metadata, HDF5 corruption, missing datasets, and run-status transitions. [`tests/`](../tests/).
- [x] Keep the 3D Jupyter notebook explicitly **synthetic only**. [`notebooks/02_3d_delta_pdf.ipynb`](../notebooks/02_3d_delta_pdf.ipynb).

### Engineering and documentation

- [x] Add GitHub Actions Python unit tests (latest checked `main` run was successful). **Standard CI does not install or run DISCUS.** [workflow](../.github/workflows/python-tests.yml).
- [x] Update the 1D quickstart to execute the already inspected prepared job, fix notebooks, ignore `runs_3d/`, and repair tracked/local-only setup document references. [PR #4](https://github.com/jmsweng/discus-python/pull/4).
- [x] Record accurate ancestry of PRs #1–#3 without rewriting shared Git history. [`docs/INTEGRATION_HISTORY.md`](../docs/INTEGRATION_HISTORY.md).
- [x] Write the four Javelin-inspired [package briefs](README.md) as future requirements. **Writing plans does not complete their requested features.**

## Next steps — prioritize these before claiming real 3D-ΔPDF capability

### Reference and scientific reliability

- [ ] **P0 — Reconcile the GaAs reference fixture with its opt-in integration test.** The archived [`config.json`](../examples/run_sample/config.json) and [`run.mac`](../examples/run_sample/run.mac) specify `qmax=30`, but [`tests/test_integration_discus.py`](../tests/test_integration_discus.py) runs `qmax=25`; the archived macro uses `outf result.csv` while the retained artifact is `result.pdf`. Establish which output is authoritative, regenerate or correct the fixture/provenance, and run the opt-in comparison under DISCUS 6.23.01; **do not change numerical tolerances just to make the test pass**.
- [ ] **P1 — Verify powder settings against real DISCUS outputs**, including `set delta`/Gaussian width, Qmin, Qmax, DQ, wavelength and termination behavior. Record test inputs and numeric differences.
- [ ] **P1 — Validate FFT normalization and real-space scale** against independent analytical and trusted diffuse-scattering/3D-ΔPDF references, not only internal NumPy consistency tests. Follow [the reproducible scientific benchmark plan](../docs/SCIENTIFIC_BENCHMARK_REPRODUCTION.md); the current C3/C4 real fixtures are Bragg-present Patterson tests.
- [ ] **P2 — Keep smoke tests reproducible** in supported WSL and report engine version, executable, fixture checksums, logs, output hashes, numerical tolerances and observed outcomes.

### Native runner integration (merged PR #7)

- [x] **P1 — Integrate the new common/NeXus reader with the 3D runner** so a configured output filename is loaded with verified source metadata and invalid geometry is rejected. [Merged PR #7](https://github.com/jmsweng/discus-python/pull/7) includes the reviewed real-engine smoke test (C1 Bragg-subtracted, 65³, `native_nexus=True` → `complete`); this is engineering evidence, **not** independent physical 3D-ΔPDF calibration.

### PR #8 independent scientific audit (in progress, not physically calibrated)

- [x] **Gate 0 groundwork — Correct independent 3×3 reciprocal grid and direct-average oracle.** [Draft PR #8](https://github.com/jmsweng/discus-python/pull/8) adds synthetic skew/2×2×2/origin-shift tests, opt-in frozen-C1 oracle test, and true 33³ archived HDF5 observed-versus-predicted metrics. **Only this reference-code and archived-evidence groundwork is complete**, not Gate 0 scientific acceptance.
- [ ] **Gate 0/1 — Resolve C1 independent ensemble, Waasmaier–Kirfel factors and DISCUS 6.23.01 `set aver` full-grid discrepancy.** Archived empirical IT92 diffuse candidate relative L2 is ≈1.77%; several low-signal voxels differ greatly. Do not mark calibrated.
- [ ] **Gate 2/3 — Correlated-disorder signed off-origin 3D-ΔPDF validation and absolute-scale audit.** Require independent predicted correlations, new real-engine evidence and coordinate/normalization checks; see [PR #8 report](../docs/SCIENTIFIC_VALIDATION_REPORT.md).

### Real 3D DISCUS pipeline (critical path)

- [ ] **P0 — Create and verify a genuine DISCUS_SUITE 6.23.01 macro** that calculates a 3D diffuse-scattering reciprocal-space dataset from a known disordered supercell. Retain input structure, exact macro, engine log and output.
- [x] **P0 — Identify actual native output layout and metadata.** C3/C4 real DISCUS 6.23.01 NeXus fixtures establish the `/entry/data` layout, `data_values`, `data_corner`, and `data_increment_vector` semantics. The [native reader](../src/discuspdf/discus_output.py) validates these values. This checks **file interpretation**, not scientific absolute normalization.
- [ ] **P0 — Calibrate Q coordinates** from lattice/reciprocal basis or trustworthy file metadata; correctly handle r.l.u. ↔ Cartesian-Q, (2\pi), axis orientation, centering and non-orthogonal cells. Reject unsupported geometry.
- [ ] **P0 — Produce physically characterized diffuse-only scattering** using a verified Bragg exclusion/average-structure handling workflow; document whether correction acts on amplitudes or intensities.
- [ ] **P0 — Demonstrate end-to-end `structure → real DISCUS 3D intensity → calibrated Python volume → 3D-ΔPDF`**, record numerical evidence and compare to at least one independently trusted reference.
- [ ] **P1 — Handle masks/windowing/background and missing Bragg voxels** without silently distorting intensity or Fourier normalization.
- [ ] **P1 — Extend real 3D integration and numerical benchmark tests**, including a near-perfect crystal, a simple known disorder model, nonorthogonal lattice geometry and failure cases.
- [ ] **P1 — Document definitive scientific acceptance criteria** (Q/r axis positions, symmetry, intensity/peak amplitudes, normalization, uncertainty and tolerances).

## Javelin-inspired expansion — 21 features, grouped by package

The existing DISCUS engine may implement some of these physical methods already. The unchecked boxes refer to **the new Python-facing workflow and its acceptance criteria**, not the Fortran feature itself. Existing prototypes are noted where applicable.

### Package A — [Structure and interoperability](01_structure_and_interoperability.md)

- [x] **A1 / Feature 1:** Editable, validated `Structure` (atomic species, sites, fractional/cartesian positions, occupancy, cell, PBC and supercell metadata). [`src/discuspdf/structure.py`](../src/discuspdf/structure.py); [PR A-2](https://github.com/jmsweng/discus-python/pull/11) · [`tests/test_structure.py`](../tests/test_structure.py).
  - *Open items: magnetic moments (Package D11), anisotropic ADPs, mixed-PBC distances.*
- [ ] **A2 / Feature 2:** CIF, ASE and optional DiffPy imports/exports with occupancy/PBC fidelity and round-trip tests.
- [ ] **A3 / Feature 3:** Real-format, version-aware DISCUS `.cll` / `.stru` import/export and round-trip tests. *Current wrapper stages existing `.cll` but does not parse or author structures.*
- [ ] **A4 / Feature 4:** Atom selection, supercell generation, substitution and reproducible displacement/disorder editing.
- [x] **A5 / Feature 5:** General direct/reciprocal lattice metrics, validated fractional↔Cartesian conversions, triclinic Q and PBC distance utilities. [`src/discuspdf/lattice.py`](../src/discuspdf/lattice.py); [PR #10](https://github.com/jmsweng/discus-python/pull/10) · [`tests/test_lattice.py`](../tests/test_lattice.py).
  - *Consolidate `DiscusVolume` / `discus_output` onto `Lattice` (separate PR, with Gate 2/3 evidence parity).*

### Package B — [Reciprocal-space grids and data](02_reciprocal_space_and_data.md)

- [ ] **B6 / Feature 6:** Flexible reciprocal-space grid objects for HKL/oblique planes or volumes with validated Q transforms. *Existing FFT supports only centered Cartesian-Q grids.*
- [ ] **B7 / Feature 7:** Labeled `xarray` scattering data, physical coordinates, uncertainties, masks, lattice and provenance; integration with existing `ReciprocalVolume`.
- [ ] **B8 / Feature 8:** Python control of **verified** average-structure/Bragg separation plus benchmarked diffuse output. *Manual preprocessing labels exist, but are not proof of correction.*
- [ ] **B9 / Feature 9:** Python controls and tests for **verified** DISCUS LOTS averaging, including lot size, sampling, seed and normalization.
- [ ] **B18 / Feature 18:** Mantid MDHisto/NeXus reciprocal data import/export with projection matrix, uncertainties, masks and frame verification.
- [ ] **B19 / Feature 19:** Versioned xarray/HDF5/Zarr persistence with validated schema and scientific provenance. *A basic HDF5 loader already exists; this task is broader.*
- [ ] **B20 / Feature 20:** Calibrated Q/3D-ΔPDF visualization, isosurfaces and VTK/ParaView output. *A basic 2D FFT slice plot already exists.*

### Package C — [Disorder and correlations](03_disorder_and_correlations.md)

- [ ] **C12 / Feature 12:** Site/species/distance/shell-resolved periodic neighbor lists on a general lattice.
- [ ] **C13 / Feature 13:** Clearly defined occupational pair probabilities and short-range-order coefficients with reference tests.
- [ ] **C14 / Feature 14:** Reference-structure displacement covariances and normalized, directional correlations.

### Package D — [Monte Carlo and advanced physics](04_monte_carlo_and_advanced_physics.md)

- [ ] **D10 / Feature 10:** Limited independent Python scattering reference calculator for physics verification, not a replacement for DISCUS.
- [ ] **D11 / Feature 11:** Magnetic neutron diffuse-scattering support with validated form factors, moments and Q-dependent projection.
- [ ] **D15 / Feature 15:** Typed Python Monte Carlo/disorder configuration compiled to validated DISCUS MMC macros.
- [ ] **D16 / Feature 16:** Multiple targeted occupational/displacement correlations with convergence and achievement reporting.
- [ ] **D17 / Feature 17:** Site-restricted occupancy/displacement move and energy configuration; verified detailed balance, units and seeded reproducibility.
- [ ] **D21 / Feature 21:** Notebook-friendly atomistic structure and selected-neighbor visualization, with optional exports.

## Maintenance rules

1. Update this file **in the same PR** that delivers a feature; link the commit/PR and test or reference dataset.
2. Use `[x]` only when the **specific checkbox's scope** is demonstrably satisfied; separate any remaining work into a new `[ ]` entry rather than overstating completion.
3. Report software-level tests and real scientific validation as distinct evidence. Passing CI or generating a synthetic FFT **does not establish** valid 3D-ΔPDF physics.
4. Keep package specs stable; use this page to track progress instead of transient references to people, agents, or work currently assigned.
5. Revisit P0/P1 priorities after every validated scientific milestone; dependencies and numerical units take precedence over visualization or optimization.
