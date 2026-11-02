# Package D — Monte Carlo models, physics options and structure visualization

**Status:** Done (Features 10, 15–17, 21; Feature 11 investigated, experimental reference only) — see [TASKS.md](TASKS.md) · **Priority:** After A–C and real 3D validation · **Suggested PR family:** `feat/models-*`

**Reference in Javelin:** [mc.py](https://github.com/rosswhitfield/javelin/blob/master/javelin/mc.py), [modifier.pyx](https://github.com/rosswhitfield/javelin/blob/master/javelin/modifier.pyx), [energies.pyx](https://github.com/rosswhitfield/javelin/blob/master/javelin/energies.pyx), [fourier.py](https://github.com/rosswhitfield/javelin/blob/master/javelin/fourier.py), [Ising tutorial](https://github.com/rosswhitfield/javelin/blob/master/docs/tutorials/ising.rst), [size-effect tutorial](https://github.com/rosswhitfield/javelin/blob/master/docs/tutorials/size_effect.rst).

## Objective

Let users define and compare disordered-structure generation workflows through Python and notebooks while using validated DISCUS Monte Carlo and diffuse-scattering backends. Add advanced physics only with explicit, independently verified conventions.

## Planned features

### 10. Optional independent Python scattering reference
- [x] Create a **small, clearly limited** structure-factor/intensity calculator for analytical fixtures and cross-checking DISCUS (not a competing full scattering engine).
- [x] Support a stated choice of neutron nuclear scattering lengths or X-ray form factors, explicit Q convention, atomic positions and sample normalization.
- [x] Benchmark a single atom, two-atom interference, perfect crystal and simple correlated disorder; document performance scaling and numerical tolerances.
- [ ] Optimize with NumPy/Numba/Cython only if reference correctness is established. *(Not needed: ~30 ns per atom per Q point with NumPy.)*

### 11. Magnetic neutron scattering (later)
- [x] Investigate what DISCUS and Javelin actually support for magnetic moments, form factors and magnetic diffuse scattering.
- [x] Specify moments, orientation relative to Q, transverse projection, scattering units and Q=0 handling.
- [x] Build tiny analytical examples and compare independently before exposing as a scientific production feature.

### 15. Python-configured Monte Carlo simulations
- [x] Define typed `MCConfig`, site selections, cycles, temperature schedules, random seeds, moves and constraints. *(One temperature per run; no schedules yet.)*
- [x] Compile to **verified** DISCUS MMC commands; save executable version, generated macro and physical parameters.
- [x] Support inspection/preparation before execution and per-run directories consistent with current runners.

### 16. Target-correlation model generation
- [x] Expose multiple neighbor-shell occupational and displacement targets via Package C's documented correlation definitions.
- [x] Produce convergence diagnostics: target vs achieved values, energies, acceptance history and seed-dependent uncertainty.
- [x] Allow failure to converge or unreachable targets to be reported explicitly rather than suggesting success.

### 17. Energy and structural move libraries
- [x] Explore occupancy swaps, displacement swaps/shifts, random Gaussian/range moves, Ising-like interactions, harmonic springs and Lennard-Jones terms.
- [x] Determine which are already present in DISCUS/MMC and **wrap verified existing behavior first**.
- [x] Document units, energy sign convention, detailed balance, periodic boundaries, site restrictions and acceptance probabilities.
- [x] Test constrained compositions and zero-energy acceptance; verify reproducibility at fixed random seed.

### 21. Atomic structure visualization
- [x] Add notebook-friendly plotting using ASE or other maintained tools for supercell snapshots, site coloring, displacement vectors and selected neighbor shells.
- [x] Allow export to established viewers without altering cell orientation or coordinates.
- [x] Make this an optional package extra; visualization must not become a dependency of headless scientific runs. *(matplotlib was already a base dependency; it and ASE are imported lazily, tested.)*

## Suggested PRs

1. Minimal tested Python configuration/compiler for already verified DISCUS MMC features.
2. MC result parsing, log and convergence plots; fixed-seed reproducibility and correlation targets.
3. Optional small independent scattering reference implementation and numerical benchmarks.
4. Optional advanced magnetic scattering and structure viewers.

## Acceptance criteria

- [x] A verified DISCUS MC macro can be generated, inspected, executed and its results loaded via a documented Python API.
- [x] Identical seeds and inputs give reproducible results within documented platform/version limits.
- [x] Achieved SRO/displacement correlations are measured independently (Package C) and compare to requested targets.
- [x] Model changes and conserved composition are validated; failure/nonconvergence is surfaced.
- [x] Magnetic and direct Python scattering (if implemented) meet independent analytical benchmarks.
- [x] Existing scientific 1D PDF/3D ΔPDF workflows remain compatible.

## Dependencies / cautions

Depends on Package A structures, Package C correlations and Claude's scientifically verified DISCUS 3D workflow. DISCUS MMC may already implement many moves/energies; do not duplicate Fortran methods by default. Javelin's aging Cython extension and public API are inspiration only; do not make them a runtime requirement. Preserve MIT notices if source is copied.
