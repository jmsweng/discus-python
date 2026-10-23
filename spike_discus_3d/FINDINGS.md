# DISCUS 6.23.01 3D diffuse scattering → 3D-ΔPDF: findings so far (2026-10-22)

Repo: `main` = `ef06490` (PR #4 merge; handoff baseline was `b5f98d3`). Engine: `/usr/local/bin/discus_suite`,
banner "D I S C U S - S U I T E Version 6.23.01", sha256 `d4274f13…45706f3`, WSL2 Ubuntu 22.04, linked to libhdf5.
Nothing below has been committed or pushed.

## Observed (real DISCUS runs; macros in `macros/`, inspector JSON/logs in the result tarballs)
- Command chain that works: `fourier` → `xray`, `wvle moa1`, corners `ll/lr/ul/tl` **in hkl (rlu)**, `na/no/nt`, `abs/ord/top`,
  `temp use`, `disp off`, `lots off`, `set technique:turbo`, `set aver,<percent>`, `run`; then `output` → `value int`,
  `form nexus`, `outf <name>`, `run`. `value 3DPDF` + `form nexus` writes DISCUS's own 3D-PDF/ΔPDF. Runs take ~1 s, ≤65 MB.
- `form hdf5` (documented obsolete) failed here with "I/O Error opening file". Not investigated further.
- File = HDF5 `/entry/data/*` (datasets only, no attributes). Intensity `data_values`, axis order = (abscissa, ordinate, top) = `data_dimension`.
  Metadata present and read by the adapter: `coordinate_unit=basecell_fractional`, `data_type_axes=hkl`, `data_type_reciprocal`, `data_type_content`,
  `data_type_with_bragg` (`bragg_present` | `bragg_subtracted`), `data_corner`, `data_increment_vector`, `unit_cell_lengths/angles`. No Cartesian-Q axes exist.
- `data_increment_vector[:, k]` (h5py read order) = hkl step along array axis k (skewed-grid case C4: column 1 = (0.025, 0.25, 0)).
- Perfect crystal (C0), `set aver,100`: max diffuse 6.4e-21 vs Bragg 1.7e10 → Bragg removed at source. Output is Friedel-merged by DISCUS
  ("Symmetry averaging: Friedel pairs merged"), so the FFT imaginary part is ~1e-16 by construction and says nothing about physics.
- Disorder model `replace GA,AS,all,0.15` (4096 atoms) gave antisite fractions 0.156 (C1) and 0.138 (C2); expected 0.15 ± 0.008.
- `seed 1` is executed but two runs gave different crystals (C1 diffuse arrays differ by about their own maximum). Cause unknown → `seed_test.sh`.

## Independently verified (Python, not DISCUS code)
- Kinematic structure factors computed from the saved atom list (`C1_model.stru`), IT92 form factors (gemmi), B=1.0, ⟨F⟩ per site × box interference
  function, same 65³ hkl grid: vs DISCUS total intensity rel. L2 = 1.2e-3 (corr 0.9999996); vs Bragg-subtracted diffuse 4.2e-3 (corr 0.99998); median ratio 1.0005.
  With DISCUS set to `xray table:inter` (C5, new random crystal, own saved model) the diffuse residual falls to rel. L2 = 1.0e-3 (corr 0.999999, median ratio 1.000000)
  and total to 1.0e-3: the default Waasmaier–Kirfel table explains the 0.05% median bias and ~4x of the diffuse residual. A ~1e-3 relative-L2 difference remains
  UNEXPLAINED (not at machine precision); agreement is therefore ~0.1%, not exact. Axis-order test repeated on C5: stored (h,k,l) rel. 1e-3, any swap ≈ 0.93–0.96.
- Array axis order verified on an asymmetric random crystal: stored (h,k,l) matches; all three axis swaps give rel. L2 ≈ 0.9.
- DISCUS 3DPDF (Bragg-subtracted, C1) equals the plain 65-point FFT of the same diffuse intensity divided by √(65³) to rel. 1.2e-15 (orthonormal scale,
  not ΔhΔkΔl/V). Its header says real-space step 0.125, but the FFT samples sit at 1/(nΔh)=0.12308 (65/64 smaller; 33 pts: 3.1%). Likely a DISCUS metadata
  inconsistency; do not use the 3DPDF header for distances. (Bragg-present 3DPDF is NOT a plain FFT: release notes describe an automatic Bragg/diffuse rescaling.)
- Non-orthogonal cell (C3: 5,6,7 Å, β=105°): the transform of DISCUS's own intensities puts the Patterson peak at u=(0.30,0.20,0.10) (exact on a 0.005 grid) →
  1.907 Å with the cell metric; ignoring β would give 2.045 Å. The metric itself comes from the file's cell parameters (the length check is not independent of the cell definition).

## Assumed / not yet tested
- Handedness/sign of hkl axes (intensity is centrosymmetric, so sign is untestable here). Cartesian orientation of the cell (only lengths are validated).
- Off-origin ΔPDF structure: uncorrelated antisite disorder has none; a correlated model (e.g. `mmc`) is needed for a quantitative peak-ratio benchmark.
- Absolute normalization of ΔPDF (relative units only). Edge/Nyquist/windowing effects at hkl = ±4 (truncation at finite Q) are not quantified.
- `lots`, `nufft`, larger/odd/even grids, other cells, neutron/electron radiation, anomalous dispersion.
- 1D (S4): `set delta` is deprecated per manual; the committed 1D reference was made with qmax=30 while README/test use qmax=25 — unresolved.

## Deliverables staged (untracked, in working tree)
`src/discuspdf/discus_output.py` (reader with source-metadata validation, fractional-coordinate ΔPDF with cell metric, exact evaluator),
`tests/test_discus_output.py` (15 tests, 2 real-file), `tests/fixtures/discus_6.23.01/` (two 0.5 MB real DISCUS files + macros), `spike_discus_3d/*`.
Full suite in the Linux VM: 69 passed, 1 skipped.
