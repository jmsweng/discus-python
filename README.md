# discuspdf — a Python layer for DISCUS: disorder models, diffuse scattering, 3D-ΔPDF and PDF

# See notebooks for usage examples, documentation is currently a work in progress.

`discuspdf` drives [DISCUS](https://github.com/tproffen/DiffuseCode) (DISCUS_SUITE 6.23.01) from Python. It
generates DISCUS macros from typed configurations (inspectable before execution), runs each job in an isolated,
hashed run directory, and loads native output (NeXus, `.stru`, logs) back as labelled data with explicit metric,
units and Bragg treatment. Around that it provides structure handling and I/O, reciprocal-space grids and data
containers, the 3D-ΔPDF transform on hkl grids, short-range-order analysis, a Monte Carlo (MMC) front end and a
small independent structure-factor reference.

DISCUS remains the scattering and modelling engine. Engine behaviour that the package relies on was taken from
the 6.23.01 source and confirmed in pre-registered comparisons with the engine (see
[Validation against DISCUS 6.23.01](#validation-against-discus-62301)); scope limits are stated where they apply.
The package is a research prototype.

## Capabilities

| Area | Main API | Docs |
|---|---|---|
| **Lattice**: direct/reciprocal metric, fractional ↔ Cartesian, Q = 2π·hkl·B, d-spacings, exact triclinic minimum image | `Lattice` | [LATTICE](docs/LATTICE.md) |
| **Structures**: average (site-occupancy) and configuration models; selection, editing, supercells, substitution, seeded displacements and configuration sampling | `Structure`, `structure_ops` | [STRUCTURE](docs/STRUCTURE.md), [STRUCTURE_OPS](docs/STRUCTURE_OPS.md) |
| **I/O**: DISCUS `.cll`/`.stru` (version-aware, B-value policy of 6.23.01), CIF, ASE, DiffPy | `read_stru`/`write_stru`, `read_cif`, `to_ase`, `to_diffpy` | [DISCUS_IO](docs/DISCUS_IO.md), [INTEROP](docs/INTEROP.md) |
| **Powder PDF** G(r) via DISCUS (Qmin/Qmax/dQ, `qsig`, wavelength; the 6.23.01 Qmax limit is flagged) | `PDFConfig`, `PDFRunner` | [VERIFIED_DISCUS_MACRO](docs/VERIFIED_DISCUS_MACRO.md) |
| **Single-crystal Fourier** macros: oblique planes/lines/volumes, `set aver`, LOTS, 8-value seeds, X-ray or neutron, wavelength, NeXus output | `ReciprocalGrid`, `FourierConfig`, `Diffuse3DRunner` | [RECIPROCAL_GRID](docs/RECIPROCAL_GRID.md), [FOURIER_CONTROLS](docs/FOURIER_CONTROLS.md) |
| **Scattering data**: `ScatteringVolume` with variance, masks, units, Bragg label and provenance; xarray, versioned HDF5, Mantid MDHisto (HKL frame, `W_MATRIX`) | `ScatteringVolume`, `read_mantid_mdhisto` | [SCATTERING_DATA](docs/SCATTERING_DATA.md), [MANTID_IO](docs/MANTID_IO.md) |
| **3D-ΔPDF**: transform on fractional hkl grids with absolute scaling; slices in true geometry, isosurfaces, `.vti`/`.vts` for ParaView | `delta_pdf_fractional`, `SampledVolume`, `plot_slice`, `write_vtk` | [3D_INTEGRATION](docs/3D_INTEGRATION.md), [VISUALIZATION](docs/VISUALIZATION.md) |
| **Short-range order**: periodic neighbour lists (`set vec`-style vectors, shells), Warren–Cowley α and correlation coefficients with null-model errors, displacement covariances/correlations (Pearson or DISCUS convention) | `neighbor_list`, `correlations` | [NEIGHBORS](docs/NEIGHBORS.md), [CORRELATIONS](docs/CORRELATIONS.md) |
| **MMC**: occupational/displacement (CORR/ENER), spring and Lennard-Jones targets; swap, swdisp and shift moves; feedback history, acceptance, target vs achieved, independent re-measurement, replica spread | `discuspdf.mc` | [MONTE_CARLO](docs/MONTE_CARLO.md), [MC_FACTS](docs/MC_FACTS.md) |
| **Reference F(Q)**: direct-sum amplitudes, \|F\|² and \|F − ⟨F⟩\|², stated X-ray (Waasmaier–Kirfel) or neutron factors, DW, occupancy | `discuspdf.reference_scattering` | [REFERENCE_SCATTERING](docs/REFERENCE_SCATTERING.md) |
| **Structure views**: atoms by species/site, displacement fields, neighbour bonds; export (extxyz, POSCAR, xsf, CIF) with read-back check | `discuspdf.structure_viz` | [STRUCTURE_VIZ](docs/STRUCTURE_VIZ.md) |
| **Disorder builder back end**: model JSON from the discus-disorder-builder page → configuration with VOID padding, ordering MC for molecules/vacancies, DISCUS 3D/1D/mmc runs, selectable output formats; DISCUS atom-name table and guessing rules; extended XYZ without ASE | `discuspdf.builder`, `discuspdf.atom_names`, `write_extxyz` | [DISORDER_BUILDER](docs/DISORDER_BUILDER.md) |
| **Magnetic scattering** (experimental): unpolarised dipole-approximation reference; DISCUS 6.23.01 computes none | `discuspdf.magnetic_reference` | [MAGNETIC_SCATTERING](docs/MAGNETIC_SCATTERING.md) |

## Notebooks

Worked examples; except where noted they run on DISCUS output committed to this repository.

| Notebook | What it shows | Needs DISCUS? |
|---|---|---|
| [`01_powder_pdf_with_discus`](notebooks/01_powder_pdf_with_discus.ipynb) | Powder PDF macro generation, execution and a Qmax sweep (incl. the 6.23.01 limit) | yes (to run) |
| [`02_3d_delta_pdf_synthetic_example`](notebooks/02_3d_delta_pdf_synthetic_example.ipynb) | 3D-ΔPDF transform of a synthetic Cartesian-Q volume; diagnostics; X-ray/neutron radiation and wavelength in DISCUS 3D macros (recorded runs) | no |
| [`03_lattices_and_coordinates`](notebooks/03_lattices_and_coordinates.ipynb) | Lattice metric, coordinate transforms, Q and d-spacings, minimum image in skewed cells | no |
| [`04_building_and_editing_structures`](notebooks/04_building_and_editing_structures.ipynb) | Average vs configuration models, selection, editing, wrapping, validation | no |
| [`05_reading_and_writing_discus_files`](notebooks/05_reading_and_writing_discus_files.ipynb) | `.cll` vs `.stru`, unit-cell ↔ supercell fractions, B-value handling | no (optional live check) |
| [`06_cif_ase_and_diffpy_conversion`](notebooks/06_cif_ase_and_diffpy_conversion.ipynb) | COD GaAs CIF → supercell → ASE/DiffPy, with loss reporting | no |
| [`07_diffuse_scattering_data_and_3d_delta_pdf`](notebooks/07_diffuse_scattering_data_and_3d_delta_pdf.ipynb) | hkl grids, DISCUS NeXus loading, Fourier macros, Bragg subtraction vs masking, 3D-ΔPDF, HDF5/xarray, VTK, Mantid | no |
| [`08_disordered_supercells`](notebooks/08_disordered_supercells.ipynb) | Sampling a configuration from site occupancies, substitution, seeded displacements, composition checks | no |
| [`09_neighbors_and_short_range_order`](notebooks/09_neighbors_and_short_range_order.ipynb) | Neighbour shells, occupational and displacement SRO, comparison with `chem corr` | no |
| [`10_monte_carlo_scattering_checks_and_magnetism`](notebooks/10_monte_carlo_scattering_checks_and_magnetism.ipynb) | Reference F(Q) vs DISCUS, structure views, analysis of an MMC run, magnetic reference | no |

Started from `notebooks/`, they switch to the repository root; elsewhere, adjust the path in the first cell.

## Installation

1. DISCUS_SUITE from the [DiffuseCode releases](https://github.com/tproffen/DiffuseCode/releases) (needed only to
   run new calculations); the runners call `discus_suite -noautorun -macro <file>` (configurable `launcher`, e.g.
   through WSL on Windows). Validation used 6.23.01 under WSL Ubuntu 22.04.
2. Python ≥ 3.10:

```bash
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -e '.[notebook,dev]'
jupyter lab
```

Optional extras: `data` (xarray, h5py), `viz` (scikit-image), `cif` (gemmi), `ase`, `diffpy`, `interop` (all three),
`reference` (xraydb, periodictable). The core import needs none of them.

## Example: powder PDF

```python
from discuspdf import PDFConfig, PDFRunner

cfg = PDFConfig(structure="examples/gaas.cll", supercell=(5, 5, 5), rmax=12.0, dr=0.02,
                qmax=7.5, qsig=0.0, radiation="xray")
runner = PDFRunner(runs_directory="runs")
prepared = runner.prepare(cfg)                       # stage only
print((prepared.directory / "run.mac").read_text())  # inspect the macro

result = runner.execute(prepared)                    # run exactly that macro
r, gr = result.load()   # reads <run>/result.csv (r (Angstrom),G(r)); DISCUS's raw text is <run>/result.pdf
result.plot()
```

Run directories hold the staged inputs with SHA-256 hashes, `run.mac`, `config.json`, the engine log and a status
record; the runners check engine error lines and the expected output grid before reporting success.

**Qmax.** DISCUS 6.23.01 silently caps the powder Qmax at ≈ 12.3/λ (≈ 8.0 Å⁻¹ for 1.5406 Å; `powder.f90` clamps
`pow_qmax` before its own check). `PDFConfig` raises `QmaxLimitWarning` when `qmax·λ > 12`. See
[VERIFIED_DISCUS_MACRO](docs/VERIFIED_DISCUS_MACRO.md).

## Validation against DISCUS 6.23.01

Predictions and pass/fail criteria were committed before each engine run; inputs, macros, logs, evidence and
analysis scripts are in [`benchmarks/discus_6_23_01/`](benchmarks/discus_6_23_01/).

* **3D-ΔPDF** — commensurate grids, one site per cell, B = 0: `set aver,100` intensities transformed with
  `delta_pdf_fractional` reproduce an independently predicted signed off-origin 3D-ΔPDF (ρ/ρ₀ to 4.5 × 10⁻⁶) and its
  absolute scale (e² Å⁻³, whole crystal).
  [Gate 2](benchmarks/discus_6_23_01/gate2_prereg/RESULT_v4.md), [Gate 3](benchmarks/discus_6_23_01/gate3/RESULT.md),
  [validation report](docs/SCIENTIFIC_VALIDATION_REPORT.md).
* **Single-crystal intensities** — X-ray and neutron, total and `set aver,100`, B ≠ 0: the reference F(Q) agrees
  to 2.5 × 10⁻¹⁵–1.8 × 10⁻¹³ (rel. L2) on commensurate grids and 9 × 10⁻⁶ for a displaced monoclinic crystal on a
  non-commensurate grid, with DISCUS's s-rounding of f and DW emulated. [Result](benchmarks/discus_6_23_01/d_reference/RESULT.md).
* **Fourier controls** — generated macros reproduce the probe runs byte for byte (`set aver`, LOTS, seeds,
  planes). [Probe result](benchmarks/discus_6_23_01/b_probes/RESULT.md).
* **SRO** — `chem corr occ`/`disp` reproduced (29/29 values to 5 × 10⁻¹³, identical pair counts).
  [Probe result](benchmarks/discus_6_23_01/c_probes/RESULT.md).
* **MMC** — generated macros reproduce hand-written runs exactly; printed achieved values are reproduced by an
  independent re-measurement with DISCUS's definitions. [Probe result](benchmarks/discus_6_23_01/d_mc_probes/RESULT.md).

DISCUS 6.23.01 behaviour worth knowing:

* Native `value 3DPDF` is a unitary FFT (÷√N, no physical units); its written axis step is off by n/(n−1) and the
  origin by one voxel on even grids. Use `value int` + `delta_pdf_fractional`.
* f and the isotropic Debye–Waller factor are tabulated at s rounded to 0.001 Å⁻¹ (≈10⁻³ relative); B = 0 in a
  structure file is read as 1.0 Å².
* MMC acceptance is heat-bath for uphill and ΔE = 0 moves, `1/(1+e^{ΔE/kT})` with ΔE rounded to 10⁻⁴ — not
  Metropolis; OpenMP runs are not reproducible. `discuspdf.mc` always writes an 8-value seed and runs serially.
  More in [MC_FACTS](docs/MC_FACTS.md).
* Magnetic moments are stored per atom but used by no scattering routine, and `.stru` I/O drops them.

**Open:** 3D-ΔPDF validation for multi-site cells, B ≠ 0 and non-commensurate grids; background/absorption
corrections; molecules; magnetic scattering beyond the reference. See the [task tracker](todo/TASKS.md).

## Notes on 3D data

* hkl grids are never relabelled as Cartesian Q. DISCUS NeXus output (`form nexus`) is read with
  `read_discus_volume()` / `Diffuse3DRun.load_discus_volume()` and transformed with `delta_pdf_fractional()`, keeping
  the hkl axes and cell metric. Other HDF5 volumes need verified centred Cartesian-Q axes (`read_hdf5_volume`,
  `ReciprocalVolume` with explicit frame, units and Bragg treatment).
* The Cartesian-Q FFT returns an `FFTResult` (complex volume, r axes, diagnostics) and raises `FFTValidationError`
  (carrying the result) when ‖Im‖/‖Re‖ exceeds `imaginary_tolerance` (default 10⁻⁴); `raise_on_imaginary=False`
  downgrades this to a warning.
* Masking/windowing of experimental data and its effect on the transform are not automated.

## Repository layout

| Path | Contents |
|---|---|
| `src/discuspdf/` | The Python package. |
| `notebooks/` | The notebooks listed above. |
| `examples/` | Example scripts, the GaAs unit cell and a recorded 1D PDF run. |
| `docs/` | Per-feature documentation, DISCUS 6.23.01 source/probe facts, validation report. |
| `benchmarks/discus_6_23_01/` | Pre-registered engine comparisons: inputs, predictions, analysis, results and evidence (large volumes are regenerated by `regenerate_evidence.py`). |
| `tests/` | `python -m pytest -q` (CI on every push); engine tests are opt-in with `DISCUS_LIVE=1`. |
| `spike_discus_3d/` | Original 3D probe kit (verified macros used by the 3D runner example and tests). |
| `todo/` | Roadmap and [task tracker](todo/TASKS.md). |

## Notes

* The package talks to DISCUS only through macros; it does not use DiffuseCode's experimental F2PY bridge
  (`python/suite_python.f90`).
* DISCUS licence: free for non-commercial use; commercial use requires the authors' permission.
