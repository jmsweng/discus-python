# discuspdf — Python tools for disorder, diffuse scattering and PDF with DISCUS

**What this is.** Many crystals are not perfectly ordered: atoms swap places, sit slightly off their ideal
positions, or prefer certain neighbours. This disorder shows up in diffraction as weak **diffuse scattering**
between the sharp Bragg peaks, and in real space as the **pair distribution function** (PDF, the 1D histogram of
interatomic distances) and the **3D-ΔPDF** (a 3D map of how pair correlations differ from the average crystal).

[DISCUS](https://github.com/tproffen/DiffuseCode) (part of DiffuseCode) is an established program that builds
disordered crystal models and calculates their scattering. It is driven by text command files (macros).
**discuspdf** is a Python package around it: you describe structures and calculations in Python, read the exact
DISCUS macro before anything runs, run DISCUS in a separate folder for every job, and load the results back as
labelled data for analysis in notebooks. It does **not** replace DISCUS's own calculations; where it computes
something itself (3D-ΔPDF transforms, correlations, a small scattering check), that code is tested against
DISCUS and against textbook results.

The package is a working prototype. Everything that is described as checked below has a recorded test or a
comparison with real DISCUS 6.23.01 output in this repository; limits are stated where they apply.

## What you can do with it

| Task | Main tools | Docs |
|---|---|---|
| **Describe crystals**: unit cells, coordinates, reciprocal space, periodic distances | `Lattice` | [LATTICE](docs/LATTICE.md) |
| **Build and edit atomic structures**: average (partial occupancy) or realised models, select/edit atoms, supercells, substitution, displacements | `Structure`, `structure_ops` | [STRUCTURE](docs/STRUCTURE.md), [STRUCTURE_OPS](docs/STRUCTURE_OPS.md) |
| **Exchange structures** with DISCUS (`.cll`, `.stru`), CIF, ASE and DiffPy | `read_stru`/`write_stru`, `read_cif`, `to_ase`, `to_diffpy` | [DISCUS_IO](docs/DISCUS_IO.md), [INTEROP](docs/INTEROP.md) |
| **Calculate a powder PDF** G(r) with DISCUS | `PDFConfig`, `PDFRunner` | [VERIFIED_DISCUS_MACRO](docs/VERIFIED_DISCUS_MACRO.md) |
| **Calculate single-crystal diffuse scattering** with DISCUS, with Bragg-peak subtraction (`set aver`) and LOTS averaging | `ReciprocalGrid`, `FourierConfig`, `Diffuse3DRunner` | [RECIPROCAL_GRID](docs/RECIPROCAL_GRID.md), [FOURIER_CONTROLS](docs/FOURIER_CONTROLS.md) |
| **Work with scattering data**: labelled volumes with masks and units, xarray/HDF5, Mantid MDHisto files | `ScatteringVolume`, `read_mantid_mdhisto` | [SCATTERING_DATA](docs/SCATTERING_DATA.md), [MANTID_IO](docs/MANTID_IO.md) |
| **Compute and view a 3D-ΔPDF**, slices in true geometry, isosurfaces, ParaView (VTK) files | `delta_pdf_fractional`, `SampledVolume`, `plot_slice`, `write_vtk` | [3D_INTEGRATION](docs/3D_INTEGRATION.md), [VISUALIZATION](docs/VISUALIZATION.md) |
| **Measure local order**: neighbour shells, Warren–Cowley short-range order, displacement correlations | `neighbor_list`, `correlations` | [NEIGHBORS](docs/NEIGHBORS.md), [CORRELATIONS](docs/CORRELATIONS.md) |
| **Make disordered models with DISCUS Monte Carlo**: target correlations, energies, swap/shift moves, with a report of target vs achieved | `discuspdf.mc` | [MONTE_CARLO](docs/MONTE_CARLO.md), [MC_FACTS](docs/MC_FACTS.md) |
| **Check scattering independently** with a small direct-sum calculator | `discuspdf.reference_scattering` | [REFERENCE_SCATTERING](docs/REFERENCE_SCATTERING.md) |
| **Draw structures**: atoms, displacement arrows, neighbour bonds; export for VESTA/OVITO | `discuspdf.structure_viz` | [STRUCTURE_VIZ](docs/STRUCTURE_VIZ.md) |
| **Magnetic neutron scattering** (experimental; DISCUS has none) | `discuspdf.magnetic_reference` | [MAGNETIC_SCATTERING](docs/MAGNETIC_SCATTERING.md) |

## Notebooks

Each notebook is a guided walkthrough with explanations; most run on data stored in this repository.

| Notebook | What it shows | Needs DISCUS? |
|---|---|---|
| [`01_powder_pdf_with_discus`](notebooks/01_powder_pdf_with_discus.ipynb) | Set up a powder PDF G(r) calculation, read the generated macro, run it and compare Qmax values | yes (to run) |
| [`02_3d_delta_pdf_synthetic_example`](notebooks/02_3d_delta_pdf_synthetic_example.ipynb) | Fourier transform of a made-up diffuse-scattering volume into a 3D-ΔPDF | no |
| [`03_lattices_and_coordinates`](notebooks/03_lattices_and_coordinates.ipynb) | Unit cells, fractional ↔ Cartesian coordinates, Q and d-spacings, periodic distances | no |
| [`04_building_and_editing_structures`](notebooks/04_building_and_editing_structures.ipynb) | Build an atomic structure, select and edit atoms, wrap coordinates, validation | no |
| [`05_reading_and_writing_discus_files`](notebooks/05_reading_and_writing_discus_files.ipynb) | DISCUS `.cll` (unit cell) and `.stru` (crystal) files | no (optional live check) |
| [`06_cif_ase_and_diffpy_conversion`](notebooks/06_cif_ase_and_diffpy_conversion.ipynb) | Read a CIF, convert to and from ASE and DiffPy, write `.cll` and CIF | no |
| [`07_diffuse_scattering_data_and_3d_delta_pdf`](notebooks/07_diffuse_scattering_data_and_3d_delta_pdf.ipynb) | Reciprocal-space grids, DISCUS output, Bragg subtraction or masking, 3D-ΔPDF, xarray/HDF5, VTK and Mantid export | no (uses stored DISCUS data) |
| [`08_disordered_supercells`](notebooks/08_disordered_supercells.ipynb) | From an average structure to a disordered supercell: sampling, substitution, displacements | no |
| [`09_neighbors_and_short_range_order`](notebooks/09_neighbors_and_short_range_order.ipynb) | Neighbour shells and how ordered a crystal is: occupational and displacement correlations, compared with DISCUS | no |
| [`10_monte_carlo_scattering_checks_and_magnetism`](notebooks/10_monte_carlo_scattering_checks_and_magnetism.ipynb) | Independent scattering check against DISCUS, structure pictures, analysing a DISCUS Monte Carlo run, magnetic scattering basics | no (uses stored DISCUS data) |

Notebooks set the working directory to the repository root when started from `notebooks/`; if you start them
elsewhere, edit the path near the top.

## Installation

1. **DISCUS** (only needed to run new calculations): install DISCUS_SUITE from the
   [DiffuseCode releases](https://github.com/tproffen/DiffuseCode/releases) and check that
   `discus_suite -noautorun -macro <file>` works on your system. The checks in this repository used
   DISCUS_SUITE 6.23.01 (on Windows through WSL Ubuntu 22.04).
2. **Python 3.10 or newer**:

```bash
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -e '.[notebook,dev]'
jupyter lab
```

Optional extras (`pip install -e '.[name]'`): `data` (xarray, HDF5), `viz` (isosurfaces), `cif` (CIF reading),
`ase`, `diffpy`, `interop` (all three converters), `reference` (X-ray form factors and neutron lengths for the
scattering checks). `import discuspdf` works without them; a feature that needs one says which.

## A first calculation: powder PDF

```python
from discuspdf import PDFConfig, PDFRunner

cfg = PDFConfig(structure="examples/gaas.cll", supercell=(5, 5, 5), rmax=12.0, dr=0.02,
                qmax=7.5, qsig=0.0, radiation="xray")
runner = PDFRunner(runs_directory="runs")
prepared = runner.prepare(cfg)                     # writes the job folder; nothing runs yet
print((prepared.directory / "run.mac").read_text())  # read the exact DISCUS macro

result = runner.execute(prepared)                  # runs exactly that macro
r, gr = result.load()
result.plot()
```

Every job folder keeps a copy of the input with its SHA-256 hash, the macro (`run.mac`), the settings
(`config.json`), the DISCUS log and a status file, so a result can always be traced back. The runner checks the
expected r-grid and DISCUS error messages; still look at the log for new kinds of calculations. To run DISCUS
some other way (for example through WSL on Windows), pass a different `launcher`.

**Qmax limit.** DISCUS 6.23.01 cannot use Q above about 12.3/λ — about 8.0 Å⁻¹ at the default X-ray wavelength
of 1.5406 Å. A larger `qmax` still runs, but gives the same G(r) as Qmax ≈ 8; `PDFConfig` warns
(`QmaxLimitWarning`) when `qmax × wavelength > 12`. Details: [VERIFIED_DISCUS_MACRO](docs/VERIFIED_DISCUS_MACRO.md).

## What has been checked against DISCUS

The comparisons are pre-registered (expected results and pass/fail limits committed before DISCUS was run) and
their inputs, macros, logs and analysis scripts are in [`benchmarks/discus_6_23_01/`](benchmarks/discus_6_23_01/).

* **3D-ΔPDF from DISCUS diffuse scattering** — for commensurate grids, one atom per cell and B = 0, DISCUS's
  Bragg-subtracted intensity (`set aver,100`) transformed with `delta_pdf_fractional` reproduces an independently
  predicted 3D-ΔPDF, in relative and in absolute units (e² Å⁻³ for the whole crystal).
  [Gate 2](benchmarks/discus_6_23_01/gate2_prereg/RESULT_v4.md), [Gate 3](benchmarks/discus_6_23_01/gate3/RESULT.md),
  [validation report](docs/SCIENTIFIC_VALIDATION_REPORT.md).
* **Single-crystal intensities** — the independent scattering calculator reproduces DISCUS X-ray and neutron
  intensities, with and without Bragg subtraction, to ~10⁻¹⁵ once DISCUS's form-factor table rounding is
  reproduced. [Result](benchmarks/discus_6_23_01/d_reference/RESULT.md).
* **Fourier options** — macros written by `FourierConfig` reproduce DISCUS probe runs byte for byte
  (`set aver`, LOTS, seeds). [Probe result](benchmarks/discus_6_23_01/b_probes/RESULT.md).
* **Short-range order** — DISCUS `chem corr` values are reproduced (29/29 to 5 × 10⁻¹³).
  [Probe result](benchmarks/discus_6_23_01/c_probes/RESULT.md).
* **Monte Carlo** — the generated MMC macros reproduce hand-written DISCUS runs exactly, and the achieved
  correlations are confirmed by an independent measurement. [Probe result](benchmarks/discus_6_23_01/d_mc_probes/RESULT.md).

Things found along the way that matter when using DISCUS 6.23.01:

* Its native `value 3DPDF` output is a unitary FFT (divided by √N, no physical units) and its written axis step is
  wrong by n/(n−1), with the origin one voxel off on even grids. Use `value int` + `delta_pdf_fractional` instead.
* Its Monte Carlo keeps uphill and zero-energy moves with probability `1/(1+e^{ΔE/kT})` (not the Metropolis rule),
  and multi-threaded runs are not reproducible; `discuspdf.mc` always runs single-threaded with a full seed.
* It stores magnetic moments but calculates no magnetic scattering.

**Not covered yet:** 3D-ΔPDF validation beyond the scope above (several atoms per cell, B ≠ 0, non-commensurate
grids), automatic background and absorption handling, molecules, and a production magnetic workflow. The
[task tracker](todo/TASKS.md) lists open items.

## Working with 3D data: points to know

* A grid in h, k, l (reciprocal-lattice units) is not a Cartesian Q grid. For DISCUS NeXus output
  (`form nexus`), use `read_discus_volume()` / `Diffuse3DRun.load_discus_volume()` and `delta_pdf_fractional()`,
  which keep the hkl axes and the cell metric. For other HDF5 files with checked, centred Cartesian-Q axes, use
  `read_hdf5_volume(...)` and a `ReciprocalVolume` with explicit coordinate system, units and Bragg treatment.
* The Cartesian-Q FFT returns an `FFTResult` (complex volume, real-space axes, diagnostics). If the imaginary part
  is larger than `imaginary_tolerance` (default 10⁻⁴ of the real part) it raises `FFTValidationError`, which still
  carries the computed result; pass `raise_on_imaginary=False` to get a warning instead.
* Notebook `02` uses a synthetic volume only. Quantitative work needs masking/windowing, correct intensity
  normalisation and a comparison with a trusted reference.

## Repository layout

| Path | Contents |
|---|---|
| `src/discuspdf/` | The Python package. |
| `notebooks/` | The walkthrough notebooks listed above. |
| `examples/` | Small scripts, the GaAs unit cell and a recorded real 1D PDF run. |
| `docs/` | User guides for each feature, DISCUS facts found in the source and in test runs, and the validation report. |
| `benchmarks/discus_6_23_01/` | The pre-registered comparisons with DISCUS 6.23.01: inputs, predictions, analysis scripts, results and small evidence files. Large engine volumes are not committed; `regenerate_evidence.py` recreates them. |
| `tests/` | Unit tests (`python -m pytest -q`, run by GitHub Actions on every push). Tests that need a real DISCUS installation are opt-in (`DISCUS_LIVE=1`). |
| `spike_discus_3d/` | The original 3D probe kit (verified macros used by the 3D runner example and its tests) and early findings. |
| `todo/` | Roadmap and [task tracker](todo/TASKS.md). |

## Notes

* DiffuseCode contains an early, explicitly not release-ready Python (F2PY) bridge in `python/suite_python.f90`;
  this package does not use it and talks to DISCUS only through macros.
* DISCUS is intended as a public-domain program; commercial use needs the authors' permission. Check DiffuseCode's
  licence for your use.
