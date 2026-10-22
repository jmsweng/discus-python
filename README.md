# discuspdf starter — PDF and 3D-ΔPDF notebooks

A **prototype** Python package around the existing [DiffuseCode](https://github.com/tproffen/DiffuseCode) executable. It does **not** reimplement atomic pair-distribution calculations.

## Scope

- Implemented: prepare reproducible per-run directories, compile a basic 1D atomic PDF macro from a `.cll` unit-cell file, launch DISCUS, validate and plot its ASCII output, sweep parameters.
- Implemented independently: FFT of an **already Bragg-subtracted 3D diffuse intensity volume** on a uniform *Cartesian* Q grid in Å⁻¹, and a 2D slice viewer.
- **Not implemented yet:** automatic generation of the 3D diffuse intensity from a structure, Bragg removal, physical normalization, direct `.cif` import, advanced atom-pair selection, or 3D volume exports from native DISCUS formats.

## Installation

1. Install the binary for DISCUS_SUITE separately from the [DiffuseCode releases](https://github.com/tproffen/DiffuseCode/releases). Confirm `discus_suite` works interactively and check its `-macro` command-line invocation on your platform.
2. Set up Python in this repository:

```bash
python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate
pip install -e '.[notebook,dev]'
jupyter lab
```

Open `notebooks/01_pdf_quickstart.ipynb`.

## First run

```python
from discuspdf import PDFConfig, PDFRunner

cfg = PDFConfig(
    structure="examples/gaas.cll",
    supercell=(5, 5, 5), rmax=12.0, dr=0.02,
    qmax=25.0, qsig=0.0, radiation="xray",
)
runner = PDFRunner(runs_directory="runs")
prepared = runner.prepare(cfg)  # no execution; inspect generated run.mac
print((prepared.directory / "run.mac").read_text())

# After confirming the macro runs in your DISCUS installation:
result = runner.execute(prepared)  # runs the exact inspected macro
r, gr = result.load()
result.plot()
```

**Important:** The macro syntax follows [documented DISCUS 1D PDF examples](https://www.xray.cz/ecm-cd/soft/xray/general/discus/discus/disman/node57.html), which are older than the latest release; validate with your local version. `result.pdf` is a **plain-text numeric PDF G(r) file**, *not* an Adobe document. This prototype uses the conventional program name/flags `discus_suite -noautorun -macro run.mac`; override the runner `launcher` if needed.

## Notebook working-directory convention

The example notebook sets the working directory to the project root if launched from `notebooks/`. For other launch locations, edit the path near the top.

## Safety/reproducibility

Each run directory records the copied `.cll` input, a SHA-256 hash, the generated `run.mac`, parameters (`config.json`), log (`run.log`) and status. Do not delete these when comparing calculations. The wrapper also verifies the expected r-grid and checks basic engine errors; always inspect DISCUS logs and benchmark one result manually.

## 3D-ΔPDF analysis

See `notebooks/02_3d_delta_pdf.ipynb` for a **synthetic demonstration only**. For real 3D data, start with a measured or simulated volume after Bragg/background removal, on a centered, uniform Cartesian Q grid in Å⁻¹. Inspect HDF5 datasets and use `read_hdf5_volume(...)` with *explicit* dataset paths, or load NumPy arrays yourself. A grid in h,k,l reciprocal-lattice units is not interchangeable with Cartesian Q for real-space axis calibration.

For quantitative work, add mask/window handling, correct reciprocal-space intensity normalization, and compare one FFT result to a trusted 3D-ΔPDF workflow. No 3D-ΔPDF calculations are run automatically inside DISCUS by this starter.

## Next stages

1. Run the unit tests and validate the 1D PDF output against your local DISCUS version.
2. Add your frequently changed structure/disorder command blocks and partial PDFs.
3. Write a 3D reciprocal-volume output adapter for an existing validated DISCUS macro.
4. Standardize Q-grid metadata and mask/Bragg subtraction for 3D-ΔPDF analysis.
5. Add optimized parameter sweeps and experimental comparison metrics.

## Notes

The DiffuseCode repository contains an initial, explicitly *not release-ready* F2PY bridge in `python/suite_python.f90`. This MVP does not depend on it. Inspect DiffuseCode's license and check permission for intended commercial use.
## Explicit powder settings

`PDFConfig` exposes `qmin`, `qmax`, `dq`, `qsig` (DISCUS `set delta`), and optional `wavelength` (Å), in addition to `rmax` and `dr`. They are saved with every run. `qmax=0` is rejected: omitting the command does not establish that the termination correction is disabled. The real-space range must be divisible by `dr`.

Run `python -m pytest -q` for unit tests; see the separately documented manual DISCUS integration check. The local installer previously added `/etc/sudoers.d/90-discus` for passwordless sudo; review and remove it with administrator approval if no longer required.
