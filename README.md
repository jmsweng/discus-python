# discuspdf starter — PDF and 3D-ΔPDF notebooks

A **prototype** Python package around the existing [DiffuseCode](https://github.com/tproffen/DiffuseCode) executable. It does **not** reimplement atomic pair-distribution calculations.

## Scope

- Implemented: prepare reproducible per-run directories, compile a basic 1D atomic PDF macro from a `.cll` unit-cell file, launch DISCUS, validate and plot its ASCII output, sweep parameters.
- Implemented independently: FFT of an **already Bragg-subtracted 3D diffuse intensity volume** on a uniform *Cartesian* Q grid in Å⁻¹, and a 2D slice viewer.
- Implemented for **DISCUS_SUITE 6.23.01 native common/NeXus**: a source-metadata-checked `read_discus_volume` adapter, physical cell metric and a fractional hkl Fourier transform; the 3D runner can load the native file from a configurable output name. See [3D integration instructions](docs/3D_INTEGRATION.md).
- **Not yet complete:** automatic 3D macro generation from arbitrary input structures, independently benchmarked correlated-disorder 3D-ΔPDF and absolute normalization, robust Bragg/background/mask processing, direct `.cif` import, advanced atom-pair selection, or general 3D visualization/export.

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

See `notebooks/02_3d_delta_pdf.ipynb` for a **synthetic demonstration only**. For DISCUS 6.23.01 `form nexus` outputs, use `read_discus_volume()` or `Diffuse3DRun.load_discus_volume()` and `delta_pdf_fractional()`, which preserve hkl axes and the cell metric. For other HDF5 files with independently verified centered uniform Cartesian-Q arrays, use `read_hdf5_volume(...)` or load NumPy arrays. A grid in h,k,l reciprocal-lattice units must never be relabeled as Cartesian Q. See [3D integration](docs/3D_INTEGRATION.md) and the [reproducible-science benchmark plan](docs/SCIENTIFIC_BENCHMARK_REPRODUCTION.md).

For quantitative work, add mask/window handling, correct reciprocal-space intensity normalization, and compare one FFT result to a trusted 3D-ΔPDF workflow. **Breaking change:** the 3D FFT no longer returns a tuple — it returns an `FFTResult` containing `complex_volume`, `r_axes`, and `diagnostics`; `real_volume` is available for plotting. Do not unpack the result as a tuple. If the imaginary/real norm ratio exceeds `imaginary_tolerance` (default `1e-4`), the transform raises `FFTValidationError` (a `ValueError` subclass) that carries the computed `FFTResult` on `exc.result`; pass `raise_on_imaginary=False` to get a warning and the same computed values instead. No 3D-ΔPDF calculations are run automatically inside DISCUS by this starter; integration history is recorded in [docs/INTEGRATION_HISTORY.md](docs/INTEGRATION_HISTORY.md).

## Next stages

1. Run the unit tests and validate the 1D PDF output against your local DISCUS version.
2. Add your frequently changed structure/disorder command blocks and partial PDFs.
3. Integrate a validated real-DISCUS 3D intensity macro with the native adapter and record full model-level benchmarking.
4. Standardize Q-grid metadata and mask/Bragg subtraction for 3D-ΔPDF analysis.
5. Add optimized parameter sweeps and experimental comparison metrics.

## Notes

The DiffuseCode repository contains an initial, explicitly *not release-ready* F2PY bridge in `python/suite_python.f90`. This MVP does not depend on it. Inspect DiffuseCode's license and check permission for intended commercial use.
## Explicit powder settings

`PDFConfig` exposes `qmin`, `qmax`, `dq`, `qsig` (DISCUS `set delta`), and optional `wavelength` (Å), in addition to `rmax` and `dr`. They are saved with every run. `qmax=0` is rejected: omitting the command does not establish that the termination correction is disabled. The real-space range must be divisible by `dr`.

Run `python -m pytest -q` for unit tests; see the separately documented manual DISCUS integration check.

For new work, pass a ReciprocalVolume carrying explicit coordinate_system='cartesian_q', units='1/angstrom', and confirmed Bragg removal metadata. See [3D integration notes](docs/3D_INTEGRATION.md). This does not by itself validate quantitative 3D-ΔPDF interpretation.
