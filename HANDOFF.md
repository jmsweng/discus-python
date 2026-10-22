# HANDOFF — discus-python (DISCUS 1D PDF automation)

**Purpose of this file:** a self-contained brief for an orchestrating model/agent. It states what
exists, what is verified, how to run things, and what remains. Read this before acting.

---

## 1. One-paragraph status

A Python package `discuspdf` wraps the DISCUS_SUITE Fortran executable to run **1D atomic
pair-distribution-function (PDF)** calculations from a `.cll` unit-cell file, plus a
**synthetic-only** 3D-ΔPDF FFT helper. The full Windows 10 + WSL2 setup was completed and
**verified end-to-end**: DISCUS_SUITE 6.23.01 is installed in Ubuntu (WSL2), the Python
package is installed in a venv, all tests pass, a **real** DISCUS 1D PDF job runs from
Python and from Jupyter, and the project is pushed to a private GitHub repo.

**Stage 1 (1D PDF): COMPLETE and verified.** 3D-ΔPDF: synthetic demonstration only (not a
real DISCUS 3D simulation).

## 2. Repository

- **GitHub:** `https://github.com/jmsweng/discus-python` (private), branch `main`
- **Local (Windows):** `C:\Users\jmswe\Documents\Programming projects\discus-python`
- **Local (WSL, where DISCUS + venv live):** `~/projects/discuspdf_starter`
  (the WSL working copy is still named `discuspdf_starter`; only the Windows folder was
  renamed to match the repo. Rename the WSL copy too if you want the names to agree.)
- Keep DISCUS and the Python venv in the **same OS context** (WSL Ubuntu). Do not mix native
  Windows Python with the Linux `discus_suite`.

### Layout
```
discus-python/
├─ pyproject.toml            # package: discuspdf-starter 0.1.0; needs Python >=3.10
├─ src/discuspdf/
│  ├─ __init__.py            # public API
│  ├─ pdf.py                # PDFConfig, PDFRunner, PDFRun, read_pdf, render_pdf_macro
│  └─ volume.py             # delta_pdf_from_diffuse, read_hdf5_volume, plot_delta_pdf_slice (SYNTHETIC)
├─ examples/
│  ├─ gaas.cll              # GaAs unit cell (F-43m) — the standard test structure
│  └─ run_sample/           # ONE real DISCUS run (reference output; see its README)
├─ notebooks/
│  ├─ 01_pdf_quickstart.ipynb   # REAL 1D PDF (RUN_DISCUS=True)
│  └─ 02_3d_delta_pdf.ipynb     # SYNTHETIC 3D-ΔPDF FFT only
├─ tests/                    # test_pdf.py, test_volume.py  (7 tests)
├─ docs/VERIFIED_DISCUS_MACRO.md   # verified macro syntax + CLI for the installed version
├─ SETUP_REPORT.md           # full setup log, versions, results, blockers
└─ .gitignore               # excludes .venv/, runs/, caches
```

## 3. Environment (verified)

| Item | Value |
|---|---|
| OS | Windows 10 22H2, build 19045 |
| WSL distro | `Ubuntu-22.04`, WSL **version 2**, user `user` |
| DISCUS | `/usr/local/bin/discus_suite`, **version 6.23.01** (installed via official `bbb_install_script.sh`) |
| Python (WSL) | 3.10.12 |
| venv | `~/projects/discuspdf_starter/.venv` (pip 26.2.2) |
| Package | `discuspdf-starter 0.1.0` (editable, extras `notebook,dev`) |
| Jupyter kernel | `discuspdf` (registered via ipykernel) |

## 4. How to run (in the WSL Ubuntu shell)

```bash
cd ~/projects/discuspdf_starter
source .venv/bin/activate
python -m pytest -q                 # 7 passed
```

### One real 1D PDF job from Python
```python
from discuspdf import PDFConfig, PDFRunner
cfg = PDFConfig(structure="examples/gaas.cll", supercell=(5, 5, 5),
                rmax=12.0, dr=0.02, qmax=25.0, qsig=0.0, radiation="xray")
job = PDFRunner(runs_directory="runs").run(cfg)
r, gr = job.load()                # r: Å (increasing), gr: G(r)
job.plot()                        # matplotlib
```
Each run writes an isolated dir `runs/<id>/` with `run.mac`, `model.cll`, `config.json`,
`run.log`, `status.json`, `result.pdf` (two-column ASCII: `r` Å, `G(r)` — **not** an Adobe PDF).

### Jupyter
```bash
jupyter lab --no-browser --ip=127.0.0.1 --port=8888
# open http://localhost:8888/lab?token=... in a Windows browser; select the
# "Python 3 (discuspdf .venv)" kernel
```

## 5. Verified DISCUS macro + CLI (IMPORTANT — do not use the old syntax)

The starter's **original** candidate macro used the OLD DISCUS PDF syntax
(`pdf` / `set range` / `calc` / `save pdf,result.pdf`). That **FAILS** on v6.23.01.
The **verified working** invocation and macro are:

```bash
discus_suite -noautorun -macro run.mac        # exit 0; writes result.pdf
```
```text
discus
read
cell model.cll,5,5,5
powder
reset
xray
set axis,q
set calc,debye
set qmin,0.5
set qmax,25.0
set dq,0.01
set delta,0.0
set wvle,1.5406
run
exit
output
outf result.pdf
form pdf, r, 0.0, 12.0, 0.02
value PDF
run
exit
exit
```
Full parameter mapping and verification evidence: `docs/VERIFIED_DISCUS_MACRO.md`.

**Headless note:** `PDFRunner.run()` launches the suite with `DISPLAY` cleared and
`PGPLOT=none` so no blank PGPLOT X-window opens (the suite is linked to `libX11`). This is
verified — the run still completes (exit 0) with only harmless `xdpyinfo`/`pgxwin_server`
"cannot connect to X" messages.

## 6. What is / is not true (scope guardrails)

- **True & verified:** 1D atomic PDF from a `.cll` file, via the real DISCUS engine, from
  Python and Jupyter.
- **Synthetic only (NOT a real DISCUS 3D simulation):** `volume.py` / notebook 2 FFT an
  *already Bragg-subtracted, uniformly spaced, zero-centered Cartesian Q-grid in Å⁻¹*.
  It does **not** simulate 3D diffuse intensity, remove Bragg peaks, or handle missing
  pixels, symmetry, resolution, background, absolute normalization, or masks.
- Passing the Python unit tests does **not** prove Fortran correctness (the end-to-end test
  uses a fake executable). Real correctness was established by the actual DISCUS runs above.
- Do **not** claim quantitative 3D-ΔPDF accuracy; benchmark against a trusted independent
  workflow before comparing to experimental data.

## 7. Known issues / caveats

- **Version drift:** the DISCUS installer tracks the *latest* release (was 6.23.01 at install
  time). If a future release changes PDF command syntax, re-validate the macro (only
  `render_pdf_macro()` / the launcher may need updating).
- **Sudo drop-in:** `/etc/sudoers.d/90-discus` grants `user` passwordless sudo (added to
  unblock the installer). Remove it to restore normal password-protected sudo.
- **PGPLOT window on interactive launch:** double-clicking the `DiscusSuite_WSL` desktop icon
  opens the interactive GUI + a PGPLOT window. Use the command line / `PDFRunner` for batch.
- **Credentials:** a GitHub SSH key (`~/.ssh/id_ed25519`) and `gh` CLI auth are set up for
  `jmsweng`. The repo is private.

## 8. Suggested next tasks (priority order)

1. **HDF5 adapter (version-aware):** v6.22+ writes unified HDF5 structure/diffraction formats
   (with a `coordinate_unit` field in 6.22.01+). Inspect a *real* DISCUS HDF5 output first
   (do not guess dataset paths), then add a reader that maps it into `volume.py`.
2. **3D-ΔPDF from DISCUS:** drive a real DISCUS 3D single-crystal → PDF/ΔPDF macro and feed
   its Bragg-subtracted volume into `delta_pdf_from_diffuse`; validate against an
   independent workflow before trusting numbers.
3. **Parameter sweeps / disorder:** extend `PDFConfig` for the disorder/defect command blocks
   the user changes most often; add sweep helpers.
4. **`.cif` import** and advanced atom-pair selection (currently `.cll` only).
5. **Packaging:** consider a real release (the PyPI name `discus` is taken by something else;
   do not `pip install discuspdf` — it's a local editable project).

## 9. Key files to read first

1. `docs/VERIFIED_DISCUS_MACRO.md` — the ground truth for macro/CLI on the installed version.
2. `SETUP_REPORT.md` — full setup log, versions, test results, code changes, blockers.
3. `src/discuspdf/pdf.py` — the runner + macro template (the part most likely to need
   updating for a new DISCUS version).
4. `examples/run_sample/README.md` — what a real run's outputs look like.
