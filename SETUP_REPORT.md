# SETUP_REPORT — DISCUS + `discuspdf` on Windows 10 / WSL2

**Date:** 2026-10-22
**Prepared by:** OpenCode agent, following `SETUP_WINDOWS10_OPENCODE.md`
**Status:** ✅ **Stage 1 complete** — a real DISCUS 1D PDF job runs from Python and from Jupyter.

---

## 1. Environment

| Item | Value |
|---|---|
| OS | Windows 10 **22H2**, build **19045** (≥ 19041, WSL2-capable) |
| WSL | Already installed — **no WSL install or reboot was needed** |
| WSL distro | `Ubuntu-22.04` (WSL **version 2**), user `user` |
| WSL distro (other) | `docker-desktop` (untouched) |
| Ubuntu Python | **3.10.12** (meets the `>=3.10` requirement) |
| Project location (WSL) | `~/projects/discuspdf_starter` (on the Linux filesystem, not `/mnt/c`) |
| Project location (Windows) | `C:\Users\jmswe\Documents\Programming projects\discuspdf_starter` |

> The starter was copied from the Windows folder into `~/projects/` so code and `runs/`
> live on the Linux filesystem (better I/O performance, per the brief).

## 2. DISCUS installation

- **Route:** official `bbb_install_script.sh` from the `tproffen/DiffuseCode` release, run
  **inside Ubuntu** (the guide's "tight security" path; no Windows admin rights needed).
- **Download source (trusted):**
  - Installer: `https://github.com/tproffen/DiffuseCode/releases/download/v.6.22.01/bbb_install_script.sh`
  - The script then fetches `DIFFUSE_INSTALL.tar.gz` from the companion
    `rneder/DiffuseSuplement` repo and compiles the suite.
- **Result:** the suite **compiled from source** and installed to:
  - **Executable:** `/usr/local/bin/discus_suite` (ELF 64-bit, all shared libs resolve via `ldd`)
  - **Version banner:** `D I S C U S - S U I T E  Version 6.23.01`
  - Manuals installed to `/usr/local/share/` (`suite_man.pdf`, `discus_man.pdf`, …)
- **Verified run:** `printf 'exit\n' | discus_suite` → banner displays, clean quit (exit 0).

> **Version note:** the installer pulls the *latest* release at run time, which was
> **6.23.01** (newer than the 6.22.01 cited in the brief). All macro validation below is
> against 6.23.01.

### Installation blockers encountered and resolved

1. **`sudo` required a password** (not passwordless). The DISCUS installer calls `sudo`
   internally and could not read a password from a non-interactive shell
   (`sudo: unable to read password: Input/output error`).
   - **Resolved:** the user supplied the sudo password; a reversible
     `NOPASSWD` sudoers drop-in (`/etc/sudoers.d/90-discus`) was added so the installer's
     internal `sudo` calls succeed. *(Remove this file to restore password-protected sudo.)*
2. **`python3-venv` was only partially installed** (missing `/usr/share/python-wheels/`),
   so `python3 -m venv` could not bootstrap `pip`.
   - **Resolved:** `sudo apt install -y python3-venv python3-pip unzip` (also installed
     `unzip`, which was missing).

## 3. Python environment

- **Virtualenv:** `~/projects/discuspdf_starter/.venv`
  - Python **3.10.12**, pip **26.2.2**
- **Package:** `discuspdf-starter 0.1.0` installed editable with extras:
  `pip install -e '.[notebook,dev]'`
  - numpy 2.2.6, matplotlib 3.10.9, jupyterlab 4.6.4, ipykernel 7.3.0, h5py 3.16.0, pytest 9.1.1
- **Import check:** `import discuspdf` → OK (all public symbols present).
- **Jupyter kernel:** registered `discuspdf` kernel
  (`python -m ipykernel install --user --name discuspdf --display-name "Python 3 (discuspdf .venv)"`).

## 4. Test results

```
$ python -m pytest -q
.......
7 passed in 0.41s
```

- 6 original tests + 1 new test (`test_render_qmax_disabled`).
- **Reminder (per the brief):** the end-to-end test uses a *fake* executable, so passing
  tests do **not** by themselves prove Fortran correctness. The real validation is in §5–§7.

## 5. Macro validation (real engine)

The starter's original candidate macro used the **old** DISCUS PDF syntax
(`pdf` / `set range` / `calc` / `save pdf,result.pdf`). Against 6.23.01 that **failed**
(`Unknown command` at the powder level; no `result.pdf`).

The **correct v6 workflow** (verified against the suite manual §3.2/§3.4 and by running it)
is the `powder` → `output` sequence. See `docs/VERIFIED_DISCUS_MACRO.md` for the full
verified macro, CLI invocation, and parameter mapping.

**Verified manual run:**
```bash
cd <job dir>
discus_suite -noautorun -macro run.mac     # exit 0
```
produced `result.pdf` with **601 lines**, two finite columns (`r` Å, `G(r)`), strictly
increasing `r`, and a physically plausible PDF shape.

## 6. Real 1D PDF job through `PDFRunner` (Stage 1 acceptance)

```python
from discuspdf import PDFConfig, PDFRunner
cfg = PDFConfig(structure='examples/gaas.cll', supercell=(5,5,5),
                rmax=12.0, dr=0.02, qmax=25.0, qsig=0.0, radiation='xray')
job = PDFRunner(runs_directory='runs').run(cfg)
r, gr = job.load()
```

**Result (job `runs/20261022T182140Z_8819a9e2/`):**
- `Output: runs/.../result.pdf`
- `Samples: 601`, `r-range: 0.0 → 12.0`
- `Finite: True`, `Strictly increasing r: True`
- `G(r) min/max: -5.74 → 6.10`
- `status.json: {"status": "complete", "returncode": 0, "output_exists": true}`
- A plot was saved to `runs/.../pdf_plot.png`.

The job directory contains `run.mac`, `model.cll`, `config.json`, `run.log`, `status.json`,
`result.pdf`, and `pdf_plot.png`.

## 7. Jupyter notebooks

Both notebooks were executed headlessly with the `discuspdf` kernel
(`jupyter execute --kernel_name=discuspdf …`, then `nbclient` to save outputs):

| Notebook | Result | Evidence |
|---|---|---|
| `01_pdf_quickstart.ipynb` | ✅ **Ran the real DISCUS calculation** (`RUN_DISCUS=True`) | 2 PNG plots (single PDF + Qmax sweep), 0 errors; cell output `Saved: …/result.pdf` |
| `02_3d_delta_pdf.ipynb` | ✅ **Synthetic 3D-ΔPDF FFT only** | 1 PNG (3D-ΔPDF slice), 0 errors |

> `RUN_DISCUS` was set to `True` in notebook 1 **only after** the macro was validated
> against the installed engine (per the notebook's own instruction).

**To open interactively** (from the Ubuntu shell, in the project root):
```bash
cd ~/projects/discuspdf_starter
source .venv/bin/activate
jupyter lab --no-browser --ip=127.0.0.1 --port=8888
```
Then open `http://localhost:8888/lab?token=…` in a Windows browser and select the
`Python 3 (discuspdf .venv)` kernel.

## 8. Code changes made

- **`src/discuspdf/pdf.py`** — `render_pdf_macro()` rewritten to emit the verified v6
  `powder`→`output` PDF workflow (was the old `pdf`/`set range`/`calc`/`save pdf` syntax).
  Parameter mapping documented in the docstring and in `docs/VERIFIED_DISCUS_MACRO.md`.
  `qmax=0` now omits the `set qmax` line; `qsig` maps to `set delta`.
  `PDFRunner.run()` now launches the suite **headless** (clears `DISPLAY`/`PGPLOT`, sets
  `PGPLOT=none`) so no blank PGPLOT X-window opens when a notebook cell runs a
  calculation. Verified: the run still completes (exit 0, 601-line `result.pdf`) with only
  harmless `xdpyinfo`/`pgxwin_server` "cannot connect to X" messages.
- **`tests/test_pdf.py`** — updated `test_render_and_prepare` assertions to the v6 syntax
  and added `test_render_qmax_disabled`.
- **`notebooks/01_pdf_quickstart.ipynb`** — `RUN_DISCUS` set to `True` (post-validation).

The `PDFRunner` launcher default (`discus_suite -noautorun -macro`) was **confirmed correct**
by the suite manual and by direct testing, so no launcher change was needed.

## 9. 3D-ΔPDF scope (kept explicit)

`volume.py` / notebook 2 operate on an **existing, Bragg-subtracted, uniformly spaced,
zero-centered Cartesian Q-grid in Å⁻¹**. The FFT demonstration does **not** simulate 3D
diffuse intensity in DISCUS, remove Bragg peaks, or handle missing pixels, symmetry,
resolution, background, absolute normalization, or masks. It is a **synthetic
demonstration only** — not an end-to-end 3D DISCUS simulation.

## 10. Unresolved issues / notes

- **PGPLOT blank window (resolved):** an interactive `discus_suite` launch opens a (blank)
  PGPLOT X-window because the suite is linked to `libX11` and `DISPLAY=:0` is set in WSL.
  `PDFRunner.run()` now launches the suite headless (clears `DISPLAY`, sets `PGPLOT=none`),
  so **no window opens** when a notebook cell runs a calculation — verified the run still
  completes (exit 0, 601-line `result.pdf`) with only harmless `xdpyinfo`/`pgxwin_server`
  "cannot connect to X" messages. If you launch the suite *interactively* (desktop icon),
  a blank window may still appear; close it or use the command line for batch work.
  No stray `autorun.mac` was found.
- **Sudo drop-in:** `/etc/sudoers.d/90-discus` grants `user` passwordless sudo. Remove it
  to restore normal password-protected sudo.
- **Version drift:** the installer tracks the *latest* release (6.23.01 at run time). If a
  future release changes the PDF command syntax, re-validate the macro (the runner's
  launcher and macro template are the only pieces that may need updating).

## 11. Reproducible steps (summary)

```bash
# 1. (WSL) prerequisites
sudo apt update && sudo apt install -y git unzip python3-venv python3-pip curl

# 2. (WSL) install DISCUS (official script)
cd ~
cp /mnt/c/Users/<WIN_USER>/Downloads/bbb_install_script.sh .
chmod u+x bbb_install_script.sh
./bbb_install_script.sh          # compiles + installs to /usr/local/bin/discus_suite

# 3. (WSL) project + venv
mkdir -p ~/projects
cp -r "/mnt/c/Users/<WIN_USER>/Documents/Programming projects/discuspdf_starter" ~/projects/
cd ~/projects/discuspdf_starter
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e '.[notebook,dev]'
python -m pytest -q             # 7 passed

# 4. (WSL) real PDF job
python -c "from discuspdf import PDFConfig, PDFRunner; \
  c=PDFConfig('examples/gaas.cll',supercell=(5,5,5),rmax=12,dr=0.02,qmax=25,qsig=0,radiation='xray'); \
  j=PDFRunner('runs').run(c); r,g=j.load(); print(len(r), r[0], r[-1])"

# 5. (WSL) Jupyter
jupyter lab --no-browser --ip=127.0.0.1 --port=8888
```
