# DISCUS + `discuspdf` on Windows 10 — local OpenCode implementation brief

**Document type:** executable setup runbook / task specification for a local coding agent  
**Prepared:** 2026-10-22  
**Project:** Python/Jupyter automation of DISCUS 1D pair distribution functions, with an exploratory 3D-ΔPDF analysis module  
**Upstream engine:** [tproffen/DiffuseCode](https://github.com/tproffen/DiffuseCode)  
**Existing prototype:** `discuspdf_starter.zip` / `discuspdf_starter/` supplied with this brief  
**Target:** 64-bit Windows 10; **preferred execution environment:** Ubuntu on WSL2

> **Agent: read this entire brief before installing anything.** Distinguish VERIFIED upstream facts from commands that must be tested on the local computer. Do not report success solely because Python unit tests pass; a real Fortran/DISCUS run is required.

## 0. Mission and boundaries

Your job is to deliver a *working, reproducible local setup* for this user:

1. Get an official version of DISCUS_SUITE running on Windows 10 (prefer Ubuntu on WSL2).
2. Install the pre-existing `discuspdf` Python package into a virtual environment in the *same OS context* as the DISCUS executable.
3. Make one real 1D PDF simulation run from the Python API, save a numeric `G(r)` series, and display a plot in Jupyter.
4. Confirm that notebook 2 runs the **synthetic** 3D-ΔPDF FFT example; do not misrepresent this as an end-to-end 3D DISCUS simulation.
5. Produce `SETUP_REPORT.md` with what was installed, actual command paths, version information, verified results, changes made, and remaining problems.

**Human approval required:** before enabling Windows features, invoking elevated PowerShell, rebooting, executing newly downloaded installer scripts/binaries, making system-wide PATH changes, changing firewall/security policy, or deleting/overwriting pre-existing data. You may inspect/download official releases and work within the new project directory as allowed. Do not execute opaque remote scripts (`curl ... | bash`, `irm ... | iex`) without source review and permission. Do not disable endpoint security to get an installer running. Maintain a change log.

**Non-goals:** recompiling the Fortran engine unless necessary; claiming quantitative accuracy of the 3D-ΔPDF FFT; or building a new GUI. Avoid mixing native Windows Python with Linux `discus_suite` until a deliberate adapter has been implemented and tested.

## 1. Why WSL2 is the recommended layout

```
Windows 10 desktop
  ├─ Windows browser → Jupyter at http://localhost:8888
  └─ Ubuntu (WSL2)
       ├─ ~/projects/discuspdf_starter/
       │    ├─ .venv/              Python + Jupyter + numpy etc.
       │    ├─ src/discuspdf/      Python runner + analysis
       │    ├─ notebooks/         user notebooks
       │    └─ runs/              isolated simulation outputs
       └─ installed DISCUS_SUITE binary
            └─ invoked locally from the same WSL Python process
```

Keeping both programs inside Ubuntu avoids Windows/Unix executable/path translation issues and is aligned with the upstream project's Windows-WSL installation scripts. WSL2 support for the standard `wsl --install` workflow requires Windows 10 version 2004, build 19041 or later; see Microsoft documentation. If this PC cannot support WSL2, follow the **native-Windows fallback** near the end instead of inventing an unsupported install path.

## 2. Prepare Windows and identify the execution context

**In a normal Windows PowerShell terminal** (NOT inside Ubuntu):

```powershell
winver
wsl --status
wsl -l -v
```

Record Windows version/build and whether Ubuntu/WSL is already installed. `wsl --status` may not be supported on older builds; handle this gracefully. If Ubuntu is already installed and WSL version is 2, **reuse it** instead of reinstalling.

If WSL2 is absent and the user approves administrative changes, open **PowerShell as Administrator** and run:

```powershell
wsl --install -d Ubuntu
```

Restart Windows **only after user approval**, launch Ubuntu from the Start menu, and create the Ubuntu user when prompted. If Ubuntu is already installed as WSL1, verify before converting it:

```powershell
wsl -l -v
wsl --set-version Ubuntu 2
```

If installation reports unavailable virtualization or unsupported Windows build, STOP and report the specific prerequisite with links to the Microsoft manual installation guide. Do not silently switch to experimental approaches.

**In an Ubuntu / WSL terminal** (all commands below use bash unless labeled PowerShell):

```bash
sudo apt update
sudo apt install -y git unzip python3 python3-venv python3-pip
python3 --version
uname -a
mkdir -p ~/projects
```

Only add compilers/system packages if a source build is explicitly needed. Keep active code and `runs/` on the Linux filesystem (e.g., `~/projects`), rather than under `/mnt/c`, for more predictable performance.

## 3. Install DISCUS_SUITE from an official upstream release

The DiffuseCode README recommends the GitHub release packages for Windows. At preparation time, GitHub identifies **v.6.22.01** as the latest release; its release notes describe HDF5 coordinate unit metadata. The release list and asset names may change in the future, so inspect the current release and installation guide rather than hard-coding an asset filename.

Official release page: https://github.com/tproffen/DiffuseCode/releases/latest  
Specific verified version: https://github.com/tproffen/DiffuseCode/releases/tag/v.6.22.01

**Agent procedure:**

1. Open the release page and expand **Assets**. Identify the Windows-specific installation guide, and note whether the supported route is Windows executable, Windows+WSL helper scripts, or a Linux package to run inside WSL.
2. Read the install guide for that **same release**, including its prerequisite list. Upstream previously mentioned Windows-WSL scripts named `bbb_install_suite_Windows_WSL.bat` and `bbb_install_Windows_wsl.ps` (v.6.18.02). **Do not assume these are present in v.6.22.01** or copy commands blindly from an older release.
3. With user approval, install through the supported method. Prefer the release's WSL-compatible instructions for this setup. If the upstream script does not support Windows 10 / the installed WSL environment, document the blocker and use an official supported fallback.
4. Do **not** assume that a generic ZIP of source files is a runnable Windows binary. Do **not** download binaries from third-party mirrors.
5. Discover the actual executable name and directory; installation scripts may use paths or program names different from `discus_suite`.

In Ubuntu, test executable discovery:

```bash
command -v discus_suite || true
command -v discus_suite.exe || true
find "$HOME" /usr/local /opt -type f \
    \( -iname 'discus_suite' -o -iname 'discus_suite.exe' \) \
    2>/dev/null | head -25
```

If the launcher lives in `~/.local/bin` and that directory is not on PATH, add an appropriate PATH entry to the **user's** shell configuration, after permission; avoid system-wide changes. If the upstream installed a wrapper script, retain it and use its documented invocation.

Once an executable is identified, test it interactively first. For a Linux executable found as `discus_suite`, for example:

```bash
command -v discus_suite
printf 'exit\n' | timeout 30s discus_suite
```

An interactive suite prompt/banner followed by a clean quit is the expected sign. Consult the release's guide if the executable requires a different exit command or cannot read piped input. Save its banner and version to `SETUP_REPORT.md`. **Do not assume `--version` or `--help` flags exist.**

**Source-build fallback, only if official binaries/installer cannot work:** follow the upstream Linux installation guide distributed with the same release. Its CMake build requires a Fortran compiler specified by the `FC` environment variable and third-party libraries. Possible Ubuntu build prerequisites include `gfortran`, CMake, FFTW, Readline, and HDF5; exact versions/flags must come from the current guide. Do not guess a full build command or present a successful source build without testing it.

## 4. Bring the pre-existing Python starter into WSL

Two files should have been downloaded on Windows:

- This task document (`SETUP_WINDOWS10_OPENCODE.md`) or the combined bundle (`discuspdf_windows10_opencode_bundle.zip`).
- `discuspdf_starter.zip` if **not** using the combined bundle.

**Preferred:** extract the combined bundle into the WSL filesystem. First locate the Windows account folder under `/mnt/c/Users/`:

```bash
ls /mnt/c/Users
```

Substitute your real Windows account folder for `<WIN_USER>` in ONE of the following:

```bash
# Combined bundle (preferred):
mkdir -p ~/projects
unzip "/mnt/c/Users/<WIN_USER>/Downloads/discuspdf_windows10_opencode_bundle.zip" \
      -d ~/projects

# OR the original starter ZIP (then copy this Markdown into the project):
unzip "/mnt/c/Users/<WIN_USER>/Downloads/discuspdf_starter.zip" \
      -d ~/projects
```

Both ZIPs contain the `discuspdf_starter/` root directory. Agent: handle cases where the browser saved into Desktop or another folder; search existing files before asking the user to re-download. Do not overwrite an existing working directory without a backup and approval.

```bash
cd ~/projects/discuspdf_starter
ls -la
cat pyproject.toml
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e '.[notebook,dev]'
python -m pytest -q
```

Expected Python version: **3.10+**. `pyproject.toml` defines runtime `numpy`, `matplotlib` and extras for Jupyter, `h5py`, and `pytest`. Do not install the PyPI package named `discus` in place of the upstream Fortran suite. Do not attempt `pip install discuspdf`; this is a *local editable project*, installed with `-e .` plus its extras.

**Important:** passing Python tests only proves code-level behavior; the existing tests use a fake DISCUS executable in the subprocess test. They do NOT prove the upstream `.cll` format, CLI options, PDF macro commands, or numeric output are compatible with the installed DISCUS release.

## 5. Inspect the generated 1D PDF macro before executing

The existing starter's 1D runner is in `src/discuspdf/pdf.py`. It builds `run.mac` with commands similar to:

```text
# Generated by discuspdf: 1D atomic PDF calculation
discus
read
cell model.cll,5,5,5
pdf
set range,12,0.02
set qmax,25
set qsig,0
set rad,xray
calc
save pdf,result.pdf
exit
exit
exit
```

**Treat this as a candidate macro**, based on documented older PDF examples, not as proven syntax for v.6.22.01. Also validate the example structure `examples/gaas.cll` (and its atom/space-group conventions) against the installed parser. The extension `result.pdf` here means *numeric pair-distribution output*, not a printable Adobe PDF document.

To generate a preview without executing DISCUS:

```bash
source .venv/bin/activate
python - <<'PY'
from discuspdf import PDFConfig, PDFRunner
cfg = PDFConfig(
    structure='examples/gaas.cll',
    supercell=(5, 5, 5),
    rmax=12.0, dr=0.02,
    qmax=25.0, qsig=0.0,
    radiation='xray',
)
job = PDFRunner(runs_directory='runs').prepare(cfg)
print('Run directory:', job.directory)
print((job.directory / 'run.mac').read_text())
PY
```

Manually run that generated macro with the installed suite, using the current release's supported invocation. The starter defaults to the argument sequence:

```bash
# Run from the generated job directory, not the project root:
# `cd` into the printed Run directory first.
printf 'exit\n' | timeout 120s discus_suite -noautorun -macro run.mac
```

**These CLI flags require local verification.** If the installed suite does not accept them, find the supported macro-launch syntax in its docs/help and record it. A manual interactive command such as `@run.mac` *may* work in the suite, but test it rather than assuming. Capture stdout/stderr, inspect output and command errors, and make the runner match the real invocation.

After the macro completes, check that a nonempty `result.pdf` exists in the job directory and contains two or more finite numeric columns. Independently review the output to determine which column is `r` and which is `G(r)`; do not silently treat every `.pdf` output as equivalent.

## 6. Run a real PDF calculation from Python

Use the **same Ubuntu shell** and activated virtual environment. Only after the manual macro succeeds:

```bash
cd ~/projects/discuspdf_starter
source .venv/bin/activate
python - <<'PY'
from discuspdf import PDFConfig, PDFRunner
from pathlib import Path
import numpy as np

cfg = PDFConfig(
    structure=Path('examples/gaas.cll'),
    supercell=(5, 5, 5),
    rmax=12.0, dr=0.02, qmax=25.0,
    qsig=0.0, radiation='xray',
)
runner = PDFRunner(runs_directory='runs')
# If discovery found a different executable/flags, pass launcher=[...]
job = runner.run(cfg)
r, gr = job.load()
print('Output:', job.output)
print('Samples:', len(r), 'r-range:', r[0], r[-1])
print('Finite:', np.isfinite(r).all() and np.isfinite(gr).all())
print('Status:', (job.directory/'status.json').read_text())
PY
```

If necessary, modify only the appropriate launcher/template parts of `src/discuspdf/pdf.py` and add/update tests so that the corrected invocation is reproducible. Do not hard-code the original author's Windows username or install location. Prefer a discoverable/overridable launcher, for example `PDFRunner(launcher=['/actual/path/to/discus_suite', '-noautorun', '-macro'])`, after validating which arguments the current executable accepts.

**Acceptance for stage 1:** an actual DISCUS-generated `result.pdf`, with a log and recorded configuration, can be parsed by `PDFRun.load()`; numeric data are finite, the radius axis is increasing, and the user can see a plausible PDF plot. Report any discrepancies with an independently produced/reference curve rather than silently guessing a physical explanation.

## 7. Open JupyterLab and validate both notebooks

```bash
cd ~/projects/discuspdf_starter
source .venv/bin/activate
jupyter lab --no-browser --ip=127.0.0.1 --port=8888
```

Copy the URL with its authentication token from the Ubuntu terminal and open it in a **Windows** browser (normally `http://localhost:8888/lab?...`). Do not disable authentication or bind Jupyter publicly to `0.0.0.0` merely for convenience. If `localhost` forwarding fails on this WSL version, diagnose network forwarding rather than changing firewall rules without approval.

- `notebooks/01_pdf_quickstart.ipynb`: must run the real calculation only after the runner is validated; confirm the plot and paths.
- `notebooks/02_3d_delta_pdf.ipynb`: synthetic FFT/viewer validation only. This does **not** simulate 3D diffuse intensity in DISCUS nor remove Bragg peaks from experimental data.

If the notebook working directory is `notebooks/` instead of the project root, follow its documented directory adjustment. Use the Jupyter kernel belonging to `.venv`, not an unrelated Windows Python interpreter.

## 8. 3D-ΔPDF scientific limitations — keep explicit

The prototype's `volume.py` operates on an existing **Bragg-subtracted**, uniformly spaced, zero-centered **Cartesian Q-grid** in Å⁻¹. Do not treat crystallographic `(h,k,l)` reciprocal-lattice units as physical Cartesian Q without a metric conversion. The FFT demonstration does not automatically account for missing pixels, symmetry, resolution effects, instrumental/background corrections, absolute normalization, masks, or proper Bragg removal.

Version v.6.22.00 introduced unified HDF5 structure/diffuse formats and expanded KUPLOT's 3D capabilities; v.6.22.01 added a coordinate-unit field. An eventual next milestone is a *version-aware HDF5 adapter*, but **do not invent HDF5 dataset paths**; inspect real files first.

Before comparing to experimental 3D-ΔPDF data, require scientifically validated Q-space processing and benchmark against a trusted independent workflow.

## 9. Troubleshooting guide

| Symptom | Action |
|---|---|
| `wsl --install` unsupported | Check Windows version/build; follow Microsoft's older/manual WSL instructions; require approval for system features/reboot. |
| WSL reports virtualization disabled | Report BIOS/UEFI virtualization requirement to user; do not alter firmware settings automatically. |
| `discus_suite: command not found` | Locate executable; verify installer completed and current Ubuntu PATH; do not guess binary name. |
| `cannot execute binary file` | Check `file <path>` and whether it is a native Windows executable or Linux ELF; choose matching OS environment. |
| Missing `.so`/DLL | Run `ldd` for ELF libraries or use upstream supported installer; avoid random DLL downloads. |
| `Unknown option -macro` | Test release-specific CLI; update launcher only after validated manual execution. |
| `result.pdf` absent | Inspect generated `run.mac`, `run.log`, and DISCUS error messages; validate structure and PDF command syntax. |
| File exists but parser fails | Inspect headers/columns, numeric formats, and units; make parsing schema explicit. |
| Python tests pass but real job fails | Expected possibility: simulated tests cannot validate the Fortran executable; debug real integration. |
| Jupyter can't import `discuspdf` | Select the `.venv` kernel and confirm `python -m pip show discuspdf-starter`. |
| Jupyter URL doesn't open in Windows | Check WSL localhost forwarding, server bind address, port, and token. |
| Windows PowerShell can't run Ubuntu `discus_suite` | Prefer OpenCode and Python in Ubuntu; do not mix environments through relative paths. |

## 10. Native Windows fallback (only if WSL2 is unavailable)

The DiffuseCode README mentions Windows binary releases. If the release **actually provides** a compatible Windows 10 x64 executable/installer, an all-native workflow is possible. Validate installer, bitness, dependencies and supported version before adopting it.

1. Install DISCUS using the *official* Windows binary and user-approved installer steps.
2. Install Python **3.10+** for Windows from a trusted distributor; verify `py -3 --version`.
3. Extract `discuspdf_starter.zip` to a clean path such as `%USERPROFILE%\Documents\discuspdf_starter`.
4. In native PowerShell, inside that directory:

```powershell
py -3 -m venv .venv
.\.venv\Scripts\python.exe -m pip install --upgrade pip
.\.venv\Scripts\python.exe -m pip install -e '.[notebook,dev]'
.\.venv\Scripts\python.exe -m pytest -q
```

5. Set `PDFRunner(launcher=[r'C:\actual\location\discus_suite.exe', ...])` with the **verified supported CLI arguments**; run the real PDF test and notebooks on the *Windows* Python kernel.
6. In PowerShell, launch Jupyter with `& .\.venv\Scripts\jupyter-lab.exe` (or `& .\.venv\Scripts\python.exe -m jupyterlab` if supported). Never use `source .venv/bin/activate` in PowerShell.

If there is no functioning native binary, do not represent this fallback as available. Report that the official WSL workflow is required or propose a user-approved alternative.

## 11. Expected agent outputs and stop conditions

Produce these files under the project root:

- `SETUP_REPORT.md`: Windows/WSL version; Python version; DISCUS version; install route and trustworthy download URL; executable path; environment directory; original commands used; tests; output file path; screenshot/plot location if generated; unresolved issues.
- `docs/VERIFIED_DISCUS_MACRO.md`: actual working macro syntax/flags and how they were verified (create `docs/` if needed).
- Changes to `src/discuspdf/pdf.py` and its tests **only if needed** for the real installation. Describe diffs and preserve parameter tracking.
- A real `runs/<job_id>/` output containing `run.mac`, `config.json`, `run.log`, `status.json`, and numeric `result.pdf` if stage 1 succeeds.

**Explicit pass/fail checklist:**

- [ ] Windows build and WSL status recorded.
- [ ] Official DISCUS executable installed and run at least once.
- [ ] Local Python `discuspdf` installed in dedicated `.venv`.
- [ ] `python -m pytest -q` passes (number reported).
- [ ] Generated macro inspected, release-compatible syntax confirmed.
- [ ] One genuine DISCUS job launched from Python successfully.
- [ ] `result.pdf` numeric series parsed and 1D PDF plotted.
- [ ] Jupyter notebook 1 verified with real engine.
- [ ] Notebook 2 verified as synthetic-only.
- [ ] Final `SETUP_REPORT.md` gives reproducible steps and blockers.

If the executable cannot be installed or another required authorization is missing, stop at the last demonstrated checkpoint and provide a concrete report. **Do not fabricate a successful run, output, or version.**

## 12. Prompt to give OpenCode

Open OpenCode in `~/projects/discuspdf_starter` (for a WSL install, start OpenCode **inside Ubuntu** if possible). OpenCode recognizes a project-root `AGENTS.md` as project instructions; this bundle includes one. Paste the following as the initial task:

> Read `SETUP_WINDOWS10_OPENCODE.md` and `AGENTS.md` fully. Execute the documented Windows 10 + WSL2 DISCUS_SUITE and Python/Jupyter integration plan step by step. First perform read-only environment checks. Ask my approval before admin tasks, installs, reboots, or running downloaded installers. Prefer official DiffuseCode v.6.22.01 release assets/current Windows-WSL guide, but verify the actual available asset names. Reuse an existing WSL environment and existing installations. Run the project's Python tests, inspect and validate the generated DISCUS macro, then demonstrate one **real** 1D PDF job through `PDFRunner` and notebook 1. Verify notebook 2 only as a synthetic 3D-ΔPDF FFT example. Fix integration problems in code and add regression tests if required. Save all evidence, exact paths/commands, changes, and unresolved issues in `SETUP_REPORT.md`. Never claim end-to-end 3D-ΔPDF simulation is implemented. Stop and report blockers rather than inventing results.

## 13. References (official primary sources)

- DiffuseCode repository / installation overview: https://github.com/tproffen/DiffuseCode
- Releases (choose correct platform/install guide): https://github.com/tproffen/DiffuseCode/releases/latest
- v.6.22.01 release: https://github.com/tproffen/DiffuseCode/releases/tag/v.6.22.01
- v.6.22.00 release (HDF5 changes): https://github.com/tproffen/DiffuseCode/releases/tag/v.6.22.00
- Microsoft WSL install: https://learn.microsoft.com/en-us/windows/wsl/install
- Microsoft manual WSL install: https://learn.microsoft.com/en-us/windows/wsl/install-manual
- OpenCode project rules (`AGENTS.md`): https://docs.opencode.ai/docs/rules/
- OpenCode on Windows (WSL recommendation): https://docs.opencode.ai/docs/

**Verification status of this brief:** upstream GitHub README/release details, Microsoft WSL documentation and OpenCode instructions were checked online. The included Python starter's local files were inspected. DISCUS was **not** installed or tested on the user's Windows machine from this environment. Actual executable arguments, `.cll` compatibility, Windows release asset names and real output still require local agent validation.
