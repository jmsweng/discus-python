# Baseline audit (read-only) — 2026-10-22

Starting points: handoff baseline `main` = `b5f98d350ca8eb7aeb7775b979cfd2164b2fc2ed`.
During the session `main` advanced to `ef06490` (PR #4 merge of OpenCode cleanup
`abd52f0`); the code inspected here is that merged content.

## Observed
- `python -m pytest -q` (Linux VM, Python 3.10, numpy/h5py/scipy/matplotlib installed
  in a throwaway venv): **54 passed, 1 skipped** (the skip is the opt-in real-DISCUS test).
  These are unit/fake-binary tests only; they say nothing about a real 3D run.
- `discus_suite` is **not** present on the Linux VM available to Claude; no DISCUS
  command was executed. The WSL2 install in the handoff was not reachable.
- Working tree shows 15 tracked files as modified; `git diff --ignore-cr-at-eol` is empty,
  so this is CRLF-vs-LF only (Windows checkout), not a content change.

## Findings from reading the code/docs (not executed against DISCUS)
1. **1D reference/test mismatch (S4).** `examples/run_sample/config.json` and `run.mac`
   record `qmax=30`, but `README.md` and `tests/test_integration_discus.py` use `qmax=25.0`
   and compare to that committed `result.pdf` at relative L2 < 1e-6. If the reference was made
   with qmax=30, the opt-in integration test should fail against a real engine. Untested.
2. **Doc drift.** `docs/VERIFIED_DISCUS_MACRO.md` says `set qmax` is omitted when `qmax == 0`;
   `pdf.py` now rejects `qmax=0`. `examples/run_sample/run.mac` writes `outf result.csv`
   while the committed file is `result.pdf` (renamed in `59cdb7b`), and its `config.json`
   contains a `/home/user/...` structure path.
3. `Diffuse3DRunner` hard-codes the output name `diffuse.h5`; nothing shows the installed
   release writes that file (handoff §1).
4. `ReciprocalVolume` accepts caller-supplied `coordinate_system`, `units`, and
   `preprocessing` strings with no source-metadata check (handoff §2). It also requires an axis
   centered by `(arange(n) - n//2) * dq`, which a native hkl or skewed grid will not satisfy.
5. `delta_pdf_from_diffuse`: r-axis = `2π·fftfreq(n, dq)` (spacing `2π/(n·dq)`, limit `π/dq`) and
   normalization `dQx·dQy·dQz/(2π)^3` are internally consistent for a centered grid; a quantitative
   check against an analytic transform pair is still to be added.
6. `plot_delta_pdf_slice` passes `extent=(a[0], a[-1], b[0], b[-1])` — pixel *centers* — so
   displayed positions are shifted by half a voxel (`dr/2`) and the full span is under-drawn.

## Not tested / not claimed
Anything about real DISCUS 3D syntax, output format, Bragg handling, Q calibration, or
benchmark agreement. S1–S3 acceptance criteria are unmet; see `PROBE_PROTOCOL.md`.

## Added after reading the installed manual (`discus_man.pdf`, v6.23, pages cited by manual page)
Observed environment (WSL2, from `probe_environment.sh`): `/usr/local/bin/discus_suite`,
sha256 `d4274f13f938e8fda040a2a3636f1accfc18ec3c10ca8a1b37069ac4d45706f3`, banner
"D I S C U S - S U I T E Version 6.23.01", linked against `libhdf5_serial*` (HDF5 output is built in).
Headless `help` prints nothing, so the PDF manual (copied from `/usr/local/share/`) is the source.

From the manual (documentation claims, NOT yet confirmed by a run):
- p.31-32: 3D single-crystal Fourier = `four` / `ll,lr,ul,tl` corners / `na,no,nt`; the grid is in
  **hkl (rlu), not Å⁻¹**; use grid step 1/N_cells; `set aver,<per>` subtracts the average structure
  factor <F> (per=0 disables), `lots off` disables lots; with lots only intensity can be saved.
- p.43-44: outputs `stan, hdf5 (obsolete, Yell-compatible), nexus/common (recommended), vesta, mrc`;
  `value int|amp|real|imag|pha|3DPDF`. DISCUS can write its own 3D-PDF (a possible independent reference).
  `3DPDF` zeros intensity beyond an optional limit and truncates later output: write it last.
- p.50-53 (ch.4, "To be completed" in the manual): common/nexus HDF5 data lives under `/entry/data`
  and is meant to record unit-cell lengths/angles, `data_type_axes` ("hkl"/"Q"/...), `data_type_content`
  ("3d-delta-pdf", ...), `data_type_with_bragg` ("bragg_present"/"bragg_subtracted"). Exact dataset
  paths and whether these attributes are really written are unverified until a file is inspected.
- p.38-39: "The `set delta` command should no longer be used except with 0.00" — recommended is
  `set profile`. `PDFConfig.qsig` maps to `set delta`, so it is a deprecated parameter (S4).
- `replace <at1>,<at2>,all,<prob>` randomly substitutes atoms; keep VOIDs, and do not `purge`,
  before `set aver` (manual warns the average structure is then wrong).

## First real DISCUS 3D output (run 1, WSL2, 2026-10-22 23:30 UTC) — OBSERVED
Macro `C0.mac` (GaAs 8x8x8, perfect crystal, `set aver,0.0`, `lots off`, turbo, 65^3 hkl grid -4..4).
Wall time 0.66 s, max RSS 54 MB. `form nexus` + `value int` wrote `C0_total_nexus.h5`
(2,221,944 B, sha256 `ea744bea3680bae2049698c1cf4e16c81ff2e2da7c54f898b24eda372ba5a253`).
- Layout (inspector JSON): everything under `/entry/data`; intensity = `/entry/data/data_values`,
  shape (65,65,65), float64, no NaN/Inf, min 1.5e-25, max 1.716e10 at index [32,32,32] (= hkl 000,
  consistent with corner -4 and step 0.125 → center index n//2).
- Metadata actually present: `coordinate_unit = "basecell_fractional"`, `data_type_axes = "hkl"`,
  `data_type_with_bragg = "bragg_present"`, `data_type_content = "intensity"`, `data_type_experiment =
  "calculated"`, `data_corner = (-4,-4,-4)`, `data_dimension = (65,65,65)`, `data_increment_vector`
  (3x3, 0.125 on the diagonal), `data_axes = (1,2,3)`, `unit_cell_lengths/angles`, space-group fields.
  `h5py` attrs are empty; all metadata are datasets. There are NO Cartesian-Q axes: Q must be derived.
- NOT determined (isotropic cubic grid cannot tell): axis order of `data_values` vs. `data_axes`,
  and whether `data_increment_vector` rows or columns are the increments. Needs an anisotropic-grid
  (e.g. na != no != nt) and a non-orthogonal-cell run.
- `form hdf5` (documented as obsolete) failed in this install: "I/O Error opening file" (cause unknown).
- My template bug (disorder lines injected into a comment line) made C1/C2 fail at macro line 2;
  template fixed, not yet re-run. C0's later passes (Bragg-subtracted, 3DPDF) did not execute.
