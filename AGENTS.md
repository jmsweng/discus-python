# OpenCode project instructions — DISCUS PDF automation

Read the tracked docs before working: `README.md`, `docs/VERIFIED_DISCUS_MACRO.md`, `docs/3D_INTEGRATION.md`, and `docs/INTEGRATION_HISTORY.md`. DISCUS_SUITE 6.23.01 is already installed and verified (see `docs/VERIFIED_DISCUS_MACRO.md`); do **not** reinstall it by default. Local-only setup logs (`SETUP_WINDOWS10_OPENCODE.md`, `SETUP_REPORT.md`) are gitignored and optional: consult them only if present on the machine.

## Working conventions

- Keep DISCUS and the Python virtual environment in the **same** OS context; WSL2 Ubuntu is preferred.
- Only use official DiffuseCode release assets; verify the asset/guide for the selected version.
- Ask human approval before elevated permissions, installations, reboots, executing unknown installers, or modifying pre-existing files.
- Never claim a real scientific simulation passed based only on Python unit tests or synthetic volume plots.
- Do not break scientific correctness: validate `.cll` parsing, command syntax, axes, output columns, units, and FFT assumptions.
- Preserve a separate, reproducible working directory for each simulation; keep macros, parameters, logs and output data.

## Project commands (Ubuntu WSL2)

```bash
cd ~/projects/discuspdf_starter
source .venv/bin/activate
python -m pytest -q
jupyter lab --no-browser --ip=127.0.0.1 --port=8888
```

## Deliverables

- `docs/VERIFIED_DISCUS_MACRO.md`: verified engine CLI and 1D PDF macro details (tracked in Git).
- At least one real job through `PDFRunner` with nonempty, validated numeric output (recorded in `examples/run_sample/`).
- Document 3D-ΔPDF notebook as synthetic demonstration, not end-to-end simulation.
- Local-only setup reports (e.g. `SETUP_REPORT.md`) are optional and gitignored; keep them out of Git.

If prerequisite access or required software is unavailable, report the blocker rather than bypassing permissions or inventing a result.

## 1D PDF output is `result.csv` — do not change it back to `.pdf`

This is the repository owner's decision (2026-11-02). It has been reverted once already
(merge `f8b2ba1`, 2026-10-22) on the mistaken belief that DISCUS can only write `result.pdf`.

- DISCUS writes whatever file name `outf` gives it. The validated macro keeps
  `outf result.pdf` (DISCUS's raw space-delimited text) only because archived runs and probe
  evidence use that exact macro.
- `PDFRunner` converts the raw file to **`result.csv`** (`r (Angstrom),G(r)` header,
  comma-separated). `PDFRun.output` and `PDFRun.load()` use the CSV; `PDFRun.raw_output` is the
  raw file. The name is `discuspdf.pdf.PDF_CSV_NAME`, guarded by `tests/test_pdf_csv_output.py`.
- Do not rename the output, the constant or the test expectations to `.pdf`. If a change seems to
  need that, ask the owner first.
- `examples/run_sample/` is an archived run from before this change and correctly has only
  `result.pdf`; leave it as it is.
