# OpenCode project instructions — DISCUS PDF automation

Read `SETUP_WINDOWS10_OPENCODE.md` before any system change. Implement its staged Windows 10 / Ubuntu WSL2 installation and verification workflow.

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

- `SETUP_REPORT.md`: exact commands, versions, executables, tests, output paths and blockers.
- `docs/VERIFIED_DISCUS_MACRO.md`: verified engine CLI and 1D PDF macro details.
- At least one real job through `PDFRunner` with nonempty, validated numeric output when installation succeeds.
- Document 3D-ΔPDF notebook as synthetic demonstration, not end-to-end simulation.

If prerequisite access or required software is unavailable, report the blocker rather than bypassing permissions or inventing a result.
