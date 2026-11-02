# Package D reference-scattering probe (DISCUS_SUITE 6.23.01)

Compares `discuspdf.reference_scattering` with real DISCUS intensities. Criteria: [PREREG.md](PREREG.md).

| File | Purpose |
|---|---|
| `make_inputs.py` | writes `inputs/*.stru`, `inputs/*.mac`, `runs.json` (committed; do not edit) |
| `criteria.json` | pass/fail limits fixed before the run |
| `run_probe.sh` | WSL runner: `bash -l run_probe.sh <inputs> <out>` → `<out>/work/*.h5`, logs, hashes |
| `analyze.py` | `python analyze.py <out>` prints the comparison and exits 1 if a criterion fails |
| `evidence/` | committed engine outputs of the recorded run |
| `RESULT.md` | the outcome |
