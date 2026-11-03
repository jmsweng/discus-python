# Builder-backend probe (DISCUS_SUITE 6.23.01)

Checks three things the disorder builder needs and that had never been run in DISCUS: own neutron scattering
lengths (`scat`) in Fourier and powder runs, `VOID` padding atoms, and the powder PDF from a `.stru`
(`read stru`). Criteria: [PREREG.md](PREREG.md), fixed before the run.

| File | Purpose |
|---|---|
| `make_inputs.py` | writes `inputs/<run>/` and `runs.json` (committed; do not edit, do not re-run after registration) |
| `criteria.json` | pass/fail limits fixed before the run |
| `run_probe.sh` | WSL runner: `bash -l run_probe.sh <inputs dir> <output dir>` |
| `analyze.py` | `python analyze.py <output dir> [--json out.json]`; exit 1 if a criterion fails; `--selftest` |
| `evidence/` | engine outputs of the recorded run (added after the run) |
| `RESULT.md` | the outcome (added after the run) |

## Running it (owner's machine, WSL Ubuntu-22.04)

From a checkout of the branch (Windows paths in the agent handoff):

```bash
cd <checkout>/benchmarks/discus_6_23_01/e_builder_backend
RUN=~/e_probe/$(date -u +%Y%m%dT%H%M%SZ)_$(hostname)
bash -l run_probe.sh inputs "$RUN"        # 12 runs, a few seconds to minutes
python3 analyze.py "$RUN" --json "$RUN/analysis.json"
```

`exit_codes.txt` must list 12 runs with `exit=0`. Copy the whole `$RUN` directory to `evidence/<RUN name>/`
(the `.h5` and `result.pdf` files are small: 17³ voxels and 501 points).
