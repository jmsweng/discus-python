# Redactions before publication (2026-11-03)

Local absolute paths of the machines the runs were made on were replaced in files that are **not** part of any
hash registration (no MANIFEST, predictions or evidence hash lists them). Nothing else in these files changed.

| File | Replaced | With |
|---|---|---|
| `discus_6_23_01/correlated_run/run_corr.sh`, `run_corr_diffuse.sh`, `short_range_run/run_sr.sh` | `cd` to the absolute checkout path | `cd "$(dirname "$0")"` |
| `discus_6_23_01/correlated_run/corr_run.log`, `corr_diffuse_run.log`, `short_range_run/sr_run.log` | the Windows checkout path in DISCUS's `Start directory` line | `<checkout>` |
| `discusvolume_lattice/inventory.txt`, `compare_bases.log` | the local agent's worktree path | `<worktree>` |
| `discus_6_23_01/LOCAL_AGENT_HANDOFF.md` | the Windows checkout path | `/mnt/c/.../discus-python` |

Left unchanged on purpose: hash-registered evidence (run directories named `<UTC time>_<host name>` and the logs
inside them, `tests/fixtures/structure_ops/live_evidence/` with its `MANIFEST.sha256`), because editing it would
break the registrations that make it evidence.
