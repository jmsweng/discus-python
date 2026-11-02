# Gate 0 v2.1 pre-registered experiment: RESULT

- Verdict: **PASS** (0/18 criteria failed)
- RUN_ID: `20261102T171021_LAPTOP-NGU3BT0N`
- Evidence dir: `benchmarks/discus_6_23_01/gate0_prereg/evidence/20261102T171021_LAPTOP-NGU3BT0N/`
- Registration commit: `3ffc4e6d6458a3f9c208b9b38a55079890883c45` (branch `science/gate0-prereg`, v2.1)
- predictions.json sha256: `7c45945c12e84ac89774e3b1785608e8d714e2924c681ddca78080af15774be4`
- DISCUS binary sha256: `d4274f13f938e8fda040a2a3636f1accfc18ec3c10ca8a1b37069ac4d45706f3` (`/usr/local/bin/discus_suite`, WSL2 Ubuntu 22.04.5 LTS, 16 cores)
- Engine: 256 runs (arms A/B x total/diffuse x kk=00-63), 0 non-zero exits; all 256 `.h5` + 256 `.log` present.
- Analyzer: `benchmarks/discus_6_23_01/gate0_prereg/gate0_analyze.py` (Python 3.10.12, numpy 2.2.6, h5py 3.16.0), run without `--cache`

## Criterion table (verbatim from GATE0_RESULTS_PREREG.md)

| criterion | value | limit | pass |
|---|---|---|---|
| E1 A: max total vs realization engine oracle (relL2) | 8.620e-15 | <= 1.000e-09 | PASS |
| E2 A: max diffuse vs realization engine oracle (relL2) | 7.251e-15 | <= 1.000e-06 | PASS |
| E3 A: max diffuse at integer h / max diffuse | 0.000e+00 | <= 1.000e-06 | PASS |
| E4 A total: log shows 'Temp. factors : ignored' | 1.000e+00 | >= 1.000e+00 | PASS |
| E4 A diffuse: log shows 'Temp. factors : ignored' | 1.000e+00 | >= 1.000e+00 | PASS |
| S1 A: |mean x - 1| | 1.851e-03 | <= 3.937e-02 | PASS |
| S2 A: |var x - 1| | 4.067e-03 | <= 1.114e-01 | PASS |
| S3 A: n_voxels outside Gamma interval (15561 tested) | 0.000e+00 | <= 0.000e+00 | PASS |
| S4 A: max |E_phys - E_eng| / E_eng | 1.038e-03 | <= 5.000e-03 | PASS |
| E1 B: max total vs realization engine oracle (relL2) | 6.727e-15 | <= 1.000e-09 | PASS |
| E2 B: max diffuse vs realization engine oracle (relL2) | 7.765e-15 | <= 1.000e-06 | PASS |
| E3 B: max diffuse at integer h / max diffuse | 0.000e+00 | <= 1.000e-06 | PASS |
| E4 B total: log shows 'Temp. factors : used' | 1.000e+00 | >= 1.000e+00 | PASS |
| E4 B diffuse: log shows 'Temp. factors : used' | 1.000e+00 | >= 1.000e+00 | PASS |
| S1 B: |mean x - 1| | 1.851e-03 | <= 3.937e-02 | PASS |
| S2 B: |var x - 1| | 4.067e-03 | <= 1.114e-01 | PASS |
| S3 B: n_voxels outside Gamma interval (15561 tested) | 0.000e+00 | <= 0.000e+00 | PASS |
| S4 B: max |E_phys - E_eng| / E_eng | 1.899e-03 | <= 5.000e-03 | PASS |

## Run notes

- WSL kept alive for the whole run (`wsl.exe -d Ubuntu-22.04 -- sleep 21600` in a separate background process); the 256 macros ran sequentially in ONE WSL call (bash loop), each exit code logged to `RUN_STATUS.txt` (256 x `exit=0`, ends `ALL_DONE`).
- `MANIFEST.sha256` lists only the files committed in the evidence dir (646 lines); `sha256sum -c MANIFEST.sha256` exit 0. The 256 `.h5` hashes are in `EVIDENCE_CONTENT_HASHES.json` only.

## Prior attempts

- **v1.5**, RUN_ID `20261102T063824_LAPTOP-NGU3BT0N`, registration `2132146`, evidence commit `321a123`: **FAIL 2/18** — only E4 B total/diffuse failed, because the registered regex looked for per-atom B values that DISCUS 6.23.01 never prints; all E1–E3 and S1–S4 passed. Diagnosis: `DIAGNOSIS_v1_5.md`. Re-registered as v2.1 with fresh seeds before any v2.1 engine run.
- Arms A and B use the same seeds, so their S1–S3 values are identical by construction; arm B adds engine-precision coverage with B ≠ 0 (E1/E2), not independent statistics.
