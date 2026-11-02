# Gate 0 v2 pre-registered experiment: RESULT

- Verdict: **FAIL** (2/18 criteria failed)
- RUN_ID: `20261102T063824_LAPTOP-NGU3BT0N`
- Evidence dir: `benchmarks/discus_6_23_01/gate0_prereg/evidence/20261102T063824_LAPTOP-NGU3BT0N/`
- Registration commit: `2132146c6b53f79f0b2ffb6dfa51bd0d005e814e` (branch `science/gate0-prereg`, v1.5)
- predictions.json sha256: `72129b03b6293afd60356ba5eb8165a88b07f9e7ef31c621ede55911cbd1ee3c`
- DISCUS binary sha256: `d4274f13f938e8fda040a2a3636f1accfc18ec3c10ca8a1b37069ac4d45706f3` (`/usr/local/bin/discus_suite`, WSL2 Ubuntu 22.04.5 LTS, 16 cores)
- Engine: 256 runs (arms A/B x total/diffuse x kk=00-63), 0 non-zero exits; all 256 `.h5` + 256 `.log` present.
- Analyzer: `benchmarks/discus_6_23_01/gate0_prereg/gate0_analyze.py` (Python 3.10.12, numpy 2.2.6, h5py 3.16.0)

## Criterion table (verbatim from GATE0_RESULTS_PREREG.md)

| criterion | value | limit | pass |
|---|---|---|---|
| E1 A: max total vs realization engine oracle (relL2) | 8.502e-15 | <= 1.000e-09 | PASS |
| E2 A: max diffuse vs realization engine oracle (relL2) | 7.833e-15 | <= 1.000e-06 | PASS |
| E3 A: max diffuse at integer h / max diffuse | 0.000e+00 | <= 1.000e-06 | PASS |
| E4 A total: log shows 'temp ignore' | 1.000e+00 | >= 1.000e+00 | PASS |
| E4 A diffuse: log shows 'temp ignore' | 1.000e+00 | >= 1.000e+00 | PASS |
| S1 A: |mean x - 1| | 7.115e-03 | <= 3.937e-02 | PASS |
| S2 A: |var x - 1| | 1.501e-02 | <= 1.114e-01 | PASS |
| S3 A: n_voxels outside Gamma interval (15561 tested) | 0.000e+00 | <= 0.000e+00 | PASS |
| S4 A: max |E_phys - E_eng| / E_eng | 1.038e-03 | <= 5.000e-03 | PASS |
| E1 B: max total vs realization engine oracle (relL2) | 6.662e-15 | <= 1.000e-09 | PASS |
| E2 B: max diffuse vs realization engine oracle (relL2) | 8.471e-15 | <= 1.000e-06 | PASS |
| E3 B: max diffuse at integer h / max diffuse | 0.000e+00 | <= 1.000e-06 | PASS |
| E4 B total: log shows B values | 0.000e+00 | >= 1.000e+00 | **FAIL** |
| E4 B diffuse: log shows B values | 0.000e+00 | >= 1.000e+00 | **FAIL** |
| S1 B: |mean x - 1| | 7.115e-03 | <= 3.937e-02 | PASS |
| S2 B: |var x - 1| | 1.501e-02 | <= 1.114e-01 | PASS |
| S3 B: n_voxels outside Gamma interval (15561 tested) | 0.000e+00 | <= 0.000e+00 | PASS |
| S4 B: max |E_phys - E_eng| / E_eng | 1.899e-03 | <= 5.000e-03 | PASS |

## Failed criteria

- E4 B total: log shows B values — value 0.000e+00 (limit >= 1.000e+00)
- E4 B diffuse: log shows B values — value 0.000e+00 (limit >= 1.000e+00)

No interpretation is offered; per the pre-registered stop rules this is a reportable outcome.
