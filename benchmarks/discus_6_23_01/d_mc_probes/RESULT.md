# D MMC probe — result (DISCUS_SUITE 6.23.01)

**10 of 13 predictions held; 2 informational; 2 were wrong as written and are explained below — both point
to how DISCUS's MMC actually works, which the Python API now documents and guards.**

| | |
|---|---|
| Registration | inputs, macros, `PREDICTIONS.md` in `fff6e69`, before the run |
| Engine run | owner's laptop, WSL Ubuntu-22.04, `discus_suite` SHA-256 `d4274f13…45706f3`, 16 logical CPUs, 2026-11-01 22:47 UTC, local agent task 002 |
| Evidence | `evidence/run1/` (11 logs, saved structures, exit codes, hashes); every run exited 0 |
| Analysis | `python benchmarks/discus_6_23_01/d_mc_probes/analyze.py` → `evidence/run1/analysis.json` |

| # | Prediction | Outcome |
|---|---|---|
| P1 | same 8-value seed, serial → identical | **PASS**: saved structures byte-identical (SHA-256 `20a42b46…`), identical tables |
| P2 | 4 OpenMP threads, same seed | INFO: **not** reproducible (different outputs, achieved −0.335 / −0.334) and **45 056 moves instead of 40 960** (a serial tail is added) |
| P3 | swchem conserves composition | **PASS** (2048 Cu, 2048 Au in every occupational run) |
| P4 | `set cyc` = trial moves | **PASS**: final `Gen: 40960` |
| P5 | printed `corr` achieved = Package C α on the saved structure | **PASS**: −0.306 vs −0.306315; 24 576 pairs both |
| P6 | CORR feedback reaches −0.3 | INFO: −0.306 after 10 moves/atom (feedback trajectory −0.27, −0.35, −0.37, … −0.30) |
| P7 | J = 0: every move neutral, accepted with ½ | **PASS**: 10 001 of 20 000, all neutral |
| P8 | J = 1, kT = 10⁴: accepted fraction ½ | **FAIL as written**: 0.639. The prediction forgot that **downhill moves are always accepted** in the source rule; the moves that are not downhill were accepted with 0.491 ≈ ½, exactly the heat-bath `1/(1+e^{ΔE/kT})` at ΔE/kT → 0. Metropolis would accept ≈ 100 %. The source reading (`mmc_basic.f90:152-198`) is confirmed. |
| P9 | J = +0.5 orders | **PASS**: −0.658 |
| P10 | `shift` without `set move` does nothing | **PASS**: positions unchanged (default σ = 0, not the help's 0.2) |
| P11 | spring / LJ with σ = 0.01 cell narrows the NN distances | **FAIL**: the spread grew 0.072 → 0.096 (spring) / 0.094 (LJ) Å. With k = 1 kT/Å² and moves of 0.03 Å, ΔE ≈ 10⁻³ kT, so almost every move is accepted with ≈ ½ or 1 and the atoms random-walk. DISCUS also rounds ΔE to 10⁻⁴ before the test (|ΔE| < 5·10⁻⁵ counts as neutral). Energy scales must be chosen so that typical ΔE ≫ 10⁻⁴ and ≳ kT. |
| P12 | printed `cd` achieved = Package C uncentred correlation from the site averages | **PASS**: 0.052 vs 0.052257 (first analysis used the displaced input as reference — an analysis error, fixed before RESULT) |
| P13 | exit 0 | **PASS** |

## Further observations used by the API
* `save` from the main level writes **no `ncell` line** and **upper-case atom names** (`CU`, `AU`); positions of
  shifted atoms are not re-wrapped. `MCRunner.load` reads the output with the input's supercell and restores
  the input's spelling of the species.
* Spring / Lennard-Jones tables print the number of pairs counting (a,b) and (b,a) separately: 20 736 for the
  6-vector neighbourhood of 1 728 Cu atoms (10 368 directed pairs).
* The "Initial" table is printed again after the run with the initial values; the final values follow
  "--- Final multiple energy configuration ---". The final `try:` count is `cyc + 1`.
