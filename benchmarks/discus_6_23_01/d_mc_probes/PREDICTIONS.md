# D MMC probe — questions and predictions (fixed before the engine run)

Source reading of DISCUS 6.23.01 (`discus/prog/mmc*.f90`, summarised in `docs/MC_FACTS.md`)
gives the predictions below. The probe records engine behaviour; a wrong prediction is a finding about
DISCUS or about our reading, and decides what the Python API may offer.

Structures (`make_inputs.py`): `cuau.stru` simple cubic a = 3 Å, 16³ cells, 2048 Cu + 2048 Au at random;
`cu_disp.stru` simple cubic Cu, 12³ cells, Gaussian displacements σ = 0.05 Å. All runs: `seed` with 8 values,
neighbourhood 1 = six vectors ±x, ±y, ±z (site 1 → 1), `save` from the main level after `exit` from mmc.

| # | Run(s) | Question | Prediction (source) |
|---|---|---|---|
| P1 | r1, r1b | Same 8-value seed, `set parallel, useomp:off`: identical result? | saved structures byte-identical; logs identical apart from timing lines |
| P2 | r1c, r1d | 4 OpenMP threads, same seed: identical? | **not** guaranteed: threads share the crystal and RNG (`mmc.f90:3736-3771`); recorded, no pass/fail |
| P3 | r1, r2, r3 | swchem conserves composition? | Cu and Au counts unchanged in every saved structure |
| P4 | r1 | `set cyc, N` = number of valid trial moves | final `Gen:` = 40960 |
| P5 | r1 | printed final achieved `corr` vs Package C on the saved structure (binary, ±d vectors → symmetric counts) | equal to the printed 3 decimals: achieved = Warren–Cowley α_CuAu = (p_CuCu − θ²)/(θ(1−θ)) over the 6-vector pair set; number of pairs = 6 × 4096 |
| P6 | r1 | CORR feedback reaches the target? | achieved within 0.05 of −0.3 (informational) |
| P7 | r2 | J = 0 (ENER): every ΔE = 0 → acceptance 1/(1+e⁰) = ½ (`mmc_basic.f90:152-198`) | accepted (g+n+b)/Gen = 0.50 ± 0.02, all in the `n` (neutral) column |
| P8 | r3 | J = 1, kT = 10⁴: heat-bath gives ≈ ½ for every move; Metropolis would give ≈ 1 | accepted fraction 0.50 ± 0.02 |
| P9 | r3b | J = +0.5 (ENER, kT = 1): positive J favours unlike neighbours | final achieved < 0 (ordering) |
| P10 | r4b | `shift` without `set move`: default σ = 0 (`discus_allocate_appl_mod.f90:1897`) | positions unchanged (to the 6 written decimals) |
| P11 | r4, r5 | `shift` with σ = 0.01 cell and spring / Lennard-Jones d0 = 3 Å | atoms move; positions change; spread of NN distances around 3 Å (printed sigma) decreases from its initial value |
| P12 | r6 | `swdisp` + `cd` CORR 0.3 along [100]: printed achieved vs Package C `displacement_correlation(reference="site_average", centered=False, direction=(1,0,0), fractional)` on the ±x pairs | equal to 3 decimals when Package C uses all six directed vectors ±x, ±y, ±z (the whole neighbourhood), x-components at both ends |
| P13 | all | engine exit code 0 except where DISCUS is known to abort at exit | exit 0 |
