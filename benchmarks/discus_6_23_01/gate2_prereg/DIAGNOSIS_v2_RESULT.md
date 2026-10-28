# Gate 2 (gate2-prereg-v2) — written diagnosis of the FAIL

Run `20261028T234015Z_LAPTOP-NGU3BT0N`, engine DISCUS_SUITE 6.23.01 (binary SHA-256 `d4274f13…45706f3`).
Registered result: **FAIL, 16 of 29 criteria** (`evidence/…/GATE2_RESULTS_PREREG.md`, unchanged). Stop rule 1 applies
(C0 failed), so this is a diagnosis. It is not a re-analysis, and nothing registered was edited.

## Finding 1 — the engine applied an unrequested Debye–Waller factor with B = 1.000 Å² (cause of C0, A, B1)

The structure files carry `Biso = 0.0` and the macros say `temp use`. The engine output is nevertheless
reduced by exp(−2·B·s²) with B = 1 relative to the B = 0 oracle:

- Fit of ln(I_engine / I_oracle) against s² over all voxels: slope −2.000012 (total) and −2.000012 (diffuse).
- Oracle × exp(−2 s²), with s rounded to 0.001 Å⁻¹ (NINT, as the engine does), agrees with the engine **total** to
  3.6×10⁻⁷ (block) and 2.7×10⁻⁷ (random), against the registered engine-precision limit of 1×10⁻⁶.
- The best-fit B is 1.00000000 for the total and 1.0000003 for the diffuse. The engine's DW is B = 1 to eight digits.
- C0's relL2 is **identical for block and random (0.3001)** because the total is dominated by the Bragg peaks,
  which do not depend on the arrangement.
  This signature points to a global, arrangement-independent factor and not to a problem with the order model.
- The ratio engine/oracle at the 8 Bragg points falls from 0.37 to 0.22 with |H|, as exp(−2 s²) predicts.

Why B = 1: not established. DISCUS apparently does not take the `Biso` column of a `format numbers` atom line
as the Debye–Waller value for `temp use`, or it treats 0.0 as unset and uses a default of 1.0.
The DISCUS-written C1 file has Biso = 1.0, so the C1 runs could not show this. Cause to be identified in the
DISCUS source or by a probe (see the v3 plan).

Not a defect of the oracle's physics, and not of the pipeline: see Finding 2.

## Finding 2 — exploratory, post hoc: with the B = 1 factor divided out, the registered predictions are reproduced

This was **done after seeing the v2 results. It is not a pass.** Engine diffuse (block) divided by exp(−2 s_NINT²),
then run through the unmodified `delta_pdf_fractional`:

| τ (cells) | registered prediction | engine, DW removed |
|---|---|---|
| (1,0,0) | −0.4886 | −0.4886 |
| (0,1,0) | +0.5503 | +0.5503 |
| (0,0,1) | −0.0242 | −0.0242 |
| (1,1,0) | −0.2636 | −0.2636 |
| (1,−1,0) | −0.2665 | −0.2665 |
| (0,2,0) | +0.1008 | +0.1008 |

- B1 against the physics oracle: 5.1×10⁻⁴ (limit 1.5×10⁻³ — would pass).
- B1 against the engine-precision oracle: 4.3×10⁻⁵ (limit 1×10⁻⁶ — would **fail**).
- **Sign and position are correct.** Amplitudes agree to four digits at every predeclared τ. The registered B2 criterion (0.01) had already passed
  in the v2 run, since the DW is a smooth factor that mostly cancels in ratios.

## Finding 3 — an unexplained small diffuse residual (open)

After removing the DW, the diffuse volume differs from the engine-precision oracle by 4.1×10⁻⁵ (block) and 4.8×10⁻⁵ (random),
about 40× the registered 1×10⁻⁶. The worst voxels are of order 0.3 % at |I| ~ 10⁷. The total agrees to 3×10⁻⁷, and the
dyadic Gate 1 run agreed to 7×10⁻¹² with the same machinery, so DISCUS's phase/exp-table arithmetic is not the obvious suspect.
The fit excludes a different B (best-fit B = 1.0000003). Not diagnosed. It is below the physics tolerance, and no claim is made about its origin.

## Finding 4 — evaluator defect in the L criteria (corrected 2026-10-28)

*The first version of this note blamed the atom-count regex. That was wrong: the regex tolerates the padding.* The cause is
the error-line check, which matched `***` and so fired on DISCUS's start-up banner (a row of bare asterisks). The engine
logs are clean: counts are 2048 Cu and 2048 Au, and the only DISCUS error markers look like `***APPL***`, which are absent.
The v2 selftest did not catch it because its synthetic logs had no banner. v3 adds the banner to the synthetic logs.
It does not affect any physical criterion.

## Consequences for the gate status

- **Gate 2: FAILED as registered (v2); not passed.** The cause is an unrequested engine setting, not a demonstrated
  failure of DISCUS's off-origin ΔPDF physics. The exploratory result is strong evidence for the latter but is not a gate pass.
- **Gate 1 (dyadic C1): PASSED** — total 7.7×10⁻¹⁵, diffuse 6.9×10⁻¹², archived-macro rerun bit-identical (see `gate1_c1_dyadic/evidence/…`).
  This supersedes the "C1 diffuse residual 1.77 %" item: it was not intrinsic to the engine on a dyadic grid.

## Proposed gate2-prereg-v3 (needs owner approval before any commit of registered files)

1. Macros: add `temp ignore` so that B = 0 is what the engine actually uses, as designed. A probe must first confirm the
   keyword is accepted and the log's `Temp. factors` line changes from `used`.
   Fallback if it is not accepted: keep `temp use` and set the oracle to B = 1 explicitly (the engine behaviour is now measured).
2. Evaluator: error check matches `***WORD***` markers, not the banner; new T criterion on the `Temp. factors` line; selftest gains the banner, an unrequested-DW case and a temp-used-log case.
3. Everything else unchanged: seeds, grid, predictions, thresholds. In particular the 1×10⁻⁶ engine-precision limits stay, so Finding 3 will show up honestly if it persists.
4. Registered in the record as a deviation made **after seeing v2 results**; v2 results and this diagnosis are kept alongside.

## Status: v3 built (2026-10-28, owner approved)

`predictions.json` models, design and thresholds are identical to v2 (checked); only the four data macros, `generate.py` and
`analyze.py` (and their hashes) changed. v3 `predictions.json` SHA-256 `40dd2b14…31cb488d`. See `LOCAL_AGENT_HANDOFF.md` ("Registration v3").
