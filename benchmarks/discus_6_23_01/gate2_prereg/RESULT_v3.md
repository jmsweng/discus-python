# Gate 2 (gate2-prereg-v3) — result and interpretation

Run `20261028T235829Z_LAPTOP-NGU3BT0N`, DISCUS_SUITE 6.23.01 (binary SHA-256 `d4274f13…45706f3`), registered commit `a243b4e`.
Evidence and the evaluator's own report: `evidence/20261028T235829Z_LAPTOP-NGU3BT0N/GATE2_RESULTS_PREREG.md` (unchanged).

## Registered verdict: **FAIL, 2 of 33 criteria** (one criterion, in both models)

Failed: *B1, pipeline 3D-ΔPDF of engine diffuse vs oracle at engine precision* — block 2.0×10⁻⁵, random 2.3×10⁻⁵, limit 1×10⁻⁶.
This verdict stands. No threshold was changed and nothing was rerun.

## Everything else passed, for both models

| Criterion | Block | Random | Limit |
|---|---|---|---|
| C0 total vs engine-precision oracle | 1.8×10⁻⁷ | 1.4×10⁻⁷ | 1×10⁻⁶ |
| C0 total vs physics (DABAX) oracle | 8.0×10⁻⁴ | 8.0×10⁻⁴ | 2.4×10⁻³ |
| C1 diffuse = 0 at Bragg positions; diffuse = total elsewhere | 0; 4×10⁻¹⁵ | 0; 7×10⁻¹⁵ | 1×10⁻⁹ |
| A correlations recovered (physics g / engine g) | 9.8×10⁻⁵ / 4×10⁻¹⁶ | 8.1×10⁻⁵ / 6×10⁻¹⁶ | 1×10⁻³ / 1×10⁻⁶ |
| B1 ΔPDF vs physics oracle | 5.1×10⁻⁴ | 5.1×10⁻⁴ | 1.5×10⁻³ |
| **B2 ρ/ρ₀ at every predeclared τ vs registered prediction** | **3.6×10⁻⁶** | 2.3×10⁻⁶ | 1×10⁻² |
| B2 signs where required | 100 % | 100 % | 100 % |
| B3 imag/real; Friedel | 7×10⁻¹⁶; 1×10⁻¹⁵ | 1×10⁻¹⁵; 2×10⁻¹⁵ | 1×10⁻¹⁰; 1×10⁻¹² |
| N negative control, max nearest-neighbour \|ρ/ρ₀\| | – | 0.015 | 0.05 |
| L, T (logs clean; `Temp. factors : ignored`) | pass | pass | required |

Block-model peaks, engine vs registered prediction: (1,0,0) −0.4886 / −0.4886, (0,1,0) +0.5503, (0,0,1) −0.0242, (1,1,0) −0.2636, (1,−1,0) −0.2665, (0,2,0) +0.1008.

## Why the one criterion fails (post hoc analysis, read-only)

Engine diffuse ÷ engine-precision oracle, voxel by voxel:
- Median ratio 1.00000000; in all but a few voxels the two agree to ~10⁻¹⁴.
- **14 voxels (block) and 44 (random) out of ~2×10⁶ differ, by 0.17–0.21 % each. Every one lies exactly on a rounding tie, s/(0.001 Å⁻¹) = n + 0.5.**
- 72 (block) and 216 (random) voxels sit on such ties; the engine and the oracle agree on the rest and disagree on these.
- The A criterion with engine g (one Brillouin zone only) is exact to 10⁻¹⁶, so the arrangement-dependent content is reproduced exactly.

Interpretation: the oracle's "engine-precision" mode emulates DISCUS's `NINT(s/0.001)` lookup. At exact ties the two programs
round some voxels to different neighbouring table entries (floating-point noise in s² and sqrt, or a different tie convention).
This is an emulation detail at the 10⁻⁵ level, not a physical disagreement. **I have not shown which program's rounding is "right", and
I do not claim the engine tie-breaking is understood.** The v2 residual of 4×10⁻⁵ is most likely the same effect (same size, same location pattern), not checked voxel by voxel.

## Consequences

- **Gate 2 as registered: FAILED (strict).** The only failing criterion tests agreement with a rounding emulation to 10⁻⁶.
- **The scientific question Gate 2 was built to answer is answered: yes.** With the instrument configured as designed (B = 0), DISCUS 6.23.01
  produces a Bragg-subtracted diffuse volume whose pipeline 3D-ΔPDF has the predicted off-origin peaks at the right positions, with the right
  signs and with amplitudes within 4×10⁻⁶ of the independent prediction, and nothing at nearest neighbours for the random control.
- I recommend **not** registering a v4 to tune the oracle's tie-breaking to the engine: that would make the oracle less independent to rescue a 10⁻⁶ limit that carries no physics.
  If the owner wants a clean registered PASS, the honest route is a v4 whose engine-precision B1 limit is stated in advance for tie voxels, labelled as a post hoc change. Not done.
- Not covered (Gate 3): absolute normalization and the native `value 3DPDF` convention.
