# Gate 2 (gate2-prereg-v4) — result: **PASS**

Run `20261029T010305Z_LAPTOP-NGU3BT0N` · DISCUS_SUITE 6.23.01 (binary SHA-256 `d4274f13…45706f3`) · registered commit `87a382d`
· `predictions.json` SHA-256 `bebdc86a…acba311` · evaluator report: `evidence/20261029T010305Z_LAPTOP-NGU3BT0N/GATE2_RESULTS_PREREG.md`.

**0 of 37 registered criteria failed.** The reviewer independently re-ran `analyze.py` on the committed evidence, on a different machine: PASS, all manifest
hashes verified, and the values agree with the operator's to the 10⁻¹⁶ level (cross-machine floating point). The operator's report files were left as committed.

## What was predicted blind, and what the engine produced

The structures used fresh seeds (block 20261101, random 20261102), and the predictions were computed and committed before any engine run of them.

| Quantity | Block | Random (control) | Limit |
|---|---|---|---|
| Total intensity vs engine-precision oracle (C0) | 2.0×10⁻¹⁶ | 2.0×10⁻¹⁶ | 1×10⁻⁶ |
| Total vs independent DABAX-physics oracle | 8.0×10⁻⁴ | 8.0×10⁻⁴ | 2.4×10⁻³ (3 ε_tab) |
| `set aver,100` diffuse: 0 at Bragg; = total elsewhere (C1) | 0; 4.5×10⁻¹⁵ | 0; 6.8×10⁻¹⁵ | 1×10⁻⁹ |
| Pair correlations from one Brillouin zone, engine g (A) | 4.4×10⁻¹⁶ | 5.0×10⁻¹⁶ | 1×10⁻⁶ |
| 3D-ΔPDF vs engine-precision oracle (B1) | 3.8×10⁻¹⁵ | 5.6×10⁻¹⁵ | 1×10⁻⁶ |
| 3D-ΔPDF vs physics oracle (B1) | 5.1×10⁻⁴ | 5.1×10⁻⁴ | 1.5×10⁻³ |
| **max \|ρ/ρ₀ − prediction\| over predeclared τ (B2)** | **4.5×10⁻⁶** | 1.9×10⁻⁶ | 1×10⁻² |
| Required signs (B2) | 100 % | – | 100 % |
| Nearest-neighbour \|ρ/ρ₀\| of control (N) | – | 0.019 | 0.05 |
| Tie voxels equal a candidate lookup (R) | 216/216 | 216/216 | 100 % at 1×10⁻¹² |
| imag/real; Friedel (B3); logs (L, T) | pass | pass | |

Example, block (1,0,0): engine −0.5022, prediction −0.5022, realized pair sum C/N −0.5020.

## Rounding ties (supporting observation, post hoc)

The same 44 of the 216 tie voxels are rounded down by the engine in all four v3/v4 runs, whatever the structure. In the block model only 72 tie voxels have
distinguishable candidates, and its 14 round-downs are exactly the 44-set restricted to those. So the engine's tie choice depends on the grid only. It is
deterministic floating-point behaviour of DISCUS's s computation, not noise. Without the R treatment the v4 data would have given B1 = 1.7×10⁻⁶ (block)
and 2.8×10⁻⁵ (random) against the NINT-only emulation (`tie_rounding` in `gate2_results.json`).

## What this establishes

With the instrument configured as registered (B = 0, `set aver,100`, commensurate 16-cell grid, |H| ≤ 4):
- DISCUS_SUITE 6.23.01's total intensity is the textbook kinematic sum with Waasmaier–Kirfel f at s rounded to 0.001 Å⁻¹. It agrees to machine precision.
- Its Bragg-subtracted diffuse volume is exactly |F − ⟨F⟩·L|².
- That volume, transformed by `discuspdf.delta_pdf_fractional`, gives a 3D-ΔPDF whose off-origin peaks have the predicted positions, signs and amplitudes
  for a structure nobody had seen, and that is featureless for the random control.

## What it does not establish

- Absolute normalization, and the native `value 3DPDF` convention (Gate 3).
- Behaviour with B ≠ 0 under this design (v2 showed `temp use` applies B = 1 Å² although the files say 0; cause not identified), with multi-site cells,
  with non-commensurate grids, or with non-dyadic steps (Gate 1's 1.77 % residual on step 0.01 is attributed, not proven).
- Anything about other DISCUS versions or about `four_aver` with partial averaging (`set aver` < 100).

## History of this gate (all evidence kept)

v1: the structure files were rejected by the engine (−49). v2: FAIL, an unrequested Debye–Waller factor. v3: FAIL 2/33 on rounding ties. v4: PASS with
fresh blind seeds. Each deviation is recorded in `predictions.json` and in `PREREGISTRATION.md`. Each version change after v2 was made after seeing data, as stated there.
