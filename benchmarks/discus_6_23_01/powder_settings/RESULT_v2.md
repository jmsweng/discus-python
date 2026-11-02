# Powder settings probe v2: result

- **Registered verdict: FAIL, 1 of 15 criteria.**
- Run `20261102T221050Z_LAPTOP-NGU3BT0N`, DISCUS_SUITE 6.23.01, binary SHA-256 `d4274f13…45706f3` (the same as Gates 0–3), registered commit `6c58688`.
- `probe_v2.py`, `probe.py` and all 8 macros match `predictions_v2.json`, and `MANIFEST.sha256` verifies.
- The single failure is a `set delta` prediction. **The F(Q)/G(r) model diagnosed from v1.1 passed every criterion on new settings and with neutrons.**

## Registered criteria (`evidence_v2/…/RESULT_TABLE.md`)

| Criterion | Cases | Result |
|---|---|---|
| P1: first written F point is the S = 0 placeholder | b, b_qmin, b_qmax, b_dq, n | PASS (≤ 2.6×10⁻³, limit 1×10⁻²) |
| P2: G(r) = Σ F sin(Qr) ΔQ with F from the second line, on the multiples of ΔQ starting at nint(Qmin/ΔQ)·ΔQ (as implemented in the registered analyzer; the docstring's "label − ΔQ" is exact only when Qmin is a multiple of ΔQ), zeroed from the last sign change, ramp Q·F₁/Q_label₁; no fit | same | PASS (1.4–2.2×10⁻⁴, limit 1×10⁻³); free scale fits 1.00000 |
| P3: F at the true Q independent of qmin / qmax | b_qmin, b_qmax | PASS (6.6×10⁻¹⁶, 5.4×10⁻¹⁶) |
| P4: `set delta` 0.03 aborts | d, d_omp1 (one thread) | PASS: `double free`, exit −6 |
| P4: `set delta` 0.03 aborts | d_grid (qmin 0.5, qmax 7.0, ΔQ 0.02) | **FAIL**: ran and wrote output |
| I1: ΔQ dependence of F (0.008 vs 0.012) | b_dq | 2.3×10⁻³ (informational) |

## What this establishes for DISCUS_SUITE 6.23.01 powder PDF (`set calc,debye`, `set axis,q`)

1. **F(Q) output labels are one ΔQ too high.** The value written at Q was computed at Q − ΔQ, and the first line is a placeholder (S = 0, F = −Q). Read DISCUS F(Q) files accordingly.
2. **G(r) = Σ F(Q) sin(Qr) ΔQ**, which is π/2 times the textbook (2/π)∫F sin(Qr) dQ of DISCUS's own F(Q). Other conventions also enter, such as how S(Q) is normalised, so this alone does not fix the absolute scale compared with other programs. But the factor is there.
3. **Termination:** F is set to zero from its last sign change up to Qmax. Below Qmin it is a straight line to F(0) = 0.
3a. **The transform puts F on the multiples of ΔQ** (first real point at nint(Qmin/ΔQ)·ΔQ, source: `iqmin = NINT(qmin/deltaq)` in `fft_fq`). F was computed at Qmin, Qmin + ΔQ, …. So when Qmin is not a multiple of ΔQ, DISCUS transforms F at Q values up to ΔQ/2 away from where it was computed (v2 base: 0.348 instead of 0.35). The registered analyzer implemented this placement. The independent check found that the docstring wording ("label − ΔQ") is 2×10⁻² off for these cases, and that the snapped placement gives 1.4–1.7×10⁻⁴.
4. **Qmin and Qmax** change only which part of F enters (F itself is unchanged to 10⁻¹⁵). ΔQ also changes F slightly: 0.2 % at 0.008 vs 0.012, up to 1.4 % at 0.05 vs 0.005 (v1.1). The cause is not identified.
5. **`set delta` > 0 is not usable.** It aborted on 5 of 6 grid/thread combinations (v1.1 and v2), with heap corruption in `powder_convolute`, and did not abort on one. The likely cause from the source is an out-of-bounds write one point past the end of the pattern (`powder_conv_res` runs to `xmax + xdel`). That would explain why the outcome depends on the grid, and it means a run that does not crash cannot be trusted either.

## Post hoc: the one non-crashing `set delta` run

Compared with the v1.1 `dq020` case, which has the same grid and gives the zero-width reference via P3, these numbers are not registered.
- `d_grid` is internally consistent: its G(r) is the P2 transform of its own F(Q) to 1.9×10⁻⁴.
- `set delta,0.03` changes G(r) by 2.8 %. F(Q) is smoothed: it is 8×10⁻³ from a Gaussian convolution (σ = 0.03/√2) of the zero-width F, and 4.8×10⁻² from the unsmoothed F.
- The predicted damping G₀·exp(−δ²r²/4) leaves 1.4 % of the 2.8 % change unexplained.
- So the width convention is only roughly confirmed. This is expected in part, because DISCUS convolves the intensity, not F, and adds a mirror term and an integral rescaling.

## Independent check

A separate agent recomputed claims 1–3 from the raw `run.mac` / `fq.dat` / `result.pdf` / `run.log` files with its own code, without reading the analyzers or these reports. Its results:
- **Confirmed:** placeholder line, one-ΔQ label offset (best free shift exactly ΔQ), F independent of Qmin/Qmax (≤ 7×10⁻¹⁶), ΔQ dependence ≈ 0.3·ΔQ, G = Σ F sin(Qr) ΔQ (π/2 vs textbook in every case), `set delta` aborts.
- **Model accuracy:** every case to ≤ 2.2×10⁻⁴ with the ΔQ-multiple placement. It also found the docstring wording issue in 3a.
- **Alternatives tested:** no label shift (6×10⁻² … 0.61), no zeroing (1.6×10⁻² … 0.26), no ramp below Qmin (4×10⁻³ … 0.36). All are clearly worse.
