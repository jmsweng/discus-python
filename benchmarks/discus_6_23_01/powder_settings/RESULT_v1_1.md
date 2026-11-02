# Powder settings probe v1.1: result

- **Registered verdict: FAIL.** Run `20261102T220205Z_LAPTOP-NGU3BT0N`, DISCUS_SUITE 6.23.01, binary SHA-256 `d4274f13…45706f3` (the same as Gates 0–3), registered commit `67bfd62`. Evidence: `evidence/20261102T220205Z_LAPTOP-NGU3BT0N/` (`MANIFEST.sha256` verified).
- The engine ran 10 of the 13 cases cleanly. **All three `set delta` cases aborted.**
- For the 10 that ran, two of the registered predictions about the output (T1, T3) were wrong. T4 passed.
- Everything below "Diagnosis" was worked out **after seeing the data**. It is confirmed only by the v2 registration (`README.md`, "Registration history").

## Registered criteria (`RESULT_TABLE.md`, `results.json`)

| Criterion | Outcome |
|---|---|
| R: `set delta` 0.02 / 0.05 / 0.10 produce G(r) | **FAIL**: exit −6 (SIGABRT), no output, all three |
| T1: G(r) = (2/π) Σ F sin(Qr) ΔQ of the run's F(Q), best placement | **FAIL**: off by 36–56 %; best-fit scale ≈ 0.43–0.61 |
| T3: F(Q) at common Q independent of qmin / dq / qmax | **FAIL** for qmin and dq (10⁻² … 0.47); PASS for qmax (≤ 7×10⁻¹⁶) |
| T4: the F(Q) output block leaves G(r) unchanged | PASS, 4.4×10⁻¹⁶ |

To make the analyzer report engine failures instead of crashing on missing files, one change was made after the run (`analyze` only). It turns a case without `result.pdf` into an "R … FAIL" row. Criteria and limits are unchanged.

## Diagnosis (post hoc)

1. **`set delta` > 0 aborts DISCUS 6.23.01** with `double free or corruption (!prev)`, inside `powder_write_mod::powder_convolute` (backtrace in `delta*/run.log`). This happens for every width tried. Likely cause from the source (not confirmed): `powder_conv_res` runs to `xmax + xdel`, one point past the end of `pow_conv(0:num)`, and writes `dat(i)` there.
2. **F(Q) is written one Q step too high.** The value written at label Q was computed at Q − ΔQ, and the first written point is a placeholder with S = 0 (F = −Q). The fitted shifts between DQ cases are 0.005 / 0.015 / 0.045, which is ΔQ − 0.005 each time. With the labels corrected and the placeholder dropped, F(Q) is identical across Qmin and Qmax changes to ≤ 7×10⁻¹⁶. So the apparent Qmin effect in T3 was the placeholder alone.
3. **DISCUS's G(r) is Σ F(Q) sin(Qr) ΔQ, without the textbook factor 2/π.** It is π/2 = 1.5708 times (2/π)∫F sin(Qr) dQ of its own F(Q). This matches the source: `fft_fq` multiplies by 2/π·ΔQ, and the caller then multiplies by π/2. Other conventions also enter, such as how S(Q) is normalised, so this does not by itself say how DISCUS's G(r) compares with other programs. But users comparing DISCUS G(r) with another program should expect this factor.
4. **Termination:** F(Q) is set to zero from its last sign change up to Qmax. Below Qmin, F is a straight line to 0 with slope F₁ / (label Q of the first real point), i.e. using the shifted label. With items 2–4 the model reproduces **all 10 runs to ≤ 2.2×10⁻⁴ (max rel), with no fitted parameter**. Without the zeroing step the mismatch is 15–18 %.
5. **DQ changes F(Q) itself slightly**, beyond the label shift: 1.6×10⁻³ (0.01 vs 0.005), 4.8×10⁻³ (0.02), 1.4×10⁻² (0.05). The cause has not been identified.

## Measured effect of each setting on G(r) (max |G − G_base| / max |G_base|, base qmin 0.5, qmax 7.5, dq 0.01)

| Setting | Values | Effect on G(r) |
|---|---|---|
| Qmin | 0.2 / 1.0 / 2.0 | 0.14 / 0.020 / 0.44 |
| DQ | 0.005 / 0.02 / 0.05 | 1.3×10⁻³ / 2.6×10⁻³ / 1.0×10⁻² |
| Qmax (termination) | 5.0 / 6.0 | 0.58 / 0.38 |
| `set delta` | 0.02 / 0.05 / 0.10 | engine aborts |
