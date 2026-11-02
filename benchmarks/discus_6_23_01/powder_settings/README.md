# Powder-PDF settings probe (powder-settings-v1)

What `set delta`, Qmin, DQ and Qmax termination do in DISCUS_SUITE 6.23.01's 1D powder PDF. This is for the TASKS.md P1 item "Verify powder settings against real DISCUS outputs". Qmax and wavelength were already measured on 2026-10-31 (`docs/VERIFIED_DISCUS_MACRO.md`).

**Registered before any engine run:** `inputs/` (13 macros), `predictions.json` (macro, structure and `probe.py` hashes, limits) and the analyzer in `probe.py`. `run` refuses to start if any of them differs.

## Design

Each case is the notebook-01 macro (`render_pdf_macro`): GaAs 5×5×5, X-ray, λ = 1.5406 Å, r = 0–12 Å in 0.02 Å steps. Exactly one setting is changed from the base (qmin 0.5, qmax 7.5, dq 0.01, delta 0). Qmax 7.5 is below the 8 Å⁻¹ limit, so termination is real. Each macro also writes F(Q) (`form powder, q` / `value F(Q)`), placed before the unchanged PDF block. `base_plain` is the base macro without that block.

| Cases | Setting |
|---|---|
| `delta002`, `delta005`, `delta010` | `set delta` 0.02, 0.05, 0.10 |
| `qmin02`, `qmin10`, `qmin20` | Qmin 0.2, 1.0, 2.0 (2.0 is above GaAs (111), Q = 1.925) |
| `dq005`, `dq020`, `dq050` | DQ 0.005, 0.02, 0.05 |
| `qmax50`, `qmax60` | Qmax 5.0, 6.0 (termination at other points) |

Source reading (DiffuseCode tag `v.6.23.01`), which the criteria test:
- `set delta,d` convolves the pattern along Q with exp(−x²/d²)/(√π d), so σ = d/√2 (`powder_write_mod.f90` `powder_conv_res`). This should damp G(r) by about exp(−d²r²/4) (informational).
- Before the transform, F(Q) is set to zero from its last sign change up to Qmax (`zero_last`). Below Qmin it is replaced by a straight line to F(0) = 0 (`discus_powder_fft.f90` `fft_fq`).

| Criterion | Limit |
|---|---|
| T1: each G(r) = (2/π) Σ F sin(Qr) ΔQ of the same run's F(Q), with the rules above; best of two Q placements (written Q, or one step lower) | max rel 2×10⁻³; scale within 2×10⁻³ of 1 |
| T3: F(Q) at common Q does not depend on qmin, dq or qmax | 1×10⁻⁶ |
| T4: the F(Q) block does not change G(r) (`base` vs `base_plain`) | 1×10⁻¹² |
| I1: `set delta` damping exp(−d²r²/4) | informational |

The placement variant that matches is reported, and it is labelled post hoc: the source alone does not settle it.

## Running it (owner's machine, WSL)

Use a fresh LF clone of the branch on the Linux side. The commands are in the PR / conversation. In short:

```bash
python benchmarks/discus_6_23_01/powder_settings/probe.py selftest    # 6 lines starting "ok"
python benchmarks/discus_6_23_01/powder_settings/probe.py run --work ~/discus_powder/work
```

`run` copies the macros and `examples/gaas.cll` into a new work folder and runs each macro once. It then copies everything to `evidence/<RUN_ID>/`, together with `ENGINE_ENV.txt` (binary SHA-256), `RUN_STATUS.json` and `MANIFEST.sha256`, and runs the analysis (`results.json`, `RESULT_TABLE.md`). Exit code 0 means PASS. It never edits inputs or deletes anything.
