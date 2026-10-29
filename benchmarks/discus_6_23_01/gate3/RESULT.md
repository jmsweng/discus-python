# Gate 3 (gate3-prereg-v2) — result: **PASS**

Run `20261029T105111Z_LAPTOP-NGU3BT0N` · DISCUS_SUITE 6.23.01 (binary SHA-256 `d4274f13…45706f3`) · registered commit `31b1740`
· `predictions.json` SHA-256 `69b7f475…8ecadfe` · evaluator report: `evidence/20261029T105111Z_LAPTOP-NGU3BT0N/GATE3_RESULTS_PREREG.md`.

**0 of 15 registered criteria failed.** v1 never ran: its clean selftest failed on the operator machine because of an evaluator tolerance defect (see `PREREGISTRATION.md`).
The `.h5` outputs were not committed. Their content hashes are in `../EVIDENCE_CONTENT_HASHES.json`, and they can be recreated with `../regenerate_evidence.py`.

| Criterion (both `set aver,0` and `set aver,100`) | Result |
|---|---|
| R0: `value int` reproduces Gate 2 v4 exactly (content hash) | yes, both, on a fresh run |
| M1: native corner / increment equal the source prediction (−7.937008 / 0.1259843) | exact (deviation 0) |
| E1: native = emulation applied to the same run's `value int` | 3.0×10⁻¹⁶, 2.9×10⁻¹⁶ |
| E2: native = blind registered prediction | yes |
| P1: native maximum at index (64, 64, 64) — the true origin, which DISCUS labels u = +0.126 | yes |
| C1: `delta_pdf_fractional` = k × native, k = 5.524×10⁻³ | 5.1×10⁻¹⁶ |
| L: logs clean, `Temp. factors : ignored` | yes |

## What Gate 3 establishes

- **Absolute normalization (part A).** `delta_pdf_fractional(value int)` gives ρ in e² Å⁻³ for the whole crystal, with no hidden scale. At every
  predeclared lattice vector it agrees with the independent physics prediction to 1.1×10⁻⁵ ρ(0). Divide by N_cells for a per-cell value.
- **DISCUS's native `value 3DPDF` is fully understood (part B).** It is Re[FFT_forward(I′)]/√(n₁n₂n₃): a unitary-DFT quantity with no
  physical units, which depends on the grid size. I′ = I, except that for `set aver 0` on non-commensurate grids the Bragg voxels are re-weighted by
  1/(N_cells |det vi|). Convert with ρ = √(n₁n₂n₃) |det vi| / V_cell × native.
- **The native axis metadata are mislabelled.** The written step is 1/(vi (n−1)); the true step is 1/(vi n). That is a stretch of n/(n−1): 0.8 % here, 3 % on 33³.
  For even n, the labelled origin is also off by one voxel. Users of native 3DPDF files must relabel the axes. `value int` + `delta_pdf_fractional`
  avoids both problems and is the recommended path.
- The same formula reproduced the archived 33³ native files from PR #8 to 2×10⁻¹⁶ (post hoc). That resolves the earlier "5.517×10⁻³ × direct + constant" observation.

## Not covered

The `sharp` option, `output range` (d*-limited), Lanczos filtering, `lots`, skewed grids, B ≠ 0, multi-site cells, other DISCUS versions.
