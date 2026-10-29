# Gate 3: absolute normalization and DISCUS's native `value 3DPDF` (pre-registered)

**Registered:** 2026-10-28, before any native 3DPDF engine output on this grid existed. `predictions.json` SHA-256 `4b878496c559a190cf5792ab50506ad8e065fd85ca3ac920db9372e7e9644e8b`.
Source read: DiffuseCode tag `v.6.23.01` (`19a214a`).

## Part A: absolute normalization of the package pipeline (established; no new run needed)

`delta_pdf_fractional` returns

  ρ(u) = (ΔH ΔK ΔL / V_cell) Σ_H I(H) e^{−2πi H·u},

with H in reciprocal-lattice units and I the DISCUS `value int` intensity in electron units, for the whole crystal. Gate 2 C0
showed that the intensity is the unscaled |Σ f e^{2πiH·x}|²: the agreement is 2×10⁻¹⁶ at engine precision, so no hidden scale factor exists.
ρ is therefore in **e² Å⁻³ for the whole crystal**; divide by N_cells for a per-cell value. At lattice vectors n:

  ρ(n) = (1/V_cell) Σ_τ C(τ) K(n − τ),  where K(n) = ΔH³ Σ_H g(H)² e^{−2πi H·n},

C(τ) is the pair sum Σ_j σ_j σ_{j+τ} and g = (f_Au − f_Cu)/2. K is the form-factor kernel truncated to the measured box.
Gate 2 v4 registered and passed this in absolute units, not only as ratios (B1 relL2 against the independent oracle).
Restated here for the v4 block diffuse at every predeclared τ: the engine-derived ρ matches the physics oracle to 1.1×10⁻⁵ ρ(0), and the
engine-precision oracle to 9×10⁻⁹ ρ(0), where the remainder is the known rounding ties. ρ(0) = 9.691×10⁶ e² Å⁻³.

**Scope:** commensurate grids, B = 0, one site per cell (the Gate 2 conditions).

## Part B: the native `value 3DPDF` (this experiment)

From the source (details in `native3dpdf.py`), DISCUS writes

  native(p) = Re[ FFTW_FORWARD(I′) ] / √(n₁n₂n₃)

where I′ is the intensity array shifted so that H = 0 is the FFT origin. Taking the real part makes the axis reversal in the
reverse map irrelevant. I′ = I for `set aver ≠ 0`. For `set aver 0`, Bragg voxels are divided by scalef = N_cells·|det vi| (`four_weight`).

The written axes are increment 1/(vi (n−1)) and corner −inc·INT((n−1)/2). The true axes of the same values are increment 1/(vi n),
with u = 0 at output index INT(n/2) for even n.

**Consequences, all of them predictions:**
1. **The native output has no physical units.** It is a unitary-DFT quantity that depends on the grid size. It converts to the package's
   absolute ρ as ρ = k·native, with k = √(n₁n₂n₃)·|det vi|/V_cell = 5.524272e-03 here. This holds only when scalef = 1, or for `set aver ≠ 0`.
2. **The axis labels are stretched by n/(n−1).** That is 3 % on 33³ grids and 0.8 % here. For even n, the label also misplaces the origin by one voxel:
   DISCUS labels the true u = 0 voxel as u = +0.126. Odd n is centred correctly.
3. **Total-intensity runs:** on non-commensurate grids, DISCUS re-weights the Bragg voxels before the FFT. On this commensurate grid,
   scalef = 4096·(1/16)³ = 1 exactly, so no re-weighting.

**Post hoc evidence (before registration):** applied to the archived 33³ runs `correlated_run/corr_total` and `short_range_run/sr_total`
(N_cells = 1, scalef = 0.125³), the emulation reproduces DISCUS's native files to 2×10⁻¹⁶, with both corner and increment exact. Other N_cells values fail (relL2 > 0.87).
This explains the earlier "5.517×10⁻³ × direct + constant, 4–5 % residual" observation: the unitary scaling, the Bragg re-weighting
and the n/(n−1) axis labels together produce it. These checks run in CI (`tests/test_gate3_prereg.py`).

## Design

The Gate 2 v4 block structure, grid and `temp ignore` are unchanged. Two macros, `set aver,0` (total) and `set aver,100` (diffuse), each write
`value int` and then `value 3DPDF` from one Fourier run. Predicted native fingerprints were computed from the Gate 2 v4 engine volumes,
recovered from git history and checked against their content hashes.

## Criteria (both macros unless noted)

| ID | Criterion | Limit |
|---|---|---|
| R0 | `value int` content hash equals Gate 2 v4 (engine reproducibility) | required |
| M1 | native corner and increment equal the source prediction (−7.937008, 0.1259843); style `single_pdf` | 1×10⁻⁹ |
| E1 | native values vs the emulation applied to the same run's `value int` | relL2 ≤ 1×10⁻¹² |
| E2 | native values vs the registered prediction fingerprint | rtol 1×10⁻⁹ |
| P1 | native maximum at index (64, 64, 64), the true origin | required |
| C1 | diffuse only: package `delta_pdf_fractional` equals k × native on the true axes | relL2 ≤ 1×10⁻¹² |
| L | log: 2048 Cu / 2048 Au, no errors, `Temp. factors` not `used` | required |

**Gate 3 PASSES only if every criterion passes.** A selftest (`selftest.py`) passes clean data and fails each of these: missing 1/√N, origin
shift, corrected metadata, wrong Bragg re-weight, changed intensities, temp used, and an edited macro (INVALID).

**Not covered:** the `sharp` option, `output range` (d*-limited), Lanczos filtering, `lots`, non-cubic or skewed grids, and odd grids beyond the archived 33³ runs.
