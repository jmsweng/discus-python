# Gate 2 — pre-registered signed off-origin 3D-ΔPDF experiment (DISCUS_SUITE 6.23.01)

**Registered:** 2026-10-28, in the git commit that adds this directory. That commit predates any engine run of these inputs.
**Owner approval:** design approved by the project owner on 2026-10-28.
**Replaces:** the retracted Gate 2 comparison in `../GATE2_RESULTS.md` (see `../REVIEW_2026-10-28.md`).

Everything that decides PASS/FAIL is fixed here and in `predictions.json` before the engine runs:
- structures
- macros
- oracle
- predicted numbers
- tolerances
- analysis code

`analyze.py` refuses to evaluate (status INVALID) if any registered file, or any input copied into the run directory, differs from its registered SHA-256.

## Question

Does DISCUS 6.23.01, run on a frozen structure with known, anisotropic, signed short-range occupational correlations, produce a Bragg-subtracted diffuse volume whose 3D-ΔPDF shows the predicted off-origin peaks? The test covers the right **positions**, **signs** and **relative amplitudes**, to stated tolerances, and the pipeline being tested is `src/discuspdf.read_discus_volume` followed by `delta_pdf_fractional`.

## Design and why each choice removes a confound

| Choice | Value | Confound it removes |
|---|---|---|
| Lattice | simple cubic, a = 4 Å, **one atom per primitive cell** | `set aver` averages per site over unit cells. With one site per cell, ⟨F⟩ is fixed by composition alone. |
| Crystal | 16 × 16 × 16 cells (N = 4096), `ncell 16,16,16,1` written explicitly | The PR #8 runs had no `ncell` and one 20 Å P1 "cell", so DISCUS used 1 cell and ⟨F⟩ = F, giving exactly zero diffuse. |
| Atom order | x fastest, then y, then z (DISCUS `celltoindex`) | `four_aver` assigns atoms to cells by index, not by position. |
| Species | Cu (σ = −1) / Au (σ = +1), exactly 2048 / 2048 | High contrast. Σσ = 0 exactly, so ⟨F⟩ has no composition-fluctuation term. |
| B | 0 | Removes Debye–Waller conventions from Gate 2. The convention, exp(−B s²) on F, is separately confirmed on the PR #8 total-intensity runs. |
| hkl grid | H = −4 + k/16, k = 0…127 per axis (`lr = 3.9375`, 128 points, even) | **Commensurate**: the finite-crystal sum equals the periodic sum exactly, so the open-vs-periodic ambiguity is gone. All phases are exact in DISCUS's 2¹⁶ exp table. With an even point count Friedel merging is *not* used, so Friedel symmetry is a real diagnostic. |
| Q range | \|H\| ≤ 4 rlu per axis, s ≤ 0.87 Å⁻¹ | The real-space kernel's leakage one lattice vector away is K(a)/K(0) = −7.8×10⁻⁴. The PR #8 window (d*max 0.1 / 0.05 Å⁻¹) gave 62 % / 89 % leakage at 2.5 Å. |
| Real-space step | 1/(128 · 1/16) = 1/8 cell = 0.5 Å | Every lattice vector is an exact grid point. |

### Order model: `block` (seed 20261028)

Each z layer is tiled by independent 2 × 2 (x, y) blocks with σ(x, y) = s_block · (−1)ˣ. Within a block, x-neighbours are unlike and y-neighbours are like, and s_block = ±1 with probability ½. Analytic ensemble values of C(τ)/N:
- −½ at ±a x̂
- +½ at ±a ŷ
- 0 at ±a ẑ
- −¼ at (±a, ±a, 0)
- 0 at every other vector in the predeclared list

### Negative control: `random` (seed 20261029)

A random permutation of 2048 Cu and 2048 Au. It has no correlation beyond sampling noise.

## Oracle (`oracle.py`, numpy only; nothing imported from `src/discuspdf`)

- **Total intensity:** I_tot = |Σⱼ fⱼ(s) e^{2πi H·xⱼ}|², with s = |H| / 2a.
- **`set aver,100`**, as implemented in `four_aver` (discus/prog/fourier_sup.f90, tag v.6.23.01): ⟨F⟩ = (mean cell structure factor, positions reduced into the first cell) × L(H), where L(H) = Σ_cells e^{2πi H·n}. Here this gives I_diff = g(s)² S(H), with g = (f_Au − f_Cu)/2 and S = |Σ σⱼ e^{2πi H·nⱼ}|².
- **Two precision modes:**
  - **physics**: DABAX W–K coefficients, exact s. This is the independent reference.
  - **engine**: DISCUS's own coefficient table, with s rounded to 0.001 Å⁻¹ using Fortran NINT, exactly as `four_formtab`/`istl` do. On the PR #8 runs this emulation reproduces the engine's total intensity to 6×10⁻¹⁶ and 7×10⁻¹⁵ relative L2.
- **ε_tab:** the physics-vs-engine difference, computed now, before any run. It is 8.0×10⁻⁴ in total, 5.1×10⁻⁴ in diffuse and 5.1×10⁻⁴ in ΔPDF. Every "physics" tolerance is 3 ε_tab.
- **Pair sums** C(τ) = Σⱼ σⱼ σ_{j+τ} are computed by explicit shifted products and cross-checked against FFT.
- **Independent real-space prediction:** ρ(τ) = Σ_τ′ C(τ′) K(τ − τ′). K is computed by zone folding plus explicit 16-point DFT matrices, which is a different code path from the 128³ FFT in the pipeline. Applied to oracle data, `delta_pdf_fractional` reproduces this prediction to 1.1×10⁻¹⁵ of ρ(0). That closes the PR #8 "direct vs FFT" discrepancy: boundaries, grid and normalization now match.

## Registered predictions (block model; physics oracle; ρ(τ)/ρ(0))

| τ (cells) | τ (Å) | realized C(τ)/N | **predicted ρ/ρ₀** | sign required | ensemble (MC mean ± sd) |
|---|---|---|---|---|---|
| ±(1,0,0) | ±4 x̂ | −0.4883 | **−0.4886** | negative | −0.500 ± 0.016 |
| ±(0,1,0) | ±4 ŷ | +0.5508 | **+0.5503** | positive | +0.500 ± 0.015 |
| ±(0,0,1) | ±4 ẑ | −0.0234 | **−0.0242** | – | 0.000 ± 0.032 |
| (1,1,0), (−1,−1,0) | | −0.2637 | **−0.2636** | negative | −0.250 ± 0.013 |
| (1,−1,0), (−1,1,0) | | −0.2666 | **−0.2665** | negative | −0.250 ± 0.013 |
| (1,0,1) | | +0.0078 | +0.0082 | – | 0.000 ± 0.023 |
| (0,1,1) | | +0.0117 | +0.0113 | – | 0.000 ± 0.023 |
| (2,0,0) | | −0.0234 | −0.0232 | – | 0.000 ± 0.032 |
| (0,2,0) | | +0.1016 | **+0.1008** | positive | 0.000 ± 0.031 |
| (0,0,2) | | −0.0312 | −0.0314 | – | 0.000 ± 0.032 |

**Honest note on the realization.** With the declared seed, the y-direction cross-block sum is a +3.25σ draw. That is why C(0,1,0)/N = +0.551 rather than ≈ +0.500, and C(0,2,0)/N = +0.102 rather than ≈ 0. The gate compares the engine with the **realized** pair sums of this frozen atom list, not with the ensemble, so this does not affect the test. The seed was not changed after seeing it. The negative control's realized nearest-neighbour ratios are all |ρ/ρ₀| ≤ 0.015.

## Acceptance criteria

These are evaluated for *both* models by `analyze.py`. Values come from `predictions.json → thresholds`.

| ID | Criterion | Limit |
|---|---|---|
| C0 | engine total vs engine-precision oracle, rel L2 | ≤ 1×10⁻⁶ |
| C0 | engine total vs physics oracle, rel L2 | ≤ 3 ε_tab (2.4×10⁻³) |
| C1 | engine diffuse at integer H (Bragg positions) / max diffuse | ≤ 1×10⁻⁹ |
| C1 | engine (diffuse − total) at non-integer H, rel L2 | ≤ 1×10⁻⁹ |
| A | lattice correlations recovered from one Brillouin zone (C = DFT of I_diff/g²) vs oracle pair sums, max\|ΔC\|/N, physics g | ≤ 1×10⁻³ |
| A | same, engine-precision g | ≤ 1×10⁻⁶ |
| B1 | pipeline 3D-ΔPDF of engine diffuse vs pipeline 3D-ΔPDF of oracle diffuse, rel L2 (physics / engine precision) | ≤ 3 ε_tab / ≤ 1×10⁻⁶ |
| B2 | at every predeclared τ: \|engine ρ/ρ₀ − predicted\| | ≤ 0.01 |
| B2 | sign agreement wherever \|predicted ratio\| ≥ 0.05 | 100 % |
| B3 | imag/real of the transform; Friedel relative error of the engine diffuse | ≤ 1×10⁻¹⁰; ≤ 1×10⁻¹² |
| N | negative control: max \|ρ/ρ₀\| over the six nearest neighbours | ≤ 0.05 |
| L | each engine log: `Atom typ = CU (# 2048)`, `AU (# 2048)`, and no error lines | required |

Metadata must also match before any criterion is evaluated: shape 128³, corner −4, increments diag(1/16), cell 4 Å / 90°, and `bragg_present` / `bragg_subtracted` as appropriate. Otherwise the run is INVALID.

**Gate 2 PASSES only if every criterion passes for both models on real engine output.** Gate 3 (absolute normalization) is explicitly *not* tested here.

## Stop rules and deviations

1. **C0 fails:** stop. The forward model or the inputs are not what was registered, so do not interpret B.
2. **C0 passes but C1 fails:** the `set aver` semantics differ from `four_aver` as read. Gate 2 is UNRESOLVED. Report the failing voxels; do not reinterpret.
3. **C0/C1 pass but A/B fail:** Gate 2 FAILS with quantitative evidence. The next step is a written diagnosis, not a parameter change.
4. **No edits after the run.** Registered files, thresholds and the analysis must not be edited after engine data exist. If a genuine bug is found in `analyze.py`:
   - preserve the original result;
   - write a deviation note;
   - register a new version (`gate2-prereg-v2`) with its own commit;
   - report both results.

## Files

- `oracle.py`, `generate.py`: oracle and generator.
- `inputs/`: registered structures and macros.
- `predictions.json`: numbers, thresholds and hashes.
- `analyze.py`: the only evaluator.
- `selftest.py`: proves the evaluator passes clean synthetic data and fails seven kinds of corruption (no averaging, zero diffuse, transposed axes, 1 % scale, sign-flipped y pattern, engine-without-rounding, edited macro).
- `compare_form_factors.py`, `wk_coefficients.json`: the W–K transcription comparison.
- `tests/test_gate2_prereg.py`: CI checks of all of the above. These are synthetic only and do not show that DISCUS passes.
