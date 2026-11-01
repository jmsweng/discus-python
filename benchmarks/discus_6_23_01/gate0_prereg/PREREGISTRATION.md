# Gate 0 v2 — pre-registered prescribed-ensemble diffuse scattering experiment (DISCUS_SUITE 6.23.01)

**Registered:** 2026-11-01. The current version is **v1**.

**Why:** Gate 0 is PARTIAL: the forward model is verified, but C1's diffuse was only compared with occupancies **fitted from the sample itself**; the independently prescribed C1 ensemble is not documented. Gate 0 v2 replaces that with an ensemble that **is** prescribed before any engine run, so DISCUS's `set aver` diffuse can be compared with an *analytic ensemble expectation* that uses no sample information. It also extends coverage beyond Gates 2/3 (which were one-site cells, B = 0) to a **multi-site cell** and **B ≠ 0**. **C1 itself stays sample-conditional; the report must say so.**

Everything that decides PASS/FAIL is fixed here and in `predictions.json` before the engine runs:
- structures
- macros
- oracle
- predicted numbers
- tolerances
- analysis code

`gate0_analyze.py` refuses to evaluate (status INVALID) if any registered file, or any input copied into the run directory, differs from its registered SHA-256.

## Question

Does DISCUS 6.23.01, run on a prescribed multi-site fcc crystal with B ≠ 0, produce a `set aver` diffuse volume whose ensemble statistics (mean, variance, per-voxel distribution) agree with the analytic expectation to the registered tolerances? The test covers engine precision (E1–E4) and ensemble statistics (S1–S4), for both Arm A (B = 0, `temp ignore`) and Arm B (B_Cu = 0.5, B_Au = 0.8, `temp use`).

## Design and why each choice removes a confound

| Choice | Value | Confound it removes |
|---|---|---|
| Lattice | fcc, conventional cubic cell a = 3.75 Å, **4 sites per cell** at (0,0,0), (½,½,0), (½,0,½), (0,½,½) | Multi-site cell: `set aver` averages per site over unit cells. With 4 sites, ⟨F⟩ is fixed by composition and site occupancies alone. |
| Crystal | repeat(average, (8, 8, 8)) → N_c = 512 cells, N = 2048 sites | Sufficient size for the central limit theorem to make x = \|D\|²/E ≈ Exp(1). |
| Species | Cu (p = 0.75) / Au (p = 0.25), every site mixed | High contrast; uncorrelated occupancies by construction (independent Bernoulli draws). |
| B | Arm A: 0 (`temp ignore`); Arm B: B_Cu = 0.5, B_Au = 0.8 (`temp use`) | Covers both the B = 0 convention (Gate 2 v3+) and the B ≠ 0 case (new in Gate 0 v2). The log must show the B values DISCUS used. |
| hkl grid | dyadic, corner (0,0,0), step 1/8, 25 points per axis (h, k, l ∈ {0, 1/8, …, 3}) | **Commensurate** with the 8-cell crystal: the finite-crystal sum equals the periodic sum exactly. L(h) = 0 off integer h, so E = 4 N_c p(1-p)\|g\|² off integer h and E = 0 at integer h. |
| Realizations | K = 64 configurations `sample_configuration(crystal, seed=s_k)` with 64 **fresh** seeds registered in `predictions.json` | Independent, uncorrelated occupancies by construction. Seeds drawn once with `secrets.randbits(63)` and recorded; never run anything on them before registration. |

## Oracle (`oracle.py`, numpy only; nothing imported from `src/discuspdf`)

- **Analytic ensemble expectation:**
  E\|D(h)\|² = p(1-p) \|g(s)\|² · Σ_j [N_c - \|L(h)\|²/N_c]
  where g(s) = f_Au(s) e^{-B_Au s²} - f_Cu(s) e^{-B_Cu s²}, L(h) = Σ_R e^{2πi h·R}, s = \|h\|/(2a).
  On this grid L(h) = 0 unless h is integer, so E = 4 N_c p(1-p)\|g\|² off integer h and E = 0 at integer h.
  This is the **independent** oracle: it uses no sample information, only the prescribed p and B.

- **Realization-exact oracle:** D_k(h) = Σ_j g_j Σ_R δc_{jR,k} e^{2πi h·(R+x_j)} from the actual k-th structure file. This separates engine precision from ensemble statistics.

- **Two precision modes:**
  - **physics**: DABAX W–K coefficients, exact s. This is the independent reference.
  - **engine**: DISCUS's own coefficient table, with s rounded to 0.001 Å⁻¹ using Fortran NINT, exactly as `four_formtab`/`istl` do.

- **δc_{jR}** is the deviation of site j in cell R from that site's average over cells (what DISCUS `set aver` subtracts — confirmed against the Gate 1 `four_aver` oracle).

## Registered predictions

- **Seeds:** 64 fresh seeds drawn with `secrets.randbits(63)`, registered in `predictions.json`.
- **Thresholds:**
  - E1: total vs realization engine oracle, rel L2 ≤ 1×10⁻⁹ per realization.
  - E2: diffuse vs realization engine oracle, rel L2 ≤ 1×10⁻⁶ per realization (looser than Gate 1's 6.9×10⁻¹² because B ≠ 0 is new; record the achieved maximum).
  - E3: diffuse at integer h ≤ 1×10⁻⁶ × max(diffuse) per realization.
  - E4: Arm B log shows B_Cu = 0.5 and B_Au = 0.8 used; Arm A log shows `temp ignore`.
  - S1: over S, mean over voxels and realizations of x = D/E_eng within 1 ± 5/√(K·|S|).
  - S2: over S, the sample variance of x within 1 ± 5·√(8/(K·|S|)) (Exp(1): var 1, fourth central moment 9).
  - S3: every voxel of the **whole grid** off integer h: x̄_v = mean_k x lies inside the two-sided Gamma(K, 1/K) interval — **except voxels with 2h integer (q ≡ −q mod 1), where D is real and x ~ χ²₁, so x̄_v ~ Gamma(K/2, 2/K)** — with family-wise false-alarm 1×10⁻³ (Bonferroni over the number of tested voxels; quantiles computed with `scipy.stats.gamma` at registration time and stored in `predictions.json`).
  - S4: physics vs engine expectation agree to ≤ 1×10⁻⁶ relative (documents form-factor/table consistency).

- **Independent evaluation set S:** h in the first zone [0,1)³ of the grid, q ≠ 0, 2q not integer, and keep only one of each pair {q, −q mod 1}. Across realizations every voxel is independent.

- **Gamma quantiles:** computed with `scipy.stats.gamma.ppf` at registration time and stored in `predictions.json`. For standard voxels: Gamma(K, 1/K). For 2h-integer voxels: Gamma(K/2, 2/K). Bonferroni alpha = 1×10⁻³ / n_tested.

## Acceptance criteria

These are evaluated for *both* arms by `gate0_analyze.py`. Values come from `predictions.json → thresholds`.

| ID | Criterion | Limit |
|---|---|---|
| E1 | engine total vs realization engine oracle (rel L2), max over realizations | ≤ 1×10⁻⁹ |
| E2 | engine diffuse vs realization engine oracle (rel L2), max over realizations | ≤ 1×10⁻⁶ |
| E3 | engine diffuse at integer h / max diffuse, max over realizations | ≤ 1×10⁻⁶ |
| E4 | Arm B: log shows B_Cu = 0.5 and B_Au = 0.8; Arm A: log shows `temp ignore` | required |
| S1 | over S: \|mean x − 1\| | ≤ 5/√(K·\|S\|) |
| S2 | over S: \|var x − 1\| | ≤ 5·√(8/(K·\|S\|)) |
| S3 | every voxel off integer h: x̄_v inside two-sided Gamma interval (Bonferroni FWER 1×10⁻³) | 0 violations |
| S4 | max \|E_physics − E_engine\| / E_engine | ≤ 1×10⁻⁶ |

Metadata must also match before any criterion is evaluated: shape 25³, corner 0, increments diag(1/8), cell 3.75 Å / 90°, and `bragg_present` / `bragg_subtracted` as appropriate. Otherwise the run is INVALID.

**Gate 0 v2 PASSES only if every criterion passes for both arms on real engine output.**

**Scope/limits:** C1 stays sample-conditional. The analytic ensemble expectation uses the prescribed p and B, not occupancies fitted from the sample. The report must state this explicitly.

## Stop rules and deviations

1. **E1 fails:** stop. The forward model or the inputs are not what was registered, so do not interpret S.
2. **E1/E2 pass but S1/S2 fail:** the ensemble statistics differ from the analytic expectation. Gate 0 v2 is UNRESOLVED. Report the failing voxels; do not reinterpret.
3. **E1–E4 and S1–S4 pass:** Gate 0 v2 PASSES. The prescribed ensemble is documented and verified.
4. **No edits after the run.** Registered files, thresholds and the analysis must not be edited after engine data exist. If a genuine bug is found in `gate0_analyze.py`:
   - preserve the original result;
   - write a deviation note;
   - register a new version (`gate0-v2-prereg-v2`) with its own commit;
   - report both results.

## Files

- `oracle.py`: independent oracle (numpy only, no `discuspdf` imports).
- `gate0_generate.py`: generator (writes structures, macros, predictions.json).
- `inputs/`: registered structures and macros.
- `predictions.json`: seeds, thresholds, Gamma quantiles, SHA-256 hashes.
- `gate0_analyze.py`: the only evaluator.
- `gate0_selftest.py`: proves the evaluator passes clean synthetic data and fails nine kinds of corruption (p = 0.30, correlated, DW convention, ×1.01 scale, sample-fitted, species swapped, missing `set aver`, shuffled seeds, Friedel noise, duplicated realization).

## Selftest (offline; must be run before registration)

Clean synthetic data built from the oracle (with complex-Gaussian sampling for the ensemble part) → PASS. Each corruption must FAIL the criterion named:
- p = 0.30 instead of 0.25 (S1)
- a correlated configuration with Warren–Cowley α₁ ≈ −0.2 (S2/S3)
- DW convention e^{−2Bs²} in amplitude (E2 and S1, Arm B)
- intensities × 1.01 (E1/E2/S1)
- sample-fitted occupancies used in the expectation (must make S-criteria uninformative → selftest asserts the evaluator never reads sample occupancies; static check)
- species swapped Cu↔Au in one realization (E1)
- missing `set aver` (E3)
- realizations shuffled against their seeds (E1)
- Friedel-violating noise 1e-4 (E1)
- one realization duplicated (S2 or a duplicate-hash check)
