# Gate 1 — pre-registered C1 dyadic-grid experiment

**Question.** Is the archived C1 `set aver,100` residual (1.77 % relative L2) caused by DISCUS's integer phase table at the non-dyadic hkl step 0.01? See `../REVIEW_2026-10-28.md`, section "Gate 1".

**Prediction**, registered in `predictions.json` before any engine run:
- On the dyadic grid (corner −20/128, step 1/128, 41³; almost the archived ±0.16 window), the engine matches the DISCUS-semantics oracle (`four_aver` as read from the 6.23.01 source; engine-precision form factors):
  - total ≤ 1×10⁻⁹ relative L2;
  - diffuse ≤ 1×10⁻⁶ relative L2.
- Re-running the archived macro reproduces the archived volumes to ≤ 1×10⁻¹² of max.

**Interpretation, fixed in advance:**
- *Dyadic diffuse passes:* the averaging algorithm is identified exactly, and the archived residual is the engine's phase-table precision on non-dyadic grids. Gate 1 can be marked RESOLVED.
- *Dyadic diffuse misses at the percent level:* the phase-table hypothesis is refuted. The averaging model is incomplete and Gate 1 stays DIAGNOSED.
- *Archived rerun fails:* the engine or binary on this machine does not reproduce the archived evidence. Investigate before anything else.

**Files**
- `c1_dyadic.py`: oracle, `predict`, `analyze RUN_DIR`.
- `inputs/`: `C1_frozen.stru` (SHA-256 `ceefca57…`), `c1_dyadic.mac`, `c1_repro.mac`.
- `predictions.json`

It reads the archived volumes from `structure_factor_benchmark_v4.2.tar.gz` at the repo root, SHA-256 checked. Run instructions: `../LOCAL_AGENT_HANDOFF.md`, Step 5.
