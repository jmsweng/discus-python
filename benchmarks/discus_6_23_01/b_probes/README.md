# Package B engine probes (DISCUS_SUITE 6.23.01)

These probes record how DISCUS 6.23.01 behaves for the reciprocal-space controls that Package B
wants to expose in Python (Features 6, 8 and 9). They are **observations, not pass/fail tests of the
engine**: `runs.json` only states whether a run is expected to finish (`ok`), to fail (`error`) or
either. What the numbers mean is analysed separately, by the reviewer, from the committed evidence.

Source reading (upstream DiffuseCode `discus/prog/fourier.f90`, `fourier_sup.f90`, help text) gave the
hypotheses below. Upstream `master` is newer than 6.23.01, so each one needs the engine.

| Probe | Question | Hypothesis from source |
|---|---|---|
| P1 planes/lines (`p1*`) | What do `data_dimension`, `data_corner`, `data_increment_vector` hold for an axis with one point? Is `tl ≠ ll` with `nt = 1` accepted? | 1-point axes are written with their own column; value unknown |
| P2 `set aver` (`p2*`) | Is the subtraction on amplitudes, `|F − ⟨F⟩|²`? Does `value real/imag` give `F − ⟨F⟩`? Is `set aver 50` seed-dependent and reproducible with an 8-value seed, a 1-value seed? | ⟨F⟩ = mean cell factor × lot interference function; `p < 100` selects random cells via `ran1` |
| P3 `lots` (`p3*`) | Normalisation (`N_cells/(cells per lot · nlots)`), periodic vs open lots, box vs ellipsoid, seed reproducibility, `lots` + `set aver 100`, `value amp` with lots, lot larger than the crystal | intensities summed over lots then scaled to the whole crystal; `box,16,16,16,2,y` should equal `lots off`; size > crystal → error −101 |
| P4 technique (`p4*`) | `nufft` vs `turbo` on a commensurate grid (step × cells integer) and a non-commensurate one, with and without `set aver 100` | 6.23.01 NUFFT is wrong when step × cells is not an integer (found in the DiffuseCode performance review) |
| P5 symmetry (`p5*`) | Default vs `set symm:ignore / friedel / apply` on an odd centred grid | default uses Friedel merging; `apply` adds point-group averaging |

Structures: `gate2_random.stru` and `gate2_block.stru` (Gate 2 v4 inputs: simple cubic a = 4 Å,
16³ cells, one atom per cell, Cu/Au 2048/2048, `ncell` written) and `gaas.cll` (`read cell …,4,4,4`).
All macros use `temp ignore` (DISCUS replaces B = 0 by 1.0 on read) and `set technique:turbo` unless
the probe is about the technique.

## Files

- `make_macros.py` — generates `macros/*.mac` and `runs.json` (committed; do not edit them by hand).
- `run_b_probes.py` — runs the macros (WSL from Windows, or a local `discus_suite`), writes
  `evidence/<stamp>/` with logs, HDF5 outputs, `runs_status.json`, `help_sections_b.txt` and
  `MANIFEST.sha256`. Inputs are not copied into the evidence (they are tracked); their hashes are.
- `summarize_b_probes.py <evidence dir>` — read-only HDF5/log summary: `summary.json`, `SUMMARY.md`.

## Running (owner's machine, PowerShell, repository root)

```powershell
.\.venv\Scripts\python.exe benchmarks\discus_6_23_01\b_probes\run_b_probes.py --work C:\tmp\bprobe
.\.venv\Scripts\python.exe benchmarks\discus_6_23_01\b_probes\summarize_b_probes.py benchmarks\discus_6_23_01\b_probes\evidence\<stamp>
```

The kit was dry-run in the cloud with a stub `discus_suite` (no engine); that checks the plumbing only.
