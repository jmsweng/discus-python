# DISCUS Fourier controls: Bragg subtraction and LOTS (Package B, Features 8 and 9)

`discuspdf.FourierConfig` writes the `fourier` and `output` blocks of a DISCUS_SUITE 6.23.01 macro
from a `ReciprocalGrid`. It offers only commands whose behaviour was observed on 6.23.01; the
observations and numbers are in
[`benchmarks/discus_6_23_01/b_probes/RESULT.md`](../benchmarks/discus_6_23_01/b_probes/RESULT.md).
The macro is plain text: render it, read it, then run it with `Diffuse3DRunner`.

```python
from discuspdf import (Lattice, ReciprocalGrid, FourierConfig, Lots, Seed, StructureSource,
                       Diffuse3DRunner, prepare_fourier_run)

lat = Lattice.from_parameters(4, 4, 4, 90, 90, 90)
grid = ReciprocalGrid.centered(lat, (128, 128, 128), (1/16, 1/16, 1/16))   # commensurate with 16 cells

cfg = FourierConfig(grid, outputs=(("int", "diffuse.h5"), ("3DPDF", "diffuse_3dpdf.h5")),
                    aver_percent=100)                       # subtract <F> from all cells
structure = StructureSource.stru("model.stru", crystal_cells=(16, 16, 16))
print(cfg.render(structure))                                # inspect before running

runner = Diffuse3DRunner(output_filename="diffuse.h5")
job = runner.execute(prepare_fourier_run(runner, cfg, structure, "path/to/model.stru"), native_nexus=True)
vol = job.load_discus_volume()                              # with_bragg == 'bragg_subtracted'
```

## What each option does in DISCUS 6.23.01

| Option | Macro | Observed behaviour |
|---|---|---|
| `aver_percent=0` | `set aver,0.0` | total scattering; label `bragg_present` |
| `aver_percent=100` | `set aver,100.0` | intensity = \|F − ⟨F⟩\|², **amplitudes** subtracted; on a commensurate grid the diffuse part off Bragg is unchanged and Bragg voxels become 0 (one undisplaced site per cell) |
| `0 < aver_percent < 100` | `set aver,p` | ⟨F⟩ from a random subset of cells ("Used 50.3 % …"); **needs `seed=Seed(8 values)`** |
| `lots=Lots("box"/"ellipsoid", (lx,ly,lz), n or "all", periodic)` | `lots box,4,4,4,64,y` | per-lot intensities summed and scaled to the whole crystal; `all` = one lot per origin cell (N³ periodic, (N − l + 1)³ open); random lots need a `Seed` |
| `symmetry=None/"friedel"/"ignore"/"apply"` | `set symm:…` | default merges Friedel pairs; `ignore` differs by ~1e-15 |
| `temperature="ignore"/"use"` | `temp …` | with `use`, DISCUS replaces B < 0.001 Å² in the file by 1.0 |
| `outputs` | `value int/real/imag/3DPDF`, `form nexus` | `real`/`imag` are F − ⟨F⟩, rounded to single precision |

Bragg label: `cfg.bragg` is `bragg_subtracted` when `aver_percent > 0`, else `bragg_present`.
A label is what the configuration requested; whether the subtraction suits a given model is a
scientific question (see "Not covered").

## Refused, with the reason

| Request | Why |
|---|---|
| `nufft` (no option exists) | 6.23.01 wrote **all-zero** files with exit code 0 (`FINUFFT … setpts error (ier=4)`) |
| `set aver` < 100 or random lots without a `Seed` | the result depends on random draws |
| a 1-value seed | not reproducible in 6.23.01 (two runs, different files); `Seed` needs 8 values in 1…99999 |
| any output other than `int` with LOTS | DISCUS error −60 |
| lot larger than the crystal | DISCUS error −101 (checked in Python when `crystal_cells` is known; LOTS require it) |
| `value amp`, other wavelengths, neutrons, `disp` | not verified on 6.23.01 |
| axis labels without a component along their index | DISCUS error −4 ("Component of increment vector is zero") |
| file names with spaces, commas or paths | DISCUS macro arguments |

**Axis labels.** DISCUS needs each active axis to have a nonzero hkl component along its `abs/ord/top`
label. `FourierConfig` picks `(h, k, l)` when that works and otherwise the first permutation that
does (an `(h h l)` plane gets `abs h, ord l, top k`). `cfg.labels_verified` is True only for
`(h, k, l)` on an hkl-aligned grid; the relabelled oblique case was run live for one `(h h l)` plane
(below), so other oblique planes follow the same rule but have not each been run.

## Runner changes

`Diffuse3DRunner` now also treats as failed:
- any DISCUS error line (`***APPL***`, `***FOUR***`, `***MAC ***`, … — the start-up banner of bare
  asterisks is not matched) and `FINUFFT … error` lines;
- native NeXus output whose values are all zero.
On the 35 probe logs the pattern flags exactly the 7 failed runs, and on the earlier Gate logs only
the v1 run that DISCUS rejected.

`prepare_fourier_run` stages the macro, the structure file (with LF line endings; a Windows
checkout has CRLF) and `fourier_config.json` (config + structure, hashed in `inputs.json`).

## Evidence

- `tests/test_fourier_macro.py` (offline, in CI): 25 configs render **byte for byte** the macros
  DISCUS ran in the probes (every successful turbo probe except the `tl ≠ ll` plane and the 1-value
  seeds, which are deliberately not offered); refusals; the engine facts above asserted on the
  committed evidence files; runner failure detection with a fake engine. Ten mutations of the new
  code (labels, seed rule, periodic flag, `aver` format, lot count, int-only, lot size, both error
  patterns, zero check) are each caught.
- `tests/test_fourier_live.py` (`DISCUS_LIVE=1`), **run on DISCUS 6.23.01 on 2026-10-31, 7/7 passed**
  ([`live_b3/results.json`](../benchmarks/discus_6_23_01/b_probes/live_b3/results.json)):
  - generated macros run through `prepare_fourier_run` + `Diffuse3DRunner` reproduced five probe
    outputs **byte for byte** (`set aver 100`, `set aver 50` with seed, 64 random lots, all lots +
    `aver 100`, `symm:ignore`; 9 files);
  - the `(h h l)` plane with the generated labels `abs h, ord l, top k` ran (the probe with `ord k` had
    been refused); its l = 0 row equals the h = k diagonal of the probe's hk0 plane to 9.7e-8
    (relative L2), and DISCUS wrote the increment columns (0.1, 0.1, 0), (0, 0, 0.1), (0, 0, 0);
  - a lot larger than the crystal gave status `failed` with `engine_error = ***MAC ***`.

## Not covered

Displacements and nonzero B (Bragg subtraction then also changes off-Bragg intensity),
multi-site cells, `set aver` on non-commensurate grids (Gate 1's 1.8 % residual), `set mode:accumulate`,
`value 3DPDF` with LOTS, `nufft`, neutron/electron radiation. Loading a run: `job.load_scattering_volume()`
(see [SCATTERING_DATA.md](SCATTERING_DATA.md)) checks the output against the staged config.
