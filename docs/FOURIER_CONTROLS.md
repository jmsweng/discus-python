# DISCUS Fourier macros: Bragg subtraction and LOTS

`FourierConfig` writes the `fourier` and `output` part of a DISCUS_SUITE 6.23.01 macro for you,
from a [`ReciprocalGrid`](RECIPROCAL_GRID.md) and a few options:

- **Bragg subtraction** (`set aver`), to get diffuse scattering without the Bragg peaks;
- **LOTS**, to average the scattering of many small pieces of a large crystal;
- a fixed random **seed**, so runs can be repeated exactly;
- which **outputs** to write (intensity, real/imaginary part, 3D-PDF) as NeXus files;
- the **radiation** (X-ray or neutron) and the **wavelength**.

It only offers commands whose behaviour was checked on DISCUS 6.23.01, and it refuses settings
that were seen to give wrong or unrepeatable results. The macro is plain text: you can print it and
read it before anything runs.

## Quick start

```python
from discuspdf import (Lattice, ReciprocalGrid, FourierConfig, Lots, Seed, StructureSource,
                       Diffuse3DRunner, prepare_fourier_run)

lat = Lattice.from_parameters(4, 4, 4, 90, 90, 90)
grid = ReciprocalGrid.centered(lat, (128, 128, 128), (1/16, 1/16, 1/16))   # fits a 16×16×16-cell crystal

cfg = FourierConfig(grid, outputs=(("int", "diffuse.h5"),), aver_percent=100)   # subtract Bragg
structure = StructureSource.stru("model.stru", crystal_cells=(16, 16, 16))
print(cfg.render(structure))          # read the macro first

runner = Diffuse3DRunner(output_filename="diffuse.h5")
job = runner.execute(prepare_fourier_run(runner, cfg, structure, "path/to/model.stru"), native_nexus=True)
sv = job.load_scattering_volume()     # labelled 'bragg_subtracted'; checked against cfg
```

The printed macro contains, among other lines:

```
fourier
  ...
  set aver,100.0
  run
  exit
```

## Options

| Option | DISCUS command | What you get |
|---|---|---|
| `aver_percent=0` (default) | `set aver,0.0` | total scattering (Bragg peaks included) |
| `aver_percent=100` | `set aver,100.0` | diffuse scattering: the average structure's amplitude is subtracted before squaring, \|F − ⟨F⟩\|². On a grid that fits the crystal, Bragg points become 0 and the rest is unchanged (for a model without displacements) |
| `0 < aver_percent < 100` | `set aver,p` | ⟨F⟩ is estimated from a random p % of the unit cells; **needs a `seed`** |
| `lots=Lots("box", (4, 4, 4), 64, periodic=True)` | `lots box,4,4,4,64,y` | 64 random 4×4×4-cell pieces, summed and scaled to the whole crystal; **needs a `seed`** |
| `lots=Lots("box", (4, 4, 4), "all", periodic=True)` | `lots box,4,4,4,all,y` | one piece starting at every unit cell (no seed needed) |
| `seed=Seed((1, 2, 3, 4, 5, 6, 7, 8))` | `seed 1,2,…,8` | repeatable random choices; exactly 8 numbers, each 1–99999 |
| `symmetry="ignore"` / `"apply"` | `set symm:…` | default (`None`) lets DISCUS merge Friedel pairs, as it does anyway |
| `temperature="use"` | `temp use` | apply the structure's B values (DISCUS replaces B < 0.001 Å² by 1.0) |
| `outputs=(("int", f), ("real", f2), ("imag", f3), ("3DPDF", f4))` | `value …`, `form nexus` | one NeXus file per entry; `real`/`imag` are of F − ⟨F⟩, in single precision |

`cfg.bragg` tells you how the result will be labelled: `bragg_subtracted` when `aver_percent > 0`,
otherwise `bragg_present`. That is what was requested. Whether subtraction is appropriate for your
model is a scientific judgement (see "Not covered" below).

**"Fits the crystal"** means the grid step is 1/N in each direction for an N-cell crystal (as with
1/16 above). Then Bragg positions fall exactly on grid points.

## What it refuses, and why

| Request | Reason |
|---|---|
| the `nufft` method | DISCUS 6.23.01 wrote files of all zeros and still reported success |
| `set aver` below 100, or random LOTS, without a `Seed` | the result would change from run to run |
| a seed with one number | not repeatable in 6.23.01 (two runs gave different files) |
| any output other than `int` with LOTS | DISCUS stops with error −60 |
| a piece (lot) larger than the crystal | DISCUS stops with error −101 (checked in Python; give `crystal_cells`) |
| `value amp`, electrons, `disp anom` (anomalous dispersion) | not checked on 6.23.01 yet |
| a wavelength symbol DISCUS does not know (e.g. `cuka`) | DISCUS would stop; use `cua1` or a number in Å |
| file names with spaces, commas or folders | they break DISCUS macro arguments |

**Axis labels.** For each grid axis DISCUS needs an `abs`/`ord`/`top` label (h, k or l) along which
the axis actually moves; otherwise it stops with error −4. `FourierConfig` uses `h, k, l` when that
works and otherwise picks a working order (for an h h l plane: `abs h, ord l, top k`).
`cfg.labels_verified` is True when the labels are the plain `h, k, l` on a grid along h, k, l. The
reordered case was run for one h h l plane, and other oblique planes follow the same rule.

## Radiation and wavelength

```python
cfg = FourierConfig(grid, outputs=(("int", "diffuse.h5"),), aver_percent=100,
                    radiation="neutron", wavelength=1.54)      # or wavelength="cua1"
```

| Option | Values | Macro |
|---|---|---|
| `radiation` | `"xray"` (default), `"neutron"` | `xray` / `neut` |
| `wavelength` | a DISCUS symbol (`"moa1"` default, `"cua1"`, `"aga1"`, `"mob"`, … — `WAVELENGTHS`) or a number in Å | `wvle moa1` / `wvle 1.54` |

The radiation line is always written before `wvle`, because `neut` clears a wavelength symbol set
earlier. X-ray and neutron runs (with `moa1`) are checked against an independent calculation in the
[D reference probe](../benchmarks/discus_6_23_01/d_reference/RESULT.md) (PASS, ≤ 1.8e-13); the
generated `fourier` block equals the macros that probe ran, byte for byte.

**The wavelength does not change the intensities** in these macros. They always contain `disp off`
(no anomalous dispersion), and then DISCUS 6.23.01 uses the wavelength only for the file's metadata
(`data_rad_symbol`, `data_rad_length`): X-ray form factors depend on sin θ/λ = |Q|/4π, not on λ, and
neutron scattering lengths are constants. In the 6.23.01 source the Fourier calculation uses λ only for
anomalous terms (`disp anom`) and extinction, neither of which these macros switch on. Setting it
still documents the experiment, and `load_scattering_volume()` checks that the file's radiation,
symbol or length (to 1e-4) match the configuration.

## Running and checking

- `prepare_fourier_run(runner, cfg, structure, path)` puts the macro, a copy of the structure file
  (with Unix line endings) and `fourier_config.json` into the run folder, with hashes in `inputs.json`.
- `Diffuse3DRunner` marks a run **failed** when DISCUS prints an error line (`***APPL***`,
  `***FOUR***`, `***MAC ***`, …), when FINUFFT reports an error, or when an output file is all zeros.
- `job.load_scattering_volume()` reads the result and checks it against the saved configuration
  (see [SCATTERING_DATA.md](SCATTERING_DATA.md)).

The measurements behind every row above are in
[`benchmarks/discus_6_23_01/b_probes/RESULT.md`](../benchmarks/discus_6_23_01/b_probes/RESULT.md).

## Not covered

- Models with atomic displacements or nonzero B values (Bragg subtraction then also changes the
  scattering away from Bragg points), and cells with several sites.
- `set aver` on grids that do not fit the crystal (an earlier test left a 1.8 % residual).
- `set mode:accumulate`, `value 3DPDF` together with LOTS, `nufft`, electron radiation, anomalous
  dispersion. A numeric `wvle` has not yet been run in DISCUS (the source accepts it; the symbol `moa1`
  has, for X-rays and neutrons).

<details>
<summary>How this was checked</summary>

- `tests/test_fourier_macro.py` (runs in CI, no DISCUS needed): 25 configurations produce, **byte
  for byte**, the macros that DISCUS ran in the engine tests; all refusals; the DISCUS behaviour in
  the tables, checked against the committed output files; the runner's failure detection, using a
  fake DISCUS. Ten deliberate bugs in the new code are each caught.
- `tests/test_fourier_live.py` (needs DISCUS; `DISCUS_LIVE=1`) **passed 7/7 on DISCUS 6.23.01 on
  2026-10-31** ([`live_b3/results.json`](../benchmarks/discus_6_23_01/b_probes/live_b3/results.json)):
  - macros generated here reproduced five engine-test outputs byte for byte (`set aver 100`,
    `set aver 50` with a seed, 64 random lots, all lots with `set aver 100`, `symm:ignore`);
  - the h h l plane with the reordered labels ran, and its l = 0 row equals the h = k diagonal of
    the hk0 plane (relative difference 9.7e-8);
  - a lot larger than the crystal ended with status `failed`.
- On the 35 engine-test logs the error-line check flags exactly the 7 runs that failed.
</details>
