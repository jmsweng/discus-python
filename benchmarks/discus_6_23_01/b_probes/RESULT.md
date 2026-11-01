# Package B engine probes — result (DISCUS_SUITE 6.23.01)

**Evidence:** `evidence/20261031T182936Z_LAPTOP-NGU3BT0N/` (commit `43ee250`, 35 runs, 116 files, 13 MB,
`check_evidence.py` PASS). Owner's laptop, WSL Ubuntu-22.04, `/usr/local/bin/discus_suite`
SHA-256 `d4274f13f938e8fda040a2a3636f1accfc18ec3c10ca8a1b37069ac4d45706f3` (recorded in the B-P1
report; the runner's own lookup through `wsl.exe` returned nothing, so `runs_status.json` has `null`).
Every number below is printed by `analyze_b_probes.py`; that script and this page were written after
the data were seen. Analysis was done by the reviewer in the cloud, not on the engine machine.

## Summary

| Question | Answer (6.23.01) |
|---|---|
| Plane/line metadata | single-point axis → zero column in `data_increment_vector`, or `tl − ll` if given; values unaffected |
| Oblique plane `hhl` with `abs h, ord k, top l` | **refused** (`Component of increment vector is zero`, −4) |
| `set aver` subtracts | **amplitudes**: `int = |F − ⟨F⟩|²`; `value real/imag` give `F − ⟨F⟩` |
| `set aver 100` on a commensurate grid | ⟨F⟩ = 0 off Bragg, ⟨F⟩ = F on Bragg for this model; diffuse is unchanged off Bragg and exactly 0 on Bragg |
| `set aver 50` | random cells; reproducible with an 8-value `seed`, **not** with a 1-value `seed` |
| LOTS normalisation | intensity scaled to the whole crystal: `box,16,16,16,2,y` = `lots off` to 6e-18 |
| LOTS counts | `all,y` → 4096 lots (one per cell); `all,n` → 2197 = (16 − 4 + 1)³ |
| Random LOTS | reproducible with an 8-value `seed` |
| `value amp` with LOTS | **refused** (`Output value NOT allowed using lots`, −60) |
| Lot larger than crystal | refused (−101) |
| `set technique:nufft` | **broken without failing**: FINUFFT `setpts error (ier=4)`, exit 0, all-zero output file |
| `set symm` | default = `friedel` = `apply` (P1 structure); `ignore` differs only at 1e-15 |

## P1 — planes and lines (GaAs `read cell …,4,4,4`)

`data_increment_vector` columns (= `ReciprocalGrid.axis_vectors` rows) as written by the engine:

| Run | shape | columns |
|---|---|---|
| `p1a_hk0` (`nt=1`, `tl=ll`) | 41×41×1 | (0.1,0,0), (0,0.1,0), **(0,0,0)** |
| `p1b_h00` | 41×1×1 | (0.1,0,0), **(0,0,0)**, **(0,0,0)** |
| `p1d_tloff` (`nt=1`, `tl = ll + (0,0,1)`) | 41×41×1 | (0.1,0,0), (0,0.1,0), **(0,0,1)** — values identical to `p1a` |
| `p1e_h0l` (`no=1`) | 41×1×41 | (0.1,0,0), **(0,0,0)**, (0,0,0.1) |

- `ReciprocalGrid` already ignores the vector of a one-point axis, so it reads all four. `DiscusVolume`
  rejects them (singular increment matrix); plane files must be read through `ReciprocalGrid`.
- `ReciprocalGrid.from_corners` currently **refuses** `tl ≠ ll` with `nt = 1`; DISCUS accepts it and
  ignores it. Keep the refusal (it catches a mistyped corner) but document it.
- **`p1c_hhl` was refused** at `run`: the ordinate step (0,0,0.1) has a zero *k* component while the
  axis is labelled `ord k`. DISCUS requires each step to have a nonzero component along its labelled
  index. A macro writer must choose `abs/ord/top` accordingly (`ord l` here). Not yet run with the
  fix → follow-up probe.

## P2 — `set aver` (gate2_random, 16³ cells, commensurate 32³ grid, 64 integer-hkl voxels)

- `value int` = |`value real` + i·`value imag`|² to 6e-11 (aver 0) and 3e-8 (aver 100): DISCUS
  subtracts **amplitudes** and the real/imag outputs are `F − ⟨F⟩`.
  `value real/imag` are stored in float64 datasets but rounded to single precision (relative
  ≤ 1.2e-7 per voxel, within the float32 bound everywhere); `value int` is double precision.
- ⟨F⟩ = F(aver 0) − F(aver 100): max |⟨F⟩| off Bragg 1.5e-11 (|F| up to 2.2e5); on Bragg ⟨F⟩ equals F
  exactly (one undisplaced site per cell). So `aver 100` leaves the off-Bragg intensity unchanged
  (relative 5.8e-15) and sets the 64 Bragg voxels to exactly 0.
- `aver 50`: "Used 50.3 %" (seed A, twice: byte-identical), "Used 48.7 %" (seed B), "Used 50.3 %" /
  "49.9 %" for two runs with the same **single-value** seed (different files). Off-Bragg intensity is
  unchanged (5.8e-15); the Bragg voxels keep a residual |F − ⟨F⟩₅₀|² of up to 1.5e-4 of the total
  Bragg maximum.

## P3 — LOTS (gate2_block, 16³ cells, same grid)

| Run | lots reported | vs `lots off` | Bragg-voxel sum ratio | total ratio |
|---|---|---|---|---|
| `box,16,16,16,2,y` | 2 box | **6e-18** | 1.000 | 1.000 |
| `box,4,4,4,all,y` | 4096 box | 0.985 | 0.01566 | 0.153 |
| `box,4,4,4,all,n` | 2197 box | 0.985 | 0.01566 | 0.153 |
| `ell,4,4,4,all,y` | 4096 ellipsoid | 0.992 | 0.00784 | 0.153 |
| `box,4,4,4,64,y` seed A (×2) | 64 box | 0.985 | 0.01537 | 0.150 |
| `box,4,4,4,64,y` seed B | 64 box | 0.985 | 0.01564 | 0.153 |
| `box,4,4,4,all,y` + `aver 100` | 4096 box | 1.0 | 3e-5 | 0.030 |

- Normalisation: per-lot intensities are summed and scaled by N_cells/(cells per lot · n_lots); the
  whole-crystal lot reproduces `lots off`.
- With 4³ lots the Bragg-voxel ratio 0.01566 ≈ (4/16)³ = 0.015625: on this grid (step 1/8) a lot of
  4 cells spreads each reflection over neighbouring voxels, so grid sums of `lots off` and LOTS are
  not comparable without accounting for the sampling. LOTS changes what "diffuse" contains (Bragg
  tails of finite lots), which is why the spec warns not to conflate finite-domain averaging with
  Bragg removal.
- Seed A twice: byte-identical. `value amp` with LOTS: refused (−60). Lot 17³ in a 16³ crystal:
  refused (−101).

## P4 — technique

- `set technique:nufft` printed `FINUFFT invokeGuru: setpts error (ier=4)!` in every run that got
  past setup, **exited 0 and wrote an HDF5 file of zeros** (commensurate and non-commensurate grids,
  with and without `set aver`). The runner's error pattern did not match that line, so those runs
  show `OK` in `runs_status.json`. FINUFFT error 4 means non-uniform points outside the allowed range;
  these structures have coordinates 0…15 (not centred on the origin).
- `nufft` + `set aver 100` on the non-commensurate grid was refused (`Reciprocal grid does not
  include zero point`, `Shift corners by h = 0.2`, −188), although h = 0 is a grid point.
- `turbo` with `lots off`/`aver 0` reproduced `p3_off` byte for byte.

## P5 — symmetry (odd 33³ grid, aver 100)

Default log: "Friedel pairs merged". `symm:friedel` and `symm:apply` are byte-identical to the
default (P1 structure, so `apply` adds nothing); `symm:ignore` ("Not applied") differs by 1.2e-15.
Even grids also report "Friedel pairs merged".

## Consequences for the Python API (B-3)

1. Generate only `set technique:turbo`. Refuse `nufft` until a probe shows correct output for a
   given structure, and treat `FINUFFT … error` in a log as a failed run in `Diffuse3DRunner`
   (also check for all-zero output).
2. `set aver` is amplitude subtraction; with `p < 100` the result depends on the random cell
   selection, so require an explicit 8-value seed and record it. Never offer a 1-value seed.
3. LOTS: offer `box`/`ellipsoid`, size ≤ crystal, `n` or `all`, periodic flag; only `value int` is
   allowed with LOTS; require an 8-value seed for random lots; record that the result is scaled to
   the whole crystal.
4. Planes: read through `ReciprocalGrid`; macro writer must pick `abs/ord/top` so that each step has a
   nonzero component along its labelled index (needs one follow-up probe).
5. Status `complete` from a run must also require a nonzero output (the `nufft` failure mode).

## Not covered

Displacement disorder, anisotropic/nonzero B, multi-site cells, `set aver` on non-commensurate grids
(the Gate 1 residual), `set mode:accumulate`, neutron/electron radiation, `value 3DPDF` with LOTS.
