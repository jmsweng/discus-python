# DISCUS 6.23.01 `rmc` menu — source facts used by `discuspdf.discus_rmc`

Files: `discus/prog/rmc.f90` (menu), `rmc_sup_mod.f90` (everything else), `rmc_mod.f90`
(variables), `discus_blk_appl.f90:137-160` (start-up values), `discus_allocate_appl_mod.f90`
(`alloc_rmc`, `alloc_rmc_planes` defaults). Line numbers are for tag `v.6.23.01` (`19a214a`). Every
fact below was read in the code; items not checked by a run are marked *(not run)*.

## Menu (`rmc.f90:145-213`)

| Command | Compared | Calls |
|---|---|---|
| `reset` | 3 of 5 | `rmc_reset` (`:267-282`): clears planes, **deselects every atom type** (`rmc_allowed = .FALSE.`) |
| `select` / `deselect` | 3 of 6 / 2 of 8 | `atom_select` (names, numbers, `all`; `site:` optional) |
| `mselect` / `mdeselect` | 2 of 7 / 2 of 8 | molecule mode (not wrapped) |
| `data` | 3 of 4 | `rmc_readdata` |
| `run` | 2 of 3 | `rmc_run` |
| `save` | 2 of 4 | `rmc_save` (not wrapped; the job saves from the main menu) |
| `set` | 3 of 3 | `rmc_set` |
| `show` | 3 of 4 | `rmc_show`, `show rv` → `rmc_rvalue` |

## `data` (`rmc_readdata`, `rmc_sup_mod.f90:1224-1670`)

`data {xray|neutron}, {symbol|λ}, file, weights, e11,e12,e13, e21,e22,e23, e31,e32,e33` — exactly 13
parameters needed (`:1303`).

- Radiation: first letter `N` or `X` after upper-casing (`:1298-1309`).
- Wavelength: a parameter starting with a letter is a symbol (stored in 4 characters); otherwise a
  number in Å (`:1310-1318`).
- e1, e2, e3: hkl of the lower-left, lower-right and upper-left corners. Steps
  `vi1 = (e2−e1)/(nx−1)`, `vi2 = (e3−e1)/(ny−1)` (`:1500-1503`, `:1564-1569`).
- Weights: first four characters, upper case: `ONE`, `SQRT`, `LOG`, `LIN`, `SQUA`, `INV`, `ISQ`,
  otherwise a weight file (`:1369-1392`). With a weight file the size test compares `wx, wy`, which
  are never set (`:1420`), so weight files are not offered.
- NIPL file (`:1361-1392`): leading lines starting with `#` skipped; `nx ny`; four numbers (ranges,
  stored in `rmc_xy`, not used for the q points); then `DO j = 1, ny: READ (rmc_int(i,j), i=1,nx)`
  (list-directed, so line breaks and blank lines do not matter).
- PGM (`set data, pgm`): the read branch is commented out (`:1459-1493`); not offered.
- After reading, the Laue-equivalent planes are generated (`rmc_symmetry`); more than
  `RMC_MAX_SYM` disables `set sym` for good (`rmc_nosym`, `:1524-1533`).
- The `fourier` menu's NIPL output (`graphic.f90:1786-1802`, `discus_out_file.f90:37-71`) writes
  `nx ny`, the ranges, then `DO j=1,npkt2: (z(i,j), i=1,npkt1)` with `#` header lines: the same
  layout and orientation. `ll/lr/ul/na/no` of `fourier` therefore map to e1/e2/e3/nx/ny.

## `set` (`rmc_set`, `rmc_sup_mod.f90:7-450`; needs ≥ 2 parameters, else error −6)

| Keyword (compared) | Effect | Start value |
|---|---|---|
| `av` | `rmc_ave = value/100`, value in [0, 100] | 0 |
| `ba`, `on\|off`[, b[, plane]] | fit background; b sets `rmc_back` (all planes or one); no b → 0 | off, 0 |
| `co`, p1, p2 | `rmc_constrain(p1) = p2`, needs 0 < p2 ≤ p1 | **0** (allocated as 0) |
| `cy`, n | `rmc_maxcyc` = valid tries | 5000 |
| `da`, `n\|p` | NIPL or PGM | NIPL |
| `db`, `on\|off`[, planes] | Debye–Waller in the form factors, per plane read so far | off |
| `di`, n | log every n generated moves | 500 |
| `log`, `on\|off` | append the log to `rmc.log` | off |
| `lot`, `o…` / `box\|ell`, lx, ly, lz, n, `y\|n`[, file] | lots off / on | off |
| `md`, A, B, d | minimum distance (Å) between types A and B | 0.5 for all |
| `mod`, p, mode[, pairing] | `rmc_set_mode` | all probabilities **0** |
| `mov`, A, sx, sy, sz | Gaussian σ per axis, fractional coordinates | 0.2 |
| `ran`, `all` / qmin, qmax / hkl min, hkl max | limits stored as **half** of 1/d (`:340-361`) | data range |
| `sca`, `on\|off`[, s[, plane]] | fit scale; s sets `rmc_skal = 1/s` | off, 1 |
| `si`, σ | acceptance width | 0.01 |
| `sym`, `on\|off` | use equivalent planes | off |

`set mode` (`rmc_set_mode`, `:565-650`): mode `SHI`, `SWD`, `SWC`, `EXT`, `ROT` (first 3 letters);
pairing `A`, `L`, `SL`, `SI`, `CO`. The probabilities are normalised into cumulative values.

## Defaults that matter

- **No `set mode` → no move.** `rmc_move_prob` and `rmc_move_cprob` start at 0 (`rmc_mod.f90:50-51`).
  `rmc_genmove` (`:3300-3306`) then walks to mode 5 (rotate), which has no branch in the atom
  generator, and `laccept` stays true with `natoms` unset. The wrapper always writes `set mode`.
- **Selection:** `reset` deselects all types; `alloc_rmc` allocates `rmc_allowed` as true. The wrapper
  writes `dese all` then `sele …`.
- **Scale constraint:** `rmc_calcskal` (`:2581-2584`): `IF (rmc_constrain(ip) .ne. ip)` copies scale
  and background from plane `rmc_constrain(ip)`. The start value is 0 for every plane
  (`alloc_rmc_planes`), so without `set const` it reads element 0 of a 1-based array. The wrapper
  writes `set const, p, p` for every plane. *(Effect without it: not run yet.)*

## `run` (`rmc_run`, `:1765-2176`)

- Checks (`rmc_check_input`, `:2424-2542`): atoms present, ≥ 1 plane, ≥ 1 selected type, ≥ 2 for swaps
  (tested against the last mode drawn, `rmc_mode`, start value shift), a displacement spread for
  `swdisp` (error −12).
- Per plane: `chi²_p = Σ w (I − s·dsi − b)² / Σ w`, `dsi = |F|²` (sum over the lots), summed over
  the planes (`:1925-1941`, `:2009-2033`). Error −14 if the starting chi² < 1e-6 (a perfect start).
- Scale and background (`rmc_calcskal`): closed-form least squares when fitted; a result ≤ 0 is
  replaced by the fixed value.
- Acceptance (`:2029-2052`): Δχ² < 0 always; else `exp(−Δχ²/sig2) > RANDOM_NUMBER`, `sig2 = σ²/2`;
  σ = 0 rejects every uphill move. (σ = −9999 accepts everything; not offered.)
- Moves (`rmc_genmove`, `:3263-3416`): `shift` adds `gasdev(σ_axis)` to the fractional position;
  `swchem` exchanges the types of two selected atoms of different type (same-type draws are not
  tries); `swdisp` exchanges displacements from the site average. Each candidate is refused if
  closer than the minimum distance to any other atom (`rmc_check_dist`, `:3098-3143`).
- Loop ends after `set cycl` valid tries; error −20 when 1000 × `set disp` moves are generated with
  no valid try.
- Log (format 1300/1400, `:2166-2169`):
  ` Gen: I9 try: I9 acc: (good/bad): I6 / I6 s2x2: G15.8` then per plane
  `   Plane I2: scal: G11.4 / back: G11.4     /        s2x2: G15.8`; the last block follows
  ` ---- Final configuration ---- `; then ` Delta Chi    : ave: … sig: … max: …` when an uphill move
  was tried. "scal" is 1/s.
- `res[1..8]`: chi², tries, good, bad, ⟨Δχ²⟩, σ(Δχ²), max Δχ², time per try.
- Random numbers: `RANDOM_NUMBER` for every draw (`ran1` calls are commented out); `gasdev` for shifts.
  The `seed` command sets that stream (`lib_f90/random_state_mod.f90`).

## What the fork changes (`rmc_sup_mod.f90`, commit `a299cfa`)

| Line | Stock | Used when |
|---|---|---|
| 2300 | average subtraction loops `j = 1, nx` instead of `ny` | always, but the subtracted array is 0 unless `set aver` > 0 |
| 2371, 2397 | sums lots from `il = 1` after copying lot 1 (lot 1 twice) | more than one lot |
| 2924, 2934 | updates lot 1 instead of lot `il` | lots |
| 2965-3010 | lot membership tests `ilots` (Fourier's) instead of `rmc_ilots`; wrong wrap | lots |

`discuspdf.engine.fork_reasons` sends RMC `set lots` (not off) and `set aver` (not 0) to the fork.
