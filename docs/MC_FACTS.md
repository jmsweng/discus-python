# DISCUS 6.23.01 MMC menu — source-cited facts (Package D)

Source tree: `/home/claude/dc6231` (read-only). All paths below are relative to `discus/prog/` unless they start with `lib_f90/`.
Abbreviations: `mmc.f90:N` = line N. "UNCONFIRMED" = inferred, not verified by running the code.

---

## 0. Top-level structure

- Menu loop: `mmc.f90:66-346`. Commands recognised: `reset` (`str_comp(befehl,'reset',2,...)`, :227 → `mmc_init`), `run` (:232, min 2 chars → `mmc_run_multi(.TRUE.)` for atoms, `mmc_run_multi_mole` for molecules), `calc [corr]` (:261-275, prints current correlations with detail), `set` (:279, needs 3 chars `set`), `show` (:284), `grow` (:289), `apply_symm` (:294), plus generic `@`, `=`, `continue`, `echo`, `eval`, `exit`, `help`, `system`, `wait`.
- **`save` is NOT an mmc-menu command**: the code is commented out (`mmc.f90:241-247`), so typing `save` in the mmc menu falls through to `ier_num=-8, ER_COMM` (unknown command, `mmc.f90:299-302`). The help (`appl_dis.hlp:4874`) still lists it. Save from the main `discus` level.
- Any `set` parse error inside a macro aborts the menu with `ier_num=-10` 'Error occured in mmc menu' (`mmc.f90:307-317`).

## 1. Command parsing (`mmc_set`, `mmc.f90:1002-1352`)

General: `get_params` splits on top-level commas (brackets `()`/`[]` are protected, `lib_f90/get_params_mod.f90:58-85`). `get_optional` then strips `name:value` tokens (names must be lowercase a-z/0-9, else `-12 ER_COMM`, `lib_f90/take_param_mod.f90:64-80`). **Only `cpara(1)` (the sub-command) is upper-cased** (`mmc.f90:1075`); all later tokens keep their case. `set` requires ≥2 parameters except `set output` (`mmc.f90:1076`), else `-6 ER_COMM` (:1342-1344). Unknown sub-command → `-6 ER_COMM` (:1338-1340).
Optional parameters accepted on *any* `set` line: `type:` (default `atoms`), `feed:`, `final:`, `site:` (`mmc.f90:1054-1058`).

| Sub-command | Match rule (after upper-casing) | Effect | Source |
|---|---|---|---|
| `set cyc, N` | `str_comp(...,'CYCLE',2,...)` → `cy`, `cyc`, `cycle`, `cycles`(compares first 5) | `mo_cyc = NINT(N)` (int64) — number of **valid trial moves**, not cycles/atom | `mmc.f90:1140-1143` |
| `set feed, N` / `set feed, auto` | `'FEED'`, ≥2 chars; `auto` ≥3 chars lowercase | `mo_feed = NINT(N)` or `mmc_feed_auto=.TRUE.` | :1148-1155 |
| `set temp, kT` | `cpara(1)(1:2)=='TE'` | `mo_kt = value` | :1281-1284 |
| `set mode, p, type [, site [, atoms...]]` | `(1:3)=='MOD'` | see below | :1188-1196, :3278-3410 |
| `set move, name, sx, sy, sz [, sig]` | `(1:3)=='MOV'` | `rmc_set_move(mo_maxmove,...,dimen=4)` | :1200-1213, `rmc_sup_mod.f90:454-522` |
| `set vec, iv, site1, site2, dx, dy, dz` | `(1:2)=='VE'` | `chem_set_vec` | :1296-1298, `chem.f90:989-1073` |
| `set neig, ...` | `(1:2)=='NE'` | `chem_set_neig` | :1253-1255, `chem.f90:1650-2057` |
| `set target, ic, type, ...` | `(1:2)=='TA'` | `mmc_set_target` | :1276-1277, :1356-2203 |
| `set finish, ...` | `'FINISH'`, ≥3 chars | `get_finish` | :1156-1158, :2559-2631 |
| `set output, feed:on/off, final:on/off` | `'OUTPUT'`, ≥2 chars (`ou`) | `mmc_out_feed`, `mmc_out_final` | :1217-1238 |
| `set valid, N` | `(1:2)=='VA'` | `mmc_no_valid = NINT(N)` | :1289-1292 |
| `set allowed, at...` | `(1:2)=='AL'` | `mmc_allowed(type)=.TRUE.`; `all` → all true | :1080-1090 |
| `set disallowed, at...` | `(1:2)=='DI'` | `mmc_allowed(type)=.FALSE.`; `all` → all false | :1107-1117 |
| others | `AN` angle, `CON` conn, `ENV` envir, `LIM` limited, `RA` range, `ST` style, `ROT` rotation, `FIX` fixed, `RBN` (undocumented algorithm switch) | | :1091-1135, :1162-1184, :1242-1272, :1302-1334 |

Note: test order matters: `'AL'` before `'AN'`; `DI` = disallowed; any word starting with `VA` is `valid`, `VE` is `vec`, `TE` temp, `TA` target, `NE` neig.

### Defaults (`discus_blk_appl.f90`, `mmc_init` at :200-297, called at startup :133 and by `reset`)
- `mo_cyc = 100`, `mo_feed = 10`, `mo_kt = 1.0` (`discus_blk_appl.f90:91-93`). These are **not** reset by the mmc `reset` command (`mmc_init` does not touch them, :200-297).
- `mmc_no_valid = 1000` (`discus_blk_appl.f90:286`, `mmc_mod.f90:175`).
- Shift sigma `mo_maxmove`: allocated with default **0.0** (`discus_allocate_appl_mod.f90:1897`); the 0.2 default from the help (`appl_dis.hlp:5150`) is commented out (`discus_blk_appl.f90:107-109`). **A `shift` mode without `set move` therefore gives zero displacements.** (UNCONFIRMED that nothing else initialises it; grep shows no other assignment.)
- Finish: `mmc_h_stop=.TRUE.`, `MMC_H_NNNN=3`, all thresholds `0.0` (`discus_blk_appl.f90:288-293`). Because the test is strict `<` (`mmc_basic.f90:1462-1464`), 0.0 thresholds can never be satisfied ⇒ by default the run always uses all cycles.
- `mmc_out_feed = mmc_out_final = .true.` (`mmc_mod.f90:69-70`) — **not** reset by `mmc_init`.
- `mmc_algo = MMC_CLASSIC` (`mmc_mod.f90:61`).
- `mmc_feed_auto = .FALSE.` (`discus_blk_appl.f90:295`).

### `set mode` (`mmc_set_mode`, `mmc.f90:3278-3410`)
Syntax as parsed: `set mode, <prob>, <type> [, <site> [, <atom1>, <atom2>, ...]]` (`mmc.f90:1188-1196`).
- `<type>` (upper-cased, first 3 chars): `SHI`→shift (`MC_MOVE_DISP`), `SWD`→swdisp, `SWC`→swchem, `SWV`→swvalue, `SWN`→swneig (+ optional `site:[..]`), `EXN`, `INV`, `ROT`; else `-9 ER_RMC` 'Invalid RMC/MC mode selected' (:3312-3346).
- `<site>` (upper-cased): first char `A`→all, `L`→local, `SL`→slocal, `SI`→site, `CO`→connectivity; else `-9 ER_RMC` (:3350-3365). Note `L` is tested before `SL`, `SI`, so `local`/`L…` works and `slocal` is distinguished only by `S`. Omitted ⇒ `all` (:3348-3349).
- Atom list (only if ≥3 params): `get_iscat`; sets `mmc_allowed_site(imode,type)=.TRUE.` (`all`→all types); type out of range → `-27 ER_APPL` (:3366-3388). `mmc_allowed_site` is **only consulted for SWCHEM with `site`/`slocal`** (`mmc.f90:5624-5629`). It does not make atoms movable in other modes (selection uses `mmc_allowed`, filled by `set target`/`set allowed`).
- Probabilities: `mmc_move_prob(imode)=werte(1)`; cumulative normalised `mmc_move_cprob` (:3393-3404).
- **BUG (code-evident):** when an atom list is given, `get_iscat` (:3368) overwrites `werte`, so `mmc_move_prob(imode)` becomes the *first atom-type number* (e.g. `1` for Cu, `-1` for `all`, `0` for VOID), not `<prob>`. With one mode and a positive type index this is harmless after normalisation; with `all` (−1) or VOID (0) `sump<=0` and `mmc_move_cprob` is left at its old value (:3398), and `mmc_select_atoms` then walks to `i=MC_N_MOVE` (=8, `swvalue`) (`mmc.f90:5526-5531`). With several modes the relative probabilities are wrong.
- Bad `<type>` still executes `mmc_local(imode)` with the previous `imode` (SAVE'd local, initial `MC_MOVE_NONE=0`) — out-of-range index (UNCONFIRMED effect).

### `set move` (`rmc_set_move`, `rmc_sup_mod.f90:454-522`, with `dimen=4`)
- Accepts exactly `name, sx, sy, sz` (ianz=4) or `name, u, v, w, sig` (ianz=5); otherwise `-6 ER_COMM` (:492-522). Only ONE atom name (`ii=1`, :497). `all` → types `1..cr_nscat` (not VOID 0) (:511-516).
- 4-value form: per-axis Gaussian sigma in **fractional (lattice) units**: `disp(i)=gasdev(mo_maxmove(i,type))` (`mmc.f90:6140-6145`). `mo_maxmove(4,·)` keeps its old value (0 unless set before).
- 5-value form: `disp = gasdev(sig) * [u,v,w]` (`mmc.f90:6146-6153`), i.e. 1-D move along `[u,v,w]` (fractional).
- `type:mole` routes to molecule arrays (`mmc.f90:1208-1212`).

### `set vec` (`chem_set_vec`, `chem.f90:989-1073`)
- `set vec, iv, is1, is2, dx, dy, dz` — exactly 6 numbers else `-6 ER_COMM` (:1035, :1066-1069). `is1,is2` are **site numbers in the unit cell** (1..`cr_ncatoms`), not atom types; out of range → `-10 ER_CHEM` 'Invalid crystal site or atom index given' (:1054-1057). `iv<=0` → `-9 ER_CHEM` (:1050-1052). `dx,dy,dz` are integer cell offsets (`int()`, :1061-1063). `set vec, reset` clears all (:1017-1031).

### `set neig` (`chem_set_neig`, `chem.f90:1650-2057`)
- Optional `number:N` sets `chem_ncor=N`; `number:next` increments; absent → `chem_ncor=max(1,chem_ncor)` i.e. overwrites the current neighbourhood (:1708-1721).
- `set neig, vec, iv1, iv2, ...` or `set neig, vec, iv1, iv2, range` (contiguous; `range` ≥2 chars) (:1813-1832, `chem_set_nei_range` :2059-2102). Undefined vector (`chem_cvec(1,iv)==-9999`) → `-9 ER_CHEM` (:1819-1828). Sets `chem_ctyp=CHEM_VEC`.
- `set neig, add` (deprecated) increments `chem_ncor` (:1758-1798); `set neig, rese` → `chem_ncor=0` (:1802-1809).
- **`set neig` must precede `set target`**: `set target` with `ic<=0` or `ic>chem_ncor` → `-14 ER_CHEM` 'Invalid neighbour definition selected' (`mmc.f90:1454-1458`).

### `set target` (`mmc_set_target`, `mmc.f90:1356-2203`)
Parsed list after removing `TARGET`: `[ic, type, ...]`. `ianz<3` → `-6 ER_COMM` (:1425-1429). Type keywords (lowercase, `str_comp` min. lengths): `corr`(2) :1475, `unid`(2) :1477, `group`(2) :1570, `pref`(2) :1572, `cd`(2) :1574, `heli`(4) :1703, `cn`(2) :1802, `spring`(2) :1859, `angle`(2) :1904, `lennard`(2) :1975, `repulsive`(2) :2040, `bucking`(2) :2112, `value`(2) :2176; else `-6 ER_COMM` (:2178-2180). Atom type outside `0..cr_nscat` → `-97 ER_APPL` 'Atom type …' (e.g. :1602-1610). Every successful `set target` increments `mmc_h_number` (:2192). Every target sets `mmc_allowed(type)=.TRUE.` for its atom types (e.g. :1611-1612, :2796, :2807) — hence `set disallowed` must come after the last target.

**corr** (`set_target_occ`, :2207-2296): `set target, ic, corr, A, B, <corr>, <energy>, CORR|ENER`. `A`/`B` may be single names or `(a,b,..)` groups (`get_iscat` handles brackets, `get_iscat_mod.f90:54-60`).
- The keyword test is `cpara(ianz)(1:2)=='CO'` / `=='EN'` — **upper-case only** (:2265-2273), because only `cpara(1)` was capitalised. A lower-case `corr`/`ener` is not recognised; it then remains as a third value for `ber_params` (UNCONFIRMED: likely an expression/variable error).
- If neither keyword matches, `mmc_cfac(ic,MC_OCC)` keeps its previous value (default 1.0 ⇒ CORR).
- CORR: `mmc_cfac=1`, `mmc_lfeed=.true.`; the user `<energy>` is **ignored**; initial depth `winit = 0` if `|corr|<0.1`, else `winit = -corr - 3*corr**5` (:2277-2286).
- ENER: `mmc_cfac=0`, `mmc_lfeed=.false.`, depth = `<energy>` fixed (:2287-2290).
- `mmc_set_disp_occ` (:2765-2829) builds the pair table: `mmc_pair=+1` (A–A), `+2` (B–B), `-1` (A→B), `-2` (B→A); `mmc_left(A)=1`, `mmc_right(B)=1`; same depth for all pairs. Atoms in neither group have `mmc_pair=0`.

**cd** (:1574-1702): `set target, ic, cd, A, B, <corr>, <energy>, CORR|ENER, dir:[u,v,w]` or `dir:[u1,v1,w1,u2,v2,w2]` (or `dir:[all]`). `dir:` is **mandatory** (else `-6` with message '"set target, ic, cd" needs "dir:[u,v,w]"', :1663-1669); 3 or 6 components else `-6` (:1653-1660). Only single atom names (`iianz=jjanz=1`). CORR → `mmc_cfac=5.0`, `mmc_lfeed=.true.`, initial depth `= -2*corr` (user energy ignored) (:1585-1597); ENER → depth = `<energy>`. `ber_params(2,...)` reads exactly 2 numbers.

**spring** (:1859-1903): `set target, ic, spring, A|all, B|all, <d0>, <k>`; reads exactly 2 numbers; `mmc_target_corr=d0`, `mmc_depth=k` for (A,B) and (B,A) (`mmc_set_disp`, :2635-2695). No CORR/ENER keyword, no feedback.

**lennard** (:1975-2039): `set target, ic, lennard, A|all, B|all, <d0>, <E> [, <pr>, <pa>]`.
`a = -|E| * pa/(pa-pr) * d0**pr`, `b = -|E| * pr/(pa-pr) * d0**pa` (:2030-2033), stored with `m=pr, n=pa` (`mmc_set_lenn`, :3072-3138). Defaults are applied only `IF (ianz == 3)` and write `werte(4)=12, werte(5)=6` (:1998-2001) — i.e. with **2 values** (`d0, E`) `pr=pa=0` (werte zeroed at :1468) ⇒ `0/0` = NaN coefficients (UNCONFIRMED at runtime; code-evident); with 3 values `pa` is forced to 12. **Always give all four values**, e.g. `..., 2.55, 1.0, 12, 6` (then E(d0) = −|E|).

Other energy types that exist: `unid`, `group`, `pref`, `heli`, `cn`, `angle`, `repulsive`, `bucking`, `value` (`mc_mod.f90:13-28`: MC_OCC=1, DISP=2, SPRING=3, ANGLE=4, VECTOR=5, BLEN=6, LENNARD=7, BUCKING=8, REPULSIVE=9, COORDNUM=10, UNI=11, GROUP=12, PREF=13, VALUE=14, HELI=15).

### `set finish` (`get_finish`, :2559-2631)
Optional-only: `stop:converge|cycles`, `log:none|screen`, `feed:N` (MMC_H_NNNN), `diff:`, `rdiff:`, `change:`, `aver:`. Defaults **when the command is given**: `stop:cycles`, `log:none`, `feed:3`, all thresholds `99999.` (:2594-2596). Bad `stop`/`log` value → `-6 ER_COMM` (:2602-2623).
⚠ `set finish, stop:converge` without explicit thresholds uses 99999 ⇒ converges as soon as `mmc_h_nfeed>=MMC_H_NNNN` (≈ third feedback call, the initial call counts) (`mmc_basic.f90:1462-1464`, :392).
⚠ `feed:N` with N>3 sizes the history arrays but the PID integral/derivative arrays are fixed `0:2` (`discus_allocate_appl_mod.f90:1840-1851`) and are indexed with `mmc_h_index` ∈ `0..N-1` (`mmc_basic.f90:2216-2219`) ⇒ out-of-bounds (UNCONFIRMED crash/corruption).

### `set output` (:1217-1238)
`set output, feed:on|off, final:on|off` (each defaults to `on` *on that line*, so omitting one re-enables it). Values ≥2 chars (`of`, `off`); else `-6 ER_COMM` 'Feedback parameter must be ''on'' or ''off'''.

### `set valid` (:1289-1292)
Stores `mmc_no_valid`; it is printed by `show` (`mmc.f90:415`) but **never used anywhere in the run** (grep: only `mmc.f90:415,1292`, `discus_blk_appl.f90:286`, `mmc_mod.f90:175`). The help's "terminates mmc once the number of invalid tried configurations exceeds <invalid_moves>" (`appl_dis.hlp:5447-5459`) is not implemented. The real limit is a hard-coded 10,000,000 re-draws in `mmc_select_atoms` → `-22 ER_RMC` 'RMC did not find a valid pair' (`mmc.f90:5510-5520`; also `rmc_sup_mod.f90:3181`).

### `reset` → `mmc_init` (`discus_blk_appl.f90:200-297`)
Clears targets, depths, pairs, move probabilities, allowed flags, angles, finish settings; does **not** clear `mo_cyc/mo_feed/mo_kt/mo_maxmove`, `chem_*` vectors/neighbourhoods, `mmc_out_*`.

### `show` → `mmc_show` (`mmc.f90:350-1000`)
Prints e.g. `'   Max. number of MC cycles     : ',f9.1,' per atom'` = `mo_cyc/cr_natoms` (:948, :411), feedback interval per atom, `mmc_no_valid`, `kT` (formats :947-997).

## 2. Algorithm of `run` (`mmc_run_multi`, `mmc.f90:3414-3991`, classic algorithm)

Setup: `mmc_initial` (`mmc_basic.f90:16-122`) disables `chem_period(i)` for any direction with `cr_icc(i)==1` (:47-50; restored at end, `mmc.f90:3978`), zeroes counters, PID gains `Kp=0.30, Ki=0.05, Kd=-0.10` (:116-118). If >1 OCC target exists: `mmc_cfac(ic,MC_OCC)=2.5/nocc` (`mmc.f90:3572-3578`; cfac is not used by the OCC feedback formula, only as the CORR/ENER flag at parse time). Initial correlations computed with `rel_cycl=0` (no feedback, :3613) and copied to `mmc_ini_*` (:3617-3621).

**Loop (serial path, `nthreads==1`)** (`mmc.f90:3859-3880`): `DO itry=1,mo_cyc`: `mmc_select_atoms` → `mmc_run_loop` → `mmc_run_loop_final`. So **`set cyc` = number of trial moves** (each a *valid* selection), not cycles per atom; the help says the same (`appl_dis.hlp:4993-4997`).

**Move selection** (`mmc_select_atoms`, `mmc.f90:5457-5708`): draw `z=RANDOM_NUMBER`; move type `i` = first with `z < mmc_move_cprob(i)` (:5526-5532). Loop is repeated (up to 1e7) until `laccept`:
- shift/invdisp: one uniformly random atom (`INT(r*cr_natoms)+1`), valid if `mmc_allowed(type)` and property selection OK (:5535-5547).
- swchem: `rmc_select` (`rmc_sup_mod.f90:3145-3261`): atom 1 uniform among allowed types; atom 2 uniform among all atoms (`all`), or same site any cell (`site`), or ±1 cell & random site (`local`), or ±1 cell & same site (`slocal`) — for local modes atom 1 is clamped to cells 2..icc−1 (:3219-3255). Valid only if **types differ** and both allowed and selection status OK (`mmc.f90:5619-5629`). Same-type pairs are silently re-drawn and **not counted** as tries.
- swdisp: as above, valid only if both allowed and `u·v<0` (displacements from average opposite) (:5585-5595).
- Then geometric constraints `check_geometry_multi` (:5698-5701).

**Modification** (`mmc_modify`, :6083-6224): shift adds `gasdev(sigma)` to `cr_pos` (fractional); swchem swaps `cr_iscat` and `cr_prop` of the two atoms; swdisp swaps the displacements relative to `chem_ave_pos`. Undo in `mmc_unmodify` (:6521-6604).

**Energy** (`mmc_energies`, :5742-6079): for each selected atom `ia` (1 or 2), for every neighbourhood `ic`, for every centre returned by `chem_neighbour_multi(isel(ia),ic,...)`, sum the contributions; called before (`e_old`) and after (`e_new`) the modification (`mmc_run_loop`, :4141-4182). Only the environment of the moved atoms is evaluated. `valid_all` is reset per `ia` (:5815), so effectively only the last atom's validity counts.

- **corr** (`mmc_energy_occ`, :6816-6951), only when the move is swchem (:5844):
  `E += mmc_depth(ic,MC_OCC,is,js) * sign(1, mmc_pair(ic,MC_OCC,is,js))` over neighbours (is=centre type, js=neighbour type, orientation via `tatom`), i.e. `E = J·(N_like − N_unlike)` (dimensionless, units of kT). Negative J favours like pairs ⇒ positive correlation (consistent with `winit=-corr-…`).
  ⚠ If a neighbour pair has `mmc_pair==0` (an atom type in neither group), the code sets `mmc_energy_occ = 0.0` (**resets the running sum**, :6885-6888) instead of skipping — in ternary systems this discards contributions already accumulated for that centre.
- **cd** (`mmc_energy_dis`, :7581-7716): `delta = (u_sel · dir1)/|dir1|` (`mmc.f90:5917-5921`), `dx = (u_neigh · dir2)/|dir2|` (:7660), where `u = cr_pos − chem_ave_pos − (cell−1) − cr_dim0` (fractional), and dot products/lengths use the metric tensor `cr_gten` ⇒ **Å**; `E += mmc_depth(ic,MC_DISP,is,js)*delta*dx` (Å²·depth) (:7671). If the selected atom is the *neighbour* (vector end) it uses `mmc_depth(ic,MC_DISP,0,0)` (:7703-7704), which is only set when a target involves type 0/`all` ⇒ normally 0, so only the "forward" direction contributes. The `loop_scat1` over `jjs` (:7645-7677) repeats the neighbour sum once for every type `jjs` paired with `is` ⇒ multiple counting when `is` has several partners (e.g. `all`).
- **spring** (`mmc_energy_spr`, :8034-8131): `E = (1/ncalc) Σ_neigh k_{is,js} (d − d0_{is,js})²`, `d = do_blen(.TRUE.,u,v)` = `sqrt(w·G·w)` in **Å** (`metric.f90:6-31`), only pairs with `d0≠0`; **averaged over the ncalc neighbours** (:8123-8125), not summed. `k` = `<energy>`.
- **lennard** (`mmc_energy_len`, :8225-8333): `E = (1/ncalc) Σ_neigh (a/d**pr − b/d**pa)`, d in Å, also **averaged** per centre (:8326-8328); `mmc_depth` is not used (commented, :8311-8314).

**Acceptance** (`mmc_test_multi`, `mmc_basic.f90:126-227`):
```
e_delta = NINT((e_new(i)-e_old(i))*1.0E4)*1.0E-4   ! rounded to 1e-4, summed over active energy types
if(mo_kt>0) then
   IF (e_del < 0) THEN laccept = .TRUE.
   ELSE IF (mo_kt < 1.0e-10) laccept = .FALSE.
   ELSE e_ran = exp(-e_del/mo_kt); e_ran = e_ran/(1+e_ran); laccept = (e_ran > r1)
else   ! mo_kt<=0: special Lennard-only rule using e_new(7) and mmc_depth(1,7,1,1)
```
(:152-198). So downhill moves are always accepted, uphill (and ΔE=0) moves with the **Glauber/heat-bath probability** `1/(1+exp(ΔE/kT))`; ΔE=0 (incl. |ΔE|<5e-5) ⇒ p=0.5. This hybrid is neither Metropolis nor pure Glauber: `P(+Δ)/P(−Δ) = e^{−Δ/kT}/(1+e^{−Δ/kT}) ≠ e^{−Δ/kT}` ⇒ **detailed balance with respect to exp(−E/kT) is not satisfied** (code-evident). Temperature enters only through `ΔE/mo_kt`. Counters: `good` (ΔE<0), `neut` (ΔE==0, accepted), `bad` (ΔE>0, accepted) (:200-210); rejected moves are not counted anywhere. Average ΔE per type accumulated for accepted moves only (:211-224).
Rejected → `mmc_unmodify`; accepted → if `.NOT.chem_quick` and any periodic direction → `chem_apply_period` for isel(1), isel(2) (`mmc.f90:4222-4231`).

**Lines 4430-4700 ("closer to target" rule)**: that is `mmc_run_grow`, the `MMC_GROWTH` algorithm, used only if `mmc_algo==MMC_GROWTH` (`mmc.f90:3882-3899`), which is set only by the undocumented `set rbn, growth, rate, bad, neut` (`mmc.f90:1091-1103`). Default is `MMC_CLASSIC` (`mmc_mod.f90:61`), so this code is **not executed** in normal runs. Its rule (swchem only, :4459): for each correlation, vote −1 if the local correlation sum moves toward the target, +1 if away (:4630-4691); accept if net <0, else with fixed probability `mmc_g_neut=0.25` (net 0) or `mmc_g_bad=0.05` (net>0) (:4693-4718; defaults `mmc_mod.f90:62-64`). No Boltzmann factor, no temperature, no detailed balance. Its local `damp = 0.01+0.99*exp(-4.0*rel_cycl)` (:4456) is computed but not used for anything in that routine.

**Feedback (CORR mode, OCC)** — performed in `mmc_correlations` at each feedback step:
- Damping actually used: `damp = 0.01 + 0.99*exp(-4.0*rel_cycl)` immediately **overwritten** by `damp = 0.01 + 0.99*exp(-4.0*rel_cycl**2)` (`mmc_basic.f90:350-351`). `rel_cycl = igen/mo_cyc` (`mmc.f90:4320`).
- Only if `lfeed` and `mmc_lfeed(ic,MC_OCC)` (CORR) and `0<rel_cycl<=1` (`mmc_basic.f90:2388-2390`). PID step (`calc_change_pid`, :2191-2238):
  ```
  diff = target - achieved
  inte(:, h) = diff ;  deri(:, h) = achieved - previous_achieved      (h = mmc_h_index, 3 slots 0:2)
  change = -( Kp*diff + Ki*SUM(inte(0:2)) + Kd*SUM(deri(0:2)) ) * damp
  ```
  then scaling (:2397-2414): `|target|<0.001` ⇒ depth:=0, change:=0; `|target|<0.10` ⇒ change*0.005; `|target|<0.2` ⇒ change*0.05; else if `target*achieved>=0` and the step would flip the sign of depth ⇒ `change = -0.1*depth`. Finally `depth += change` for (is,js) and (js,is) (:2416-2417).
  Gains are adapted after every feedback (`mmc_basic.f90:1336-1367`): Kp×1.05 (≤1) if mean ΣΔachieved>0 else ×0.90 (≥0.001); Ki×1.025 or ×0.95; Kd×0.95 if |mean change| grew, else ×1.025.
- **cd feedback** (proportional only): `depth -= cfac*(target-achieved)/2*|target|*damp`, cfac=5 (:1000-1004). ENER (cfac=0) ⇒ no change.
- spring, lennard: no feedback.
- Small CORR targets (|c|<0.1) start at J=0 and adapt with 0.5% of the PID step ⇒ very slow; |c|<0.001 ⇒ J forced to 0.

## 3. Achieved correlation and pair set

**Pair enumeration** (`mmc_correlations`, `mmc_basic.f90:438-875`): loop over all atoms `i` (every atom, no `mmc_allowed`/property filtering), neighbours from `chem_neighbour_multi`, and only pairs with `tatom(j)=.TRUE.` are counted (:575).
For `CHEM_VEC` (`chem_neig_multi.f90:287-355`): if the atom's site == vector start site → neighbour at `cell + (dx,dy,dz)`, `tatom=.TRUE.`; **ELSE IF** site == vector end site → neighbour at `cell − d`, `tatom=.FALSE.`. Hence each vector contributes each pair **once, directed start→end**; a vector with start site == end site never generates the backward neighbour (the ELSEIF is not reached). For symmetric counting **and** for symmetric energies you must define both `+d` and `−d` vectors. The same neighbour routine is used for the energies, so with only `+d` the energy of a swapped atom ignores its `−d` bonds (code-evident; UNCONFIRMED practical impact).
Periodic boundaries: `check_bound` wraps the cell index with `chem_period` and returns an `offset` (`check_bound_mod.f90:7-67`); positions `patom = cr_pos + offset` (`chem_neig_multi.f90:310-312`). Non-periodic directions drop the pair (`lok=.false.`). `chem_period(i)` is forced off for `cr_icc(i)==1` during the run (`mmc_basic.f90:48-50`).

**corr achieved value** (`mmc_correlations_occ`, `mmc_basic.f90:2242-2456`): counts are classified by group: `pneig(a,b)` with a,b ∈ {1=left group, 2=right group, 3=other} (:583-616).
```
r_nneigh = SUM(pneig(:,:))
m_left   = SUM(pneig(1,:))/r_nneigh       ! fraction of pairs whose CENTRE is in group 1
m_right  = SUM(pneig(2,:))/r_nneigh
achieved = 0.5*( (1 - pneig(2,1)/SUM(pneig(2,:))/m_left)  +  (1 - pneig(1,2)/SUM(pneig(1,:))/m_right) )
```
(:2319-2333); if a group is empty → 1.0 (:2332). This is the average of the two Warren-Cowley α's using conditional probabilities and centre-based compositions; the denominators include "other" neighbours.
`chem corr occ` uses instead the pooled `θ = (2p11+p12+p21)/(2N)` and `c = (P11 − θ²)/(θ(1−θ))` over the 4 A/B pair classes only (`chem.f90:4566-4576`). The two are identical only for a binary system with symmetric pair counts (n12=n21) and no third species (algebra: α=1−n12/(θN(1−θ)) = (p11−θ²)/(θ(1−θ))). A Welberry value is computed (:2335-2352) but only shown in detailed output (`mmc_write_detail`).
Printed "Number of pairs" for the feedback lines = `pair11+pair12+pair21+pair22` (:2339, :2439); in the initial/summary table (`mmc_correlation_write`) it is `mmc_ach_pairs = pneig(1,2)` (A→B count) (:2359, :2061/2070) — different meanings.

**cd achieved** (:509-530, :741-759, :977-1033): per pair `dpi = (u_c·dir1)/|dir1|`, `dpj = (u_n·dir2)/|dir2|` (Å, metric `cr_gten`; displacements from `chem_ave_pos`, recomputed by `chem_aver` at every feedback, :406-411); `achieved = <dpi·dpj>/sqrt(<dpi²><dpj²>)` with (is,js) and (js,is) sums pooled (:981-996), i.e. a Pearson-like correlation without mean subtraction. With `dir:[all]` (`chem_ldall`): `dpi=1`, `jdir = u_c` ⇒ `dpj = u_n·u_c/|u_c|` (:516-526).

**spring/lennard achieved**: mean pair distance (Å) and its standard deviation over (is,js)+(js,is) pairs (:1101-1142, :1205-1245).

At the end (`lfinished`) achieved values then sigmas are written to `res[]` (`mmc_basic.f90:1526-1582`): for OCC only the first unlike pair (`js>is`) per `ic`.

## 4. Screen output (all to `output_io`, `prompt_mod`)

Sequence for `run` (`mmc_run_multi`):
1. If `lout_feed = .TRUE. .and. mmc_out_feed` (:3608): `mmc_correlation_write` prints `' --- Initial multiple energy configuration ---'`, blank, header 410, one line per target with **initial** values (`mmc_basic.f90:1875-2011`), then `' ---------------------------------------------'`.
2. At each feedback (`MOD(igen, imodulus*nthreads)==0`, `mmc.f90:4318`), if `lout_feed`: format 2000 line, then `mmc_correlations` prints header 410 and one line per target.
3. Final, if `lout_feed_in .and. mmc_out_final` (:3913): blank + 3 banner lines, **`mmc_correlation_write` again (initial values, still titled "Initial")**, format 3000, final format 2000 (without rel_cycl), header 410 + final achieved lines (`lfinished`, no feedback), convergence block (only if `mmc_h_stop`), format 5000 + 5010 lines, format 4000.

`set output, feed:off` suppresses items 1–2 (initial table and all feedback lines); `final:off` suppresses item 3 (including the final achieved lines, energies table and timing). Feedback/convergence computations still run.

Verbatim formats:
```
! per-feedback progress line (mmc.f90:4341-4342)
 2000 FORMAT (/,' Gen: ',I10,' try: ',I10,' acc: (g/n/b): ',I8,' / ',I8,' / ',I8,'  MC moves ', f6.4)
!   igen (moves done), itry*nthreads, iacc_good, iacc_neut, iacc_bad (cumulative), rel_cycl=igen/mo_cyc
! final line (mmc.f90:3980-3981) — same without f6.4
 2000 FORMAT (/,' Gen: ',I10,' try: ',I10,' acc: (g/n/b): ',I8,' / ',I8,' / ',I8,'  MC moves ')
 3000 FORMAT (/,' --- Final multiple energy configuration ---')
! header (mmc_basic.f90:1585-1589 and 2002-2006)
  410 FORMAT ( 45x,'Correlations/',/ ' Neig.- Energy-',7x,'Atoms',11x,'Target',2x,'Distance/',4x, 'Sigma',5x,'Diff',6x,'Diff/',5x,'Number',/ ' Def.   Type    central  Neighbors',13x,'Angle',27x,'Target  of pairs')
! corr line during feedback/final (mmc_basic.f90:2453-2454)
 3100 FORMAT (1x,i3,3x,'Occupancy',a5,3x,a5, 8x,2(f7.3,3x), 10x,f7.3,3x,f7.3,3x,i8)
!   ic, at(is), at(js), target, achieved, [sigma blank], target-achieved, (target-achieved)/|target| (|target|=0 -> /1), npairs
! cd (mmc_basic.f90:1595-1596)
 3200 FORMAT (1x,i3,3x,'Disp.Cor.',a5,3x,a5, 8x,2(f7.3,3x), 10x,f7.3,3x,f7.3,3x,i8)
!   ic, at(is), at(js), target, achieved, diff, diff/target (signed divider), nneigh
! spring (mmc_basic.f90:1599-1600) / lennard (1605-1606)
 3300 FORMAT (1x,i3,3x,'Hooke    ',a5,3x,a5, 8x,5(f7.3,3x) ,i8)
 3700 FORMAT (1x,i3,3x,'Lennard  ',a5,3x,a5, 8x,5(f7.3,3x) ,i8)
!   ic, at(is), at(js), target d0, mean distance, sigma, d0-mean, (d0-mean)/d0, npairs
! initial/summary table (mmc_basic.f90:2081-2084, 2138-2139): title a9 passed as 'Occupancy','Disp.Cor.','Hooke    ','Lennard  '
 3100 FORMAT (1x,i3,3x,a9,a5,3x,a5, 8x,2(f7.3,3x), 10x,f7.3,3x,f7.3,3x,i8)   ! corr/cd: ic,title,at,at,target,initial,diff,diff/|t|,ini_pairs
 3100 FORMAT (1x,i3,3x,a9,a5,3x,a5, 8x,5(f7.3,3x) ,i8)                      ! spring/LJ: target,initial,sigma,diff,diff/|t|,pairs
! energy-change summary (mmc.f90:3986-3989), one 5010 line per energy type that had accepted moves
 5000 FORMAT (/,' --- Average energy changes ---',//, ' Energy type',14x,'number',17x,'number',2x,'number',/, 29x,'E < 0',5x,' <E>',11x,'E = 0',5x,'E > 0',4x,' <E>')
 5010 FORMAT(a24, 2x,i8,f15.6,2x,i8,2x,i8,f15.6)   ! name, n(dE<0), <dE<0>, n(dE=0), n(dE>0), <dE>0>
 4000 FORMAT (/,' Elapsed time : ',I4,' h ',I2,' min ',I2,' sec ',/ ' Time/cycle   : ',F9.3,' sec',/)
```
Atom names `cr_at_lis` are `LEN=4` (`crystal_mod.f90:89`) written with `a5` ⇒ one leading blank. `f7.3` overflows to `*******` for |x|≥1000 (e.g. relative diff with tiny targets). There is **no "rejected" or "invalid" count** in the output: rejected ≈ Gen − (g+n+b) (minus moves with no valid energy); invalid (same-type) selections are invisible. In the serial path, after normal completion the final `try:` prints `mo_cyc+1` (Fortran DO index after loop, `mmc.f90:3859`, passed as `itry*nthreads`) — UNCONFIRMED but standard Fortran semantics. Convergence messages: `mmc_basic.f90:1482-1523` (`'(a,f12.4,a,f12.4)'`), printed only when `mmc_h_stop`.

## 5. Random numbers, seeds, OpenMP

- All MMC randomness uses the Fortran intrinsic `RANDOM_NUMBER` (e.g. `mmc.f90:5526,5539`, `rmc_sup_mod.f90:3183-3233`, `mmc_basic.f90:171`) and `gasdev` (Box–Muller/Marsaglia polar from `RANDOM_NUMBER`, `lib_f90/random.f90:11-49`); `ran1` is a wrapper of `RANDOM_NUMBER` (:108-122).
- `seed` (suite/discus top level, `lib_f90/kdo_all_mod.f90:179-180` → `do_seed`, `lib_f90/random_state_mod.f90:141-237`): values (optionally grouped via `group:`) → `ini_ran_ix`: copies `min(np, nseeds)` integers into `seed_vals` and `RANDOM_SEED(PUT=...)`; `seed` / `seed 0` → time-based (:95-110); also resets `gasdev` state `iset=0` (:118). `nseeds` = compiler `RANDOM_SEED(SIZE)` (8 for recent gfortran — UNCONFIRMED for this build). If fewer than `nseeds` values are given, the remaining slots keep the previous `seed_vals` (possibly time-derived) ⇒ give all `nseeds` values for reproducibility.
- **Active OpenMP in MMC**: `!$OMP PARALLEL` to query thread count (`mmc.f90:3630-3639`) and the parallel trial loop (`mmc.f90:3736-3771`, `!$OMP DO SCHEDULE(static)`), used when `nthreads>1`, `par_omp_use` and no CON/ENVIR neighbourhood (`lserial`, :3598-3604, :3629). `nthreads` = all OMP threads unless `set parallel, nthread:N` (`lib_f90/do_set_mod.f90:411-440`; default `par_omp_maxthreads=-1`, `lib_f90/parallel_mod.f90:5-6`). Parallel scheme: master serially draws `nthreads` pairs (excluding shared atoms for swchem/swdisp only, :3711-3732), then threads evaluate/modify **shared** `cr_iscat/cr_pos` concurrently, call `RANDOM_NUMBER` inside threads (`mmc_test_multi`), and update shared `n_e_av_*`/`e_aver_*` without atomics. Afterwards an extra serial tail of `max(1,min(nthreads/4*cr_natoms, 2*mo_feed, mo_cyc))` moves (:3832-3855). ⇒ With >1 thread results depend on thread count and timing (races between neighbouring pairs; per-thread RNG streams) — **not reproducible**; no explicit parallel RNG streams exist. The number of moves also changes (`mo_cyc/nthreads*nthreads` + tail). Reproducible runs need `seed` with all values **and** serial execution (`set parallel, useomp:off` or `nthread:1`, or a CON/ENV neighbourhood) and the same binary. (UNCONFIRMED by execution.)
- `!$OMP CRITICAL` in `mmc_run_grow` (:4804-4820) is in the unused growth path.

## 6. Periodic boundaries, positions, save

- Neighbour search wraps cell indices (`check_bound`), periodic only in directions with `chem_period` true and `cr_icc>1` during the run.
- Shift moves change `cr_pos` directly in fractional coordinates (`mmc.f90:6144,6152`); accepted shifts are kept. Wrapping of a moved atom (`chem_apply_period`, `chem.f90:5825-5885`: shift by ±`cr_icc` if outside `[NINT(cr_dim(i,1))−0.1, NINT(cr_dim(i,1))+cr_icc+0.1]`) is applied **only if `.NOT.chem_quick`** (`mmc.f90:4225-4230`); `chem_quick` defaults to `.true.` (`chem_mod.f90:38`) ⇒ by default shifted atoms are never re-wrapped (positions may drift; neighbour lookup is index-based in quick/vec mode so unaffected).
- swchem changes only `cr_iscat` and `cr_prop` (positions/displacements stay at the sites) (:6196-6201).
- `save` (main menu) writes the current `cr_pos`/`cr_iscat` (UNCONFIRMED in detail; not inspected).

## 7. Things that differ from naive expectations (summary)

1. `set cyc` = total valid trial moves (not per atom); same-type swchem draws are silently redrawn and not counted.
2. Acceptance = always-downhill + Glauber `1/(1+e^{ΔE/kT})` uphill; ΔE=0 ⇒ 50 %; ΔE rounded to 1e-4; no detailed balance.
3. `set valid` is a no-op; the only abort is 1e7 failed draws (`-22 ER_RMC`).
4. `damp` in the OCC/cd feedback is `0.01+0.99*exp(-4*rel_cycl**2)`, not `exp(-4*rel_cycl)`.
5. CORR targets: user energy ignored; initial J=−c−3c⁵ (0 if |c|<0.1); PID with adaptive gains; |c|<0.1 adapts ×0.005, |c|<0.001 fixes J=0. Several CORR targets each get their own independent feedback loop (no coupling); the help recommends ENER when combining targets (`appl_dis.hlp:5316-5321`).
6. `CORR`/`ENER` keywords must be upper case.
7. Lennard needs all 4 numbers; spring/LJ energies are per-centre **averages**.
8. Vector neighbourhoods are directed; define ±d for symmetric pair counts and energies.
9. Achieved `corr` is the average of two Warren-Cowley α's (centre-based compositions), equal to `chem corr occ` only for binary, symmetric counts.
10. `set mode` with an atom list sets the probability to the first atom-type index (bug).
11. Default shift sigma is 0 (help says 0.2).
12. `set finish, stop:converge` alone (thresholds 99999) stops after ~3 feedbacks; default (no `set finish`) never converges early.
13. OCC history (`mmc_h_diff`) is filled only inside `IF(lout ...)` (`mmc_basic.f90:2425-2432`), so with `feed:off` OCC targets do not enter the convergence test; then `conv_val` (a SAVE variable, `mmc_basic.f90:290`) keeps old/zero values ⇒ with `stop:converge` and output off, OCC-only runs may stop early (UNCONFIRMED).
14. OpenMP runs are non-deterministic.
15. The "Initial" table is reprinted at the end with initial values; final achieved values follow under "--- Final multiple energy configuration ---".
