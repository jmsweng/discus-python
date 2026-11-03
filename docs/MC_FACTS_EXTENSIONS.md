# DISCUS 6.23.01 MMC: extension facts (neighbourhoods, other targets, finish, misc)

Companion to `docs/MC_FACTS.md` (covers `set vec`, targets `corr`/`cd`/`spring`/`lennard`, swchem/swdisp/shift,
temp/cyc/feed/output, acceptance, log format). Not repeated here.

Source (read-only): `/home/claude/dc_6.23.01/discus/prog/` (`lib_f90/` paths: `/home/claude/dc_6.23.01/lib_f90/`).
"UNCONFIRMED" = read from code, not run. "Code-evident" = follows from the code, not executed.

---

## 0. Case sensitivity

- `mmc_set` upper-cases only the sub-command (`mmc.f90:1075`). `str_comp` (`lib_f90/str_comp_mod.f90`) is a
  minimum-length prefix match on `.eq.`: **case-sensitive**.
- Case-insensitive (`do_cap` runs first): `set mode` type and site (`mmc.f90:3311`, `:3351`); `set neig` type (`chem.f90:1754`).
- Lowercase required: target keywords (`mmc.f90:1475-2176`); `range` (`chem.f90:2087`); `reset`/`all`/`direc`/
  `central`/`neig`/`short` (`chem.f90:1113`, `:1121`, `:1250`, `:1261`, `:1345`, `:1381`, `:1417`, `:1490`, `:1588`);
  `set limited` `cell`/`atom` (`mmc.f90:1164`, `:1174`); `set style` (`:1266`, `:1269`); `set output` `on`/`off`
  (`:1218`, `:1229`); `type:` (`:1203`); `set finish` `converge`/`cycles`/`screen`/`none` (`:2602-2619`).
- Uppercase required: `CO`/`EN` (e.g. `mmc.f90:1516`, `:1714`, `:2265`); `OFF` in `set limited` (`:1182`).
- A wrong-case keyword in a branch with no `ELSE` is **silently ignored** (`set limited, off` does nothing);
  `set target, 1, CORR, …` gives `-6` (`mmc.f90:2178-2180`).

---

## 1. Neighbourhoods other than `vec`

### 1.1 Types (`chem_neig_multi.f90:14-19`)

| Type | Code | Defined by | Search (`chem_neig_multi.f90`) |
|---|---|---|---|
| DIST | 1 | `set neig, dis` | sphere shell `[|r|-σ/2, |r|+σ/2]` via `do_find_env`, all types; optional direction cone (`:75-140`) |
| VEC | 2 | `set vec` + `set neig, vec` | index based, directed (`:287-355`) |
| ANG | 3 | `set angle` + `set neig, ang` | three-site geometry (`:402-544`) |
| ENVIR | 4 | `set env` + `set neig, env` | sphere `[rmin,rmax]`, typed (`:545-…`) |
| RANGE | 5 | `set range` + `set neig, ran` | direction windows (`:141-286`) |
| CON | 6 | `set con` + `set neig, con` | connectivity list of the centre (`:361-401`) |

### 1.2 `set neig` dispatcher (`chem_set_neig`, `chem.f90:1650-2057`)

- `number:` (`chem.f90:1702-1721`): absent → `chem_ncor = max(1,chem_ncor)` (`:1720`), so the **current slot is
  overwritten**. `number:N` → `N` (`:1717`). `number:next` → `+1`.
- `rese` (`:1802-1809`): `chem_ncor=0`; resets only element 1 of the per-type counters. Stale higher counters remain
  (code-evident; impact UNCONFIRMED).
- `add` (`:1758-1798`) is deprecated. Dispatch after `do_cap` (`:1754`): `VE` (`:1813`), `CO` (`:1836`), `RA` (`:1882`),
  `AN` (`:1905`), `ENV` (`:1928`), `DIS` (`:1951`), `DIR` (`:2016`). Undefined ids for vec/ran/ang/env → `-9` (`ER_CHEM`).
- `range` keyword (`chem_set_nei_range`, `chem.f90:2059-2102`): lowercase only (`:2087`); `a, b, range` expands to `a..b`.
- `con, <names>` (`:1836-1857`): exact, case-sensitive match against connectivity names; falls back to numbers
  `1..CHEM_MAX_CON` if a name fails.

### 1.3 Definition commands

| Command | Routine | Syntax as parsed | Notes |
|---|---|---|---|
| `set con` | `chem_set_con` `chem.f90:1075-1206` | `reset` (`:1113`) · `all` (`:1121`) · `iv, <atomtype>, <connname or no.>` (`:1147`, 3 tokens) | `iv<=0` → `-28` (`:1183`); bad type → `-10` (`:1186-1187`). Connectivities come from the CONNECT menu (`appl_dis.hlp:34`). |
| `set angle` | `chem_set_angle` `chem.f90:1462-1555` | `iv, siteA, siteB, du, dv, dw, siteC, eu, ev, ew`: exactly 10 numbers (`:1509`) | Sites are unit-cell site numbers, not types. Integer offsets. Centre = `siteA`. `chem_cwin(1..9,iv)` set at `:1535-1543`. `reset` at `:1490`. |
| `set env` | `chem_set_envir` `chem.f90:1557-1649` | `iv, rmin, rmax, <centre type>, <neighbour types>…` (`:1609`, `del_params(3)`) | First type = centre (`chem_cenv(0,iv)`, `-1` = any). `chem_env_neig` = #types−1 (`:1635`). Distances in Å. `reset` at `:1588`. |
| `set range` | `chem_set_ranges` `chem.f90:1207-1461` | `reset` (`:1250`), `direc, iv, u, v, w, sig [, wsig] [, sym\|nosym]` (`:1261`), `central` (`:1345`), `neig` (`:1381`), `short` (`:1417`) | Only `direc` read in full. Other sub-keywords UNCONFIRMED. |
| `set neig, dir` | `chem.f90:2016-2055` | `all` · `u,v,w` · `u1,v1,w1,u2,v2,w2` | Not used by `heli`, which takes `dir:[…]` (`mmc.f90:1419-1421`). |

### 1.4 `set neig, dis` (`chem.f90:1951-2015`)

- Syntax: `DIS, u, v, w, sigma [, angsig] [, S]`. Branch test `ianz.ge.5` (`:1951`) is before `del_params`.
- `angsig` is read **only if `S` is also given**; otherwise the fifth value is silently dropped (code-evident,
  `:1959-1975`).
- `chem_rmax = |vector|`; `chem_rmin/rmax = |vector| ∓ σ/2` (`:1981-1985`).
- Search (`chem_neig_multi.f90:75-140`): all atom types in the shell; with a cone, only within `chem_wink_sigma` of a
  listed direction (`:87-111`). `S` expands the vector with `rmc_symmetry`.

### 1.5 Target / neighbourhood compatibility

From the `chem_ctyp` test in each energy routine:

| Energy routine | VEC | ENVIR | RANGE | DIST | CON | ANG |
|---|---|---|---|---|---|---|
| occupational `corr` `mmc_energy_occ` (`mmc.f90:6816`; test `:6855`) | yes | yes | yes | yes | yes | no |
| `unid` `mmc_energy_uni` (`:6955`; test `:7045`) | yes | yes | yes | yes | yes | no |
| `group`, `pref` (`:7238`, `:7320`; test `:7277`) | yes | yes | yes | yes | yes | no |
| `cd` (`:7581`; not re-verified) | yes | yes | yes | yes | yes | no |
| `heli` `mmc_energy_heli` (`:7720`; test `:7791`) | yes | yes | yes | yes | yes | no |
| `cn` `mmc_energy_cn` (`:7476`; test ≈`:7490`) | yes | yes | yes | yes | yes | **no** |
| `spring` (`:8034`; test `:8074`) | yes | yes | yes | yes | yes | no |
| `lennard` (`:8225`; test `:8266`) | yes | yes | yes | yes | yes | no |
| `repulsive` `mmc_energy_rep` (`:8338`; test `:8378`) | yes | yes | yes | yes | yes | no |
| `bucking` `mmc_energy_buck` (`:8453`; test `:8493`) | yes | yes | yes | yes | yes | no |
| `angle` `mmc_energy_angle` (`:8560`; test `:8599-8600`) | yes | yes | **no** | **no** | yes | yes |

A mismatch gives **zero energy silently** (code-evident).

---

## 2. Targets

Common: `set target, ic, <type>, …` (`mmc.f90:1356-2203`). `ic` must be `1..chem_ncor` (`:1454-1458`). Success
increments `mmc_h_number` (`:2192`). `get_iscat` splits a single token `(A,B,C)` (`discus/prog/get_iscat_mod.f90:54-62`)
and returns `-1` for `all` (`:93-94`).

### 2.1 `corr` with multi-type lists (`set_target_occ`, `mmc.f90:2207-2296`)

- Syntax: `set target, ic, corr, <grp1>, <grp2>, c [, CO | EN [, E]]`. `<grp>` = `A`, `(A,B,…)` or `all`
  (`get_iscat` at `:2249`, `:2258`).
- `CO` (`:2265`, `cfac=1`, feedback) or `EN` (`cfac=0`, fixed `depth=E`, `:2290`). **No keyword**: `mmc_cfac(ic,OCC)`
  keeps its previous value (initially `1.0`, `discus_blk_appl.f90:225`). Always write `CO` or `EN`.
- CO initial depth: `winit = -c − 3c⁵`, zero if `|c|<0.1` (`:2278-2282`).
- Pairs (`mmc_set_disp_occ`, `mmc.f90:2765-2829`): grp1×grp1 `+1` (`:2793-2803`); grp2×grp2 `+2` (`:2804-2814`);
  grp1×grp2 and reverse `-1`/`-2` (`:2815-2825`). All carry target `c`. Energy: MC_FACTS §2; achieved value: MC_FACTS §3.
- **Clears `(ic,OCC)` first** (`:2787-2791`: `mmc_left`, `mmc_target_corr`, `mmc_depth`, `mmc_pair` zeroed). A second
  `corr` on the same `ic` **erases the first** (code-evident). Use one neighbourhood per independent OCC target.

### 2.2 `unid` (`mmc.f90:1477-1569`)

- Syntax: `set target, ic, unid, <grp1>, <grp2>, c [, CO | EN]`.
- The `IF(is_corr)` branch (`:1480-1528`) is dead (`cpara(2)=='corr'` is false for `unid`).
- CO: `winit = -2c − 3c⁵` (`:1556`; no `|c|<0.1` rule). EN: `depth = E`.
- Uses `mmc_set_unid_occ` (`:2833-2903`), which **does not clear arrays** (contrast `:2787-2791`).
- Energy `mmc_energy_uni` (`:6955-7234`): directed. Depth × sign of `mmc_pair(is,js)`; start/end roles from
  `mmc_pair(is,is)` = ±1 (start) or ±2 (end) (`:7068-7138`). Sign for the `ABS(depth)` branches: UNCONFIRMED.

### 2.3 `cn` coordination number (`mmc.f90:1802-1858`; `mmc_set_cn` `:2699-2761`; energy `:7476-7580`)

- Syntax: `set target, ic, cn, (A), (B), target [, depth]`. `depth` defaults to `1.0` (`:1837`).
- `mmc_cfac(ic,CN)=1.0` (`:1858`), but the feedback code is commented out (`mmc_basic.f90:945-950`): **no feedback**.
- `mmc_pair(ic,CN,is,js) = -1` for `is∈grp1`, `js∈grp2`, **one direction only** (`:2724-2732`).
  `mmc_target_corr = target`, `mmc_depth = depth·target` (`:2731`).
- Energy only when the **moved atom is the centre** (`:7526`, `isel(ia)==iatom(0,icent)`):
  `E = Σ_js depth(is,js)·(n_js − target)²` over `js` with `pair==-1` (`:7545`). `n_js` counts all `iatom` entries of
  type `js`; the `tatom` flag is not tested (code-evident).
- `all` as centre or neighbour gives index `-1`, which never matches (`cr_iscat≥1`): no energy, silently.
- Achieved: `mmc_ach_corr = n_cn/ncentral(is)` (`mmc_basic.f90:944`). Meaning of `n_cn`: UNCONFIRMED.

### 2.4 `angle` (`mmc.f90:1904-1950`; energy `mmc_energy_angle` `:8560-8656`)

- Syntax: `set target, ic, angle, <centre>, <end1>, <end2>, <target_deg>, <depth>` (`ber_params(2)` at `:1943`).
- Stored at `mmc_target_corr(ic,ANGLE,js,ls)` and `(ls,js)` (`:1945-1946`) for the two **end** types. **The centre
  type is not part of the key**; it only sets `mmc_allowed` (`:1928`, `:1937-1938`).
- Energy: for each neighbour pair `ii<jj` of the centre with non-zero end-type target (`:8612`), `a` = angle at the
  centre in degrees (`do_bang`, `metric.f90`), `b` = target:
  - `b ≤ 90` and `a > b`: `depth·(mod(a+b/2, b) − b/2)²` (`:8625-8631`);
  - otherwise: `depth·(a−b)²` (`:8634-8636`).
  Summed over pairs.
- **A target of exactly 0° is never evaluated** (`:8612`).
- Neighbourhoods accepted: VEC, CON, ENVIR, ANG only (`:8599-8600`).

### 2.5 `repulsive` (`mmc.f90:2040-2111`; energy `mmc_energy_rep` `:8338-8448`)

- Syntax: `set target, ic, repulsive, <T1>, <T2> [, A [, B [, C [, M]]]]`. Defaults: `A=0, B=1, C=0, M=1`, filled from
  the right as arguments are omitted (`mmc.f90` after `ber_params`).
- `mmc_set_disp(ic,REP,is,js,100.0,|A|)` (`:2106`): gate `100`, depth `|A|`. `mmc_set_rep` (`:2107`, `:3142-3211`)
  stores `a=|A|, b=B, c=C, m=M`; `rep_low = min|A|`.
- Per neighbour distance `d` (Å), gate ≠ 0 (`:8409`):
  - `d > C`: `−|A| + 1/((d−C)/B)^M` (`:8418-8425`);
  - `d ≤ C`: `−|A| + 1/(1e-6)^M` (`:8425-8431`).
  Average `E = (1/ncalc)Σ` (`:8441`). No neighbours: `E = −|rep_low|` (`:8444`), a constant per centre.

### 2.6 `bucking` (`mmc.f90:2112-2175`; energy `mmc_energy_buck` `:8453-8552`)

- Syntax: `set target, ic, bucking, <T1>, <T2>, <p1>, A, rho, B, <dist>, <depth>, rmin, atmin`: **eight numbers**
  after the types. `werte(1)` (`p1`) is never used.
- `A=werte(2)`, `rho=werte(3)`, `B=werte(4)`, `rmin=werte(7)`, `atmin=werte(8)` (`:2158-2166`).
- `mmc_set_disp(ic,BUCK,is,js,werte(5),|werte(6)|)` (`:2171`). **Energy gated on `werte(5) ≠ 0`** (`:8523`).
- Per neighbour: `d > rmin` → `depth·(A·exp(−d/rho) − B/d⁶)` (`:8533-8537`); else `atmin` (`:8538`). Averaged (`:8549`).

### 2.7 `heli` (`mmc.f90:1703-1801`; energy `mmc_energy_heli` `:7720-7910`)

- Syntax: `set target, ic, heli, <T1>, <T2>, <value>, <depth> [, CO | EN] [, dir:[u,v,w]]`. `dir` default
  `[1.00, 0.00, 0.00]` (`:1419-1421`, `get_optional` at `:1704`).
- `CO` → `cfac=5.0`, `lfeed=.true.` (`:1715-1716`); `EN` → `cfac=0` (`:1719-1720`).
- **Pitfall**: `werte(2) = −2·value` whenever `cfac≠0` (`:1724-1725`). With no `CO`/`EN`, `cfac` keeps its default
  `1.0` and the user's depth is overwritten. Always give `EN` or `CO`.
- `mmc_set_disp(ic,HELI,is,js,value,depth)` (`:1742`).
- Energy per neighbour: `depth·½·e1·e2·(1+e3)`. Centre branch (`:7826`) divides by `sqrt`-normalised lengths; the
  neighbour-as-centre branch (`:7882`) divides by `skalpro(ww,ww)` with **no `sqrt`** (code-evident; intended form UNCONFIRMED).

Not covered here: `group` (`:1570`), `pref` (`:1572`), `value` (`:2176`); `cd`, `spring`, `lennard` are in MC_FACTS.

---

## 3. `set finish` and the convergence test

### 3.1 Syntax and defaults

- Dispatch at `mmc.f90:1156-1158`; parser `get_finish` `mmc.f90:2559-2631`. Keywords (`:2591`): `feed`, `diff`,
  `rdiff`, `change`, `aver`, `stop`, `log`.
- Defaults **when the command is given** (`:2594-2596`): `feed=3`, four thresholds `99999.`, `stop=cycles`, `log=none`.
- `stop:converge` → `mmc_h_stop=.TRUE.` (`:2602-2603`); `log:screen` → `mmc_h_log` (`:2613-2614`).
- Keyword to variable (`:2624-2628`): `diff`→`mmc_h_conv_m`, `rdiff`→`_r`, `change`→`_c`, `aver`→`_a`.
- Without the command (`discus_blk_appl.f90:288-293`): `mmc_h_stop=.TRUE.`, `MMC_H_NNNN=3`, thresholds `0.0`. With strict
  `<`, `0.0` is never met, so the default never converges early.

### 3.2 Quantities (`mmc_correlations`, `mmc_basic.f90:1395-1464`)

Per feedback, `mmc_h_diff(i) = target_i − achieved_i`, `i = 1..mmc_h_number`:

| conv_val | Definition | Source | Test (strict `<`) |
|---|---|---|---|
| (1) | `max_i |diff_i|` | `:1415` | `diff` |
| (2) | `max_i |diff_i|/|target_i|`, targets ≠ 0 | `:1412` | `rdiff` |
| (3) | `max_i` of mean signed relative change over the history | `:1402-1410`, `:1418` | `|·| < aver` |
| (4) | `max_i` of largest absolute relative change between feedbacks | `:1419` | `|·| < change` |

- `done = (mmc_h_nfeed ≥ MMC_H_NNNN) .AND. (1)<diff .AND. (2)<rdiff .AND. |(3)|<aver .AND. |(4)|<change`
  (`mmc_basic.f90:1462-1464`). With `mmc_m_index=0`, `conv_val = 1e20` (`:1401`, `:1421`): no pass at the first feedback.
- `mmc_h_nfeed` is reset at `mmc.f90:3685` and incremented at `mmc_basic.f90:392` on **every** `mmc_correlations` call.
  With `NNNN=3` the earliest pass is at the **3rd feedback**.
- **Bug (code-evident)**: `mmc_h_aver = mmc_h_aver / mmc_m_index` (`mmc_basic.f90:1410`) sits inside
  `do i = 1, mmc_h_number`, so the whole array is divided once per target. Entry `k` is divided `N−k+1` times.
  `conv_val(3)` is too small for `N ≥ 2`; exact for `N = 1`. `conv_val` is `SAVE` (`mmc_basic.f90:290`); the final call
  (`mmc.f90:3934-3935`) recomputes `done`.

### 3.3 When the test runs

- Serial: `igen = itry` (`mmc.f90:3863`), `imodulus = mo_feed/nnthreads` (`:3829`), feedback when
  `MOD(igen, imodulus·nthreads)==0` (`mmc.f90:4318-4340`). So the test runs every `set feed` iterations (`mo_feed`,
  `discus_blk_appl.f90:92`).
- `done` also ends the loop at `igen > mo_cyc` (`:4336`). The initial call (`mmc.f90:3613`) precedes the reset at
  `:3685` and is not counted.

### 3.4 Does `stop:converge` work for OCC-only runs?

- **Output on** (default): history recorded only inside `IF (lout .AND. mmc_pair<0 .AND. lfirst)` (`mmc_basic.f90:2425-2433`,
  one entry per centre). For feedback calls `lout = lout_feed_in .AND. mmc_out_feed` (`mmc.f90:3608`;
  `mmc_mod.f90:69`, default `.true.`). A **single OCC target** is consistent with the code. Several targets:
  `conv_val(3)` is underestimated (§3.2).
- **Output feed off** (`set output, feed:off`): `mmc_h_ctarg` stays 0, the block at `:1396` is skipped, `conv_val` stays
  at its zero initial value, and `done` (`:1462`) is still evaluated. With thresholds `> 0` the run stops at feedback
  `NNNN`. UNCONFIRMED by execution (MC_FACTS item 13 agrees).
- `log:screen` prints only inside `IF(lout)` (`mmc_basic.f90:1482-1520`), so it also needs feed output on.

### 3.5 Parsing collision on `feed:`

> **Correction (review, 2026-11-03):** `get_optional` in `mmc_set` scans from the last parameter
> backwards and stops with error −12 at the first unknown `name:` (`take_param_mod.f90:131-136`);
> `mmc_set` does not check it, and `get_finish`'s `ber_params` resets it (`berechne_mod.f90:50`). So
> `feed:` is stripped only when it comes after every other `name:` option; written first, it reaches
> `get_finish`. The wrapper writes no `feed:`. The text below describes the simpler (wrong) picture.

- `mmc_set` runs `get_optional` for `type`, `feed`, `final`, `site` (`mmc.f90:1054-1058`, `:1073-1074`) on every `set`
  line, removing matches before dispatch.
- `set finish, feed:N` alone leaves `ianz=1`, so `IF (ianz >= 2 …)` (`mmc.f90:1076`) skips it silently. With another
  keyword (e.g. `stop:converge, feed:5`), `get_finish` never sees `feed:` and `MMC_H_NNNN` stays 3 (code-evident).

---

## 4. Other commands

### 4.1 `set allowed` / `set disallowed` (`mmc.f90:1080-1090`, `:1107-1117`)

- `all` → `mmc_allowed = .TRUE.` / `.FALSE.`; otherwise the listed types (`get_iscat`).
- Every target also sets `mmc_allowed` for its types (`:2796`, `:2807`, `:1928`, `:1937-1938`): place `set disallowed`
  after the last `set target`.

### 4.2 `set limited` (`mmc.f90:1162-1184`)

- `cell, c1,c2,c3, e1,e2,e3` (`:1164-1173`); `atom, lo, hi` (`:1174-1181`); `OFF` (`:1182`, uppercase).
- Affects only shift/invdisp selection (`mmc.f90:5536-5541`; `mmc_limit_selection` `:8932-8949`, which draws `RANDOM_NUMBER`).

### 4.3 `set style` (`mmc.f90:1265-1271`)

- `set style, atom` / `set style, mole` on raw `cpara(2)`, no `ELSE`. `mole` (also `set rota`, `:1211`, and `type:mole`
  on `set move`, `:1200-1215`) switches `mmc_style` to molecule mode; atom-mode wrappers must not emit them.

### 4.4 `set move` (`mmc.f90:1200-1215`; `rmc_set_move`, `rmc_sup_mod.f90:454-560`)

- Atom mode: `set move, <atom|all>, u, v, w, sig` (ianz 5 after `del_params`, `:492-497`) or `set move, <atom|all>, u, v, w`.
- **Atom name must be the first token.** A numeric first token is read as a type index (code-evident, `:497`).
- The 3-number form writes `move(1..3)` only (`:507-508`), so `sig` keeps its previous value. Setting one type twice with
  different sigma leaves a stale sigma (code-evident).
- Units: fractional cell coordinates (`cr_pos`). Shift (`mmc.f90:6136-6152`): `sig==0` → per-axis `gasdev(u)`,
  `gasdev(v)`, `gasdev(w)`; `sig≠0` → one `gasdev(sig)` scaling `(u,v,w)`.
- `mc.py` `Shift` (`src/discuspdf/mc.py:286-310`) writes the 3-number form, so `sig` stays 0 and the per-axis path runs.

### 4.5 `set mode` (`mmc_set_mode`, `mmc.f90:3278-3410`)

- **Bug** (MC_FACTS item 10): `get_iscat(…, werte, …)` at `mmc.f90:3368` overwrites `werte` after `werte(1)` held the
  probability; `mmc_move_prob(imode) = werte(1)` (`:3393`) then stores the **first atom-type index** (`-1` for `all`,
  `0` for void) instead of `<prob>`. `mmc_move_cprob` is rebuilt only when `sump>0` (`:3398-3402`). Only for lists
  (`ianz>2`, `:3366`).
- Arguments (`mmc.f90:1188-1196`): `cpara(1)`=type, `cpara(2)`=site, `cpara(3…)`=atoms (`del_params(2)` at `:3367`).
  Site codes (`:3352-3364`): `A`→all, `L`→local, `SL`→slocal, `SI`→site, `CO`→connectivity; omitted → all (`:3348`).
- `set mode, p, shift, all` has `ianz=2`: no atom list, no bug (`mc.py` relies on this).

### 4.6 `apply_symm` (`mmc.f90:294-296` → `chem_symm_mod.f90:11-243`)

- Uses the average structure `chem_aver` and every space-group operator on each defined vector (`chem_cvec`).
- New vectors are appended to `chem_cvec` and to each CHEM_VEC neighbourhood that used the original (`chem_use_vec`,
  `chem_nvec`, `:200-206`). **It changes the neighbourhoods the next run uses.**
- Only new vectors are appended (`lisnew`, `:170-181`), so a repeat call should add nothing (code-evident; not run).
  An optional file argument writes `set neig, rese` and `set neig, number:next, vec, …` (`:217-220`).

### 4.7 `calc corr` (`mmc.f90:261-275`)

- Sets `lout`, `done`, `lfinished`, `ldetail` = true and `rel_cycl=2.0`, then `mmc_correlations` (`:262-270`).
  `rel_cycl>1` disables feedback (`mmc_basic.f90:2390`, active for `0<rel_cycl<=1`): depths unchanged.
- Side effect: `mmc_h_nfeed`/`mmc_h_index` advance (`mmc_basic.f90:390-392`), reset by a run (`mmc.f90:3685`). `lfeed`
  at that point: UNCONFIRMED.

### 4.8 `grow` (`mmc.f90:289-291`; `mmc_grow` `mmc.f90:9092-9174`)

- Not Metropolis, no energy test. For `i = 1..mo_cyc`: random atom (`:9140`); count neighbours of every type within
  `[0, a0·√3+0.1]` Å (`:9133-9134`); pick the most frequent type, random start index (`:9151-9154`); write `cr_iscat` (`:9161`).
- **Hard-coded to 6 types**: `nneig(6)` (`:9120`), `istart = INT(r1·6)`, `DO j=0,5` (`:9154-9160`). Types above 6 are out of
  bounds. No neighbours → `is_max=0` → invalid type 0 (code-evident).
- Calls `mmc_correlations(…, lfeed=.TRUE., …)` every `mo_feed` iterations (`:9162-9163`), `rel_cycl = i/mo_cyc`:
  **feedback changes depths during `grow`**. `done` is ignored.

---

## 5. Random numbers and reproducibility

- Source: intrinsic `RANDOM_NUMBER`, seeded through `ini_ran_ix` (`lib_f90/random_state_mod.f90`).
  - Seed values given → `RANDOM_SEED(PUT=…)` (`:96-118`) and `iset=0` (`:118`).
  - Single value `0` → default `RANDOM_SEED()` plus time of day (`:95-96`). **Not reproducible.**
- Draws: MMC atom selection and acceptance (`mmc.f90:5014-5541`, selected lines), `mmc_limit_selection` (`:8935-8949`),
  `grow` (`:9140`, `:9151`), shift via `gasdev` (`:6142`, `:6148`).
- `gasdev` (`lib_f90/random.f90:16-45`): Marsaglia polar method on `RANDOM_NUMBER` (`:31`, `:33`). Rejects `r≥1` (`:38`),
  so the number of uniforms per call varies. Caches the second deviate in `SAVE` variables (`:20`, `:24`, `:40-45`).
  `iset` is in `random_mod.f90:13`, reset only by seeding (`random_state_mod.f90:118`).
- Consequence: a shift changes the RNG stream length, so different move mixes do not share a stream even with the same
  seed. Deterministic for a fixed macro and seed.
- No RNG use: `set neig`, `set target`, `set finish`, `set valid` (no-op, MC_FACTS item 3), `calc corr`, `apply_symm`,
  `set allowed`. OpenMP: MC_FACTS §5 (not reproducible with more than one thread).

---

## 6. Wrapper safety (atom mode), closing table

| Feature | Safe to wrap in atom mode? | Why |
|---|---|---|
| `set neig, vec` (+ `range`, `number:`) | yes | MC_FACTS §1; `range` lowercase (`chem.f90:2087`) |
| `set neig, con` + `set con` | with caveat | needs CONNECT data (`appl_dis.hlp:34`); swchem effect on connectivity lists UNCONFIRMED |
| `set angle` + `set neig, ang` | with caveat | `angle` target only; site numbers, not types (`chem.f90:1509`) |
| `set env` + `set neig, env` | with caveat | `-1` centre = any; distances in Å |
| `set neig, dis` | with caveat | never for `angle` (`mmc.f90:8599-8600`); `angsig` needs `S` (`chem.f90:1959-1972`) |
| `set range` (`direc` only) | with caveat | other sub-keywords UNCONFIRMED |
| `corr`, one per `ic`, `CO`/`EN` explicit | yes | `mmc_set_disp_occ` (`mmc.f90:2765-2829`) |
| `corr` with `(A,B),(C)` lists | with caveat | a second `corr` on the same `ic` erases the first (`:2787-2791`) |
| `unid` | with caveat | no array clear (`:2833`); directed; `ABS(depth)` sign UNCONFIRMED |
| `cn` | with caveat | one-way pairs (`:2724-2732`); `all` dropped; no feedback; centre-only energy (`:7526`) |
| `angle` | with caveat | centre type ignored (`:1945`); 0° target inert (`:8612`); no DIST/RANGE (`:8599-8600`) |
| `bucking` | with caveat | eight numbers, `werte(1)` unused; `werte(5)` must be non-zero (`:8523`) |
| `repulsive` | with caveat | averaged; constant when no neighbours (`:8444`) |
| `heli` | with caveat | always `CO` or `EN`, else depth overwritten (`:1724-1725`); normalisation mismatch (`:7826` vs `:7882`) |
| `set finish, stop:converge` | with caveat | thresholds explicit; feed output on (§3.4); `feed:` collision (§3.5); `aver` bug for N≥2 (§3.2) |
| `set finish` defaults only | yes | thresholds 0 never converge (`discus_blk_appl.f90:290-293`) |
| `set output, feed:off` with converge | no | breaks the convergence history (§3.4) |
| `set allowed` / `set disallowed` | yes | `disallowed` after the last target (§4.1) |
| `set limited` | yes, shift only | lowercase `cell`/`atom`, uppercase `OFF` (`mmc.f90:1164`, `:1174`, `:1182`) |
| `set move, <atom>, u,v,w[,sig]` | yes | atom first; 3-number form keeps old sigma (`rmc_sup_mod.f90:492-508`) |
| `set mode, p, shift, all` / `swchem, all\|local\|slocal\|site` | yes | no atom list (`mmc.f90:3366`) |
| `set mode, p, <type>, <site>, <atoms…>` | **no** | probability overwritten (`mmc.f90:3368`, `:3393`) |
| `apply_symm` | no | appends vectors and rewires neighbourhoods (§4.6) |
| `calc corr` | yes | read-only w.r.t. depths (`rel_cycl=2.0`) |
| `grow` | no | not Metropolis; 6-type hard-coding; feedback side effect (§4.8) |
| `set rbn, growth` | no | undocumented; not run by default (MC_FACTS §2) |
| `set style, atom` / `set style, mole`, `set rota`, `type:mole` | atom: yes; the rest: no | `mole` and `rota` switch to molecule mode (§4.3) |

---

## 7. UNCONFIRMED (not verified by execution)

1. Sign and `ABS(depth)` conventions in `mmc_energy_uni` for each start/end case (`mmc.f90:7068-7138`).
2. `cn` achieved `n_cn` definition (`mmc_basic.f90:917-944`), and whether `tatom=.FALSE.` neighbours enter `n_neig` (§2.3).
3. `set range` `central`/`neig`/`short` behaviour (`chem.f90:1345-1460`).
4. Whether the partial counter reset in `set neig, rese` causes wrong results later (§1.2).
5. Whether connectivity lists follow species after `swchem` (§1.5, `chem_neig_multi.f90:361-401`).
6. Output-off convergence: early stop after `NNNN` feedbacks with thresholds `> 0` (§3.4).
7. Whether the division bug at `mmc_basic.f90:1410` changes published results (§3.2, N≥2).
8. Whether `mmc_h_number` and `mmc_h_ctarg` stay aligned for multi-pair `corr` (§3.2).
9. Which `heli` normalisation is intended, `mmc.f90:7826` or `:7882` (§2.7).
10. `apply_symm` idempotence (§4.6).
11. Effect of `lfeed` on a `calc corr` made before any run (§4.7).
12. Whether invalid trial draws advance `igen = itry`, and so count toward `set feed` (§3.3).
13. Mixing `unid` (no clear) and `corr` (clears) on one `ic` (§2.1, §2.2).
14. Error path for unknown names inside `(A,B)` lists (`get_iscat`).
