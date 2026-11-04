# DISCUS 6.23.01: powder menu and powder output facts

Scope: the `powder` menu (`set`, `run`, `show`, `reset`, `xray`, `neutron`, `electron`, `asym`, `chemistry`)
and the powder parts of the `output` menu (`form powd`, `form pdf`, `value`, file writing).
Written 2026-11-03 for the Python wrapper (`src/discuspdf/pdf.py`). Source reading only; **no DISCUS run was
made for this document**. Behaviour inferred from code is marked "derived". Anything not confirmed is marked
**UNCONFIRMED** and listed in section 10.

## 0. Conventions

- `prog/` = `/home/claude/dc_6.23.01/discus/prog/`; `lib/` = `/home/claude/dc_6.23.01/lib_f90/`.
- Line numbers refer to the tree as shipped. The command catalogue (`src/discuspdf/data/suite_catalogue/discus_suite_6.23.01.json`,
  source commit `19a214a`) lists the same `powder.f90` dispatch lines (126-233) and `set` help at `appl_dis.hlp:7435`.
- Keyword "j" = minimum abbreviation length passed to `str_comp` (section 1.3).
- Numbers below are at the wrapper's usual settings unless stated.

## 1. Parsing rules (all powder and output commands)

### 1.1 Command word
- `get_cmd` (`lib/sup.f90:12`) takes the first word up to a blank or tab as `befehl` (`sup.f90:189-195`).
- No case folding is performed in `get_cmd` or in `do_powder`. Whether the macro/line reader folds case upstream: **UNCONFIRMED**.

### 1.2 Arguments
- `get_params` (`lib/get_params_mod.f90:14`) splits on commas that are not inside quotes or `()`/`[]`.
  Leading and trailing blanks are removed. Case is not changed. Empty fields are an error (`ier -2`, `:132-138`).
- Count limit: more fields than the array allows gives `ier -17` (`:73-79`).

### 1.3 Keyword matching (`str_comp`, `lib/str_comp_mod.f90:6-57`)
`str_comp(token, keyword, j, la, lb)` returns true iff `la >= j` and the first `min(la, lb)` characters are equal (`:33-52`).
Consequences:
- Prefix match. `qmin`, `q` and `qmax` differ only through `j`.
- An over-long token matches a keyword it starts with: `qmaxX` matches `qmax` (`lb` caps the compared length).
- The first matching `ELSEIF` wins. Order matters: `pr` matches `preferred` (`powder.f90:1060`) and not `profile` (`:1112`).
- Comparison is Fortran `.eq.`, so it is case-sensitive. Upper-case `Q` will not match `q` in `form powd`.

### 1.4 Errors
- Unknown `set` keyword: `ier_num = -8` (`powder.f90:1478-1479`).
- Wrong argument count for most `set` keywords: `ier -6` (e.g. `dtth`, `:910-932`).
- Several options accept an unknown second argument **without any error and without any change** (`calc`, `fourier`, `temp`, `axis`; section 3).
- `set refine` appears in `appl_dis.hlp` but has no branch in `do_pow_set` (`powder.f90:741-1482`); it gives `-8`.

## 2. Menu commands (`do_powder`, `powder.f90:14-297`)

| Command (j) | Dispatch | What it does | Notes |
|---|---|---|---|
| `asym` (2) | `:126-128` | lists asymmetric unit | body not traced (UNCONFIRMED details) |
| `chemistry` (2) | `:136-138` | lists atoms in the crystal | body not traced |
| `continue` (2) | `:131` | continues a macro | |
| `echo`, `evaluate`, `exit`, `help`, `?`, `system`, `wait` | `:139-230` | standard | `help` maps to `discus powd` (`:163`) |
| `electron` (1) | `:168-172` | `lxray=.true.`, `diff_radiation=RAD_ELEC`, `diff_table=RAD_INTER` | electron form-factor path not traced (UNCONFIRMED) |
| `neutron` (1) | `:175-179` | `lxray=.false.`, `RAD_NEUT`, `RAD_INTER` | neutron scatter = `b_coh` (`fourier_sup.f90:2976+`, `dlink`); `form()` returns `scat(1)` only (`fourier_sup.f90:2947-2971`). Does **not** change `lpcor` (stays Bragg). |
| `reset` (2) | `:182-184` | `CALL powder_reset` only | see section 4 for what is not restored |
| `run` (2) | `:187-190` | `guess_atom_all`, then `powder_run(zeile,lp)` | `zeile` (`mode:`, `current:`, `fraction:`) goes to `phases_set` (`powder.f90:335`); help `appl_dis.hlp:7407` |
| `show` (2) | `:193-210` | `dlink`, then `pow_conv_limits`, `pow_show` (`powder.f90:451`) | `show scat,<atom>` prints the scattering curve (`:197-208`) |
| `set` (2) | `:213-215` | `do_pow_set` (`:741`) | section 3 |
| `xray` (1) | `:233-247` | `lxray=.true.`, `RAD_XRAY`, `diff_table=RAD_INTER`; then table | default table `opara(O_TABLE)='waas'` is set on each menu entry (`:80`); `if(opara=='waas')` gives `RAD_WAAS` (`:242-243`); `'discamb'` gives `RAD_DISC` (`:244-246`). The help spells `waass` (`appl_dis.hlp:7847`); the code compares `'waas'` exactly, so `table:waass` leaves `RAD_INTER`. |

Unknown command in the powder menu: `ier -8` (`powder.f90:252`).

## 3. `set` subcommands (`do_pow_set`, `powder.f90:741-1482`)

Start value = value after `reset` (section 4), unless stated. "Silent" = accepted, no error, no change.

| Keyword (j) | Syntax as parsed | Start value | Effect in code | Refs | Pitfalls |
|---|---|---|---|---|---|
| `axis` (2) | `set axis,<x>` | Q (`powder_mod.f90:46`) | sets `POW_AXIS_Q` **regardless of `<x>`** | `:796-797` | `set axis,tth` silently gives Q. Help says `"q"\|"tth"` (`appl_dis.hlp:7440`). |
| `back` (2) | `set back,p0[,p1..]` | 0 terms | `pow_back(0..n)`; `pow_back(-1)` origin | `:798-821` | Applied only to `inte` (`powder_write_mod.f90:468-478`). `pow_back(-1)` is not reset. |
| `bragg` (2) | `set bragg,incl\|excl` | include | `pow_l_all` | `:822-835` | Used only in complete mode (`powder.f90:2347`). Prefixes `incl`/`excl` accepted. |
| `dispersion` (2) | `set dispersion,anomalous\|off` | off (`fourier_reset.f90:70`) | `ano`; anomalous f' f'' from `cr_delfr/cr_delfi` | `:836-849`; `fourier_sup.f90:2730-2736` | Not reset by `powder reset`. Neutron unaffected. |
| `corrlin` (5) | `set corrlin,<a>[,<rcut>]` | 0 | `pdf_clin_a`, `pdf_rcut` | `:850-866` | Used only by the correlated-motion PDF branch (`powder_write_mod.f90:868-1063`). |
| `corrquad` (5) | `set corrquad,<a>[,<rcut>]` | 0 | `pdf_cquad_a` | `:867-884` | As above. Also doubles Qmax in `powder_run` (`powder.f90:352`). |
| `delta` (2) | `set delta,<res>` | 0.0 (`powder.f90:3589` reset) | `pow_delta=res`, **and forces `pow_profile=GAUSS`** (`:904`) | `:897-909` | Gaussian convolution only if `res>0` (`powder_write_mod.f90:1375`). `res>0` aborts on most grids (issue #65). With `res=0` the profile stays Gaussian and no convolution is done: this is why the wrapper's `set delta,0.0` disables the PSVGT default (section 4). |
| `dtth` (2) | `set dtth,<v>` | 0.05 deg | `pow_deltatth`, `pow_deltaqtth=.false.`, `pow_constlam=.true.`; range 0..180 | `:910-933` | Angle-axis instrument mode. |
| `dq` (2) | `set dq,<v>` | 0.001 Å⁻¹ | `pow_deltaq = NINT(v·1e5)/1e5`: **rounded to 1e-5 Å⁻¹**; `pow_deltaqtth=.true.` | `:934-947` | `0.0123456` is stored as `0.01235`. |
| `dh`, `dk`, `dl` (2) | `set dh,<v>` etc. | 0.05 r.l.u. (`powder_mod.f90:103`) | `pow_hkl_del(1..3)` | `:948-986` | Used by `calc complete` only (`powder_complete`, `powder.f90:1936-2030`). |
| `energy` (2) | `set energy,<E>` | 17.480782 keV (`fourier_reset.f90:83`), `l_energy=.false.` | `renergy=E`, `lambda=' '`, `l_energy=.true.`. Units: keV (X-ray, electron), meV (neutron). | `:987-999`; `lib/element_data_mod.f90:3005-3060` | X-ray: λ(Å)=12.398419/E(keV) (derived; reproduces default 0.709260 for 17.480782). Neutron: λ(Å)=9.048/√E(meV) (derived). |
| `ka21` (4) | `set ka21,<ratio>` | 0.0, off | `pow_ka21`, `pow_ka21_u=.true.`; applied in `pow_k12` on output (`powder_write_mod.f90:221`, `:747`) | `:1002-1015` | Output only. Formula in `pow_k12` (`powder_write_mod.f90:3049`) not traced. |
| `partial` (3) | `set partial,<i>,<j>` or `set partial,off` | off | `pow_l_partial`, `pow_ipartial(1..2)`; Debye histogram restricted to that pair (`powder_pdf_hist.f90:349-352`) | `:1016-1031` | PDF output routed to `powder_out_partial` (`graphic.f90:568-569`). Normalisation UNCONFIRMED (section 6.7). |
| `period` (4) | `set period,<d>` or `set period,"off"` | off (`pow_lperiod=.false.`, `powder_mod.f90:93`) | Debye: finite sphere cut of radius `d/(2·0.8)` (`powder_pdf_hist.f90:80-101`, `powder_mod.f90:95`); PDF divided by shape function `1-1.5x+0.5x³` (`powder_write_mod.f90` default PDF branch, `:1123-1131`) | `:1034-1059` | **Not reset by `powder reset`** (section 4). Structure restore after the sphere cut UNCONFIRMED (`powder_pdf_hist.f90:91-92` saves `internal.powderback.stru`). Value `0` turns it off (`:1044`). |
| `preferred` (2) | `set preferred,off\|rietveld\|march\|g1,v\|damping,v\|g2,v\|portion,v\|hkl,h,k,l` | off | `pow_pref`; applied to `|csf|²` in complete mode only (`powder.f90:644`) | `:1060-1111` | `pr` is taken by this keyword, not by `profile`. |
| `profile` (2) | `set profile,off\|gauss\|pearson\|pseudo\|eta,..\|uvw,u,v,w[,fwhm:..]\|asym..\|tof..\|bank..` | PSVGT, `u=v=0`, `w=0.05`, `eta=0.5` (`powder.f90:3609-3617` reset) | `pow_profile`, `pow_u/v/w`, `pow_eta` | `:1112-1113` and `pow_set_profile` `:1486-1780` | `off` sets `pow_profile=0` and `pow_delta=0` (`:1575-1578`); this is the only switch that turns off the PSVGT broadening (section 4). |
| `scale` (2) | `set scale,<s>` | 1.0 | `pow_scale`: multiplies `inte`, `F(Q)`, `PDF` (not `S(Q)`) | `:1114-1125`; `powder_write_mod.f90:477, 489, 497, 1135` | |
| `sh`, `sk`, `sl` (2) | `set sh,<v>` etc. | 0 | `pow_hkl_shift(1..3)` | `:1126-1163` | Complete mode only. |
| `stepr` (2) | `set stepr,<v>` | 0.01 Å (`debye_mod.f90:20`) | `pow_del_hist`: pair-distance bin width of the Debye histogram | `:1150-1161`, `powder_pdf_hist.f90:418, 443` | **Not reset by `powder reset`**. Sets the r-accuracy of every Debye run (section 5.2). |
| `calc` (1) | `set calc,complete\|debye\|nufft\|grid` | complete (`powder_mod.f90:50`) | `pow_four_type` = 0/1/2/3 (`powder_mod.f90:15-18`) | `:1177-1195` | One-letter abbreviations work (`d`=debye, `c`=complete, `n`=nufft, `g`=grid). Unknown second argument: **silent**. `fast` does not exist (no branch). |
| `fourier` (1) | `set fourier,fourier\|stack` | fourier (`pow_four_mode=0`) | `POW_STACK` for stacking faults | `:1196-1210` | Unknown second argument: silent. |
| `lpcor` (1) | `set lpcor,none\|bragg[,ang]\|neutron\|synchrotron,f[,ang]\|corre[,ang]\|tof` | bragg, `fac=0.88`, `ang=20` (`powder_mod.f90:53, 88`; reset) | `pow_lp`, `pow_lp_fac`, `pow_lp_ang`; `polarisation()` (`phases_set_mod.f90:1008`) | `:1212-1283` | `bragg` with no angle gives `fac=0.5`, `ang=0` (`:1225-1230`). Polarisation enters `inte` (`phases_set_mod.f90:456, 466`), and the `F(Q)`/`PDF` path does not use it. Unknown option: `-90`. |
| `temp` (1) | `set temp,use\|ignore` | ignore (`fourier_reset.f90:71`, `ldbw=.false.`) | `ldbw`: Debye–Waller `exp(-B·s²)` in `cfact` (`fourier_sup.f90:2730`) | `:1286-1298` | Unknown second argument: silent. **Not reset by `powder reset`**. |
| `tthmin` (5), `tthmax` (5), `tthzero` (4) | `set tthmin,<deg>` | 0.1, 40.0, 0 | `pow_tthmin/max`, `pow_qtthmin/max=.false.`, `pow_constlam=.true.`; converted by `pow_conv_limits` | `:1299-1358`; `powder.f90:3423-3492` | Angle-axis instrument mode. |
| `tmin` (4), `tmax` (4), `tstep` (5) | TOF time limits | | `pow_constlam=.false.` | `:1361-1406` | TOF mode. |
| `qmin` (4) | `set qmin,<v>` | 0.2 Å⁻¹ | `pow_qmin`, `pow_qtthmin=.true.`, `pow_constlam=.true.` | `:1421-1434` | Use multiples of `dq` (section 6.4). |
| `qmax` (4) | `set qmax,<v>` | 7.0 Å⁻¹ | `pow_qmax`, `pow_qtthmax=.true.` | `:1407-1420` | Must satisfy `qmax·λ/(4π) <= 1` (`powder.f90:1816-1819`). |
| `qzero` (4) | `set qzero,<v>` | 0 | `pow_qzero`, `pow_qtthzero=.true.` | `:1435-1449` | Output x-array `xpl = q - qzero` (`powder_write_mod.f90:247-249`). Effect on written `F(Q)` UNCONFIRMED. |
| `wvle`/`wvlength` (1) | `set wvle,<λ>` or `set wvle,<symbol>` | `MOA1`, 0.709260 Å (`fourier_reset.f90:66, 84`) | numeric: `rlambda=λ`, `lambda=' '`, `l_energy=.false.`; a symbol is upper-cased and stored in `lambda` | `:1452-1480` | Symbol table `per_symwl` values not extracted (UNCONFIRMED). `wvle` alone with prefix `w` also matches. |

Additional `set` notes:
- `set calc,debye` is the only way the wrapper switches the algorithm: `pow_four_type` defaults to `POW_COMPL` (`powder_mod.f90:50`, reset `powder.f90:3552+`).
- `set disp` with `off` is `ano=.false.` (`:836-849`).
- `show` prints the settings (`pow_show`, `powder.f90:451-737`), including `DELTA Q` (`:1291`) and `Profile` lines.

## 4. Defaults and what `reset` restores

Start-up chain: `discus.f90:13` → `discus_setup.f90:49` → `discus_reset_all.f90:54` (`fourier_reset`, which sets
`lambda`, `ldbw`, `ano`, `renergy`, `rlambda`, `diff_power=4`) and `:60` (`powder_reset`).

The powder-menu `reset` calls `powder_reset` only (`powder.f90:182-184`, `:3552-3645`).

| Setting | Restored by menu `reset`? | Where it is set | Note |
|---|---|---|---|
| Q-limits, dq, axis, calc, profile (→PSVGT, w=0.05), hkl steps, scale, back, preferred, lpcor (fac 0.88), partial | Yes (except below) | `powder.f90:3552-3645` | `pow_profile` returns to PSVGT, so **a `set delta` before `reset` does not persist** |
| Wavelength, energy/λ, radiation, table | **No** | `fourier_reset.f90:66-85` | `xray`/`neutron` must be re-issued after `reset` if the radiation changed |
| `ldbw` (`set temp`) | **No** | `fourier_reset.f90:71` | persists |
| `ano` (`set dispersion`) | **No** | `fourier_reset.f90:70` | persists |
| `stepr` (`pow_del_hist`) | **No** | `debye_mod.f90:20` | persists (section 7, L10) |
| `period` (`pow_lperiod`, `pow_period`) | **No** (not in reset block) | `powder_mod.f90:93-94` | persists |
| `back_q`, `pow_back(-1)` | **No** | `powder_mod.f90:98, 99` | `pow_nback` and `pow_back(:)` are reset |
| `pow_pr_fwhm` | **No** | `powder_mod.f90:113` | `pow_pr_par` is reset |

Module versus reset values that differ: `pow_ds_max` `0.0001` (`powder_mod.f90:85`) versus `0.001` after reset (`powder.f90:3587`); both are used only in complete mode.

## 5. `run`

### 5.1 Common steps (`powder_run`, `powder.f90:301-411`)
1. `dlink` → `rlambda`, `lambda`, form-factor table (`powder_run`, `:301-310`); `phases_set(zeile)` (`:335`).
2. `pow_conv_limits` (`:336`) converts angle limits to Q. Saves user values to `pow_*_u`.
3. Widens Qmax for the convolution: `pow_qmax += min(5, 2·pow_width·fwhm + 2|qzero|)` (`:350`), doubled if `corrlin/corrquad` is set.
4. Caps Qmax: `pow_qmax = min(pow_qmax, 4π/λ − 1e-4)` (`:356`). Note this is `4π/λ`; the hard error in `powder_qcheck` is `qmax·λ/(4π) > 1` (`:1816`).
5. `check_dq` (`:361`): **doubles `pow_deltaq` silently** until `(qmax−qmin)/dq+1 <= 2^18`.
6. `powder_qcheck` (`powder.f90:1784-1825`): errors if `λ < tiny`, `qmax <= qmin`, `dq <= 0`, or `qmax·λ/(4π) > 1` (`-108`).
7. Dispatch on `pow_four_type` (`:383-390`): debye → `pow_pdf_hist`; complete → `powder_complete`; nufft/grid → `powder_nufft`.
8. `powder_run_post` (`:415-447`): `powder_convolute` (`powder_write_mod.f90:1271`) then `phases_place` (`powder.f90:437`).
9. Restores user `pow_qmin/qmax/deltaq/tth*` (`:404-409`), so output uses the user Q-range.

### 5.2 Debye path (`pow_pdf_hist`, `powder_pdf_hist.f90:14-154`; `powder_debye_hist_cart`, `:204-697`)

Structure used: **all atoms of the crystal as read** (the supercell), in Cartesian coordinates, as a **finite cluster**
with no periodic images (derived from the `powder_trans_atoms_tocart` call in `pow_pdf_hist`; `cr_natoms` loops at `:427` (parallel) and `:469` (serial)). Atoms with
`cr_iscat(1,j) <= 0` are skipped. Periodic PDF only with `set period` (sphere cut, section 3).

Q grid (derived):
- `xm(1)=pow_qmin/2π` (`:294`), `uin(1)=pow_deltaq/2π`, `ss=(pow_qmax+pow_qmax_buf)/2π` with `pow_qmax_buf=0.5` (`:295`, `powder_mod.f90:84`).
- `num(1) = NINT((ss-xm)/uin)+1` (`:298`). Grid point `k` (1-based) is `Q_k = pow_qmin + (k-1)·pow_deltaq`.
- Here `pow_qmax` is the widened value from 5.1 step 3.

Pair histogram (derived):
- `distance = |udist|`, `n_hist = NINT(distance/pow_del_hist)+2` (`:305`).
- For every unordered atom pair `j<l`: `bin = INT((r + 0.5·Δ)/Δ)`, `Δ = pow_del_hist` (`:443`, `:485`). So `bin = round(r/Δ)` and the distance used is `r̃ = bin·Δ`.
- Histogram is kept per species pair `look(i,j)` (`:364`). Bin 0 must be empty: if any entry is in bin 0 (r < Δ/2) the run stops with `ier -123` (`:533-534`).

Sum (derived; `:587`, `:606`, `:623-665`):
- `sinc_table(x) = sin(2π·INT(x·I2PI/2π) mod 2^16 / 2^16)`: a 65536-entry sine table (`powder_tables.f90:13-43`, `diffuse_mod.f90:20-21`).
  Its argument is `r̃·Q`, truncated to multiples of 2π/65536 (floor).
- `partial(k,pair) = Σ_bins H_pair(bin) · sin_table(r̃·Q_k) / (r̃·Q_k)` (`:587-589`), using `Q_k` as above.
- `rsf(k) = Σ_pairs 2·partial(k,pair)·[Re f_i Re f_j + Im f_i Im f_j] + Σ_species natom_s · |f_s(Q_k)|²` (`:655`, diagonal term in the same routine), with `f` from `cfact_pure` (no DW) when `ldbw` is off, or from `cfact` (with `exp(-B s²)`) when `ldbw` is on (`:623`).
  `cfact_pure = f + f' + i f''`; `cfact = occ·exp(-B·s²)·(f + f' + i f'')` for isotropic B (`fourier_sup.f90:2730-2736`). `s² = (Q/4π)²` (`four_formtab`, `:2691`; table index `nint(s/CFINC)`, `CFINC=0.001`, `discus_config_mod.f90:28`).
- The factor 2 turns the unordered-pair histogram into ordered pairs. The diagonal (self) term is separate.
- Form factor: `f(s) = c + Σ_{i=1..4} a_i·exp(-b_i·s²)` for X-ray (`fourier_sup.f90:2947-2971`, `diff_power=4`); `f = b_coh` for neutron.
- `pow_nreal = Σ NINT(N_species·occ)` (`:516`).

Result arrays used later (derived):
- `rsf(1:num(1))` → `pow_conv(1:npkt)` (`powder_write_mod.f90:1351`). `pow_conv(0)` is not set in the Debye branch (it stays 0): see the one-ΔQ shift in 6.4.
- `pow_faver2`, `pow_f2aver`, `pow_f2(i,species)`, `pow_u2aver`, `pow_nreal` from `powder_f2aver` (`powder_write_mod.f90:2969`).
- `pow_nn_partial = natom` (`:682`).

### 5.3 Complete path (`powder_complete`, `powder.f90:1829-2593`)
- Loops over an hkl grid with steps `pow_hkl_del` (default 0.05 r.l.u.) inside `|h_i| <= a_i·pow_ds_max` (`:1936-1938`, `:2011-2015`).
- For each point with `q = 2π·|h*|` in `[pow_qmin, pow_qmax+pow_deltaq]` (`:2377`): `inten = |csf|²·scale_diffuse` (`:2357`), optional `calc_preferred` (`:644`, `:2360-2363`), then `powder_split(... pow_qsp)` into the Q bins (call sites `:2391`, `:2527`, `:2954`; definition not read).
- `pow_conv(k) = 2·pow_qsp(k)·(1/xdel)` (`powder_write_mod.f90:1355`, `:1436-1437`).
- `csf` construction, `scale_diffuse`, and the weights of `powder_split` are **UNCONFIRMED** (not read).
- `pow_sq` for complete mode includes a Lorentz factor (`phases_set_mod.f90:441-455`, `lorentz(ttheta,q,...)`) and is normalised by `sum(pow_sq)/sum(pow_fu)` (`:511`). The Debye and complete S(Q) definitions therefore differ.

### 5.4 Convolution (`powder_convolute`, `powder_write_mod.f90:1271-1451`)
- `pow_profile == 0` → no convolution (`:1362`).
- `GAUSS` with `pow_delta > 0` → `powder_conv_res` (`:1375`, `:1599`) (issue #65 region).
- `PSVGT` (default: `w=0.05`) → `powder_conv_psvgt_uvw` (`:1945`), **including the Debye branch** (`pow_type==POW_DEBYE` at `:2025`). The Debye pattern is therefore broadened by default. Size of the effect: UNCONFIRMED (section 7, L11).
- Debye convolution step is `xdel = (qmax−qmin)/num(1)` (`:1329`), which is not `pow_deltaq` in general (UNCONFIRMED impact).
- Debye scale: `scalef = Σpow_tmp/Σpow_uuu` (`:1436`), i.e. the convolution preserves the sum.

## 6. Output (`output` menu, `graphic.f90`)

### 6.1 `form` (`graphic.f90:218`, dispatched by `format`, min 1 character)
- `form powd,{tth|q|tof}[,start,end,step]` (`:245-265`). Needs 5 fields for limits (`:245-265`). Step rounded to 1e-5 (`:274`).
  Token `q` matches any token starting with `q` (`str_comp(...,'q',1,...)`, `:251-253`).
- `form pdf,r[,rmin,rmax,dr]` (`:283-300`). Needs 5 fields for limits (`:283-300`); no rounding. Without limits the `pdf_rminu/rmaxu/deltaru` module values are used (`pdf_mod`, UNCONFIRMED values).

### 6.2 `value` (`graphic.f90:782-840`, match is case-sensitive)
| Token | Code | Meaning (derived) |
|---|---|---|
| `inte` / `in` (`:796`) | `val_inten`=1 (`qval.f90:3`) | intensity: `pow_conv` (with polarisation), times `pow_scale`, plus background |
| `S(Q)` (`:819`) | `val_sq`=7 (`qval.f90:9`) | `pow_sq` (section 6.3) |
| `F(Q)` (`:822`) | `val_fq`=8 | `Q·(S(Q)−1)` |
| `I(Q)` (`:816`) | `val_iq`=13 | `pol·rsf/N + lpv·(⟨f⟩²−⟨f²⟩)` (6.3) |
| `f2aver` (`:829`) / `faver2` (`:831`) | 9 / 10 | `⟨f²⟩` / `⟨f⟩²` |
| `PDF` (`:837-838`) | `val_pdf`=14 (`qval.f90:16`) | G(r) (6.5) |
| `<` prefix | `laver=.true.` | averaged variant |
| `f2averb`, `faver2b`, `faverb` | 17, 18, 19 | UNCONFIRMED (not traced) |
| `pdf` (lower case) | falls into `'p'` branch: `ianz==1` gives `val_phase` (`:800-807`) | **a lower-case `pdf` silently selects the phase** (derived; not run) |

### 6.3 Formulas for S(Q), F(Q), I(Q) (Debye, default `ldbw=.false.`, `deb_conv=.false.`)
Per-atom weights (`powder_f2aver`, `powder_write_mod.f90:2969-3047`): `w_s = N_s·occ_s / N_real`, `N_real = Σ NINT(N_s·occ_s)`.
- `⟨f²⟩(Q) = Σ_s w_s·|f_s(Q)|²` (from `cfact_pure`, no DW). `pow_f2aver`.
- `⟨f⟩(Q) = Σ_s w_s·sgn·|f_s(Q)|`, with `sgn` = sign of `Re f_s(0)`; `pow_faver2 = ⟨f⟩²` (`:3023`, squared at `:3035`). This uses `|f|`, not the complex `f`, so it is not exact for anomalous scattering.
- Per-atom pattern: `pha_powder = pow_conv/pha_nreal` (`phases_set_mod.f90:254`).
- `pow_conv` (inte path) = `Σ pha_scale·pha_powder·polarisation(2θ)` (`phases_average`, `:466-471`).
- `pow_sq` = `Σ pha_scale·pha_powder` (no polarisation) (`phases_average`, `:476-478`).

Then (`phases_average`, `:540`, `:559-574`):
- **S(Q) = pow_sq/⟨f⟩² + 1 − ⟨f²⟩/⟨f⟩²**, then divided by `sq_aver`, the mean of S over `Q > 18 Å⁻¹` if more than 50 points exist, else 1 (`:559-574`).
- **F(Q) = Q·(S(Q) − 1)** (`powder_write_mod.f90:270-273`).
- **I(Q) = pow_conv + lpv·(⟨f⟩² − ⟨f²⟩)** for the branch `deb_conv .or. .not.ldbw` (`:276-283`); `pow_conv` already carries the polarisation factor; `normalizer=1` (`:263`). With `ldbw` on, the second term is `lpv·(⟨f⟩² − exp(−q²·u²)·⟨f²⟩)` with `u² = Σ w_s·B_s/(8π²)` (`:282`, `:290-293`; `powder_f2aver:3033-3036`).
- Single species (derived): `S(Q) = 1 + (2/N)·Σ_{i<j} sin(Q r̃_ij)/(Q r̃_ij)` (the `⟨f⟩²` cancels), `I(Q) = pol(2θ)·f²·S(Q)`.

Complete path S(Q): Lorentz-weighted (`phases_set_mod.f90:441-455`) and normalised differently (`:511`). Not analysed further (UNCONFIRMED).

Scale and background: `pow_scale` multiplies `inte`, `F(Q)`, `PDF`; `pow_back` only `inte`; neither applies to `S(Q)` (`powder_write_mod.f90:468-497`).

### 6.4 Grid and the one-ΔQ label shift
- Output grid (`powder_out`, `:34-581`): `xmin=pow_qmin_c`, `xdel=pow_deltaq_c`, `npkt=NINT((xmax+xdel−xmin)/xdel)` (`:154`), with `xmax = pow_qmax_c` (widened). Labels `xpl(j) = xmin + j·xdel` (`:243-249`).
- Then a cubic spline (`spline`/`splint`, `:427`) onto the user grid `qmin + i·deltaq` (`:429`). If the user grid equals the calc grid and is aligned (`qmin` multiple of `dq`), the spline returns the nodes.
- **Label shift (derived):** the Debye value `rsf(k)` is computed at `Q_k = qmin + (k−1)·dq`, but `powder_convolute` stores it at index `k` (`:1351`), which is labelled `qmin + k·dq`. Index 0 is a placeholder (S=0, so F=−Q). **Written S(Q), F(Q), I(Q) at label `Q_j` are the values for `Q_j − dq`.** The `F(Q)` labels being one ΔQ high (VERIFIED_DISCUS_MACRO.md) follows from this. The PDF transform compensates by dropping `ywrt(0)` (`discus_powder_fft.f90:68`).
- `form powd,q,start,end,step`: the written grid is the user grid (`powder_out` Q branch, `:388-430`), not the calc grid.

### 6.5 PDF (`powder_out_pdf`, `powder_write_mod.f90:585-1227`; default branch `:1064-1113`)
- `F(Q) = Q·(S−1)` at the output grid (`:792-797`, derived).
- Termination (`zero_last`, `:1081`): from the top of the grid downward, zero every point until the last sign change; that point is also zeroed.
- Transform: `fft_fq` (`discus_powder_fft.f90:7-112`): a DST-II (`ddst`, FFTPACK) of `F(Q_i)` on a zero-padded grid of `npkt_fft = 2^18` points (`powder_write_mod.f90:122`). Points below `qmin` are filled by `F = F_1·Q/Q_1` (`:61-66`). Samples `(2/π)·dq·F` are cubic-spline interpolated to the r grid (`:90`, `spline_prep`).
- `yfour = yfour·π/2` (`:1102`). Net: **G(r) = dq·Σ_i F(Q_i)·sin(Q_i·r)** over the zeroed grid. This is π/2 times the textbook `(2/π)∫F sin(Qr) dQ`, which matches VERIFIED_DISCUS_MACRO.md.
- r grid: `npkt = NINT((rmax−rmin)/dr)+1` (`:1074`); `rstep = (rmax−rmin)/(npkt−1)` (`:1076`), so the written step differs from `dr` if `(rmax−rmin)/dr` is not an integer.
- Periodic division by the shape function (`:1123-1131`) only with `set period`.
- Scale `pow_scale` (`:1135`).
- Lorch or other window functions: none found in `powder_write_mod.f90` (grep). The only damping is the Gaussian tail damping inside the correlated-motion branch (`:912-966`, `corr_if`).

### 6.6 File written (`powder_do_write`, `discus_powder_file.f90:5-27`)
- One file, no header, no comments. Each line is `WRITE(iff,*) x, y` (list-directed, default precision): two columns separated by spaces.
- Number of lines = `npkt_wrt+1` (PDF: `npkt_pdf` lines). The verified macro gives 601 lines for `r` 0–12 Å at 0.02 Å (VERIFIED_DISCUS_MACRO.md).
- File opened with `oeffne(iff,outfile,'unknown')`; overwrite behaviour UNCONFIRMED. The wrapper refuses to overwrite beforehand (`pdf.py`).
- Writer dispatch: `graphic.f90:568-580`: `pow_l_partial` → `powder_out_partial`; `cpow_form=='r'` → `powder_out_pdf`; `tof` → `powder_out_tof`; else `powder_out`.

### 6.7 Partial PDF (`powder_out_partial.f90:7-187`)
- Saves the structure (`internal.partial.temp`, `:40-41`), writes a pseudo-intensity into kuplot (`kuplot.inte_temp`, `:57-59`), then does further steps (`exp2pdf_reset`, `:60`, and beyond). **Normalisation and the pair-partial G(r) formula: UNCONFIRMED** beyond line 60.

## 7. Live checks that could be run against an independent Python calculation

Setup for L1-L9 and L13: a single-species cubic cell (one atom, `a=3.0 Å`, supercell 6×6×6 = 216 atoms, or
smaller for speed), `read`/`cell` as in the wrapper, `powder; reset; neutron; set lpcor,none; set profile,off;
set calc,debye; set axis,q; set qmin,0.5; set dq,0.01; set qmax,6.0; set stepr,0.01; set wvle,<λ>; run;
output; outf check.raw; form powd,q,0.5,6.0,0.01; value S(Q); run`. `S(Q)` for one species does not depend on `b`,
so any element will do. Python side: N atoms, all pairs, `r̃ = round(r/Δ)·Δ`, `sinc(x)=sin(x)/x`.

| ID | Check | Expected relation (derived) | Tolerance |
|---|---|---|---|
| L1 | `set axis,tth` then `show` | Still reports `Q` (`powder.f90:796-797`, `:517`) | exact |
| L2 | Label shift | `S_out[j] = S_py(qmin + (j−1)·dq)` for `j>=1`; `S_out[0] = 0`. The unshifted alignment `S_py(qmin + j·dq)` should fail by about `dq·dS/dQ` | `1e-4` for the shifted alignment |
| L3 | Debye S(Q), single species neutron | `S(Q) = 1 + (2/N)·Σ_{i<j} sin(Q r̃)/(Q r̃)`, `N = pow_nreal` (the `S` path is `phases_average:540`) | `1e-4` absolute |
| L4 | `stepr` sensitivity | Repeat L3 with `stepr` = 0.001, 0.01, 0.05. DISCUS should match the binned Python at each step (within `1e-4`) and deviate from unbinned Python by an amount that grows with `Δ` (phase error up to `Q·Δ/2`, i.e. 0.15 rad at Q=6, Δ=0.05) | `1e-4` binned |
| L5 | `I(Q) = b²·S(Q)` for neutron | `value I(Q) / value S(Q) = |b|²` (neutron, `lpcor none`: pol=1). Use `value f2aver` to confirm `⟨f²⟩=b²` at every Q | `1e-6` relative |
| L6 | X-ray polarisation | `I(Q)/(S(Q)·f2aver(Q)) = pol(2θ) = (1+0.88·cos²2θ)/2`, with `sinθ = Qλ/4π` (default `lpcor` bragg, fac 0.88). Derived from `phases_average:466` (pol in `inte` only) | `1e-6` |
| L7 | F(Q) identity | `value F(Q)` at each label equals `Q_label·(S_out−1)` from the same run | `1e-10` |
| L8 | Dimer | 2 atoms at `d=1.5 Å`: `S(Q) = 1 + sinc(Q·r̃)`, `r̃=1.50` | `1e-4` |
| L9 | `dq` rounding | `set dq,0.0123456` then `show` reports `0.01235` (`powder.f90:940`) | exact |
| L10 | `stepr` not reset | `set stepr,0.05; reset; set calc,debye; ... ` gives the same S(Q) as a fresh session with `stepr 0.05` and different from `stepr 0.01` (`reset` does not touch `pow_del_hist`) | exact |
| L11 | PSVGT broadening | Compare `reset` (no `set profile`) with `set profile,off` on the same structure. Derived expectation: differences are nonzero, of size set by FWHM from `w=0.05`. If they are exactly zero, the reading of `powder_convolute` for Debye is wrong | size: UNCONFIRMED |
| L12 | S(Q) renormalisation | With `λ=0.5` and `qmax=20` (`qmax·λ/4π` ≤ 1): mean of `S_out` over `Q>18` equals 1 (`phases_average:559-574`) | `1e-10` |
| L13 | `qmax` error | `set wvle,0.5` with `set qmax,26` gives `-108` ("Qmax is too large"); `qmax=25` runs (`powder.f90:1816`) | exact |
| L14 | `temp,use` (neutron, isotropic B) | Derived: `I(Q) = b²·[1 + e^{−2Bs²}·(S_pure(Q)−1)]`, `s=Q/4π`, `S_pure` = the no-DW Debye S from L3. Source: `fourier_sup.f90:2730`, `powder_pdf_hist.f90:623`, `powder_write_mod.f90:288-294`. The B-value used must match the structure's `cr_dw` (UNCONFIRMED: source of `cr_dw`) | `1e-4` |

## 8. Pitfalls

1. `set axis,tth` silently does Q (`powder.f90:796-797`), while the help documents `tth`.
2. `set delta,x` replaces the profile with GAUSS (`:904`). Only `set profile,off` (`:1575-1578`) or `set delta,0` with a following Gaussian profile removes the PSVGT broadening that `reset` restores (`:3609-3617`). The verified macro works only because `set delta,0.0` follows `reset`.
3. `set delta,res>0` aborts or gives wrong values on most grids (issue #65). Keep blocked.
4. `reset` in the powder menu does not restore `lambda`, radiation, `ldbw`, `ano`, `stepr`, `period` (section 4). A session-level `set temp` or `set stepr` leaks into later runs.
5. `set calc,<bad>`, `set fourier,<bad>`, `set temp,<bad>`: accepted silently with no change (`:1177-1195`, `:1196-1210`, `:1286-1298`).
6. `set dq` is rounded to 1e-5 (`:940`). `check_dq` (`:361`) doubles `dq` silently when there are more than 2^18 points.
7. Debye r-accuracy is set by `stepr` (default 0.01 Å). Phase error up to `Q·Δ/2` (0.035 rad at Q=7). Python checks must use the same binning (L4).
8. Atoms closer than `stepr/2` (0.005 Å) make the run fail with `-123` (`powder_pdf_hist.f90:533-534`).
9. Output label shift of one `dq` for S(Q), F(Q), I(Q) (6.4). The Python side must match labels, or read `Q − dq`.
10. `value pdf` (lower case) selects the phase, not the PDF (6.2). Use `value PDF`.
11. `value r` is the real part of the structure factor (`graphic.f90:810-811`), not a PDF token.
12. `form powd,q...`: any token beginning with `q` matches (`str_comp` quirk, 1.3).
13. The `form pdf` r-step is recalculated as `(rmax−rmin)/(npkt−1)` (`powder_write_mod.f90:1076`). Use `(rmax−rmin)/dr` that is an integer.
14. S(Q) is divided by its mean over `Q>18 Å⁻¹` when the grid has more than 50 such points (`phases_set_mod:559-574`). S(Q) absolute scale therefore depends on `qmax`.
15. `set scale` does not change S(Q) (`powder_write_mod.f90:477-497`).
16. `lpcor` changes `inte`, not S(Q), F(Q), PDF (derived). `set lpcor,bragg` without an angle gives `fac=0.5`, not 0.88 (`:1225-1230`).
17. Neutron runs keep Bragg polarisation unless `set lpcor,neutron` or `none` is given (`lpcor` is not changed by `neutron`).
18. `set period` survives `reset` and changes the PDF (section 3).
19. Occupancy: the Debye pair sum uses `cfact_pure` without `occ` (`powder_pdf_hist.f90:655`). Occupancy enters only through `pow_nreal` and the `⟨f⟩` weights. For partial occupancy this is UNCONFIRMED.
20. The Debye sum is over the **finite** cluster of the cell as read, not a periodic crystal (section 5.2). `cell model.cll,5,5,5` gives a finite 5×5×5 cube. Python must use the same atoms with no images.
21. Qmax saturates empirically near `12.3/λ` (VERIFIED_DISCUS_MACRO.md), while the code cap is `4π/λ − 1e-4` (`:356`) and the error is at `4π/λ` (`:1816`). The cause of the saturation is not identified (UNCONFIRMED).
22. The default X-ray table is Waasmaier–Kirfel only if `xray` is typed without `table:` (`:242-243`). `table:waass` (help spelling) silently keeps International Tables.
23. `dispersion anomalous` takes `|f|` in `⟨f⟩` (6.3): not exact for anomalous scattering.
24. `set preferred` and `set bragg` act only in complete mode.

## 9. Feature table (for the wrapper)

| Feature | Worth wrapping? | Why |
|---|---|---|
| `set qmin`, `qmax`, `dq` | Yes (done) | Required for every run; keep `dq` a multiple of 1e-5 and `qmin` a multiple of `dq`. |
| `set wvle` (numeric) | Yes (done) | Required; `energy` is an alternative but its units differ by radiation. |
| `set calc,debye` | Yes (done) | The per-atom Debye path is what the wrapper relies on. Add a check that `show` reports Debye. |
| `set profile,off` | **Yes, add explicitly** | Removes the PSVGT default (L11) without relying on `set delta,0`. |
| `set delta` (`res>0`) | No | Issue #65. |
| `set stepr` | **Yes** | Controls Debye accuracy (section 7, L4). Expose it, default 0.01, and report it with results. |
| `set lpcor` | Maybe | Matters only for `inte` and `I(Q)`. Expose for `I(Q)` work. |
| `set temp,use` | Later | Changes `I(Q)`/`S(Q)` through B; needs the L14 check first. |
| `set dispersion anomalous` | Later | `|f|` approximation in `⟨f⟩` (pitfall 23). |
| `set scale`, `set back` | No | Post-scaling only; no need for G(r). |
| `set axis,tth`, `tthmin/max`, `dtth` | No | Changes nothing for the Q-based wrapper; `axis,tth` is misleading (pitfall 1). |
| `set calc,complete` | No (for now) | Different S(Q) definition (Lorentz, `:511`); details UNCONFIRMED. |
| `set calc,nufft/grid` | No | Not read. |
| `set period` | No (for now) | Periodic correction and structure restore UNCONFIRMED. |
| `set partial` | Later | Valuable, but normalisation UNCONFIRMED (6.7). |
| `set preferred`, `set bragg`, `set sh/sk/sl`, `set ka21`, `set fourier stack`, `set profile uvw/asym/tof` | No | Not part of the 1D isotropic powder PDF. |
| `set qzero`, `set tthzero` | No | Shift the x-axis only; effect on F(Q) UNCONFIRMED. |
| `show` | Yes | Record the settings in each run for reproducibility. Confirm the `set` values took effect (L1, L9). |
| `form pdf,r,rmin,rmax,dr` | Yes (done) | The PDF grid. |
| `form powd,q,start,end,step` | Yes | Needed for S(Q), F(Q), I(Q) checks. Use an aligned grid. |
| `value S(Q)`, `F(Q)`, `I(Q)`, `f2aver`, `faver2` | Yes, for checks | Give the intermediate values needed for L2-L8 and L14. |
| `value inte` | Maybe | Only for the polarisation check (L6). |
| `value PDF` | Yes (done) | The output. |
| `value PDF` with `set period` | No (for now) | Periodic correction UNCONFIRMED. |
| `form tth`, `tof`, `value 3DPDF`, `value 3DBETA` | No | Not part of this wrapper. |
| `powder_out_partial` (`set partial` + PDF) | Later | See 6.7. |

## 10. UNCONFIRMED list

1. Case folding of commands and keywords (`get_cmd` and `do_powder` do not fold; any upstream folding not found).
2. Runtime values of `ldbw`, `ano`, `lambda`, `diff_table` after a macro start (call chain confirmed, values not run).
3. Size of the PSVGT (`w=0.05`) broadening of the Debye pattern under the default profile (L11).
4. Cause of the empirical Qmax saturation near `12.3/λ` (VERIFIED_DISCUS_MACRO.md).
5. `cr_dw` meaning: inferred to be `B` from `powder_f2aver` (`u² = ΣB/8π²`); the reader that sets `cr_dw` was not traced.
6. Occupancy weighting in the Debye pair sum (`cfact_pure` has no `occ`).
7. `phases_place` (`powder.f90:437`) and the exact copy path from convoluted `pow_conv` into `pha_powder`.
8. `powder_complete` internals: `csf`, `scale_diffuse`, `powder_split` weights, `calc_preferred`.
9. `powder_nufft` and `calc grid`.
10. `set period`: the structure restore after the sphere cut, and the exact PDF correction for finite `period`.
11. `set partial`: normalisation and the formula in `powder_out_partial` beyond its line 60.
12. Electron form factors and electron energy units (`dlink`, `get_wave` electron branch).
13. Symbolic wavelengths (`per_symwl`) values.
14. Bodies of `asym` and `chemistry` (`powder.f90:126-138`).
15. Effect of `set qzero` on written `F(Q)` values.
16. Module values `pdf_rminu`, `pdf_rmaxu`, `pdf_deltaru` used when `form pdf` gives no limits (`pdf_mod`).
17. Behaviour of `value` with lower-case or unknown tokens at run time (inferred from code, pitfalls 10-11).
18. Slope of the low-Q extrapolation in `fft_fq` uses `xwrt(1)`, which is `qmin+dq` (`discus_powder_fft.f90:65`). Effect size not measured.
19. Debye convolution step `xdel=(qmax−qmin)/num(1)` versus `pow_deltaq` (`powder_write_mod.f90:1329`). Effect not measured.
20. Overwrite behaviour of `oeffne(...,'unknown')` in `powder_do_write`.
21. The `ldbw` I(Q) closed form (L14) and the `ldbw` S(Q) branch (uses `pow_fu`, not traced).
22. `f2averb`, `faver2b`, `faverb` value codes.
23. `pow_k12` (Kα2) formula.
24. Whether the `xray table:discamb` setting persists in a later plain `xray` in the same menu session (`opara` is set on menu entry only, `:80`).
25. `powder_calc_fwhm_symm` (used in the Qmax widening, `:350`).
