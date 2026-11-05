# DISCUS 6.23.01: `pdf` and `exp2pdf` menu facts for the Python wrapper

Read-only audit of the Fortran. No source, build file or test was changed.
Every fact cites a file and line. Anything not established by reading the code
is marked **UNCONFIRMED**. "Code-read defect" means the code looks wrong from
reading; its runtime effect is UNCONFIRMED until it is run.

Citation keys (paths relative to the roots shown):

| Key | File |
|---|---|
| P | `/home/claude/dc_6.23.01/discus/prog/pdf.f90` |
| PM | `.../discus/prog/pdf_mod.f90` (module defaults) |
| AL | `.../discus/prog/discus_allocate_appl_mod.f90` |
| MS | `.../discus/prog/modify_mod.f90` (`atom_select`, `mole_select`) |
| FS | `.../discus/prog/fourier_sup.f90` (`form`, `dlink`) |
| RS | `.../discus/prog/rmc_sup_mod.f90` (`rmc_set`, `rmc_set_mode`, `rmc_genmove`) |
| RM | `.../discus/prog/rmc_mod.f90` |
| RB | `.../discus/prog/discus_blk_appl.f90` (RMC defaults) |
| PS | `.../discus/prog/discus_s_pdf_file.f90` (compiled PDF writer) |
| X | `.../discus/prog/discus_exp2pdf.f90` (exp2pdf dispatch) |
| XD | `.../discus/prog/discus_exp2pdf_data.f90` |
| XL | `.../discus/prog/discus_exp2pdf_load.f90` |
| XR | `.../discus/prog/discus_exp2pdf_run.f90` |
| XS | `.../discus/prog/discus_exp2pdf_supp.f90` |
| XW | `.../discus/prog/discus_exp2pdf_write.f90` |
| PF | `.../discus/prog/discus_powder_fft.f90` (`fft_fq`) |
| QV | `.../discus/prog/qval.f90` (`val_faver2`) |
| OF | `.../discus/prog/output_format_mod.f90` |
| UW | `/home/claude/dc_6.23.01/lib_f90/unified_write.f90` |
| SC | `/home/claude/dc_6.23.01/lib_f90/str_comp_mod.f90` |
| RN | `/home/claude/dc_6.23.01/lib_f90/random.f90` (`ran1`, `gasdev`) |
| RMOD | `/home/claude/dc_6.23.01/lib_f90/random_mod.f90` |
| RSM | `/home/claude/dc_6.23.01/lib_f90/random_state_mod.f90` |
| KD | `/home/claude/dc_6.23.01/lib_f90/kdo_all_mod.f90` (`seed`) |
| W | `/home/claude/dc_6.23.01/lib_f90/wink_mod.f90` (`zpi`) |
| H | `.../discus/prog/appl_dis.hlp` (help text) |
| CAT | `/home/claude/discus-python/src/discuspdf/data/suite_catalogue/discus_suite_6.23.01.json` (cited as JSON path `menus/<menu>/commands[i]`) |

---

## 0. Fork check and build selection

- `diff -q -r /home/claude/dc_6.23.01/discus/prog /home/claude/dc_fork/discus/prog`
  lists these differing files: `CMakeLists.txt`, `appl_dis.hlp`, `diffuse_mod.f90`,
  `discus_allocate_appl_mod.f90`, `four_strucf_OMP_mod.f90`, `fourier_finufft.f90`,
  `fourier_sup.f90`, `graphic.f90`, `rmc_sup_mod.f90`; the fork adds
  `four_blas_threads.c` and `four_strucf_gemm_mod.f90`.
- **`pdf.f90`, `pdf_mod.f90`, all six `discus_exp2pdf*.f90` files,
  `discus_s_pdf_file.f90` and `discus_powder_fft.f90` are identical in the fork.**
- Fork changes on the same route, checked by hunk:
  - `rmc_sup_mod.f90`: hunks in `rmc_calc_scat`, `rmc_inten`, `rmc_strucf` and
    `rmc_inlot` (the LOT filter called from `rmc_genmove`). The `rmc_genmove` body
    itself has no hunk.
  - `fourier_sup.f90`: hunks only in `four_run_nufft` and `four_dbwtab`.
    `form()` (FS:2947-2971) is identical in the fork.
  - `discus_allocate_appl_mod.f90`: hunk in `alloc_diffuse_dbw`, not `alloc_pdf`.
  - `diffuse_mod.f90`: one added flag `four_dbw_valid` (Debye-Waller cache).
  - Effect on the PDF and exp2pdf paths: UNCONFIRMED (only `rmc_inlot` with LOT
    enabled and the DW cache are plausibly reachable).
- `lib_f90` lives at `/home/claude/dc_6.23.01/lib_f90`, not under `discus/prog`.
- Build: `CMakeLists.txt:37-43` compiles the exp2pdf files and
  `discus_s_pdf_file.f90`. `discus_pdf_file.f90` is marked "OBSOLETE; standalone no
  longer supported" (`discus_pdf_file.f90:2`) and appears only in a commented
  `#set` (`CMakeLists.txt:79`).

---

## 1. `pdf` menu: commands and syntax as parsed

### 1.1 Parser rules

- `str_comp(a, b, j, la, lb)` (SC:6-57) compares at most `lb` characters of the
  keyword, and requires at least `j` characters typed. Keywords therefore match by
  prefix: `calculate` matches `calc`, `runaway` matches `run`.
- The menu runs `get_cmd` then a fixed `ELSEIF` ladder (P:67-286). A line with `=`
  is evaluated as maths unless it starts with `echo`, `system`, `help` or `?`
  (P:80-86).
- The menu does not inherit `common` (CAT `menus/discus/pdf/inherits` = `[]`), so
  `seed` and other common commands are not available inside `pdf` (UNCONFIRMED at
  runtime; from catalogue and dispatch).

### 1.2 Top-level dispatch (P:98-285)

| Command (min abbrev) | Line | Syntax as parsed | Effect |
|---|---|---|---|
| `continue` (3) | 98 | `continue` | macro continue |
| `echo` (2) | 103 | `echo <text>` | echo |
| `evaluate` (2) | 108 | `evaluate <expr>` | `do_eval` |
| `exit` (2) | 113 | `exit` | `GOTO 9999` (P:328) |
| `help` (2), `?` (1) | 118-126 | `help [errors]` | `do_hel` |
| `system` (2) | 130-138 | `system <cmd>` | `do_operating` |
| `wait` (3) | 142 | `wait` | `do_input` |
| `reset` (2) | 147-162 | `reset all` only | `pdf_reset` (P:3558); any other form gives error -6 (P:155-160) |
| `calc` (2) | 182-190 | `calc` (no args) | `pdf_setup(CALC)`, `pdf_skal = 1/rmc_skal(1)` (P:187), `pdf_determine(.true.)`; sets `pdf_success` |
| `data` (2) | 194-196 | `data <file>` | `pdf_mode = FIT`, `pdf_readdata` (P:889) |
| `run` (2) | 200-206 | `run` | `pdf_mode = FIT`, `pdf_setup(FIT)`, `pdf_run` (P:1738) |
| `save` (2) | 210-211 | `save {stru\|pdf\|mark}, <file>` | `pdf_save` (P:779) |
| `set` (2) | 215-216 | `set <subcommand> ...` | `pdf_set` (P:1140) |
| `show` (2) | 220-232 | `show [pdf\|all\|<token>]` | 0 args: `pdf_show('PDF')`; 1 arg: first 3 chars upper-cased; other tokens go to `rmc_show` (P:720-724) |
| `iselect`/`ideselect` (3) | 236-242 | `isele`/`idese` `{list}` | `atom_select` on `pdf_allowed_i` (P:239) |
| `jselect`/`jdeselect` (3) | 244-250 | `jsele`/`jdese` `{list}` | `atom_select` on `pdf_allowed_j` (P:247) |
| `seleect`/`deselect` (3/2) | 255-270 | `sele`/`dese` `{list}` | `atom_select` on RMC arrays `rmc_allowed` (P:267) |
| `mselect`/`mdeselect` (2) | 274-279 | `msel`/`mdes` `{list}` | `mole_select` on `rmc_allowed` (P:277) |
| unknown | 283-285 | | error -8 |

- Before most commands, `pdf` grows its arrays if `cr_nscat > PDF_MAXSCAT`
  (P:169-178).
- The catalogue lists the same 23 names under `menus/discus/pdf/commands[0..22]`,
  with `seleect` and `deselect` as written in code. The help file documents `sele`,
  `dese`, `msel` and `mdes` (H:6411-6428).

### 1.3 `set` subcommands (`pdf_set`, P:1140-1736)

Match is by the first 3 characters for most keywords (`cpara(1)(1:3)`), or by
`str_comp` with the minimum length shown. Arguments after the keyword are the
parameter count checked by the code.

| Keyword (code match) | Arguments | Sets | Lines | Default (PM / reset) |
|---|---|---|---|---|
| `BOUNDARY` (3, 8) | `periodic` or `crystal` `[,2D\|3D]` `[,exact]` | `chem_period(1..3)`, `pdf_2d`, `pdf_lexact` | 1189-1243 | `pdf_2d=F` (PM:94), `pdf_lexact=F` (PM:109) |
| `CALC` (3, 4) | `exact` | `pdf_lexact` | 1247-1253 | F |
| `FINITE` (3, 6) | `periodic` / `sphere,D` / `polygon,D` / `tanh` | `pdf_finite` 0/1/2/3; `pdf_sphere` or `pdf_diam_poly` | 1257-1284 | `PDF_BACK_PERIOD` (PM:82, PM:10) |
| `QALPHA`/`QBROAD` (3, 6) | `<a>` | `pdf_qalp` (width ∝ r) | 1296-1306 | 0 (PM:70) |
| `GAMMA`/`CORRLINEAR` (3, 5/10) | `<g>` | `pdf_clin_a` | 1310-1320 | 0 |
| `DELTA`/`CORRQUADRATIC` (3, 5/13) | `<d>` | `pdf_cquad_a` | 1324-1334 | 0 |
| `DEN*`/`RDE*` (3 chars) | `auto` or `<rho0>` | `pdf_lrho0` (T = auto); `pdf_rho0`; RDE sets `pdf_lrho0_rel` | 1338-1354 | auto (PM:108) |
| `DIAMETER` (3, 8) | `<d>` | `pdf_diam = abs(d)` | 1358-1367 | 0 |
| `FRANGE` (3, 6) | `<rmin>,<rmax>` | `pdf_rfminu/rfmaxu`; mode FIT | 1371-1384 | sentinel 0.01 (PM:61-63) |
| `QMAX` (3, 4) | `<qmax>` | `pdf_qmax` (0 disables termination) | 1388-1397 | 30 (PM:50) |
| `QSIGMA`/`QDAMP` (3, 6/5) | `<sigma_q>` | `pdf_sigmaq` | 1401-1411 | 0 (PM:58) |
| `PARTIAL` (3, 7) | `internal` or `a-list, b-list, w` | `pdf_lweights`; `pdf_weight(a,b)=w` both orders | 1415-1454 | `pdf_lweights=F` (PM:107) |
| `POLYGON` (3, 7) | `<c1> [,c2 ...]` | `pdf_poly_n`, `pdf_poly(1..n)` | 1458-1470 | `pdf_poly_n=0` (PM:83) |
| `RADIATION` (3, 9) | `neutron` / `xray [,qx]` / `electron [,qx]` | `pdf_lxray`, `pdf_radiation`, `pdf_xq` | 1474-1511 | see 2.3 |
| `RANGE` (3, 5) | `<rmax>,<dr>` | refused if data loaded (1516-1520); `pdf_rmaxu`, `pdf_deltaru`, `pdf_us_int`; mode CALC | 1515-1558 | rmax 50, dr 0.01 (PM:49, PM:54) |
| `REFINE` (3, 6) | `none`/`scale`/`density`/`lattice [a\|b\|c\|alpha\|beta\|gamma]` | `pdf_refine_*` flags (PM:112-114) | 1562-1655 | all F |
| `WEIGHT` (3, 6) | `<s>` | `pdf_scale = abs(s)` | 1656-1665 | 1 (PM:77) |
| `SHAPE` (3, 5) | `<s>` | `pdf_shape = abs(s)` | 1669-1678 | 0 |
| `SRATIO` (3, 6) | `<srat>,<rcut>` | `pdf_srat`, `pdf_rcut` | 1682-1692 | 1, 0 |
| `THERMAL` (3, 7) | `gauss` or `off` | `pdf_gauss = (arg starts 'gauss')` | 1696-1702 | F (PM:91) |
| `CYC*`,`DIS*`,`SIG*`,`MOD*`,`MDI*`,`SCA*`,`MOV*` | forwarded to `rmc_set` | see 1.4 | 1706-1722 | RB:137-183 |
| other | | error -6 | 1726-1728 | |

- `set` with fewer than 2 tokens gives error -6 (P:1730-1732).
- The help file's keyword list (H:6490-6520) matches these branches. Help also
  documents `set "bound","crystal"` with `quick`/`exact`; the code reads `exact`
  only (P:1208), and `quick` has no branch (UNCONFIRMED whether it matters).

### 1.4 RMC keywords forwarded from `pdf set`

| `pdf set` keyword | `rmc_set` branch | Meaning | Lines | Default |
|---|---|---|---|---|
| `cycles` (`CY`) | RS:138-146 | `rmc_maxcyc` | | 5000 (RB:138) |
| `display` (`DI`) | RS:205-212 | `rmc_display` (print interval) | | 500 (RB:137) |
| `sigma` (`SI`) | RS:413-420 | `rmc_sigma`, the **acceptance temperature** for CHI2 (not broadening) | | 0.01 (RB:148) |
| `mode` (`MOD`) | RS:321-333 -> `rmc_set_mode` RS:565-650 | move mix: `<prob>, <SHI\|SWD\|SWC\|EXT>` (RS:593-599) | | `rmc_mode` = shift (RB:155) |
| `mdis` (`MD`) | RS:287 | `rmc_mindist(i,j)` | | |
| `scal` (`SCA`) | RS:381-391 | `rmc_doskal` on/off (RS:383), `rmc_skal = 1/value` | | off (RB:141); skal 1.0 (AL:2993) |
| `move` (`MOV`) | RS:334-341 -> `rmc_set_move` RS:454 | `rmc_maxmove` (gasdev step sizes) | | 0 (RB:161-163) |

- The exact argument order of `set mode` is read from the code: the first numeric
  argument is the probability and the next is the mode name (RS:321-333).
  UNCONFIRMED at runtime.
- `rmc_move_cprob` is filled only inside `rmc_set_mode` (RS:639-642). Its default is
  zero (RM:50-51). See pitfall 3.

### 1.5 Start values

Module defaults are in `pdf_mod.f90` (PM:10-116). `pdf_reset` (P:3558-3660)
restores the same values, with these exceptions noted:

- `pdf_reset` assigns `pdf_poly(5) = 0` only (P:3620-3630 region; the array is not
  cleared as a whole).
- `pdf_reset` does not touch `chem_period`, `pdf_2d` is reset (P:3636), `pdf_lexact`
  is reset (P:3644).
- Key defaults: `pdf_qmax=30` (PM:50), `pdf_rmaxu=50` (PM:49), `pdf_deltaru=0.01`
  (PM:54), `pdf_deltari=0.001` (PM:51), `pdf_us_int=1` (PM:55), `pdf_skal=1` (PM:57),
  `pdf_scale=1` (PM:77), `pdf_rho0=0` with `pdf_lrho0=T` (PM:72, PM:108),
  `pdf_sigmaq=0` (PM:58), `pdf_xq=0` (PM:59), `pdf_qalp=0` (PM:70),
  `pdf_gauss=F` (PM:91), `pdf_finite=0` (PM:82), `pdf_radiation=XRAY` (PM:86),
  `pdf_lxray=F` (PM:90), `pdf_power=4` (PM:87), `pdf_mode=CALC` (PM:56),
  `pdf_success=F` (PM:89).

### 1.6 `select`, `deselect`, `show`, `save`, `reset` details

- `atom_select` (MS:1764-1886): arguments are atom types by name or number. `-1`
  selects all (MS:1849-1857 region). The optional `site:` list is parsed
  (MS:1860-1885) and stored in `pdf_lsite_i/j`. **`pdf_lsite_i/j` are never read by
  the PDF calculation** (grep: only allocation and the `atom_select` argument), so
  `site:` has no effect on G(r) (code-read).
- `mole_select` (MS:1889-1954): first token starting with `a` selects all molecule
  types; otherwise integer types are used (MS:1941-1947).
- `show`: `pdf_show` (P:568-777). `show all`/`show pdf` print the PDF settings
  (P:593-716). `show mode` and `show atom` call `rmc_show('MOD')` and
  `rmc_show('ATO')` via the general branch (P:720-724).
- `save` options are described in 3.4.
- `reset` requires the keyword `all` (P:147-162). Catalogue syntax says `reset`
  alone (CAT `menus/discus/pdf/commands[8]`).

---

## 2. `calc`: the G(r) formula

### 2.1 Pipeline (`pdf_determine`, P:2347-2626)

1. Quick (unit-cell) mode needs `cr_ncatoms /= 0` (P:2383-2392).
2. Molecule B-values and look-up tables are built (P:2402-2465).
3. Exact mode (`pdf_lexact`) or unit-cell mode, with the accumulation routine chosen
   at P:2543-2587.
4. `pdf_convtherm` (P:2609) turns integer pair counts into weighted histograms, with
   optional Gaussian broadening.
5. `pdf_convert` (P:2613) forms G(r): normalisation, baseline, finite-size factor,
   Q-damping and termination.
6. The sum of weights is printed (P:2620; format P:2624).

### 2.2 Pair enumeration and binning

| Path | Routine (lines) | Binning and cut | Notes |
|---|---|---|---|
| Unit cell, all atom types (default) | `pdf_addcorr_n_fast` (P:2927-3038) | integer counts per species pair in `pdf_temp`; `ibin = int((d+dr/2)/dr)` (P:3027); accept `d <= rmax` (P:3024) | self pair falls in bin 0 and is skipped (P:3469) |
| Unit cell, partial or Gaussian | `pdf_addcorr` (P:2761-2925) | `ibin = nint(d/dr)` (P:2854); accept `dr < d <= rmax` (P:2853); weight added per pair (P:2912-2913) | used by `pdf_run` atom moves (P:1910-1914) |
| Exact, all types | `pdf_addcorr_e_all` (P:3205-3267) | pairs `i<j` counted in both orders; `ibin <= int(rmax/dr)+1` (P:3205 region) | |
| Exact, periodic | `pdf_addcorr_ep_all` (P:3269-3355) | unit-cell cell list | |
| Exact, molecules | `pdf_addcorr_e_all_mol` (P:3356) | | |
| Exact, partial | `pdf_addcorr_e` (P:3137-3204) | | |

- Periodic images come from `pdf_bnd` (P:433-443) with cell loops of width
  `int(rmax/a)+1` (P:2799-2800). For non-periodic boundaries the loops are clamped to
  the cell (P:2809-2815).
- `set boundary periodic` sets `chem_period(1..3)` (P:1199-1201). `2D` sets `pdf_2d`
  (P:1206), which limits loops to one cell in directions with `icc <= 1` (P:2801-2804).

### 2.3 Weights

- `<f>` and weights (P:469-485):
  - `pdf_natoms(is) = NINT(count_is * occ_is)` (P:471).
  - `bave = sum_is f_is * pdf_natoms(is) / N_real` (P:473-475), `N_real = cr_n_real_atoms` (P:472).
  - `w_ij = f_i * f_j / bave**2 * occ_i * occ_j` (P:480-482).
  - User partials replace `w` (P:1444-1445).
- `f` is `form(is, cr_scat, pdf_lxray, h2, pdf_power)` with `h2 = pdf_xq**2` (P:450,
  P:473, P:480). `form` (FS:2947-2971):
  - X-ray (`lxray=T`): `f = scat(1) + sum_{i=1..power} scat(2i) * exp(-scat(2i+1) * h2)`.
  - Neutron (`lxray=F`): `f = scat(1)`.
  - `power = pdf_power = 4` (PM:87).
- The coefficient table `cr_scat` is filled by `dlink(pdf_radiation, ...)` (P:447-448):
  - Neutron: `cr_scat(:,i)=0` then `cr_scat(1,i) = temp_bcoh` (FS:3045-3046).
  - X-ray: `cr_scat(:,i) = temp_scat` from `get_scat_xray` (FS:3069-3074). That
    routine fills `temp_scat(1) = per_fc` (the constant c), then `per_fa1, per_fb1, ...,
    per_fa4, per_fb4` at indices 2-9 (`lib_f90/element_data_mod.f90:3244-3262`).
    So index 1 of the X-ray table is the constant term, not a coherent b value.
    The comment at FS:3013 ("a1,b1,---a4,b4,c") lists c last and is inconsistent with
    the fill order.
- **Default weights (code-read defect; see pitfall 1).** `pdf_radiation` defaults to
  X-ray (PM:86), so `dlink` loads the X-ray table, but `pdf_lxray` defaults to
  `.false.` (PM:90; reset P:3634). Only `set radiation xray|electron` sets
  `pdf_lxray=T` (P:1483, P:1494). Without that call, `form` takes the neutron branch
  and returns `scat(1)` = `per_fc` alone (FS:2947-2971). The X-ray term sum
  `sum scat(2i)*exp(...)` is dropped, so the weight is not f(0) = c + sum a_i.
  `show` reports "X-Rays" by default (P:602-612), so the displayed and used radiation
  disagree. Runtime size of the error UNCONFIRMED.

### 2.4 Normalisation and baseline (`pdf_convert`, P:2628-2759)

- Number density: default `pdf_lrho0=T` gives `r0 = N_real/V/ncc` (P:2664-2665), with
  `ncc = icc1*icc2*icc3` (P:2655). `set density <v>` gives `r0 = v` (P:2661).
  `set density rde <v>` gives `r0 = rho0*N_real/V/ncc` (P:2659). The line before it
  (P:2658) is overwritten and is dead code.
- Scale: `norm = pdf_scale / N_real` (P:2667). It is also divided by `dr` unless
  Gaussian broadening or `qalp` is active (P:2668-2670).
- Baseline: `rr = 2*zpi*r*r0*dnorm = 4*pi*r*rho0*dnorm` (P:2681), with `zpi = 2*pi`
  (W:11). The help definition is `G(r) = 4 pi r (rho(r) - rho(0))` (H:6395).
- G at bin i (P:2717): `G_i = norm * corr_i / r - rr`. For periodic sphere
  (P:2706-2715) the same form is multiplied by the sphere factor.
- Finite-size options (P:2680-2703):
  - `PDF_BACK_POLY` (P:2682-2691): `rr = 4*pi*r*rho0*dnorm - sum_k p_k r^k` for
    `r < diam_poly`, clamped at zero; zero beyond.
  - `PDF_BACK_SPHERE` (P:2692-2699): `rr *= (1 - 1.5 x + 0.5 x^3)` for `r < D`, where
    `x = r/D` and `D = pdf_sphere` (a diameter). The same factor multiplies G in the
    periodic branch (P:2708-2711).
  - `PDF_BACK_TANH` (P:2700-2703): `rr = max(0, -4*pi*r*r0*dnorm*tanh(shape*(r-diam)))`.
- `dnorm` for partials: `pdf_dnorm = (sum of selected w) / (sum of all w)` (P:494-511).

### 2.5 Q-resolution damping (`set qsig`, `set qdamp`)

- Applied only when `pdf_sigmaq > 0` (P:2723-2733): `pdf_calc(i) *= pdf_exp(NINT((i+1)*fac4))`
  with `fac4 = dr*sigma_q/step` and `step = 0.0005` (PM:93). By the code this is
  `exp(-(sigma_q*r)^2/2)` with `r = (i+1)*dr`, an off-by-one bin against the G index
  (UNCONFIRMED against a run).
- **Code-read defect (pitfall 2):** `pdf_exp` is initialised only inside `pdf_setup`,
  and only if `pdf_gauss` or `pdf_qalp > 0` (P:552-562). `pdf_reset` zeroes it
  (P:3587). With `set qsig` alone, `pdf_exp` stays zero, so G is damped to zero.

### 2.6 Termination (Qmax) and ripple

- Default `pdf_qmax = 30` (PM:50), so termination is on unless `set qmax 0`.
- Kernel (P:535-545): `K(0) = Qmax` and `K(+-k) = sin(Qmax*k*dr)/(k*dr)` for `k` up to
  `2048` (PM:28 gives `PDF_MAXSINCC = 2**12+1`; P:541-544 fills the kernel).
- Convolution (P:2749-2755): `CONVLV_SUB` (utilities.f90:4657-4696, NR `realft`
  based), then `G_term(r_i) = (dr/pi) * sum_k G_k * K(i-k)`. `factor = dr/zpi*2`
  (P:2751) equals `dr/pi`.
- `CONVLV_SUB` needs a power-of-two `n` and an odd kernel length (utilities.f90:4670-4671).
  `n` comes from `alloc_pdf` as `2**(INT(LOG(n_dat)/LOG(2))+2)` (AL:2381-2384).
- No Lorch or other window was found in `pdf_setup` or `pdf_convert`.
- The kernel is truncated at `+-2048` steps: about 2.05 A at `dr = 0.001`, and wider
  for larger `dr`. Effect UNCONFIRMED.

### 2.7 Thermal Gaussian and `set delta`/`gamma`/`qalpha`/`srat`

- Enabled by `set therm gauss` (P:1698) or `qalp > 0` (P:3476).
- Per-pair width (P:3479-3485):
  `sigma^2 = fac*(cr_dw_i + cr_dw_j + B_mole) - cquad/r^2 - cqua_mole/r^2 - clin/r - clin_mole/r`,
  clamped to `>= 0`, with `fac = 1/(2*zpi^2) = 1/(8 pi^2)` (P:3459).
  `cr_dw` is the structure's per-type Debye-Waller term (`class_crystal.f90`
  `cr_dw(0:nscat)`). **UNCONFIRMED** whether `cr_dw` is B or U^2; the arithmetic
  reads it as B.
- Then `sigma^2 += (qalp*r)^2`, `sigma = sqrt(sigma)`, `sigma *= srat` if `r <= rcut`
  (P:3484-3489). Window `igaus = 1 + nint(5*sigma/dr)` (P:3491).
- Kernel (P:3499-3519): `gnorm = 1/(sqrt(2 pi)*sigma)`, with asymmetry factor
  `(1 + ig*dr/r)` (P:3518).
- `set delta`/`set gamma` change `sigma^2` only when Gaussian is on (P:3480-3481).
- `set sigma` is **not** a broadening parameter (RS:413-420; see 1.4).

### 2.8 Partial PDFs and i/j selections

- `iselect`/`jselect` set `pdf_allowed_i/j` (P:236-250). A pair `(i,j)` counts if
  `(allowed_i(is) and allowed_j(js)) or (allowed_j(is) and allowed_i(js))`
  (P:2840-2841, P:3465-3466).
- `all_atoms` skips the partial path when every type is allowed (P:2531-2534).
- Partial weights from `set partial` (P:1420-1449) override `pdf_weight`.

---

## 3. `data`, `run`, `save`

### 3.1 `data` file format (`pdf_readdata`, P:889-1047)

- One filename argument only; any other form gives error -6 (P:1028-1031).
- Optional header:
  - A first line starting `History` opens a block read by `extract_hist`
    (P:1049-1087). Keys `temp=` and `Qmax=` are read until a line beginning
    `##### start data`. `pdf_qmax` takes `res_para(2)` (P:1060-1063).
    UNCONFIRMED which key lands in `res_para(2)`.
  - Lines starting with `#` are skipped (`skip_spec`, P:1088-1105).
- Data records: four free-format columns `r  G  dr(ignored)  sigma` (P:959-993).
  - `pdf_wic = 1/sigma^2`, or 1 if `sigma = 0` (P:968-972, P:989-993).
  - **Conflict:** help says the last value "gives the weight w(G) ... e.g.
    1/(error of G)" (H:6440-6452). Code squares it, so **column 4 must be the
    error**. Code wins; see pitfall 9.
- Grid: `deltaru = (re - ra)/(N-1)` (P:1006); `us_int = nint(deltaru/deltari)`
  (P:1008); `dr = deltaru/us_int` (P:1009); `bin = nint((re-ra)/dr)` (P:1014);
  `rmax = re` (P:1015); fit bounds `rfminf = ra`, `rfmaxf = re` (P:1017-1018).
- **Requirement:** `|ra - deltaru| <= 1e-6`, otherwise error -5 (P:1020-1024).
  The first `r` must equal the grid step (help H:6446-6447 says the same).
- Success sets `pdf_ldata = T` (P:1027), which blocks `set range` (P:1516).
- `data` sets `pdf_mode = FIT` (P:195).

### 3.2 `run` (`pdf_run`, P:1738-2050)

Setup and guards:

- `pdf_setup(FIT)` (P:202, P:334-397). Fit window: `rfmin = max(rfminu, rfminf)`,
  `rfmax = max(rfmaxu, rfmaxf)` (P:393-394). Errors -9 and -10 for out-of-file user
  limits (P:381-392).
- `pdf_run` guard `pdf_obs(1) == -9999` gives error -4 (P:1781). The default of
  `pdf_obs` is 0 (AL:2409-2411 region), so this guard does not catch a missing data
  file (UNCONFIRMED impact; see pitfall 4).

Fit window and chi-squared (P:1801-1848):

- `nmi = nint(rfmin/dr)`, `nma = nint(rfmax/dr)` (P:1801-1802).
- Sums over `ip = nmi .. nma/us_int`, comparing `pdf_obs(ip)` with
  `pdf_calc((ip-1)*us_int + 1)` (P:1833-1841), with weight `pdf_wic(ip)`.
- Scale `s = pdf_skal`:
  - `rmc_doskal` off (default): `s = 1/rmc_skal(1)` (P:1845).
  - `rmc_doskal` on: least squares `s = sum(w*calc*obs)/sum(w*calc^2)` (P:1843).
- `chi = (ee + s^2*cc - 2*s*ce) / wtot` (P:1847-1848), with
  `ee = sum w*obs^2`, `cc = sum w*calc^2`, `ce = sum w*calc*obs`, `wtot = sum w`.
- Initial state from `pdf_determine(.false.)` (P:1823).

Main loop (P:1867-1999):

- `sig2 = rmc_sigma^2 / 2` (P:1865). Default `rmc_sigma = 0.01` (RB:148).
- Step schedule per `igen` (P:1876-1881):
  - scale/density step when `MOD(igen,10) == 2` (the early-phase clause is redundant);
  - lattice step when `MOD(igen,200) == 5`, or `MOD(igen,10) == 5` when
    `igen < 0.05*rmc_maxcyc`;
  - otherwise an atom move (P:1886-1890).
- Atom move (`rmc_genmove` RS:3263-3416; `rmc_genmove_mol` RS:3420 when not
  `rmc_sel_atom`):
  - move type from `RANDOM_NUMBER` against cumulative `rmc_move_cprob` (RS:3300-3303);
  - SHIFT: random atom, `p_new = p + gasdev(rmc_maxmove)` (RS:3329 onwards);
  - SWCHEM and SWDISP: swaps or displacements of two atoms;
  - EXTERNAL: `rmc_genmove_ext` (`extrmc.f90:112`).
  - ROTATE (`rmc_mode = 5`) has no branch in `rmc_genmove` (see pitfall 3).
- Accept test (P:1910-1955):
  1. Subtract old contributions (`pdf_addcorr(..., -2)`, P:1910), apply the move
     (`pdf_makemove`, P:2314-2345), add new contributions (`+2`, P:1914), then
     `pdf_convert` (P:1916).
  2. `cnew` from the same chi expression; `prob = cnew - cold` (P:1940).
  3. If `prob < 0`: accept (P:1942-1943).
  4. Else if `sig2 > 0`: `p = exp(-prob/sig2)`, accept when `p > ran1` (P:1950-1951).
     `ran1` is `RANDOM_NUMBER` (RN:108-122).
  5. Else reject (P:1953).
  6. `rmc_sigma == -9999` accepts all moves (P:1957).
- Reject: positions restored (P:1971), `pdf_corr` restored from `pdf_old`
  (P:1972-1974), `pdf_convert` (P:1975).
- Counters: `itry` counts atom attempts (P:1918); `iacc_good` for downhill and
  `iacc_bad` for accepted uphill moves (P:1961-1967).
- Scale/density step, `pdf_rmc_scale` (P:2052-2152):
  - `pdf_scale += gasdev(ref_maxpdfsd(1))` if `refine scale` (P:2086-2088);
    `pdf_rho0 += gasdev(ref_maxpdfsd(2))` if `refine density` (P:2090-2092).
  - Same Metropolis test (P:2122-2141; `ran1` at P:2131). Reject reverts scale and
    rho0 (P:2147-2148).
- Lattice step, `pdf_rmc_lattice` (P:2154-2313):
  - `cr_a0` and `cr_win` changed by `gasdev` for the flagged parameters
    (P:2180-2250 region); `setup_lattice` rebuilds the cell.
  - Full recalculation: `pdf_corr = 0` and `pdf_addcorr` over all atoms
    (P:2243-2247), then `pdf_convert`.
  - Same Metropolis test; reject restores the lattice and recalculates.
- Loop control: `loop = itry < rmc_maxcyc` (P:1987). Error -20 if
  `igen > 1000*rmc_display` with `itry = 0` (P:1989-1993).
- Output:
  - every `rmc_display` generations (P:1995-1998), format at P:1300:
    `Gen: <igen> try: <itry> acc: (good/bad): <g> / <b> s2x2: <cold>`;
  - "Final configuration (s= ...)" at P:1250;
  - "Delta Chi : ave: ... sig: ... max: ..." at P:3000 (printed at P:2006);
  - timing at P:4000.
- After the run: `rmc_skal(1) = 1/pdf_skal` (P:2015), so the fitted scale persists
  and is picked up by the next `calc` (P:187). `update_cr_dim` (P:2032).
- Results are stored in `res_para(1..8)`: final chi, tries, good, bad, average
  delta, sigma, max delta, time per cycle (P:2019-2028).

### 3.3 Things the fit does not do

- No compensation of `wic` in the chi, beyond the `1/sigma^2` weights.
- `pdf_addcorr` in atom moves uses per-pair Gaussians (P:2858-2906), while
  `pdf_determine` uses the integer-count route plus `pdf_convtherm`
  (P:3027, P:3476-3530). Both must agree; they are two code paths (pitfall 18).

### 3.4 `save` (`pdf_save`, P:779-888)

- `save stru, <file>`: `save_struc` (P:821-825). Format not checked here
  (UNCONFIRMED).
- `save pdf, <file>` (P:829-846):
  - requires `pdf_success` (P:830), else error -12 (P:834-836);
  - calls `pdf_save_file` (P:837-838) (PS:9-98).
  - Rows: `i = nmi..nma step us_int`, `r = i*dr`, `G = pdf_skal*calc(i)` (PS:81-88).
  - Format `F9.4,3X,F21.10,5X,2(F6.2,1X)` (PS:94), with two constant columns
    `0.0` and `1.0`.
  - A file name beginning `kuplot` writes into KUPLOT memory instead (PS:47-60
    region).
  - Calc-mode window: `rfmin = deltaru`, `rfmax = rmaxu` (P:374-377).
- `save mark, <file>` (P:847-866): error -7 if `pdf_gauss` (P:849-851). Otherwise
  writes `r, 0.00` (format P:886, `F9.4,3X,F6.2`) for bins with `calc > 0`.
- CAT lists a `format` optional for `save` (`menus/discus/pdf/commands[12]`), but
  `pdf_save` parses no format keyword (UNCONFIRMED: ignored).

### 3.5 Probe results: fit and partials (2026-11-05, `tools/probes/pdf_fit_probes.py`, tasks 019-021)

- **`data` after `reset all` (or in a fresh session) fails** with "Too many points in PDF" (-1): `reset` sets
  PDF_MAXDAT = 1 (P:3568), `data` allocates exactly n lines (P:944-951) and the read loop tests `ip > PDF_MAXDAT`
  before the read that meets the end of the file (P:978-981). A `calc` over a larger range first sizes the arrays.
- A `save pdf` file of `calc` starts at r = dr and is valid `data` input as it stands.
- The fit log has the RMC format (`Gen: … try: … acc: (good/bad): … s2x2:`, "Final configuration", "Delta Chi").
  Swaps of random Cu0.75Au0.25 (108 atoms, 2000 tries, σ 0) against the L1_2 G(r): χ² 0.576 → 0.206, never rising,
  α1 0.036 → −0.318, Au count kept.
- `save pdf` after `run` and `calc` + `save pdf` on the fitted crystal differ by 8e-4 of max |G| for r ≤ 7.5
  (pitfall 18: two accumulation paths). Near rmax every Gaussian lacks its partners beyond rmax, so a model on a
  range ending at rmax falls off within ~0.5 Å of it (0.2–0.6 of max |G| at 7.8–8.0 Å for B = 0.5).
- **Partials** (neutron, periodic, no Gaussian; binary Cu3Au and ternary Cu2AuAg; residuals ~1e-8):
  - (a, a): G = w_aa P_aa/(N dr r) − c_a² w_aa 4πrρ0 (P over both orders);
  - `isele a; jsele b`, a ≠ b: G = w_ab P_ab/(2 N dr r) − 2 c_a c_b w_ab 4πrρ0, half the pairs its baseline assumes,
    so G_total ≠ G_aa + G_bb + G_ab (the difference is exactly the missing half);
  - {a, b} on both sides counts the cross pairs in full: G_{a,b} − G_aa − G_bb is the correct cross partial;
  - no selection state leaks between `calc`s in one session.

---

## 4. Random numbers

- `ran1(idum)` ignores `idum` and returns `RANDOM_NUMBER` (RN:108-122; comment at
  RN:110 says it replaces an old Numerical Recipes generator).
- `gasdev` is a polar method on `RANDOM_NUMBER` (RN:11-49).
- `idum = -182783467` is declared but unused by `ran1` (RMOD:12).
- Seeding: the `seed` command is in the common menu (KD:179-180). With `0` it calls
  `RANDOM_SEED()` and seeds from the time (RSM:85-99).
- The `pdf` menu does not inherit `common` (CAT `menus/discus/pdf/inherits` = `[]`),
  so the seed must be set in the parent menu before `pdf` (UNCONFIRMED at runtime).
- The generator is the compiler's `RANDOM_NUMBER`. Which compiler and version is used
  for the build is UNCONFIRMED, so runs are reproducible only for a fixed build and
  seed.

---

## 5. Fork (repeat of section 0)

- Source files of the PDF and exp2pdf route are identical in the fork.
- The fork differs in `rmc_inlot` (used by RMC moves when the LOT filter is on) and in
  `four_run_nufft`/`four_dbwtab`. `form()` is unchanged.
- Consequence for the wrapper: none in the PDF and exp2pdf sources. The RMC LOT filter
  and Debye-Waller cache differences can affect fits or powder inputs (UNCONFIRMED).

---

## 6. `exp2pdf` menu

### 6.1 Dispatch (X:116-185)

| Command | Line | Handler | Syntax as parsed | Defaults and effect |
|---|---|---|---|---|
| `back` (2) | 116-117 | `exp2pdf_load(2)` (XL:9) | `back xy,<file> [,scale:<s>]`; `back csv,<file>`; `back "kuplot",<n>` | background; `scale:` default 1.0 (XL:52); sets `exp_bscale` (XL:138) |
| `composition` (2) | 121-122 | `exp2pdf_composition` (XS:9-232) | `composition [comp:current \| comp:<string>]` or bare `<string>` | `exp_comp_current` default T (XD:29); bare string sets F; string gives element names and occupancies |
| `data` (2) | 126-127 | `exp2pdf_load(1)` | `data xy,<file>`; `data csv,<file>`; `data "kuplot",<n>` | data set; x, y and sigma = abs(dy) (XL:237-246) |
| `exit` (2) | 131 | | | `lend = T` |
| `faver` (5) | 136-137 | `exp2pdf_load(4)` | `faver xy,<file>`; `faver csv,<file>`; `faver "kuplot",<n>` | `<f>` vs Q; squared at load (XL:299-310); file form bug (pitfall 8) |
| `outfile`/`outputfile` (2/10) | 141-143 | `exp2pdf_outfile` (XS:343-441) | see 6.2 | |
| `poly` (2) | 147-148 | `exp2pdf_poly` (XS:651-704) | `poly order:<n>` or `poly <n>` | default in data module 7 (XD:28, XR:1737 region); `poly` default opara 9.000 (XS:687-689) |
| `radiation` (2) | 152-153 | `exp2pdf_radiation` (XS:236-338) | `radiation "xray"\|"neutron"\|"electron" [lambda:<A> \| energy:<v>]` | lambda 0.10 A (XS:279, XD:63); energy 123.9842004 (XS:279, XD:64) |
| `reset` (2) | 157-158 | `exp2pdf_reset` (XR:1702-1763) | `reset` | see pitfall 21 |
| `run` (2) | 162-164 | `exp2pdf_run` (XR:11-90), then `exp2pdf_init` | `run [mode:silent\|inter]` | default `silent` (XR:37-38 region, XR:53-63) |
| `show` (2) | 168-169 | `exp2pdf_show` (XS:780-811) | `show` | prints settings |
| `qscale` (2) | 173-174 | `exp2pdf_qscale` (XS:708-776) | `qscale qcrystal:<v> [,qobs:<v>]` or `dcrystal:<v> [,qobs:<v>]` | see pitfall 6 |
| `limits` (2) | 178-179 | `exp2pdf_qmax` (XS:445-647) | `limits inst:[<qmin>,<qmax>] [,fourier:[<qmin>,<qmax>]]` | see pitfall 9 |
| other | 184 | `kdo_all` | | generic commands such as `set` (KD:~179-190) |

- Commands from `discus_kdo_common` are tried first (X:108-110). `exp2pdf` inherits
  `discus/common` and `common` (CAT `menus/discus/exp2pdf/inherits`), so `seed` is
  available here.
- `faver` requires the full word (min abbreviation 5, CAT `commands[4]`).
- Catalogue `run` lists options `nthread`, `useomp`, `plot`, `scale`, `damp`,
  `width`, `filter`, `range`, `dchi`, `pshift`, `conf`, `chi`, `status`
  (CAT `menus/discus/exp2pdf/commands[10]`). The code reads only `mode`
  (XR:31-40). The Lanczos parameters are hard-coded (XR:639). See pitfall 23.

### 6.2 `outfile` (XS:343-441)

- Keywords: `rmin`, `rmax`, `rstep`, `scale`, `gr`, `iq`, `sq`, `fq`, `format`
  (XS:390). Defaults: names `discus.grobs`, `discus.iqobs`, `discus.sqobs`,
  `discus.fqobs` (XS:392-396); `rmin 0.01`, `rmax 100.01`, `rstep 0.01`, `scale 1`
  (XS:392-396); `format stan` as the printed default (XS:395).
- `rmin/rmax/rstep/scale/format` take effect only when `gr` is present (XS:405-418).
- `format:nexus|common` gives 16 and `format:stan` gives 0 (XS:408-415). The data
  module default is 16 (XD:74, XD:75-77), so **the default GR output is the HDF5
  "common" format**. See pitfall 5.
- `iq`, `sq`, `fq` set the file name and enable the output (XS:421-441).
- Initial `exp_outgr_l = T`, `exp_outiq_l = F`, `exp_outfq_l = F`, `exp_outsq_l = F`
  (XD:70-73; reset XR:1738-1762 region).

### 6.3 Input files

- `data`, `back`, `faver` go through KUPLOT's loader (`do_load`, XL:130). Formats
  `xy`, `csv` and `kuplot` are handled there. Their details are UNCONFIRMED here,
  because the loader is outside this scope.
- After load, `exp_x = x`, `exp_data = y`, `exp_sigma = abs(dy)` (XL:237-246).
- Background arrays: `exp_back = y`, `exp_sigmab = abs(dy)`, `exp_xb = x`
  (XL:272-279).
- `faver`: `exp_faver2 = y**2` (XL:310). The file therefore holds `<f>` (average
  form factor), not `<f>^2`. Help does not say this (H:2583-2586). UNCONFIRMED from
  help, confirmed from code.
- Data must be on a Q scale (H:2570-2572 summary text).

### 6.4 `run`: algorithm, in order (XR:11-90 and onwards)

1. Mode parse (XR:59-67). `silent` sets `prompt, off, off, save` (XR:63).
2. Spline onto an equidistant Q grid with step `exp_qstep = 0.001` (XD:60,
   XR:91-224; natural cubic spline via `spline_mod`, XR:225-272):
   - Grid nodes: `q1 = (nint(xmin/dq)+1)*dq`, `q2 = (nint(xmax/dq)-1)*dq` (XR:225-233 region).
   - Effective range: `qmin = max(first data Q, (int(inst_qmin/dq)+1)*dq)` and
     `qmax = min(last data Q, (int(inst_qmax/dq)-1)*dq)` (XR:150-200 region).
   - With background: `I_pure = I - bscale*B(Q)` (XR:166); sigma sum linear:
     `sigma = sigma_data + bscale*sigma_bg` (XR:167).
   - `exp_qmax_f = min(exp_qmax_f, exp_qmax)` (XR:177, XR:204).
3. `<f>^2` (XR:314-441):
   - If no `faver` file (`exp_kfaver == 0`, XR:314), and composition is not
     `current`, a synthetic structure is built: cubic cell a = 1.000 A, angles 90 deg,
     one atom per element at `x = i`, occupancy `atocc/sum` (XR:314-355 region).
     With `comp:current` the loaded crystal is used.
   - Radiation: x-ray lambda 0.1 A (XR:~340-355 region); electron energy forced to
     200 (XR:375); neutron lambda 0.1 (XR:~380 region).
   - `exp_is_rad` selects the user lambda or energy (XR:~385-395).
   - `powder_run` with output `val_faver2` (XR:429, XR:434). `val_faver2 = 10`
     (QV:12) computes `<f>^2 = (sum_l f_l(Q) * amount_l / N)^2` (QV:216-232).
   - Read back from KUPLOT dataset `kuplot.exp2pdf.faver2` (XR:433-441).
   - After `powder_run`, `exp_renergy = renergy` (XR:430-431), which overwrites an
     electron energy the user set (code-read; pitfall 22).
   - With a `faver` file, `exp_faver2` is taken from the file (XL:310).
4. Normalise and multiply by Q (XR:487-489): `y = I_pure/<f>^2 * Q`,
   sigma likewise. Sigmas below `TOL` become 9999, or all 1 if every sigma is below
   `TOL` (XR:452-486 region).
5. Polynomial background (XR:501-555):
   - KUPLOT `poly, N` with `N = exp_npoly` (XR:501), unit weights (`one`, XR:515),
     80 cycles (XR:517), status off (XR:519).
   - Fit by `do_fit` (XR:542); residual kept as `F_init - poly` (XR:555 comment:
     "Get F(Q)-poly").
   - Help formula: `F = F_init - sum_{i=1..N} p_i Q^i` (H:2530-2537), which has no
     constant term. Whether the KUPLOT `poly` fit includes one is UNCONFIRMED.
6. Lanczos smoothing for Q limits (XR:639-641): `filter:lanczos, width:151, damp:0.005`
   on a copy (dataset `ik`). Only used to find Qmax and Qmin (below).
7. Fourier Qmax (XR:645-691): search for a zero crossing of the smoothed F(Q) near
   `exp_qmax_f`, then take the nearer sample (XR:683-691). Rule UNCONFIRMED in detail.
8. Fourier Qmin (XR:693-762): if `fourier:` lower limit given (XR:741-744), use it;
   otherwise choose among minima of a 51-point smoothed curve (XR:725 `, 51`) below
   the first-peak Q (XR:746-762). Rule UNCONFIRMED in detail.
9. Q-scale correction (XR:741-800): only if `qscale ... qobs:` was given
   (`exp_qfirst_o > 0`, XR:765). Find the nearest maximum to `qobs`
   (or to `qcrystal`/`2pi/dcrystal` if `qobs < 0.01`), then
   `scalex = qcrystal / Q_peak` (XR:791). All Q values and `exp_qmax_f` are scaled by
   `scalex` (XR:~792-800).
10. Empirical scale (XR:826-827):
    `scalef = (0.40 / (sum(|y|)*dq / sqrt(max(y)) / Qmax))^2`, then `y *= scalef`.
    **The G(r) absolute scale is data-dependent (help H:2538-2541 says so).** `scalef`
    is applied before the FQ and SQ outputs, so those files are scaled too.
11. Transform `fft_fq` (XR:843; PF:7-112):
    - Linear ramp from `Q = 0` to the first point: `temp(i) = i*dq*y1/x1` (PF:~40).
    - Discrete sine transform (`ddst`) on `2**18` points.
    - R grid `r_i = (i+0.5)*pi/((N-1)*dq)`, weights `(2/pi)*dq` (PF:~65-72).
    - Spline onto the output grid `rmin..rmax` with step `rstep` (PF:~75-110).
    - Result: `G(r) = (2/pi) * sum_k dq * F_k * sin(Q_k r)` (code). The help formula
      omits `2/pi` (H:2534-2536). No Lorch window.
    - Output points: `npdf = int((rmax-rmin)/rstep)+1` (XR:837); default 10001.
12. Scale output: `yfour *= exp_scale` (XR:878). `exp_scale` is the `outfile scale:`
    value (XS:418).
13. Outputs (XR:879-955):
    - GR (XR:879-900): `ncols = 2`, labels `r[Angstroem]` and `G(r)[Angstroem^-2]`
      (XR:870-873 region).
      - format 16: `exp2pdf_write_unified` (XW:10-199) writes HDF5 via `h5fortran`
        (UW:27-33, UW:110-114); `output_format_common = 16` (OF:24).
      - format 0 (`stan`): ASCII with `# HEADER`/`# DATA`/`# EXPER`/`# SAMPLE`/
        `# CORREC`/`# OUTPUT` lines (XR:1599-1700 header writer; key lines at
        XR:1666-1697), then `powder_do_write` (XR:900).
      - Header keys useful for parsing: `# OUTPUT : n_points`, `# OUTPUT : min,max,step`,
        `# OUTPUT : value ;2C::` (XR:1694-1697 region).
    - FQ (XR:909-914): `exp_outfq_l` gives `Q[Angstroem^-1]`, `F(Q)[Angstroem^-1]`,
      written from the scaled `y_wrt` (XR:913-914).
    - SQ (XR:916-932): `y_wrt = y_wrt/Q + 1` (XR:919-921), i.e. `S(Q) = 1 + F(Q)/Q`
      using the scaled F. Written at XR:932.
    - IQ (XR:938-953): `exp_x` and `exp_data` after spline and background subtraction,
      not divided by `<f>^2` (XR:943-944, XR:953).
14. Cleanup (XR:72-83; `exp2pdf_init` at X:164): `exp2pdf_init` clears the IQ, SQ,
    FQ flags and Q-scale flags (XS:815-837) after every `run`, so `outfile iq:` etc.
    must be repeated.

### 6.5 Probe results on the stock engine (2026-11-05, `tools/probes/exp2pdf_probes.py`)

Three rounds (mailbox tasks 015-017). They settle what the code reading left open:

- **`run` needs a set-up the help does not mention.** It finishes only after a crystal is read, a `powder`
  calculation runs on it in the same session, and `composition comp:current` is used.
  - `composition <string>` builds an internal crystal (XR:315-355) whose form factors are not found
    ("Unknown element, no Fourier calculated", -20).
  - Without a preceding powder run every `run` reports "Fractions add up to 0.000" (-169): `pha_frac(1)` is set
    to 1 only inside a powder run that reaches its calculation (powder.f90:376-377), and exp2pdf does not check
    the error code after its own `powder_run` (XR:429-437), so an earlier error surfaces as -169.
  - `comp:current` without a preceding powder run aborts DISCUS (SIGABRT). A powder run on one crystal followed
    by `read` of another, then exp2pdf, also aborts.
  - The composition crystal may give the composition by whole atoms (3 Cu + 1 Au) or by occupancies (Cu 0.75,
    Au 0.25): identical G(r).
- **Composition strings**: `Cu`, `Cu3 Au1` parse; `Cu:3 Au:1` is an error (-12, read as an optional parameter
  name), `Cu3Au` is an error (-1). (Only relevant to the broken string route.)
- **`limits`**: `inst:` alone fails ("Lower fourier limit is less than lower instrumental limit"); both given
  works. The Fourier Qmax is moved to a nearby zero of F(Q): 16.5 asked, 16.473 used (header `lim_four_app`).
- **Format**: without `format:` the GR file is HDF5 (`\x89HDF`); `format:stan` gives the text header of 6.4 step 13.
- **Transform**: G(r) equals (2/π) Σ dq F sin(Qr) over the written F(Q) plus the linear ramp, to 2.8e-4 of max |G|.
- **Electron energy** 200 vs 300 keV, **sigma column** present or not: identical G(r). **Background** with
  `scale:2` on data + 2B: equal to the background-free run to 1e-11.

The wrapper (`discuspdf.discus_exp2pdf`) writes this set-up for every run.

---

## 7. Proposed deterministic live checks

Keep structures small. Use cubic cells with a few hundred atoms at most. Test all
relations with a Python reimplementation that reads only the DISCUS output files.

**T1. `calc` vs brute-force periodic pair histogram (unit-cell path).**
- Structure S1: simple cubic, a = 2.5 A, one atom, `icc = 4 4 4` (64 atoms).
  Structure S2: rock-salt 1:1, a = 5.64 A, `icc = 3 3 3` (216 atoms).
- Commands: `set boundary periodic`, `set qmax 0`, `set therm off` (default), `set
  finite periodic` (default), `set density auto` (default), `set radiation neutron`
  (so `f = b`; for S1 the weight is 1), `set range 10, 0.01`, `calc`, `save pdf, t1.dat`.
- Python: for every ordered pair of atoms and every lattice image with
  `0.005 < d <= 10`, bin with `k = round(d/0.01)` and count `H_k` weighted by
  `w_ij = b_i b_j / <b>^2 * occ_i occ_j` (`<b>^2` as in P:473-482).
  Expected: `G(r_k) = H_k/(N * 0.01 * r_k) - 4*pi*r_k*rho0`, `rho0 = N/V_total`.
  Here `N` = atoms in the supercell (`cr_n_real_atoms`) and `V_total` = supercell volume.
  Tolerance: 1e-8 relative (integer counts; only floating-point rounding).
  - UNCONFIRMED: whether `cr_v` is the unit-cell volume (then `rho0` as stated).
  - Expect the `k = 0` bin to be absent in the output (P:3469).
- Repeat with `set calc exact` and `set boundary crystal`, against a non-periodic
  Python cluster (no images, ordered pairs, `k >= 1`).

**T2. Termination kernel.** With `set qmax 30` (default) and the T1 structure, compare
`calc` to the Python convolution of the T1 `G` with `K(k) = sin(Qmax*k*dr)/(k*dr)`,
`K(0) = Qmax`, times `dr/pi`, for `|k| <= 2048`. Tolerance set by the kernel truncation.

**T3. Q-damping: exact relation, and a test for pitfall 2.**
- `set qsig 0.5`, then `calc`. Expected (by the code): `G_d(r_i) = G(r_i) * exp(-(0.5*r_{i+1})^2/2)`,
  with the one-bin shift noted in 2.5 (UNCONFIRMED).
- Run once with only `set qsig` and once with `set therm gauss` active. The first
  run is expected to show zeros (pitfall 2). That result is itself the check.

**T4. Thermal width.** With `set therm gauss` and `cr_dw` set in the structure file
(units UNCONFIRMED), the peak sigma should follow `sqrt((B_i+B_j)/(8 pi^2))`. Check the
broadening trend across two B values rather than an absolute value.

**T5. `fft_fq` and output transform (exp2pdf).** From a synthetic I(Q) (DISCUS powder
of S1 or S2, or any smooth function), run `exp2pdf` with `format:stan` and `fq` output.
Read the FQ file. Expected: `G_out(r_k) = c * (2/pi) * dq * sum_k F_k sin(Q_k r_k)`,
with one constant `c = scalef * exp_scale` for all `r_k`. Test the constancy of
`G_out / Python_transform` across `r`, not its value. Tolerance 1e-3 relative after
spline interpolation.

**T6. `<f>^2` from `val_faver2`.** For the composition-only path, Python
`(sum_l c_l f_l(Q) / N)^2` with the same form-factor coefficients should match the
`kuplot.exp2pdf.faver2` values (QV:216-232). Requires access to KUPLOT memory or
an equivalent output. UNCONFIRMED availability.

**T7. Q-scale correction.** Scale a synthetic I(Q) by a known factor `f`. With
`qscale qcrystal:<q0> qobs:<q0/f>`, the output peak positions should move by
`scalex = q0/(q0/f)` (XR:791). Test the peak-position ratio.

**T8. Peak positions only.** Compare `exp2pdf` G(r) peak positions with the `pdf calc`
G(r) of the same structure. Expect agreement to within `rstep` (0.01 A). Absolute G
values are not expected to match (empirical scale, pitfall 7).

---

## 8. Pitfalls and feature table

### 8.1 Pitfalls (ordered by severity)

1. **Default `calc` uses the neutron branch on the X-ray table.** `pdf_radiation`
   defaults to X-ray (PM:86) and `dlink` loads the X-ray table (P:447-448), but
   `pdf_lxray` defaults to `.false.` (PM:90; P:3634), so `form()` returns `scat(1)` (= `per_fc`, element_data_mod.f90:3244)
   (FS:2947-2971). Always set `set radiation xray` or `neutron` explicitly before
   `calc` or `run`. The code-level error is confirmed; runtime size UNCONFIRMED. `show` says X-Rays (P:602-612).
2. **`set qsig` alone zeroes G(r).** `pdf_exp` is initialised only under Gaussian or
   `qalp` (P:552-562), but the damping reads it (P:2723-2731). UNCONFIRMED at runtime.
3. **`run` without `set mode` may hit an unimplemented move.** `rmc_move_cprob` is
   zero by default (RM:50-51) and set only in `rmc_set_mode` (RS:639-642).
   `rmc_genmove` then selects mode 5 (RS:3300-3303), which has no branch there
   (RS:3329-3360 covers modes 1-4). `natoms` is then undefined.
   UNCONFIRMED at runtime.
4. **`run` without data.** `pdf_wic` defaults to 0 (AL:~2408-2412) and the guard
   `pdf_obs(1) == -9999` (P:1781) never fires by default. `wtot = 0` then divides
   (P:1848). UNCONFIRMED.
5. **Default exp2pdf GR file is HDF5.** `format` defaults to 16 (XD:74), giving
   `exp2pdf_write_unified` (XR:879-900, XW:10). A text reader needs
   `outfile gr:<f>, format:stan` (format applies only with `gr`, XS:405-418).
6. **`qscale ... dcrystal:` ignores the `dcrystal` value.** Lines XS:763, XS:771-772
   read `owerte(O_QVAL)` inside the `O_DVAL` branch. A user `dcrystal` gives a zero
   divisor. Code-read defect.
7. **exp2pdf G(r) and S(Q) are empirically scaled.** `scalef` (XR:826-827) and
   `exp_scale` (XR:878) are data-dependent. FQ and SQ outputs carry `scalef` (XR:913-921).
   The help says so (H:2538-2541). Do not treat absolute values as physical.
8. **`faver` file branch never records its name or index.** `elseif(itype==3)` is
   repeated (XL:122, XL:126), so `exp_cfaver`/`exp_kfaver` are not set for `itype==4`
   (XL:106-109). The squared `<f>` is still set (XL:299-310). Effect on the fit
   UNCONFIRMED.
9. **`limits` default bounds fail validation.** The default `inst` range `[0, 1e9]`
   (XS:483-485) is checked unconditionally (XS:539-557, XS:599-637). The code path
   suggests error -61 when `inst:` is omitted, and it compares against
   `exp_x(exp_dim(1))`, so `limits` before `data` reads unallocated arrays.
   UNCONFIRMED at runtime.
10. **`data` column 4 is an error, not a weight.** Code uses `1/sigma^2` (P:968-993),
    while help says "weight ... e.g. 1/(error of G)" (H:6440-6452). Supply sigma.
11. **`data` first r must equal Δr** (P:1020-1024; H:6446-6447).
12. **Keyword prefix matching is loose** (SC:6-57). `calculate` acts as `calc`.
13. **Full-word `select` is not accepted.** The branch compares `'seleect'` for 6
    characters (P:255-257), so the word `select` does not match. Use `sele`
    (H:6418-6420). UNCONFIRMED at runtime (code-read -8).
14. **Full-word `mselect` deselects.** `lselect` uses `str_comp(befehl,'mselelect',2,lbef,9)`
    (P:279). For `mselect` the comparison is FALSE, so molecules are deselected. `msel`
    (4 chars) gives TRUE. UNCONFIRMED at runtime.
15. **Uninitialised `ldummy` into `atom_select`.** `ldummy` (P:56) is passed as
    `sel_atom` (P:239-242, P:247-250), which is read at MS:1832. Whether the previous
    i/j selection is cleared depends on undefined memory. UNCONFIRMED.
16. **`site:` has no effect** on the PDF (P:236-250, MS:1860-1885; no reads of
    `pdf_lsite_i/j`).
17. **`set polygon` with more than 5 coefficients** writes past `pdf_poly(5)`
    (P:1462-1466; PM:78). No bound check. Code-read defect.
18. **Two accumulation paths.** Per-pair Gaussians in `pdf_addcorr` (P:2858-2906) and
    the integer route plus `pdf_convtherm` (P:3027, P:3476-3530) should agree. Test
    both with T1 and a Gaussian case.
19. **Occupancy enters twice.** `NINT(count*occ)` (P:471) and `occ_i*occ_j` in
    `w` (P:482). Correct only for full occupancy UNCONFIRMED.
20. **Sinc kernel truncated at 2048 steps** (PM:28; P:541-544). Effect UNCONFIRMED.
21. **`exp2pdf_reset` keeps wavelength, energy and `exp_is_rad`** (XR:1702-1763;
    XD:61-64 defaults not restored). Set `radiation` again after `reset`.
22. **Electron energy overwritten.** `renergy = 200` (XR:375) is copied back to
    `exp_renergy` (XR:430-431). A user-set electron energy is lost. Code-read.
23. **Catalogue `run` options are not implemented.** Only `mode` is read (XR:31-40).
    The Lanczos width and damping are fixed at 151 and 0.005 (XR:639).
24. **Help and code disagree on the exp2pdf formula.** Help omits `2/pi` (H:2534-2536)
    and states a polynomial without a constant term (H:2530-2537). UNCONFIRMED whether
    the fitted polynomial has a constant term.
25. **`save mark` fails with Gaussian broadening** (P:849-851, error -7).
26. **RNG seeding is outside the `pdf` menu** (CAT `inherits` = `[]`; KD:179-180).
27. **`reset` needs the word `all`** (P:147-162), though CAT says `reset`.
28. **`exp2pdf` uses KUPLOT global state.** It writes and reads datasets (XR:433,
    `kuplot.exp2pdf.faver2`) and clears them (XR:836, XR:834-837). Do not mix with an
    open KUPLOT session.
29. **Exact-mode normalisation of `rho0`** depends on `cr_v` meaning (UNCONFIRMED).

### 8.2 Feature table

| Feature | Worth wrapping? | Why |
|---|---|---|
| `pdf calc` G(r), unit-cell and exact | Yes | Deterministic. T1 checks it exactly against Python. |
| `pdf set range`, `qmax`, `density`, `finite sphere`, `radiation`, `boundary` | Yes | Simple scalars. Set radiation explicitly (pitfall 1). |
| `pdf set qsig`, `qalpha`, `therm`, `delta`, `gamma` | Later | Pitfall 2 (qsig), and Gaussian path has two code routes (18). |
| `pdf save pdf` (ASCII G(r)) | Yes | Simple fixed format (PS:94). |
| `pdf save stru` | No (for now) | Structure output belongs to the structure layer. |
| `pdf data` + `run` (RMC fit) | Later | Stochastic. Needs `set mode` (3), data (4), column convention (10), seed (4). |
| `pdf set partial`, `iselect`/`jselect` | Later | Pitfalls 13-16. Test after fixes. |
| `pdf show` | Yes | Read-only. |
| `exp2pdf data`, `back`, `composition`, `radiation` | Yes | Inputs with clear defaults. |
| `exp2pdf run` to G(r) | Yes, with `format:stan` | Empirical scale. Wrapped as shape only (`discus_exp2pdf`), with the set-up of 6.5. |
| `exp2pdf outfile fq/sq/iq` | Yes, with caveats | Scaled (pitfall 7); reset flags after each run (6.4 step 14). |
| `exp2pdf limits`, `qscale`, `poly` | Later | Defects 6 and 9; fix before use. |
| `exp2pdf faver` | Later | File-branch bug (8). |
| KUPLOT-based loading (`"kuplot",n`) | No | Couples to global KUPLOT state (28). |
| 3D-PDF menu (`3dpdf`) | Out of scope | Separate menu. |

---

## 9. UNCONFIRMED index (summary)

Each item is in the body with its citation. This index lists them together.

- Default-weight numerics on the X-ray table with `pdf_lxray = F` (pitfall 1; P:90, P:447-448, FS:2947-2971).
- Runtime size of the X-ray/neutron mismatch in `form()` (confirmed in code at FS:2947-2971 and element_data_mod.f90:3244; magnitude not run).
- `pdf_exp` zeroing with `set qsig` only (pitfall 2; P:552-562, P:2723-2731).
- `set mode` requirement before `run` (pitfall 3; RS:639-642, RS:3300-3303).
- Missing-data `run` behaviour (pitfall 4; P:1781, P:1848).
- Units of `cr_dw` (B or U^2) (P:3479; `class_crystal.f90`).
- `cr_v` meaning for `rho0` in exact mode (T1).
- Off-by-one of the Q-damping index (P:2731).
- Sinc kernel truncation effect (pitfall 20).
- Occupancy double counting (pitfall 19).
- Which key feeds `res_para(2)` in a History block (P:1060-1063).
- `set bound crystal, quick` (help; code has no `quick` branch, P:1208).
- `set mode` argument order and syntax (RS:321-333).
- `save` `format` keyword (CAT only; pitfall; P:779-888).
- KUPLOT `xy`/`csv` loader specifics (XL:130).
- `faver` file branch effect on fit (pitfall 8).
- `limits` default validation (pitfall 9).
- Fourier Qmax and Qmin selection rules (XR:645-762).
- Whether the exp2pdf `poly` fit includes a constant term (XR:501-555; H:2530-2537).
- Help vs code `2/pi` factor (H:2534-2536; PF:~40-72).
- Electron energy overwrite (pitfall 22).
- `poly` without order argument (XS:651-704; `ber_params` on no arguments).
- Effect of fork differences in `rmc_inlot` and DW cache (section 0).
- Compiler and version behind `RANDOM_NUMBER` (section 4).
- `seed` availability inside `pdf` (CAT inherits `[]`).
- Whether `pdf_obs`, `pdf_wic` defaults are zero at runtime (AL).
- `rmc_genmove` behaviour for `rmc_mode = 5` at runtime (RS:3300-3360).
- Whether `exp_sigma` from `dy` is zero for two-column files (XL:246).
- Behaviour of `iselect` ldummy (pitfall 15).
- Exact-mode `cr_v` usage in `rho0` (T1).
