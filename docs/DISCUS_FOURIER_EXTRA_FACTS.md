# DISCUS 6.23.01: Fourier, Output and 3DPDF Menu Facts

> **Checked by the overseer against the source (2026-11-03) and used by `FourierConfig`:**
> `electron` (fourier.f90:387-391) and the 5-Gaussian table `get_scat_elec` (lib_f90/element_data_mod.f90:3359);
> `energy` units from `get_wave` (element_data_mod.f90:3005-3075): X-ray λ = 12.398/E[keV], neutron
> λ = 9.045/√E[meV], electron relativistic with E in keV; `dispersion a|o` (fourier.f90:351-367);
> `delf NAME,f′,f″` / `delf NAME,internal` (fourier.f90:295-349, first atom only, `all` refused with the
> 3-value form); internal f′, f″ only for the Kα1 symbols (`get_scat_ano`, :3309-3355, any other gives 0);
> `value 3DPDF` replaces dsi by dsi3d and frees dsi3d (graphic.f90:520-545); **bugs**: `output_limit`
> (graphic.f90:2861-2889) zeroes `dsi(i,j,k)` with the component index k inside the loop over h's
> components (only layers 1-3), and `qval` evaluates ⟨f⟩-type values at `istl(i,j,1)` (qval.f90:188-253),
> which `type:sharp` uses (graphic.f90:530-540). Both unchanged on the fork and upstream master.
> Everything below this box was written by an agent and is only partly verified.


Scope: single-crystal scattering menus of DISCUS 6.23.01 that discus-python does
not wrap yet. Read-only audit of the Fortran; no code was changed.

Conventions used below:

- `prog/` = `/home/claude/dc_6.23.01/discus/prog/`
- `lib/` = `/home/claude/dc_6.23.01/lib_f90/`
- Every fact carries a `file:line` citation. Line numbers were checked against
  the source in this pass.
- **UNCONFIRMED** marks anything that is not verified in code, or whose runtime
  behaviour was not checked by running DISCUS. Those items are collected in the
  final section.
- "Wrapped" means the item is already written by
  `src/discuspdf/fourier_macro.py` (FourierConfig).

## 0. Shared facts

### 0.1 Menu dispatch and parsing

- The Fourier menu dispatches on abbreviated keywords with `str_comp(befehl, name, min_len, lbef, full_len)`.
  Examples: `electron` `prog/fourier.f90:387`, `energy` `:373`, `scat` `:751`,
  `set` `:819`, `show` `:824`.
- Parameters are split by `get_params` and numbers by `ber_params`
  (`prog/fourier.f90:353-355` is the pattern; the helper is in `lib/take_param_mod.f90:10-152`).
- Keyword prefix matching is defined in `lib/str_comp_mod.f90`. A keyword is
  accepted when its typed prefix is at least `min_len` characters long.
- The menu command list used by the suite is in
  `src/discuspdf/data/suite_catalogue/discus_suite_6.23.01.json` under
  `/menus/discus/fourier` (command list, `min_abbrev`, `compare_len`).

### 0.2 Wavelength, energy and the default table

- `get_wave` in `lib/element_data_mod.f90:3005-3075` converts the energy or
  wavelength. The physical constants are stored as mantissas only (see the
  PER_RAD_POWER definition near `lib/element_data_mod.f90:37`).
- Conversions as implemented (numerically checked in Python against the physical
  formulas, agreement to about 1e-7 relative):
  - X-ray: lambda (Angstrom) = 12.398 / E (keV)
  - Neutron: lambda (Angstrom) = 9.045 / sqrt(E in meV) (mantissa form 9.0448)
  - Electron: relativistic formula. At 200 keV this gives 0.025079 Angstrom.
- The energy/wavelength symbol table `per_symwl` / `per_wavel` is at `lib/element_data_mod.f90:2591-2610`.
- Default X-ray table before any `xray` command: International Tables. This follows
  from the `diff_table` initialisation (see `prog/diffuse_mod.f90:97-146`
  defaults and `prog/fourier_reset.f90`). Switching to Waasmaier happens through
  `xray` (see section 6). Runtime confirmation: UNCONFIRMED.

### 0.3 Scattering convention

- `s = sin(theta)/lambda`. Grid constants: `CFINC = 0.001`, `CFPKT = 9999`
  (`prog/discus_config_mod.f90:27-28`).
- Form factor is a Gaussian sum: `f(s) = c + sum_i a_i * exp(-b_i * s^2)`, with
  the Debye-Waller term `exp(-B * s^2)`. Implemented in `dlink`
  (`prog/fourier_sup.f90:2976-3145`); `form` at `prog/fourier_sup.f90:2947-2971`;
  `four_formtab` at `prog/fourier_sup.f90:2691-2743`.
- The help text for `scat` describes the sign as `exp(+b/4 h^2)`. The code uses
  `exp(-b s^2)`. The code is taken as authoritative (see pitfalls, section 10).

### 0.4 Radiation flags

- `neut`, `xray`, `electron` change `diff_radiation` (`RAD_NEUT`, `RAD_XRAY`,
  `RAD_ELEC`). The menu names are listed at `prog/fourier.f90:1355`.
- None of `neut`, `xray`, `electron` resets `l_energy`. An energy set earlier
  therefore persists across radiation switches (confirmed in the branches listed
  in section 1 and section 6).

## 1. `electron`

- Branch: `prog/fourier.f90:387-391`.
- Effect (code):
  - `lxray = .true.`
  - `diff_radiation = RAD_ELEC`
  - `diff_table = RAD_INTER`
  - `lambda = ' '`
- Does not set `l_energy`. `electron` without a prior `energy` therefore uses the
  wavelength that `get_wave` derives from the current settings. For the default
  electron wavelength 0.70926 Angstrom, the derived energy is about 0.2989 keV
  (derived from the relativistic conversion above, not run in DISCUS;
  UNCONFIRMED at runtime).
- Electron scattering-factor table: the help (`prog/appl_dis.hlp:850-860`, release
  note) cites "Table 4.3.2.3 Int. table Vol C (2006) respectively" and
  Peng, Ren, Dudarev, Whelan, Acta Cryst A52 (1996) 257. Which table is
  actually loaded by `get_scat_elec` (`lib/element_data_mod.f90:3359`) is
  UNCONFIRMED.
- `diff_table = RAD_INTER` (International Tables) is set at this branch; `xray`
  can change it later (section 6).
- Status: not wrapped by `FourierConfig`. The wrapper writes `xray` or `neut` only.

## 2. `energy`

- Branch: `prog/fourier.f90:373-383`.
- Takes exactly one parameter (`ianz == 1`), read with `ber_params` into
  `renergy`. Sets `lambda = ' '` and `l_energy = .true.`. Any other argument count
  gives error -6 (`ier_num = -6`, `ier_typ = ER_COMM`).
- Unit of the value: keV, implied by the X-ray conversion in section 0.2
  (`get_wave`). The help text unit was not checked; UNCONFIRMED against the help.
- Relation to `wvlength`: `wvlength` is at `prog/fourier.f90:922-941`. Setting
  one clears the other (the `lambda = ' '` on the energy branch, the `l_energy`
  flag on the wavelength branch). The exact clearing logic in `wvlength` was not
  re-read here; UNCONFIRMED.
- Status: not wrapped. The wrapper writes `wvle` (wavelength) before `xray` in
  some paths (see section 9 and `docs/FOURIER_CONTROLS.md`).

## 3. `dispersion` and `delf`

### 3.1 `dispersion`

- Branch: `prog/fourier.f90:353-369`.
- Accepts exactly one argument. The first character decides the mode:
  - starts with `a`: `ano = .true.` (anomalous dispersion on)
  - starts with `o`: `ano = .false.` (off)
  - any other value, or a bare `dispersion` with no argument, gives error -6.
- The flag `ano` is consumed by the anomalous table, `get_scat_ano`
  (`lib/element_data_mod.f90:3309`).
- Pitfall: the help text says a bare `disp` is allowed. The code gives error -6
  (see section 10).
- The wrapper writes `disp off` (`src/discuspdf/fourier_macro.py`). Per
  `docs/FOURIER_CONTROLS.md`, this makes the wavelength metadata only. Runtime
  effect: UNCONFIRMED.

### 3.2 `delf`

- Branch: `prog/fourier.f90:298-349`.
- The branch exists and sets the anomalous correction for an element. The full
  argument syntax and meaning of the stored values were not traced in this pass.
  Semantics: UNCONFIRMED.
- Status: not wrapped.

## 4. `scat`

- Branch: `prog/fourier.f90:751-815`.
- Element is resolved through `get_iscat` (`prog/get_iscat_mod.f90`). The lookup
  is an exact match after `do_cap`. If the name is not found, the error is
  `ier_num = -27` (`prog/get_iscat_mod.f90`, not-found branch). Absent-name
  behaviour in the Fourier path: UNCONFIRMED at runtime.

### 4.1 Form with 10 or 12 parameters (user-supplied Gaussian curve)

- Condition: `ianz == 10` or `ianz == 12` (`prog/fourier.f90:752-753`).
- Target index: `jj(k) = nint(werte(1))` for the element index.
- Coefficients written:
  - `cr_scat(1, jj) = werte(ianz)`, the constant term `c`
  - `cr_scat(i, jj) = werte(i)` for `i = 2 .. ianz-1` (the `a_i` and `b_i` values)
- Sets `cr_scat_int(jj) = .false.`, so the element no longer uses the internal
  table.
- The `ianz == 12` branch reuses the same loop with a larger count. The mapping
  of values to Gaussian terms for the 12-value form (number of Gaussian terms,
  and whether `cr_scat` layout is `(index, element)`) is UNCONFIRMED.
- Code observation: the inner loop is `DO i = 2, ianz-1` nested in `DO k = 1, i`
  at `prog/fourier.f90:761-765` and reuses `i`. With the current count of one
  element (`i = 1`) the outer loop body runs once, so the shadowing has no effect
  today. Recorded for the record; not a bug at present.
- Negative element index (`werte(1) <= 0`) gives error -27 (`prog/fourier.f90:776-779`).

### 4.2 Form with 2 parameters (internal or equivalence)

- Condition: `ianz == 2` (`prog/fourier.f90:785-803`).
- `scat <elem> internal`: sets `cr_scat_int(k) = .true.` and
  `cr_scat_equ(k) = .false.`. A value of `-1` applies to all elements
  (`prog/fourier.f90:791-797`).
- `scat <elem> <other_elem>`: sets `cr_scat_equ(k) = .true.` and
  `cr_at_equ(k)` to the other element name (`prog/fourier.f90:798-803`).
  This makes one element use another element's curve. Whether the equivalence
  lookup is performed at run time is UNCONFIRMED.
- Any other argument count gives error -6 (`prog/fourier.f90:809-812`).

### 4.3 Reset and show

- `prog/fourier_reset.f90` resets `cr_scat_int` to TRUE, so the internal table
  is the default after a reset (confirmed in code; the reset does not clear
  `diff_table` or `diff_exti`).
- `show scat` (`prog/fourier.f90:828`) calls `do_show_scat`
  (`prog/discus_show.f90:814-950`). It prints `a(i)`, `b(i)`, `c`, `f'`, `f''`.
  For `diff_power` cases 0, 4 and 5 the printed values come from `res_para`
  (a's, b's and c).

### 4.4 Status

- Not wrapped. The wrapper writes `scat neutron b` only (per
  `src/discuspdf/fourier_macro.py`). The 10/12-parameter form is the one useful
  for user-defined form factors.

## 5. `calculate`, `hkl`, `zone`, `layer`, `lots`, `ll`/`lr`/`ul`/`tl`

### 5.1 `calculate` (single reciprocal point)

- Branch: `prog/fourier.f90:243-294`.
- Evaluates the structure factor at one reciprocal point. Two helpers are
  involved:
  - `calc_000` (`prog/fourier_sup.f90:3149-3226`): handles the 000 reflection,
    uses `res_para` 1-4 and normalises by `cr_icc`.
  - `calc_hkl` (`prog/fourier_sup.f90:3230-3549`): writes an output file and
    handles the `.fcf` list. Extinction is applied at `prog/fourier_sup.f90:3520`.
    `calc_hkl` sets the globals `eck`, `vi` and `inc`.
- Which helper is called for a given `calculate` point: UNCONFIRMED (dispatch
  not traced line by line).
- Status: not wrapped.

### 5.2 `hkl` (SHELXL list of reciprocal points)

- Branch: `prog/fourier.f90:395-400`. The handler body is at
  `prog/fourier.f90:1253` (`process 'hkl' command`).
- Calculates at each point listed in a SHELXL-style list. Extinction correction
  from `set exti` is applied in `calc_hkl` (see 5.1) and in the SHELX / list5 /
  list9 writers (see section 7). It is not applied to the `four_run` amplitudes.
- Input file format: UNCONFIRMED (not traced).
- Status: not wrapped.

### 5.3 `zone`

- Branch: `prog/fourier.f90:976-1009`.
- `zone OFF`: `l_zone = .false.`, `inc(3) = 1`, `ltop_c = .true.`
  (`prog/fourier.f90:977-981`).
- `zone h k l res` (four parameters, `ianz == 4`): sets the following
  (`prog/fourier.f90:983-1005`):
  - `zone_uvw(1..3)` = h, k, l
  - `zone_res` = res
  - `lxray = .true.` and `diff_radiation = RAD_ELEC` (electron diffraction is
    switched on by the zone command itself)
  - `fave = 1.0` (100 % average)
  - `ilots = LOT_BOX`, `lperiod = .true.`
  - `ls_xyz(i) = MIN(10, cr_icc(i))`
  - `inc = (512, 512, 1)`
  - `renergy = 200.00` keV, `l_energy = .true.`, `lambda = ' '`
  - `l_zone = .true.`
- Other argument counts: error -6 (`prog/fourier.f90:1006-1008`).
- Implementation in `prog/zone.f90`:
  - `zone_setup` (`prog/zone.f90:14-179`) uses `rlambda` and `zone_ewald`.
  - `zone_project` (`prog/zone.f90:181-391`) sets `fave = 0`, `ilots = LOT_OFF`
    and `nlots = 1` (`:275-277`), then calls `four_run` (`:304`).
- Whether the zone pattern depends on a random lot choice: UNCONFIRMED.
- Pitfall: `zone` changes radiation to electron and energy to 200 keV, silently.
  Any earlier `energy` or `xray` setting is overwritten.
- Status: not wrapped.

### 5.4 `layer`

- Branch: `prog/fourier.f90:402-434`.
- Defines a single layer of the reciprocal grid. Argument syntax and effect on
  `inc` and `ls_xyz`: UNCONFIRMED (branch located, body not traced in this pass).
- Status: not wrapped.

### 5.5 `ll`, `lr`, `ul`, `tl`, `lots`

- Corner settings `ll`/`lr`/`ul`/`tl` at `prog/fourier.f90:436-463`. The wrapper
  already writes these (`src/discuspdf/fourier_macro.py`).
- `lots` at `prog/fourier.f90:467-528`: controls the random-lot (LOT_*) sampling.
  Not wrapped. Option semantics beyond LOT_BOX / LOT_OFF: UNCONFIRMED.

### 5.6 `run`

- Branch: `prog/fourier.f90:639-747`.
- Requires `natoms` to be set (error otherwise). Calls `four_run`
  (`prog/fourier_sup.f90:11-323`). Runs `zone_setup` / `zone_project` when
  `zone` is active.
- `four_run` sequence (code): `csf-acsf`, `dsi`, `four_weight`, `diff_scale`
  (`prog/fourier_sup.f90:281-282`), filter, `four_conv`, `four_qinfo`.
- Status: wrapped (`run`), see `src/discuspdf/fourier_macro.py`.

## 6. Other `set` options and related menu items

### 6.1 `set` dispatch

- `set` at `prog/fourier.f90:819-820` calls `fourier_set`, whose body is at
  `prog/fourier.f90:1066-1247`.
- `top` at `prog/fourier.f90:889-902` sets the ordinate. Wrapped (abs/ord/top).
- `temp` at `prog/fourier.f90:872-885` switches temperature coefficients on or
  off. Wrapped.

### 6.2 Options already used by the wrapper

- `set technique turbo`: wrapped.
- `set aver`: wrapped. `four_aver` is at `prog/fourier_sup.f90:1194-1357`. It
  returns immediately when `ave == 0`. Site selection uses `ran1 .le. ave`.
- `set symm`: wrapped as optional.

### 6.3 Options not wrapped

- `set technique` values other than `turbo`: the accepted list was not traced in
  this pass. UNCONFIRMED.
- Filter settings: `four_run` applies a filter (`prog/fourier_sup.f90:11-323`);
  defaults are set in `prog/fourier_reset.f90` (`four_*` filter defaults).
  Syntax of the `set` keywords for filter and conv: UNCONFIRMED (not traced).
- Convolution: `four_conv` in `prog/fourier_sup.f90:11-323` path. Syntax:
  UNCONFIRMED.

### 6.4 `xray` and the table choice

- Branch: `prog/fourier.f90:949-972`. Accepts a table selector: `inter`,
  `waas`, or `disc`.
  - `inter`: International Tables (the default, see 0.2)
  - `waas`: Waasmaier-Kirfel (`get_scat_xray_waas`, `lib/element_data_mod.f90:3269`)
  - `disc`: not supported by the form-factor path. `dlink`
    (`prog/fourier_sup.f90:2976-3145`) has no `RAD_DISC` branch. Runtime
    behaviour (error or silent fallback): UNCONFIRMED.
- Neutron scattering table: `get_scat_neut` (`lib/element_data_mod.f90:3294`).
- Anomalous table: `get_scat_ano` (`lib/element_data_mod.f90:3309`).
- Electron table: `get_scat_elec` (`lib/element_data_mod.f90:3359`).

### 6.5 `show`

- `show` at `prog/fourier.f90:824-868`. The `scat` subcommand is at `:828`.
- Prints current settings. Other subcommands: UNCONFIRMED (not traced).

## 7. Output menu (`graphic.f90`)

The output menu is handled in `prog/graphic.f90`. The writer core is
`do_output` (`prog/graphic.f90:1518-2128`) and `write_unified_data`
(`prog/graphic.f90:2435-2812`).

### 7.1 `format`

- Branch: `prog/graphic.f90:218-440`. Keywords (with prefix checks) and their
  `ityp` codes (`prog/output_format_mod.f90`, codes stan 0 ... common 16):
  - `3d` (`:227`)
  - `gnup` (`:232`)
  - `pgm ` (`:237`): writer disabled in 6.23.01, `do_pgm` at `prog/graphic.f90:1279`
    returns error -6.
  - `post` (`:242`): writer disabled, `do_post` at `prog/graphic.f90:1179`
    returns error -6.
  - `powd` (`:247`) with sub-options `tth` (`:250`), `q` (`:252`), `tof` (`:254`).
    The `stl`, `dst`, `lop` options are commented out (`:256-260`).
  - `pdf` (`:283`) with sub-option `r` (`:286`).
  - `ppm ` (`:306`): writer disabled, `do_ppm` at `prog/graphic.f90:1350`
    returns error -6.
  - `stan` (`:311`): 1D, 2D and 3D; `PART_` names at `prog/graphic.f90:1836`.
  - `shelx` (`:316`); `hklf4` (`:323`) maps to the same SHELX writer.
  - `list5` (`:333`), `list9` (`:343`).
  - `nexus` (`:348`): commented out. `hdf5`, `common` and `nexus` are handled in
    the branch at `:363-389`; `common` maps to OUTPUT_FORMAT_COMMON at `:389`.
    The `scale` option is read at `:366`.
  - `vtk` (`:353`), `mrc` (`:358`), `xplor` (`:396`), `vesta` (`:401`).
- The 3d NIPL writer is at `prog/graphic.f90` (`do_output`, 3d branch).
  Behaviour when `out_inc(3) > 1`: the file unit is closed inside the loop
  (UNCONFIRMED, see section 10).

### 7.2 `input`

- Branch: `prog/graphic.f90:452-484`.
- Reads an existing grid file: two integers (`out_inc(1)`, `out_inc(2)`), then
  four reals (`xmin`, `xmax`, `ymin`, `ymax`), then one real `zmax`. The data
  are then read row by row into `dsi(ix, iy, 1)`.
- The STAN files written by this same program start with an NIPL header. Whether
  `input` can read them back is UNCONFIRMED; the header layout suggests it is
  not directly compatible.

### 7.3 `outfile`

- Branch: `prog/graphic.f90:485-504`.
- Uses `get_optional` with `O_MODE` options `new`, `old`, `add`
  (`out_mode = NEW | OLD | ADD`). The default output name is `fcalc.dat`
  (`prog/output_mod.f90`).
- The output file name is taken as typed. No extension is appended
  (`lib/support_mod.f90:656-747`, `oeffne`).

### 7.4 `range`

- Branch: `prog/graphic.f90:505-509`. Calls `output_range`
  (`prog/graphic.f90:948-1047`).
- Only the XPLOR writer uses the range (`output_range`). Other writers ignore it.

### 7.5 `run`

- Branch: `prog/graphic.f90:515-670`.
- Requires that a Fourier was computed (`four_was_run`). Otherwise the error is
  -118 (`prog/graphic.f90`, `run` branch).
- Calls `chem_elem`, `set_output`, then the writer selected by `ityp`.
- For `value == val_3DPDF` the writer uses `dsi3d`. If `dsi3d` is not allocated,
  error -180 is returned (`prog/graphic.f90:521-527`). `dsi3d` is cleared after
  use.
- `set_output` is at `prog/graphic.f90:2131-2183`.

### 7.6 `show`

- Branch: `prog/graphic.f90:670-695`. Prints the current output settings.

### 7.7 `threshold`

- Branch: `prog/graphic.f90:706-782`. Sets `ps_high` and `ps_low` (percent of
  `diffumax`) and derives `zmax` / `zmin`.
- Inert for the writers in 6.23.01 because pgm, ppm and post are disabled
  (see 7.1). Not useful to wrap.

### 7.8 `value`

- Branch: `prog/graphic.f90:782-874`. Keywords are case-sensitive.
- Leading `<` selects the average: `laver = .true.`, `ix = 2`
  (`prog/graphic.f90:782-794`).
- Value codes (`prog/qval.f90:1-20`):

| code | constant | keyword family |
|---|---|---|
| 1 | val_inten | intensity |
| 2 | val_ampli | amplitude |
| 3 | val_phase | phase |
| 4 | val_real | real part |
| 5 | val_imag | imaginary part |
| 6 | val_ranph | random phase |
| 7 | val_sq | S(Q) |
| 8 | val_fq | F(Q) |
| 9 | val_f2aver | <F^2> |
| 10 | val_faver2 | <F>^2 |
| 11 | val_faver | <F> |
| 12 | val_norm | commented out (`prog/qval.f90`, line 12) |
| 13 | val_iq | I(Q) |
| 14 | val_pdf | PDF |
| 15 | val_3DPDF | 3D PDF |
| 16 | val_3DBETA | 3D beta |
| 17 | val_f2averb | <F^2> (variant b) |
| 18 | val_faver2b | <F>^2 (variant b) |
| 19 | val_faverb | <F> (variant b) |

- Formula selection is in the `qval` function (`prog/qval.f90`). The intensity
  uses `dsi` unless `laver` is set. Amplitude, phase, real and imaginary parts
  use `csf`.
- Value 15 (3DPDF) is the only path that uses `dsi3d` (see 7.5).

### 7.9 `reset`, `set_output`, and the 3D prep path

- `reset` at `prog/graphic.f90:513-515` region calls `output_reset`
  (`prog/graphic.f90:2816-2850`).
- `out_prep_3dpdf` (`prog/graphic.f90:2187-2301`) prepares the 3D PDF output.
  `out_prep_3dpdf_3d` (`prog/graphic.f90:2393-2430`) performs the forward FFTW
  transform of `dsi` and keeps the real part. Line mapping of complex to real
  data: `lib/map_1dtofield.f90` (`fftfdtoline_33D_C_R_DP`).

### 7.10 Other writers

- SHELX / list5 / list9 writers apply extinction (`prog/graphic.f90:1911`,
  `:1942`, `:1970`). The four_run amplitudes do not get extinction applied.
- `prog/discus_s_out_file.f90` holds `output_save_file_1d` and
  `output_save_file_2d` (header lines written around `:121-129`).
- `prog/discus_out_file.f90` is obsolete and commented out in `CMakeLists`
  (line 79).

## 8. 3DPDF menu (`discus_3dpdf.f90`)

### 8.1 Command dispatch

- `three_main` (`prog/discus_3dpdf.f90:55-208`) dispatches on:
  - `exit` (`:167`)
  - `help` / `?` (`:172-173`)
  - `reset` (`:184-186`), calls `three_reset` (`prog/discus_3dpdf.f90:1605`)
  - `data` (`:190-197`), stores the file name in `three_data`
  - `run` (`:200-202`), calls `three_run`
- The suite catalogue agrees: `/menus/discus/3dpdf/commands` lists
  `exit, help, ?, reset, data, run` (`discus_suite_6.23.01.json`).
- `three_load`, `three_find`, `three_show` and `three_inter` are not dispatched
  from the menu. They are called only from `three_run`.

### 8.2 `data`

- Stores `cpara(1)` in `three_data` after `do_build_name` (`:190-197`).
- Input file layout (1D, 2D or 3D data, header format): UNCONFIRMED. The
  branch at `prog/discus_3dpdf.f90:491` (`ndims == 2`) indicates a 2D case
  exists.

### 8.3 `run`

- `three_run` (`prog/discus_3dpdf.f90:210-227`):
  1. Deallocates and reallocates `three_peaks(5, three_max_peaks, 2)` and zeroes it.
  2. `three_load(three_data)` (`:229-287`): reads the kuplot HDF5 data through
     `do_load`, node 1.
  3. `three_find` (`:351-461`): peak search. Criteria in code; maximum 50 peaks.
  4. `three_show` (`:516-561`): writes to the output file. This is the only
     active write in the 3DPDF path.
  5. `three_inter` (`:614-782`): interpretation (inter-atomic vector matching,
     `three_inter_match`, `three_inter_rigid`, `three_inter_rigid_shift`, called
     at `:698`, `:774`, `:776`). No active file output; results stay in memory.

### 8.4 `reset`

- `three_reset` at `prog/discus_3dpdf.f90:1605`. Sets `linit = .true.`
  (`:162` and `:186`).

## 9. Feature table

Columns: feature | worth wrapping? | why | small live check.

| feature | worth wrapping? | why | small live check |
|---|---|---|---|
| `electron` (no energy) | Yes, with `energy` | Electron diffraction is a separate path; the energy default is not obvious (section 1) | Run `electron` then `show`; compare the derived energy with 0.2989 keV |
| `energy` | Yes | Needed to set a specific X-ray or electron energy; clears wavelength | Set `energy 200`, check the wavelength against 0.025079 A |
| `dispersion` | Yes, as a boolean | Already partly covered by `disp off`; bare form errors | `dispersion a` then `dispersion o` should succeed; bare `dispersion` should give -6 |
| `delf` | Maybe, after reading the body | Anomalous correction per element; syntax not traced | Run a `delf` on one element, then `show scat` |
| `scat` 10/12-param form | Yes | Only way to supply a user Gaussian form factor | Set one element with `scat <elem> ...` (10 params), then `show scat` and compare printed a, b, c |
| `scat` internal / equivalence | Yes | Lets one element reuse another element's curve | `scat Fe internal`, then `show scat` |
| `calculate` | Maybe | Single-point structure factor; useful for checks | `calculate 1 1 1` after `run`, compare with `value amplitude` |
| `hkl` | Maybe | Batch calculation from a list; extinction applies here | Create a 2-line hkl list and run; check the output file |
| `zone` | Yes | Strong defaults (electron, 200 keV, 512x512); silently changes state | `zone 0 0 1 0.1`, then `show`; check renergy 200 and radiation electron |
| `zone OFF` | Yes, with `zone` | Needed to return to normal mode | `zone OFF` after `zone`, check l_zone false |
| `layer` | Unknown | Body not traced | Read `prog/fourier.f90:402-434` first |
| `lots` | No, unless needed | Random-lot control; fine default | Leave default; not needed for current macros |
| `set technique` (not turbo) | Maybe | Alternate algorithms; syntax unverified | Try `set technique` with `help` output and inspect the list |
| `set filter`, `set conv` | Maybe | Post-processing of the transform | Run with a filter set and compare `value intensity` |
| `xray waas` / `xray inter` | Yes | Table choice changes results; default is International Tables | Run `xray waas` then `show scat` and compare f values for one element |
| `xray disc` | No, until checked | No RAD_DISC branch in `dlink` | Run `xray disc`, look for an error message |
| `show scat` | Yes | Confirms form-factor coefficients for the record | `show scat` after `scat`, compare |
| `output format` (stan, 3d, gnup, shelx, list5, list9, common, hdf5, xplor, vesta, vtk, mrc) | Yes for the used ones | Each writer has its own layout; several already used | Write one file per format, check the header lines |
| `output format pgm/ppm/post` | No | Disabled in 6.23.01 (error -6) | `output format pgm` should give -6 |
| `output input` | Maybe | Reads an existing grid; header layout looks incompatible with STAN files | Feed a STAN file to `input` and check the error |
| `output outfile` with mode | Yes | Needed to control output names and append | `outfile name add` then run twice, check file contents |
| `output range` | Only for XPLOR | Other writers ignore it | Set range, write XPLOR, compare |
| `output run` | Already used | Core write step | Covered by existing macros |
| `output show` | Low | Print only | Not needed |
| `output threshold` | No | Inert for active writers | Set threshold, write stan, confirm file is unchanged |
| `output value` (codes 1-19) | Yes | Selects which quantity is written; case-sensitive | `value amplitude` vs `value intensity`, compare output |
| `output value 3DPDF` | Yes | Needs `dsi3d`; error -180 if missing | Run without a 3D transform, expect -180 |
| `3dpdf data` | Yes, with `run` | Only input command | Set `data <file>` then `run` on a sample file |
| `3dpdf run` | Yes | Runs load, find, show, inter in one step | Run on a sample and check the peak list in the output |
| `3dpdf reset` | Yes, trivial | Re-initialises | `reset` then `run` |

## 10. UNCONFIRMED list

Items below are not verified in code, or their runtime behaviour was not
checked by running DISCUS.

1. Energy unit in `energy` (keV assumed from `get_wave`); compare with help text.
2. Electron table actually loaded by `get_scat_elec` (Peng et al. 1996 vs
   International Tables, as cited in `prog/appl_dis.hlp:850-860`).
3. Derived energy for `electron` without `energy` (about 0.2989 keV, computed
   from the formula, not run in DISCUS).
4. `delf` argument syntax and meaning of stored values (`prog/fourier.f90:298-349`).
5. `scat` 12-parameter layout: number of Gaussian terms and index layout of
   `cr_scat`.
6. `scat` absent-name behaviour in the Fourier path (`ier_num = -27` expected).
7. `scat` equivalence lookup at run time.
8. `calculate` dispatch between `calc_000` and `calc_hkl`.
9. `hkl` input file format.
10. `zone` pattern dependence on the random lot.
11. `layer` syntax and effect (`prog/fourier.f90:402-434`, body not traced).
12. `lots` options beyond LOT_BOX / LOT_OFF.
13. `set technique` values other than `turbo`; filter and conv `set` keywords.
14. `xray disc`: no RAD_DISC branch in `dlink`; runtime behaviour unknown.
15. Default table before any `xray` command (International Tables expected from
    `diff_table` initialisation); runtime confirmation.
16. `dispersion` effect on wavelength (metadata only per `docs/FOURIER_CONTROLS.md`).
17. `energy` / `wvlength` clearing logic (`prog/fourier.f90:922-941`).
18. `input` compatibility with STAN files that contain a NIPL header.
19. 3d NIPL writer when `out_inc(3) > 1` (file unit closed inside the loop).
20. Scale meaning inside `do_lanczos`.
21. 3DPDF input file layout for `data` (1D, 2D, 3D; header format).
22. Whether `three_load`, `three_find`, `three_show`, `three_inter` have any
    path from the menu other than `three_run` (none found by grep).
23. `help` text sign for `scat` (`exp(+b/4 h^2)`) versus code (`exp(-b s^2)`).
24. Runtime checks of all pitfalls in section 11.

## 11. Pitfalls

- `dispersion` bare form: help says allowed, code gives error -6
  (`prog/fourier.f90:353-369`).
- `scat` help sign versus code (section 0.3). Code is authoritative.
- `zone` silently sets electron radiation, 200 keV, and 512x512 grid
  (`prog/fourier.f90:983-1005`).
- `electron` does not reset `l_energy`; an earlier `energy` can persist.
- `neut`, `xray` and `electron` do not reset `l_energy`.
- `output format pgm`, `ppm` and `post` return error -6 in 6.23.01. Do not
  wrap them.
- `output threshold` is inert for writers in 6.23.01.
- `output range` affects only XPLOR.
- Extinction (`set exti`) applies to `calc_hkl` and to SHELX / list5 / list9
  writers only, not to `four_run` amplitudes.
- `output value 3DPDF` needs `dsi3d` (error -180 if not allocated).
- `output run` requires a previous Fourier (error -118).
- `output value` keywords are case-sensitive.
- `prog/fourier.f90:761-765`: the inner loop reuses the index `i` of the outer
  loop. No effect with one element; recorded for future edits.
- `prog/fourier_reset.f90` does not reset `diff_table` or `diff_exti`.
- `prog/discus_out_file.f90` is obsolete and commented out (`CMakeLists`
  line 79); do not use it as a reference.
