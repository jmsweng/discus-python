# DISCUS 6.23.01: `patterson`, `inverse` and `diff-four` menu facts

Agent-written from the source (2026-11-05); the items used by `discuspdf.discus_patterson` were checked against
`patters.f90` and on the engine (section 9). Line numbers refer to `discus/prog` of DiffuseCode 6.23.01.


Source root: /home/claude/dc_6.23.01/discus/prog (unless noted). Line numbers are from the files as read.
CODE-READ = read from source. UNCONFIRMED = not verified by running code.

## 0. Shared code

- All three menus run the same subroutine `patterson(inverse_type)` (patters.f90:7), with the menu
  mode in the INOUT argument `inverse_type` (patters.f90:46).
- Mode constants: `INV_INV=0, INV_DIFF=1, INV_PATT=2, INV_FCALC=3` (inverse_mod.f90:7-10; INV_FCALC unused in the paths read).
- Value-type constants (intens_mod.f90:7-12): INTENSITY=1, AMPLITUDE=2, PHASE_ANG=3, REAL_PART=4, IMAG_PART=5, FCALC=6.
- Module state for the menu: patters_mod.f90 (defaults at lines 10-33). Reciprocal layer (eck, vi, inc) is shared with the
  main `fourier` menu (fourier_reset.f90:68,75-81).
- Run pipeline: `patterson_check_values` (920) -> `patterson_run` (1175) -> `do_patters` (1305) -> `calc_patters` (1813)
  with `set_patt_value` (1905) or `set_patt_value_s` (1994, SHELXL and difference path).
- Command-matching: `str_comp` (lib_f90/str_comp_mod.f90:6-57) compares the first min(len) characters, requires the
  typed length >= j, and accepts a typed word that is longer than the key if the key is its prefix.
  Consequences: 'diff-four' is accepted by the 'diff' branch (discus_kdo.f90:244, key 'diff' min 2);
  bare 't' hits 'top' (patters.f90:174) before 'type' (691); bare 's' hits 'statistics' (660) before 'scale','set','show','sign'.

## 1. How each menu is entered

| Menu | Main-menu command | Code | Mode value | Notes |
|---|---|---|---|---|
| inverse | `inverse` (min 3 chars) | discus_kdo.f90:312-315 | INV_INV = 0 | sets lout_rho=.true. (313-314); gated by linteractive/lblock/lmakro (312) |
| diff-four | `diff` (min 2 chars); 'diff-four' also matches by prefix | discus_kdo.f90:244-247 | INV_DIFF = 1 | gated the same way (244); help documents 'diff-four' (appl_dis.hlp:1078, 3143) |
| patterson | `patterson` (min 2 chars) | discus_kdo.f90:335-338 | INV_PATT = 2 | lout_rho=.true. (337) |

- The mode variable is passed by reference, so `format hklf4` (patters.f90:277) changes it to INV_PATT for the rest of
  the session in that menu, even inside `inverse` or `diff`. CODE-READ, see Bug B4.
- Prompt suffix: '/inverse', '/diff-four', '/patterson' (patters.f90:77-83); restored at 914.
- `help` dispatch: 289-305 (`discus inverse`, `discus diff-four`, `discus patterson`; the offsets 15/17/17 match the strings).
- At entry: `diff_lsingle = .TRUE.` (67) and `patt_scale` = cr_ncatoms/cr_natoms when both counts > 0 (70-71); else 1.0 if it is 0 (72-73).

## 2. Commands and their arguments (patterson(), patters.f90 lines 154-908)

Abbreviation minimums are the `j` argument of str_comp; the key length is the third argument.

| Command (key, min) | Lines | Arguments as parsed | Effect |
|---|---|---|---|
| abscissa (1) | 154-170 | h/k/l -> extr_abs=1/2/3; else -6 | reciprocal axis for output |
| top (1) | 174-190 | h/k/l -> extr_top | |
| continue (1) | 194 | macro_continue | |
| echo (2) | 199 | echo | |
| evaluate (2) | 204 | do_eval | |
| exit (3) | 209 | leaves loop | |
| file (1) | 214-237 | `file a|b,<name>` (ianz<=2, else -6 at 220); a/b selector; do_build_name (226/232 use cpara(2)) | rho_file(1|2); else -41 (235) |
| format (1) | 241-285 | `gnup` -> ftyp=3 (249); `stan` -> ftyp=0 (254); `shelx[,fobs|fcalc]` -> ftyp=4, rho_type(1)=AMPLITUDE (260), second arg fobs->AMPLITUDE, fcalc->FCALC (265), else -52 (267); `hklf4` -> ftyp=5, rho_type(1)=INTENSITY (276), inverse_type=INV_PATT (277); else -9 (279) | |
| help / ? | 289-305 | `help errors` or menu help | |
| layer (2) | 309-327 | 11 values: 3x3 eck columns (ll,lr,ul), inc(1)=w(10), inc(2)=w(11) | |
| ll (2) | 331-346 | 3 reals -> eck(:,1) | |
| lr (2) | 350-365 | 3 reals -> eck(:,2) | |
| nabs (2) | 369-384 | 1 int > 0 -> inc(1) | shared inc |
| nord (2) | 388-405 | 1 int > 0 -> inc(2) | |
| ordinate (2) | 428-441 | h/k/l -> extr_ord; bad letter -6; **wrong count: silently ignored** (no ELSE) | |
| ntop (2) | 409-424 | 1 int > 0 -> inc(3) | |
| rhoabs (4) | 445-461 | x/y/z -> rho_extr_abs | |
| rhoord (5) | 465-478 | x/y/z -> rho_extr_ord; wrong count silently ignored | |
| rhotop (4) | 482-498 | x/y/z -> rho_extr_top | |
| rholayer (5) | 502-520 | 11 values -> rho_eck(:,1..3), rho_inc(1)=w(10), rho_inc(2)=w(11) | |
| rholl (5) | 524-535 | 3 reals -> rho_eck(:,1) | |
| rholr (5) | 540-551 | 3 reals -> rho_eck(:,2) | |
| rhotl (4) | 555-566 | 3 reals -> rho_eck(:,4) (not in help) | |
| rhoul (4) | 570-581 | 3 reals -> rho_eck(:,3) | |
| rhona (5) | 586-602 | 1 int > 0 -> rho_inc(1) | |
| rhono (5) | 607-623 | 1 int > 0 -> rho_inc(2) | |
| rhont (5) | 628-644 | 1 int > 0 -> rho_inc(3) | |
| reset (2) | 648-649 | fourier_reset | see section 6 |
| run (2, key 'run ') | 653-656 | patterson_check_values, then patterson_run | |
| statistics (1) | 660-687 | needs cr_v>0 (661); `screen` (e_io=6) or `file <name>` (e_io=38, do_build_name); anything else -6 (676-678); calls e_create(e_io,cpara(2),lpara(2),.TRUE.) | writes <name>.ehkl etc. |
| type (1) | 691-728 | exactly 2 params (694); a->i=1, b->i=2 else -41 and **RETURN from patterson()** (704-706); types: intensity, amplitude, phase, real, imaginary, fobs (->AMPLITUDE), fcalc (->FCALC); else rho_type(i)=0 and -42 (723-726) | |
| scale (2) | 733-745 | 1 real -> patt_scale | |
| set (2) | 749-820 | see below | |
| sign (2) | 824-835 | `inverse` -> patt_sign=-1; `fourier` -> +1; else -6 | **no effect**, see Bug B1 |
| show (2) | 839-841 | patterson_check_values + patterson_show | |
| tl (1) | 845-857 | 3 reals -> eck(:,4) | |
| ul (1) | 861-874 | 3 reals -> eck(:,3) | |
| wilson (1) | 878-886 | needs cr_v>0; wilson_calc | |
| system (2) | 890-896 | OS command | |
| wait (3) | 901-903 | do_input | |
| other | 906-908 | -8 | |

`set` sub-commands (patters.f90:749-820):
- `set accu init|add` (757-767) -> patt_accu = PATT_INIT(0) / PATT_ADD(1). Default 0 = init.
- `set excl none|<value>` (768-781) -> patt_excl9999 false / true with patt_excl_val = second value (the code passes cpara(1)='0' so ber_params returns the value in werte(2)).
- `set mode normal|sharp|super` (782-795) -> patt_mode 0/1/2.
- `set origin normal|subtract` (796-804) -> patt_origin 0/1.
- `set rsym apply|ignore` (805-811) -> patt_rsym .true./.false. Any other value: **silently ignored** (no ELSE).
- other first word -> -6 (812-815). Wrong count -> -6 (816-818).
- Prefix checks: 'apply' min 1 length 5 (809) matches 'appl', as help says.

## 3. Input files

### 3a. Which files each menu needs (patterson_run, patters.f90:1175-1301)
- Common: `cr_v <= 0` -> -35 (1206-1212). Array allocation (1214-1227). `MAXQXY(1..2)` check on rho_inc -> -71 (1229-1233; rho_inc(3) not checked).
- Swap (1235-1244): if rho_type(2) is INTENSITY, AMPLITUDE or REAL_PART, files and types of a and b are swapped. This runs for all three menus.
- ftyp=4 (shelxl): one file; rho_type(1) must be AMPLITUDE or FCALC, else -51 (1248-1253). ltwo_files=.false.
- ftyp!=4 for INV_INV (1256-1265): allowed pairs (INTENSITY,PHASE_ANG), (AMPLITUDE,PHASE_ANG), (REAL_PART,IMAG_PART); else -36. Two files.
- ftyp!=4 for INV_PATT (1266-1274): allowed rho_type(1) = INTENSITY, or AMPLITUDE, or (REAL_PART with rho_type(2)=IMAG_PART); else -38. ltwo_files stays TRUE unless rho_type(1)=INTENSITY (1275-1276).
  - For INTENSITY with sharp/super: if ftyp=5 -> `e_create(36, rho_file(1), .FALSE.)` (1278-1282) to build e_aver_f2; else -38 (1283-1286).
- ftyp!=4 for INV_DIFF (1290-1292): one file, no type check.

### 3b. Formats (do_patters, patters.f90:1305-1809; reading is by unit ifa=21, ifb=22, line 1335)

- **stan (ftyp=0), 2D (inc(1)>1 and inc(2)>1)**: KUPL-style. Header line 1: `nx ny` (1697); line 2: `xmin xmax ymin ymax` (1698-1699). Second file must have identical header (1714-1716), else -22 (1744). Then ny rows of nx values (1718-1723). Header min/max values are read but used only for the equality test. Loops use inc(1) (1730) and ny(1) (1718), not nx(1) or inc(2) consistently.
- **stan (ftyp=0), 1D (inc(1)=1 or inc(2)=1)**: lines `h value` (1670), one per point, abscissa value read and discarded (h(extr_abs) overwritten at 1680). Second file same (1672).
- **gnup (ftyp=3)**: no header. Blocks of inc(1) lines, each `x y value third` (1754-1755), with the value in column 3 (z1(1)) and the third index in column 4 (zz1) (1762-1764 assign h(extr_ima)=zz1). Blocks separated by one blank line (1777-1787). Two-file match requires x, y and column 4 equal (1761).
- **shelxl (ftyp=4)**: `h k l Fobs Fcalc phase` free format (1413, phase in degrees); fallback fixed format `(3I4,2F10.2,F7.2)` (1415, 1417). Column use: zz1=Fobs, zz2=Fcalc, zz3=phase (1426, 2016-2018, 2029-2030).
- **hklf4 (ftyp=5)**: `h k l I` read as integers and one real (1496-1498), fallback `(3I4,F8.2)` (1498). Sigma is ignored. Read stops at EOF (END=40/60). Reflection list defines the data; no grid is used.
- **Reciprocal grid**: for stan/gnup it is the reciprocal layer: h = eck(:,1) + (i-1)*vi(:,1) + (j-1)*vi(:,2) (1680, 1735-1736), where vi is computed in patterson_check_values (949-951) from eck and inc. For shelxl/hklf4 the grid is the reflection list itself.

### 3c. Input-type meaning (`type a|b`, `format shelx,...`)
- Intensity: |F|^2 (used directly for Patterson, or as I for SHARP/SUPER with hklf4).
- Amplitude: |F|; Patterson squares it (1978, 1987, 2023-2024); inverse combines with phase.
- Real/imag: Patterson uses Re^2+Im^2 (1972, 1985); inverse uses Re + i Im directly (1940).
- Phase: degrees, used with cosd/sind (1950-1951, 1954-1955).
- fobs/fcalc: map to AMPLITUDE/FCALC (719-722).
- The help table (appl_dis.hlp:6092-6096) says Patterson 'amplitude' needs only file a. The code still requires file b for AMPLITUDE (see Bug B12).

### 3d. Files per menu
- inverse: stan/gnup need 2 files (INTENSITY+PHASE, AMPLITUDE+PHASE or REAL+IMAG); shelxl needs 1 file (AMPLITUDE or FCALC).
- diff-four: stan/gnup nominally 1 file (no INV_DIFF branch, see Bug B3); shelxl 1 file (Fobs-Fcalc).
- patterson: intensity 1 file (stan/gnup/hklf4); amplitude 2 files in code (help says 1); real+imag 2 files; shelxl 1 file; hklf4 1 file.

## 4. Computation in run

### 4a. Real-space grid (calc_patters, 1852-1860, 1872-1885)
- r(j,i,k) = rho_eck(:,1) + (j-1)*rho_vi(:,1) + (i-1)*rho_vi(:,2) + (k-1)*rho_vi(:,3).
- rho_vi(:,a) = (rho_eck(:,a+1) - rho_eck(:,1)) / (rho_inc(a)-1), for a=1,2,3 with corner columns 2 (lr), 3 (ul), 4 (tl) (patterson_check_values 952-954).
- rho_eck column 1=rholl, 2=rholr, 3=rhoul, 4=rhotl.

### 4b. Summation (calc_patters)
- csf(j,i,k) = csf(j,i,k) + exp(+2*pi*i*h.r(j,i,k)) * A(h), with the phase handled as an integer table index (1857-1880). cex(i)=exp(+i*2*pi*i/I2PI) (fourier_sup.f90:2602-2632), so the exponent sign is **+** (CODE-READ; see Bug B1).
- A(h) is the per-reflection value times patt_scale * vr (patters.f90:1428, 1463, 1475, 1544, 1571, 1648, 1686, 1738, 1768). vr = cr_vr (1297), the reciprocal cell volume (class_crystal.f90:59 comment; crystal_mod.f90:49 default 1.0D0). So the 1/V factor appears as multiplication by the reciprocal volume.
- Normalisation therefore: rho = cr_vr * patt_scale * sum_h A(h) exp(+2 pi i h.r) (UNCONFIRMED against a reference result).

### 4c. Per-reflection values
- INV_INV (set_patt_value 1930-1956): (REAL,IMAG): A = Re + i Im. (INTENSITY,PHASE): A = sqrt(I)(cosd phi + i sind phi), negative I skipped (1946). (AMPLITUDE,PHASE): A = |F|(cosd phi + i sind phi).
- INV_PATT (set_patt_value 1957-1989): INTENSITY -> A=I; AMPLITUDE -> A=|F|^2; (REAL,IMAG) -> A=Re^2+Im^2.
- INV_PATT shelxl (set_patt_value_s 2020-2027): INTENSITY and AMPLITUDE -> zz1^2 (Fobs^2); FCALC -> zz2^2.
- INV_INV shelxl (2014-2018): AMPLITUDE -> Fobs e^{i phase}; FCALC -> Fcalc e^{i phase}.
- INV_DIFF shelxl (2028-2030): A = (Fobs - Fcalc) e^{i phase}. **No INV_DIFF branch in set_patt_value** (1905-1992), see Bug B3.
- hklf4 (1508-1544 rsym, 1608-1648 non-rsym), using intensity zz1:
  - NORMAL: A = I*wilson_scale (1515, 1617); with origin subtract: A = I*wilson_scale - fj2(d*^2) (1513, 1613).
  - SHARP: E = sqrt(|I|)/e_aver_f2(j) (1624 non-rsym; 1522 rsym, see B5); A = sqrt|I|*E*sqrt(wilson_scale) (1526/1529, 1626/1629); with origin subtract minus sqrt(fj2) (1526-1527, 1626-1627). Guard e_aver_f2>0 only in non-rsym (1623).
  - SUPER: E^2 (1536/1540, 1639/1643), with origin subtract E^2-1 (1538, 1641).
  - d* for bins: dstar from do_blen; bins from e_hist (2898-2939), 20 bins, limits 0.2 ... 3.8, 1e6.
- fj2(d*^2) = sum over atoms of form(...)^2 (patters.f90:3225-3248).

### 4d. Symmetry and Friedel (rsym, set rsym)
- shelxl, rsym=.FALSE. (1405-1430): one term per reflection, no symmetry, no Friedel.
- shelxl, rsym=.TRUE. (1432-1480): loop over rec_n_sym operations (rec_sym), y(i)=sum_j R(i,j)*h(j) (1448-1453); phase = sum_i y(i)*t(i) (1455); zz4 = zz3 + 360*phase (1457); if acentric, also -y and -zz4 (1468-1477).
- hklf4, rsym=.TRUE. (1481-1587): loop over ops, y = R h (1548-1553); no phase used (intensity); if acentric also -y (1565-1572).
- hklf4, rsym=.FALSE. (1589-1659): one term per reflection, no expansion, no Friedel pair (help 6146-6152 says 'all' means no expansion).
- 1D/2D/gnup (1665-1791): no symmetry at all.

### 4e. Exclusion
- set_patt_value (INV_INV 1931-1937; INV_PATT 1958-1979) exclude points equal to patt_excl_val (exact equality, -9999 default when set).
- set_patt_value_s (SHELXL) and the hklf4 branches (1495-1649) do not test patt_excl9999 (see B10).

### 4f. Accumulation
- patt_accu=PATT_INIT: csf(1:rho_inc(1),1:rho_inc(2),1:rho_inc(3)) zeroed (1389-1390).
- PATT_ADD: no zeroing; the previous csf is used. Requires same rho_inc (not checked).

### 4g. Origin subtraction and sharpening summary
- `set origin subtract` is used only in the hklf4 branches (1510, 1523, 1537, 1610, 1625, 1640). It has no effect for stan/gnup/shelxl (B9).
- `set mode sharp|super` is used only in the hklf4 branches (1509, 1518, 1532, 1609, 1618, 1633). For INTENSITY+stan/gnup with sharp/super it errors (-38, 1283-1286); for other types, or shelxl, it is ignored.

## 5. Output

- `run` writes nothing itself. It sets four_was_run=.TRUE. (patters.f90:1298).
- Output is produced by the `output` command (discus_kdo.f90:330-331 -> do_niplps(lout_rho) at graphic.f90:7). The `run` branch there is gated by four_was_run (graphic.f90:515-516).
- set_output(linverse) (graphic.f90:2131-2170) copies rho_eck, rho_vi, rho_inc and rho_extr_* into out_eck, out_vi, out_inc, out_extr_*. For linverse it copies eck columns 1-3 only (2151-2155), not column 4, and does not set out_extr_top (2148-2149).
- do_output (graphic.f90:1518) reads qval(...) (e.g. 1644, 1842) and writes to `outfile` opened at 1856 (name from the output menu; exact default not traced, UNCONFIRMED).
- `wilson` (hklf4 only) writes `<name>.wilson` (fourier name without '.hkl' + 'wilson'; patters.f90:3135-3144) and prints the scale to the output unit (3212-3217). It is also called from do_patters for hklf4 input (1374-1376).
- `stat` writes, if 'file <name>' is used (e_io=38): <name>.ehkl (2520-2521, 2560-2563), <name>.statistics (2682), <name>.short (2740), <name>.hkl.histogram etc. (2798-2875). Screen (e_io=6): no file.
- Progress counters are written to output_io during long loops (1669, 1719, 1752).

## 6. Start values and what reset restores

### 6a. Program start (patters_mod.f90:10-33)
- rho_file = blanks; patt_accu=0 (INIT); rho_inc=(121,121,1); rho_type=(4,5) = (REAL_PART, IMAG_PART); patt_sign=-1; patt_mode=0 (NORMAL); patt_origin=0; ftyp=0 (stan); patt_excl9999=.false.; patt_rsym=.false.
- rho_eck = (0,0,0 | 1,0,0 | 0,1,0 | 0,0,1) (column = corner 1..4).
- rho_vi = diag(0.02, 0.02, 0.02) (not recomputed until run/show).
- patt_scale=1.0 (overwritten at each patterson entry, see 1.); patt_excl_val=-9999.0; e_aver_f2=0; wilson_scale=1.0; wilson_b=0.0.
- Reciprocal layer default (from fourier_reset, also used at first use): eck=(0,0,0 | 5,0,0 | 0,5,0 | 0,0,0), vi=(0.05,0,0 | 0,0.05,0 | 0,0,0), inc=(121,121,1), lambda='MOA1', rlambda=0.709260, renergy=17.480782.

### 6b. `reset` (patters.f90:648-649 -> fourier_reset.f90:7-116)
- Restores the shared Fourier state: csf/tcsf/acsf zeroed (30-32); alloc to size 1 (18-20); eck, vi, inc, lambda, four_* options, four_was_run=.FALSE. (53), zone and radiation settings, cr_delf* (108-114).
- Does NOT touch patters_mod: rho_eck, rho_vi, rho_inc, rho_type, rho_file, ftyp, patt_*, patt_scale, wilson_scale/wilson_b, e_aver_f2 keep their values. (CODE-READ; this is a trap.)
- Does not reset the mode or the prompt.

## 7. Bugs and traps (CODE-READ unless marked)

B1. `sign` has no effect. calc_patters sets `patt_sign = 1.0D0` at patters.f90:1840 on every call, so the inverse (-1) setting is lost. patt_sign is used at 1858-1862 only. The exponent sign is therefore fixed at + (see 4b). Sign matters for inverse and diff-four maps (Patterson is even, so unaffected). Severity: high for inverse/diff-four.

B2. Phase base not reset between rows (calc_patters 1872-1885). The base is reset only inside the i loop at 1883: `iarg = iarg0 + (iincu*(j-1) + iincv*i)`. When the j loop advances, iarg is not set to iarg0 + iincu*j; it stays at iarg0 + iincu*(j-1) + iincv*rho_inc(2). The error is iincv*rho_inc(2) - iincu (mod the table period). The commented old code (1891-1899) reset with `iarg0 + iincu*j*patt_sign` at the j level, which is the intended value. UNCONFIRMED numerically; a test with rho_inc(1)>=2 and rho_inc(2)>=2 against a direct sum would confirm. Severity: high if confirmed (affects every j>=2 column).

B3. diff-four has no branch for non-shelxl input. set_patt_value (1905-1992) has INV_INV and INV_PATT only, so for stan/gnup input (1683, 1731, 1765) a_b is never assigned in the INV_DIFF case. a_b is initialised once to (0,0) at 1401 and remains so, giving a zero map (UNCONFIRMED by run). Only shelxl reaches the INV_DIFF formula (2028-2030). Also, format hklf4 converts the menu to INV_PATT (277), so diff-four with hklf4 is not reachable.

B4. `format hklf4` changes inverse_type to INV_PATT (patters.f90:277), silently, in any menu. The help for hklf4 says the optional 'all' parameter exists (help 6135-6152). The code does not parse 'all' (274-277); the second parameter is ignored and INV_PATT is always set.

B5. Inconsistent SHARP/SUPER normalisation for hklf4 with rsym. Non-rsym (1624, 1639) and e_create's own E (2554) use E = sqrt(|I|)/e_aver_f2 (e_aver_f2 = sqrt(mean I), set at 2479). The rsym branch uses e_f = sqrt(|I|/e_aver_f2) (1522, 1536), which is not E. The rsym branch also has no guard against e_aver_f2 = 0 (guard exists at 1623, 1638). Severity: high for sharp/super with rsym.

B6. zz4 in the hklf4 rsym branch is computed from zz3 (1557, 1569), but zz3 is never read in the hklf4 path (only zz1 is read at 1496-1498). zz4 is not used in any value. Dead, misleading.

B7. `lll` (1341, 1586) is uninitialised and used only in a dead check. Harmless.

B8. 1D input with inc(1)=1 and inc(2)>1: the loop index i runs over inc(2) but h uses vi(:,1) (1680). vi(:,1) = (eck(:,2)-eck(:,1))/(inc(1)-1) is divided by zero (949, Inf/NaN). Should use vi(:,2). CODE-READ; severity medium for that geometry.

B9. `set origin subtract` and `set mode sharp|super` are silently ignored for stan/gnup/shelxl input (no error), except INTENSITY+stan/gnup with sharp/super, which errors (1283-1286). Help 6298 and 6308-6309 says they apply only to HKLF4, so the code matches the help except that it does not warn.

B10. `set excl` is not applied to shelxl (set_patt_value_s has no exclusion argument, 1994-2033) or to any hklf4 branch (1495-1649). Help 6282-6287 describes it generically.

B11. Silent truncation. Every read error or EOF in do_patters jumps to label 900 with ier_num unchanged (1670, 1697, 1722, 1754, 1777, 1795), so a short file gives no error. SHELXL and HKLF4 loops use END= jumps (1412, 1439, 1495, 1596) with the same effect.

B12. AMPLITUDE Patterson requires file b in the code. For INV_PATT, only INTENSITY sets ltwo_files=.FALSE. (1275-1276); AMPLITUDE (1266-1274) still opens and reads file b (1383, 1671-1673). The help (6092-6095) lists amplitude with file a only. Severity: medium (user must supply a dummy b).

B13. `type` with a bad selector returns from patterson() (patters.f90:706 RETURN, not cycle). The menu is left and the prompt is not restored (914 is skipped). Other error paths use cycle loop_main.

B14. `file` uses the length of the second word (lpara(2)) in the str_comp test (222, 228) instead of lpara(1). Low severity: the match is capped at the file-name length.

B15. `ordinate` (428-441) and `rhoord` (465-478) silently ignore a wrong parameter count (no ELSE). `abscissa` and `top` report -6.

B16. `set rsym` with any value other than apply/ignore is silently ignored (805-811). Other set values report -6 (812-815).

B17. Help vs code, statistics: help 6349 `stat {"screen" | <outputfile>}` but code accepts only the keywords screen or file (671-678); a bare file name gives -6.

B18. Help vs code, shelxl and hklf4 symmetry: help 6128-6129 and 6143-6144 say all reflections are expanded by symmetry. The code expands only when `set rsym apply` (1405 vs 1432; 1481 vs 1589), and patt_rsym defaults to .FALSE. (patters_mod.f90:26), which format does not change.

B19. Help vs code, gnup: help 6109-6112 describes 'gnu' as an output file type. The code reads it as input (248-249, 1747-1791).

B20. Help vs code, wilson: help 6375-6376 says Icalc = k**2*Fobs. The code prints Fcalc = sqrt(wilson_scale)*Fobs and Icalc = wilson_scale*Iobs (3213-3215).

B21. wilson_scale factor (patters.f90:3202). The regression intercept is (m*sumy - p*sumx)/(number*m) with q = ln(fj2/<I>) (3189, 3200-3202), i.e. ln of an intensity ratio. wilson_scale is then set to exp(0.5*intercept), which is the square root of the intensity scale, while the print (3213-3215) treats it as the intensity scale. Possible factor-2 error in the log scale. UNCONFIRMED numerically; a test on a known file would confirm.

B22. wilson_calc uses a list-directed read `READ(ifa,*) h, rint, sigma` (3151) with no fixed-width fallback (e_create 2311-2317 and do_patters 1413-1416 have one). Fixed-format HKLF4 files stop silently at err=902 (3151) and wilson_scale keeps its old value. wilson_calc errors (-98/-99 at 3102, 3110) are not checked by do_patters before it continues (1375-1378); whether oeffne at 1378 resets ier_num is UNCONFIRMED.

B23. e_create pass 2 (patters.f90:2391-2466) runs `DO i=1,n` where n counts only non-blank records (pass 1 cycles blank lines at 2305-2307 before n=n+1 at 2318), but blank lines in pass 2 also consume an iteration (2407). A file with blank lines loses its last records silently (UNCONFIRMED impact). The early RETURN at 2401-2404 on read error leaves e_aver_f2 unnormalised (sums, not RMS) with no ier_num.

B24. e_create writes column zz2 (2562) that is never read for hklf4 (2540 reads h and zz1; zz2 stays 0.0 from 2252). The division zz2/(e_aver_f2(j))**2 (2562) has no guard, so an empty bin gives Inf/NaN in <name>.ehkl.

B25. Degenerate grid. patterson_check_values divides by (inc-1) and (rho_inc-1) (949-954). With the default reset, inc(3)=1 and eck(:,4)-eck(:,1) = 0 gives vi(:,3) = 0/0 = NaN (fourier_reset.f90:68, 75-78; printed by `show`, 1058-1088). With rho_inc(3)=1, rho_vi(:,3) = (0,0,1)/0 = Inf (patters_mod.f90:19, 27; 954), used in calc_patters at 1855, 1860 through nint/int of Inf/NaN (UB). The k loop runs once so the effect is UNCONFIRMED.

B26. KUPL 2D branch (1694-1746): header xmin/xmax/ymin/ymax read (1698-1703) and not used beyond the equality test; loops use inc(1) (1730, 1735) while reading nx(1) values per row (1722-1723). If inc(1)!=nx(1), stale z1 values are used (UNCONFIRMED impact). The row loop uses ny(1) (1718) while the itic uses inc(2) (1717).

B27. GNU branch (1747-1791): the blank-line handling at 1777-1787 has an unreachable check at 1782-1786 (inside the blank-line IF). A non-blank line after a block is consumed silently.

B28. hkl column for stan/gnu 1D: abscissa value read and discarded (1670, overwritten at 1680). Not a bug, but the input's abscissa column is ignored.

B29. e_create and do_patters use a mixture of h as REAL (e_create 2070) and h as INTEGER (wilson_calc 3081). Consistent with formats; noted only.

B30. patterson_run swap (1235-1244) runs for every mode. For diff-four (only file a is read) a file b with type intensity/amplitude/real would silently replace file a. CODE-READ.

B31. The help `rhotl` command is not documented (code 555-566). Minor.

B32. Non-rsym hklf4 does not add Friedel mates (1589-1659), while the help implies the file is a full set (6146-6152). Consistent with 'all', noted only.

## 8. Quick reference of the three menus

- inverse (INV_INV): 2 files (stan/gnup) or 1 file with amplitude/fcalc (shelxl). Value: complex Fourier from F with phases. Output via `output`. Sign: no effect (B1). Phase-base bug (B2) affects all rows j>=2.
- diff-four (INV_DIFF): 1 file, shelxl only in practice (B3). Value: Fobs - Fcalc with phase. Sign no effect (B1).
- patterson (INV_PATT): intensity (1 file), amplitude (code wants 2, B12), real+imag (2), shelxl (1), hklf4 (1; forces INV_PATT, sets intensity). Value: |F|^2 or E-based sharpened values. Output via `output`. Sign: no effect (B1).

## 9. Probe results on the stock engine (2026-11-05, `tools/probes/patterson_probes.py`, tasks 025-028)

- **The stan reader cannot skip `#` lines.** A NIPL file written by DISCUS's own `output` starts with 26 `#` header
  lines; `patterson`/`inverse`/`diff-four` `run` then fails with "Error reading file" (-3), because the reader starts
  with the `nx ny` line (patters.f90:1689-1723). Copies without the `#` lines are read.
- **Patterson and inverse follow** rho = s / V_cell * sum_h A(h) exp(+2 pi i h.r), s = cr_ncatoms / cr_natoms (1/8 for
  2x2x2 cells), A = I (Patterson) or Re F + i Im F (inverse), to 3e-5 of max |rho| at every point except Bug B2's.
  The map file has rows = ordinate points, as the fourier output.
- **Bug B2 (row phase) confirmed**: the first ordinate point of every abscissa column after the first is wrong
  (errors up to 0.95 of max |rho|); every other point is right.
- **diff-four with stan input writes an all-zero map** (B3 confirmed).
- `sign -1` is refused ("Missing or wrong parameters", -6); its value is overwritten in every run anyway.
