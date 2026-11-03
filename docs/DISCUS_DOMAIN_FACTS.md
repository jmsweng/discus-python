# DISCUS 6.23.01 domain menu: facts behind `discuspdf.discus_domain.Domains`

Source: `/home/claude/dc_6.23.01` (read-only). Paths below are relative to `discus/prog/`
unless they start with `lib_f90/`. Line numbers are those of the 6.23.01 files. Nothing here
was confirmed by running the DISCUS binary; "verified" means read from the source.

## 1. Pseudo atoms

- A pseudo atom is **not** an atom of the host crystal and needs no special name. It is an atom
  line of the input file (`inputfile`, a DISCUS `.stru`). The name is the first 4 characters
  (`line(1:4)`), the position is free-format in columns 5-52 (`domain.f90:2894-2913`).
- The name is matched against names given by `assign` (`domain.f90:2897`, `str_comp` with
  4 characters; names are upper-cased, `domain.f90:2894`). The name is stored in
  `CHARACTER(LEN=4)` (`domain_mod.f90:28`), so names are at most 4 characters.
- An input line whose name was never assigned stops the run with error -91 "Unknown pseudo atom
  name" (`domain.f90:2903-2906`). The input file therefore must hold only names assigned in the
  same domain menu session.
- The help text's "special names that define a domain or an object" (`appl_dis.hlp`, section
  `2 domain`, instructions) refers to `domain`/`molecule`/`object` header lines inside a
  content file (`domain.f90:2435` and `micro_read_atom`), not to pseudo atoms.
- The input file header is read with `stru_readheader` before the pseudo-atom lines
  (`domain.f90:1068`), so it is a normal `.stru` header followed by the pseudo-atom lines.

### Modes (`mode`)

| Word (first 4 chars compared) | Value | Source |
|---|---|---|
| `pseudo` / `pseudoatoms` | `CLU_IN_PSEUDO` = 1, the default (`domain_reset`, `domain.f90:3227`) | `domain.f90:227-233`, `domain_mod.f90:12` |
| `domain` | `CLU_IN_CLUSTER` = 0: explicit `domain ... domain end` blocks in the input file | `domain.f90:231`, `domain_mod.f90:11`, `domain.f90:740-762` |
| `irregular` | `CLU_IN_IRREG` = 2: all pseudo atoms together define one irregular domain | `domain.f90:235`, `1206-1208`, `domain_mod.f90:13` |

`mode` compares only 4 characters (`str_comp(cpara(1), ..., 1, lbef, ...)` with `lbef` the length of
the word `mode`, `domain.f90:227-236`), so `pseudoatoms` and `pseudo` are the same. The step
implements `pseudoatoms` only. `domain` and `irregular` are not implemented (their parsing is
in `micro_filereading`, `domain.f90:740-790`, and `domain_irregular`, `domain.f90:1206-1300`).

## 2. `assign` (domain menu)

Dispatch: `domain.f90:187` (`domain_assign`, `domain.f90:332-462`). The first parameter is the
keyword; the second is the pseudo-atom name, which is created on first use
(`domain.f90:364-387`).

**Quotes matter.** `get_params` keeps quote characters in the parameter (`lib_f90/get_params_mod.f90:67-68`,
both quote characters are copied), and `str_comp` then fails to match `"character"` against
`character` (`lib_f90/str_comp_mod.f90:6-57`). A quoted keyword is **silently ignored**, with no
error. The step therefore writes unquoted keywords (`assign character,DOM1,sphere`), while the
help shows quotes.

Keyword matching is by prefix: `str_comp(..., 'contentfile', 2, ...)` accepts `content`
(`domain.f90:401`), `orientation` accepts `orient` (`domain.f90:417`). The step writes the short
forms.

- **character** (`domain.f90:389-399`): `cube` -> -1, `cylinder` -> -2, `sphere` -> -3,
  `fuzzy` -> -4, `irregular` -> -5 (`domaindis_mod.f90:11-15`). Default -4 (`domain.f90:3219`).
  Other words are ignored silently.
- **content** (`domain.f90:401-408`): the content file name. It is read as a `.stru` by
  `micro_read_atom` (`domain.f90:2194`) with the same header checks as the host
  (`domain.f90:1083-1103`):
  - lattice constants must match the host within 1e-4 A and angles within 1e-3
    (`lmetric`, `domain.f90:1083`), otherwise error -87;
  - the first character of the space group must match the host (`domain.f90:1098`), otherwise -88.
  - Content coordinates are **fractional** coordinates of the content cell, with `u(4) = 1`
    (`domain.f90:2468-2471`). The atoms are placed by the orientation and translation
    (see 4 below).
  - A name starting with `internal` is internal storage (`domain.f90:221` for the input file,
    `domain.f90:2042` for content); the step refuses such names.
- **fuzzy** (`domain.f90:409-416`): the removal distance in A. Default **2.55 A**
  (`domain.f90:3220`), which applies if `assign fuzzy` is never given; the step always writes it.
- **shape** (`domain.f90:443-453`): `assign shape,<name>,<row>,<t1>,<t2>,<t3>,<t4>[,<sigma>]`.
  `row = werte(1)`; `t_j = werte(j+1)` for j = 1..4; `sigma = |werte(6)|` (`domain.f90:449-453`).
  Default all zero (`domain.f90:3222`).
  - Units: host fractional (unit-cell) coordinates. `micro_read_atom` multiplies the host
    fractional position by the inverse of the shape matrix (`domain.f90:2474-2475`), so a shape
    entry of 1 is one unit cell.
  - Clamp: in the 3x3 part, every entry is `MAX(shape * size_sigma, 0.001)` (`domain.f90:2918`).
    Negative entries cannot be represented, and zeros are stored as 0.001. The 4th column
    (translation) is not clamped (`domain.f90:2922`).
  - Sigma: `size_sigma = MAX(1 + gasdev(sigma), 0.01)` per row (`domain.f90:2916`), so sigma
    is a size distribution. Not exposed by the step.
- **orient** (`domain.f90:417-438`): `assign orient,<name>,<row>,<o1>,<o2>,<o3>,<t>`, the same
  row/value layout. Copied without clamping (`domain.f90:432`, `domain.f90:2920`, `2924`).
  Default all zero (`domain.f90:3221`). An all-zero 3x3 part collapses every content atom onto
  the pseudo-atom position, so the step always writes an orientation (identity by default).

## 3. The shape frame and the inside test

Per pseudo-atom line (`micro_read_simple`, `domain.f90:2826-2956`):

- Matrix `M` (4x4) has `M[i][j] = S[i][j]` (shape 3x3 part, clamped as above) and
  `M[i][4] = t_i + xyz_i` (shape translation plus pseudo-atom position) (`domain.f90:2918-2922`).
  `mc_idimen = inverse(M)` (`domain.f90:2933`).
- A content atom at fractional `u` is first placed at `v = O u + (t_o + xyz)` with
  `O` = orient matrix (`domain.f90:2474`: `v = matmul(mc_matrix, u)`), then the shape
  coordinates are `w = inverse(M) v` (`domain.f90:2475`).
- Inside tests on `w` (`domain.f90:2482-2518`):
  - `sphere`: metric length of `w` < |b| (`radius(2)`, `domain.f90:2342`);
  - `cube`: `|w_1|, |w_2|, |w_3| <= 1` (`domain.f90:2489`);
  - `cylinder`: `|w_3| <= 1` and the metric length of `(w_1, w_2, 0)` < |b| (`domain.f90:2503-2510`);
  - `fuzzy`: all atoms (`domain.f90:2517-2518`).
- So the default identity shape is a cube of half-width one cell, and a sphere of radius |b|
  (matches the help, `appl_dis.hlp` section `3 assign`).

## 4. What `run` does

`run` (`domain.f90:240-256`) calls `micro_filereading` (`domain.f90:605`). In pseudo mode
(`domain.f90:1112-1190`) it loops over the pseudo-atom lines and, for each one, inserts the
content atoms (`micro_read_atom`, `domain.f90:2194-2822`).

Per content atom:

1. Placed at `v = O u + (t_o + xyz)` in host fractional coordinates (`domain.f90:2474`); the
   position is stored as is, not wrapped into the cell (`domain.f90:2545`).
2. Inserted only if inside the shape (`domain.f90:2489-2518`). Inserted atoms are **appended**
   at the end of the atom list: `cr_natoms = cr_natoms + 1` (`domain.f90:2542`).
3. Type: an existing type is reused only if the name **and** the Debye-Waller factor are equal
   (`domain.f90:2555`, `dw1.eq.cr_dw(j)`). Otherwise a new type is appended to the host's list
   (`domain.f90:2587-2590`). Content atoms of the host elements with another B therefore add
   types.
4. Flags: `PROP_NORMAL` and `PROP_DOMAIN` set (`domain.f90:2559`, `2590`;
   `prop_para_mod.f90:11`, `13`).

Host atoms (`remove_strict` branch, `domain.f90:2776-2800`; runs only when `fuzzy > 0`):

5. Removal is by **distance to the inserted atoms**, not by shape membership. A host atom
   (indices `1 .. natoms_old`) within `fuzzy` of any inserted atom is set to type 0 (vacancy),
   `cr_iscat(1,i) = 0` (`domain.f90:3014` for fuzzy, `domain.f90:3075` for other shapes).
   The other shapes use a bounding-box prefilter (`domain.f90:3060`) and skip host atoms that
   belong to molecules (`domain.f90:3044-3052`).
6. The help says "All host atoms that are within the domain are removed" (`appl_dis.hlp`,
   `assign` section, line 3332). The source does not do shape-membership removal in this path.
   The shape-membership routine `micro_remove_old` (`domain.f90:2105-2190`) is never called from
   the source tree (grep of `discus/prog` and the rest of `dc_6.23.01` finds only its definition).
7. A warning (`ier_num = +1`) is raised if the shortest host-to-domain distance exceeds
   1.5 x `fuzzy` (`domain.f90:2793-2796`).
8. With `fuzzy = 0` nothing is removed, so host atoms overlap the domain.

Afterwards (`domain.f90:1314-1317`): `chem_purge = .TRUE.`, `chem_quick = .FALSE.`,
`chem_period = .FALSE.` (the crystal is no longer periodic; the help says the same).

### Result: atom count and ordering

- Vacancies stay in the list. A vacancy is an atom with `cr_iscat = 0`, and its type name is
  `VOID` (`discus_allocate_appl_mod.f90:973`). They are removed only by `purge`
  (`modify_mod.f90:1048`; help `appl_dis.hlp` line 8017 "Deletes all vacancies"). Purge is a
  discus-menu command, so a `Purge()` step after `Domains` removes them.
- Atom count: `N_result = N_host + N_inserted` before any purge. After a purge it is
  `N_host - N_removed + N_inserted`. Neither is in general a multiple of the host cell count.
- Order: inserted atoms follow the host list, so the cell-major order of `read_stru` is broken
  when any atom is inserted.
- Read-back (`discus_io.py`): `read_stru` requires `N % prod(n) == 0` (`discus_io.py:758-770`)
  and the cell-major order check (`discus_io.py:791-815`). Both are needed for
  `StructureJob.run` to return a `Structure`. Expected result, not verified by running DISCUS:
  `StructureJob.run` raises `ResultNotAStructure` for any `Domains` job, because domain atoms are
  appended after the host list (so order fails) and the count is not a multiple of the cells.
  `check_order=False` only skips the order check, not the divisibility check.
- Save: `save` writes every atom whose type passes `check_select_status` (`discus_save.f90:959`),
  and `sav_latom` is initialised to TRUE (`discus_save.f90:2024`). So vacancies are probably
  written as `VOID` lines. Not verified.
- New atom types can be added (see 4.3). Property is non-constant after removal (vacancies
  clear `PROP_NORMAL`), which the job tolerates (`discus_ops.py`, `run`: only `Property` may vary).

## 5. Random numbers

`micro_read_simple` calls `gasdev(sigma)` three times per matching pseudo-atom line
(`domain.f90:2916`). `gasdev` calls `RANDOM_NUMBER` unconditionally, including for sigma = 0
(`lib_f90/random.f90:11-49`; the cache `iset` only halves the draws). `Domains.uses_random_numbers`
is therefore True, and `StructureJob` requires a seed.

## 6. Other observations

- `reset` (`domain.f90:168`) calls `domain_reset`, which the catalogue summary describes as
  "Reset all Fourier definitions" (`suite_catalogue` entry for `discus/domain` `reset`). The
  summary is wrong; the command resets the domain definitions only.
- `inputfile` is a file name only (`domain.f90:217-222`); the file is read at `run`.
- The step writes one `reset` at the start, so one `Domains` step does not inherit an earlier
  one's definitions. Several `Domains` steps in one job have not been checked against the source;
  the domain tables persist between `run`s because `clu_number` is only reset in the
  cluster mode (`domain.f90:248`).

## 7. Not determined

- Whether a DISCUS 6.23.01 run writes the `VOID` atoms into `result.stru` (see section 4).
- Whether `StructureJob.run` raises for the result of a real `Domains` job (needs a DISCUS run).
- The exact content of `appl_dis.hlp` for the `shape` sigma and the `set` removal modes; not
  pinned down, not implemented.
- Whether DISCUS draws the RANDOM_NUMBER stream in the same order on every platform (only that the
  step needs a seed).
