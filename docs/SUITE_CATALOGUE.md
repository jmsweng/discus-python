# DISCUS_SUITE command catalogue

`discuspdf.suite_catalogue` lists every command the DISCUS_SUITE programs accept — the
suite level, DISCUS, KUPLOT, REFINE, DIFFEV, MIXSCAT and the common command language — and
checks macros against it. It is the inventory for extending discus-python to the whole suite:
which commands exist, which menu each one opens, and which ones discus-python can already
produce ([`SUITE_COVERAGE.md`](SUITE_COVERAGE.md)).

The catalogue says what the engine **accepts**. It does not say what a command does
physically; that needs a live engine check (see the verification levels below).

```python
from discuspdf.suite_catalogue import load_catalogue

cat = load_catalogue()                 # DISCUS_SUITE 6.23.01, 44 menus
cat.resolve("discus", "four").enters   # 'discus/fourier'
cat.menu("discus/fourier").command("set").syntax[:2]
# ['set "aver" [, <per>]', 'set {"external"|"internal"}']
for step in cat.check(macro_text):     # lines the engine would reject as unknown
    print(step.line, step.menu, step.text)
```

## Where it comes from

`tools/suite_catalogue/build_catalogue.py` reads a DiffuseCode checkout:

- **Fortran dispatchers** (the source of truth). Each menu routine compares the command word
  with `str_comp(befehl, 'name', min, lbef, len)` in an IF/ELSEIF chain; a sub-menu appends a
  suffix to the prompt. The scanner records each branch (name, minimum abbreviation, the
  `len` argument, source line), the keywords tested inside it (`str_comp(cpara(n), ...)`), the
  `name:value` options of `get_optional`, the routines it calls, and the menus a routine
  falls through to (`kdo_all`, `discus_kdo_common`, `kuplot_kdo_common`).
- **Help files** (`*.hlp`) add the syntax lines and a summary. Help topics with no matching
  command are kept as `help_only` (e.g. DISCUS's `mc`, which has help but no command in
  6.23.01).

```
python tools/suite_catalogue/build_catalogue.py <DiffuseCode> src/discuspdf/data/suite_catalogue/discus_suite_6.23.01.json \
    --version 6.23.01 --commit 19a214a659f1e4e3be0b4540f8a3bcb79ff5548e --label stock
```

The bundled file was built from tag `v.6.23.01` (`19a214a`). Built from the owner's fork
`master` (`a85c9de`), it lists the same menus and commands with the same abbreviations,
keywords and options: the fork changes behaviour (see below), not the command set.

## How the engine reads a command

These rules are implemented in `str_comp` and `Catalogue.resolve`, and checked by
`tests/test_suite_catalogue.py`:

- **Matching** (`lib_f90/str_comp_mod.f90`): the first `min(la, lb)` characters are compared,
  and only if that is at least the minimum. `la` is the typed length, `lb` the call's length
  argument, usually the name's length. Case-sensitive. A typed word longer than the name
  still matches (`fourierx` runs `fourier`).
- **Word length**: the dispatcher keeps only part of the word. KUPLOT and MIXSCAT keep 4
  characters, DIFFEV 9, and most DISCUS menus 8–14 (`word_length` per menu). In KUPLOT,
  `mergeg` therefore runs `merge`.
- **Order**: the first matching branch in source order wins, so an earlier, shorter branch
  can shadow a later one. In the KUPLOT fit menu, `con` runs `continue`, and `co` is the
  shortest form of `convergence`. `Catalogue.shortest_abbreviation` gives the real minimum.
- **Fall-through**: an unknown word goes to the routines in `inherits`. A few menus (`exp2pdf`,
  `period`, `super`, `symm`, KUPLOT `adt`) ask `discus_kdo_common` / `kuplot_kdo_common`
  first (`inherits_first`), so that routine's `exit` and `show` take precedence there.
- **Menus**: `exit` returns to the parent menu. `read` handles one command and returns,
  which is why `read` / `cell x.cll` needs no `exit`. In a macro, every other menu reads until
  `exit`.
- **Conditional entry** (`enters_when`, verified in the source):
  - `save <file>` writes at once, and only a bare `save` opens the save menu.
  - `insert` opens its menu only for `domain`, `molecule` or `object`.
  - In `output`, `run` opens the transformation menu only with `trans:yes`.
  - `save` in the `pdf` and `rmc` menus never opens a menu.
- **Not executed**: the walker reads each line once. It does not run `do`/`if` blocks or
  `@macro` calls, evaluate variables, or follow `branch`.

## Verification levels

The plan for the suite interface gives every typed wrapper one of three levels:

| Level | Meaning |
|---|---|
| generated | The macro the wrapper writes passes `Catalogue.check` and its own unit tests. |
| engine-ran | A live run on DISCUS 6.23.01 (owner's machine) finished without engine errors and produced the expected files. |
| verified | A pre-registered check against an independent reference passed. Used only for scattering intensities, PDFs and refinement results, with the smallest model that tests the claim (2–3 structures for a deterministic check). |

## Stock engine or fork

Commands run on stock DISCUS 6.23.01 unless they touch something the owner's fork
(`jmsweng/DiffuseCode` `master` `a85c9de`, [FORK_CHANGES.md](https://github.com/jmsweng/DiffuseCode/blob/master/FORK_CHANGES.md))
fixed. Those run on the fork:

| Area | Catalogue commands | Fork change |
|---|---|---|
| REFINE derivatives | `refine` `newparam` (`points:5` under MPI, `points:2`), `run` | MPI step sizes; `points:2` is new |
| REFINE convergence | `refine` `run` | `MAX(dP/sig)` from parameter sigmas |
| Seed restore with Gaussian noise | `common` `seed`, used by REFINE | Box–Muller cache reset |
| RMC with lots | `discus/rmc` `set` (lots), `run` | lot selection, updates, averages |
| REFINE under MPI | `refine` `run` with KUPLOT data | `kuplot_exit` saves the whole data set; HDF5 `hdferr` initialised |
| NUFFT | `discus/fourier` `set technique:nufft`, `run` | coordinate scaling on non-integer grids |
| Debye–Waller cache | `discus/fourier` `temp use`, `run` | table size and cache key |

The fork's speed-ups (ZGEMM kernel, serial kernel, 3D-PDF FFT, REFINE run counts) do not
change results, so they do not decide the engine.

## Quirks of 6.23.01 found while building the catalogue

These were found by reading the source. Except where noted, they have not been confirmed
on a running engine. They are candidates for reports to the DiffuseCode maintainers.

- **`pdf` menu `select`:** the branch literal is `'seleect'` with length 6 (`pdf.f90:255`), so
  a typed `select` is an unknown command and only `seleec`… reaches the branch.
  `mselect` has a similar typo (`'mselelect'`, `pdf.f90:279`).
- **`shear` `mselect`:** the inner test uses `'mseliect'` (`shear.f90:481`), so the molecule
  mode is never set for `mselect`.
- **Typo spellings that are the only accepted ones:** `chemstry` in the `tran` menu
  (`transfrm.f90:129`; `chemistry` is unknown there) and `accoustic` in `wave`.
- **KUPLOT `mergeg`:** unreachable, because KUPLOT keeps 4 characters of a command.
- **No error for unknown commands:** the `3dpdf` and `plot` menus appear to ignore them
  (no `-8` error branch; in `plot` it is commented out, `discus_plot.f90:793`). A generated
  macro must not rely on the engine to report a misspelt command there; `check` does.
- **KUPLOT `adt` `show <args>`:** swallowed by the common `show`, which prints nothing.
- **Dead code:** `kuplot_3dm_menu` has no caller and is not built. DIFFEV `refine` and
  `restrial` only print an "obsolete" error.
- **Interactive-only commands:** `help` works only interactively (`?` always works). KUPLOT
  `mouse` is handled by its main loop.

## Files

- `src/discuspdf/suite_catalogue.py` — loader, `str_comp`, `resolve`, `walk`, `check`.
- `src/discuspdf/data/suite_catalogue/discus_suite_6.23.01.json` — the catalogue (schema
  `discuspdf.suite_catalogue` v1).
- `tools/suite_catalogue/` — `fortran_scan.py`, `help_parse.py`, `build_catalogue.py`,
  `render_samples.py` (macros from the package's generators), `coverage.py`.
- `tests/test_suite_catalogue.py` — matching rules, menu tracking, conditional entry, and
  that every committed macro and every generated macro uses only known commands. With
  `DIFFUSECODE_SRC=<checkout of v.6.23.01>` it also rebuilds the catalogue and compares it
  with the bundled file.
