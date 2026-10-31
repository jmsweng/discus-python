# DISCUS 6.23.01 manual summary — `.cll` / `.stru` interchange

Source: `/usr/local/share/suite.hlp` (extracted to `help_grep.txt` / `help_sections.txt`).
This records **what the manual says**, separately from what the Q1–Q7 probes show.

## File layout (keyword-controlled)

Both the unit-cell (`.cll`) and saved-structure (`.stru`) files use the same
keyword-controlled style. Except for the atom lines, each line starts with a
keyword in the first few characters, followed by comma-separated parameters.

- The **first line MUST contain `title`**; the sequence of the other keywords is
  not relevant.
- The **last keyword MUST be `atoms`**, which starts the list of atoms.
- Comments and empty lines may appear anywhere between `title` and `atoms`.
- The keywords **`title`, `spcgr`, `cell` and `atoms` are required**; all others
  are optional.
- The first character after a keyword must be a blank or TAB. Parameters are
  comma-separated; whitespace between parameters may be blanks or TABs.

## `spcgr`

`spcgr {<symbol> | no.} [, "2"] [, setting:{abc|bac|cab|cba|bca|acb}]`

Hermann-Mauguin symbol (or number). Used to **generate all atoms within the unit
cell** from the listed asymmetric-unit atoms. Origin choice 2 via the `"2"`
parameter; non-standard orthorhombic settings via `setting:`. The symbol is
checked against the lattice constants; on contradiction the unit cell is not read.

## `cell`

`cell <a>,<b>,<c>,<alpha>,<beta>,<gamma>` (comma form encouraged; free-form
`a b c alpha beta gamma` also accepted). All six constants required.

## `ncells` / `ncell`

`ncells  nx, ny, nz, ncatoms` — tells DISCUS how large the crystal is in unit
cells along x, y, z. `ncatoms` = number of atoms in **each unit cell** (also
obtainable via `eval n[3]`). Useful when continuing work on a crystal saved by
a previous simulation. When saving in keyword-controlled format, DISCUS can write
the number of unit cells and atoms-per-unit-cell for you.

## `atoms` — allowed columns

```
<name>  <x>, <y>, <z>, <B> [, <propertyflag>] [, <MoleNo>, <MoleAt>, <occ>]
<name>  <x>  <y>  <z>  <B>
<name>  <x>  <y>  <z>
```

- **Mandatory:** atom name and fractional coordinates x, y, z.
- **B** (isotropic thermal coefficient) is the next column.
- Optional trailing columns: property flag, MoleNo, MoleAt, occ (occupancy).
- The optional parameters should be listed on the `atoms` keyword line (header).
- Coordinates and B may be free format (no comma) — in which case the property
  flag must be omitted — or arithmetic expressions.
- For a unit-cell file, **only the asymmetric-unit atoms** should be listed.
- Atom name is left-bound in the first columns; charged atoms e.g. `NA1+`.
- Empty lines are ignored.

## `save` (writing a `.stru`)

- `save <name>` — saves the current structure to `<name>` (probes Q1/Q3).
  Format is identical to the unit-cell/structure file. Each atom is written in
  format `A4,3(2X,F11.6),5X,F6.4`.
  - **Do not use** `save stru,<name>`: it is accepted but writes **no file**.
  - **Do not use** `save "stru",<name>`: the quoted keyword is treated as the
    file name, so a file literally named `stru` is written.
  - The PDF/RMC context also lists `save {"pdf" | "stru" | "mark"}, <file>`;
    for saving the structure from the `discus>` menu, use bare `save <name>`.

## `read`

- `read cell <file>,n1,n2,n3` — reads a unit cell and expands it by the space
  group, repeated n1×n2×n3.
- `read stru <file>` — reads a saved structure (whole crystal, atom by atom; no
  expansion).

## Tracked-fixture observations (grounding, not manual text)

- `examples/gaas.cll`: 4-line header (`title`, `spcgr F-43m`, `cell`, `atoms`),
  2 atom lines with 4 columns (`name x y z B`), no `ncell`, no `format`.
- `gate2_block.stru`: header has `title`, `spcgr P1`, `cell`, `ncell 16,16,16,1`,
  `format numbers,XYZBPMMOS`, `atoms x,y,z,Biso,Property,MoleNo,MoleAt,Occ,St,Sh,Sk,Sl`.
  Atom coordinates run 0…15 (per-unit-cell integers), 2 atoms/cell.
- `C1_frozen.stru`: `title`, `spcgr F -4 3 m , setting:abc`, `cell`,
  `format numbers,XYZBPMMOS`, `atoms …` (no `ncell`). Coordinates are centred
  (e.g. −4.0 … −3.0), B = 1.0, occ = 1.0.
