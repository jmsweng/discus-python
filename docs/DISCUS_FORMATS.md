# DISCUS 6.23.01 `.cll` / `.stru` format — verified facts

Every fact below is backed by a probe run against the real DISCUS 6.23.01
(`/usr/local/bin/discus_suite`, WSL `Ubuntu-22.04`). Evidence files are in
`tests/fixtures/discus_6_23_01/evidence/`; the macros/inputs are in
`benchmarks/discus_6_23_01/a3_format_probes/`.

## `.cll` (unit-cell source, `read cell file,n1,n2,n3`)

| Fact | Evidence file | Verbatim excerpt |
|---|---|---|
| A `.cll` needs `title`, `spcgr`, `cell`, `atoms`; only `title` is hard-required (see `.stru` Q7). | `Q7_v1.log` | `***APPL*** Structure q7_v1.stru *** -48 ***` (title removed → fail) |
| `read cell` expands the listed atoms by the space group and repeats n1×n2×n3. | `Q2.log` | `cell gaas.cll,2,2,2` → `Value of n[1] = 64.00000000000`; `…,1,1,1` → `8.00000000000` |
| F-43m GaAs (2 listed atoms) expands to 8 atoms/cell. | `Q2.log` | `Value of n[1] = 8.00000000000` for `…,1,1,1` |
| Atom line `name x y z` (no B) is accepted. | `Q1.log` | `cell q1_xyz.cll,1,1,1` → `Value of n[1] = 8.00000000000` |
| Atom line `name x y z B` is accepted. | `Q1.log` | `cell q1_xyzB.cll,1,1,1` → `Value of n[1] = 8.00000000000` |
| Atom line `name x y z B occ` (5 cols) is **rejected**. | `Q1.log` | `***FORT*** Error reading parameters for *** -1 ***` |
| A missing B silently becomes **Biso = 1.0** on save. | `Q1.log` + `q1_xyz.stru` | saved line `GA … 1.000000, 1, 0, 0, 1.000000, _, 0, 0, 0` from a no-B input |

## `.stru` (saved structure, `read stru file`)

| Fact | Evidence file | Verbatim excerpt |
|---|---|---|
| `title` is the only hard-required header line. | `Q7_v1.log` | `***APPL*** Structure q7_v1.stru *** -48 ***` |
| `spcgr`, `cell`, `format` may be absent (defaults used). | `Q7_v2.log`, `Q7_v3.log`, `Q7_v4.log` | each → `Value of n[1] = 8.00000000000` |
| `#` and `!` comment lines and blank lines are tolerated. | `Q7_v7.log`, `Q7_v8.log`, `Q7_v9.log` | each → `Value of n[1] = 8.00000000000` |
| An `ncell` line is **not** required to load. | `Q5.log` | with `ncell`: n[1]=4096; without: n[1]=4096 (identical) |
| `ncell nx,ny,nz,ncatoms` is informational for `read stru`. | `Q5.log` | `gate2` (ncell 16,16,16,1): n[1]=4096, n[2]=2, n[3]=1, n[4]=0 |
| Saved atom format is `A4,3(2X,F11.6),5X,F6.4`; 11 numeric fields after the name. | `gaas_222.stru` | `GA -1.000000, -1.000000, -1.000000, 1.000000, 1, 0, 0, 1.000000, _, 0, 0, 0` |
| `format numbers,XYZBPMMOS` is written by `save`. | `gaas_222.stru` | `format numbers,XYZBPMMOS` |
| `save` writes **no `ncell` line** (it is dropped on save), so a saved `.stru` must be read with `supercell=`. | `gaas_222.stru`, `C1_frozen.stru` | neither file has an `ncell` line |
| 5-column `.stru` (no `format`): col 4 = B, col 5 = ignored. | `Q8.log` + `corr_model_mod_out.stru` | input col4=0.5, col5=0.8 → saved `Biso 0.500000 … Occ 1.000000` |
| A `.stru` atom line with **fewer** columns than the `format` line promises **fails**. | `Q10.log` | `***APPL*** File format specified as numbers only *** -49 ***` / `Error reading atom coordinates` |
| Coordinates must be the fixed 6-decimal `F11.6` form; 10-decimals, scientific, and leading `+` all **fail**. | `Q12.log`, `Q12b.log`, `Q12c.log` | `***APPL*** File format specified as numbers only *** -49 ***` (all three) |
| CRLF line endings and trailing spaces load; **tab separators fail**. | `Q11.log`, `Q11b.log` | CRLF `.stru`/`.cll` → `Value of n[1] = 8.00000000000`; tab → `***APPL*** … -49 ***` |
| Species names are uppercased on save; charges kept; 9-char names fail (truncated to 4). | `Q9.log` + saved files | `ga`→`GA`, `Fe3+`→`FE3+`, `O2-`→`O2-`; `ABCDEFGHI` → `***APPL*** Atom name: ABCD` / `No chemical element name detected` |
| Occupancy < 1 is representable and round-trips. | `Q13.log` + `q13_out.stru` | input Occ 0.5 on 2 atoms → saved keeps `0.500000` on those 2, `1.000000` on the rest |
| For a plain crystal DISCUS writes `Property=1, MoleNo=0, MoleAt=0, Occ=1.000000, St=_, Sh=0, Sk=0`; only Biso varies. | `q14_analyze.py` output | `C1_frozen/gate2/gaas_222`: Biso ∈ {1.0, 0.0}, all other columns single-valued as listed |

## Coordinate convention and atom order

| Fact | Evidence file | Verbatim excerpt |
|---|---|---|
| `read cell …,2,2,2` then `save` yields **centred** coordinates (−1…+0.75 for a 2×2×2 GaAs), not 0…1 supercell fractions. | `gaas_222.stru` | x,y,z min=−1, max=+0.75 (64 atoms) |
| `read cell …,1,1,1` then `save` yields 0…1 coordinates. | `gaas_111.stru` | x,y,z in 0…1 (8 atoms) |
| x varies fastest across cells; 8 atoms per unit cell for F-43m GaAs. | `gaas_222.stru` | first 16 rows: x cycles −1,−1,−0.5,−0.5,−0.75,−0.25,−0.75,−0.25 then 0,0,0.5,0.5,… |
| A hand-written `.stru` may use 0…(n−1) integer-style per-unit-cell coords (does **not** match the centred convention). | `gate2_block.stru` (tracked) | `ncell 16,16,16,1`, coords 0…15 |

## The double-expansion hazard (Q6)

Reading a saved `.stru` with `read cell` re-expands it, silently multiplying the
crystal:

| Command | Atoms |
|---|---|
| `read stru gaas_222.stru` | 64 |
| `read cell gaas_222.stru,1,1,1` | **256** (4×) |
| `read stru C1_frozen.stru` | 4096 |
| `read cell C1_frozen.stru,1,1,1` | **16384** (4×) |

Evidence: `Q6.log`. DISCUS prints no warning. This is the PR #8 failure mode.

## Defaults that bite (B)

A `.cll` atom line with no B column is read with **B = 1.0** (DISCUS default), and
`save` writes `Biso 1.000000`. Evidence: `Q1.log` / `q1_xyz.stru`. PR #8 failed
once because a missing B silently became DISCUS's default. The Python writer must
therefore always write B explicitly (`0.0` when `b_iso` is `None`).

## Decisions for the Python I/O

1. **Writer header for `.stru`** — write `title`, `spcgr <sym>, setting:abc`,
   `cell a,b,c,90,90,90` (comma form), `format numbers,XYZBPMMOS`, and an
   `atoms x,y,z,Biso,Property,MoleNo,MoleAt,Occ,St,Sh,Sk,Sl` header. Include an
   `ncell n1,n2,n3,m` line (m = atoms per unit cell) so the file round-trips with
   `read stru` and carries the supercell info. Rationale: Q3/Q5 show DISCUS writes
   exactly this header and reads it back; `ncell` is harmless and informative.
2. **Writer decimal places** — **6** decimals, `F11.6`-style, no scientific
   notation, no leading `+`. Rationale: Q12 shows DISCUS's `numbers` format accepts
   only the fixed 6-decimal form it writes; anything else fails.
3. **Deriving `cell_index` and `site_label`** — from `ncell n1,n2,n3,m` and atom
   order: `cell_index = i // m`, and within a cell the site index `i % m`; the
   (ix,iy,iz) cell position is the mixed-radix decomposition of `cell_index` in
   (n1,n2,n3) with **x fastest** (Q4: x varies fastest). `site_label` = the species
   name at that site. Rationale: Q4/Q5 establish m = atoms/cell and the x-fastest
   ordering.
4. **Handling a `.stru` without `ncell`** — a `.stru` without an `ncell` line
   (the common case: DISCUS `save` drops the `ncell` line) is read with the
   `supercell=` argument of `read_stru`/`parse_stru`:

   * `supercell=None` (default): read as **one** cell (`supercell=(1,1,1)`,
     `m = N`) only if every axis has `hi − lo < 1 + 1e-6` (lo/hi = min/max of
     the file coordinates per axis).  Otherwise raise `DiscusFormatError`
     naming the per-axis spans: the atoms span more than one unit cell.
   * `supercell=(n1,n2,n3)`: use it (`m = N / prod(n)`); the cell-major order
     check then runs.
   * `supercell="infer"`: per axis `n = round(−2·lo)`; accepted only for a
     centred crystal — `n ≥ 1`, `|−2·lo − n| ≤ 1e-5`, `hi < n/2 + 1e-6`, and
     `N % prod(n) == 0` — i.e. the layout `read cell …,n1,n2,n3` followed by
     `save` produces (Q4: −1…+0.75 for n = 2).  `gaas_222.stru` → (2,2,2);
     `C1_frozen.stru` → (8,8,8); `gaas_111.stru` (0…0.75) reads with the
     default.

   Rationale: Q5/Q7 show `ncell` is optional and its absence does not change
   the atom count, but Q4 shows a saved supercell is **centred** and spans n
   unit cells, so reading it as one cell puts the atoms on top of each other
   (the A3-6a failure).  The default must therefore refuse multi-cell files.
5. **Telling a `.cll` from a `.stru` by content** — a file is a `.clu`/`.cll`
   (unit-cell source) iff it has **no** `format` line **and** **no** `ncell` line
   **and** its atom lines carry **≤ 4** numeric columns (name,x,y,z[,B]); a
   `.stru` has a `format numbers,…` line (and usually `ncell`) with 11-column atom
   lines. Rule true for all fixtures: `gaas.cll` (no format/ncell, 4-col) vs
   `gaas_111/222.stru` (format present, 11-col) and `p1_2x2x2_a4.stru` (format +
   ncell, 11-col). Rationale: Q1/Q3/Q8 show the two kinds differ exactly in these
   markers.
6. **Species name policy** — the Python writer writes species names
   **verbatim** (it does not upper-case them); DISCUS upper-cases them on save
   (`Ga` → `GA`, `Fe3+` → `FE3+`, `O2-` unchanged).  Names longer than 4
   characters are **rejected on read and on write** (DISCUS truncates to 4 and
   then fails to find the element).  An optional trailing charge token (`2-`,
   `3+`) is kept. Rationale: Q9.
7. **Columns droppable only with `allow_loss=True`** — the trailing
   `Property, MoleNo, MoleAt, Occ, St, Sh, Sk` columns (DISCUS writes them as
   constants for a plain crystal; dropping them loses no per-atom data for the
   fixtures). `Biso` is **never** droppable (a missing B silently becomes 1.0).
   Rationale: Q14 (constants) + Q1/Q10 (B default).
