# Package A PR A3 — DISCUS 6.23.01 file I/O

## 1. Base / head / PR / environment

- **Base:** `origin/main` at `b4f7721` (Package A PR A-2: Editable Structure model)
- **Head:** see PR head (A3-6b: correct example and notebook, live physics checks, committed live evidence, report)
- **PR:** [#12](https://github.com/jmsweng/discus-python/pull/12) (draft)
- **Environment:** Windows 10, Python 3.14, DISCUS_SUITE 6.23.01 in WSL `Ubuntu-22.04`

## 2. What is delivered

- `src/discuspdf/discus_io.py`: `.cll` and `.stru` reader/writer for DISCUS_SUITE 6.23.01.
- `tests/test_discus_io_cll.py`: 22 offline tests for `.cll` parsing and writing.
- `tests/test_discus_io_stru.py`: 34 offline tests for `.stru` parsing and writing (15 original + 19 from A3-6a).
- `tests/test_discus_io_live.py`: 7 opt-in live tests (L0–L5) against the real DISCUS engine, with physics (neighbour) checks added in A3-6b.
- `tests/discus_live.py`: helper for running DISCUS macros via WSL.
- `docs/DISCUS_FORMATS.md`: verified facts and 7 design decisions.
- `docs/DISCUS_IO.md`: user documentation.
- `examples/discus_io_basics.py`: offline example.
- `notebooks/05_discus_io.ipynb`: notebook with a guarded live cell.

## 3. Design decisions table

| Decision | Choice | Rationale |
|---|---|---|
| `.cll` atom columns | ≤ 4 (name,x,y,z[,B]) | Q1: 5-column form is rejected by DISCUS |
| `.stru` atom columns | 12 (x,y,z,Biso,Property,MoleNo,MoleAt,Occ,St,Sh,Sk) | Q14: DISCUS writes 12 columns |
| Missing B in `.cll` | Written as explicit `0.0` | Q1: DISCUS defaults missing B to 1.0 |
| `.stru` without `ncell` | `supercell=` selects the supercell: `None` reads one cell only if every axis spans < 1 unit cell (else raise); `(n1,n2,n3)` uses it; `"infer"` recovers a centred crystal (`n = round(−2·lo)`) | Q4/Q5/Q7: a saved supercell is centred and spans n cells, so one-cell is wrong (A3-6a) |
| Species names | Written **verbatim** (no upper-casing); DISCUS upper-cases on save; >4 chars rejected on read and write | Q9: DISCUS uppercases and truncates |
| Decimal places | 6 (F11.6-style) | Q12: DISCUS accepts only fixed 6-decimal form |
| Atomic write | Temp file + `os.replace` | Safety: no partial writes |
| Overwrite | Refused by default | Safety: no accidental data loss |

## 4. Verified facts table

Copied from `docs/DISCUS_FORMATS.md` "Decisions":

1. **Writer header for `.stru`** — write `title`, `spcgr`, `cell`, `format numbers,XYZBPMMOS`, `atoms` header, and `ncell` line.
2. **Writer decimal places** — 6 decimals, F11.6-style, no scientific notation, no leading `+`.
3. **Deriving `cell_index`** — from `ncell n1,n2,n3,m` and atom order: `cell_index = i // m`, x-fastest.
4. **Handling a `.stru` without `ncell`** — `supercell=` selects the supercell: `None` reads one cell only if every axis spans < 1 unit cell (else raise); `(n1,n2,n3)` uses it; `"infer"` recovers a centred crystal (`n = round(−2·lo)`). (A3-6a; overrides the original "treat as one cell".)
5. **Telling a `.cll` from a `.stru`** — `.cll` has no `format`/`ncell` and ≤ 4 numeric columns; `.stru` has `format` and 12 columns.
6. **Species name policy** — written verbatim (no upper-casing); DISCUS upper-cases on save; >4 chars rejected on read and write.
7. **Columns droppable only with `allow_loss=True`** — `Property`, `MoleNo`, `MoleAt`, `Occ`, `St`, `Sh`, `Sk`; `Biso` is never droppable.

## 5. Tests

- **Offline:** 56 tests (22 `.cll` + 34 `.stru`), all passing.
- **Live:** 7 tests (L0–L5), all passing with `DISCUS_LIVE=1` (~5 s); all skipped without the env var.
- **Mutation (A3-3):** 5 breaks caught (ncell count, occupancy, order, ncell line, `/n` → `*n`).
- **Mutation (A3-4):** 3 breaks caught (B default, ncell line, order check).
- **Mutation (A3-6a):** 4 breaks caught — (a) `min_separation=None` back → `test_overlap_check_is_active`; (b) span guard removed → `test_no_ncell_policy`; (c) `infer` → `round(−lo)` → `test_no_ncell_policy`; (d) order check skips `m=1` → new `test_order_check_runs_for_m1`.
- **Mutation (A3-6b):** 3 breaks caught — (a) `format_stru` writes supercell lengths → L2 + `test_roundtrip_gaas_222`; (b) `format_stru` writes `frac` unscaled → L2 + `test_roundtrip_gaas_222`; (c) `infer` off by 2 → L1/L3 + `test_no_ncell_policy`/`test_gaas_222_is_zinc_blende`.

## 6. Evidence

- **Fixtures:** `tests/fixtures/discus_6_23_01/` (gaas.cll, gaas_111.stru, gaas_222.stru, p1_2x2x2_a4.stru).
- **Probe logs:** `tests/fixtures/discus_6_23_01/evidence/` (Q1–Q14, smoke).
- **Live evidence:** `tests/fixtures/discus_6_23_01/live_evidence/20261031T0403Z/` (7 passing-run logs + `MANIFEST.sha256`, 157 419 bytes total). **Superseded in part:** this run predates the A3-6c correction of the live-test crystals (L2 used unit-cell-sized supercell frames and a 32-atom Ga-only "GaAs") and of the B assertion (see below). A new run with the corrected tests must be committed as an additional folder; it is pending on the owner's DISCUS machine.

## 7. Caveats / not done

- Symmetry expansion of `.cll` (the Python API does not expand; DISCUS does).
- Occupancy < 1 in `.stru` (rejected by `parse_stru`).
- Molecule/property columns (dropped only with `allow_loss=True`).
- CIF/ASE formats (next PR).

## 8. Changed files

```
88 files changed, 2917 insertions(+), 4 deletions(-)
```

Key files:
- `src/discuspdf/discus_io.py` (917 lines)
- `src/discuspdf/__init__.py` (+21 lines)
- `pyproject.toml` (+5 lines, `discus_live` marker)
- `tests/test_discus_io_cll.py` (292 lines)
- `tests/test_discus_io_stru.py` (243 lines)
- `tests/test_discus_io_live.py` (412 lines)
- `tests/discus_live.py` (170 lines)
- `docs/DISCUS_FORMATS.md` (105 lines)
- `docs/DISCUS_IO.md` (new)
- `examples/discus_io_basics.py` (new)
- `notebooks/05_discus_io.ipynb` (new)
