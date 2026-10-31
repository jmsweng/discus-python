# Package A PR A3 — DISCUS 6.23.01 file I/O

## 1. Base / head / PR / environment

- **Base:** `origin/main` at `b4f7721` (Package A PR A-2: Editable Structure model)
- **Head:** `cb7ccdf` (A3-4: opt-in live DISCUS engine tests)
- **PR:** [#12](https://github.com/jmsweng/discus-python/pull/12) (draft)
- **Environment:** Windows 10, Python 3.14, DISCUS_SUITE 6.23.01 in WSL `Ubuntu-22.04`

## 2. What is delivered

- `src/discuspdf/discus_io.py`: `.cll` and `.stru` reader/writer for DISCUS_SUITE 6.23.01.
- `tests/test_discus_io_cll.py`: 22 offline tests for `.cll` parsing and writing.
- `tests/test_discus_io_stru.py`: 15 offline tests for `.stru` parsing and writing.
- `tests/test_discus_io_live.py`: 7 opt-in live tests (L0–L5) against the real DISCUS engine.
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
| `.stru` without `ncell` | Treated as a single unit cell | Q5/Q7: `ncell` is optional |
| Species names | 1–4 chars, uppercased on save, charges kept | Q9: DISCUS uppercases and truncates |
| Decimal places | 6 (F11.6-style) | Q12: DISCUS accepts only fixed 6-decimal form |
| Atomic write | Temp file + `os.replace` | Safety: no partial writes |
| Overwrite | Refused by default | Safety: no accidental data loss |

## 4. Verified facts table

Copied from `docs/DISCUS_FORMATS.md` "Decisions":

1. **Writer header for `.stru`** — write `title`, `spcgr`, `cell`, `format numbers,XYZBPMMOS`, `atoms` header, and `ncell` line.
2. **Writer decimal places** — 6 decimals, F11.6-style, no scientific notation, no leading `+`.
3. **Deriving `cell_index`** — from `ncell n1,n2,n3,m` and atom order: `cell_index = i // m`, x-fastest.
4. **Handling a `.stru` without `ncell`** — treat as a single unit cell.
5. **Telling a `.cll` from a `.stru`** — `.cll` has no `format`/`ncell` and ≤ 4 numeric columns; `.stru` has `format` and 12 columns.
6. **Species name policy** — 1–4 chars, uppercase on write, optional trailing charge.
7. **Columns droppable only with `allow_loss=True`** — `Property`, `MoleNo`, `MoleAt`, `Occ`, `St`, `Sh`, `Sk`; `Biso` is never droppable.

## 5. Tests

- **Offline:** 37 tests (22 `.cll` + 15 `.stru`), all passing.
- **Live:** 7 tests (L0–L5), all passing with `DISCUS_LIVE=1` (~5 s).
- **Mutation (A3-3):** 5 breaks caught (ncell count, occupancy, order, ncell line, `/n` → `*n`).
- **Mutation (A3-4):** 3 breaks caught (B default, ncell line, order check).

## 6. Evidence

- **Fixtures:** `tests/fixtures/discus_6_23_01/` (gaas.cll, gaas_111.stru, gaas_222.stru, p1_2x2x2_a4.stru).
- **Probe logs:** `tests/fixtures/discus_6_23_01/evidence/` (Q1–Q14, smoke).
- **Live evidence:** `C:\tmp\a3live\` (macro outputs, saved `.stru` files).

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
