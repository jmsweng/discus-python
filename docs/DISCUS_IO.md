# DISCUS 6.23.01 file I/O

## What is supported

`.cll` unit-cell files and `.stru` saved structures of DISCUS_SUITE 6.23.01,
exactly the variants in [DISCUS_FORMATS.md](DISCUS_FORMATS.md); nothing else.
Other DISCUS versions or format variants are not claimed.

## The two kinds and the hazard

| Kind | Command | Expands? | Example |
|---|---|---|---|
| `.cll` (unit-cell source) | `read cell file,n1,n2,n3` | Yes — space group × n1×n2×n3 | `gaas.cll` (2 atoms) → 8/cell → 64 for 2×2×2 |
| `.stru` (saved structure) | `read stru file` | No — loads as-is | `gaas_222.stru` (64 atoms) → 64 |

**The hazard (Q6):** reading a saved `.stru` with `read cell` re-expands it,
silently multiplying the crystal:

| Command | Atoms |
|---|---|
| `read stru gaas_222.stru` | 64 |
| `read cell gaas_222.stru,1,1,1` | **256** (4×) |

DISCUS prints no warning. The Python API refuses the wrong kind:
`parse_cll` raises `DiscusFormatError` if the file has a `format` line or
`ncell` line (a `.stru` marker); `parse_stru` raises if the file has
≤ 4 numeric columns per atom line (a `.cll` marker).

## Coordinate conversion

Files store coordinates in **unit-cell** form; `Structure` stores
**supercell** fractions.

- **Reading a `.stru` with `ncell n1,n2,n3,m`:**
  `frac = x_file / n` (element-wise, no shift, no wrap).
- **Reading a `.stru` without `ncell`** (the common case: DISCUS `save` drops
  the `ncell` line) — selected by the `supercell=` argument:
  - `supercell=None` (default): read as one cell (`n = (1,1,1)`, `frac = x_file`)
    **only if every axis spans less than one unit cell**; otherwise
    `DiscusFormatError` (pass `supercell=(n1,n2,n3)` or `supercell="infer"`).
  - `supercell=(n1,n2,n3)`: `frac = x_file / n`; the cell-major order check runs.
  - `supercell="infer"`: per axis `n = round(−2·lo)` for a centred crystal
    (`n ≥ 1`, `|−2·lo − n| ≤ 1e-5`, `hi < n/2 + 1e-6`, `N % prod(n) == 0`);
    `frac = x_file / n`.
- **Writing a `.stru`:** `x_file = frac * n` (element-wise).

**Worked example** (cell 4 Å, `ncell 2,2,2,1`):

```
Unit cell: a = b = c = 4 Å
Supercell: 2×2×2, m = 1 (one atom per unit cell)
Atom at supercell frac (0.25, 0.5, 0.0):
  x_file = 0.25 × 2 = 0.5
  y_file = 0.5 × 2 = 1.0
  z_file = 0.0 × 2 = 0.0
Reading back:
  frac = (0.5/2, 1.0/2, 0.0/2) = (0.25, 0.5, 0.0) ✓
```

## API

### Reading

```python
from discuspdf import read_cll, read_stru, parse_cll, parse_stru

# Read a .cll file (unit-cell source).
s = read_cll("examples/gaas.cll")          # 8 atoms (F-43m GaAs)
s = read_cll("examples/gaas.cll", supercell=(2, 2, 2))  # 64 atoms

# Read a .stru file (saved structure).  DISCUS `save` drops the ncell line,
# so a saved multi-cell file needs supercell=:
s = read_stru("tests/fixtures/discus_6_23_01/gaas_222.stru", supercell="infer")  # 64 atoms, 2x2x2
s = read_stru("some_file.stru", supercell=(8, 8, 8))  # 4096 atoms (no ncell)
s = read_stru("tests/fixtures/discus_6_23_01/gaas_111.stru")  # one cell: 0..0.75
```

### Writing

```python
from discuspdf import write_cll, write_stru, format_cll, format_stru

# Write a Structure as .cll (unit-cell source).
write_cll("out.cll", structure)

# Write a Structure as .stru (saved structure).
write_stru("out.stru", structure)

# Get the text without writing to a file.
text = format_cll(structure)
text = format_stru(structure)
```

### Detection

```python
from discuspdf import detect_discus_kind

kind = detect_discus_kind("file.txt")  # "cll" or "stru"
```

## Loss and identity policy

- **Not stored:** atom ids, site labels, lattice orientation, symmetry expansion
  (not performed), `Property`/`MoleNo`/`MoleAt`/`St` columns (only with
  `allow_loss=True`).
- **`b_iso=None`** → explicit `0.0` in the file (DISCUS defaults missing B to 1.0).
- **Species names** are written **verbatim** (the writer does not upper-case
  them); DISCUS upper-cases on save (`Ga` → `GA`).  Names longer than 4
  characters are rejected on read and on write.
- **`defaults_applied`:** recorded when `allow_loss=True` drops non-constant columns.

## Safety

- **Atomic write:** writes to a temporary file then `os.replace`.
- **No overwrite by default:** `FileExistsError` unless `overwrite=True`.
- **Read-back check:** the written text is parsed back and must `isclose` the input.
- **`\n` line endings:** all output uses `\n` (not CRLF).

## Testing

- **Offline tests:** `tests/test_discus_io_cll.py` (22 tests),
  `tests/test_discus_io_stru.py` (15 tests).
- **Live tests:** set `DISCUS_LIVE=1` and run:
  ```powershell
  $env:DISCUS_LIVE="1"; .\.venv\Scripts\python.exe -m pytest -q -m discus_live
  ```
  7 tests (L0–L5) against the real DISCUS 6.23.01 engine via WSL.
- **Committed evidence:** `tests/fixtures/discus_6_23_01/evidence/` (probe logs),
  `tests/fixtures/discus_6_23_01/` (fixture files).

## Not done

- Symmetry expansion of `.cll` (the Python API does not expand; DISCUS does).
- Occupancy < 1 in `.stru` (rejected by `parse_stru`).
- Molecule/property columns (dropped only with `allow_loss=True`).
- CIF/ASE formats (next PR).
