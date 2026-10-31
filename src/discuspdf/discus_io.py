"""DISCUS 6.23.01 ``.cll`` / ``.stru`` interchange (Package A, Feature A3).

This module reads and writes the two DISCUS structure file kinds and tells them
apart:

* ``.cll`` — a **unit-cell source** (DISCUS ``read cell file,n1,n2,n3``).  The
  listed atoms are the symmetry-independent (asymmetric-unit) atoms; DISCUS
  expands them by the space group and repeats them n1×n2×n3.
* ``.stru`` — a **saved structure** (DISCUS ``read stru file``): the whole crystal,
  atom by atom, with **no** expansion.

Mixing the two up silently multiplies the crystal (the PR #8 failure mode), so
:func:`parse_cll` rejects any ``.stru``-only marker (``ncell``, ``format``, or an
atom line with the extended columns) with a message that says to use
``read_stru`` (DISCUS ``read stru``) instead.

All behaviour is checked against files produced by the real DISCUS 6.23.01; see
``docs/DISCUS_FORMATS.md`` for the verified facts.  This module is **offline**:
it never starts DISCUS.
"""
from __future__ import annotations

import os
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .lattice import Lattice
from .structure import Structure

# DISCUS's default isotropic B (Å²) when a ``.cll`` atom line has no B column
# (verified by probe Q1/Q10: a missing B silently becomes 1.0 on save).
_DISCUS_B_DEFAULT = 1.0

# Number of decimal places the writer uses (DISCUS's ``numbers`` format is
# ``F11.6``; probe Q12 shows DISCUS accepts only the fixed 6-decimal form).
_DECIMALS = 6

# Header keywords a ``.cll`` may contain.  Anything else is an unknown directive.
_CLL_KEYWORDS = frozenset({"title", "spcgr", "cell", "atoms"})

# ``.stru``-only markers: their presence means the file is a saved structure, not
# a unit-cell source.
_STRU_MARKERS = frozenset({"ncell", "format"})


class DiscusFormatError(ValueError):
    """A DISCUS file could not be parsed.

    The message always names the file (or ``<string>``), the line number, and the
    offending text.
    """


@dataclass(frozen=True)
class DiscusCell:
    """A parsed ``.cll`` unit-cell source.

    ``structure`` holds the atoms **as listed** (the asymmetric unit) — it is
    *not* expanded by the space group.  For ``spacegroup != "P1"`` the structure
    is therefore not the whole unit cell.  ``defaults_applied`` records which
    fields were filled in with DISCUS defaults (e.g. ``("b_iso",)`` when the file
    had no B column).
    """

    title: str
    spacegroup: str
    structure: Structure
    defaults_applied: tuple[str, ...]


def _fail(source: str, line_no: int, text: str, msg: str) -> None:
    raise DiscusFormatError(f"{source}:{line_no}: {msg}: {text!r}")


def detect_discus_kind(text: str, *, source: str = "<string>") -> str:
    """Classify DISCUS file text as ``"cll"`` or ``"stru"``.

    Rule (facts decision #5): a file is a ``.cll`` (unit-cell source) iff it has
    **no** ``format`` line **and** **no** ``ncell`` line **and** its atom lines
    carry at most 4 numeric columns (name,x,y,z[,B]).  A ``.stru`` (saved
    structure) has a ``format numbers,…`` line (and usually ``ncell``) with
    11-column atom lines.  Raises ``DiscusFormatError`` if neither.
    """
    has_format = False
    has_ncell = False
    max_atom_cols = 0
    in_atoms = False
    for raw in text.splitlines():
        line = raw.rstrip("\r").strip()
        if not line or line.startswith("#") or line.startswith("!"):
            continue
        keyword, rest = _split_directive(line)
        if not keyword:
            continue
        if keyword == "format":
            has_format = True
        elif keyword == "ncell":
            has_ncell = True
        elif keyword == "atoms":
            in_atoms = True
        elif in_atoms:
            n_values = len(line.replace(",", " ").split()) - 1
            max_atom_cols = max(max_atom_cols, n_values)
    if has_format or has_ncell:
        return "stru"
    if max_atom_cols == 0:
        # No atom lines: cannot decide by content; treat as .cll (parse_cll
        # will reject the missing atom lines).
        return "cll"
    if max_atom_cols <= 4:
        # A .cll has at most 4 numeric columns (name,x,y,z[,B]).
        return "cll"
    if max_atom_cols == 12:
        # A .stru has 12 fields after the name.
        return "stru"
    _fail(
        source, 0, "",
        f"cannot determine DISCUS file kind (atom lines have {max_atom_cols} "
        "fields; expected <=4 for .cll or 12 for .stru)",
    )


def _split_directive(line: str) -> tuple[str, str]:
    """Split a line into (keyword, rest).  The keyword is the first token; the
    first character after it must be a blank or TAB (DISCUS rule)."""
    stripped = line.lstrip()
    if not stripped:
        return "", ""
    # keyword = leading run of non-space, non-tab chars
    i = 0
    while i < len(stripped) and not stripped[i].isspace():
        i += 1
    keyword = stripped[:i].lower()
    rest = stripped[i:]
    return keyword, rest


def _parse_cell_line(rest: str, source: str, line_no: int) -> Lattice:
    """Parse a ``cell`` line into a Lattice (canonical orientation).

    Accepts the comma form (encouraged) and the free-form ``a b c alpha beta
    gamma``.  Requires exactly six finite numbers.
    """
    # Remove any comma, then split on whitespace.
    tokens = rest.replace(",", " ").split()
    if len(tokens) != 6:
        _fail(source, line_no, rest, f"cell line must have exactly six numbers, got {len(tokens)}")
    try:
        values = [float(t) for t in tokens]
    except ValueError:
        _fail(source, line_no, rest, f"cell line has a non-numeric value: {tokens!r}")
    if not all(np.isfinite(v) for v in values):
        _fail(source, line_no, rest, f"cell line has a non-finite value: {tokens!r}")
    a, b, c, alpha, beta, gamma = values
    try:
        return Lattice.from_parameters(a, b, c, alpha, beta, gamma)
    except ValueError as e:
        _fail(source, line_no, rest, f"invalid cell parameters: {e}")


def _parse_atom_line(
    line: str, source: str, line_no: int
) -> tuple[str, list[float]]:
    """Parse an atom line into (species, [x, y, z] [, B] [, occ]).

    A ``.cll`` atom line has 3, 4, or 5 numeric fields after the name.  More than
    5 is a ``.stru``-only marker (the extended columns) and is rejected.  The
    species is the first token of the line (left-bound).
    """
    tokens = line.replace(",", " ").split()
    n_values = len(tokens) - 1  # numeric fields after the species name
    if n_values < 3:
        _fail(source, line_no, line, f"atom line needs name + x,y,z, got {n_values} numeric fields")
    if n_values > 4:
        _fail(
            source, line_no, line,
            "atom line has more than 4 numeric columns; a .cll has at most 4 "
            "numeric columns (name,x,y,z[,B])",
        )
    # n_values is 3 (x,y,z) or 4 (x,y,z,B).  A .cll has no occupancy column
    # (probe Q1: the 5-column B+occ form is rejected by DISCUS), so there is no
    # occupancy to parse here.
    species = tokens[0]
    if not species or any(c.isspace() for c in species):
        _fail(source, line_no, line, f"bad species name: {species!r}")
    if len(species) > 4:
        _fail(source, line_no, line, f"species name exceeds 4 characters: {species!r}")
    try:
        values = [float(t) for t in tokens[1:]]
    except ValueError:
        _fail(source, line_no, line, f"atom line has a non-numeric value: {tokens[1:]!r}")
    if not all(np.isfinite(v) for v in values):
        _fail(source, line_no, line, f"atom line has a non-finite value: {tokens[1:]!r}")
    return species, values


def parse_cll(text: str, *, source: str = "<string>") -> DiscusCell:
    """Parse DISCUS ``.cll`` text into a :class:`DiscusCell`.

    The returned ``structure`` holds the atoms **as listed** (the asymmetric
    unit); it is never expanded by the space group.  For ``spacegroup != "P1"``
    the structure is not the whole unit cell.  A ``.cll`` stores lengths and
    angles, so a non-canonical Cartesian orientation is not kept — fractional
    coordinates are kept as read.
    """
    kind = detect_discus_kind(text, source=source)
    if kind == "stru":
        _fail(
            source, 0, "",
            "this looks like a saved structure; use read_stru (DISCUS `read stru`), "
            "not read cell",
        )
    title: str | None = None
    spacegroup: str | None = None
    lattice: Lattice | None = None
    atoms: list[tuple[str, list[float]]] = []
    seen: set[str] = set()
    in_atoms = False

    for line_no, raw in enumerate(text.splitlines(), start=1):
        line = raw.rstrip("\r")
        stripped = line.strip()
        if not stripped:
            continue  # blank lines are ignored
        if stripped.startswith("#") or stripped.startswith("!"):
            continue  # comment lines are ignored
        keyword, rest = _split_directive(line)
        if not keyword:
            continue

        if keyword in _STRU_MARKERS:
            _fail(
                source, line_no, stripped,
                f"'{keyword}' is a .stru-only marker; this looks like a saved "
                "structure; use read_stru (DISCUS `read stru`), not read cell",
            )

        if keyword == "atoms":
            if "atoms" in seen:
                _fail(source, line_no, stripped, "repeated 'atoms' line")
            seen.add("atoms")
            in_atoms = True
            continue

        if not in_atoms:
            if keyword not in _CLL_KEYWORDS:
                _fail(source, line_no, stripped, f"unknown directive: {keyword!r}")
            if keyword in seen:
                _fail(source, line_no, stripped, f"repeated {keyword!r} line")
            seen.add(keyword)
            if keyword == "title":
                title = rest.strip()
                if not title:
                    _fail(source, line_no, stripped, "empty title")
            elif keyword == "spcgr":
                spacegroup = rest.strip()
                if not spacegroup:
                    _fail(source, line_no, stripped, "empty spcgr")
            elif keyword == "cell":
                lattice = _parse_cell_line(rest, source, line_no)
            continue

        # in_atoms: every non-blank, non-comment line is an atom line.
        species, values = _parse_atom_line(line, source, line_no)
        atoms.append((species, values))

    # Required keywords.
    if title is None:
        _fail(source, 0, "", "missing required 'title' line")
    if spacegroup is None:
        _fail(source, 0, "", "missing required 'spcgr' line")
    if lattice is None:
        _fail(source, 0, "", "missing required 'cell' line")
    if "atoms" not in seen:
        _fail(source, 0, "", "missing required 'atoms' line")
    if not atoms:
        _fail(source, 0, "", "no atom lines found")

    n = len(atoms)
    species_arr = np.array([a[0] for a in atoms], dtype=object)
    frac = np.array([a[1][:3] for a in atoms], dtype=float)

    defaults: list[str] = []
    b_iso: np.ndarray | None
    if all(len(a[1]) == 4 for a in atoms):
        b_iso = np.array([a[1][3] for a in atoms], dtype=float)
        if np.any(b_iso < 0):
            _fail(source, 0, "", "negative B value")
    elif all(len(a[1]) == 3 for a in atoms):
        b_iso = np.full(n, _DISCUS_B_DEFAULT, dtype=float)
        defaults.append("b_iso")
    else:
        _fail(source, 0, "", "mixed atom-line column counts (some with B, some without)")

    # A .cll has no occupancy column (probe Q1: the 5-column B+occ form is
    # rejected), so occupancy is always the Structure default (1.0) and is not
    # a default-applied field.
    occupancy = None

    structure = Structure(
        lattice,
        species_arr,
        frac,
        model="average",
        occupancy=occupancy,
        b_iso=b_iso,
        supercell=(1, 1, 1),
        pbc=(True, True, True),
        min_separation=None,
    )
    return DiscusCell(
        title=title,
        spacegroup=spacegroup,
        structure=structure,
        defaults_applied=tuple(defaults),
    )


def read_cll(path) -> DiscusCell:
    """Read a DISCUS ``.cll`` file from ``path``."""
    p = Path(path)
    text = p.read_text(encoding="ascii")
    return parse_cll(text, source=str(p))


def format_cll(obj, *, title: str | None = None, spacegroup: str | None = None) -> str:
    """Format a :class:`DiscusCell` or :class:`Structure` as ``.cll`` text.

    A ``Structure`` must be a unit cell (``supercell == (1, 1, 1)``, all ``pbc``
    True); a supercell raises ``ValueError`` (use ``write_stru``).  The
    spacegroup defaults to ``"P1"`` for a bare ``Structure`` (the atoms are
    written explicitly, so no symmetry is applied) and to the cell's own
    spacegroup for a ``DiscusCell``.  B is always written explicitly (``0.0``
    when ``b_iso`` is ``None``).  Lines end with ``\\n``.
    """
    if isinstance(obj, DiscusCell):
        structure = obj.structure
        if title is None:
            title = obj.title
        if spacegroup is None:
            spacegroup = obj.spacegroup
    else:
        structure = obj
        if title is None:
            title = "structure"
        if spacegroup is None:
            spacegroup = "P1"

    if tuple(structure.supercell) != (1, 1, 1):
        raise ValueError("a .cll is a unit cell; use write_stru for supercells")
    if not all(structure.pbc):
        raise ValueError("a .cll requires all pbc True")

    a, b, c, alpha, beta, gamma = structure.lattice.parameters
    d = _DECIMALS
    lines = [
        f"title {title}",
        f"spcgr {spacegroup}",
        f"cell {a:.{d}f} {b:.{d}f} {c:.{d}f} "
        f"{alpha:.{d}f} {beta:.{d}f} {gamma:.{d}f}",
        "atoms",
    ]
    b_iso = structure.b_iso
    for i in range(structure.n_atoms):
        x, y, z = (f"{v:.{d}f}" for v in structure.frac[i])
        b_val = 0.0 if b_iso is None else float(b_iso[i])
        lines.append(f"{structure.species[i]:<4} {x} {y} {z} {b_val:.{d}f}")
    return "\n".join(lines) + "\n"


def write_cll(
    path,
    obj,
    *,
    title: str | None = None,
    spacegroup: str | None = None,
    overwrite: bool = False,
) -> Path:
    """Write ``obj`` to ``path`` as a ``.cll`` file, atomically.

    Refuses an existing path unless ``overwrite=True``.  Writes to a temporary
    file in the same directory, then ``os.replace``.  The written text is parsed
    back and must agree with the input (lattice ``isclose`` rtol=1e-8, frac to
    the written precision, exact species, b_iso); on any failure the temporary
    file is deleted, the target is left untouched, and an error is raised.
    """
    target = Path(path)
    if target.exists() and not overwrite:
        raise FileExistsError(f"refusing to overwrite existing file: {target}")

    text = format_cll(obj, title=title, spacegroup=spacegroup)

    # ASCII check (species/titles that are not ASCII -> error).
    try:
        text.encode("ascii")
    except UnicodeEncodeError as e:
        raise ValueError(f"non-ASCII text in .cll output: {e}") from e

    directory = target.parent if str(target.parent) else Path(".")
    fd, tmp_name = tempfile.mkstemp(dir=str(directory), prefix=".cll_", suffix=".tmp")
    os.close(fd)
    tmp = Path(tmp_name)
    try:
        with open(tmp, "w", newline="\n", encoding="ascii") as f:
            f.write(text)
        # Read back and verify agreement with the input.
        written = parse_cll(tmp.read_text(encoding="ascii"), source=str(tmp))
        source = obj if isinstance(obj, DiscusCell) else obj
        if not _cll_roundtrip_agrees(source, written, text):
            raise ValueError("written .cll does not round-trip to the input")
        os.replace(tmp, target)
    except Exception:
        if tmp.exists():
            os.unlink(tmp)
        raise
    return target


@dataclass(frozen=True)
class DiscusStru:
    """A parsed DISCUS ``.stru`` saved structure.

    ``structure`` is the whole crystal (``model="configuration"``): ``lattice``
    is the whole periodic frame and ``frac`` are fractions of it.  ``spacegroup``
    is verbatim and informational — ``read stru`` does **not** expand by it.
    ``ncell`` is the file's ``ncell n1,n2,n3,m`` line (``None`` if absent).
    ``atom_id`` is ``0 … N−1`` in file order (atom ids are *not* stored in the
    file).  ``dropped`` records non-default columns dropped under
    ``allow_loss=True``.
    """

    title: str
    spacegroup: str
    ncell: tuple[int, int, int, int] | None
    structure: Structure
    dropped: tuple[str, ...]
    defaults_applied: tuple[str, ...]


def _parse_ncell(rest: str, source: str, line_no: int) -> tuple[int, int, int, int]:
    tokens = rest.replace(",", " ").split()
    if len(tokens) != 4:
        _fail(source, line_no, rest, f"ncell line must have four integers, got {len(tokens)}")
    try:
        values = [int(t) for t in tokens]
    except ValueError:
        _fail(source, line_no, rest, f"ncell line has a non-integer value: {tokens!r}")
    if any(v <= 0 for v in values):
        _fail(source, line_no, rest, f"ncell line has a non-positive value: {values!r}")
    return tuple(values)  # type: ignore[return-value]


def _parse_stru_atoms(
    text: str, source: str
) -> tuple[list[tuple[str, list[float], dict[str, str]]], int]:
    """Parse the atom lines of a .stru.

    Returns (atoms, max_cols) where each atom is (species, [x,y,z,Biso,Occ,
    Property, MoleNo, MoleAt, St, Sh, Sk], extras).  The .stru atom line has 11
    numeric fields after the name.
    """
    atoms: list[tuple[str, list[float], dict[str, str]]] = []
    in_atoms = False
    max_cols = 0
    for line_no, raw in enumerate(text.splitlines(), start=1):
        line = raw.rstrip("\r")
        stripped = line.strip()
        if not stripped or stripped.startswith("#") or stripped.startswith("!"):
            continue
        keyword, rest = _split_directive(line)
        if not keyword:
            continue
        if keyword == "atoms":
            in_atoms = True
            continue
        if not in_atoms:
            continue
        tokens = line.replace(",", " ").split()
        if len(tokens) < 4:
            _fail(source, line_no, line, f"atom line needs name + x,y,z, got {len(tokens)} fields")
        n_values = len(tokens) - 1
        max_cols = max(max_cols, n_values)
        if n_values != 12:
            _fail(
                source, line_no, line,
                f".stru atom line must have 12 fields after the name, got {n_values}",
            )
        species = tokens[0]
        if not species or any(c.isspace() for c in species):
            _fail(source, line_no, line, f"bad species name: {species!r}")
        if len(species) > 4:
            _fail(source, line_no, line, f"species name exceeds 4 characters: {species!r}")
        try:
            values = [float(t) for t in tokens[1:5]]  # x,y,z,Biso
        except ValueError:
            _fail(source, line_no, line, f"atom line has a non-numeric value: {tokens[1:5]!r}")
        if not all(np.isfinite(v) for v in values):
            _fail(source, line_no, line, f"atom line has a non-finite value: {tokens[1:5]!r}")
        # Property, MoleNo, MoleAt, Occ, St, Sh, Sk are tokens[5:12]
        extras = {
            "Property": tokens[5],
            "MoleNo": tokens[6],
            "MoleAt": tokens[7],
            "Occ": tokens[8],
            "St": tokens[9],
            "Sh": tokens[10],
            "Sk": tokens[11],
        }
        atoms.append((species, values, extras))
    return atoms, max_cols


def parse_stru(
    text: str,
    *,
    source: str = "<string>",
    supercell: tuple[int, int, int] | None = None,
    allow_loss: bool = False,
) -> DiscusStru:
    """Parse DISCUS ``.stru`` text into a :class:`DiscusStru`.

    The returned ``structure`` is the whole crystal (``model="configuration"``):
    ``lattice`` is the whole periodic frame and ``frac`` are fractions of it.
    Atom coordinates in the file are measured in unit cells; they are divided by
    ``ncell`` (no origin shift, no wrapping).  See the module docstring for the
    conversion.
    """
    kind = detect_discus_kind(text, source=source)
    if kind == "cll":
        _fail(
            source, 0, "",
            "this looks like a unit-cell source; use read_cll (DISCUS `read cell`), "
            "not read stru",
        )

    title: str | None = None
    spacegroup: str | None = None
    lattice_params: tuple[float, ...] | None = None
    ncell: tuple[int, int, int, int] | None = None
    seen: set[str] = set()

    # First pass: header lines.
    in_atoms = False
    for line_no, raw in enumerate(text.splitlines(), start=1):
        line = raw.rstrip("\r")
        stripped = line.strip()
        if not stripped or stripped.startswith("#") or stripped.startswith("!"):
            continue
        keyword, rest = _split_directive(line)
        if not keyword:
            continue
        if keyword == "atoms":
            if "atoms" in seen:
                _fail(source, line_no, stripped, "repeated 'atoms' line")
            seen.add("atoms")
            in_atoms = True
            continue
        if in_atoms:
            continue  # atom lines handled by _parse_stru_atoms
        if keyword in ("title", "spcgr", "cell", "ncell", "format"):
            if keyword in seen:
                _fail(source, line_no, stripped, f"repeated {keyword!r} line")
            seen.add(keyword)
            if keyword == "title":
                title = rest.strip()
                if not title:
                    _fail(source, line_no, stripped, "empty title")
            elif keyword == "spcgr":
                spacegroup = rest.strip()
            elif keyword == "cell":
                tokens = rest.replace(",", " ").split()
                if len(tokens) != 6:
                    _fail(source, line_no, rest, f"cell line must have six numbers, got {len(tokens)}")
                try:
                    lattice_params = tuple(float(t) for t in tokens)
                except ValueError:
                    _fail(source, line_no, rest, f"cell line has a non-numeric value: {tokens!r}")
                if not all(np.isfinite(v) for v in lattice_params):
                    _fail(source, line_no, rest, f"cell line has a non-finite value: {tokens!r}")
            elif keyword == "ncell":
                ncell = _parse_ncell(rest, source, line_no)
            elif keyword == "format":
                # accepted; content not validated further
                pass
        else:
            _fail(source, line_no, stripped, f"unknown directive: {keyword!r}")

    if title is None:
        _fail(source, 0, "", "missing required 'title' line")
    if spacegroup is None:
        _fail(source, 0, "", "missing required 'spcgr' line")
    if lattice_params is None:
        _fail(source, 0, "", "missing required 'cell' line")
    if "atoms" not in seen:
        _fail(source, 0, "", "missing required 'atoms' line")

    atoms, _ = _parse_stru_atoms(text, source)
    if not atoms:
        _fail(source, 0, "", "no atom lines found")

    a, b, c, alpha, beta, gamma = lattice_params
    N = len(atoms)

    # Determine n and m.
    if ncell is not None:
        n1, n2, n3, m = ncell
        if n1 * n2 * n3 * m != N:
            _fail(
                source, 0, "",
                f"atom count {N} != n1*n2*n3*m = {n1 * n2 * n3 * m}",
            )
        n = (n1, n2, n3)
    else:
        if supercell is not None:
            n = supercell
            prod_n = n[0] * n[1] * n[2]
            if N % prod_n != 0:
                _fail(
                    source, 0, "",
                    f"atom count {N} not divisible by supercell volume {prod_n}",
                )
            m = N // prod_n
        else:
            # No ncell, no supercell: treat as a single unit cell (decision #4).
            n = (1, 1, 1)
            m = N

    # If both ncell and supercell= are given and differ -> error.
    if ncell is not None and supercell is not None and tuple(ncell[:3]) != tuple(supercell):
        _fail(
            source, 0, "",
            f"ncell {ncell[:3]} and supercell= {supercell} disagree",
        )

    # Build the whole-frame lattice and fractional coordinates.
    unit_lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)
    full_lattice = Lattice.from_parameters(a * n[0], b * n[1], c * n[2], alpha, beta, gamma)

    species_arr = np.array([at[0] for at in atoms], dtype=object)
    b_iso = np.array([at[1][3] for at in atoms], dtype=float)
    if np.any(b_iso < 0):
        _fail(source, 0, "", "negative B value")
    # Occ is the 8th value (index 8 in tokens, index 4 in the values list is
    # Biso; the Occ column is tokens[8]).
    occ = np.array([float(at[2]["Occ"]) for at in atoms], dtype=float)
    if np.any(occ != 1.0):
        _fail(source, 0, "", "occupancy must be exactly 1 in a .stru")

    # frac = x_file / n (element-wise, no shift, no wrap).
    frac = np.empty((N, 3), dtype=float)
    for i, at in enumerate(atoms):
        frac[i, 0] = at[1][0] / n[0]
        frac[i, 1] = at[1][1] / n[1]
        frac[i, 2] = at[1][2] / n[2]

    # cell_index: row r -> cell number r // m, slot r % m; unravel with x fastest.
    cell_index = np.zeros((N, 3), dtype=np.int64)
    site_label = np.empty(N, dtype=object)
    for r in range(N):
        cell_num = r // m
        slot = r % m
        iz = cell_num // (n[0] * n[1])
        rem = cell_num % (n[0] * n[1])
        iy = rem // n[0]
        ix = rem % n[0]
        cell_index[r] = (ix, iy, iz)
        site_label[r] = f"site{slot + 1}"

    # Order check: for every slot and axis, ref = median over cells of
    # (x_file - cell_index); every atom must satisfy |x_file - cell_index - ref| <= 0.5.
    if m > 1:
        for slot in range(m):
            rows = list(range(slot, N, m))
            for axis in range(3):
                vals = np.array([atoms[r][1][axis] - cell_index[r, axis] for r in rows], dtype=float)
                ref = float(np.median(vals))
                for r in rows:
                    if abs(atoms[r][1][axis] - cell_index[r, axis] - ref) > 0.5:
                        _fail(
                            source, 0, "",
                            "atom order does not match ncell cell-major order "
                            f"(row {r}, slot {slot}, axis {axis})",
                        )

    # Column policy (decision #7): the trailing Property, MoleNo, MoleAt, Occ,
    # St, Sh, Sk columns are droppable only with allow_loss=True.  DISCUS writes
    # them as constants for a plain crystal (Q14); a non-constant column is an
    # error unless allow_loss (then dropped + recorded).  Biso is never
    # droppable.
    dropped: list[str] = []
    for name in ("Property", "MoleNo", "MoleAt", "St", "Sh", "Sk"):
        values = {at[2][name] for at in atoms}
        if len(values) > 1:
            if not allow_loss:
                _fail(
                    source, 0, "",
                    f"non-constant {name} column {sorted(values)}; "
                    "pass allow_loss=True to drop it",
                )
            if name not in dropped:
                dropped.append(name)

    structure = Structure(
        full_lattice,
        species_arr,
        frac,
        model="configuration",
        b_iso=b_iso,
        cell_index=cell_index,
        supercell=n,
        pbc=(True, True, True),
        min_separation=None,
    )
    return DiscusStru(
        title=title,
        spacegroup=spacegroup,
        ncell=ncell,
        structure=structure,
        dropped=tuple(dropped),
        defaults_applied=(),
    )


def read_stru(
    path,
    *,
    supercell: tuple[int, int, int] | None = None,
    allow_loss: bool = False,
) -> DiscusStru:
    """Read a DISCUS ``.stru`` file from ``path``."""
    p = Path(path)
    text = p.read_text(encoding="ascii")
    return parse_stru(text, source=str(p), supercell=supercell, allow_loss=allow_loss)


def format_stru(obj, *, title: str | None = None, spacegroup: str | None = None) -> str:
    """Format a :class:`DiscusStru` or :class:`Structure` as ``.stru`` text.

    Only ``model="configuration"`` (a realised crystal) is written; an
    ``average`` model raises ``ValueError`` (a seeded sampling step is needed
    first — Package A, A5).  Rows must already be in DISCUS cell-major order
    consistent with ``cell_index`` and ``supercell``; otherwise
    ``DiscusFormatError`` (the writer never reorders).  B is always written
    explicitly (``0.0`` when ``b_iso is None``).  ``ncell = (*supercell, m)``
    with ``m = N / prod(supercell)``.  Lines end with ``\\n``.
    """
    if isinstance(obj, DiscusStru):
        structure = obj.structure
        if title is None:
            title = obj.title
        if spacegroup is None:
            spacegroup = obj.spacegroup
    else:
        structure = obj
        if title is None:
            title = "discuspdf structure"
        if spacegroup is None:
            spacegroup = "P1"

    if structure.model != "configuration":
        raise ValueError(
            "a .stru is a realised crystal; an average model needs a seeded "
            "sampling step first (Package A, A5)"
        )
    if not all(structure.pbc):
        raise ValueError("a .stru requires all pbc True")

    n = tuple(structure.supercell)
    N = structure.n_atoms
    prod_n = n[0] * n[1] * n[2]
    if N % prod_n != 0:
        raise ValueError(f"atom count {N} not divisible by supercell volume {prod_n}")
    m = N // prod_n

    # Rows must be in DISCUS cell-major order: row r -> cell r//m, slot r%m.
    for r in range(N):
        expected_cell = r // m
        ix, iy, iz = (int(v) for v in structure.cell_index[r])
        if iz * (n[0] * n[1]) + iy * n[0] + ix != expected_cell:
            raise DiscusFormatError(
                f"rows are not in DISCUS order; the writer never reorders "
                f"(row {r} is in cell {ix},{iy},{iz}, expected cell {expected_cell})"
            )

    a, b, c, alpha, beta, gamma = structure.lattice.parameters
    unit_a, unit_b, unit_c = a / n[0], b / n[1], c / n[2]
    d = _DECIMALS
    lines = [
        f"title {title}",
        f"spcgr {spacegroup}",
        f"cell    {unit_a:.{d}f},  {unit_b:.{d}f},  {unit_c:.{d}f}, "
        f"{alpha:.{d}f}, {beta:.{d}f}, {gamma:.{d}f}",
        f"ncell       {n[0]},       {n[1]},       {n[2]},       {m}",
        "format numbers,XYZBPMMOS",
        "atoms      x,              y,              z,             Biso,    "
        "Property,  MoleNo,  MoleAt,   Occ,     St,  Sh,  Sk,  Sl",
    ]
    b_iso = structure.b_iso
    for i in range(N):
        x = structure.frac[i, 0] * n[0]
        y = structure.frac[i, 1] * n[1]
        z = structure.frac[i, 2] * n[2]
        b_val = 0.0 if b_iso is None else float(b_iso[i])
        lines.append(
            f"{structure.species[i]:<4} {x:.{d}f}, {y:.{d}f}, {z:.{d}f}, "
            f"{b_val:.{d}f},       1,       0,       0,   1.000000, _,   0,   0,   0"
        )
    return "\n".join(lines) + "\n"


def write_stru(
    path,
    obj,
    *,
    title: str | None = None,
    spacegroup: str | None = None,
    overwrite: bool = False,
) -> Path:
    """Write ``obj`` to ``path`` as a ``.stru`` file, atomically.

    Same safety as :func:`write_cll`: refuse an existing path unless
    ``overwrite=True``; write to a temporary file in the same directory then
    ``os.replace``; parse the written text back and require ``isclose`` agreement
    with ordering preserved; on failure delete the temp file, leave the target
    untouched, and raise.
    """
    target = Path(path)
    if target.exists() and not overwrite:
        raise FileExistsError(f"refusing to overwrite existing file: {target}")

    text = format_stru(obj, title=title, spacegroup=spacegroup)
    try:
        text.encode("ascii")
    except UnicodeEncodeError as e:
        raise ValueError(f"non-ASCII text in .stru output: {e}") from e

    directory = target.parent if str(target.parent) else Path(".")
    fd, tmp_name = tempfile.mkstemp(dir=str(directory), prefix=".stru_", suffix=".tmp")
    os.close(fd)
    tmp = Path(tmp_name)
    try:
        with open(tmp, "w", newline="\n", encoding="ascii") as f:
            f.write(text)
        written = parse_stru(tmp.read_text(encoding="ascii"), source=str(tmp))
        ref = obj if isinstance(obj, DiscusStru) else obj
        if not _stru_roundtrip_agrees(ref, written):
            raise ValueError("written .stru does not round-trip to the input")
        os.replace(tmp, target)
    except Exception:
        if tmp.exists():
            os.unlink(tmp)
        raise
    return target


def _stru_roundtrip_agrees(source, written: DiscusStru) -> bool:
    """Check that a written .stru parses back to the input structure."""
    ref = source if isinstance(source, Structure) else source.structure
    if written.structure.n_atoms != ref.n_atoms:
        return False
    if not np.array_equal(written.structure.species, ref.species):
        return False
    if not written.structure.lattice.isclose(ref.lattice, rtol=1e-8):
        return False
    if not np.allclose(written.structure.frac, ref.frac, atol=0.5 * 10 ** -_DECIMALS):
        return False
    # ordering preserved: cell_index identical
    if not np.array_equal(written.structure.cell_index, ref.cell_index):
        return False
    ref_b = np.zeros(ref.n_atoms) if ref.b_iso is None else ref.b_iso
    if not np.allclose(ref_b, written.structure.b_iso, atol=0.5 * 10 ** -_DECIMALS):
        return False
    return True


def _cll_roundtrip_agrees(source, written: DiscusCell, text: str) -> bool:
    """Check that a written ``.cll`` parses back to the input structure."""
    if isinstance(source, DiscusCell):
        ref = source.structure
        if written.title != source.title:
            return False
        if written.spacegroup != source.spacegroup:
            return False
    else:
        ref = source
    if not ref.lattice.isclose(written.structure.lattice, rtol=1e-8):
        return False
    if ref.n_atoms != written.structure.n_atoms:
        return False
    if not np.array_equal(ref.species, written.structure.species):
        return False
    # frac to the written precision (6 decimals).
    if not np.allclose(ref.frac, written.structure.frac, atol=0.5 * 10 ** -_DECIMALS):
        return False
    # b_iso: the writer always writes it explicitly, so compare to the written
    # value (0.0 when the input had b_iso=None).
    ref_b = np.zeros(ref.n_atoms) if ref.b_iso is None else ref.b_iso
    if not np.allclose(ref_b, written.structure.b_iso, atol=0.5 * 10 ** -_DECIMALS):
        return False
    return True
