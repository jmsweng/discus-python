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
            "atom line has the extended .stru columns; this looks like a saved "
            "structure; use read_stru (DISCUS `read stru`), not read cell",
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
