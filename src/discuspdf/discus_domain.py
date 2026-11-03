"""DISCUS ``domain`` menu as a typed step: insert domains of a second structure into the host crystal.

A domain is a region of a *content* structure (another phase or orientation) that replaces part
of the host. DISCUS places each domain at a *pseudo-atom* position: an atom of the input file
whose name is the pseudo-atom name. The input file is a ``.stru``, read into internal storage
with the full structure reader (``read_to_internal``, ``domain.f90`` 733-746), so the name
comes from its atom type. Its lattice must equal the host's (within 1e-4 A and 1e-3 deg) and its
space-group letter must match (``domain.f90`` 1083-1103, errors -87 / -88).

Positions are in the host's frame: after ``read cell ...,n1,n2,n3`` that frame is centred,
running from about -n/2 to n/2 cells. Content atoms are placed at ``O u + position``, so the
content file must be a crystal centred on its own origin that spans the domain shape. The shape, orientation, character and fuzzy distance of each pseudo atom are set with
``assign`` at the domain level. The file-level facts behind every choice here are in
``docs/DISCUS_DOMAIN_FACTS.md`` (with file:line references into DISCUS 6.23.01).

::

    from discuspdf.discus_domain import Domains
    from discuspdf.discus_ops import StructureJob
    from discuspdf.fourier_macro import Seed, StructureSource

    step = Domains(content="dom_cu.stru", inputfile="pseudo.stru", name="DOM1", character="sphere",
                   shape=((2.0, 0.0, 0.0, 0.0), (0.0, 2.0, 0.0, 0.0), (0.0, 0.0, 2.0, 0.0)),
                   fuzzy=2.0)
    job = StructureJob(StructureSource.cell("cu_b05.cll", (6, 6, 6)), [step],
                       seed=Seed((11, 22, 33, 44, 55, 66, 77, 88)))
    result = job.run(runner, {"cu_b05.cll": ..., "dom_cu.stru": ..., "pseudo.stru": ...},
    #                with StructureJob(..., readback="box"): the result is not whole unit cells

Semantics used (all from ``discus/prog/domain.f90`` unless stated):

* ``shape`` row i is ``(S[i][0], S[i][1], S[i][2], t_i)``. A point ``s`` of the unit shape maps to
  host fractional coordinates ``h = S s + (t + position)`` (``micro_read_simple`` 2918-2924,
  ``micro_read_atom`` 2474-2475). The unit shape is a cube of half-width 1 (``|s_i| <= 1``), a
  cylinder along the third axis (``|s_3| <= 1`` and radial metric length < |b|), or a sphere of
  metric radius |b| (2482-2513). The identity S gives a half-width of one unit cell (two cells
  across).
* ``orient`` row i is ``(O[i][0], O[i][1], O[i][2], t_i)``. Each content atom with fractional
  coordinates ``u`` is placed at ``O u + (t + position)`` (``micro_read_simple`` 2920-2924,
  ``micro_read_atom`` 2474). The default in DISCUS is an all-zero matrix (``domain_reset``
  3221), which would collapse every content atom onto the pseudo-atom position, so this step
  always writes the orientation and uses the identity when none is given.
* Content atoms outside the shape are not inserted; atoms inside are appended to the crystal
  (``micro_read_atom`` 2542-2556, 2587).
* Host atoms are removed by distance, not by shape: a host atom becomes a vacancy (type 0, kept
  in the list) when it lies within ``fuzzy`` of an inserted atom (``micro_fuzzy_rem`` 3014,
  3075). The removal runs only when ``fuzzy > 0`` (``micro_read_atom`` 2776). DISCUS's default
  fuzzy distance is 2.55 A (``domain_reset`` 3220), so this step always writes it.
* ``run`` draws random numbers (``gasdev`` three times on every pseudo-atom line, ``domain.f90``
  2916), even though with no size spread the result does not depend on them. The job still asks
  for a seed because these draws change the random state that later steps see.

Left out on purpose (not pinned down for this step): mode ``domain`` (explicit domain lists,
``CLU_IN_CLUSTER``), mode ``irregular`` (``CLU_IN_IRREG``, ``domain_irregular``), the size
distribution ``sigma`` of ``assign shape``, the ``set`` menu (removal modes), and ``symm``.
"""

from __future__ import annotations

import math
import numbers
import re
from dataclasses import dataclass
from typing import Sequence

from .discus_ops import Step
from .macro import Macro

CHARACTERS = ("cube", "cylinder", "sphere", "fuzzy")  # domain.f90 389-399; irregular is out of scope
MODES = ("pseudoatoms",)  # domain.f90 227-236; 'domain' and 'irregular' are out of scope
IDENTITY_ORIENT = ((1.0, 0.0, 0.0, 0.0), (0.0, 1.0, 0.0, 0.0), (0.0, 0.0, 1.0, 0.0))

_NAME = re.compile(r"[A-Za-z0-9_]{1,4}")  # clu_name is CHARACTER(LEN=4) (domain_mod.f90:28)
_FORBIDDEN = frozenset(",\"'#!&/\\")
_MAX_FILE = 200  # clu_content is CHARACTER(LEN=200) (domain_mod.f90:27)


def _file(label: str, v) -> str:
    if not isinstance(v, str) or not v or v != v.strip():
        raise ValueError(f"{label} must be a file name without surrounding blanks, got {v!r}")
    if any(c.isspace() or c in _FORBIDDEN or ord(c) < 32 for c in v):
        raise ValueError(f"{label} must be a plain file name (no blanks, commas, quotes, slashes, "
                         f"or comment characters), got {v!r}")
    if len(v) > _MAX_FILE:
        raise ValueError(f"{label} is longer than {_MAX_FILE} characters")
    if v.lower().startswith("internal"):
        raise ValueError(f"{label} must not start with 'internal' (DISCUS reads that as internal storage, "
                         "domain.f90:221, 2042)")
    return v


def _name(v) -> str:
    if not isinstance(v, str) or not _NAME.fullmatch(v):
        raise ValueError(f"pseudo-atom name must be 1 to 4 letters, digits or underscores, got {v!r}")
    return v.upper()  # DISCUS upper-cases names (domain.f90 2894, 2896)


def _real(label: str, v) -> float:
    if isinstance(v, bool) or not isinstance(v, numbers.Real) or not math.isfinite(float(v)):
        raise ValueError(f"{label} must be a finite number, got {v!r}")
    return float(v)


def _matrix(label: str, v, cols: int) -> tuple[tuple[float, ...], ...]:
    try:
        rows = tuple(tuple(r) for r in v)
    except TypeError:
        raise ValueError(f"{label} must be 3 rows of {cols} numbers, got {v!r}") from None
    if len(rows) != 3 or any(len(r) != cols for r in rows):
        raise ValueError(f"{label} must be 3 rows of {cols} numbers, got {v!r}")
    return tuple(tuple(_real(f"{label}[{i}][{j}]", x) for j, x in enumerate(r)) for i, r in enumerate(rows))


@dataclass(frozen=True)
class Domains(Step):
    """Insert the domains of ``content`` at the pseudo-atom positions of ``inputfile``.

    * ``content``: DISCUS ``.stru`` of the domain atoms, fractional coordinates of the host cell,
      same cell and space group as the host (``domain.f90`` 1083-1103).
    * ``inputfile``: DISCUS ``.stru`` whose atom lines name the pseudo atom (``name``, columns 1-4)
      and give its position (a normal ``.stru`` atom line). It must contain only atoms of this step's ``name``:
      an unknown name is an error (``domain.f90`` 2903-2906).
    * ``name``: pseudo-atom name, 1 to 4 characters.
    * ``character``: ``"cube"``, ``"cylinder"``, ``"sphere"`` or ``"fuzzy"``.
    * ``shape``: 3 rows ``(s_i1, s_i2, s_i3, t_i)`` in host cell units (see the module docstring).
      The 3x3 part must not be negative: DISCUS floors entries below 0.001 to 0.001
      (``domain.f90`` 2918), so a zero entry is stored as 0.001.
    * ``orient``: optional 3 rows ``(o_i1, o_i2, o_i3, t_i)``; ``None`` writes the identity.
    * ``fuzzy``: distance in Angstrom; host atoms closer than this to an inserted atom become
      vacancies. ``0`` disables the removal, so host atoms remain (overlap).
    * ``mode``: ``"pseudoatoms"`` (the only mode implemented here).

    ``apply`` runs ``reset`` first, so the step does not inherit definitions of an earlier
    domain step. The result is not a periodic host crystal (``domain.f90`` 1314-1316): see
    ``docs/DISCUS_DOMAIN_FACTS.md`` for the atom count and ordering, which decide whether
    ``StructureJob.run`` can read it back.
    """

    uses_random_numbers = True

    content: str
    inputfile: str
    name: str
    character: str
    shape: Sequence[Sequence[float]]
    orient: Sequence[Sequence[float]] | None = None
    fuzzy: float = 2.55  # DISCUS's own default (domain_reset, domain.f90 3220); 0 removes nothing
    mode: str = "pseudoatoms"

    def __post_init__(self):
        object.__setattr__(self, "content", _file("content", self.content))
        object.__setattr__(self, "inputfile", _file("inputfile", self.inputfile))
        object.__setattr__(self, "name", _name(self.name))
        if not isinstance(self.character, str) or self.character.lower() not in CHARACTERS:
            raise ValueError(f"character must be one of {CHARACTERS}, got {self.character!r}")
        object.__setattr__(self, "character", self.character.lower())
        shape = _matrix("shape", self.shape, 4)
        for i, row in enumerate(shape):
            for j in range(3):
                if row[j] < 0.0:
                    raise ValueError(f"shape[{i}][{j}] is negative ({row[j]!r}); DISCUS cannot store negative "
                                     "shape entries (domain.f90:2918)")
        object.__setattr__(self, "shape", shape)
        if self.orient is None:
            object.__setattr__(self, "orient", IDENTITY_ORIENT)
        else:
            object.__setattr__(self, "orient", _matrix("orient", self.orient, 4))
        fuzzy = _real("fuzzy", self.fuzzy)
        if fuzzy < 0.0:
            raise ValueError(f"fuzzy must be >= 0 Angstrom, got {self.fuzzy!r}")
        object.__setattr__(self, "fuzzy", fuzzy)
        if self.mode not in MODES:
            raise ValueError(f"mode must be one of {MODES} in this step, got {self.mode!r}")

    def input_names(self) -> tuple[str, ...]:
        """The files the job must stage: the domain content and the pseudo-atom input file."""
        return (self.content, self.inputfile)

    def apply(self, m: Macro) -> None:
        with m.menu("domain"):
            m.command("reset")  # domain.f90 168 (domain_reset): start from the system defaults
            m.command("mode", self.mode)
            m.command("inputfile", self.inputfile)
            m.command("assign", "character", self.name, self.character)
            m.command("assign", "content", self.name, self.content)
            m.command("assign", "fuzzy", self.name, self.fuzzy)
            for row, values in enumerate(self.shape, start=1):
                m.command("assign", "shape", self.name, row, *values)
            for row, values in enumerate(self.orient, start=1):
                m.command("assign", "orient", self.name, row, *values)
            m.command("run")


__all__ = ["CHARACTERS", "MODES", "IDENTITY_ORIENT", "Domains"]
