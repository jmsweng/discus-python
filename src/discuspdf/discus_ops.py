"""Apply DISCUS structure operations to a crystal and get the result back as a Structure.

A :class:`StructureJob` reads a crystal into DISCUS (a ``.cll`` unit cell repeated, a ``.stru``
file, or a :class:`~discuspdf.structure.Structure`), seeds the random generator, applies a list
of *steps*, saves the crystal and reads it back::

    from discuspdf.discus_ops import StructureJob, Thermal, ReplaceRandom
    from discuspdf.fourier_macro import Seed, StructureSource

    job = StructureJob(StructureSource.cell("cu.cll", (6, 6, 6)),
                       [ReplaceRandom("CU", "AU", 0.25), Thermal()],
                       seed=Seed((11, 22, 33, 44, 55, 66, 77, 88)))
    print(job.macro().render())                 # inspect
    result = job.run(runner, {"cu.cll": "cu.cll"})  # a SuiteRunner (stock or fork)
    result.structure                            # discuspdf Structure, configuration model

Each step is a small frozen dataclass whose ``apply(macro)`` writes its DISCUS lines into the
DISCUS main menu (or opens and closes its own sub-menu). The steps check their arguments; the
macro builder checks every command against the catalogue. Which physics a step produces is
documented with the source it comes from; the wrappers here are at the *engine-ran* level of
``docs/SUITE_CATALOGUE.md`` unless a step's docstring says otherwise.
"""

from __future__ import annotations

import math
import numbers
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Mapping, Sequence

from .fourier_macro import Seed, StructureSource
from .macro import Macro, MacroError, format_value

RESULT_NAME = "result.stru"
INPUT_NAME = "input.stru"


class Step:
    """One DISCUS operation. Subclasses implement :meth:`apply` (called in the DISCUS main menu)."""

    uses_random_numbers = False  # True when the step draws from DISCUS's random generator

    def input_names(self) -> tuple[str, ...]:
        """Files the step reads; the job's ``inputs`` must contain them."""
        return ()

    def apply(self, m: Macro) -> None:  # pragma: no cover - interface
        raise NotImplementedError


def _int(name: str, v, minimum: int = 1) -> int:
    if isinstance(v, bool) or not isinstance(v, numbers.Integral) or v < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}, got {v!r}")
    return int(v)


def _prob(name: str, v) -> float:
    if isinstance(v, bool) or not isinstance(v, numbers.Real) or not 0.0 <= float(v) <= 1.0:
        raise ValueError(f"{name} must be a probability in [0, 1], got {v!r}")
    return float(v)


def _atom(name: str, v) -> str | int:
    """An atom name (as DISCUS reads it) or a scattering-type number."""
    if isinstance(v, numbers.Integral) and not isinstance(v, bool):
        return _int(name, v, 0)
    if not isinstance(v, str) or not v.strip() or any(c in v for c in " ,:#!"):
        raise ValueError(f"{name} must be an atom name or type number, got {v!r}")
    return v.strip()


@dataclass(frozen=True)
class Raw(Step):
    """Lines written unchanged in the DISCUS main menu (checked against the catalogue on render).
    Set ``random=True`` if they draw random numbers, so the job insists on a seed."""

    lines: tuple
    random: bool = False

    @property
    def uses_random_numbers(self) -> bool:
        return self.random

    def apply(self, m: Macro) -> None:
        for line in self.lines:
            m.raw(line)


@dataclass(frozen=True)
class Thermal(Step):
    """Random Gaussian displacements from each atom's displacement parameters (``thermal``).

    ``discus/prog/thermal.f90`` (``ther_displ``): every atom is shifted by a Gaussian vector built
    from its type's ADPs (isotropic: sigma = sqrt(B / 8 pi^2) per Cartesian direction). ``two_d``
    leaves the directions in which the crystal is one cell thick unchanged (``ther 2d``);
    ``molecules`` moves whole molecules rigidly with the molecule B (``ther mol``). The random
    numbers come from the job's seed. Note DISCUS 6.23.01 reads an atom B below 0.001 as 1.0
    (``docs/FOURIER_CONTROLS.md``)."""

    uses_random_numbers = True

    two_d: bool = False
    molecules: bool = False

    def apply(self, m: Macro) -> None:
        params = (["mol"] if self.molecules else []) + (["2d"] if self.two_d else [])
        m.command("thermal", *params)


@dataclass(frozen=True)
class ReplaceAtom(Step):
    """Replace the atom with index ``index`` (1-based, DISCUS order) by ``atom``
    (``replace <index>,<at>[,<dbw>]``); ``b`` is required if the type is new to the crystal."""

    index: int
    atom: str | int
    b: float | None = None

    def __post_init__(self):
        _int("index", self.index)
        _atom("atom", self.atom)
        if self.b is not None and (isinstance(self.b, bool) or not isinstance(self.b, numbers.Real) or self.b < 0):
            raise ValueError("b must be a non-negative number")

    def apply(self, m: Macro) -> None:
        params = [self.index, self.atom] + ([self.b] if self.b is not None else [])
        m.command("replace", *params)


@dataclass(frozen=True)
class ReplaceRandom(Step):
    """Replace atoms of type ``old`` by ``new`` with probability ``probability`` on all sites or
    on one site of the unit cell (``replace <at1>,<at2>,{all|<site>},<prob>``). Uses the seed."""

    uses_random_numbers = True

    old: str | int
    new: str | int
    probability: float
    site: int | str = "all"

    def __post_init__(self):
        _atom("old", self.old)
        _atom("new", self.new)
        if _prob("probability", self.probability) == 0.0:
            raise ValueError("probability must be > 0 (modify_mod.f90 do_replace refuses 0)")
        if self.site != "all":
            _int("site", self.site)

    def apply(self, m: Macro) -> None:
        m.command("replace", self.old, self.new, self.site, self.probability)


@dataclass(frozen=True)
class Switch(Step):
    """Exchange the chemistry of atoms ``first`` and ``second`` (1-based indices; ``switch``)."""

    first: int
    second: int

    def __post_init__(self):
        _int("first", self.first)
        _int("second", self.second)

    def apply(self, m: Macro) -> None:
        m.command("switch", self.first, self.second)


@dataclass(frozen=True)
class Remove(Step):
    """Turn atoms ``start``..``end`` (1-based; ``start`` may be "all" or "last") into vacancies
    (``remove``). Follow with :class:`Purge` to drop them from the list."""

    start: int | str
    end: int | None = None

    def __post_init__(self):
        if self.start not in ("all", "last"):
            _int("start", self.start)
        if self.end is not None:
            _int("end", self.end)
            if self.start in ("all", "last"):
                raise ValueError("end is not used with start 'all' or 'last'")
            if self.end < self.start:
                raise ValueError("end must be >= start (DISCUS error -19)")

    def apply(self, m: Macro) -> None:
        m.command("remove", *([self.start] + ([self.end] if self.end is not None else [])))


@dataclass(frozen=True)
class Purge(Step):
    """Delete vacancies from the atom list (``purge``); ``types=True`` also drops unused atom types."""

    types: bool = False

    def apply(self, m: Macro) -> None:
        if self.types:
            m.command("purge", type="yes")
        else:
            m.command("purge")


WAVE_KINDS = ("longitudinal", "transverse", "density")
WAVE_FUNCTIONS = ("sinus", "box", "triangle")
WAVE_BRANCHES = ("acoustic", "optical")


def _real(name: str, v) -> float:
    if isinstance(v, bool) or not isinstance(v, numbers.Real) or not math.isfinite(float(v)):
        raise ValueError(f"{name} must be a finite number, got {v!r}")
    return float(v)


def _vec3(name: str, v) -> tuple[float, float, float]:
    try:
        items = tuple(v)
    except TypeError:
        raise ValueError(f"{name} must be three numbers, got {v!r}") from None
    if len(items) != 3:
        raise ValueError(f"{name} must be three numbers, got {v!r}")
    vec = tuple(_real(f"{name}[{i}]", x) for i, x in enumerate(items))
    if all(x == 0.0 for x in vec):
        raise ValueError(f"{name} must not be the zero vector")
    return vec


@dataclass(frozen=True)
class Wave(Step):
    """A displacement or density wave of the DISCUS ``waves`` menu (``discus/wave``).

    Atom-based modes only. Molecule selection (``mselect``, ``mreplace``) and rotation
    (``rotation``) are left out on purpose; so are the ``site:`` option of ``select``: the
    source parses it into ``wv_lsite`` (``modify_mod.f90`` 1862-1884) but ``wave_run`` never reads
    that array (``waves.f90`` 715-1036), so it would have no effect.

    * ``kind``: ``"longitudinal"`` (``long``), ``"transverse"`` (``tran``, needs
      ``oscillation``) or ``"density"`` (``dens``, needs ``replace``, ``plow``, ``phigh``).
    * ``vector`` (``vect``): wave vector in lattice units, not all zero. Only its direction is
      used: ``wave_run`` normalises it (``waves.f90`` 815-817).
    * ``length`` (``len``): wavelength in Angstrom, > 0 (source default 50).
    * ``function`` (``func``): ``"sinus"``, ``"box"`` or ``"triangle"``; ``box_width`` is the
      optional ``box,<delta>`` with 0 < delta < 1 (``waves.f90`` 242-243, default 0.5).
    * ``amplitude`` (``amp``): displacement amplitude in Angstrom, required for the displacement
      modes. ``sinus`` is ``amp*cos(2 pi arg) + shift`` (``waves.f90`` sinus function).
    * ``shift`` (``shift``): constant displacement in Angstrom, optional, displacement modes only.
    * ``phase`` (``phase``): degrees, or ``"random"`` (written unquoted, see below).
    * ``branch`` (``accoustic`` / ``optical``): displacement modes only, default acoustic. An
      optical wave reverses the displacement of atoms whose name contains '-' (``waves.f90``
      1019-1022). The density branch ignores it, so it must be None there.
    * ``select`` (``select``): displacement modes only: ``"all"`` (default) or names/type numbers.
    * ``replace`` (``replace <from>,<to>[,viceversa]``), ``plow``, ``phigh``: density modes
      only, 0 <= plow <= phigh <= 1. The wave profile f runs between ``plow`` and ``phigh``
      (``waves.f90`` 825-836), and an atom of type ``from`` is replaced by ``to`` when
      ``ran1 > f`` (``waves.f90`` 960): its replacement probability is ``1 - f``, i.e. between
      ``1 - phigh`` and ``1 - plow``. With ``viceversa`` an atom of type ``to`` becomes ``from``
      with probability ``f`` (``waves.f90`` 981-984).

    Density runs overwrite ``amp`` and ``shift`` (``waves.f90`` 825-836), so they must be None
    there. The macro is ``waves``: optional ``reset`` (only when the menu has that command), the
    settings, then ``run`` (``waves.f90`` 470-471, 715). ``select`` is needed because a reset
    leaves every atom unselected (``waves.f90`` 1310); without it ``run`` moves nothing.

    The source reads quoted values literally (``get_params_mod.f90`` 67-68 keeps the quotes and
    ``str_comp`` does not strip them), so ``func "box"``, ``phase "random"`` and a quoted
    ``"viceversa"`` would not be recognised. The step writes them unquoted.
    """

    kind: str
    vector: Sequence[float]
    length: float

    @property
    def uses_random_numbers(self) -> bool:
        return self.kind == "density" or self.phase == "random"

    function: str = "sinus"
    box_width: float | None = None
    amplitude: float | None = None
    shift: float | None = None
    phase: float | str = 0.0
    branch: str | None = None
    oscillation: Sequence[float] | None = None
    select: object = None
    replace: Sequence | None = None
    viceversa: bool = False
    plow: float | None = None
    phigh: float | None = None

    def __post_init__(self):
        if self.kind not in WAVE_KINDS:
            raise ValueError(f"kind must be one of {WAVE_KINDS}, got {self.kind!r}")
        density = self.kind == "density"
        object.__setattr__(self, "vector", _vec3("vector", self.vector))
        length = _real("length", self.length)
        if length <= 0:
            raise ValueError(f"length must be > 0 Angstrom, got {self.length!r}")
        object.__setattr__(self, "length", length)

        if self.function not in WAVE_FUNCTIONS:
            raise ValueError(f"function must be one of {WAVE_FUNCTIONS}, got {self.function!r}")
        if self.function == "box" and self.box_width is not None:
            width = _real("box_width", self.box_width)
            if not 0.0 < width < 1.0:
                raise ValueError(f"box_width must be in (0, 1), got {self.box_width!r}")
            object.__setattr__(self, "box_width", width)
        elif self.box_width is not None:
            raise ValueError("box_width is only accepted with function 'box'")

        if isinstance(self.phase, str):
            if self.phase != "random":
                raise ValueError(f"phase must be a number or 'random', got {self.phase!r}")
        else:
            object.__setattr__(self, "phase", _real("phase", self.phase))

        if density:
            self._check_density()
        else:
            self._check_displacement()

    def _check_density(self) -> None:
        for name in ("amplitude", "shift", "branch", "select", "oscillation"):
            if getattr(self, name) is not None:
                raise ValueError(f"{name} is not used by density waves (it is overwritten at run time)")
        if not isinstance(self.viceversa, bool):
            raise ValueError("viceversa must be True or False")
        if self.replace is None:
            raise ValueError("a density wave needs replace=(from, to)")
        pair = tuple(self.replace)
        if len(pair) != 2:
            raise ValueError(f"replace must be (from, to), got {self.replace!r}")
        object.__setattr__(self, "replace", (_atom("replace from", pair[0]), _atom("replace to", pair[1])))
        if self.plow is None or self.phigh is None:
            raise ValueError("a density wave needs plow and phigh")
        plow, phigh = _real("plow", self.plow), _real("phigh", self.phigh)
        if not 0.0 <= plow <= phigh <= 1.0:
            raise ValueError(f"need 0 <= plow <= phigh <= 1, got plow={plow!r}, phigh={phigh!r}")
        object.__setattr__(self, "plow", plow)
        object.__setattr__(self, "phigh", phigh)

    def _check_displacement(self) -> None:
        if self.amplitude is None:
            raise ValueError(f"a {self.kind} wave needs amplitude (Angstrom)")
        object.__setattr__(self, "amplitude", _real("amplitude", self.amplitude))
        if self.shift is not None:
            object.__setattr__(self, "shift", _real("shift", self.shift))
        branch = "acoustic" if self.branch is None else self.branch
        if branch not in WAVE_BRANCHES:
            raise ValueError(f"branch must be one of {WAVE_BRANCHES}, got {self.branch!r}")
        object.__setattr__(self, "branch", branch)

        select = self.select if self.select is not None else "all"
        if isinstance(select, (str, numbers.Integral)):
            select = (select,)
        items = tuple(_atom("select", x) for x in select)
        if not items:
            raise ValueError("select must name at least one atom")
        if "all" in items and len(items) > 1:
            raise ValueError("select 'all' cannot be combined with other atoms")
        object.__setattr__(self, "select", items)

        if self.kind == "transverse":
            if self.oscillation is None:
                raise ValueError("a transverse wave needs oscillation=(x, y, z)")
            object.__setattr__(self, "oscillation", _vec3("oscillation", self.oscillation))
        elif self.oscillation is not None:
            raise ValueError("oscillation is only used by transverse waves")

        if self.viceversa or self.replace is not None or self.plow is not None or self.phigh is not None:
            raise ValueError("replace, viceversa, plow and phigh are only used by density waves")

    def apply(self, m: Macro) -> None:
        density = self.kind == "density"
        with m.menu("waves"):
            if m.catalogue.menu("discus/wave").command("reset") is not None:
                m.command("reset")  # waves.f90 163-164; the menu has it, so the step starts from defaults
            if density:
                m.command("density")
            elif self.kind == "transverse":
                m.command("transverse")
            else:
                m.command("long")
            if not density:
                # waves.f90 182 spells it 'accoustic'; 'optical' at 294
                m.command("accoustic" if self.branch == "acoustic" else "optical")
            if self.box_width is not None:
                m.command("function", self.function, self.box_width)
            else:
                m.command("function", self.function)
            m.command("length", self.length)
            if not density:
                m.command("amplitude", self.amplitude)
            m.command("vector", *self.vector)
            m.command("phase", self.phase)
            if not density and self.shift is not None:
                m.command("shift", self.shift)
            if self.kind == "transverse":
                m.command("osci", *self.oscillation)
            if density:
                m.command("plow", self.plow)
                m.command("phigh", self.phigh)
                m.command("replace", *(list(self.replace) + (["viceversa"] if self.viceversa else [])))
            else:
                m.command("select", *self.select)
            m.command("run")


class ResultNotAStructure(ValueError):
    """The saved crystal could not be read back as a Structure (e.g. after ``purge`` the atom
    count is no longer a whole number per unit cell). ``run`` is the finished SuiteRun, so the
    file is still available as ``run.output("result.stru")``."""

    def __init__(self, message, run):
        super().__init__(message)
        self.run = run


def estimate_cells(xyz) -> tuple[int, int, int]:
    """DISCUS's own crystal size estimate from atom coordinates in unit cells: floor(span) + 1
    per axis (``estimate_ncells``, ``discus_estimate.f90`` 57-63), with a 1e-6 tolerance."""
    import numpy as np

    span = np.max(xyz, axis=0) - np.min(xyz, axis=0)
    return tuple(int(np.floor(float(v) + 1e-6)) + 1 for v in span)


def read_box(path, cells=None) -> tuple:
    """Read a DISCUS-saved crystal as ONE periodic box of ``cells`` unit cells, in any atom order
    and with any atom count (after ``purge``, domains, stacking...). ``cells=None`` uses DISCUS's
    own estimate from the atom extents (:func:`estimate_cells`). The lattice is the whole box
    (cell x cells), ``frac`` are fractions of it, ``supercell`` is (1, 1, 1) and no cell labels
    are kept. B values are taken as saved (DISCUS writes the B it used); occupancies are not
    kept. Returns ``(Structure, None)``."""
    from .discus_io import _parse_cell_line, _parse_stru_atoms, _split_directive, DiscusFormatError
    from .lattice import Lattice
    from .structure import Structure
    import numpy as np

    text = Path(path).read_text(errors="replace")
    cell = None
    for line_no, line in enumerate(text.splitlines(), start=1):
        key, rest = _split_directive(line)
        if key == "cell":
            cell = _parse_cell_line(rest, str(path), line_no)
        if key == "atoms":
            break
    if cell is None:
        raise DiscusFormatError(f"{path}: no cell line")
    atoms, _ = _parse_stru_atoms(text, str(path))
    if not atoms:
        raise DiscusFormatError(f"{path}: no atoms")
    species = np.array([at[0] for at in atoms], dtype=object)
    xyz = np.array([at[1][:3] for at in atoms], dtype=float)
    b_iso = np.array([at[1][3] for at in atoms], dtype=float)
    if np.any(b_iso < 0):
        raise DiscusFormatError(f"{path}: negative B values")
    n = np.asarray([int(c) for c in (cells if cells is not None else estimate_cells(xyz))], dtype=float)
    a, b, c, al, be, ga = cell.parameters
    box = Lattice.from_parameters(a * n[0], b * n[1], c * n[2], al, be, ga)
    structure = Structure(box, species, xyz / n, model="configuration", b_iso=b_iso, min_separation=None)
    return structure, None


class CrystalBuilder:
    """A source that *builds* the crystal inside DISCUS instead of reading one file (for example
    the stacking-fault generator). Subclasses implement:

    * ``build(m)``: write the DISCUS lines that create the crystal (called in the DISCUS main
      menu, before the seed is written, so a builder that draws random numbers writes its own
      seed via ``seed``);
    * ``input_names()``: the file names it reads, which the job's ``inputs`` must contain;
    * ``read_result(path, check_order)``: read the saved crystal back as ``(Structure, DiscusStru)``,
      raising :class:`~discuspdf.discus_io.DiscusFormatError` when it cannot.
    """

    uses_random_numbers = False

    def build(self, m: Macro, seed: Seed | None) -> None:  # pragma: no cover - interface
        raise NotImplementedError

    def input_names(self) -> tuple[str, ...]:  # pragma: no cover - interface
        raise NotImplementedError

    def read_result(self, path: Path, check_order: bool):  # pragma: no cover - interface
        raise NotImplementedError


@dataclass
class JobResult:
    structure: object  # discuspdf Structure
    run: object  # SuiteRun
    stru: object  # DiscusStru as read


@dataclass(frozen=True)
class StructureJob:
    """Read a crystal, apply ``steps`` in order, save it as ``result.stru``.

    ``source`` is a :class:`StructureSource` (the file is given to :meth:`run` as an input) or a
    :class:`~discuspdf.structure.Structure` (written as ``input.stru`` at run time). ``seed`` is
    written right after reading the crystal; every step that draws random numbers needs one for
    a reproducible result."""

    source: object
    steps: Sequence[Step] = field(default_factory=tuple)
    seed: Seed | None = None
    title: str | None = None
    allow_unseeded: bool = False
    readback: str = "cells"  # "cells": a Structure of unit cells; "box": one box (read_box; for a
    # CrystalBuilder the cell count is DISCUS's estimate from the atom extents); "none"

    def __post_init__(self):
        from .structure import Structure

        if not isinstance(self.source, (StructureSource, Structure, CrystalBuilder)):
            raise TypeError("source must be a StructureSource, a Structure or a CrystalBuilder")
        if isinstance(self.source, Structure) and self.source.model != "configuration":
            raise ValueError("a Structure source must be a configuration (a realised crystal)")
        steps = tuple(self.steps)
        for s in steps:
            if not isinstance(s, Step):
                raise TypeError(f"not a Step: {s!r}")
        object.__setattr__(self, "steps", steps)
        if self.seed is not None and not isinstance(self.seed, Seed):
            raise TypeError("seed must be a Seed (8 values)")
        if self.readback not in ("cells", "box", "none"):
            raise ValueError("readback must be 'cells', 'box' or 'none'")
        random_steps = [s for s in steps if s.uses_random_numbers]
        if isinstance(self.source, CrystalBuilder) and self.source.uses_random_numbers:
            random_steps.insert(0, self.source)
        if random_steps and self.seed is None and not self.allow_unseeded:
            raise ValueError(f"{type(random_steps[0]).__name__} draws random numbers; give seed=Seed(...) "
                             "(without one DISCUS seeds from the clock), or allow_unseeded=True")

    @property
    def _source(self) -> StructureSource:
        from .structure import Structure

        if isinstance(self.source, Structure):
            return StructureSource.stru(INPUT_NAME, tuple(int(n) for n in self.source.supercell))
        return self.source

    def macro(self) -> Macro:
        m = Macro(self.title or "discuspdf StructureJob")
        if isinstance(self.source, CrystalBuilder):
            with m.program("discus"):
                # the builder writes the seed itself where it needs it (before drawing)
                self.source.build(m, self.seed)
                for step in self.steps:
                    step.apply(m)
                m.command("save", RESULT_NAME)
            return m
        src = self._source
        with m.program("discus"):
            with m.menu("read"):
                if src.kind == "stru":
                    m.command("structure", src.file)
                else:
                    n = src.crystal_cells
                    m.command("cell", src.file, *n)
            if self.seed is not None:
                m.command("seed", *self.seed.values)
            for step in self.steps:
                step.apply(m)
            m.command("save", RESULT_NAME)
        return m

    def crystal_cells(self, inputs: Mapping[str, Path | str] = ()) -> tuple[int, int, int]:
        src = self._source
        if src.crystal_cells is not None:
            return src.crystal_cells
        from .discus_io import read_stru

        path = Path(dict(inputs)[src.file])
        for line in path.read_text(errors="replace").splitlines():
            words = line.replace(",", " ").split()
            if words and words[0].lower() == "ncell" and len(words) >= 4:
                return tuple(int(float(w)) for w in words[1:4])
            if words and words[0].lower() == "atoms":
                break
        raise ValueError(f"{src.file} has no ncell line; give StructureSource.stru(..., crystal_cells=(n1, n2, n3))")

    def run(self, runner, inputs: Mapping[str, Path | str] | None = None, *,
            engine: str | None = None, check_order: bool = True) -> JobResult:
        """Run with a :class:`~discuspdf.suite_runner.SuiteRunner` and read ``result.stru`` back
        (DISCUS ``save`` writes no ``ncell`` line, so the crystal size comes from the source).
        ``check_order=False`` skips the cell-major order check, e.g. after displacements larger
        than half a cell; the cell labels are then assigned by row and are not checked."""
        from .discus_io import read_stru, write_stru
        from .structure import Structure

        inputs = dict(inputs or {})
        needed = [n for step in self.steps for n in step.input_names() if n not in inputs]
        if needed:
            raise ValueError(f"inputs must contain {sorted(set(needed))} (files the steps read)")
        with tempfile.TemporaryDirectory() as tmp:
            if isinstance(self.source, CrystalBuilder):
                missing = [n for n in self.source.input_names() if n not in inputs]
                if missing:
                    raise ValueError(f"inputs must contain {missing}")
                cells = None
            elif isinstance(self.source, Structure):
                path = Path(tmp) / INPUT_NAME
                write_stru(path, self.source, check_order=False, min_separation=None)
                inputs[INPUT_NAME] = path
                cells = self.crystal_cells(inputs)
            elif self.source.file not in inputs:
                raise ValueError(f"inputs must contain {self.source.file!r}")
            else:
                cells = self.crystal_cells(inputs)
            run = runner.run(self.macro(), inputs, [RESULT_NAME], engine=engine)
        from .discus_io import DiscusFormatError

        if self.readback == "none":
            return JobResult(None, run, None)
        if self.readback == "box":
            try:
                structure, _ = read_box(run.output(RESULT_NAME), cells)
            except DiscusFormatError as exc:
                raise ResultNotAStructure(f"DISCUS ran, but {RESULT_NAME} could not be read: {exc}", run) from exc
            return JobResult(structure, run, None)
        if isinstance(self.source, CrystalBuilder):
            try:
                structure, stru = self.source.read_result(run.output(RESULT_NAME), check_order)
            except DiscusFormatError as exc:
                raise ResultNotAStructure(f"DISCUS ran, but {RESULT_NAME} could not be read back: {exc}", run) from exc
            return JobResult(structure, run, stru)

        try:
            # read cell / read stru keep DISCUS's cell-major order, so the order check holds unless a
            # purge shifted rows; vacancies make the Property column non-constant (do_remove clears
            # PROP_NORMAL), which is the only column allowed to vary
            stru = read_stru(run.output(RESULT_NAME), supercell=cells, check_order=check_order,
                             min_separation=None, allow_loss=True)
        except DiscusFormatError as exc:
            raise ResultNotAStructure(f"DISCUS ran, but {RESULT_NAME} is not a {cells} crystal: {exc}; "
                                      "readback='box' reads it as one periodic box", run) from exc
        extra = [c for c in stru.dropped if c != "Property"]
        if extra:
            raise ResultNotAStructure(f"{RESULT_NAME} has non-constant {extra} columns the Structure cannot keep", run)
        return JobResult(stru.structure, run, stru)


__all__ = ["estimate_cells", "read_box", "CrystalBuilder", "ResultNotAStructure", "Step", "Raw", "Thermal", "ReplaceAtom", "ReplaceRandom", "Switch", "Remove", "Purge",
           "Wave", "StructureJob", "JobResult", "MacroError", "format_value"]
