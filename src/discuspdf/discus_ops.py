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
    """Lines written unchanged in the DISCUS main menu (checked against the catalogue on render)."""

    lines: tuple

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

    old: str | int
    new: str | int
    probability: float
    site: int | str = "all"

    def __post_init__(self):
        _atom("old", self.old)
        _atom("new", self.new)
        _prob("probability", self.probability)
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

    def __post_init__(self):
        from .structure import Structure

        if not isinstance(self.source, (StructureSource, Structure)):
            raise TypeError("source must be a StructureSource or a Structure")
        if isinstance(self.source, Structure) and self.source.model != "configuration":
            raise ValueError("a Structure source must be a configuration (a realised crystal)")
        steps = tuple(self.steps)
        for s in steps:
            if not isinstance(s, Step):
                raise TypeError(f"not a Step: {s!r}")
        object.__setattr__(self, "steps", steps)
        if self.seed is not None and not isinstance(self.seed, Seed):
            raise TypeError("seed must be a Seed (8 values)")

    @property
    def _source(self) -> StructureSource:
        from .structure import Structure

        if isinstance(self.source, Structure):
            return StructureSource.stru(INPUT_NAME, tuple(int(n) for n in self.source.supercell))
        return self.source

    def macro(self) -> Macro:
        m = Macro(self.title or "discuspdf StructureJob")
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
        ncell = read_stru(path, check_order=False, min_separation=None).ncell
        if ncell is None:
            raise ValueError(f"{src.file} has no ncell line; give StructureSource.stru(..., crystal_cells)")
        return tuple(int(n) for n in ncell[:3])

    def run(self, runner, inputs: Mapping[str, Path | str] | None = None, *,
            engine: str | None = None) -> JobResult:
        """Run with a :class:`~discuspdf.suite_runner.SuiteRunner` and read ``result.stru`` back
        (DISCUS ``save`` writes no ``ncell`` line, so the crystal size comes from the source)."""
        from .discus_io import read_stru, write_stru
        from .structure import Structure

        inputs = dict(inputs or {})
        with tempfile.TemporaryDirectory() as tmp:
            if isinstance(self.source, Structure):
                path = Path(tmp) / INPUT_NAME
                write_stru(path, self.source, check_order=False, min_separation=None)
                inputs[INPUT_NAME] = path
            elif self.source.file not in inputs:
                raise ValueError(f"inputs must contain {self.source.file!r}")
            cells = self.crystal_cells(inputs)
            run = runner.run(self.macro(), inputs, [RESULT_NAME], engine=engine)
        stru = read_stru(run.output(RESULT_NAME), supercell=cells, check_order=False, min_separation=None)
        return JobResult(stru.structure, run, stru)


__all__ = ["Step", "Raw", "Thermal", "ReplaceAtom", "ReplaceRandom", "Switch", "Remove", "Purge",
           "StructureJob", "JobResult", "MacroError", "format_value"]
