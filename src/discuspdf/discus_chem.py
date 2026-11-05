"""DISCUS's own correlation analysis (the ``chem`` menu), run on a crystal and read back as data.

The ``chem corr`` calculations are the engine's side of the comparison documented in
``benchmarks/discus_6_23_01/c_probes/RESULT.md``: :class:`ChemAnalysis` writes the neighbour
definitions and the correlations of a probe macro, :func:`parse_chem_log` reads every correlation
the log prints, and :func:`run_chem` runs a :class:`~discuspdf.discus_ops.StructureJob` and returns
the numbers together with the kept run::

    from discuspdf.discus_chem import ChemAnalysis, ChemOccupational, Neighbourhood, run_chem
    from discuspdf.mc import Neighborhood, Vector
    from discuspdf.discus_ops import StructureJob
    from discuspdf.fourier_macro import StructureSource

    analysis = ChemAnalysis(
        neighbourhoods=(Neighbourhood(Neighborhood((Vector(1, 2, (0, 0, 0)),), symmetric=False)),),
        calculations=(ChemOccupational("Au", "Cu"),),
    )
    job = StructureJob(StructureSource.stru("occ_bcc.stru"), [analysis], readback="none")
    outcome = run_chem(job, runner, {"occ_bcc.stru": "occ_bcc.stru"})
    outcome.results[0].value          # DISCUS value (12 significant digits)
    outcome.engine_crashed_at_exit    # True when DISCUS aborted after printing everything

Commands are those of the probe macros (``inputs/*.mac``), in the same order: ``set mode, quick,
period``; ``set vec`` for each vector of the vector neighbourhoods (numbered 1, 2, ... in order); ``set neig, rese``; per neighbourhood
``set neig, vec`` or ``set neig, dis`` (and ``set neig, dir`` for a displacement direction); ``show``;
then per calculation ``corr occ`` or ``corr disp`` followed by ``eval res[k]`` for k = 0..N
(``res[0]`` is the number of neighbourhoods, ``res[k]`` the value of neighbourhood k). The menu is
entered with ``chem``, the spelling the probes use (``chemistry`` is the same command).

What the commands compute, as checked against the probe runs (see ``RESULT.md``):

* ``corr occ, a, b``: occupational correlation of two different species over the pairs; DISCUS
  pools the two species' concentrations and the result does not depend on the order of a and b.
  Python side: :meth:`discuspdf.correlations.OccupationalCorrelation.discus`.
* ``corr disp, a, b`` (``a``/``b`` may be ``all``): the **uncentred** directional correlation of
  the displacements of the pairs, with displacements measured from the per-site average.
  Python side: :func:`discuspdf.correlations.displacement_correlation` with
  ``reference="site_average", centered=False, direction_frame="fractional"``.

Known engine issue: DISCUS 6.23.01 aborted with ``free(): invalid pointer`` (SIGABRT, exit 134) after
``chem corr occ`` followed by ``exit``, after printing every result (``occ_bcc`` in the probes).
:class:`~discuspdf.suite_runner.SuiteRunner` treats a nonzero exit as a failed run and is not
changed for this. :func:`run_chem` accepts that one case only when the log holds every requested
correlation, the signature follows the last result line, the status is exit 134 and the saved
``result.stru`` is present; it then returns with ``engine_crashed_at_exit=True``. Whether the crash
also happens for a macro that ends without the suite-level ``exit`` (as :class:`StructureJob`
writes it) has not been tested.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Sequence

from .discus_ops import RESULT_NAME, Step, StructureJob
from .engine import _ANSI
from .fourier_macro import StructureSource
from .macro import Macro
from .suite_runner import SuiteRun
from .atom_names import check_discus_name
from .mc import Neighborhood, Vector
from ._validate import integer as _whole, real as _real, vec3 as _triple_real

CRASH_SIGNATURE = "free(): invalid pointer"
CRASH_RETURN_CODE = 134  # 128 + SIGABRT, as the shell reports it


class ChemIncomplete(RuntimeError):
    """The engine finished, but its log does not hold every requested correlation."""


# ----------------------------------------------------------------------------- validation helpers

def _species(name: str, v, allow_all: bool) -> str:
    """An atom type name, or a type number (written as its digits), or 'all' where DISCUS allows it."""
    if allow_all is False and isinstance(v, str) and v.strip().lower() == "all":
        raise ValueError(f"{name}: 'all' is only allowed for displacement correlations")
    if isinstance(v, int) and not isinstance(v, bool):
        return str(_whole(name, v))
    return str(check_discus_name(name, v, allow_all=allow_all))


# ------------------------------------------------------------------------------- the definitions

@dataclass(frozen=True)
class DistanceNeighbourhood:
    """``set neig, dis, <x>, <y>, <z>, <sigma>``: the pairs whose distance lies within
    ``sigma`` of the length of the lattice vector (x, y, z) (DISCUS ``fsig``; the probes use 0.1)."""

    direction: tuple[float, float, float]
    sigma: float

    def __post_init__(self):
        object.__setattr__(self, "direction", _triple_real("direction", self.direction))
        sigma = _real("sigma", self.sigma)
        if sigma <= 0:
            raise ValueError(f"sigma must be > 0, got {self.sigma!r}")
        object.__setattr__(self, "sigma", sigma)


@dataclass(frozen=True)
class Neighbourhood:
    """One numbered neighbourhood (``number:k`` is its 1-based position in the analysis).

    ``definition`` is either a :class:`discuspdf.mc.Neighborhood` (the same vectors the ``mmc`` menu
    uses; written as one ``set vec`` per vector and ``set neig, vec, <ids>``, numbered as
    :meth:`discuspdf.mc.MCConfig.mmc_lines` numbers them) or a :class:`DistanceNeighbourhood`.
    ``mc.Neighborhood`` adds the reversed vectors by default (``symmetric=True``); the probe macros
    use one direction only, so pass ``symmetric=False`` to reproduce them. In a periodic crystal the
    reversed vector gives the same pairs the other way round, so an occupational value does not change;
    the pair count doubles.

    ``direction`` is the displacement direction pair ``(a, b)`` of ``set neig, dir, a..., b...``
    used by displacement correlations; it is only written for a vector neighbourhood, the only form
    the probes use with it. ``None`` leaves DISCUS's default ("all directions")."""

    definition: Neighborhood | DistanceNeighbourhood
    direction: tuple[tuple[float, float, float], tuple[float, float, float]] | None = None

    def __post_init__(self):
        if not isinstance(self.definition, (Neighborhood, DistanceNeighbourhood)):
            raise TypeError("definition must be a discuspdf.mc.Neighborhood or a DistanceNeighbourhood")
        if self.direction is not None:
            if not isinstance(self.definition, Neighborhood):
                raise ValueError("a displacement direction is only written for a vector neighbourhood")
            pair = tuple(self.direction)
            if len(pair) != 2:
                raise ValueError(f"direction must be (a, b), got {self.direction!r}")
            object.__setattr__(self, "direction", (_triple_real("direction a", pair[0]),
                                                   _triple_real("direction b", pair[1])))


@dataclass(frozen=True)
class ChemOccupational:
    """``corr occ, <a>, <b>``: occupational correlation between two different atom types."""

    species_a: str
    species_b: str
    word = "occ"

    def __post_init__(self):
        a = _species("species_a", self.species_a, allow_all=False)
        b = _species("species_b", self.species_b, allow_all=False)
        if a.lower() == b.lower():
            raise ValueError("corr occ needs two different atom types")
        object.__setattr__(self, "species_a", a)
        object.__setattr__(self, "species_b", b)


@dataclass(frozen=True)
class ChemDisplacement:
    """``corr disp, <a>, <b>``: displacement correlation; ``a`` or ``b`` may be ``"all"``."""

    species_a: str
    species_b: str
    word = "disp"

    def __post_init__(self):
        object.__setattr__(self, "species_a", _species("species_a", self.species_a, allow_all=True))
        object.__setattr__(self, "species_b", _species("species_b", self.species_b, allow_all=True))


# ----------------------------------------------------------------------------------------- the step

@dataclass(frozen=True)
class ChemAnalysis(Step):
    """The ``chem`` analysis of a crystal: numbered neighbourhoods and the correlations to calculate.
    It writes the lines into the DISCUS main menu (``chem`` ... ``exit``) and does not change the
    crystal, so it can sit between the crystal steps and ``save``. The correlation vectors are numbered
    1, 2, ... in the order the vector neighbourhoods list them."""

    neighbourhoods: Sequence[Neighbourhood] = ()
    calculations: Sequence[ChemOccupational | ChemDisplacement] = ()

    def __post_init__(self):
        nbs = tuple(self.neighbourhoods)
        for nb in nbs:
            if not isinstance(nb, Neighbourhood):
                raise TypeError(f"not a Neighbourhood: {nb!r}")
        calcs = tuple(self.calculations)
        for c in calcs:
            if not isinstance(c, (ChemOccupational, ChemDisplacement)):
                raise TypeError(f"not a calculation: {c!r}")
        if not calcs:
            raise ValueError("ChemAnalysis needs at least one calculation (ChemOccupational or ChemDisplacement)")
        if not nbs:
            raise ValueError("a correlation needs at least one neighbourhood")
        object.__setattr__(self, "neighbourhoods", nbs)
        object.__setattr__(self, "calculations", calcs)

    def vector_numbers(self) -> list[tuple[int, ...]]:
        """The ``set vec`` numbers of each neighbourhood (empty for a distance neighbourhood)."""
        out, iv = [], 0
        for nb in self.neighbourhoods:
            if isinstance(nb.definition, Neighborhood):
                n = len(nb.definition.vectors)
                out.append(tuple(range(iv + 1, iv + n + 1)))
                iv += n
            else:
                out.append(())
        return out

    def apply(self, m: Macro) -> None:
        numbers = self.vector_numbers()
        with m.menu("chem"):
            m.command("set", "mode", "quick", "period")
            for nb, ids in zip(self.neighbourhoods, numbers):
                if isinstance(nb.definition, Neighborhood):
                    for i, v in zip(ids, nb.definition.vectors):
                        m.command("set", "vec", i, v.origin_site, v.target_site, *v.translation)
            m.command("set", "neig", "rese")
            for k, (nb, ids) in enumerate(zip(self.neighbourhoods, numbers), start=1):
                d = nb.definition
                if isinstance(d, Neighborhood):
                    m.command("set", "neig", "vec", *ids, number=k)
                else:
                    m.command("set", "neig", "dis", *d.direction, d.sigma, number=k)
                if nb.direction is not None:
                    m.command("set", "neig", "dir", *nb.direction[0], *nb.direction[1], number=k)
            m.command("show")
            for calc in self.calculations:
                m.command("corr", calc.word, calc.species_a, calc.species_b)
                for k in range(len(self.neighbourhoods) + 1):
                    m.command("eval", f"res[{k}]")


# ----------------------------------------------------------------------------------------- parsing

@dataclass(frozen=True)
class ChemResult:
    """One correlation the log prints: one neighbourhood of one ``corr`` command."""

    block: int  # index of the corr command in the log (0-based)
    kind: str  # "occ" or "disp"
    species_a: str  # as typed in the corr command
    species_b: str
    neighbour: int  # 1-based neighbourhood number
    definition: str | None  # from the show table: "vec 1: site 1 -> site 2, neig 0 0 0" or "dis 0.50 0.50 0.50 fsig 0.10"
    pairs: int  # '# pairs' column
    value: float  # DISCUS's 12-digit eval value when present, else the table value
    value_exact: bool  # True when value is the eval value (12 significant digits)
    table_value: float  # the value in the table (4 decimals)
    occupation_percent: tuple[float, float, float] | None  # (AA, AB, BB) in percent, occ only
    direction_a: tuple[float, float, float] | None  # disp only
    direction_b: tuple[float, float, float] | None


_CORR_CMD = re.compile(r"discus/chem>\s*corr\s+(occ|disp)\s*,\s*([^,\s]+)\s*,\s*([^,\s]+)", re.IGNORECASE)
_CORR_VECTOR = re.compile(r"Correlation vector\s+(\d+)\s*:\s*Site\s+(\d+)\s*->\s*site\s+(\d+),"
                          r"\s*neig\s*=\s*(-?\d+)\s+(-?\d+)\s+(-?\d+)")
_DEF_ROW = re.compile(r"^\s*(\d+)\s+(vec|neig)\s+(\S.*?)\s*$")
_RES = re.compile(r"Value of res\[(\d+)\]\s*=\s*(\S+)")
_OCC_ROW = re.compile(r"^\s*(\d+)\s+(\d+\.\d+)\s*%\s+(\d+\.\d+)\s*%\s+(\d+\.\d+)\s*%\s+(\d+)\s+(-?\d+\.\d+)\s*$")
_F = r"\s+(-?\d+\.\d+)"
_DISP_ROW = re.compile(r"^\s*(\d+)" + _F * 6 + r"\s+(\d+)" + _F + r"\s*$")


@dataclass
class _Row:
    pairs: int
    table: float
    pct: tuple[float, float, float] | None = None
    dir_a: tuple[float, float, float] | None = None
    dir_b: tuple[float, float, float] | None = None


@dataclass
class _Block:
    index: int
    kind: str
    a: str
    b: str
    defs: dict  # neighbour number -> ("vec", n) or ("neig", "x y z fsig") from the last show table
    vectors: dict  # correlation vector number -> "site 1 -> site 2, neig 0 0 0"
    rows: dict  # neighbour number -> _Row
    res: dict  # neighbour number -> eval value


def _describe(defn, vectors: dict) -> str | None:
    if defn is None:
        return None
    kind, text = defn
    if kind == "vec":
        n = int(text)
        return f"vec {n}: {vectors[n]}" if n in vectors else f"vec {n}"
    return f"dis {text}"


def parse_chem_log(text: str) -> list[ChemResult]:
    """Every correlation a DISCUS log prints, in the order of the ``corr`` commands.

    Reads the ``corr occ``/``corr disp`` echo, the table rows (neighbour, percentages for occ,
    directions for disp, pair count, correlation) and the ``Value of res[k]`` lines that follow
    each ``corr``; the neighbour definitions come from the last ``show`` table before the ``corr``.
    A ``corr`` whose table is missing still counts as a block, so a truncated log is detected by
    the caller, not silently shortened. ANSI colour codes are removed first.
    """
    text = _ANSI.sub("", text)
    defs: dict = {}
    vectors: dict = {}
    blocks: list[_Block] = []
    current: _Block | None = None
    for line in text.splitlines():
        m = _CORR_CMD.search(line)
        if m:
            current = _Block(len(blocks), m.group(1).lower(), m.group(2), m.group(3),
                             dict(defs), dict(vectors), {}, {})
            blocks.append(current)
            continue
        m = _CORR_VECTOR.search(line)
        if m:
            n, a, b, u, v, w = (int(x) for x in m.groups())
            vectors[n] = f"site {a} -> site {b}, neig {u} {v} {w}"
            continue
        m = _DEF_ROW.match(line)
        if m:
            k, mode, rest = int(m.group(1)), m.group(2), m.group(3)
            if mode == "vec":
                defs[k] = ("vec", rest.split()[0])
            else:
                tokens = rest.split()
                defs[k] = ("neig", f"{tokens[0]} {tokens[1]} {tokens[2]} fsig {tokens[3]}") \
                    if len(tokens) >= 4 else None
            continue
        if current is None:
            continue
        m = _RES.search(line)
        if m:
            try:
                current.res[int(m.group(1))] = float(m.group(2))
            except ValueError:
                pass
            continue
        if current.kind == "occ":
            m = _OCC_ROW.match(line)
            if m:
                k = int(m.group(1))
                current.rows[k] = _Row(int(m.group(5)), float(m.group(6)),
                                       pct=(float(m.group(2)), float(m.group(3)), float(m.group(4))))
        else:
            m = _DISP_ROW.match(line)
            if m:
                k = int(m.group(1))
                f = [float(x) for x in m.groups()[1:7]]
                current.rows[k] = _Row(int(m.group(8)), float(m.group(9)),
                                       dir_a=tuple(f[0:3]), dir_b=tuple(f[3:6]))

    results: list[ChemResult] = []
    for blk in blocks:
        for k in sorted(blk.rows):
            row = blk.rows[k]
            exact = blk.res.get(k)
            results.append(ChemResult(
                block=blk.index, kind=blk.kind, species_a=blk.a, species_b=blk.b, neighbour=k,
                definition=_describe(blk.defs.get(k), blk.vectors), pairs=row.pairs,
                value=row.table if exact is None else exact, value_exact=exact is not None,
                table_value=row.table, occupation_percent=row.pct,
                direction_a=row.dir_a, direction_b=row.dir_b))
    return results


# ----------------------------------------------------------------------------------------- running

@dataclass(frozen=True)
class ChemOutcome:
    results: tuple  # of ChemResult, in log order
    run: SuiteRun  # the kept run directory (run.mac, run.log, result.stru, status.json)
    engine_crashed_at_exit: bool = False


def _single_chem(job: StructureJob) -> ChemAnalysis:
    found = [s for s in job.steps if isinstance(s, ChemAnalysis)]
    if len(found) != 1:
        raise ValueError(f"the job needs exactly one ChemAnalysis step, found {len(found)}")
    return found[0]


def _problems(results, chem: ChemAnalysis) -> list[str]:
    """Why ``results`` does not hold every correlation of ``chem`` (empty list = complete)."""
    by_block: dict[int, list[ChemResult]] = {}
    for r in results:
        by_block.setdefault(r.block, []).append(r)
    out: list[str] = []
    expected_blocks = set(range(len(chem.calculations)))
    if set(by_block) - expected_blocks:
        out.append(f"unexpected corr blocks {sorted(set(by_block) - expected_blocks)}")
    want_k = list(range(1, len(chem.neighbourhoods) + 1))
    for i, calc in enumerate(chem.calculations):
        rows = sorted(by_block.get(i, []), key=lambda r: r.neighbour)
        label = f"corr {calc.word} {calc.species_a},{calc.species_b}"
        if [r.neighbour for r in rows] != want_k:
            out.append(f"{label}: neighbourhoods {[r.neighbour for r in rows]}, expected {want_k}")
        elif any((r.kind, r.species_a.lower(), r.species_b.lower())
                 != (calc.word, calc.species_a.lower(), calc.species_b.lower()) for r in rows):
            out.append(f"{label}: printed for other species")
        elif not all(r.value_exact for r in rows):
            out.append(f"{label}: missing eval value")
    return out


def _log_text(run: SuiteRun) -> str:
    return run.log.read_text(errors="replace") if run.log.is_file() else ""


def _crashed_after_results(run: SuiteRun, results, chem: ChemAnalysis) -> bool:
    try:
        status = run.status()
    except (OSError, ValueError):
        return False
    if status.get("status") != "failed" or status.get("returncode") != CRASH_RETURN_CODE:
        return False
    text = _ANSI.sub("", _log_text(run))
    at = text.find(CRASH_SIGNATURE)
    if at < 0 or text.rfind("Value of res[", 0, at) < 0:
        return False
    if _problems(results, chem):
        return False
    out = run.directory / RESULT_NAME
    return out.is_file() and out.stat().st_size > 0


def run_chem(job: StructureJob, runner, inputs=(), *, engine: str | None = None) -> ChemOutcome:
    """Run ``job`` (a :class:`StructureJob` on a file source with one :class:`ChemAnalysis` step)
    with a :class:`~discuspdf.suite_runner.SuiteRunner` and return the parsed correlations.

    ``inputs`` maps the crystal file name to its path. The job's ``readback`` is not used: the
    crystal is saved as ``result.stru`` in the kept run, and nothing is read back.

    Raises :class:`RuntimeError` when the run fails, unless the engine exited with
    :data:`CRASH_RETURN_CODE` after printing every correlation (then the outcome is returned with
    ``engine_crashed_at_exit=True``). Raises :class:`ChemIncomplete` when a run that exited cleanly
    does not hold every correlation.
    """
    if not isinstance(job, StructureJob):
        raise TypeError("run_chem takes a StructureJob")
    chem = _single_chem(job)
    if not isinstance(job.source, StructureSource):
        raise TypeError("run_chem needs a file source (StructureSource); write the crystal to a file first")
    inputs = dict(inputs)
    needed = {job.source.file} | {n for s in job.steps for n in s.input_names()}
    missing = sorted(n for n in needed if n not in inputs)
    if missing:
        raise ValueError(f"inputs must contain {missing}")
    run = runner.prepare(job.macro(), inputs, [RESULT_NAME], engine=engine)
    try:
        runner.execute(run)
    except RuntimeError:
        results = parse_chem_log(_log_text(run))
        if _crashed_after_results(run, results, chem):
            return ChemOutcome(tuple(results), run, engine_crashed_at_exit=True)
        raise
    results = parse_chem_log(_log_text(run))
    problems = _problems(results, chem)
    if problems:
        raise ChemIncomplete(f"DISCUS exited cleanly but the log is incomplete ({'; '.join(problems)}); "
                             f"see {run.directory}")
    return ChemOutcome(tuple(results), run, engine_crashed_at_exit=False)


__all__ = [
    "CRASH_RETURN_CODE", "CRASH_SIGNATURE", "ChemAnalysis", "ChemDisplacement", "ChemIncomplete",
    "ChemOccupational", "ChemOutcome", "ChemResult", "DistanceNeighbourhood", "Neighborhood",
    "Neighbourhood", "Vector", "parse_chem_log", "run_chem",
]
