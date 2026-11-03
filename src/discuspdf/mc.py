"""Python-configured DISCUS Monte Carlo (MMC) runs (Package D, Features 15–17).

``MCConfig`` describes one DISCUS 6.23.01 ``mmc`` run — neighbourhoods built from ``set vec`` vectors, targets
(energies), moves, number of trial moves, temperature and an 8-value seed — and renders the exact macro text.
``MCRunner`` stages the structure, macro and configuration in a fresh run directory, executes DISCUS, keeps the
log and loads the final structure; ``MCResult`` holds the parsed feedback history and the final table.

Only engine behaviour that was read in the source *and* checked by the engine probe
(``benchmarks/discus_6_23_01/d_mc_probes``) is offered; see ``docs/MONTE_CARLO.md`` and ``docs/MC_FACTS.md``.
In particular DISCUS's MMC acceptance rule is **not Metropolis**: downhill moves are always accepted and
uphill or neutral moves with probability ``1/(1 + exp(ΔE/kT))`` (heat bath; ΔE = 0 → ½). Energies are in
units of kT.
"""
from __future__ import annotations

import hashlib
import json
import re
import shutil
import subprocess
import tempfile
import uuid
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

from .engine import current_engine, Engine, as_engine
from .fourier_macro import Seed
from .structure import Structure

__all__ = [
    "Vector", "Neighborhood", "OccupationalTarget", "DisplacementTarget", "SpringTarget",
    "LennardJonesTarget", "Swap", "Shift", "SwapDisplacement", "MCConfig", "MCRunner", "MCRun", "MCResult",
    "parse_mmc_log", "match_row", "replica_summary",
]

_NAME = re.compile(r"^[A-Za-z][A-Za-z0-9]{0,3}$")
_PAIRING = ("all", "local", "slocal", "site")


def _species(name, what: str) -> str:
    if not isinstance(name, str) or not _NAME.match(name):
        raise ValueError(f"{what} must be a DISCUS atom name (letter + up to 3 letters/digits), got {name!r}")
    return name


def _finite(x, what: str) -> float:
    try:
        v = float(x)
    except (TypeError, ValueError):
        raise ValueError(f"{what} must be a number, got {x!r}") from None
    if not np.isfinite(v):
        raise ValueError(f"{what} must be finite")
    return v


def _fmt(v: float) -> str:
    """Plain decimal (no exponent), shortest round-trip digits; e.g. 1e-05 -> '0.00001'."""
    s = np.format_float_positional(float(v), trim="-")
    return s + ".0" if "." not in s else s


# ------------------------------------------------------------------------------------------- neighbourhoods
@dataclass(frozen=True)
class Vector:
    """``set vec``: from site ``origin_site`` (1-based slot in the unit cell, DISCUS numbering) to site
    ``target_site`` in the cell displaced by the integer ``translation``."""
    origin_site: int
    target_site: int
    translation: tuple[int, int, int]

    def __post_init__(self) -> None:
        for name in ("origin_site", "target_site"):
            v = getattr(self, name)
            if isinstance(v, bool) or not isinstance(v, (int, np.integer)) or v < 1:
                raise ValueError(f"{name} must be a 1-based site number, got {v!r}")
        t = tuple(self.translation)
        if len(t) != 3 or any(isinstance(x, bool) or not isinstance(x, (int, np.integer)) for x in t):
            raise ValueError("translation must be three integers (unit cells)")
        object.__setattr__(self, "origin_site", int(self.origin_site))
        object.__setattr__(self, "target_site", int(self.target_site))
        object.__setattr__(self, "translation", tuple(int(x) for x in t))

    def reversed(self) -> "Vector":
        return Vector(self.target_site, self.origin_site, tuple(-x for x in self.translation))


@dataclass(frozen=True)
class Neighborhood:
    """A group of vectors (``set neig, vec, …``). DISCUS counts each vector **in one direction only**
    (origin → target); with ``symmetric=True`` (default) the reversed vector of each is added unless it is
    already listed, so pair counts and energies are symmetric."""
    vectors: tuple[Vector, ...]
    symmetric: bool = True

    def __post_init__(self) -> None:
        vecs = tuple(self.vectors)
        if not vecs or not all(isinstance(v, Vector) for v in vecs):
            raise ValueError("a neighbourhood needs at least one Vector")
        out: list[Vector] = []
        for v in vecs:
            if v.origin_site == v.target_site and v.translation == (0, 0, 0):
                raise ValueError("a vector from a site to itself in the same cell is not a pair")
            if v not in out:
                out.append(v)
        if self.symmetric:
            for v in list(out):
                if v.reversed() not in out:
                    out.append(v.reversed())
        object.__setattr__(self, "vectors", tuple(out))

    @classmethod
    def from_translations(cls, translations, site: int = 1, symmetric: bool = True) -> "Neighborhood":
        """Same-site vectors ``site → site`` for each integer translation."""
        return cls(tuple(Vector(site, site, tuple(t)) for t in translations), symmetric)


# ------------------------------------------------------------------------------------------- targets
def _one_of(correlation, energy) -> tuple[str, float, float]:
    if (correlation is None) == (energy is None):
        raise ValueError("give exactly one of correlation= (DISCUS adjusts the energy: CORR) "
                         "or energy= (fixed energy: ENER)")
    if correlation is not None:
        c = _finite(correlation, "correlation")
        if not -1.0 <= c <= 1.0:
            raise ValueError("correlation must be in [-1, 1]")
        return "CORR", c, 0.0
    return "ENER", 0.0, _finite(energy, "energy")


@dataclass(frozen=True)
class OccupationalTarget:
    """``set target, ic, corr, A, B, c, E, CORR|ENER`` — chemical short-range order between species ``a``
    and ``b`` on neighbourhood ``neighborhood`` (0-based index into ``MCConfig.neighborhoods``).

    ``correlation=c``: DISCUS's feedback adjusts the pair energy J towards the target; the energy written
    is ignored by DISCUS. ``energy=J``: fixed J in kT, energy = J·(N_like − N_unlike) per centre; J < 0
    favours like pairs (clustering), J > 0 unlike pairs (ordering)."""
    neighborhood: int
    a: str
    b: str
    correlation: float | None = None
    energy: float | None = None

    def __post_init__(self) -> None:
        _species(self.a, "a"), _species(self.b, "b")
        if self.a.upper() == self.b.upper():
            raise ValueError("an occupational target needs two different species")
        _one_of(self.correlation, self.energy)

    @property
    def mode(self) -> str:
        return _one_of(self.correlation, self.energy)[0]

    def line(self, ic: int) -> str:
        mode, c, e = _one_of(self.correlation, self.energy)
        return f"set target, {ic}, corr, {self.a}, {self.b}, {_fmt(c)}, {_fmt(e)}, {mode}"


@dataclass(frozen=True)
class DisplacementTarget:
    """``set target, ic, cd, A, B, c, E, CORR|ENER, dir:[u,v,w]`` — correlation of the displacement
    components along ``direction`` ([uvw] of the unit cell) of the two ends of each pair. DISCUS's energy
    only counts the pair when the moved atom is the vector's origin (``docs/MC_FACTS.md`` §2)."""
    neighborhood: int
    a: str
    b: str
    direction: tuple[float, float, float]
    correlation: float | None = None
    energy: float | None = None

    def __post_init__(self) -> None:
        _species(self.a, "a"), _species(self.b, "b")
        d = tuple(_finite(x, "direction") for x in self.direction)
        if len(d) != 3 or not any(d):
            raise ValueError("direction must be a nonzero [u, v, w]")
        object.__setattr__(self, "direction", d)
        _one_of(self.correlation, self.energy)

    @property
    def mode(self) -> str:
        return _one_of(self.correlation, self.energy)[0]

    def line(self, ic: int) -> str:
        mode, c, e = _one_of(self.correlation, self.energy)
        u = ",".join(_fmt(x) for x in self.direction)
        return f"set target, {ic}, cd, {self.a}, {self.b}, {_fmt(c)}, {_fmt(e)}, {mode}, dir:[{u}]"


@dataclass(frozen=True)
class SpringTarget:
    """``set target, ic, spring, A, B, d0, k`` — harmonic energy k (d − d0)² (d in Å, k in kT/Å²),
    **averaged** over the neighbours of each centre (DISCUS ``mmc_energy_spr``). No feedback."""
    neighborhood: int
    a: str
    b: str
    distance: float
    k: float

    def __post_init__(self) -> None:
        _species(self.a, "a"), _species(self.b, "b")
        if _finite(self.distance, "distance") <= 0:
            raise ValueError("distance must be positive (DISCUS ignores pairs with d0 = 0)")
        _finite(self.k, "k")

    mode = "ENER"

    def line(self, ic: int) -> str:
        return f"set target, {ic}, spring, {self.a}, {self.b}, {_fmt(self.distance)}, {_fmt(self.k)}"


@dataclass(frozen=True)
class LennardJonesTarget:
    """``set target, ic, lennard, A, B, d0, E, n_rep, n_att`` — a/d^n_rep − b/d^n_att with its minimum −|E|
    at d0 (Å), averaged over the neighbours of each centre. All four numbers are always written (with two,
    DISCUS 6.23.01 uses exponents 0 → NaN; with three it forces n_att = 12)."""
    neighborhood: int
    a: str
    b: str
    distance: float
    depth: float
    n_rep: float = 12.0
    n_att: float = 6.0

    def __post_init__(self) -> None:
        _species(self.a, "a"), _species(self.b, "b")
        if _finite(self.distance, "distance") <= 0:
            raise ValueError("distance must be positive")
        _finite(self.depth, "depth")
        if not _finite(self.n_rep, "n_rep") > _finite(self.n_att, "n_att") > 0:
            raise ValueError("need n_rep > n_att > 0")

    mode = "ENER"

    def line(self, ic: int) -> str:
        return (f"set target, {ic}, lennard, {self.a}, {self.b}, {_fmt(self.distance)}, {_fmt(self.depth)}, "
                f"{_fmt(self.n_rep)}, {_fmt(self.n_att)}")


_TARGETS = (OccupationalTarget, DisplacementTarget, SpringTarget, LennardJonesTarget)


# ------------------------------------------------------------------------------------------- moves
def _prob(p) -> float:
    v = _finite(p, "probability")
    if v <= 0:
        raise ValueError("probability must be positive")
    return v


@dataclass(frozen=True)
class Swap:
    """``set mode, p, swchem, pairing``: exchange the species of two atoms of different type (positions
    stay). ``pairing``: ``all`` (any second atom), ``local`` (±1 cell), ``slocal`` (±1 cell, same site),
    ``site`` (same site, any cell). Same-species pairs are redrawn and not counted as moves."""
    probability: float = 1.0
    pairing: str = "all"

    def __post_init__(self) -> None:
        _prob(self.probability)
        if self.pairing not in _PAIRING:
            raise ValueError(f"pairing must be one of {_PAIRING}")

    def lines(self) -> list[str]:
        return [f"set mode, {_fmt(self.probability)}, swchem, {self.pairing}"]


@dataclass(frozen=True)
class SwapDisplacement:
    """``set mode, p, swdisp, pairing``: exchange the displacements (from the site averages) of two atoms
    whose displacements point in opposite directions."""
    probability: float = 1.0
    pairing: str = "all"

    def __post_init__(self) -> None:
        _prob(self.probability)
        if self.pairing not in _PAIRING:
            raise ValueError(f"pairing must be one of {_PAIRING}")

    def lines(self) -> list[str]:
        return [f"set mode, {_fmt(self.probability)}, swdisp, {self.pairing}"]


@dataclass(frozen=True)
class Shift:
    """``set mode, p, shift, all`` + ``set move, A, sx, sy, sz``: Gaussian displacement of one atom with
    standard deviation ``sigma[A]`` per axis in **unit-cell fractions**. Every species to be moved needs a
    sigma: DISCUS 6.23.01's default is 0 (the help's 0.2 is not in the code)."""
    sigma: dict
    probability: float = 1.0

    def __post_init__(self) -> None:
        _prob(self.probability)
        if not self.sigma:
            raise ValueError("Shift needs sigma per species, e.g. {'Cu': (0.01, 0.01, 0.01)}")
        clean = {}
        for sp, s in dict(self.sigma).items():
            _species(sp, "species")
            s = tuple(_finite(x, "sigma") for x in s)
            if len(s) != 3 or any(x < 0 for x in s) or not any(s):
                raise ValueError(f"sigma for {sp!r} must be three non-negative numbers, not all zero")
            clean[sp] = s
        object.__setattr__(self, "sigma", clean)

    def lines(self) -> list[str]:
        out = [f"set mode, {_fmt(self.probability)}, shift, all"]
        out += [f"set move, {sp}, " + ", ".join(_fmt(x) for x in s) for sp, s in self.sigma.items()]
        return out


_MOVES = (Swap, SwapDisplacement, Shift)


# ------------------------------------------------------------------------------------------- configuration
@dataclass(frozen=True)
class MCConfig:
    """One DISCUS ``mmc`` run.

    ``moves_per_atom`` × number of atoms = ``set cyc`` (DISCUS counts *valid trial moves*, not cycles);
    ``feedback_per_atom`` × atoms = ``set feed`` (feedback and log interval). ``temperature`` is kT in the
    energy units of the targets (keep 1.0 when a target uses ``correlation=``: the feedback works in kT).
    Runs are always serial (``set parallel, useomp:off``): DISCUS's OpenMP MMC loop is not reproducible.
    """
    neighborhoods: tuple[Neighborhood, ...]
    targets: tuple
    moves: tuple
    seed: Seed
    moves_per_atom: float = 10.0
    feedback_per_atom: float = 1.0
    temperature: float = 1.0
    output_file: str = "mc_out.stru"

    def __post_init__(self) -> None:
        nbs = tuple(self.neighborhoods)
        if not nbs or not all(isinstance(n, Neighborhood) for n in nbs):
            raise ValueError("neighborhoods must be a non-empty tuple of Neighborhood")
        tg = tuple(self.targets)
        if not tg or not all(isinstance(t, _TARGETS) for t in tg):
            raise ValueError(f"targets must be a non-empty tuple of {[t.__name__ for t in _TARGETS]}")
        for t in tg:
            if isinstance(t.neighborhood, bool) or not isinstance(t.neighborhood, (int, np.integer)) \
                    or not 0 <= t.neighborhood < len(nbs):
                raise ValueError(f"target neighborhood index {t.neighborhood!r} out of range 0..{len(nbs) - 1}")
        for kind, why in ((OccupationalTarget, "DISCUS keeps only the last occupational target of a neighbourhood"),
                          (DisplacementTarget, "DISCUS keeps one displacement direction per neighbourhood")):
            used = [t.neighborhood for t in tg if isinstance(t, kind)]
            if len(used) != len(set(used)):
                raise ValueError(f"two {kind.__name__}s on one neighbourhood: {why}; use separate neighbourhoods")
        n_feedback = sum(1 for t in tg if t.mode == "CORR")
        if n_feedback > 1:
            raise ValueError("at most one target may use correlation= (feedback); DISCUS runs one independent "
                             "feedback loop per target and recommends fixed energies when combining targets")
        mv = tuple(self.moves)
        if not mv or not all(isinstance(m, _MOVES) for m in mv):
            raise ValueError("moves must be a non-empty tuple of Swap / SwapDisplacement / Shift")
        if len({type(m) for m in mv}) != len(mv):
            raise ValueError("give each move type once")
        if not isinstance(self.seed, Seed):
            raise ValueError("seed must be an 8-value Seed: a run without one cannot be repeated")
        for name in ("moves_per_atom", "feedback_per_atom"):
            if _finite(getattr(self, name), name) <= 0:
                raise ValueError(f"{name} must be positive")
        if self.feedback_per_atom > self.moves_per_atom:
            raise ValueError("feedback_per_atom must not exceed moves_per_atom")
        if _finite(self.temperature, "temperature") <= 0:
            raise ValueError("temperature (kT) must be positive")
        if not (isinstance(self.output_file, str) and Path(self.output_file).name == self.output_file
                and self.output_file.endswith(".stru") and " " not in self.output_file):
            raise ValueError("output_file must be a plain '.stru' file name without spaces")
        object.__setattr__(self, "neighborhoods", nbs)
        object.__setattr__(self, "targets", tg)
        object.__setattr__(self, "moves", mv)

    # -- facts ---------------------------------------------------------------------------------------
    def cautions(self) -> list[str]:
        """Known DISCUS behaviours that affect this configuration (stored with every run)."""
        out = ["acceptance is heat-bath for uphill and neutral moves (ΔE = 0 accepted with probability ½), "
               "not Metropolis"]
        for t in self.targets:
            if t.mode == "CORR" and isinstance(t, OccupationalTarget):
                c = t.correlation
                if abs(c) < 0.001:
                    out.append(f"target {c:g}: DISCUS fixes J = 0 for |target| < 0.001 (random swaps)")
                elif abs(c) < 0.1:
                    out.append(f"target {c:g}: DISCUS scales the feedback by 0.005 for |target| < 0.1; "
                               "convergence is very slow")
                elif abs(c) < 0.2:
                    out.append(f"target {c:g}: DISCUS scales the feedback by 0.05 for |target| < 0.2")
            if isinstance(t, DisplacementTarget):
                out.append("cd energy counts a pair only when the moved atom is the vector origin")
            if isinstance(t, (SpringTarget, LennardJonesTarget)):
                out.append(f"{type(t).__name__}: energy is averaged over each centre's neighbours")
        if self.temperature != 1.0 and any(t.mode == "CORR" for t in self.targets):
            out.append("feedback targets work in units of kT; temperature != 1 rescales them")
        return out

    def counts(self, n_atoms: int) -> tuple[int, int]:
        """(``set cyc``, ``set feed``) for a crystal of ``n_atoms``."""
        if not isinstance(n_atoms, (int, np.integer)) or n_atoms < 1:
            raise ValueError("n_atoms must be a positive integer")
        cyc = max(1, int(round(self.moves_per_atom * n_atoms)))
        feed = max(1, int(round(self.feedback_per_atom * n_atoms)))
        return cyc, min(feed, cyc)

    # -- rendering -----------------------------------------------------------------------------------
    def mmc_lines(self, n_atoms: int) -> list[str]:
        cyc, feed = self.counts(n_atoms)
        vec_lines, neig_lines, iv = [], [], 0
        for k, nb in enumerate(self.neighborhoods, start=1):
            ids = []
            for v in nb.vectors:
                iv += 1
                ids.append(iv)
                t = v.translation
                vec_lines.append(f"set vec, {iv}, {v.origin_site}, {v.target_site}, {t[0]}, {t[1]}, {t[2]}")
            neig_lines.append("set neig, vec, " + ", ".join(str(i) for i in ids) + f", number:{k}")
        lines = ["set vec, reset", *vec_lines, "set neig, rese", *neig_lines]
        for m in self.moves:
            lines += m.lines()
        for t in self.targets:
            lines.append(t.line(t.neighborhood + 1))
        lines += [f"set temp, {_fmt(self.temperature)}", f"set cyc, {cyc}", f"set feed, {feed}",
                  "set output, feed:on, final:on"]
        return lines

    def render(self, structure_file: str, n_atoms: int, *, title: str | None = None) -> str:
        """Complete macro text (LF): read the structure, run mmc, save ``output_file`` from the main level."""
        if Path(structure_file).name != structure_file or not structure_file.endswith(".stru") \
                or " " in structure_file:
            raise ValueError("structure_file must be a plain '.stru' file name without spaces")
        if structure_file == self.output_file:
            raise ValueError("structure_file and output_file must differ")
        head = [f"# {title}"] if title else []
        if title and "\n" in title:
            raise ValueError("title must be one line")
        body = head + ["discus", self.seed.command(), "set parallel, useomp:off", "read",
                       f"  stru {structure_file}", "mmc"]
        body += self.mmc_lines(n_atoms) + ["show", "run", "exit", f"save {self.output_file}", "exit", "exit", ""]
        return "\n".join(body)

    def check_structure(self, structure: Structure) -> None:
        """Refuse structures this configuration cannot describe: sites out of range, unknown species."""
        n_cells = int(np.prod(structure.supercell))
        if structure.n_atoms % n_cells:
            raise ValueError("the structure does not have the same number of atoms in every cell")
        m = structure.n_atoms // n_cells
        for nb in self.neighborhoods:
            for v in nb.vectors:
                if max(v.origin_site, v.target_site) > m:
                    raise ValueError(f"vector uses site {max(v.origin_site, v.target_site)} but a cell has {m} sites")
        present = {s.upper() for s in structure.species}
        for t in self.targets:
            for sp in (t.a, t.b):
                if sp.upper() not in present:
                    raise ValueError(f"target species {sp!r} is not in the structure")
        for mv in self.moves:
            if isinstance(mv, Shift):
                for sp in mv.sigma:
                    if sp.upper() not in present:
                        raise ValueError(f"Shift sigma given for {sp!r}, which is not in the structure")
        if any(n == 1 for n in structure.supercell):
            raise ValueError("DISCUS switches off periodic boundaries along an axis with one cell; "
                             "use at least 2 cells per axis")

    # -- persistence ---------------------------------------------------------------------------------
    def to_dict(self) -> dict[str, Any]:
        def target(t):
            d = {"type": type(t).__name__, **{k: getattr(t, k) for k in t.__dataclass_fields__}}
            if "direction" in d:
                d["direction"] = list(d["direction"])
            return d

        def move(m):
            d = {"type": type(m).__name__, **{k: getattr(m, k) for k in m.__dataclass_fields__}}
            if "sigma" in d:
                d["sigma"] = {k: list(v) for k, v in d["sigma"].items()}
            return d

        return {
            "schema": "discuspdf.MCConfig", "version": 1,
            "neighborhoods": [{"vectors": [[v.origin_site, v.target_site, list(v.translation)] for v in nb.vectors],
                               "symmetric": nb.symmetric} for nb in self.neighborhoods],
            "targets": [target(t) for t in self.targets],
            "moves": [move(m) for m in self.moves],
            "seed": list(self.seed.values), "moves_per_atom": self.moves_per_atom,
            "feedback_per_atom": self.feedback_per_atom, "temperature": self.temperature,
            "output_file": self.output_file,
        }

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> "MCConfig":
        if d.get("schema") != "discuspdf.MCConfig" or d.get("version") != 1:
            raise ValueError("not a discuspdf.MCConfig version 1 dictionary")
        types = {c.__name__: c for c in (*_TARGETS, *_MOVES)}
        nbs = tuple(Neighborhood(tuple(Vector(a, b, tuple(t)) for a, b, t in nb["vectors"]), nb["symmetric"])
                    for nb in d["neighborhoods"])

        def build(item):
            item = dict(item)
            c = types[item.pop("type")]
            if "direction" in item:
                item["direction"] = tuple(item["direction"])
            if "sigma" in item:
                item["sigma"] = {k: tuple(v) for k, v in item["sigma"].items()}
            return c(**item)

        return cls(nbs, tuple(build(t) for t in d["targets"]), tuple(build(m) for m in d["moves"]),
                   Seed(tuple(d["seed"])), d["moves_per_atom"], d["feedback_per_atom"], d["temperature"],
                   d["output_file"])


# ------------------------------------------------------------------------------------------- log parsing
_GEN = re.compile(r"Gen:\s*(\d+)\s+try:\s*(\d+)\s+acc:\s*\(g/n/b\):\s*(\d+)\s*/\s*(\d+)\s*/\s*(\d+)\s+MC moves\s*([0-9.]*)")
_NUM = r"([-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?|\*+)"
_ROW_CORR = re.compile(r"^\s*(\d+)\s+(Occupancy|Disp\.Cor\.)\s*(\S+)\s+(\S+)\s+" + r"\s+".join([_NUM] * 4) + r"\s+(\d+|\*+)\s*$")
_ROW_DIST = re.compile(r"^\s*(\d+)\s+(Hooke|Lennard)\s+(\S+)\s+(\S+)\s+" + r"\s+".join([_NUM] * 5) + r"\s+(\d+|\*+)\s*$")


def _val(s: str) -> float:
    return float("nan") if s.startswith("*") else float(s)


def _count(s: str) -> int | None:
    return None if s.startswith("*") else int(s)


_KIND = {OccupationalTarget: "occupational", DisplacementTarget: "displacement", SpringTarget: "spring",
         LennardJonesTarget: "lennard_jones"}


def match_row(rows: list[dict], target, *, required: bool = True) -> dict | None:
    """The row of ``target`` in one printed table: same neighbourhood, energy type and species pair.

    DISCUS prints rows by neighbourhood, then energy type, then atom-type index — not in the order the
    targets were given — and skips rows without pairs; so rows are matched, never taken by position."""
    want = (target.neighborhood + 1, _KIND[type(target)], frozenset({target.a.upper(), target.b.upper()}))
    hits = [r for r in rows if (r["neighborhood"], r["kind"], frozenset({r["a"].upper(), r["b"].upper()})) == want]
    if len(hits) > 1:
        raise ValueError(f"{len(hits)} rows match {type(target).__name__} on neighbourhood {target.neighborhood}")
    if not hits:
        if required:
            raise ValueError(f"no printed row for {type(target).__name__} {target.a}-{target.b} on neighbourhood "
                             f"{target.neighborhood} (DISCUS skips a target without pairs)")
        return None
    return hits[0]


def _row(line: str) -> dict | None:
    m = _ROW_CORR.match(line)
    if m:
        ic, kind, a, b, tgt, ach, diff, rel, n = m.groups()
        return {"neighborhood": int(ic), "kind": "occupational" if kind == "Occupancy" else "displacement",
                "a": a, "b": b, "target": _val(tgt), "achieved": _val(ach), "sigma": None,
                "diff": _val(diff), "rel_diff": _val(rel), "pairs": _count(n)}
    m = _ROW_DIST.match(line)
    if m:
        ic, kind, a, b, tgt, ach, sig, diff, rel, n = m.groups()
        return {"neighborhood": int(ic), "kind": "spring" if kind == "Hooke" else "lennard_jones",
                "a": a, "b": b, "target": _val(tgt), "achieved": _val(ach), "sigma": _val(sig),
                "diff": _val(diff), "rel_diff": _val(rel), "pairs": _count(n)}
    return None


def parse_mmc_log(text: str) -> dict:
    """Split a DISCUS 6.23.01 ``mmc run`` log into ``initial`` rows, ``feedback`` records (``gen``, ``tries``,
    ``accepted`` = (good, neutral, bad), ``fraction``, ``rows``), the ``final`` record and ``energies``.

    Formats: ``docs/MC_FACTS.md`` §4. The "Initial" table that DISCUS prints again at the end is ignored."""
    text = re.sub(r"\x1b\[[0-9;]*m", "", text)
    lines = text.splitlines()
    out: dict = {"initial": [], "feedback": [], "final": None, "energies": [], "errors": []}
    state, current = None, None
    for line in lines:
        if re.search(r"\*\*\*[A-Z ]{3,4}\*\*\*", line) or "did not find a valid pair" in line:
            out["errors"].append(line.strip())
        if "--- Initial multiple energy configuration ---" in line:
            state = "initial" if not out["initial"] and not out["feedback"] else "reprint"
            continue
        if "--- Final multiple energy configuration ---" in line:
            state = "final"
            continue
        if "--- Average energy changes ---" in line:
            state = "energies"
            continue
        m = _GEN.search(line)
        if m:
            g, tr, a1, a2, a3, frac = m.groups()
            current = {"gen": int(g), "tries": int(tr), "accepted": (int(a1), int(a2), int(a3)),
                       "fraction": float(frac) if frac else None, "rows": []}
            if state == "final":
                out["final"] = current
            else:
                out["feedback"].append(current)
                state = "feedback"
            continue
        r = _row(line)
        if r is not None:
            if state == "initial":
                out["initial"].append(r)
            elif state in ("feedback", "final") and current is not None:
                current["rows"].append(r)
            continue
        if state == "energies":
            m = re.match(r"^\s*(.+?)\s{2,}(\d+)\s+" + _NUM + r"\s+(\d+)\s+(\d+)\s+" + _NUM + r"\s*$", line)
            if m:
                name, n_down, e_down, n_zero, n_up, e_up = m.groups()
                out["energies"].append({"energy": name.strip(), "n_down": int(n_down), "mean_down": _val(e_down),
                                        "n_zero": int(n_zero), "n_up": int(n_up), "mean_up": _val(e_up)})
    return out


# ------------------------------------------------------------------------------------------- results
@dataclass
class MCResult:
    """Parsed outcome of one run: ``structure`` (final configuration, sites labelled by slot), the log
    tables and the configuration. ``status`` is ``complete`` when DISCUS finished all moves without an
    error; it says nothing about whether the targets were reached (see :meth:`target_report`)."""
    directory: Path
    config: MCConfig
    structure: Structure
    log: dict
    status: str
    n_atoms: int

    def target_report(self, tol: float = 0.02) -> list[dict]:
        """Target vs achieved (final table) for every feedback target, with ``reached`` = |diff| ≤ tol.
        Fixed-energy targets have no target value to reach; they are reported with ``reached=None``."""
        final = self.log["final"]
        if final is None:
            raise ValueError("the log has no final table")
        rows = []
        for t in self.config.targets:
            r = match_row(final["rows"], t)
            rec = dict(r, mode=t.mode)
            rec["reached"] = bool(abs(r["achieved"] - r["target"]) <= tol) if t.mode == "CORR" else None
            rows.append(rec)
        return rows

    def acceptance(self) -> np.ndarray:
        """(n_feedback, 5): gen, tries, good, neutral, bad (cumulative counts)."""
        fb = self.log["feedback"]
        return np.array([[f["gen"], f["tries"], *f["accepted"]] for f in fb], dtype=np.int64).reshape(-1, 5)

    def measure(self, reference: Structure | None = None) -> list[dict]:
        """Independent measurement of every target on the final structure with Package C.

        Pairs are exactly the neighbourhood's directed vectors (as DISCUS counts them); ``reference`` gives
        the undisplaced site positions (default: the per-site average cell positions of the final structure,
        as DISCUS's ``chem_aver``). Each value follows DISCUS's own definition (``docs/MC_FACTS.md`` §3), computed
        independently: occupational — the mean of the two Warren–Cowley α with centre-based compositions
        (``pooled_alpha`` gives Package C's pooled-θ α as well; the two agree for a binary system with
        symmetric pair counts); displacement — uncentred ⟨x_i x_j⟩/√(⟨x_i²⟩⟨x_j²⟩) of the components along the
        direction, with (a→b) and (b→a) pairs pooled; spring / Lennard-Jones — mean and standard deviation of
        the a–b pair distances (Å) and their number."""
        from .correlations import displacements, node_species, occupational_correlation
        from .neighbors import neighbor_list
        st = self.structure
        ref = _site_average_reference(st) if reference is None else reference
        out = []
        for t in self.config.targets:
            nb = self.config.neighborhoods[t.neighborhood]
            pairs = _neighborhood_pairs(neighbor_list, ref, nb)
            a, b = _spelling(st, t.a), _spelling(st, t.b)
            spec = node_species(pairs, st)
            so, sj = spec[pairs.i], spec[pairs.j]
            rec = {"type": type(t).__name__, "neighborhood": t.neighborhood, "a": t.a, "b": t.b,
                   "pairs": pairs.n_pairs}
            if isinstance(t, OccupationalTarget):
                rec["value"] = _discus_occupational(so, sj, a, b)
                rec["pooled_alpha"] = float(occupational_correlation(pairs, st).discus(a, b))
            elif isinstance(t, DisplacementTarget):
                u = displacements(pairs, st, reference="site_average")
                d = np.asarray(t.direction, float) @ (st.lattice.matrix / st.supercell[:, None])
                x = u @ (d / np.linalg.norm(d))
                rec["value"] = _discus_displacement(x[pairs.i], x[pairs.j], so, sj, a, b)
            else:
                dist = _pair_distances(pairs, ref, st, a, b)
                rec["value"] = float(dist.mean()) if len(dist) else float("nan")
                rec["std"] = float(dist.std()) if len(dist) else float("nan")
                rec["pairs"] = int(len(dist))
            out.append(rec)
        return out

    def history(self, target_index: int = 0) -> np.ndarray:
        """(n_feedback, 3): gen, target, achieved of ``config.targets[target_index]`` (feedback lines)."""
        t = self.config.targets[target_index]
        rows = []
        for f in self.log["feedback"]:
            r = match_row(f["rows"], t, required=False)
            if r is not None:
                rows.append((f["gen"], r["target"], r["achieved"]))
        return np.array(rows, dtype=float).reshape(-1, 3)


def _discus_occupational(so, sj, a: str, b: str) -> float:
    """DISCUS ``mmc_correlations_occ``: groups 1 = a, 2 = b, 3 = other; mean of the two α."""
    g = lambda s: np.where(s == a, 0, np.where(s == b, 1, 2))   # noqa: E731
    p = np.zeros((3, 3))
    np.add.at(p, (g(so), g(sj)), 1.0)
    r = p.sum()
    if r == 0 or p[0].sum() == 0 or p[1].sum() == 0:
        return 1.0                                             # DISCUS's value for an empty group
    m_left, m_right = p[0].sum() / r, p[1].sum() / r
    return float(0.5 * ((1 - p[1, 0] / p[1].sum() / m_left) + (1 - p[0, 1] / p[0].sum() / m_right)))


def _discus_displacement(xi, xj, so, sj, a: str, b: str) -> float:
    """DISCUS cd achieved: Σ x_i x_j / √(Σ x_i² Σ x_j²) over (a→b) pairs, plus (b→a) pairs when a ≠ b."""
    keep = (so == a) & (sj == b)
    if a != b:
        keep |= (so == b) & (sj == a)
    xi, xj = xi[keep], xj[keep]
    den = np.sqrt(np.sum(xi * xi) * np.sum(xj * xj))
    return float(np.sum(xi * xj) / den) if den > 0 else float("nan")


def _site_average_reference(st: Structure) -> Structure:
    """Undisplaced reference: cell origin + per-site average position within the cell."""
    n = st.supercell
    local = st.frac * n - st.cell_index
    avg = {lab: local[st.site_label == lab].mean(axis=0) for lab in set(st.site_label)}
    frac = (st.cell_index + np.array([avg[lab] for lab in st.site_label])) / n
    return st.replace(frac=frac, min_separation=None)


def _neighborhood_pairs(neighbor_list, ref: Structure, nb: Neighborhood):
    lengths = [np.linalg.norm((np.array(v.translation, float)) @ (ref.lattice.matrix / ref.supercell[:, None]))
               for v in nb.vectors]
    cell = ref.lattice.matrix / ref.supercell[:, None]
    cutoff = max(lengths) + float(np.linalg.norm(cell.sum(axis=0))) + float(np.linalg.norm(cell[0] - cell[1] + cell[2]))
    nl = neighbor_list(ref, cutoff)
    mask = np.zeros(nl.n_pairs, bool)
    for v in nb.vectors:
        mask |= ((nl.origin_site == f"site{v.origin_site}") & (nl.target_site == f"site{v.target_site}")
                 & np.all(nl.translation == np.array(v.translation), axis=1))
    return nl._subset(mask)


def _spelling(st: Structure, name: str) -> str:
    """The structure's own spelling of a species name (DISCUS names are case-insensitive)."""
    for s in set(st.species):
        if s.upper() == name.upper():
            return s
    raise ValueError(f"species {name!r} is not in the structure")


def _pair_distances(pairs, ref: Structure, st: Structure, a: str, b: str) -> np.ndarray:
    """Distances (Å) of the a–b pairs in ``st``; the pair image comes from ``ref`` (minimum-image
    correction of each end's displacement, so atoms that crossed a boundary are handled)."""
    from .correlations import node_species
    spec = node_species(pairs, st)
    rows = pairs.node_row
    so, sj = spec[pairs.i], spec[pairs.j]
    keep = ((so == a) & (sj == b)) | ((so == b) & (sj == a))
    io, it = rows[pairs.i[keep]], rows[pairs.j[keep]]
    d_ref = ref.frac[it] + pairs.shift[keep] - ref.frac[io]
    du = st.lattice.minimum_image(st.frac[it] - ref.frac[it])[0] - st.lattice.minimum_image(st.frac[io] - ref.frac[io])[0]
    return np.linalg.norm((d_ref + du) @ st.lattice.matrix, axis=1)


def replica_summary(results: list["MCResult"]) -> list[dict]:
    """Seed-to-seed spread: mean and standard deviation over replicas of DISCUS's final achieved value and
    of the independent measurement, per target. All results must come from the same configuration
    apart from the seed."""
    if len(results) < 2:
        raise ValueError("give at least two replicas")
    base = results[0].config.to_dict()
    for r in results[1:]:
        d = r.config.to_dict()
        if {**d, "seed": None} != {**base, "seed": None}:
            raise ValueError("replicas differ in more than the seed")
    if len({tuple(r.config.seed.values) for r in results}) != len(results):
        raise ValueError("replicas must use different seeds")
    out = []
    for k, t in enumerate(results[0].config.targets):
        ach = np.array([match_row(r.log["final"]["rows"], t)["achieved"] for r in results], float)
        mea = np.array([r.measure()[k]["value"] for r in results], float)
        out.append({"target_index": k, "type": type(t).__name__, "mode": t.mode,
                    "target": getattr(t, "correlation", None), "n": len(results),
                    "achieved_mean": float(ach.mean()), "achieved_std": float(ach.std(ddof=1)),
                    "measured_mean": float(mea.mean()), "measured_std": float(mea.std(ddof=1))})
    return out


# ------------------------------------------------------------------------------------------- runner
@dataclass(frozen=True)
class MCRun:
    directory: Path
    config: MCConfig
    n_atoms: int


_ENGINE_ERROR = re.compile(r"(?im)\b(?:unknown command|fatal error|segmentation fault)"
                           r"|^(?!.*\bwarning\b).*\*\*\*(?:[A-Z]{3,4} ?)\*\*\*.*$|did not find a valid pair")


class MCRunner:
    """Stage and execute :class:`MCConfig` runs, one fresh directory per run.

    ``launcher`` is the command that runs ``discus_suite`` with a macro in the run directory, e.g. the
    default ``("discus_suite", "-noautorun", "-macro")`` inside WSL/Linux."""

    def __init__(self, runs_directory="runs_mc", launcher=("discus_suite", "-noautorun", "-macro"),
                 timeout_seconds: float = 3600, *, engine: Engine | None = None):
        """``engine`` (a :class:`discuspdf.engine.Engine`) replaces ``launcher`` and ``timeout_seconds``."""
        self.runs_directory = Path(runs_directory)
        launcher = tuple(launcher)
        if engine is None and (not launcher or timeout_seconds <= 0):
            raise ValueError("invalid launcher or timeout")
        self.engine = as_engine(engine, launcher, timeout_seconds)
        self.launcher = self.engine.launcher
        self.timeout_seconds = self.engine.timeout_seconds

    def prepare(self, config: MCConfig, structure: Structure, *, title: str | None = None) -> MCRun:
        from .discus_io import format_stru
        if not isinstance(config, MCConfig):
            raise TypeError("config must be an MCConfig")
        config.check_structure(structure)
        text_stru = format_stru(structure, title=title or "discuspdf MC input")
        macro = config.render("mc_in.stru", structure.n_atoms, title=title)
        self.runs_directory.mkdir(parents=True, exist_ok=True)
        directory = self.runs_directory / uuid.uuid4().hex
        directory.mkdir()
        files = {"mc_in.stru": text_stru, "run.mac": macro,
                 "mc_config.json": json.dumps({**config.to_dict(), "cautions": config.cautions(),
                                               "n_atoms": structure.n_atoms,
                                               "supercell": [int(x) for x in structure.supercell],
                                               "species": sorted(set(structure.species)),
                                               "counts": config.counts(structure.n_atoms)}, indent=1) + "\n"}
        hashes = {}
        for name, text in files.items():
            (directory / name).write_bytes(text.encode("utf-8"))
            hashes[name] = hashlib.sha256(text.encode("utf-8")).hexdigest()
        (directory / "inputs.json").write_text(json.dumps(hashes, indent=1) + "\n")
        (directory / "status.json").write_text(json.dumps({"status": "prepared"}) + "\n")
        return MCRun(directory, config, structure.n_atoms)

    def _status(self, run: MCRun, state: str, **kw) -> None:
        (run.directory / "status.json").write_text(json.dumps({"status": state, **kw}, indent=1) + "\n")

    def execute(self, run: MCRun) -> MCResult:
        """Run DISCUS on the staged macro; raise ``RuntimeError`` (status ``failed``) on any engine error,
        missing output or incomplete move count."""
        hashes = json.loads((run.directory / "inputs.json").read_text())
        for name, digest in hashes.items():
            if hashlib.sha256((run.directory / name).read_bytes()).hexdigest() != digest:
                self._status(run, "failed", reason=f"staged {name} was modified", engine_started=False)
                raise RuntimeError(f"staged {name} was modified after preparation")
        out = run.directory / run.config.output_file
        if out.exists():
            raise FileExistsError(f"{out} exists; prepare a new run")
        try:
            engine = current_engine(self)
        except ValueError as exc:
            self._status(run, "failed", reason=f"invalid engine settings: {exc}", engine_started=False)
            raise
        self._status(run, "running")
        try:
            proc = engine.run(run.directory, "run.mac", error_pattern=None)  # load() checks the log
        except (FileNotFoundError, subprocess.TimeoutExpired) as exc:
            self._status(run, "failed", reason=str(exc))
            raise RuntimeError(f"DISCUS could not complete; see {run.directory}") from exc
        return self.load(run, returncode=proc.returncode)

    def load(self, run: MCRun, *, returncode: int | None = None) -> MCResult:
        """Validate and load an executed run directory. ``returncode`` defaults to the one recorded in
        ``status.json`` by :meth:`execute`; a run whose exit code is unknown is not loaded."""
        from .discus_io import read_stru
        from .neighbors import label_sites_by_slot
        problems = []
        try:
            hashes = json.loads((run.directory / "inputs.json").read_text())
            for name, digest in hashes.items():
                if hashlib.sha256((run.directory / name).read_bytes()).hexdigest() != digest:
                    problems.append(f"staged {name} was modified")
            staged = json.loads((run.directory / "mc_config.json").read_text())
        except (OSError, json.JSONDecodeError) as exc:
            problems.append(f"run directory incomplete: {exc}")
            staged = None
        if returncode is None:
            try:
                returncode = json.loads((run.directory / "status.json").read_text()).get("returncode")
            except (OSError, json.JSONDecodeError):
                returncode = None
        if returncode is None:
            problems.append("exit code unknown")
        elif returncode != 0:
            problems.append(f"exit code {returncode}")
        log_path = run.directory / "run.log"
        text = log_path.read_text(encoding="utf-8", errors="replace") if log_path.is_file() else ""
        if not text:
            problems.append("no run.log")
        log = parse_mmc_log(text)
        err = _ENGINE_ERROR.search(re.sub(r"\x1b\[[0-9;]*m", "", text))
        if err:
            problems.append(f"engine error: {err.group(0).strip()}")
        cyc = staged["counts"][0] if staged else None
        if log["final"] is None:
            problems.append("no final table in the log")
        elif cyc is not None and log["final"]["gen"] != cyc:
            problems.append(f"{log['final']['gen']} moves done, {cyc} requested")
        out = run.directory / run.config.output_file
        if not out.is_file():
            problems.append(f"{run.config.output_file} was not written")
        structure = None
        if not problems:
            try:
                # DISCUS 'save' writes no ncell line (the crystal size is the input's) and upper-cases names
                raw = read_stru(out, supercell=tuple(staged["supercell"])).structure
                names = {s.upper(): s for s in staged["species"]}
                unknown = sorted({s for s in raw.species if s.upper() not in names})
                if unknown:
                    raise ValueError(f"species {unknown} were not in the input")
                if raw.n_atoms != staged["n_atoms"]:
                    raise ValueError(f"{raw.n_atoms} atoms in the output, {staged['n_atoms']} in the input")
                raw = raw.replace(species=np.array([names[s.upper()] for s in raw.species], dtype=object),
                                  min_separation=None)
                structure = label_sites_by_slot(raw)
                for t in run.config.targets:
                    match_row(log["final"]["rows"], t)
            except Exception as exc:                    # any reading problem is a failed run, recorded
                problems.append(f"output not usable: {exc}")
        if problems:
            self._status(run, "failed", returncode=returncode, reasons=problems)
            raise RuntimeError("DISCUS MC run failed: " + "; ".join(problems) + f" (see {run.directory})")
        self._status(run, "complete", returncode=returncode, moves=cyc,
                     output_sha256=hashlib.sha256(out.read_bytes()).hexdigest())
        return MCResult(run.directory, run.config, structure, log, "complete", run.n_atoms)
