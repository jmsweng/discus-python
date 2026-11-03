"""Reverse Monte Carlo fits to single-crystal diffuse scattering planes (the DISCUS ``rmc`` menu).

:class:`RmcFit` is a :class:`~discuspdf.discus_ops.Step`: it reads one or more planes of intensity
(NIPL files), moves or swaps atoms, and keeps a move with DISCUS's acceptance rule, so the crystal
saved by the :class:`~discuspdf.discus_ops.StructureJob` is the fitted model. :class:`FourierPlane`
computes a plane from a model in the same geometry, which is how a test target is made, and
:func:`read_nipl` / :func:`write_nipl` move planes between numpy and DISCUS::

    from discuspdf.discus_rmc import FourierPlane, PlaneGrid, RmcFit, RmcPlane, RmcSwap, run_rmc

    grid = PlaneGrid(ll=(0, -1, 0), lr=(3, -1, 0), ul=(0, 1, 0), nx=31, ny=21)
    fit = RmcFit(planes=(RmcPlane("target.nipl", grid),), moves=(RmcSwap(),), select=("CU", "AU"),
                 tries=2000, sigma=0.0)
    outcome = run_rmc(StructureJob(source, [fit], seed=seed), runner, {..., "target.nipl": path})
    outcome.log.final.chi2, outcome.result.structure

Facts from the 6.23.01 source (``docs/DISCUS_RMC_FACTS.md`` has the details):

* ``data`` (``rmc_readdata``, ``rmc_sup_mod.f90:1224-1670``) takes the radiation, the wavelength
  (a symbol or Å), the file, the weighting scheme and the hkl of three corners: lower left, lower
  right, upper left. Point (i, j) of the file is ``ll + (i-1)(lr-ll)/(nx-1) + (j-1)(ul-ll)/(ny-1)``,
  the same as the ``fourier`` menu's ``ll/lr/ul/na/no`` and its NIPL output, rows j = 1..ny of nx
  values. The ranges in the file header are read but not used. Only NIPL data are offered: the
  PGM reader is commented out. Weight files are not offered: their size check compares two
  variables that are never set (``:1420``).
* Cost: per plane chi² = Σ w (I_obs − s·I_calc − b)² / Σ w with I_calc = |F|² (not normalised),
  summed over planes (``:1925-1941``). The scale s and background b are fitted in closed form when
  switched on, otherwise fixed (``rmc_calcskal``, ``:2544-2591``). DISCUS prints and stores 1/s.
* Acceptance (``:2029-2052``): downhill always; uphill with probability exp(−Δχ²/(σ²/2)),
  ``σ = set sigm`` (default 0.01); σ = 0 accepts downhill moves only.
* ``set cycl`` counts valid tries (moves that passed the type and distance checks), not cycles per
  atom; the log line is printed every ``set disp`` generated moves (``:2081-2100``).
* **Defaults that break a fit:** the move probabilities start at 0 for every mode, and with no
  ``set mode`` the move generator falls through to a mode with no branch (``:3300-3378``); the
  atom selection is cleared by ``reset`` (``rmc.f90:274``). The step therefore always writes
  ``set mode`` for each move and ``dese all`` then ``sele``.
* **Scale constraint:** ``rmc_constrain`` is allocated as 0 for every plane and only ``set const``
  changes it, but ``rmc_calcskal`` replaces the scale and background of plane p by those of plane
  ``rmc_constrain(p)`` when it differs from p (``:2581-2584``) — index 0, outside the array. The
  step writes ``set const, p, p`` for every plane, which makes that branch inactive.
* Random numbers come from the Fortran ``RANDOM_NUMBER`` stream that the ``seed`` command sets.
* Lots and ``set aver`` use code the fork corrects (lot index, lot membership, the row loop of
  the average subtraction); :func:`discuspdf.engine.fork_reasons` routes them to the fork.
"""
from __future__ import annotations

import math
import numbers
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np

from .discus_ops import JobResult, Step, StructureJob
from .fourier_macro import _RADIATION_KEYWORD, RADIATIONS, _check_wavelength
from .macro import Macro

__all__ = [
    "PAIRINGS", "WEIGHTS", "PlaneGrid", "RmcPlane", "RmcShift", "RmcSwap", "RmcSwapDisplacement", "RmcLots",
    "RmcFit", "FourierPlane", "RmcProgress", "RmcLog", "RmcOutcome", "parse_rmc_log", "run_rmc",
    "read_nipl", "write_nipl",
]

# rmc_set_mode (rmc_sup_mod.f90:606-622): A, L, SL, SI (CO, connectivity, is not offered)
PAIRINGS = ("all", "local", "slocal", "site")
# rmc_readdata (rmc_sup_mod.f90:1369-1392): the first four characters, upper case
WEIGHTS = {"one": "one", "sqrt": "sqrt", "log": "log", "linear": "lin", "square": "squa",
           "inverse": "inv", "inverse_sqrt": "isq"}
_NAME = re.compile(r"^[A-Za-z0-9_.\-]+$")


def _real(name: str, v) -> float:
    if isinstance(v, bool) or not isinstance(v, numbers.Real) or not math.isfinite(float(v)):
        raise ValueError(f"{name} must be a finite number, got {v!r}")
    return float(v)


def _int(name: str, v, minimum: int = 1) -> int:
    if isinstance(v, bool) or not isinstance(v, numbers.Integral) or v < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}, got {v!r}")
    return int(v)


def _vec3(name: str, v) -> tuple[float, float, float]:
    t = tuple(v)
    if len(t) != 3:
        raise ValueError(f"{name} must have three components, got {v!r}")
    return tuple(_real(name, x) for x in t)


def _atom(name: str, v) -> str:
    if not isinstance(v, str) or not v.strip() or any(c in v for c in " ,:#!()[]"):
        raise ValueError(f"{name} must be a DISCUS atom name, got {v!r}")
    return v.strip().upper()


def _file(name: str, v) -> str:
    if not isinstance(v, str) or not _NAME.match(v) or v in (".", ".."):
        raise ValueError(f"{name} must be a plain file name (letters, digits, '_', '-', '.'), got {v!r}")
    return v


# ----------------------------------------------------------------------------------------- geometry
@dataclass(frozen=True)
class PlaneGrid:
    """A plane of nx × ny points in reciprocal space, given by the hkl of three corners: ``ll``
    (first point of the first row), ``lr`` (last point of the first row) and ``ul`` (first point of
    the last row). The same numbers go to ``fourier`` (``ll/lr/ul/na/no``) and to ``rmc data``."""
    ll: tuple
    lr: tuple
    ul: tuple
    nx: int
    ny: int

    def __post_init__(self) -> None:
        for k in ("ll", "lr", "ul"):
            object.__setattr__(self, k, _vec3(k, getattr(self, k)))
        _int("nx", self.nx, 2), _int("ny", self.ny, 2)
        a = np.subtract(self.lr, self.ll)
        b = np.subtract(self.ul, self.ll)
        if not np.any(a) or not np.any(b) or np.linalg.norm(np.cross(a, b)) < 1e-12 * (
                np.linalg.norm(a) * np.linalg.norm(b)):
            raise ValueError("the corners must span a plane (lr - ll and ul - ll not parallel or zero)")

    def hkl(self) -> np.ndarray:
        """hkl of every point, shape (ny, nx, 3): ``hkl()[j, i]`` is file row j, column i."""
        a = (np.subtract(self.lr, self.ll)) / (self.nx - 1)
        b = (np.subtract(self.ul, self.ll)) / (self.ny - 1)
        i = np.arange(self.nx)[None, :, None]
        j = np.arange(self.ny)[:, None, None]
        return np.asarray(self.ll)[None, None, :] + i * a + j * b

    def corner_values(self) -> list[float]:
        return [*self.ll, *self.lr, *self.ul]


# ----------------------------------------------------------------------------------------- NIPL files
def write_nipl(path, z, ranges=None) -> Path:
    """Write a 2D intensity array as a DISCUS NIPL file: ``nx ny``, ``xmin xmax ymin ymax``, then
    ny rows of nx values. ``z`` has shape (ny, nx), row j being the j-th line of the plane. The
    ranges are only metadata for DISCUS (``rmc data`` reads but does not use them)."""
    z = np.asarray(z, dtype=float)
    if z.ndim != 2 or min(z.shape) < 2 or not np.isfinite(z).all():
        raise ValueError("z must be a finite 2D array with at least 2 x 2 points")
    ny, nx = z.shape
    r = (0.0, float(nx - 1), 0.0, float(ny - 1)) if ranges is None else tuple(_real("ranges", x) for x in ranges)
    if len(r) != 4:
        raise ValueError("ranges must be (xmin, xmax, ymin, ymax)")
    lines = [f"{nx} {ny}", " ".join(f"{x:.8g}" for x in r)]
    lines += [" ".join(f"{v:.9e}" for v in row) for row in z]
    p = Path(path)
    p.write_text("\n".join(lines) + "\n")
    return p


def read_nipl(path) -> tuple[np.ndarray, tuple[float, float, float, float]]:
    """Read a NIPL file as DISCUS does (leading ``#`` lines skipped, free format): returns
    ``(z, ranges)`` with z of shape (ny, nx)."""
    lines = Path(path).read_text().splitlines()
    k = 0
    while k < len(lines) and lines[k].startswith("#"):
        k += 1
    if k >= len(lines):
        raise ValueError(f"{path}: no NIPL header")
    head = lines[k].split()
    nx, ny = int(head[0]), int(head[1])
    rest = " ".join(lines[k + 1:]).split()
    if len(rest) < 4 + nx * ny:
        raise ValueError(f"{path}: expected 4 + {nx}*{ny} numbers, found {len(rest)}")
    ranges = tuple(float(x) for x in rest[:4])
    z = np.array([float(x) for x in rest[4:4 + nx * ny]]).reshape(ny, nx)
    return z, ranges


# ----------------------------------------------------------------------------------------- data planes
@dataclass(frozen=True)
class RmcPlane:
    """One ``data`` line: a NIPL file with its plane geometry, radiation and weights.

    ``weights``: ``one`` (default), ``sqrt``, ``log``, ``linear``, ``square``, ``inverse``,
    ``inverse_sqrt`` — DISCUS's ``rmc_dowic`` schemes applied to the observed intensity."""
    file: str
    grid: PlaneGrid
    radiation: str = "xray"
    wavelength: str | float = "moa1"
    weights: str = "one"

    def __post_init__(self) -> None:
        _file("file", self.file)
        if not isinstance(self.grid, PlaneGrid):
            raise TypeError("grid must be a PlaneGrid")
        if self.radiation not in RADIATIONS:
            raise ValueError(f"radiation must be one of {RADIATIONS}")
        object.__setattr__(self, "wavelength", _check_wavelength(self.wavelength))
        if self.weights not in WEIGHTS:
            raise ValueError(f"weights must be one of {tuple(WEIGHTS)} (weight files are not offered)")

    def params(self) -> list:
        rad = "xray" if self.radiation == "xray" else "neutron"  # first letter X / N (:1300-1309)
        return [rad, self.wavelength, self.file, WEIGHTS[self.weights], *self.grid.corner_values()]

    def check_file(self, path) -> None:
        """Refuse a file whose size differs from the grid (DISCUS would read it with other axes)."""
        z, _ = read_nipl(path)
        if z.shape != (self.grid.ny, self.grid.nx):
            raise ValueError(f"{path}: {z.shape[1]} x {z.shape[0]} points, the grid says "
                             f"{self.grid.nx} x {self.grid.ny}")


# ----------------------------------------------------------------------------------------- moves
def _prob(p) -> float:
    v = _real("probability", p)
    if v <= 0:
        raise ValueError("probability must be positive")
    return v


@dataclass(frozen=True)
class RmcShift:
    """``set mode, p, shift`` and ``set move, A, sx, sy, sz``: a Gaussian displacement of one atom,
    sigma per axis in **unit-cell fractions** (``gasdev`` added to the fractional position,
    ``rmc_sup_mod.f90:3338-3342``). Every species to be moved needs a sigma (the default is 0.2)."""
    sigma: dict
    probability: float = 1.0

    def __post_init__(self) -> None:
        _prob(self.probability)
        if not self.sigma:
            raise ValueError("RmcShift needs a sigma per species, e.g. {'CU': (0.01, 0.01, 0.01)}")
        clean = {}
        for sp, s in dict(self.sigma).items():
            s = _vec3("sigma", s)
            if any(x < 0 for x in s) or not any(s):
                raise ValueError(f"sigma for {sp!r} must be three non-negative numbers, not all zero")
            clean[_atom("species", sp)] = s
        object.__setattr__(self, "sigma", clean)

    def apply(self, m: Macro) -> None:
        m.command("set", "mode", self.probability, "shift", "all")
        for sp, s in self.sigma.items():
            m.command("set", "move", sp, *s)


@dataclass(frozen=True)
class RmcSwap:
    """``set mode, p, swchem, pairing``: exchange the types of two selected atoms of different type.
    ``pairing``: ``all`` (any two atoms), ``local`` (neighbouring cells, any site), ``slocal``
    (neighbouring cells, same site), ``site`` (same site, any cell) (``rmc_select``, ``:3145-3261``)."""
    probability: float = 1.0
    pairing: str = "all"

    def __post_init__(self) -> None:
        _prob(self.probability)
        if self.pairing not in PAIRINGS:
            raise ValueError(f"pairing must be one of {PAIRINGS}")

    def apply(self, m: Macro) -> None:
        m.command("set", "mode", self.probability, "swchem", self.pairing)


@dataclass(frozen=True)
class RmcSwapDisplacement:
    """``set mode, p, swdisp, pairing``: exchange the displacements (from the site averages) of two
    atoms. DISCUS refuses it (error -12) when no site has a displacement spread."""
    probability: float = 1.0
    pairing: str = "all"

    def __post_init__(self) -> None:
        _prob(self.probability)
        if self.pairing not in PAIRINGS:
            raise ValueError(f"pairing must be one of {PAIRINGS}")

    def apply(self, m: Macro) -> None:
        m.command("set", "mode", self.probability, "swdisp", self.pairing)


_MOVES = (RmcShift, RmcSwap, RmcSwapDisplacement)


@dataclass(frozen=True)
class RmcLots:
    """``set lots, box|ellipsoid, lx, ly, lz, n, yes|no``: compute the intensities from ``count``
    randomly placed lots of lx × ly × lz cells. Runs on the fork (the stock code updates the wrong
    lot and counts lot 1 twice)."""
    size: tuple
    count: int
    shape: str = "box"
    periodic: bool = True

    def __post_init__(self) -> None:
        size = tuple(_int("lot size", n) for n in self.size)
        if len(size) != 3:
            raise ValueError("lot size must be (lx, ly, lz)")
        object.__setattr__(self, "size", size)
        _int("count", self.count)
        if self.shape not in ("box", "ellipsoid"):
            raise ValueError("shape must be 'box' or 'ellipsoid'")
        if not isinstance(self.periodic, bool):
            raise TypeError("periodic must be True or False")

    def params(self) -> list:
        return ["lots", self.shape, *self.size, self.count, "yes" if self.periodic else "no"]


# ----------------------------------------------------------------------------------------- the fit
@dataclass(frozen=True)
class RmcFit(Step):
    """One RMC fit in the ``rmc`` menu, applied to the job's crystal.

    Parameters
    ----------
    planes: the data planes (files must be given in the job's inputs).
    moves: :class:`RmcShift`, :class:`RmcSwap`, :class:`RmcSwapDisplacement`, each type once.
    select: atom names that may be moved or swapped (``"all"`` for every type).
    tries: valid tries (``set cycl``).
    sigma: σ of the acceptance rule (``set sigm``); 0 keeps downhill moves only. Δχ² is in the
        units of the intensities squared, so σ must be chosen relative to the data.
    display: generated moves between log lines (``set disp``).
    fit_scale / scale: fit the scale in closed form, or keep it fixed at ``scale`` (DISCUS's 1/s,
        the number its log prints).
    fit_background / background: the same for a constant background (a fitted background that
        comes out ≤ 0 is replaced by the fixed one).
    min_distances: ``{(A, B): d}`` (Å, ``set mdis``); a move closer than d to another atom is
        refused. DISCUS's default is 0.5 Å for every pair.
    debye_waller: use the atoms' B in the structure factors (``set dbw``, default off). Must match
        how the data were computed (``FourierPlane(debye_waller=...)``).
    symmetry: also fit the symmetry-equivalent planes of the Laue group (``set sym``).
    average_percent: ``set aver`` (subtract the average structure); > 0 runs on the fork.
    lots: :class:`RmcLots` or None.
    q_range: ``(qmin, qmax)`` in 1/d (Å⁻¹, = 2 sinθ/λ) for ``set range`` (DISCUS stores half).
    """
    planes: tuple
    moves: tuple
    select: tuple | str
    tries: int
    sigma: float = 0.01
    display: int = 500
    fit_scale: bool = True
    scale: float = 1.0
    fit_background: bool = False
    background: float = 0.0
    min_distances: dict = field(default_factory=dict)
    debye_waller: bool = False
    symmetry: bool = False
    average_percent: float = 0.0
    lots: RmcLots | None = None
    q_range: tuple | None = None

    uses_random_numbers = True

    def __post_init__(self) -> None:
        planes = tuple(self.planes)
        if not planes or not all(isinstance(p, RmcPlane) for p in planes):
            raise ValueError("planes must be a non-empty tuple of RmcPlane")
        if len({p.file for p in planes}) != len(planes):
            raise ValueError("two planes read the same file")
        object.__setattr__(self, "planes", planes)
        moves = tuple(self.moves)
        if not moves or not all(isinstance(mv, _MOVES) for mv in moves):
            raise ValueError("moves must be a non-empty tuple of RmcShift / RmcSwap / RmcSwapDisplacement")
        if len({type(mv) for mv in moves}) != len(moves):
            raise ValueError("give each move type once")
        object.__setattr__(self, "moves", moves)
        if isinstance(self.select, str):
            if self.select != "all":
                raise ValueError("select must be 'all' or a tuple of atom names")
        else:
            sel = tuple(_atom("select", s) for s in self.select)
            if not sel:
                raise ValueError("select must name at least one atom type")
            object.__setattr__(self, "select", sel)
            if any(isinstance(mv, RmcSwap) for mv in moves) and len(set(sel)) < 2:
                raise ValueError("swaps need at least two selected types (DISCUS error -11)")
            for mv in moves:
                if isinstance(mv, RmcShift):
                    extra = sorted(set(mv.sigma) - set(sel))
                    if extra:
                        raise ValueError(f"RmcShift sigma for types not selected: {extra}")
        _int("tries", self.tries)
        if _real("sigma", self.sigma) < 0:
            raise ValueError("sigma must be >= 0 (0 accepts downhill moves only)")
        _int("display", self.display)
        for k in ("fit_scale", "fit_background", "debye_waller", "symmetry"):
            if not isinstance(getattr(self, k), bool):
                raise TypeError(f"{k} must be True or False")
        if _real("scale", self.scale) <= 0:
            raise ValueError("scale must be positive")
        if _real("background", self.background) < 0:
            raise ValueError("background must be >= 0 (DISCUS replaces a fitted value <= 0 by it)")
        md = {}
        for pair, d in dict(self.min_distances).items():
            a, b = pair
            if _real("min distance", d) < 0:
                raise ValueError("min distances must be >= 0")
            md[(_atom("min distance type", a), _atom("min distance type", b))] = float(d)
        object.__setattr__(self, "min_distances", md)
        if not 0 <= _real("average_percent", self.average_percent) <= 100:
            raise ValueError("average_percent must be in [0, 100]")
        if self.lots is not None and not isinstance(self.lots, RmcLots):
            raise TypeError("lots must be an RmcLots or None")
        if self.q_range is not None:
            lo, hi = (_real("q_range", x) for x in self.q_range)
            if not 0 <= lo < hi:
                raise ValueError("q_range must be (qmin, qmax) with 0 <= qmin < qmax")
            object.__setattr__(self, "q_range", (lo, hi))

    def input_names(self) -> tuple[str, ...]:
        return tuple(p.file for p in self.planes)

    def apply(self, m: Macro) -> None:
        with m.menu("rmc"):
            for p in self.planes:
                m.command("data", *p.params())
            for i in range(1, len(self.planes) + 1):
                m.command("set", "const", i, i)  # rmc_constrain starts at 0 (see module docstring)
            if self.debye_waller:
                m.command("set", "dbw", "on")
            if self.symmetry:
                m.command("set", "sym", "on")
            m.command("dese", "all")
            m.command("sele", *(("all",) if self.select == "all" else self.select))
            for mv in self.moves:
                mv.apply(m)
            for (a, b), d in self.min_distances.items():
                m.command("set", "mdis", a, b, d)
            m.command("set", "scal", "on" if self.fit_scale else "off", self.scale)
            m.command("set", "back", "on" if self.fit_background else "off", self.background)
            if self.average_percent:
                m.command("set", "aver", self.average_percent)
            if self.lots is not None:
                m.command("set", *self.lots.params())
            if self.q_range is not None:
                m.command("set", "range", *self.q_range)
            m.command("set", "sigm", self.sigma)
            m.command("set", "cycl", self.tries)
            m.command("set", "disp", self.display)
            m.command("run")

    def check_inputs(self, paths: Mapping[str, Path | str]) -> None:
        """Check every plane file against its grid."""
        for p in self.planes:
            p.check_file(paths[p.file])


# ----------------------------------------------------------------------------------------- data from a model
@dataclass(frozen=True)
class FourierPlane(Step):
    """Compute the intensity on a :class:`PlaneGrid` with the ``fourier`` menu and write it as a NIPL
    file (``output``: ``value int``, ``form stan``), in the orientation ``rmc data`` reads. Use it
    on a model to make a target plane; give the same radiation, wavelength and Debye–Waller choice
    to the :class:`RmcPlane` and :class:`RmcFit`."""
    grid: PlaneGrid
    output: str
    radiation: str = "xray"
    wavelength: str | float = "moa1"
    debye_waller: bool = False

    def __post_init__(self) -> None:
        if not isinstance(self.grid, PlaneGrid):
            raise TypeError("grid must be a PlaneGrid")
        _file("output", self.output)
        if self.radiation not in RADIATIONS:
            raise ValueError(f"radiation must be one of {RADIATIONS}")
        object.__setattr__(self, "wavelength", _check_wavelength(self.wavelength))
        if not isinstance(self.debye_waller, bool):
            raise TypeError("debye_waller must be True or False")

    def apply(self, m: Macro) -> None:
        g = self.grid
        with m.menu("fourier"):
            m.command(_RADIATION_KEYWORD[self.radiation])
            m.command("wvle", self.wavelength)
            m.command("ll", *g.ll)
            m.command("lr", *g.lr)
            m.command("ul", *g.ul)
            m.command("tl", *g.ll)
            m.command("na", g.nx)
            m.command("no", g.ny)
            m.command("nt", 1)
            m.command("temp", "use" if self.debye_waller else "ignore")
            m.command("disp", "off")
            m.command("lots", "off")
            m.command("set", "aver", 0.0)
            m.command("run")
        with m.menu("output"):
            m.command("value", "int")
            m.command("form", "stan")
            m.command("outf", self.output)
            m.command("run")


# ----------------------------------------------------------------------------------------- the log
@dataclass(frozen=True)
class RmcProgress:
    """One progress block of the log (format 1300/1400, ``rmc_sup_mod.f90:2166-2169``)."""
    generated: int
    tries: int
    accepted_downhill: int | None
    accepted_uphill: int | None
    chi2: float
    planes: tuple  # ((scale, background, chi2), ...) with DISCUS's scale 1/s
    final: bool = False


@dataclass(frozen=True)
class RmcLog:
    progress: tuple
    delta_chi2: tuple | None  # (average, sigma, max) of the uphill Δχ², if any were tried

    @property
    def final(self) -> RmcProgress | None:
        for p in reversed(self.progress):
            if p.final:
                return p
        return None


_INT = r"(\d+|\*+)"
_FLT = r"([-+]?(?:\d+\.?\d*|\.\d+)(?:[eEdD][-+]?\d+)?|NaN|Infinity|\*+)"
_GEN = re.compile(r"Gen:\s*" + _INT + r"\s+try:\s*" + _INT + r"\s+acc:\s*\(good/bad\):\s*" + _INT
                  + r"\s*/\s*" + _INT + r"\s+s2x2:\s*" + _FLT)
_PLANE = re.compile(r"Plane\s*(\d+):\s*scal:\s*" + _FLT + r"\s*/\s*back:\s*" + _FLT + r"\s*/\s*s2x2:\s*" + _FLT)
_DELTA = re.compile(r"Delta Chi\s*:\s*ave:\s*" + _FLT + r"\s+sig:\s*" + _FLT + r"\s+max:\s*" + _FLT)


def _i(s: str) -> int | None:
    return None if s.startswith("*") else int(s)


def _f(s: str) -> float:
    if s.startswith("*"):
        return float("nan")
    return float(s.replace("D", "E").replace("d", "e"))


def parse_rmc_log(text: str) -> RmcLog:
    """Every progress block in a DISCUS log, in order; the last one after
    ``---- Final configuration ----`` is marked ``final``. Counts printed as ``****`` (overflow of
    I6) are None."""
    progress: list[RmcProgress] = []
    cur = None
    planes: list = []
    final_next = False
    delta = None

    def close():
        nonlocal cur, planes
        if cur is not None:
            progress.append(RmcProgress(*cur, planes=tuple(planes), final=cur_final))
        cur, planes = None, []

    cur_final = False
    for line in text.splitlines():
        if "Final configuration" in line:
            final_next = True
            continue
        g = _GEN.search(line)
        if g:
            close()
            cur = (_i(g[1]), _i(g[2]), _i(g[3]), _i(g[4]), _f(g[5]))
            cur_final, final_next = final_next, False
            continue
        p = _PLANE.search(line)
        if p and cur is not None:
            planes.append((_f(p[2]), _f(p[3]), _f(p[4])))
            continue
        d = _DELTA.search(line)
        if d:
            delta = (_f(d[1]), _f(d[2]), _f(d[3]))
    close()
    return RmcLog(tuple(progress), delta)


@dataclass(frozen=True)
class RmcOutcome:
    result: JobResult
    log: RmcLog


def run_rmc(job: StructureJob, runner, inputs: Mapping[str, Path | str] | None = None, **kw) -> RmcOutcome:
    """Run a job containing one :class:`RmcFit` and parse its log. The plane files are checked
    against their grids first. Raises RuntimeError when the log has no final block."""
    fits = [s for s in job.steps if isinstance(s, RmcFit)]
    if len(fits) != 1:
        raise ValueError("run_rmc needs a job with exactly one RmcFit step")
    paths = dict(inputs or {})
    missing = [p.file for p in fits[0].planes if p.file not in paths]
    if missing:
        raise ValueError(f"inputs must contain {missing}")
    fits[0].check_inputs(paths)
    result = job.run(runner, paths, **kw)
    log = parse_rmc_log(result.run.log.read_text(errors="replace") if result.run.log.is_file() else "")
    if log.final is None:
        raise RuntimeError(f"the RMC log has no final configuration; see {result.run.directory}")
    return RmcOutcome(result, log)
