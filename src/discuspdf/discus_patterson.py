"""Patterson and inverse Fourier maps of a reciprocal plane (the DISCUS ``patterson`` and ``inverse`` menus).

Both menus run ``patters.f90``: they read values on a reciprocal plane (a NIPL file written by
:class:`~discuspdf.discus_rmc.FourierPlane`) and sum them over a real-space plane::

    rho(r) = s / V_cell * sum_h A(h) exp(+2 pi i h.r),   s = (atoms per unit cell) / (atoms in the crystal)

with A = I for :class:`PattersonMap` and A = Re F + i Im F for :class:`InverseMap`; r in unit-cell fractions. The
normalisation s is set from the crystal loaded when the menu is entered (patters.f90:70-71), so the job must hold the
crystal the data were computed from (or one with the same counts). The map is written by the ``output`` menu.

Facts (stock 6.23.01, probed 2026-11-05 with ``tools/probes/patterson_probes.py``, mailbox tasks 025-028; checks to
3e-5 of max |rho|):

* **The reader cannot skip the ``#`` header lines** that DISCUS's own ``output`` writes at the top of a NIPL file
  ("Error reading file", patters.f90:1689-1723). :func:`run_map` gives the menus copies without them.
* **The first ordinate point of every abscissa column after the first is wrong**: the phase is not reset when the
  outer loop advances (patters.f90:1872-1885; the point gets the phase of the previous column's end). The steps
  therefore ask DISCUS for one extra ordinate row below the plane, and :func:`read_map` drops it.
* The sign of the exponent is +; the Patterson and the inverse map follow the formula above. ``sign`` is overwritten
  by +1 in every run (patters.f90:1840) and is not offered; ``sign -1`` is a syntax error anyway.
* Not offered: ``diff-four`` (with NIPL input it writes an all-zero map: ``set_patt_value`` has no branch for it,
  patters.f90:1905-1992); SHELXL/HKLF4 input (symmetry expansion, sharpening and Wilson scaling are untested); the
  ``set`` options (accumulation, exclusion, sharpening, origin), which a single job uses at their start values.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Mapping

import numpy as np

from .discus_ops import Step
from .discus_rmc import PlaneGrid, _plane_grid, read_nipl
from .macro import Macro
from ._validate import integer as _int, plain_file_name as _file, vec3 as _vec3

__all__ = ["InverseMap", "PattersonMap", "RealPlane", "read_map", "run_map", "strip_nipl_header"]


@dataclass(frozen=True)
class RealPlane:
    """A real-space plane in unit-cell fractions: point (a, o) = ll + a (lr - ll)/(na - 1) + o (ul - ll)/(no - 1)."""
    ll: tuple
    lr: tuple
    ul: tuple
    na: int
    no: int

    def __post_init__(self) -> None:
        for k in ("ll", "lr", "ul"):
            object.__setattr__(self, k, _vec3(k, getattr(self, k)))
        for k in ("na", "no"):
            if _int(k, getattr(self, k)) < 2:
                raise ValueError(f"{k} must be at least 2")
        if np.allclose(np.subtract(self.lr, self.ll), 0) or np.allclose(np.subtract(self.ul, self.ll), 0):
            raise ValueError("lr and ul must differ from ll")

    def points(self) -> np.ndarray:
        """(no, na, 3) array of the plane's points (rows = ordinate), as :func:`read_map` returns the map."""
        ll, lr, ul = (np.asarray(v, float) for v in (self.ll, self.lr, self.ul))
        a = np.arange(self.na)[None, :, None]
        o = np.arange(self.no)[:, None, None]
        return ll + a * (lr - ll) / (self.na - 1) + o * (ul - ll) / (self.no - 1)

    def extended(self) -> "RealPlane":
        """The plane with one more ordinate row below ll (the row DISCUS gets wrong)."""
        d = (np.asarray(self.ul, float) - np.asarray(self.ll, float)) / (self.no - 1)
        clean = lambda v: tuple(float(x) + 0.0 for x in np.round(v, 12))  # noqa: E731
        return RealPlane(clean(np.asarray(self.ll) - d), clean(np.asarray(self.lr) - d), tuple(self.ul),
                         self.na, self.no + 1)


def _menu_lines(m: Macro, menu: str, files, grid: PlaneGrid, plane: RealPlane) -> None:
    ext = plane.extended()
    with m.menu(menu):
        m.command("reset")
        m.command("format", "stan")
        for ab, (name, kind) in zip("ab", files):
            m.command("file", ab, name)
            m.command("type", ab, kind)
        m.command("ll", *grid.ll)
        m.command("lr", *grid.lr)
        m.command("ul", *grid.ul)
        m.command("na", grid.nx)
        m.command("no", grid.ny)
        m.command("rholl", *ext.ll)
        m.command("rholr", *ext.lr)
        m.command("rhoul", *ext.ul)
        m.command("rhona", ext.na)
        m.command("rhono", ext.no)
        m.command("run")


def _output(m: Macro, value: str, name: str) -> None:
    with m.menu("output"):
        m.command("value", value)
        m.command("form", "stan")
        m.command("outf", name)
        m.command("run")


@dataclass(frozen=True)
class PattersonMap(Step):
    """The Patterson map of an intensity plane: ``intensity`` is a NIPL file on ``grid`` (for example from
    ``FourierPlane(grid, file)``), ``plane`` the real-space plane, ``output`` the NIPL map file DISCUS writes (read it
    with :func:`read_map`, which drops the extra row). Run the job with :func:`run_map`."""
    intensity: str
    grid: PlaneGrid
    plane: RealPlane
    output: str

    def __post_init__(self) -> None:
        _file("intensity", self.intensity)
        _file("output", self.output)
        object.__setattr__(self, "grid", _plane_grid(self.grid))
        if not isinstance(self.plane, RealPlane):
            raise TypeError("plane must be a RealPlane")
        if self.output == self.intensity:
            raise ValueError("output must differ from the input")

    def input_names(self) -> tuple[str, ...]:
        return (self.intensity,)

    def output_names(self) -> tuple[str, ...]:
        return (self.output,)

    def apply(self, m: Macro) -> None:
        _menu_lines(m, "patterson", [(self.intensity, "intensity")], self.grid, self.plane)
        _output(m, "real", self.output)


@dataclass(frozen=True)
class InverseMap(Step):
    """The inverse Fourier transform of a complex amplitude plane: ``real`` and ``imag`` are NIPL files of Re F and
    Im F on ``grid`` (``FourierPlane(grid, file, value="real")`` / ``"imag"``). ``output`` gets the real part of the
    map, ``imag_output`` (optional) the imaginary part."""
    real: str
    imag: str
    grid: PlaneGrid
    plane: RealPlane
    output: str
    imag_output: str | None = None

    def __post_init__(self) -> None:
        for k in ("real", "imag", "output"):
            _file(k, getattr(self, k))
        if self.imag_output is not None:
            _file("imag_output", self.imag_output)
        object.__setattr__(self, "grid", _plane_grid(self.grid))
        if not isinstance(self.plane, RealPlane):
            raise TypeError("plane must be a RealPlane")
        names = [self.real, self.imag, self.output] + ([self.imag_output] if self.imag_output else [])
        if len(set(names)) != len(names):
            raise ValueError("the input and output files must all differ")

    def input_names(self) -> tuple[str, ...]:
        return (self.real, self.imag)

    def output_names(self) -> tuple[str, ...]:
        return (self.output,) + ((self.imag_output,) if self.imag_output else ())

    def apply(self, m: Macro) -> None:
        _menu_lines(m, "inverse", [(self.real, "real"), (self.imag, "imag")], self.grid, self.plane)
        _output(m, "real", self.output)
        if self.imag_output:
            _output(m, "imag", self.imag_output)


def strip_nipl_header(src, dst) -> Path:
    """Copy a NIPL file without its leading ``#`` lines (the patterson/inverse reader cannot skip them)."""
    lines = Path(src).read_text(encoding="utf-8", errors="replace").splitlines()
    dst = Path(dst)
    dst.write_text("\n".join(l for l in lines if not l.lstrip().startswith("#")) + "\n", encoding="ascii")
    return dst


def read_map(path, plane: RealPlane) -> np.ndarray:
    """The map of a :class:`PattersonMap` / :class:`InverseMap` output on ``plane``: shape (no, na), rows =
    ordinate, as :meth:`RealPlane.points`; the extra row DISCUS computes wrongly is dropped."""
    z, _ = read_nipl(path)
    if z.shape != (plane.no + 1, plane.na):
        raise ValueError(f"{path}: map of shape {z.shape}, expected {(plane.no + 1, plane.na)} (plane plus one row)")
    return z[1:]


def run_map(job, runner, inputs: Mapping[str, Path | str], work: Path | str):
    """Run a job containing :class:`PattersonMap` / :class:`InverseMap` steps. Their input files (looked up in
    ``inputs``) are copied to ``work`` without the ``#`` header lines first; other inputs are passed as they are.
    Returns the :class:`~discuspdf.discus_ops.JobResult`."""
    work = Path(work)
    work.mkdir(parents=True, exist_ok=True)
    paths = dict(inputs)
    for step in job.steps:
        if isinstance(step, (PattersonMap, InverseMap)):
            for name in step.input_names():
                if name not in paths:
                    raise ValueError(f"inputs must contain {name!r}")
                if Path(paths[name]).resolve() != (work / name).resolve():
                    paths[name] = strip_nipl_header(paths[name], work / name)
    return job.run(runner, paths)
