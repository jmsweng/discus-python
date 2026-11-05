"""Experimental PDF from a measured powder pattern with DISCUS's ``exp2pdf`` menu (stock DISCUS 6.23.01).

``exp2pdf`` subtracts an optional scaled background, divides by the mean form factor squared <f>^2,
multiplies by Q, subtracts a fitted polynomial, and sine-transforms F(Q) to G(r). **The scale of G(r) is
empirical** (a data-dependent factor, discus_exp2pdf_run.f90:826-827): use the *shape*; a fit against a model
needs a free scale factor. Facts below are from the source (``docs/DISCUS_PDFMENU_EXP2PDF_FACTS.md`` section
6) and from three probe rounds on the real engine (``tools/probes/exp2pdf_probes.py``, 2026-11-05):

* ``run`` finishes only after a crystal is read, a ``powder`` calculation is run on it in the same session, and
  ``composition comp:current`` is used (round 2). ``composition <string>`` builds an internal crystal that fails
  ("Unknown element"), ``comp:current`` without the powder run aborts DISCUS, and the other failures show up as
  "Fractions add up to 0.000". The wrapper therefore writes a small composition crystal (one atom per species,
  occupancy = its fraction; round 3: occupancies give the same G(r) as whole atoms), runs a short powder
  calculation on it and then ``exp2pdf``. A powder run on one crystal followed by ``read`` of another aborts.
* ``limits inst:`` alone fails (the Fourier default [0, 1e9] is still checked, XS:445-647); the wrapper writes
  ``limits`` only when a limit is asked for, and then both ranges in full, inside the data range.
* DISCUS moves the Fourier Qmax to a nearby zero of F(Q) (16.5 asked -> 16.473 used, round 3); the limits used
  are in :attr:`ExpPdfResult.header` (``lim_four_app``).
* The default G(r) file is HDF5; the wrapper always asks for the text format (``format:stan``).
* G(r) = (2/pi) sum dq F(Q) sin(Qr) over the written F(Q) with a linear ramp from Q = 0 (round 3: 2.8e-4
  of max |G|); no Lorch window.
* The electron energy and the wavelength do not change G(r) when the data are on a Q scale (round 3:
  identical files); the data must be on a Q scale, so 2-theta data are converted here with the wavelength.
* Background with ``scale:`` and a third (sigma) column: as documented (round 3: 1e-11 and 0 differences).
* Not offered: ``qscale`` (ignores ``dcrystal``), ``faver`` (the file branch never records the file),
  KUPLOT data sets.
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Mapping

import numpy as np

from ._validate import positive as _pos, real as _real
from .atom_names import check_discus_name, discus_atom_name

__all__ = ["ExpPdf", "ExpPdfResult", "read_pattern", "read_exp2pdf_file", "run_exp2pdf", "two_theta_to_q"]

RADIATIONS = ("xray", "neutron", "electron")
DATA_NAME, BACK_NAME, COMP_NAME = "e2p_data.xy", "e2p_back.xy", "e2p_comp.stru"
GR_NAME, FQ_NAME, SQ_NAME, IQ_NAME = "e2p_gr.dat", "e2p_fq.dat", "e2p_sq.dat", "e2p_iq.dat"
# The powder calculation that sets DISCUS up for exp2pdf (round 2/3 set-up; its own settings are overwritten
# by exp2pdf, XR:395-427).
_POWDER = ("powder", "reset", "xray", "set wvle,0.5", "set calc,debye", "set qmin,0.5", "set qmax,17.9",
           "set dq,0.01", "set stepr,0.001", "set profile,off", "set lpcor,none", "set temp,ignore",
           "set disp,off", "set period,off", "run", "exit")


def two_theta_to_q(two_theta_deg, wavelength: float) -> np.ndarray:
    """Q = 4 pi sin(theta) / lambda (A^-1) for 2-theta in degrees and lambda in A."""
    lam = _pos("wavelength", wavelength)
    tth = np.asarray(two_theta_deg, dtype=float)
    if np.any((tth <= 0) | (tth >= 180)):
        raise ValueError("2-theta must be between 0 and 180 degrees")
    return 4 * math.pi * np.sin(np.radians(tth) / 2) / lam


def read_pattern(path) -> np.ndarray:
    """Columns (x, y[, dy]) of a powder pattern file: whitespace or comma separated; lines that do not
    start with two numbers (headers, ``#`` comments, .chi preambles) are skipped. Returns an (n, 2|3) array
    sorted by x, with repeated x refused."""
    rows = []
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        parts = line.replace(",", " ").split()
        try:
            vals = [float(p) for p in parts[:3]]
        except ValueError:
            continue
        if len(vals) >= 2 and all(math.isfinite(v) for v in vals):
            rows.append(vals)
    if len(rows) < 10:
        raise ValueError(f"no usable (x, y) columns in {path}")
    width = min(len(r) for r in rows)
    a = np.array([r[:width] for r in rows])
    a = a[np.argsort(a[:, 0], kind="stable")]
    if np.any(np.diff(a[:, 0]) <= 0):
        raise ValueError(f"repeated x values in {path}")
    return a


def _composition(comp) -> tuple[tuple[str, float], ...]:
    if isinstance(comp, str):
        items = []
        for tok in comp.replace(",", " ").split():
            m = re.fullmatch(r"([A-Z][a-z]?)(\d*\.?\d*)", tok)
            if not m:
                raise ValueError(f"composition token {tok!r}: write e.g. 'Cu3 Au1' or a dict {{'Cu': 3, 'Au': 1}}")
            items.append((m.group(1), float(m.group(2)) if m.group(2) else 1.0))
    elif isinstance(comp, Mapping):
        items = list(comp.items())
    else:
        raise TypeError("composition must be a string such as 'Fe1 O0.978' or a mapping {'Fe': 1, 'O': 0.978}")
    out, seen = [], set()
    for name, amount in items:
        name = check_discus_name("composition atom", str(name))
        if discus_atom_name(name) is None:
            raise ValueError(f"composition atom {name!r} is not an element DISCUS knows")
        if name.upper() in seen:
            raise ValueError(f"composition lists {name!r} twice")
        seen.add(name.upper())
        out.append((name, _pos(f"amount of {name}", amount)))
    if not out:
        raise ValueError("composition is empty")
    if len(out) > 99:
        raise ValueError("at most 99 species")
    total = sum(a for _, a in out)
    return tuple((n, a / total) for n, a in out)


def _comp_stru(fractions) -> str:
    """A P1 crystal (cell 20 A) with one atom per species, occupancy = its fraction (round 3 layout)."""
    lines = ["title discuspdf exp2pdf composition crystal", "spcgr P1",
             "cell    20.000000,  20.000000,  20.000000, 90.000000, 90.000000, 90.000000",
             f"ncell        1,        1,        1,       {len(fractions)}", "format numbers,XYZBPMMOS",
             "atoms      x,              y,              z,             Biso,    Property,  MoleNo,  MoleAt,   Occ,"
             "     St,  Sh,  Sk,  Sl"]
    for k, (name, occ) in enumerate(fractions):
        x = 0.2 * (k % 5)
        y = 0.2 * (k // 5 % 5)
        z = 0.2 * (k // 25)
        lines.append(f"{name:<10} {x:.6f},       {y:.6f},       {z:.6f},      0.000000,       1,       0,"
                     f"       0,   {occ:.6f}, _,   0,   0,   0")
    return "\n".join(lines) + "\n"


@dataclass(frozen=True)
class ExpPdf:
    """One ``exp2pdf`` reduction of a measured powder pattern to G(r) (shape only: the scale is empirical).

    ``data`` and ``background`` are file paths (x, y[, sigma]); ``axis`` says whether x is Q (A^-1) or
    2-theta (degrees; then ``wavelength`` in A is required and both files are converted to Q). The background
    must cover the data's Q range. ``composition`` is a string ``"Cu3 Au1"`` / ``"Fe1 O0.978"`` or a mapping.
    ``qmin``/``qmax`` are the instrument limits and ``fourier_qmin``/``fourier_qmax`` the transform limits
    (A^-1; any one given makes the wrapper write all four, the missing ones from the data range or the
    instrument limits). ``poly_order`` is the polynomial subtracted from F(Q) (DISCUS default 7).
    ``rmin``, ``rmax``, ``rstep`` set the G(r) grid. ``outputs`` may add ``"fq"``, ``"sq"``, ``"iq"``."""

    data: Path | str
    composition: str | Mapping[str, float]
    radiation: str = "xray"
    axis: str = "q"
    wavelength: float | None = None
    energy: float | None = None
    background: Path | str | None = None
    background_scale: float = 1.0
    qmin: float | None = None
    qmax: float | None = None
    fourier_qmin: float | None = None
    fourier_qmax: float | None = None
    poly_order: int = 7
    rmin: float = 0.01
    rmax: float = 30.0
    rstep: float = 0.01
    outputs: tuple[str, ...] = ("fq",)
    fractions: tuple = field(init=False, repr=False)

    def __post_init__(self):
        if self.radiation not in RADIATIONS:
            raise ValueError(f"radiation must be one of {RADIATIONS}")
        if self.axis not in ("q", "tth"):
            raise ValueError("axis must be 'q' or 'tth' (2-theta in degrees)")
        if self.wavelength is not None and self.energy is not None:
            raise ValueError("give the wavelength or the energy, not both")
        if self.wavelength is not None:
            _pos("wavelength", self.wavelength)
        if self.energy is not None:
            _pos("energy", self.energy)
        if self.axis == "tth" and self.wavelength is None:
            raise ValueError("2-theta data need the wavelength (A)")
        _pos("background_scale", self.background_scale)
        if not isinstance(self.poly_order, int) or isinstance(self.poly_order, bool) or not 1 <= self.poly_order <= 20:
            raise ValueError("poly_order must be an integer 1..20")
        for name in ("qmin", "qmax", "fourier_qmin", "fourier_qmax"):
            if getattr(self, name) is not None:
                _real(name, getattr(self, name))
        _pos("rmin", self.rmin), _pos("rstep", self.rstep)
        if not self.rmax > self.rmin:
            raise ValueError("rmax must be larger than rmin")
        if int((self.rmax - self.rmin) / self.rstep) + 1 > 200001:
            raise ValueError("more than 200001 G(r) points")
        outs = tuple(self.outputs)
        if not set(outs) <= {"fq", "sq", "iq"}:
            raise ValueError("outputs may contain 'fq', 'sq', 'iq'")
        object.__setattr__(self, "outputs", outs)
        object.__setattr__(self, "fractions", _composition(self.composition))

    # ------------------------------------------------------------------ inputs
    def _q_columns(self, path) -> np.ndarray:
        a = read_pattern(path)
        if self.axis == "tth":
            a[:, 0] = two_theta_to_q(a[:, 0], self.wavelength)
        return a

    def limits(self, q_data) -> tuple[tuple[float, float], tuple[float, float]] | None:
        """(inst, fourier) ranges written to ``limits``, or None when none was asked for."""
        asked = (self.qmin, self.qmax, self.fourier_qmin, self.fourier_qmax)
        if all(v is None for v in asked):
            return None
        lo, hi = float(q_data[0]), float(q_data[-1])
        iq = (lo if self.qmin is None else float(self.qmin), hi if self.qmax is None else float(self.qmax))
        fq = (iq[0] if self.fourier_qmin is None else float(self.fourier_qmin),
              iq[1] if self.fourier_qmax is None else float(self.fourier_qmax))
        # DISCUS's own checks (XS:576-647), here with a readable message
        if not lo <= iq[0] < iq[1] <= hi:
            raise ValueError(f"instrument limits {iq} must lie inside the data range {lo:.4f} .. {hi:.4f}")
        if not iq[0] <= fq[0] < fq[1] <= iq[1]:
            raise ValueError(f"Fourier limits {fq} must lie inside the instrument limits {iq}")
        return iq, fq

    def macro_lines(self, q_data, with_background: bool) -> list[str]:
        rad = f"radiation {self.radiation}"
        if self.energy is not None:
            rad += f", energy:{self.energy:g}"
        elif self.wavelength is not None:
            rad += f", lambda:{self.wavelength:g}"
        else:
            rad += ", lambda:0.1"  # DISCUS's default; irrelevant on a Q scale (round 3)
        lines = ["discus", "read", f"stru {COMP_NAME}", *_POWDER, "exp2pdf", "reset", f"data xy, {DATA_NAME}"]
        if with_background:
            lines.append(f"back xy, {BACK_NAME}, scale:{self.background_scale:.10g}")
        lines += ["composition comp:current", rad]
        lim = self.limits(q_data)
        if lim is not None:
            (a, b), (c, d) = lim
            up, down = (lambda v: math.ceil(v * 1e6 - 1e-6) / 1e6), (lambda v: math.floor(v * 1e6 + 1e-6) / 1e6)
            a, c = up(a), up(c)            # rounded inwards, so the written values stay inside DISCUS's checks
            b, d = down(b), down(d)
            lines.append(f"limits inst:[{a:.6f},{b:.6f}], fourier:[{c:.6f},{d:.6f}]")
        if self.poly_order != 7:
            lines.append(f"poly order:{self.poly_order}")
        out = (f"outfile gr:{GR_NAME}, rmin:{self.rmin:g}, rmax:{self.rmax:g}, rstep:{self.rstep:g}, format:stan")
        for key, name in (("fq", FQ_NAME), ("sq", SQ_NAME), ("iq", IQ_NAME)):
            if key in self.outputs:
                out += f", {key}:{name}"
        lines += [out, "run mode:silent", "exit", "exit"]
        return lines

    def output_names(self) -> list[str]:
        names = {"fq": FQ_NAME, "sq": SQ_NAME, "iq": IQ_NAME}
        return [GR_NAME] + [names[k] for k in self.outputs]


@dataclass(frozen=True)
class ExpPdfResult:
    r: np.ndarray
    g: np.ndarray
    header: dict
    fq: tuple[np.ndarray, np.ndarray] | None
    sq: tuple[np.ndarray, np.ndarray] | None
    iq: tuple[np.ndarray, np.ndarray] | None
    run: object

    @property
    def fourier_limits(self) -> tuple[float, float] | None:
        """The Q range DISCUS transformed (``lim_four_app`` of the header)."""
        v = self.header.get("lim_four_app")
        return (float(v[0]), float(v[1])) if v else None


def read_exp2pdf_file(path) -> tuple[np.ndarray, np.ndarray, dict]:
    """(x, y, header) of an exp2pdf text file (``format:stan``); header values split at commas."""
    header, rows = {}, []
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        s = line.strip()
        if s.startswith("#"):
            m = re.match(r"#\s*\w+\s*:\s*(.+?)\s*;\w+::\s*(.*)$", s)
            if m:
                header[m.group(1).strip()] = [v.strip() for v in m.group(2).split(",")]
            continue
        parts = s.replace("D+", "E+").replace("D-", "E-").split()
        if len(parts) >= 2:
            rows.append((float(parts[0]), float(parts[1])))
    if not rows:
        raise ValueError(f"no data in {path}")
    a = np.array(rows)
    return a[:, 0], a[:, 1], header


def run_exp2pdf(config: ExpPdf, runner, work_inputs: Path | str | None = None) -> ExpPdfResult:
    """Convert the files, run ``exp2pdf`` with ``runner`` (a :class:`~discuspdf.suite_runner.SuiteRunner`) and
    read the result. The converted input files are written to ``work_inputs`` (default: a folder
    ``exp2pdf_inputs`` next to the runner's run directories) and staged with their hashes."""
    data = config._q_columns(config.data)
    back = config._q_columns(config.background) if config.background is not None else None
    if back is not None and not (back[0, 0] <= data[0, 0] and back[-1, 0] >= data[-1, 0]):
        raise ValueError(f"the background ({back[0, 0]:.4f} .. {back[-1, 0]:.4f} A^-1) must cover the data's Q range "
                         f"({data[0, 0]:.4f} .. {data[-1, 0]:.4f})")
    folder = Path(work_inputs) if work_inputs is not None else Path(runner.runs_directory) / "exp2pdf_inputs"
    folder.mkdir(parents=True, exist_ok=True)
    inputs = {DATA_NAME: folder / DATA_NAME, COMP_NAME: folder / COMP_NAME}
    np.savetxt(inputs[DATA_NAME], data, fmt="%.10g")
    if back is not None:
        inputs[BACK_NAME] = folder / BACK_NAME
        np.savetxt(inputs[BACK_NAME], back, fmt="%.10g")
    inputs[COMP_NAME].write_text(_comp_stru(config.fractions), encoding="ascii")
    text = "\n".join(config.macro_lines(data[:, 0], back is not None)) + "\n"
    run = runner.run(text, inputs, config.output_names())
    r, g, header = read_exp2pdf_file(run.output(GR_NAME))
    extra = {}
    for key, name in (("fq", FQ_NAME), ("sq", SQ_NAME), ("iq", IQ_NAME)):
        extra[key] = read_exp2pdf_file(run.output(name))[:2] if key in config.outputs else None
    return ExpPdfResult(r, g, header, extra["fq"], extra["sq"], extra["iq"], run)
