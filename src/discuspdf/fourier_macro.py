"""Python control of DISCUS 6.23.01 single-crystal Fourier runs: Bragg subtraction and LOTS (Package B, 8–9).

:class:`FourierConfig` writes the ``fourier`` and ``output`` blocks of a DISCUS macro using only
commands whose behaviour was observed on DISCUS_SUITE 6.23.01
(``benchmarks/discus_6_23_01/b_probes/RESULT.md``, Gates 2–3):

* ``set aver,<p>`` subtracts the **average amplitude**: intensity = |F − ⟨F⟩|²; ``value real/imag``
  are F − ⟨F⟩. For ``0 < p < 100`` the unit cells used for ⟨F⟩ are drawn at random, so an explicit
  8-value seed is required (a 1-value ``seed`` was not reproducible).
* ``lots {box|ell},lx,ly,lz,{n|all},{y|n}``: per-lot intensities are summed and scaled to the whole
  crystal; random lot origins need an 8-value seed; only ``value int`` may be written.
* ``set technique:nufft`` is **not offered**: in 6.23.01 it wrote all-zero files with exit code 0
  (FINUFFT ``setpts error``). Only ``turbo`` is generated.
* ``set symm:{friedel|ignore|apply}``; ``friedel`` is the engine default.
* Radiation ``xray``, ``neutron`` or ``electron`` (macro keywords ``xray`` / ``neut`` / ``electron``). The
  radiation line is written **before** ``wvle`` / ``energy``: ``neut`` and ``electron`` clear a previously
  set wavelength symbol.
* ``wvle`` takes a DISCUS tube symbol (``moa1``, ``cua1``, …) or a number in Å; ``energy`` an energy (keV
  for X-rays and electrons, meV for neutrons). Without anomalous dispersion (``disp off``, the default)
  DISCUS 6.23.01 uses the wavelength only for metadata: f₀(sinθ/λ) and neutron b do not depend on it.
  With ``dispersion=True`` (``disp anom``) the internal f′, f″ exist only for the Kα1 symbols.

The macro is plain text you can read before running it (``render``). Nothing here runs DISCUS;
:func:`prepare_fourier_run` stages the macro for :class:`~discuspdf.diffuse3d.Diffuse3DRunner`.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from itertools import permutations
from pathlib import Path
from typing import Any
import json
import re
import tempfile

import numpy as np

from .reciprocal_grid import ReciprocalGrid

SCHEMA_NAME = "discuspdf.FourierConfig"
SCHEMA_VERSION = 2          # v2 added 'radiation' and numeric wavelengths; v1 descriptions still read
OUTPUT_VALUES = ("int", "real", "imag", "3DPDF")
SYMMETRY = ("friedel", "ignore", "apply")
RADIATIONS = ("xray", "neutron", "electron")
_RADIATION_KEYWORD = {"xray": "xray", "neutron": "neut", "electron": "electron"}   # the keywords run in the probes
# Internal f', f'' exist only for these K-alpha1 symbols (lib_f90/element_data_mod.f90 get_scat_ano, :3309-3355);
# any other wavelength (number, energy, a2 or b line) gives f' = f'' = 0 without a message.
ANOMALOUS_SYMBOLS = ("tia1", "cra1", "fea1", "coa1", "cua1", "moa1", "aga1", "taa1", "wa1", "aua1")
ENERGY_UNITS = {"xray": "keV", "electron": "keV", "neutron": "meV"}  # get_wave, element_data_mod.f90:3027-3036
# 'value 3DPDF, qmax:/dsmax:/hklmax:' is not offered: output_limit (graphic.f90:2861-2889) tests |h| inside the
# loop over the components of h and zeroes dsi(i,j,k) with the component index k, i.e. only layers 1-3.
# DISCUS 6.23.01 wavelength symbols (lib_f90/element_data_mod.f90, per_symwl): K-alpha1, K-alpha2 and
# K-beta of the usual anodes. 'moa1' is the default and the one run in all probes.
_ANODES = ("ti", "cr", "fe", "co", "cu", "mo", "ag", "ta", "w", "au")
WAVELENGTHS = tuple(f"{a}{line}" for a in _ANODES for line in ("a1", "a2", "b"))
# per_wavel of the same table (Angstrom), for every symbol in WAVELENGTHS
WAVELENGTH_ANGSTROM = {
    "tia1": 2.74841, "tia2": 2.75207, "tib": 2.51381,
    "cra1": 2.28962, "cra2": 2.29351, "crb": 2.0848,
    "fea1": 1.93597, "fea2": 1.93991, "feb": 1.75653,
    "coa1": 1.78892, "coa2": 1.79278, "cob": 1.62075,
    "cua1": 1.540592, "cua2": 1.544343, "cub": 1.39217,
    "moa1": 0.70926, "moa2": 0.713543, "mob": 0.632253,
    "aga1": 0.559363, "aga2": 0.563775, "agb": 0.49701,
    "taa1": 0.215484, "taa2": 0.22029, "tab": 0.190076,
    "wa1": 0.208992, "wa2": 0.213813, "wb": 0.184363,
    "aua1": 0.180185, "aua2": 0.185064, "aub": 0.158971,
}
_NAME_RE = re.compile(r"^[A-Za-z0-9_.\-]+$")
_LABELS = ("h", "k", "l")


def _fmt(v) -> str:
    return ",".join(f"{float(x):.6g}" for x in v)


def _check_name(kind: str, name: str, suffixes: tuple[str, ...]) -> str:
    if not isinstance(name, str) or not _NAME_RE.match(name) or name in (".", ".."):
        raise ValueError(f"{kind} must be a plain file name (letters, digits, '_', '-', '.'), got {name!r}")
    if not name.lower().endswith(suffixes):
        raise ValueError(f"{kind} {name!r} must end with one of {suffixes}")
    return name


@dataclass(frozen=True)
class Seed:
    """Eight-value DISCUS seed (``seed v1,…,v8``). A single value was not reproducible in 6.23.01."""
    values: tuple[int, int, int, int, int, int, int, int]

    def __post_init__(self) -> None:
        vals = tuple(self.values)
        if len(vals) != 8:
            raise ValueError("Seed needs exactly 8 values; DISCUS 6.23.01 fills missing words from the clock")
        for v in vals:
            if isinstance(v, (bool, np.bool_)) or not isinstance(v, (int, np.integer)) or not 1 <= v <= 99999:
                raise ValueError(f"seed values must be integers in 1..99999 (exact in DISCUS's real parser), got {v!r}")
        object.__setattr__(self, "values", tuple(int(v) for v in vals))

    def command(self) -> str:
        return "seed " + ",".join(str(v) for v in self.values)


@dataclass(frozen=True)
class Lots:
    """LOTS sampling: ``shape`` 'box' or 'ellipsoid', ``size`` in unit cells, ``count`` an integer or
    ``'all'`` (one lot per origin cell), ``periodic`` boundary conditions."""
    shape: str
    size: tuple[int, int, int]
    count: int | str
    periodic: bool

    def __post_init__(self) -> None:
        if self.shape not in ("box", "ellipsoid"):
            raise ValueError("lots shape must be 'box' or 'ellipsoid'")
        size = tuple(self.size)
        if len(size) != 3 or any(isinstance(n, bool) or not isinstance(n, (int, np.integer)) or n < 1 for n in size):
            raise ValueError("lots size must be three integers >= 1 (unit cells)")
        object.__setattr__(self, "size", tuple(int(n) for n in size))
        if self.count != "all" and (isinstance(self.count, bool) or not isinstance(self.count, (int, np.integer))
                                    or self.count < 1):
            raise ValueError("lots count must be a positive integer or 'all'")
        if not isinstance(self.periodic, bool):
            raise ValueError("lots periodic must be True or False")

    @property
    def random(self) -> bool:
        return self.count != "all"

    def command(self) -> str:
        shape = "box" if self.shape == "box" else "ell"
        count = "all" if self.count == "all" else str(int(self.count))
        return f"lots {shape},{self.size[0]},{self.size[1]},{self.size[2]},{count},{'y' if self.periodic else 'n'}"

    def n_lots(self, crystal_cells) -> int:
        """Number of lots DISCUS uses: ``count``, or for 'all' the number of possible origins
        (all cells if periodic, else (N − l + 1) per axis; both observed in 6.23.01)."""
        if self.random:
            return int(self.count)
        cells = np.asarray(crystal_cells, int)
        per = cells if self.periodic else cells - np.asarray(self.size) + 1
        return int(np.prod(per))


@dataclass(frozen=True)
class StructureSource:
    """How the macro loads the crystal: ``read stru <file>`` or ``read cell <file>,n1,n2,n3``.

    ``crystal_cells`` is the crystal size in unit cells; it is required for 'stru' if LOTS are used
    (DISCUS refuses lots larger than the crystal). For a ``.stru`` written by DISCUS it is the
    ``ncell`` line."""
    kind: str
    file: str
    crystal_cells: tuple[int, int, int] | None = None

    def __post_init__(self) -> None:
        if self.kind not in ("stru", "cell"):
            raise ValueError("kind must be 'stru' or 'cell'")
        _check_name("structure file", self.file, (".stru",) if self.kind == "stru" else (".cll", ".cell"))
        if self.kind == "cell" and self.crystal_cells is None:
            raise ValueError("read cell needs crystal_cells (n1, n2, n3)")
        if self.crystal_cells is not None:
            c = tuple(self.crystal_cells)
            if len(c) != 3 or any(isinstance(n, bool) or not isinstance(n, (int, np.integer)) or n < 1 for n in c):
                raise ValueError("crystal_cells must be three integers >= 1")
            object.__setattr__(self, "crystal_cells", tuple(int(n) for n in c))

    @classmethod
    def stru(cls, file: str, crystal_cells=None) -> "StructureSource":
        return cls("stru", file, crystal_cells)

    @classmethod
    def cell(cls, file: str, crystal_cells) -> "StructureSource":
        return cls("cell", file, crystal_cells)

    def lines(self) -> list[str]:
        if self.kind == "stru":
            arg = f"stru {self.file}"
        else:
            n = self.crystal_cells
            arg = f"cell {self.file},{n[0]},{n[1]},{n[2]}"
        return ["read", f"  {arg}", "eval n[1]", "eval n[2]", "eval n[3]"]


def _check_wavelength(w) -> str | float:
    """A known DISCUS symbol (stored lower case) or a finite positive number of Å."""
    if isinstance(w, str):
        sym = w.strip().lower()
        if sym not in WAVELENGTHS:
            raise ValueError(f"wavelength symbol {w!r} is not a DISCUS 6.23.01 symbol {WAVELENGTHS}; "
                             "or give a number in Å")
        return sym
    if isinstance(w, (bool, np.bool_)) or not isinstance(w, (int, float, np.integer, np.floating)):
        raise ValueError(f"wavelength must be a DISCUS symbol or a number in Å, got {w!r}")
    w = float(w)
    if not np.isfinite(w) or w <= 0:
        raise ValueError(f"wavelength must be finite and > 0 Å, got {w!r}")
    return w


# DISCUS atom names: up to 4 characters (longer names break the .stru reader), compared upper case.
_ATOM_NAME_RE = re.compile(r"^[A-Z][A-Z0-9+\-.]{0,3}$")


def normalize_scattering_lengths(value) -> tuple[tuple[str, float], ...]:
    """Own coherent neutron scattering lengths (fm) per DISCUS atom name, as sorted ``(NAME, b)`` pairs.

    Accepts a mapping ``{"NI62": 10.3}`` or pairs. Names are upper-cased (DISCUS ignores case) and must
    be DISCUS atom names (1–4 characters, starting with a letter); names equal up to case are refused.
    """
    if value is None:
        return ()
    items = list(value.items()) if hasattr(value, "items") else [tuple(p) for p in value]
    out: dict[str, float] = {}
    for item in items:
        if len(item) != 2:
            raise ValueError("scattering_lengths entries are (atom name, b in fm)")
        name, b = item
        if not isinstance(name, str):
            raise ValueError(f"scattering_lengths: atom name must be a string, got {name!r}")
        up = name.strip().upper()
        if not _ATOM_NAME_RE.match(up):
            raise ValueError(f"scattering_lengths: {name!r} is not a DISCUS atom name (1-4 characters, "
                             "starting with a letter)")
        if up == "VOID":
            raise ValueError("scattering_lengths: VOID does not scatter and cannot get a scattering length")
        if up in out:
            raise ValueError(f"scattering_lengths: {up!r} given twice (DISCUS ignores upper/lower case)")
        if isinstance(b, (bool, np.bool_)) or not isinstance(b, (int, float, np.integer, np.floating)) \
                or not np.isfinite(float(b)):
            raise ValueError(f"scattering_lengths: b of {up} must be a finite number of fm, got {b!r}")
        out[up] = float(b)
    return tuple(sorted(out.items()))


def _normalize_delf(value) -> tuple[tuple[str, tuple[float, float]], ...]:
    items = value.items() if isinstance(value, dict) else tuple(value)
    out = []
    for item in items:
        name, fv = item
        if not isinstance(name, str) or not _ATOM_NAME_RE.match(name.strip().upper()):
            raise ValueError(f"delf atom name {name!r} is not a DISCUS atom name")
        fv = tuple(fv)
        if len(fv) != 2 or any(isinstance(x, (bool, np.bool_)) or not isinstance(x, (int, float, np.integer, np.floating))
                               or not np.isfinite(x) for x in fv):
            raise ValueError(f"delf for {name!r} must be (f', f'') numbers")
        out.append((name.strip().upper(), (float(fv[0]), float(fv[1]))))
    names = [n for n, _ in out]
    if len(set(names)) != len(names):
        raise ValueError("delf names one atom twice")
    return tuple(out)


def scat_commands(scattering_lengths) -> list[str]:
    """``scat NAME, 0,0,0,0,0,0,0,0, b`` lines (fourier menu). The eight form-factor coefficients are zero
    and the constant term is b; for neutrons DISCUS 6.23.01 uses only that term (fourier_sup.f90 ``form``),
    so b is the scattering length in fm. The setting is kept for a later ``powder`` run."""
    return [f"  scat {name}, 0,0,0,0,0,0,0,0, {b:.8g}" for name, b in scattering_lengths]


def axis_labels(grid: ReciprocalGrid, tol: float = 1e-12) -> tuple[str, str, str]:
    """``abs/ord/top`` labels for DISCUS: each active axis must have a nonzero hkl component along its
    label (6.23.01 refused an hhl plane labelled ``ord k``: 'Component of increment vector is zero').
    Prefers (h, k, l); otherwise the first permutation that satisfies the rule."""
    V = grid.axis_vectors
    scale = float(np.max(np.abs(V))) or 1.0
    active = grid.active_axes
    candidates = [_LABELS] + [p for p in permutations(_LABELS) if p != _LABELS]
    for labels in candidates:
        if all(abs(V[k, _LABELS.index(labels[k])]) > tol * scale for k in active):
            return tuple(labels)  # type: ignore[return-value]
    raise ValueError("no abs/ord/top labelling gives every active axis a nonzero labelled hkl component")


@dataclass(frozen=True)
class FourierConfig:
    """One DISCUS single-crystal Fourier calculation and its outputs.

    Parameters
    ----------
    grid: the reciprocal-space grid (``ll/lr/ul/tl/na/no/nt`` come from it).
    outputs: ``((value, file), ...)`` with value in ``int``, ``real``, ``imag``, ``3DPDF`` and an
        ``.h5`` file name; written in that order with ``form nexus``.
    aver_percent: ``set aver`` percentage, 0 (Bragg present) … 100.
    lots: :class:`Lots` or ``None`` (``lots off``).
    seed: :class:`Seed`; required when ``0 < aver_percent < 100`` or for random lots.
    temperature: ``'ignore'`` or ``'use'`` (``temp``). DISCUS 6.23.01 replaces B < 0.001 by 1.0 when
        reading, so ``'use'`` with B = 0 in the file gives B = 1 Å².
    symmetry: ``set symm:`` value, or ``None`` to keep the engine default (Friedel merging).
    radiation: ``'xray'`` (default) or ``'neutron'``.
    wavelength: a DISCUS symbol from :data:`WAVELENGTHS` (default ``'moa1'``, 0.70926 Å in DISCUS) or a
        number in Å (> 0). With ``disp off`` it changes only the file metadata, not the intensities.
    energy: the beam energy instead of a wavelength (``energy`` command): keV for X-rays and electrons
        (relativistic for electrons), meV for neutrons (``get_wave``). Leave ``wavelength`` at its default.
        Intensities on an hkl grid do not depend on the wavelength except through anomalous dispersion.
    dispersion: anomalous dispersion for X-rays (``disp anom``). DISCUS has internal f′, f″ only for the
        Kα1 symbols in :data:`ANOMALOUS_SYMBOLS`; for any other wavelength or an energy, every species needs
        its own values in ``delf``, otherwise it silently gets f′ = f″ = 0.
    delf: ``{atom name: (f′, f″)}`` X-ray anomalous terms in units of e, written as ``delf NAME,f′,f″``;
        needs ``dispersion=True``. Every name must be in the structure (DISCUS error -27 otherwise); species
        without an entry keep the internal values (0 for a non-Kα1 wavelength). Only atoms with DISCUS's
        internal form factors use them.
    scattering_lengths: own coherent neutron scattering lengths (fm) by DISCUS atom name, e.g.
        ``{"NI62": 10.3}`` for an isotope label; neutron runs only. Each becomes a
        ``scat NAME, 0,0,0,0,0,0,0,0, b`` line before ``show``. Give only names that occur in the
        structure (what DISCUS does with ``scat`` for an absent name was not checked).
        Empty (the default) leaves the macro unchanged.
    """
    grid: ReciprocalGrid
    outputs: tuple[tuple[str, str], ...]
    aver_percent: float = 0.0
    lots: Lots | None = None
    seed: Seed | None = None
    temperature: str = "ignore"
    symmetry: str | None = None
    wavelength: str | float = "moa1"
    radiation: str = "xray"
    labels: tuple[str, str, str] | None = field(default=None)
    scattering_lengths: Any = ()
    energy: float | None = None
    dispersion: bool = False
    delf: Any = ()

    def __post_init__(self) -> None:
        if not isinstance(self.grid, ReciprocalGrid):
            raise TypeError("grid must be a ReciprocalGrid")
        p = float(self.aver_percent)
        if not np.isfinite(p) or not 0.0 <= p <= 100.0:
            raise ValueError("aver_percent must be in [0, 100]")
        object.__setattr__(self, "aver_percent", p)
        if self.lots is not None and not isinstance(self.lots, Lots):
            raise TypeError("lots must be a Lots instance or None")
        if self.seed is not None and not isinstance(self.seed, Seed):
            raise TypeError("seed must be a Seed (8 values) or None")
        if 0.0 < p < 100.0 and self.seed is None:
            raise ValueError(f"set aver,{p:g} draws unit cells at random; give an explicit 8-value Seed")
        if self.lots is not None and self.lots.random and self.seed is None:
            raise ValueError("random lot origins need an explicit 8-value Seed")
        if self.temperature not in ("ignore", "use"):
            raise ValueError("temperature must be 'ignore' or 'use'")
        if self.symmetry is not None and self.symmetry not in SYMMETRY:
            raise ValueError(f"symmetry must be one of {SYMMETRY} or None")
        if self.radiation not in RADIATIONS:
            raise ValueError(f"radiation must be one of {RADIATIONS}")
        if self.wavelength is None and self.energy is None:
            raise ValueError("wavelength=None needs energy=")
        object.__setattr__(self, "wavelength", "moa1" if self.wavelength is None else _check_wavelength(self.wavelength))
        if self.energy is not None:
            e = self.energy
            if isinstance(e, (bool, np.bool_)) or not isinstance(e, (int, float, np.integer, np.floating)) \
                    or not np.isfinite(e) or e <= 0:
                raise ValueError(f"energy must be a positive number ({ENERGY_UNITS[self.radiation]}), got {e!r}")
            if self.wavelength not in ("moa1", None):
                raise ValueError("give either wavelength or energy, not both")
            object.__setattr__(self, "energy", float(e))
        if self.radiation == "electron" and self.energy is None and isinstance(self.wavelength, str):
            raise ValueError("electrons need energy= (keV) or a wavelength in Å; an X-ray tube symbol would "
                             "be recorded as the electron wavelength")
        if not isinstance(self.dispersion, bool):
            raise TypeError("dispersion must be True or False")
        delf = _normalize_delf(self.delf)
        if self.dispersion and self.radiation != "xray":
            raise ValueError("anomalous dispersion applies to X-rays only (DISCUS sets f' = f'' = 0 otherwise)")
        if delf and not self.dispersion:
            raise ValueError("delf values are used only with dispersion=True")
        if self.dispersion and not delf and (self.energy is not None or self.wavelength not in ANOMALOUS_SYMBOLS):
            raise ValueError("DISCUS 6.23.01 has internal f', f'' only for the K-alpha1 symbols "
                             f"{ANOMALOUS_SYMBOLS}; give delf for every species (any species without one gets 0)")
        object.__setattr__(self, "delf", delf)
        scat = normalize_scattering_lengths(self.scattering_lengths)
        if scat and self.radiation != "neutron":
            raise ValueError("scattering_lengths are neutron scattering lengths; use radiation='neutron' "
                             "(X-rays always use the element's form factor)")
        object.__setattr__(self, "scattering_lengths", scat)
        outs = tuple(tuple(o) for o in self.outputs)
        if not outs:
            raise ValueError("at least one output is required")
        names = set()
        for o in outs:
            if len(o) != 2:
                raise ValueError("each output is (value, file)")
            value, name = o
            if value not in OUTPUT_VALUES:
                hint = " ('amp' was not verified in 6.23.01)" if value == "amp" else ""
                raise ValueError(f"output value {value!r} not offered; use one of {OUTPUT_VALUES}{hint}")
            _check_name("output file", name, (".h5", ".hdf5"))
            if name in names:
                raise ValueError(f"output file {name!r} used twice")
            names.add(name)
            if self.lots is not None and value != "int":
                raise ValueError(f"with lots DISCUS writes only 'int' (6.23.01 refused {value!r}, error -60)")
        values = [o[0] for o in outs]
        if values.count("3DPDF") > 1 or ("3DPDF" in values and values[-1] != "3DPDF"):
            raise ValueError("write at most one 3DPDF output and write it last: the 3D-PDF replaces the intensity "
                             "in memory (dsi = dsi3d, then dsi3d is freed, graphic.f90:520-545)")
        object.__setattr__(self, "outputs", outs)
        labels = axis_labels(self.grid) if self.labels is None else tuple(self.labels)
        if sorted(labels) != sorted(_LABELS):
            raise ValueError("labels must be a permutation of ('h', 'k', 'l')")
        V = self.grid.axis_vectors
        for k in self.grid.active_axes:
            if V[k, _LABELS.index(labels[k])] == 0:
                raise ValueError(f"axis {k} is labelled {labels[k]!r} but has no {labels[k]} component; "
                                 "DISCUS refuses this ('Component of increment vector is zero')")
        object.__setattr__(self, "labels", labels)

    # -- facts derived from the configuration ---------------------------------------
    @property
    def bragg(self) -> str:
        """Label for the written intensity: 'bragg_subtracted' if ⟨F⟩ is subtracted, else 'bragg_present'."""
        return "bragg_subtracted" if self.aver_percent > 0 else "bragg_present"

    @property
    def wvle_argument(self) -> str:
        """The ``wvle`` argument: the symbol, or the wavelength in Å with up to 12 significant digits."""
        return self.wavelength if isinstance(self.wavelength, str) else f"{self.wavelength:.12g}"

    @property
    def labels_verified(self) -> bool:
        """True when the labels are (h, k, l) on an hkl-aligned grid — the case run in the probes."""
        return self.labels == _LABELS and self.grid.is_hkl_axis_aligned()

    def check_crystal(self, crystal_cells) -> None:
        """Refuse lots larger than the crystal (DISCUS error -101)."""
        if self.lots is None or crystal_cells is None:
            return
        cells = np.asarray(crystal_cells, int)
        if np.any(np.asarray(self.lots.size) > cells):
            raise ValueError(f"lot size {self.lots.size} exceeds the crystal {tuple(int(n) for n in cells)} "
                             "(DISCUS error -101)")

    def fourier_lines(self) -> list[str]:
        c = self.grid.discus_corners()
        wave = f"  energy {self.energy:.12g}" if self.energy is not None else f"  wvle {self.wvle_argument}"
        lines = ["fourier", f"  {_RADIATION_KEYWORD[self.radiation]}", wave,
                 f"  ll {_fmt(c['ll'])}", f"  lr {_fmt(c['lr'])}", f"  ul {_fmt(c['ul'])}", f"  tl {_fmt(c['tl'])}",
                 f"  na {c['na']}", f"  no {c['no']}", f"  nt {c['nt']}",
                 f"  abs {self.labels[0]}", f"  ord {self.labels[1]}", f"  top {self.labels[2]}",
                 f"  temp {self.temperature}", f"  disp {'anom' if self.dispersion else 'off'}",
                 f"  {self.lots.command() if self.lots else 'lots off'}",
                 "  set technique:turbo", f"  set aver,{self.aver_percent:.1f}"]
        if self.symmetry is not None:
            lines.append(f"  set symm:{self.symmetry}")
        lines += scat_commands(self.scattering_lengths)
        lines += [f"  delf {name},{fp:.12g},{fpp:.12g}" for name, (fp, fpp) in self.delf]
        lines += ["  show", "  run", "exit"]
        return lines

    def output_lines(self) -> list[str]:
        lines = ["output"]
        for value, name in self.outputs:
            lines += [f"  value {value}", "  form nexus", f"  outf {name}", "  run"]
        lines.append("exit")
        return lines

    def render(self, structure: StructureSource, *, title: str | None = None) -> str:
        """Complete DISCUS macro text (LF line endings), ready to inspect and stage."""
        if not isinstance(structure, StructureSource):
            raise TypeError("structure must be a StructureSource")
        if self.lots is not None:
            if structure.crystal_cells is None:
                raise ValueError("lots need the crystal size: give StructureSource(crystal_cells=...)")
            self.check_crystal(structure.crystal_cells)
        head = []
        if title is not None:
            if "\n" in title:
                raise ValueError("title must be one line")
            head.append(f"# {title}")
        body = head + ["discus"]
        if self.seed is not None:
            body.append(self.seed.command())
        body += structure.lines() + self.fourier_lines() + self.output_lines() + ["exit", "exit", ""]
        return "\n".join(body)

    # -- persistence / provenance ---------------------------------------------------
    def to_dict(self) -> dict[str, Any]:
        return {"schema": SCHEMA_NAME, "version": SCHEMA_VERSION, "grid": self.grid.to_dict(),
                "outputs": [list(o) for o in self.outputs], "aver_percent": self.aver_percent,
                "lots": None if self.lots is None else {"shape": self.lots.shape, "size": list(self.lots.size),
                                                       "count": self.lots.count, "periodic": self.lots.periodic},
                "seed": None if self.seed is None else list(self.seed.values),
                "temperature": self.temperature, "symmetry": self.symmetry, "wavelength": None if self.energy is not None else self.wavelength,
                "radiation": self.radiation, "labels": list(self.labels), "technique": "turbo", "engine": "DISCUS_SUITE 6.23.01",
                "bragg": self.bragg, "labels_verified": self.labels_verified,
                # only written when used, so descriptions of existing runs stay unchanged
                **({"scattering_lengths": {n: b for n, b in self.scattering_lengths}}
                   if self.scattering_lengths else {}),
                **({"energy": self.energy} if self.energy is not None else {}),
                **({"dispersion": True} if self.dispersion else {}),
                **({"delf": {n: list(v) for n, v in self.delf}} if self.delf else {})}

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> "FourierConfig":
        if d.get("schema") != SCHEMA_NAME or d.get("version") not in (1, SCHEMA_VERSION):
            raise ValueError(f"not a {SCHEMA_NAME} v1/v{SCHEMA_VERSION} description")
        # v1 had no 'radiation' key: every v1 configuration was X-ray.
        radiation = d["radiation"] if d["version"] >= 2 else "xray"
        lots = None if d["lots"] is None else Lots(d["lots"]["shape"], tuple(d["lots"]["size"]),
                                                   d["lots"]["count"], d["lots"]["periodic"])
        return cls(ReciprocalGrid.from_dict(d["grid"]), tuple(tuple(o) for o in d["outputs"]),
                   aver_percent=d["aver_percent"], lots=lots,
                   seed=None if d["seed"] is None else Seed(tuple(d["seed"])),
                   temperature=d["temperature"], symmetry=d["symmetry"], wavelength=d["wavelength"],
                   radiation=radiation, labels=tuple(d["labels"]),
                   scattering_lengths=d.get("scattering_lengths") or (), energy=d.get("energy"),
                   dispersion=bool(d.get("dispersion", False)),
                   delf={n: tuple(v) for n, v in (d.get("delf") or {}).items()})


def prepare_fourier_run(runner, config: FourierConfig, structure: StructureSource, structure_path,
                        *, title: str | None = None):
    """Stage a generated macro, the structure file and ``fourier_config.json`` with
    :meth:`Diffuse3DRunner.prepare`. The runner's ``output_filename`` must be one of the config's
    output files (it is the file the runner validates)."""
    structure_path = Path(structure_path)
    if structure_path.name != structure.file:
        raise ValueError(f"structure_path {structure_path.name!r} does not match the macro's {structure.file!r}")
    if runner.output_filename not in {name for _, name in config.outputs}:
        raise ValueError(f"runner.output_filename {runner.output_filename!r} is not an output of this config")
    text = config.render(structure, title=title)
    raw = structure_path.read_bytes()
    if b"\0" in raw[:4096]:
        raise ValueError(f"{structure_path} is not a text file")
    with tempfile.TemporaryDirectory() as tmp:
        macro = Path(tmp) / "run.mac"
        macro.write_bytes(text.encode("ascii"))
        # A Windows checkout (core.autocrlf) has CRLF; DISCUS would read the CR as part of each line.
        staged = Path(tmp) / structure.file
        staged.write_bytes(raw.replace(b"\r\n", b"\n"))
        cfg = Path(tmp) / "fourier_config.json"
        cfg.write_text(json.dumps({**config.to_dict(), "structure": {
            "kind": structure.kind, "file": structure.file, "crystal_cells": structure.crystal_cells,
            "staged_with_lf": raw != staged.read_bytes()}}, indent=1) + "\n", encoding="utf-8")
        return runner.prepare(macro, assets=(staged, cfg))
