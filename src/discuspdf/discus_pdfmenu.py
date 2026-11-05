"""Real-space PDFs with DISCUS's ``pdf`` menu (``calc``), as a :class:`~discuspdf.discus_ops.Step`.

This is DISCUS's direct route: a histogram of pair distances in the crystal, weighted and normalised,
with an optional sinc termination for a finite Qmax. (The powder route, :mod:`discuspdf.pdf`, goes
through the Debye powder pattern instead.) Facts are from the 6.23.01 source, see
``docs/DISCUS_PDFMENU_EXP2PDF_FACTS.md`` (agent-written; the items used here were checked):

* DISCUS works on an internal step dr_i = dr / NINT(dr/0.001) (0.001 Å or finer, ``set range``,
  pdf.f90:1536-1546) and ``save pdf`` writes every NINT(dr/0.001)-th point, r = dr, 2 dr, … (discus_s_pdf_file.f90:81-88).
  Without Gaussian peaks or termination the file is therefore a 0.001-Å histogram *sampled*, not averaged, at
  the output step, so the step refuses dr > 0.001 unless ``thermal``, ``qalpha`` or ``qmax`` smooth G.
* G(r_i) = pdf_scale · H_i / (N · dr_i · r_i) − 4π r_i ρ0 (``pdf_convert``, :2654-2717), r_i = i·dr_i,
  H_i = weighted ordered pair count (over periodic images) in bin i = INT((d + dr_i/2)/dr_i) for d ≤ rmax
  (:3015-3028), N the number of real atoms, ρ0 = N / (V_cell · n_cells) unless given (:2656-2665).
* Weights w_ij = f_i f_j / ⟨f⟩² · occ_i occ_j (:445-485), X-ray f from the Waasmaier 5-Gaussian table (the
  powder menu uses the International 4-Gaussian table) at s = sinθ/λ = ``weight_s`` (the number given to
  ``set radiation`` enters exp(−b h2) as h2 = value², i.e. it is s, not Q, whatever the help says).
  **Pitfall fixed here:** the start value of the radiation is X-ray but the flag that adds the Gaussian
  terms of the X-ray form factor (``pdf_lxray``) starts false, so without ``set radiation`` DISCUS weights by
  the constant term c of the table only (``form``, fourier_sup.f90:2947-2971). The step always writes
  ``set radiation``.
* ``qmax`` > 0 convolves G with sin(Qmax r)/r (kernel truncated at ±2048 steps, :535-545, :2749-2755);
  ``qmax=0`` switches it off (DISCUS's start value is 30).
* ``set qsig`` (Q-resolution damping) is not offered: it multiplies by ``pdf_exp``, which is filled only when
  the Gaussian thermal broadening is on (:552-562, :2723-2733); probed in ``tools/probes``.
* ``save pdf`` writes r, G and two constant columns (format F9.4,3X,F21.10, discus_s_pdf_file.f90:94).
"""
from __future__ import annotations

import math
import numbers
from dataclasses import dataclass

import numpy as np

from .discus_ops import Step
from .fourier_macro import RADIATIONS
from .macro import Macro
from ._validate import PLAIN_FILE_RE as _NAME, positive as _pos

__all__ = ["PdfCalc", "read_pdf_menu_file"]



@dataclass(frozen=True)
class PdfCalc(Step):
    """``pdf`` … ``calc`` and ``save pdf, <output>`` on the job's crystal.

    Parameters
    ----------
    output: file name for ``save pdf`` (declared as a run output).
    rmax, dr: ``set range, rmax, dr`` (Å); see the module notes on the internal step.
    radiation: ``xray``, ``neutron`` or ``electron`` (always written; see the module notes).
    weight_s: s = sinθ/λ (Å⁻¹; Q = 4π s) at which X-ray/electron form factors are taken for the weights.
    qmax: termination (``set qmax``); 0 = none. DISCUS's start value is 30.
    periodic: ``set boundary, periodic`` (images through the crystal's periodic boundaries) or ``crystal``.
    exact: the all-atoms path instead of the unit-cell path (``set boundary, …, exact``).
    sphere: ``set finite, sphere, D`` (a spherical particle of diameter D Å) instead of a periodic background.
    density: ρ0 in atoms/Å³, or None for DISCUS's own N/V (``set density``).
    thermal: ``set therm, gauss`` — Gaussian peaks from the atoms' B values; ``delta``/``gamma`` (correlated
        motion, ``set delta``/``set gamma``) then sharpen them. ``qalpha`` (``set qalpha``) adds a width ∝ r and
        broadens on its own. All three are always written (0 when unused), because DISCUS keeps them in a session.
    """
    output: str
    rmax: float = 20.0
    dr: float = 0.001
    radiation: str = "xray"
    weight_s: float = 0.0
    qmax: float = 0.0
    periodic: bool = True
    exact: bool = False
    sphere: float | None = None
    density: float | None = None
    thermal: bool = False
    delta: float = 0.0
    gamma: float = 0.0
    qalpha: float = 0.0

    def __post_init__(self) -> None:
        if not isinstance(self.output, str) or not _NAME.match(self.output):
            raise ValueError(f"output must be a plain file name, got {self.output!r}")
        _pos("rmax", self.rmax), _pos("dr", self.dr)
        if self.dr >= self.rmax:
            raise ValueError("dr must be smaller than rmax")
        if self.radiation not in RADIATIONS:
            raise ValueError(f"radiation must be one of {RADIATIONS}")
        _pos("weight_s", self.weight_s, allow_zero=True)
        if self.weight_s and self.radiation == "neutron":
            raise ValueError("weight_s applies to X-ray and electron form factors only")
        _pos("qmax", self.qmax, allow_zero=True)
        for k in ("periodic", "exact", "thermal"):
            if not isinstance(getattr(self, k), bool):
                raise TypeError(f"{k} must be True or False")
        if self.sphere is not None:
            _pos("sphere", self.sphere)
        if self.density is not None:
            _pos("density", self.density)
        for k in ("delta", "gamma", "qalpha"):
            v = getattr(self, k)
            if isinstance(v, bool) or not isinstance(v, numbers.Real) or not math.isfinite(float(v)):
                raise ValueError(f"{k} must be a finite number")
            if v and k != "qalpha" and not self.thermal:
                raise ValueError(f"{k} acts only with thermal=True (Gaussian peaks)")
        smooth = self.thermal or self.qalpha > 0 or self.qmax > 0
        if self.dr > 0.001 + 1e-12 and not smooth:
            raise ValueError("dr > 0.001 without thermal/qalpha/qmax: DISCUS would write a 0.001-A histogram sampled "
                             "every dr (not averaged); use dr <= 0.001 or smooth G")

    def output_names(self) -> tuple[str, ...]:
        return (self.output,)

    def apply(self, m: Macro) -> None:
        with m.menu("pdf"):
            if self.radiation == "neutron":
                m.command("set", "radiation", "neutron")
            else:
                m.command("set", "radiation", self.radiation, self.weight_s)
            # always three parameters: with two, pdf.f90:1206-1209 tests cpara(3)/(4) left from earlier lines,
            # and 'set boundary' also sets the exact flag, so it is given here ('3D' is just not '2D')
            m.command("set", "boundary", "periodic" if self.periodic else "crystal",
                      "exact" if self.exact else "3D")
            m.command("set", "range", self.rmax, self.dr)
            m.command("set", "qmax", self.qmax)
            if self.sphere is not None:
                m.command("set", "finite", "sphere", self.sphere)
            else:
                m.command("set", "finite", "periodic")
            m.command("set", "density", "auto" if self.density is None else self.density)
            m.command("set", "therm", "gauss" if self.thermal else "off")
            for k in ("delta", "gamma", "qalpha"):
                m.command("set", k, getattr(self, k))
            m.command("calc")
            m.command("save", "pdf", self.output)


def read_pdf_menu_file(path) -> tuple[np.ndarray, np.ndarray]:
    """(r, G) from a ``save pdf`` file (r, G and two constant columns per line), via :func:`discuspdf.pdf.read_pdf`."""
    from .pdf import read_pdf

    return read_pdf(path)
