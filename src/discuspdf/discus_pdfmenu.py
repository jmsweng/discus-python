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
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .discus_ops import Step
from .fourier_macro import RADIATIONS
from .macro import Macro
from ._validate import PLAIN_FILE_RE as _NAME, positive as _pos

__all__ = ["PdfCalc", "PdfData", "PdfFit", "PdfPartials", "read_pdf_menu_file", "read_pdf_partials", "run_pdf_fit"]



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
            self.settings(m)
            m.command("calc")
            m.command("save", "pdf", self.output)

    def settings(self, m: Macro, *, rmax: float | None = None, dr: float | None = None) -> None:
        """The ``set`` lines of this calculation inside an open ``pdf`` menu (``rmax``/``dr`` override the range;
        :class:`PdfPartials` and :class:`PdfFit` reuse them)."""
        if self.radiation == "neutron":
            m.command("set", "radiation", "neutron")
        else:
            m.command("set", "radiation", self.radiation, self.weight_s)
        # always three parameters: with two, pdf.f90:1206-1209 tests cpara(3)/(4) left from earlier lines,
        # and 'set boundary' also sets the exact flag, so it is given here ('3D' is just not '2D')
        m.command("set", "boundary", "periodic" if self.periodic else "crystal",
                  "exact" if self.exact else "3D")
        m.command("set", "range", self.rmax if rmax is None else rmax, self.dr if dr is None else dr)
        m.command("set", "qmax", self.qmax)
        if self.sphere is not None:
            m.command("set", "finite", "sphere", self.sphere)
        else:
            m.command("set", "finite", "periodic")
        m.command("set", "density", "auto" if self.density is None else self.density)
        m.command("set", "therm", "gauss" if self.thermal else "off")
        for k in ("delta", "gamma", "qalpha"):
            m.command("set", k, getattr(self, k))


def read_pdf_menu_file(path) -> tuple[np.ndarray, np.ndarray]:
    """(r, G) from a ``save pdf`` file (r, G and two constant columns per line), via :func:`discuspdf.pdf.read_pdf`."""
    from .pdf import read_pdf

    return read_pdf(path)


# ---------------------------------------------------------------------------------------------- partial PDFs
def _species(name) -> str:
    from .atom_names import check_discus_name

    return check_discus_name("species", name).upper()


@dataclass(frozen=True)
class PdfPartials(Step):
    """Partial PDFs G_ab(r) of the job's crystal, with the settings of ``model`` (a :class:`PdfCalc`, whose
    ``output`` is not written).

    DISCUS selects the pair types with ``isele``/``jsele``. Probed on the stock engine (2026-11-05,
    ``tools/probes/pdf_fit_probes.py``, rounds 1-3, binary and ternary crystals, residuals 1e-8):

    * same-type partial (a, a): G_aa = w_aa P_aa / (N dr r) - c_a^2 w_aa 4 pi r rho0, as expected
      (w = f_a f_b / <f>^2, P the pair count over both orders);
    * **cross partial selected directly is wrong**: ``isele a; jsele b`` counts half of the a-b pairs,
      G = w_ab P_ab / (2 N dr r) - 2 c_a c_b w_ab 4 pi r rho0, so the three partials of a binary do not add up to the
      total (issue #73);
    * selecting {a, b} on both sides counts the cross pairs in full, so the step writes that and the two same-type
      partials, and :func:`read_pdf_partials` returns G_ab = G_{a,b} - G_aa - G_bb (exact for the ternary probe).

    The selections are reset to all atoms at the end (``isele all``, ``jsele all``). Files: ``{prefix}_{A}_{A}.dat`` for a
    same-type selection and ``{prefix}_{A}_{B}_union.dat`` for {A, B} on both sides (the union, not G_AB; see
    :meth:`output_names`)."""

    model: PdfCalc
    pairs: tuple
    prefix: str = "partial"

    def __post_init__(self) -> None:
        if not isinstance(self.model, PdfCalc):
            raise TypeError("model must be a PdfCalc (its settings are used)")
        pairs = []
        for pair in tuple(self.pairs):
            if len(tuple(pair)) != 2:
                raise ValueError(f"pairs must be (A, B) tuples, got {pair!r}")
            a, b = sorted(_species(x) for x in pair)
            if (a, b) in pairs:
                raise ValueError(f"pair {pair!r} given twice")
            pairs.append((a, b))
        if not pairs:
            raise ValueError("give at least one pair")
        object.__setattr__(self, "pairs", tuple(pairs))
        if not isinstance(self.prefix, str) or not _NAME.match(self.prefix + "_x.dat"):
            raise ValueError(f"prefix must make plain file names, got {self.prefix!r}")

    def selections(self) -> list[tuple[str, ...]]:
        """The species sets selected on both sides: {a} for (a, a); {a}, {b} and {a, b} for (a, b)."""
        out: list[tuple[str, ...]] = []
        for a, b in self.pairs:
            for sel in (((a,),) if a == b else ((a,), (b,), (a, b))):
                if sel not in out:
                    out.append(sel)
        return out

    def file_of(self, sel: tuple[str, ...]) -> str:
        return f"{self.prefix}_{sel[0]}_{sel[0]}.dat" if len(sel) == 1 else f"{self.prefix}_{sel[0]}_{sel[1]}_union.dat"

    def output_names(self) -> tuple[str, ...]:
        return tuple(self.file_of(sel) for sel in self.selections())

    def apply(self, m: Macro) -> None:
        with m.menu("pdf"):
            self.model.settings(m)
            for sel in self.selections():
                m.command("idese", "all")
                m.command("isele", *sel)
                m.command("jdese", "all")
                m.command("jsele", *sel)
                m.command("calc")
                m.command("save", "pdf", self.file_of(sel))
            m.command("isele", "all")
            m.command("jsele", "all")


def read_pdf_partials(step: PdfPartials, run) -> dict:
    """``{(A, B): (r, G_AB)}`` from a run that contains ``step`` (cross partials combined, see :class:`PdfPartials`)."""
    g = {sel: read_pdf_menu_file(run.output(step.file_of(sel))) for sel in step.selections()}
    out = {}
    for a, b in step.pairs:
        if a == b:
            out[(a, b)] = g[(a,)]
        else:
            r, gu = g[(a, b)]
            out[(a, b)] = (r, gu - g[(a,)][1] - g[(b,)][1])
    return out


# ---------------------------------------------------------------------------------------------- fitting G(r)
@dataclass(frozen=True)
class PdfData:
    """A G(r) data file for ``pdf data``: four columns r, G, 0, sigma on a uniform grid r = dr, 2 dr, ..., rmax
    (``pdf_readdata``, pdf.f90:889-1047: the first r must equal the step, column 4 is the error, weight 1/sigma^2,
    sigma 0 means weight 1). Write one with :meth:`write` (which can also exclude the top of the range from the fit),
    or check an existing one with :meth:`from_file`."""

    file: str
    dr: float
    rmax: float
    points: int

    @classmethod
    def write(cls, path, r, g, sigma=None, *, fit_rmax: float | None = None) -> "PdfData":
        """Write r, G, 0, sigma. ``fit_rmax`` leaves the points above it in the file (DISCUS's model range ends at
        the data's rmax) but gives them sigma 1e6, i.e. weight 1e-12: ``set frange`` cannot lower the upper end of
        the fit window (pdf.f90:393-394 takes the larger of the user's and the file's rmax), so this is how the
        window is shortened. Recommended: rmax - 0.5 with Gaussian peaks (see :class:`PdfFit`)."""
        r = np.asarray(r, dtype=float)
        g = np.asarray(g, dtype=float)
        sig = np.zeros_like(r) if sigma is None else np.broadcast_to(np.asarray(sigma, dtype=float), r.shape)
        if r.ndim != 1 or r.shape != g.shape or len(r) < 3:
            raise ValueError("r and g must be 1-D arrays of the same length (>= 3)")
        if not (np.isfinite(r).all() and np.isfinite(g).all() and np.isfinite(sig).all()) or (sig < 0).any():
            raise ValueError("r, g and sigma must be finite, sigma >= 0")
        dr = (r[-1] - r[0]) / (len(r) - 1)
        if dr <= 0 or np.max(np.abs(np.diff(r) - dr)) > 1e-6 * max(1.0, dr) or abs(r[0] - dr) > 1e-6:
            raise ValueError("r must be the uniform grid dr, 2 dr, ..., rmax (pdf data needs the first r = dr); "
                             "interpolate the data onto it first")
        if fit_rmax is not None:
            fr = float(fit_rmax)
            if not r[0] < fr < r[-1]:
                raise ValueError(f"fit_rmax must lie inside the data ({r[0]} .. {r[-1]})")
            sig = np.where(r > fr + 1e-9, 1e6, sig)
        path = Path(path)
        with path.open("w", encoding="ascii") as f:
            for x, y, s in zip(r, g, sig):
                f.write(f"{x:.6f} {y:.10g} 0.0 {s:.10g}\n")
        return cls(path.name, float(dr), float(r[-1]), len(r))

    @classmethod
    def from_file(cls, path) -> "PdfData":
        rows = [l.split() for l in Path(path).read_text(encoding="utf-8", errors="replace").splitlines()
                if l.strip() and not l.lstrip().startswith("#")]
        if len(rows) < 3 or any(len(row) < 4 for row in rows):
            raise ValueError(f"{path}: pdf data needs four columns r, G, dr, sigma")
        r = np.array([float(row[0].replace("D", "E")) for row in rows])
        dr = (r[-1] - r[0]) / (len(r) - 1)
        if dr <= 0 or np.max(np.abs(np.diff(r) - dr)) > 1e-5 or abs(r[0] - dr) > 1e-6:
            raise ValueError(f"{path}: r must be the uniform grid dr, 2 dr, ..., rmax")
        return cls(Path(path).name, float(dr), float(r[-1]), len(r))

    def __post_init__(self) -> None:
        if not isinstance(self.file, str) or not _NAME.match(self.file):
            raise ValueError(f"file must be a plain file name, got {self.file!r}")
        _pos("dr", self.dr), _pos("rmax", self.rmax)


@dataclass(frozen=True)
class PdfFit(Step):
    """A reverse Monte Carlo fit of the job's crystal to G(r) data in the ``pdf`` menu (``data`` + ``run``).

    ``model`` (a :class:`PdfCalc`) gives the calculation settings (radiation, weights, qmax, boundary, density,
    Gaussian peaks); its ``rmax``/``dr`` are replaced by the data's and its ``output`` is the file the fitted G(r)
    is saved to (``save pdf`` after ``run``, from ``fit_rmin`` to the data's rmax). ``moves``, ``select``, ``tries``,
    ``sigma``, ``display``, ``fit_scale``/``scale`` and ``min_distances`` are those of
    :class:`~discuspdf.discus_rmc.RmcFit` and are checked and written by the same code (``pdf set`` forwards them to
    the RMC settings). ``fit_rmin`` is the lower end of the fit window (``set frange``).

    Facts (stock engine, probed 2026-11-05 with ``tools/probes/pdf_fit_probes.py``, and the source):

    * after ``reset all`` (or in a fresh session) ``data`` fails with "Too many points in PDF": the arrays are
      sized to exactly the number of lines and the read loop needs one more (pdf.f90:944-981). The step therefore
      runs one ``calc`` over ``rmax + 1`` first, which sizes them;
    * ``set frange`` cannot lower the upper end of the window: pdf_setup takes the larger of the user's and the
      file's rmax (pdf.f90:393-394). To leave out the top of the range, give those points a huge sigma
      (``PdfData.write(..., fit_rmax=...)``). Near rmax every Gaussian lacks its partners beyond rmax, so the model
      falls off within about 0.5 A of rmax: use ``fit_rmax = rmax - 0.5`` with Gaussian peaks;
    * the log has the RMC format (:func:`discuspdf.discus_rmc.parse_rmc_log`); a swap fit of random Cu0.75Au0.25
      (108 atoms, 2000 tries) to the L1_2 G(r) reached alpha1 = -0.318;
    * the G(r) saved after ``run`` (DISCUS's per-pair route) and ``calc`` on the fitted crystal differ by 8e-4 of
      max |G| (two code paths);
    * ``run`` stores the fitted scale, which a later ``calc`` in the session would apply (pdf.f90:187, :2015), and
      minimum distances persist from earlier RMC settings; the step resets the RMC settings first (``rmc``,
      ``reset``) and sets the scale back to fixed 1 after saving.
    """

    model: PdfCalc
    data: PdfData
    moves: tuple
    select: tuple | str
    tries: int
    sigma: float = 0.01
    display: int = 500
    fit_scale: bool = True
    scale: float = 1.0
    min_distances: dict = field(default_factory=dict)
    fit_rmin: float | None = None

    uses_random_numbers = True

    def __post_init__(self) -> None:
        from .discus_rmc import _real, check_fit_settings

        if not isinstance(self.model, PdfCalc):
            raise TypeError("model must be a PdfCalc (its settings are used)")
        if not isinstance(self.data, PdfData):
            raise TypeError("data must be a PdfData (PdfData.write or PdfData.from_file)")
        if self.data.file == self.model.output:
            raise ValueError("the model's output file must differ from the data file")
        smooth = self.model.thermal or self.model.qalpha > 0 or self.model.qmax > 0
        if self.data.dr > 0.001 + 1e-12 and not smooth:
            raise ValueError("data step > 0.001 A needs thermal, qalpha or qmax in the model (see PdfCalc)")
        check_fit_settings(self)
        if self.fit_rmin is not None:
            v = float(_real("fit_rmin", self.fit_rmin))
            if not self.data.dr - 1e-9 <= v < self.data.rmax:
                raise ValueError(f"fit_rmin must lie inside the data ({self.data.dr} .. {self.data.rmax})")
            object.__setattr__(self, "fit_rmin", v)

    def input_names(self) -> tuple[str, ...]:
        return (self.data.file,)

    def output_names(self) -> tuple[str, ...]:
        return (self.model.output,)

    def check_inputs(self, paths) -> None:
        found = PdfData.from_file(paths[self.data.file])
        if abs(found.dr - self.data.dr) > 1e-6 or abs(found.rmax - self.data.rmax) > 1e-6:
            raise ValueError(f"{self.data.file}: grid dr {found.dr}, rmax {found.rmax}; the step expects "
                             f"{self.data.dr}, {self.data.rmax}")

    def apply(self, m: Macro) -> None:
        from .discus_rmc import write_fit_settings

        with m.menu("rmc"):
            m.command("reset")  # clears selections, move probabilities and minimum distances left in the session
        with m.menu("pdf"):
            self.model.settings(m, rmax=round(self.data.rmax + 1.0, 6), dr=self.data.dr)
            m.command("calc")  # sizes the arrays so that 'data' can read the file (see the class notes)
            m.command("data", self.data.file)
            if self.fit_rmin is not None:
                m.command("set", "frange", self.fit_rmin, self.data.rmax)
            write_fit_settings(m, self)
            m.command("run")
            m.command("save", "pdf", self.model.output)
            m.command("set", "scal", "off", 1)  # 'run' keeps the fitted scale for later 'calc's (pdf.f90:2015, :187)


def run_pdf_fit(job, runner, inputs=None, **kw):
    """Run a job containing one :class:`PdfFit` and parse its log (an :class:`~discuspdf.discus_rmc.RmcOutcome`).
    The data file is checked against the step first; RuntimeError when the log has no final block."""
    from .discus_rmc import run_logged_fit

    return run_logged_fit(job, runner, inputs, PdfFit, **kw)
