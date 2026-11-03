"""Minimal, isolated DISCUS 1D PDF simulation runner.

This adapter uses documented legacy DISCUS PDF commands. Validate once against
YOUR installed version of discus_suite before running a large sweep.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Sequence
import hashlib
import json
import math
import re
import shutil
import subprocess
import uuid
import warnings

import numpy as np

from .engine import current_engine, Engine, as_engine
from .fourier_macro import normalize_scattering_lengths, scat_commands


# DISCUS_SUITE 6.23.01 cannot use Q above about 12.3/wavelength (Å⁻¹); larger qmax values give the same
# G(r). Measured: G(r) stops changing between Qmax 8.0 and 8.1 at λ = 1.5406 Å and between 17.0 and 17.5
# at λ = 0.7107 Å. Above Q ≈ 8 Å⁻¹ the DISCUS F(Q) also depends on the wavelength by up to several percent
# (debye mode, X-ray and neutron). See docs/VERIFIED_DISCUS_MACRO.md, "Qmax limit and wavelength".
QMAX_WAVELENGTH_WARN = 12.0      # warn when qmax * wavelength exceeds this (margin below the measured limit)
QMAX_WAVELENGTH_CEILING = 12.3   # approximate measured limit, used only in the warning text


class QmaxLimitWarning(UserWarning):
    """qmax is above what DISCUS can use at this wavelength; the PDF is computed unchanged."""


SET_DELTA_ISSUE = "https://github.com/jmsweng/discus-python/issues/65"
_SET_DELTA_PROBLEM = (
    "DISCUS_SUITE 6.23.01 aborts on most Q grids when 'set delta' > 0 (heap corruption in "
    "powder_convolute), and a run that does not abort may still be wrong. See "
    f"{SET_DELTA_ISSUE} and docs/VERIFIED_DISCUS_MACRO.md."
)


class SetDeltaError(ValueError):
    """qsig > 0 (DISCUS 'set delta') was requested without allow_qsig=True."""


class SetDeltaWarning(UserWarning):
    """qsig > 0 runs 'set delta', which DISCUS 6.23.01 cannot do reliably (allow_qsig=True was given)."""


def _caller_stacklevel() -> int:
    """stacklevel for a warning issued in PDFConfig.__post_init__ that points at the user's code.

    Skips the dataclass-generated ``__init__`` and ``dataclasses.replace()`` frames, so the
    warning is attributed to the line that created the config, however it was created.
    """
    import dataclasses
    import sys
    frame = sys._getframe(2)  # 0: this helper, 1: __post_init__, 2: its caller
    level = 2
    while frame is not None and (
        frame.f_code.co_filename == dataclasses.__file__
        or (frame.f_code.co_name == "__init__" and frame.f_code.co_filename.startswith("<"))
    ):
        frame = frame.f_back
        level += 1
    return level


@dataclass(frozen=True)
class PDFConfig:
    """1D powder PDF calculation from a DISCUS .cll unit cell or a .stru supercell.

    ``structure`` is a ``.cll`` unit cell, repeated ``supercell`` times by ``read cell``, or a
    ``.stru`` configuration (for example a disordered supercell), read whole by ``read stru``; a
    ``.stru`` needs ``supercell=(1, 1, 1)``. ``.cll`` macros are unchanged by the ``.stru`` option.

    ``scattering_lengths`` gives atom names their own coherent neutron scattering length (fm), e.g.
    ``{"NI62": 10.3}`` for an isotope label (neutron only). They are set with ``scat`` in a
    ``fourier`` block before ``powder``; DISCUS keeps them for the powder run. Empty (the default)
    leaves the macro unchanged.

    All Q values are inverse angstroms. qsig is the DISCUS powder
    "set delta" convolution width; it is not assumed to be a universal
    instrument-resolution parameter. qmax must be explicitly positive:
    omitting "set qmax" is not a verified way to disable termination.

    qsig > 0 raises :class:`SetDeltaError`, because DISCUS 6.23.01 aborts on most Q grids
    with 'set delta' > 0 (issue #65). Pass ``allow_qsig=True`` to run it anyway, e.g. with a
    DISCUS version in which this is fixed; a :class:`SetDeltaWarning` is then issued.
    """
    structure: Path | str
    supercell: tuple[int, int, int] = (5, 5, 5)
    rmax: float = 15.0
    dr: float = 0.02
    qmax: float = 25.0
    qsig: float = 0.0
    radiation: str = "xray"
    qmin: float = 0.5
    dq: float = 0.01
    wavelength: float | None = None
    allow_qsig: bool = False
    scattering_lengths: Any = ()

    def __post_init__(self) -> None:
        if len(self.supercell) != 3 or any(type(n) is not int or n < 1 for n in self.supercell):
            raise ValueError("supercell must contain three positive integers")
        for name in ("rmax", "dr", "qmax", "qsig", "qmin", "dq"):
            v = getattr(self, name)
            if isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v):
                raise ValueError(f"{name} must be a finite number")
        if not (0 < self.dr < self.rmax):
            raise ValueError("Require 0 < dr < rmax")
        if not math.isclose(self.rmax / self.dr, round(self.rmax / self.dr), abs_tol=1e-7):
            raise ValueError("rmax must be a whole number of dr intervals")
        if not (0 < self.qmin < self.qmax) or self.dq <= 0:
            raise ValueError("Require 0 < qmin < qmax and dq > 0 (qmax=0 is not supported)")
        if self.qsig < 0:
            raise ValueError("qsig must be >= 0")
        if not isinstance(self.allow_qsig, bool):
            raise ValueError("allow_qsig must be True or False")
        if self.qsig > 0:
            if not self.allow_qsig:
                raise SetDeltaError(
                    f"qsig = {self.qsig:g} needs 'set delta', and {_SET_DELTA_PROBLEM} "
                    "Use qsig = 0, or pass allow_qsig=True to run it anyway (for example with a fixed DISCUS).")
            warnings.warn(f"qsig = {self.qsig:g} with allow_qsig=True: {_SET_DELTA_PROBLEM}",
                          SetDeltaWarning, stacklevel=_caller_stacklevel())
        if self.wavelength is not None:
            if isinstance(self.wavelength, bool) or not isinstance(self.wavelength, (int, float)) or not math.isfinite(self.wavelength) or self.wavelength <= 0:
                raise ValueError("wavelength must be finite and > 0")
        if self.radiation not in {"xray", "neutron"}:
            raise ValueError("radiation must be 'xray' or 'neutron'")
        wavelength = default_wavelength(self)
        if self.qmax * wavelength > QMAX_WAVELENGTH_WARN:
            ceiling = QMAX_WAVELENGTH_CEILING / wavelength
            warnings.warn(
                f"qmax = {self.qmax:g} Å⁻¹ is above what DISCUS 6.23.01 can use at wavelength {wavelength:g} Å "
                f"(about {ceiling:.1f} Å⁻¹). The calculation runs unchanged, but its G(r) is the same as for "
                f"Qmax ≈ {ceiling:.1f} Å⁻¹. A shorter explicit wavelength raises the limit, but above Q ≈ 8 Å⁻¹ "
                "the DISCUS result then depends on the wavelength by up to several percent "
                "(see docs/VERIFIED_DISCUS_MACRO.md, 'Qmax limit and wavelength').",
                QmaxLimitWarning, stacklevel=_caller_stacklevel())
        suffix = Path(self.structure).suffix.lower()
        if suffix not in (".cll", ".stru"):
            raise ValueError("structure must be a DISCUS .cll unit cell or a .stru supercell")
        if suffix == ".stru" and tuple(self.supercell) != (1, 1, 1):
            raise ValueError("a .stru is read whole (read stru) and cannot be repeated; use supercell=(1, 1, 1)")
        scat = normalize_scattering_lengths(self.scattering_lengths)
        if scat and self.radiation != "neutron":
            raise ValueError("scattering_lengths are neutron scattering lengths; use radiation='neutron'")
        object.__setattr__(self, "scattering_lengths", scat)

    @property
    def structure_kind(self) -> str:
        """'cell' for a .cll unit cell, 'stru' for a .stru supercell."""
        return "stru" if Path(self.structure).suffix.lower() == ".stru" else "cell"

    @property
    def staged_name(self) -> str:
        """Name of the structure file inside the run directory (and in the macro)."""
        return "model.stru" if self.structure_kind == "stru" else "model.cll"

def default_wavelength(config: PDFConfig) -> float:
    """Wavelength (Å) written to the macro: config.wavelength, else 1.5406 (xray) or 1.5400 (neutron)."""
    if config.wavelength is not None:
        return float(config.wavelength)
    return 1.5406 if config.radiation == "xray" else 1.5400


def render_pdf_macro(config: PDFConfig) -> str:
    """Compile typed parameters to a DISCUS macro using a staged local file.

    Uses the DISCUS v6.x powder->PDF workflow, verified against the installed
    DISCUS_SUITE 6.23.01 (see docs/VERIFIED_DISCUS_MACRO.md):

      discus
      read
      cell <file>,<nx>,<ny>,<nz>
      powder
      reset
      <radiation>
      set axis,q
      set calc,debye
      set qmin,<qmin>
      set qmax,<qmax>
      set dq,<dq>
      set wvle,<wavelength>
      run
      exit
      output
      outf result.pdf
      form pdf, r, 0.0, <rmax>, <dr>
      value PDF
      run
      exit
      exit

    Notes:
    * In v6 the 1D powder PDF is computed from the powder pattern, so the
      calculation is driven from the ``powder`` menu (``set calc,debye``) and
      the PDF is written from the ``output`` menu (``form pdf, r, ...``).
    * ``qsig`` maps to the powder Gaussian convolution width (``set delta``);
      a value of 0.0 disables the convolution.
    * qmax=0 is rejected because the omission semantics are not verified.
    * ``result.pdf`` is DISCUS's raw space-delimited text (not an Adobe PDF).
      It stays in the macro because the archived runs and probe evidence use
      exactly this macro. ``PDFRunner`` converts it to ``result.csv``, the
      run's output (see ``PDF_CSV_NAME``); do not change that back to .pdf.
    """
    nx, ny, nz = config.supercell
    # Defaults are explicit in PDFConfig and can be overridden per run.
    wavelength = default_wavelength(config)
    read = "stru model.stru" if config.structure_kind == "stru" else f"cell model.cll,{nx},{ny},{nz}"
    scat = ["fourier", *scat_commands(config.scattering_lengths), "exit"] if config.scattering_lengths else []
    lines = [
        "# Generated by discuspdf: 1D atomic PDF calculation (DISCUS v6.x)",
        "discus",
        "read",
        read,
        *scat,
        "powder",
        "reset",
        config.radiation,
        "set axis,q",
        "set calc,debye",
        f"set qmin,{config.qmin:.12g}",
    ]
    lines += [
        f"set qmax,{config.qmax:.12g}",
        f"set dq,{config.dq:.12g}",
        f"set delta,{config.qsig:.12g}",
        f"set wvle,{wavelength:.12g}",
        "run",
        "exit",  # back to DISCUS main menu
        "output",
        "outf result.pdf",
        f"form pdf, r, 0.0, {config.rmax:.12g}, {config.dr:.12g}",
        "value PDF",
        "run",
        "exit",  # output menu
        "exit",  # DISCUS section
        "",
    ]
    return "\n".join(lines)


def read_pdf(path: Path | str) -> tuple[np.ndarray, np.ndarray]:
    """Read a DISCUS ASCII PDF file (r, G(r), possibly other columns)."""
    rows: list[tuple[float, float]] = []
    with Path(path).open(encoding="utf-8", errors="replace") as stream:
        for line in stream:
            stripped = line.strip()
            if not stripped or stripped.startswith(("#", "!", "%", ";")):
                continue
            parts = stripped.replace("D+", "E+").replace("D-", "E-").split()
            if len(parts) < 2:
                continue
            try:
                r, gr = float(parts[0]), float(parts[1])
            except ValueError:
                continue  # Text header
            if math.isfinite(r) and math.isfinite(gr):
                rows.append((r, gr))
    if len(rows) < 2:
        raise ValueError(f"No usable (r, G(r)) series found in {path}")
    arr = np.asarray(rows, dtype=float)
    if not np.all(np.diff(arr[:, 0]) > 0):
        raise ValueError("PDF r coordinates must be strictly increasing")
    return arr[:, 0], arr[:, 1]


# The 1D PDF output of a run is a comma-separated file. This is the owner's decision
# (2026-11-02): do not rename it back to result.pdf. DISCUS itself still writes its raw
# space-delimited text to RAW_PDF_NAME, which is kept next to it for provenance.
PDF_CSV_NAME = "result.csv"
RAW_PDF_NAME = "result.pdf"
PDF_CSV_HEADER = "r (Angstrom),G(r)"


def write_pdf_csv(path: Path | str, r: np.ndarray, gr: np.ndarray) -> Path:
    """Write (r, G(r)) as CSV: header ``r (Angstrom),G(r)``, then one row per point.

    Values are written with ``repr`` so they read back as exactly the same floats.
    """
    r = np.asarray(r, dtype=float)
    gr = np.asarray(gr, dtype=float)
    if r.shape != gr.shape or r.ndim != 1:
        raise ValueError("r and G(r) must be 1D arrays of the same length")
    lines = [PDF_CSV_HEADER] + [f"{a!r},{b!r}" for a, b in zip(r.tolist(), gr.tolist())]
    path = Path(path)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    return path


def read_pdf_csv(path: Path | str) -> tuple[np.ndarray, np.ndarray]:
    """Read a CSV written by ``write_pdf_csv`` (or any ``r,G(r)`` CSV with one header line)."""
    with Path(path).open(encoding="utf-8") as stream:
        header = stream.readline().strip()
        if header != PDF_CSV_HEADER:
            raise ValueError(f"{path}: expected the header {PDF_CSV_HEADER!r}, found {header!r}")
        arr = np.loadtxt(stream, delimiter=",", ndmin=2)
    if arr.shape[0] < 2 or arr.shape[1] != 2:
        raise ValueError(f"No usable (r, G(r)) series found in {path}")
    if not np.all(np.diff(arr[:, 0]) > 0):
        raise ValueError("PDF r coordinates must be strictly increasing")
    return arr[:, 0], arr[:, 1]


@dataclass(frozen=True)
class PDFRun:
    directory: Path

    @property
    def output(self) -> Path:
        """The run's PDF output, ``<directory>/result.csv`` (comma-separated, with a header)."""
        return self.directory / PDF_CSV_NAME

    @property
    def raw_output(self) -> Path:
        """The file DISCUS wrote, ``<directory>/result.pdf`` (space-delimited text, not Adobe PDF)."""
        return self.directory / RAW_PDF_NAME

    def load(self) -> tuple[np.ndarray, np.ndarray]:
        return read_pdf_csv(self.output)

    def plot(self, *, ax=None, label=None):
        import matplotlib.pyplot as plt
        r, gr = self.load()
        if ax is None:
            _, ax = plt.subplots(figsize=(8, 4))
        ax.plot(r, gr, label=label)
        ax.set(xlabel=r"$r$ (Å)", ylabel=r"$G(r)$")
        if label is not None:
            ax.legend()
        return ax


# The 1D runner's engine-error markers (narrower than engine.ENGINE_ERROR_PATTERN; unchanged).
_PDF_ENGINE_ERROR = re.compile(r"(?im)\b(?:unknown command|fatal error|segmentation fault|error\s*:)")


class PDFRunner:
    def __init__(
        self,
        runs_directory: Path | str = "runs",
        launcher: Sequence[str] = ("discus_suite", "-noautorun", "-macro"),
        timeout_seconds: int = 3600,
        *,
        engine: "Engine | None" = None,
    ) -> None:
        """launcher is the executable and args preceding 'run.mac'; ``engine`` (a
        :class:`discuspdf.engine.Engine`) replaces ``launcher`` and ``timeout_seconds``."""
        self.runs_directory = Path(runs_directory)
        if engine is None:
            if not tuple(launcher):
                raise ValueError("launcher must not be empty")
            if timeout_seconds <= 0:
                raise ValueError("timeout_seconds must be positive")
        self.engine = as_engine(engine, launcher, timeout_seconds)
        self.launcher = self.engine.launcher
        self.timeout_seconds = self.engine.timeout_seconds

    def prepare(self, config: PDFConfig) -> PDFRun:
        """Create inspectable files without executing DISCUS."""
        source = Path(config.structure).expanduser().resolve()
        if not source.is_file():
            raise FileNotFoundError(source)
        self.runs_directory.mkdir(parents=True, exist_ok=True)
        run_dir = self.runs_directory / (
            datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "_" + uuid.uuid4().hex[:8]
        )
        run_dir.mkdir(exist_ok=False)
        raw = source.read_bytes()
        staged = run_dir / config.staged_name
        # A Windows checkout (core.autocrlf) gives CRLF; DISCUS would read the CR as part of each line.
        if config.structure_kind == "stru":
            staged.write_bytes(raw.replace(b"\r\n", b"\n"))
        else:
            shutil.copy2(source, staged)   # unchanged: archived .cll runs record this file's hash
        (run_dir / "run.mac").write_text(render_pdf_macro(config), encoding="utf-8")
        params = asdict(config)
        if not params["scattering_lengths"]:
            del params["scattering_lengths"]   # descriptions of plain runs stay as before
        else:
            params["scattering_lengths"] = [list(p) for p in params["scattering_lengths"]]
        params["structure"] = str(source)
        params["input_sha256"] = hashlib.sha256(staged.read_bytes()).hexdigest()
        # The exact staged macro is the exact executed macro: record its hash.
        params["macro_sha256"] = hashlib.sha256((run_dir / "run.mac").read_bytes()).hexdigest()
        params["generated_utc"] = datetime.now(timezone.utc).isoformat()
        params["launcher"] = [*self.launcher, "run.mac"]
        (run_dir / "config.json").write_text(json.dumps(params, indent=2), encoding="utf-8")
        (run_dir / "status.json").write_text('{"status": "prepared"}\n', encoding="utf-8")
        return PDFRun(run_dir)

    def execute(self, prepared: PDFRun) -> PDFRun:
        """Execute the exact macro that was inspected via prepare().

        Status transitions: ``prepared`` → ``running`` → ``complete`` /
        ``failed``.  Any preflight failure (invalid config, tampered staged
        files, missing macro) is recorded as ``failed`` *before* the exception
        is raised, and notes that no engine process was started, so a run
        directory never stays ``prepared`` after an attempted execution.
        """
        params = json.loads((prepared.directory / "config.json").read_text(encoding="utf-8"))
        config_keys = ("structure", "supercell", "rmax", "dr", "qmax", "qsig",
                       "radiation", "qmin", "dq", "wavelength", "allow_qsig", "scattering_lengths")

        def status(state: str, **details: object) -> None:
            (prepared.directory / "status.json").write_text(
                json.dumps({"status": state, **details}, indent=2), encoding="utf-8"
            )

        def preflight_fail(exc: Exception) -> None:
            status("failed", reason=f"Preflight check failed: {exc}", engine_started=False)
            raise exc

        with warnings.catch_warnings():  # already warned when the config was created
            warnings.simplefilter("ignore", QmaxLimitWarning)
            warnings.simplefilter("ignore", SetDeltaWarning)
            config = PDFConfig(**{key: params[key] for key in config_keys if key in params})
        staged = prepared.directory / config.staged_name
        if not staged.is_file() or hashlib.sha256(staged.read_bytes()).hexdigest() != params["input_sha256"]:
            preflight_fail(ValueError("Staged structure is missing or differs from the recorded hash"))
        macro = prepared.directory / "run.mac"
        if not macro.is_file():
            preflight_fail(FileNotFoundError("Missing prepared macro"))
        if "macro_sha256" in params and hashlib.sha256(macro.read_bytes()).hexdigest() != params["macro_sha256"]:
            preflight_fail(ValueError("Staged macro differs from the recorded hash; the inspected macro is not the executed macro"))
        for existing in (prepared.raw_output, prepared.output):
            if existing.exists():
                preflight_fail(FileExistsError(f"Refusing to overwrite existing output: {existing}; "
                                              "delete it or use a new run directory to rerun"))

        status("running")
        try:
            completed = current_engine(self).run(prepared.directory, "run.mac", error_pattern=_PDF_ENGINE_ERROR)
        except (FileNotFoundError, subprocess.TimeoutExpired) as exc:
            status("failed", reason=str(exc), engine_started=False)
            raise RuntimeError(f"DISCUS could not complete; see {prepared.directory}") from exc

        engine_error = completed.engine_error
        if completed.returncode != 0 or engine_error or not prepared.raw_output.is_file():
            status("failed", returncode=completed.returncode,
                   reason="Engine error, nonzero exit, or missing output",
                   engine_error=engine_error)
            raise RuntimeError(f"DISCUS run failed; inspect {prepared.directory / 'run.log'}")

        try:
            r, gr = read_pdf(prepared.raw_output)
            n = round(config.rmax / config.dr) + 1
            expected = np.arange(n) * config.dr
            if len(r) != n or not np.allclose(r, expected, rtol=0,
                                              atol=max(1e-6, config.dr * 1e-5)):
                raise ValueError(f"Unexpected r grid: expected {n} samples from 0 to {config.rmax} Å")
            if not np.all(np.isfinite(gr)):
                raise ValueError("Non-finite G(r) values")
        except (ValueError, OSError) as exc:
            status("failed", returncode=completed.returncode, reason=f"Invalid PDF data: {exc}")
            raise RuntimeError(f"Invalid PDF output; inspect {prepared.raw_output}") from exc

        write_pdf_csv(prepared.output, r, gr)
        status("complete", returncode=completed.returncode, output_exists=True,
               samples=len(r), output=PDF_CSV_NAME, raw_output=RAW_PDF_NAME)
        return prepared

    def run(self, config: PDFConfig) -> PDFRun:
        return self.execute(self.prepare(config))
