"""Minimal, isolated DISCUS 1D PDF simulation runner.

This adapter uses documented legacy DISCUS PDF commands. Validate once against
YOUR installed version of discus_suite before running a large sweep.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Sequence
import hashlib
import json
import math
import os
import shutil
import subprocess
import uuid
import warnings

import numpy as np


# DISCUS_SUITE 6.23.01 cannot use Q above about 12.3/wavelength (Å⁻¹); larger qmax values give the same
# G(r). Measured: G(r) stops changing between Qmax 8.0 and 8.1 at λ = 1.5406 Å and between 17.0 and 17.5
# at λ = 0.7107 Å. Above Q ≈ 8 Å⁻¹ the DISCUS F(Q) also depends on the wavelength by up to several percent
# (debye mode, X-ray and neutron). See docs/VERIFIED_DISCUS_MACRO.md, "Qmax limit and wavelength".
QMAX_WAVELENGTH_WARN = 12.0      # warn when qmax * wavelength exceeds this (margin below the measured limit)
QMAX_WAVELENGTH_CEILING = 12.3   # approximate measured limit, used only in the warning text


class QmaxLimitWarning(UserWarning):
    """qmax is above what DISCUS can use at this wavelength; the PDF is computed unchanged."""


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
    """1D powder PDF calculation from a DISCUS .cll unit cell.

    All Q values are inverse angstroms. qsig is the DISCUS powder
    "set delta" convolution width; it is not assumed to be a universal
    instrument-resolution parameter. qmax must be explicitly positive:
    omitting "set qmax" is not a verified way to disable termination.
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
        if Path(self.structure).suffix.lower() != ".cll":
            raise ValueError("MVP supports DISCUS .cll unit-cell files only")

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
    """
    nx, ny, nz = config.supercell
    # Defaults are explicit in PDFConfig and can be overridden per run.
    wavelength = default_wavelength(config)
    lines = [
        "# Generated by discuspdf: 1D atomic PDF calculation (DISCUS v6.x)",
        "discus",
        "read",
        f"cell model.cll,{nx},{ny},{nz}",
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


@dataclass(frozen=True)
class PDFRun:
    directory: Path

    @property
    def output(self) -> Path:
        return self.directory / "result.pdf"

    def load(self) -> tuple[np.ndarray, np.ndarray]:
        return read_pdf(self.output)

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


class PDFRunner:
    def __init__(
        self,
        runs_directory: Path | str = "runs",
        launcher: Sequence[str] = ("discus_suite", "-noautorun", "-macro"),
        timeout_seconds: int = 3600,
    ) -> None:
        """launcher is the executable and args preceding 'run.mac'."""
        self.runs_directory = Path(runs_directory)
        self.launcher = tuple(str(x) for x in launcher)
        if not self.launcher:
            raise ValueError("launcher must not be empty")
        if timeout_seconds <= 0:
            raise ValueError("timeout_seconds must be positive")
        self.timeout_seconds = timeout_seconds

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
        shutil.copy2(source, run_dir / "model.cll")
        (run_dir / "run.mac").write_text(render_pdf_macro(config), encoding="utf-8")
        params = asdict(config)
        params["structure"] = str(source)
        params["input_sha256"] = hashlib.sha256(source.read_bytes()).hexdigest()
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
                       "radiation", "qmin", "dq", "wavelength")

        def status(state: str, **details: object) -> None:
            (prepared.directory / "status.json").write_text(
                json.dumps({"status": state, **details}, indent=2), encoding="utf-8"
            )

        def preflight_fail(exc: Exception) -> None:
            status("failed", reason=f"Preflight check failed: {exc}", engine_started=False)
            raise exc

        with warnings.catch_warnings():  # already warned when the config was created
            warnings.simplefilter("ignore", QmaxLimitWarning)
            config = PDFConfig(**{key: params[key] for key in config_keys if key in params})
        staged = prepared.directory / "model.cll"
        if not staged.is_file() or hashlib.sha256(staged.read_bytes()).hexdigest() != params["input_sha256"]:
            preflight_fail(ValueError("Staged structure is missing or differs from the recorded hash"))
        macro = prepared.directory / "run.mac"
        if not macro.is_file():
            preflight_fail(FileNotFoundError("Missing prepared macro"))
        if "macro_sha256" in params and hashlib.sha256(macro.read_bytes()).hexdigest() != params["macro_sha256"]:
            preflight_fail(ValueError("Staged macro differs from the recorded hash; the inspected macro is not the executed macro"))
        if prepared.output.exists():
            preflight_fail(FileExistsError(f"Refusing to overwrite existing output: {prepared.output}; "
                                          "delete it or use a new run directory to rerun"))

        command = [*self.launcher, "run.mac"]
        env = os.environ.copy()
        env.pop("DISPLAY", None)
        env["PGPLOT"] = "none"

        status("running")
        try:
            completed = subprocess.run(
                command, cwd=prepared.directory, input="exit\n", text=True,
                capture_output=True, timeout=self.timeout_seconds, check=False, env=env,
            )
        except (FileNotFoundError, subprocess.TimeoutExpired) as exc:
            status("failed", reason=str(exc), engine_started=False)
            raise RuntimeError(f"DISCUS could not complete; see {prepared.directory}") from exc

        log = completed.stdout + "\n--- STDERR ---\n" + completed.stderr
        (prepared.directory / "run.log").write_text(log, encoding="utf-8")
        import re
        engine_error = re.search(r"(?im)\b(?:unknown command|fatal error|segmentation fault|error\s*:)", log)
        if completed.returncode != 0 or engine_error or not prepared.output.is_file():
            status("failed", returncode=completed.returncode,
                   reason="Engine error, nonzero exit, or missing output",
                   engine_error=engine_error.group(0) if engine_error else None)
            raise RuntimeError(f"DISCUS run failed; inspect {prepared.directory / 'run.log'}")

        try:
            r, gr = prepared.load()
            n = round(config.rmax / config.dr) + 1
            expected = np.arange(n) * config.dr
            if len(r) != n or not np.allclose(r, expected, rtol=0,
                                              atol=max(1e-6, config.dr * 1e-5)):
                raise ValueError(f"Unexpected r grid: expected {n} samples from 0 to {config.rmax} Å")
            if not np.all(np.isfinite(gr)):
                raise ValueError("Non-finite G(r) values")
        except (ValueError, OSError) as exc:
            status("failed", returncode=completed.returncode, reason=f"Invalid PDF data: {exc}")
            raise RuntimeError(f"Invalid PDF output; inspect {prepared.output}") from exc

        status("complete", returncode=completed.returncode, output_exists=True,
               samples=len(r))
        return prepared

    def run(self, config: PDFConfig) -> PDFRun:
        return self.execute(self.prepare(config))
