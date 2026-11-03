"""One way to start DISCUS_SUITE, shared by every runner in the package.

An :class:`Engine` is a DISCUS_SUITE installation: the command that runs a macro (e.g.
``("discus_suite", "-noautorun", "-macro")`` inside WSL/Linux, or
``("wsl.exe", "-d", "Ubuntu-22.04", "--", "discus_suite", "-noautorun", "-macro")`` from
Windows), and a label saying which build it is:

* ``"stock"`` — DISCUS_SUITE 6.23.01 as released (the engine of all validation evidence);
* ``"fork"`` — the owner's fork (``jmsweng/DiffuseCode`` ``master``), which fixes bugs in
  REFINE, RMC with lots, NUFFT and the Debye–Waller cache (``docs/SUITE_CATALOGUE.md``).

:meth:`Engine.run` starts the engine the way every validated run in this repository did: in
the run directory, with ``exit`` on stdin, no ``DISPLAY`` and ``PGPLOT=none`` (otherwise
``pgxwin_server`` can hold the output pipe open after DISCUS exits), and the log written to
``run.log`` as stdout, a ``--- STDERR ---`` line, then stderr.

:func:`fork_reasons` says whether a macro uses something only the fork does right;
:func:`choose_engine` picks the engine accordingly.
"""

from __future__ import annotations

import os
import re
import shlex
import subprocess
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Sequence

DEFAULT_LAUNCHER = ("discus_suite", "-noautorun", "-macro")
ENGINE_LABELS = ("stock", "fork")
# Environment variables that name an engine's launcher (split on whitespace).
LAUNCHER_ENV = {"stock": "DISCUSPDF_ENGINE", "fork": "DISCUSPDF_ENGINE_FORK"}

# Definitive engine-failure markers. PGPLOT/X11 noise ("xdpyinfo: unable to open display",
# "pgxwin_server: cannot connect") is harmless for headless runs and is not matched.
ENGINE_ERROR_PATTERN = re.compile(
    r"(?im)\b(?:unknown command|fatal error|segmentation fault|error\s*:)"
    # DISCUS error lines carry a tag such as ***APPL***, ***FOUR***, ***MAC *** (the start-up
    # banner is a row of bare asterisks and does not match).
    r"|\*\*\*(?:[A-Z]{3,4} ?)\*\*\*"
    # DISCUS 6.23.01 nufft: FINUFFT reports 'setpts error (ier=4)', exits 0 and writes zeros.
    r"|\bFINUFFT\b[^\n]*\berror\b",
)
_ANSI = re.compile(r"\x1b\[[0-9;]*m")


@dataclass(frozen=True)
class EngineResult:
    """What one engine process did."""

    returncode: int
    log: str  # stdout + "\n--- STDERR ---\n" + stderr, as written to run.log
    engine_error: str | None  # first match of the error pattern, or None
    seconds: float
    command: tuple

    @property
    def ok(self) -> bool:
        return self.returncode == 0 and self.engine_error is None


@dataclass(frozen=True)
class Engine:
    """A DISCUS_SUITE installation and how to start it."""

    launcher: tuple = DEFAULT_LAUNCHER
    label: str = "stock"
    version: str = "6.23.01"
    timeout_seconds: float | None = 3600
    env: dict = field(default_factory=dict, hash=False, compare=False)  # extra variables

    def __post_init__(self):
        launcher = tuple(str(x) for x in self.launcher)
        object.__setattr__(self, "launcher", launcher)
        if not launcher:
            raise ValueError("launcher must not be empty")
        if self.label not in ENGINE_LABELS:
            raise ValueError(f"label must be one of {ENGINE_LABELS}, got {self.label!r}")
        if self.timeout_seconds is not None and not self.timeout_seconds > 0:
            raise ValueError("timeout_seconds must be positive (or None for no limit)")

    @classmethod
    def from_env(cls, label: str = "stock", **kwargs) -> "Engine":
        """The engine named by ``DISCUSPDF_ENGINE`` (stock) or ``DISCUSPDF_ENGINE_FORK`` (fork),
        e.g. ``DISCUSPDF_ENGINE_FORK="/opt/fork/bin/discus_suite -noautorun -macro"``. Without the
        variable, the stock engine is ``discus_suite -noautorun -macro``; the fork has no default."""
        value = shlex.split(os.environ.get(LAUNCHER_ENV[label], ""), posix=os.name != "nt")
        if not value:
            if label == "fork":
                raise LookupError(f"no fork engine configured; set {LAUNCHER_ENV['fork']}")
            value = list(DEFAULT_LAUNCHER)
        return cls(tuple(value), label, **kwargs)

    def command(self, macro: str = "run.mac") -> list[str]:
        return [*self.launcher, macro]

    def environment(self) -> dict:
        env = os.environ.copy()
        env.pop("DISPLAY", None)
        env["PGPLOT"] = "none"
        env.update(self.env)
        return env

    def run(self, directory: Path, macro: str = "run.mac", *, log_name: str | None = "run.log",
            error_pattern: re.Pattern | None = ENGINE_ERROR_PATTERN, strip_ansi: bool = False,
            timeout_seconds: float | None = None) -> EngineResult:
        """Run ``macro`` in ``directory`` and write the log.

        ``FileNotFoundError`` (launcher missing) and ``subprocess.TimeoutExpired`` propagate:
        callers record them as "engine not started / did not finish". ``error_pattern`` is
        searched in the log (after removing ANSI colour codes if ``strip_ansi``)."""
        directory = Path(directory)
        timeout = self.timeout_seconds if timeout_seconds is None else timeout_seconds
        if timeout is not None and not timeout > 0:
            raise ValueError("timeout_seconds must be positive")
        command = self.command(macro)
        start = time.monotonic()
        try:
            proc = subprocess.run(command, cwd=directory, input="exit\n", text=True, capture_output=True,
                                  timeout=timeout, env=self.environment(), check=False)
        except subprocess.TimeoutExpired as exc:
            if log_name:  # keep what the engine printed before it was stopped
                def text(x):
                    return x.decode(errors="replace") if isinstance(x, bytes) else (x or "")
                (directory / log_name).write_text(text(exc.stdout) + "\n--- STDERR ---\n" + text(exc.stderr)
                                                  + f"\n--- TIMEOUT after {timeout} s ---\n", encoding="utf-8")
            raise
        seconds = time.monotonic() - start
        log = proc.stdout + "\n--- STDERR ---\n" + proc.stderr
        if log_name:
            (directory / log_name).write_text(log, encoding="utf-8")
        found = None
        if error_pattern is not None:
            m = error_pattern.search(_ANSI.sub("", log) if strip_ansi else log)
            found = m.group(0) if m else None
        return EngineResult(proc.returncode, log, found, seconds, tuple(command))

    def describe(self) -> dict:
        """JSON-serialisable identity, recorded next to every run."""
        return {"label": self.label, "version": self.version, "launcher": list(self.launcher)}


def as_engine(engine: "Engine | None", launcher: Sequence[str] | None, timeout_seconds) -> Engine:
    """The engine a runner uses: ``engine`` if given, else one built from the runner's legacy
    ``launcher`` / ``timeout_seconds`` arguments."""
    if engine is not None:
        if not isinstance(engine, Engine):
            raise TypeError(f"engine must be an Engine, got {type(engine).__name__}")
        return engine
    return Engine(tuple(launcher if launcher is not None else DEFAULT_LAUNCHER), "stock",
                  timeout_seconds=timeout_seconds)


def current_engine(runner) -> Engine:
    """A runner's engine, honouring later assignments to its legacy ``launcher`` and
    ``timeout_seconds`` attributes (older code and tests set them after construction)."""
    from dataclasses import replace

    engine = runner.engine
    launcher = tuple(str(x) for x in runner.launcher)
    if launcher != engine.launcher or runner.timeout_seconds != engine.timeout_seconds:
        engine = replace(engine, launcher=launcher, timeout_seconds=runner.timeout_seconds)
    return engine


# ---------------------------------------------------------------------------------------------
# Which macros need the fork


def _options(params: list[str], declared) -> dict[str, str]:
    """``name:value`` parameters as ``get_optional`` reads them (lib_f90/take_param_mod.f90): an
    abbreviated name selects the first declared option it matches."""
    from .suite_catalogue import str_comp

    out: dict[str, str] = {}
    for p in params:
        name, sep, value = p.partition(":")
        if not sep or not name:
            continue
        for opt in declared:
            if str_comp(name, opt, min(len(name), len(opt)), len(opt)):
                out.setdefault(opt, value)
                break
    return out


def _number(text: str) -> float | None:
    try:
        return float(text)
    except ValueError:
        return None


def fork_reasons(macro: str, catalogue=None) -> list[str]:
    """Reasons a macro should run on the fork rather than stock 6.23.01 (empty = stock is fine).

    The rules follow ``docs/SUITE_ENGINE.md`` ("Stock or fork"). They flag features whose
    results the fork corrects, not its speed-ups, and read options and keywords with DISCUS's
    own abbreviation rules. A Fourier ``run`` inside a ``do``/``if`` block or after an
    ``@macro`` call may run more than once, so it counts as a repeated run.
    """
    from .suite_catalogue import load_catalogue

    cat = catalogue or load_catalogue()
    reasons: list[str] = []
    temp_use = False
    runs_with_temp_use = 0
    depth = 0
    macro_called = False
    for step in cat.walk(macro):
        if step.kind == "control":
            depth += 1 if step.command in ("do", "if") else -1 if step.command in ("enddo", "endif") else 0
            continue
        if step.kind == "macro":
            macro_called = True
            continue
        if step.kind != "command":
            continue
        name = step.command
        cmd = cat.resolve(step.menu, step.text.split()[0])
        args = step.text.split(None, 1)[1] if len(step.text.split(None, 1)) > 1 else ""
        params = [p.strip().replace(" ", "") for p in args.split(",")] if args else []
        if step.menu == "refine" and name == "run":
            reasons.append(f"line {step.line}: REFINE 'run' (fork fixes MPI derivatives, the dP/sigma "
                           "criterion, seed restore with Gaussian noise and MPI result transfer)")
        elif step.menu == "refine" and name == "newparam":
            points = _options(params, cmd.optional).get("points")
            if points is not None and _number(points) == 2:
                reasons.append(f"line {step.line}: REFINE 'points:2' exists only on the fork")
        elif step.menu == "discus/rmc" and name == "set" and params and params[0].upper()[:3] == "LOT":
            # rmc_sup_mod.f90: 'set lots,O...' (one parameter starting with O) switches lots off
            off = len(params) == 2 and params[1].upper()[:1] == "O"
            if not off:
                reasons.append(f"line {step.line}: RMC with lots (fork fixes lot moves and updates)")
        elif step.menu == "discus/rmc" and name == "set" and params and params[0].upper()[:2] == "AV":
            # rmc_sup_mod.f90:62-67 ('AV', one value in percent); rmc_calc_scat subtracts the average
            # with a row loop bounded by nx instead of ny (:2300), fixed on the fork
            value = _number(params[1]) if len(params) == 2 else None
            if value is None or value != 0:
                reasons.append(f"line {step.line}: RMC 'set aver' (fork fixes the average subtraction "
                               "on planes with nx != ny)")
        elif step.menu == "discus/fourier" and name == "set":
            if _options(params, cmd.optional).get("technique") == "nufft":
                reasons.append(f"line {step.line}: NUFFT (fork fixes the scaling on non-integer grids)")
        elif step.menu == "discus/fourier" and name == "temp" and params:
            # fourier.f90:876-878 test the first two characters
            temp_use = params[0][:2] == "us" if params[0][:2] in ("us", "ig") else temp_use
        elif step.menu == "discus/fourier" and name == "reset":
            temp_use = False  # fourier_reset clears the Debye-Waller flag
        elif step.menu == "discus/fourier" and name == "run" and temp_use:
            runs_with_temp_use += 2 if (depth or macro_called) else 1
            if runs_with_temp_use >= 2 and not any("Debye-Waller" in r for r in reasons):
                reasons.append(f"line {step.line}: more than one Fourier run with 'temp use' in one session "
                               "(fork fixes the Debye-Waller cache key)")
    return reasons


class EngineChoiceError(RuntimeError):
    """A macro needs the fork but no fork engine was given."""


def choose_engine(macro: str, stock: Engine, fork: Engine | None = None, catalogue=None) -> Engine:
    """``fork`` if :func:`fork_reasons` finds a reason, otherwise ``stock``."""
    if fork is not None and fork.label != "fork":
        raise EngineChoiceError(f"engine given as fork is labelled {fork.label!r}")
    if stock.label != "stock":
        raise EngineChoiceError(f"engine given as stock is labelled {stock.label!r}")
    reasons = fork_reasons(macro, catalogue)
    if not reasons:
        return stock
    if fork is None:
        raise EngineChoiceError("this macro needs the fork engine: " + "; ".join(reasons))
    return fork
