"""Run any DISCUS_SUITE macro in an isolated, inspectable run directory.

:class:`SuiteRunner` is the general runner behind the suite interface: it accepts a
:class:`~discuspdf.macro.Macro` (or macro text), the input files it reads and the output
files it must produce, and follows the same rules as the package's specialised runners
(``PDFRunner``, ``Diffuse3DRunner``, ``MCRunner``):

* ``prepare`` writes ``run.mac``, copies the inputs and records their SHA-256 in
  ``run.json`` without starting DISCUS, so the macro can be inspected first;
* ``execute`` runs exactly that macro (hashes re-checked), writes ``run.log`` and
  ``status.json`` (``prepared`` → ``running`` → ``complete`` / ``failed``), and fails on a
  nonzero exit, an engine error in the log, or a missing or empty output;
* the engine is chosen per macro: the fork when :func:`discuspdf.engine.fork_reasons` finds
  something only the fork does right, otherwise stock DISCUS 6.23.01.

``complete`` means the engine ran cleanly and wrote the files. It says nothing about whether
the numbers are physically right.
"""

from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
import uuid
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Mapping

from .engine import Engine, EngineChoiceError, fork_reasons
from .macro import Macro
from .suite_catalogue import load_catalogue

RUN_FORMAT = "discuspdf.suite_run"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _plain_name(name: str, what: str) -> str:
    if (not isinstance(name, str) or not name or name in (".", "..") or Path(name).name != name
            or "\\" in name or any(c.isspace() for c in name)):
        raise ValueError(f"{what} must be a plain file name without spaces, got {name!r}")
    return name


@dataclass
class SuiteRun:
    directory: Path

    @property
    def macro(self) -> Path:
        return self.directory / "run.mac"

    @property
    def log(self) -> Path:
        return self.directory / "run.log"

    @property
    def record(self) -> dict:
        return json.loads((self.directory / "run.json").read_text(encoding="utf-8"))

    def status(self) -> dict:
        return json.loads((self.directory / "status.json").read_text(encoding="utf-8"))

    @property
    def outputs(self) -> list[Path]:
        return [self.directory / name for name in self.record["outputs"]]

    def output(self, name: str) -> Path:
        if name not in self.record["outputs"]:
            raise KeyError(f"{name!r} is not a declared output of this run")
        return self.directory / name


class SuiteRunner:
    """Prepare and run macros with the stock engine, or the fork when the macro needs it."""

    def __init__(self, runs_directory: Path | str = "runs_suite", engine: Engine | None = None,
                 fork: Engine | None = None):
        self.runs_directory = Path(runs_directory)
        self.engine = engine or Engine()
        self.fork = fork
        if self.engine.label != "stock":
            raise ValueError("engine must be the stock engine; pass the fork as fork=")
        if fork is not None and fork.label != "fork":
            raise ValueError("fork must be an Engine labelled 'fork'")

    def prepare(self, macro: Macro | str, inputs: Mapping[str, Path | str] | Iterable[Path | str] = (),
                outputs: Iterable[str] = (), *, engine: str | None = None) -> SuiteRun:
        """Stage a run without starting DISCUS.

        ``inputs`` maps the names the macro uses to files (a list keeps each file's own name).
        ``outputs`` are file names the macro must write. ``engine`` forces ``"stock"`` or
        ``"fork"``; by default it is chosen from the macro (stock unless the fork is needed;
        a macro that needs the fork fails here if the runner has none)."""
        text = macro.render() if isinstance(macro, Macro) else str(macro)
        if not text.endswith("\n"):
            text += "\n"
        unknown = load_catalogue().check(text)
        if unknown:
            raise ValueError("macro has lines DISCUS_SUITE 6.23.01 does not know: "
                             + "; ".join(f"line {s.line} in {s.menu}: {s.text!r}" for s in unknown))
        reasons = fork_reasons(text)
        chosen = self._choose(engine, reasons)
        if not isinstance(inputs, Mapping):
            inputs = {Path(p).name: p for p in inputs}
        names = {_plain_name(n, "input name"): Path(p).expanduser() for n, p in inputs.items()}
        outs = [_plain_name(n, "output name") for n in outputs]
        reserved = {"run.mac", "run.log", "run.json", "status.json"}
        clash = (set(names) | set(outs)) & reserved or set(names) & set(outs)
        if clash or len(set(outs)) != len(outs):
            raise ValueError(f"file names clash or are reserved: {sorted(clash) or outs}")
        for name, path in names.items():
            if not path.is_file():
                raise FileNotFoundError(f"input {name!r}: {path}")
        self.runs_directory.mkdir(parents=True, exist_ok=True)
        directory = self.runs_directory / (datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
                                           + "_" + uuid.uuid4().hex[:8])
        directory.mkdir()
        (directory / "run.mac").write_text(text, encoding="utf-8", newline="\n")
        hashes = {"run.mac": _sha256(directory / "run.mac")}
        for name, path in names.items():
            shutil.copyfile(path, directory / name)
            hashes[name] = _sha256(directory / name)
        record = {"format": RUN_FORMAT, "format_version": 1, "engine": chosen.describe(),
                  "fork_reasons": reasons,
                  "stock_despite_fork_reasons": bool(reasons) and chosen.label == "stock", "inputs": sorted(names), "outputs": outs,
                  "sha256": hashes}
        (directory / "run.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
        run = SuiteRun(directory)
        self._status(run, "prepared")
        return run

    def execute(self, run: SuiteRun) -> SuiteRun:
        record = run.record
        if record.get("format") != RUN_FORMAT:
            raise ValueError(f"{run.directory} is not a suite run")
        engine = self.engine if record["engine"]["label"] == "stock" else self.fork
        if engine is None or list(engine.launcher) != record["engine"]["launcher"]:
            self._fail(run, "the engine recorded at prepare() is not this runner's", engine_started=False)
        for name, digest in record["sha256"].items():
            path = run.directory / name
            if not path.is_file() or _sha256(path) != digest:
                self._fail(run, f"staged {name} is missing or was modified after prepare()",
                           engine_started=False)
        for path in run.outputs:
            if path.exists():
                self._fail(run, f"output {path.name} exists already; prepare a new run", engine_started=False)
        self._status(run, "running")
        try:
            result = engine.run(run.directory, "run.mac", strip_ansi=True)
        except FileNotFoundError as exc:
            self._fail(run, str(exc), engine_started=False)
        except subprocess.TimeoutExpired as exc:
            self._fail(run, f"timed out after {exc.timeout} s; the engine was stopped", engine_started=True)
        missing = [p.name for p in run.outputs if not p.is_file() or p.stat().st_size == 0]
        if result.returncode != 0 or result.engine_error or missing:
            reason = ("nonzero exit" if result.returncode != 0 else
                      f"engine error in log: {result.engine_error!r}" if result.engine_error else
                      f"missing or empty output: {', '.join(missing)}")
            self._fail(run, reason, returncode=result.returncode, engine_error=result.engine_error)
        warnings = (["run on stock although the fork corrects: " + "; ".join(record["fork_reasons"])]
                    if record.get("stock_despite_fork_reasons") else [])
        self._status(run, "complete", returncode=result.returncode, seconds=round(result.seconds, 3),
                     outputs={p.name: _sha256(p) for p in run.outputs}, warnings=warnings)
        return run

    def run(self, macro: Macro | str, inputs=(), outputs=(), *, engine: str | None = None) -> SuiteRun:
        return self.execute(self.prepare(macro, inputs, outputs, engine=engine))

    # -- helpers ---------------------------------------------------------------------------

    def _choose(self, forced: str | None, reasons: list[str]) -> Engine:
        if forced not in (None, "stock", "fork"):
            raise ValueError("engine must be 'stock', 'fork' or None")
        if forced == "stock":
            return self.engine
        if forced == "fork" or reasons:
            if self.fork is None:
                raise EngineChoiceError("this run needs the fork engine"
                                        + (": " + "; ".join(reasons) if reasons else ""))
            return self.fork
        return self.engine

    @staticmethod
    def _status(run: SuiteRun, state: str, **details) -> None:
        (run.directory / "status.json").write_text(json.dumps({"status": state, **details}, indent=2) + "\n",
                                                   encoding="utf-8")

    def _fail(self, run: SuiteRun, reason: str, **details):
        self._status(run, "failed", reason=reason, **details)
        raise RuntimeError(f"DISCUS_SUITE run failed: {reason} (see {run.directory})")
