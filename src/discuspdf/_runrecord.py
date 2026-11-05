"""Run-directory bookkeeping shared by the runners (private).

PDFRunner, Diffuse3DRunner, MCRunner and SuiteRunner each keep a run directory with the staged
inputs, their SHA-256 hashes and a ``status.json`` (``prepared`` -> ``running`` -> ``complete`` /
``failed``). The file names and the JSON layout of each runner are unchanged (archived run
directories are read back by ``load``); only the code that writes and checks them is shared here.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Mapping


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_status(directory: Path, state: str, details: Mapping | None = None, *,
                 indent: int | None = 2, newline: bool = True) -> None:
    """Write ``{"status": state, **details}`` to ``directory/status.json``. ``indent`` and ``newline``
    keep each runner's existing layout byte for byte."""
    text = json.dumps({"status": state, **(details or {})}, indent=indent)
    (Path(directory) / "status.json").write_text(text + ("\n" if newline else ""), encoding="utf-8")


def changed_inputs(directory: Path, hashes: Mapping[str, str]) -> list[tuple[str, str]]:
    """``(name, "missing" | "modified")`` for every recorded input that is no longer as staged."""
    out = []
    for name, digest in hashes.items():
        path = Path(directory) / name
        if not path.is_file():
            out.append((name, "missing"))
        elif sha256_file(path) != digest:
            out.append((name, "modified"))
    return out
