"""Shared set-up for the opt-in live DISCUS tests (``DISCUS_LIVE=1``): one launcher, one Engine, one runner.

On Windows DISCUS runs through WSL (``DISCUS_WSL_DISTRO``, default Ubuntu-22.04; ``wsl.exe`` translates the
Windows working directory to /mnt/<drive>/...); elsewhere ``discus_suite`` must be on PATH.
"""
from __future__ import annotations

import os
import platform
from pathlib import Path

import pytest

from discuspdf.engine import Engine
from discuspdf.suite_runner import SuiteRunner

DISCUS_ARGS = ("discus_suite", "-noautorun", "-macro")


def live_launcher() -> tuple[str, ...]:
    if platform.system() == "Windows":
        return ("wsl.exe", "-d", os.environ.get("DISCUS_WSL_DISTRO", "Ubuntu-22.04"), "--", *DISCUS_ARGS)
    return DISCUS_ARGS


def live_engine(timeout_seconds: float = 600, name: str | None = None) -> Engine:
    if platform.system() == "Windows":
        return Engine(live_launcher(), *(() if name is None else (name,)), timeout_seconds=timeout_seconds)
    return Engine(timeout_seconds=timeout_seconds)


def live_work(default_name: str) -> Path:
    """``DISCUS_LIVE_WORK`` or C:\\tmp\\<name> (Windows) / /tmp/<name>."""
    return Path(os.environ.get("DISCUS_LIVE_WORK",
                               rf"C:\tmp\{default_name}" if platform.system() == "Windows" else f"/tmp/{default_name}"))


def live_runner(work: Path, timeout_seconds: float = 600, **kw) -> SuiteRunner:
    """A SuiteRunner on the stock engine; skips when the work directory has spaces (DISCUS cannot use it)."""
    if " " in str(work):
        pytest.skip("DISCUS_LIVE_WORK must not contain spaces")
    return SuiteRunner(work, live_engine(timeout_seconds), **kw)
