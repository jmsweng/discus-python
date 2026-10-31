"""Helper for opt-in live DISCUS 6.23.01 tests.

Run with ``DISCUS_LIVE=1`` to enable; tests are skipped otherwise.
On Windows, DISCUS is invoked through WSL (``wsl.exe -d <distro> -- bash -lc …``).
On Linux, ``discus_suite`` is called directly.
"""
from __future__ import annotations

import os
import platform
import re
import subprocess
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path


@dataclass
class LiveResult:
    stdout: str
    returncode: int
    banner: str
    macro: str
    workdir: str


def _wsl_path(win_path: str) -> str:
    """Convert a Windows path (e.g. ``C:\\x\\y``) to a WSL path (``/mnt/c/x/y``)."""
    p = win_path.replace("\\", "/")
    # p is now e.g. "C:/x/y"
    drive = p[0].lower()
    rest = p[2:]  # e.g. "x/y"
    return f"/mnt/{drive}/{rest}"


def run_discus(
    macro_text: str,
    workdir: str,
    *,
    files: dict[str, bytes | str] | None = None,
    timeout: int = 300,
) -> LiveResult:
    """Run a DISCUS macro in ``workdir`` and return the result.

    ``files`` is a mapping of filename → content (bytes or str) to write into
    ``workdir`` before running the macro.  The macro is written as
    ``<workdir>/macro.mac``.

    Raises ``RuntimeError`` if the return code is non-zero, the banner is
    missing, or the banner shows a version other than 6.23.01.
    """
    if " " in workdir:
        raise ValueError(f"workdir must not contain spaces: {workdir!r}")

    workdir_path = Path(workdir)
    workdir_path.mkdir(parents=True, exist_ok=True)

    # Write input files.
    if files:
        for name, content in files.items():
            p = workdir_path / name
            if isinstance(content, bytes):
                p.write_bytes(content)
            else:
                p.write_text(content, encoding="ascii")

    macro_name = "macro.mac"
    macro_path = workdir_path / macro_name
    macro_path.write_text(macro_text, encoding="ascii")

    system = platform.system()
    if system == "Windows":
        distro = os.environ.get("DISCUS_WSL_DISTRO", "Ubuntu-22.04")
        wsl_workdir = _wsl_path(workdir)
        cmd = [
            "wsl.exe", "-d", distro, "--",
            "bash", "-lc",
            f"cd {wsl_workdir} && discus_suite -noautorun -macro {macro_name}",
        ]
    else:
        cmd = [
            "discus_suite", "-noautorun", "-macro", macro_name,
        ]

    try:
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout,
            cwd=str(workdir_path) if system != "Windows" else None,
        )
    except subprocess.TimeoutExpired:
        raise RuntimeError(f"DISCUS timed out after {timeout}s")

    # DISCUS writes to stderr on Windows/WSL; merge stdout and stderr.
    stdout = (proc.stdout or "") + (proc.stderr or "")
    returncode = proc.returncode

    # Find the banner line.
    banner = ""
    for line in stdout.splitlines():
        if "Version 6.23.01" in line:
            banner = line.strip()
            break

    if not banner:
        raise RuntimeError(
            f"DISCUS banner 'Version 6.23.01' not found in output. "
            f"Return code: {returncode}. "
            f"First 500 chars: {stdout[:500]!r}"
        )

    if returncode != 0:
        raise RuntimeError(
            f"DISCUS returned non-zero code {returncode}. "
            f"First 500 chars: {stdout[:500]!r}"
        )

    # Check for DISCUS error patterns.
    for pattern in (r"\*\*\*APPL\*\*\*", r"\*\*\*FORT\*\*\*", r"\*\*\*MAC\*\*\*"):
        if re.search(pattern, stdout):
            raise RuntimeError(
                f"DISCUS error pattern {pattern} found in output. "
                f"First 500 chars: {stdout[:500]!r}"
            )

    return LiveResult(
        stdout=stdout,
        returncode=returncode,
        banner=banner,
        macro=macro_name,
        workdir=workdir,
    )


def evidence_log(
    test_name: str,
    macro_text: str,
    files: dict[str, bytes | str] | None,
    result: LiveResult,
) -> None:
    """Append a log entry to the evidence directory if ``DISCUS_LIVE_EVIDENCE`` is set."""
    evidence_dir = os.environ.get("DISCUS_LIVE_EVIDENCE")
    if not evidence_dir:
        return
    Path(evidence_dir).mkdir(parents=True, exist_ok=True)
    log_path = Path(evidence_dir) / f"{test_name}.log"
    ts = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    lines = [
        f"=== {test_name} @ {ts} ===",
        f"command: discus_suite -noautorun -macro {result.macro}",
        f"workdir: {result.workdir}",
        f"returncode: {result.returncode}",
        f"banner: {result.banner}",
        "--- macro ---",
        macro_text,
        "--- files ---",
    ]
    if files:
        import hashlib
        for name, content in files.items():
            data = content if isinstance(content, bytes) else content.encode("ascii")
            sha = hashlib.sha256(data).hexdigest()
            lines.append(f"{name}: sha256={sha}")
    lines.append("--- stdout ---")
    lines.append(result.stdout)
    with open(log_path, "a", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


def evidence_note(test_name: str, text: str) -> None:
    """Append free text (e.g. the head of a DISCUS-saved file) to the test's evidence log."""
    evidence_dir = os.environ.get("DISCUS_LIVE_EVIDENCE")
    if not evidence_dir:
        return
    Path(evidence_dir).mkdir(parents=True, exist_ok=True)
    with open(Path(evidence_dir) / f"{test_name}.log", "a", encoding="utf-8") as f:
        f.write(text.rstrip("\n") + "\n")
