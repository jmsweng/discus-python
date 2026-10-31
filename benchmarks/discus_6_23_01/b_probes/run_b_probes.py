"""Run the Package B DISCUS probes and collect evidence (logs, HDF5, hashes).

Owner's machine (PowerShell, repository root, DISCUS in WSL Ubuntu-22.04):
    .\\.venv\\Scripts\\python.exe benchmarks\\discus_6_23_01\\b_probes\\run_b_probes.py --work C:\\tmp\\bprobe

Each run gets its own directory under <work>/<stamp>/<id>/ holding run.mac, the input file, run.log
and the HDF5 outputs. The whole <stamp> directory is then copied to
benchmarks/discus_6_23_01/b_probes/evidence/<stamp>/ with MANIFEST.sha256 and runs_status.json.
The script never edits macros and never deletes anything. Exit code 0 means every run behaved as
runs.json expects (ok / error / either); 1 otherwise. Evidence is written in both cases.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import platform
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path, PureWindowsPath

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
# DISCUS errors look like "***APPL*** ..."; the banner is a row of bare asterisks (same rule as gate2 analyze.py).
ERROR_RE = re.compile(r"(?m)^.*(?:\*\*\*[A-Z]+\*\*\*|\bERROR\b|\bError\b|[Ss]egmentation fault|[Uu]nknown command).*$")
HELP_KEYS = ("lots", "set", "value", "seed", "technique")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def stage_lf(src: Path, dest: Path) -> None:
    """Copy a text file with LF line endings. A Windows checkout with core.autocrlf=true has CRLF,
    and DISCUS would read the carriage return as part of each command or atom line. The staged
    copy then has the same bytes (and SHA-256) as the committed git blob."""
    data = src.read_bytes()
    if b"\0" in data[:4096]:
        raise ValueError(f"{src} is not a text file")
    dest.write_bytes(data.replace(b"\r\n", b"\n"))


def to_wsl(path: Path) -> str:
    p = PureWindowsPath(str(path))
    if not p.drive:
        raise ValueError(f"--work must be an absolute Windows path with a drive letter, got {path}")
    rest = "/".join(p.parts[1:])
    return f"/mnt/{p.drive[0].lower()}/{rest}"


def shell(cmd: str, cwd: Path, launcher: str, distro: str, timeout: int) -> subprocess.CompletedProcess:
    if launcher == "wsl":
        full = f"cd '{to_wsl(cwd)}' && {cmd}"
        return subprocess.run(["wsl.exe", "-d", distro, "--", "bash", "-lc", full], capture_output=True,
                              text=True, timeout=timeout, encoding="utf-8", errors="replace")
    return subprocess.run(["bash", "-lc", cmd], cwd=cwd, capture_output=True, text=True, timeout=timeout,
                          encoding="utf-8", errors="replace")


def engine_info(work: Path, launcher: str, distro: str) -> dict:
    info = {}
    r = shell("p=$(command -v discus_suite); echo \"$p\"; sha256sum \"$p\"", work, launcher, distro, 120)
    lines = r.stdout.strip().splitlines()
    info["executable"] = lines[0] if lines else None
    info["executable_sha256"] = lines[1].split()[0] if len(lines) > 1 else None
    return info


def help_sections(work: Path, launcher: str, distro: str) -> str:
    r = shell("cat /usr/local/share/suite.hlp 2>/dev/null || true", work, launcher, distro, 120)
    out, keep = [], False
    for line in r.stdout.splitlines():
        m = re.match(r"^(\d)\s+(\S+)", line)
        if m:
            keep = m.group(2).lower() in HELP_KEYS
        if keep:
            out.append(line)
    return "\n".join(out) + "\n"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--work", required=True, type=Path, help="scratch directory without spaces, e.g. C:\\tmp\\bprobe")
    ap.add_argument("--launcher", choices=("wsl", "local"), default="wsl" if sys.platform == "win32" else "local")
    ap.add_argument("--distro", default="Ubuntu-22.04")
    ap.add_argument("--only", default="", help="comma-separated run ids (default: all)")
    ap.add_argument("--timeout", type=int, default=3600, help="seconds per run")
    ap.add_argument("--evidence", type=Path, default=HERE / "evidence")
    args = ap.parse_args(argv)
    if " " in str(args.work):
        ap.error("--work must not contain spaces (DISCUS macro paths)")

    spec = json.loads((HERE / "runs.json").read_text(encoding="utf-8"))
    runs = spec["runs"]
    if args.only:
        wanted = [s.strip() for s in args.only.split(",") if s.strip()]
        unknown = sorted(set(wanted) - {r["id"] for r in runs})
        if unknown:
            ap.error(f"unknown run ids: {unknown}")
        runs = [r for r in runs if r["id"] in wanted]

    stamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "_" + platform.node()
    base = args.work / stamp
    base.mkdir(parents=True)
    status = {"stamp": stamp, "launcher": args.launcher, "distro": args.distro,
              "python": sys.version.split()[0], "engine": engine_info(base, args.launcher, args.distro),
              "runs": []}
    all_as_expected = True
    for r in runs:
        d = base / r["id"]
        d.mkdir()
        stage_lf(HERE / r["macro"], d / "run.mac")
        input_hashes = {}
        for inp in r["inputs"]:
            stage_lf(REPO / inp, d / Path(inp).name)
            input_hashes[inp] = sha256(d / Path(inp).name)
        t0 = time.time()
        try:
            proc = shell("unset DISPLAY; export PGPLOT_DEV=/null; printf 'exit\\n' | discus_suite -noautorun -macro run.mac",
                         d, args.launcher, args.distro, args.timeout)
            log, rc = proc.stdout + "\n--- STDERR ---\n" + proc.stderr, proc.returncode
        except subprocess.TimeoutExpired as exc:
            log, rc = f"TIMEOUT after {args.timeout} s\n{exc}", None
        seconds = round(time.time() - t0, 2)
        (d / "run.log").write_text(log, encoding="utf-8", newline="\n")
        errors = [m.group(0).strip() for m in ERROR_RE.finditer(log)][:20]
        present = [o for o in r["outputs"] if (d / o).is_file() and (d / o).stat().st_size > 0]
        engine_ok = rc == 0 and not errors and len(present) == len(r["outputs"])
        if r["expect"] == "ok":
            as_expected = engine_ok
        elif r["expect"] == "error":
            as_expected = bool(errors) or rc not in (0,)
        else:
            as_expected = True
        all_as_expected &= as_expected
        status["runs"].append({"id": r["id"], "expect": r["expect"], "returncode": rc, "seconds": seconds,
                               "error_lines": errors, "outputs_present": present,
                               "outputs_missing": sorted(set(r["outputs"]) - set(present)),
                               "inputs_sha256": input_hashes,
                               "as_expected": as_expected})
        print(f"{r['id']:<28} rc={rc} {seconds:>8.1f}s outputs {len(present)}/{len(r['outputs'])} "
              f"errors {len(errors)} {'OK' if as_expected else 'UNEXPECTED'}", flush=True)

    (base / "help_sections_b.txt").write_text(help_sections(base, args.launcher, args.distro),
                                              encoding="utf-8", newline="\n")
    status["all_as_expected"] = all_as_expected
    (base / "runs_status.json").write_text(json.dumps(status, indent=1) + "\n", encoding="utf-8", newline="\n")

    # Inputs are tracked in the repository; the evidence keeps their hashes, not second copies.
    input_names = {Path(i).name for r in runs for i in r["inputs"]}
    dest = args.evidence / stamp
    shutil.copytree(base, dest, ignore=lambda folder, names: [n for n in names if n in input_names])
    files = sorted(p for p in dest.rglob("*") if p.is_file())
    manifest = "".join(f"{sha256(p)}  {p.relative_to(dest).as_posix()}\n" for p in files)
    (dest / "MANIFEST.sha256").write_text(manifest, encoding="utf-8", newline="\n")
    print(f"evidence: {dest}")
    print("ALL AS EXPECTED" if all_as_expected else "SOME RUNS UNEXPECTED (see runs_status.json)")
    return 0 if all_as_expected else 1


if __name__ == "__main__":
    sys.exit(main())
