"""Read-only check of one probe evidence directory before it is committed.

    python benchmarks/discus_6_23_01/b_probes/check_evidence.py benchmarks/discus_6_23_01/b_probes/evidence/<stamp>
    python benchmarks/discus_6_23_01/b_probes/check_evidence.py <smoke evidence dir> --partial

Checks: MANIFEST.sha256 matches every file (and lists every file except itself), runs_status.json,
summary.json and SUMMARY.md exist, every run in runs.json has a directory with run.mac and run.log,
staged macros equal the committed ones byte for byte, no tracked input was copied in, no file is
larger than 5 MB and the directory is under 40 MB. Prints PASS or the list of problems (exit 1).
Run it before `git add`: after a fresh Windows checkout (autocrlf) text files change on disk.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    partial = "--partial" in argv
    if partial:
        argv.remove("--partial")
    if len(argv) != 1:
        print(__doc__)
        return 2
    ev = Path(argv[0])
    problems: list[str] = []
    files = {p.relative_to(ev).as_posix(): p for p in ev.rglob("*") if p.is_file()}
    manifest = {}
    if "MANIFEST.sha256" not in files:
        problems.append("MANIFEST.sha256 missing")
    else:
        for line in files["MANIFEST.sha256"].read_text(encoding="utf-8").splitlines():
            digest, name = line.split("  ", 1)
            manifest[name] = digest
    late = {"MANIFEST.sha256", "summary.json", "SUMMARY.md"}       # summary is written after the manifest
    for name, path in files.items():
        if name in late:
            continue
        if name not in manifest:
            problems.append(f"not in manifest: {name}")
        elif hashlib.sha256(path.read_bytes()).hexdigest() != manifest[name]:
            problems.append(f"hash mismatch: {name}")
    for name in manifest:
        if name not in files:
            problems.append(f"in manifest but missing: {name}")
    for need in ("runs_status.json", "summary.json", "SUMMARY.md", "help_sections_b.txt"):
        if need not in files:
            problems.append(f"missing {need}")
    spec = json.loads((HERE / "runs.json").read_text(encoding="utf-8"))
    status_ids = set()
    if "runs_status.json" in files:
        status_ids = {r["id"] for r in json.loads(files["runs_status.json"].read_text(encoding="utf-8"))["runs"]}
    inputs = {Path(i).name for r in spec["runs"] for i in r["inputs"]}
    for r in spec["runs"]:
        if r["id"] not in status_ids:
            continue                                  # partial run (--only): checked runs only
        for need in ("run.mac", "run.log"):
            if f"{r['id']}/{need}" not in files:
                problems.append(f"missing {r['id']}/{need}")
        staged = files.get(f"{r['id']}/run.mac")
        committed = (HERE / r["macro"]).read_bytes().replace(b"\r\n", b"\n")
        if staged is not None and staged.read_bytes() != committed:
            problems.append(f"{r['id']}/run.mac differs from {r['macro']}")
    missing_runs = sorted({r["id"] for r in spec["runs"]} - status_ids)
    if missing_runs and not partial:
        problems.append(f"runs not executed: {missing_runs}")
    for name, path in files.items():
        if Path(name).name in inputs:
            problems.append(f"tracked input copied into evidence: {name}")
        if path.stat().st_size > 5_000_000:
            problems.append(f"file larger than 5 MB: {name}")
    total = sum(p.stat().st_size for p in files.values())
    if total > 40_000_000:
        problems.append(f"evidence is {total/1e6:.1f} MB (> 40 MB)")
    print(f"{len(files)} files, {total/1e6:.2f} MB, {len(status_ids)} runs")
    if problems:
        print("FAIL")
        for p in problems:
            print(" -", p)
        return 1
    print("PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
