"""Re-create the DISCUS engine HDF5 volumes of the Gate 1 / Gate 2 evidence directories.

The volumes (~200 MB) are not kept in the repository. Everything needed to recreate them is: each evidence directory
keeps the exact .stru and .mac files that were run, the logs, ENGINE_ENV.txt (engine binary SHA-256, registered commit)
and MANIFEST.sha256. This script reruns those macros with DISCUS_SUITE 6.23.01 and checks the result against the
content hashes in EVIDENCE_CONTENT_HASHES.json, which were computed from the original files before they were removed.

Why a *content* hash: DISCUS writes the run date (audit_creation_date) into every HDF5 file, so a regenerated file can
never have the original file's SHA-256 (listed in MANIFEST.sha256 for the record). The content hash covers the
intensity array (float64, little-endian, C order) and the grid/cell/Bragg metadata, and nothing date- or host-dependent.

Usage
  python regenerate_evidence.py list
  python regenerate_evidence.py regenerate <evidence_dir> <new_empty_out_dir> [--discus /usr/local/bin/discus_suite]
                                [--analyze]
  python regenerate_evidence.py hash <dir_with_h5>          # print content hashes (used to build the JSON)
  python regenerate_evidence.py fetch-history               # fetch the PR #8 commits (needed after a squash merge)

The registered commits named in ENGINE_ENV.txt, and the commit 8cd59ac whose original Gate 2 v4 volumes Gate 3 reads, are
not in main's history: PR #8 was squash-merged. GitHub keeps them under refs/pull/8/head even after the branch is deleted,
and `fetch-history` (run automatically by --analyze) fetches that ref.

`regenerate` refuses to run unless every input it copies matches MANIFEST.sha256 of the evidence directory, never
writes into the evidence directory, and warns if the engine binary differs from ENGINE_ENV.txt.

Possible outcomes per volume
  MATCH          content identical to the original (expected on the same binary; the archived Gate 1 rerun was
                 bit-identical in its data).
  DIFFERS        content differs. Rounding order can change with another binary, compiler or OpenMP thread count. Use
                 --analyze: it reruns the evaluator *of the registered commit* (from a temporary git worktree) on the
                 new volumes. For Gate 2 v4 that evaluator compares against an independent oracle that reproduces the
                 engine to ~1e-16, so its verdict and numbers show whether the regenerated data are scientifically the same.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
HASHES = HERE / "EVIDENCE_CONTENT_HASHES.json"
META_KEYS = ("data_corner", "data_increment_vector", "data_dimension", "data_axes", "unit_cell_lengths",
             "unit_cell_angles", "data_type_with_bragg", "data_type_content", "data_radiation", "data_rad_symbol")
SCHEMA = "discus-h5-content-v1"


def content_hash(path: Path) -> str:
    """SHA-256 over data_values (float64 LE, C order) + selected metadata; excludes audit_* fields."""
    import h5py
    h = hashlib.sha256(SCHEMA.encode())
    with h5py.File(path, "r") as f:
        g = f["entry/data"]
        h.update(np.ascontiguousarray(g["data_values"][()], dtype="<f8").tobytes())
        meta = {}
        for k in META_KEYS:
            v = g[k][()]
            meta[k] = v.decode().strip() if isinstance(v, bytes) else np.asarray(v).tolist()
        h.update(json.dumps(meta, sort_keys=True).encode())
    return h.hexdigest()


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_manifest(d: Path) -> dict:
    out = {}
    for line in (d / "MANIFEST.sha256").read_text().splitlines():
        if line.strip():
            digest, name = line.split(maxsplit=1)
            out[name.lstrip("*")] = digest
    return out


def engine_env(d: Path) -> dict:
    txt = (d / "ENGINE_ENV.txt").read_text()
    commit = re.search(r"registered_commit:\s*([0-9a-f]{7,40})", txt)
    binary = re.search(r"^([0-9a-f]{64})\s+\S*discus_suite", txt, re.M)
    return {"registered_commit": commit.group(1) if commit else None, "binary_sha256": binary.group(1) if binary else None}


def macro_outputs(mac: Path) -> list:
    return re.findall(r"^\s*outf\s+(\S+\.h5)\s*$", mac.read_text(), re.M)


def load_registry() -> dict:
    return json.loads(HASHES.read_text())


def cmd_hash(args):
    d = Path(args.dir)
    for p in sorted(d.glob("*.h5")):
        print(content_hash(p), p.name)


def cmd_list(args):
    reg = load_registry()
    for rel, ent in reg["evidence"].items():
        present = sum((REPO / rel / n).exists() for n in ent["volumes"])
        print(f"{rel}\n   {len(ent['volumes'])} volumes ({present} present locally), macros: {', '.join(ent['macros'])}"
              f"\n   registered commit {ent['registered_commit']}, engine sha256 {ent['binary_sha256'][:12]}…")


def cmd_regenerate(args):
    reg = load_registry()
    ev = Path(args.evidence_dir).resolve()
    rel = str(ev.relative_to(REPO))
    if rel not in reg["evidence"]:
        sys.exit(f"{rel} is not a registered evidence directory (see `list`)")
    ent = reg["evidence"][rel]
    out = Path(args.out_dir).resolve()
    if out == ev or ev in out.parents or out in ev.parents:
        sys.exit("out_dir must be outside the evidence directory")
    if out.exists() and any(out.iterdir()):
        sys.exit(f"{out} exists and is not empty; use a new directory")
    out.mkdir(parents=True, exist_ok=True)

    man = read_manifest(ev)
    for name in ent["inputs"]:
        if sha256_file(ev / name) != man.get(name):
            sys.exit(f"input {name} does not match MANIFEST.sha256; refusing to run")
        shutil.copy2(ev / name, out / name)

    discus = Path(args.discus)
    bsha = sha256_file(discus) if discus.exists() else None
    if bsha != ent["binary_sha256"]:
        print(f"WARNING: engine {discus} sha256 {bsha} differs from the recorded {ent['binary_sha256']}; "
              "content may differ (use --analyze)")
    status = []
    for m in ent["macros"]:
        with open(out / f"{m}.mac") as fin, open(out / f"{m}.log", "w") as flog:
            rc = subprocess.run([str(discus)], stdin=fin, stdout=flog, stderr=subprocess.STDOUT, cwd=out,
                                timeout=3600).returncode
        status.append(f"{m} exit={rc}")
        print(status[-1])
    (out / "REGENERATION.txt").write_text(
        f"source_evidence: {rel}\nengine: {discus}\nengine_sha256: {bsha}\n" + "\n".join(status) + "\n")

    ok = True
    for name, want in ent["volumes"].items():
        p = out / name
        got = content_hash(p) if p.exists() else None
        verdict = "MATCH" if got == want else ("MISSING" if got is None else "DIFFERS")
        ok &= verdict == "MATCH"
        print(f"{verdict:8s} {name}")
    if args.analyze:
        run_registered_analysis(ent, out)
    print("ALL VOLUMES MATCH the original content" if ok else
          "Not all volumes match the original content; see above (and --analyze)")
    sys.exit(0 if ok else 1)


PR_REF = "refs/pull/8/head"


def ensure_commit(commit: str):
    """Make `commit` available locally, fetching the PR #8 history from GitHub if needed."""
    have = lambda: subprocess.run(["git", "-C", str(REPO), "cat-file", "-e", f"{commit}^{{commit}}"],
                                  capture_output=True).returncode == 0
    if have():
        return
    shallow = subprocess.run(["git", "-C", str(REPO), "rev-parse", "--is-shallow-repository"],
                             capture_output=True, text=True).stdout.strip() == "true"
    cmd = ["git", "-C", str(REPO), "fetch", *(["--unshallow"] if shallow else []), "origin",
           f"+{PR_REF}:refs/remotes/origin/pr8-history"]
    print("fetching PR #8 history:", " ".join(cmd[3:]))
    subprocess.run(cmd, check=True)
    if not have():
        sys.exit(f"commit {commit} is not available even after fetching {PR_REF}")


def cmd_fetch_history(args):
    ensure_commit("8cd59ac")
    print("PR #8 history available (e.g. 8cd59ac, and every registered commit in EVIDENCE_CONTENT_HASHES.json)")


def run_registered_analysis(ent: dict, out: Path):
    """Run the evaluator of the registered commit, from a temporary detached worktree, on the regenerated dir."""
    commit = ent["registered_commit"]
    ensure_commit(commit)
    wt = Path(tempfile.mkdtemp(prefix="evidence_wt_"))
    subprocess.run(["git", "-C", str(REPO), "worktree", "add", "--detach", str(wt), commit], check=True)
    try:
        script = wt / ent["analysis"][0]
        rc = subprocess.run([sys.executable, str(script), *ent["analysis"][1:], str(out)], cwd=wt).returncode
        print(f"registered-commit analysis ({commit[:7]}) exit={rc}; reports written into {out}")
    finally:
        subprocess.run(["git", "-C", str(REPO), "worktree", "remove", "--force", str(wt)], check=False)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("list").set_defaults(fn=cmd_list)
    sub.add_parser("fetch-history").set_defaults(fn=cmd_fetch_history)
    p = sub.add_parser("hash"); p.add_argument("dir"); p.set_defaults(fn=cmd_hash)
    p = sub.add_parser("regenerate"); p.add_argument("evidence_dir"); p.add_argument("out_dir")
    p.add_argument("--discus", default="/usr/local/bin/discus_suite")
    p.add_argument("--analyze", action="store_true"); p.set_defaults(fn=cmd_regenerate)
    a = ap.parse_args(); a.fn(a)


if __name__ == "__main__":
    main()
