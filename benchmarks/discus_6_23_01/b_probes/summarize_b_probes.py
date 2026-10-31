"""Summarise one Package B probe evidence directory (read-only on the HDF5 files).

    python benchmarks/discus_6_23_01/b_probes/summarize_b_probes.py benchmarks/discus_6_23_01/b_probes/evidence/<stamp>

Writes summary.json and SUMMARY.md into that directory. Reads /entry/data with raw h5py (planes have a
singular increment matrix that DiscusVolume rejects by design). Interprets nothing: the analysis of
what the numbers mean is done separately, by the reviewer.
"""
from __future__ import annotations

import hashlib
import json
import re
import sys
from pathlib import Path

import numpy as np

LOG_PATTERNS = re.compile(r"(?i)(<F>|aver|used .*% of|lots?\b|lot +\d|friedel|symm|technique|turbo|nufft|"
                          r"temp\.|number of atoms|eval|\*\*\*[A-Z]+\*\*\*|\berror\b|warning)")


def _text(v):
    return v.decode("utf-8", "replace") if isinstance(v, (bytes, np.bytes_)) else v


def describe_h5(path: Path) -> dict:
    import h5py
    out = {"file": path.name, "bytes": path.stat().st_size}
    with h5py.File(path, "r") as f:
        if "/entry/data" not in f:
            out["error"] = "no /entry/data group"
            out["groups"] = list(f.keys())
            return out
        g = f["/entry/data"]
        meta = {}
        for k in g:
            v = g[k][()]
            if np.ndim(v) == 0:
                meta[k] = _text(v.item() if hasattr(v, "item") else v)
            elif k != "data_values":
                meta[k] = np.asarray(v).tolist()
        vals = np.asarray(g["data_values"][()], dtype=float)
    out["meta"] = meta
    out["values_shape"] = list(vals.shape)
    out["values_sha256"] = hashlib.sha256(np.ascontiguousarray(vals).tobytes()).hexdigest()
    finite = np.isfinite(vals)
    out["n_nonfinite"] = int(vals.size - finite.sum())
    if finite.any():
        v = vals[finite]
        out["stats"] = {"min": float(v.min()), "max": float(v.max()), "mean": float(v.mean()),
                        "sum": float(v.sum()), "n_negative": int((v < 0).sum()), "n_zero": int((v == 0).sum())}
    return out


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 1:
        print(__doc__)
        return 2
    ev = Path(argv[0])
    status = json.loads((ev / "runs_status.json").read_text(encoding="utf-8"))
    summary = {"stamp": status["stamp"], "engine": status["engine"], "runs": []}
    for r in status["runs"]:
        d = ev / r["id"]
        log = (d / "run.log").read_text(encoding="utf-8", errors="replace") if (d / "run.log").is_file() else ""
        lines = [l.rstrip() for l in log.splitlines() if LOG_PATTERNS.search(l) and not set(l.strip()) <= {"*"}]
        files = [describe_h5(p) for p in sorted(d.glob("*.h5"))]
        summary["runs"].append({**r, "log_lines": lines[:80], "h5": files})
    (ev / "summary.json").write_text(json.dumps(summary, indent=1) + "\n", encoding="utf-8", newline="\n")

    md = [f"# Package B probe summary — {summary['stamp']}", "",
          f"Engine: `{summary['engine'].get('executable')}` sha256 `{summary['engine'].get('executable_sha256')}`", "",
          "| run | expect | as expected | rc | s | file | shape | corner | increments | bragg | min | max | sum |",
          "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in summary["runs"]:
        if not r["h5"]:
            md.append(f"| {r['id']} | {r['expect']} | {r['as_expected']} | {r['returncode']} | {r['seconds']} "
                      f"| — | | | | | | | |")
        for h in r["h5"]:
            m, s = h.get("meta", {}), h.get("stats", {})
            md.append(f"| {r['id']} | {r['expect']} | {r['as_expected']} | {r['returncode']} | {r['seconds']} | "
                      f"{h['file']} | {h.get('values_shape')} | {m.get('data_corner')} | "
                      f"{m.get('data_increment_vector')} | {m.get('data_type_with_bragg')} | "
                      f"{s.get('min', '')!s:.10} | {s.get('max', '')!s:.10} | {s.get('sum', '')!s:.12} |")
    md += ["", "Error lines and selected log lines are in summary.json."]
    (ev / "SUMMARY.md").write_text("\n".join(md) + "\n", encoding="utf-8", newline="\n")
    print(f"wrote {ev / 'summary.json'} and SUMMARY.md ({len(summary['runs'])} runs)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
