"""Compare a DISCUS ``chem corr`` probe run with the registered predictions.

Usage: python benchmarks/discus_6_23_01/c_probes/analyze.py [EVIDENCE_DIR]

Reads ``<run>.log`` for occ_bcc, disp_sc and mixed_bcc, takes every ``corr`` block
(the table rows "Neig. ... # pairs  correlation" and the following ``eval res[k]``
values, which DISCUS prints with 12 significant digits), checks the input SHA-256
values against ``predictions.json`` and prints one table row per neighborhood.
Exit status 0 if every correlation agrees to 1e-9 (absolute) and every pair count
agrees exactly, 1 otherwise.
"""
from __future__ import annotations

import json
import math
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
TOL = 1e-9

_CORR = re.compile(r"discus/chem> corr (occ|disp), *(\w+), *(\w+)")
_ROW = re.compile(r"^\s+(\d+)\s.*?\s(\d+)\s+(-?\d+\.\d+)\s*$")
_RES = re.compile(r"Value of res\[(\d+)\] =\s+(\S+)")


def parse_log(text: str) -> list[dict]:
    """[{kind, a, b, rows: {k: n_pairs}, res: {k: value}}] in order of the corr commands."""
    blocks = []
    for line in text.splitlines():
        m = _CORR.search(line)
        if m:
            blocks.append({"kind": m.group(1), "a": m.group(2), "b": m.group(3), "rows": {}, "res": {}})
            continue
        if not blocks:
            continue
        r = _RES.search(line)
        if r:
            blocks[-1]["res"][int(r.group(1))] = float(r.group(2))
            continue
        row = _ROW.match(line)
        if row and ("%" in line or blocks[-1]["kind"] == "disp"):
            blocks[-1]["rows"][int(row.group(1))] = int(row.group(2))
    return blocks


def main(evidence: Path) -> int:
    preds = json.loads((HERE / "predictions.json").read_text(encoding="ascii"))
    ok = True
    lines = []
    sha = {}
    for ln in (evidence / "input_sha256.txt").read_text().splitlines():
        h, name = ln.split()
        sha[name] = h
    bad_hash = [k for k, v in preds["input_sha256"].items() if sha.get(k) != v]
    lines.append(f"input SHA-256: {len(preds['input_sha256']) - len(bad_hash)}/{len(preds['input_sha256'])} match"
                 + (f" (mismatch: {bad_hash})" if bad_hash else ""))
    ok &= not bad_hash

    worst = 0.0

    def check_blocks(run, blocks, expected, n_neigh):
        """Every expected corr block present, in order, with every neighborhood's table row and res value."""
        got = [(b["kind"], b["a"].upper(), b["b"].upper()) for b in blocks]
        want = [(k, a.upper(), b.upper()) for k, a, b in expected]
        good = got == want and all(set(b["rows"]) >= set(range(1, n_neigh + 1)) and
                                   set(b["res"]) >= set(range(1, n_neigh + 1)) for b in blocks)
        lines.append(f"{run}: corr blocks {got} (expected {want}), "
                     f"{'complete' if good else 'INCOMPLETE'}")
        return good

    # occupational
    blocks = parse_log((evidence / "occ_bcc.log").read_text(errors="replace"))
    nbh = preds["occ_bcc"]["neighborhoods"]
    ok &= check_blocks("occ_bcc", blocks, [("occ", "Au", "Cu"), ("occ", "Cu", "Au")], len(nbh))
    for blk, key in zip(blocks, ("corr_AuCu", "corr_CuAu")):
        for k, p in enumerate(nbh, start=1):
            want = p[key]
            got = blk["res"].get(k, math.nan)
            d = abs(got - want)
            worst = max(worst, d)
            n_want = p.get("n_pairs_discus", p.get("n_pairs_ours"))
            n_got = blk["rows"].get(k)
            good = d <= TOL and n_got == n_want
            ok &= good
            lines.append(f"occ_bcc  corr occ {blk['a']},{blk['b']}  #{k} {p['kind']:3s} {str(p['def']):32s} "
                         f"pairs {n_got}/{n_want}  DISCUS {got:+.12f}  ours {want:+.12f}  |d| {d:.1e}  "
                         f"{'OK' if good else 'FAIL'}")
    # displacement
    for run, pairs in (("disp_sc", ["Cu_Cu"]), ("mixed_bcc", ["Au_Cu", "Cu_Cu", "all_all"])):
        blocks = parse_log((evidence / f"{run}.log").read_text(errors="replace"))
        nbh = preds[run]["neighborhoods"]
        ok &= check_blocks(run, blocks, [("disp",) + tuple(p.split("_")) for p in pairs], len(nbh))
        for blk, key in zip(blocks, pairs):
            for k, p in enumerate(nbh, start=1):
                want = p[key]["discus_convention"]
                got = blk["res"].get(k, math.nan)
                d = abs(got - want)
                worst = max(worst, d)
                n_got = blk["rows"].get(k)
                good = d <= TOL and n_got == p[key]["n_pairs"]
                ok &= good
                lines.append(f"{run:9s} corr disp {blk['a']},{blk['b']}  #{k} {str(p['def']):48s} "
                             f"pairs {n_got}/{p[key]['n_pairs']}  DISCUS {got:+.12f}  ours {want:+.12f}  "
                             f"(Pearson {p[key]['pearson']:+.6f})  |d| {d:.1e}  {'OK' if good else 'FAIL'}")
    lines.append(f"largest |DISCUS - ours| = {worst:.2e} (tolerance {TOL:g}; DISCUS prints 12 significant digits)")
    lines.append("PASS" if ok else "FAIL")
    print("\n".join(lines))
    return 0 if ok else 1


if __name__ == "__main__":
    default = sorted((HERE / "evidence").iterdir())[-1]
    sys.exit(main(Path(sys.argv[1]) if len(sys.argv) > 1 else default))
