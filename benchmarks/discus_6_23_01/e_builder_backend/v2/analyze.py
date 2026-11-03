"""Evaluate the builder-backend probe v2 (PREREG.md, criteria.json).

    python benchmarks/discus_6_23_01/e_builder_backend/v2/analyze.py <evidence dir> [--json out.json] [--v1 <v1 run dir>]
    python benchmarks/discus_6_23_01/e_builder_backend/v2/analyze.py --selftest

<evidence dir> as written by ../run_probe.sh (work/<run>/result.pdf, exit_codes.txt, input_sha256.txt, logs).
--v1: v1 evidence run directory; then I1 reports whether p_mix_n reproduces v1's p_mix byte for byte.
Exit code 0 = all criteria passed, 1 = a criterion failed, 2 = selftest failed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
CRITERIA = json.loads((HERE / "criteria.json").read_text(encoding="utf-8"))
SPEC = json.loads((HERE / "runs.json").read_text(encoding="utf-8"))
R_EXPECTED = np.arange(round(SPEC["pdf"]["rmax"] / SPEC["pdf"]["dr"]) + 1) * SPEC["pdf"]["dr"]


def max_rel(a, b) -> float:
    return float(np.abs(np.asarray(a) - np.asarray(b)).max() / np.abs(b).max())


def evaluate(g: dict) -> tuple[dict, list[str]]:
    """g: run id -> (r, G)."""
    c = SPEC["counts"]
    s_vac = c["vac"]["atoms"] / c["vac"]["rows"]
    s_void = c["void"]["atoms"] / c["void"]["rows"]
    G = {k: v[1] for k, v in g.items()}
    vals = {
        "V0": float(sum(len(r) != len(R_EXPECTED) or not np.allclose(r, R_EXPECTED, atol=1e-6) or not np.isfinite(y).all()
                        for r, y in g.values())),
        "V1": max_rel(G["p_mix_big_n"], G["p_mix_n"]),
        "V2": max_rel(G["p_mix_big_x"], G["p_mix_x"]),
        "V3": max_rel(s_vac * G["p_vac_void_n"], G["p_vac_big_n"]),
        "V4": max_rel(s_vac * G["p_vac_void_x"], G["p_vac_big_x"]),
        "V5": max_rel(s_void * G["p_void_x"], G["p_mix_x"]),
        "V6": max_rel(G["p_vac_big_n"], G["p_mix_n"]),
    }
    res, failed = {}, []
    for key, crit in CRITERIA.items():
        v = vals[key]
        ok = v <= crit["max"] if "max" in crit else v >= crit["min"]
        res[key] = dict(crit, value=v, passed=bool(ok))
        if not ok:
            failed.append(key)
    return res, failed


def check_integrity(ev: Path) -> tuple[list[str], list[str]]:
    from discuspdf.diffuse3d import _ENGINE_ERROR_PATTERN
    problems, notes = [], []
    codes = dict(l.split(" exit=") for l in (ev / "exit_codes.txt").read_text().splitlines() if " exit=" in l)
    ids = {r["id"] for r in SPEC["runs"]}
    if set(codes) != ids or any(v.strip() != "0" for v in codes.values()):
        problems.append(f"exit codes: {codes}")
    committed = {p.relative_to(HERE / "inputs").as_posix() for p in (HERE / "inputs").rglob("*") if p.is_file()}
    hashed = set()
    for line in (ev / "input_sha256.txt").read_text().splitlines():
        h, name = line.split(maxsplit=1)
        name = name.strip().removeprefix("./")
        hashed.add(name)
        f = HERE / "inputs" / name
        if not f.is_file():
            problems.append(f"input {name} is not a committed file")
        elif hashlib.sha256(f.read_bytes().replace(b"\r\n", b"\n")).hexdigest() != h:
            problems.append(f"input {name} differs from the committed file")
    if committed - hashed:
        problems.append(f"committed inputs not used: {sorted(committed - hashed)}")
    for r in SPEC["runs"]:
        log = ev / f"{r['id']}.log"
        if not log.is_file():
            problems.append(f"{r['id']}: no engine log")
            continue
        notes += [f"{r['id']}: {l.strip()}" for l in log.read_text(errors="replace").splitlines()
                  if _ENGINE_ERROR_PATTERN.search(l)]
    return problems, notes


def load(ev: Path) -> dict:
    from discuspdf.pdf import read_pdf
    return {r["id"]: read_pdf(ev / "work" / r["id"] / "result.pdf") for r in SPEC["runs"]}


def selftest() -> int:
    ok = True
    c = SPEC["counts"]
    r = R_EXPECTED
    gm = np.sin(3 * r) * np.exp(-r / 5) - 0.97 * r * (r < 1.5)
    gv = 0.94 * gm + 0.01 * np.cos(r)                        # vacancies change the shape
    s_vac, s_void = c["vac"]["atoms"] / c["vac"]["rows"], c["void"]["atoms"] / c["void"]["rows"]
    good = {"p_mix_n": gm, "p_mix_x": 1.1 * gm, "p_mix_big_n": gm, "p_mix_big_x": 1.1 * gm, "p_void_x": 1.1 * gm / s_void,
            "p_vac_void_n": gv / s_vac, "p_vac_void_x": 1.1 * gv / s_vac, "p_vac_big_n": gv, "p_vac_big_x": 1.1 * gv}
    pack = lambda d: {k: (r, v) for k, v in d.items()}   # noqa: E731
    _, failed = evaluate(pack(good))
    if failed:
        print("selftest: consistent data failed", failed)
        ok = False
    cases = {"V1": {"p_mix_big_n": gm * 1.25}, "V2": {"p_mix_big_x": 1.1 * gm * 0.8},
             "V3": {"p_vac_void_n": gv}, "V4": {"p_vac_void_x": 1.1 * gv}, "V5": {"p_void_x": 1.1 * gm},
             "V6": {"p_vac_big_n": gm, "p_vac_void_n": gm / s_vac}}
    for key, change in cases.items():
        _, failed = evaluate(pack({**good, **change}))
        if key not in failed:
            print(f"selftest: case {key} not caught ({failed})")
            ok = False
    bad = pack(good)
    bad["p_mix_n"] = (r[:-1], gm[:-1])
    try:
        _, failed = evaluate(bad)
        caught = "V0" in failed
    except ValueError:
        caught = True                                          # shape mismatch is caught too
    if not caught:
        print("selftest: short r grid not caught")
        ok = False
    with tempfile.TemporaryDirectory() as tmp:                # integrity on a fake evidence folder
        ev = Path(tmp)
        shutil.copytree(HERE / "inputs", ev / "work")
        lines = [f"{hashlib.sha256(p.read_bytes()).hexdigest()}  ./{p.relative_to(ev / 'work').as_posix()}"
                 for p in sorted((ev / "work").rglob("*")) if p.is_file()]
        (ev / "input_sha256.txt").write_text("\n".join(lines) + "\n")
        (ev / "exit_codes.txt").write_text("".join(f"{x['id']} exit=0\n" for x in SPEC["runs"]))
        for x in SPEC["runs"]:
            (ev / f"{x['id']}.log").write_text("done\n")
        p, n = check_integrity(ev)
        if p or n:
            print("selftest: clean fake evidence gave", p, n)
            ok = False
        (ev / "input_sha256.txt").write_text("\n".join(lines[1:]) + "\n")
        (ev / "exit_codes.txt").write_text("p_mix_n exit=134\n")
        p, _ = check_integrity(ev)
        if not any("exit codes" in x for x in p) or not any("not used" in x for x in p):
            print("selftest: integrity problems not caught", p)
            ok = False
    print("selftest", "PASS" if ok else "FAIL")
    return 0 if ok else 2


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("evidence", type=Path, nargs="?")
    ap.add_argument("--json", type=Path)
    ap.add_argument("--v1", type=Path)
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args(argv)
    if a.selftest:
        return selftest()
    if a.evidence is None:
        ap.error("give the evidence directory or --selftest")
    problems, notes = check_integrity(a.evidence)
    res, failed = evaluate(load(a.evidence))
    if problems:
        failed = ["E0"] + failed
        res["E0"] = {"what": "integrity", "problems": problems, "passed": False}
        print(f"E0 FAIL integrity: {problems}")
    for key, row in res.items():
        if key != "E0":
            bound = f"<= {row['max']:g}" if "max" in row else f">= {row['min']:g}"
            print(f"{key} {'PASS' if row['passed'] else 'FAIL'} {row['value']:.3e} {bound}  {row['what']}")
    res["I0"] = {"what": "engine-marker lines (informational)", "lines": notes}
    print(f"I0 informational: {len(notes)} engine-marker line(s)" + "".join(f"\n   {l}" for l in notes[:20]))
    if a.v1:
        same = (a.v1 / "work" / "p_mix" / "result.pdf").read_bytes() == (a.evidence / "work" / "p_mix_n" / "result.pdf").read_bytes()
        res["I1"] = {"what": "p_mix_n reproduces v1 p_mix byte for byte (informational)", "value": same}
        print(f"I1 informational: p_mix_n == v1 p_mix byte for byte: {same}")
    print(f"criteria failed: {len(failed)} {failed}")
    if a.json:
        a.json.write_text(json.dumps(res, indent=1, default=float) + "\n", encoding="utf-8")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
