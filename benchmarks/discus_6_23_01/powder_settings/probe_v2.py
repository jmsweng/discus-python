"""Powder settings probe v2: pre-registered test of the model diagnosed (post hoc) from the v1.1 run.

See RESULT_v1_1.md. Registered predictions for DISCUS_SUITE 6.23.01, on settings not used in v1.1:

  P1  F(Q) output: the value written at label Q was computed at Q - dQ; the first written point is a
      placeholder with S = 0, i.e. F = -Q (|F0/Q0 + 1| <= 1e-2).
  P2  G(r) = sum_j F(Q_j) sin(Q_j r) dQ (no 2/pi) on the grid j*dQ, where F is the written F from the
      second point on, placed at label - dQ, set to zero from its last sign change, and below Qmin a
      straight line F(Q) = Q * F1 / Qlabel1. No fitted parameter. Max rel deviation <= 1e-3.
  P3  F at the true Q does not depend on qmin or qmax: <= 1e-12 at common Q.
  P4  set delta > 0 aborts the engine (no result.pdf), also with OMP_NUM_THREADS=1 and on another grid.
  I1  (informational) dQ dependence of F at the true Q.
  I2  (informational) neutron radiation follows P1-P3 as well (graded as P1-P2 like the others).

Usage (WSL, repository root, venv):
    python benchmarks/discus_6_23_01/powder_settings/probe_v2.py generate
    python benchmarks/discus_6_23_01/powder_settings/probe_v2.py selftest
    python benchmarks/discus_6_23_01/powder_settings/probe_v2.py run --work ~/discus_powder/work2
    python benchmarks/discus_6_23_01/powder_settings/probe_v2.py analyze <evidence_dir>
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import platform
import re
import shutil
import signal
import subprocess
import sys
import warnings
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import probe as V1  # noqa: E402  (shared helpers: read_xy, zero_last, sha256, render_pdf_macro)

INPUTS = HERE / "inputs_v2"
PREDICTIONS = HERE / "predictions_v2.json"
STRUCTURE = V1.STRUCTURE

BASE = dict(qmin=0.35, qmax=6.83, dq=0.012, qsig=0.0, radiation="xray")  # 540 steps; all v2 values on whole steps
# id, changes, extra environment, purpose
CASES = [
    ("b", {}, {}, "new base (qmin 0.35, qmax 6.83, dq 0.012)"),
    ("b_qmin", {"qmin": 1.31}, {}, "P3: qmin (base grid + 80 steps)"),
    ("b_qmax", {"qmax": 5.51}, {}, "P3: qmax / termination (430 steps)"),
    ("b_dq", {"dq": 0.008}, {}, "I1: dq"),
    ("n", {"qmin": 0.5, "qmax": 7.5, "dq": 0.01, "radiation": "neutron"}, {}, "I2: neutron"),
    ("d", {"qsig": 0.03}, {}, "P4: set delta"),
    ("d_omp1", {"qsig": 0.03}, {"OMP_NUM_THREADS": "1"}, "P4: set delta, one thread"),
    ("d_grid", {"qsig": 0.03, "qmin": 0.5, "qmax": 7.0, "dq": 0.02}, {}, "P4: set delta, other grid"),
]
LIMITS = {"P1_placeholder": 1e-2, "P2_transform_rel_max": 1e-3, "P3_invariance_rel": 1e-12}


def case_config(over: dict):
    s = {**BASE, **over}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", V1.QmaxLimitWarning)
        return V1.PDFConfig(STRUCTURE, supercell=(5, 5, 5), rmax=12.0, dr=0.02,
                            wavelength=1.5406 if s["radiation"] == "xray" else 1.54, **s)


def case_macro(over: dict) -> str:
    m = V1.render_pdf_macro(case_config(over))
    assert m.count("\noutput\n") == 1
    return m.replace("\noutput\n", "\noutput\n" + V1.FQ_BLOCK, 1)


def cmd_generate(_):
    INPUTS.mkdir(exist_ok=True)
    macros = {}
    for cid, over, _e, _p in CASES:
        p = INPUTS / f"{cid}.mac"
        p.write_bytes(case_macro(over).encode("ascii"))
        macros[cid] = V1.sha256(p)
    PREDICTIONS.write_text(json.dumps({
        "probe": "powder-settings-v2", "engine": "DISCUS_SUITE 6.23.01", "base": BASE,
        "cases": [{"id": c, "changed": o, "env": e, "purpose": p} for c, o, e, p in CASES],
        "structure_sha256": V1.sha256(STRUCTURE), "macro_sha256": macros,
        "code_sha256": {"probe_v2.py": V1.sha256(Path(__file__)), "probe.py": V1.sha256(HERE / "probe.py")},
        "limits": LIMITS, "predictions": __doc__.split("Usage")[0].strip(),
    }, indent=1) + "\n", encoding="utf-8")
    print(f"wrote {len(macros)} macros and {PREDICTIONS.name}")


# ---------------------------------------------------------------- model
def true_f(q_label: np.ndarray, f: np.ndarray, dq: float):
    """P1: drop the placeholder, move each value from its label to label - dq."""
    return q_label[1:] - dq, f[1:]


def g_model(q_label: np.ndarray, f: np.ndarray, r: np.ndarray, dq: float) -> np.ndarray:
    qt, ft = true_f(q_label, f, dq)
    ft = V1.zero_last(ft)
    j0 = int(round(qt[0] / dq))
    qq = np.arange(j0 + len(qt)) * dq
    ff = np.empty_like(qq)
    ff[j0:] = ft
    ff[:j0] = qq[:j0] * ft[0] / q_label[1]
    return dq * np.sin(np.outer(r, qq)) @ ff


def analyze(ev: Path) -> dict:
    pred = json.loads(PREDICTIONS.read_text(encoding="utf-8"))
    crit = []

    def add(name, value, limit, ok):
        crit.append({"criterion": name, "value": None if value is None else float(value), "limit": limit,
                     "pass": None if ok is None else bool(ok)})

    data = {}
    for cid, over, _e, _p in CASES:
        d = ev / cid
        e = {"set": {**BASE, **over}, "has_pdf": (d / "result.pdf").exists(), "has_fq": (d / "fq.dat").exists()}
        if e["has_pdf"]:
            e["r"], e["g"] = V1.read_xy(d / "result.pdf")
        if e["has_fq"]:
            e["q"], e["f"] = V1.read_xy(d / "fq.dat")
        data[cid] = e
    for cid, e in data.items():
        if cid.startswith("d"):
            add(f"P4 {cid}: engine aborts, no G(r)", None, None, not e["has_pdf"])
            continue
        if not (e["has_pdf"] and e["has_fq"]):
            add(f"R {cid}: engine produced F(Q) and G(r)", None, None, False)
            continue
        dq = e["set"]["dq"]
        v = abs(e["f"][0] / e["q"][0] + 1)
        add(f"P1 {cid}: first F point is the S = 0 placeholder, |F0/Q0 + 1|", v, LIMITS["P1_placeholder"],
            v <= LIMITS["P1_placeholder"])
        gm = g_model(e["q"], e["f"], e["r"], dq)
        v = np.max(np.abs(e["g"] - gm)) / np.max(np.abs(e["g"]))
        add(f"P2 {cid}: G vs model (no fit)", v, LIMITS["P2_transform_rel_max"], v <= LIMITS["P2_transform_rel_max"])
        e["model_scale_fit"] = float(gm @ e["g"] / (gm @ gm))
    b = data["b"]
    if b["has_fq"]:
        qb, fb = true_f(b["q"], b["f"], b["set"]["dq"])
        for cid, kind in (("b_qmin", "P3"), ("b_qmax", "P3"), ("b_dq", "I1")):
            e = data[cid]
            if not e["has_fq"]:
                continue
            q, f = true_f(e["q"], e["f"], e["set"]["dq"])
            step = 0.002  # true-Q grids 0.338 + 0.012k and 0.342 + 0.008k are both on 0.002
            kb = {int(round(x / step)): y for x, y in zip(qb, fb)}
            pairs = np.array([(kb[int(round(x / step))], y) for x, y in zip(q, f) if int(round(x / step)) in kb])
            v = np.max(np.abs(pairs[:, 0] - pairs[:, 1])) / np.max(np.abs(fb)) if len(pairs) else float("inf")
            if kind == "P3":
                add(f"P3 {cid}: F at true Q vs b, {len(pairs)} common points", v, LIMITS["P3_invariance_rel"],
                    v <= LIMITS["P3_invariance_rel"])
            else:
                add(f"I1 {cid}: F at true Q vs b (dq 0.008 vs 0.012), {len(pairs)} common points (informational)",
                    v, None, None)
    graded = [c for c in crit if c["pass"] is not None]
    res = {"probe": pred["probe"], "criteria": crit, "n_graded": len(graded),
           "n_failed": sum(not c["pass"] for c in graded),
           "model_scale_fit": {c: e.get("model_scale_fit") for c, e in data.items() if "model_scale_fit" in e}}
    res["verdict"] = "PASS" if res["n_failed"] == 0 else "FAIL"
    return res


def table(res: dict) -> str:
    out = [f"# Powder settings probe ({res['probe']}): {res['verdict']} ({res['n_failed']}/{res['n_graded']} failed)", "",
           "| criterion | value | limit | pass |", "|---|---|---|---|"]
    for c in res["criteria"]:
        val = "–" if c["value"] is None else f"{c['value']:.3e}"
        lim = "–" if c["limit"] is None else f"{c['limit']:.0e}"
        ok = "info" if c["pass"] is None else ("PASS" if c["pass"] else "FAIL")
        out.append(f"| {c['criterion']} | {val} | {lim} | {ok} |")
    out += ["", "Least-squares scale of the model onto G(r) (1 expected, not graded): " +
            ", ".join(f"{k} {v:.5f}" for k, v in res["model_scale_fit"].items())]
    return "\n".join(out) + "\n"


def cmd_analyze(a):
    ev = Path(a.evidence).resolve()
    res = analyze(ev)
    (ev / "results.json").write_text(json.dumps(res, indent=1) + "\n", encoding="utf-8")
    (ev / "RESULT_TABLE.md").write_text(table(res), encoding="utf-8")
    print(table(res))
    return 0 if res["verdict"] == "PASS" else 1


# ---------------------------------------------------------------- run (owner's machine)
def cmd_run(a):
    pred = json.loads(PREDICTIONS.read_text(encoding="utf-8"))
    for cid, h in pred["macro_sha256"].items():
        if V1.sha256(INPUTS / f"{cid}.mac") != h:
            sys.exit(f"STOP: inputs_v2/{cid}.mac differs from predictions_v2.json")
    for name, h in pred["code_sha256"].items():
        if V1.sha256(HERE / name) != h:
            sys.exit(f"STOP: {name} differs from the registered version")
    if V1.sha256(STRUCTURE) != pred["structure_sha256"]:
        sys.exit("STOP: examples/gaas.cll differs from predictions_v2.json")
    exe = shutil.which("discus_suite")
    if not exe:
        sys.exit("STOP: discus_suite not on PATH")
    run_id = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "_" + platform.node()
    work = Path(a.work).expanduser().resolve() / run_id
    ev = HERE / "evidence_v2" / run_id
    if work.exists() or ev.exists():
        sys.exit(f"STOP: {work} or {ev} exists")
    work.mkdir(parents=True)
    status = []
    for cid, _o, extra, _p in CASES:
        d = work / cid
        d.mkdir()
        shutil.copyfile(INPUTS / f"{cid}.mac", d / "run.mac")
        shutil.copyfile(STRUCTURE, d / "model.cll")
        env = {k: v for k, v in os.environ.items() if k != "DISPLAY"}
        env["PGPLOT"] = "none"
        env.update(extra)
        with open(d / "run.log", "w", encoding="utf-8") as log:
            p = subprocess.Popen([exe, "-noautorun", "-macro", "run.mac"], cwd=d, stdin=subprocess.PIPE,
                                 stdout=log, stderr=subprocess.STDOUT, text=True, env=env, start_new_session=True)
            try:
                p.communicate("exit\n", timeout=600)
            finally:
                try:
                    os.killpg(p.pid, signal.SIGTERM)
                except (ProcessLookupError, PermissionError):
                    pass
        text = (d / "run.log").read_text(encoding="utf-8", errors="replace")
        err = re.search(r"\*\*\*[A-Z]+\*\*\*|[Ss]egmentation fault|double free|SIGABRT", text)
        outs = [n for n in ("fq.dat", "result.pdf") if (d / n).exists()]
        status.append({"id": cid, "env": extra, "exit": p.returncode, "engine_error": err.group(0) if err else None,
                       "outputs": outs})
        print(cid, "exit", p.returncode, err.group(0) if err else "", outs)
    shutil.copytree(work, ev)
    head = subprocess.run(["git", "-C", str(HERE), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    env_txt = [f"run_id: {run_id}", f"git_head: {head}", f"discus_path: {exe}",
               f"discus_sha256: {V1.sha256(Path(exe).resolve())}", f"platform: {platform.platform()}",
               f"nproc: {os.cpu_count()}", f"python: {sys.version.split()[0]}  numpy {np.__version__}"]
    (ev / "ENGINE_ENV.txt").write_text("\n".join(env_txt) + "\n", encoding="utf-8")
    (ev / "RUN_STATUS.json").write_text(json.dumps(status, indent=1) + "\n", encoding="utf-8")
    files = sorted(p for p in ev.rglob("*") if p.is_file())
    (ev / "MANIFEST.sha256").write_text(
        "".join(f"{V1.sha256(p)}  {p.relative_to(ev).as_posix()}\n" for p in files), encoding="ascii")
    print("\n".join(env_txt))
    print(f"evidence: {ev}")
    return cmd_analyze(argparse.Namespace(evidence=str(ev)))  # P4 expects some runs to abort: always analyse


# ---------------------------------------------------------------- selftest
def _synthetic(ev: Path, *, scale=1.0, no_shift=False, aborts=True, qmin_changes_f=False):
    shells = np.array([2.448, 3.998, 4.688, 5.654, 6.161, 6.925, 7.996, 8.481])
    amps = np.array([4, 12, 12, 6, 12, 24, 6, 12]) / shells
    r = np.round(np.arange(601) * 0.02, 10)
    for cid, over, _e, _p in CASES:
        s = {**BASE, **over}
        d = ev / cid
        d.mkdir(parents=True)
        if cid.startswith("d") and aborts:
            continue
        dq = s["dq"]
        qlab = np.round(s["qmin"] + np.arange(int(round((s["qmax"] - s["qmin"]) / dq)) + 1) * dq, 10)
        qtrue = qlab - (0 if no_shift else dq)
        f = (amps * np.sin(np.outer(qtrue, shells)) * np.exp(-0.004 * qtrue[:, None] ** 2)).sum(1)
        if qmin_changes_f and cid == "b_qmin":
            f = f * 1.001
        f[0] = -qlab[0]
        if no_shift:  # an engine without the label offset: G is the transform at the written Q
            ft = V1.zero_last(f[1:])
            j0 = int(round(qlab[1] / dq))
            qq = np.arange(j0 + len(ft)) * dq
            ff = np.empty_like(qq)
            ff[j0:] = ft
            ff[:j0] = qq[:j0] * ft[0] / qlab[1]
            g = dq * np.sin(np.outer(r, qq)) @ ff
        else:
            g = scale * g_model(qlab, f, r, dq)
        np.savetxt(d / "result.pdf", np.c_[r, g])
        np.savetxt(d / "fq.dat", np.c_[qlab, f])


def cmd_selftest(_):
    import tempfile
    cases = {"clean": ({}, "PASS"), "scale_pi_over_2_missing": ({"scale": 2 / np.pi}, "FAIL"),
             "no_label_offset": ({"no_shift": True}, "FAIL"), "delta_does_not_abort": ({"aborts": False}, "FAIL"),
             "qmin_changes_F": ({"qmin_changes_f": True}, "FAIL")}
    ok = True
    for name, (kw, want) in cases.items():
        with tempfile.TemporaryDirectory() as t:
            _synthetic(Path(t), **kw)
            res = analyze(Path(t))
        good = res["verdict"] == want
        ok &= good
        bad = [c["criterion"] for c in res["criteria"] if c["pass"] is False][:2]
        print(f"{'ok ' if good else 'BAD'} {name}: {res['verdict']} (expected {want}) {bad}")
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)
    sp.add_parser("generate").set_defaults(fn=cmd_generate)
    sp.add_parser("selftest").set_defaults(fn=cmd_selftest)
    p = sp.add_parser("run"); p.add_argument("--work", required=True); p.set_defaults(fn=cmd_run)
    p = sp.add_parser("analyze"); p.add_argument("evidence"); p.set_defaults(fn=cmd_analyze)
    a = ap.parse_args()
    sys.exit(a.fn(a) or 0)


if __name__ == "__main__":
    main()
