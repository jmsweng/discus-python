"""Powder-PDF settings probe: what `set delta`, Qmin, DQ and Qmax termination do in DISCUS 6.23.01.

TASKS.md P1 "Verify powder settings against real DISCUS outputs". Qmax and wavelength were measured on
2026-10-31 (docs/VERIFIED_DISCUS_MACRO.md); this probe covers the rest on the same GaAs macro.

Every case is the notebook-01 macro (``render_pdf_macro``; GaAs 5x5x5, X-ray, lambda 1.5406 A, rmax 12,
dr 0.02) with one setting changed from the base (qmin 0.5, qmax 7.5, dq 0.01, qsig 0). Qmax 7.5 is below
the 8 A^-1 limit, so termination is real. Each macro also writes DISCUS's F(Q) (``form powder, q`` /
``value F(Q)``) before the unchanged PDF block, so every G(r) can be checked against an independent
sine transform of the F(Q) of the same run:

    G_ref(r) = (2/pi) * sum_j F(Q_j) sin(Q_j r) dQ

with the two rules read from the 6.23.01 source (powder_write_mod.f90 ``zero_last``, discus_powder_fft.f90
``fft_fq``): F(Q) is set to zero from its last sign change up to Qmax, and below Qmin it is replaced by a
straight line to F(0) = 0. Whether the transform places F at its written Q or one step lower (the source
copies ywrt(1:n) to index nint(qmin/dq)) is not clear from the source alone, so both placements are
evaluated and the result states which one the engine follows (post hoc, labelled as such).

Usage (WSL, repository root, venv with numpy):
    python benchmarks/discus_6_23_01/powder_settings/probe.py generate      # rewrites inputs/ + predictions.json
    python benchmarks/discus_6_23_01/powder_settings/probe.py selftest      # analyzer on synthetic data
    python benchmarks/discus_6_23_01/powder_settings/probe.py run --work ~/discus_powder/work
    python benchmarks/discus_6_23_01/powder_settings/probe.py analyze <evidence_dir>

``run`` never edits inputs and never deletes anything; it refuses an existing work or evidence directory.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
import platform
import signal
import re
import shutil
import subprocess
import sys
import warnings
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO / "src"))
from discuspdf import PDFConfig, QmaxLimitWarning, SetDeltaWarning, render_pdf_macro  # noqa: E402

STRUCTURE = REPO / "examples" / "gaas.cll"
INPUTS = HERE / "inputs"
PREDICTIONS = HERE / "predictions.json"

BASE = dict(qmin=0.5, qmax=7.5, dq=0.01, qsig=0.0)
# id, changed setting(s), purpose
CASES = [
    ("base", {}, "reference case"),
    ("base_plain", {}, "same as base, without the F(Q) block: the extra output must not change G(r)"),
    ("delta002", {"qsig": 0.02}, "set delta"),
    ("delta005", {"qsig": 0.05}, "set delta"),
    ("delta010", {"qsig": 0.10}, "set delta"),
    ("qmin02", {"qmin": 0.2}, "Qmin"),
    ("qmin10", {"qmin": 1.0}, "Qmin"),
    ("qmin20", {"qmin": 2.0}, "Qmin above the first GaAs reflection (111), Q = 1.925"),
    ("dq005", {"dq": 0.005}, "DQ"),
    ("dq020", {"dq": 0.02}, "DQ"),
    ("dq050", {"dq": 0.05}, "DQ"),
    ("qmax50", {"qmax": 5.0}, "termination at another Qmax"),
    ("qmax60", {"qmax": 6.0}, "termination at another Qmax"),
]

FQ_BLOCK = "outf fq.dat\nform powder, q\nvalue F(Q)\nrun\n"

# Registered limits (before any engine run of this probe).
LIMITS = {
    "T1_transform_rel_max": 2e-3,   # max|G - G_ref| / max|G|, best placement variant, with zeroing
    "T1_scale_abs": 2e-3,           # |c - 1|, c = least-squares scale of G_ref onto G
    "T3_fq_invariance_rel": 1e-6,   # F(Q) at common Q unchanged by qmin, dq, qmax
    "T4_parity_rel": 1e-12,         # base vs base_plain G(r)
}
# Informational (approximate theory, no pass/fail): Gaussian `set delta` damps G(r) by exp(-delta^2 r^2 / 4).
INFO = {"I1_delta_damping": "max|G_delta - G_0 exp(-delta^2 r^2/4)| / max|G_0|"}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def case_config(over: dict) -> PDFConfig:
    s = {**BASE, **over}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", QmaxLimitWarning)
        warnings.simplefilter("ignore", SetDeltaWarning)  # the set delta cases are the point of the probe
        return PDFConfig(STRUCTURE, supercell=(5, 5, 5), rmax=12.0, dr=0.02, radiation="xray",
                         wavelength=1.5406, allow_qsig=s["qsig"] > 0, **s)


def case_macro(cid: str, over: dict) -> str:
    m = render_pdf_macro(case_config(over))
    if cid == "base_plain":
        return m
    assert m.count("\noutput\n") == 1
    return m.replace("\noutput\n", "\noutput\n" + FQ_BLOCK, 1)


# ---------------------------------------------------------------- generate
def cmd_generate(_):
    INPUTS.mkdir(exist_ok=True)
    macros = {}
    for cid, over, _p in CASES:
        p = INPUTS / f"{cid}.mac"
        p.write_bytes(case_macro(cid, over).encode("ascii"))
        macros[cid] = sha256(p)
    pred = {
        "probe": "powder-settings-v1",
        "engine": "DISCUS_SUITE 6.23.01",
        "base": BASE,
        "cases": [{"id": c, "changed": o, "purpose": p} for c, o, p in CASES],
        "structure_sha256": sha256(STRUCTURE),
        "macro_sha256": macros,
        "code_sha256": {"probe.py": sha256(Path(__file__))},
        "limits": LIMITS,
        "informational": INFO,
        "expectations": {
            "T1": "every case: G(r) = (2/pi) sum F sin(Qr) dQ of its own F(Q), F zeroed from its last sign change, "
                  "linear below Qmin; at least one Q placement variant within the limits",
            "T3": "F(Q) is a pointwise Debye sum: at common Q it does not depend on qmin, dq or qmax",
            "T4": "the F(Q) output block does not change G(r)",
            "I1": "set delta,d convolves along Q with exp(-x^2/d^2)/(sqrt(pi) d), i.e. G(r) ~ G_0(r) exp(-d^2 r^2/4)",
        },
    }
    PREDICTIONS.write_text(json.dumps(pred, indent=1) + "\n", encoding="ascii")
    print(f"wrote {len(macros)} macros and {PREDICTIONS.name}")


# ---------------------------------------------------------------- reading
def read_xy(path: Path) -> tuple[np.ndarray, np.ndarray]:
    rows = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        parts = line.replace("D+", "E+").replace("D-", "E-").split()
        if len(parts) < 2:
            continue
        try:
            rows.append((float(parts[0]), float(parts[1])))
        except ValueError:
            continue
    a = np.asarray(rows, float)
    if a.ndim != 2 or len(a) < 2:
        raise ValueError(f"no (x, y) data in {path}")
    return a[:, 0], a[:, 1]


# ---------------------------------------------------------------- reference transform
def zero_last(y: np.ndarray) -> np.ndarray:
    """DISCUS 6.23.01 powder_write_mod.f90 ``zero_last``: walk down from the last point, zero each point
    until a sign change (the upper point of the pair is zeroed too) or an exact zero is met."""
    y = y.copy()
    for i in range(len(y) - 1, 0, -1):
        if y[i] * y[i - 1] < 0:
            y[i] = 0.0
            break
        if y[i] == 0.0:
            break
        y[i] = 0.0
    return y


def g_ref(q: np.ndarray, f: np.ndarray, r: np.ndarray, dq: float, *, shift: int, zero: bool) -> np.ndarray:
    """(2/pi) sum F sin(Q r) dQ on the grid j*dq, F linear to 0 below its first point.

    shift=0: F placed at its written Q. shift=1: placed one step lower (Q - dq)."""
    f = zero_last(f) if zero else f.copy()
    q = q - shift * dq
    j0 = int(round(q[0] / dq))
    qq = np.arange(j0 + len(q)) * dq
    ff = np.empty_like(qq)
    ff[j0:] = f
    ff[:j0] = qq[:j0] * f[0] / q[0] if q[0] > 0 else 0.0
    return (2.0 / np.pi) * dq * np.sin(np.outer(r, qq)) @ ff


def transform_check(q, f, r, g, dq) -> dict:
    out = {}
    for shift in (0, 1):
        for zero in (True, False):
            ref = g_ref(q, f, r, dq, shift=shift, zero=zero)
            c = float(ref @ g / (ref @ ref))
            out[f"shift{shift}_{'zero' if zero else 'nozero'}"] = {
                "rel_max": float(np.max(np.abs(g - ref)) / np.max(np.abs(g))),
                "scale": c,
            }
    best = min((k for k in out if k.endswith("_zero")), key=lambda k: out[k]["rel_max"])
    return {"variants": out, "best_with_zeroing": best}


def common(q1, f1, q2, f2, step=0.005):
    k1 = {int(round(x / step)): v for x, v in zip(q1, f1)}
    k2 = {int(round(x / step)): v for x, v in zip(q2, f2)}
    keys = sorted(set(k1) & set(k2))
    a = np.array([k1[k] for k in keys])
    b = np.array([k2[k] for k in keys])
    return len(keys), a, b


# ---------------------------------------------------------------- analyze
def analyze(ev: Path) -> dict:
    pred = json.loads(PREDICTIONS.read_text())
    data = {}
    failed_runs = []
    for cid, over, _p in CASES:
        d = ev / cid
        if not (d / "result.pdf").exists():  # engine produced no output (v1.1 run: set delta aborts)
            failed_runs.append(cid)
            continue
        r, g = read_xy(d / "result.pdf")
        entry = {"r": r, "g": g, "set": {**BASE, **over}}
        if (d / "fq.dat").exists():
            entry["q"], entry["f"] = read_xy(d / "fq.dat")
        data[cid] = entry
    res = {"probe": pred["probe"], "limits": LIMITS, "cases": {}, "criteria": [], "effects": {},
           "engine_failures": failed_runs}
    crit = res["criteria"]
    for cid in failed_runs:
        crit.append({"criterion": f"R {cid}: engine produced G(r)", "value": 0.0, "limit": None, "pass": False})

    def add(name, value, limit, ok, info=False):
        crit.append({"criterion": name, "value": value, "limit": limit,
                     "pass": None if info else bool(ok)})

    base = data["base"]
    for cid, e in data.items():
        cres = {"settings": e["set"], "n_r": int(len(e["r"]))}
        if "q" in e:
            dq = e["set"]["dq"]
            cres["fq_first_last"] = [float(e["q"][0]), float(e["q"][-1])]
            cres["fq_points"] = int(len(e["q"]))
            cres["fq_written_step"] = float(np.median(np.diff(e["q"])))
            t = transform_check(e["q"], e["f"], e["r"], e["g"], dq)
            cres["transform"] = t
            b = t["variants"][t["best_with_zeroing"]]
            add(f"T1 {cid}: G vs sine transform of own F(Q) ({t['best_with_zeroing']})", b["rel_max"],
                LIMITS["T1_transform_rel_max"], b["rel_max"] <= LIMITS["T1_transform_rel_max"])
            add(f"T1 {cid}: |scale - 1|", abs(b["scale"] - 1), LIMITS["T1_scale_abs"],
                abs(b["scale"] - 1) <= LIMITS["T1_scale_abs"])
        if cid not in ("base", "base_plain"):
            cres["effect_on_G_rel_max"] = float(np.max(np.abs(e["g"] - base["g"])) / np.max(np.abs(base["g"])))
            res["effects"][cid] = cres["effect_on_G_rel_max"]
        res["cases"][cid] = cres
    # T3: F(Q) invariance
    for cid in ("qmin02", "qmin10", "qmin20", "dq005", "dq020", "dq050", "qmax50", "qmax60"):
        e = data[cid]
        n, a, b = common(e["q"], e["f"], base["q"], base["f"])
        v = float(np.max(np.abs(a - b)) / np.max(np.abs(base["f"]))) if n else float("inf")
        res["cases"][cid]["fq_common_points_with_base"] = n
        add(f"T3 {cid}: F(Q) at {n} common Q vs base", v, LIMITS["T3_fq_invariance_rel"],
            n > 0 and v <= LIMITS["T3_fq_invariance_rel"])
    # T4: parity
    v = float(np.max(np.abs(data["base_plain"]["g"] - base["g"])) / np.max(np.abs(base["g"])))
    add("T4 base_plain vs base: G(r)", v, LIMITS["T4_parity_rel"], v <= LIMITS["T4_parity_rel"])
    # I1: delta damping (informational)
    for cid in [c for c in ("delta002", "delta005", "delta010") if c in data]:
        d = data[cid]["set"]["qsig"]
        pred_g = base["g"] * np.exp(-(d * base["r"]) ** 2 / 4)
        v = float(np.max(np.abs(data[cid]["g"] - pred_g)) / np.max(np.abs(base["g"])))
        res["cases"][cid]["damping_model_rel_max"] = v
        res["cases"][cid]["damping_at_rmax"] = float(np.exp(-(d * base["r"][-1]) ** 2 / 4))
        add(f"I1 {cid}: G vs G_0 exp(-d^2 r^2/4) (informational)", v, None, None, info=True)
    graded = [c for c in crit if c["pass"] is not None]
    res["n_failed"] = sum(not c["pass"] for c in graded)
    res["n_graded"] = len(graded)
    res["verdict"] = "PASS" if res["n_failed"] == 0 else "FAIL"
    return res


def table(res: dict) -> str:
    lines = [f"# Powder settings probe ({res['probe']}): {res['verdict']} ({res['n_failed']}/{res['n_graded']} failed)", "",
             "| criterion | value | limit | pass |", "|---|---|---|---|"]
    for c in res["criteria"]:
        lim = "–" if c["limit"] is None else f"{c['limit']:.0e}"
        if c["criterion"].startswith("R "):
            lines.append(f"| {c['criterion']} | no output | – | FAIL |")
            continue
        ok = "info" if c["pass"] is None else ("PASS" if c["pass"] else "FAIL")
        lines.append(f"| {c['criterion']} | {c['value']:.3e} | {lim} | {ok} |")
    lines += ["", "## Effect of each setting on G(r): max|G - G_base| / max|G_base|", "", "| case | settings | effect |", "|---|---|---|"]
    for cid, v in res["effects"].items():
        s = res["cases"][cid]["settings"]
        lines.append(f"| {cid} | qmin {s['qmin']}, qmax {s['qmax']}, dq {s['dq']}, delta {s['qsig']} | {v:.3e} |")
    return "\n".join(lines) + "\n"


def cmd_analyze(a):
    ev = Path(a.evidence).resolve()
    res = analyze(ev)
    (ev / "results.json").write_text(json.dumps(res, indent=1) + "\n", encoding="utf-8")
    (ev / "RESULT_TABLE.md").write_text(table(res), encoding="utf-8")
    print(table(res))
    return 0 if res["verdict"] == "PASS" else 1


# ---------------------------------------------------------------- run (owner's machine)
def cmd_run(a):
    pred = json.loads(PREDICTIONS.read_text())
    for cid, h in pred["macro_sha256"].items():
        if sha256(INPUTS / f"{cid}.mac") != h:
            sys.exit(f"STOP: inputs/{cid}.mac differs from predictions.json")
    if sha256(Path(__file__)) != pred["code_sha256"]["probe.py"]:
        sys.exit("STOP: probe.py differs from the registered version in predictions.json")
    if sha256(STRUCTURE) != pred["structure_sha256"]:
        sys.exit("STOP: examples/gaas.cll differs from predictions.json")
    exe = shutil.which("discus_suite")
    if not exe:
        sys.exit("STOP: discus_suite not on PATH")
    run_id = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "_" + platform.node()
    work = Path(a.work).expanduser().resolve() / run_id
    ev = HERE / "evidence" / run_id
    if work.exists() or ev.exists():
        sys.exit(f"STOP: {work} or {ev} exists")
    work.mkdir(parents=True)
    status = []
    for cid, _o, _p in CASES:
        d = work / cid
        d.mkdir()
        shutil.copyfile(INPUTS / f"{cid}.mac", d / "run.mac")
        shutil.copyfile(STRUCTURE, d / "model.cll")
        # Same environment as PDFRunner (no DISPLAY, PGPLOT=none). The log goes to a file, not a pipe:
        # v1 hung after base_plain because a process left behind by DISCUS kept the output pipe open.
        env = {k: v for k, v in os.environ.items() if k != "DISPLAY"}
        env["PGPLOT"] = "none"
        with open(d / "run.log", "w", encoding="utf-8") as log:
            p = subprocess.Popen([exe, "-noautorun", "-macro", "run.mac"], cwd=d, stdin=subprocess.PIPE,
                                 stdout=log, stderr=subprocess.STDOUT, text=True, env=env,
                                 start_new_session=True)
            try:
                p.communicate("exit\n", timeout=600)
            finally:
                try:  # stop anything this run left behind in its own session
                    os.killpg(p.pid, signal.SIGTERM)
                except (ProcessLookupError, PermissionError):
                    pass
        text = (d / "run.log").read_text(encoding="utf-8", errors="replace")
        err = re.search(r"\*\*\*[A-Z]+\*\*\*|[Ss]egmentation fault", text)
        outs = [n for n in ("fq.dat", "result.pdf") if (d / n).exists()]
        status.append({"id": cid, "exit": p.returncode, "engine_error": err.group(0) if err else None, "outputs": outs})
        print(cid, "exit", p.returncode, "error" if err else "", outs)
    shutil.copytree(work, ev)
    env = [f"run_id: {run_id}", f"git_head: {subprocess.run(['git', '-C', str(REPO), 'rev-parse', 'HEAD'], capture_output=True, text=True).stdout.strip()}",
           f"discus_path: {exe}", f"discus_sha256: {sha256(Path(exe).resolve())}", f"platform: {platform.platform()}",
           f"python: {sys.version.split()[0]}  numpy {np.__version__}"]
    (ev / "ENGINE_ENV.txt").write_text("\n".join(env) + "\n", encoding="utf-8")
    (ev / "RUN_STATUS.json").write_text(json.dumps(status, indent=1) + "\n", encoding="utf-8")
    files = sorted(p for p in ev.rglob("*") if p.is_file())
    (ev / "MANIFEST.sha256").write_text("".join(f"{sha256(p)}  {p.relative_to(ev).as_posix()}\n" for p in files), encoding="ascii")
    print("\n".join(env))
    print(f"evidence: {ev}")
    need = {cid: ({"result.pdf"} if cid == "base_plain" else {"result.pdf", "fq.dat"}) for cid, _o, _p in CASES}
    if any(s["exit"] != 0 or s["engine_error"] or not need[s["id"]] <= set(s["outputs"]) for s in status):
        print("STOP: at least one run failed; see RUN_STATUS.json. Do not analyze; report.")
        return 1
    return cmd_analyze(argparse.Namespace(evidence=str(ev)))


# ---------------------------------------------------------------- selftest (no engine)
def _synthetic(ev: Path, *, scale=1.0, shift=1, zero=True, perturb_plain=0.0, fq_noise=0.0):
    """Fake evidence: F(Q) = damped sines (GaAs-like shells) on each case's grid; 'DISCUS' G = g_ref(F).

    For the delta cases each shell is damped by exp(-d^2 r_k^2 / 4), so their own F transforms to a
    damped G, as with the real Gaussian convolution."""
    rng = np.random.default_rng(7)
    shells = np.array([2.448, 3.998, 4.688, 5.654, 6.161, 6.925, 7.996, 8.481])
    amps = np.array([4, 12, 12, 6, 12, 24, 6, 12]) / shells
    r = np.round(np.arange(601) * 0.02, 10)
    for cid, over, _p in CASES:
        s = {**BASE, **over}
        q = np.round(s["qmin"] + np.arange(int(round((s["qmax"] - s["qmin"]) / s["dq"])) + 1) * s["dq"], 10)
        damp = np.exp(-(s["qsig"] * shells) ** 2 / 4)
        f = (amps * damp * np.sin(np.outer(q, shells)) * np.exp(-0.004 * q[:, None] ** 2)).sum(1)
        if cid.startswith("dq"):
            f = f + fq_noise * rng.standard_normal(len(f))
        g = scale * g_ref(q, f, r, s["dq"], shift=shift, zero=zero)
        if cid == "base_plain":
            g = g * (1 + perturb_plain)
        d = ev / cid
        d.mkdir(parents=True)
        np.savetxt(d / "result.pdf", np.c_[r, g])
        if cid != "base_plain":
            np.savetxt(d / "fq.dat", np.c_[q, f])


def cmd_selftest(_):
    import tempfile
    cases = {
        "clean_shift1": ({}, "PASS"),
        "clean_shift0": ({"shift": 0}, "PASS"),
        "scale_1.01": ({"scale": 1.01}, "FAIL"),
        "no_zeroing": ({"zero": False}, "FAIL"),
        "plain_differs": ({"perturb_plain": 1e-6}, "FAIL"),
        "fq_changes_with_dq": ({"fq_noise": 1e-3}, "FAIL"),
    }
    ok = True
    for name, (kw, want) in cases.items():
        with tempfile.TemporaryDirectory() as t:
            _synthetic(Path(t), **kw)
            res = analyze(Path(t))
        good = res["verdict"] == want
        ok &= good
        failed = [c["criterion"] for c in res["criteria"] if c["pass"] is False][:3]
        print(f"{'ok ' if good else 'BAD'} {name}: {res['verdict']} (expected {want}) {failed}")
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
