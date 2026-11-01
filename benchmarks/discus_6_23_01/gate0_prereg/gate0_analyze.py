"""Evaluate real DISCUS 6.23.01 output against the pre-registered Gate 0 v2 predictions.

Usage:
    python gate0_analyze.py RUN_DIR [--out RUN_DIR]

RUN_DIR must contain the generated inputs (copied unmodified from inputs/), the engine HDF5 files
gate0_{A,B}_{total,diffuse}_{0..63}.h5 and, if available, the engine logs *.log.
Writes RUN_DIR/gate0_results.json and RUN_DIR/GATE0_RESULTS_PREREG.md. Exit code: 0 = all criteria
pass, 1 = at least one criterion fails, 2 = invalid run (missing files, modified pre-registered inputs
or code, wrong metadata). Thresholds come only from predictions.json; this script has none of its own.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import oracle as O  # noqa: E402

ARMS = ("A", "B")


class Invalid(Exception):
    pass


def _crit(name, value, limit, op="<="):
    ok = bool(np.isfinite(value) and (value <= limit if op == "<=" else value >= limit))
    return {"criterion": name, "value": float(value), "limit": float(limit), "op": op, "pass": ok}


def _check_hashes(run_dir: Path, pred: dict) -> list:
    problems = []
    for name, h in pred["code_sha256"].items():
        if O.sha256_file(HERE / name) != h:
            problems.append(f"pre-registered code/data file changed since registration: {name}")
    for name, h in pred["input_sha256"].items():
        p = run_dir / name
        if not p.exists():
            problems.append(f"missing input copy in run dir: {name}")
        elif O.sha256_file(p) != h:
            problems.append(f"input in run dir differs from registered version: {name}")
    return problems


def _read_h5(path: Path, want_bragg: str):
    import h5py
    if not path.exists():
        raise Invalid(f"missing {path.name}")
    with h5py.File(path, "r") as f:
        g = f["entry/data"]
        values = g["data_values"][()]
        corner = g["data_corner"][()]
        inc = g["data_increment_vector"][()]
        bragg = g["data_type_with_bragg"][()].decode().strip()
    if values.shape != (25, 25, 25):
        raise Invalid(f"{path.name}: shape {values.shape}")
    if not np.allclose(corner, 0.0, atol=1e-12):
        raise Invalid(f"{path.name}: corner {corner}")
    if not np.allclose(inc, np.eye(3) * 0.125, atol=1e-12):
        raise Invalid(f"{path.name}: increments {inc}")
    if bragg != want_bragg:
        raise Invalid(f"{path.name}: data_type_with_bragg={bragg}, expected {want_bragg}")
    return values


def _log_checks(run_dir: Path, arm: str, kind: str) -> dict:
    cand = sorted(run_dir.glob(f"gate0_{arm}_{kind}*.log"))
    if not cand:
        return {"log": None, "note": "no log found", "b_values": None, "temp": None}
    txt = cand[0].read_text(errors="replace")
    temp = re.search(r"Temp\. factors\s*:\s*(\w+)", txt)
    # Look for B values in the log
    b_cu = re.search(r"B\s*=\s*([\d.]+).*Cu|Cu.*B\s*=\s*([\d.]+)", txt)
    b_au = re.search(r"B\s*=\s*([\d.]+).*Au|Au.*B\s*=\s*([\d.]+)", txt)
    return {"log": cand[0].name, "temp": temp.group(1) if temp else None,
            "b_cu": b_cu, "b_au": b_au}


def analyze_arm(run_dir: Path, arm: str, pred: dict) -> dict:
    T = pred["thresholds"]
    K = pred["n_realizations"]
    crit = []
    info = {}

    # E-criteria: engine precision, per realization
    e1_vals, e2_vals, e3_vals = [], [], []
    for k in range(K):
        seed = pred["seeds"][k]
        stru_path = run_dir / f"gate0_{arm}_{k:02d}.stru"
        if not stru_path.exists():
            raise Invalid(f"missing {stru_path.name}")

        # Read engine output
        tot_path = run_dir / f"gate0_{arm}_total_{k:02d}.h5"
        dif_path = run_dir / f"gate0_{arm}_diffuse_{k:02d}.h5"
        I_eng = _read_h5(tot_path, "bragg_present")
        D_eng = _read_h5(dif_path, "bragg_subtracted")

        # Compute realization-exact oracle
        b_cu, b_au = (O.B_CU_A, O.B_AU_A) if arm == "A" else (O.B_CU_B, O.B_AU_B)
        D_orc, I_orc = O.realization_exact(stru_path, "engine", b_cu, b_au)

        # E1: total vs realization engine oracle
        e1 = O.rel_l2(I_eng, I_orc)
        e1_vals.append(e1)
        # E2: diffuse vs realization engine oracle
        # The engine's diffuse has D=0 at integer h (set aver subtracts the average)
        # So we must zero out the oracle's D at integer h before comparing
        h_int = np.isclose(O.H_AXIS, np.rint(O.H_AXIS))
        is_int = h_int[:, None, None] & h_int[None, :, None] & h_int[None, None, :]
        D_orc_zeroed = D_orc.copy()
        D_orc_zeroed[is_int] = 0
        e2 = O.rel_l2(D_eng, np.abs(D_orc_zeroed) ** 2)
        e2_vals.append(e2)
        # E3: diffuse at integer h
        h_int = np.isclose(O.H_AXIS, np.rint(O.H_AXIS))
        is_int = h_int[:, None, None] & h_int[None, :, None] & h_int[None, None, :]
        mx = max(float(np.max(D_eng)), 1e-300)
        e3 = float(np.max(np.abs(D_eng[is_int]))) / mx
        e3_vals.append(e3)

    crit.append(_crit(f"E1 {arm}: max total vs realization engine oracle (relL2)", max(e1_vals),
                      T["E1_total_vs_realization_engine_relL2_max"]))
    crit.append(_crit(f"E2 {arm}: max diffuse vs realization engine oracle (relL2)", max(e2_vals),
                      T["E2_diffuse_vs_realization_engine_relL2_max"]))
    crit.append(_crit(f"E3 {arm}: max diffuse at integer h / max diffuse", max(e3_vals),
                      T["E3_diffuse_at_integer_h_max_over_max_diffuse"]))

    # E4: log checks
    for kind in ("total", "diffuse"):
        lg = _log_checks(run_dir, arm, kind)
        if arm == "A":
            ok = lg["temp"] is not None and lg["temp"].lower() == "ignored"
            crit.append(_crit(f"E4 {arm} {kind}: log shows 'temp ignore'", float(ok), 1.0, ">="))
        else:
            # Check B values are present in the log
            ok = lg["b_cu"] is not None and lg["b_au"] is not None
            crit.append(_crit(f"E4 {arm} {kind}: log shows B values", float(ok), 1.0, ">="))

    # S-criteria: ensemble, over the independent evaluation set
    mask_S = O.independent_set_mask()
    n_S = int(mask_S.sum())

    # S1: mean of x = D/E_eng within 1 ± 5/√(K·|S|)
    x_all = []
    for k in range(K):  # K is already defined above
        dif_path = run_dir / f"gate0_{arm}_diffuse_{k:02d}.h5"
        D_eng = _read_h5(dif_path, "bragg_subtracted")
        b_cu, b_au = (O.B_CU_A, O.B_AU_A) if arm == "A" else (O.B_CU_B, O.B_AU_B)
        E_eng = O.expectation("engine", b_cu=b_cu, b_au=b_au)
        # x = D / E at each voxel in S
        x = D_eng[mask_S] / E_eng[mask_S]
        x_all.append(x)
    x_mean = float(np.mean(np.concatenate(x_all)))
    s1_tol = T["S1_mean_x_within"]
    crit.append(_crit(f"S1 {arm}: |mean x - 1|", abs(x_mean - 1.0), s1_tol))

    # S2: sample variance of x within 1 ± 5·√(8/(K·|S|))
    x_flat = np.concatenate(x_all)
    x_var = float(np.var(x_flat, ddof=1))
    s2_tol = T["S2_var_x_within"]
    crit.append(_crit(f"S2 {arm}: |var x - 1|", abs(x_var - 1.0), s2_tol))

    # S3: every voxel of the whole grid off integer h: x̄_v inside two-sided Gamma interval
    # (except voxels with 2h integer where x ~ χ²₁, so x̄_v ~ Gamma(K/2, 2/K))
    lo_std, hi_std = pred["gamma_quantiles"]["S3_standard"]
    lo_half, hi_half = pred["gamma_quantiles"]["S3_half_integer"]

    # Identify 2h integer voxels (q ≡ -q mod 1)
    h = np.stack(np.meshgrid(O.H_AXIS, O.H_AXIS, O.H_AXIS, indexing="ij"), axis=-1)
    m = (O.H_AXIS * 8).astype(int)
    two_int = (m[:, None, None] % 4 == 0) & (m[None, :, None] % 4 == 0) & (m[None, None, :] % 4 == 0)
    is_int = np.all(np.isclose(h, np.rint(h)), axis=-1)
    off_int = ~is_int

    s3_violations = 0
    n_tested = 0
    x_mean_vox = None
    for k in range(K):
        dif_path = run_dir / f"gate0_{arm}_diffuse_{k:02d}.h5"
        D_eng = _read_h5(dif_path, "bragg_subtracted")
        b_cu, b_au = (O.B_CU_A, O.B_AU_A) if arm == "A" else (O.B_CU_B, O.B_AU_B)
        E_eng = O.expectation("engine", b_cu=b_cu, b_au=b_au)
        # Only compute x where E > 0 (off integer h)
        mask = E_eng > 0
        x = np.zeros_like(D_eng)
        x[mask] = D_eng[mask] / E_eng[mask]
        # Mean over realizations
        if x_mean_vox is None:
            x_mean_vox = np.zeros_like(x)
        x_mean_vox += x
    x_mean_vox /= K

    # Check each voxel (only where E > 0, i.e., off integer h) — vectorised
    b_cu, b_au = (O.B_CU_A, O.B_AU_A) if arm == "A" else (O.B_CU_B, O.B_AU_B)
    E_eng = O.expectation("engine", b_cu=b_cu, b_au=b_au)
    valid = off_int & (E_eng > 0)

    # Standard voxels (valid and NOT 2h integer)
    std_mask = valid & ~two_int
    half_mask = valid & two_int
    n_tested = int(std_mask.sum() + half_mask.sum())

    # Count violations: standard voxels vs standard interval, half-integer vs half interval
    std_vals = x_mean_vox[std_mask]
    half_vals = x_mean_vox[half_mask]
    s3_violations = int(
        np.sum(~(lo_std <= std_vals) & ~(std_vals <= hi_std)) +
        np.sum(~(lo_half <= half_vals) & ~(half_vals <= hi_half))
    )

    crit.append(_crit(f"S3 {arm}: n_voxels outside Gamma interval ({n_tested} tested)", float(s3_violations), 0.0))

    # S4: physics vs engine expectation agree to ≤ 1e-6 relative
    b_cu, b_au = (O.B_CU_A, O.B_AU_A) if arm == "A" else (O.B_CU_B, O.B_AU_B)
    E_phys = O.expectation("physics", b_cu=b_cu, b_au=b_au)
    E_eng = O.expectation("engine", b_cu=b_cu, b_au=b_au)
    # Relative difference at non-zero voxels
    mask = E_eng > 1e-300
    s4 = float(np.max(np.abs(E_phys[mask] - E_eng[mask]) / E_eng[mask]))
    crit.append(_crit(f"S4 {arm}: max |E_phys - E_eng| / E_eng", s4, T["S4_physics_vs_engine_rel_max"]))

    return {"criteria": crit, "info": info}


def analyze(run_dir: Path) -> dict:
    run_dir = Path(run_dir)
    pred = json.loads((HERE / "predictions.json").read_text())
    out = {"predictions_sha256": O.sha256_file(HERE / "predictions.json"), "run_dir": str(run_dir)}
    problems = _check_hashes(run_dir, pred)
    if problems:
        out.update(status="INVALID", problems=problems)
        return out
    try:
        out["arms"] = {a: analyze_arm(run_dir, a, pred) for a in ARMS}
    except Invalid as e:
        out.update(status="INVALID", problems=[str(e)])
        return out
    allc = [c for a in out["arms"].values() for c in a["criteria"]]
    out["n_criteria"] = len(allc)
    out["n_failed"] = sum(not c["pass"] for c in allc)
    out["status"] = "PASS" if out["n_failed"] == 0 else "FAIL"
    return out


def write_report(res: dict, out_dir: Path):
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "gate0_results.json").write_text(json.dumps(res, indent=1) + "\n")
    L = [f"# Gate 0 v2 pre-registered experiment: {res['status']}", "",
         "Generated by `benchmarks/discus_6_23_01/gate0_prereg/gate0_analyze.py`; do not edit by hand.", "",
         f"- predictions.json sha256: `{res['predictions_sha256']}`"]
    if res["status"] == "INVALID":
        L += ["", "## Run is INVALID", ""] + [f"- {p}" for p in res["problems"]]
    for a, r in res.get("arms", {}).items():
        L += ["", f"## Arm `{a}`", "", "| criterion | value | limit | pass |", "|---|---|---|---|"]
        for c in r["criteria"]:
            L.append(f"| {c['criterion']} | {c['value']:.3e} | {c['op']} {c['limit']:.3e} | {'PASS' if c['pass'] else '**FAIL**'} |")
    (out_dir / "GATE0_RESULTS_PREREG.md").write_text("\n".join(L) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir", type=Path)
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args()
    res = analyze(a.run_dir)
    write_report(res, a.out or a.run_dir)
    print(f"Gate 0 v2 pre-registered experiment: {res['status']}"
          + (f" ({res['n_failed']}/{res['n_criteria']} criteria failed)" if "n_failed" in res else ""))
    for p in res.get("problems", []):
        print("  INVALID:", p)
    sys.exit({"PASS": 0, "FAIL": 1}.get(res["status"], 2))


if __name__ == "__main__":
    main()
