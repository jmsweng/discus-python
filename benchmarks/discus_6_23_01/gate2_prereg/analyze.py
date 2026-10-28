"""Evaluate real DISCUS 6.23.01 output against the pre-registered Gate 2 predictions.

Usage:
    python analyze.py RUN_DIR [--out RUN_DIR]

RUN_DIR must contain the generated inputs (copied unmodified from inputs/), the four engine HDF5 files
gate2_{block,random}_{total,diffuse}.h5 and, if available, the engine logs *.log.
Writes RUN_DIR/gate2_results.json and RUN_DIR/GATE2_RESULTS_PREREG.md. Exit code: 0 = all criteria
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
sys.path.insert(0, str(HERE.parents[2] / "src"))
import oracle as O  # noqa: E402

MODELS = ("block", "random")


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


def _log_checks(run_dir: Path, model: str, kind: str) -> dict:
    cand = sorted(run_dir.glob(f"gate2_{model}_{kind}*.log"))
    if not cand:
        return {"log": None, "note": "no log found (gating log checks skipped, recorded as missing)",
                "atom_counts_ok": False, "errors_found": None}
    txt = cand[0].read_text(errors="replace")
    counts = {m.group(1): int(m.group(2)) for m in
              re.finditer(r"Atom typ\s*=\s*(\w+)\s*\(#\s*(\d+)\s*\)", txt)}
    used = re.findall(r"Used\s+([\d.]+)\s*% of the crystal", txt)
    err = [l.strip() for l in txt.splitlines() if re.search(r"\*\*\*|ERROR|Error", l)]
    return {"log": cand[0].name, "sha256": O.sha256_file(cand[0]), "atom_type_counts": counts,
            "atom_counts_ok": counts.get("CU") == 2048 and counts.get("AU") == 2048,
            "used_percent_lines": used, "form_factor_line": "Waasmeier & Kirfel" in txt,
            "error_lines": err[:20], "errors_found": bool(err)}


def _read(path: Path, want_bragg: str):
    from discuspdf.discus_output import read_discus_volume
    vol = read_discus_volume(path)
    problems = []
    if vol.values.shape != (O.N_H,) * 3:
        problems.append(f"{path.name}: shape {vol.values.shape}")
    if not np.allclose(vol.corner, -O.H_MAX, atol=1e-12):
        problems.append(f"{path.name}: corner {vol.corner}")
    if not np.allclose(vol.increments, np.eye(3) * O.DH, atol=1e-12):
        problems.append(f"{path.name}: increments {vol.increments.tolist()}")
    if not (np.allclose(vol.cell_lengths, O.A_CELL) and np.allclose(vol.cell_angles, 90.0)):
        problems.append(f"{path.name}: cell {vol.cell_lengths} {vol.cell_angles}")
    if vol.with_bragg != want_bragg:
        problems.append(f"{path.name}: data_type_with_bragg={vol.with_bragg}, expected {want_bragg}")
    return vol, problems


def analyze_model(run_dir: Path, model: str, pred: dict) -> dict:
    from discuspdf.discus_output import delta_pdf_fractional
    T = pred["thresholds"]; P = pred["models"][model]
    tot, p1 = _read(run_dir / f"gate2_{model}_total.h5", "bragg_present")
    dif, p2 = _read(run_dir / f"gate2_{model}_diffuse.h5", "bragg_subtracted")
    if p1 or p2:
        raise Invalid("; ".join(p1 + p2))
    sig = O.parse_stru(run_dir / f"gate2_{model}.stru")
    if not np.array_equal(sig, O.build_sigma(model)):
        raise Invalid(f"{model}: structure file does not match the registered realization")
    orc, bitwise = {}, {}
    for mode in ("physics", "engine"):
        It, Id, g2 = O.forward(sig, mode)
        # Validity: the recomputed oracle must equal the registered one numerically (1e-9). Bitwise
        # identity is only reported, since libm/FFT results may differ by an ulp across machines.
        if not (O.fingerprint_matches(It, P[mode]["I_total_fingerprint"]) and
                O.fingerprint_matches(Id, P[mode]["I_diffuse_fingerprint"])):
            raise Invalid(f"{model}/{mode}: recomputed oracle differs from the registered fingerprint")
        bitwise[mode] = (O.sha256_array(It) == P[mode]["I_total_sha256"] and
                         O.sha256_array(Id) == P[mode]["I_diffuse_sha256"])
        orc[mode] = (It, Id, g2)
    It_p, Id_p, g2_p = orc["physics"]; It_e, Id_e, g2_e = orc["engine"]
    eps = P["eps_tab"]; floor = T["floor_for_factor_thresholds"]
    It_eng, Id_eng = tot.values, dif.values
    crit = []
    # C0 forward model of total intensity
    crit.append(_crit("C0 total vs engine-precision oracle (relL2)", O.rel_l2(It_eng, It_e),
                      T["C0_total_vs_engine_precision_relL2_max"]))
    crit.append(_crit("C0 total vs physics oracle (relL2)", O.rel_l2(It_eng, It_p),
                      max(T["C0_total_vs_physics_relL2_factor_of_eps_tab"] * eps["I_total_relL2"], floor)))
    # C1 averaging (set aver,100) on the same frozen list
    is_int = np.isclose(O.H_AXIS, np.rint(O.H_AXIS))
    mask_int = np.zeros(It_eng.shape, bool); mask_int[np.ix_(is_int, is_int, is_int)] = True
    mx = max(float(np.max(Id_eng)), 1e-300)
    crit.append(_crit("C1 diffuse at integer H / max diffuse", float(np.max(np.abs(Id_eng[mask_int]))) / mx,
                      T["C1_diffuse_at_integer_H_max_over_max_diffuse"]))
    crit.append(_crit("C1 (diffuse - total) off integer H (relL2)",
                      O.rel_l2(Id_eng[~mask_int], It_eng[~mask_int]),
                      T["C1_diffuse_minus_total_off_integer_relL2_max"]))
    # A exact lattice correlations from one Brillouin zone, H in [-1/2, 1/2)
    i0 = int(np.argmin(np.abs(O.H_AXIS + 0.5))); zone = slice(i0, i0 + O.N_CELL)
    kz = np.mod(np.rint(O.H_AXIS[zone] * O.N_CELL).astype(int), O.N_CELL)
    C_or = O.pair_sums_direct(sig).astype(float); Nn = sig.size
    zone_dev = {}
    for mode, g2 in (("physics", g2_p), ("engine", g2_e)):
        S = np.empty((O.N_CELL,) * 3)
        S[np.ix_(kz, kz, kz)] = (Id_eng / g2)[zone, zone, zone]
        C_eng = np.real(np.fft.fftn(S)) / Nn
        err = float(np.max(np.abs(C_eng - C_or)) / Nn)
        key = "physics" if mode == "physics" else "engine_precision"
        crit.append(_crit(f"A correlation recovery, {mode} g (max|dC|/N)", err,
                          T[f"A_corr_recovery_{key}_max_abs_err_over_N"]))
        Sz = (Id_eng / g2).reshape(8, 16, 8, 16, 8, 16)
        zone_dev[mode] = float(np.max(np.abs(Sz - Sz.mean(axis=(0, 2, 4), keepdims=True))) / max(Sz.max(), 1e-300))
    # B pipeline 3D-dPDF of the engine diffuse volume
    r = delta_pdf_fractional(dif, imaginary_tolerance=1.0, raise_on_imaginary=False)
    from discuspdf.discus_output import DiscusVolume
    rp = {m: delta_pdf_fractional(DiscusVolume(orc[m][1], dif.corner, dif.increments, dif.cell_lengths,
                                               dif.cell_angles, "bragg_subtracted"),
                                  imaginary_tolerance=1.0, raise_on_imaginary=False) for m in orc}
    crit.append(_crit("B1 dPDF vs physics-oracle dPDF (relL2)", O.rel_l2(r.real_volume, rp["physics"].real_volume),
                      max(T["B1_dpdf_vs_oracle_physics_relL2_factor_of_eps_tab"] * eps["dpdf_relL2"], floor)))
    crit.append(_crit("B1 dPDF vs engine-precision-oracle dPDF (relL2)",
                      O.rel_l2(r.real_volume, rp["engine"].real_volume),
                      T["B1_dpdf_vs_oracle_engine_precision_relL2_max"]))
    c0 = O.N_H // 2; st = O.N_H // O.N_CELL
    rho = {",".join(map(str, t)): float(r.real_volume[c0 + t[0] * st, c0 + t[1] * st, c0 + t[2] * st])
           for t in O.TAU_POINTS}
    rho0 = rho["0,0,0"]
    table = []
    for t in O.TAU_POINTS:
        k = ",".join(map(str, t)); pr = P["physics"]["rho_pred_ratio"][k]
        ratio = rho[k] / rho0 if rho0 != 0 else float("nan")
        row = {"tau_cells": k, "tau_A": [v * O.A_CELL for v in t], "engine_rho": rho[k],
               "engine_ratio": ratio, "predicted_ratio": pr,
               "realized_C_over_N": P["C_over_N_realized"][k],
               "ratio_ok": bool(abs(ratio - pr) <= T["B2_ratio_abs_tolerance"])}
        if abs(pr) >= T["B2_sign_required_if_abs_predicted_ratio_ge"]:
            row["sign_required"] = True
            row["sign_ok"] = bool(np.sign(ratio) == np.sign(pr))
        table.append(row)
    worst = max(abs(row["engine_ratio"] - row["predicted_ratio"]) for row in table)
    crit.append(_crit("B2 max |ratio - predicted| at predeclared tau", worst, T["B2_ratio_abs_tolerance"]))
    nsign = [row for row in table if row.get("sign_required")]
    crit.append(_crit("B2 sign agreement where required (fraction)",
                      sum(row["sign_ok"] for row in nsign) / max(len(nsign), 1), 1.0, ">="))
    crit.append(_crit("B3 imag/real of transform", r.imaginary_ratio, T["B3_imag_over_real_max"]))
    crit.append(_crit("B3 Friedel relative error of engine diffuse", r.friedel_relative_error,
                      T["B3_friedel_rel_error_max"]))
    if model == "random":
        nn = [abs(row["engine_ratio"]) for row in table if sum(abs(v) for v in map(int, row["tau_cells"].split(","))) == 1]
        crit.append(_crit("N control max |rho(NN)/rho(0)|", max(nn), T["N_control_abs_nn_ratio_max"]))
    logs = {k: _log_checks(run_dir, model, k) for k in ("total", "diffuse")}
    for k, lg in logs.items():
        crit.append(_crit(f"L {k} log: atom counts CU=2048, AU=2048 and no error lines",
                          float(bool(lg["atom_counts_ok"]) and lg["errors_found"] is False), 1.0, ">="))
    return {"criteria": crit, "tau_table": table, "zone_consistency_max_rel_dev": zone_dev,
            "oracle_bitwise_identical_to_registration": bitwise,
            "hdf5_sha256": {"total": tot.provenance["sha256"], "diffuse": dif.provenance["sha256"]},
            "max_I_total_engine": float(It_eng.max()), "max_I_diffuse_engine": float(Id_eng.max()),
            "logs": logs}


def analyze(run_dir: Path) -> dict:
    run_dir = Path(run_dir)
    pred = json.loads((HERE / "predictions.json").read_text())
    out = {"predictions_sha256": O.sha256_file(HERE / "predictions.json"), "run_dir": str(run_dir)}
    problems = _check_hashes(run_dir, pred)
    if problems:
        out.update(status="INVALID", problems=problems)
        return out
    try:
        out["models"] = {m: analyze_model(run_dir, m, pred) for m in MODELS}
    except Invalid as e:
        out.update(status="INVALID", problems=[str(e)])
        return out
    allc = [c for m in out["models"].values() for c in m["criteria"]]
    out["n_criteria"] = len(allc)
    out["n_failed"] = sum(not c["pass"] for c in allc)
    out["status"] = "PASS" if out["n_failed"] == 0 else "FAIL"
    return out


def write_report(res: dict, out_dir: Path):
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "gate2_results.json").write_text(json.dumps(res, indent=1) + "\n")
    L = [f"# Gate 2 pre-registered experiment: {res['status']}", "",
         "Generated by `benchmarks/discus_6_23_01/gate2_prereg/analyze.py`; do not edit by hand.", "",
         f"- predictions.json sha256: `{res['predictions_sha256']}`"]
    if res["status"] == "INVALID":
        L += ["", "## Run is INVALID", ""] + [f"- {p}" for p in res["problems"]]
    for m, r in res.get("models", {}).items():
        L += ["", f"## Model `{m}`", "", f"- HDF5 sha256: total `{r['hdf5_sha256']['total']}`, "
              f"diffuse `{r['hdf5_sha256']['diffuse']}`", "", "| criterion | value | limit | pass |", "|---|---|---|---|"]
        for c in r["criteria"]:
            L.append(f"| {c['criterion']} | {c['value']:.3e} | {c['op']} {c['limit']:.3e} | {'PASS' if c['pass'] else '**FAIL**'} |")
        L += ["", "| tau (cells) | engine rho/rho0 | predicted | realized C/N | sign req. |", "|---|---|---|---|---|"]
        for t in r["tau_table"]:
            sg = ("ok" if t["sign_ok"] else "**WRONG**") if t.get("sign_required") else "-"
            L.append(f"| {t['tau_cells']} | {t['engine_ratio']:+.4f} | {t['predicted_ratio']:+.4f} | "
                     f"{t['realized_C_over_N']:+.4f} | {sg} |")
    (out_dir / "GATE2_RESULTS_PREREG.md").write_text("\n".join(L) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir", type=Path)
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args()
    res = analyze(a.run_dir)
    write_report(res, a.out or a.run_dir)
    print(f"Gate 2 pre-registered experiment: {res['status']}"
          + (f" ({res['n_failed']}/{res['n_criteria']} criteria failed)" if "n_failed" in res else ""))
    for p in res.get("problems", []):
        print("  INVALID:", p)
    sys.exit({"PASS": 0, "FAIL": 1}.get(res["status"], 2))


if __name__ == "__main__":
    main()
