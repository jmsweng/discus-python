"""Evaluate a Gate 3 engine run against predictions.json.

Usage: python analyze.py RUN_DIR
RUN_DIR holds the registered inputs (copied unmodified) and the engine outputs
gate3_block_{total,diffuse}_{int,3dpdf}.h5 and logs. Exit: 0 PASS, 1 FAIL, 2 INVALID.
Writes RUN_DIR/gate3_results.json and RUN_DIR/GATE3_RESULTS_PREREG.md.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
BENCH = HERE.parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(BENCH)); sys.path.insert(0, str(BENCH / "gate2_prereg"))
sys.path.insert(0, str(BENCH.parents[1] / "src"))
import native3dpdf as N  # noqa: E402
import oracle as O  # noqa: E402
import regenerate_evidence as R  # noqa: E402

PRED_PATH = HERE / "predictions.json"


class Invalid(Exception):
    pass


def _crit(name, value, limit, op="<="):
    ok = bool(value <= limit) if op == "<=" else bool(value >= limit)
    return {"criterion": name, "value": float(value), "limit": float(limit), "op": op, "pass": ok}


def _h5(path: Path) -> dict:
    import h5py
    if not path.exists():
        raise Invalid(f"missing {path.name}")
    with h5py.File(path, "r") as f:
        g = f["entry/data"]
        dec = lambda k: g[k][()].decode().strip() if k in g else None
        return {"values": g["data_values"][()], "corner": g["data_corner"][()],
                "inc": g["data_increment_vector"][()], "style": dec("data_type_style"),
                "bragg": dec("data_type_with_bragg"), "cell": g["unit_cell_lengths"][()]}


def _log(run_dir: Path, kind: str) -> dict:
    p = run_dir / f"gate3_block_{kind}.log"
    if not p.exists():
        return {"ok": False, "note": "no log"}
    t = p.read_text(errors="replace")
    counts = {m.group(1): int(m.group(2)) for m in re.finditer(r"Atom typ\s*=\s*(\w+)\s*\(#\s*(\d+)\s*\)", t)}
    err = [l.strip() for l in t.splitlines() if re.search(r"\*\*\*[A-Z]+\*\*\*|\bERROR\b|\bError\b", l)]
    temp = re.search(r"Temp\. factors\s*:\s*(\w+)", t)
    ok = counts.get("CU") == 2048 and counts.get("AU") == 2048 and not err and bool(temp) and temp.group(1) != "used"
    return {"ok": ok, "counts": counts, "errors": err[:10], "temp": temp.group(1) if temp else None}


def analyze(run_dir: Path) -> dict:
    from discuspdf.discus_output import DiscusVolume, delta_pdf_fractional
    pred = json.loads(PRED_PATH.read_text())
    out = {"schema": pred["schema"], "predictions_sha256": R.sha256_file(PRED_PATH), "run_dir": str(run_dir)}
    try:
        for n, h in pred["code_sha256"].items():
            if R.sha256_file(HERE / n) != h:
                raise Invalid(f"registered code {n} changed")
        for n, h in pred["input_sha256"].items():
            if not (run_dir / n).exists() or R.sha256_file(run_dir / n) != h:
                raise Invalid(f"input {n} missing or differs from registration")
        T = pred["thresholds"]; k = pred["conversion_factor_k"]
        crit, info = [], {}
        mc = np.array(pred["native_metadata_predicted"]["data_corner"])
        mi = np.array(pred["native_metadata_predicted"]["increment"])
        for kind, P in pred["kinds"].items():
            vi = _h5(run_dir / f"gate3_block_{kind}_int.h5")
            vn = _h5(run_dir / f"gate3_block_{kind}_3dpdf.h5")
            if vi["values"].shape != (O.N_H,) * 3 or vn["values"].shape != (O.N_H,) * 3:
                raise Invalid(f"{kind}: unexpected shapes")
            if not (np.allclose(vi["corner"], -O.H_MAX) and np.allclose(vi["inc"], np.eye(3) * O.DH)):
                raise Invalid(f"{kind}: int volume grid is not the registered grid")
            if vi["bragg"] != P["with_bragg"]:
                raise Invalid(f"{kind}: int volume data_type_with_bragg={vi['bragg']}")
            h = R.content_hash(run_dir / f"gate3_block_{kind}_int.h5")
            crit.append(_crit(f"R0 {kind}: int volume content hash equals Gate 2 v4 (fraction)",
                              float(h == P["int_content_hash_gate2_v4"]), 1.0, ">="))
            dev = max(float(np.max(np.abs(vn["corner"] - mc))), float(np.max(np.abs(vn["inc"] - np.diag(mi)))))
            crit.append(_crit(f"M1 {kind}: native corner/increment equal source prediction (max abs dev)",
                              dev, T["M1_native_metadata_abs_tol"]))
            crit.append(_crit(f"M1 {kind}: native data_type_style == single_pdf",
                              float(vn["style"] == pred["native_metadata_predicted"]["data_type_style"]), 1.0, ">="))
            I = vi["values"]
            J, scalef = (N.bragg_reweight(I, np.full(3, -O.H_MAX), np.full(3, O.DH), pred["design"]["ncells"])
                         if kind == "total" else (I, None))
            emu = N.native_3dpdf(J)
            crit.append(_crit(f"E1 {kind}: native vs emulation of same-run int (relL2)",
                              O.rel_l2(vn["values"], emu), T["E1_native_vs_emulation_of_same_run_int_relL2_max"]))
            crit.append(_crit(f"E2 {kind}: native matches registered fingerprint (fraction)",
                              float(O.fingerprint_matches(vn["values"], P["native_fingerprint"],
                                                          rtol=T["E2_native_fingerprint_rtol"])), 1.0, ">="))
            am = [int(v) for v in np.unravel_index(np.argmax(vn["values"]), vn["values"].shape)]
            crit.append(_crit(f"P1 {kind}: native argmax at predicted index {P['native_argmax_index']} (fraction)",
                              float(am == P["native_argmax_index"]), 1.0, ">="))
            info[kind] = {"int_content_hash": h, "four_weight_scalef_used": scalef, "native_argmax_index": am,
                          "native_label_at_argmax": (vn["corner"] + np.diag(vn["inc"]) * np.array(am)).tolist()}
            if kind == "diffuse":
                vol = DiscusVolume(I, np.full(3, -O.H_MAX), np.eye(3) * O.DH, vi["cell"], np.full(3, 90.0),
                                   "bragg_subtracted")
                r = delta_pdf_fractional(vol, imaginary_tolerance=1.0, raise_on_imaginary=False)
                assert np.isclose(r.u_axes[0][64], 0.0) and np.isclose(r.u_axes[0][65] - r.u_axes[0][64], 0.125)
                crit.append(_crit("C1 diffuse: package delta_pdf_fractional == k * native on true axes (relL2)",
                                  O.rel_l2(k * vn["values"], r.real_volume), T["C1_pipeline_vs_converted_native_relL2_max"]))
            lg = _log(run_dir, kind); info[kind]["log"] = lg
            crit.append(_crit(f"L {kind}: log atom counts, no errors, Temp. factors not 'used'", float(lg["ok"]), 1.0, ">="))
        out.update(criteria=crit, info=info, n_failed=sum(not c["pass"] for c in crit))
        out["status"] = "PASS" if out["n_failed"] == 0 else "FAIL"
    except Invalid as e:
        out.update(status="INVALID", problems=[str(e)])
    return out


def write_report(res: dict, run_dir: Path):
    (run_dir / "gate3_results.json").write_text(json.dumps(res, indent=1, default=str) + "\n")
    L = [f"# Gate 3 pre-registered experiment: {res['status']}", "",
         f"predictions.json sha256 {res['predictions_sha256']}", ""]
    if res["status"] == "INVALID":
        L += ["Problems:"] + [f"- {p}" for p in res["problems"]]
    else:
        L += [f"{res['n_failed']}/{len(res['criteria'])} criteria failed.", "",
              "| criterion | value | limit | pass |", "|---|---|---|---|"]
        L += [f"| {c['criterion']} | {c['value']:.3e} | {c['op']} {c['limit']:.3e} | {'PASS' if c['pass'] else 'FAIL'} |"
              for c in res["criteria"]]
        L += ["", "```json", json.dumps(res["info"], indent=1, default=str), "```"]
    (run_dir / "GATE3_RESULTS_PREREG.md").write_text("\n".join(L) + "\n")


def main():
    run_dir = Path(sys.argv[1]).resolve()
    res = analyze(run_dir)
    write_report(res, run_dir)
    print(f"Gate 3 pre-registered experiment: {res['status']}"
          + (f" ({res['n_failed']}/{len(res['criteria'])} criteria failed)" if "criteria" in res else ""))
    sys.exit({"PASS": 0, "FAIL": 1}.get(res["status"], 2))


if __name__ == "__main__":
    main()
