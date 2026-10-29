"""Waasmaier & Kirfel (1995) X-ray form factors: transcription and comparison with DISCUS 6.23.01.

Usage:
    python compare_form_factors.py --discus-src <DiffuseCode checkout at tag v.6.23.01>

Sources compared (Cu, Ga, As, Au; neutral atoms):
  dabax   ESRF DABAX f0_WaasKirf.dat (data supplied by Waasmaier & Kirfel), copied verbatim below.
          This is the reference used by the oracle.
  discus  lib_f90/element_data_mod.f90 arrays per_wa1..5, per_wb1..5, per_wc (parsed from source).
  xraydb  xraydb package 'Waasmaier' table (optional, if installed).
  ase     ase.utils.xrdebye.waasmaier, which cites "Table (1)" of the paper (optional, Cu and Au only).

The printed paper (Acta Cryst. A51, 416-431, doi:10.1107/S0108767394013292) was not reachable from
the build environment (paywalled; the IUCr pages refused the connection), so no value here was read
off the PDF. Agreement between independently maintained transcriptions is the evidence instead.
f(s) = c + sum_i a_i exp(-b_i s^2), s = sin(theta)/lambda in 1/Å.
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
Z = {"CU": 29, "GA": 31, "AS": 33, "AU": 79}

# Verbatim DABAX rows, column order "a1 a2 a3 a4 a5 c b1 b2 b3 b4 b5"
DABAX = {
    "CU": "14.014192  4.784577  5.056806  1.457971  6.932996 -3.254477   3.738280  0.003744 13.034982 72.554794  0.265666",
    "GA": "15.758946  6.841123  4.121016  2.714681  2.395246 -0.847395   3.121754  0.226057 12.482196 66.203621  0.007238",
    "AS": "17.025642  4.503441  3.715904  3.937200  6.790175 -2.984117   2.597739  0.003012 14.272119 50.437996  0.193015",
    "AU": "16.777390 19.317156 32.979683  5.595453 10.576854 -6.279078   0.122737  8.621570  1.256902 38.008820  0.000601",
}


def from_dabax(row: str) -> dict:
    v = [float(x) for x in row.split()]
    return {"a": v[0:5], "b": v[6:11], "c": v[5]}


def from_discus(src_root: Path) -> dict:
    src = (src_root / "lib_f90" / "element_data_mod.f90").read_text()

    def arr(name, kind="real"):
        m = re.search(r"::\s*%s\s*=\s*\(/(.*?)/\)" % re.escape(name), src, re.S | re.I)
        body = re.sub(r"!.*", "", m.group(1))
        if kind == "real":
            return [float(x) for x in re.findall(r"([-+]?\d+\.\d*)_PREC_DP", body)]
        return re.findall(r"'([^']{4})'", body)

    names = arr("per_name", "str")
    A = {k: arr(k) for k in ["per_wc"] + [f"per_w{t}{i}" for i in range(1, 6) for t in "ab"]}
    out = {}
    for el in Z:
        j = names.index(f"{el:<4s}")
        out[el] = {"a": [A[f"per_wa{i}"][j] for i in range(1, 6)],
                   "b": [A[f"per_wb{i}"][j] for i in range(1, 6)], "c": A["per_wc"][j],
                   "array_index_1based": j + 1}
    return out


def from_xraydb() -> dict | None:
    try:
        import xraydb
    except ImportError:
        return None
    db = xraydb.get_xraydb(); t = db.tables["Waasmaier"]
    out = {}
    for el in Z:
        row = db.query(t).filter(t.c.ion == el.capitalize()).one()
        out[el] = {"a": json.loads(row.scale), "b": json.loads(row.exponents), "c": row.offset}
    out["_version"] = xraydb.__version__
    return out


def from_ase() -> dict | None:
    try:
        import ase
        from ase.utils.xrdebye import waasmaier
    except ImportError:
        return None
    out = {"_version": ase.__version__}
    for el in Z:
        v = waasmaier.get(el.capitalize())
        if v:
            out[el] = {"a": v[0:10:2], "b": v[1:10:2], "c": v[10]}
    return out


def f_of_s(c, s):
    return c["c"] + sum(a * np.exp(-b * s * s) for a, b in zip(c["a"], c["b"]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--discus-src", type=Path, required=True)
    args = ap.parse_args()
    commit = subprocess.run(["git", "-C", str(args.discus_src), "rev-parse", "HEAD"],
                            capture_output=True, text=True).stdout.strip()
    tag = subprocess.run(["git", "-C", str(args.discus_src), "describe", "--tags", "--exact-match"],
                         capture_output=True, text=True).stdout.strip()
    srcs = {"dabax": {el: from_dabax(r) for el, r in DABAX.items()}, "discus": from_discus(args.discus_src)}
    xr, asd = from_xraydb(), from_ase()
    if xr: srcs["xraydb"] = xr
    if asd: srcs["ase"] = asd
    s = np.linspace(0.0, 2.0, 4001)
    elements, comparison = {}, {}
    for el in Z:
        elements[el] = {k: {kk: vv for kk, vv in v[el].items()} for k, v in srcs.items() if el in v}
        ref = f_of_s(srcs["dabax"][el], s)
        comparison[el] = {"Z": Z[el]}
        for k, v in srcs.items():
            if el not in v:
                continue
            fk = f_of_s(v[el], s)
            comparison[el][k] = {"f0": float(fk[0]),
                                 "max_rel_diff_vs_dabax_s_le_2": float(np.max(np.abs(fk - ref) / np.abs(ref))),
                                 "max_abs_coeff_diff_vs_dabax": float(np.max(np.abs(
                                     np.r_[v[el]["a"], v[el]["b"], v[el]["c"]] -
                                     np.r_[srcs["dabax"][el]["a"], srcs["dabax"][el]["b"], srcs["dabax"][el]["c"]])))}
    out = {"reference_for_oracle": "dabax",
           "formula": "f(s) = c + sum_i a_i exp(-b_i s^2), s = sin(theta)/lambda [1/A]",
           "paper": "D. Waasmaier & A. Kirfel, Acta Cryst. A51 (1995) 416-431, doi:10.1107/S0108767394013292",
           "paper_pdf_consulted": False,
           "sources": {"dabax": "https://ftp.esrf.fr/pub/scisoft/DabaxFiles/f0_WaasKirf.dat (verbatim rows)",
                       "discus": f"tproffen/DiffuseCode {tag or '?'} ({commit}) lib_f90/element_data_mod.f90",
                       "xraydb": f"xraydb {xr['_version']}" if xr else "not installed",
                       "ase": f"ase {asd['_version']} ase/utils/xrdebye.py" if asd else "not installed"},
           "elements": elements, "comparison": comparison}
    (HERE / "wk_coefficients.json").write_text(json.dumps(out, indent=1) + "\n")
    for el, c in comparison.items():
        print(el, {k: (round(v["f0"], 6), f"{v['max_rel_diff_vs_dabax_s_le_2']:.1e}")
                   for k, v in c.items() if k != "Z"})


if __name__ == "__main__":
    main()
