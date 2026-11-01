"""Compare the D reference probe's DISCUS outputs with discuspdf.reference_scattering (criteria: PREREG.md).

    python benchmarks/discus_6_23_01/d_reference/analyze.py <evidence dir> [--json out.json]

<evidence dir> holds work/<run>.h5 and work/*.stru as produced by run_probe.sh.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

from discuspdf.neighbors import label_sites_by_slot
from discuspdf.discus_io import read_stru
from discuspdf.discus_output import read_discus_volume
from discuspdf.reference_scattering import (ScatteringFactors, diffuse_intensity, intensity)

HERE = Path(__file__).resolve().parent
WK = json.loads((HERE.parent / "gate2_prereg" / "wk_coefficients.json").read_text())["elements"]
# DISCUS 6.23.01 coherent lengths (fm), lib_f90/element_data_mod.f90 per_bcoh (Sears 1992 values)
B_DISCUS = {"Cu": 7.718, "Au": 7.63, "Ga": 7.288, "As": 6.58}
CFINC = 0.001
CRITERIA = json.loads((HERE / "criteria.json").read_text())


def wk_f(coef, s):
    s2 = np.asarray(s, float) ** 2
    return coef["c"] + sum(a * np.exp(-b * s2) for a, b in zip(coef["a"], coef["b"]))


def factors(radiation: str, mode: str, species) -> ScatteringFactors:
    if radiation == "neut":
        return ScatteringFactors.neutron({s: B_DISCUS[s] for s in species}, source="DISCUS 6.23.01 per_bcoh")
    if mode == "physics":       # independent: DABAX coefficients, exact s
        return ScatteringFactors.xray_wk({s: WK[s.upper()]["dabax"] for s in species}, source="DABAX W-K")
    # engine emulation: DISCUS's own coefficients, f at s rounded to 0.001 with NINT (four_formtab)
    return ScatteringFactors.xray_functions(
        {s: (lambda x, c=WK[s.upper()]["discus"]: wk_f(c, np.floor(np.asarray(x) / CFINC + 0.5) * CFINC))
         for s in species}, source="DISCUS W-K table, s rounded to 0.001 (emulation)")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("evidence", type=Path)
    ap.add_argument("--json", type=Path)
    a = ap.parse_args(argv)
    runs = json.loads((HERE / "runs.json").read_text())["runs"]
    work = a.evidence / "work"
    out, failed = {}, []
    for r in runs:
        st = label_sites_by_slot(read_stru(work / r["structure"]).structure)
        vol = read_discus_volume(work / r["output"])
        hkl = np.einsum("ab,ijkb->ijka", vol.increments, np.indices(vol.values.shape).transpose(1, 2, 3, 0)) \
            + vol.corner
        fn = diffuse_intensity if r["aver"] > 0 else intensity
        row = {"bragg_meta": vol.with_bragg, "max_discus": float(vol.values.max())}
        modes = ("physics",) if r["radiation"] == "neut" else ("physics", "engine")
        for mode in modes:
            ref = fn(st, factors(r["radiation"], mode, set(st.species)), hkl=hkl, debye_waller=True)
            d = vol.values - ref
            row[mode] = {"rel_l2": float(np.linalg.norm(d) / np.linalg.norm(ref)),
                         "max_abs_over_max": float(np.abs(d).max() / np.abs(ref).max())}
        crit = CRITERIA.get(r["id"])
        if crit:
            value = row[crit["mode"]][crit["metric"]]
            row["criterion"] = dict(crit, value=value, passed=bool(value <= crit["max"]))
            if not row["criterion"]["passed"]:
                failed.append(r["id"])
        out[r["id"]] = row
        print(f"{r['id']:16s} " + "  ".join(f"{m}: L2 {row[m]['rel_l2']:.3e} max {row[m]['max_abs_over_max']:.3e}"
                                          for m in modes)
              + (f"  [{'PASS' if row['criterion']['passed'] else 'FAIL'} {crit['mode']} {crit['metric']} <= {crit['max']:g}]"
                 if crit else "  [informational]"))
    print(f"criteria failed: {len(failed)} {failed}")
    if a.json:
        a.json.write_text(json.dumps(out, indent=1) + "\n")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
