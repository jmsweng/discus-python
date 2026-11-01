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


def _round_s(s):
    return np.floor(np.asarray(s, float) / CFINC + 0.5) * CFINC     # Fortran NINT for s >= 0


def factors(radiation: str, mode: str, b_iso: dict) -> tuple[ScatteringFactors, bool]:
    """(factors, debye_waller flag) for one reference mode.

    physics: X-ray DABAX W-K at the exact s, neutron DISCUS lengths; DW exp(-B s^2) at the exact s.
    engine:  as DISCUS 6.23.01 four_formtab for isotropic B: the whole factor exp(-B s_r^2) * f(s_r) is
             tabulated at s_r = NINT(s / 0.001) * 0.001 (fourier_sup.f90, cfact; four_dbwtab runs only for
             anisotropic ADPs). X-ray uses DISCUS's own W-K coefficients."""
    species = sorted(b_iso)
    if mode == "physics":
        if radiation == "neut":
            return ScatteringFactors.neutron({s: B_DISCUS[s] for s in species}, "DISCUS 6.23.01 per_bcoh"), True
        return ScatteringFactors.xray_wk({s: WK[s.upper()]["dabax"] for s in species}, "DABAX W-K"), True
    funcs = {}
    for s in species:
        if radiation == "neut":
            funcs[s] = lambda x, b=B_DISCUS[s], B=b_iso[s]: b * np.exp(-B * _round_s(x) ** 2)
        else:
            funcs[s] = lambda x, c=WK[s.upper()]["discus"], B=b_iso[s]: (
                wk_f(c, _round_s(x)) * np.exp(-B * _round_s(x) ** 2))
    # container label 'xray' only means "s-dependent factor"; the source string says what it is
    return ScatteringFactors.xray_functions(funcs, f"DISCUS {radiation} table emulation: f(s_r)·exp(-B s_r^2)"), False


def species_b(st) -> dict:
    out = {}
    for sp in set(st.species):
        b = np.unique(st.b_iso[np.asarray(st.species) == sp])
        if len(b) != 1:
            raise ValueError(f"{sp}: B differs between atoms; DISCUS has one B per atom type")
        out[sp] = float(b[0])
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("evidence", type=Path)
    ap.add_argument("--json", type=Path)
    a = ap.parse_args(argv)
    spec = json.loads((HERE / "runs.json").read_text())
    runs, GRIDS = spec["runs"], spec["grids"]
    work = a.evidence / "work"
    out, failed = {}, []
    # integrity: every run exited 0, inputs are the committed ones, outputs exist
    codes = dict(l.split(" exit=") for l in (a.evidence / "exit_codes.txt").read_text().split("\n") if " exit=" in l)
    bad = {k: v for k, v in codes.items() if v.strip() != "0"}
    if bad or set(codes) != {r["id"] for r in runs}:
        print(f"engine exit codes not all 0 / runs missing: {codes}")
        failed.append("exit_codes")
    import hashlib
    for line in (a.evidence / "input_sha256.txt").read_text().splitlines():
        h, name = line.split()
        if hashlib.sha256((HERE / "inputs" / name).read_bytes().replace(b"\r\n", b"\n")).hexdigest() != h:
            print(f"input {name} differs from the committed file")
            failed.append(f"input:{name}")
    for r in runs:
        st = label_sites_by_slot(read_stru(work / r["structure"]).structure)
        vol = read_discus_volume(work / r["output"])
        hkl = np.einsum("ab,ijkb->ijka", vol.increments, np.indices(vol.values.shape).transpose(1, 2, 3, 0)) \
            + vol.corner
        fn = diffuse_intensity if r["aver"] > 0 else intensity
        row = {"bragg_meta": vol.with_bragg, "max_discus": float(vol.values.max())}
        want_bragg = "bragg_subtracted" if r["aver"] > 0 else "bragg_present"
        if vol.with_bragg != want_bragg or list(vol.values.shape) != list(GRIDS[r["grid"]]["n"]):
            print(f"{r['id']}: metadata {vol.with_bragg} / shape {vol.values.shape} not as requested")
            failed.append(f"meta:{r['id']}")
        modes = ("physics", "engine")
        bmap = species_b(st)
        for mode in modes:
            fac, dw = factors(r["radiation"], mode, bmap)
            ref = fn(st, fac, hkl=hkl, debye_waller=dw)
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
