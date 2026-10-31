"""Reviewer analysis of the Package B probe evidence (read-only). Prints the numbers quoted in RESULT.md.

    python benchmarks/discus_6_23_01/b_probes/analyze_b_probes.py [evidence dir]
"""
from __future__ import annotations

import hashlib
import json
import re
import sys
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
DEFAULT = HERE / "evidence" / "20261031T182936Z_LAPTOP-NGU3BT0N"


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    E = Path(argv[0]) if argv else DEFAULT
    out: dict = {}

    def ld(run, v="int"):
        with h5py.File(E / run / f"{run}_{v}.h5", "r") as f:
            g = f["/entry/data"]
            return np.asarray(g["data_values"][()], float), {k: g[k][()] for k in g if k != "data_values"}

    def log(run):
        return re.sub(r"\x1b\[[0-9;]*m", "", (E / run / "run.log").read_text(errors="replace"))

    def rel(a, b):
        return float(np.linalg.norm(a - b) / max(np.linalg.norm(b), 1e-300))

    def digest(a):
        return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()[:16]

    def errors(run):
        return [l.strip() for l in log(run).splitlines() if re.search(r"\*\*\*[A-Z ]+\*\*\*.*\*\*\*|[Ee]rror \(ier", l)]

    # P1 planes and lines
    p1 = {}
    for r in ("p1a_hk0", "p1b_h00", "p1d_tloff", "p1e_h0l"):
        v, m = ld(r)
        p1[r] = {"shape": list(v.shape), "corner": m["data_corner"].tolist(),
                 "increment_columns": m["data_increment_vector"].T.tolist()}
    p1["p1a_equals_p1d_values"] = bool(np.array_equal(ld("p1a_hk0")[0], ld("p1d_tloff")[0]))
    p1["p1c_hhl_errors"] = errors("p1c_hhl")
    out["P1"] = p1

    # P2 set aver
    I0, m = ld("p2_aver0")
    F0 = ld("p2_aver0", "real")[0] + 1j * ld("p2_aver0", "imag")[0]
    I1 = ld("p2_aver100")[0]
    F1 = ld("p2_aver100", "real")[0] + 1j * ld("p2_aver100", "imag")[0]
    axes = [m["data_corner"][k] + m["data_increment_vector"][k, k] * np.arange(n) for k, n in enumerate(I0.shape)]
    G = np.meshgrid(*axes, indexing="ij")
    bragg = np.all([np.abs(g - np.round(g)) < 1e-9 for g in G], axis=0)
    Fav = F0 - F1
    p2 = {"n_bragg_voxels": int(bragg.sum()),
          "int0_vs_absF0sq_rel": rel(I0, abs(F0) ** 2),
          "int100_vs_absF100sq_rel": rel(I1, abs(F1) ** 2),
          "avgF_offbragg_max_abs": float(np.abs(Fav[~bragg]).max()),
          "avgF_vs_F0_on_bragg_rel": rel(Fav[bragg], F0[bragg]),
          "int100_on_bragg_max": float(I1[bragg].max()), "int0_on_bragg_max": float(I0[bragg].max()),
          "int100_vs_int0_offbragg_rel": rel(I1[~bragg], I0[~bragg])}
    for r in ("p2_aver50_s8a", "p2_aver50_s8b", "p2_aver50_s8c", "p2_aver50_s1a", "p2_aver50_s1b"):
        I = ld(r)[0]
        p2[r] = {"sha": digest(I), "used_percent": re.findall(r"Used\s+([\d.]+)\s*%", log(r)),
                 "offbragg_vs_int0_rel": rel(I[~bragg], I0[~bragg]),
                 "on_bragg_max_over_int0_max": float(I[bragg].max() / I0[bragg].max())}
    out["P2"] = p2

    # P3 lots
    off, m3 = ld("p3_off")
    p3 = {}
    for r in ("p3_box16_n2_y", "p3_box4_all_y", "p3_box4_all_n", "p3_ell4_all_y", "p3_box4_n64_s8a",
              "p3_box4_n64_s8b", "p3_box4_n64_s8c", "p3_box4_all_y_aver100"):
        I = ld(r)[0]
        nl = re.findall(r"Fourier volume\s*:\s*(\d+)\s+(\w+) shaped lots", log(r))
        p3[r] = {"sha": digest(I), "lots_reported": nl, "rel_vs_off": rel(I, off),
                 "bragg_sum_ratio": float(I[bragg].sum() / off[bragg].sum()),
                 "total_sum_ratio": float(I.sum() / off.sum())}
    p3["p3_box4_amp_errors"] = errors("p3_box4_amp")
    p3["p3_box17_err_errors"] = errors("p3_box17_err")
    out["P3"] = p3

    # P4 technique
    p4 = {}
    for g in ("g32", "gnc"):
        for a in ("a0", "a100"):
            t = ld(f"p4_turbo_{g}_{a}")[0]
            entry = {"turbo_max": float(t.max()), "nufft_errors": errors(f"p4_nufft_{g}_{a}")}
            path = E / f"p4_nufft_{g}_{a}" / f"p4_nufft_{g}_{a}_int.h5"
            if path.is_file():
                n = ld(f"p4_nufft_{g}_{a}")[0]
                entry.update(nufft_nonzero=int(np.count_nonzero(n)), nufft_vs_turbo_rel=rel(n, t))
            p4[f"{g}_{a}"] = entry
    p4["turbo_g32_a0_equals_p3_off"] = bool(np.array_equal(ld("p4_turbo_g32_a0")[0], off))
    out["P4"] = p4

    # P5 symmetry
    d = ld("p5_default")[0]
    p5 = {}
    for r in ("p5_symm_ignore", "p5_symm_friedel", "p5_symm_apply"):
        I = ld(r)[0]
        p5[r] = {"identical_to_default": bool(np.array_equal(I, d)), "rel_vs_default": rel(I, d),
                 "log": [l.strip() for l in log(r).splitlines() if "Symmetry averaging" in l][:1]}
    out["P5"] = p5

    print(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
