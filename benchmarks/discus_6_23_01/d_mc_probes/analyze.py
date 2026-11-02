"""Check the D MMC probe evidence against PREDICTIONS.md (prints one line per prediction).

    python benchmarks/discus_6_23_01/d_mc_probes/analyze.py [evidence dir] [--json out.json]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

from discuspdf.correlations import displacement_correlation, occupational_correlation
from discuspdf.discus_io import read_stru
from discuspdf.mc import parse_mmc_log
from discuspdf.neighbors import label_sites_by_slot, neighbor_list

HERE = Path(__file__).resolve().parent
SIX = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("evidence", type=Path, nargs="?", default=HERE / "evidence" / "run1")
    ap.add_argument("--json", type=Path)
    a = ap.parse_args(argv)
    E, W = a.evidence, a.evidence / "work"
    runs = [r["id"] for r in json.loads((HERE / "runs.json").read_text())["runs"]]
    log = {r: parse_mmc_log((E / f"{r}.log").read_text(errors="replace")) for r in runs}
    inp = {n: label_sites_by_slot(read_stru(W / n).structure) for n in ("cuau.stru", "cu_disp.stru")}
    src = {r["id"]: r["structure"] for r in json.loads((HERE / "runs.json").read_text())["runs"]}
    # DISCUS 'save' drops the ncell line; the crystal size is that of the input
    st = {r: label_sites_by_slot(read_stru(W / f"{r}_out.stru",
                                           supercell=tuple(inp[src[r]].supercell)).structure) for r in runs}
    res: dict = {}

    def report(key, ok, text):
        res[key] = {"ok": ok, "text": text}
        print(f"{key:4s} {'PASS' if ok else ('INFO' if ok is None else 'FAIL')}  {text}")

    codes = dict(l.split(" exit=")[0:2] for l in (E / "exit_codes.txt").read_text().splitlines())
    sha = {l.split()[1]: l.split()[0] for l in (E / "output_sha256.txt").read_text().splitlines()}
    report("P1", sha["r1_occ_corr_out.stru"] == sha["r1b_occ_corr_repeat_out.stru"]
           and log["r1_occ_corr"]["final"] == log["r1b_occ_corr_repeat"]["final"],
           "serial runs with the same seed: saved structures and final tables identical")
    f_c, f_d = log["r1c_occ_corr_omp4"]["final"], log["r1d_occ_corr_omp4_repeat"]["final"]
    report("P2", None, f"4 threads: outputs identical = {sha['r1c_occ_corr_omp4_out.stru'] == sha['r1d_occ_corr_omp4_repeat_out.stru']}; "
           f"moves {f_c['gen']} / {f_d['gen']} (requested 40960); achieved {f_c['rows'][0]['achieved']} / {f_d['rows'][0]['achieved']}")
    comp0 = {s: int(np.sum(np.asarray(inp['cuau.stru'].species) == s)) for s in ("Cu", "Au")}
    ok3 = all({s: int(np.sum(np.char.upper(np.asarray(st[r].species, str)) == s.upper())) for s in comp0} == comp0
              for r in runs if r.startswith(("r1", "r2", "r3")))
    report("P3", ok3, f"swchem conserves composition {comp0} in all occupational runs")
    report("P4", log["r1_occ_corr"]["final"]["gen"] == 40960, f"final Gen = {log['r1_occ_corr']['final']['gen']} for set cyc 40960")

    # P5: printed achieved vs Package C on the saved structure (six directed vectors)
    s1 = st["r1_occ_corr"]
    nl = neighbor_list(s1, 3.1)
    mask = np.zeros(nl.n_pairs, bool)
    for t in SIX:
        mask |= np.all(nl.translation == np.array(t), axis=1)
    pairs = nl._subset(mask)
    occ = occupational_correlation(pairs, s1)
    alpha = occ.discus("CU", "AU")   # DISCUS saves names in upper case
    printed = log["r1_occ_corr"]["final"]["rows"][0]
    report("P5", abs(alpha - printed["achieved"]) <= 0.0005 + 1e-12 and printed["pairs"] == pairs.n_pairs,
           f"printed achieved {printed['achieved']} (pairs {printed['pairs']}) vs Package C {alpha:.6f} (pairs {pairs.n_pairs})")
    report("P6", None, f"CORR target -0.3, final achieved {printed['achieved']}")
    for r, key in (("r2_occ_ener_zero", "P7"), ("r3_occ_ener_hot", "P8")):
        f = log[r]["final"]
        g, n, b = f["accepted"]
        frac = (g + n + b) / f["gen"]
        nondown = (n + b) / (f["gen"] - g)
        ok = abs(frac - 0.5) <= 0.02 if key == "P7" else abs(frac - 0.5) <= 0.02
        report(key, ok, f"{r}: accepted {g}/{n}/{b} of {f['gen']} = {frac:.4f}; accepted share of non-downhill "
                         f"moves {nondown:.4f}")
    report("P9", log["r3b_occ_ener_cold"]["final"]["rows"][0]["achieved"] < 0,
           f"J=+0.5: achieved {log['r3b_occ_ener_cold']['final']['rows'][0]['achieved']}")
    d0 = inp["cu_disp.stru"].frac
    dmax = {r: float(np.abs(st[r].frac - d0).max() * st[r].supercell.max()) for r in
            ("r4_shift_spring", "r4b_shift_no_move", "r5_shift_lennard", "r6_swdisp_cd")}
    report("P10", dmax["r4b_shift_no_move"] < 2e-6, f"no set move: max position change {dmax['r4b_shift_no_move']:.2e} cell")
    ini = {r: log[r]["initial"][0]["sigma"] for r in ("r4_shift_spring", "r5_shift_lennard")}
    fin = {r: log[r]["final"]["rows"][0]["sigma"] for r in ini}
    report("P11", dmax["r4_shift_spring"] > 1e-3 and all(fin[r] < ini[r] for r in ini),
           f"atoms moved (max {dmax['r4_shift_spring']:.3f} cell); NN distance sigma initial -> final: "
           + ", ".join(f"{r}: {ini[r]} -> {fin[r]}" for r in ini))
    s6 = st["r6_swdisp_cd"]
    ideal = inp["cu_disp.stru"].replace(frac=inp["cu_disp.stru"].cell_index / 12.0)   # undisplaced sites
    nl6 = neighbor_list(ideal, 3.1)
    m6 = np.zeros(nl6.n_pairs, bool)
    for t in SIX:
        m6 |= np.all(nl6.translation == np.array(t), axis=1)
    dc = displacement_correlation(nl6._subset(m6), s6, direction=(1, 0, 0), direction_frame="fractional",
                                  reference="site_average", centered=False)
    pr6 = log["r6_swdisp_cd"]["final"]["rows"][0]
    report("P12", abs(dc.correlation - pr6["achieved"]) <= 0.0005 + 1e-12,
           f"printed cd achieved {pr6['achieved']} (pairs {pr6['pairs']}) vs Package C uncentred {dc.correlation:.6f}")
    report("P13", all(v.split()[0] == "0" for v in codes.values()), f"exit codes {sorted(set(v.split()[0] for v in codes.values()))}")
    if a.json:
        a.json.write_text(json.dumps(res, indent=1) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
