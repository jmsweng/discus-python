"""Package C engine probe: DISCUS ``chem corr`` versus discuspdf.correlations.

Writes, into ``inputs/`` next to this file:

* three ``.stru`` crystals (written with ``discuspdf.write_stru``, so in DISCUS order);
* one macro per crystal that defines neighbor vectors / distances and runs
  ``chem corr occ`` or ``chem corr disp`` and prints every ``res[i]``;
* ``predictions.json``: our values for the same neighborhoods, computed from the
  **written files read back** (so file rounding is included), plus SHA-256 of every input.

The predictions are fixed before DISCUS is run.  Run on the engine machine with
``run_probes.sh``; compare with ``analyze.py``.  No discuspdf code runs on the
engine side.

Usage: python benchmarks/discus_6_23_01/c_probes/generate.py
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np

from discuspdf import Lattice, Structure, neighbor_list
from discuspdf.correlations import displacement_correlation, occupational_correlation
from discuspdf.discus_io import read_stru, write_stru

HERE = Path(__file__).resolve().parent
INPUTS = HERE / "inputs"


def _cell(base, species, a, labels):
    return Structure(Lattice.from_parameters(a, a, a, 90, 90, 90), np.array(species, object), base,
                     model="configuration", site_label=labels)


def _write(name: str, structure: Structure) -> Structure:
    path = INPUTS / f"{name}.stru"
    write_stru(path, structure, title=f"c_probe {name}", overwrite=True)
    back = read_stru(path).structure
    # read_stru labels sites by species; DISCUS numbers sites by their slot in the cell.
    m = back.n_atoms // int(np.prod(back.supercell))
    return back.replace(site_label=np.array([f"site{r % m + 1}" for r in range(back.n_atoms)], dtype=object))


# ------------------------------------------------------------------ crystals
def crystal_occ() -> Structure:
    """bcc 8x8x8, Au/Cu with partial B2 order and site-dependent concentrations (seed 20261101)."""
    rng = np.random.default_rng(20261101)
    ref = _cell([[0, 0, 0], [0.5, 0.5, 0.5]], ["Au", "Cu"], 3.0, ["M1", "M2"]).repeat((8, 8, 8))
    on1 = ref.site_label == "M1"
    p_au = np.where(on1, 0.75, 0.30)                    # site 1 mostly Au, site 2 mostly Cu
    spec = np.where(rng.random(ref.n_atoms) < p_au, "Au", "Cu")
    return ref.replace(species=spec.astype(object))


def crystal_disp() -> Structure:
    """sc 12x6x6, Cu, longitudinal wave h=2 along a plus transverse sine along b, plus noise (seed 7)."""
    rng = np.random.default_rng(7)
    ref = _cell([[0, 0, 0]], ["Cu"], 3.0, ["C1"]).repeat((12, 6, 6))
    phase = 2 * np.pi * 2 * ref.cell_index[:, 0] / 12
    u = np.zeros((ref.n_atoms, 3))
    u[:, 0] = 0.08 * np.cos(phase)
    u[:, 1] = 0.05 * np.sin(phase)
    u += rng.normal(scale=0.01, size=u.shape)
    return ref.with_cartesian(ref.cartesian + u)


def crystal_mixed() -> Structure:
    """bcc 6x6x6, random Au/Cu (50 %), species-dependent random displacements (seed 99)."""
    rng = np.random.default_rng(99)
    ref = _cell([[0, 0, 0], [0.5, 0.5, 0.5]], ["Au", "Au"], 3.2, ["M1", "M2"]).repeat((6, 6, 6))
    spec = np.where(rng.random(ref.n_atoms) < 0.5, "Au", "Cu")
    ref = ref.replace(species=spec.astype(object))
    u = rng.normal(scale=0.04, size=(ref.n_atoms, 3))
    u[spec == "Cu", 2] += 0.03                       # Cu shifted along +c on average
    return ref.with_cartesian(ref.cartesian + u)


# ------------------------------------------------------------------ neighborhoods
# DISCUS: set vec, iv, is1, is2, dx, dy, dz  (sites numbered 1.. in the cell)
OCC_VECS = [  # (is1, is2, translation)
    (1, 2, (0, 0, 0)),
    (1, 2, (-1, -1, -1)),
    (1, 1, (1, 0, 0)),
    (2, 2, (0, 1, 0)),
    (2, 1, (1, 0, 1)),
]
OCC_DIST = [((0.5, 0.5, 0.5), 0.1), ((1.0, 0.0, 0.0), 0.1)]      # set neig, dis, u, v, w, fsig

DISP_VECS = [  # (is1, is2, translation, dir1, dir2)
    (1, 1, (1, 0, 0), (1, 0, 0), (1, 0, 0)),
    (1, 1, (2, 0, 0), (1, 0, 0), (1, 0, 0)),
    (1, 1, (3, 0, 0), (1, 0, 0), (1, 0, 0)),
    (1, 1, (0, 1, 0), (1, 0, 0), (1, 0, 0)),
    (1, 1, (1, 0, 0), (1, 0, 0), (0, 1, 0)),
    (1, 1, (1, 0, 0), (0, 1, 0), (0, 1, 0)),
]

MIXED_VECS = [
    (1, 2, (0, 0, 0), (0, 0, 1), (0, 0, 1)),
    (1, 1, (0, 0, 1), (0, 0, 1), (0, 0, 1)),
    (1, 2, (0, 0, 0), (1, 1, 0), (1, 1, 0)),
]


def _vec_lines(vecs):
    lines = []
    for iv, v in enumerate(vecs, start=1):
        is1, is2, t = v[0], v[1], v[2]
        lines.append(f"set vec, {iv}, {is1}, {is2}, {t[0]}, {t[1]}, {t[2]}")
    return lines


def _macro_occ(name):
    lines = ["discus", "read", f"  stru {name}.stru", "chem", "set mode, quick, period"]
    lines += _vec_lines(OCC_VECS)
    lines.append("set neig, rese")
    nn = 0
    for iv in range(1, len(OCC_VECS) + 1):
        nn += 1
        lines.append(f"set neig, vec, {iv}, number:{nn}")
    for (u, fsig) in OCC_DIST:
        nn += 1
        lines.append(f"set neig, dis, {u[0]}, {u[1]}, {u[2]}, {fsig}, number:{nn}")
    lines += ["show", "corr occ, Au, Cu"]
    lines += [f"eval res[{k}]" for k in range(0, nn + 1)]
    lines += ["corr occ, Cu, Au"]
    lines += [f"eval res[{k}]" for k in range(0, nn + 1)]
    lines += ["exit", "exit", "exit"]
    return "\n".join(lines) + "\n"


def _macro_disp(name, vecs, pairs):
    lines = ["discus", "read", f"  stru {name}.stru", "chem", "set mode, quick, period"]
    lines += _vec_lines(vecs)
    lines.append("set neig, rese")
    for iv, v in enumerate(vecs, start=1):
        d1, d2 = v[3], v[4]
        lines.append(f"set neig, vec, {iv}, number:{iv}")
        lines.append(f"set neig, dir, {d1[0]}, {d1[1]}, {d1[2]}, {d2[0]}, {d2[1]}, {d2[2]}, number:{iv}")
    lines.append("show")
    for a, b in pairs:
        lines.append(f"corr disp, {a}, {b}")
        lines += [f"eval res[{k}]" for k in range(0, len(vecs) + 1)]
    lines += ["exit", "exit", "exit"]
    return "\n".join(lines) + "\n"


def _labels(structure):
    """DISCUS site number (1-based, order in the cell) -> our site label."""
    m = structure.n_atoms // int(np.prod(structure.supercell))
    return {k + 1: str(structure.site_label[k]) for k in range(m)}


def _vec_selection(nl, labels, is1, is2, t):
    """DISCUS ``set vec``: both directions when is1 != is2; only +t for every atom when is1 == is2."""
    return nl.select(origin_site=labels[is1], target_site=labels[is2], translation=t, symmetric=is1 != is2)


def predictions_occ(s: Structure) -> dict:
    labels = _labels(s)
    nl = neighbor_list(s, 6.0)
    out = []
    for is1, is2, t in OCC_VECS:
        sel = _vec_selection(nl, labels, is1, is2, t)
        res = occupational_correlation(sel)
        n_discus = sel.n_pairs
        out.append({"kind": "vec", "def": [is1, is2, list(t)], "n_pairs_discus": n_discus,
                    "corr_AuCu": res.discus("Au", "Cu"), "corr_CuAu": res.discus("Cu", "Au"),
                    "pearson_AuAu": res.value("pearson", "Au", "Au"),
                    "warren_cowley_AuCu": res.value("warren_cowley", "Au", "Cu")})
    for u, fsig in OCC_DIST:
        d = math.sqrt(sum((x * 3.0) ** 2 for x in u))
        sel = nl.select(distance=(d - fsig / 2, d + fsig / 2))   # shells are > 1 Å apart here
        res = occupational_correlation(sel)
        out.append({"kind": "dis", "def": [list(u), fsig], "distance": d, "n_pairs_ours": sel.n_pairs,
                    "corr_AuCu": res.discus("Au", "Cu"), "corr_CuAu": res.discus("Cu", "Au")})
    return {"neighborhoods": out}


def predictions_disp(s: Structure, vecs, pairs) -> dict:
    labels = _labels(s)
    nl = neighbor_list(_undisplaced(s), 10.0)
    out = []
    for is1, is2, t, d1, d2 in vecs:
        sel = _vec_selection(nl, labels, is1, is2, t)
        row = {"def": [is1, is2, list(t), list(d1), list(d2)]}
        for a, b in pairs:
            kw = {} if (a, b) == ("all", "all") else {"origin_species": a, "target_species": b}
            res = displacement_correlation(sel, s, direction=d1, target_direction=d2,
                                           direction_frame="fractional", reference="site_average",
                                           centered=False, **kw)
            pear = displacement_correlation(sel, s, direction=d1, target_direction=d2,
                                            direction_frame="fractional", **kw)
            row[f"{a}_{b}"] = {"discus_convention": res.correlation, "pearson": pear.correlation,
                               "n_pairs": res.n_pairs}
        out.append(row)
    return {"neighborhoods": out}


def _undisplaced(s: Structure) -> Structure:
    """Ideal positions: site position of the parent cell + cell_index (same sites and cells)."""
    m = s.n_atoms // int(np.prod(s.supercell))
    sc = s.supercell.astype(float)
    parent = (s.frac * sc - s.cell_index)[:m]
    ideal = np.round(parent * 8) / 8          # the probe crystals have sites at multiples of 1/2
    frac = (ideal[np.arange(s.n_atoms) % m] + s.cell_index) / sc
    return s.replace(frac=frac)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    INPUTS.mkdir(exist_ok=True)
    preds = {"discus_version": "6.23.01", "note": "values from the written files read back"}
    occ = _write("occ_bcc", crystal_occ())
    (INPUTS / "occ_bcc.mac").write_text(_macro_occ("occ_bcc"), encoding="ascii", newline="\n")
    preds["occ_bcc"] = predictions_occ(occ)

    disp = _write("disp_sc", crystal_disp())
    (INPUTS / "disp_sc.mac").write_text(_macro_disp("disp_sc", DISP_VECS, [("Cu", "Cu")]),
                                        encoding="ascii", newline="\n")
    preds["disp_sc"] = predictions_disp(disp, DISP_VECS, [("Cu", "Cu")])

    mixed = _write("mixed_bcc", crystal_mixed())
    pairs = [("Au", "Cu"), ("Cu", "Cu"), ("all", "all")]
    (INPUTS / "mixed_bcc.mac").write_text(_macro_disp("mixed_bcc", MIXED_VECS, pairs),
                                          encoding="ascii", newline="\n")
    preds["mixed_bcc"] = predictions_disp(mixed, MIXED_VECS, pairs)

    preds["input_sha256"] = {p.name: sha256(p) for p in sorted(INPUTS.iterdir()) if p.suffix in (".stru", ".mac")}
    (HERE / "predictions.json").write_text(json.dumps(preds, indent=1, sort_keys=True) + "\n",
                                           encoding="ascii", newline="\n")
    print(json.dumps({k: v for k, v in preds.items() if k != "input_sha256"}, indent=1)[:3000])


if __name__ == "__main__":
    main()
