"""Generate the builder-backend probe v2 (DISCUS_SUITE 6.23.01): VOID atoms in the 1D powder PDF.

Follows v1's finding P3 (../RESULT.md): with VOID padding the powder G(r) is scaled by rows / atoms. v2 tests
(a) the crystal written as one cell without VOID rows (``ncell 1,1,1,N``) against the verified layout, and
(b) the scale law with another VOID fraction (vacancies) and with X-rays.

Run from the repository root:  python benchmarks/discus_6_23_01/e_builder_backend/v2/make_inputs.py
Writes inputs/<run>/ and runs.json. Committed before the engine run; do not edit by hand (PREREG.md).
"""
from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import make_inputs as v1  # noqa: E402  (v1 structures, lattice, PDF settings)

from discuspdf import Lattice, PDFConfig  # noqa: E402
from discuspdf.discus_io import format_stru  # noqa: E402
from discuspdf.pdf import render_pdf_macro  # noqa: E402
from discuspdf.structure import Structure  # noqa: E402

INPUTS = HERE / "inputs"
N_VAC = 16                                 # vacancies in the 'vac' crystal


def vac_rows() -> np.ndarray:
    """Which of the 256 fcc rows of 'mix' become vacancies (drawn without replacement, seed 20261104)."""
    return np.sort(np.random.default_rng(20261104).permutation(256)[:N_VAC])


def structure(kind: str) -> Structure:
    """'mix', 'void': as v1. 'vac': 'mix' with 16 rows renamed VOID (padded layout, 256 rows, 240 atoms).
    '*_big': the real atoms of 'mix' / 'vac' as ONE cell of the whole crystal (no VOID rows, ncell 1,1,1,N)."""
    if kind in ("mix", "void"):
        return v1.structure(kind)
    base = v1.structure("mix")
    species = np.array(base.species, dtype=object)
    if kind.startswith("vac"):
        species[vac_rows()] = "VOID"
    st = base.replace(species=species)
    if not kind.endswith("_big"):
        return st
    keep = species != "VOID"
    big = Lattice.from_parameters(*st.lattice.parameters)          # the whole 4x4x4 crystal as one cell
    return Structure(big, species[keep], st.frac[keep], model="configuration",
                     site_label=[f"a{i}" for i in range(int(keep.sum()))], b_iso=np.ones(int(keep.sum())))


RUNS = [
    {"id": "p_mix_n", "structure": "mix", "radiation": "neutron"},
    {"id": "p_mix_x", "structure": "mix", "radiation": "xray"},
    {"id": "p_mix_big_n", "structure": "mix_big", "radiation": "neutron"},
    {"id": "p_mix_big_x", "structure": "mix_big", "radiation": "xray"},
    {"id": "p_void_x", "structure": "void", "radiation": "xray"},
    {"id": "p_vac_void_n", "structure": "vac", "radiation": "neutron"},
    {"id": "p_vac_void_x", "structure": "vac", "radiation": "xray"},
    {"id": "p_vac_big_n", "structure": "vac_big", "radiation": "neutron"},
    {"id": "p_vac_big_x", "structure": "vac_big", "radiation": "xray"},
]


def main() -> None:
    if INPUTS.exists():
        shutil.rmtree(INPUTS)
    counts = {}
    for run in RUNS:
        d = INPUTS / run["id"]
        d.mkdir(parents=True)
        st = structure(run["structure"])
        counts[run["structure"]] = {"rows": int(st.n_atoms), "atoms": int((np.asarray(st.species) != "VOID").sum()),
                                    "ncell": [int(x) for x in st.supercell]}
        cfg = PDFConfig("model.stru", supercell=(1, 1, 1), radiation=run["radiation"], **v1.PDF)
        (d / "run.mac").write_text(render_pdf_macro(cfg), encoding="ascii", newline="\n")
        (d / "model.stru").write_text(format_stru(st, title=run["structure"]), encoding="ascii", newline="\n")
    spec = {"engine": "DISCUS_SUITE 6.23.01", "follows": "v1 P3 (../RESULT.md)", "pdf": v1.PDF, "counts": counts,
            "runs": RUNS}
    (HERE / "runs.json").write_text(json.dumps(spec, indent=1) + "\n", encoding="utf-8")
    print(f"wrote {len(RUNS)} runs to {INPUTS}; counts {counts}")


if __name__ == "__main__":
    main()
