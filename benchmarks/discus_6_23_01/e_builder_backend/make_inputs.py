"""Generate the builder-backend probe (DISCUS_SUITE 6.23.01): own neutron scattering lengths (``scat``),
VOID padding atoms and the ``.stru`` source of the powder PDF.

Run from the repository root:  python benchmarks/discus_6_23_01/e_builder_backend/make_inputs.py
Writes inputs/<run>/ (run.mac and the structure file each macro reads) and runs.json. The macros come from
discuspdf itself (FourierConfig.render, render_pdf_macro), so the probe tests the generated text. The
generated files are committed before the engine run and must not be edited by hand; PREREG.md fixes the
criteria in advance.
"""
from __future__ import annotations

import json
import shutil
from pathlib import Path

import numpy as np

from discuspdf import FourierConfig, Lattice, PDFConfig, ReciprocalGrid, Seed, StructureSource
from discuspdf.discus_io import format_cll, format_stru
from discuspdf.pdf import render_pdf_macro
from discuspdf.structure import Structure
from discuspdf.structure_ops import repeat

HERE = Path(__file__).resolve().parent
INPUTS = HERE / "inputs"

A = 3.52                      # fcc Ni, Å
N = (4, 4, 4)                 # cells
FCC = [[0, 0, 0], [.5, .5, 0], [.5, 0, .5], [0, .5, .5]]
VOID_POS = [.25, .25, .25]    # one padding atom per cell, on an empty interstitial position
B_DISCUS = {"NI": 10.3, "CU": 7.718}   # DISCUS 6.23.01 per_bcoh (fm)
SEED = Seed((20261, 1103, 4242, 777, 31415, 2718, 1618, 99))
PDF = dict(rmax=10.0, dr=0.02, qmin=0.5, qmax=7.0, dq=0.01)


def lattice() -> Lattice:
    return Lattice.from_parameters(A, A, A, 90, 90, 90)


def _cu_cells() -> np.ndarray:
    """Which of the 256 fcc rows (cell-major, 4 per cell) are Cu: 64 drawn without replacement."""
    rng = np.random.default_rng(20261103)
    return np.sort(rng.permutation(4 * int(np.prod(N)))[:64])


def structure(kind: str) -> Structure:
    """'mix' Ni with 25 % Cu; 'relabel' the same with Cu renamed NI63 (DISCUS reads it as Ni);
    'void' = 'mix' plus one VOID per cell; 'ni' pure Ni."""
    sites = FCC + ([VOID_POS] if kind == "void" else [])
    cell = Structure(lattice(), ["NI"] * len(sites), sites, model="configuration",
                     site_label=[f"s{i + 1}" for i in range(len(sites))])
    st = repeat(cell, N).structure
    species = np.array(st.species, dtype=object)
    per = len(sites)
    real = np.array([i for i in range(st.n_atoms) if i % per < 4])      # fcc rows, cell-major
    if kind in ("mix", "relabel", "void"):
        species[real[_cu_cells()]] = "NI63" if kind == "relabel" else "CU"
    if kind == "void":
        species[np.arange(st.n_atoms) % per == 4] = "VOID"
    return st.replace(species=species, b_iso=np.ones(st.n_atoms))


def unit_cell() -> Structure:
    return Structure(lattice(), ["NI"] * 4, FCC, model="configuration", site_label=["s1", "s2", "s3", "s4"],
                     b_iso=np.ones(4))


def grid() -> ReciprocalGrid:
    return ReciprocalGrid.centered(lattice(), (17, 17, 17), (0.25, 0.25, 0.25))


RUNS = [
    # 3D neutron, Bragg subtracted (set aver,100)
    {"id": "f_mix", "type": "fourier", "structure": "mix", "scat": {}},
    {"id": "f_relabel_scat", "type": "fourier", "structure": "relabel", "scat": {"NI63": B_DISCUS["CU"]}},
    {"id": "f_relabel_noscat", "type": "fourier", "structure": "relabel", "scat": {}},
    {"id": "f_void", "type": "fourier", "structure": "void", "scat": {}},
    # 1D powder PDF from a .stru (read stru)
    {"id": "p_mix", "type": "pdf", "structure": "mix", "radiation": "neutron", "scat": {}},
    {"id": "p_relabel_scat", "type": "pdf", "structure": "relabel", "radiation": "neutron",
     "scat": {"NI63": B_DISCUS["CU"]}},
    {"id": "p_relabel_noscat", "type": "pdf", "structure": "relabel", "radiation": "neutron", "scat": {}},
    {"id": "p_void", "type": "pdf", "structure": "void", "radiation": "neutron", "scat": {}},
    {"id": "p_ni_stru_n", "type": "pdf", "structure": "ni", "radiation": "neutron", "scat": {}},
    {"id": "p_ni_stru_x", "type": "pdf", "structure": "ni", "radiation": "xray", "scat": {}},
    # the verified .cll path (read cell, 4x4x4) on the same crystal
    {"id": "p_ni_cll_n", "type": "pdf", "structure": "ni_cll", "radiation": "neutron", "scat": {}},
    {"id": "p_ni_cll_x", "type": "pdf", "structure": "ni_cll", "radiation": "xray", "scat": {}},
]


def macro_and_files(run: dict) -> tuple[str, dict[str, str]]:
    kind = run["structure"]
    if run["type"] == "fourier":
        name = f"{kind}.stru"
        cfg = FourierConfig(grid(), outputs=(("int", f"{run['id']}.h5"),), aver_percent=100.0, seed=SEED,
                            radiation="neutron", scattering_lengths=run["scat"])
        text = cfg.render(StructureSource.stru(name, crystal_cells=N), title=f"e_builder_backend {run['id']}")
        return text, {name: format_stru(structure(kind), title=kind)}
    if kind == "ni_cll":
        cfg = PDFConfig("model.cll", supercell=N, radiation=run["radiation"], **PDF)
        return render_pdf_macro(cfg), {"model.cll": format_cll(unit_cell(), title="ni")}
    cfg = PDFConfig("model.stru", supercell=(1, 1, 1), radiation=run["radiation"],
                    scattering_lengths=run["scat"], **PDF)
    return render_pdf_macro(cfg), {"model.stru": format_stru(structure(kind), title=kind)}


def main() -> None:
    if INPUTS.exists():
        shutil.rmtree(INPUTS)
    for run in RUNS:
        d = INPUTS / run["id"]
        d.mkdir(parents=True)
        text, files = macro_and_files(run)
        (d / "run.mac").write_text(text, encoding="ascii", newline="\n")
        for name, body in files.items():
            (d / name).write_text(body, encoding="ascii", newline="\n")
    spec = {"engine": "DISCUS_SUITE 6.23.01", "lattice_a": A, "cells": list(N), "b_discus": B_DISCUS,
            "pdf": PDF, "grid": grid().to_dict(), "runs": RUNS}
    (HERE / "runs.json").write_text(json.dumps(spec, indent=1) + "\n", encoding="utf-8")
    print(f"wrote {len(RUNS)} runs to {INPUTS}")


if __name__ == "__main__":
    main()
