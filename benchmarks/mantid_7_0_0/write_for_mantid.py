"""Write MDHisto files with discuspdf for the Mantid check (run with discuspdf's Python).

    python benchmarks/mantid_7_0_0/write_for_mantid.py tests/fixtures/mantid_7_0_0 <out dir>

Each Mantid fixture is read with read_mantid_mdhisto and written back with write_mantid_mdhisto
(<name>.discuspdf.nxs). A synthetic volume on a skewed grid in a monoclinic cell (values
100*i0 + 10*i1 + i2, variance 0.5*value, one masked voxel) is added as skewed.discuspdf.nxs together
with skewed.expected.json, so the Mantid side can check positions and values it did not produce.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

from discuspdf import Lattice, ReciprocalGrid, ScatteringVolume
from discuspdf.mantid_io import read_mantid_mdhisto, write_mantid_mdhisto

SOURCES = ("aligned_hkl.nxs", "projected_hh0.nxs", "w_matrix_after.nxs", "histo_formula.nxs")


def skewed_volume() -> ScatteringVolume:
    lat = Lattice.from_parameters(5, 6, 7, 90, 105, 90)
    V = np.array([[0.25, 0.0, 0.125], [0.0, 0.25, 0.0], [-0.125, 0.0, 0.5]])
    grid = ReciprocalGrid(lat, (-1.0, -1.0, -1.0), V, (6, 5, 4))
    i0, i1, i2 = np.meshgrid(*[np.arange(n) for n in grid.shape], indexing="ij")
    values = 100.0 * i0 + 10.0 * i1 + i2
    mask = np.zeros(grid.shape, bool)
    mask[1, 2, 3] = True
    return ScatteringVolume(grid, values, "arbitrary", variance=0.5 * values, mask=mask, mask_reasons=("test",))


def main(src: Path, out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    for name in SOURCES:
        sv = read_mantid_mdhisto(src / name)
        write_mantid_mdhisto(sv, out / name.replace(".nxs", ".discuspdf.nxs"), overwrite=True)
    sv = skewed_volume()
    write_mantid_mdhisto(sv, out / "skewed.discuspdf.nxs", overwrite=True)
    (out / "skewed.expected.json").write_text(json.dumps({
        "hkl_of_index": {f"{i},{j},{k}": sv.grid.index_to_hkl((i, j, k)).tolist()
                         for i, j, k in [(0, 0, 0), (5, 0, 0), (0, 4, 0), (0, 0, 3), (3, 2, 1)]},
        "value_formula": "100*i0 + 10*i1 + i2", "variance_formula": "0.5*value", "masked_index": [1, 2, 3],
        "cell": list(sv.grid.lattice.parameters)}, indent=1) + "\n", encoding="utf-8")
    print("wrote", sorted(p.name for p in out.iterdir()))


if __name__ == "__main__":
    main(Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve())
