"""DISCUS file I/O basics — offline example.

Reads ``examples/gaas.cll``, prints the listed atoms; builds a 2×2×2 P1 GaAs
``Structure`` with numpy; writes ``.cll`` and ``.stru`` into a temporary
directory; reads them back and prints atom counts; prints the first 6 lines
of each written file.  Runs offline in under 5 seconds.
"""
from __future__ import annotations

import tempfile
from pathlib import Path

import numpy as np

from discuspdf import (
    Lattice,
    Structure,
    read_cll,
    read_stru,
    write_cll,
    write_stru,
)

# 1. Read the committed .cll fixture.
print("=== Read examples/gaas.cll ===")
s = read_cll("examples/gaas.cll")
print(f"atoms: {s.structure.n_atoms}")
for i in range(s.structure.n_atoms):
    x, y, z = s.structure.frac[i]
    b = 0.0 if s.structure.b_iso is None else s.structure.b_iso[i]
    print(f"  {s.structure.species[i]:<4} {x:.6f} {y:.6f} {z:.6f} B={b:.6f}")

# 2. Build a 2×2×2 P1 GaAs Structure.
print("\n=== Build 2×2×2 P1 GaAs ===")
lat = Lattice.from_parameters(5.6537, 5.6537, 5.6537, 90.0, 90.0, 90.0)
unit_frac = np.array(
    [
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.0],
        [0.0, 0.5, 0.5],
        [0.5, 0.0, 0.5],
        [0.25, 0.25, 0.25],
        [0.75, 0.75, 0.25],
        [0.25, 0.75, 0.75],
        [0.75, 0.25, 0.75],
    ]
)
unit_species = np.array(
    ["GA", "GA", "GA", "GA", "AS", "AS", "AS", "AS"],
    dtype=object,
)
n = (2, 2, 2)
m = 4  # N / prod(n) = 32 / 8 = 4
N = 32
species = np.empty(N, dtype=object)
frac = np.empty((N, 3))
cell_index = np.zeros((N, 3), dtype=np.int64)
for r in range(N):
    cell_num = r // m
    slot = r % m
    iz = cell_num // 4
    rem = cell_num % 4
    iy = rem // 2
    ix = rem % 2
    cell_index[r] = (ix, iy, iz)
    species[r] = unit_species[slot]
    frac[r] = (unit_frac[slot] + np.array([ix, iy, iz])) / np.array(n)
s2 = Structure(
    lat, species, frac,
    model="configuration",
    b_iso=np.full(N, 0.35),
    cell_index=cell_index,
    supercell=n,
    pbc=(True, True, True),
    min_separation=None,
)
print(f"atoms: {s2.n_atoms}")
from collections import Counter
print(f"species: {dict(Counter(s2.species))}")

# 3. Write .cll (unit cell) and .stru (supercell) into a temporary directory.
print("\n=== Write .cll and .stru ===")
with tempfile.TemporaryDirectory() as tmpdir:
    # Write the unit cell as .cll (8 atoms, supercell=(1,1,1)).
    unit_cell = Structure(
        lat, unit_species, unit_frac,
        model="configuration",
        b_iso=np.full(8, 0.35),
        supercell=(1, 1, 1),
        pbc=(True, True, True),
        min_separation=None,
    )
    cll_path = Path(tmpdir) / "gaas.cll"
    write_cll(cll_path, unit_cell)

    # Write the supercell as .stru (32 atoms, supercell=(2,2,2)).
    stru_path = Path(tmpdir) / "gaas.stru"
    write_stru(stru_path, s2)

    # 4. Read them back and print atom counts.
    # A .cll stores the unit cell (8 atoms for P1 GaAs); DISCUS would expand
    # it to 32 for a 2×2×2 supercell, but the Python API does not expand.
    s_cll = read_cll(cll_path)
    s_stru = read_stru(stru_path)
    print(f".cll read-back (unit cell): {s_cll.structure.n_atoms} atoms")
    print(f".stru read-back (2×2×2): {s_stru.structure.n_atoms} atoms")

    # 5. Print the first 6 lines of each written file.
    print("\n--- gaas.cll (first 6 lines) ---")
    for line in cll_path.read_text().splitlines()[:6]:
        print(line)
    print("\n--- gaas.stru (first 6 lines) ---")
    for line in stru_path.read_text().splitlines()[:6]:
        print(line)

print("\nDone.")
