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

# 2. Build a 2×2×2 P1 GaAs Structure (8 atoms per cell: 4 Ga + 4 As).
print("\n=== Build 2×2×2 P1 GaAs ===")
lat = Lattice.from_parameters(5.6537, 5.6537, 5.6537, 90.0, 90.0, 90.0)
# 8 atoms per unit cell: the four fcc positions for Ga and the same shifted
# by (1/4, 1/4, 1/4) for As.
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
m = len(unit_frac)  # atoms per unit cell (8)
N = m * n[0] * n[1] * n[2]  # 64
# cell-major, x fastest (vectorised)
cell_num = np.arange(N) // m
ix = (cell_num % n[0]).astype(np.int64)
iy = ((cell_num // n[0]) % n[1]).astype(np.int64)
iz = (cell_num // (n[0] * n[1])).astype(np.int64)
cell_index = np.stack([ix, iy, iz], axis=1)
slot = np.arange(N) % m
species = unit_species[slot]
frac = (unit_frac[slot] + cell_index) / np.array(n, dtype=float)
s2 = Structure(
    lat, species, frac,
    model="configuration",
    b_iso=np.full(N, 0.35),
    cell_index=cell_index,
    supercell=n,
    pbc=(True, True, True),
)
print(f"atoms: {s2.n_atoms}")
from collections import Counter
print(f"species: {dict(Counter(s2.species))}")
# Neighbour check: in the P1 8-atom cell the Ga (and As) sit on the fcc
# sites, so the nearest-neighbour distance is a/sqrt(2) = 3.9978 A and every
# atom has 4 neighbours at that distance (a Ga-As pair is closer, at
# a*sqrt(3)/4 = 2.4480 A, but same-species neighbours are 3.9978 A apart).
# (The 4-neighbour / 2.4480 A result is for the F-43m conventional cell in
# gaas_222.stru, which the live tests check.)
D = s2.distance_matrix()
np.fill_diagonal(D, 99.0)
nn = (D < 1.3).sum(axis=1)
print(f"nearest-neighbour check: min distance {D.min():.4f} A "
      f"(expected {5.6537 / 2 ** 0.5:.4f}), "
      f"neighbours at that distance per atom: {sorted(set(nn.tolist()))}")

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
    # it to 64 for a 2×2×2 supercell, but the Python API does not expand.
    s_cll = read_cll(cll_path)
    # The .stru we wrote carries an ncell line, so it reads back as the
    # 2×2×2 supercell (64 atoms).
    s_stru = read_stru(stru_path)
    print(f".cll read-back (unit cell): {s_cll.structure.n_atoms} atoms")
    print(f".stru read-back (2×2×2): {s_stru.structure.n_atoms} atoms")

    # A DISCUS-saved .stru has no ncell line (DISCUS `save` drops it), so the
    # 64-atom centred gaas_222 fixture is read with supercell="infer".
    s222 = read_stru("tests/fixtures/discus_6_23_01/gaas_222.stru", supercell="infer")
    print(f"gaas_222.stru (supercell='infer'): {s222.structure.n_atoms} atoms, "
          f"supercell {tuple(s222.structure.supercell)}")

    # 5. Print the first 6 lines of each written file.
    print("\n--- gaas.cll (first 6 lines) ---")
    for line in cll_path.read_text().splitlines()[:6]:
        print(line)
    print("\n--- gaas.stru (first 6 lines) ---")
    for line in stru_path.read_text().splitlines()[:6]:
        print(line)

print("\nDone.")
