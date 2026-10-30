"""Structure basics: GaAs conventional cell, a Cu/Au average model, select/edit, distances.

Run:  python examples/structure_basics.py
"""
import numpy as np

from discuspdf import Lattice, Structure

# --- GaAs conventional cell (zinc blende, a = 4.53 Å) ---
a = 4.53
lat = Lattice.from_parameters(a, a, a, 90.0, 90.0, 90.0)
print(f"GaAs cubic cell: a = {a} Å")

# Ga at (0, 0, 0), As at (0.25, 0.25, 0.25)
gaas = Structure(
    lat,
    np.array(["Ga", "As"]),
    np.array([[0.0, 0.0, 0.0], [0.25, 0.25, 0.25]]),
    model="configuration",
)
print(f"  {gaas}")
print(f"  counts: {gaas.counts()}")
print(f"  composition: {gaas.composition()}")
print(f"  formula: {gaas.formula()}")
print(f"  atomic numbers: {gaas.atomic_numbers()}")

# Ga-As distance
d = gaas.distance(0, 1)
print(f"  Ga-As distance: {d:.6f} Å  (expected {a * np.sqrt(3) / 4:.6f})")

# --- Cu/Au average model (FCC, a = 4.0 Å) ---
a_cu = 4.0
lat_cu = Lattice.from_parameters(a_cu, a_cu, a_cu, 90.0, 90.0, 90.0)
cupau = Structure(
    lat_cu,
    np.array(["Cu", "Au"]),
    np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]),
    model="average",
    occupancy=np.array([0.75, 0.25]),
    site_label=np.array(["S1", "S1"]),
)
print(f"\nCu/Au average model:")
print(f"  {cupau}")
print(f"  composition: {cupau.composition()}")
print(f"  formula: {cupau.formula()}")

# --- Select and edit ---
ga_only = gaas.select(species="Ga")
print(f"\nSelect Ga only: {ga_only}")

# Edit: move As
gaas2 = gaas.with_frac(np.array([[0.0, 0.0, 0.0], [0.3, 0.25, 0.25]]))
d2 = gaas2.distance(0, 1)
print(f"After moving As to (0.3, 0.25, 0.25): Ga-As distance = {d2:.6f} Å")

# Original unchanged
d_orig = gaas.distance(0, 1)
print(f"Original Ga-As distance: {d_orig:.6f} Å (unchanged)")

# --- Wrap ---
gaas_wrapped = gaas.wrapped()
print(f"\nWrapped: {gaas_wrapped.frac}")

print("\nDone.")
