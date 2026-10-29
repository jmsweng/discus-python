"""Lattice basics: cubic and triclinic Q / d-spacings and a PBC distance.

Run:  python examples/lattice_basics.py
"""
import math

import numpy as np

from discuspdf import Lattice

# --- cubic a = 5 Å ---
cubic = Lattice.from_parameters(5.0, 5.0, 5.0, 90.0, 90.0, 90.0)
print("cubic cell:", cubic)
q100 = cubic.hkl_to_q((1, 0, 0))
print(f"  |Q(100)| = {np.linalg.norm(q100):.6f} Å^-1  (expected {2 * math.pi / 5:.6f})")
print(f"  d(100) = {cubic.d_spacing((1, 0, 0)):.6f} Å,  d(111) = {cubic.d_spacing((1, 1, 1)):.6f} Å")
print(f"  volume = {cubic.volume:.4f} Å³")

# --- triclinic ---
tric = Lattice.from_parameters(4.1, 5.3, 6.7, 78.0, 95.0, 112.0)
print("\ntriclinic cell:", tric)
q111 = tric.hkl_to_q((1, 1, 1))
print(f"  Q(111) = [{', '.join(f'{x:.6f}' for x in q111)}] Å^-1")
print(f"  d(111) = {tric.d_spacing((1, 1, 1)):.6f} Å")
print(f"  recovered parameters: {tuple(round(x, 6) for x in tric.parameters)}")

# --- PBC distance: (0,0,0) to (0.9, 0, 0) in the cubic cell ---
d = cubic.distance_pbc((0.0, 0.0, 0.0), (0.9, 0.0, 0.0))
print(f"\nPBC distance (0,0,0) -> (0.9,0,0) = {d:.6f} Å  (expected 0.5)")
