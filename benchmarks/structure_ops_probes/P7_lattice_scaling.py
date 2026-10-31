"""P7: Lattice scaling by (2,3,4)."""
import numpy as np
from discuspdf import Lattice

# Start with a cubic lattice
lat = Lattice.from_parameters(5.0, 5.0, 5.0, 90.0, 90.0, 90.0)
print(f"Original lattice.parameters = {lat.parameters}")
print(f"Original lattice.matrix =\n{lat.matrix}")

# Scale rows by (2, 3, 4)
A = lat.matrix
A_scaled = A * np.array([2.0, 3.0, 4.0])[:, None]  # scale rows
print(f"\nScaled matrix (rows * [2,3,4]):\n{A_scaled}")

# Build new lattice from scaled matrix
lat2 = Lattice.from_matrix(A_scaled)
print(f"\nScaled lattice.parameters = {lat2.parameters}")
print(f"Scaled lattice.matrix =\n{lat2.matrix}")

# Verify: a=10, b=15, c=20, angles unchanged
a, b, c, alpha, beta, gamma = lat2.parameters
print(f"\nCheck: a={a:.4f} (expect 10), b={b:.4f} (expect 15), c={c:.4f} (expect 20)")
print(f"Check: alpha={alpha:.4f} (expect 90), beta={beta:.4f} (expect 90), gamma={gamma:.4f} (expect 90)")
