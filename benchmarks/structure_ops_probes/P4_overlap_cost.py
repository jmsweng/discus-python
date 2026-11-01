"""P4: Overlap check cost for 4096-atom and 64000-atom structures."""
import time
import numpy as np
from discuspdf import Lattice, Structure

def _gaas_supercell(N):
    """GaAs zinc-blende N x N x N supercell as a configuration (8 atoms per cell)."""
    a = 5.6537
    fcc = np.array([[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]])
    base = np.vstack([fcc, fcc + 0.25])
    cells = np.array(np.meshgrid(range(N), range(N), range(N), indexing="ij")).reshape(3, -1).T
    frac = ((base[None, :, :] + cells[:, None, :]) / N).reshape(-1, 3)
    labels = [f"Ga{i}" for i in range(1, 5)] + [f"As{i}" for i in range(1, 5)]
    return dict(
        lattice=Lattice.from_parameters(a * N, a * N, a * N, 90.0, 90.0, 90.0),
        species=np.tile(["Ga"] * 4 + ["As"] * 4, len(cells)).astype(object),
        frac=frac,
        site_label=np.tile(labels, len(cells)).astype(object),
        cell_index=np.repeat(cells, 8, axis=0),
        supercell=(N, N, N),
    )

# 4096 atoms: GaAs 8x8x8 (8*8*8=512 cells * 8 atoms = 4096)
d = _gaas_supercell(8)
print(f"4096-atom structure: {d['lattice'].parameters}")
t0 = time.perf_counter()
s = Structure(d["lattice"], d["species"], d["frac"], model="configuration",
              site_label=d["site_label"], cell_index=d["cell_index"], supercell=d["supercell"])
elapsed = time.perf_counter() - t0
print(f"4096 atoms: {elapsed:.3f} s")

# 64000 atoms: use a larger supercell. 8 atoms per cell, need 8000 cells.
# 20x20x20 = 8000 cells * 8 = 64000 atoms
d2 = _gaas_supercell(20)
print(f"64000-atom structure: {d2['lattice'].parameters}")
t0 = time.perf_counter()
s2 = Structure(d2["lattice"], d2["species"], d2["frac"], model="configuration",
               site_label=d2["site_label"], cell_index=d2["cell_index"], supercell=d2["supercell"])
elapsed2 = time.perf_counter() - t0
print(f"64000 atoms: {elapsed2:.3f} s")
