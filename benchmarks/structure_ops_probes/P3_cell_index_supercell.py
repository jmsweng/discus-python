"""P3: Cell-index/supercell semantics for 2x2x2 GaAs supercell."""
import numpy as np
from discuspdf import Lattice, Structure

# Copy the _gaas_supercell helper from tests/test_structure.py
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

d = _gaas_supercell(2)
s = Structure(d["lattice"], d["species"], d["frac"], model="configuration",
              site_label=d["site_label"], cell_index=d["cell_index"], supercell=d["supercell"])

print(f"lattice.parameters = {s.lattice.parameters}")
print(f"supercell = {s.supercell}")
print(f"n_atoms = {s.n_atoms}")
print(f"cell_index[:20] = {s.cell_index[:20].tolist()}")
print(f"frac min = {s.frac.min(axis=0)}")
print(f"frac max = {s.frac.max(axis=0)}")
print(f"atom_id[:20] = {s.atom_id[:20].tolist()}")
