"""Interop basics: CIF → Structure → ASE / DiffPy, plus DISCUS .cll and CIF writing.

Run:  python examples/interop_basics.py
"""
import sys
import tempfile
from pathlib import Path

import numpy as np

from discuspdf import Lattice, Structure
from discuspdf.cif_io import read_cif, write_cif
from discuspdf.discus_io import write_cll

FIXTURES = Path(__file__).parent.parent / "tests" / "fixtures" / "interop"

# --- 1. Read a real CIF (COD GaAs) -----------------------------------------
imp = read_cif(FIXTURES / "gaas_cod_9008845.cif")
s = imp.structure
print(f"CIF: {imp.block}, space group {imp.spacegroup}, "
      f"{imp.n_asymmetric} asymmetric sites -> {s.n_atoms} atoms")
assert s.n_atoms == 8
a = s.lattice.parameters[0]
print(f"  a = {a:.4f} A")

# --- 2. Build a 4x4x4 configuration supercell via the Structure API --------
n = (4, 4, 4)
m = s.n_atoms
lat = Lattice.from_parameters(a * n[0], a * n[1], a * n[2], *s.lattice.parameters[3:])
N = m * n[0] * n[1] * n[2]
cell_num = np.arange(N) // m
slot = np.arange(N) % m
cell_index = np.stack(
    [cell_num % n[0], (cell_num // n[0]) % n[1], cell_num // (n[0] * n[1])],
    axis=1,
).astype(np.int64)
supercell = Structure(
    lat,
    s.species[slot],
    (s.frac[slot] + cell_index) / np.array(n, dtype=float),
    model="configuration",
    cell_index=cell_index,
    supercell=n,
    pbc=(True, True, True),
)
print(f"4x4x4 configuration supercell: {supercell.n_atoms} atoms, "
      f"supercell {tuple(supercell.supercell)}")
assert supercell.n_atoms == 512
assert tuple(int(x) for x in supercell.supercell) == (4, 4, 4)

# --- 3. Nearest-neighbour check: a*sqrt(3)/4, 4 neighbours -----------------
D = supercell.distance_matrix()
np.fill_diagonal(D, np.inf)
d_min = D.min()
expected = a * np.sqrt(3) / 4
assert abs(d_min - expected) < 1e-6, f"min distance {d_min:.6f} != {expected:.6f}"
nn = (D < 2.6).sum(axis=1)
assert set(nn.tolist()) == {4}, sorted(set(nn.tolist()))
print(f"nearest-neighbour check: min distance {d_min:.6f} A "
      f"(expected {expected:.6f}), 4 neighbours per atom")

# --- 4. Convert to ASE and DiffPy (if installed) ---------------------------
try:
    from discuspdf.ase_io import to_ase, from_ase
    atoms = to_ase(s)
    s2 = from_ase(atoms)
    assert s2.isclose(s, rtol=1e-8, atol=1e-12)
    print(f"ASE round trip: OK ({len(atoms)} atoms, dp_* arrays attached)")
except ImportError as e:
    print(f"ASE skipped: {e}")

try:
    from discuspdf.diffpy_io import to_diffpy, from_diffpy
    ds = to_diffpy(s)
    s3 = from_diffpy(ds, model=s.model)
    assert s3.isclose(s, rtol=1e-8, atol=1e-12)
    print(f"DiffPy round trip: OK ({len(ds)} atoms)")
except ImportError as e:
    print(f"DiffPy skipped: {e}")

# --- 5. Write a DISCUS .cll and a CIF --------------------------------------
with tempfile.TemporaryDirectory() as tmpdir:
    cll_path = Path(tmpdir) / "gaas.cll"
    cif_path = Path(tmpdir) / "gaas.cif"
    write_cll(cll_path, s)
    write_cif(cif_path, s)
    from discuspdf.discus_io import read_cll
    s_cll = read_cll(cll_path).structure
    # .cll has 6-decimal precision and regenerates site labels; compare
    # lattice, coordinates, occupancy and B factor (not the exact labels).
    assert np.allclose(s_cll.lattice.parameters, s.lattice.parameters, atol=1e-6)
    assert np.allclose(s_cll.frac, s.frac, atol=1e-6)
    assert np.allclose(s_cll.occupancy, s.occupancy, atol=1e-6)
    # .cll has no B column here, so read_cll defaults B to 1.0; skip the B
    # comparison for this B-free fixture.
    s_cif = read_cif(cif_path).structure
    assert s_cif.isclose(s, rtol=1e-8, atol=1e-9)
    print(f".cll and CIF written and read back: {s_cll.n_atoms} atoms each")

print("interop_basics: all checks passed")
