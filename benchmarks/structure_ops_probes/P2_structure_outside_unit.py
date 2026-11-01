"""P2: Structure with positions outside [0,1)."""
import numpy as np
from discuspdf import Lattice, Structure
from discuspdf.cif_io import read_cif

# Build the GaAs 8-atom cell via read_cif
cif = read_cif("tests/fixtures/interop/gaas_cod_9008845.cif")
s = cif.structure
print(f"Original: n_atoms={s.n_atoms}, frac range=[{s.frac.min():.4f}, {s.frac.max():.4f}]")

# Try with_frac with values outside [0,1)
frac_shifted = s.frac + np.array([1.2, -0.3, 0.05])
print(f"frac_shifted range=[{frac_shifted.min():.4f}, {frac_shifted.max():.4f}]")
try:
    s2 = s.with_frac(frac_shifted)
    print(f"with_frac accepted: n_atoms={s2.n_atoms}, frac range=[{s2.frac.min():.4f}, {s2.frac.max():.4f}]")
except Exception as e:
    print(f"with_frac raised: {type(e).__name__}: {e}")

# Try a single atom moved to -0.01
frac_neg = s.frac.copy()
frac_neg[0] = -0.01
print(f"frac_neg[0] = {frac_neg[0]}")
try:
    s3 = s.with_frac(frac_neg)
    print(f"with_frac(-0.01) accepted: n_atoms={s3.n_atoms}, frac[0]={s3.frac[0]}")
except Exception as e:
    print(f"with_frac(-0.01) raised: {type(e).__name__}: {e}")

# If accepted, check wrapped()
if 's2' in dir():
    try:
        w = s2.wrapped()
        print(f"wrapped() frac range=[{w.frac.min():.4f}, {w.frac.max():.4f}]")
    except Exception as e:
        print(f"wrapped() raised: {type(e).__name__}: {e}")

# distance_matrix and overlap check with shifted atoms
if 's2' in dir():
    try:
        dm = s2.distance_matrix()
        print(f"distance_matrix shape={dm.shape}, min={dm.min():.4f}")
    except Exception as e:
        print(f"distance_matrix raised: {type(e).__name__}: {e}")
