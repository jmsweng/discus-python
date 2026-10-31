"""P5: DISCUS I/O ordering for 4x4x4 GaAs supercell."""
import numpy as np
from discuspdf import Lattice, Structure
from discuspdf.discus_io import format_stru, parse_stru, format_cll

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

d = _gaas_supercell(4)
s = Structure(d["lattice"], d["species"], d["frac"], model="configuration",
              site_label=d["site_label"], cell_index=d["cell_index"], supercell=d["supercell"])
print(f"4x4x4 supercell: n_atoms={s.n_atoms}, supercell={s.supercell}")

# The helper uses np.meshgrid(indexing="ij") which gives cells in order:
# (0,0,0), (0,0,1), ..., (0,1,0), ... (iz fastest, ix slowest)
# DISCUS expects: ix fastest, iy, iz slowest: (0,0,0), (1,0,0), ..., (0,1,0), ...
# So the rows are NOT in DISCUS order.

print(f"cell_index[:8] = {s.cell_index[:8].tolist()}")
print(f"cell_index[8:16] = {s.cell_index[8:16].tolist()}")

# Try format_stru — should raise DiscusFormatError
try:
    text = format_stru(s)
    print(f"format_stru accepted (unexpected)")
    s_rt = parse_stru(text)
    print(f"Round trip equal: {s == s_rt}")
except Exception as e:
    print(f"format_stru raised: {type(e).__name__}: {e}")

# Note: Structure.select(rows=...) does NOT reorder rows — it uses a boolean
# mask and returns rows in ascending order. So select cannot fix the ordering.
# The only way to get DISCUS-order rows is to build the Structure with rows
# already in DISCUS order (ix fastest).

# Build a DISCUS-ordered version manually:
N = 4
a = 5.6537
fcc = np.array([[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]])
base = np.vstack([fcc, fcc + 0.25])

# DISCUS order: ix fastest, iy, iz slowest
cells_disco = np.array([(ix, iy, iz) for iz in range(N) for iy in range(N) for ix in range(N)])
frac_disco = ((base[None, :, :] + cells_disco[:, None, :]) / N).reshape(-1, 3)
labels = [f"Ga{i}" for i in range(1, 5)] + [f"As{i}" for i in range(1, 5)]
s_disco = Structure(
    Lattice.from_parameters(a * N, a * N, a * N, 90.0, 90.0, 90.0),
    np.tile(["Ga"] * 4 + ["As"] * 4, len(cells_disco)).astype(object),
    frac_disco,
    model="configuration",
    site_label=np.tile(labels, len(cells_disco)).astype(object),
    cell_index=np.repeat(cells_disco, 8, axis=0),
    supercell=(N, N, N),
)
print(f"\nDISCUS-ordered structure: n_atoms={s_disco.n_atoms}")
print(f"cell_index[:8] = {s_disco.cell_index[:8].tolist()}")
print(f"cell_index[8:16] = {s_disco.cell_index[8:16].tolist()}")

try:
    text = format_stru(s_disco)
    print(f"format_stru accepted (DISCUS-ordered)")
    ds_rt = parse_stru(text)
    s_rt = ds_rt.structure
    print(f"Round trip equal: {s_disco == s_rt}")
    if s_disco != s_rt:
        print(f"  frac diff max: {np.abs(s_disco.frac - s_rt.frac).max():.2e}")
        print(f"  species equal: {np.array_equal(s_disco.species, s_rt.species)}")
        print(f"  site_label equal: {np.array_equal(s_disco.site_label, s_rt.site_label)}")
        print(f"  cell_index equal: {np.array_equal(s_disco.cell_index, s_rt.cell_index)}")
        print(f"  supercell equal: {s_disco.supercell == s_rt.supercell}")
        print(f"  atom_id equal: {np.array_equal(s_disco.atom_id, s_rt.atom_id)}")
except Exception as e:
    print(f"format_stru raised (DISCUS-ordered): {type(e).__name__}: {e}")

# Displace 10 random atoms by 0.05 Å and 0.4 Å using with_cartesian
if 's_disco' in dir():
    rng = np.random.default_rng(42)
    idx = rng.choice(s_disco.n_atoms, size=10, replace=False)
    cart = s_disco.cartesian.copy()
    cart[idx] += 0.05  # 0.05 Å displacement
    s_disp = s_disco.with_cartesian(cart)
    try:
        text2 = format_stru(s_disp)
        ds_rt2 = parse_stru(text2)
        s_rt2 = ds_rt2.structure
        print(f"Displaced 0.05 Å round trip equal: {s_disp == s_rt2}")
    except Exception as e:
        print(f"Displaced 0.05 Å parse_stru raised: {type(e).__name__}: {e}")

    cart3 = s_disco.cartesian.copy()
    cart3[idx] += 0.4  # 0.4 Å displacement
    s_disp3 = s_disco.with_cartesian(cart3)
    try:
        text3 = format_stru(s_disp3)
        ds_rt3 = parse_stru(text3)
        s_rt3 = ds_rt3.structure
        print(f"Displaced 0.4 Å round trip equal: {s_disp3 == s_rt3}")
    except Exception as e:
        print(f"Displaced 0.4 Å parse_stru raised: {type(e).__name__}: {e}")

# format_cll with supercell != (1,1,1)
try:
    cll_text = format_cll(s)
    print(f"\nformat_cll accepted supercell={s.supercell}, output length={len(cll_text)}")
except Exception as e:
    print(f"\nformat_cll raised: {type(e).__name__}: {e}")
