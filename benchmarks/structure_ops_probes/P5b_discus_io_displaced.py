"""P5b: DISCUS I/O round-trip with displaced atoms (Gaussian noise)."""
import numpy as np
from discuspdf import Lattice, Structure
from discuspdf.discus_io import format_stru, parse_stru

def _gaas_disco(N):
    """GaAs zinc-blende N x N x N supercell in DISCUS order (ix fastest)."""
    a = 5.6537
    fcc = np.array([[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]])
    base = np.vstack([fcc, fcc + 0.25])
    cells = np.array([(ix, iy, iz) for iz in range(N) for iy in range(N) for ix in range(N)])
    frac = ((base[None, :, :] + cells[:, None, :]) / N).reshape(-1, 3)
    labels = [f"Ga{i}" for i in range(1, 5)] + [f"As{i}" for i in range(1, 5)]
    return Structure(
        Lattice.from_parameters(a * N, a * N, a * N, 90.0, 90.0, 90.0),
        np.tile(["Ga"] * 4 + ["As"] * 4, len(cells)).astype(object),
        frac,
        model="configuration",
        site_label=np.tile(labels, len(cells)).astype(object),
        cell_index=np.repeat(cells, 8, axis=0),
        supercell=(N, N, N),
    )

s = _gaas_disco(4)
print(f"4x4x4 DISCUS-ordered GaAs: n_atoms={s.n_atoms}")

# Displace ALL atoms with independent Gaussian Cartesian noise
for sigma in [0.05, 0.4, 1.0, 2.0]:
    rng = np.random.default_rng(1)
    noise = rng.normal(0, sigma, size=(s.n_atoms, 3))
    cart = s.cartesian + noise
    s_disp = s.with_cartesian(cart)
    try:
        text = format_stru(s_disp)
        ds = parse_stru(text)
        s_rt = ds.structure
        max_diff = np.abs(s_disp.frac - s_rt.frac).max()
        print(f"sigma={sigma}: round trip ok, max |df| = {max_diff:.2e}")
    except Exception as e:
        print(f"sigma={sigma}: {type(e).__name__}: {e}")

# One run where only ONE atom is moved by 3.0 Å
cart_one = s.cartesian.copy()
cart_one[0] += 3.0
s_one = s.with_cartesian(cart_one)
try:
    text = format_stru(s_one)
    ds = parse_stru(text)
    s_rt = ds.structure
    max_diff = np.abs(s_one.frac - s_rt.frac).max()
    print(f"single atom +3.0 A: round trip ok, max |df| = {max_diff:.2e}")
except Exception as e:
    print(f"single atom +3.0 Å: {type(e).__name__}: {e}")
