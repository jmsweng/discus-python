"""Independent oracle for the pre-registered Gate 0 v2 experiment (DISCUS_SUITE 6.23.01).

Nothing here imports the package under test (``src/discuspdf``). numpy only.

Physical system
---------------
FCC lattice, conventional cubic cell a = 3.75 Å, 4 sites per cell at
(0,0,0), (½,½,0), (½,0,½), (0,½,½); every site Cu with p_Cu = 0.75,
Au with p_Au = p = 0.25.  Crystal: repeat(average, (8, 8, 8)) →
N_c = 512 cells, N = 2048 sites.

Reciprocal grid: dyadic, corner (0,0,0), step 1/8 along each reciprocal
axis, 25 points per axis (h, k, l ∈ {0, 1/8, …, 3}).

Analytic ensemble expectation (the independent oracle):
    E|D(h)|² = p(1-p) |g(s)|² · Σ_j [N_c - |L(h)|²/N_c]
    g(s) = f_Au(s) e^{-B_Au s²} - f_Cu(s) e^{-B_Cu s²}
    L(h) = Σ_R e^{2πi h·R},  s = |h|/(2a)

On this grid L(h) = 0 unless h is integer, so
    E = 4 N_c p(1-p)|g|²  off integer h,  E = 0  at integer h.

Two precision modes:
    "physics": DABAX W-K coefficients, f at the exact s.
    "engine" : DISCUS's own W-K coefficient table, f at s rounded to
               0.001 Å⁻¹ with Fortran NINT (half away from zero).

The realization-exact oracle (D_k(h) = Σ_j g_j Σ_R δc_{jR,k} e^{2πi h·(R+x_j)}
from the actual k-th structure file) is also computed, to separate engine
precision from ensemble statistics.

Coefficients: wk_coefficients.json (written by compare_form_factors.py in gate2_prereg).
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np

# ----------------------------------------------------------------------------- import guard
def _assert_no_discuspdf() -> None:
    """Assert at import that no ``discuspdf`` module is loaded."""
    loaded = [m for m in sys.modules if m == "discuspdf" or m.startswith("discuspdf.")]
    if loaded:
        raise ImportError(
            f"oracle.py must not import discuspdf; found loaded modules: {loaded}"
        )

_assert_no_discuspdf()

HERE = Path(__file__).resolve().parent
WK_PATH = HERE.parent / "gate2_prereg" / "wk_coefficients.json"

# ----------------------------------------------------------------------------- constants
A_CELL = 3.75
N_CELL = 8
N_C = N_CELL ** 3          # 512 cells
N_SITES = 4 * N_C          # 2048 sites
P_AU = 0.25
P_CU = 0.75
B_CU_A = 0.0
B_AU_A = 0.0               # Arm A
B_CU_B = 0.5              # Arm B
B_AU_B = 0.8
CFINC = 0.001
SPECIES = ("CU", "AU")

# FCC sites in fractional coordinates of the conventional cell
X_SITES = np.array([
    [0.0, 0.0, 0.0],
    [0.5, 0.5, 0.0],
    [0.5, 0.0, 0.5],
    [0.0, 0.5, 0.5],
])

# Reciprocal grid: dyadic, corner (0,0,0), step 1/8, 25 points per axis
H_AXIS = np.arange(25) / 8.0   # 0, 1/8, ..., 3

# ----------------------------------------------------------------------------- form factors
def load_wk(source: str = "dabax") -> dict:
    data = json.loads(WK_PATH.read_text())
    return {el: data["elements"][el][source] for el in SPECIES}


def coefficients_for(mode: str) -> dict:
    return load_wk({"physics": "dabax", "engine": "discus"}[mode])


def wk_f(coef: dict, s) -> np.ndarray:
    s2 = np.asarray(s, dtype=float) ** 2
    return coef["c"] + sum(a * np.exp(-b * s2) for a, b in zip(coef["a"], coef["b"]))


def s_of_grid(mode: str) -> np.ndarray:
    """sin(theta)/lambda on the full H grid, (25, 25, 25); 'engine' mimics DISCUS rounding."""
    h2 = H_AXIS ** 2
    q2 = (h2[:, None, None] + h2[None, :, None] + h2[None, None, :]) * (1.0 / A_CELL ** 2) / 4.0
    s = np.sqrt(q2)
    if mode == "physics":
        return s
    if mode == "engine":
        return np.floor(s * (1.0 / CFINC) + 0.5) * CFINC     # Fortran NINT for s >= 0
    raise ValueError(mode)


# ----------------------------------------------------------------------------- analytic expectation
def g_squared(mode: str, s: np.ndarray | None = None,
             b_cu: float = B_CU_A, b_au: float = B_AU_A) -> np.ndarray:
    """|g(s)|² = |f_Au(s) e^{-B_Au s²} - f_Cu(s) e^{-B_Cu s²}|²."""
    wk = coefficients_for(mode)
    s = s_of_grid(mode) if s is None else s
    f_au = wk_f(wk["AU"], s)
    f_cu = wk_f(wk["CU"], s)
    g = f_au * np.exp(-b_au * s ** 2) - f_cu * np.exp(-b_cu * s ** 2)
    return np.abs(g) ** 2


def lattice_sum_sq(h: np.ndarray) -> np.ndarray:
    """|L(h)|² = |Σ_R e^{2πi h·R}|² for R = 0..7 in each axis (8³ cells).

    L(h) = Π_axis l(h_axis) where l(h) = Σ_{r=0}^{7} e^{2πi h r}.
    For h = m/8: l = 0 if m not divisible by 8, else 8.
    |L|² = 8^6 if h integer, else 0.
    """
    is_int = np.all(np.isclose(h, np.rint(h)), axis=-1)
    return np.where(is_int, float(N_C ** 2), 0.0)


def expectation(mode: str, s: np.ndarray | None = None,
               b_cu: float = B_CU_A, b_au: float = B_AU_A) -> np.ndarray:
    """E|D(h)|² = p(1-p)|g(s)|² · Σ_j [N_c - |L(h)|²/N_c] on the full grid.

    On this grid L(h) = 0 off integer h, so
        E = 4 N_c p(1-p)|g|²  off integer h,  E = 0  at integer h.
    """
    g2 = g_squared(mode, s, b_cu, b_au)
    h = np.stack(np.meshgrid(H_AXIS, H_AXIS, H_AXIS, indexing="ij"), axis=-1)
    l2 = lattice_sum_sq(h)
    # Σ_j [N_c - |L|²/N_c] = 4 * (N_c - |L|²/N_c)
    # Off integer h: |L|² = 0 → 4 N_c
    # At integer h: |L|² = N_c² → 0
    sum_j = 4.0 * (N_C - l2 / N_C)
    return P_AU * (1.0 - P_AU) * g2 * sum_j


# ----------------------------------------------------------------------------- independent evaluation set S
def independent_set_mask() -> np.ndarray:
    """Boolean mask for the independent evaluation set S.

    h in the first zone [0,1)³ of the grid, q ≠ 0, 2q not integer,
    and keep only one of each pair {q, −q mod 1}.
    """
    h = np.stack(np.meshgrid(H_AXIS, H_AXIS, H_AXIS, indexing="ij"), axis=-1)
    # First zone: h in [0, 1) on each axis
    in_zone = np.all(h < 1.0, axis=-1)
    # q ≠ 0
    q_nonzero = ~np.all(h == 0.0, axis=-1)
    # 2q not integer: for h = m/8, 2h = m/4; integer iff m divisible by 4
    # h = m/8, 2h = m/4; m/4 integer iff m % 4 == 0
    m = (H_AXIS * 8).astype(int)
    two_int = (m[:, None, None] % 4 == 0) & (m[None, :, None] % 4 == 0) & (m[None, None, :] % 4 == 0)
    # Keep one of each {q, -q mod 1} pair: keep q where q_x < 0.5 OR (q_x == 0.5 AND q_y < 0.5) OR ...
    # Simpler: keep q where the first nonzero component is positive (in [0,1))
    # Since we're already in [0,1), -q mod 1 = 1-q. Keep q if q != 1-q, i.e. 2q != 0 mod 1.
    # 2q != 0 mod 1 means 2q not integer, which we already require.
    # For q in [0,1) with 2q not integer: q and -q mod 1 = 1-q are distinct.
    # Keep q if q < 1-q in lexicographic order, i.e. q is "closer to 0".
    # Simplest: keep q where q_x < 0.5 OR (q_x == 0.5 AND q_y < 0.5) OR ...
    # But since 2q not integer, no component is exactly 0.5 (that would make 2q have that component = 1, integer).
    # So: keep q where q_x < 0.5, OR (q_x > 0.5 AND ...) — no, just keep the half-space.
    # For a pair {q, 1-q}: exactly one has each component < 0.5 (since no component is 0.5).
    # Keep q if q_x < 0.5 (the "positive" half).
    keep = h[..., 0] < 0.5
    return in_zone & q_nonzero & ~two_int & keep


def n_independent_set() -> int:
    return int(independent_set_mask().sum())


# ----------------------------------------------------------------------------- realization-exact oracle
def parse_structure(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Parse a .stru file written by format_stru.

    Returns (species, cell_index, b_iso) where:
      species: (N,) string array
      cell_index: (N, 3) int array (unit-cell indices)
      b_iso: (N,) float array
    """
    species, ci, bi = [], [], []
    for line in Path(path).read_text().splitlines():
        t = line.replace(",", " ").split()
        if t and t[0] in SPECIES:
            species.append(t[0])
            x, y, z = float(t[1]), float(t[2]), float(t[3])
            b = float(t[4])
            # x = ox + x_site where x_site ∈ {0, 0.5}
            # Determine cell index from fractional part
            ox = int(np.round(x)) if abs(x - np.round(x)) < 0.1 else int(np.floor(x))
            oy = int(np.round(y)) if abs(y - np.round(y)) < 0.1 else int(np.floor(y))
            oz = int(np.round(z)) if abs(z - np.round(z)) < 0.1 else int(np.floor(z))
            ci.append((ox, oy, oz))
            bi.append(b)
    return (np.array(species, dtype=object),
            np.array(ci, dtype=np.int64),
            np.array(bi, dtype=float))


def _occupancy_grid(species: np.ndarray, ci: np.ndarray) -> np.ndarray:
    """Build the 8³ occupancy grid: occ[ox,oy,oz,j] = 1 if site j in cell (ox,oy,oz) is Au."""
    n = len(species)
    occ = np.zeros((N_CELL, N_CELL, N_CELL, 4), dtype=float)
    is_au = (species == "AU")
    if not is_au.any():
        return occ
    idx = ci[is_au]
    # Determine which of the 4 FCC sites each Au atom belongs to
    # x_file = ox + x_site, so x_site = x_file - ox
    # We need to recover x_site from the file coordinates
    # Re-read the file to get exact coordinates
    # Actually, we can compute: for each Au atom, x_site = x - ox
    # But we only have ci here. Let's re-parse to get coordinates.
    # Simpler: use the fact that x_site ∈ {0, 0.5} and determine from the fractional part
    # We already have ci, so let's re-read to get the fractional parts
    return occ


def _occupancy_grid_from_file(path: Path) -> np.ndarray:
    """Build the 8³ × 4 occupancy grid from a .stru file.

    occ[ox, oy, oz, j] = 1 if site j in cell (ox, oy, oz) is Au in this realization.
    """
    occ = np.zeros((N_CELL, N_CELL, N_CELL, 4), dtype=float)
    for line in Path(path).read_text().splitlines():
        t = line.replace(",", " ").split()
        if not t or t[0] not in SPECIES:
            continue
        sp = t[0]
        x, y, z = float(t[1]), float(t[2]), float(t[3])
        # Determine cell index and site
        # x = ox + x_site, x_site ∈ {0, 0.5}
        ox = int(np.round(x)) if abs(x - np.round(x)) < 0.1 else int(np.floor(x))
        oy = int(np.round(y)) if abs(y - np.round(y)) < 0.1 else int(np.floor(y))
        oz = int(np.round(z)) if abs(z - np.round(z)) < 0.1 else int(np.floor(z))
        # Determine which FCC site: x_site = x - ox
        x_site = x - ox
        y_site = y - oy
        z_site = z - oz
        # Match to one of the 4 FCC sites
        for j, (xs, ys, zs) in enumerate(X_SITES):
            if abs(x_site - xs) < 0.1 and abs(y_site - ys) < 0.1 and abs(z_site - zs) < 0.1:
                if sp == "AU":
                    occ[ox, oy, oz, j] = 1.0
                break
    return occ


def total_direct(path: Path, mode: str = "engine") -> np.ndarray:
    """Compute I(h) = |F(h)|² directly from the .stru file, with per-atom
    f_el(s) and B from the file (no averaging, no site-level assumption).

    F(h) = Σ_n f_n e^{-B_n s²} e^{2πi h·x_n}
    I(h) = |F(h)|²

    This is the ground-truth total intensity for a single realization.
    """
    s = s_of_grid(mode)
    wk = coefficients_for(mode)
    f_au = wk_f(wk["AU"], s)
    f_cu = wk_f(wk["CU"], s)

    # Parse all atoms from the .stru file (skip header lines)
    atoms = []
    for line in Path(path).read_text().splitlines():
        t = line.replace(",", " ").split()
        if not t or t[0] not in SPECIES:
            continue
        sp = t[0]
        x, y, z = float(t[1]), float(t[2]), float(t[3])
        B = float(t[4])
        atoms.append((sp, x, y, z, B))

    # Vectorized computation: process all atoms at once
    # F = Σ_n f_n e^{-B_n s²} e^{2πi h·x_n}
    n = len(atoms)
    if n == 0:
        return np.zeros(s.shape)

    # Build numpy arrays for all atoms
    sp_arr = np.array([a[0] for a in atoms])
    x_arr = np.array([a[1] for a in atoms])
    y_arr = np.array([a[2] for a in atoms])
    z_arr = np.array([a[3] for a in atoms])
    B_arr = np.array([a[4] for a in atoms])

    # Process in chunks to limit memory: (chunk, 25, 25, 25)
    F = np.zeros(s.shape, complex)
    chunk = 256
    for i in range(0, n, chunk):
        sl = slice(i, i + chunk)
        x_c = x_arr[sl]  # (chunk,)
        y_c = y_arr[sl]  # (chunk,)
        z_c = z_arr[sl]  # (chunk,)
        B_c = B_arr[sl]  # (chunk,)
        sp_c = sp_arr[sl]  # (chunk,)

        # Phase: (chunk, 25, 25, 25)
        phase = np.exp(2j * np.pi * (
            H_AXIS[None, :, None, None] * x_c[:, None, None, None] +
            H_AXIS[None, None, :, None] * y_c[:, None, None, None] +
            H_AXIS[None, None, None, :] * z_c[:, None, None, None]
        ))
        # Amplitude: (chunk, 25, 25, 25)
        # f_el(s) is (25,25,25), B_c is (chunk,)
        # amp[n, ix, iy, iz] = f_el[ix, iy, iz] * exp(-B_n * s[ix, iy, iz]²)
        is_au = (sp_c == "AU")[:, None, None, None]  # (chunk, 1, 1, 1)
        f_c = np.where(is_au,
                       np.broadcast_to(f_au, (len(sp_c), 25, 25, 25)),
                       np.broadcast_to(f_cu, (len(sp_c), 25, 25, 25)))
        amp = f_c * np.exp(-B_c[:, None, None, None] * s[None, :, :, :] ** 2)
        F += np.sum(amp * phase, axis=0)
    return np.abs(F) ** 2


def realization_exact(path: Path, mode: str = "engine",
                     b_cu: float = B_CU_A, b_au: float = B_AU_A) -> tuple[np.ndarray, np.ndarray]:
    """Compute the realization-exact D(h) and I(h) on the full grid.

    D_k(h) = Σ_j g_j Σ_R δc_{jR,k} e^{2πi h·(R+x_j)}
    where δc_{jR,k} = 1 if site j in cell R is Au in realization k, else 0.

    Returns (D, I) where:
      D: (25, 25, 25) complex array, the diffuse structure factor
      I: (25, 25, 25) float array, the total intensity |F|²
    """
    s = s_of_grid(mode)
    wk = coefficients_for(mode)
    f_au = wk_f(wk["AU"], s)
    f_cu = wk_f(wk["CU"], s)

    # Build the occupancy grid: occ[ox,oy,oz,j] = 1 if Au
    occ = _occupancy_grid_from_file(path)

    # Subtract the per-site mean over cells: δc_{jR} = c_{jR} - c̄_j
    c_mean = occ.mean(axis=(0, 1, 2), keepdims=True)  # (1, 1, 1, 4)
    delta_c = occ - c_mean  # (8, 8, 8, 4)

    # Compute D via FFT:
    # D(h) = Σ_j g_j Σ_R δc_{jR} e^{2πi h·R} e^{2πi h·x_j}
    #      = Σ_j g_j e^{2πi h·x_j} · DFT_j(h)
    # where DFT_j(h) = Σ_R δc_{jR} e^{2πi h·R}

    # For h = m/8 (m = 0..24), the DFT index is m mod 8 (since the grid is 8³)
    # numpy ifftn gives: ifft(c)[k] = (1/N³) Σ_n c[n] e^{2πi k n / N}
    # So Σ_n c[n] e^{2πi h n} = N³ · ifft(c)[m] where h = m/8

    # For the 25-point grid: h = m/8 for m = 0..24
    # The DFT index is m mod 8 (periodic with period 8)
    m = np.arange(25)
    k_idx = m % N_CELL  # 0..7, 0..7, 0..7, 0..1

    D = np.zeros(s.shape, complex)
    for j in range(4):
        delta_j = delta_c[:, :, :, j]
        # DFT: Σ_R δc_{jR} e^{2πi h·R} for h = m/8
        # numpy ifftn uses e^{-2πi k n / N}, so we need ifft
        # Σ_R δc e^{2πi h·R} = N³ · ifft(δc) at the right index
        dft_full = np.fft.ifftn(delta_j) * (N_CELL ** 3)  # (8, 8, 8)
        # Evaluate at the 25 grid points: index m mod 8
        dft = dft_full[np.ix_(k_idx, k_idx, k_idx)]  # (25, 25, 25)
        # Phase factor: e^{2πi h·x_j}
        xj = X_SITES[j]
        phase = np.exp(2j * np.pi * (
            H_AXIS[:, None, None] * xj[0] +
            H_AXIS[None, :, None] * xj[1] +
            H_AXIS[None, None, :] * xj[2]
        ))
        g_j = f_au * np.exp(-b_au * s ** 2) - f_cu * np.exp(-b_cu * s ** 2)
        D += g_j * dft * phase

    # Total intensity: direct sum over atoms (per-atom f and per-atom B from the file)
    I = total_direct(path, mode)

    return D, I


def realization_exact_diffuse(path: Path, mode: str = "engine",
                              b_cu: float = B_CU_A, b_au: float = B_AU_A) -> np.ndarray:
    """|D(h)|² from the realization-exact oracle."""
    D, _ = realization_exact(path, mode, b_cu, b_au)
    return np.abs(D) ** 2


# ----------------------------------------------------------------------------- utilities
def sha256_file(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def sha256_array(a: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(a, dtype=np.float64).tobytes()).hexdigest()


def rel_l2(a, b) -> float:
    return float(np.linalg.norm(np.asarray(a) - np.asarray(b)) / max(np.linalg.norm(b), 1e-300))
