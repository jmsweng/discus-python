"""Independent oracle for the pre-registered Gate 2 experiment (DISCUS_SUITE 6.23.01).

Nothing here imports the package under test (``src/discuspdf``). numpy only.

Physical system
---------------
Simple-cubic lattice, a = 4 Å, one atom per primitive cell, 16 x 16 x 16 cells (N = 4096),
binary Cu/Au, B = 0. Occupation variable sigma_j = +1 (Au) / -1 (Cu), sum(sigma) = 0 exactly.

Reciprocal grid: H = -4 + k/16, k = 0..127 on each axis (basecell rlu), i.e. commensurate with the
16-cell crystal. Because every atom sits on an integer lattice vector n_j, exp(2 pi i H.n_j) is
periodic in k with period 16, so the finite-crystal sum computed by DISCUS equals the periodic
(circular) sum exactly. Open vs periodic boundary conventions therefore cannot differ.

Forward model (verified on the PR #8 engine runs to 6e-16 / 7e-15 relative L2 for total intensity):
    F(H)     = sum_j f_j(s) exp(+2 pi i H.x_j),   s = sin(theta)/lambda = |H| / (2a)
    I_tot    = |F|^2
DISCUS 6.23.01 `set aver,100` (four_aver in discus/prog/fourier_sup.f90): <F>(H) = (mean over unit
cells of the cell structure factor, positions reduced into the first cell) x (lattice interference
function L(H) = sum_cells exp(2 pi i H.n)). With one site per cell this is fbar(s) L(H), fbar = mean f.
    I_diff   = |F - fbar L|^2 = g(s)^2 S(H),  g = (f_Au - f_Cu)/2,  S(H) = |sum_j sigma_j e^{2 pi i H.n_j}|^2
Two precision modes:
    "physics": independent reference. W-K coefficients from the DABAX transcription, f at the exact s.
    "engine" : emulates the engine's inputs exactly. DISCUS 6.23.01's own W-K coefficient table
               (lib_f90/element_data_mod.f90; agrees with DABAX to <= 4.3e-8 relative in f), f at s
               rounded to the nearest 0.001 Å^-1 with Fortran NINT (half away from zero), as DISCUS
               does (four_formtab / CFINC = 0.001, istl = nint(sqrt(q2)/CFINC)).
The physics-vs-engine difference (eps_tab, ~1e-3) is computed before any engine run and sets the
tolerance of every "physics" criterion; "engine" criteria test the algorithm with no approximation.
Coefficients: wk_coefficients.json (written by compare_form_factors.py).
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
N_CELL = 16
A_CELL = 4.0
H_MAX = 4.0
DH = 1.0 / N_CELL
N_H = int(round(2 * H_MAX / DH))           # 128
H_AXIS = -H_MAX + DH * np.arange(N_H)      # -4 .. 3.9375, exact in binary
CFINC = 0.001
SPECIES = ("CU", "AU")                     # sigma = -1, +1
SEEDS = {"block": 20261028, "random": 20261029}

# Predeclared real-space evaluation points, in lattice vectors (cells)
TAU_POINTS = [(0, 0, 0), (1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1),
              (1, 1, 0), (1, -1, 0), (-1, 1, 0), (-1, -1, 0), (1, 0, 1), (0, 1, 1),
              (2, 0, 0), (0, 2, 0), (0, 0, 2)]

ENSEMBLE = {   # analytic ensemble values of C(tau)/N for the block model
    "block": {(0, 0, 0): 1.0, (1, 0, 0): -0.5, (-1, 0, 0): -0.5, (0, 1, 0): 0.5, (0, -1, 0): 0.5,
              (0, 0, 1): 0.0, (0, 0, -1): 0.0, (1, 1, 0): -0.25, (1, -1, 0): -0.25,
              (-1, 1, 0): -0.25, (-1, -1, 0): -0.25, (1, 0, 1): 0.0, (0, 1, 1): 0.0,
              (2, 0, 0): 0.0, (0, 2, 0): 0.0, (0, 0, 2): 0.0},
}


# ----------------------------------------------------------------------------- form factors
def load_wk(source: str = "dabax") -> dict:
    data = json.loads((HERE / "wk_coefficients.json").read_text())
    return {el: data["elements"][el][source] for el in SPECIES}


def coefficients_for(mode: str) -> dict:
    return load_wk({"physics": "dabax", "engine": "discus"}[mode])


def wk_f(coef: dict, s) -> np.ndarray:
    s2 = np.asarray(s, dtype=float) ** 2
    return coef["c"] + sum(a * np.exp(-b * s2) for a, b in zip(coef["a"], coef["b"]))


def s_of_grid(mode: str) -> np.ndarray:
    """sin(theta)/lambda on the full H grid, (N_H, N_H, N_H); 'engine' mimics DISCUS rounding."""
    h2 = H_AXIS ** 2
    q2 = (h2[:, None, None] + h2[None, :, None] + h2[None, None, :]) * (1.0 / A_CELL ** 2) / 4.0
    s = np.sqrt(q2)
    if mode == "physics":
        return s
    if mode == "engine":
        return np.floor(s * (1.0 / CFINC) + 0.5) * CFINC     # Fortran NINT for s >= 0
    raise ValueError(mode)


# ----------------------------------------------------------------------------- structures
def build_sigma(model: str, seed: int | None = None) -> np.ndarray:
    """sigma[ix, iy, iz] in {-1, +1}, sum exactly zero."""
    rng = np.random.default_rng(SEEDS[model] if seed is None else seed)
    n = N_CELL
    x = np.arange(n)
    if model == "block":
        # independent 2x2 (x,y) blocks in each z layer: x-neighbours unlike, y-neighbours like
        s_blk = rng.choice(np.array([-1, 1]), size=(n // 2, n // 2, n))
        sig = s_blk[x[:, None, None] // 2, x[None, :, None] // 2, x[None, None, :]] \
            * ((-1) ** x)[:, None, None]
    elif model == "random":
        flat = np.array([1] * (n ** 3 // 2) + [-1] * (n ** 3 // 2))
        sig = rng.permutation(flat).reshape(n, n, n)
    else:
        raise ValueError(model)
    sig = sig.astype(np.int8)
    assert sig.sum() == 0
    return sig


def atom_list(sig: np.ndarray):
    """Atoms in DISCUS celltoindex order (x fastest, then y, then z), one site per cell."""
    names, pos = [], []
    n = sig.shape[0]
    for iz in range(n):
        for iy in range(n):
            for ix in range(n):
                names.append(SPECIES[1] if sig[ix, iy, iz] > 0 else SPECIES[0])
                pos.append((ix, iy, iz))
    return names, np.asarray(pos, dtype=float)


def stru_text(sig: np.ndarray, title: str) -> str:
    """DISCUS structure file in the layout DISCUS itself writes (as the archived C1 file)."""
    names, pos = atom_list(sig)
    n = sig.shape[0]
    lines = [f"title {title}",
             "spcgr P1",
             f"cell  {A_CELL:10.6f},{A_CELL:10.6f},{A_CELL:10.6f},{90:10.6f},{90:10.6f},{90:10.6f}",
             f"ncell {n:8d},{n:8d},{n:8d},{1:10d}",
             "format numbers,XYZBPMMOS",
             "atoms      x,              y,              z,             Biso,    Property,"
             "  MoleNo,  MoleAt,   Occ,     St,  Sh,  Sk,  Sl"]
    for nm, p in zip(names, pos):
        lines.append(discus_atom_line(nm, p, biso=0.0))
    return "\n".join(lines) + "\n"


def discus_atom_line(name: str, xyz, biso: float, prop: int = 1, occ: float = 1.0) -> str:
    """One atom line exactly as DISCUS 6.23.01 writes it (discus_save.f90, FORMAT 4):
    (a4,3(1x,f14.6,','),4x,f10.6,',',i8,',',I8,',',I8,', ',F10.6,', ',A1,3(', ',I3))
    gate2-prereg-v1 wrote two blanks after the name, shifting every field one column; DISCUS reads
    'format numbers' lines with a FIXED-column format and rejected it (error -49, engine probe).
    """
    s = f"{name.upper():<4.4s}" + "".join(f" {v:14.6f}," for v in xyz)
    s += f"    {biso:10.6f},{prop:8d},{0:8d},{0:8d}, {occ:10.6f}, _" + "".join(f", {0:3d}" for _ in range(3))
    return s


def discus_numbers_read(line: str):
    """Emulate DISCUS 6.23.01 read_atom_line (atom_line.f90) for 'format numbers,XYZBPMMOS':
    read(line(5:), '(3(1x,f14.6,1x),4x,f10.6,1x,i8,1x,I8,1x,I8,2x,F10.6,2x,A1,3(2x,I3))').
    Returns (name, x, y, z, B) or raises ValueError where DISCUS would report error -49.
    Fortran blank fields read as zero; any other non-numeric character in a field is an error.
    """
    rest = line[4:]
    spec = [("x", 1), ("f", 14), ("x", 1)] * 3 + [("x", 4), ("f", 10), ("x", 1), ("i", 8), ("x", 1),
                                                   ("i", 8), ("x", 1), ("i", 8), ("x", 2), ("f", 10),
                                                   ("x", 2), ("a", 1)] + [("x", 2), ("i", 3)] * 3
    vals, pos = [], 0
    for kind, w in spec:
        field = rest[pos:pos + w].ljust(w); pos += w
        if kind == "x":
            continue
        txt = field.strip()
        if kind == "a":
            vals.append(field); continue
        if txt == "":
            vals.append(0); continue
        if kind == "i":
            if not txt.lstrip("+-").isdigit():
                raise ValueError(f"integer field {field!r}")
            vals.append(int(txt))
        else:
            try:
                vals.append(float(txt))
            except ValueError:
                raise ValueError(f"real field {field!r}") from None
    return line[:4].strip(), vals[0], vals[1], vals[2], vals[3]


def parse_stru(path: Path):
    """Read back a structure file written by stru_text (used to prove the file is what we analyse)."""
    names, pos = [], []
    for line in Path(path).read_text().splitlines():
        t = line.replace(",", " ").split()
        if t and t[0] in SPECIES:
            names.append(t[0]); pos.append([float(v) for v in t[1:4]])
    pos = np.asarray(pos)
    n = int(round(pos.max())) + 1
    sig = np.zeros((n, n, n), dtype=np.int8)
    for nm, p in zip(names, pos.astype(int)):
        sig[tuple(p)] = 1 if nm == SPECIES[1] else -1
    return sig


# ----------------------------------------------------------------------------- correlations
def pair_sums_direct(sig: np.ndarray) -> np.ndarray:
    """C(tau) = sum_j sigma_j sigma_{j+tau} (periodic), by explicit shifted products (no FFT)."""
    s = sig.astype(np.int64)
    n = s.shape[0]
    C = np.empty((n, n, n), dtype=np.int64)
    for tx in range(n):
        sx = np.roll(s, -tx, axis=0)
        for ty in range(n):
            sxy = np.roll(sx, -ty, axis=1)
            for tz in range(n):
                C[tx, ty, tz] = int(np.sum(s * np.roll(sxy, -tz, axis=2)))
    return C


def ensemble_sd(model: str, n_real: int = 2000, seed: int = 7) -> dict:
    """Monte-Carlo SD of C(tau)/N over independent realizations (FFT is fine here: sanity only)."""
    rng = np.random.default_rng(seed)
    vals = {t: [] for t in TAU_POINTS}
    for _ in range(n_real):
        sig = build_sigma(model, seed=int(rng.integers(2 ** 31)))
        C = np.rint(np.real(np.fft.ifftn(np.abs(np.fft.fftn(sig)) ** 2)) * 1).astype(np.int64)
        for t in TAU_POINTS:
            vals[t].append(C[tuple(np.mod(t, N_CELL))] / sig.size)
    return {t: (float(np.mean(v)), float(np.std(v, ddof=1))) for t, v in vals.items()}


# ----------------------------------------------------------------------------- intensities
def lattice_sums(sig: np.ndarray):
    """A_sp(k) = sum_{j in sp} exp(+2 pi i k.n_j / 16) for k = 0..15 (one Brillouin zone)."""
    n = sig.shape[0]
    occ_au = (sig > 0).astype(float)
    occ_cu = (sig < 0).astype(float)
    # numpy ifftn uses exp(+2 pi i k n / N) / N^3
    return {"AU": np.fft.ifftn(occ_au) * n ** 3, "CU": np.fft.ifftn(occ_cu) * n ** 3}


def forward(sig: np.ndarray, mode: str, wk: dict | None = None):
    """Return (I_total, I_diffuse, g2) on the full (128,128,128) grid, axes (h, k, l)."""
    wk = coefficients_for(mode) if wk is None else wk
    s = s_of_grid(mode)
    f = {sp: wk_f(wk[sp], s) for sp in SPECIES}
    A = lattice_sums(sig)
    k = np.mod(np.rint(H_AXIS * N_CELL).astype(int), N_CELL)
    ix = np.ix_(k, k, k)
    F = f["AU"] * A["AU"][ix] + f["CU"] * A["CU"][ix]
    n_au = int((sig > 0).sum()); n_cu = int((sig < 0).sum())
    fbar = (n_au * f["AU"] + n_cu * f["CU"]) / sig.size
    is_int = np.isclose(H_AXIS, np.rint(H_AXIS))
    L = np.zeros(s.shape); L[np.ix_(is_int, is_int, is_int)] = sig.size
    I_tot = np.abs(F) ** 2
    I_diff = np.abs(F - fbar * L) ** 2
    g2 = ((f["AU"] - f["CU"]) / 2.0) ** 2
    return I_tot, I_diff, g2


def direct_F(sig: np.ndarray, hkl: np.ndarray, mode: str = "physics", wk: dict | None = None):
    """Brute-force structure factor at arbitrary H points (m,3); used to verify forward()."""
    wk = coefficients_for(mode) if wk is None else wk
    names, pos = atom_list(sig)
    names = np.asarray(names)
    s = np.linalg.norm(hkl, axis=1) / (2 * A_CELL)
    if mode == "engine":
        s = np.floor(s * (1.0 / CFINC) + 0.5) * CFINC
    F = np.zeros(len(hkl), complex)
    for sp in SPECIES:
        ph = np.exp(2j * np.pi * hkl @ pos[names == sp].T).sum(1)
        F += wk_f(wk[sp], s) * ph
    return F


# ----------------------------------------------------------------------------- real space
def kernel_on_lattice(g2: np.ndarray) -> np.ndarray:
    """K(n) = (dH^3 / V) sum_H g2(H) exp(-2 pi i H.n) at integer lattice vectors n (16^3, periodic).

    Implemented by zone folding + explicit 16-point DFT matrices (not fftn), so it is an
    independent code path from the 128^3 transform used by the pipeline.
    """
    n = N_CELL
    k = np.mod(np.rint(H_AXIS * n).astype(int), n)
    # H = -4 + i/16 and 64 = 0 (mod 16): inside each 16-point zone block, grid index r IS k = r.
    assert np.array_equal(k[:n], np.arange(n))
    G = g2.reshape(N_H // n, n, N_H // n, n, N_H // n, n).sum(axis=(0, 2, 4))   # fold 8^3 zones
    m = np.arange(n)
    E = np.exp(-2j * np.pi * np.outer(m, m) / n)      # E[k, n]
    K = np.einsum("ijk,ia,jb,kc->abc", G, E, E, E, optimize=True)
    return (K * DH ** 3 / A_CELL ** 3).real


def predicted_rho(C: np.ndarray, K: np.ndarray, taus=TAU_POINTS) -> dict:
    """rho(tau) = sum_tau' C(tau') K(tau - tau'), periodic on the 16^3 lattice (direct sum)."""
    n = N_CELL
    idx = np.indices((n, n, n)).reshape(3, -1).T
    Cf = C.reshape(-1).astype(float)
    out = {}
    for t in taus:
        d = np.mod(np.asarray(t) - idx, n)
        out[t] = float(np.sum(Cf * K[d[:, 0], d[:, 1], d[:, 2]]))
    return out


# ----------------------------------------------------------------------------- utilities
def sha256_file(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def sha256_array(a: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(a, dtype=np.float64).tobytes()).hexdigest()


FINGERPRINT_VOXELS = [(0, 0, 0), (64, 64, 64), (72, 64, 64), (64, 72, 64), (64, 64, 72), (5, 77, 101),
                      (127, 3, 60), (40, 90, 13), (66, 70, 74), (100, 100, 100)]


def fingerprint(a: np.ndarray) -> dict:
    """Platform-robust summary of an oracle array (libm/FFT may differ by an ulp between machines)."""
    a = np.asarray(a, dtype=np.float64)
    return {"sum": float(a.sum()), "sum_sq": float(np.sum(a * a)),
            "voxels": [float(a[v]) for v in FINGERPRINT_VOXELS]}


def fingerprint_matches(a: np.ndarray, ref: dict, rtol: float = 1e-9) -> bool:
    f = fingerprint(a)
    scale = max(abs(x) for x in ref["voxels"]) or 1.0
    return (np.isclose(f["sum"], ref["sum"], rtol=rtol, atol=0) and
            np.isclose(f["sum_sq"], ref["sum_sq"], rtol=rtol, atol=0) and
            np.allclose(f["voxels"], ref["voxels"], rtol=rtol, atol=rtol * scale))


def rel_l2(a, b) -> float:
    return float(np.linalg.norm(np.asarray(a) - np.asarray(b)) / max(np.linalg.norm(b), 1e-300))
