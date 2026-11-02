"""Brute-force anchor check for the Gate 0 oracle.

For realizations A_00, A_31, B_63, compare the oracle I with an independent
brute-force sum (no FFT, no site-level assumption) at 20 seeded random voxels.
Assert max relative difference <= 1e-10. Print the three anchor values for A_00.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import oracle


def brute_force_I(path: Path, mode: str, b_cu: float = 0.0, b_au: float = 0.0) -> np.ndarray:
    """Independent brute-force I(h) = |F(h)|², no FFT, no site-level assumption.

    F(h) = Σ_n f_el(s) exp(-B_n s²) exp(2πi h·x_n)
    B_n: from the file in arm B, 0 in arm A (temp ignore).
    """
    s = oracle.s_of_grid(mode)
    wk = oracle.coefficients_for(mode)
    f_au = oracle.wk_f(wk["AU"], s)
    f_cu = oracle.wk_f(wk["CU"], s)

    use_file_B = (b_cu != 0.0 or b_au != 0.0)

    atoms = []
    for line in path.read_text().splitlines():
        t = line.replace(",", " ").split()
        if not t or t[0] not in oracle.SPECIES:
            continue
        sp = t[0]
        x, y, z = float(t[1]), float(t[2]), float(t[3])
        B_file = float(t[4])
        B = B_file if use_file_B else 0.0
        atoms.append((sp, x, y, z, B))

    F = np.zeros(s.shape, complex)
    for sp, x, y, z, B in atoms:
        f = f_au if sp == "AU" else f_cu
        amp = f * np.exp(-B * s ** 2)
        phase = np.exp(2j * np.pi * (
            oracle.H_AXIS[:, None, None] * x +
            oracle.H_AXIS[None, :, None] * y +
            oracle.H_AXIS[None, None, :] * z
        ))
        F += amp * phase
    return np.abs(F) ** 2


def main() -> int:
    rng = np.random.default_rng(42)
    tests = [
        ("A_00", "engine", 0.0, 0.0),
        ("A_31", "engine", 0.0, 0.0),
        ("B_63", "engine", oracle.B_CU_B, oracle.B_AU_B),
    ]
    max_rel = 0.0
    for name, mode, b_cu, b_au in tests:
        path = HERE / "inputs" / f"gate0_{name}.stru"
        _, I_oracle = oracle.realization_exact(path, mode, b_cu=b_cu, b_au=b_au)
        I_brute = brute_force_I(path, mode, b_cu=b_cu, b_au=b_au)
        # 20 seeded random voxels
        idx = rng.integers(0, 25, size=(20, 3))
        for i, j, k in idx:
            rel = abs(I_oracle[i, j, k] - I_brute[i, j, k]) / max(I_brute[i, j, k], 1e-300)
            max_rel = max(max_rel, rel)
        print(f"{name}: max rel diff over 20 voxels = {max_rel:.2e}")

    assert max_rel <= 1e-10, f"Max relative difference {max_rel:.2e} > 1e-10"

    # Print the three anchor values for A_00
    path = HERE / "inputs" / "gate0_A_00.stru"
    _, I_oracle = oracle.realization_exact(path, "engine", b_cu=0.0, b_au=0.0)
    print(f"\nA_00 anchor values (arm A, engine mode):")
    print(f"  h=(1/8,0,0):    {I_oracle[1, 0, 0]:.2f}")
    print(f"  h=(3/8,2/8,1/8): {I_oracle[3, 2, 1]:.2f}")
    print(f"  h=(0,0,1/8):    {I_oracle[0, 0, 1]:.2f}")

    print(f"\nPASS: max relative difference = {max_rel:.2e} <= 1e-10")
    return 0


if __name__ == "__main__":
    sys.exit(main())
