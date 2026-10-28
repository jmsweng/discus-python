"""Correlated-disorder 3D-ΔPDF model for DISCUS 6.23.01 validation.

This module defines a minimal, fully controlled correlated-disorder model with
a prescribed nearest-neighbor correlation, and provides independent references
(pair-correlation function and PDF peak) for validating the DISCUS 3D-ΔPDF
output.

Model:
    - Cubic supercell: n=4, a=5.0 Å (supercell = 20 Å)
    - 4-site basis per unit cell: (0,0,0), (0.5,0,0), (0,0.5,0), (0,0,0.5)
    - Binary species: Ga (A) and As (B)
    - Species assignment: checkerboard pattern at the basis-site level

The frozen structure is a 4×4×4 supercell with 256 sites (128 Ga, 128 As).
The species are assigned in a checkerboard pattern at the basis-site level:
    - Ga where (i + j + k + site_idx) % 2 == 0
    - As where (i + j + k + site_idx) % 2 == 1

This gives 50% Ga and 50% As overall, with a specific correlation that can
be computed exactly from the frozen atom list.

The independent references are:
    1. Directional pair correlations: c(tau) computed at the 2.5 Å NN distance
       from the frozen atom list.
    2. Pair-shell statistics: g(r) computed from the frozen atom list with
       proper per-shell normalization.
    3. PDF peak predictions: derived from actual directional pair counts,
       scattering contrast, and boundary conditions.
    4. Structure factor: F(hkl) computed directly from the frozen atom list
       with full complex phase factors.
    5. Diffuse intensity: |F(hkl) - F_avg(hkl)|² where F_avg(hkl) is the
       full complex average structure factor.

IMPORTANT: The scattering factors used here are SYNTHETIC (atomic numbers as
proxy). For real X-ray scattering validation, use documented Q-dependent
form factors (e.g., IT92 or Waasmaier–Kirfel) and Debye–Waller conventions.

This model is designed to be:
    - Fully deterministic (no random sampling)
    - Independent of the DISCUS engine (analytical references)
    - Minimal (4 sites per unit cell, 4×4×4 supercell)
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np


def build_frozen_structure(n: int = 4, a: float = 5.0) -> dict:
    """Build a frozen n×n×n crystal with realized binary species assignments.

    The structure has a 4-site basis per unit cell:
        (0,0,0), (0.5,0,0), (0,0.5,0), (0,0,0.5)

    Species are assigned in a checkerboard pattern at the basis-site level:
        - Ga where (i + j + k + site_idx) % 2 == 0
        - As where (i + j + k + site_idx) % 2 == 1

    This gives 50% Ga and 50% As overall.

    The positions are in CARTESIAN coordinates (Å) for the full supercell.
    The supercell has dimensions n*a in each direction.

    Parameters
    ----------
    n : int
        Supercell size (default 4).
    a : float
        Unit cell length in Å (default 5.0).

    Returns
    -------
    dict with keys:
        n, a, basis, positions (Cartesian, Å), species
    """
    # 4-site basis in unit cell fractional coordinates
    basis = [
        (0.0, 0.0, 0.0),
        (0.5, 0.0, 0.0),
        (0.0, 0.5, 0.0),
        (0.0, 0.0, 0.5),
    ]

    # Generate all site positions in the supercell
    # Positions are in Cartesian coordinates (Å)
    positions = []
    species = []
    for i in range(n):
        for j in range(n):
            for k in range(n):
                for site_idx, (bx, by, bz) in enumerate(basis):
                    # Convert to Cartesian: unit cell index + basis offset
                    x = (i + bx) * a
                    y = (j + by) * a
                    z = (k + bz) * a
                    positions.append((x, y, z))
                    # Checkerboard pattern at basis-site level
                    sp = "Ga" if (i + j + k + site_idx) % 2 == 0 else "As"
                    species.append(sp)

    return {
        "n": n,
        "a": a,
        "basis": basis,
        "positions": positions,  # Cartesian, Å
        "species": species,
    }


def write_stru(path: Path, struct: dict) -> None:
    """Write a DISCUS .stru file for the frozen correlated model.

    The .stru file contains the full frozen crystal with realized species assignments.
    Each line is: atom_type x y z occupancy Biso

    IMPORTANT: DISCUS .stru files use FRACTIONAL coordinates relative to the
    declared cell. Since we're declaring a 5 Å cell but have a 20 Å supercell,
    we need to convert Cartesian positions to fractional coordinates of the
    SUPERCELL (20 Å), not the unit cell.

    The convention is:
    - Declare the supercell dimensions (20 Å) in the cell line
    - Use fractional coordinates where 0 to 1 spans the full supercell
    """
    n = struct["n"]
    a = struct["a"]
    positions = struct["positions"]  # Cartesian, Å
    species = struct["species"]

    # Supercell dimensions
    supercell_a = n * a

    lines = [
        "# Frozen correlated-disorder model for 3D-dPDF validation",
        f"# Supercell: n={n}, a={a} Å (supercell = {supercell_a} Å)",
        f"# Total sites: {len(species)}",
        f"spcgr 1",
        f"cell {supercell_a} {supercell_a} {supercell_a} 90 90 90",
    ]

    for (x, y, z), sp in zip(positions, species):
        # Convert Cartesian to fractional (relative to supercell)
        x_frac = x / supercell_a
        y_frac = y / supercell_a
        z_frac = z / supercell_a
        # Occupancy: 1.0 for the assigned species
        occ = 1.0
        lines.append(f"{sp} {x_frac:.6f} {y_frac:.6f} {z_frac:.6f} {occ:.6f} 1.0")

    path.write_text("\n".join(lines) + "\n")


def compute_directional_pair_correlations(struct: dict) -> dict:
    """Compute the actual signed connected pair correlations at the NN distance.

    For a binary system with species A and B, the occupation variable is:
        η_i = +1 if site i is A, -1 if site i is B

    The connected covariance at displacement tau is:
        c(tau) = <η_i * η_{i+tau}> - <η_i> * <η_{i+tau}>

    For p=0.5, <η_i> = 0, so c(tau) = <η_i * η_{i+tau}>.

    The normalized correlation coefficient is:
        ρ(tau) = c(tau) / sqrt(<η_i²> * <η_{i+tau}²>)

    For p=0.5, <η_i²> = 1, so ρ(tau) = c(tau).

    This function computes c(tau) at the 2.5 Å NN distance (a/2) along x, y, z
    and their negatives, directly from the frozen atom list.

    Returns
    -------
    dict with keys:
        c_tau : dict mapping displacement vector (Å) to covariance
        rho_tau : dict mapping displacement vector (Å) to normalized correlation
        n_pairs : dict mapping displacement vector (Å) to pair count
    """
    n = struct["n"]
    a = struct["a"]
    positions = np.array(struct["positions"])  # Cartesian, Å
    species = struct["species"]

    # Convert species to occupation variables: +1 for Ga, -1 for As
    eta = np.array([1.0 if sp == "Ga" else -1.0 for sp in species])

    # Compute <η_i> (should be 0 for p=0.5)
    eta_mean = np.mean(eta)

    # NN displacements: ±a/2 along x, y, z
    tau_vectors = [
        (a / 2.0, 0.0, 0.0),
        (-a / 2.0, 0.0, 0.0),
        (0.0, a / 2.0, 0.0),
        (0.0, -a / 2.0, 0.0),
        (0.0, 0.0, a / 2.0),
        (0.0, 0.0, -a / 2.0),
    ]

    c_tau = {}
    rho_tau = {}
    n_pairs = {}

    for tau in tau_vectors:
        count = 0
        product_sum = 0.0

        # Compute all pairwise displacements
        diff = positions[:, np.newaxis, :] - positions[np.newaxis, :, :]  # (N, N, 3)

        for i in range(len(eta)):
            # Find all sites j such that position_j - position_i ≈ tau
            dist_to_tau = np.sqrt(np.sum((diff[i, :, :] - np.array(tau)) ** 2, axis=1))
            mask = np.isclose(dist_to_tau, 0.0, atol=1e-6)
            mask[i] = False  # Exclude self

            if np.any(mask):
                count += np.sum(mask)
                product_sum += np.sum(eta[i] * eta[mask])

        if count > 0:
            c_val = product_sum / count - eta_mean**2
            c_tau[tau] = c_val
            rho_tau[tau] = c_val  # For p=0.5, ρ = c
            n_pairs[tau] = count

    return {
        "c_tau": c_tau,
        "rho_tau": rho_tau,
        "n_pairs": n_pairs,
        "eta_mean": float(eta_mean),
    }


def compute_pair_shell_statistics(struct: dict) -> dict:
    """Compute properly normalized pair-shell statistics from the frozen atom list.

    For a binary system, the pair-correlation function is:
        g_AB(r) = <n_A(r) * n_B(r + tau)> / (<n_A> * <n_B> / V)

    where n_A(r) is the number density of A atoms at position r, and V is
    the volume.

    For a frozen structure, we compute g(r) by counting pairs at each distance
    and normalizing by the EXPECTED number of pairs in that specific distance
    shell for a random distribution.

    The expected number of A-B pairs in a shell of volume ΔV at distance r is:
        E[N_AB] = n_A * n_B * ΔV / V

    where ΔV is the volume of the shell.

    Returns
    -------
    dict with keys:
        r : list of distances (Å)
        g_ab : list of g_AB(r) values (cross-correlation)
        g_aa : list of g_AA(r) values (auto-correlation)
        g_bb : list of g_BB(r) values (auto-correlation)
        shell_volumes : list of shell volumes (Å³)
    """
    n = struct["n"]
    a = struct["a"]
    positions = np.array(struct["positions"])  # Cartesian, Å
    species = struct["species"]

    n_total = len(positions)
    n_ga = sum(1 for sp in species if sp == "Ga")
    n_as = sum(1 for sp in species if sp == "As")

    # Volume of the supercell
    V = (n * a) ** 3

    # Create species array: 1 for Ga, 0 for As
    sp_array = np.array([1.0 if sp == "Ga" else 0.0 for sp in species])

    # Distances to consider: 0, a/2, a, 1.5a, 2a, 2.5a, 3a
    r_values = [0.0, a / 2.0, a, 1.5 * a, 2.0 * a, 2.5 * a, 3.0 * a]
    g_ab = []
    g_aa = []
    g_bb = []
    shell_volumes = []

    # Shell width (Å)
    dr = 0.1

    for r in r_values:
        if r == 0:
            g_ab.append(1.0)
            g_aa.append(1.0)
            g_bb.append(1.0)
            shell_volumes.append(0.0)
        else:
            # Compute all pairwise distances (vectorized)
            diff = positions[:, np.newaxis, :] - positions[np.newaxis, :, :]
            dist = np.sqrt(np.sum(diff**2, axis=2))

            # Shell: r - dr/2 < dist < r + dr/2
            mask = (dist > r - dr / 2.0) & (dist < r + dr / 2.0)
            np.fill_diagonal(mask, False)

            # Count A-B, A-A, B-B pairs
            count_ab = np.sum((sp_array[:, np.newaxis] == 1) & (sp_array[np.newaxis, :] == 0) & mask)
            count_ab += np.sum((sp_array[:, np.newaxis] == 0) & (sp_array[np.newaxis, :] == 1) & mask)
            count_aa = np.sum((sp_array[:, np.newaxis] == 1) & (sp_array[np.newaxis, :] == 1) & mask)
            count_bb = np.sum((sp_array[:, np.newaxis] == 0) & (sp_array[np.newaxis, :] == 0) & mask)

            # Shell volume: 4πr² * dr (for a thin shell)
            shell_vol = 4.0 * np.pi * r**2 * dr
            shell_volumes.append(shell_vol)

            # Expected pairs in this shell for random distribution
            # E[N_AB] = n_A * n_B * shell_vol / V
            expected_ab = n_ga * n_as * shell_vol / V
            # E[N_AA] = n_A * (n_A - 1) * shell_vol / V
            expected_aa = n_ga * (n_ga - 1) * shell_vol / V
            # E[N_BB] = n_B * (n_B - 1) * shell_vol / V
            expected_bb = n_as * (n_as - 1) * shell_vol / V

            g_ab.append(count_ab / expected_ab if expected_ab > 0 else 0.0)
            g_aa.append(count_aa / expected_aa if expected_aa > 0 else 0.0)
            g_bb.append(count_bb / expected_bb if expected_bb > 0 else 0.0)

    return {
        "r": r_values,
        "g_ab": g_ab,
        "g_aa": g_aa,
        "g_bb": g_bb,
        "shell_volumes": shell_volumes,
        "n_ga": n_ga,
        "n_as": n_as,
        "n_total": n_total,
        "V": V,
    }


def predict_pdf_peaks(struct: dict, pair_corr: dict) -> dict:
    """Predict signed 3D-ΔPDF peak positions and amplitudes from the frozen structure.

    The 3D-ΔPDF peak at displacement tau is proportional to:
        Δρ(tau) = (f_A - f_B)² * c(tau) * N_pairs(tau) / V

    where:
    - f_A, f_B are the scattering factors for species A and B
    - c(tau) is the connected occupation covariance
    - N_pairs(tau) is the number of pairs at displacement tau
    - V is the volume of the supercell

    For a binary system with p=0.5, the contrast is (f_A - f_B)².
    The peak amplitude includes the pair multiplicity and is normalized by volume.

    Returns
    -------
    dict with keys:
        tau : list of peak positions (Cartesian, Å)
        amplitude : list of expected peak amplitudes (normalized by volume)
        c_tau : list of covariance values at each peak
        n_pairs : list of pair multiplicities at each peak
        contrast : (f_A - f_B)²
    """
    a = struct["a"]
    n = struct["n"]
    V = (n * a) ** 3

    # Scattering contrast: (f_Ga - f_As)²
    # SYNTHETIC: using atomic numbers as proxy (Ga: Z=31, As: Z=33)
    f_ga = 31.0
    f_as = 33.0
    contrast = (f_ga - f_as) ** 2  # = 4.0

    # NN displacements: ±a/2 along x, y, z
    tau_vectors = [
        (a / 2.0, 0.0, 0.0),
        (-a / 2.0, 0.0, 0.0),
        (0.0, a / 2.0, 0.0),
        (0.0, -a / 2.0, 0.0),
        (0.0, 0.0, a / 2.0),
        (0.0, 0.0, -a / 2.0),
    ]

    c_tau = []
    n_pairs = []
    amplitude = []

    for tau in tau_vectors:
        # Get covariance and pair count from the computed pair correlations
        c_val = pair_corr["c_tau"].get(tau, 0.0)
        n_pairs_val = pair_corr["n_pairs"].get(tau, 0)

        c_tau.append(c_val)
        n_pairs.append(n_pairs_val)

        # Amplitude: contrast * c(tau) * N_pairs(tau) / V
        amp = contrast * c_val * n_pairs_val / V
        amplitude.append(amp)

    return {
        "tau": tau_vectors,
        "amplitude": amplitude,
        "c_tau": c_tau,
        "n_pairs": n_pairs,
        "contrast": contrast,
        "V": V,
    }


def compute_structure_factor(struct: dict, hkl_list: list) -> np.ndarray:
    """Compute the structure factor F(hkl) directly from the frozen atom list.

    F(hkl) = sum_i f_i * exp(2πi * (h*x_i + k*y_i + l*z_i))

    where f_i is the scattering factor for atom i, and (x_i, y_i, z_i)
    are the FRACTIONAL coordinates (0 to 1 spans the supercell).

    IMPORTANT: The positions in struct are in Cartesian (Å). We need to
    convert to fractional coordinates relative to the supercell.

    Parameters
    ----------
    struct : dict
        Structure dictionary with positions (Cartesian, Å) and species.
    hkl_list : list
        List of (h, k, l) tuples to compute F for.

    Returns
    -------
    np.ndarray of complex structure factors, shape (len(hkl_list),)
    """
    n = struct["n"]
    a = struct["a"]
    supercell_a = n * a
    positions_cart = np.array(struct["positions"])  # Cartesian, Å
    species = struct["species"]

    # Convert to fractional coordinates (relative to supercell)
    positions_frac = positions_cart / supercell_a

    # Scattering factors (SYNTHETIC: atomic numbers as proxy)
    # Ga: Z=31, As: Z=33
    f_ga = 31.0
    f_as = 33.0

    f = np.array([f_ga if sp == "Ga" else f_as for sp in species])

    # Compute F(hkl) for each hkl
    F = np.empty(len(hkl_list), complex)
    for idx, (h, k, l) in enumerate(hkl_list):
        phase = 2j * np.pi * (
            h * positions_frac[:, 0]
            + k * positions_frac[:, 1]
            + l * positions_frac[:, 2]
        )
        F[idx] = np.sum(f * np.exp(phase))

    return F


def compute_average_structure_factor(struct: dict, hkl_list: list) -> np.ndarray:
    """Compute the average structure factor F_avg(hkl) for the disordered system.

    For a binary system with occupation probabilities p_A and p_B at each site,
    the average scattering factor at site i is:
        f_avg_i = p_A * f_A + p_B * f_B

    The average structure factor is:
        F_avg(hkl) = sum_i f_avg_i * exp(2πi * (h*x_i + k*y_i + l*z_i))

    For p_A = p_B = 0.5, f_avg_i = 0.5 * (f_A + f_B) for all sites.

    Returns
    -------
    np.ndarray of complex average structure factors, shape (len(hkl_list),)
    """
    n = struct["n"]
    a = struct["a"]
    supercell_a = n * a
    positions_cart = np.array(struct["positions"])  # Cartesian, Å

    # Convert to fractional coordinates (relative to supercell)
    positions_frac = positions_cart / supercell_a

    # Scattering factors (SYNTHETIC: atomic numbers as proxy)
    f_ga = 31.0
    f_as = 33.0

    # Average scattering factor: p_A * f_A + p_B * f_B = 0.5 * (f_A + f_B)
    f_avg = 0.5 * (f_ga + f_as)

    # Compute F_avg(hkl) for each hkl
    F_avg = np.empty(len(hkl_list), complex)
    for idx, (h, k, l) in enumerate(hkl_list):
        phase = 2j * np.pi * (
            h * positions_frac[:, 0]
            + k * positions_frac[:, 1]
            + l * positions_frac[:, 2]
        )
        F_avg[idx] = np.sum(f_avg * np.exp(phase))

    return F_avg


def compute_diffuse_intensity(struct: dict, hkl_list: list) -> np.ndarray:
    """Compute the diffuse intensity I_diff(hkl) = |F(hkl) - F_avg(hkl)|².

    F(hkl) is the structure factor of the frozen (realized) structure.
    F_avg(hkl) is the average structure factor (full complex sum with phase factors).

    Returns
    -------
    np.ndarray of diffuse intensities, shape (len(hkl_list),)
    """
    F = compute_structure_factor(struct, hkl_list)
    F_avg = compute_average_structure_factor(struct, hkl_list)

    # Diffuse intensity: |F - F_avg|²
    I_diff = np.abs(F - F_avg) ** 2

    return I_diff


def main() -> None:
    """Build the frozen correlated model and write all outputs."""
    out_dir = Path(__file__).parent / "correlated_run"
    out_dir.mkdir(parents=True, exist_ok=True)

    struct = build_frozen_structure(n=4, a=5.0)

    # Write .stru file
    stru_path = out_dir / "corr_model.stru"
    write_stru(stru_path, struct)
    print(f"Wrote {stru_path}")

    # Compute independent references
    pair_corr = compute_directional_pair_correlations(struct)
    pcf = compute_pair_shell_statistics(struct)
    peak = predict_pdf_peaks(struct, pair_corr)

    # Compute structure factor and diffuse intensity for validation
    hkl_test = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (0.5, 0, 0)]
    F_test = compute_structure_factor(struct, hkl_test)
    F_avg_test = compute_average_structure_factor(struct, hkl_test)
    I_diff_test = compute_diffuse_intensity(struct, hkl_test)

    # Write references
    refs_path = out_dir / "references.json"
    refs_data = {
        "structure": {
            "n": int(struct["n"]),
            "a": float(struct["a"]),
            "basis": struct["basis"],
            "n_sites": int(len(struct["positions"])),
            "n_ga": int(sum(1 for sp in struct["species"] if sp == "Ga")),
            "n_as": int(sum(1 for sp in struct["species"] if sp == "As")),
            "note": "Positions are in Cartesian (Å). Supercell = 20 Å.",
        },
        "directional_pair_correlations": {
            "eta_mean": float(pair_corr["eta_mean"]),
            "c_tau": {
                f"{k[0]:.1f},{k[1]:.1f},{k[2]:.1f}": float(v)
                for k, v in pair_corr["c_tau"].items()
            },
            "rho_tau": {
                f"{k[0]:.1f},{k[1]:.1f},{k[2]:.1f}": float(v)
                for k, v in pair_corr["rho_tau"].items()
            },
            "n_pairs": {
                f"{k[0]:.1f},{k[1]:.1f},{k[2]:.1f}": int(v)
                for k, v in pair_corr["n_pairs"].items()
            },
        },
        "pair_shell_statistics": {
            "r": [float(r) for r in pcf["r"]],
            "g_ab": [float(g) for g in pcf["g_ab"]],
            "g_aa": [float(g) for g in pcf["g_aa"]],
            "g_bb": [float(g) for g in pcf["g_bb"]],
            "shell_volumes": [float(v) for v in pcf["shell_volumes"]],
            "n_ga": int(pcf["n_ga"]),
            "n_as": int(pcf["n_as"]),
            "n_total": int(pcf["n_total"]),
            "V": float(pcf["V"]),
        },
        "pdf_peaks": {
            "tau": peak["tau"],
            "amplitude": [float(a) for a in peak["amplitude"]],
            "c_tau": [float(c) for c in peak["c_tau"]],
            "n_pairs": [int(n) for n in peak["n_pairs"]],
            "contrast": float(peak["contrast"]),
            "V": float(peak["V"]),
            "note": "SYNTHETIC: contrast uses atomic numbers as proxy (Ga:31, As:33).",
        },
        "validation": {
            "hkl_test": hkl_test,
            "F_test": [f"{f.real:.6f}+{f.imag:.6f}j" for f in F_test],
            "F_avg_test": [f"{f.real:.6f}+{f.imag:.6f}j" for f in F_avg_test],
            "I_diff_test": [f"{v:.6f}" for v in I_diff_test],
            "note": "SYNTHETIC: structure factors use atomic numbers as proxy.",
        },
    }
    refs_path.write_text(json.dumps(refs_data, indent=2))
    print(f"Wrote {refs_path}")

    # Print summary
    n_ga = sum(1 for sp in struct["species"] if sp == "Ga")
    n_as = sum(1 for sp in struct["species"] if sp == "As")
    print(f"\nStructure: n={struct['n']}, a={struct['a']} Å (supercell = {struct['n'] * struct['a']} Å)")
    print(f"Basis: {struct['basis']}")
    print(f"Total sites: {len(struct['positions'])}")
    print(f"Ga sites: {n_ga}, As sites: {n_as}")
    print(f"eta_mean: {pair_corr['eta_mean']:.6f}")
    print(f"\nDirectional pair correlations (at ±2.5 Å):")
    for tau, c in pair_corr["c_tau"].items():
        n_pairs = pair_corr["n_pairs"][tau]
        print(f"  tau={tau} Å: c={c:.6f}, N_pairs={n_pairs}")
    print(f"\nPair-shell statistics:")
    for r, g_ab, g_aa, g_bb, shell_vol in zip(
        pcf["r"], pcf["g_ab"], pcf["g_aa"], pcf["g_bb"], pcf["shell_volumes"]
    ):
        print(f"  r={r:.2f} Å: g_AB={g_ab:.4f}, g_AA={g_aa:.4f}, g_BB={g_bb:.4f} (shell_vol={shell_vol:.2f} Å³)")
    print(f"\nExpected PDF peaks (SYNTHETIC contrast):")
    for tau, amp, c, np in zip(peak["tau"], peak["amplitude"], peak["c_tau"], peak["n_pairs"]):
        print(f"  tau={tau} Å: c={c:.4f}, N_pairs={np}, amplitude={amp:.6f}")
    print(f"\nValidation (SYNTHETIC structure factors):")
    for hkl, F, F_avg, I in zip(hkl_test, F_test, F_avg_test, I_diff_test):
        print(f"  hkl={hkl}: F={F:.4f}, F_avg={F_avg:.4f}, I_diff={I:.4f}")


if __name__ == "__main__":
    main()
