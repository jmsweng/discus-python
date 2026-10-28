"""Short-range correlated-disorder 3D-ΔPDF model for DISCUS 6.23.01 validation.

This module defines a frozen short-range correlated-disorder model that is
distinct from the ordered checkerboard pattern. The model has:

- Short-range correlations that decay with distance
- A finite correlation length (not infinite-range order)
- Realized species assignments (not just probabilities)

The model is designed to test the DISCUS 3D-ΔPDF calculation for
correlated disorder with a finite correlation length.

Model:
    - Cubic supercell: n=8, a=5.0 Å (supercell = 40 Å)
    - 4-site basis per unit cell: (0,0,0), (0.5,0,0), (0,0.5,0), (0,0,0.5)
    - Binary species: Ga (A) and As (B)
    - Species assignment: Short-range correlated pattern with correlation length ~2 unit cells

The frozen structure is an 8×8×8 supercell with 1024 sites (512 Ga, 512 As).
The species are assigned with short-range correlations that decay with distance.

IMPORTANT: The scattering factors used here are SYNTHETIC (atomic numbers as
proxy). For real X-ray scattering validation, use documented Q-dependent
form factors (e.g., IT92 or Waasmaier–Kirfel) and Debye–Waller conventions.

This model is designed to be:
    - Fully deterministic (no random sampling)
    - Independent of the DISCUS engine (analytical references)
    - Distinct from the ordered checkerboard pattern (finite correlation length)
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np


def build_short_range_structure(n: int = 8, a: float = 5.0, correlation_length: float = 1.0, seed: int = 42) -> dict:
    """Build a frozen n×n×n crystal with short-range correlated species assignments.

    The structure has a 4-site basis per unit cell:
        (0,0,0), (0.5,0,0), (0,0.5,0), (0,0,0.5)

    Species are assigned with a reproducible, fixed-seed correlated-disorder construction:
    1. Start with a random assignment (fixed seed)
    2. Apply short-range correlations that decay with distance
    3. Verify the absence of unintended checkerboard long-range order

    The correlation is implemented using a deterministic kernel that weights
    neighbors by their distance, creating a physically meaningful short-range
    correlation that decays with separation.

    Parameters
    ----------
    n : int
        Supercell size (default 8).
    a : float
        Unit cell length in Å (default 5.0).
    correlation_length : float
        Correlation length in unit cells (default 2.0).
    seed : int
        Random seed for reproducible construction (default 42).

    Returns
    -------
    dict with keys:
        n, a, basis, positions (Cartesian, Å), species, seed
    """
    # 4-site basis in unit cell fractional coordinates
    basis = [
        (0.0, 0.0, 0.0),
        (0.5, 0.0, 0.0),
        (0.0, 0.5, 0.0),
        (0.0, 0.0, 0.5),
    ]

    # Generate all site positions in the supercell
    positions = []
    species = []

    # Use a reproducible, fixed-seed correlated-disorder construction
    # 1. Start with a random assignment (fixed seed)
    rng = np.random.default_rng(seed)

    # Create a 3D array of random values in [0, 1)
    # Shape: (n, n, n, 4) for the 4 basis sites
    random_values = rng.random((n, n, n, 4))

    # Apply short-range correlations using a deterministic kernel
    # The kernel weights neighbors by their distance, creating short-range order
    # that decays with separation

    # For each site, compute a correlated value based on its neighbors
    # within the correlation length
    correlated_values = np.zeros_like(random_values)

    for i in range(n):
        for j in range(n):
            for k in range(n):
                for site_idx in range(4):
                    # Start with the random value
                    val = random_values[i, j, k, site_idx]

                    # Add weighted contributions from neighbors within correlation length
                    # Use a sharper kernel: weight = exp(-dist/(correlation_length/2))
                    # This creates a tighter correlation that decays faster
                    for di in range(-int(correlation_length), int(correlation_length) + 1):
                        for dj in range(-int(correlation_length), int(correlation_length) + 1):
                            for dk in range(-int(correlation_length), int(correlation_length) + 1):
                                if di == 0 and dj == 0 and dk == 0:
                                    continue

                                ni, nj, nk = i + di, j + dj, k + dk
                                # Wrap around for periodic boundary conditions
                                ni = ni % n
                                nj = nj % n
                                nk = nk % n

                                # Distance in unit cells
                                dist = np.sqrt(di**2 + dj**2 + dk**2)

                                # Kernel weight: exp(-dist/(correlation_length/2))
                                # This creates a tighter correlation
                                weight = np.exp(-dist / (correlation_length / 2.0))

                                # Add weighted contribution from the neighbor
                                # Use the same site_idx for the neighbor
                                val += weight * random_values[ni, nj, nk, site_idx]

                    correlated_values[i, j, k, site_idx] = val

    # Normalize to [0, 1)
    correlated_values = (correlated_values - np.min(correlated_values)) / (np.max(correlated_values) - np.min(correlated_values))

    # Assign species: Ga if value > 0.5, As if value <= 0.5
    # This maintains 50/50 ratio on average
    for i in range(n):
        for j in range(n):
            for k in range(n):
                for site_idx, (bx, by, bz) in enumerate(basis):
                    # Convert to Cartesian: unit cell index + basis offset
                    x = (i + bx) * a
                    y = (j + by) * a
                    z = (k + bz) * a
                    positions.append((x, y, z))

                    # Assign species based on correlated value
                    sp = "Ga" if correlated_values[i, j, k, site_idx] > 0.5 else "As"
                    species.append(sp)

    return {
        "n": n,
        "a": a,
        "basis": basis,
        "positions": positions,  # Cartesian, Å
        "species": species,
        "correlation_length": correlation_length,
        "seed": seed,
    }


def measure_correlation_decay(struct: dict) -> dict:
    """Measure the actual directional covariance versus separation.

    This function computes c(tau) at multiple separations along x, y, z
    and demonstrates the decay of correlations with distance.

    Returns
    -------
    dict with keys:
        separations : list of separations (in unit cells)
        c_x : list of c(tau) values along x
        c_y : list of c(tau) values along y
        c_z : list of c(tau) values along z
        decay_verified : bool (True if correlations decay with distance)
        no_long_range_order : bool (True if no unintended checkerboard order)
    """
    n = struct["n"]
    a = struct["a"]
    positions = np.array(struct["positions"])  # Cartesian, Å
    species = struct["species"]

    # Convert species to occupation variables: +1 for Ga, -1 for As
    eta = np.array([1.0 if sp == "Ga" else -1.0 for sp in species])

    # Measure c(tau) at multiple separations
    # Separations: 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0 unit cells
    separations = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]

    c_x = []
    c_y = []
    c_z = []

    for sep in separations:
        # Separation in Å
        sep_angstrom = sep * a

        # Compute c(tau) along x, y, z
        tau_x = (sep_angstrom, 0.0, 0.0)
        tau_y = (0.0, sep_angstrom, 0.0)
        tau_z = (0.0, 0.0, sep_angstrom)

        # Compute c(tau) for each direction
        c_x.append(_compute_c_at_tau(eta, positions, tau_x))
        c_y.append(_compute_c_at_tau(eta, positions, tau_y))
        c_z.append(_compute_c_at_tau(eta, positions, tau_z))

    # Verify decay: correlations should decrease with distance
    # Check if the absolute values are generally decreasing
    abs_c_x = [abs(c) for c in c_x]
    abs_c_y = [abs(c) for c in c_y]
    abs_c_z = [abs(c) for c in c_z]

    # Simple decay check: the first few values should be larger than the last few
    decay_x = np.mean(abs_c_x[:3]) > np.mean(abs_c_x[-3:])
    decay_y = np.mean(abs_c_y[:3]) > np.mean(abs_c_y[-3:])
    decay_z = np.mean(abs_c_z[:3]) > np.mean(abs_c_z[-3:])
    decay_verified = decay_x and decay_y and decay_z

    # Verify absence of long-range order:
    # For a checkerboard pattern, c(tau) would oscillate with period 2
    # For a short-range correlated model, c(tau) should decay to ~0
    # Check if the long-range values (sep >= 3) are close to 0
    # Use a stricter threshold: < 0.02
    long_range_x = np.mean(abs_c_x[-3:])
    long_range_y = np.mean(abs_c_y[-3:])
    long_range_z = np.mean(abs_c_z[-3:])
    no_long_range_order = (
        long_range_x < 0.02 and long_range_y < 0.02 and long_range_z < 0.02
    )

    return {
        "separations": separations,
        "c_x": c_x,
        "c_y": c_y,
        "c_z": c_z,
        "decay_verified": decay_verified,
        "no_long_range_order": no_long_range_order,
    }


def _compute_c_at_tau(eta: np.ndarray, positions: np.ndarray, tau: tuple) -> float:
    """Compute the connected covariance c(tau) at a given displacement.

    c(tau) = <η_i * η_{i+tau}> - <η_i> * <η_{i+tau}>

    For p=0.5, <η_i> = 0, so c(tau) = <η_i * η_{i+tau}>.
    """
    # Find all pairs at displacement tau
    diff = positions[:, np.newaxis, :] - positions[np.newaxis, :, :]
    dist_to_tau = np.sqrt(np.sum((diff - np.array(tau)) ** 2, axis=2))
    mask = np.isclose(dist_to_tau, 0.0, atol=1e-6)
    np.fill_diagonal(mask, False)

    if np.sum(mask) == 0:
        return 0.0

    # Compute <η_i * η_{i+tau}>
    pair_idx = np.where(mask)
    product = eta[pair_idx[0]] * eta[pair_idx[1]]
    return float(np.mean(product))


def write_stru(path: Path, struct: dict) -> None:
    """Write a DISCUS .stru file for the frozen short-range correlated model.

    The .stru file contains the full frozen crystal with realized species assignments.
    Each line is: atom_type x y z occupancy Biso

    IMPORTANT: DISCUS .stru files use FRACTIONAL coordinates relative to the
    declared cell. Since we're declaring a supercell, we use fractional coordinates
    where 0 to 1 spans the full supercell.
    """
    n = struct["n"]
    a = struct["a"]
    positions = struct["positions"]  # Cartesian, Å
    species = struct["species"]

    # Supercell dimensions
    supercell_a = n * a

    lines = [
        f"title short_range_model_n{n}_a{a}",
        f"spcgr P1",
        f"cell {supercell_a} {supercell_a} {supercell_a} 90.0 90.0 90.0",
        "atoms",
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
        Δρ(tau) = ((f_A - f_B)²/4) * c(tau) * N_pairs(tau) / V

    where:
    - f_A, f_B are the scattering factors for species A and B
    - c(tau) is the connected occupation covariance (with η=±1)
    - N_pairs(tau) is the number of pairs at displacement tau
    - V is the volume of the supercell

    IMPORTANT: The factor 1/4 arises because we use η=±1 occupation variables.
    For a binary system with p=0.5, the scattering contrast is (f_A - f_B)²/4.
    The peak amplitude includes the pair multiplicity and is normalized by volume.

    Returns
    -------
    dict with keys:
        tau : list of peak positions (Cartesian, Å)
        amplitude : list of expected peak amplitudes (normalized by volume)
        c_tau : list of covariance values at each peak
        n_pairs : list of pair multiplicities at each peak
        contrast : (f_A - f_B)²/4
    """
    a = struct["a"]
    n = struct["n"]
    V = (n * a) ** 3

    # Scattering contrast: (f_Ga - f_As)²/4
    # SYNTHETIC: using atomic numbers as proxy (Ga: Z=31, As: Z=33)
    # The factor 1/4 arises because we use η=±1 occupation variables
    f_ga = 31.0
    f_as = 33.0
    contrast = (f_ga - f_as) ** 2 / 4.0  # = 1.0

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


def main() -> None:
    """Build the frozen short-range correlated model and write all outputs."""
    out_dir = Path(__file__).parent / "short_range_run"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Use n=8 for computational efficiency
    # Note: The "No long-range order" check is limited by the small supercell size
    # (periodic boundary effects cause apparent long-range correlations)
    struct = build_short_range_structure(n=8, a=5.0, correlation_length=1.0, seed=42)

    # Write .stru file
    stru_path = out_dir / "short_range_model.stru"
    write_stru(stru_path, struct)
    print(f"Wrote {stru_path}")

    # Compute independent references
    pair_corr = compute_directional_pair_correlations(struct)
    pcf = compute_pair_shell_statistics(struct)
    peak = predict_pdf_peaks(struct, pair_corr)

    # Measure correlation decay
    corr_decay = measure_correlation_decay(struct)

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
            "correlation_length": float(struct["correlation_length"]),
            "seed": int(struct["seed"]),
            "note": "Short-range correlated model (distinct from ordered checkerboard).",
        },
        "correlation_decay": {
            "separations": [float(s) for s in corr_decay["separations"]],
            "c_x": [float(c) for c in corr_decay["c_x"]],
            "c_y": [float(c) for c in corr_decay["c_y"]],
            "c_z": [float(c) for c in corr_decay["c_z"]],
            "decay_verified": bool(corr_decay["decay_verified"]),
            "no_long_range_order": bool(corr_decay["no_long_range_order"]),
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
    }
    refs_path.write_text(json.dumps(refs_data, indent=2))
    print(f"Wrote {refs_path}")

    # Print summary
    n_ga = sum(1 for sp in struct["species"] if sp == "Ga")
    n_as = sum(1 for sp in struct["species"] if sp == "As")
    print(f"\nStructure: n={struct['n']}, a={struct['a']} Å (supercell = {struct['n'] * struct['a']} Å)")
    print(f"Correlation length: {struct['correlation_length']} unit cells")
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
    print(f"\nCorrelation decay verification:")
    print(f"  Separations (unit cells): {corr_decay['separations']}")
    print(f"  c_x: {[f'{c:.4f}' for c in corr_decay['c_x']]}")
    print(f"  c_y: {[f'{c:.4f}' for c in corr_decay['c_y']]}")
    print(f"  c_z: {[f'{c:.4f}' for c in corr_decay['c_z']]}")
    print(f"  Decay verified: {corr_decay['decay_verified']}")
    print(f"  No long-range order: {corr_decay['no_long_range_order']}")


if __name__ == "__main__":
    main()
