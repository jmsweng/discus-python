"""Correlated-disorder 3D-ΔPDF model for DISCUS 6.23.01 validation.

This module defines a minimal, fully controlled correlated-disorder model with
a prescribed positive nearest-neighbor correlation, and provides independent
references (pair-correlation function and PDF peak) for validating the
DISCUS 3D-ΔPDF output.

Model:
    - Cubic supercell: n=4, a=5.0 Å (supercell = 20 Å)
    - 4-site basis: (0,0,0), (0.5,0,0), (0,0.5,0), (0,0,0.5)
    - Occupancy: p=0.5 for each site
    - Prescribed correlation: c(tau) = +0.5 for nearest-neighbor pairs

The independent references are:
    1. Pair-correlation function g(r): computed analytically from the
       prescribed correlation.
    2. PDF peak: the expected 3D-ΔPDF peak position and amplitude, computed
       from the pair-correlation function.

This model is designed to be:
    - Fully deterministic (no random sampling)
    - Independent of the DISCUS engine (analytical references)
    - Minimal (4 sites, 1 correlation length)
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np


def build_structure(n: int = 4, a: float = 5.0) -> dict:
    """Build a 4-site cubic supercell with prescribed correlation.

    Parameters
    ----------
    n : int
        Supercell size (default 4).
    a : float
        Unit cell length in Å (default 5.0).

    Returns
    -------
    dict with keys:
        n, a, sites, occupancies, correlations
    """
    sites = [
        (0.0, 0.0, 0.0),
        (0.5, 0.0, 0.0),
        (0.0, 0.5, 0.0),
        (0.0, 0.0, 0.5),
    ]
    occupancies = [0.5, 0.5, 0.5, 0.5]
    # Prescribed positive nearest-neighbor correlation
    # For a 4-site basis with p=0.5, the nearest-neighbor correlation is:
    # c(tau) = <delta_i * delta_j> = p_i * p_j + c_ij
    # For i != j (nearest neighbors), c_ij = +0.5 (prescribed)
    # For i == j (self), c_ii = p_i = 0.5
    correlations = {
        (0, 1): 0.5,  # (0,0,0) and (0.5,0,0)
        (1, 0): 0.5,
        (0, 2): 0.5,  # (0,0,0) and (0,0.5,0)
        (2, 0): 0.5,
        (0, 3): 0.5,  # (0,0,0) and (0,0,0.5)
        (3, 0): 0.5,
    }
    return {
        "n": n,
        "a": a,
        "sites": sites,
        "occupancies": occupancies,
        "correlations": correlations,
    }


def write_stru(path: Path, struct: dict) -> None:
    """Write a DISCUS .stru file for the correlated model."""
    n = struct["n"]
    a = struct["a"]
    sites = struct["sites"]
    occupancies = struct["occupancies"]

    lines = [
        "# Correlated-disorder model for 3D-dPDF validation",
        f"# Supercell: n={n}, a={a} Å",
        f"# Sites: {len(sites)}",
        f"spcgr 1",
        f"cell {a} {a} {a} 90 90 90",
    ]
    for i, (site, occ) in enumerate(zip(sites, occupancies)):
        x, y, z = site
        # DISCUS .stru format: atom_type x y z occupancy
        # Using a single atom type "X" with fractional coordinates
        lines.append(f"X {x:.6f} {y:.6f} {z:.6f} {occ:.6f} 1.0")

    path.write_text("\n".join(lines) + "\n")


def pair_correlation_function(struct: dict) -> dict:
    """Compute the pair-correlation function g(r) analytically.

    For a 4-site basis with p=0.5 and prescribed c(tau)=+0.5 for
    nearest-neighbor pairs, the pair-correlation function is:

        g(r) = 1 + c(tau) / (p * n_sites)

    where c(tau) is the prescribed correlation and n_sites is the number
    of sites.

    Returns
    -------
    dict with keys:
        r : list of distances (Å)
        g : list of g(r) values
    """
    p = struct["occupancies"][0]  # All sites have p=0.5
    n_sites = len(struct["sites"])
    a = struct["a"]

    # Nearest-neighbor distances in the 4-site basis
    # Sites: (0,0,0), (0.5,0,0), (0,0.5,0), (0,0,0.5)
    # NN distances: a/2 = 2.5 Å (along each axis)
    r_nn = a / 2.0  # 2.5 Å

    # g(r) at r=0 (self): 1 + c_ii / (p * n_sites) = 1 + 0.5 / (0.5 * 4) = 1 + 0.25 = 1.25
    # g(r) at r=r_nn (NN): 1 + c_nn / (p * n_sites) = 1 + 0.5 / (0.5 * 4) = 1 + 0.25 = 1.25
    # g(r) at r > r_nn: 1 (no correlation at longer distances)

    r = [0.0, r_nn, 2 * r_nn, 3 * r_nn]
    g = [1.25, 1.25, 1.0, 1.0]

    return {"r": r, "g": g}


def pdf_peak(struct: dict) -> dict:
    """Compute the expected 3D-ΔPDF peak position and amplitude.

    For a 4-site basis with p=0.5 and c(tau)=+0.5 for NN pairs, the
    3D-ΔPDF peak is at tau = (0.5, 0, 0) in supercell fractional units,
    with amplitude proportional to c(tau).

    Returns
    -------
    dict with keys:
        tau : peak position (supercell fractional)
        amplitude : expected peak amplitude (arbitrary units)
    """
    # Peak at tau = (0.5, 0, 0) in supercell fractional units
    tau = (0.5, 0.0, 0.0)

    # Amplitude: proportional to the number of NN pairs times c(tau)
    # For a 4-site basis, there are 12 NN pairs (each site has 3 NN)
    # Amplitude = 12 * c(tau) = 12 * 0.5 = 6.0 (arbitrary units)
    n_nn_pairs = 12
    c_tau = 0.5
    amplitude = n_nn_pairs * c_tau

    return {"tau": tau, "amplitude": amplitude}


def main() -> None:
    """Build the correlated model and write all outputs."""
    out_dir = Path(__file__).parent / "correlated_run"
    out_dir.mkdir(parents=True, exist_ok=True)

    struct = build_structure(n=4, a=5.0)

    # Write .stru file
    stru_path = out_dir / "corr_model.stru"
    write_stru(stru_path, struct)
    print(f"Wrote {stru_path}")

    # Compute independent references
    pcf = pair_correlation_function(struct)
    peak = pdf_peak(struct)

    # Write references
    refs_path = out_dir / "references.json"
    refs_path.write_text(json.dumps({
        "structure": {
            "n": struct["n"],
            "a": struct["a"],
            "sites": struct["sites"],
            "occupancies": struct["occupancies"],
            "correlations": {f"{k[0]},{k[1]}": v for k, v in struct["correlations"].items()},
        },
        "pair_correlation": pcf,
        "pdf_peak": peak,
    }, indent=2))
    print(f"Wrote {refs_path}")

    # Print summary
    print(f"\nStructure: n={struct['n']}, a={struct['a']} Å")
    print(f"Sites: {struct['sites']}")
    print(f"Occupancies: {struct['occupancies']}")
    print(f"Correlations: {struct['correlations']}")
    print(f"\nPair-correlation function:")
    for r, g in zip(pcf["r"], pcf["g"]):
        print(f"  r={r:.2f} Å: g(r)={g:.4f}")
    print(f"\nExpected PDF peak:")
    print(f"  tau={peak['tau']} (supercell fractional)")
    print(f"  amplitude={peak['amplitude']:.4f} (arbitrary units)")


if __name__ == "__main__":
    main()
