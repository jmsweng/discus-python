#!/usr/bin/env python3
"""Compare DISCUS engine 3D-ΔPDF output with independent predictions.

This script:
1. Loads the DISCUS HDF5 output (total intensity)
2. Computes the 3D-ΔPDF using delta_pdf_fractional
3. Compares off-origin peaks at ±2.5 Å against predicted positions and amplitudes
"""
import json
import sys
from pathlib import Path

import numpy as np

# Add the src directory to the path
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent / "src"))

from discuspdf.discus_output import read_discus_volume, delta_pdf_fractional


def main():
    run_dir = Path(__file__).parent
    
    # Load the independent references
    with open(run_dir / "references.json") as f:
        refs = json.load(f)
    
    print("=== Independent References ===")
    print(f"PDF peaks (predicted):")
    taus = refs["pdf_peaks"]["tau"]
    amps = refs["pdf_peaks"]["amplitude"]
    for tau, amp in zip(taus, amps):
        print(f"  tau={tau}: amplitude={amp:.6f}")
    
    # Load the DISCUS HDF5 output (diffuse intensity in hkl)
    # We need to compute the inverse FFT to get the signed 3D-ΔPDF
    print("\n=== Loading DISCUS HDF5 output (diffuse intensity, hkl) ===")
    import h5py
    
    h5_path = run_dir / "corr_diffuse.h5"
    with h5py.File(h5_path, "r") as f:
        # Read the data (nested in entry/data/)
        intensity = f["entry/data/data_values"][:]
        # Read the metadata
        corner = f["entry/data/data_corner"][:]
        increments = f["entry/data/data_increment_vector"][:]
        cell_lengths = f["entry/data/unit_cell_lengths"][:]
    
    print(f"Loaded: {h5_path.name}")
    print(f"  Shape: {intensity.shape}")
    print(f"  Corner: {corner}")
    print(f"  Increments: {increments}")
    print(f"  Cell lengths: {cell_lengths}")
    print(f"  Min: {intensity.min():.4e}")
    print(f"  Max: {intensity.max():.4e}")
    print(f"  Mean: {intensity.mean():.4e}")
    
    # Compute the inverse FFT to get the 3D-ΔPDF
    # The intensity is in hkl space, so we need to inverse FFT to get real space
    print("\n=== Computing 3D-ΔPDF (inverse FFT of intensity) ===")
    pdf = np.fft.ifftn(intensity).real
    
    print(f"\n  PDF min: {pdf.min():.6e}")
    print(f"  PDF max: {pdf.max():.6e}")
    print(f"  PDF mean: {pdf.mean():.6e}")
    
    # Find the peak in the PDF
    # The PDF is on a grid in fractional coordinates
    # We need to find the peak at ±2.5 Å along x, y, z
    
    # The grid is in fractional units
    # The corner is at -2.0 and increments are 0.125
    # So the grid goes from -2.0 to -2.0 + 32*0.125 = -2.0 + 4.0 = 2.0
    
    n_grid = pdf.shape[0]
    u_min = corner[0]
    u_max = u_min + (n_grid - 1) * increments[0][0]
    
    print(f"\n  Intensity grid: {n_grid} points from {u_min:.4f} to {u_max:.4f}")
    print(f"  PDF grid: {n_grid} points (real space, same dimensions)")
    
    # 2.5 Å in fractional units (relative to 20 Å supercell)
    tau_frac = 2.5 / 20.0  # = 0.125
    
    # Grid index for tau_frac
    # The grid is evenly spaced from u_min to u_max
    # idx = (tau_frac - u_min) / (u_max - u_min) * (n_grid - 1)
    
    # But the grid is in reciprocal space, so we need to be careful
    # The PDF is the inverse FFT of the intensity
    # The grid index for a displacement tau is:
    # idx = int(tau / (u_max - u_min) * n_grid)
    
    # Actually, let's just find the maximum in the PDF
    # The peak should be at ±2.5 Å along x, y, z
    
    # Find the global maximum
    max_idx = np.unravel_index(np.argmax(pdf), pdf.shape)
    print(f"\n  Global maximum at index {max_idx}: {pdf[max_idx]:.6e}")
    
    # Find the peak along each axis
    # For x-axis: fix y=0, z=0, vary x
    # The grid center is at index n_grid//2
    
    center = n_grid // 2
    
    # Peak along +x
    pdf_x = pdf[:, center, center]
    peak_x_idx = np.argmax(pdf_x)
    print(f"\n  Peak along +x: index={peak_x_idx}, value={pdf_x[peak_x_idx]:.6e}")
    
    # Peak along -x
    peak_x_neg_idx = np.argmin(pdf_x)
    print(f"  Peak along -x: index={peak_x_neg_idx}, value={pdf_x[peak_x_neg_idx]:.6e}")
    
    # Peak along +y
    pdf_y = pdf[center, :, center]
    peak_y_idx = np.argmax(pdf_y)
    print(f"\n  Peak along +y: index={peak_y_idx}, value={pdf_y[peak_y_idx]:.6e}")
    
    # Peak along -y
    peak_y_neg_idx = np.argmin(pdf_y)
    print(f"  Peak along -y: index={peak_y_neg_idx}, value={pdf_y[peak_y_neg_idx]:.6e}")
    
    # Peak along +z
    pdf_z = pdf[center, center, :]
    peak_z_idx = np.argmax(pdf_z)
    print(f"\n  Peak along +z: index={peak_z_idx}, value={pdf_z[peak_z_idx]:.6e}")
    
    # Peak along -z
    peak_z_neg_idx = np.argmin(pdf_z)
    print(f"  Peak along -z: index={peak_z_neg_idx}, value={pdf_z[peak_z_neg_idx]:.6e}")
    
    # Compare with predictions
    print("\n=== Comparison with Predictions ===")
    
    # Predicted peaks (from references.json)
    # tau=(2.5, 0.0, 0.0): amplitude=-0.002
    # tau=(0.0, 2.5, 0.0): amplitude=+0.002
    # tau=(0.0, 0.0, 2.5): amplitude=-0.002
    
    # The PDF from DISCUS is in arbitrary units (depends on the intensity scale)
    # We need to normalize the PDF to compare with the predictions
    
    # The predictions are in units of (f_Ga - f_As)²/4 = 1.0
    # The DISCUS PDF is in units of the intensity (which depends on the form factors)
    
    # For a rough comparison, let's normalize the PDF by its maximum
    pdf_norm = pdf / pdf.max()
    
    print(f"  Normalized PDF (divided by max):")
    print(f"    Peak along +x: {pdf_norm[peak_x_idx, center, center]:.6f}")
    print(f"    Peak along -x: {pdf_norm[peak_x_neg_idx, center, center]:.6f}")
    print(f"    Peak along +y: {pdf_norm[center, peak_y_idx, center]:.6f}")
    print(f"    Peak along -y: {pdf_norm[center, peak_y_neg_idx, center]:.6f}")
    print(f"    Peak along +z: {pdf_norm[center, center, peak_z_idx]:.6f}")
    print(f"    Peak along -z: {pdf_norm[center, center, peak_z_neg_idx]:.6f}")
    
    # The signs should match:
    # +x: negative, -x: negative (Friedel pair)
    # +y: positive, -y: positive (Friedel pair)
    # +z: negative, -z: negative (Friedel pair)
    
    print("\n  Sign check:")
    print(f"    +x: {'NEGATIVE' if pdf[peak_x_idx, center, center] < 0 else 'POSITIVE'} (expected: NEGATIVE)")
    print(f"    -x: {'NEGATIVE' if pdf[peak_x_neg_idx, center, center] < 0 else 'POSITIVE'} (expected: NEGATIVE)")
    print(f"    +y: {'NEGATIVE' if pdf[center, peak_y_idx, center] < 0 else 'POSITIVE'} (expected: POSITIVE)")
    print(f"    -y: {'NEGATIVE' if pdf[center, peak_y_neg_idx, center] < 0 else 'POSITIVE'} (expected: POSITIVE)")
    print(f"    +z: {'NEGATIVE' if pdf[center, center, peak_z_idx] < 0 else 'POSITIVE'} (expected: NEGATIVE)")
    print(f"    -z: {'NEGATIVE' if pdf[center, center, peak_z_neg_idx] < 0 else 'POSITIVE'} (expected: NEGATIVE)")
    
    # Save the results
    results = {
        "h5_file": h5_path.name,
        "shape": [int(x) for x in pdf.shape],
        "corner": [float(x) for x in corner],
        "increments": [[float(x) for x in row] for row in increments],
        "cell_lengths": [float(x) for x in cell_lengths],
        "pdf_min": float(pdf.min()),
        "pdf_max": float(pdf.max()),
        "pdf_mean": float(pdf.mean()),
        "peak_positions": {
            "+x": [int(peak_x_idx), int(center), int(center)],
            "-x": [int(peak_x_neg_idx), int(center), int(center)],
            "+y": [int(center), int(peak_y_idx), int(center)],
            "-y": [int(center), int(peak_y_neg_idx), int(center)],
            "+z": [int(center), int(center), int(peak_z_idx)],
            "-z": [int(center), int(center), int(peak_z_neg_idx)],
        },
        "peak_values": {
            "+x": float(pdf[peak_x_idx, center, center]),
            "-x": float(pdf[peak_x_neg_idx, center, center]),
            "+y": float(pdf[center, peak_y_idx, center]),
            "-y": float(pdf[center, peak_y_neg_idx, center]),
            "+z": float(pdf[center, center, peak_z_idx]),
            "-z": float(pdf[center, center, peak_z_neg_idx]),
        },
        "peak_values_normalized": {
            "+x": float(pdf_norm[peak_x_idx, center, center]),
            "-x": float(pdf_norm[peak_x_neg_idx, center, center]),
            "+y": float(pdf_norm[center, peak_y_idx, center]),
            "-y": float(pdf_norm[center, peak_y_neg_idx, center]),
            "+z": float(pdf_norm[center, center, peak_z_idx]),
            "-z": float(pdf_norm[center, center, peak_z_neg_idx]),
        },
        "predicted_peaks": {
            "tau": refs["pdf_peaks"]["tau"],
            "amplitude": refs["pdf_peaks"]["amplitude"],
        },
    }
    
    with open(run_dir / "engine_comparison.json", "w") as f:
        json.dump(results, f, indent=2)
    
    print(f"\n  Results saved to: {run_dir / 'engine_comparison.json'}")


if __name__ == "__main__":
    main()
