#!/usr/bin/env python3
"""Compare DISCUS engine 3D-ΔPDF output with independent predictions for the short-range model."""
import json
import sys
from pathlib import Path

import numpy as np
import h5py


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
    
    # Load the DISCUS HDF5 output (total intensity in hkl)
    print("\n=== Loading DISCUS HDF5 output (total intensity, hkl) ===")
    h5_path = run_dir / "sr_total.h5"
    with h5py.File(h5_path, "r") as f:
        intensity = f["entry/data/data_values"][:]
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
    print("\n=== Computing 3D-ΔPDF (inverse FFT of intensity) ===")
    pdf = np.fft.ifftn(intensity).real
    
    print(f"\n  PDF min: {pdf.min():.6e}")
    print(f"  PDF max: {pdf.max():.6e}")
    print(f"  PDF mean: {pdf.mean():.6e}")
    
    # Find the peaks along each axis
    n_grid = pdf.shape[0]
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
    
    # The PDF from DISCUS is in arbitrary units (depends on the intensity scale)
    # We need to normalize the PDF to compare with the predictions
    
    # For a rough comparison, let's normalize the PDF by its maximum
    if pdf.max() != 0:
        pdf_norm = pdf / pdf.max()
    else:
        pdf_norm = pdf
    
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
    # +z: positive, -z: positive (Friedel pair)
    
    print("\n  Sign check:")
    print(f"    +x: {'NEGATIVE' if pdf[peak_x_idx, center, center] < 0 else 'POSITIVE'} (expected: NEGATIVE)")
    print(f"    -x: {'NEGATIVE' if pdf[peak_x_neg_idx, center, center] < 0 else 'POSITIVE'} (expected: NEGATIVE)")
    print(f"    +y: {'NEGATIVE' if pdf[center, peak_y_idx, center] < 0 else 'POSITIVE'} (expected: POSITIVE)")
    print(f"    -y: {'NEGATIVE' if pdf[center, peak_y_neg_idx, center] < 0 else 'POSITIVE'} (expected: POSITIVE)")
    print(f"    +z: {'NEGATIVE' if pdf[center, center, peak_z_idx] < 0 else 'POSITIVE'} (expected: POSITIVE)")
    print(f"    -z: {'NEGATIVE' if pdf[center, center, peak_z_neg_idx] < 0 else 'POSITIVE'} (expected: POSITIVE)")
    
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
