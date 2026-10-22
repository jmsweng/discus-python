"""3D-Delta PDF analysis helpers. Not a replacement for diffuse-scattering calculation."""
from __future__ import annotations

from pathlib import Path
from typing import Sequence
import numpy as np


def read_hdf5_volume(
    path: str | Path,
    intensity_dataset: str,
    q_datasets: Sequence[str],
) -> tuple[np.ndarray, tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """Read an HDF5 reciprocal-space volume with explicitly selected datasets.

    Each coordinate dataset must be a 1D axis in physical Q units (Å^-1).
    Inspect the HDF5 file before choosing paths; layouts vary by source/version.
    """
    import h5py
    if len(q_datasets) != 3:
        raise ValueError("Specify three HDF5 coordinate dataset paths")
    with h5py.File(path, "r") as file:
        intensity = np.asarray(file[intensity_dataset], dtype=float)
        q_axes = tuple(np.asarray(file[key], dtype=float) for key in q_datasets)
    return intensity, q_axes


def delta_pdf_from_diffuse(
    diffuse_intensity: np.ndarray,
    q_axes: Sequence[np.ndarray],
) -> tuple[np.ndarray, tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """Illustrative 3D Fourier transform of BRAGG-SUBTRACTED diffuse intensity.

    Requires a finite cubic/rectangular reciprocal-space lattice sampled on
    *uniform centered Cartesian* Q axes (Å^-1), not raw h,k,l axes. Returns
    the real component of 1/(2π)^3 ∫ I_diff(Q) exp(-i Q.r) dQ on r axes (Å).
    Absolute experimental normalization, form factors, masking, apodization,
    missing voxels and background corrections are not handled here.
    """
    intensity = np.asarray(diffuse_intensity, dtype=float)
    if intensity.ndim != 3 or len(q_axes) != 3:
        raise ValueError("Provide a 3D intensity array and three 1D Q axes")
    if not np.all(np.isfinite(intensity)):
        raise ValueError("Intensity contains missing/NaN/Inf values; preprocess first")

    spacings = []
    r_axes = []
    for axis_index, (n, raw) in enumerate(zip(intensity.shape, q_axes)):
        q = np.asarray(raw, dtype=float)
        if q.ndim != 1 or len(q) != n or n < 2 or not np.all(np.isfinite(q)):
            raise ValueError(f"Q axis {axis_index} must be finite, 1D, length {n}")
        diffs = np.diff(q)
        dq = float(np.mean(diffs))
        if dq <= 0 or not np.allclose(diffs, dq, rtol=1e-5, atol=1e-10):
            raise ValueError(f"Q axis {axis_index} must be uniformly increasing")
        expected = (np.arange(n) - (n // 2)) * dq
        if not np.allclose(q, expected, rtol=1e-5, atol=max(1e-8, 1e-5*dq)):
            raise ValueError(f"Q axis {axis_index} must be centered on zero")
        spacings.append(dq)
        r_axes.append(np.fft.fftshift(np.fft.fftfreq(n, d=dq)) * 2*np.pi)

    # fftn implements exp(-i Q.r) when Q-grid's origin has been shifted to index 0.
    volume = np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(intensity)))
    volume = volume * (np.prod(spacings) / (2*np.pi)**3)
    # A physical centrosymmetric I(Q) should give a nearly real ΔPDF.
    imag_ratio = np.linalg.norm(volume.imag) / max(np.linalg.norm(volume.real), 1e-30)
    if imag_ratio > 1e-4:
        import warnings
        warnings.warn("Large imaginary component: check Friedel symmetry, background, and Q axes")
    return volume.real, tuple(r_axes)


def plot_delta_pdf_slice(volume: np.ndarray, r_axes: Sequence[np.ndarray], axis: int = 2, index: int | None = None):
    """Plot a central 2D slice of a 3D ΔPDF volume."""
    import matplotlib.pyplot as plt
    volume = np.asarray(volume)
    if volume.ndim != 3 or axis not in (0, 1, 2) or len(r_axes) != 3:
        raise ValueError("Expected a 3D volume and axis 0, 1, or 2")
    if index is None:
        index = volume.shape[axis] // 2
    if not 0 <= index < volume.shape[axis]:
        raise ValueError("slice index outside volume")
    coords = [i for i in range(3) if i != axis]
    sliced = np.take(volume, index, axis=axis).T
    a, b = np.asarray(r_axes[coords[0]]), np.asarray(r_axes[coords[1]])
    vmax = float(np.max(np.abs(sliced))) or 1.0
    fig, ax = plt.subplots(figsize=(6, 5))
    image = ax.imshow(sliced, origin="lower", aspect="equal", cmap="RdBu_r",
                      vmin=-vmax, vmax=vmax,
                      extent=(a[0], a[-1], b[0], b[-1]))
    ax.set(xlabel=f"r{coords[0]+1} (Å)", ylabel=f"r{coords[1]+1} (Å)")
    fig.colorbar(image, ax=ax, label="ΔPDF (relative units)")
    return ax
