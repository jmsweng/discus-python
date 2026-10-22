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


from dataclasses import dataclass


@dataclass(frozen=True)
class FFTDiagnostics:
    imaginary_ratio: float
    max_imaginary: float
    max_real: float
    friedel_relative_error: float
    dq: tuple[float, float, float]
    reciprocal_voxel_volume: float
    normalization: str


@dataclass(frozen=True)
class FFTResult:
    complex_volume: np.ndarray
    r_axes: tuple[np.ndarray, np.ndarray, np.ndarray]
    diagnostics: FFTDiagnostics

    @property
    def real_volume(self) -> np.ndarray:
        return self.complex_volume.real


def delta_pdf_from_diffuse(
    diffuse_intensity: np.ndarray,
    q_axes: Sequence[np.ndarray] | None = None,
    *,
    imaginary_tolerance: float = 1e-4,
    raise_on_imaginary: bool = True,
) -> FFTResult:
    """Transform a Bragg-subtracted diffuse volume to complex real-space data.

    Inputs must be finite and sampled on centered uniform Cartesian Q axes
    expressed in Å^-1 (rad/Å). The FFT uses exp(-i Q.r) with
    1/(2π)^3 times the reciprocal voxel volume normalization.

    When supplied a ReciprocalVolume, requires explicit Bragg preprocessing
    metadata. Raw arrays are accepted for backwards-compatible analysis, but
    their units/preprocessing cannot be verified automatically.
    """
    from .reciprocal import ReciprocalVolume
    if isinstance(diffuse_intensity, ReciprocalVolume):
        volume_input = diffuse_intensity
        volume_input.require_bragg_subtracted()
        if q_axes is not None:
            raise ValueError("Do not supply q_axes separately with ReciprocalVolume")
        q_axes = volume_input.q_axes
        intensity = np.asarray(volume_input.intensity, dtype=float)
    else:
        intensity = np.asarray(diffuse_intensity, dtype=float)
    if not np.isfinite(imaginary_tolerance) or imaginary_tolerance < 0:
        raise ValueError("imaginary_tolerance must be finite and >= 0")
    if intensity.ndim != 3 or q_axes is None or len(q_axes) != 3:
        raise ValueError("Provide 3D intensity and three centered Cartesian Q axes")
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

    voxel = float(np.prod(spacings))
    complex_volume = np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(intensity)))
    complex_volume = complex_volume * (voxel / (2*np.pi)**3)
    real_norm = float(np.linalg.norm(complex_volume.real))
    imag_norm = float(np.linalg.norm(complex_volume.imag))
    ratio = imag_norm / max(real_norm, 1e-30)

    # I(Q) and I(-Q) are related through periodic FFT-index inversion.
    # This diagnostic assumes the discrete periodic representation and is
    # meaningful independently of array parity.
    unshifted = np.fft.ifftshift(intensity)
    negate = tuple((-np.arange(n)) % n for n in intensity.shape)
    opposite = unshifted[np.ix_(*negate)]
    friedel = float(np.linalg.norm(unshifted - opposite) /
                    max(np.linalg.norm(unshifted), 1e-30))
    diagnostics = FFTDiagnostics(
        imaginary_ratio=ratio,
        max_imaginary=float(np.max(np.abs(complex_volume.imag))),
        max_real=float(np.max(np.abs(complex_volume.real))),
        friedel_relative_error=friedel,
        dq=tuple(spacings),
        reciprocal_voxel_volume=voxel,
        normalization="dQx*dQy*dQz/(2*pi)^3; exp(-i Q.r)",
    )
    if ratio > imaginary_tolerance:
        import warnings
        message = (f"FFT imaginary/real norm ratio {ratio:.3g} exceeds "
                   f"tolerance {imaginary_tolerance:.3g}; Friedel error {friedel:.3g}")
        if raise_on_imaginary:
            raise ValueError(message)
        warnings.warn(message, RuntimeWarning, stacklevel=2)
    return FFTResult(complex_volume, tuple(r_axes), diagnostics)


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
