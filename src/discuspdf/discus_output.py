"""Reader and Fourier transform for DISCUS_SUITE 6.23.01 single-crystal "common/nexus" HDF5 output.

Layout facts below were OBSERVED from real files written by DISCUS_SUITE 6.23.01
(``form nexus`` + ``value int|3DPDF``); see ``tests/fixtures/discus_6.23.01``:

* everything lives under ``/entry/data`` as datasets (no HDF5 attributes);
* intensity is ``data_values`` with axes ordered (abscissa, ordinate, top) = ``data_dimension``;
* the grid is given in *basecell-fractional* units (``coordinate_unit``), i.e. hkl in rlu, NOT Å^-1;
* ``data_corner`` is the hkl of element [0,0,0]; ``data_increment_vector[:, k]`` (as read by h5py,
  column k) is the hkl step along array axis k (verified with a deliberately skewed grid);
* ``data_type_with_bragg`` is ``bragg_present`` or ``bragg_subtracted``.

Nothing here converts hkl to Cartesian Q by relabelling: for non-orthogonal cells use
:func:`delta_pdf_fractional`, which transforms in fractional coordinates and keeps the metric.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import hashlib
import warnings
import numpy as np

_GROUP = "/entry/data"
_REQUIRED = ("data_values", "data_dimension", "data_axes", "data_corner", "data_increment_vector",
             "coordinate_unit", "data_type_axes", "data_type_reciprocal", "data_type_content",
             "data_type_with_bragg", "unit_cell_lengths", "unit_cell_angles",
             "audit_conform_dict_name")


def _text(value) -> str:
    return value.decode("utf-8", "replace") if isinstance(value, (bytes, np.bytes_)) else str(value)


def _validated_direct_basis(lengths, angles) -> np.ndarray:
    """Construct a cell basis without concealing degenerate/impossible metrics."""
    lengths = np.asarray(lengths, dtype=float)
    angles = np.asarray(angles, dtype=float)
    if lengths.shape != (3,) or not np.isfinite(lengths).all() or np.any(lengths <= 0):
        raise ValueError("unit_cell_lengths must be three finite positive values in Å")
    if angles.shape != (3,) or not np.isfinite(angles).all() or np.any(angles <= 0) or np.any(angles >= 180):
        raise ValueError("unit_cell_angles must be finite degrees strictly between 0 and 180")
    a, b, c = lengths
    alpha, beta, gamma = np.deg2rad(angles)
    ca, cb, cg = np.cos((alpha, beta, gamma))
    metric = np.array([[1, cg, cb], [cg, 1, ca], [cb, ca, 1]], dtype=float)
    eigenvalues = np.linalg.eigvalsh(metric)
    if eigenvalues[0] <= 1e-12 * eigenvalues[-1]:
        raise ValueError("Invalid or degenerate unit cell: metric is not positive definite")
    sg = np.sin(gamma)
    bx, by = b * cg, b * sg
    cx = c * cb
    cy = c * (ca - cb * cg) / sg
    cz_sq = c*c - cx*cx - cy*cy
    if not np.isfinite(cz_sq) or cz_sq <= 0:
        raise ValueError("Invalid unit cell: nonpositive squared c_z component")
    basis = np.array([[a, 0.0, 0.0], [bx, by, 0.0], [cx, cy, np.sqrt(cz_sq)]])
    if not np.isfinite(basis).all() or np.linalg.det(basis) <= 0:
        raise ValueError("Invalid unit cell basis")
    return basis


@dataclass(frozen=True)
class DiscusVolume:
    values: np.ndarray
    corner: np.ndarray                 # hkl of element [0,0,0]
    increments: np.ndarray             # increments[:, k] = hkl step along array axis k
    cell_lengths: np.ndarray           # Å
    cell_angles: np.ndarray            # degrees
    with_bragg: str
    meta: dict = field(default_factory=dict)
    provenance: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Validate physical cell and reciprocal-grid geometry, including direct construction."""
        v = np.asarray(self.values)
        if v.ndim != 3 or any(n < 1 for n in v.shape) or not np.issubdtype(v.dtype, np.number):
            raise ValueError("data_values must be a numeric 3D array")
        if not np.all(np.isfinite(v)):
            raise ValueError("data_values must be finite")
        shapes = {"data_corner": ((3,), self.corner),
                  "data_increment_vector": ((3, 3), self.increments),
                  "unit_cell_lengths": ((3,), self.cell_lengths),
                  "unit_cell_angles": ((3,), self.cell_angles)}
        for name, (shape, raw) in shapes.items():
            a = np.asarray(raw)
            if a.shape != shape or not np.issubdtype(a.dtype, np.number) or not np.isfinite(a).all():
                raise ValueError(f"{name} must be a finite numeric array with shape {shape}")
        inc = np.asarray(self.increments, dtype=float)
        if abs(np.linalg.det(inc)) < 1e-14 * (np.linalg.norm(inc, ord=2) ** 3) or np.linalg.cond(inc) > 1e10:
            raise ValueError("data_increment_vector is singular or nearly singular")
        if self.with_bragg not in ("bragg_present", "bragg_subtracted"):
            raise ValueError("Unknown Bragg metadata")
        _validated_direct_basis(self.cell_lengths, self.cell_angles)

    # -- geometry ---------------------------------------------------------
    def is_axis_aligned(self, tol: float = 1e-9) -> bool:
        off = self.increments - np.diag(np.diag(self.increments))
        return bool(np.all(np.abs(off) <= tol * np.max(np.abs(self.increments))))

    def hkl_axes(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        if not self.is_axis_aligned():
            raise ValueError("Grid is skewed in hkl (non-diagonal data_increment_vector); independent "
                             "1D hkl axes do not exist. Regrid explicitly before using this volume.")
        return tuple(self.corner[k] + np.diag(self.increments)[k] * np.arange(n)
                     for k, n in enumerate(self.values.shape))

    def direct_basis(self) -> np.ndarray:
        """Rows a, b, c in Cartesian Å; orientation is a convention."""
        return _validated_direct_basis(self.cell_lengths, self.cell_angles)

    def reciprocal_basis(self, two_pi: bool = True) -> np.ndarray:
        """Rows a*, b*, c*. Q (Å^-1) = h @ reciprocal_basis(two_pi=True); two_pi=False gives d* (Å^-1)."""
        b = np.linalg.inv(self.direct_basis()).T
        return 2 * np.pi * b if two_pi else b

    @property
    def cell_volume(self) -> float:
        return float(abs(np.linalg.det(self.direct_basis())))

    def is_orthogonal(self, tol: float = 1e-6) -> bool:
        return bool(np.allclose(self.cell_angles, 90.0, atol=tol))

    def to_reciprocal_volume(self):
        """Cartesian-Q view, only when this is exact (axis-aligned grid AND orthogonal cell)."""
        from .reciprocal import ReciprocalVolume
        if not (self.is_axis_aligned() and self.is_orthogonal()):
            raise ValueError("Exact Cartesian Q axes exist only for an axis-aligned grid in an orthogonal "
                             "cell. Use delta_pdf_fractional() instead of relabelling hkl as Q.")
        axes = tuple(2 * np.pi * h / a for h, a in zip(self.hkl_axes(), self.cell_lengths))
        pre = ("bragg_subtracted",) if self.with_bragg == "bragg_subtracted" else ()
        return ReciprocalVolume(self.values, axes, "cartesian_q", "1/angstrom", preprocessing=pre,
                                metadata={"source": self.provenance.get("path"),
                                          "q_from": "2*pi*hkl/a (orthogonal cell, computed here)"})


def read_discus_volume(path: str | Path, *, expect_content: str = "intensity") -> DiscusVolume:
    """Read a DISCUS 6.23.01 reciprocal-space (hkl) volume; fail with actionable messages otherwise."""
    import h5py
    p = Path(path)
    if not p.is_file() or not h5py.is_hdf5(p):
        raise ValueError(f"Not an HDF5 file: {p}")
    with h5py.File(p, "r") as f:
        if _GROUP not in f:
            raise ValueError(f"{p.name}: no {_GROUP} group; this is not DISCUS common/nexus output "
                             "(write it with 'form nexus' in the output menu).")
        g = f[_GROUP]
        missing = [k for k in _REQUIRED if k not in g]
        if missing:
            raise ValueError(f"{p.name}: missing datasets {missing}")
        raw = {k: g[k][()] for k in g if g[k].ndim <= 3}
    meta = {k: _text(v) for k, v in raw.items() if np.ndim(v) == 0}
    if meta["audit_conform_dict_name"] != "Disorder unified data":
        raise ValueError(f"Unrecognised file family {meta['audit_conform_dict_name']!r}")
    checks = {
        "coordinate_unit": ("basecell_fractional", "coordinates are not hkl in rlu; refusing to guess units"),
        "data_type_axes": ("hkl", "axes are not hkl (a real-space 'uvw' file is a 3D-PDF, not intensity)"),
        "data_type_reciprocal": ("reciprocal", "data are not in reciprocal space"),
        "data_type_content": (expect_content, f"content is not {expect_content!r}"),
    }
    for key, (want, why) in checks.items():
        if meta[key] != want:
            raise ValueError(f"{p.name}: {key}={meta[key]!r}, expected {want!r}: {why}.")
    if meta["data_type_with_bragg"] not in ("bragg_present", "bragg_subtracted"):
        raise ValueError(f"Unknown data_type_with_bragg={meta['data_type_with_bragg']!r}")
    values = np.asarray(raw["data_values"], dtype=float)
    if values.ndim != 3 or not np.isfinite(values).all():
        raise ValueError(f"{p.name}: data_values must be a finite 3D array")
    if tuple(np.asarray(raw["data_dimension"]).tolist()) != values.shape:
        raise ValueError(f"{p.name}: data_dimension {raw['data_dimension']} != array shape {values.shape}")
    if tuple(np.asarray(raw["data_axes"]).tolist()) != (1, 2, 3):
        raise ValueError(f"{p.name}: data_axes={raw['data_axes']}; only (1,2,3) was verified against DISCUS output")
    inc = np.asarray(raw["data_increment_vector"], dtype=float)
    if inc.shape != (3, 3):
        raise ValueError(f"{p.name}: data_increment_vector must be 3x3")
    prov = {"path": str(p), "sha256": hashlib.sha256(p.read_bytes()).hexdigest(),
            "creation_method": meta.get("audit_creation_method"),
            "metadata_asserted": {k: meta[k] for k in checks} | {"data_type_with_bragg": meta["data_type_with_bragg"]},
            "independently_checked": ["finite 3D array", "shape == data_dimension", "data_axes == (1,2,3)"]}
    return DiscusVolume(values, np.asarray(raw["data_corner"], float), inc,
                        np.asarray(raw["unit_cell_lengths"], float), np.asarray(raw["unit_cell_angles"], float),
                        meta["data_type_with_bragg"], meta, prov)


@dataclass(frozen=True)
class FractionalFFTResult:
    complex_volume: np.ndarray
    u_axes: tuple                      # fractional real-space coordinates (cell units) of each axis
    direct_basis: np.ndarray           # rows a,b,c (Å)
    imaginary_ratio: float
    friedel_relative_error: float
    normalization: str
    provenance: dict

    @property
    def real_volume(self) -> np.ndarray:
        return self.complex_volume.real

    def cartesian_length(self, u) -> float:
        """Length in Å of a fractional vector u, using the cell metric (not axis labels)."""
        return float(np.linalg.norm(np.asarray(u, float) @ self.direct_basis))


def _centered_steps(vol: DiscusVolume) -> np.ndarray:
    steps = np.diag(vol.increments)
    if not np.isfinite(steps).all() or np.any(steps <= 0):
        raise ValueError("Fractional FFT requires positive, increasing hkl grid steps")
    for k, n in enumerate(vol.values.shape):
        if not np.isclose(vol.corner[k], -(n // 2) * steps[k], rtol=1e-6, atol=1e-9):
            raise ValueError(f"hkl axis {k} is not centered on 0 at index n//2 (corner={vol.corner[k]}, "
                             f"step={steps[k]}, n={n}); the centered FFT would be mis-phased.")
    return steps


def delta_pdf_fractional(vol: DiscusVolume, *, require_diffuse: bool = True,
                         imaginary_tolerance: float = 1e-4,
                         raise_on_imaginary: bool = True) -> FractionalFFTResult:
    """Transform in fractional coordinates: rho(u) = dh dk dl / V_cell * sum_h I(h) exp(-2 pi i h.u).

    Equals the Cartesian formula dQx dQy dQz/(2 pi)^3 * sum I exp(-i Q.r) for any cell (no regridding,
    no interpolation). u spacing is 1/(n*dh) per axis. Raises for a skewed or uncentered grid and,
    by default, for Bragg-present data (the result would be a Patterson function, not a ΔPDF).
    """
    if require_diffuse and vol.with_bragg != "bragg_subtracted":
        raise ValueError("Volume is 'bragg_present'; a 3D-ΔPDF needs bragg_subtracted input "
                         "(or pass require_diffuse=False to obtain a Patterson function).")
    if not np.isfinite(imaginary_tolerance) or imaginary_tolerance < 0:
        raise ValueError("imaginary_tolerance must be finite and >= 0")
    vol.hkl_axes()                      # raises for skewed grids
    steps = _centered_steps(vol)
    I = vol.values
    c = np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(I))) * (np.prod(steps) / vol.cell_volume)
    u_axes = tuple(np.fft.fftshift(np.fft.fftfreq(n, d=d)) for n, d in zip(I.shape, steps))
    ratio = float(np.linalg.norm(c.imag) / max(np.linalg.norm(c.real), 1e-30))
    un = np.fft.ifftshift(I)
    opp = un[np.ix_(*[(-np.arange(n)) % n for n in I.shape])]
    fr = float(np.linalg.norm(un - opp) / max(np.linalg.norm(un), 1e-30))
    result = FractionalFFTResult(c, u_axes, vol.direct_basis(), ratio, fr,
                                 "dh*dk*dl/V_cell; exp(-2*pi*i*h.u); u step = 1/(n*dh)", dict(vol.provenance))
    if ratio > imaginary_tolerance:
        from .volume import FFTValidationError
        message = (f"imaginary/real ratio {ratio:.3g} > {imaginary_tolerance:g}; "
                   f"Friedel error {fr:.3g}")
        if raise_on_imaginary:
            raise FFTValidationError(message, result)
        warnings.warn(message, RuntimeWarning, stacklevel=2)
    return result


def evaluate_real_space(vol: DiscusVolume, u_points) -> np.ndarray:
    """Exact (not gridded) rho at fractional points u (m,3); same normalization as delta_pdf_fractional."""
    h = vol.hkl_axes(); steps = np.diag(vol.increments)
    u = np.atleast_2d(np.asarray(u_points, float))
    E = [np.exp(-2j * np.pi * np.outer(u[:, k], h[k])) for k in range(3)]
    rho = np.einsum("pi,pj,pk,ijk->p", *E, vol.values, optimize=True)
    return rho * (np.prod(steps) / vol.cell_volume)
