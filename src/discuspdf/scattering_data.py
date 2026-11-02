"""Labelled reciprocal-space scattering data (Package B, Features 7 and 19).

:class:`ScatteringVolume` couples values on a :class:`~discuspdf.reciprocal_grid.ReciprocalGrid` with
their uncertainty, an explicit mask, units, the Bragg treatment and provenance. It is a plain NumPy
container; :meth:`ScatteringVolume.to_xarray` / :meth:`ScatteringVolume.from_xarray` convert to and
from ``xarray`` (optional dependency), and :meth:`ScatteringVolume.save_hdf5` /
:meth:`ScatteringVolume.load_hdf5` give a versioned HDF5 file (schema ``discuspdf.ScatteringVolume``).

Rules that are enforced, not just documented:

* ``mask`` is boolean, **True = excluded**. Missing data must be masked: NaN or Inf in ``intensity``
  (or ``variance``) is accepted only on masked voxels.
* ``variance`` is the variance (σ²) of ``intensity`` in ``units²``; it must be >= 0 where unmasked.
* ``bragg`` is one of ``bragg_present`` (total scattering), ``bragg_subtracted`` (e.g. DISCUS
  ``set aver``), ``bragg_masked`` (Bragg regions removed by the mask; needs a mask) or ``unknown``.
  A label is a declaration by the source; it is not a proof that the correction was right.
* Negative intensities are kept (Bragg-subtracted data can be negative); nothing is clipped.
* Every transform to real space refuses masked voxels unless an explicit fill value is given, and
  records that fill in ``corrections``.
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any
import datetime as _dt
import hashlib
import json

import numpy as np

from .reciprocal_grid import ReciprocalGrid

SCHEMA_NAME = "discuspdf.ScatteringVolume"
SCHEMA_VERSION = 1
BRAGG_LABELS = ("bragg_present", "bragg_subtracted", "bragg_masked", "unknown")
_DIMS = ("i0", "i1", "i2")


def _json_safe(name: str, value) -> str:
    try:
        return json.dumps(value, sort_keys=True, allow_nan=False)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be JSON-serialisable (str, numbers, lists, dicts; no NaN)") from exc


def _readonly(a: np.ndarray) -> np.ndarray:
    a = np.array(a, copy=True)
    a.flags.writeable = False
    return a


def _digest(a: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()


def _package_version() -> str:
    try:
        from importlib.metadata import version
        return version("discuspdf-starter")
    except Exception:  # pragma: no cover - not installed as a distribution
        return "unknown"


def mask_near_bragg(grid: ReciprocalGrid, radius: float) -> np.ndarray:
    """Boolean mask, True where |Q − Q(nearest integer hkl)| <= ``radius`` (Å⁻¹, 2π convention).

    The nearest reflection is found among the 27 integer points around ``round(hkl)``, so the
    Euclidean Q metric is used even for oblique reciprocal lattices. ``radius`` must be smaller than
    half the shortest reciprocal-lattice vector, otherwise masks of neighbouring reflections merge
    and the "nearest" choice is no longer meaningful."""
    radius = float(radius)
    B = 2 * np.pi * grid.lattice.reciprocal_matrix
    shortest = min(np.linalg.norm(v) for v in (*B, B[0] + B[1], B[0] - B[1], B[1] + B[2], B[1] - B[2],
                                                 B[0] + B[2], B[0] - B[2]))
    if not np.isfinite(radius) or radius <= 0 or radius >= shortest / 2:
        raise ValueError(f"radius must be in (0, {shortest / 2:.6g}) Å⁻¹ for this lattice")
    hkl = grid.hkl()
    base = np.round(hkl)
    best = np.full(grid.shape, np.inf)
    for off in np.ndindex(3, 3, 3):
        d = (hkl - (base + np.asarray(off) - 1)) @ B
        best = np.minimum(best, np.einsum("...i,...i->...", d, d))
    return best <= radius * radius


@dataclass(frozen=True, eq=False)
class ScatteringVolume:
    """Values on a reciprocal grid with uncertainty, mask, units, Bragg label and provenance."""

    grid: ReciprocalGrid
    intensity: np.ndarray
    units: str
    bragg: str = "unknown"
    variance: np.ndarray | None = None
    mask: np.ndarray | None = None
    mask_reasons: tuple[str, ...] = ()
    corrections: tuple[str, ...] = ()
    provenance: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not isinstance(self.grid, ReciprocalGrid):
            raise TypeError("grid must be a ReciprocalGrid")
        shape = self.grid.shape
        I = np.asarray(self.intensity)
        if not np.issubdtype(I.dtype, np.floating) and not np.issubdtype(I.dtype, np.integer):
            raise ValueError(f"intensity must be real numeric, got dtype {I.dtype}")
        if I.shape != shape:
            raise ValueError(f"intensity shape {I.shape} != grid shape {shape}")
        I = I.astype(float)
        if not isinstance(self.units, str) or not self.units.strip():
            raise ValueError("units must be a non-empty string (e.g. 'arbitrary', 'e^2 per crystal')")
        if self.bragg not in BRAGG_LABELS:
            raise ValueError(f"bragg must be one of {BRAGG_LABELS}, got {self.bragg!r}")
        M = None
        if self.mask is not None:
            M = np.asarray(self.mask)
            if M.dtype != bool:
                raise ValueError(f"mask must be a boolean array (True = excluded), got dtype {M.dtype}")
            if M.shape != shape:
                raise ValueError(f"mask shape {M.shape} != grid shape {shape}")
        unmasked = np.ones(shape, bool) if M is None else ~M
        if not np.isfinite(I[unmasked]).all():
            raise ValueError("intensity has NaN/Inf on unmasked voxels; mask missing data explicitly")
        V = None
        if self.variance is not None:
            V = np.asarray(self.variance)
            if V.shape != shape or not (np.issubdtype(V.dtype, np.floating) or np.issubdtype(V.dtype, np.integer)):
                raise ValueError(f"variance must be a real array with the grid shape {shape}")
            V = V.astype(float)
            vu = V[unmasked]
            if not np.isfinite(vu).all() or np.any(vu < 0):
                raise ValueError("variance must be finite and >= 0 on unmasked voxels")
        if self.bragg == "bragg_masked" and (M is None or not M.any()):
            raise ValueError("bragg='bragg_masked' needs a mask that excludes the Bragg regions")
        if M is None and self.mask_reasons:
            raise ValueError("mask_reasons given without a mask")
        for name in ("mask_reasons", "corrections"):
            items = tuple(getattr(self, name))
            if not all(isinstance(s, str) and s for s in items):
                raise ValueError(f"{name} must be non-empty strings")
            object.__setattr__(self, name, items)
        if not isinstance(self.provenance, dict):
            raise ValueError("provenance must be a dict")
        prov = json.loads(_json_safe("provenance", self.provenance))
        object.__setattr__(self, "provenance", prov)
        object.__setattr__(self, "intensity", _readonly(I))
        object.__setattr__(self, "variance", None if V is None else _readonly(V))
        object.__setattr__(self, "mask", None if M is None else _readonly(M))

    # -- constructors -----------------------------------------------------------
    @classmethod
    def from_discus_volume(cls, volume, *, units: str = "DISCUS intensity (whole crystal)") -> "ScatteringVolume":
        """Wrap a verified DISCUS 6.23.01 ``DiscusVolume``; keeps its file hash and asserted metadata."""
        prov = {"source": "DISCUS common/NeXus", "reader": "discuspdf.read_discus_volume",
                **{k: v for k, v in volume.provenance.items()},
                "discus_metadata": {k: volume.meta[k] for k in sorted(volume.meta)}}
        return cls(ReciprocalGrid.from_discus_volume(volume), volume.values, units, bragg=volume.with_bragg,
                   provenance=prov)

    @classmethod
    def from_discus_file(cls, path, *, config=None, units: str | None = None) -> "ScatteringVolume":
        """Read a DISCUS 6.23.01 ``form nexus`` reciprocal-space file, including **planes and lines**.

        Unlike :func:`~discuspdf.discus_output.read_discus_volume` (volumes only), the grid is built with
        :class:`ReciprocalGrid`, which ignores the step vector of a one-point axis (DISCUS writes zeros
        there, or ``tl − ll`` if given). The file must be hkl / reciprocal / basecell-fractional and its
        content ``intensity``, ``real part`` or ``imaginary part`` (all observed in 6.23.01).

        With ``config`` (the :class:`~discuspdf.fourier_macro.FourierConfig` that wrote the file), the
        file's grid, Bragg label and content are checked against it and the configuration (seed,
        ``set aver``, lots, labels) is stored under ``provenance['fourier_config']``."""
        import h5py
        from .discus_output import _GROUP, _REQUIRED, _text
        from .lattice import Lattice
        p = Path(path)
        if not p.is_file() or not h5py.is_hdf5(p):
            raise ValueError(f"not an HDF5 file: {p}")
        with h5py.File(p, "r") as f:
            if _GROUP not in f:
                raise ValueError(f"{p.name}: no {_GROUP} group; not DISCUS common/nexus output")
            g = f[_GROUP]
            missing = [k for k in _REQUIRED if k not in g]
            if missing:
                raise ValueError(f"{p.name}: missing datasets {missing}")
            raw = {k: g[k][()] for k in g}
        meta = {k: _text(v) for k, v in raw.items() if np.ndim(v) == 0}
        want = {"audit_conform_dict_name": "Disorder unified data", "coordinate_unit": "basecell_fractional",
                "data_type_axes": "hkl", "data_type_reciprocal": "reciprocal"}
        for key, value in want.items():
            if meta.get(key) != value:
                raise ValueError(f"{p.name}: {key}={meta.get(key)!r}, expected {value!r}")
        content = meta["data_type_content"]
        if content not in _DISCUS_CONTENT.values():
            raise ValueError(f"{p.name}: content {content!r} is not a reciprocal-space value this reader knows")
        bragg = meta["data_type_with_bragg"]
        if bragg not in ("bragg_present", "bragg_subtracted"):
            raise ValueError(f"{p.name}: unknown data_type_with_bragg={bragg!r}")
        values = np.asarray(raw["data_values"], dtype=float)
        if values.ndim != 3 or tuple(np.asarray(raw["data_dimension"]).tolist()) != values.shape:
            raise ValueError(f"{p.name}: data_dimension {raw['data_dimension']} != array shape {values.shape}")
        if tuple(np.asarray(raw["data_axes"]).tolist()) != (1, 2, 3):
            raise ValueError(f"{p.name}: data_axes={raw['data_axes']}; only (1,2,3) was verified")
        lattice = Lattice.from_parameters(*np.asarray(raw["unit_cell_lengths"], float),
                                          *np.asarray(raw["unit_cell_angles"], float))
        grid = ReciprocalGrid(lattice, np.asarray(raw["data_corner"], float),
                              np.asarray(raw["data_increment_vector"], float).T, values.shape)
        prov: dict[str, Any] = {"source": "DISCUS common/NeXus", "reader": "ScatteringVolume.from_discus_file",
                                "path": str(p), "sha256": hashlib.sha256(p.read_bytes()).hexdigest(),
                                "content": content, "discus_metadata": {k: meta[k] for k in sorted(meta)}}
        if config is not None:
            prov["fourier_config"] = _check_against_config(config, p.name, grid, bragg, content,
                                                           radiation=_file_radiation(meta, raw))
        return cls(grid, values, units or f"DISCUS {content} (whole crystal)", bragg=bragg, provenance=prov)

    # -- derived volumes --------------------------------------------------------
    def with_mask(self, mask, reason: str) -> "ScatteringVolume":
        """New volume whose mask is the union of the current mask and ``mask``; ``reason`` is recorded."""
        m = np.asarray(mask)
        if m.dtype != bool or m.shape != self.grid.shape:
            raise ValueError("mask must be a boolean array with the grid shape")
        if not isinstance(reason, str) or not reason:
            raise ValueError("a non-empty reason is required")
        new = m if self.mask is None else (self.mask | m)
        return replace(self, mask=new, mask_reasons=(*self.mask_reasons, reason))

    def mask_bragg(self, radius: float) -> "ScatteringVolume":
        """Mask voxels within ``radius`` Å⁻¹ of integer hkl and label the result ``bragg_masked``.

        Refused for data that are already ``bragg_subtracted`` (masking would hide the subtraction's
        residuals without saying so); mask those explicitly with :meth:`with_mask` instead."""
        if self.bragg == "bragg_subtracted":
            raise ValueError("data are already bragg_subtracted; use with_mask() with an explicit reason")
        m = mask_near_bragg(self.grid, radius)
        out = self.with_mask(m, f"bragg: |Q - Q_hkl| <= {radius:g} 1/angstrom")
        return replace(out, bragg="bragg_masked")

    def with_correction(self, intensity, description: str, *, variance=None) -> "ScatteringVolume":
        """New volume with corrected values (e.g. background subtracted) and the step recorded."""
        if not isinstance(description, str) or not description:
            raise ValueError("a non-empty description is required")
        return replace(self, intensity=intensity, variance=self.variance if variance is None else variance,
                       corrections=(*self.corrections, description))

    def filled(self, fill_value: float) -> np.ndarray:
        """Intensity with masked voxels replaced by ``fill_value`` (a copy)."""
        out = np.array(self.intensity, copy=True)
        if self.mask is not None:
            out[self.mask] = float(fill_value)
        return out

    # -- transforms -------------------------------------------------------------
    def _transform_input(self, fill_masked, require_diffuse: bool):
        if require_diffuse and self.bragg not in ("bragg_subtracted", "bragg_masked"):
            raise ValueError(f"bragg={self.bragg!r}: a 3D-ΔPDF needs bragg_subtracted or bragg_masked data "
                             "(pass require_diffuse=False for a Patterson function)")
        if self.mask is not None and self.mask.any():
            if fill_masked is None:
                raise ValueError(f"{int(self.mask.sum())} voxels are masked; pass fill_masked=<value> to state "
                                 "what the transform should use there")
            if not np.isfinite(float(fill_masked)):
                raise ValueError("fill_masked must be finite")
            return self.filled(fill_masked), f"masked voxels filled with {float(fill_masked):g}"
        return np.array(self.intensity), None

    def to_reciprocal_volume(self, *, fill_masked: float | None = None):
        """Cartesian-Q :class:`ReciprocalVolume`, only when the grid's Cartesian transform is exact."""
        from .reciprocal import ReciprocalVolume
        comp = self.grid.fft_compatibility()
        if not comp.cartesian:
            raise ValueError("no exact Cartesian Q axes for this grid: " + "; ".join(comp.reasons))
        values, fill = self._transform_input(fill_masked, require_diffuse=False)
        pre = {"bragg_subtracted": ("bragg_subtracted",), "bragg_masked": ("bragg_excluded",)}.get(self.bragg, ())
        meta = {"units": self.units, "corrections": list(self.corrections) + ([fill] if fill else []),
                "provenance": self.provenance}
        return ReciprocalVolume(values, self.grid.q_axes(), "cartesian_q", "1/angstrom", preprocessing=pre,
                                metadata=meta)

    def delta_pdf(self, *, fill_masked: float | None = None, require_diffuse: bool = True, **kwargs):
        """3D-ΔPDF via :func:`delta_pdf_fractional` (any cell metric; hkl-aligned, centred volume).

        The result's ``direct_basis`` is in the canonical orientation (a ∥ x); fractional ``u``
        coordinates do not depend on orientation."""
        from .discus_output import DiscusVolume, delta_pdf_fractional
        comp = self.grid.fft_compatibility()
        if not comp.fractional:
            raise ValueError("delta_pdf_fractional is not exact for this grid: " + "; ".join(comp.reasons))
        values, fill = self._transform_input(fill_masked, require_diffuse)
        prov = {**self.provenance, "units": self.units, "bragg": self.bragg,
                "corrections": list(self.corrections) + ([fill] if fill else [])}
        a, b, c, al, be, ga = self.grid.lattice.parameters
        vol = DiscusVolume(values, self.grid.corner, self.grid.axis_vectors.T, np.array([a, b, c]),
                           np.array([al, be, ga]),
                           "bragg_subtracted" if self.bragg != "bragg_present" else "bragg_present",
                           {}, prov)
        return delta_pdf_fractional(vol, require_diffuse=False, **kwargs)

    # -- xarray -----------------------------------------------------------------
    def to_xarray(self):
        """``xarray.Dataset`` with dims ``(i0, i1, i2)``, integer index coords, and h/k/l (hkl-aligned grids)
        and qx/qy/qz (Cartesian-aligned grids) as extra coords. Metadata are JSON strings in ``attrs``."""
        import xarray as xr
        coords: dict[str, Any] = {d: np.arange(n) for d, n in zip(_DIMS, self.grid.shape)}
        if self.grid.is_hkl_axis_aligned():
            for d, name, axis in zip(_DIMS, "hkl", self.grid.hkl_axes()):
                coords[name] = (d, axis, {"units": "r.l.u."})
        if self.grid.is_q_axis_aligned():
            for d, name, axis in zip(_DIMS, ("qx", "qy", "qz"), self.grid.q_axes()):
                coords[name] = (d, axis, {"units": "1/angstrom", "convention": "Q = 2*pi*hkl@B"})
        data = {"intensity": (_DIMS, np.array(self.intensity), {"units": self.units})}
        if self.variance is not None:
            data["variance"] = (_DIMS, np.array(self.variance), {"units": f"({self.units})^2"})
        if self.mask is not None:
            data["mask"] = (_DIMS, np.array(self.mask), {"meaning": "True = excluded"})
        attrs = {"discuspdf_schema": SCHEMA_NAME, "discuspdf_schema_version": SCHEMA_VERSION,
                 "grid": self.grid.to_json(), "intensity_units": self.units, "bragg": self.bragg,
                 "mask_reasons": _json_safe("mask_reasons", list(self.mask_reasons)),
                 "corrections": _json_safe("corrections", list(self.corrections)),
                 "provenance": _json_safe("provenance", self.provenance)}
        return xr.Dataset(data, coords=coords, attrs=attrs)

    @classmethod
    def from_xarray(cls, ds) -> "ScatteringVolume":
        """Inverse of :meth:`to_xarray`. Slicing with ``isel`` along i0/i1/i2 (any constant integer step)
        is supported: the grid is rebuilt from the index coordinates. Stale h/k/l or Q coordinates,
        reordered dims, dropped attributes and unknown schemas are refused."""
        attrs = dict(ds.attrs)
        if attrs.get("discuspdf_schema") != SCHEMA_NAME:
            raise ValueError("dataset has no discuspdf schema attribute; arithmetic in xarray drops attrs "
                             "unless xr.set_options(keep_attrs=True) is used")
        if attrs.get("discuspdf_schema_version") != SCHEMA_VERSION:
            raise ValueError(f"unsupported schema version {attrs.get('discuspdf_schema_version')!r}")
        if "intensity" not in ds:
            raise ValueError("dataset has no 'intensity' variable")
        for name in ("intensity", "variance", "mask"):
            if name in ds and tuple(ds[name].dims) != _DIMS:
                raise ValueError(f"{name} dims {tuple(ds[name].dims)} != {_DIMS}; transpose back first")
        parent = ReciprocalGrid.from_json(attrs["grid"])
        start, steps, shape = [], [], []
        for d in _DIMS:
            idx = np.asarray(ds[d].values)
            if idx.ndim != 1 or idx.size == 0 or not np.issubdtype(idx.dtype, np.integer):
                raise ValueError(f"index coordinate {d} must be a non-empty 1D integer array")
            step = 1
            if idx.size > 1:
                diffs = np.diff(idx)
                if np.any(diffs != diffs[0]) or diffs[0] == 0:
                    raise ValueError(f"index coordinate {d} is not an evenly spaced selection")
                step = int(diffs[0])
            start.append(int(idx[0]))
            steps.append(step)
            shape.append(int(idx.size))
        V = parent.axis_vectors * np.asarray(steps, float)[:, None]
        corner = parent.corner + np.asarray(start, float) @ np.where(
            np.asarray(parent.shape)[:, None] > 1, parent.axis_vectors, 0.0)
        grid = ReciprocalGrid(parent.lattice, corner, V, tuple(shape))
        if grid.is_hkl_axis_aligned():
            for d, name, axis in zip(_DIMS, "hkl", grid.hkl_axes()):
                if name in ds.coords and not np.allclose(ds[name].values, axis, rtol=0, atol=1e-9):
                    raise ValueError(f"coordinate {name} does not match the grid; was it edited?")
        S, q0 = grid.q_step_vectors(), grid.q_corner()
        for k, (d, name) in enumerate(zip(_DIMS, ("qx", "qy", "qz"))):
            if name not in ds.coords:
                continue
            # A reversed selection is still parallel to Cartesian axis k, only with a negative step.
            expected = q0[k] + S[k, k] * np.arange(grid.shape[k])
            off = np.delete(S[k], k)
            if (np.any(np.abs(off) > 1e-9 * max(np.linalg.norm(S[k]), 1e-300))
                    or not np.allclose(ds[name].values, expected, rtol=1e-12, atol=1e-12)):
                raise ValueError(f"coordinate {name} does not match the grid; was it edited?")
        return cls(grid, np.asarray(ds["intensity"].values), attrs["intensity_units"], bragg=attrs["bragg"],
                   variance=np.asarray(ds["variance"].values) if "variance" in ds else None,
                   mask=np.asarray(ds["mask"].values) if "mask" in ds else None,
                   mask_reasons=tuple(json.loads(attrs.get("mask_reasons", "[]"))),
                   corrections=tuple(json.loads(attrs.get("corrections", "[]"))),
                   provenance=json.loads(attrs.get("provenance", "{}")))

    # -- HDF5 persistence -------------------------------------------------------
    def save_hdf5(self, path: str | Path, *, overwrite: bool = False, compression: str | None = None) -> Path:
        """Write schema ``discuspdf.ScatteringVolume`` v1. Refuses to overwrite unless asked."""
        import h5py
        p = Path(path)
        if p.exists() and not overwrite:
            raise FileExistsError(f"{p} exists; pass overwrite=True to replace it")
        if compression not in (None, "gzip"):
            raise ValueError("compression must be None or 'gzip'")
        with h5py.File(p, "w") as f:
            f.attrs["schema"] = SCHEMA_NAME
            f.attrs["schema_version"] = SCHEMA_VERSION
            f.attrs["written_by"] = f"discuspdf {_package_version()}"
            f.attrs["written_utc"] = _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds")
            f.attrs["grid"] = self.grid.to_json()
            f.attrs["intensity_units"] = self.units
            f.attrs["bragg"] = self.bragg
            f.attrs["mask_reasons"] = _json_safe("mask_reasons", list(self.mask_reasons))
            f.attrs["corrections"] = _json_safe("corrections", list(self.corrections))
            f.attrs["provenance"] = _json_safe("provenance", self.provenance)
            f.attrs["intensity_sha256"] = _digest(self.intensity)
            if self.variance is not None:
                f.attrs["variance_sha256"] = _digest(self.variance)
            if self.mask is not None:
                f.attrs["mask_sha256"] = _digest(self.mask.astype(np.uint8))
            f.create_dataset("intensity", data=self.intensity, dtype="f8", compression=compression)
            if self.variance is not None:
                f.create_dataset("variance", data=self.variance, dtype="f8", compression=compression)
            if self.mask is not None:
                f.create_dataset("mask", data=self.mask.astype(np.uint8), dtype="u1", compression=compression)
        return p

    @classmethod
    def load_hdf5(cls, path: str | Path) -> "ScatteringVolume":
        """Read a file written by :meth:`save_hdf5`; refuses other schemas/versions and corrupted data."""
        import h5py
        p = Path(path)
        if not p.is_file() or not h5py.is_hdf5(p):
            raise ValueError(f"not an HDF5 file: {p}")
        with h5py.File(p, "r") as f:
            a = {k: (v.decode() if isinstance(v, bytes) else v) for k, v in f.attrs.items()}
            if a.get("schema") != SCHEMA_NAME:
                raise ValueError(f"{p.name}: schema {a.get('schema')!r} is not {SCHEMA_NAME!r}")
            version = int(a.get("schema_version", -1))
            if version not in _READERS:
                raise ValueError(f"{p.name}: schema version {version} is not supported "
                                 f"(this code reads {sorted(_READERS)})")
            return _READERS[version](cls, f, a, p)


# DISCUS 6.23.01 'output value' -> data_type_content written in the file (observed in the B probes).
_DISCUS_CONTENT = {"int": "intensity", "real": "real part", "imag": "imaginary part"}


# DISCUS 6.23.01 data_radiation values (observed in the D reference probe).
_DISCUS_RADIATION = {"xray": "xray", "neutron": "neutron"}


def _file_radiation(meta: dict, raw: dict) -> dict:
    """Radiation metadata of a DISCUS file: data_radiation, data_rad_symbol, data_rad_length[0] (Å)."""
    out = {"radiation": meta.get("data_radiation"), "symbol": meta.get("data_rad_symbol"), "length": None}
    if "data_rad_length" in raw:
        length = np.ravel(np.asarray(raw["data_rad_length"], float))
        if length.size and np.isfinite(length[0]) and length[0] > 0:
            out["length"] = float(length[0])
    return out


def _check_against_config(config, name: str, grid: ReciprocalGrid, bragg: str, content: str,
                          radiation: dict | None = None) -> dict:
    """Refuse a file that the given FourierConfig cannot have written; return the config as a dict.

    ``radiation`` (from :func:`_file_radiation`) is compared where the file records it: the radiation
    always, the symbol for a symbolic wavelength, the length (rel. 1e-4) for a numeric one."""
    from .fourier_macro import FourierConfig
    if not isinstance(config, FourierConfig):
        raise TypeError("config must be a FourierConfig")
    values = {v for v, n in config.outputs if n == name}
    if not values:
        raise ValueError(f"{name} is not one of this configuration's outputs {[n for _, n in config.outputs]}")
    value = values.pop()
    if _DISCUS_CONTENT.get(value) != content:
        raise ValueError(f"{name}: file content {content!r} but the configuration wrote 'value {value}'")
    if bragg != config.bragg:
        raise ValueError(f"{name}: file says {bragg!r}, configuration (set aver,{config.aver_percent:g}) "
                         f"implies {config.bragg!r}")
    if radiation is not None:
        rad = radiation.get("radiation")
        if rad is not None and rad.strip().lower() != _DISCUS_RADIATION[config.radiation]:
            raise ValueError(f"{name}: file radiation {rad!r} but the configuration used {config.radiation!r}")
        sym, length = radiation.get("symbol"), radiation.get("length")
        if isinstance(config.wavelength, str):
            if sym and sym.strip() and sym.strip().lower() != config.wavelength:
                raise ValueError(f"{name}: file wavelength symbol {sym!r} but the configuration used "
                                 f"{config.wavelength!r}")
        elif length is not None and abs(length - config.wavelength) > 1e-4 * config.wavelength:
            raise ValueError(f"{name}: file wavelength {length:g} Å but the configuration used "
                             f"{config.wavelength:g} Å")
    cg = config.grid
    if cg.shape != grid.shape:
        raise ValueError(f"{name}: grid shape {grid.shape} != configured {cg.shape}")
    # Corners are written to the macro with 6 significant digits.
    scale = max(1.0, float(np.max(np.abs(cg.corner))), float(np.max(np.abs(cg.axis_vectors))))
    V_file, V_cfg = grid.axis_vectors, cg.axis_vectors
    for k in grid.active_axes:
        if not np.allclose(V_file[k], V_cfg[k], rtol=0, atol=1e-6 * scale):
            raise ValueError(f"{name}: step of axis {k} {V_file[k].tolist()} != configured {V_cfg[k].tolist()}")
    if not np.allclose(grid.corner, cg.corner, rtol=0, atol=1e-6 * scale):
        raise ValueError(f"{name}: corner {grid.corner.tolist()} != configured {cg.corner.tolist()}")
    return config.to_dict()


def _read_v1(cls, f, a, p: Path) -> ScatteringVolume:
    required = ("grid", "intensity_units", "bragg", "mask_reasons", "corrections", "provenance", "intensity_sha256")
    missing = [k for k in required if k not in a]
    if missing:
        raise ValueError(f"{p.name}: missing attributes {missing}")
    for name in ("intensity", "variance", "mask"):
        if name in f and f[name].dtype.kind not in "fiu":
            raise ValueError(f"{p.name}: dataset {name} has non-numeric dtype {f[name].dtype}")
    if "intensity" not in f:
        raise ValueError(f"{p.name}: no intensity dataset")
    I = np.asarray(f["intensity"][()], dtype=float)
    if _digest(I) != a["intensity_sha256"]:
        raise ValueError(f"{p.name}: intensity does not match its stored SHA-256 (corrupted or edited)")
    var = None
    if "variance" in f:
        var = np.asarray(f["variance"][()], dtype=float)
        if "variance_sha256" not in a or _digest(var) != a["variance_sha256"]:
            raise ValueError(f"{p.name}: variance does not match its stored SHA-256 (corrupted or edited)")
    mask = None
    if "mask" in f:
        raw = np.asarray(f["mask"][()])
        if not np.isin(raw, (0, 1)).all():
            raise ValueError(f"{p.name}: mask must contain only 0 and 1")
        if "mask_sha256" not in a or _digest(raw.astype(np.uint8)) != a["mask_sha256"]:
            raise ValueError(f"{p.name}: mask does not match its stored SHA-256 (corrupted or edited)")
        mask = raw.astype(bool)
    return cls(ReciprocalGrid.from_json(a["grid"]), I, a["intensity_units"], bragg=a["bragg"],
               variance=var,
               mask=mask, mask_reasons=tuple(json.loads(a["mask_reasons"])),
               corrections=tuple(json.loads(a["corrections"])), provenance=json.loads(a["provenance"]))


# Schema version -> reader. A future version adds its reader (or a migration to the current class) here.
_READERS = {1: _read_v1}
