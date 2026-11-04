"""Physically labelled slices, isosurfaces and VTK export of reciprocal-space and 3D-ΔPDF volumes
(Package B, Feature 20).

Everything goes through :class:`SampledVolume`: values on a regular grid whose point ``(i0, i1, i2)``
sits at the Cartesian position ``origin + i @ steps`` (rows of ``steps`` are the step vectors). Its
constructors take the calibrated objects of this package, so the geometry comes from source
metadata, never from array indices:

* :meth:`SampledVolume.from_scattering_volume` — Q in Å⁻¹ (2π convention) from the ``ReciprocalGrid``;
* :meth:`SampledVolume.from_fractional_fft` — real space in Å from ``delta_pdf_fractional``
  (fractional u axes times the direct basis, so non-orthogonal cells keep their angles);
* :meth:`SampledVolume.from_fft_result` — real space in Å from the Cartesian transform.

Slices are drawn with their true in-plane geometry (``aspect='equal'``, non-orthogonal axes appear
sheared). Signed data (Bragg-subtracted intensity, ΔPDF) use a symmetric diverging colour map;
masked voxels are drawn grey. VTK output: ``.vti`` (ImageData) only when every axis runs along +x,
+y, +z; otherwise ``.vts`` (StructuredGrid) with explicit point coordinates. Both are VTK XML with raw
appended binary data, readable by ParaView and the ``vtk`` package; no VTK installation is needed
to write them.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any
import json
import re

import numpy as np

_TOL = 1e-9


def _readonly(a, dtype=float):
    a = np.array(a, dtype=dtype, copy=True)
    a.flags.writeable = False
    return a


@dataclass(frozen=True, eq=False)
class SampledVolume:
    """Values on a regular grid with physical (Cartesian) coordinates.

    ``index_vectors`` (optional) gives each axis's step in crystal coordinates (hkl rows for
    reciprocal space, fractional uvw rows for real space); it is used only for axis labels.
    """
    values: np.ndarray
    origin: np.ndarray
    steps: np.ndarray
    space: str
    quantity: str
    signed: bool = False
    mask: np.ndarray | None = None
    index_vectors: np.ndarray | None = None
    extra: dict[str, np.ndarray] = field(default_factory=dict)
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        v = np.asarray(self.values, dtype=float)
        if v.ndim != 3 or v.size == 0:
            raise ValueError("values must be a non-empty 3D array (use length-1 axes for planes)")
        if self.space not in ("reciprocal", "real"):
            raise ValueError("space must be 'reciprocal' or 'real'")
        if not isinstance(self.quantity, str) or not re.fullmatch(r"[A-Za-z0-9_]+", self.quantity):
            raise ValueError("quantity must be a simple identifier (used as the VTK array name)")
        o = np.asarray(self.origin, float)
        S = np.asarray(self.steps, float)
        if o.shape != (3,) or S.shape != (3, 3) or not (np.isfinite(o).all() and np.isfinite(S).all()):
            raise ValueError("origin must be (3,) and steps (3, 3), finite")
        active = [k for k in range(3) if v.shape[k] > 1]
        if active and np.linalg.matrix_rank(S[active]) != len(active):
            raise ValueError("step vectors of the active axes are linearly dependent")
        m = None
        if self.mask is not None:
            m = np.asarray(self.mask)
            if m.dtype != bool or m.shape != v.shape:
                raise ValueError("mask must be a boolean array with the shape of values")
        if not np.isfinite(v[~m] if m is not None else v).all():
            raise ValueError("values must be finite where not masked")
        iv = None
        if self.index_vectors is not None:
            iv = np.asarray(self.index_vectors, float)
            if iv.shape != (3, 3):
                raise ValueError("index_vectors must be 3x3")
        extra = {}
        for name, arr in dict(self.extra).items():
            a = np.asarray(arr, dtype=float)
            if a.shape != v.shape or not re.fullmatch(r"[A-Za-z0-9_]+", name):
                raise ValueError(f"extra array {name!r} must have the values' shape and a simple name")
            extra[name] = _readonly(a)
        object.__setattr__(self, "values", _readonly(v))
        object.__setattr__(self, "origin", _readonly(o))
        object.__setattr__(self, "steps", _readonly(S))
        object.__setattr__(self, "mask", None if m is None else _readonly(m, bool))
        object.__setattr__(self, "index_vectors", None if iv is None else _readonly(iv))
        object.__setattr__(self, "extra", extra)
        object.__setattr__(self, "metadata", json.loads(json.dumps(self.metadata, default=str)))

    # -- constructors -------------------------------------------------------------
    @classmethod
    def from_scattering_volume(cls, sv, *, signed: bool | None = None) -> "SampledVolume":
        """Q (Å⁻¹, 2π convention) of a :class:`~discuspdf.scattering_data.ScatteringVolume`.

        ``signed`` (diverging colour map, no log scale) follows the data by default: True when any
        unmasked value is negative. A ``bragg_subtracted`` volume from DISCUS ``set aver`` is
        |F − ⟨F⟩|² ≥ 0, so it is unsigned and can be drawn on a log scale; an intensity difference
        computed elsewhere has negative values and is signed. Pass ``signed=`` to choose explicitly."""
        g = sv.grid
        if signed is None:
            values = np.where(sv.mask, np.nan, sv.intensity) if sv.mask is not None else sv.intensity
            signed = bool(np.isfinite(values).any() and np.nanmin(values) < 0)
        extra = {} if sv.variance is None else {"variance": np.nan_to_num(sv.variance, nan=0.0)}
        return cls(sv.intensity, g.q_corner(), g.q_step_vectors(), "reciprocal", "intensity",
                   signed=bool(signed), mask=sv.mask, index_vectors=g.axis_vectors, extra=extra,
                   metadata={"units": sv.units, "bragg": sv.bragg, "grid": g.to_dict(),
                             "provenance": sv.provenance, "coordinates": "Q = 2*pi*hkl@B, 1/angstrom"})

    @classmethod
    def from_fractional_fft(cls, result, *, quantity: str = "delta_pdf", part: str = "real") -> "SampledVolume":
        """Real space (Å) of a :class:`~discuspdf.discus_output.FractionalFFTResult`."""
        if part not in ("real", "imag"):
            raise ValueError("part must be 'real' or 'imag'")
        A = np.asarray(result.direct_basis, float)
        u0 = np.array([ax[0] for ax in result.u_axes], float)
        du = np.array([ax[1] - ax[0] if len(ax) > 1 else 0.0 for ax in result.u_axes], float)
        values = result.complex_volume.real if part == "real" else result.complex_volume.imag
        return cls(values, u0 @ A, du[:, None] * A, "real", quantity, signed=True,
                   index_vectors=np.diag(du),
                   metadata={"normalization": result.normalization, "direct_basis": A.tolist(),
                             "part": part, "provenance": result.provenance, "coordinates": "r = u@A, angstrom"})

    @classmethod
    def from_fft_result(cls, result, *, quantity: str = "delta_pdf", part: str = "real") -> "SampledVolume":
        """Real space (Å) of an :class:`~discuspdf.volume.FFTResult` (Cartesian transform)."""
        if part not in ("real", "imag"):
            raise ValueError("part must be 'real' or 'imag'")
        r = [np.asarray(ax, float) for ax in result.r_axes]
        dr = np.array([ax[1] - ax[0] if len(ax) > 1 else 0.0 for ax in r])
        values = result.complex_volume.real if part == "real" else result.complex_volume.imag
        return cls(values, np.array([ax[0] for ax in r]), np.diag(dr), "real", quantity, signed=True,
                   metadata={"normalization": result.diagnostics.normalization, "part": part,
                             "coordinates": "Cartesian r, angstrom"})

    # -- geometry -------------------------------------------------------------------
    @property
    def shape(self) -> tuple[int, int, int]:
        return tuple(self.values.shape)  # type: ignore[return-value]

    @property
    def unit(self) -> str:
        return "Å⁻¹" if self.space == "reciprocal" else "Å"

    @property
    def active_axes(self) -> tuple[int, ...]:
        return tuple(k for k in range(3) if self.shape[k] > 1)

    def points(self) -> np.ndarray:
        """Cartesian position of every point, shape ``shape + (3,)``."""
        idx = np.stack(np.meshgrid(*[np.arange(n, dtype=float) for n in self.shape], indexing="ij"), axis=-1)
        return self.origin + idx @ self.steps

    def is_axis_aligned(self, tol: float = _TOL) -> bool:
        """Every active axis k steps along +Cartesian axis k (the condition for VTK ImageData)."""
        for k in self.active_axes:
            s = self.steps[k]
            if s[k] <= 0 or np.any(np.abs(np.delete(s, k)) > tol * np.linalg.norm(s)):
                return False
        return True

    def filled(self, fill_value: float = np.nan) -> np.ndarray:
        out = np.array(self.values, copy=True)
        if self.mask is not None:
            out[self.mask] = fill_value
        return out

    def _axis_label(self, k: int) -> str:
        sym = "Q" if self.space == "reciprocal" else "r"
        if self.index_vectors is None:
            return f"{sym} along axis {k} ({self.unit})"
        v = self.index_vectors[k]
        v = v / np.max(np.abs(v))
        txt = " ".join(f"{x:.3g}" for x in v)
        bracket = f"[{txt}]*" if self.space == "reciprocal" else f"[{txt}]"
        return f"{sym} ∥ {bracket} ({self.unit})"


# ----------------------------------------------------------------------------- slices
def _plane_geometry(vol: SampledVolume, axis: int, index: int):
    p, q = [k for k in range(3) if k != axis]
    sp, sq = vol.steps[p], vol.steps[q]
    if np.linalg.norm(sp) == 0 or np.linalg.norm(sq) == 0:
        raise ValueError("the slice plane has an axis with a single point; choose another axis")
    e1 = sp / np.linalg.norm(sp)
    t = sq - (sq @ e1) * e1
    e2 = t / np.linalg.norm(t)
    base = vol.origin + index * vol.steps[axis]
    i = np.arange(vol.shape[p] + 1) - 0.5
    j = np.arange(vol.shape[q] + 1) - 0.5
    I, J = np.meshgrid(i, j, indexing="ij")
    P = base + I[..., None] * sp + J[..., None] * sq
    normal = np.cross(e1, e2)
    return p, q, P @ e1, P @ e2, float(base @ normal)


def plot_slice(vol: SampledVolume, axis: int = 2, index: int | None = None, *, ax=None, cmap=None,
               log: bool = False, vmax: float | None = None, title: str | None = None):
    """Draw one layer of ``vol`` perpendicular to array ``axis`` with its true in-plane geometry.

    Coordinates are the projections of the voxel corners on an orthonormal basis of the plane
    (first in-plane axis along x), in Å⁻¹ or Å. Signed data get a symmetric diverging map
    (``RdBu_r``); unsigned data ``viridis`` (``log=True`` for a log scale). Masked voxels are grey.
    Returns the matplotlib Axes."""
    import matplotlib.pyplot as plt
    from matplotlib import colors

    if axis not in (0, 1, 2):
        raise ValueError("axis must be 0, 1 or 2")
    n = vol.shape[axis]
    index = n // 2 if index is None else index
    if isinstance(index, bool) or not isinstance(index, (int, np.integer)) or not 0 <= index < n:
        raise ValueError(f"index must be an integer in [0, {n})")
    p, q, X, Y, offset = _plane_geometry(vol, axis, int(index))
    data = np.take(vol.filled(np.nan), index, axis=axis)
    data = np.ma.masked_invalid(data)
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 5))
    if vol.signed:
        if log:
            raise ValueError("log scale is not meaningful for signed data")
        top = vmax if vmax is not None else (float(np.max(np.abs(data))) if data.count() else 1.0) or 1.0
        cm = plt.get_cmap(cmap or "RdBu_r").copy()
        norm = colors.Normalize(-top, top)
    else:
        cm = plt.get_cmap(cmap or "viridis").copy()
        if log:
            pos = data[data > 0]
            if pos.size == 0:
                raise ValueError("log scale needs positive values")
            norm = colors.LogNorm(float(pos.min()), vmax or float(pos.max()))
            data = np.ma.masked_less_equal(data, 0)
        else:
            norm = colors.Normalize(vmax=vmax)
    cm.set_bad("0.6")
    mesh = ax.pcolormesh(X, Y, data, cmap=cm, norm=norm, shading="flat")
    ax.set_aspect("equal")
    ax.set_xlabel(vol._axis_label(p))
    S = vol.steps
    angle = np.degrees(np.arccos(np.clip(S[p] @ S[q] / np.linalg.norm(S[p]) / np.linalg.norm(S[q]), -1, 1)))
    ax.set_ylabel(f"⊥ in plane (axis {q} at {angle:.4g}° to axis {p}) ({vol.unit})")
    ax.set_title(title or f"{vol.quantity}: layer {index} of axis {axis}, plane offset {offset:.4g} {vol.unit}")
    ax.figure.colorbar(mesh, ax=ax, label=f"{vol.quantity} ({vol.metadata.get('units', 'arb.')})")
    return ax


# ----------------------------------------------------------------------------- isosurfaces
def isosurface(vol: SampledVolume, level: float, *, fill_value: float | None = None):
    """Marching-cubes isosurface at ``level`` in Cartesian coordinates (needs scikit-image).

    Returns ``(vertices (m, 3), faces (k, 3))``. Masked voxels must be given a ``fill_value``."""
    try:
        from skimage.measure import marching_cubes
    except ImportError as exc:  # pragma: no cover - depends on environment
        raise ImportError("isosurface needs scikit-image: pip install 'discuspdf-starter[viz]'") from exc
    if len(vol.active_axes) != 3:
        raise ValueError("an isosurface needs a 3D volume")
    if vol.mask is not None and vol.mask.any():
        if fill_value is None:
            raise ValueError("masked voxels present; pass fill_value=")
        data = vol.filled(fill_value)
    else:
        data = np.array(vol.values)
    if not data.min() < level < data.max():
        raise ValueError(f"level {level} is outside the data range [{data.min():.4g}, {data.max():.4g}]")
    verts, faces, _, _ = marching_cubes(data, level=level)
    return vol.origin + verts @ vol.steps, faces


def plot_isosurfaces(vol: SampledVolume, levels, *, ax=None, colors_=None, alpha: float = 0.6,
                     fill_value: float | None = None):
    """3D plot of one or more isosurfaces (for signed data, e.g. ``levels=(+c, -c)``)."""
    import matplotlib.pyplot as plt
    if ax is None:
        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_subplot(projection="3d")
    levels = list(levels)
    palette = colors_ or ["tab:red" if lv > 0 else "tab:blue" for lv in levels]
    for lv, c in zip(levels, palette):
        v, f = isosurface(vol, lv, fill_value=fill_value)
        ax.plot_trisurf(v[:, 0], v[:, 1], f, v[:, 2], color=c, alpha=alpha, linewidth=0)
    corners = vol.origin + np.array(np.meshgrid([0, 1], [0, 1], [0, 1], indexing="ij")).reshape(3, -1).T \
        * (np.array(vol.shape) - 1) @ vol.steps
    span = corners.max(0) - corners.min(0)
    ax.set_box_aspect(np.where(span > 0, span, 1))
    for setter, k in ((ax.set_xlabel, "x"), (ax.set_ylabel, "y"), (ax.set_zlabel, "z")):
        setter(f"{k} ({vol.unit})")
    ax.set_title(f"{vol.quantity} isosurfaces at {', '.join(f'{lv:.3g}' for lv in levels)}")
    return ax


# ----------------------------------------------------------------------------- VTK XML
def _appended(arrays: list[tuple[str, np.ndarray, int]]):
    """Header fragments and raw appended block (UInt64 byte-count headers, little endian)."""
    blocks, offset, refs = [], 0, []
    for name, arr, ncomp in arrays:
        if np.asarray(arr).dtype == np.uint8:
            vtk_type, payload = "UInt8", np.ascontiguousarray(arr, dtype=np.uint8).tobytes()
        else:
            vtk_type, payload = "Float64", np.ascontiguousarray(arr, dtype="<f8").tobytes()
        refs.append(f'<DataArray type="{vtk_type}" Name="{name}" NumberOfComponents="{ncomp}" '
                    f'format="appended" offset="{offset}"/>')
        block = np.uint64(len(payload)).tobytes() + payload
        blocks.append(block)
        offset += len(block)
    return refs, b"".join(blocks)


def write_vtk(vol: SampledVolume, path, *, overwrite: bool = False) -> Path:
    """Write ``vol`` as VTK XML: ``.vti`` (only if axis-aligned) or ``.vts`` (any regular grid).

    Point data: the values (named after ``quantity``; masked voxels written as NaN), ``mask``
    (UInt8, 1 = excluded) and any ``extra`` arrays. A JSON sidecar ``<file>.json`` keeps the
    geometry, units and provenance. Coordinates are Cartesian Å⁻¹ or Å."""
    p = Path(path)
    if p.suffix not in (".vti", ".vts"):
        raise ValueError("path must end in .vti (ImageData) or .vts (StructuredGrid)")
    if p.exists() and not overwrite:
        raise FileExistsError(f"{p} exists; pass overwrite=True")
    if p.suffix == ".vti" and not vol.is_axis_aligned():
        raise ValueError("ImageData (.vti) needs every axis along +x, +y, +z; this grid is oblique or "
                         "non-orthogonal in Cartesian coordinates — write .vts (explicit points) instead")
    n = vol.shape
    extent = f"0 {n[0] - 1} 0 {n[1] - 1} 0 {n[2] - 1}"
    order = lambda a: np.asarray(a).ravel(order="F")        # VTK: x (axis 0) fastest
    arrays = [(vol.quantity, order(vol.filled(np.nan)).astype("<f8"), 1)]
    if vol.mask is not None:
        arrays.append(("mask", order(vol.mask).astype(np.uint8), 1))
    for name, arr in vol.extra.items():
        arrays.append((name, order(arr).astype("<f8"), 1))
    if p.suffix == ".vts":
        pts = vol.points().transpose(2, 1, 0, 3).reshape(-1, 3).astype("<f8")
        arrays.append(("Points", pts, 3))
    refs, blob = _appended(arrays)
    point_refs = "\n".join("        " + r for r in refs[:len(refs) - (p.suffix == ".vts")])
    if p.suffix == ".vti":
        spacing = [vol.steps[k, k] if n[k] > 1 else 1.0 for k in range(3)]
        head = (f'  <ImageData WholeExtent="{extent}" Origin="{" ".join(repr(float(x)) for x in vol.origin)}" '
                f'Spacing="{" ".join(repr(float(s)) for s in spacing)}">\n'
                f'    <Piece Extent="{extent}">\n      <PointData Scalars="{vol.quantity}">\n{point_refs}\n'
                f'      </PointData>\n    </Piece>\n  </ImageData>\n')
        kind = "ImageData"
    else:
        head = (f'  <StructuredGrid WholeExtent="{extent}">\n    <Piece Extent="{extent}">\n'
                f'      <PointData Scalars="{vol.quantity}">\n{point_refs}\n      </PointData>\n'
                f'      <Points>\n        {refs[-1]}\n      </Points>\n    </Piece>\n  </StructuredGrid>\n')
        kind = "StructuredGrid"
    xml = (f'<?xml version="1.0"?>\n<VTKFile type="{kind}" version="1.0" byte_order="LittleEndian" '
           f'header_type="UInt64">\n{head}  <AppendedData encoding="raw">\n   _')
    with open(p, "wb") as fh:
        fh.write(xml.encode("ascii"))
        fh.write(blob)
        fh.write(b"\n  </AppendedData>\n</VTKFile>\n")
    side = {"schema": "discuspdf.SampledVolume.vtk", "version": 1, "space": vol.space, "unit": vol.unit,
            "quantity": vol.quantity, "signed": vol.signed, "shape": list(n), "origin": vol.origin.tolist(),
            "steps": vol.steps.tolist(), "point_order": "axis 0 fastest (VTK x)",
            "index_vectors": None if vol.index_vectors is None else vol.index_vectors.tolist(),
            "metadata": vol.metadata}
    Path(str(p) + ".json").write_text(json.dumps(side, indent=1) + "\n", encoding="utf-8")
    return p


def read_vtk(path) -> dict[str, Any]:
    """Read a file written by :func:`write_vtk` back into numpy (no VTK needed): ``{'kind', 'shape',
    'points' (shape+(3,)), 'arrays': {name: array with the original (i0, i1, i2) layout}}``."""
    raw = Path(path).read_bytes()
    marker = raw.index(b'<AppendedData encoding="raw">')
    start = raw.index(b"_", marker) + 1
    head = raw[:marker].decode("ascii")
    kind = re.search(r'<VTKFile type="(\w+)"', head).group(1)
    ext = [int(x) for x in re.search(r'WholeExtent="([^"]+)"', head).group(1).split()]
    shape = (ext[1] + 1, ext[3] + 1, ext[5] + 1)
    out: dict[str, Any] = {"kind": kind, "shape": shape, "arrays": {}}
    for m in re.finditer(r'<DataArray type="(\w+)" Name="(\w+)" NumberOfComponents="(\d)" format="appended" '
                         r'offset="(\d+)"/>', head):
        typ, name, ncomp, off = m.group(1), m.group(2), int(m.group(3)), int(m.group(4))
        pos = start + off
        nbytes = int(np.frombuffer(raw[pos:pos + 8], "<u8")[0])
        arr = np.frombuffer(raw[pos + 8:pos + 8 + nbytes], {"Float64": "<f8", "UInt8": "u1"}[typ])
        if ncomp == 3:
            out["points"] = arr.reshape(shape[2], shape[1], shape[0], 3).transpose(2, 1, 0, 3)
        else:
            out["arrays"][name] = arr.reshape(shape, order="F")
    if kind == "ImageData":
        o = np.array([float(x) for x in re.search(r'Origin="([^"]+)"', head).group(1).split()])
        s = np.array([float(x) for x in re.search(r'Spacing="([^"]+)"', head).group(1).split()])
        idx = np.stack(np.meshgrid(*[np.arange(n) for n in shape], indexing="ij"), -1)
        out["points"] = o + idx * s
    return out
