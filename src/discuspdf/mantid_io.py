"""Mantid MDHistoWorkspace NeXus files (``SaveMD`` version 2) ↔ :class:`ScatteringVolume` (Package B, 18).

Conventions observed on files written by Mantid 7.0.0 (``tests/fixtures/mantid_7_0_0``) and on two
older real files (Shiver 2023, Javelin):

* ``/MDHistoWorkspace`` (attr ``SaveMDVersion`` = 2) holds ``data/D0…D{n-1}`` (bin **edges**, attrs
  ``frame``, ``long_name``, ``units``), ``data/signal`` with attr ``axes = "D{n-1}:…:D0"`` (the first
  array axis is the **last** dimension), ``errors_squared``, ``mask`` (int8, 1 = masked) and
  ``num_events`` of the same shape.
* ``transform_to_orig`` (16 float32, row-major 4×4 affine): original coordinates = M·[x, 1]. Its
  columns are the projection basis vectors (``BinMD`` with ``[H,H,0]`` gives column (1, 1, 0)).
* ``experiment0/logs/W_MATRIX/value`` (9 values, row-major 3×3; written by ``ConvertToMD``): its
  **columns** are the projection vectors of the original workspace, hkl = W·orig. Absent → identity.
* So for an HKL-frame workspace: **hkl = W · (M₃ · x + t)**.
* ``units`` is ``r.l.u.`` or ``in <L> A^-1``, where L = 2π·|basis vector| in Å⁻¹ (checked).
* ``experiment0/sample/oriented_lattice`` gives a, b, c, α, β, γ and ``orientation_matrix``, which is
  **UB** (UBᵀ·UB equals the reciprocal metric; B in Busing–Levy form, a* along x).

The reader cross-checks three independent records of each axis (affine transform × W_MATRIX, the
axis name such as ``[H,H,0]``, and the Å⁻¹ length in ``units``) and refuses the file if they disagree.
It refuses non-HKL frames (Q_lab, Q_sample), more or fewer than three HKL dimensions, and extra
dimensions (e.g. DeltaE) unless they are integrated (one bin; the range is recorded).
"""
from __future__ import annotations

from pathlib import Path
from typing import Any
import hashlib
import re

import numpy as np

from .lattice import Lattice
from .reciprocal_grid import ReciprocalGrid
from .scattering_data import ScatteringVolume

_ROOT = "MDHistoWorkspace"
_UNITS_RE = re.compile(r"^in\s+([0-9.eE+-]+)\s*A\^-1$")


class MantidFormatError(ValueError):
    """The file is not a Mantid MDHisto file this reader understands, or its metadata disagree."""


def _text(v) -> str:
    if isinstance(v, np.ndarray):
        v = v.ravel()[0] if v.size else ""
    return v.decode("utf-8", "replace") if isinstance(v, (bytes, np.bytes_)) else str(v)


def parse_axis_name(name: str) -> np.ndarray | None:
    """Direction of a Mantid HKL axis name such as ``[H,H,0]``, ``[-H,H,0]``, ``[0,0,L]``, ``[2H,-H,0]``.

    Returns the hkl vector, or ``None`` if the name is not of that form (constant offsets such as
    ``[H,1,0]`` are not interpreted)."""
    m = re.fullmatch(r"\s*\[(.*)\]\s*", name)
    if not m:
        return None
    parts = [p.strip() for p in m.group(1).split(",")]
    if len(parts) != 3:
        return None
    out = []
    for p in parts:
        t = re.fullmatch(r"([+-]?)(\d*\.?\d*(?:[eE][+-]?\d+)?)([HKLhkl]?)", p)
        if not t:
            return None
        sign, num, letter = t.groups()
        if letter:
            value = float(num) if num not in ("", ".") else 1.0
            out.append(-value if sign == "-" else value)
        else:
            if num in ("", ".") or float(num) != 0.0:
                return None
            out.append(0.0)
    return np.array(out)


def _format_axis_name(v: np.ndarray) -> str:
    nz = [k for k in range(3) if v[k] != 0]
    letter = "HKL"[nz[0]] if nz else "H"
    toks = []
    for c in v:
        if c == 0:
            toks.append("0")
        elif c == 1:
            toks.append(letter)
        elif c == -1:
            toks.append("-" + letter)
        else:
            toks.append(f"{c:.6g}{letter}")
    return "[" + ",".join(toks) + "]"


def _busing_levy_b(lattice: Lattice) -> np.ndarray:
    """Mantid's B (Busing & Levy 1967): upper triangular, BᵀB = reciprocal metric, a* along x."""
    return np.linalg.cholesky(lattice.reciprocal_metric).T


def _u_from_ub(UB: np.ndarray, lattice: Lattice):
    U = UB @ np.linalg.inv(_busing_levy_b(lattice))
    return U.tolist() if np.allclose(U @ U.T, np.eye(3), atol=1e-6) else None


def _affine(g, name: str, n: int) -> np.ndarray | None:
    if name not in g:
        return None
    d = g[name]
    rows, cols = int(np.ravel(d.attrs.get("rows", [n + 1]))[0]), int(np.ravel(d.attrs.get("columns", [n + 1]))[0])
    a = np.asarray(d[()], dtype=float)
    if rows * cols != a.size:
        raise MantidFormatError(f"{name}: {a.size} values for a {rows}x{cols} matrix")
    return a.reshape(rows, cols)


def read_mantid_mdhisto(path, *, units: str | None = None) -> ScatteringVolume:
    """Read a Mantid ``SaveMD`` (v2) MDHistoWorkspace in the HKL frame as a :class:`ScatteringVolume`.

    NaN signal values (no coverage, e.g. from ``MDNorm``) are masked with reason
    ``'NaN in Mantid signal'``; Mantid's own mask is kept with reason ``'Mantid mask'``.
    ``errors_squared`` becomes ``variance``."""
    import h5py
    p = Path(path)
    if not p.is_file() or not h5py.is_hdf5(p):
        raise MantidFormatError(f"not an HDF5 file: {p}")
    with h5py.File(p, "r") as f:
        if _ROOT not in f:
            raise MantidFormatError(f"{p.name}: no /{_ROOT} group (an MDEventWorkspace or other file?)")
        g = f[_ROOT]
        version = int(np.ravel(g.attrs.get("SaveMDVersion", [-1]))[0])
        if version != 2:
            raise MantidFormatError(f"{p.name}: SaveMDVersion {version}; only version 2 is understood")
        d = g["data"]
        if "signal" not in d:
            raise MantidFormatError(f"{p.name}: no data/signal")
        axes_attr = _text(d["signal"].attrs.get("axes", b""))
        names_slow_first = [a for a in axes_attr.split(":") if a]
        if not names_slow_first or any(a not in d for a in names_slow_first):
            raise MantidFormatError(f"{p.name}: signal 'axes' attribute {axes_attr!r} does not name the axis datasets")
        dim_names = names_slow_first[::-1]                    # dimension order D0, D1, …
        n = len(dim_names)
        signal = np.asarray(d["signal"][()], dtype=float)
        if signal.ndim != n:
            raise MantidFormatError(f"{p.name}: signal has {signal.ndim} axes, 'axes' names {n}")
        perm = tuple(range(n))[::-1]
        signal = np.transpose(signal, perm)
        arrays = {k: np.transpose(np.asarray(d[k][()]), perm) for k in ("errors_squared", "mask", "num_events")
                  if k in d}
        dims = []
        for k, name in enumerate(dim_names):
            ds = d[name]
            edges = np.asarray(ds[()], dtype=float)
            if edges.ndim != 1 or len(edges) != signal.shape[k] + 1:
                raise MantidFormatError(f"{p.name}: {name} has {len(edges)} edges for {signal.shape[k]} bins")
            dims.append({"dataset": name, "edges": edges, "frame": _text(ds.attrs.get("frame", b"")),
                         "long_name": _text(ds.attrs.get("long_name", b"")), "units": _text(ds.attrs.get("units", b""))})
        to_orig = _affine(g, "transform_to_orig", n)
        w_path = "experiment0/logs/W_MATRIX/value"
        W = np.asarray(g[w_path][()], float).reshape(3, 3) if w_path in g else None
        ol = "experiment0/sample/oriented_lattice"
        if ol not in g:
            raise MantidFormatError(f"{p.name}: no oriented lattice (experiment0/sample/oriented_lattice); "
                                    "hkl cannot be related to Q")
        o = g[ol]
        cell = [float(np.ravel(o[k][()])[0]) for k in ("unit_cell_a", "unit_cell_b", "unit_cell_c",
                                                        "unit_cell_alpha", "unit_cell_beta", "unit_cell_gamma")]
        UB = np.asarray(o["orientation_matrix"][()], float) if "orientation_matrix" in o else None
        root_attrs = {k: _text(v) if not np.issubdtype(np.asarray(v).dtype, np.number) else np.asarray(v).tolist()
                      for k, v in g.attrs.items()}
        extra_meta = {k: np.asarray(g[k][()]).tolist() for k in ("coordinate_system", "visual_normalization") if k in g}

    # ---- dimensions: three HKL axes, anything else must be integrated -------------------
    hkl_idx = [k for k, dm in enumerate(dims) if dm["frame"] == "HKL"]
    for k, dm in enumerate(dims):
        if dm["frame"] in ("QLab", "QSample"):
            raise MantidFormatError(f"{p.name}: axis {dm['long_name']!r} is in frame {dm['frame']}; convert to the "
                                    "HKL frame in Mantid first (this reader does not apply UB/goniometer)")
    if len(hkl_idx) != 3:
        raise MantidFormatError(f"{p.name}: {len(hkl_idx)} HKL axes; a reciprocal-space volume needs exactly 3 "
                                f"(axes: {[dm['long_name'] for dm in dims]})")
    integrated = []
    for k, dm in enumerate(dims):
        if k in hkl_idx:
            continue
        if signal.shape[k] != 1:
            raise MantidFormatError(f"{p.name}: non-HKL axis {dm['long_name']!r} ({dm['frame']}) has "
                                    f"{signal.shape[k]} bins; integrate or slice it in Mantid first")
        integrated.append({"name": dm["long_name"], "frame": dm["frame"], "units": dm["units"],
                           "range": [float(dm["edges"][0]), float(dm["edges"][-1])]})
    if to_orig is not None:
        M = to_orig[:n, :n]
        t = to_orig[:n, n]
        other = [k for k in range(n) if k not in hkl_idx]
        if other and (np.any(M[np.ix_(hkl_idx, other)] != 0) or np.any(M[np.ix_(other, hkl_idx)] != 0)):
            raise MantidFormatError(f"{p.name}: HKL axes are coupled to non-HKL axes in transform_to_orig")
        M3, t3 = M[np.ix_(hkl_idx, hkl_idx)], t[hkl_idx]
    else:
        M3, t3 = np.eye(3), np.zeros(3)
    W3 = np.eye(3) if W is None else W
    P = W3 @ M3                                   # columns: hkl of a unit step of each axis
    shift = W3 @ t3

    # ---- grid and cross-checks -------------------------------------------------------------
    lattice = Lattice.from_parameters(*cell)
    B2pi = 2 * np.pi * lattice.reciprocal_matrix
    checks = []
    steps, centres0 = [], []
    for col, k in enumerate(hkl_idx):
        dm = dims[k]
        e = dm["edges"]
        w = np.diff(e)
        if np.any(w <= 0) or not np.allclose(w, w[0], rtol=1e-6, atol=0):
            raise MantidFormatError(f"{p.name}: axis {dm['long_name']!r} has non-uniform or decreasing bins")
        steps.append(float(w[0]))
        centres0.append(float(e[0] + w[0] / 2))
        basis = P[:, col]
        named = parse_axis_name(dm["long_name"])
        if named is not None:
            cos = abs(named @ basis) / (np.linalg.norm(named) * np.linalg.norm(basis))
            if not np.allclose(named, basis, rtol=1e-5, atol=1e-6):
                raise MantidFormatError(f"{p.name}: axis name {dm['long_name']!r} says {named.tolist()} but the "
                                        f"stored transform gives {basis.round(6).tolist()} (cos {cos:.6f})")
            checks.append(f"{dm['long_name']}: name matches transform")
        else:
            checks.append(f"{dm['long_name']}: name not interpretable, not checked")
        um = _UNITS_RE.match(dm["units"])
        if um:
            length = float(np.linalg.norm(basis @ B2pi))
            stated = float(um.group(1))
            if abs(length - stated) > 0.006 * stated + 0.0005:       # Mantid prints 3 significant digits
                raise MantidFormatError(f"{p.name}: axis {dm['long_name']!r} units say {stated} Å⁻¹ per unit, "
                                        f"the lattice and transform give {length:.4g}")
            checks.append(f"{dm['long_name']}: units length {stated} matches {length:.4g}")
        elif dm["units"] != "r.l.u.":
            raise MantidFormatError(f"{p.name}: axis {dm['long_name']!r} has units {dm['units']!r}; "
                                    "only 'r.l.u.' and 'in <L> A^-1' are understood")
    axis_vectors = (P * np.asarray(steps)[None, :]).T            # row = hkl step of each HKL axis
    corner = shift + P @ np.asarray(centres0)

    # ---- values ---------------------------------------------------------------------------
    sel = tuple(slice(None) if k in hkl_idx else 0 for k in range(n))
    values = signal[sel]
    variance = arrays["errors_squared"][sel].astype(float) if "errors_squared" in arrays else None
    mask = arrays["mask"][sel].astype(bool) if "mask" in arrays else np.zeros(values.shape, bool)
    reasons = ["Mantid mask"] if mask.any() else []
    nan = ~np.isfinite(values)
    if variance is not None:
        nan |= ~np.isfinite(variance)
    if nan.any():
        mask = mask | nan
        reasons.append("NaN in Mantid signal (no coverage)")
    grid = ReciprocalGrid(lattice, corner, axis_vectors, values.shape)
    prov: dict[str, Any] = {
        "source": "Mantid MDHistoWorkspace (SaveMD v2)", "reader": "discuspdf.read_mantid_mdhisto",
        "path": str(p), "sha256": hashlib.sha256(p.read_bytes()).hexdigest(),
        "axes": [{k2: v for k2, v in dm.items() if k2 != "edges"} for dm in dims],
        "integrated_axes": integrated, "transform_to_orig": None if to_orig is None else to_orig.tolist(),
        "W_MATRIX": None if W is None else W.tolist(),
        "UB": None if UB is None else UB.tolist(), "U": None if UB is None else _u_from_ub(UB, lattice),
        "cell": cell, "root_attrs": root_attrs, **extra_meta, "consistency_checks": checks,
        "num_events_total": float(np.nansum(arrays["num_events"])) if "num_events" in arrays else None,
        "lattice_orientation": "canonical (a along x); Mantid's UB/U are recorded, not applied"}
    if variance is not None:
        variance = np.where(mask, np.nan, variance)          # undefined where excluded
    return ScatteringVolume(grid, values, units or "Mantid signal", bragg="unknown", variance=variance,
                            mask=mask if mask.any() else None, mask_reasons=tuple(reasons), provenance=prov)


def _group(parent, group_name, nx_class, **attrs):
    grp = parent.create_group(group_name)
    grp.attrs["NX_class"] = np.bytes_(nx_class.encode())
    for k, v in attrs.items():
        grp.attrs[k] = np.int32(v) if isinstance(v, int) else np.bytes_(str(v).encode())
    return grp


def _write_experiment(root, lattice: Lattice) -> None:
    """experiment0 with the same skeleton Mantid 7.0.0 SaveMD writes (LoadMD needs it) and an oriented
    lattice whose UB = U·B with U = identity and Mantid's Busing–Levy B."""
    import datetime as dt
    s1 = lambda text: np.array([text.encode()])
    e = _group(root, "experiment0", "NXgroup", version=1)
    inst = _group(e, "instrument", "NXinstrument", version=1)
    inst.create_dataset("name", data=s1(" "))
    xml = _group(inst, "instrument_xml", "NXnote")
    xml.create_dataset("data", data=s1(" "))
    xml.create_dataset("description", data=s1("XML contents of the instrument IDF file."))
    xml.create_dataset("type", data=s1("text/xml"))
    pm = _group(inst, "instrument_parameter_map", "NXnote", version=1)
    pm.create_dataset("author", data=s1(" "))
    pm.create_dataset("data", data=s1(" "))
    pm.create_dataset("date", data=s1(dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%S")))
    pm.create_dataset("description", data=s1("A string representation of the parameter map."))
    pm.create_dataset("type", data=s1("text/plain"))
    logs = _group(e, "logs", "NXgroup", version=1)
    gon = _group(logs, "goniometer", "NXpositioner", version=1)
    gon.create_dataset("num_axes", data=np.array([0], np.int32))
    gon.create_dataset("rotation_matrix", data=np.eye(3).ravel())
    note = _group(logs, "written_by", "NXlog")
    note.create_dataset("value", data=s1("discuspdf write_mantid_mdhisto"))
    note["value"].attrs["units"] = np.bytes_(b" ")
    smp = _group(e, "sample", "NXsample", version=1, name=" ", name_empty=1, shape_id=" ",
                 shape_xml='<type name="userShape">  </type>')
    for k in ("geom_height", "geom_thickness", "geom_width"):
        smp.create_dataset(k, data=np.array([0.0]))
    smp.create_dataset("geom_id", data=np.array([0], np.int32))
    mat = _group(smp, "material", "NXdata", version=2, formulaStyle="empty", name=" ")
    for k, v in (("number_density", 0.0), ("packing_fraction", 1.0), ("pressure", 0.0), ("temperature", 0.0)):
        mat.create_dataset(k, data=np.array([v]))
    smp.create_dataset("num_oriented_lattice", data=np.array([1], np.int32))
    smp.create_dataset("num_other_samples", data=np.array([0], np.int32))
    ol = _group(smp, "oriented_lattice", "NXcrystal")
    ol.create_dataset("cross_term", data=np.array([0], np.int32))
    ol.create_dataset("maximum_order", data=np.array([0], np.int32))
    ol.create_dataset("modulated_hkl_error", data=np.zeros((3, 3)))
    ol.create_dataset("modulated_orientation_matrix", data=np.zeros((3, 3)))
    ol.create_dataset("orientation_matrix", data=_busing_levy_b(lattice))
    for key, val in zip(("a", "b", "c", "alpha", "beta", "gamma"), lattice.parameters):
        ol.create_dataset(f"unit_cell_{key}", data=np.array([val]))
        ol.create_dataset(f"unit_cell_{key}_error", data=np.array([0.0]))


def write_mantid_mdhisto(sv: ScatteringVolume, path, *, overwrite: bool = False) -> Path:
    """Write a 3D :class:`ScatteringVolume` as a ``SaveMD`` v2 MDHistoWorkspace that Mantid's ``LoadMD``
    reads (checked with Mantid 7.0.0).

    Each array axis becomes one HKL dimension whose basis vector is the grid step scaled so its largest
    component is ±1 (``[H,H,0]``-style names, units ``r.l.u.``); ``transform_to_orig`` holds those basis
    vectors, so :func:`read_mantid_mdhisto` gives back the same grid. Variance → ``errors_squared``
    (zeros if absent), mask → ``mask``; masked values are written as 0. Planes and lines are refused
    (a one-point axis has no direction)."""
    import h5py
    p = Path(path)
    if p.exists() and not overwrite:
        raise FileExistsError(f"{p} exists; pass overwrite=True")
    g = sv.grid
    if g.ndim != 3:
        raise ValueError("write_mantid_mdhisto needs a 3D grid (a one-point axis has no direction to store)")
    V = g.axis_vectors
    scale = np.array([V[k][np.argmax(np.abs(V[k]))] for k in range(3)])  # signed largest component
    basis = (V / scale[:, None]).T                                       # columns = basis vectors
    x0 = np.linalg.solve(basis, g.corner)                                # corner = basis @ x0
    edges = [x0[k] - scale[k] / 2 + scale[k] * np.arange(g.shape[k] + 1) for k in range(3)]
    for k in range(3):
        if scale[k] < 0:                    # keep edges increasing: flip the basis vector instead
            basis[:, k] *= -1
            edges[k] = -edges[k]
    if any(np.any(np.diff(e) <= 0) for e in edges):
        raise ValueError("internal error: bin edges are not increasing")
    to_orig = np.eye(4)
    to_orig[:3, :3] = basis
    from_orig = np.linalg.inv(to_orig)
    perm = (2, 1, 0)
    mask = np.zeros(g.shape, bool) if sv.mask is None else sv.mask
    signal = np.where(mask, 0.0, np.nan_to_num(sv.intensity))
    err2 = np.zeros(g.shape) if sv.variance is None else np.where(mask, 0.0, np.nan_to_num(sv.variance))
    with h5py.File(p, "w") as f:
        root = f.create_group(_ROOT)
        root.attrs["NX_class"] = np.bytes_(b"NXentry")
        root.attrs["SaveMDVersion"] = np.int32(2)
        root.attrs["QConvention"] = np.bytes_(b"Inelastic")
        root.create_dataset("coordinate_system", data=np.array([3], np.uint32))     # HKL
        root.create_dataset("visual_normalization", data=np.array([0], np.uint32))  # no normalisation
        for name, mat, kind in (("transform_to_orig", to_orig, b"CoordTransformAffine"),
                                ("transform_from_orig", from_orig, b"CoordTransformAffine")):
            ds = root.create_dataset(name, data=mat.ravel().astype(np.float32))
            ds.attrs["type"] = np.bytes_(kind)
            ds.attrs["rows"] = np.int32(4)
            ds.attrs["columns"] = np.int32(4)
        d = root.create_group("data")
        d.attrs["NX_class"] = np.bytes_(b"NXdata")
        for k in range(3):
            ds = d.create_dataset(f"D{k}", data=edges[k])
            ds.attrs["frame"] = np.bytes_(b"HKL")
            ds.attrs["long_name"] = np.bytes_(_format_axis_name(basis[:, k]).encode())
            ds.attrs["units"] = np.bytes_(b"r.l.u.")
        s = d.create_dataset("signal", data=np.transpose(signal, perm))
        s.attrs["signal"] = np.int32(1)
        s.attrs["axes"] = np.bytes_(b"D2:D1:D0")
        d.create_dataset("errors_squared", data=np.transpose(err2, perm))
        d.create_dataset("mask", data=np.transpose(mask, perm).astype(np.int8))
        d.create_dataset("num_events", data=np.zeros(tuple(g.shape[::-1])))
        _write_experiment(root, g.lattice)
    return p
