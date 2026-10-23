#!/usr/bin/env python3
"""Describe a DISCUS output file WITHOUT assuming any layout.

Usage:  python inspect_output.py <file> [--json out.json]

For HDF5: lists every group/dataset/attribute (shape, dtype, min/max/finite
counts, monotonicity and uniform spacing of 1D datasets). For anything else it
reports size, magic bytes and the first lines so the real format can be
identified. It never interprets units or coordinates; it only reports facts.
"""
import hashlib, json, sys
from pathlib import Path


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _attr(v):
    try:
        import numpy as np
        if isinstance(v, bytes):
            return v.decode("utf-8", "replace")
        if isinstance(v, np.ndarray):
            return v.tolist() if v.size <= 32 else f"<array shape={v.shape} dtype={v.dtype}>"
        if isinstance(v, np.generic):
            return v.item()
    except Exception:
        pass
    return str(v)


def describe_h5(path):
    import h5py, numpy as np
    out = {"format": "hdf5", "groups": [], "datasets": {}, "root_attrs": {}}
    with h5py.File(path, "r") as f:
        out["root_attrs"] = {k: _attr(v) for k, v in f.attrs.items()}

        def visit(name, obj):
            if isinstance(obj, h5py.Group):
                out["groups"].append({"name": "/" + name,
                                      "attrs": {k: _attr(v) for k, v in obj.attrs.items()}})
                return
            info = {"shape": list(obj.shape), "dtype": str(obj.dtype),
                    "attrs": {k: _attr(v) for k, v in obj.attrs.items()}}
            if obj.dtype.kind in "fiu" and obj.size > 0:
                a = np.asarray(obj[...], dtype=float)
                fin = np.isfinite(a)
                info.update(n=int(a.size), n_finite=int(fin.sum()),
                            n_nan=int(np.isnan(a).sum()), n_inf=int(np.isinf(a).sum()))
                if fin.any():
                    info.update(min=float(a[fin].min()), max=float(a[fin].max()),
                                mean=float(a[fin].mean()))
                if a.ndim == 1 and a.size > 1 and fin.all():
                    d = np.diff(a)
                    info.update(monotonic_increasing=bool(np.all(d > 0)),
                                monotonic_decreasing=bool(np.all(d < 0)),
                                step_mean=float(d.mean()),
                                step_relspread=float((d.max() - d.min()) / abs(d.mean())) if d.mean() else None,
                                first=float(a[0]), last=float(a[-1]))
                if a.ndim == 3 and fin.all():
                    info["max_at_index"] = [int(i) for i in np.unravel_index(np.argmax(a), a.shape)]
                    info["center_index"] = [s // 2 for s in a.shape]
            elif obj.dtype.kind in "SO" and obj.size <= 16:
                info["values"] = [_attr(x) for x in obj[...].ravel()]
            out["datasets"]["/" + name] = info

        f.visititems(visit)
    return out


def describe_other(path):
    p = Path(path)
    raw = p.read_bytes()[:4096]
    return {"format": "not-hdf5", "magic_hex": raw[:16].hex(),
            "printable_head": raw[:600].decode("utf-8", "replace"),
            "size": p.stat().st_size}


def main(argv):
    if len(argv) < 2:
        print(__doc__); return 2
    path = Path(argv[1])
    rep = {"file": str(path), "size_bytes": path.stat().st_size, "sha256": sha256(path)}
    try:
        import h5py
        is_h5 = h5py.is_hdf5(path)
    except ImportError:
        is_h5 = False
        rep["warning"] = "h5py not importable; HDF5 not inspected"
    rep.update(describe_h5(path) if is_h5 else describe_other(path))
    text = json.dumps(rep, indent=2, default=str)
    if "--json" in argv:
        Path(argv[argv.index("--json") + 1]).write_text(text)
    print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
