"""Result files of a disorder-builder run, in the formats the model selects (:data:`.model.EXPORT_DEFAULTS`).

========  ===========================================================================================
format    files
========  ===========================================================================================
hdf5      ``diffuse.h5`` / ``total.h5`` (:meth:`ScatteringVolume.save_hdf5`), ``delta_pdf.h5`` / ``pdf3d.h5``
csv       ``pdf_1d.csv`` (``r (Angstrom),G(r)``, the run's ``result.csv``)
png       ``*_slices.png`` (three orthogonal slices through the origin), ``pdf_1d.png``
vtk       ``delta_pdf.vts``, ``pdf3d_patterson.vts``
npz       ``diffuse.npz`` / ``total.npz`` / ``delta_pdf.npz`` (numpy only; the files of builder engine 1.0)
txt       ``pdf_1d.txt`` (``r(A) G(r)``, space-separated)
========  ===========================================================================================
"""
from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np

__all__ = ["save_volume_npz", "save_real_space_hdf5", "save_slices", "write_volume_outputs",
           "write_pdf_outputs", "preview_factors", "numpy_preview"]


def save_volume_npz(vol, path: Path) -> None:
    """Intensity with its hkl axes (r.l.u.) and the unit cell, readable with numpy alone."""
    h, k, l = vol.grid.hkl_axes()
    np.savez_compressed(path, intensity=np.asarray(vol.intensity), h=h, k=k, l=l,
                        cell=np.array(vol.grid.lattice.parameters), bragg=str(vol.bragg))


def save_real_space_hdf5(res, path: Path, *, quantity: str) -> None:
    """A fractional-FFT result (3D-ΔPDF or Patterson): real part, fractional axes u, v, w, direct basis (Å)."""
    import h5py
    with h5py.File(path, "w") as f:
        f.attrs["quantity"] = quantity
        f.attrs["normalization"] = str(res.normalization)
        f.attrs["imaginary_ratio"] = float(res.imaginary_ratio)
        f.attrs["generator"] = "discuspdf.builder"
        f.create_dataset("values", data=np.asarray(res.complex_volume.real), compression="gzip")
        for name, ax in zip(("u", "v", "w"), res.u_axes):
            f.create_dataset(name, data=np.asarray(ax, float))
        f.create_dataset("direct_basis", data=np.asarray(res.direct_basis, float))


def save_slices(vol, title: str, png: Path, log) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from ..viz import plot_slice
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
    for ax_i, (axis, name) in enumerate([(2, "l = 0 / w = 0"), (1, "k = 0 / v = 0"), (0, "h = 0 / u = 0")]):
        plot_slice(vol, axis=axis, ax=axes[ax_i], title=f"{title}: {name}")
    fig.tight_layout()
    fig.savefig(png, dpi=130)
    plt.close(fig)
    log(f"  wrote {png.name}")


def write_volume_outputs(vol, out: Path, name: str, formats: dict, log, *, delta_pdf: bool, real_name: str,
                         real_quantity: str, require_diffuse: bool = True) -> None:
    """Files for one reciprocal volume (``diffuse`` or ``total``) and its real-space transform."""
    from ..viz import SampledVolume, write_vtk
    if formats["hdf5"]:
        vol.save_hdf5(out / f"{name}.h5", overwrite=True)
        log(f"  wrote {name}.h5")
    if formats["npz"]:
        save_volume_npz(vol, out / f"{name}.npz")
    if formats["png"]:
        save_slices(SampledVolume.from_scattering_volume(vol), name, out / f"{name}_slices.png", log)
    if not delta_pdf:
        return
    res = vol.delta_pdf(require_diffuse=require_diffuse)
    rho = SampledVolume.from_fractional_fft(res, quantity=real_quantity)
    log(f"  {real_name}: imaginary/real = {res.imaginary_ratio:.1e}")
    if formats["hdf5"]:
        stem = "delta_pdf" if real_quantity == "delta_pdf" else "pdf3d"
        save_real_space_hdf5(res, out / f"{stem}.h5", quantity=real_quantity)
    if formats["npz"] and real_quantity == "delta_pdf":
        np.savez_compressed(out / "delta_pdf.npz", values=res.complex_volume.real,
                            u=res.u_axes[0], v=res.u_axes[1], w=res.u_axes[2],
                            direct_basis=np.asarray(res.direct_basis))
    if formats["vtk"]:
        write_vtk(rho, out / ("delta_pdf.vts" if real_quantity == "delta_pdf" else "pdf3d_patterson.vts"),
                  overwrite=True)
    if formats["png"]:
        stem = "delta_pdf" if real_quantity == "delta_pdf" else "pdf3d"
        save_slices(rho, real_name, out / f"{stem}_slices.png", log)


def write_pdf_outputs(run, out: Path, formats: dict, title: str, log) -> tuple[np.ndarray, np.ndarray]:
    r, g = run.load()
    written = []
    if formats["csv"]:
        shutil.copy2(run.output, out / "pdf_1d.csv")
        written.append("pdf_1d.csv")
    if formats["txt"]:
        np.savetxt(out / "pdf_1d.txt", np.column_stack([r, g]), header="r(A) G(r)")
        written.append("pdf_1d.txt")
    if formats["png"]:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.plot(r, g, lw=1)
        ax.set_xlabel("r (Å)")
        ax.set_ylabel("G(r)")
        ax.set_title(title)
        fig.tight_layout()
        fig.savefig(out / "pdf_1d.png", dpi=130)
        plt.close(fig)
        written.append("pdf_1d.png")
    if written:
        log("  wrote " + ", ".join(written))
    return r, g


def preview_factors(cfg: dict, species, radiation: str):
    """Scattering factors for the numpy preview and a description: neutron b from DISCUS's table (or the
    model's own lengths); X-ray f₀(s) from ``xraydb`` if installed, else f = Z (stated in the result)."""
    from ..atom_names import atomic_number, discus_atom_name, element_of
    from ..reference_scattering import ScatteringFactors
    from .runs import scattering_lengths_for
    from .model import ModelError
    names = sorted({str(s) for s in species})
    real = [s for s in names if s.upper() != "VOID"]
    unreadable = [s for s in real if element_of(s) is None]
    if unreadable:
        raise ModelError(f"atom name(s) {unreadable} are not elements DISCUS can read; no scattering factor")
    if radiation == "neutron":
        own = scattering_lengths_for(cfg, real, "neutron")
        table = {s: own.get(s.upper(), discus_atom_name(s).b) for s in real}
        table.update({s: 0.0 for s in names if s.upper() == "VOID"})
        return ScatteringFactors.neutron(table, "DISCUS 6.23.01 per_bcoh (+ model's own b)"), "neutron b (DISCUS table)"
    zero = lambda s: np.zeros(np.shape(s))   # noqa: E731
    try:
        import xraydb
        funcs = {s: (lambda q, el=element_of(s): np.asarray(xraydb.f0(el, np.asarray(q, float).ravel()),
                                                          float).reshape(np.shape(q))) for s in real}
        label = "X-ray f0(s) (xraydb, Waasmaier-Kirfel)"
    except ImportError:
        funcs = {s: (lambda q, z=float(atomic_number(s)): np.full(np.shape(q), z)) for s in real}
        label = "f = Z (install xraydb for f0(s))"
    funcs.update({s: zero for s in names if s.upper() == "VOID"})
    return ScatteringFactors.xray_functions(funcs, label), label


def numpy_preview(cfg, structure, info, out: Path, log, plane=None):
    """Diffuse intensity |F − ⟨F⟩|² per unit cell in one reciprocal plane, from
    :mod:`discuspdf.reference_scattering` (no DISCUS, no Debye–Waller). Radiation as in the 3D settings."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from ..reference_scattering import diffuse_intensity
    plane = plane or cfg.get("calc", {}).get("preview", {"u": [1, 0, 0], "v": [0, 1, 0], "range": 6})
    u, v, R = np.array(plane["u"], float), np.array(plane["v"], float), float(plane["range"])
    o = np.array(plane.get("origin") or [0, 0, 0], float)      # plane offset, e.g. (0, 0, 0.5) for hk½
    n = np.array(structure.supercell)
    radiation = (cfg.get("calc", {}).get("fourier") or {}).get("radiation") or "xray"
    factors, label = preview_factors(cfg, structure.species, radiation)
    step = 1.0 / float(max(n))                         # commensurate with the supercell: no size ripples
    t = np.arange(-round(R / step), round(R / step) + 1) * step
    H = o[None, None, :] + t[:, None, None] * u[None, None, :] + t[None, :, None] * v[None, None, :]
    intensity = diffuse_intensity(structure, factors, hkl=H, normalize="per_cell")
    np.savez(out / "numpy_preview_plane.npz", t=t, u=u, v=v, origin=o, intensity=intensity, factors=label)
    B = 2 * np.pi * np.linalg.inv(info["unit_lattice"]).T
    gu, gv = u @ B, v @ B
    e1 = gu / np.linalg.norm(gu)
    e2 = gv - (gv @ e1) * e1
    e2 /= np.linalg.norm(e2)
    TA, TB = np.meshgrid(t, t, indexing="ij")
    X = TA * (gu @ e1) + TB * (gv @ e1)
    Y = TA * (gu @ e2) + TB * (gv @ e2)
    fig, ax = plt.subplots(figsize=(6, 5.2))
    vmax = np.percentile(np.abs(intensity), 99.5) or 1.0
    im = ax.pcolormesh(X, Y, intensity, shading="nearest", cmap="magma", vmin=0, vmax=vmax)
    ax.set_aspect("equal")
    lab = lambda w: "[" + " ".join(f"{x:g}" for x in w) + "]"   # noqa: E731
    ax.set_xlabel(f"Q along u = {lab(u)} (Å⁻¹)")
    ax.set_ylabel(f"Q in the plane, normal to u (Å⁻¹); v = {lab(v)}")
    ax.set_title(f"numpy preview: diffuse intensity ({label})" + (f", offset {lab(o)}" if np.any(o) else ""))
    fig.colorbar(im, ax=ax)
    fig.tight_layout()
    fig.savefig(out / "numpy_preview_plane.png", dpi=130)
    plt.close(fig)
    log(f"  wrote numpy_preview_plane.png ({label}; true plane geometry)")
    return intensity
