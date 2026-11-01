"""Atomic-structure pictures and viewer export (Package D, Feature 21).

Notebook-friendly matplotlib views of a :class:`~discuspdf.structure.Structure`: atoms coloured by species or
site, displacement arrows against a reference structure, and bonds of selected neighbor pairs. Export to
established viewers (VESTA, OVITO, Jmol, ASE's GUI…) goes through ASE (``pip install 'discuspdf-starter[ase]'``).

Nothing here changes coordinates: positions are drawn and written exactly as stored (Cartesian =
``frac @ lattice.matrix``, no wrapping, no re-orientation of the cell). matplotlib and ASE are imported
inside the functions, so headless scientific runs never load them.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from .structure import Structure

__all__ = ["species_colors", "plot_structure", "plot_displacements", "plot_bonds", "export_structure",
           "displacement_vectors"]

# Fixed, colour-blind-friendly cycle (Okabe–Ito) so the same species keeps its colour across figures.
_CYCLE = ("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00", "#F0E442", "#000000")


def species_colors(labels) -> dict[str, str]:
    """Deterministic colour per label (sorted order, Okabe–Ito cycle)."""
    return {lab: _CYCLE[k % len(_CYCLE)] for k, lab in enumerate(sorted(set(labels)))}


def _mask(structure: Structure, select) -> np.ndarray:
    n = structure.n_atoms
    if select is None:
        return np.ones(n, bool)
    sel = np.asarray(select)
    if sel.dtype == bool:
        if sel.shape != (n,):
            raise ValueError(f"boolean select must have shape ({n},)")
        return sel
    ids = np.asarray(sel, np.int64).ravel()
    known = np.isin(ids, structure.atom_id)
    if not known.all():
        raise ValueError(f"select contains unknown atom_id values {ids[~known][:5].tolist()}")
    return np.isin(structure.atom_id, ids)


def _axes3d(ax):
    import matplotlib.pyplot as plt
    if ax is None:
        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_subplot(projection="3d")
    elif getattr(ax, "name", "") != "3d":
        raise ValueError("ax must be a 3D axes (fig.add_subplot(projection='3d'))")
    return ax


def _cell_edges(matrix: np.ndarray, origin=np.zeros(3)) -> list[np.ndarray]:
    corners = np.array([[i, j, k] for i in (0, 1) for j in (0, 1) for k in (0, 1)], float)
    edges = []
    for a in range(8):
        for b in range(a + 1, 8):
            if np.abs(corners[a] - corners[b]).sum() == 1:
                edges.append(np.array([corners[a], corners[b]]) @ matrix + origin)
    return edges


def _equal_aspect(ax, points: np.ndarray) -> None:
    lo, hi = points.min(axis=0), points.max(axis=0)
    span = max(float((hi - lo).max()), 1e-9)
    mid = (lo + hi) / 2
    ax.set_xlim(mid[0] - span / 2, mid[0] + span / 2)
    ax.set_ylim(mid[1] - span / 2, mid[1] + span / 2)
    ax.set_zlim(mid[2] - span / 2, mid[2] + span / 2)
    try:
        ax.set_box_aspect((1, 1, 1))
    except AttributeError:                                         # pragma: no cover - old matplotlib
        pass


def plot_structure(structure: Structure, *, ax=None, color_by: str = "species", select=None,
                   cell: str | None = "supercell", colors: dict | None = None, size: float = 40.0,
                   legend: bool = True, alpha: float = 0.9):
    """Scatter the atoms in Cartesian Å, coloured by ``"species"`` or ``"site"`` (site label).

    ``select``: boolean mask or ``atom_id`` values to draw (others are skipped). ``cell``: draw the
    ``"supercell"`` frame, the ``"parent"`` cell at the origin, or ``None``. Marker area is scaled by
    occupancy, so partially occupied sites of an average model look lighter. Returns the 3D axes."""
    if not isinstance(structure, Structure):
        raise TypeError("structure must be a Structure")
    if color_by not in ("species", "site"):
        raise ValueError("color_by must be 'species' or 'site'")
    if cell not in ("supercell", "parent", None):
        raise ValueError("cell must be 'supercell', 'parent' or None")
    ax = _axes3d(ax)
    m = _mask(structure, select)
    xyz = structure.cartesian
    labels = np.asarray(structure.species if color_by == "species" else structure.site_label, dtype=object)
    cmap = species_colors(labels) if colors is None else dict(colors)
    missing = sorted(set(labels[m]) - set(cmap))
    if missing:
        raise ValueError(f"no colour given for {missing}")
    occ = structure.occupancy
    for lab in sorted(set(labels[m])):
        k = m & (labels == lab)
        ax.scatter(*xyz[k].T, s=size * occ[k], c=cmap[lab], label=str(lab), alpha=alpha,
                   depthshade=True, edgecolors="none")
    pts = [xyz[m]] if m.any() else [np.zeros((1, 3))]
    if cell is not None:
        A = structure.lattice.matrix if cell == "supercell" else structure.lattice.matrix / structure.supercell[:, None]
        for e in _cell_edges(A):
            ax.plot(*e.T, color="0.5", lw=0.8)
            pts.append(e)
    _equal_aspect(ax, np.vstack(pts))
    ax.set_xlabel("x (Å)")
    ax.set_ylabel("y (Å)")
    ax.set_zlabel("z (Å)")
    if legend and m.any():
        ax.legend(loc="upper left", fontsize="small", markerscale=1.0)
    return ax


def displacement_vectors(configuration: Structure, reference: Structure) -> tuple[np.ndarray, np.ndarray]:
    """(reference Cartesian positions, Cartesian displacements) of ``configuration`` against ``reference``.

    Atoms are matched by ``(site_label, cell_index)``; both structures need the same lattice, supercell
    and pbc. The displacement is the minimum image in the supercell when all pbc are True, so an atom
    that crossed the boundary does not show a supercell-long arrow."""
    for name, s in (("configuration", configuration), ("reference", reference)):
        if not isinstance(s, Structure):
            raise TypeError(f"{name} must be a Structure")
    if not configuration.lattice.isclose(reference.lattice):
        raise ValueError("configuration and reference have different lattices")
    if not np.array_equal(configuration.supercell, reference.supercell) or configuration.pbc != reference.pbc:
        raise ValueError("configuration and reference have different supercells or pbc")
    key = {}
    for row, (site, ci) in enumerate(zip(reference.site_label, map(tuple, reference.cell_index))):
        if (site, ci) in key:
            raise ValueError(f"reference has two atoms on site {site!r} in cell {ci}")
        key[(site, ci)] = row
    rows = []
    for site, ci in zip(configuration.site_label, map(tuple, configuration.cell_index)):
        if (site, ci) not in key:
            raise ValueError(f"configuration atom on site {site!r} in cell {ci} has no reference site")
        rows.append(key[(site, ci)])
    rows = np.asarray(rows, np.int64)
    if len(np.unique(rows)) != len(rows):
        raise ValueError("two configuration atoms map to the same reference site")
    df = configuration.frac - reference.frac[rows]
    if all(configuration.pbc):
        df = df - np.round(df)
    disp = df @ configuration.lattice.matrix
    return reference.cartesian[rows], disp


def plot_displacements(configuration: Structure, reference: Structure, *, ax=None, scale: float = 1.0,
                       select=None, color_by: str = "species", colors: dict | None = None,
                       draw_atoms: bool = True, cell: str | None = "supercell"):
    """Arrows from each reference site to its atom in ``configuration``, multiplied by ``scale``.

    The scale is stated in the axes title so a magnified picture is never mistaken for real distances."""
    if scale <= 0 or not np.isfinite(scale):
        raise ValueError("scale must be positive")
    start, disp = displacement_vectors(configuration, reference)
    ax = _axes3d(ax)
    m = _mask(configuration, select)
    labels = np.asarray(configuration.species if color_by == "species" else configuration.site_label, dtype=object)
    cmap = species_colors(labels) if colors is None else dict(colors)
    if draw_atoms:
        plot_structure(reference, ax=ax, cell=cell, size=12.0, legend=False, alpha=0.35,
                       color_by="species", colors=species_colors(reference.species))
    handles = []
    for lab in sorted(set(labels[m])):
        k = m & (labels == lab)
        handles.append(ax.quiver(*start[k].T, *(scale * disp[k]).T, color=cmap[lab], label=str(lab),
                                 arrow_length_ratio=0.25, linewidth=1.2))
    ax.set_title(f"displacements × {scale:g}")
    if handles:
        ax.legend(handles=handles, loc="upper left", fontsize="small")
    return ax


def plot_bonds(neighbors, *, ax=None, color: str | None = None, linewidth: float = 1.0, draw_atoms: bool = True):
    """Draw the pairs of a (selected) :class:`~discuspdf.neighbors.NeighborList` as segments.

    Each segment runs from the origin node to the origin + pair vector, i.e. to the periodic image that
    makes the pair; it may therefore leave the supercell frame. Pairs are drawn once per direction."""
    from .neighbors import NeighborList
    if not isinstance(neighbors, NeighborList):
        raise TypeError("neighbors must be a NeighborList (e.g. nl.select(shell=1))")
    st = neighbors.structure
    ax = _axes3d(ax)
    if draw_atoms:
        plot_structure(st, ax=ax, size=25.0, alpha=0.6)
    origin = neighbors.node_frac[neighbors.i] @ st.lattice.matrix
    end = origin + neighbors.vector
    shells = np.unique(np.round(neighbors.distance, 6))
    palette = species_colors([f"{d:.6f}" for d in shells])
    for d in shells:
        k = np.isclose(neighbors.distance, d, atol=1e-6)
        c = color or palette[f"{d:.6f}"]
        for a, b in zip(origin[k], end[k]):
            ax.plot(*np.array([a, b]).T, color=c, lw=linewidth)
        ax.plot([], [], color=c, lw=linewidth, label=f"{d:.3f} Å ({int(k.sum())} pairs)")
    _equal_aspect(ax, np.vstack([origin, end, st.cartesian]) if len(origin) else st.cartesian)
    ax.legend(loc="upper left", fontsize="small")
    return ax


def export_structure(structure: Structure, path, *, format: str | None = None, check: bool = True) -> Path:
    """Write the structure for an external viewer through ASE (``extxyz``, ``cif``, ``vasp``, ``xsf``, …).

    The cell is written as stored (rows a, b, c of ``lattice.matrix`` in Å) and positions as
    ``frac @ lattice.matrix``; nothing is rotated or wrapped. With ``check=True`` the file is read back
    with ASE and cell, positions (to the format's precision, 1e-5 Å) and element symbols are compared;
    formats ASE cannot read back are refused with ``check=True``."""
    from .ase_io import to_ase
    try:
        import ase.io
    except ImportError as exc:                                     # pragma: no cover - env specific
        raise ImportError("export_structure needs ASE: pip install 'ase>=3.29'") from exc
    path = Path(path)
    atoms = to_ase(structure)
    ase.io.write(path, atoms, format=format)
    if check:
        back = ase.io.read(path, format=format)
        if len(back) != len(atoms):
            raise ValueError(f"{path.name}: {len(back)} atoms read back, {len(atoms)} written")
        if not np.allclose(back.cell.array, atoms.cell.array, atol=1e-5):
            raise ValueError(f"{path.name}: the cell changed on write (format {format or path.suffix})")
        if not np.allclose(back.get_positions(), atoms.get_positions(), atol=1e-5):
            raise ValueError(f"{path.name}: positions changed on write (format {format or path.suffix}); "
                             "this format re-orients or wraps — choose another one")
        if list(back.get_chemical_symbols()) != list(atoms.get_chemical_symbols()):
            raise ValueError(f"{path.name}: element symbols changed on write")
    return path
