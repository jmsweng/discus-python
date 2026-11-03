"""Realise a disorder-builder model as a configuration :class:`~discuspdf.Structure`.

Steps: draw an occupant (atom, molecule orientation or vacancy) for every site of every cell (random or with
the exact composition), optionally order them by Monte Carlo (:mod:`.ordering`), place the atoms with rigid-body
translation and libration of molecules and Gaussian atomic displacements, and pad every cell with ``VOID`` atoms
so that all cells hold the same number of rows (DISCUS needs that for LOTS and cell-major files). Moved from the
disorder builder's ``engine.py`` (builder @ 92fe94e) unchanged: same model and seed, same structure.
"""
from __future__ import annotations

import math

import numpy as np

from ..lattice import Lattice
from ..structure import Structure
from .model import VOID, lattice_matrix, molecule_local_coords, occupant_classes, orientation_table
from .ordering import ordering_mc

__all__ = ["draw_occupants", "random_small_rotation", "build_structure", "clash_report"]


def draw_occupants(cfg: dict, classes: list[list[dict]], n: tuple, rng) -> np.ndarray:
    """choice[s, ix, iy, iz] = index into classes[s]."""
    n_sites = len(classes)
    choice = np.zeros((n_sites, *n), dtype=np.int64)
    n_cells = int(np.prod(n))
    for s, site in enumerate(cfg["sites"]):
        p = np.array([c["p"] for c in classes[s]])
        if site.get("draw", "random") == "exact":
            counts = np.floor(p * n_cells).astype(int)
            rest = n_cells - counts.sum()
            if rest:   # largest remainders, ties broken at random
                frac = p * n_cells - counts + rng.random(len(p)) * 1e-9
                counts[np.argsort(-frac)[:rest]] += 1
            labels = np.repeat(np.arange(len(p)), counts)
            rng.shuffle(labels)
            choice[s] = labels.reshape(n[2], n[1], n[0]).transpose(2, 1, 0)
        else:
            choice[s] = rng.choice(len(p), size=n_cells, p=p).reshape(n[2], n[1], n[0]).transpose(2, 1, 0)
    return choice


def random_small_rotation(rng, sigma_deg: float) -> np.ndarray:
    if sigma_deg <= 0:
        return np.eye(3)
    axis = rng.normal(size=3)
    axis /= np.linalg.norm(axis)
    ang = math.radians(rng.normal(scale=sigma_deg))
    K = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    return np.eye(3) + math.sin(ang) * K + (1 - math.cos(ang)) * K @ K


def build_structure(cfg: dict, seed: int, log):
    rng = np.random.default_rng(seed)
    cell = cfg["cell"]
    n = tuple(int(x) for x in cfg["supercell"])
    A = lattice_matrix(cell)                       # rows a, b, c (Å)
    Ainv = np.linalg.inv(A)
    orient = orientation_table(cfg)
    mols = {name: molecule_local_coords(m) for name, m in cfg.get("molecules", {}).items()}
    classes = [occupant_classes(site, orient, A) for site in cfg["sites"]]

    # slots per site = largest number of atoms any occupant needs
    def n_atoms_of(c):
        return len(mols[c["molecule"]][0]) if c["type"] == "molecule" else (1 if c["type"] == "atom" else 0)
    slots = [max(1, max(n_atoms_of(c) for c in cl)) for cl in classes]
    # where padding VOIDs go: the atoms of the site's largest occupant, so slots keep their usual place
    pad_local = []
    for cl in classes:
        big = max(cl, key=n_atoms_of)
        if big["type"] == "molecule":
            pad_local.append(mols[big["molecule"]][1] @ big["R"].T)
        else:
            pad_local.append(np.zeros((1, 3)))
    m = int(sum(slots))

    choice = draw_occupants(cfg, classes, n, rng)
    corr_report = ordering_mc(cfg, classes, choice, A, mols, rng, log)

    dis = cfg.get("disorder", {})
    rigid = dis.get("molecule_rigid", {}) or {}
    t_sigma = float(rigid.get("translation_sigma", 0.0))
    lib_sigma = float(rigid.get("libration_sigma_deg", 0.0))
    use_b = bool(cfg.get("calc", {}).get("debye_waller", False))

    N = m * int(np.prod(n))
    species = np.empty(N, dtype=object)
    cart = np.zeros((N, 3))
    cell_index = np.zeros((N, 3), dtype=np.int64)
    site_label = np.empty(N, dtype=object)
    b_iso = np.ones(N)
    molecule_id = np.full(N, -1, dtype=np.int64)    # which molecule an atom belongs to (-1: none)
    is_void = np.zeros(N, dtype=bool)
    sigma_row = np.full(N, np.nan)                  # per-occupant displacement σ (NaN: species default)

    row = 0
    mol_counter = 0
    stats = {}
    for iz in range(n[2]):                      # DISCUS cell-major order: x fastest
        for iy in range(n[1]):
            for ix in range(n[0]):
                origin = np.array([ix, iy, iz], float)
                for s, site in enumerate(cfg["sites"]):
                    c = classes[s][choice[s, ix, iy, iz]]
                    stats.setdefault(site["label"], {}).setdefault(c["key"], 0)
                    stats[site["label"]][c["key"]] += 1
                    centre = (origin + np.asarray(site["frac"], float)) @ A
                    b_site = float(site.get("b_iso", 1.0))
                    if c["type"] == "atom":
                        sp, pos = [c["species"]], centre[None, :]
                        bs = [float(c["b"]) if c.get("b") is not None else b_site]
                    elif c["type"] == "molecule":
                        msp, local = mols[c["molecule"]]
                        lib = lib_sigma if c.get("lib") is None else c["lib"]
                        trans = t_sigma if c.get("trans") is None else c["trans"]
                        R = c["R"] @ random_small_rotation(rng, lib)
                        shift = rng.normal(scale=trans, size=3) if trans > 0 else 0.0
                        sp, pos = msp, centre + shift + local @ R.T
                        bs = [float(c["b"]) if c.get("b") is not None else b_site] * len(msp)
                    else:
                        sp, pos, bs = [], np.zeros((0, 3)), []
                    k = len(sp)
                    for j in range(slots[s]):
                        r = row + j
                        cell_index[r] = (ix, iy, iz)
                        site_label[r] = f"{site['label'][:10]}.{j}"
                        if j < k:
                            species[r] = sp[j]
                            cart[r] = pos[j]
                            b_iso[r] = bs[j]
                            if c.get("sigma") is not None:
                                sigma_row[r] = c["sigma"]
                            if c["type"] == "molecule":
                                molecule_id[r] = mol_counter
                        else:
                            species[r] = VOID
                            cart[r] = centre + pad_local[s][j]
                            is_void[r] = True
                    if c["type"] == "molecule":
                        mol_counter += 1
                    row += slots[s]

    # static Gaussian displacements: per-occupant σ where given, else the species default
    sig = {k: float(v) for k, v in (dis.get("atomic_displacements", {}) or {}).items()}
    sigma = np.array([sig.get(s, 0.0) for s in species])
    sigma = np.where(np.isfinite(sigma_row), sigma_row, sigma)
    sigma[is_void] = 0.0
    if np.any(sigma < 0):
        raise ValueError("displacement σ must be >= 0")
    move = sigma > 0
    cart[move] += rng.normal(size=(int(move.sum()), 3)) * sigma[move, None]

    frac_unit = cart @ Ainv                           # unit-cell fractions
    if cfg.get("wrap", False):   # off by default: DISCUS lists atoms by cell, molecules stay whole
        frac_unit = np.mod(frac_unit, np.array(n, float))
    frac_super = frac_unit / np.array(n, float)
    if not use_b:
        b_iso[:] = 1.0                                 # written but not used (temp ignore)
    if np.any(b_iso[~is_void] < 0.001):
        raise ValueError("DISCUS 6.23.01 replaces B < 0.001 by 1.0; use B >= 0.001 or switch Debye–Waller off")

    super_lat = Lattice.from_parameters(cell["a"] * n[0], cell["b"] * n[1], cell["c"] * n[2],
                                        cell["alpha"], cell["beta"], cell["gamma"])
    structure = Structure(super_lat, species, frac_super, model="configuration", site_label=site_label,
                          cell_index=cell_index, supercell=n, b_iso=b_iso, min_separation=None)
    info = {"atoms_per_cell": m, "slots_per_site": slots, "n_rows": N, "n_void": int(is_void.sum()),
            "n_molecules": mol_counter, "occupant_counts": stats, "correlations": corr_report,
            "molecule_id": molecule_id, "is_void": is_void, "unit_lattice": A}
    return structure, info


def clash_report(structure, info, threshold: float, log) -> dict:
    """Shortest distances between atoms that are not part of the same molecule (periodic)."""
    real = ~info["is_void"]
    idx = np.where(real)[0]
    if len(idx) < 2:
        return {"n_close_pairs": 0}
    A = structure.lattice.matrix
    f = structure.frac[idx] % 1.0
    mol = info["molecule_id"][idx]
    # cell list in fractional space
    lengths = np.linalg.norm(A, axis=1)
    nb = np.maximum(1, np.floor(lengths / max(threshold, 0.5)).astype(int))
    nb = np.minimum(nb, 60)
    b = np.floor(f * nb).astype(int) % nb
    key = (b[:, 0] * nb[1] + b[:, 1]) * nb[2] + b[:, 2]
    order = np.argsort(key)
    start = np.searchsorted(key[order], np.arange(nb.prod()))
    end = np.searchsorted(key[order], np.arange(nb.prod()), side="right")
    close, dmin = [], np.inf
    offsets = [(i, j, k) for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)]
    for p in range(len(idx)):
        for off in offsets:
            bb = (b[p] + off) % nb
            kk = (bb[0] * nb[1] + bb[1]) * nb[2] + bb[2]
            cand = order[start[kk]:end[kk]]
            cand = cand[cand > p]
            if not len(cand):
                continue
            df = f[cand] - f[p]
            df -= np.round(df)
            d = np.linalg.norm(df @ A, axis=1)
            other = (mol[cand] != mol[p]) | (mol[p] < 0)
            d = d[other]
            if len(d):
                dmin = min(dmin, float(d.min()))
                for q in cand[other][d < threshold]:
                    close.append((int(idx[p]), int(idx[q])))
    close = sorted(set(close))
    if close:
        log(f"  WARNING: {len(close)} atom pairs from different molecules/sites closer than {threshold} Å "
            f"(shortest {dmin:.3f} Å). Check molecule sizes, orientations and site positions.")
    else:
        log(f"  no contacts between different molecules/sites below {threshold} Å")
    return {"threshold": threshold, "n_close_pairs": len(close), "shortest": dmin, "examples": close[:20]}
