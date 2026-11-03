"""The disorder-builder model: reading it and turning sites into occupant choices.

A model is the JSON the discus-disorder-builder page exports (``schema: "discuspdf-disorder-builder"``,
``version: 1``). Orientation groups arrive already expanded to explicit rotation matrices with weights, so no
symmetry parsing happens here. The functions were moved from the builder's ``engine.py`` (builder @ 92fe94e)
without changing what they compute, so a model gives the same structure for the same seed.
"""
from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

import numpy as np

from ..atom_names import atomic_number, element_of, mass_number

SCHEMA = "discuspdf-disorder-builder"
SUPPORTED_VERSIONS = (1,)
VOID = "VOID"   # DISCUS's non-scattering atom type; pads cells so every cell has the same atom count

#: Export formats of :func:`discuspdf.builder.run` and their defaults (the page's Script step).
EXPORT_DEFAULTS = {"hdf5": True, "csv": True, "png": True, "vtk": False, "npz": False, "txt": False}


class ModelError(ValueError):
    """The model cannot be built as given."""


def load_model(source) -> dict:
    """A model dict from a dict, a JSON string, or a path to a ``.json`` file."""
    if isinstance(source, dict):
        cfg = source
    elif isinstance(source, Path) or (isinstance(source, str) and not source.lstrip().startswith("{")):
        cfg = json.loads(Path(source).read_text(encoding="utf-8"))
    elif isinstance(source, str):
        cfg = json.loads(source)
    else:
        raise TypeError("model must be a dict, a JSON string or a path")
    schema = cfg.get("schema", SCHEMA)
    if schema != SCHEMA:
        raise ModelError(f"not a disorder-builder model (schema {schema!r}, expected {SCHEMA!r})")
    version = cfg.get("version", 1)
    if version not in SUPPORTED_VERSIONS:
        raise ModelError(f"model version {version!r} is newer than this discuspdf understands "
                         f"({SUPPORTED_VERSIONS}); update discus-python")
    for key in ("cell", "supercell", "sites"):
        if key not in cfg:
            raise ModelError(f"model has no {key!r}")
    return cfg


def export_formats(cfg: dict) -> dict[str, bool]:
    """Which output files to write: the model's ``outputs`` over :data:`EXPORT_DEFAULTS`."""
    given = cfg.get("outputs") or {}
    unknown = set(given) - set(EXPORT_DEFAULTS)
    if unknown:
        raise ModelError(f"unknown export format(s) {sorted(unknown)}; known: {sorted(EXPORT_DEFAULTS)}")
    return {k: bool(given.get(k, v)) for k, v in EXPORT_DEFAULTS.items()}


def lattice_matrix(cell: dict) -> np.ndarray:
    """Rows a, b, c in the discuspdf canonical orientation (a ∥ x, b in the xy plane)."""
    from ..lattice import Lattice
    return Lattice.from_parameters(cell["a"], cell["b"], cell["c"],
                                   cell["alpha"], cell["beta"], cell["gamma"]).matrix


def _check_rotation(R, where: str) -> np.ndarray:
    R = np.asarray(R, float).reshape(3, 3)
    if not np.allclose(R @ R.T, np.eye(3), atol=1e-6):
        raise ModelError(f"{where}: matrix is not orthogonal")
    return R


def orientation_table(cfg: dict) -> dict:
    """group id -> (molecule name, list of R, weights normalised to 1)."""
    out = {}
    for gid, g in cfg.get("orientations", {}).items():
        mats = [_check_rotation(m, f"orientation group {gid!r}") for m in g["matrices"]]
        w = np.asarray(g.get("weights") or [1.0] * len(mats), float)
        if len(w) != len(mats) or np.any(w < 0) or w.sum() <= 0:
            raise ModelError(f"orientation group {gid!r}: bad weights")
        out[gid] = (g["molecule"], mats, w / w.sum())
    return out


def anchor_mass(name: str) -> float:
    """Rough mass used only to find a molecule's centre of mass (anchor ``"com"``): the mass number if the
    name carries one (D: 2.014), else 1.008 for H and 2Z (+0.5 above Z = 20). Kept as in the builder so that
    existing models place their molecules exactly as before."""
    a = mass_number(name)
    if a:
        return 2.014 if a == 2 and element_of(name) == "H" else float(a)
    z = atomic_number(name)
    return 1.008 if z == 1 else 2.0 * z + (0.5 if z > 20 else 0.0)


def molecule_local_coords(mol: dict) -> tuple[list[str], np.ndarray]:
    """Species and Cartesian coordinates (Å) relative to the molecule's anchor."""
    species = [a["el"] for a in mol["atoms"]]
    xyz = np.array([a["xyz"] for a in mol["atoms"]], float).reshape(-1, 3)
    anchor = mol.get("anchor", "centroid")
    if anchor == "centroid":
        origin = xyz.mean(axis=0)
    elif anchor == "com":
        m = np.array([anchor_mass(s) for s in species])
        origin = (xyz * m[:, None]).sum(axis=0) / m.sum()
    elif isinstance(anchor, int):
        origin = xyz[anchor]
    else:
        raise ModelError(f"molecule anchor {anchor!r} not understood")
    return species, xyz - origin


def site_rotation(site: dict, A: np.ndarray) -> np.ndarray | None:
    """Cartesian matrix of the rotation part of a site's symmetry operation (Jones symbol, e.g. "-y,x,z+1/2")."""
    sym = site.get("sym")
    if not sym:
        return None
    M = np.zeros((3, 3))
    parts = sym.replace(" ", "").lower().split(",")
    if len(parts) != 3:
        raise ModelError(f"site {site['label']}: cannot read the operation {sym!r}")
    for i, part in enumerate(parts):
        for sign, coef, var in re.findall(r"([+-]?)(\d*\.?\d*(?:/\d+)?)\*?([xyz])", part):
            c = 1.0 if coef == "" else (float(coef.split("/")[0] or 1) / float(coef.split("/")[1])
                                        if "/" in coef else float(coef))
            M[i, "xyz".index(var)] += -c if sign == "-" else c
    At = A.T
    R = At @ M @ np.linalg.inv(At)
    if not np.allclose(R @ R.T, np.eye(3), atol=1e-5):
        raise ModelError(f"site {site['label']}: the operation {sym} does not fit this cell")
    return R


def occupant_classes(site: dict, orient: dict, A: np.ndarray | None = None) -> list[dict[str, Any]]:
    """Expand a site's occupants into choice classes: atoms, single molecule orientations, vacancies."""
    classes = []
    for occ in site["occupants"]:
        p = float(occ["p"])
        if p < 0:
            raise ModelError(f"site {site['label']}: negative probability")
        if p == 0:
            continue
        # per-occupant disorder overrides (None: use the defaults of the Disorder step)
        over = {k: (None if occ.get(k) in (None, "") else float(occ[k])) for k in ("sigma", "trans", "lib")}
        if occ["type"] == "atom":
            classes.append({"type": "atom", "species": occ["species"], "p": p,
                            "b": occ.get("b"), "key": f"atom:{occ['species']}", **over})
        elif occ["type"] == "vacancy":
            classes.append({"type": "vacancy", "p": p, "key": "vacancy"})
        elif occ["type"] == "molecule":
            mol, mats, w = orient[occ["group"]]
            SR = site_rotation(site, A) if A is not None else None   # symmetry-equivalent site: rotate the molecule
            for k, (R, wk) in enumerate(zip(mats, w)):
                if wk > 0:
                    classes.append({"type": "molecule", "molecule": mol, "R": R if SR is None else SR @ R,
                                    "p": p * wk, "b": occ.get("b"), "key": f"mol:{occ['group']}:{k}",
                                    "group": occ["group"], "index": k, **over})
        else:
            raise ModelError(f"unknown occupant type {occ['type']!r}")
    total = sum(c["p"] for c in classes)
    if not classes or abs(total - 1.0) > 1e-6:
        raise ModelError(f"site {site['label']}: occupant probabilities sum to {total:.6g}, not 1")
    return classes
