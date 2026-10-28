"""Standalone scientific reference for frozen-C1 DISCUS 6.23.01 comparisons.

Reference predictions are IT92-based.  |F-F_average|² is a candidate
observable, NOT a certified implementation of the engine's set aver.
No production scattering code is imported and absolute scale is unresolved.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import numpy as np

IT92 = {
    "GA": ((15.2354, 6.7006, 4.3591, 2.9623),
           (3.0669, 0.2412, 10.7805, 61.4135), 1.7189),
    "AS": ((16.6723, 6.0701, 3.4313, 4.2779),
           (2.6345, 0.2647, 12.9479, 47.7972), 2.531),
}
C1_SHA = "ceefca57e26ab40e03d2adfd443cff559bd7fc3c837402f418f6e8def35406b9"


def xray_factor(species, s):
    """IT92 Gaussian form factor; s=sin(theta)/lambda [1/Å]."""
    a, b, c = IT92[species.upper()]
    s2 = np.asarray(s, dtype=float) ** 2
    return c + sum(ai * np.exp(-bi * s2) for ai, bi in zip(a, b))


def hkl_grid(shape, corner, increments):
    """Voxel (i,j,k) maps to corner + increments @ (i,j,k)."""
    if len(shape) != 3 or any(
        not isinstance(n, (int, np.integer)) or n < 1 for n in shape
    ):
        raise ValueError("shape must contain three positive integers")
    corner = np.asarray(corner, dtype=float)
    increments = np.asarray(increments, dtype=float)
    if corner.shape != (3,) or increments.shape != (3, 3):
        raise ValueError("expected corner (3,) and increments (3,3)")
    if not np.isfinite(corner).all() or not np.isfinite(increments).all():
        raise ValueError("grid coordinates must be finite")
    if np.linalg.matrix_rank(increments) != 3:
        raise ValueError("increment matrix must be nonsingular")
    ijk = np.indices(shape, dtype=float)
    return np.einsum("ab,bijk->ijka", increments, ijk) + corner


def load_c1_structure(path):
    """Load only the verified frozen C1 structure, in fractional units."""
    species, positions, biso = [], [], []
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if not line or line.lower().startswith(
            ("title", "spcgr", "cell", "format", "atoms")
        ):
            continue
        parts = line.replace(",", " ").split()
        if not parts or parts[0].upper() not in ("GA", "AS"):
            continue
        if len(parts) < 5:
            raise ValueError("malformed atom line")
        species.append(parts[0].upper())
        positions.append([float(v) for v in parts[1:4]])
        biso.append(float(parts[4]))
    species = np.asarray(species)
    positions = np.asarray(positions, dtype=float)
    biso = np.asarray(biso, dtype=float)
    if (len(species) != 4096 or (sum(species == "GA"), sum(species == "AS"))
            != (1723, 2373)):
        raise ValueError("unexpected C1 frozen species counts")
    if positions.shape != (4096, 3) or not np.isfinite(positions).all():
        raise ValueError("invalid C1 positions")
    if not np.allclose(biso, 1.0):
        raise ValueError("unexpected Biso")
    return species, positions, biso


def empirical_basis_probabilities(species, positions):
    """Sample-fitted Ga fractions: NOT independently prescribed ensemble."""
    groups = {}
    for j, b in enumerate(np.round(positions % 1.0, 4)):
        groups.setdefault(tuple(b), []).append(j)
    if len(groups) != 8 or any(len(v) != 512 for v in groups.values()):
        raise ValueError("expected eight 512-site basis families")
    return {b: float(np.mean(species[js] == "GA")) for b, js in groups.items()}


def direct_average(hkl, positions, p_ga_per_site, *, a, biso, chunk=64):
    """Independent EXPECTED amplitude: per-site prescribed Ga probabilities.

    No actual species list is an input; no unjustified division by atom count.
    """
    hkl = np.asarray(hkl, float).reshape(-1, 3)
    positions = np.asarray(positions, float)
    p = np.asarray(p_ga_per_site, float)
    biso = np.asarray(biso, float)
    if (positions.ndim != 2 or positions.shape[1] != 3
            or p.shape != (len(positions),) or biso.shape != p.shape):
        raise ValueError("incompatible site arrays")
    if not np.isfinite(p).all() or np.any((p < 0) | (p > 1)):
        raise ValueError("site probabilities out of range")
    result = np.empty(len(hkl), complex)
    for k in range(0, len(hkl), chunk):
        h = hkl[k:k + chunk]
        s = np.linalg.norm(h, axis=1) / (2 * a)
        f_expected = (xray_factor("GA", s)[:, None] * p
                      + xray_factor("AS", s)[:, None] * (1 - p))
        result[k:k + chunk] = np.sum(
            f_expected * np.exp(-s[:, None]**2 * biso)
            * np.exp(2j * np.pi * (h @ positions.T)), axis=1)
    return result


def factored_average(hkl, basis_p_ga, *, a, n_min=-4, n_max=3, biso=1.0):
    """Separate finite-box calculation to cross-check direct_average."""
    hkl = np.asarray(hkl, float).reshape(-1, 3)
    s = np.linalg.norm(hkl, axis=1) / (2 * a)
    f_ga, f_as = xray_factor("GA", s), xray_factor("AS", s)
    basis = np.zeros(len(hkl), complex)
    for site, p in basis_p_ga.items():
        basis += (p * f_ga + (1 - p) * f_as) * np.exp(
            2j * np.pi * (hkl @ np.asarray(site)))
    factors = [
        np.exp(2j * np.pi * hkl[:, d, None]
               * np.arange(n_min, n_max + 1)[None, :]).sum(axis=1)
        for d in range(3)
    ]
    return np.exp(-biso * s**2) * basis * factors[0] * factors[1] * factors[2]


def actual_amplitude(hkl, positions, species, *, a, biso, chunk=64):
    """Direct summation of the actual species on the frozen site list."""
    hkl = np.asarray(hkl, float).reshape(-1, 3)
    ga = species == "GA"
    result = np.empty(len(hkl), complex)
    for k in range(0, len(hkl), chunk):
        h = hkl[k:k + chunk]
        s = np.linalg.norm(h, axis=1) / (2 * a)
        value = np.zeros(len(h), complex)
        for mask, name in ((ga, "GA"), (~ga, "AS")):
            value += xray_factor(name, s) * np.sum(
                np.exp(-s[:, None]**2 * biso[mask])
                * np.exp(2j * np.pi * (h @ positions[mask].T)),
                axis=1)
        result[k:k + chunk] = value
    return result


def metrics(predicted, observed, grid):
    """Never fit scale or offset; mask near-zero voxels in relative ratios."""
    d = np.asarray(predicted) - np.asarray(observed)
    norm = np.linalg.norm(observed.ravel())
    floor = 1e-6 * np.max(np.abs(observed))
    mask = np.abs(observed) > floor
    worst = np.unravel_index(np.argmax(np.abs(d)), d.shape)
    return {
        "n_voxels": int(d.size), "above_floor": int(mask.sum()),
        "ratio_floor": float(floor),
        "relative_l2": float(np.linalg.norm(d.ravel()) / norm) if norm else None,
        "mean_absolute_error": float(np.mean(np.abs(d))),
        "max_absolute_error": float(np.max(np.abs(d))),
        "pearson": float(np.corrcoef(predicted.ravel(), observed.ravel())[0, 1]),
        "median_relative_error_above_floor":
            float(np.median(np.abs(d[mask] / observed[mask]))) if mask.any() else None,
        "worst_index": list(map(int, worst)),
        "worst_hkl": grid[worst].tolist(),
        "worst_observed": float(observed[worst]),
        "worst_predicted": float(predicted[worst]),
    }


def benchmark(archive, output):
    """Recompute full 33^3 C1 benchmark from *actual* archived HDF5 voxels."""
    import h5py
    archive, output = Path(archive), Path(output)
    structure = archive / "C1_frozen.stru"
    if hashlib.sha256(structure.read_bytes()).hexdigest() != C1_SHA:
        raise ValueError("frozen structure SHA-256 does not match")
    species, positions, biso = load_c1_structure(structure)
    observed, hashes, grid = {}, {}, None
    for key, filename, label in (
        ("total", "C1_offcomm_total.h5", b"bragg_present"),
        ("diffuse", "C1_offcomm_diffuse.h5", b"bragg_subtracted"),
    ):
        path = archive / filename
        hashes[filename] = hashlib.sha256(path.read_bytes()).hexdigest()
        with h5py.File(path) as handle:
            data = handle["entry/data"]
            if data["data_type_with_bragg"][()] != label:
                raise ValueError("unexpected Bragg status")
            v = data["data_values"][()]
            if v.shape != (33, 33, 33) or tuple(data["data_dimension"][()]) != v.shape:
                raise ValueError("unexpected HDF5 volume")
            new_grid = hkl_grid(v.shape, data["data_corner"][()],
                                data["data_increment_vector"][()])
            if grid is not None and not np.array_equal(grid, new_grid):
                raise ValueError("inconsistent grids")
            grid = new_grid
            observed[key] = v
    coords = grid.reshape(-1, 3)
    actual = actual_amplitude(coords, positions, species, a=5.6537, biso=biso)
    basis = empirical_basis_probabilities(species, positions)
    average = factored_average(coords, basis, a=5.6537)
    actual = actual.reshape((33, 33, 33))
    average = average.reshape((33, 33, 33))
    total = np.abs(actual)**2
    diff = np.abs(actual - average)**2
    invalid_candidate = np.abs(actual)**2 - np.abs(average)**2
    comparisons = {
        "total_IT92": metrics(total, observed["total"], grid),
        "diffuse_empirical_IT92": metrics(diff, observed["diffuse"], grid),
        "difference_of_intensities_IT92":
            metrics(invalid_candidate, observed["diffuse"], grid),
    }
    report = {
        "source": "archived DISCUS 6.23.01; no new engine execution",
        "reference_factors": "IT92 (not engine default Waasmaier-Kirfel)",
        "average_model": "basis probabilities fitted to same C1, conditional only",
        "normalization": "NOT calibrated",
        "input_hashes": {"C1_frozen.stru": C1_SHA, **hashes},
        "volume": [33, 33, 33],
        "basis_p_ga": {str(k): v for k, v in basis.items()},
        "comparisons": comparisons,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + "\n")
    np.savez_compressed(output.with_suffix(".npz"),
                        hkl=grid, observed_total=observed["total"],
                        observed_diffuse=observed["diffuse"],
                        F_actual_IT92=actual, F_average_empirical_IT92=average,
                        signed_total_residual=total-observed["total"],
                        signed_diffuse_residual=diff-observed["diffuse"])
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("archive_directory", type=Path)
    parser.add_argument("output_json", type=Path)
    args = parser.parse_args()
    benchmark(args.archive_directory, args.output_json)
