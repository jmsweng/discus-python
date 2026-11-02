"""Compare DiscusVolume geometry with Lattice for random and special cells."""

import numpy as np
import sys

from discuspdf.discus_output import _validated_direct_basis
from discuspdf.lattice import Lattice


def compare_one(lengths, angles, label):
    """Compare direct basis, reciprocal basis (both two_pi values), and cell volume."""
    results = {}
    try:
        A_dv = _validated_direct_basis(lengths, angles)
    except Exception as e:
        results["direct_basis_dv"] = f"EXCEPTION: {type(e).__name__}: {e}"
        return results

    try:
        lat = Lattice.from_parameters(lengths[0], lengths[1], lengths[2],
                                      angles[0], angles[1], angles[2])
        A_lat = lat.matrix
    except Exception as e:
        results["direct_basis_lat"] = f"EXCEPTION: {type(e).__name__}: {e}"
        return results

    diff_direct = np.abs(A_dv - A_lat)
    results["direct_basis_max_abs"] = float(np.max(diff_direct))
    denom = np.abs(A_lat)
    denom[denom == 0] = 1.0  # avoid division by zero
    results["direct_basis_max_rel"] = float(np.max(diff_direct / denom))
    results["direct_basis_exact_zero"] = bool(np.allclose(A_dv, A_lat, atol=0, rtol=0))

    # Reciprocal basis: DiscusVolume style
    B_dv_2pi = np.linalg.inv(A_dv).T
    B_dv_no2pi = B_dv_2pi / (2 * np.pi)

    # Lattice reciprocal
    B_lat = lat.reciprocal_matrix  # no 2*pi
    B_lat_2pi = B_lat * (2 * np.pi)

    diff_recip_2pi = np.abs(B_dv_2pi - B_lat_2pi)
    results["reciprocal_2pi_max_abs"] = float(np.max(diff_recip_2pi))
    denom2 = np.abs(B_lat_2pi)
    denom2[denom2 == 0] = 1.0
    results["reciprocal_2pi_max_rel"] = float(np.max(diff_recip_2pi / denom2))
    results["reciprocal_2pi_exact_zero"] = bool(np.allclose(B_dv_2pi, B_lat_2pi, atol=0, rtol=0))

    diff_recip_no2pi = np.abs(B_dv_no2pi - B_lat)
    results["reciprocal_no2pi_max_abs"] = float(np.max(diff_recip_no2pi))
    denom3 = np.abs(B_lat)
    denom3[denom3 == 0] = 1.0
    results["reciprocal_no2pi_max_rel"] = float(np.max(diff_recip_no2pi / denom3))
    results["reciprocal_no2pi_exact_zero"] = bool(np.allclose(B_dv_no2pi, B_lat, atol=0, rtol=0))

    # Cell volume
    vol_dv = float(abs(np.linalg.det(A_dv)))
    vol_lat = float(lat.volume)
    results["cell_volume_dv"] = vol_dv
    results["cell_volume_lat"] = vol_lat
    results["cell_volume_max_abs"] = abs(vol_dv - vol_lat)
    results["cell_volume_max_rel"] = abs(vol_dv - vol_lat) / vol_lat if vol_lat != 0 else float('inf')
    results["cell_volume_exact_zero"] = (vol_dv == vol_lat)

    return results


def test_invalid_inputs():
    """Test that both reject the same invalid inputs."""
    cases = [
        ("angles_sum_360", [5.0, 5.0, 5.0], [90.0, 90.0, 180.0]),
        ("degenerate_cell", [5.0, 5.0, 5.0], [90.0, 90.0, 0.0]),
        ("zero_length", [0.0, 5.0, 5.0], [90.0, 90.0, 90.0]),
        ("negative_length", [-5.0, 5.0, 5.0], [90.0, 90.0, 90.0]),
    ]
    results = {}
    for name, lengths, angles in cases:
        dv_exc = None
        lat_exc = None
        try:
            _validated_direct_basis(lengths, angles)
        except Exception as e:
            dv_exc = f"{type(e).__name__}: {e}"
        try:
            Lattice.from_parameters(*lengths, *angles)
        except Exception as e:
            lat_exc = f"{type(e).__name__}: {e}"
        results[name] = {"discus_output": dv_exc, "lattice": lat_exc}
    return results


def main():
    rng = np.random.default_rng(0)
    all_results = []
    invalid_results = test_invalid_inputs()

    # 200 random cells
    for i in range(200):
        lengths = rng.uniform(2.0, 20.0, size=3)
        angles = rng.uniform(60.0, 120.0, size=3)
        # Filter to valid, positive-volume cells
        a, b, c = lengths
        alpha, beta, gamma = np.radians(angles)
        cos_a, cos_b, cos_c = np.cos(alpha), np.cos(beta), np.cos(gamma)
        det = 1.0 - cos_a**2 - cos_b**2 - cos_c**2 + 2.0 * cos_a * cos_b * cos_c
        if det <= 0:
            continue
        result = compare_one(lengths, angles, f"random_{i}")
        all_results.append(result)

    # C1 cubic cell (from the archived evidence)
    c1_lengths = [5.6537, 5.6537, 5.6537]
    c1_angles = [90.0, 90.0, 90.0]
    result = compare_one(c1_lengths, c1_angles, "c1_cubic")
    all_results.append(result)

    # Gate 2 4 Å cubic cell
    gate2_lengths = [4.0, 4.0, 4.0]
    gate2_angles = [90.0, 90.0, 90.0]
    result = compare_one(gate2_lengths, gate2_angles, "gate2_cubic")
    all_results.append(result)

    # Triclinic fixture
    try:
        from discuspdf.cif_io import read_cif
        cif_path = "tests/fixtures/interop/triclinic_p1bar.cif"
        if __import__("os").path.exists(cif_path):
            model = read_cif(cif_path)
            lengths = [model.cell.a, model.cell.b, model.cell.c]
            angles = [model.cell.alpha, model.cell.beta, model.cell.gamma]
            result = compare_one(lengths, angles, "triclinic_p1bar")
            all_results.append(result)
        else:
            all_results.append({"triclinic_p1bar": "FILE NOT FOUND"})
    except Exception as e:
        all_results.append({"triclinic_p1bar": f"EXCEPTION: {type(e).__name__}: {e}"})

    # Print summary
    print(f"{'Case':<25} {'Direct max_abs':>14} {'Direct max_rel':>14} {'Recip2pi max_abs':>16} {'Recip2pi max_rel':>16} {'Vol max_abs':>12} {'Vol max_rel':>12} {'Exact zero':>10}")
    print("-" * 140)
    for r in all_results:
        if "direct_basis_max_abs" in r:
            print(f"{r.get('label', '?'):<25} {r['direct_basis_max_abs']:>14.2e} {r['direct_basis_max_rel']:>14.2e} {r['reciprocal_2pi_max_abs']:>16.2e} {r['reciprocal_2pi_max_rel']:>16.2e} {r['cell_volume_max_abs']:>12.2e} {r['cell_volume_max_rel']:>12.2e} {str(r['direct_basis_exact_zero']):>10}")
        else:
            print(f"{str(r):<100}")

    print("\n=== Invalid input behaviour ===")
    for name, r in invalid_results.items():
        print(f"{name}:")
        print(f"  discus_output: {r['discus_output']}")
        print(f"  lattice:       {r['lattice']}")


if __name__ == "__main__":
    main()
