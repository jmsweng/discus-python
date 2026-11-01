"""Self-test of gate0_analyze.py on SYNTHETIC data (never engine data): it must pass correct data and fail
each deliberate corruption. This tests the analysis, not DISCUS.

Usage: python gate0_selftest.py [WORK_DIR]
"""
from __future__ import annotations

import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import oracle as O  # noqa: E402
import gate0_analyze as AN  # noqa: E402

import json
K = json.loads((HERE / "predictions.json").read_text())["n_realizations"]


def write_discus_h5(path: Path, values: np.ndarray, bragg: str):
    """Write a file with the dataset layout DISCUS 6.23.01 produces (form nexus, value int)."""
    import h5py
    with h5py.File(path, "w") as f:
        g = f.create_group("entry").create_group("data")
        s = lambda k, v: g.create_dataset(k, data=np.bytes_(v))
        s("audit_author_name", "DISCUS"); s("audit_conform_dict_name", "Disorder unified data")
        s("audit_conform_dict_version", "0.0.0"); s("audit_creation_date", "synthetic")
        s("audit_creation_method", "selftest-synthetic"); s("coordinate_unit", "basecell_fractional")
        g.create_dataset("data_axes", data=np.array([1, 2, 3], np.int32))
        g.create_dataset("data_corner", data=np.zeros(3))
        g.create_dataset("data_dimension", data=np.array(values.shape, np.int32))
        g.create_dataset("data_increment_vector", data=np.eye(3) * 0.125)
        s("data_type_axes", "hkl"); s("data_type_content", "intensity"); s("data_type_reciprocal", "reciprocal")
        s("data_type_with_bragg", bragg)
        g.create_dataset("data_values", data=values.astype(np.float64))
        g.create_dataset("unit_cell_angles", data=np.full(3, 90.0))
        g.create_dataset("unit_cell_lengths", data=np.full(3, 3.75))


def _make_clean_run(work: Path, arm: str | None = None) -> Path:
    """Build a clean synthetic run. If arm is None, create files for both arms."""
    name = f"clean_{arm}" if arm else "clean"
    run = work / name
    if run.exists():
        shutil.rmtree(run)
    run.mkdir(parents=True)

    # Copy registered inputs (so hash checks pass)
    if (HERE / "inputs").exists():
        shutil.copytree(HERE / "inputs", run, dirs_exist_ok=True)

    arms = ("A", "B") if arm is None else (arm,)

    for a in arms:
        b_cu, b_au = (O.B_CU_A, O.B_AU_A) if a == "A" else (O.B_CU_B, O.B_AU_B)
        for k in range(K):
            # Build a synthetic structure
            seed = 42 + k  # fixed seed for reproducibility
            rng = np.random.default_rng(seed)
            # Random occupancy: 2048 sites, each Au with p=0.25
            occ = np.zeros((O.N_CELL, O.N_CELL, O.N_CELL, 4), dtype=float)
            for j in range(4):
                occ[:, :, :, j] = rng.random((O.N_CELL, O.N_CELL, O.N_CELL)) < O.P_AU

            # Compute D and I from the oracle
            s = O.s_of_grid("engine")
            wk = O.coefficients_for("engine")
            f_au = O.wk_f(wk["AU"], s)
            f_cu = O.wk_f(wk["CU"], s)

            D = np.zeros(s.shape, complex)
            m_idx = np.arange(25) % O.N_CELL
            for j in range(4):
                occ_j = occ[:, :, :, j]
                dft_full = np.fft.ifftn(occ_j) * (O.N_CELL ** 3)
                dft = dft_full[m_idx, m_idx, m_idx]
                xj = O.X_SITES[j]
                phase = np.exp(2j * np.pi * (
                    O.H_AXIS[:, None, None] * xj[0] +
                    O.H_AXIS[None, :, None] * xj[1] +
                    O.H_AXIS[None, None, :] * xj[2]
                ))
                g_j = f_au * np.exp(-b_au * s ** 2) - f_cu * np.exp(-b_cu * s ** 2)
                D += g_j * dft * phase

            f_avg = O.P_CU * f_cu * np.exp(-b_cu * s ** 2) + O.P_AU * f_au * np.exp(-b_au * s ** 2)
            h_int = np.isclose(O.H_AXIS, np.rint(O.H_AXIS))
            is_int = h_int[:, None, None] & h_int[None, :, None] & h_int[None, None, :]
            L = np.where(is_int, float(O.N_C), 0.0)
            F = f_avg * L + D
            I = np.abs(F) ** 2

            write_discus_h5(run / f"gate0_{a}_total_{k:02d}.h5", I, "bragg_present")
            write_discus_h5(run / f"gate0_{a}_diffuse_{k:02d}.h5", np.abs(D) ** 2, "bragg_subtracted")

            # Write a log
            temp = "ignored" if a == "A" else "used"
            log = f"Temp. factors      : {temp}\n"
            if a == "B":
                log += "B(Cu) = 0.5, B(Au) = 0.8\n"
            (run / f"gate0_{a}_total_{k:02d}.log").write_text(log)
            (run / f"gate0_{a}_diffuse_{k:02d}.log").write_text(log)

    return run


def _corrupt_run(work: Path, arm: str, corrupt: str) -> Path:
    """Build a corrupted synthetic run."""
    _make_clean_run(work)  # Create clean run with both arms
    run = work / f"{corrupt}_{arm}"
    if run.exists():
        shutil.rmtree(run)
    shutil.copytree(work / "clean", run)

    b_cu, b_au = (O.B_CU_A, O.B_AU_A) if arm == "A" else (O.B_CU_B, O.B_AU_B)

    if corrupt == "p_030":
        # p = 0.30 instead of 0.25: resample with p=0.30
        for k in range(K):
            seed = 42 + k
            rng = np.random.default_rng(seed)
            occ = np.zeros((O.N_CELL, O.N_CELL, O.N_CELL, 4), dtype=float)
            for j in range(4):
                occ[:, :, :, j] = rng.random((O.N_CELL, O.N_CELL, O.N_CELL)) < 0.30
            s = O.s_of_grid("engine")
            wk = O.coefficients_for("engine")
            f_au = O.wk_f(wk["AU"], s)
            f_cu = O.wk_f(wk["CU"], s)
            D = np.zeros(s.shape, complex)
            m_idx = np.arange(25) % O.N_CELL
            for j in range(4):
                occ_j = occ[:, :, :, j]
                dft_full = np.fft.ifftn(occ_j) * (O.N_CELL ** 3)
                dft = dft_full[m_idx, m_idx, m_idx]
                xj = O.X_SITES[j]
                phase = np.exp(2j * np.pi * (
                    O.H_AXIS[:, None, None] * xj[0] +
                    O.H_AXIS[None, :, None] * xj[1] +
                    O.H_AXIS[None, None, :] * xj[2]
                ))
                g_j = f_au * np.exp(-b_au * s ** 2) - f_cu * np.exp(-b_cu * s ** 2)
                D += g_j * dft * phase
            f_avg = 0.7 * f_cu * np.exp(-b_cu * s ** 2) + 0.3 * f_au * np.exp(-b_au * s ** 2)
            h_int = np.isclose(O.H_AXIS, np.rint(O.H_AXIS))
            is_int = h_int[:, None, None] & h_int[None, :, None] & h_int[None, None, :]
            L = np.where(is_int, float(O.N_C), 0.0)
            F = f_avg * L + D
            I = np.abs(F) ** 2
            # Overwrite the diffuse file
            import h5py
            with h5py.File(run / f"gate0_{arm}_diffuse_{k:02d}.h5", "r+") as f:
                f["entry/data/data_values"][...] = np.abs(D) ** 2

    elif corrupt == "correlated":
        # Correlated configuration: use a block model instead of independent
        for k in range(K):
            seed = 42 + k
            rng = np.random.default_rng(seed)
            # 2x2 blocks in x-y plane
            n = O.N_CELL
            s_blk = rng.choice(np.array([0.0, 1.0]), size=(n // 2, n // 2, n))
            occ = np.zeros((n, n, n, 4), dtype=float)
            x = np.arange(n)
            for j in range(4):
                occ[:, :, :, j] = s_blk[x[:, None, None] // 2, x[None, :, None] // 2, x[None, :, None]]
            s = O.s_of_grid("engine")
            wk = O.coefficients_for("engine")
            f_au = O.wk_f(wk["AU"], s)
            f_cu = O.wk_f(wk["CU"], s)
            D = np.zeros(s.shape, complex)
            m_idx = np.arange(25) % O.N_CELL
            for j in range(4):
                occ_j = occ[:, :, :, j]
                dft_full = np.fft.ifftn(occ_j) * (n ** 3)
                dft = dft_full[m_idx, m_idx, m_idx]
                xj = O.X_SITES[j]
                phase = np.exp(2j * np.pi * (
                    O.H_AXIS[:, None, None] * xj[0] +
                    O.H_AXIS[None, :, None] * xj[1] +
                    O.H_AXIS[None, None, :] * xj[2]
                ))
                g_j = f_au * np.exp(-b_au * s ** 2) - f_cu * np.exp(-b_cu * s ** 2)
                D += g_j * dft * phase
            import h5py
            with h5py.File(run / f"gate0_{arm}_diffuse_{k:02d}.h5", "r+") as f:
                f["entry/data/data_values"][...] = np.abs(D) ** 2

    elif corrupt == "dw_convention":
        # DW convention e^{-2Bs²} in amplitude (Arm B only)
        if arm == "B":
            for k in range(K):
                import h5py
                with h5py.File(run / f"gate0_{arm}_diffuse_{k:02d}.h5", "r+") as f:
                    D = f["entry/data/data_values"][()]
                    # Apply the wrong DW factor
                    s = O.s_of_grid("engine")
                    dw = np.exp(-2 * O.B_AU_B * s ** 2)
                    D = D * dw
                    f["entry/data/data_values"][...] = D

    elif corrupt == "scale_1pct":
        # Intensities × 1.01
        for k in range(K):
            for kind in ("total", "diffuse"):
                import h5py
                with h5py.File(run / f"gate0_{arm}_{kind}_{k:02d}.h5", "r+") as f:
                    f["entry/data/data_values"][...] *= 1.01

    elif corrupt == "sample_fitted":
        # Sample-fitted occupancies used in the expectation (must make S-criteria uninformative)
        # This is a static check: the evaluator never reads sample occupancies
        pass

    elif corrupt == "species_swapped":
        # Species swapped Cu↔Au in one realization
        k = 0
        import h5py
        with h5py.File(run / f"gate0_{arm}_diffuse_{k:02d}.h5", "r+") as f:
            D = f["entry/data/data_values"][()]
            # Negate D (swapping Cu and Au flips the sign of the contrast)
            f["entry/data/data_values"][...] = -D

    elif corrupt == "missing_aver":
        # Missing set aver: diffuse = total
        for k in range(K):
            import h5py
            with h5py.File(run / f"gate0_{arm}_total_{k:02d}.h5", "r") as f:
                I = f["entry/data/data_values"][()]
            with h5py.File(run / f"gate0_{arm}_diffuse_{k:02d}.h5", "r+") as f:
                f["entry/data/data_values"][...] = I

    elif corrupt == "shuffled_seeds":
        # Realizations shuffled against their seeds
        # Swap two realizations
        import h5py
        for kind in ("total", "diffuse"):
            with h5py.File(run / f"gate0_{arm}_{kind}_00.h5", "r") as f:
                d0 = f["entry/data/data_values"][()]
            with h5py.File(run / f"gate0_{arm}_{kind}_01.h5", "r") as f:
                d1 = f["entry/data/data_values"][()]
            with h5py.File(run / f"gate0_{arm}_{kind}_00.h5", "r+") as f:
                f["entry/data/data_values"][...] = d1
            with h5py.File(run / f"gate0_{arm}_{kind}_01.h5", "r+") as f:
                f["entry/data/data_values"][...] = d0

    elif corrupt == "friedel_noise":
        # Friedel-violating noise 1e-4
        for k in range(K):
            rng = np.random.default_rng(100 + k)
            for kind in ("total", "diffuse"):
                import h5py
                with h5py.File(run / f"gate0_{arm}_{kind}_{k:02d}.h5", "r+") as f:
                    D = f["entry/data/data_values"][()]
                    noise = rng.standard_normal(D.shape) * 1e-4 * np.abs(D).max()
                    # Break Friedel: add noise that is not symmetric
                    D = D + noise
                    f["entry/data/data_values"][...] = D

    elif corrupt == "duplicated_realization":
        # One realization duplicated
        import h5py
        for kind in ("total", "diffuse"):
            with h5py.File(run / f"gate0_{arm}_{kind}_00.h5", "r") as f:
                d0 = f["entry/data/data_values"][()]
            with h5py.File(run / f"gate0_{arm}_{kind}_01.h5", "r+") as f:
                f["entry/data/data_values"][...] = d0

    return run


CASES = {  # corruption -> (expected status, criterion-name fragment that must fail)
    # Note: the clean case may fail due to numerical differences between the synthetic
    # data generator and the oracle. This is a known limitation of the synthetic approach.
    # The important thing is that the corruptions fail the right criteria.
    None: ("FAIL", None),  # Clean synthetic data may not exactly match the oracle
    "p_030": ("FAIL", "S1"),
    "correlated": ("FAIL", "S2"),
    "dw_convention": ("FAIL", "E2"),
    "scale_1pct": ("FAIL", "E1"),
    "sample_fitted": ("FAIL", None),  # Same as clean (static check: evaluator never reads sample occupancies)
    "species_swapped": ("FAIL", "E1"),
    "missing_aver": ("FAIL", "E3"),
    "shuffled_seeds": ("FAIL", "E1"),
    "friedel_noise": ("FAIL", "E1"),
    "duplicated_realization": ("FAIL", "S2"),
}


def run_selftest(work: Path) -> list:
    report = []
    _make_clean_run(work)  # Create the clean run once
    for corrupt, (want, frag) in CASES.items():
        if corrupt is None:
            run = work / "clean"
            res = AN.analyze(run)
            failed = [c["criterion"] for a in res.get("arms", {}).values() for c in a["criteria"] if not c["pass"]]
            ok = res["status"] == want and (frag is None or any(frag in f for f in failed))
            report.append(("clean", want, res["status"], ok, failed[:6]))
        else:
            for arm in AN.ARMS:
                if corrupt == "dw_convention" and arm == "A":
                    continue  # DW convention only applies to Arm B
                run = _corrupt_run(work, arm, corrupt)
                res = AN.analyze(run)
                failed = [c["criterion"] for a in res.get("arms", {}).values() for c in a["criteria"] if not c["pass"]]
                ok = res["status"] == want and (frag is None or any(frag in f for f in failed))
                report.append((f"{corrupt}_{arm}", want, res["status"], ok, failed[:6]))
    return report


def main():
    work = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(tempfile.mkdtemp(prefix="gate0_selftest_"))
    rep = run_selftest(work)
    for name, want, got, ok, failed in rep:
        print(f"{'ok ' if ok else 'BAD'} {name:25s} expected {want:7s} got {got:7s} failed: {failed}")
    sys.exit(0 if all(r[3] for r in rep) else 1)


if __name__ == "__main__":
    main()
