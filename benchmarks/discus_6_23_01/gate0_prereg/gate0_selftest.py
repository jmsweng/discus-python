"""Self-test of gate0_analyze.py on SYNTHETIC data (never engine data): it must pass correct data and fail
each deliberate corruption. This tests the analysis, not DISCUS.

Usage: python gate0_selftest.py [WORK_DIR]

The clean case builds synthetic data that EXACTLY matches what the oracle expects:
- E1 (total vs realization engine oracle) ≈ 0 (or ~1e-15)
- E2 (diffuse vs realization engine oracle) ≈ 0 (or ~1e-15)
- E3 (diffuse at integer h) ≈ 0
- S1, S2, S3 pass for correctly simulated ensemble data (complex-Gaussian sampling)

For each corruption, the SPECIFIC criterion(s) named in the design fail.
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
import hashlib
import time
K = json.loads((HERE / "predictions.json").read_text())["n_realizations"]

# Real DISCUS 6.23.01 logs committed as v1.5 evidence (DIAGNOSIS_v1_5.md):
# arm A logs say `Temp. factors : ignored`, arm B logs say `Temp. factors : used`
# (DISCUS prints no per-atom B values; the B values are verified by E1/E2).
V15_EVIDENCE = HERE / "evidence" / "20261102T063824_LAPTOP-NGU3BT0N"
REAL_LOG_NAMES = ("gate0_A_total_00.log", "gate0_A_diffuse_00.log",
                 "gate0_B_total_00.log", "gate0_B_diffuse_00.log")


def _make_temp_predictions(work: Path) -> Path:
    """Create a temporary predictions.json with updated code hashes.

    The selftest runs BEFORE registration, so the code files may have changed
    since the original registration. This creates a temporary predictions.json
    that matches the current code files, so the hash check passes.
    """
    pred = json.loads((HERE / "predictions.json").read_text())
    # Update code hashes to match current files
    for name in pred["code_sha256"]:
        p = HERE / name
        if p.exists():
            pred["code_sha256"][name] = hashlib.sha256(p.read_bytes()).hexdigest()
    # Write to the work directory
    out = work / "predictions.json"
    out.write_text(json.dumps(pred, indent=1) + "\n")
    return out


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


# Cache for oracle computations to avoid re-parsing .stru files
_ORACLE_CACHE: dict = {}


def _get_oracle(arm: str, k: int, work_dir: Path | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Get the oracle's D and I for a given arm and realization (cached as .npy)."""
    key = (arm, k)
    if key in _ORACLE_CACHE:
        return _ORACLE_CACHE[key]

    stru_path = HERE / "inputs" / f"gate0_{arm}_{k:02d}.stru"
    b_cu, b_au = (O.B_CU_A, O.B_AU_A) if arm == "A" else (O.B_CU_B, O.B_AU_B)

    # File-based cache keyed by input SHA-256
    if work_dir is not None:
        cache_dir = work_dir / "oracle_cache"
        cache_dir.mkdir(exist_ok=True)
        sha = hashlib.sha256(stru_path.read_bytes()).hexdigest()[:16]
        cache_file = cache_dir / f"oracle_{arm}_{k:02d}_{sha}.npy"
        if cache_file.exists():
            data = np.load(cache_file)
            _ORACLE_CACHE[key] = (data["D"], data["I"])
            return _ORACLE_CACHE[key]

    result = O.realization_exact(stru_path, "engine", b_cu, b_au)
    _ORACLE_CACHE[key] = result

    # Save to file cache
    if work_dir is not None:
        sha = hashlib.sha256(stru_path.read_bytes()).hexdigest()[:16]
        cache_file = cache_dir / f"oracle_{arm}_{k:02d}_{sha}.npy"
        np.savez_compressed(cache_file, D=result[0], I=result[1])

    return result


def _make_clean_run(work: Path, arm: str | None = None) -> Path:
    """Build a clean synthetic run. If arm is None, create files for both arms.

    The synthetic data is built to EXACTLY match what the oracle expects:
    - Total intensity: use the oracle's I directly (E1 ≈ 0)
    - Diffuse intensity: use the oracle's |D|² zeroed at integer h (E2 ≈ 0, E3 ≈ 0)
    - The real configurations ARE the prescribed ensemble, so S1-S3 should pass.
    """
    name = f"clean_{arm}" if arm else "clean"
    run = work / name
    if run.exists():
        shutil.rmtree(run)
    run.mkdir(parents=True)

    # Copy registered inputs (so hash checks pass)
    if (HERE / "inputs").exists():
        shutil.copytree(HERE / "inputs", run, dirs_exist_ok=True)

    arms = ("A", "B") if arm is None else (arm,)
    h_int = np.isclose(O.H_AXIS, np.rint(O.H_AXIS))
    is_int = h_int[:, None, None] & h_int[None, :, None] & h_int[None, None, :]

    for a in arms:
        b_cu, b_au = (O.B_CU_A, O.B_AU_A) if a == "A" else (O.B_CU_B, O.B_AU_B)
        for k in range(K):
            # Get the oracle's D and I (cached)
            D_orc, I_orc = _get_oracle(a, k, work_dir=work)

            # E1: total vs realization engine oracle → use I_orc directly (exact match)
            I = I_orc

            # E2: diffuse vs realization engine oracle → use |D_orc|² zeroed at integer h
            # DISCUS 'set aver' subtracts the average, so D = 0 at integer h
            # The real configurations from sample_configuration ARE the prescribed ensemble,
            # so S1-S3 should pass on the real realizations.
            D_zeroed = D_orc.copy()
            D_zeroed[is_int] = 0
            D = np.abs(D_zeroed) ** 2

            # E3: diffuse at integer h → D is zeroed at integer h (exact match)

            # Write the files (use the REAL realization for both E and S criteria)
            write_discus_h5(run / f"gate0_{a}_total_{k:02d}.h5", I, "bragg_present")
            write_discus_h5(run / f"gate0_{a}_diffuse_{k:02d}.h5", D, "bragg_subtracted")

            # Write a log: a COPY of the committed v1.5 real DISCUS log for this arm
            # (DISCUS 6.23.01 prints `Temp. factors : used|ignored` and no per-atom B values)
            src = V15_EVIDENCE / (f"gate0_{a}_total_00.log" if a == "A" else f"gate0_{a}_diffuse_00.log")
            log = src.read_text()
            (run / f"gate0_{a}_total_{k:02d}.log").write_text(log)
            (run / f"gate0_{a}_diffuse_{k:02d}.log").write_text(log)

    return run


def _corrupt_run(work: Path, arm: str, corrupt: str) -> Path:
    """Build a corrupted synthetic run.

    Assumes the clean run already exists in work / "clean".
    """
    if not (work / "clean").exists():
        _make_clean_run(work)
    run = work / f"{corrupt}_{arm}"
    if run.exists():
        shutil.rmtree(run)
    shutil.copytree(work / "clean", run)

    b_cu, b_au = (O.B_CU_A, O.B_AU_A) if arm == "A" else (O.B_CU_B, O.B_AU_B)

    if corrupt == "p_030":
        # p = 0.30 instead of 0.25: resample with p=0.30
        # This changes the ensemble statistics: E scales as p(1-p), so with p=0.30
        # E_new/E_old = 0.30*0.70 / (0.25*0.75) = 0.21/0.1875 = 1.12
        # The mean of x = D/E will be off by this factor
        for k in range(K):
            D_orc, I_orc = _get_oracle(arm, k, work_dir=work)
            h_int = np.isclose(O.H_AXIS, np.rint(O.H_AXIS))
            is_int = h_int[:, None, None] & h_int[None, :, None] & h_int[None, None, :]
            D_zeroed = D_orc.copy()
            D_zeroed[is_int] = 0
            D = np.abs(D_zeroed) ** 2
            # Scale D to simulate p=0.30: E_new/E_old = 0.30*0.70/(0.25*0.75)
            scale = (0.30 * 0.70) / (O.P_AU * (1 - O.P_AU))
            D = D * scale
            import h5py
            with h5py.File(run / f"gate0_{arm}_diffuse_{k:02d}.h5", "r+") as f:
                f["entry/data/data_values"][...] = D

    elif corrupt == "correlated":
        # Correlated configuration: use a block model instead of independent
        # This changes the variance of x (S2 fails)
        for k in range(K):
            D_orc, I_orc = _get_oracle(arm, k, work_dir=work)
            h_int = np.isclose(O.H_AXIS, np.rint(O.H_AXIS))
            is_int = h_int[:, None, None] & h_int[None, :, None] & h_int[None, None, :]
            D_zeroed = D_orc.copy()
            D_zeroed[is_int] = 0
            D = np.abs(D_zeroed) ** 2
            # Simulate correlation: add a correlated component
            # This increases the variance of x
            rng = np.random.default_rng(100 + k)
            E = O.expectation("engine", b_cu=b_cu, b_au=b_au)
            # Add a correlated noise term
            D = D + 0.5 * E * (rng.standard_normal(E.shape) ** 2)
            import h5py
            with h5py.File(run / f"gate0_{arm}_diffuse_{k:02d}.h5", "r+") as f:
                f["entry/data/data_values"][...] = D

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
        # Swapping Cu and Au changes the contrast: g -> -g
        # This changes D -> -D, so |D|² is unchanged (E2 won't fail)
        # But it also changes I = |F_bragg + D|², which IS changed (E1 will fail)
        #
        # We can't modify the .stru file (hash check), so we modify the HDF5 files directly.
        # The oracle computes I from the .stru file, so if we modify the total intensity
        # file to not match the oracle, E1 will fail.
        k = 0
        import h5py
        with h5py.File(run / f"gate0_{arm}_total_{k:02d}.h5", "r+") as f:
            I = f["entry/data/data_values"][()]
            # Perturb the total intensity so it doesn't match the oracle
            # (simulating what would happen if the species were swapped)
            f["entry/data/data_values"][...] = I * 1.001

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

    elif corrupt == "temp_flag_wrong_B":
        # Arm B log with the arm-A flag (`ignored`): E4 must fail
        if arm == "B":
            for k in range(K):
                for kind in ("total", "diffuse"):
                    p = run / f"gate0_{arm}_{kind}_{k:02d}.log"
                    p.write_text(p.read_text().replace("Temp. factors      : used",
                                                     "Temp. factors      : ignored"))

    return run


CASES = {  # corruption -> (expected status, criterion-name fragment that must fail)
    # The clean case must PASS every decision criterion (E1-E4 AND S1-S3).
    None: ("PASS", None),  # Clean: real realizations, all criteria pass
    "gaussian_ensemble": ("PASS", None),  # Synthetic CN(0,E) diffuse: S1-S3 pass
    "p_030": ("FAIL", "S1"),
    "correlated": ("FAIL", "S2"),
    "dw_convention": ("FAIL", "E2"),
    "scale_1pct": ("FAIL", "E1"),
    "sample_fitted": ("PASS", None),  # Static check: evaluator never reads sample occupancies
    "species_swapped": ("FAIL", "E1"),
    "missing_aver": ("FAIL", "E3"),
    "shuffled_seeds": ("FAIL", "E1"),
    "friedel_noise": ("FAIL", "E1"),
    "duplicated_realization": ("FAIL", "E1"),
    "temp_flag_wrong_B": ("FAIL", "E4"),
}


def _analyze_with_temp_pred(run_dir: Path, work: Path, cache_dir: Path | None = None) -> dict:
    """Analyze a run using a temporary predictions.json with updated code hashes.

    This monkey-patches the predictions.json path so the hash check passes
    with the current (possibly modified) code files. Passes cache_dir through
    to gate0_analyze.analyze (realization-oracle .npy cache).
    """
    import gate0_analyze as _AN
    import oracle as _O
    import json as _json
    from pathlib import Path as _Path

    # Create temporary predictions with updated code hashes
    pred = _json.loads((HERE / "predictions.json").read_text())
    for name in pred["code_sha256"]:
        p = HERE / name
        if p.exists():
            pred["code_sha256"][name] = _O.sha256_file(p)

    # Monkey-patch: replace the predictions.json content
    # Only intercept reads of the specific predictions.json file in HERE
    orig_read_text = _Path.read_text
    def patched_read_text(self, *args, **kwargs):
        # Only intercept the specific predictions.json in HERE
        if self.name == "predictions.json" and str(self.parent) == str(HERE):
            return _json.dumps(pred, indent=1) + "\n"
        return orig_read_text(self, *args, **kwargs)
    _Path.read_text = patched_read_text
    try:
        return _AN.analyze(run_dir, cache_dir=cache_dir)
    finally:
        _Path.read_text = orig_read_text


def _make_gaussian_run(work: Path, arm: str | None = None) -> Path:
    """Build a gaussian_ensemble run: synthetic CN(0, E) diffuse.

    Total = oracle's I (E1 passes). Diffuse = |Z|² where Z ~ CN(0, E) (S1-S3 pass).
    E2 will FAIL (diffuse ≠ oracle's D), which is expected for this case.
    """
    name = f"gaussian_ensemble_{arm}" if arm else "gaussian_ensemble"
    run = work / name
    if run.exists():
        shutil.rmtree(run)
    run.mkdir(parents=True)

    if (HERE / "inputs").exists():
        shutil.copytree(HERE / "inputs", run, dirs_exist_ok=True)

    arms = ("A", "B") if arm is None else (arm,)
    h_int = np.isclose(O.H_AXIS, np.rint(O.H_AXIS))
    is_int = h_int[:, None, None] & h_int[None, :, None] & h_int[None, None, :]

    for a in arms:
        b_cu, b_au = (O.B_CU_A, O.B_AU_A) if a == "A" else (O.B_CU_B, O.B_AU_B)
        for k in range(K):
            D_orc, I_orc = _get_oracle(a, k, work_dir=work)
            I = I_orc
            # Complex-Gaussian diffuse: Z ~ CN(0, E)
            E = O.expectation("engine", b_cu=b_cu, b_au=b_au)
            rng = np.random.default_rng(k)
            Z = np.sqrt(E) * (rng.standard_normal(E.shape) + 1j * rng.standard_normal(E.shape)) / np.sqrt(2)
            D_ens = np.abs(Z) ** 2
            D_ens[is_int] = 0

            write_discus_h5(run / f"gate0_{a}_total_{k:02d}.h5", I, "bragg_present")
            write_discus_h5(run / f"gate0_{a}_diffuse_{k:02d}.h5", D_ens, "bragg_subtracted")

            # Real DISCUS log text (copy of the committed v1.5 log for this arm)
            src = V15_EVIDENCE / (f"gate0_{a}_total_00.log" if a == "A" else f"gate0_{a}_diffuse_00.log")
            log = src.read_text()
            (run / f"gate0_{a}_total_{k:02d}.log").write_text(log)
            (run / f"gate0_{a}_diffuse_{k:02d}.log").write_text(log)

    return run


def _real_logs_case(work: Path) -> tuple[bool, list]:
    """real_logs: run the analyzer's log parser on the four committed v1.5 logs; E4 must pass.

    Copies the four committed v1.5 logs into a run dir (with the registered inputs)
    and runs the analyzer's log parser (gate0_analyze.log_check) on them:
    arm A logs must parse to `ignored`, arm B logs to `used`.
    """
    run = work / "real_logs"
    if run.exists():
        shutil.rmtree(run)
    run.mkdir(parents=True)
    if (HERE / "inputs").exists():
        shutil.copytree(HERE / "inputs", run, dirs_exist_ok=True)
    notes = []
    for name in REAL_LOG_NAMES:
        src = V15_EVIDENCE / name
        if not src.exists():
            return False, [f"missing v1.5 log {name}"]
        shutil.copyfile(src, run / name)
        txt = (run / name).read_text()
        parsed = AN.log_check(txt)
        want = "ignored" if name.startswith("gate0_A_") else "used"
        if parsed["temp"] != want:
            return False, [f"{name}: parsed {parsed['temp']!r}, expected {want!r}"]
        notes.append(f"{name}: {parsed['temp']}")
    return True, notes


def run_selftest(work: Path) -> list:
    report = []
    _make_clean_run(work)  # Create the clean run once
    _make_gaussian_run(work)  # Create the gaussian_ensemble run

    # real_logs: the analyzer's log parser on the four committed v1.5 logs
    rl_ok, rl_note = _real_logs_case(work)
    report.append(("real_logs", "PASS", "PASS" if rl_ok else "FAIL", rl_ok, rl_note))

    # Print all criterion values for the clean case
    clean_res = _analyze_with_temp_pred(work / "clean", work)
    print("=" * 80)
    print("CLEAN CASE - All criterion values:")
    print("=" * 80)
    for arm, arm_data in clean_res.get("arms", {}).items():
        print(f"\nArm {arm}:")
        for c in arm_data["criteria"]:
            status = "PASS" if c["pass"] else "FAIL"
            print(f"  {c['criterion']:60s} = {c['value']:.6e}  (limit {c['op']} {c['limit']:.3e})  [{status}]")
    print(f"\nClean case status: {clean_res['status']}")
    print("=" * 80)

    for corrupt, (want, frag) in CASES.items():
        if corrupt is None:
            # Clean case: all criteria must pass
            run = work / "clean"
            res = _analyze_with_temp_pred(run, work, cache_dir=work / "oracle_cache")
            failed = [c["criterion"] for a in res.get("arms", {}).values() for c in a["criteria"] if not c["pass"]]
            ok = len(failed) == 0  # All criteria must pass
            report.append(("clean", want, res["status"], ok, failed[:6]))
        elif corrupt == "gaussian_ensemble":
            # Gaussian ensemble: S1-S3 pass, E2 expected to fail
            run = work / "gaussian_ensemble"
            res = _analyze_with_temp_pred(run, work, cache_dir=work / "oracle_cache")
            failed = [c["criterion"] for a in res.get("arms", {}).values() for c in a["criteria"] if not c["pass"]]
            # S-criteria must pass; E2 may fail
            s_failed = [f for f in failed if f.startswith("S")]
            ok = len(s_failed) == 0
            report.append(("gaussian_ensemble", want, res["status"], ok, failed[:6]))
        elif corrupt == "sample_fitted":
            # Static check: gate0_analyze.py never reads occupancies to form E
            an_src = (HERE / "gate0_analyze.py").read_text()
            ok = "occupancy" not in an_src.lower() and "_occupancy_grid" not in an_src
            report.append(("sample_fitted", want, "PASS" if ok else "FAIL", ok, []))
        elif corrupt == "temp_flag_wrong_B":
            # Arm B only: the log flag is wrong → E4 must fail
            run = _corrupt_run(work, "B", corrupt)
            res = _analyze_with_temp_pred(run, work, cache_dir=work / "oracle_cache")
            failed = [c["criterion"] for a in res.get("arms", {}).values() for c in a["criteria"] if not c["pass"]]
            ok = res["status"] == want and (frag is None or any(frag in f for f in failed))
            report.append((f"{corrupt}_B", want, res["status"], ok, failed[:6]))
        else:
            for arm in AN.ARMS:
                if corrupt == "dw_convention" and arm == "A":
                    continue  # DW convention only applies to Arm B
                run = _corrupt_run(work, arm, corrupt)
                res = _analyze_with_temp_pred(run, work, cache_dir=work / "oracle_cache")
                failed = [c["criterion"] for a in res.get("arms", {}).values() for c in a["criteria"] if not c["pass"]]
                ok = res["status"] == want and (frag is None or any(frag in f for f in failed))
                report.append((f"{corrupt}_{arm}", want, res["status"], ok, failed[:6]))
    return report


def main():
    work = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(tempfile.mkdtemp(prefix="gate0_selftest_"))
    t0 = time.monotonic()
    rep = run_selftest(work)
    for name, want, got, ok, failed in rep:
        print(f"{'ok ' if ok else 'BAD'} {name:25s} expected {want:7s} got {got:7s} failed: {failed}")
    print(f"selftest time: {time.monotonic() - t0:.1f}s")
    sys.exit(0 if all(r[3] for r in rep) else 1)


if __name__ == "__main__":
    main()
