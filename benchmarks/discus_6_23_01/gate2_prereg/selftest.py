"""Self-test of analyze.py on SYNTHETIC data (never engine data): it must pass correct data and fail
each deliberate corruption. This tests the analysis, not DISCUS.

Usage: python selftest.py [WORK_DIR]
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
import analyze as AN  # noqa: E402

LOG = """          ***********************************************************
          *         D I S C U S - S U I T E  Version 6.23.01        *
          ***********************************************************
 Fourier technique    : turbo Fourier
 Fourier calculated by: atom form factors taken from: Waasmeier & Kirfel
{used}
   Temp. factors      : {temp}
 Lot #    1/   1 : complete crystal
                   Atom typ = CU         (#      2048 )
                   Atom typ = AU         (#      2048 )
"""


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
        g.create_dataset("data_corner", data=np.full(3, -O.H_MAX))
        g.create_dataset("data_dimension", data=np.array(values.shape, np.int32))
        g.create_dataset("data_increment_vector", data=np.eye(3) * O.DH)
        s("data_type_axes", "hkl"); s("data_type_content", "intensity"); s("data_type_reciprocal", "reciprocal")
        s("data_type_with_bragg", bragg)
        g.create_dataset("data_values", data=values.astype(np.float64))
        g.create_dataset("unit_cell_angles", data=np.full(3, 90.0))
        g.create_dataset("unit_cell_lengths", data=np.full(3, O.A_CELL))


def make_run(work: Path, corrupt: str | None = None) -> Path:
    run = work / (corrupt or "clean")
    if run.exists():
        shutil.rmtree(run)
    shutil.copytree(HERE / "inputs", run)
    for model in AN.MODELS:
        sig = O.build_sigma(model)
        It, Id, _ = O.forward(sig, "engine")
        if corrupt == "no_averaging":
            Id = It.copy()
        elif corrupt == "zero_diffuse":
            Id = np.zeros_like(It)
        elif corrupt == "transpose_xy":
            It, Id = It.transpose(1, 0, 2).copy(), Id.transpose(1, 0, 2).copy()
        elif corrupt == "scale_1pct":
            It, Id = It * 1.01, Id * 1.01
        elif corrupt == "y_sign_flip" and model == "block":
            sig2 = (sig * ((-1) ** np.arange(O.N_CELL))[None, :, None]).astype(np.int8)
            It, Id, _ = O.forward(sig2, "engine")
        elif corrupt == "unrequested_dw_B1":      # exactly what gate2-prereg-v2 saw on the real engine
            H = O.H_AXIS; q = np.sqrt(H[:, None, None]**2 + H[None, :, None]**2 + H[None, None, :]**2) / (2 * O.A_CELL)
            dw = np.exp(-2 * (np.rint(q / O.CFINC) * O.CFINC) ** 2)
            It, Id = It * dw, Id * dw
        elif corrupt in ("ties_rounded_down", "tie_voxels_wrong"):
            _, s_dn = O.s_tie_candidates()
            Itd, Idd, _ = O.forward(sig, "engine", s=s_dn)
            if corrupt == "ties_rounded_down":      # legitimate alternative rounding: must still PASS
                It, Id = Itd, Idd
            else:                                    # 0.5 % off at tie voxels only: must FAIL R
                t = O.tie_mask(); It = np.where(t, It * 1.005, It); Id = np.where(t, Id * 1.005, Id)
        elif corrupt == "physics_precision":
            It, Id, _ = O.forward(sig, "physics")
        write_discus_h5(run / f"gate2_{model}_total.h5", It, "bragg_present")
        write_discus_h5(run / f"gate2_{model}_diffuse.h5", Id, "bragg_subtracted")
        (run / f"gate2_{model}_total.log").write_text(LOG.format(used="", temp=("used" if corrupt == "temp_used_log" else "ignored")))
        (run / f"gate2_{model}_diffuse.log").write_text(
            LOG.format(used=" Used 100.0 % of the crystal to calculate <F> ...",
                       temp=("used" if corrupt == "temp_used_log" else "ignored")))
    if corrupt == "edited_macro":
        p = run / "gate2_block_diffuse.mac"; p.write_text(p.read_text().replace("set aver,100.0", "set aver,50.0"))
    return run


CASES = {  # corruption -> expected status and a criterion-name fragment that must fail
    None: ("PASS", None),
    "no_averaging": ("FAIL", "C1"),
    "zero_diffuse": ("FAIL", "C1"),
    "transpose_xy": ("FAIL", "C0"),
    "scale_1pct": ("FAIL", "C0"),
    "y_sign_flip": ("FAIL", "B2 sign"),
    "physics_precision": ("FAIL", "engine-precision"),
    "unrequested_dw_B1": ("FAIL", "C0"),
    "temp_used_log": ("FAIL", "T "),
    "ties_rounded_down": ("PASS", None),
    "tie_voxels_wrong": ("FAIL", "R "),
    "edited_macro": ("INVALID", None),
}


def run_selftest(work: Path) -> list:
    report = []
    for corrupt, (want, frag) in CASES.items():
        res = AN.analyze(make_run(work, corrupt))
        failed = [c["criterion"] for m in res.get("models", {}).values() for c in m["criteria"] if not c["pass"]]
        ok = res["status"] == want and (frag is None or any(frag in f for f in failed))
        report.append((corrupt or "clean", want, res["status"], ok, failed[:6]))
    return report


def main():
    work = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(tempfile.mkdtemp(prefix="gate2_selftest_"))
    rep = run_selftest(work)
    for name, want, got, ok, failed in rep:
        print(f"{'ok ' if ok else 'BAD'} {name:18s} expected {want:7s} got {got:7s} failed: {failed}")
    sys.exit(0 if all(r[3] for r in rep) else 1)


if __name__ == "__main__":
    main()
