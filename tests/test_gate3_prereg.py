"""Gate 3 (native `value 3DPDF`) tooling checks. The archived-pair test uses real DISCUS 6.23.01 output."""
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

h5py = pytest.importorskip("h5py")
B = Path(__file__).resolve().parents[1] / "benchmarks" / "discus_6_23_01"
G3 = B / "gate3"
sys.path.insert(0, str(G3)); sys.path.insert(0, str(B))
import native3dpdf as N  # noqa: E402
import regenerate_evidence as R  # noqa: E402

PRED = json.loads((G3 / "predictions.json").read_text())


def test_registered_files_unchanged():
    for n, h in PRED["code_sha256"].items():
        assert R.sha256_file(G3 / n) == h, f"{n} changed after registration"
    for n, h in PRED["input_sha256"].items():
        assert R.sha256_file(G3 / "inputs" / n) == h


@pytest.mark.parametrize("pair", [("correlated_run", "corr_total"), ("short_range_run", "sr_total")])
def test_emulation_reproduces_archived_discus_native_3dpdf(pair):
    """Real engine output (33^3 grid, set aver 0, structure read as one cell): values and axis metadata."""
    d, stem = pair
    gi = h5py.File(B / d / f"{stem}.h5", "r")["entry/data"]
    gn = h5py.File(B / d / f"{stem}_3dpdf.h5", "r")["entry/data"]
    I, vi = gi["data_values"][()], np.diag(gi["data_increment_vector"][()])
    J, scalef = N.bragg_reweight(I, gi["data_corner"][()], vi, ncells=1)
    assert np.isclose(scalef, 0.125 ** 3)
    nat = gn["data_values"][()]
    assert np.linalg.norm(N.native_3dpdf(J) - nat) / np.linalg.norm(nat) < 1e-13
    corner, inc = N.native_axes(I.shape, vi)
    assert np.allclose(gn["data_corner"][()], corner) and np.allclose(np.diag(gn["data_increment_vector"][()]), inc)
    # and the documented metadata defect: the true step is 1/(n dH), not 1/((n-1) dH)
    assert np.allclose(N.true_axes(I.shape, vi)[1], 1 / (33 * 0.125))


def test_selftest_passes_clean_and_rejects_corruptions(tmp_path):
    if subprocess.run(["git", "-C", str(B), "cat-file", "-e", "8cd59ac"], capture_output=True).returncode != 0:
        pytest.skip("Gate 2 v4 volumes need the git history of the PR #8 branch (commit 8cd59ac)")
    import gate3_selftest as st
    bad = [r for r in st.run_selftest(tmp_path) if not r[3]]
    assert not bad, bad
