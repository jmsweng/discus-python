"""Re-run the D reference probe analysis on the committed DISCUS evidence (no engine needed)."""
import os
from pathlib import Path
import subprocess
import sys

import pytest

pytest.importorskip("h5py")
ROOT = Path(__file__).resolve().parents[1]
KIT = ROOT / "benchmarks/discus_6_23_01/d_reference"


def test_committed_evidence_passes_all_criteria():
    ev = KIT / "evidence/20261101T221751Z_LAPTOP-NGU3BT0N"
    out = subprocess.run([sys.executable, str(KIT / "analyze.py"), str(ev)], capture_output=True, text=True,
                         env={**os.environ, "PYTHONPATH": str(ROOT / "src")}, check=False)
    assert out.returncode == 0, out.stdout + out.stderr
    assert "criteria failed: 0 []" in out.stdout
