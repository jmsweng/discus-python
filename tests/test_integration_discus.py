"""Opt-in integration test against a real DISCUS_SUITE installation.

Skipped automatically when ``discus_suite`` is not on PATH (standard CI).
Run locally in the WSL2 environment, e.g.:

    DISCUS_INTEGRATION=1 python -m pytest tests/test_integration_discus.py -q

It runs the verified 1D PDF macro (docs/VERIFIED_DISCUS_MACRO.md) on the
committed ``examples/gaas.cll`` sample and compares the produced G(r) to the
committed reference ``examples/run_sample/result.pdf`` with tolerant metrics:
same grid, finite values, and a small relative L2 deviation. Tolerances are
deliberately empirical, not bitwise identity.

This is a software/engineering check of the runner workflow against the
installed engine version; it is **not** a quantitative 3D-ΔPDF validation.
"""
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from discuspdf import PDFConfig, PDFRunner, read_pdf

REPO_ROOT = Path(__file__).resolve().parent.parent
REFERENCE = REPO_ROOT / "examples" / "run_sample" / "result.pdf"
STRUCTURE = REPO_ROOT / "examples" / "gaas.cll"

# Empirical tolerances for the same engine version and identical settings.
REL_L2_TOLERANCE = 1e-6
GRID_ATOL = 1e-8


def _discus_available() -> bool:
    if os.environ.get("DISCUS_INTEGRATION") != "1":
        return False
    if shutil.which("discus_suite") is None:
        return False
    try:
        proc = subprocess.run(["discus_suite", "-noautorun", "-macro", "-h"],
                              capture_output=True, text=True, timeout=60)
    except (OSError, subprocess.TimeoutExpired):
        return False
    # The suite does not support -h (it treats it as a macro name); a fast,
    # controlled failure is proof the binary launches.
    return True


@pytest.mark.skipif(not _discus_available(), reason="DISCUS_SUITE not available (set DISCUS_INTEGRATION=1)")
def test_real_discus_1d_pdf_matches_reference(tmp_path):
    cfg = PDFConfig(
        structure=STRUCTURE,
        supercell=(5, 5, 5),
        rmax=12.0,
        dr=0.02,
        qmax=25.0,
        qsig=0.0,
        radiation="xray",
    )
    runner = PDFRunner(runs_directory=tmp_path / "runs")
    prepared = runner.prepare(cfg)
    result = runner.execute(prepared)

    st = json.loads((result.directory / "status.json").read_text())
    assert st["status"] == "complete"
    assert (result.directory / "run.log").is_file()

    r, gr = result.load()
    ref_r, ref_gr = read_pdf(REFERENCE)

    assert len(r) == len(ref_r) == 601
    assert np.allclose(r, ref_r, rtol=0, atol=GRID_ATOL)
    assert np.all(np.isfinite(gr))
    rel = float(np.linalg.norm(gr - ref_gr) / max(np.linalg.norm(ref_gr), 1e-30))
    assert rel < REL_L2_TOLERANCE, f"Relative L2 deviation {rel:.3g} exceeds {REL_L2_TOLERANCE}"
