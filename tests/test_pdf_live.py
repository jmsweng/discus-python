"""Opt-in live tests for the DISCUS Qmax limit. Run with DISCUS_LIVE=1."""
import os
import tempfile
import warnings
from pathlib import Path

import numpy as np
import pytest

from discuspdf import PDFConfig, QmaxLimitWarning, read_pdf, render_pdf_macro
from tests.discus_live import run_discus

pytestmark = [pytest.mark.discus_live, pytest.mark.skipif(
    os.environ.get("DISCUS_LIVE") != "1", reason="live DISCUS test; set DISCUS_LIVE=1 to enable")]

GAAS = Path(__file__).resolve().parents[1] / "examples" / "gaas.cll"


def _pdf(name, qmax):
    base = os.environ.get("DISCUS_LIVE_WORKDIR", os.path.join(tempfile.gettempdir(), "discuspdf_live"))
    workdir = os.path.join(base, name)
    os.makedirs(workdir, exist_ok=True)
    result = Path(workdir) / "result.pdf"
    if result.exists():
        result.unlink()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", QmaxLimitWarning)
        cfg = PDFConfig(GAAS, supercell=(5, 5, 5), rmax=12.0, dr=0.02, qmax=qmax, qsig=0.0, radiation="xray")
    run_discus(render_pdf_macro(cfg), workdir, files={"model.cll": GAAS.read_bytes()}, timeout=600)
    return read_pdf(result)[1]


def _rel(a, b):
    return float(np.max(np.abs(a - b)) / np.max(np.abs(b)))


def test_live_qmax_above_limit_gives_same_pdf():
    # default wavelength 1.5406 A: Qmax 25 and 8.5 are both above the limit (~8.0) -> identical G(r)
    assert _rel(_pdf("q25", 25.0), _pdf("q8.5", 8.5)) < 1e-6


def test_live_qmax_below_limit_changes_pdf():
    assert _rel(_pdf("q7.5", 7.5), _pdf("q25b", 25.0)) > 1e-2
