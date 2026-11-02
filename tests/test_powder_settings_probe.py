"""Powder-settings probe (benchmarks/discus_6_23_01/powder_settings): inputs and analyzer, no DISCUS needed."""
import hashlib
import importlib.util
import json
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parents[1] / "benchmarks" / "discus_6_23_01" / "powder_settings"
spec = importlib.util.spec_from_file_location("powder_probe", HERE / "probe.py")
probe = importlib.util.module_from_spec(spec)
spec.loader.exec_module(probe)


def test_committed_macros_are_what_the_generator_writes():
    pred = json.loads((HERE / "predictions.json").read_text())
    assert set(pred["macro_sha256"]) == {c for c, _o, _p in probe.CASES}
    for cid, over, _p in probe.CASES:
        text = (HERE / "inputs" / f"{cid}.mac").read_text(encoding="ascii")
        assert text == probe.case_macro(cid, over)
        assert hashlib.sha256(text.encode()).hexdigest() == pred["macro_sha256"][cid]


def test_probe_macros_keep_the_verified_pdf_block():
    """Only the F(Q) output block is added to the notebook-01 macro; base_plain is that macro unchanged."""
    for cid, over, _p in probe.CASES:
        plain = probe.render_pdf_macro(probe.case_config(over))
        m = probe.case_macro(cid, over)
        assert m == plain if cid == "base_plain" else m.replace(probe.FQ_BLOCK, "", 1) == plain


def test_zero_last_matches_the_discus_loop():
    import numpy as np
    y = np.array([1.0, -2.0, 3.0, 2.0, 1.0])  # last sign change between index 1 and 2
    assert probe.zero_last(y).tolist() == [1.0, -2.0, 0.0, 0.0, 0.0]
    y = np.array([1.0, 0.0, 2.0, 3.0])        # exact zero stops the walk
    assert probe.zero_last(y).tolist() == [1.0, 0.0, 0.0, 0.0]


@pytest.mark.parametrize("name, kw, want", [
    ("clean", {}, "PASS"), ("clean_shift0", {"shift": 0}, "PASS"), ("scale", {"scale": 1.01}, "FAIL"),
    ("no_zeroing", {"zero": False}, "FAIL"), ("plain", {"perturb_plain": 1e-6}, "FAIL"),
    ("fq_dq", {"fq_noise": 1e-3}, "FAIL"),
])
def test_analyzer_on_synthetic_evidence(tmp_path, name, kw, want):
    probe._synthetic(tmp_path, **kw)
    assert probe.analyze(tmp_path)["verdict"] == want


spec2 = importlib.util.spec_from_file_location("powder_probe_v2", HERE / "probe_v2.py")
probe_v2 = importlib.util.module_from_spec(spec2)
spec2.loader.exec_module(probe_v2)


def test_v2_committed_macros_are_what_the_generator_writes():
    pred = json.loads((HERE / "predictions_v2.json").read_text(encoding="utf-8"))
    for cid, over, _e, _p in probe_v2.CASES:
        text = (HERE / "inputs_v2" / f"{cid}.mac").read_text(encoding="ascii")
        assert text == probe_v2.case_macro(over)
        assert hashlib.sha256(text.encode()).hexdigest() == pred["macro_sha256"][cid]


@pytest.mark.parametrize("kw, want", [
    ({}, "PASS"), ({"scale": 2 / 3.141592653589793}, "FAIL"), ({"no_shift": True}, "FAIL"),
    ({"aborts": False}, "FAIL"), ({"qmin_changes_f": True}, "FAIL"),
])
def test_v2_analyzer_on_synthetic_evidence(tmp_path, kw, want):
    probe_v2._synthetic(tmp_path, **kw)
    assert probe_v2.analyze(tmp_path)["verdict"] == want
