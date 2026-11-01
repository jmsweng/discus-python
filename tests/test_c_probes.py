"""The committed DISCUS ``chem corr`` evidence agrees with the registered predictions (offline).

The predictions in ``benchmarks/discus_6_23_01/c_probes/predictions.json`` were committed before the
engine run; this test re-derives them from the committed inputs with the current code (so a change in
``discuspdf.correlations`` that alters any value fails here) and re-checks the engine logs.
"""
from __future__ import annotations

import importlib.util
import json
import math
from pathlib import Path

import pytest

PROBE = Path(__file__).resolve().parents[1] / "benchmarks" / "discus_6_23_01" / "c_probes"


def _load(name):
    spec = importlib.util.spec_from_file_location(f"c_probe_{name}", PROBE / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_engine_logs_match_registered_predictions(capsys):
    analyze = _load("analyze")
    evidence = sorted((PROBE / "evidence").iterdir())[-1]
    assert analyze.main(evidence) == 0
    out = capsys.readouterr().out
    assert "input SHA-256: 6/6 match" in out and out.rstrip().endswith("PASS")
    assert out.count("complete") == 3 and "INCOMPLETE" not in out


@pytest.mark.parametrize("cut", ["corr occ, Cu, Au", "corr disp, all, all"])
def test_analysis_fails_when_a_corr_block_is_missing(tmp_path, cut):
    """A log truncated before a whole corr block (e.g. an earlier crash) must not pass."""
    import shutil
    analyze = _load("analyze")
    evidence = sorted((PROBE / "evidence").iterdir())[-1]
    dst = tmp_path / "ev"
    shutil.copytree(evidence, dst)
    for log in dst.glob("*.log"):
        text = log.read_text(errors="replace")
        if cut in text:
            log.write_text(text[:text.index(cut)])
    assert analyze.main(dst) == 1


def test_current_code_reproduces_the_registered_predictions():
    gen = _load("generate")
    from discuspdf.discus_io import read_stru
    import numpy as np

    preds = json.loads((PROBE / "predictions.json").read_text(encoding="ascii"))

    def labelled(name):
        s = read_stru(PROBE / "inputs" / f"{name}.stru").structure
        m = s.n_atoms // int(np.prod(s.supercell))
        return s.replace(site_label=np.array([f"site{r % m + 1}" for r in range(s.n_atoms)], dtype=object))

    def same(a, b):
        return (math.isnan(a) and math.isnan(b)) or a == pytest.approx(b, abs=1e-12)

    occ = gen.predictions_occ(labelled("occ_bcc"))
    for got, want in zip(occ["neighborhoods"], preds["occ_bcc"]["neighborhoods"]):
        for k in ("corr_AuCu", "corr_CuAu", "n_pairs_discus", "n_pairs_ours"):
            if k in want:
                assert same(got[k], want[k]), k
    disp = gen.predictions_disp(labelled("disp_sc"), gen.DISP_VECS, [("Cu", "Cu")])
    for got, want in zip(disp["neighborhoods"], preds["disp_sc"]["neighborhoods"]):
        assert same(got["Cu_Cu"]["discus_convention"], want["Cu_Cu"]["discus_convention"])
    pairs = [("Au", "Cu"), ("Cu", "Cu"), ("all", "all")]
    mixed = gen.predictions_disp(labelled("mixed_bcc"), gen.MIXED_VECS, pairs)
    for got, want in zip(mixed["neighborhoods"], preds["mixed_bcc"]["neighborhoods"]):
        for k in ("Au_Cu", "Cu_Cu", "all_all"):
            assert same(got[k]["discus_convention"], want[k]["discus_convention"])
            assert same(got[k]["pearson"], want[k]["pearson"])
            assert got[k]["n_pairs"] == want[k]["n_pairs"]
