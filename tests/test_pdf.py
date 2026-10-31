from pathlib import Path
import json
import sys

import numpy as np
import pytest
import warnings
from discuspdf import PDFConfig, PDFRunner, render_pdf_macro, read_pdf, QmaxLimitWarning


def test_render_and_prepare(tmp_path):
    sample = tmp_path / "input.cll"
    sample.write_text("demo input")
    cfg = PDFConfig(sample, supercell=(3, 4, 5), rmax=12, dr=.03,
                    qmin=.8, qmax=20, dq=.02, wavelength=1.7)
    text = render_pdf_macro(cfg)
    for command in ("cell model.cll,3,4,5", "set calc,debye",
                    "set qmin,0.8", "set qmax,20", "set dq,0.02",
                    "set wvle,1.7", "form pdf, r, 0.0, 12, 0.03"):
        assert command in text
    run = PDFRunner(runs_directory=tmp_path / "runs").prepare(cfg)
    assert (run.directory / "model.cll").read_text() == "demo input"
    assert json.loads((run.directory / "status.json").read_text())["status"] == "prepared"


@pytest.mark.parametrize("changes", [
    {"qmax": 0}, {"qmin": 20, "qmax": 10}, {"dq": 0}, {"qsig": -1},
    {"dr": -1}, {"rmax": 1.03, "dr": 0.2}, {"wavelength": -1},
    {"radiation": "electron"}
])
def test_invalid_config(changes):
    with pytest.raises(ValueError):
        PDFConfig("input.cll", **changes)


def test_numeric_output(tmp_path):
    path = tmp_path / "result.pdf"
    path.write_text("# PDF output\nr G(r)\n0.0 1.0 0.0\n0.1 1.5 0.0\n")
    r, g = read_pdf(path)
    assert np.allclose(r, [0, .1])
    assert np.allclose(g, [1, 1.5])


def fake_runner(tmp_path, output):
    script = tmp_path / "fake_discuss.py"
    script.write_text("from pathlib import Path\nPath('result.pdf').write_text(" + repr(output) + ")\n")
    cell = tmp_path / "cell.cll"
    cell.write_text("dummy")
    runner = PDFRunner(runs_directory=tmp_path / "runs", launcher=[sys.executable, str(script)])
    cfg = PDFConfig(cell, rmax=.2, dr=.1)
    return runner, cfg


def test_execute_prepared_exact_macro(tmp_path):
    runner, cfg = fake_runner(tmp_path, "0.0 1.0\n0.1 1.2\n0.2 0.8\n")
    prepared = runner.prepare(cfg)
    macro_before = (prepared.directory / "run.mac").read_bytes()
    result = runner.execute(prepared)
    assert result.directory == prepared.directory
    assert (result.directory / "run.mac").read_bytes() == macro_before
    assert len(result.load()[0]) == 3
    assert json.loads((result.directory / "status.json").read_text())["status"] == "complete"
    with pytest.raises(FileExistsError):
        runner.execute(prepared)


@pytest.mark.parametrize("output", [
    "0.0 1\n0.1 2\n",                    # truncated
    "0.0 1\n0.12 2\n0.2 3\n",         # wrong grid
    "0.0 1\n0.1 2\n0.2 NaN\n",        # bad value
    "garbage\n",
])
def test_reject_bad_pdf(tmp_path, output):
    runner, cfg = fake_runner(tmp_path, output)
    prepared = runner.prepare(cfg)
    with pytest.raises(RuntimeError):
        runner.execute(prepared)
    assert json.loads((prepared.directory / "status.json").read_text())["status"] == "failed"


def test_reject_staged_input_tampering(tmp_path):
    runner, cfg = fake_runner(tmp_path, "0.0 1\n0.1 2\n0.2 3\n")
    prepared = runner.prepare(cfg)
    (prepared.directory / "model.cll").write_text("modified")
    with pytest.raises(ValueError, match="hash"):
        runner.execute(prepared)
    # E3: the failed preflight is recorded; the directory does not stay "prepared".
    st = json.loads((prepared.directory / "status.json").read_text())
    assert st["status"] == "failed"
    assert st["engine_started"] is False
    assert "returncode" not in st
    assert not (prepared.directory / "run.log").exists()


def test_reject_tampered_macro(tmp_path):
    runner, cfg = fake_runner(tmp_path, "0.0 1\n0.1 2\n0.2 3\n")
    prepared = runner.prepare(cfg)
    (prepared.directory / "run.mac").write_text("tampered macro")
    with pytest.raises(ValueError, match="macro"):
        runner.execute(prepared)
    st = json.loads((prepared.directory / "status.json").read_text())
    assert st["status"] == "failed"
    assert st["engine_started"] is False


def test_rerun_rejected_with_failed_status(tmp_path):
    runner, cfg = fake_runner(tmp_path, "0.0 1\n0.1 2\n0.2 3\n")
    prepared = runner.prepare(cfg)
    runner.execute(prepared)
    assert json.loads((prepared.directory / "status.json").read_text())["status"] == "complete"
    with pytest.raises(FileExistsError):
        runner.execute(prepared)
    st = json.loads((prepared.directory / "status.json").read_text())
    assert st["status"] == "failed"
    assert st["engine_started"] is False


def test_missing_binary_records_failed_without_engine(tmp_path):
    cell = tmp_path / "cell.cll"
    cell.write_text("dummy")
    runner = PDFRunner(runs_directory=tmp_path / "runs",
                       launcher=["definitely_not_a_binary_xyz"])
    prepared = runner.prepare(PDFConfig(cell, rmax=.2, dr=.1))
    with pytest.raises(RuntimeError):
        runner.execute(prepared)
    st = json.loads((prepared.directory / "status.json").read_text())
    assert st["status"] == "failed"
    assert st["engine_started"] is False
    assert "returncode" not in st


def test_prepare_records_macro_hash(tmp_path):
    cell = tmp_path / "cell.cll"
    cell.write_text("dummy")
    runner = PDFRunner(runs_directory=tmp_path / "runs",
                       launcher=[sys.executable, str(tmp_path / "unused.py")])
    prepared = runner.prepare(PDFConfig(cell, rmax=.2, dr=.1))
    import hashlib
    params = json.loads((prepared.directory / "config.json").read_text())
    assert params["macro_sha256"] == hashlib.sha256(
        (prepared.directory / "run.mac").read_bytes()).hexdigest()


@pytest.mark.parametrize("changes", [
    {}, {"qmax": 30}, {"qmax": 8.0}, {"radiation": "neutron"}, {"qmax": 20, "wavelength": 1.7},
])
def test_qmax_limit_warning(changes):
    with pytest.warns(QmaxLimitWarning, match="calculation runs unchanged"):
        PDFConfig("input.cll", **changes)


@pytest.mark.parametrize("changes", [
    {"qmax": 7.5}, {"qmax": 25, "wavelength": 0.4}, {"qmax": 12, "wavelength": 1.0},
])
def test_no_qmax_limit_warning(changes):
    with warnings.catch_warnings():
        warnings.simplefilter("error", QmaxLimitWarning)
        PDFConfig("input.cll", **changes)


def test_macro_unchanged_above_limit():
    with pytest.warns(QmaxLimitWarning):
        cfg = PDFConfig("input.cll", qmax=25)
    assert render_pdf_macro(cfg).splitlines() == [
        "# Generated by discuspdf: 1D atomic PDF calculation (DISCUS v6.x)",
        "discus", "read", "cell model.cll,5,5,5", "powder", "reset", "xray",
        "set axis,q", "set calc,debye", "set qmin,0.5", "set qmax,25", "set dq,0.01",
        "set delta,0", "set wvle,1.5406", "run", "exit", "output", "outf result.pdf",
        "form pdf, r, 0.0, 15, 0.02", "value PDF", "run", "exit", "exit",
    ]


def test_execute_does_not_warn_again(tmp_path):
    with pytest.warns(QmaxLimitWarning):
        runner, cfg = fake_runner(tmp_path, "0.0 1.0\n0.1 1.2\n0.2 0.8\n")
        prepared = runner.prepare(cfg)
    with warnings.catch_warnings():
        warnings.simplefilter("error", QmaxLimitWarning)
        runner.execute(prepared)
