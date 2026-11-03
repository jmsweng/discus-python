"""The 1D PDF output of a PDFRunner run is ``result.csv``.

This is a deliberate decision by the repository owner (2026-11-02). It has been reverted
to ``result.pdf`` before (merge f8b2ba1, 2026-10-22) by an agent that assumed DISCUS
requires that name. It does not: DISCUS writes whatever ``outf`` names. If a test here
fails because the output was renamed back to .pdf, restore result.csv instead of editing
the test. See CLAUDE.md / AGENTS.md, "1D PDF output is result.csv".
"""
import json
from pathlib import Path
import sys

import numpy as np
import pytest

from discuspdf import (PDF_CSV_NAME, PDFConfig, PDFRun, PDFRunner, read_pdf, read_pdf_csv,
                       render_pdf_macro, write_pdf_csv)

RAW = "0.0 0.0\n0.1 1.2345678901234567\n0.2 -0.8e-3\n"


def _runner(tmp_path, output=RAW):
    script = tmp_path / "fake_discus.py"
    script.write_text("from pathlib import Path\nPath('result.pdf').write_text(" + repr(output) + ")\n")
    cell = tmp_path / "cell.cll"
    cell.write_text("dummy")
    runner = PDFRunner(runs_directory=tmp_path / "runs", launcher=[sys.executable, str(script)])
    return runner, PDFConfig(cell, rmax=.2, dr=.1)


def test_output_is_csv_not_pdf():
    # Owner's decision: do NOT change this back to result.pdf.
    assert PDF_CSV_NAME == "result.csv"
    run = PDFRun(directory=Path("x"))
    assert run.output.name == "result.csv"
    assert run.raw_output.name == "result.pdf"


def test_run_writes_csv_with_header_and_exact_values(tmp_path):
    runner, cfg = _runner(tmp_path)
    run = runner.run(cfg)
    text = run.output.read_text(encoding="utf-8")
    lines = text.splitlines()
    assert lines[0] == "r (Angstrom),G(r)"
    assert all(line.count(",") == 1 for line in lines)
    assert len(lines) == 4
    r_raw, g_raw = read_pdf(run.raw_output)
    r, g = run.load()
    assert np.array_equal(r, r_raw) and np.array_equal(g, g_raw)   # bit-identical to DISCUS's numbers
    status = json.loads((run.directory / "status.json").read_text())
    assert status["status"] == "complete" and status["output"] == "result.csv"


def test_csv_reads_with_plain_numpy(tmp_path):
    path = write_pdf_csv(tmp_path / "x.csv", [0.0, 0.1], [1.0, 2.5])
    arr = np.loadtxt(path, delimiter=",", skiprows=1)
    assert np.array_equal(arr, [[0.0, 1.0], [0.1, 2.5]])


def test_refuses_to_overwrite_existing_csv(tmp_path):
    runner, cfg = _runner(tmp_path)
    prepared = runner.prepare(cfg)
    prepared.output.write_text("old")
    with pytest.raises(FileExistsError):
        runner.execute(prepared)


def test_invalid_engine_output_writes_no_csv(tmp_path):
    runner, cfg = _runner(tmp_path, output="0.0 1.0\n0.1 1.0\n")   # r grid too short
    with pytest.raises(RuntimeError):
        runner.run(cfg)
    run_dir = next((tmp_path / "runs").iterdir())
    assert not (run_dir / "result.csv").exists()


def test_read_pdf_csv_checks_header(tmp_path):
    path = tmp_path / "bad.csv"
    path.write_text("0.0,1.0\n0.1,2.0\n")
    with pytest.raises(ValueError, match="header"):
        read_pdf_csv(path)


def test_macro_still_writes_discus_raw_file():
    # The validated macro is unchanged: archived runs and probe evidence use it byte for byte.
    assert "outf result.pdf\n" in render_pdf_macro(PDFConfig("input.cll"))
