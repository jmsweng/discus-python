from pathlib import Path
import numpy as np
import pytest
from discuspdf import PDFConfig, PDFRunner, render_pdf_macro, read_csv


def test_render_and_prepare(tmp_path):
    sample = tmp_path / "input.cll"
    sample.write_text("demo input")
    cfg = PDFConfig(sample, supercell=(3, 4, 5), rmax=12, dr=0.03, qmax=20)
    text = render_pdf_macro(cfg)
    assert "cell model.cll,3,4,5" in text
    # v6.x powder->PDF workflow (verified against DISCUS_SUITE 6.23.01)
    assert "powder" in text
    assert "set calc,debye" in text
    assert "set qmax,20" in text
    assert "form pdf, r, 0.0, 12, 0.03" in text
    assert "value PDF" in text
    run = PDFRunner(runs_directory=tmp_path / "runs").prepare(cfg)
    assert (run.directory / "model.cll").read_text() == "demo input"
    assert (run.directory / "config.json").exists()


def test_render_qmax_disabled(tmp_path):
    """qmax=0 disables the Qmax termination correction (no set qmax line)."""
    sample = tmp_path / "input.cll"
    sample.write_text("demo input")
    cfg = PDFConfig(sample, qmax=0.0)
    text = render_pdf_macro(cfg)
    assert "set qmax" not in text
    assert "set delta,0" in text  # qsig=0 -> no Gaussian convolution


def test_validation():
    with pytest.raises(ValueError):
        PDFConfig("input.cll", dr=-1)
    with pytest.raises(ValueError):
        PDFConfig("input.cll", radiation="electron")


def test_numeric_output(tmp_path):
    path = tmp_path / "result.csv"
    path.write_text("# PDF output\nr G(r)\n0.01 1.0 0.0\n0.02 1.5 0.0\n")
    r, g = read_csv(path)
    assert np.allclose(r, [0.01, 0.02])
    assert np.allclose(g, [1.0, 1.5])


def test_runner_end_to_end_with_fake_executable(tmp_path):
    """Tests Python subprocess integration, NOT real Fortran correctness."""
    import sys
    script = tmp_path / "fake_discuss.py"
    script.write_text("from pathlib import Path\nPath('result.csv').write_text('0.0 1.0\\n0.1 1.2\\n')\n")
    cell = tmp_path / "cell.cll"
    cell.write_text("dummy")
    runner = PDFRunner(
        runs_directory=tmp_path / "runs",
        launcher=[sys.executable, str(script)],
    )
    result = runner.run(PDFConfig(cell))
    x, y = result.load()
    assert len(x) == 2
    assert (result.directory / "status.json").read_text().find('"complete"') >= 0
    # Column headings are prepended to the CSV
    content = (result.directory / "result.csv").read_text()
    assert content.startswith("r (Angstrom)\tG(r)")
