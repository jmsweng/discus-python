"""Unit tests for discuspdf.discus_pdfmenu (macro lines, validation, file reading)."""
from __future__ import annotations

import pytest

from discuspdf.discus_ops import StructureJob
from discuspdf.discus_pdfmenu import PdfCalc, read_pdf_menu_file
from discuspdf.engine import fork_reasons
from discuspdf.fourier_macro import StructureSource
from discuspdf.suite_catalogue import load_catalogue

SRC = StructureSource.cell("cu.cll", (4, 4, 4))


def lines(step):
    text = StructureJob(SRC, [step], readback="none").macro().render()
    assert load_catalogue().check(text) == []
    assert fork_reasons(text) == []
    ls = text.splitlines()
    return ls[ls.index("pdf") + 1:ls.index("save pdf,g.dat") + 1]


def test_default_lines():
    assert lines(PdfCalc("g.dat")) == [
        "set radiation,xray,0", "set boundary,periodic,3D", "set range,20,0.001", "set qmax,0",
        "set finite,periodic", "set density,auto", "set therm,off", "set delta,0", "set gamma,0", "set qalpha,0",
        "calc", "save pdf,g.dat"]


def test_options_lines():
    ls = lines(PdfCalc("g.dat", rmax=10, dr=0.02, radiation="neutron", qmax=25, periodic=False, exact=True,
                       sphere=30.0, density=0.08, thermal=True, delta=1.5, qalpha=0.01))
    for line in ("set radiation,neutron", "set boundary,crystal,exact", "set range,10,0.02", "set qmax,25",
                 "set finite,sphere,30", "set density,0.08", "set therm,gauss", "set delta,1.5", "set qalpha,0.01"):
        assert line in ls, line
    assert "set radiation,electron,1.5" in lines(PdfCalc("g.dat", radiation="electron", weight_s=1.5))


@pytest.mark.parametrize("kw", [dict(output="a b"), dict(rmax=0), dict(dr=30.0), dict(radiation="gamma"),
                                dict(radiation="neutron", weight_s=1.0), dict(dr=0.01), dict(qmax=-1), dict(sphere=0),
                                dict(density=-0.1), dict(delta=1.0), dict(periodic=1)])
def test_invalid(kw):
    with pytest.raises((ValueError, TypeError)):
        PdfCalc(**{"output": "g.dat", **kw})


def test_output_is_declared_and_file_is_read(tmp_path):
    assert PdfCalc("g.dat").output_names() == ("g.dat",)
    f = tmp_path / "g.dat"
    f.write_text("   0.0100      -0.1234567890     0.00   1.00\n   0.0200       2.5000000000     0.00   1.00\n")
    r, g = read_pdf_menu_file(f)
    assert r.tolist() == [0.01, 0.02] and g.tolist() == [-0.123456789, 2.5]


def test_qalpha_alone_and_coarse_dr_with_smoothing():
    assert "set qalpha,0.02" in lines(PdfCalc("g.dat", qalpha=0.02, dr=0.01))
    assert "set range,20,0.01" in lines(PdfCalc("g.dat", dr=0.01, qmax=25))
