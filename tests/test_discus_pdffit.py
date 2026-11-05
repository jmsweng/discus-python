"""PdfData, PdfFit, PdfPartials: validation and the macro (no engine)."""
from __future__ import annotations

import numpy as np
import pytest

from discuspdf.discus_ops import ReplaceRandom, StructureJob
from discuspdf.discus_pdfmenu import PdfCalc, PdfData, PdfFit, PdfPartials, read_pdf_partials
from discuspdf.discus_rmc import RmcShift, RmcSwap
from discuspdf.fourier_macro import Seed, StructureSource

SRC = StructureSource.cell("cu_b05.cll", (3, 3, 3))
SEED = Seed((1, 2, 3, 4, 5, 6, 7, 8))


@pytest.fixture()
def data(tmp_path):
    r = np.arange(1, 801) * 0.01
    return PdfData.write(tmp_path / "t.dat", r, np.sin(r), sigma=0.1), tmp_path / "t.dat"


def test_pdfdata_write_and_read(data):
    d, path = data
    assert (d.file, d.dr, d.rmax, d.points) == ("t.dat", pytest.approx(0.01), pytest.approx(8.0), 800)
    back = PdfData.from_file(path)
    assert (back.file, back.points) == ("t.dat", 800) and abs(back.dr - 0.01) < 1e-9 and abs(back.rmax - 8) < 1e-9
    first = path.read_text().splitlines()[0].split()
    assert first[0] == "0.010000" and first[2] == "0.0" and first[3] == "0.1"


def test_pdfdata_needs_first_r_equal_step(tmp_path):
    with pytest.raises(ValueError, match="first r = dr"):
        PdfData.write(tmp_path / "x.dat", np.arange(0, 800) * 0.01, np.zeros(800))
    with pytest.raises(ValueError, match="sigma"):
        PdfData.write(tmp_path / "x.dat", np.arange(1, 9) * 0.01, np.zeros(8), sigma=-1)


def fit(data, **kw):
    args = dict(model=PdfCalc("fit.dat", radiation="xray", thermal=True), data=data, moves=(RmcSwap(),),
                select=("CU", "AU"), tries=100, sigma=0.0)
    args.update(kw)
    return PdfFit(**args)


def test_pdffit_macro_has_the_working_order(data):
    lines = StructureJob(SRC, [ReplaceRandom("CU", "AU", 0.25), fit(data[0])], seed=SEED).macro().render().splitlines()
    i = {k: lines.index(k) for k in ("set range,9,0.01", "calc", "data t.dat", "run", "save pdf,fit.dat")}
    assert i["set range,9,0.01"] < i["calc"] < i["data t.dat"] < i["run"] < i["save pdf,fit.dat"]
    assert "set frange,0.01,7.5" in lines and "set mode,0,shift,all" in lines and "set mode,1,swchem,all" in lines
    assert lines.index("dese all") < lines.index("sele CU,AU")
    assert "set scal,on,1" in lines and "set sigm,0" in lines and "set cycl,100" in lines


def test_pdffit_refusals(data):
    d = data[0]
    with pytest.raises(ValueError, match="needs thermal"):
        fit(d, model=PdfCalc("fit.dat", radiation="xray"))
    with pytest.raises(ValueError, match="fit_range"):
        fit(d, fit_range=(0.01, 9.0))
    with pytest.raises(ValueError, match="sigma for exactly"):
        fit(d, moves=(RmcShift({"CU": (0.01, 0.01, 0.01)}),))
    with pytest.raises(ValueError, match="two selected"):
        fit(d, select=("CU",))
    with pytest.raises(ValueError, match="differ from the data"):
        fit(d, model=PdfCalc("t.dat", radiation="xray", thermal=True))
    with pytest.raises(TypeError, match="PdfData"):
        fit("t.dat")


def test_pdffit_checks_the_file(data, tmp_path):
    d, path = data
    other = PdfData.write(tmp_path / "t2.dat", np.arange(1, 401) * 0.02, np.zeros(400))
    step = fit(d)
    step.check_inputs({"t.dat": path})
    with pytest.raises(ValueError, match="the step expects"):
        step.check_inputs({"t.dat": tmp_path / other.file})


def test_partials_selections_and_files():
    p = PdfPartials(PdfCalc("x.dat", rmax=8, radiation="neutron"), (("Cu", "Au"), ("AU", "AU"), ("AU", "AG")))
    assert p.pairs == (("AU", "CU"), ("AU", "AU"), ("AG", "AU"))
    assert p.selections() == [("AU",), ("CU",), ("AU", "CU"), ("AG",), ("AG", "AU")]
    lines = StructureJob(SRC, [p], readback="none").macro().render().splitlines()
    assert lines.count("calc") == 5 and lines[-5:-3] == ["isele all", "jsele all"]
    assert "jsele AU,CU" in lines and "save pdf,partial_AU_CU.dat" in lines
    with pytest.raises(ValueError, match="twice"):
        PdfPartials(PdfCalc("x.dat", radiation="neutron"), (("AU", "CU"), ("CU", "AU")))


def test_read_partials_combines_the_cross_term(tmp_path):
    p = PdfPartials(PdfCalc("x.dat", rmax=8, radiation="neutron"), (("AU", "CU"),))
    r = np.arange(1, 6) * 0.1
    vals = {("AU",): 1.0, ("CU",): 2.0, ("AU", "CU"): 10.0}

    class Run:
        def output(self, name):
            return tmp_path / name
    for sel, v in vals.items():
        (tmp_path / p.file_of(sel)).write_text("".join(f"{x:.4f} {v * x:.6f} 0.0 1.0\n" for x in r))
    (rr, g), = read_pdf_partials(p, Run()).values()
    assert np.allclose(g, 7.0 * r)
