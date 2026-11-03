"""Unit tests for discuspdf.discus_rmc: macro lines, validation, NIPL files, log parsing, routing."""
from __future__ import annotations

import sys

import numpy as np
import pytest

from discuspdf.discus_ops import ReplaceRandom, StructureJob
from discuspdf.discus_rmc import (FourierPlane, PlaneGrid, RmcFit, RmcLots, RmcPlane, RmcShift, RmcSwap,
                                  parse_rmc_log, read_nipl, run_rmc, write_nipl)
from discuspdf.engine import Engine, fork_reasons
from discuspdf.fourier_macro import Seed, StructureSource
from discuspdf.suite_catalogue import load_catalogue
from discuspdf.suite_runner import SuiteRunner

SEED = Seed((1, 2, 3, 4, 5, 6, 7, 8))
GRID = PlaneGrid((0, -1, 0), (0, 1, 0), (3, -1, 0), 21, 31)
SRC = StructureSource.cell("cu.cll", (4, 4, 4))


def fit(**kw):
    args = dict(planes=(RmcPlane("t.nipl", GRID),), moves=(RmcSwap(),), select=("CU", "AU"), tries=2000)
    args.update(kw)
    return RmcFit(**args)


def rmc_lines(step):
    text = StructureJob(SRC, [step], seed=SEED).macro().render()
    assert load_catalogue().check(text) == []
    lines = text.splitlines()
    return lines[lines.index("rmc") + 1:lines.index("run") + 1]


def test_default_fit_lines():
    assert rmc_lines(fit()) == [
        "reset", "data xray,moa1,t.nipl,one,0,-1,0,0,1,0,3,-1,0", "dese all", "sele CU,AU",
        "set mode,0,shift,all", "set mode,1,swchem,all", "set scal,on,1", "set back,off,0", "set sigm,0.01", "set cycl,2000",
        "set disp,500", "run"]


def test_all_options_lines():
    planes = (RmcPlane("a.nipl", GRID, radiation="neutron", wavelength=1.5, weights="sqrt"),)
    f = fit(planes=planes, moves=(RmcSwap(0.7, "site"), RmcShift({"cu": (0.01, 0.02, 0.0), "AU": (0.1, 0, 0)}, 0.3)),
            sigma=0.0, display=100, fit_scale=False, scale=2.0, fit_background=True, background=1.5,
            min_distances={("Cu", "Au"): 2.0}, debye_waller=True, symmetry=True)
    lines = rmc_lines(f)
    assert lines[:2] == ["reset", "data neutron,1.5,a.nipl,sqrt,0,-1,0,0,1,0,3,-1,0"]
    assert not any(line.startswith("set mode,0,") for line in lines)
    for line in ("set dbw,on", "set sym,on", "set mode,0.7,swchem,site", "set mode,0.3,shift,all",
                 "set move,CU,0.01,0.02,0", "set move,AU,0.1,0,0", "set mdis,CU,AU,2", "set scal,off,2",
                 "set back,on,1.5", "set sigm,0", "set disp,100"):
        assert line in lines, line
    assert fork_reasons("\n".join(["discus", "rmc", *lines, "exit", "exit"])) == []


def test_aver_and_lots_go_to_the_fork():
    for f in (fit(average_percent=10.0), fit(lots=RmcLots((2, 2, 2), 20))):
        text = StructureJob(SRC, [f], seed=SEED).macro().render()
        assert load_catalogue().check(text) == []
        assert fork_reasons(text), text
    assert fork_reasons("discus\nrmc\nset aver,0\nexit\nexit\n") == []
    assert fork_reasons("discus\nrmc\nset av,5\nexit\nexit\n")
    assert "set lots,ellipsoid,2,3,4,5,no" in StructureJob(
        SRC, [fit(lots=RmcLots((2, 3, 4), 5, "ellipsoid", False))], seed=SEED).macro().render()


@pytest.mark.parametrize("make", [
    lambda: fit(planes=()),
    lambda: fit(planes=(RmcPlane("t.nipl", GRID), RmcPlane("u.nipl", GRID))),  # one plane only
    lambda: fit(moves=()),
    lambda: fit(moves=(RmcSwap(), RmcSwap(0.5))),
    lambda: fit(select=("CU",)),  # swaps need two types
    lambda: fit(select="most"),
    lambda: fit(moves=(RmcShift({"NI": (0.1, 0.1, 0.1)}),)),  # sigma for an unselected type
    lambda: fit(moves=(RmcShift({"CU": (0.1, 0.1, 0.1)}),)),  # no sigma for AU (DISCUS default 0.2)
    lambda: fit(select="all", moves=(RmcShift({"CU": (0.1, 0.1, 0.1)}),)),
    lambda: PlaneGrid((0, 0, 0), (3, 0, 0), (0, 1, 0), 31, 21),  # nx > ny
    lambda: RmcPlane("t.nipl", GRID, weights="log"),
    lambda: fit(tries=0),
    lambda: fit(sigma=-1.0),
    lambda: fit(scale=0.0),
    lambda: fit(background=-1.0),
    lambda: fit(average_percent=120),
    lambda: fit(fit_scale=1),
    lambda: RmcPlane("t.nipl", GRID, weights="file.dat"),
    lambda: RmcPlane("t.nipl", GRID, radiation="electron"),
    lambda: RmcPlane("t.nipl", GRID, wavelength="xx"),
    lambda: RmcPlane("../t.nipl", GRID),
    lambda: PlaneGrid((0, 0, 0), (1, 0, 0), (2, 0, 0), 5, 5),  # degenerate
    lambda: PlaneGrid((0, 0, 0), (1, 0, 0), (0, 1, 0), 1, 5),
    lambda: RmcSwap(pairing="conn"),
    lambda: RmcShift({"CU": (0, 0, 0)}),
    lambda: RmcShift({}),
    lambda: RmcLots((2, 2), 5),
    lambda: RmcLots((2, 2, 2), 5, shape="sphere"),
])
def test_invalid(make):
    with pytest.raises((ValueError, TypeError)):
        make()


def test_fit_needs_a_seed_and_lists_its_files():
    with pytest.raises(ValueError, match="seed"):
        StructureJob(SRC, [fit()])
    assert fit().input_names() == ("t.nipl",)


def test_plane_hkl_follows_the_rmc_data_convention():
    h = GRID.hkl()
    assert h.shape == (31, 21, 3)
    assert np.allclose(h[0, 0], (0, -1, 0)) and np.allclose(h[0, -1], (0, 1, 0))
    assert np.allclose(h[-1, 0], (3, -1, 0)) and np.allclose(h[20, 10], (2, 0, 0))


def test_nipl_round_trip_and_discus_header(tmp_path):
    z = np.arange(6.0).reshape(2, 3)
    p = write_nipl(tmp_path / "a.nipl", z, (0, 1, 2, 3))
    back, r = read_nipl(p)
    assert np.array_equal(back, z) and r == (0, 1, 2, 3)
    # DISCUS writes '#' header lines and a blank line after each row
    q = tmp_path / "b.nipl"
    q.write_text("# HEADER  : lines\n# x\n    3    2\n 0. 1. 2. 3.\n  1 2 3\n\n  4 5 6\n\n")
    back, _ = read_nipl(q)
    assert back.tolist() == [[1, 2, 3], [4, 5, 6]]
    with pytest.raises(ValueError, match="points"):
        RmcPlane("b.nipl", PlaneGrid((0, 0, 0), (1, 0, 0), (0, 1, 0), 2, 3)).check_file(q)


def test_fourier_plane_lines():
    job = StructureJob(SRC, [ReplaceRandom("CU", "AU", 1.0, site=1), FourierPlane(GRID, "t.nipl")],
                       seed=SEED, readback="none")
    text = job.macro().render()
    assert load_catalogue().check(text) == []
    for line in ("ll 0,-1,0", "lr 0,1,0", "ul 3,-1,0", "tl 0,-1,0", "na 21", "no 31", "nt 1",
                 "temp ignore", "value int", "form stan", "outf t.nipl"):
        assert line in text.splitlines(), line
    assert fork_reasons(text) == []


LOG = """ Running RMC-Fit ...

 Size of model crystal     :   4 x   4 x   4 containing       256 atoms

 Starting main RMC loop ...
 (Mode SWCHEM)

 Gen:       500 try:       480 acc: (good/bad):     40 /      0 s2x2:  0.12345678E+08
   Plane  1: scal:  0.1000E-01 / back:   0.000         /        s2x2:  0.12345678E+08

 ---- Final configuration ----

 Gen:      2100 try:      2000 acc: (good/bad): ****** /      3 s2x2:   1234.5678
   Plane  1: scal:  0.2000E-01 / back:   1.500         /        s2x2:   1234.5678

 Delta Chi    : ave:      12.50     sig:      3.000     max:      20.00
"""


def test_parse_rmc_log():
    log = parse_rmc_log(LOG)
    assert len(log.progress) == 2
    first, last = log.progress
    assert (first.generated, first.tries, first.accepted_downhill, first.accepted_uphill) == (500, 480, 40, 0)
    assert first.chi2 == pytest.approx(1.2345678e7) and not first.final
    assert last.final and log.final is last and last.accepted_downhill is None
    assert last.planes == ((pytest.approx(0.02), pytest.approx(1.5), pytest.approx(1234.5678)),)
    assert log.delta_chi2 == (12.5, 3.0, 20.0)
    assert parse_rmc_log("nothing").final is None


FAKE = """
import shutil, sys
from pathlib import Path
text = Path("cu.cll").read_text()
lines = [l for l in text.splitlines()]
Path("result.stru").write_text(text)
print('''LOGTEXT''')
"""


def test_run_rmc_checks_planes_and_parses(tmp_path):
    cell = tmp_path / "cu.cll"
    cell.write_text("title cu\nspcgr P1\ncell 3.615,3.615,3.615,90,90,90\natoms\nCU 0 0 0 0.5\n")
    nipl = write_nipl(tmp_path / "t.nipl", np.ones((31, 21)))
    script = tmp_path / "fake.py"
    script.write_text(FAKE.replace("LOGTEXT", LOG))
    runner = SuiteRunner(tmp_path / "runs", Engine((sys.executable, str(script))))
    src = StructureSource.cell("cu.cll", (1, 1, 1))
    job = StructureJob(src, [fit()], seed=SEED, readback="none")
    out = run_rmc(job, runner, {"cu.cll": cell, "t.nipl": nipl})
    assert out.log.final.tries == 2000
    bad = write_nipl(tmp_path / "bad.nipl", np.ones((21, 31)))
    with pytest.raises(ValueError, match="points"):
        run_rmc(job, runner, {"cu.cll": cell, "t.nipl": bad})
    with pytest.raises(ValueError, match="inputs"):
        run_rmc(job, runner, {"cu.cll": cell})


def test_fourier_plane_file_is_a_declared_output(tmp_path):
    cell = tmp_path / "cu.cll"
    cell.write_text("title cu\nspcgr P1\ncell 3.615,3.615,3.615,90,90,90\natoms\nCU 0 0 0 0.5\n")
    script = tmp_path / "fake.py"
    script.write_text("from pathlib import Path\nPath('result.stru').write_text(Path('cu.cll').read_text())\n"
                      "Path('t.nipl').write_text('2 2\\n0 1 0 1\\n1 2\\n3 4\\n')\n")
    runner = SuiteRunner(tmp_path / "runs", Engine((sys.executable, str(script))))
    job = StructureJob(StructureSource.cell("cu.cll", (1, 1, 1)), [FourierPlane(GRID, "t.nipl")], readback="none")
    r = job.run(runner, {"cu.cll": cell})
    assert read_nipl(r.run.output("t.nipl"))[0].tolist() == [[1, 2], [3, 4]]
    with pytest.raises(ValueError, match="twice"):
        StructureJob(SRC, [FourierPlane(GRID, "t.nipl"), FourierPlane(GRID, "t.nipl")], readback="none").run(
            runner, {"cu.cll": cell})
