"""Live checks of SuiteRunner and Macro on a real engine (opt-in: DISCUS_LIVE=1).

Small on purpose: each case is one short macro on the committed 2-atom GaAs cell.

Environment: DISCUS_LIVE=1; DISCUS_LIVE_WORK (no spaces; default C:\\tmp\\suitelive on
Windows, /tmp/suitelive); DISCUS_WSL_DISTRO (default Ubuntu-22.04). The fork case also needs
DISCUSPDF_ENGINE_FORK (the fork's launcher, e.g. ``wsl.exe -d Ubuntu-22.04 -- /path/to/fork/discus_suite
-noautorun -macro``) and is skipped without it.
"""
from __future__ import annotations

import os
import platform
from pathlib import Path

import numpy as np
import pytest

from discuspdf.discus_io import read_stru
from discuspdf.engine import Engine
from discuspdf.macro import Macro
from discuspdf.suite_runner import SuiteRunner
from tests.live_engine import live_engine

pytestmark = pytest.mark.discus_live
if os.environ.get("DISCUS_LIVE") != "1":
    pytest.skip("live DISCUS test; set DISCUS_LIVE=1 to enable", allow_module_level=True)

ROOT = Path(__file__).resolve().parents[1]
GAAS = ROOT / "tests" / "fixtures" / "discus_6_23_01" / "gaas.cll"
WORK = Path(os.environ.get("DISCUS_LIVE_WORK",
                           r"C:\tmp\suitelive" if platform.system() == "Windows" else "/tmp/suitelive"))


def stock() -> Engine:
    return live_engine(600, "stock")


def runner(**kw) -> SuiteRunner:
    if " " in str(WORK):
        pytest.skip("DISCUS_LIVE_WORK must not contain spaces")
    return SuiteRunner(WORK, stock(), **kw)


def supercell_macro(engine_note: str) -> Macro:
    m = Macro(f"suite live check: GaAs 2x2x2 supercell ({engine_note})")
    with m.program("discus"):
        with m.menu("read"):
            m.command("cell", "gaas.cll", 2, 2, 2)
        m.command("save", "supercell.stru")
    return m


def check_supercell(run):
    # DISCUS 'save' writes no ncell line (known since PR #12), so the crystal size is given
    s = read_stru(run.output("supercell.stru"), supercell=(2, 2, 2)).structure
    assert s.n_atoms == 64  # F-43m: 8 atoms per cell, 2x2x2 cells
    assert sorted({x.upper() for x in s.species}) == ["AS", "GA"]


def test_discus_read_and_save_on_stock():
    run = runner().run(supercell_macro("stock"), {"gaas.cll": GAAS}, ["supercell.stru"])
    assert run.status()["status"] == "complete"
    assert run.record["engine"]["label"] == "stock"
    check_supercell(run)


def test_kuplot_function_and_ksave():
    m = Macro("suite live check: KUPLOT function and ksave")
    with m.program("kuplot"):
        m.command("function", "sin(r[0])", 0, 6, 0.5)
        with m.menu("ksave", 1):
            m.command("outfile", "sine.xy")
            m.command("format", "xy")
            m.command("run")
    run = runner().run(m, (), ["sine.xy"])
    data = np.loadtxt(run.output("sine.xy"))
    assert data.shape[0] == 13
    assert np.allclose(data[:, 1], np.sin(data[:, 0]), atol=1e-5)


def test_engine_error_marks_the_run_failed():
    m = Macro("suite live check: missing input file must fail")
    with m.program("discus"):
        with m.menu("read"):
            m.command("cell", "does_not_exist.cll", 1, 1, 1)
        m.command("save", "never.stru")
    r = runner()
    run = r.prepare(m, (), ["never.stru"])
    with pytest.raises(RuntimeError):
        r.execute(run)
    assert run.status()["status"] == "failed"


@pytest.mark.skipif(not os.environ.get("DISCUSPDF_ENGINE_FORK"), reason="no fork engine configured")
def test_same_macro_on_the_fork():
    fork = Engine.from_env("fork", timeout_seconds=600)
    run = runner(fork=fork).run(supercell_macro("fork"), {"gaas.cll": GAAS}, ["supercell.stru"], engine="fork")
    assert run.record["engine"]["label"] == "fork"
    check_supercell(run)
