"""Unit tests for discuspdf.discus_powder (macro lines, validation, file reading)."""
from __future__ import annotations

import numpy as np
import pytest

from discuspdf.discus_ops import StructureJob
from discuspdf.discus_powder import PowderPattern, read_powder_file
from discuspdf.engine import fork_reasons
from discuspdf.fourier_macro import StructureSource
from discuspdf.suite_catalogue import load_catalogue

SRC = StructureSource.cell("cu.cll", (3, 3, 3))


def text(step):
    t = StructureJob(SRC, [step], readback="none").macro().render()
    assert load_catalogue().check(t) == []
    assert fork_reasons(t) == []
    return t.splitlines()


def test_lines():
    ls = text(PowderPattern((("S(Q)", "s.dat"), ("inte", "i.dat")), radiation="neutron", wavelength=0.7))
    i = ls.index("powder")
    assert ls[i + 1:i + 16] == ["reset", "neut", "set wvle,0.7", "set calc,debye", "set qmin,0.5", "set qmax,15",
                                "set dq,0.01", "set stepr,0.01", "set profile,off", "set lpcor,none",
                                "set temp,ignore", "set disp,off", "set period,off", "run", "exit"]
    assert ls[i + 16:i + 22] == ["output", "outf s.dat", "form powd,q,0.5,15,0.01", "value S(Q)", "run", "exit"]
    assert "value inte" in ls
    e = text(PowderPattern((("F(Q)", "f.dat"),), energy=17.48, temperature="use", polarisation="bragg"))
    assert "set energy,17.48" in e and "set temp,use" in e and "set lpcor,bragg" in e


@pytest.mark.parametrize("kw", [
    dict(outputs=()), dict(outputs=(("PDF", "p.dat"),)), dict(outputs=(("inte", "a"), ("S(Q)", "a"))),
    dict(qmin=2.0, qmax=1.0), dict(dq=0.0123456), dict(qmin=0.505), dict(qmax=20.0),
    dict(radiation="neutron"), dict(wavelength="xx1"), dict(energy=10.0, wavelength="cua1"),
    dict(wavelength=1.5406, qmax=8.2), dict(stepr=0), dict(temperature="hot"), dict(polarisation="tof"),
    dict(dq=0.00001, qmax=15.0),
])
def test_invalid(kw):
    args = {"outputs": (("S(Q)", "s.dat"),), **kw}
    with pytest.raises(ValueError):
        PowderPattern(**args)


def test_read_and_shift(tmp_path):
    f = tmp_path / "s.dat"
    f.write_text(" 0.5  1.0\n 0.51  1.1\n")
    q, y = read_powder_file(f)
    assert q.tolist() == [0.5, 0.51] and y.tolist() == [1.0, 1.1]
    q2, _ = read_powder_file(f, dq=0.01, shift_labels=True)
    assert np.allclose(q2, [0.49, 0.50])
    with pytest.raises(ValueError):
        read_powder_file(f, shift_labels=True)
