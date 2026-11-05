"""PattersonMap, InverseMap, RealPlane, read_map (no engine)."""
from __future__ import annotations

import numpy as np
import pytest

from discuspdf.discus_ops import StructureJob
from discuspdf.discus_patterson import InverseMap, PattersonMap, RealPlane, read_map, run_map, strip_nipl_header
from discuspdf.discus_rmc import FourierPlane, PlaneGrid
from discuspdf.fourier_macro import StructureSource

SRC = StructureSource.cell("cu.cll", (2, 2, 2))
GRID = PlaneGrid(ll=(-3, -3, 0), lr=(3, -3, 0), ul=(-3, 3, 0), nx=13, ny=13)
PLANE = RealPlane((0, 0, 0), (2, 0, 0), (0, 2, 0), 21, 21)


def test_extended_plane_keeps_the_points():
    ext = PLANE.extended()
    assert ext.no == 22 and ext.ll == (0.0, -0.1, 0.0) and ext.lr == (2.0, -0.1, 0.0) and ext.ul == (0.0, 2.0, 0.0)
    assert np.allclose(ext.points()[1:], PLANE.points())


def test_patterson_macro():
    lines = StructureJob(SRC, [PattersonMap("i.nipl", GRID, PLANE, "p.nipl")], readback="none").macro().render().splitlines()
    i = lines.index("patterson")
    assert lines[i + 1:i + 5] == ["reset", "format stan", "file a,i.nipl", "type a,intensity"]
    assert "rholl 0,-0.1,0" in lines and "rhono 22" in lines and "rhona 21" in lines
    j = lines.index("output")
    assert lines[j + 1:j + 5] == ["value real", "form stan", "outf p.nipl", "run"]


def test_inverse_macro_and_outputs():
    step = InverseMap("re.nipl", "im.nipl", GRID, PLANE, "rho.nipl", imag_output="rhoi.nipl")
    lines = StructureJob(SRC, [step], readback="none").macro().render().splitlines()
    assert "type a,real" in lines and "type b,imag" in lines and lines.count("output") == 2
    assert step.output_names() == ("rho.nipl", "rhoi.nipl") and step.input_names() == ("re.nipl", "im.nipl")
    with pytest.raises(ValueError, match="differ"):
        InverseMap("re.nipl", "re.nipl", GRID, PLANE, "rho.nipl")


def test_fourier_plane_values():
    lines = StructureJob(SRC, [FourierPlane(GRID, "re.nipl", value="real")], readback="none").macro().render()
    assert "value real" in lines.splitlines()
    with pytest.raises(ValueError, match="value"):
        FourierPlane(GRID, "x.nipl", value="phase")


def test_strip_and_read_map(tmp_path):
    z = np.arange(22 * 21, dtype=float).reshape(22, 21)
    text = "# HEADER\n# more\n 21 22\n 0 2 -0.1 2\n" + "\n".join(" ".join(f"{v:.5E}" for v in row) for row in z)
    (tmp_path / "m.nipl").write_text(text)
    strip_nipl_header(tmp_path / "m.nipl", tmp_path / "s.nipl")
    assert not (tmp_path / "s.nipl").read_text().startswith("#")
    assert np.array_equal(read_map(tmp_path / "m.nipl", PLANE), z[1:])
    with pytest.raises(ValueError, match="expected"):
        read_map(tmp_path / "m.nipl", RealPlane((0, 0, 0), (1, 0, 0), (0, 1, 0), 21, 23))


def test_run_map_copies_inputs_without_headers(tmp_path):
    (tmp_path / "i.nipl").write_text("# h\n 13 13\n")
    seen = {}

    class Job:
        steps = [PattersonMap("i.nipl", GRID, PLANE, "p.nipl")]

        def run(self, runner, paths):
            seen.update(paths)
            return "ok"
    assert run_map(Job(), None, {"i.nipl": tmp_path / "i.nipl", "cu.cll": "x"}, tmp_path / "w") == "ok"
    assert seen["cu.cll"] == "x" and (tmp_path / "w" / "i.nipl").read_text() == " 13 13\n"
