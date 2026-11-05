"""ExpPdf: input conversion, composition, limits and the macro (no engine)."""
from __future__ import annotations

import numpy as np
import pytest

from discuspdf.discus_exp2pdf import (ExpPdf, _comp_stru, read_exp2pdf_file, read_pattern, two_theta_to_q)
from discuspdf.suite_catalogue import load_catalogue


@pytest.fixture()
def data(tmp_path):
    q = np.linspace(0.5, 17.9, 400)
    p = tmp_path / "d.xy"
    p.write_text("# a header line\nQ I\n" + "\n".join(f"{a:.6f}, {1 + a:.6f}" for a in q) + "\n")
    return p, q


def test_two_theta_to_q_round_trips():
    q = np.array([1.0, 5.0, 12.0])
    tth = 2 * np.degrees(np.arcsin(q * 0.5 / (4 * np.pi)))
    assert np.allclose(two_theta_to_q(tth, 0.5), q, rtol=1e-12)
    with pytest.raises(ValueError):
        two_theta_to_q([0.0, 10.0], 0.5)


def test_read_pattern_skips_headers_and_reads_commas(data):
    a = read_pattern(data[0])
    assert a.shape == (400, 2) and np.allclose(a[:, 0], data[1], atol=1e-6)


@pytest.mark.parametrize("comp, expected", [
    ("Cu3 Au1", (("Cu", 0.75), ("Au", 0.25))),
    ("Cu3,Au", (("Cu", 0.75), ("Au", 0.25))),
    ({"Fe": 1, "O": 1}, (("Fe", 0.5), ("O", 0.5))),
    ("Si O2", (("Si", 1 / 3), ("O", 2 / 3))),
])
def test_composition(data, comp, expected):
    c = ExpPdf(data[0], comp)
    assert [n for n, _ in c.fractions] == [n for n, _ in expected]
    assert np.allclose([f for _, f in c.fractions], [f for _, f in expected])


@pytest.mark.parametrize("comp, message", [
    ("Cu3Au", "composition token"), ("Cu:3", "composition token"), ({"Xx": 1}, "not an element"),
    ({"Cu": 0}, "amount of Cu"), ("Cu Cu", "twice"), ("", "empty"),
])
def test_bad_composition(data, comp, message):
    with pytest.raises((ValueError, TypeError), match=message):
        ExpPdf(data[0], comp)


def test_macro_is_known_to_the_catalogue_and_has_the_working_set_up(data):
    c = ExpPdf(data[0], "Cu", fourier_qmax=16.5, poly_order=5, outputs=("fq", "sq"))
    lines = c.macro_lines(data[1], with_background=False)
    text = "\n".join(lines) + "\n"
    assert all(s.kind != "unknown" for s in load_catalogue().walk(text, "suite"))
    i_read, i_pow, i_e2p = lines.index("stru e2p_comp.stru"), lines.index("powder"), lines.index("exp2pdf")
    assert i_read < i_pow < i_e2p
    assert "composition comp:current" in lines
    assert "limits inst:[0.500000,17.900000], fourier:[0.500000,16.500000]" in lines
    assert "poly order:5" in lines
    assert any(l.startswith("outfile gr:e2p_gr.dat") and "format:stan" in l and "sq:e2p_sq.dat" in l for l in lines)


def test_no_limits_line_unless_asked(data):
    lines = ExpPdf(data[0], "Cu").macro_lines(data[1], with_background=True)
    assert not any(l.startswith("limits") for l in lines)
    assert "back xy, e2p_back.xy, scale:1" in lines and "poly order:7" not in "\n".join(lines)


def test_limits_rounded_inwards_and_checked(data):
    q = data[1] + 1e-7
    lines = ExpPdf(data[0], "Cu", qmax=17.0).macro_lines(q, with_background=False)
    assert "limits inst:[0.500001,17.000000], fourier:[0.500001,17.000000]" in lines
    with pytest.raises(ValueError, match="inside the data range"):
        ExpPdf(data[0], "Cu", qmax=20.0).limits(data[1])
    with pytest.raises(ValueError, match="inside the instrument limits"):
        ExpPdf(data[0], "Cu", qmax=15.0, fourier_qmax=16.0).limits(data[1])


def test_options_are_checked(data):
    with pytest.raises(ValueError, match="2-theta data need"):
        ExpPdf(data[0], "Cu", axis="tth")
    with pytest.raises(ValueError, match="not both"):
        ExpPdf(data[0], "Cu", wavelength=0.5, energy=100)
    with pytest.raises(ValueError, match="radiation"):
        ExpPdf(data[0], "Cu", radiation="gamma")
    with pytest.raises(ValueError, match="outputs"):
        ExpPdf(data[0], "Cu", outputs=("gr",))


def test_comp_stru_layout():
    text = _comp_stru((("Cu", 0.75), ("Au", 0.25)))
    assert "ncell        1,        1,        1,       2" in text
    assert "0.750000, _," in text and text.count("\n") == 8


def test_read_exp2pdf_file(tmp_path):
    p = tmp_path / "g.dat"
    p.write_text("# CORREC : lim_four_app;2R::       1.199,     16.473\n# OUTPUT : n_points    ;1I::     2\n"
                 "  0.5  1.0\n  0.51  -2.0D+00\n")
    r, g, h = read_exp2pdf_file(p)
    assert np.allclose(r, [0.5, 0.51]) and np.allclose(g, [1, -2])
    assert h["lim_four_app"] == ["1.199", "16.473"] and h["n_points"] == ["2"]
