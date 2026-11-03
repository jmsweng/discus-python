"""Tests for the ``Wave`` step (DISCUS ``waves`` menu, ``discus/wave``).

Source references are to ``discus/prog/waves.f90`` in the DISCUS 6.23.01 tree; the test
expectations follow the command order that the source reads (reset, settings, run).
"""
from __future__ import annotations

import math

import pytest

from discuspdf.discus_ops import Wave, StructureJob
from discuspdf.fourier_macro import StructureSource
from discuspdf.suite_catalogue import Menu, load_catalogue

WAVE_MENU = "discus/wave"

# Full command name the source dispatches to, per word the step writes (waves.f90 dispatch).
INTENDED = {
    "reset": "reset",            # 163
    "long": "long",              # 289
    "transverse": "transverse",  # 497
    "density": "density",        # 204
    "accoustic": "accoustic",    # 182
    "optical": "optical",        # 294
    "function": "function",      # 229
    "length": "length",          # 272
    "amplitude": "amplitude",    # 187
    "vector": "vector",          # 502
    "phase": "phase",            # 318
    "shift": "shift",            # 475
    "osci": "osci",              # 299
    "select": "select",          # 218
    "plow": "plow",              # 357
    "phigh": "phigh",            # 340
    "replace": "replace",        # 408
    "run": "run",                # 470
}


@pytest.fixture(scope="module")
def cat():
    return load_catalogue("6.23.01")


def _job(*steps):
    return StructureJob(StructureSource.cell("cu.cll", (2, 2, 2)), list(steps))


def _block(text: str) -> list[str]:
    """The lines of the waves menu, between its opening line and its exit."""
    lines = text.splitlines()
    start = lines.index("waves")
    end = lines.index("exit", start)
    return lines[start + 1:end]


# -- rendered macro text ---------------------------------------------------------------------

def test_transverse_sinus_exact_lines():
    step = Wave("transverse", (1, 0, 0), 20, function="sinus", amplitude=0.05,
                oscillation=(0, 1, 0), select=("Cu",), phase=0)
    assert _block(_job(step).macro().render()) == [
        "reset",
        "transverse",
        "accoustic",
        "function sinus",
        "length 20",
        "amplitude 0.05",
        "vector 1,0,0",
        "phase 0",
        "osci 0,1,0",
        "select Cu",
        "run",
    ]


def test_longitudinal_optical_box_with_shift_exact_lines():
    step = Wave("longitudinal", (1, 1, 0), 10, function="box", box_width=0.3, amplitude=0.1,
                shift=0.02, branch="optical", select=["Cu", 3], phase=90)
    assert _block(_job(step).macro().render()) == [
        "reset",
        "long",
        "optical",
        "function box,0.3",
        "length 10",
        "amplitude 0.1",
        "vector 1,1,0",
        "phase 90",
        "shift 0.02",
        "select Cu,3",
        "run",
    ]


def test_longitudinal_defaults_select_all_acoustic_triangle_random_phase():
    step = Wave("longitudinal", [0, 0, 1], 5, function="triangle", amplitude=0.2, phase="random")
    assert _block(_job(step).macro().render()) == [
        "reset",
        "long",
        "accoustic",
        "function triangle",
        "length 5",
        "amplitude 0.2",
        "vector 0,0,1",
        "phase random",
        "select all",
        "run",
    ]


def test_density_exact_lines_with_viceversa():
    step = Wave("density", (0, 0, 1), 30, function="box", box_width=0.4, replace=("Cu", "Au"),
                plow=0.2, phigh=0.8, viceversa=True, phase="random")
    assert _block(_job(step).macro().render()) == [
        "reset",
        "density",
        "function box,0.4",
        "length 30",
        "vector 0,0,1",
        "phase random",
        "plow 0.2",
        "phigh 0.8",
        "replace Cu,Au,viceversa",
        "run",
    ]


def test_density_without_viceversa_has_two_replace_fields():
    step = Wave("density", (1, 0, 0), 8, replace=(1, "Au"), plow=0, phigh=1)
    lines = _block(_job(step).macro().render())
    assert lines[-2] == "replace 1,Au"
    assert "select all" not in lines and "accoustic" not in lines


def test_the_whole_macro_passes_the_catalogue_check(cat):
    steps = [
        Wave("transverse", (1, 0, 0), 20, amplitude=0.05, oscillation=(0, 1, 0), select="all"),
        Wave("density", (0, 0, 1), 30, replace=("Cu", "Au"), plow=0.2, phigh=0.8, viceversa=True),
    ]
    for step in steps:
        text = _job(step).macro().render()
        assert cat.check(text) == []


def test_reset_is_written_only_when_the_menu_has_it(monkeypatch):
    original = Menu.command

    def without_reset(self, name):
        if self.path == WAVE_MENU and name == "reset":
            return None
        return original(self, name)

    monkeypatch.setattr(Menu, "command", without_reset)
    step = Wave("longitudinal", (1, 0, 0), 10, amplitude=0.1)
    lines = _block(_job(step).macro().render())
    assert "reset" not in lines
    assert lines[0] == "long"


def test_wave_is_exported():
    import discuspdf.discus_ops as ops

    assert "Wave" in ops.__all__


# -- every emitted command resolves to the intended one ----------------------------------------

@pytest.mark.parametrize("step", [
    Wave("transverse", (1, 0, 0), 20, amplitude=0.05, oscillation=(0, 1, 0), select="all",
         branch="optical", shift=0.01, function="box", box_width=0.5, phase="random"),
    Wave("longitudinal", (1, 0, 0), 20, amplitude=0.05, select=("Cu", 2)),
    Wave("density", (0, 0, 1), 30, replace=("Cu", "Au"), plow=0.1, phigh=0.9, viceversa=True),
])
def test_every_emitted_command_resolves_to_the_intended_command(cat, step):
    lines = _block(_job(step).macro().render())
    for line in lines:
        word = line.split()[0]
        assert word in INTENDED, f"unexpected word {word!r}"
        cmd = cat.resolve(WAVE_MENU, word)
        assert cmd is not None, word
        assert cmd.name == INTENDED[word], (line, cmd.name)


def test_short_forms_resolve_to_the_source_commands(cat):
    expected = {"acco": "accoustic", "opti": "optical", "func": "function", "len": "length",
                "amp": "amplitude", "vect": "vector", "osci": "osci", "tran": "transverse",
                "dens": "density", "repl": "replace", "sel": "select"}
    for word, name in expected.items():
        assert cat.resolve(WAVE_MENU, word).name == name, word


# -- invalid inputs --------------------------------------------------------------------------

def _longitudinal(**over):
    kwargs = dict(kind="longitudinal", vector=(1, 0, 0), length=10.0, amplitude=0.1)
    kwargs.update(over)
    return kwargs


def _density(**over):
    kwargs = dict(kind="density", vector=(1, 0, 0), length=10.0, replace=("Cu", "Au"),
                  plow=0.0, phigh=1.0)
    kwargs.update(over)
    return kwargs


INVALID = [
    ("unknown kind (molecules are not supported)", dict(_longitudinal(kind="rotation"))),
    ("zero vector", dict(_longitudinal(vector=(0, 0, 0)))),
    ("two-component vector", dict(_longitudinal(vector=(1, 0)))),
    ("nan in vector", dict(_longitudinal(vector=(1, math.nan, 0)))),
    ("boolean vector entry", dict(_longitudinal(vector=(True, 0, 0)))),
    ("zero length", dict(_longitudinal(length=0.0))),
    ("negative length", dict(_longitudinal(length=-1.0))),
    ("infinite length", dict(_longitudinal(length=math.inf))),
    ("unknown function", dict(_longitudinal(function="gauss"))),
    ("box width at 1", dict(_longitudinal(function="box", box_width=1.0))),
    ("box width at 0", dict(_longitudinal(function="box", box_width=0.0))),
    ("box width for sinus", dict(_longitudinal(function="sinus", box_width=0.3))),
    ("phase as text", dict(_longitudinal(phase="Random"))),
    ("non-finite phase", dict(_longitudinal(phase=math.nan))),
    ("unknown branch", dict(_longitudinal(branch="chiral"))),
    ("missing amplitude", dict(_longitudinal(amplitude=None))),
    ("non-finite shift", dict(_longitudinal(shift=math.inf))),
    ("transverse without oscillation", dict(_longitudinal(kind="transverse"))),
    ("transverse with zero oscillation", dict(_longitudinal(kind="transverse", oscillation=(0, 0, 0)))),
    ("oscillation on a longitudinal wave", dict(_longitudinal(oscillation=(0, 1, 0)))),
    ("empty select", dict(_longitudinal(select=()))),
    ("all with other atoms", dict(_longitudinal(select=("all", "Cu")))),
    ("atom name with a space", dict(_longitudinal(select=("Cu O",)))),
    ("atom name with a comma", dict(_longitudinal(select=("Cu,O",)))),
    ("viceversa on a displacement wave", dict(_longitudinal(viceversa=True))),
    ("replace on a displacement wave", dict(_longitudinal(replace=("Cu", "Au")))),
    ("plow on a displacement wave", dict(_longitudinal(plow=0.1))),
    ("density without replace", dict(_density(replace=None))),
    ("density replace with three entries", dict(_density(replace=("Cu", "Au", "Ag")))),
    ("density without plow", dict(_density(plow=None))),
    ("density with phigh below plow", dict(_density(plow=0.6, phigh=0.4))),
    ("density with phigh above 1", dict(_density(phigh=1.5))),
    ("density with negative plow", dict(_density(plow=-0.1))),
    ("density with amplitude (overwritten by the source)", dict(_density(amplitude=0.1))),
    ("density with shift (overwritten by the source)", dict(_density(shift=0.1))),
    ("density with branch", dict(_density(branch="optical"))),
    ("density with select", dict(_density(select=("Cu",)))),
    ("density with oscillation", dict(_density(oscillation=(0, 1, 0)))),
    ("non-boolean viceversa", dict(_density(viceversa="yes"))),
]


@pytest.mark.parametrize("kwargs", [k for _, k in INVALID], ids=[n for n, _ in INVALID])
def test_invalid_inputs_raise_value_error(kwargs):
    with pytest.raises(ValueError):
        Wave(**kwargs)


def test_valid_defaults_are_filled_in():
    step = Wave("transverse", [1, 0, 0], 20, amplitude=0.05, oscillation=[0, 1, 0])
    assert step.branch == "acoustic"
    assert step.select == ("all",)
    assert step.vector == (1.0, 0.0, 0.0)
    assert step.phase == 0.0
    density = Wave(**_density())
    assert density.plow == 0.0 and density.phigh == 1.0 and density.viceversa is False


def test_macro_writes_no_site_or_molecule_commands():
    text = _job(Wave("transverse", (1, 0, 0), 20, amplitude=0.05, oscillation=(0, 1, 0))).macro().render()
    assert "site" not in text
    assert "msel" not in text and "mrep" not in text and "rot" not in text.split()
