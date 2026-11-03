"""Tests for discuspdf.discus_domain: the domain menu step (macro lines, validation, catalogue)."""
from __future__ import annotations

import math

import pytest

from discuspdf.discus_domain import IDENTITY_ORIENT, Domains
from discuspdf.discus_ops import StructureJob
from discuspdf.fourier_macro import Seed, StructureSource
from discuspdf.suite_catalogue import load_catalogue

SEED = Seed((11, 22, 33, 44, 55, 66, 77, 88))
HOST = StructureSource.cell("cu_b05.cll", (6, 6, 6))
DIAG = ((1.5, 0.0, 0.0, 0.0), (0.0, 1.5, 0.0, 0.0), (0.0, 0.0, 1.5, 0.0))
ROT_Z = ((0.0, -1.0, 0.0, 0.0), (1.0, 0.0, 0.0, 0.0), (0.0, 0.0, 1.0, 0.0))

SPHERE_STEP = dict(content="dom_cu.stru", inputfile="pseudo.stru", name="dom1", character="sphere",
                   shape=DIAG, fuzzy=2.0)
CUBE_STEP = dict(content="dom_fe.stru", inputfile="pseudo.stru", name="CUB1", character="cube",
                 shape=((2.0, 0.0, 0.0, 0.5), (0.0, 2.0, 0.0, 0.0), (0.0, 0.0, 2.0, -0.25)),
                 orient=ROT_Z, fuzzy=1.5)

# the lines the domain step writes inside the DISCUS main menu (between 'domain' and 'exit')
SPHERE_LINES = [
    "domain",
    "reset",
    "mode pseudoatoms",
    "inputfile pseudo.stru",
    "assign character,DOM1,sphere",
    "assign content,DOM1,dom_cu.stru",
    "assign fuzzy,DOM1,2",
    "assign shape,DOM1,1,1.5,0,0,0",
    "assign shape,DOM1,2,0,1.5,0,0",
    "assign shape,DOM1,3,0,0,1.5,0",
    "assign orient,DOM1,1,1,0,0,0",
    "assign orient,DOM1,2,0,1,0,0",
    "assign orient,DOM1,3,0,0,1,0",
    "run",
    "exit",
]

CUBE_LINES = [
    "domain",
    "reset",
    "mode pseudoatoms",
    "inputfile pseudo.stru",
    "assign character,CUB1,cube",
    "assign content,CUB1,dom_fe.stru",
    "assign fuzzy,CUB1,1.5",
    "assign shape,CUB1,1,2,0,0,0.5",
    "assign shape,CUB1,2,0,2,0,0",
    "assign shape,CUB1,3,0,0,2,-0.25",
    "assign orient,CUB1,1,0,-1,0,0",
    "assign orient,CUB1,2,1,0,0,0",
    "assign orient,CUB1,3,0,0,1,0",
    "run",
    "exit",
]


def macro_lines(step) -> list[str]:
    job = StructureJob(HOST, [step], seed=SEED, title="t")
    return job.macro().render().splitlines()


def domain_block(step) -> list[str]:
    lines = macro_lines(step)
    start = lines.index("domain")
    return lines[start:lines.index("exit", start) + 1]


def test_sphere_domain_lines():
    assert domain_block(Domains(**SPHERE_STEP)) == SPHERE_LINES


def test_cube_domain_with_orientation_lines():
    assert domain_block(Domains(**CUBE_STEP)) == CUBE_LINES


def test_whole_job_lines_and_catalogue():
    job = StructureJob(HOST, [Domains(**SPHERE_STEP)], seed=SEED, title="t")
    assert macro_lines(Domains(**SPHERE_STEP))[:5] == [
        "# t", "discus", "read", "cell cu_b05.cll,6,6,6", "seed 11,22,33,44,55,66,77,88"]
    assert macro_lines(Domains(**SPHERE_STEP))[-2:] == ["save result.stru", "exit"]
    assert load_catalogue().check(job.macro().render()) == []


@pytest.mark.parametrize("kwargs", [SPHERE_STEP, CUBE_STEP])
def test_whole_job_passes_catalogue_check(kwargs):
    job = StructureJob(HOST, [Domains(**kwargs)], seed=SEED)
    assert load_catalogue().check(job.macro().render()) == []


def test_identity_orientation_is_written_when_none_given():
    lines = domain_block(Domains(**{**SPHERE_STEP, "orient": None}))
    orient = [line for line in lines if line.startswith("assign orient")]
    assert orient == ["assign orient,DOM1,1,1,0,0,0", "assign orient,DOM1,2,0,1,0,0",
                      "assign orient,DOM1,3,0,0,1,0"]
    assert Domains(**SPHERE_STEP).orient == IDENTITY_ORIENT


def test_uses_random_numbers_and_needs_a_seed():
    step = Domains(**SPHERE_STEP)
    assert step.uses_random_numbers is True
    with pytest.raises(ValueError, match="seed"):
        StructureJob(HOST, [step])


def test_input_names_are_content_and_pseudo_file():
    assert Domains(**CUBE_STEP).input_names() == ("dom_fe.stru", "pseudo.stru")


def test_name_and_character_are_normalised():
    step = Domains(**{**SPHERE_STEP, "name": "dom1", "character": "Sphere"})
    assert (step.name, step.character) == ("DOM1", "sphere")


@pytest.mark.parametrize("word, full", [
    ("reset", "reset"),
    ("mode", "mode"),
    ("inputfile", "inputfile"),
    ("assign", "assign"),
    ("run", "run"),
])
def test_every_emitted_word_resolves_to_its_command(word, full):
    cmd = load_catalogue().resolve("discus/domain", word)
    assert cmd is not None and cmd.name == full
    emitted = {line.split()[0].split(",")[0] for line in domain_block(Domains(**CUBE_STEP))[1:-1]}
    assert word in emitted


@pytest.mark.parametrize("override, message", [
    ({"name": "TOOLONG"}, "pseudo-atom name"),
    ({"name": "A,B"}, "pseudo-atom name"),
    ({"name": ""}, "pseudo-atom name"),
    ({"character": "irregular"}, "character must be"),
    ({"character": "box"}, "character must be"),
    ({"content": "dir/dom.stru"}, "plain file name"),
    ({"content": "dom cu.stru"}, "plain file name"),
    ({"content": "dom,cu.stru"}, "plain file name"),
    ({"content": ""}, "content must be a file name"),
    ({"content": "internal_dom.stru"}, "internal"),
    ({"inputfile": "internal_domain_list.stru"}, "internal"),
    ({"inputfile": "pseudo.stru # comment"}, "plain file name"),
    ({"content": "x" * 201}, "longer than"),
    ({"shape": DIAG[:2]}, "3 rows of 4"),
    ({"shape": ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))}, "3 rows of 4"),
    ({"shape": ((1.0, 0.0, 0.0, 0.0), (0.0, 1.0, 0.0, 0.0), (0.0, 0.0, math.nan, 0.0))}, "finite"),
    ({"shape": ((1.0, 0.0, 0.0, 0.0), (0.0, 1.0, 0.0, 0.0), (0.0, 0.0, math.inf, 0.0))}, "finite"),
    ({"shape": ((1.0, 0.0, 0.0, 0.0), (0.0, -1.0, 0.0, 0.0), (0.0, 0.0, 1.0, 0.0))}, "negative"),
    ({"shape": ((1.0, True, 0.0, 0.0), (0.0, 1.0, 0.0, 0.0), (0.0, 0.0, 1.0, 0.0))}, "finite"),
    ({"orient": ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))}, "3 rows of 4"),
    ({"orient": ((1.0, 0.0, 0.0, 0.0), (0.0, 1.0, 0.0, 0.0), (0.0, 0.0, math.nan, 0.0))}, "finite"),
    ({"fuzzy": -0.5}, "fuzzy must be >= 0"),
    ({"fuzzy": math.inf}, "finite"),
    ({"mode": "domain"}, "mode must be"),
    ({"mode": "irregular"}, "mode must be"),
])
def test_validation_errors(override, message):
    with pytest.raises(ValueError, match=message):
        Domains(**{**SPHERE_STEP, **override})


def test_missing_fields_are_rejected():
    with pytest.raises(TypeError):
        Domains(content="dom_cu.stru", name="DOM1", character="sphere", shape=DIAG)  # no inputfile
