"""Tests for discuspdf.discus_ops: steps, StructureJob macros, and a run with a fake engine."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

from discuspdf.discus_io import read_stru
from discuspdf.discus_ops import (Purge, Raw, Remove, ReplaceAtom, ReplaceRandom, StructureJob, Switch,
                                  Thermal)
from discuspdf.engine import Engine
from discuspdf.fourier_macro import Seed, StructureSource
from discuspdf.suite_catalogue import load_catalogue
from discuspdf.suite_runner import SuiteRunner

ROOT = Path(__file__).resolve().parents[1]
CU = ROOT / "tests" / "fixtures" / "discus_ops" / "cu_b05.cll"
SEED = Seed((11, 22, 33, 44, 55, 66, 77, 88))


def lines(job):
    return job.macro().render().splitlines()


def test_job_macro_reads_seeds_applies_and_saves():
    job = StructureJob(StructureSource.cell("cu_b05.cll", (6, 6, 6)),
                       [ReplaceRandom("CU", "AU", 0.25), Thermal()], seed=SEED, title="t")
    assert lines(job) == ["# t", "discus", "read", "cell cu_b05.cll,6,6,6", "seed 11,22,33,44,55,66,77,88",
                          "replace CU,AU,all,0.25", "thermal", "save result.stru", "exit"]
    assert load_catalogue().check(job.macro().render()) == []


@pytest.mark.parametrize("step, line", [
    (Thermal(), "thermal"),
    (Thermal(two_d=True), "thermal 2d"),
    (Thermal(molecules=True), "thermal mol"),
    (Thermal(molecules=True, two_d=True), "thermal mol,2d"),
    (ReplaceAtom(3, "AU"), "replace 3,AU"),
    (ReplaceAtom(3, "AG", b=0.7), "replace 3,AG,0.7"),
    (ReplaceRandom("CU", "AU", 0.5, site=2), "replace CU,AU,2,0.5"),
    (ReplaceRandom(1, 2, 1.0), "replace 1,2,all,1"),
    (Switch(1, 2), "switch 1,2"),
    (Remove(5), "remove 5"),
    (Remove(5, 7), "remove 5,7"),
    (Remove("all"), "remove all"),
    (Purge(), "purge"),
    (Purge(types=True), "purge type:yes"),
    (Raw(("show atom,1",)), "show atom,1"),
])
def test_step_lines(step, line):
    job = StructureJob(StructureSource.stru("x.stru", (2, 2, 2)), [step])
    assert lines(job)[4] == line  # after the comment, discus, read, stru
    assert load_catalogue().check(job.macro().render()) == []


@pytest.mark.parametrize("make", [
    lambda: ReplaceRandom("CU", "AU", 0.0),  # replace refuses prob <= 0 (modify_mod.f90)
    lambda: ReplaceRandom("CU", "AU", 1.5),
    lambda: ReplaceRandom("CU", "AU", 0.5, site=0),
    lambda: ReplaceRandom("C U", "AU", 0.5),
    lambda: ReplaceAtom(0, "AU"),
    lambda: ReplaceAtom(1, "AU", b=-1),
    lambda: Switch(1, True),
    lambda: Remove("first"),
    lambda: Remove(3, 0),
])
def test_invalid_steps(make):
    with pytest.raises(ValueError):
        make()


def test_job_type_checks():
    with pytest.raises(TypeError):
        StructureJob("cu.cll")
    with pytest.raises(TypeError):
        StructureJob(StructureSource.cell("cu.cll", (2, 2, 2)), ["thermal"])
    with pytest.raises(TypeError):
        StructureJob(StructureSource.cell("cu.cll", (2, 2, 2)), seed=(1, 2, 3))


def test_raw_unknown_command_is_refused_on_render():
    job = StructureJob(StructureSource.cell("cu.cll", (2, 2, 2)), [Raw(("frobnicate 1",))])
    with pytest.raises(ValueError):
        job.macro().render()


GAAS_222 = ROOT / "tests" / "fixtures" / "discus_6_23_01" / "gaas_222.stru"


def _cu_structure(n=2):
    """A configuration: the committed DISCUS-saved GaAs 2x2x2 crystal (64 atoms)."""
    assert n == 2
    return read_stru(GAAS_222, supercell=(2, 2, 2)).structure


def test_structure_source_is_written_as_input_stru(tmp_path):
    s = _cu_structure()
    job = StructureJob(s, [Thermal()], seed=SEED)
    assert lines(job)[3] == "structure input.stru"


FAKE = """
import sys
from pathlib import Path
src = Path("input.stru") if Path("input.stru").exists() else None
text = src.read_text() if src else Path("cu_b05.cll").read_text()
# pretend DISCUS: write back the input crystal unchanged, without an ncell line
lines = [l for l in text.splitlines() if not l.lower().startswith("ncell")]
Path("result.stru").write_text("\\n".join(lines) + "\\n")
print("done")
"""


def test_run_reads_result_with_the_source_crystal_size(tmp_path):
    script = tmp_path / "fake.py"
    script.write_text(FAKE)
    runner = SuiteRunner(tmp_path / "runs", Engine((sys.executable, str(script))))
    s = _cu_structure(2)
    result = StructureJob(s, [Thermal()], seed=SEED).run(runner)
    assert result.run.status()["status"] == "complete"
    assert result.structure.n_atoms == s.n_atoms == 64
    assert tuple(result.structure.supercell) == (2, 2, 2)
    assert np.allclose(result.structure.frac, s.frac)


def test_run_needs_the_source_file(tmp_path):
    runner = SuiteRunner(tmp_path / "runs", Engine((sys.executable, "-c", "pass")))
    with pytest.raises(ValueError, match="inputs must contain"):
        StructureJob(StructureSource.cell("cu_b05.cll", (2, 2, 2))).run(runner, {})


def test_result_that_is_not_a_whole_crystal_raises_with_the_run(tmp_path):
    from discuspdf.discus_ops import Purge, ResultNotAStructure
    script = tmp_path / "fake.py"
    script.write_text(FAKE.replace('lines = [l for l in text.splitlines() if not l.lower().startswith("ncell")]',
                                   'lines = [l for l in text.splitlines() if not l.lower().startswith("ncell")][:-3]'))
    runner = SuiteRunner(tmp_path / "runs", Engine((sys.executable, str(script))))
    with pytest.raises(ResultNotAStructure) as exc:
        StructureJob(_cu_structure(2), [Purge()]).run(runner)
    assert exc.value.run.output("result.stru").is_file()
