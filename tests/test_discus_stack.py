"""Tests for discuspdf.discus_stack: the stacking-fault builder (DISCUS ``stack`` menu)."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

from discuspdf.discus_io import DiscusFormatError
from discuspdf.discus_ops import ResultNotAStructure, StructureJob
from discuspdf.discus_stack import Stacking
from discuspdf.engine import Engine
from discuspdf.fourier_macro import Seed
from discuspdf.suite_catalogue import load_catalogue
from discuspdf.suite_runner import SuiteRunner

FIX = Path(__file__).resolve().parent / "fixtures" / "discus_stack"
SEED = Seed((11, 22, 33, 44, 55, 66, 77, 88))

THIRD = 0.3333333333333333
FULL_NAME = {"layer": "layer", "number": "number", "crow": "crow", "first": "first",
             "trans": "translation", "sigma": "sigma", "modulus": "modulus", "aver": "average",
             "set": "set", "random": "random", "create": "create", "run": "run", "exit": "exit"}


def abab(**kw):
    """Fault-free two-layer ABAB stacking: type 1 -> 2 -> 1, both steps 0.5 cell along c."""
    half = (0, 0, 0.5)
    args = dict(layers=("cu_layer.stru", "ag_layer.stru"), correlation=((0, 1), (1, 0)),
                translations=((half, half), (half, half)), modulus=((1, 0, 0), (0, 1, 0)),
                number=8, first=1)
    args.update(kw)
    return Stacking(**args)


def fcc(**kw):
    """Close packing in hexagonal axes (a = 2.864 A, c = 7.015 A for Al; A, B, C in plane at (0,0),
    (1/3,2/3), (2/3,1/3)) as a first-order Markov sequence of positions: after A comes B with
    probability 1 - alpha and C with alpha (and cyclically). alpha = 0 is fcc ABC; alpha > 0 mixes
    in reversed steps. Translation i->j is the in-plane offset of the next position plus 1/3 along
    c; the average stacking vector is the straight (0, 0, 1/3) with 'set trans:fixed' (a sheared
    average would make the origins drift, see the Stacking docstring)."""
    pos = {0: (0.0, 0.0), 1: (1 / 3, 2 / 3), 2: (2 / 3, 1 / 3)}
    trans = tuple(tuple(((pos[j][0] - pos[i][0]) % 1.0, (pos[j][1] - pos[i][1]) % 1.0, 1 / 3)
                        for j in range(3)) for i in range(3))
    alpha = kw.pop("alpha", 0.1)
    args = dict(layers=("fcc_a.stru", "fcc_b.stru", "fcc_c.stru"),
                correlation=((0, 1 - alpha, alpha), (alpha, 0, 1 - alpha), (1 - alpha, alpha, 0)),
                translations=trans, modulus=((1, 0, 0), (0, 1, 0)), average=(0, 0, 1 / 3),
                number=30, first="random", modulus_on=True)
    args.update(kw)
    return Stacking(**args)


def stack_block(text):
    lines = text.splitlines()
    start = lines.index("stack")
    end = lines.index("exit", start)
    return lines[start + 1:end]


def test_abab_macro_is_exact():
    job = StructureJob(abab(), title="ABAB")
    assert job.macro().render().splitlines() == [
        "# ABAB", "discus", "stack",
        "layer cu_layer.stru", "layer ag_layer.stru", "number 8",
        "crow 1,0,1", "crow 2,1,0", "first 1",
        "trans 1,1,0,0,0.5", "trans 1,2,0,0,0.5", "trans 2,1,0,0,0.5", "trans 2,2,0,0,0.5",
        "modulus 1,0,0,0,1,0", "set trans:average", "set modulus:off", "create", "run", "exit",
        "save result.stru", "exit"]


def test_fcc_markov_macro_is_exact():
    job = StructureJob(fcc(), seed=SEED, title="fcc")
    assert job.macro().render().splitlines() == [
        "# fcc", "discus", "seed 11,22,33,44,55,66,77,88", "stack",
        "layer fcc_a.stru", "layer fcc_b.stru", "layer fcc_c.stru", "number 30",
        "crow 1,0,0.9,0.1", "crow 2,0.1,0,0.9", "crow 3,0.9,0.1,0", "first random",
        "trans 1,1,0,0,0.333333333333", "trans 1,2,0.333333333333,0.666666666667,0.333333333333",
        "trans 1,3,0.666666666667,0.333333333333,0.333333333333",
        "trans 2,1,0.666666666667,0.333333333333,0.333333333333", "trans 2,2,0,0,0.333333333333",
        "trans 2,3,0.333333333333,0.666666666667,0.333333333333",
        "trans 3,1,0.333333333333,0.666666666667,0.333333333333",
        "trans 3,2,0.666666666667,0.333333333333,0.333333333333", "trans 3,3,0,0,0.333333333333",
        "modulus 1,0,0,0,1,0", "aver 0,0,0.333333333333", "set trans:fixed", "set modulus:on",
        "create", "run", "exit", "save result.stru", "exit"]


@pytest.mark.parametrize("make, seeded", [(abab, False), (fcc, True)])
def test_macros_pass_the_catalogue_and_every_stack_word_resolves(make, seeded):
    job = StructureJob(make(), seed=SEED if seeded else None)
    text = job.macro().render()
    cat = load_catalogue()
    assert cat.check(text) == []
    block = stack_block(text)
    assert block, "no stack commands"
    for line in block:
        word = line.split()[0]
        cmd = cat.resolve("discus/stack", word)
        assert cmd is not None, word
        assert cmd.name == FULL_NAME[word], (line, cmd.name)


def test_random_numbers_need_a_seed_and_the_flag_follows_the_chain():
    assert abab().uses_random_numbers is False  # deterministic rows, fixed first layer
    assert abab(first="random").uses_random_numbers is True
    assert fcc().uses_random_numbers is True
    with pytest.raises(ValueError, match="draws random numbers"):
        StructureJob(fcc())
    StructureJob(abab())  # no seed needed
    assert abab(correlation=((0.5, 0.5), (1, 0)), translations=((
        (0, 0, 0.5), (0, 0, 0.5)), ((0, 0, 0.5), (0, 0, 0.5)))).uses_random_numbers is True


def test_input_names_are_the_layer_files_once_each_in_order():
    assert abab().input_names() == ("cu_layer.stru", "ag_layer.stru")
    same = abab(layers=("cu_layer.stru", "cu_layer.stru"))
    assert same.input_names() == ("cu_layer.stru",)


def test_fixed_average_writes_aver_and_set_fixed():
    text = StructureJob(abab(average=(0, 0, 0.5))).macro().render()
    block = stack_block(text)
    assert "aver 0,0,0.5" in block
    assert "set trans:fixed" in block
    assert "set trans:average" not in block
    assert load_catalogue().check(text) == []


def test_sigmas_random_offset_sigma_and_atoms_are_written_only_when_set():
    job = StructureJob(abab(sigmas=(((0, 0, 0), (0.01, 0, 0)), ((0, 0, 0), (0, 0, 0))),
                            random_probability=0.2, random_offset=(0.1, 0, 0), random_sigma=(0, 0.02, 0),
                            atoms="strict"),
                       seed=SEED)
    block = stack_block(job.macro().render())
    assert "sigma 1,2,0.01,0,0" in block
    assert "random prob,0.2" in block and "random offset,0.1,0,0" in block
    assert "random sigma,0,0.02,0" in block and "set atoms:strict" in block
    plain = stack_block(StructureJob(abab()).macro().render())
    assert not any(x.startswith(("sigma", "random", "set atoms")) for x in plain)


def test_fixture_layers_are_discus_structure_files():
    from discuspdf.discus_io import read_stru
    for name in ("cu_layer", "ag_layer", "fcc_a", "fcc_b", "fcc_c"):
        assert read_stru(FIX / f"{name}.stru", supercell=(1, 1, 1)).structure.n_atoms == 1


VALID = dict(layers=("cu_layer.stru", "ag_layer.stru"), correlation=((0, 1), (1, 0)),
             translations=(((0, 0, 0.5), (0, 0, 0.5)), ((0, 0, 0.5), (0, 0, 0.5))),
             modulus=((1, 0, 0), (0, 1, 0)), number=8, first=1)


@pytest.mark.parametrize("change, match", [
    ({"layers": ()}, "at least one"),
    ({"layers": ("a b.stru", "ag_layer.stru")}, "plain DISCUS"),
    ({"layers": ("../cu.stru", "ag_layer.stru")}, "plain DISCUS"),
    ({"layers": ("cu,x.stru", "ag_layer.stru")}, "plain DISCUS"),
    ({"layers": ("cu_layer.cll", "ag_layer.stru")}, "plain DISCUS"),
    ({"layers": ("internal.x.stru", "ag_layer.stru")}, "internal"),
    ({"correlation": ((0, 1), (1,))}, "entries"),
    ({"correlation": ((0, 1), (1, 0), (0, 0))}, "rows"),
    ({"correlation": ((0.5, 0.4), (1, 0))}, "sums to"),
    ({"correlation": ((0, 1.5), (1, -0.5))}, "probability"),
    ({"correlation": ((0, float("nan")), (1, 0))}, "probability"),
    ({"translations": (((0, 0, 0.5),),)}, "rows"),
    ({"translations": (((0, 0, float("inf")), (0, 0, 0.5)), ((0, 0, 0.5), (0, 0, 0.5)))}, "finite"),
    ({"translations": (((0, 0, 0), (0, 0, 0)), ((0, 0, 0), (0, 0, 0)))}, "is zero"),
    ({"modulus": ((1, 0, 0),)}, "two vectors"),
    ({"modulus": ((0, 0, 0), (0, 1, 0))}, "zero vector"),
    ({"modulus": ((1, 0, 0), (2, 0, 0))}, "coplanar"),
    ({"average": (2, 0, 0)}, "coplanar"),
    ({"average": (0, 0, 0)}, "zero vector"),
    ({"number": 0}, "integer"),
    ({"number": 2.5}, "integer"),
    ({"number": True}, "integer"),
    ({"first": 0}, ">= 1"),
    ({"first": 3}, "larger than"),
    ({"first": "rand"}, "random"),
    ({"modulus_on": "yes"}, "True or False"),
    ({"random_probability": 1.5}, "probability"),
    ({"random_offset": (0, 0)}, "three components"),
    ({"random_sigma": (0, -0.1, 0)}, "negative"),
    ({"sigmas": (((0, 0, -1), (0, 0, 0)), ((0, 0, 0), (0, 0, 0)))}, "negative"),
    ({"atoms": "fast"}, "atoms must be"),
])
def test_invalid_stacks_are_refused(change, match):
    args = dict(VALID)
    args.update(change)
    with pytest.raises(ValueError, match=match):
        Stacking(**args)


def test_read_result_is_one_box_of_discus_estimated_cells(tmp_path):
    # 8 layers 0.5 cells apart span z 0..3.5 -> DISCUS's estimate floor(3.5)+1 = 4 cells along c
    f = tmp_path / "result.stru"
    f.write_text(saved_stack(8))
    structure, stru = Stacking(**VALID).read_result(f, True)
    assert stru is None and tuple(structure.supercell) == (1, 1, 1)
    assert structure.n_atoms == 8
    assert np.allclose(structure.lattice.parameters[:3], (4.0, 4.0, 16.0))
    assert np.allclose(structure.frac[:, 2], np.arange(8) * 0.5 / 4)


def test_read_result_keeps_any_layer_count(tmp_path):
    # 9 layers reach z = 4.0: 5 cells, 9 atoms (not a whole number per cell) - still one box
    f = tmp_path / "result.stru"
    f.write_text(saved_stack(9))
    structure, _ = Stacking(**VALID).read_result(f, True)
    assert structure.n_atoms == 9 and np.isclose(structure.lattice.parameters[2], 20.0)


def saved_stack(n_layers, dz=0.5):
    """A DISCUS-style saved crystal (no ncell line) of alternating CU/AG layers, one atom each."""
    lines = ["title stack", "spcgr P 1", "cell 4.0, 4.0, 4.0, 90.0, 90.0, 90.0",
             "format numbers,XYZBPMMOS",
             "atoms      x,              y,              z,             Biso,    Property,  MoleNo,  "
             "MoleAt,   Occ,     St,  Sh,  Sk,  Sl"]
    for k in range(n_layers):
        name = "CU" if k % 2 == 0 else "AG"
        lines.append(f"{name} 0.000000, 0.000000, {k * dz:.6f}, 0.500000, 1, 0, 0, 1.000000, _, 0, 0, 0")
    return "\n".join(lines) + "\n"


def test_job_run_reads_the_stack_back_through_a_fake_engine(tmp_path):
    fake = tmp_path / "fake.py"
    fake.write_text("from pathlib import Path\nPath('result.stru').write_text(%r)\n" % saved_stack(8))
    runner = SuiteRunner(tmp_path / "runs", Engine((sys.executable, str(fake))))
    inputs = {"cu_layer.stru": FIX / "cu_layer.stru", "ag_layer.stru": FIX / "ag_layer.stru"}
    result = StructureJob(Stacking(**VALID), seed=SEED).run(runner, inputs)
    assert result.run.status()["status"] == "complete"
    assert tuple(result.structure.supercell) == (1, 1, 1)
    assert result.structure.n_atoms == 8


def test_job_run_missing_layer_file_is_refused(tmp_path):
    runner = SuiteRunner(tmp_path / "runs", Engine((sys.executable, "-c", "pass")))
    with pytest.raises(ValueError, match="inputs must contain"):
        StructureJob(Stacking(**VALID), seed=SEED).run(runner, {"cu_layer.stru": FIX / "cu_layer.stru"})


def test_job_run_unreadable_result_raises_with_the_run(tmp_path):
    fake = tmp_path / "fake.py"
    fake.write_text("from pathlib import Path\nPath('result.stru').write_text('title x\\natoms\\n')\n")
    runner = SuiteRunner(tmp_path / "runs", Engine((sys.executable, str(fake))))
    inputs = {"cu_layer.stru": FIX / "cu_layer.stru", "ag_layer.stru": FIX / "ag_layer.stru"}
    with pytest.raises(ResultNotAStructure) as exc:
        StructureJob(Stacking(**VALID), seed=SEED).run(runner, inputs)
    assert exc.value.run.output("result.stru").is_file()

