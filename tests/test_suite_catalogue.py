"""Tests for the bundled DISCUS_SUITE command catalogue and the macro walker."""

from pathlib import Path

import pytest

import json
import os
import subprocess
import sys

from discuspdf.suite_catalogue import Catalogue, load_catalogue, str_comp

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def cat():
    return load_catalogue("6.23.01")


def test_str_comp_matches_discus_rules():
    # lib_f90/str_comp_mod.f90: compare min(len(typed), len(name)) characters, at least `minimum`
    assert str_comp("four", "fourier", 2)
    assert str_comp("fo", "fourier", 2)
    assert not str_comp("f", "fourier", 2)
    assert str_comp("fourierx", "fourier", 2)  # longer than the name still matches
    assert not str_comp("Fourier", "fourier", 2)  # case-sensitive
    assert not str_comp("", "fourier", 1)
    # the length argument truncates the comparison: str_comp(befehl, 'seleect', 3, lbef, 6)
    assert str_comp("seleecx", "seleect", 3, 6)
    assert not str_comp("seleecx", "seleect", 3)


def test_engine_metadata(cat):
    assert cat.engine["version"] == "6.23.01"
    assert cat.engine["commit"].startswith("19a214a")
    assert cat.data["schema_version"] == 1


@pytest.mark.parametrize("menu", ["suite", "discus", "discus/fourier", "discus/powd", "discus/mmc",
                                  "discus/chem", "discus/rmc", "kuplot", "kuplot/fit", "refine",
                                  "diffev", "mixscat", "common"])
def test_expected_menus_present(cat, menu):
    assert cat.menu(menu).commands


def test_menu_entry(cat):
    assert cat.resolve("suite", "discus").enters == "discus"
    assert cat.resolve("discus", "four").enters == "discus/fourier"
    assert cat.resolve("discus", "powder").enters == "discus/powd"
    assert cat.resolve("kuplot", "fit").enters == "kuplot/fit"


def test_fallthrough_to_common(cat):
    # 'variable' is a lib_f90 command reached from the DISCUS main menu
    cmd = cat.resolve("discus", "variable")
    assert cmd is not None and cmd.menu == "common"


def test_unknown_command(cat):
    assert cat.resolve("discus/fourier", "frobnicate") is None
    assert [s.text for s in cat.check("discus\nfrobnicate 1\nexit\n")] == ["frobnicate 1"]


def test_walk_tracks_menus(cat):
    macro = "discus\nread\ncell x.cll,2,2,2\nfourier\nxray\nrun\nexit\nexit\n"
    steps = [(s.menu, s.command) for s in cat.walk(macro)]
    assert steps == [("suite", "discus"), ("discus", "read"), ("discus/read", "cell"),
                     ("discus", "fourier"), ("discus/fourier", "xray"), ("discus/fourier", "run"),
                     ("discus/fourier", "exit"), ("discus", "exit")]


def test_walk_conditional_entry(cat):
    # 'save <file>' writes at once; 'save' alone opens the save menu (discus_save.f90)
    direct = [s.menu for s in cat.walk("discus\nsave out.stru\nfourier\nexit\n")]
    assert direct == ["suite", "discus", "discus", "discus/fourier"]
    menu = [s.menu for s in cat.walk("discus\nsave\noutf out.stru\nrun\nexit\n")]
    assert menu[2:] == ["discus/save", "discus/save", "discus/save"]
    # 'insert' opens its menu only for domain/molecule/object
    assert [s.menu for s in cat.walk("discus\ninsert mol\nexit\n")][-1] == "discus/insert"
    assert [s.menu for s in cat.walk("discus\ninsert Cu,0,0,0,0.5\nfourier\n")][-1] == "discus"


def test_comments_control_and_assignments(cat):
    steps = list(cat.walk("# c\ndiscus\ndo i[1]=1,2\n r[1] = 2\nenddo\n@other.mac\n"))
    assert [s.kind for s in steps] == ["comment", "command", "control", "assignment", "control", "macro"]


def test_generated_macros_are_fully_known(cat):
    """Every macro the package's generators write uses only commands the catalogue knows."""
    sys.path.insert(0, str(ROOT / "tools" / "suite_catalogue"))
    try:
        from render_samples import sample_macros
    finally:
        sys.path.pop(0)
    samples = sample_macros()
    assert len(samples) > 40
    unknown = {label: [s.text for s in cat.check(text)] for label, text in samples}
    assert {k: v for k, v in unknown.items() if v} == {}


def test_engine_typos_and_shadowing_are_reproduced(cat):
    """Quirks of 6.23.01 found in the source (see docs/SUITE_CATALOGUE.md)."""
    assert cat.resolve("discus/pdf", "select") is None  # branch literal 'seleect', length 6
    assert cat.resolve("discus/pdf", "seleecz").name == "seleect"  # only 6 characters compared
    assert cat.resolve("kuplot", "mergeg").name == "merge"  # KUPLOT keeps 4 characters of a command
    assert cat.resolve("kuplot/fit", "con").name == "continue"  # earlier branch wins
    assert cat.resolve("kuplot/fit", "co").name == "convergence"
    assert cat.resolve("common", "fe").name == "fend"
    assert cat.resolve("discus/exp2pdf", "exit").menu == "discus/common"  # common is asked first
    assert cat.shortest_abbreviation(cat.menu("kuplot").command("mergeg")) is None
    assert cat.resolve("discus/pdf", "save").enters is None  # writes a file, never opens the save menu


def test_committed_macros_are_known(cat):
    """Every committed macro that DISCUS ran must use only commands the catalogue knows."""
    unknown = {}
    for path in sorted(ROOT.rglob("*.mac")):
        if ".git" in path.parts:
            continue
        for step in cat.check(path.read_text(errors="replace")):
            unknown.setdefault(path.relative_to(ROOT).as_posix(), []).append(f"{step.menu}: {step.text}")
    assert unknown == {}


def test_every_enters_target_exists(cat):
    for cmd in cat.commands(include_standard=True):
        if cmd.enters:
            assert cmd.enters in cat.menus, (cmd.menu, cmd.name, cmd.enters)


def test_rejects_wrong_schema():
    with pytest.raises(ValueError):
        Catalogue({"schema": "something else"})
    with pytest.raises(ValueError):
        load_catalogue("5.0")


@pytest.mark.skipif(not os.environ.get("DIFFUSECODE_SRC"),
                    reason="set DIFFUSECODE_SRC to a DiffuseCode v.6.23.01 checkout to rebuild the catalogue")
def test_bundled_catalogue_matches_a_fresh_build(cat, tmp_path):
    out = tmp_path / "cat.json"
    subprocess.run([sys.executable, str(ROOT / "tools" / "suite_catalogue" / "build_catalogue.py"),
                    os.environ["DIFFUSECODE_SRC"], str(out), "--version", cat.engine["version"],
                    "--commit", cat.engine["commit"], "--label", cat.engine["label"]], check=True)
    assert json.loads(out.read_text()) == cat.data


def _kinds(cat, macro, start="suite"):
    return [(s.menu, s.kind, s.command) for s in cat.walk(macro, start)]


def test_assignment_rules_follow_the_engine(cat):
    # trailing comment is removed first, so 'fourier # a=b' is the fourier command
    assert _kinds(cat, "fourier # a=b\nrun\n", "discus")[1] == ("discus/fourier", "command", "run")
    # '==' lines are never assignments
    assert _kinds(cat, "show a==b\n", "discus")[0][1] == "command"
    # exemptions use str_comp: 'ech a=b' is echo
    assert _kinds(cat, "ech a=b\n", "discus")[0][2] == "echo"
    # per-menu exemptions: KUPLOT titles may contain '='
    assert _kinds(cat, "tit1 a=b\n", "kuplot")[0][1] == "command"
    assert _kinds(cat, "tit1 a=b\n", "discus")[0][1] == "assignment"
    # the read menu has no do_math: 'x=1' there is just an unknown line (and read is left)
    assert _kinds(cat, "read\nx=1\nfourier\n", "discus")[1:] == [
        ("discus/read", "unknown", None), ("discus", "command", "fourier")]


def test_command_word_ends_at_blank_only(cat):
    # get_cmd splits at a blank or tab, so 'run,' is not 'run'
    assert _kinds(cat, "run,\n", "discus/fourier")[0][1] == "unknown"
    assert _kinds(cat, "run\tx\n", "discus/fourier")[0][2] == "run"


def test_blocks_only_at_program_level(cat):
    assert _kinds(cat, "do i[1]=1,3\nfourier\nrun\nexit\nenddo\n", "discus") == [
        ("discus", "control", "do"), ("discus", "command", "fourier"), ("discus/fourier", "command", "run"),
        ("discus/fourier", "command", "exit"), ("discus", "control", "enddo")]
    assert _kinds(cat, "if(a>b) then\nend if\n", "discus")[1][1] == "control"
    # inside a sub-menu a 'do' line is an assignment (do_math), not a loop
    assert _kinds(cat, "do i[1]=1,3\n", "discus/fourier")[0][1] == "assignment"


def test_continuation_lines_are_joined(cat):
    steps = list(cat.walk("discus\nfour&\nier\nrun\n"))
    assert [(s.line, s.command) for s in steps] == [(1, "discus"), (2, "fourier"), (4, "run")]


def test_read_is_left_after_any_line(cat):
    assert _kinds(cat, "read\nfoo\ninsert mol\n", "discus")[-1] == ("discus", "command", "insert")


def test_lines_after_final_exit(cat):
    assert _kinds(cat, "exit\ndiscus\n")[-1] == ("suite", "after_end", None)


def test_compare_len_and_longer_words(cat):
    assert cat.resolve("discus/fourier", "run").name == "run"
    assert cat.resolve("discus/fourier", "runx") is None  # branch 'run ' with length 4


def test_leaving_commands(cat):
    steps = [(s.menu, s.command) for s in cat.walk("kuplot\nksave 1\nrun\nload xy,a.xy\nexit\n")]
    assert steps == [("suite", "kuplot"), ("kuplot", "ksave"), ("kuplot/ksav", "run"),
                     ("kuplot", "load"), ("kuplot", "exit")]
    assert _kinds(cat, "finished\ndiscus\n")[-1][1] == "after_end"
