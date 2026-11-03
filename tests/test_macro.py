"""Tests for the checked macro builder (discuspdf.macro); no DISCUS engine is run."""
import doctest
import math
from pathlib import Path

import pytest

import discuspdf.macro as macro_module
from discuspdf.macro import Macro, MacroError, format_value
from discuspdf.suite_catalogue import load_catalogue


# --- format_value -----------------------------------------------------------

@pytest.mark.parametrize("value, text", [
    (0, "0"),
    (-4, "-4"),
    (10**9, "1000000000"),
    (0.1 + 0.2, "0.3"),
    (1.5406, "1.5406"),
    (1e-12, "1e-12"),
    (2.0, "2"),
    ("n[1]", "n[1]"),
    ("  gaas.cll ", "gaas.cll"),
    ("'a#b'", "'a#b'"),
    ('"x!y"', '"x!y"'),
])
def test_format_value_accepts(value, text):
    assert format_value(value) == text


def test_format_value_float_uses_twelve_significant_digits():
    assert format_value(1 / 3) == "0.333333333333"


@pytest.mark.parametrize("value", [True, False, math.nan, math.inf, -math.inf, "", "   ",
                                   "a\nb", "a\rb", "a&", "x#y", "x!y", None, [1, 2]])
def test_format_value_rejects(value):
    with pytest.raises(MacroError):
        format_value(value)


# --- module docstring example -----------------------------------------------

def test_module_docstring_example_renders_as_shown():
    finder = doctest.DocTestFinder()
    runner = doctest.DocTestRunner(optionflags=doctest.NORMALIZE_WHITESPACE)
    tests = [t for t in finder.find(macro_module, "discuspdf.macro") if t.examples]
    assert tests, "docstring example not found"
    for test in tests:
        runner.run(test)
    assert runner.failures == 0
    assert runner.tries >= 2


# --- menus ------------------------------------------------------------------

def test_menu_blocks_write_exit_at_end_of_each_block():
    m = Macro("nested")
    with m.program("discus"):
        with m.menu("fourier"):
            m.command("run")
    text = m.render()
    assert text == "# nested\ndiscus\nfourier\nrun\nexit\nexit\n"
    assert m.current_menu == "suite"


def test_read_menu_takes_exactly_one_command_and_writes_no_exit():
    m = Macro()
    with m.program("discus"):
        with m.menu("read"):
            m.command("cell", "model.cll", 2, 2, 2)
        m.command("save", "out.stru")
    text = m.render()
    assert text == "discus\nread\ncell model.cll,2,2,2\nsave out.stru\nexit\n"


def test_read_menu_with_zero_commands_is_refused():
    m = Macro()
    with pytest.raises(MacroError, match="exactly one"):
        with m.program("discus"):
            with m.menu("read"):
                pass


def test_read_menu_with_two_commands_is_refused():
    # read leaves the menu after its one command; a second one is refused with a clear message
    m = Macro()
    with pytest.raises(MacroError, match="takes exactly one command"):
        with m.program("discus"):
            with m.menu("read"):
                m.command("cell", "a.cll", 1, 1, 1)
                m.command("cell", "b.cll", 1, 1, 1)


def test_command_refuses_a_word_that_opens_a_menu():
    m = Macro()
    with m.program("discus"):
        with pytest.raises(MacroError, match="opens the menu"):
            m.command("fourier")


def test_menu_refuses_a_word_that_does_not_open_a_menu():
    m = Macro()
    with m.program("discus"):
        with pytest.raises(MacroError, match="does not open a menu"):
            with m.menu("save", "out.stru"):
                pass


def test_command_refuses_unknown_words():
    m = Macro()
    with m.program("discus"):
        with pytest.raises(MacroError, match="not a command"):
            m.command("definitely_not_a_command")


@pytest.mark.parametrize("word", ["set lots", "set,lots", "", "a\tb"])
def test_command_refuses_words_with_spaces_or_commas(word):
    m = Macro()
    with m.program("discus"):
        with pytest.raises(MacroError, match="invalid command word"):
            m.command(word)


def test_save_is_a_plain_command_at_discus_level():
    m = Macro()
    with m.program("discus"):
        m.command("save", "out.stru")
        assert m.current_menu == "discus"
    assert m.render() == "discus\nsave out.stru\nexit\n"


def test_options_become_name_value_with_trailing_underscore_dropped():
    m = Macro()
    with m.program("discus"):
        with m.menu("fourier"):
            m.command("set", technique="nufft", global_=1)
    text = m.render()
    assert "set technique:nufft,global:1" in text.splitlines()


def test_positional_and_option_parameters_are_both_formatted():
    m = Macro()
    with m.program("discus"):
        with m.menu("fourier"):
            m.command("set", "a", 0.1 + 0.2, x=2.5)
    assert "set a,0.3,x:2.5" in m.render().splitlines()


def test_render_raises_when_macro_ends_inside_a_menu():
    m = Macro()
    m.program("discus").__enter__()
    with pytest.raises(MacroError, match="ends inside"):
        m.render()


def test_top_level_exit_is_refused():
    m = Macro()
    with pytest.raises(MacroError, match="top level"):
        m.command("exit")


def test_unknown_start_menu_is_refused():
    with pytest.raises(MacroError, match="unknown start menu"):
        Macro(start="nowhere")


# --- raw and comment lines --------------------------------------------------

def test_raw_lines_are_accepted_when_known():
    m = Macro()
    with m.program("discus"):
        m.raw("set foo,bar")
    assert m.render() == "discus\nset foo,bar\nexit\n"


def test_raw_unknown_command_makes_render_raise():
    m = Macro()
    with m.program("discus"):
        m.raw("nonsense_command 1,2")
    with pytest.raises(MacroError, match="unknown|not"):
        m.render()


def test_raw_refuses_multiple_lines():
    m = Macro()
    with pytest.raises(MacroError, match="one line"):
        m.raw("set a,1\nset b,2")


def test_comment_lines_start_with_hash_space():
    m = Macro("title line\nsecond line")
    m.comment("third")
    assert m.render() == "# title line\n# second line\n# third\n"


# --- catalogue conformance --------------------------------------------------

def _rendered_examples():
    yield Macro("pdf-like")
    m = Macro("hk0")
    with m.program("discus"):
        with m.menu("read"):
            m.command("cell", "gaas.cll", 4, 4, 4)
        with m.menu("fourier"):
            m.command("xray").command("ll", -2, -2, 0).command("run")
    yield m
    m = Macro("refine")
    with m.program("refine"):
        m.command("run")
    yield m
    m = Macro("rmc")
    with m.program("discus"):
        with m.menu("rmc"):
            m.command("set", "lots", "off")
    yield m


def test_every_rendered_macro_passes_the_catalogue_check():
    cat = load_catalogue()
    for example in _rendered_examples():
        text = example.render()
        assert cat.check(text) == []


def test_ksave_run_leaves_its_menu_so_no_exit_is_written():
    m = Macro()
    with m.program("kuplot"):
        m.command("function", "sin(r[0])", 0, 6, 0.5)
        with m.menu("ksave", 1):
            m.command("outfile", "sine.xy").command("format", "xy").command("run")
        m.command("function", "cos(r[0])", 0, 6, 0.5)  # back at the KUPLOT prompt
    assert m.render().splitlines() == [
        "kuplot", "function sin(r[0]),0,6,0.5", "ksave 1", "outfile sine.xy", "format xy", "run",
        "function cos(r[0]),0,6,0.5", "exit"]


def test_nothing_may_follow_a_leaving_command_in_its_block():
    m = Macro()
    with pytest.raises(MacroError, match="already left"):
        with m.program("kuplot"):
            with m.menu("ksave", 1):
                m.command("run")
                m.command("format", "xy")


def test_finished_ends_a_program_like_exit():
    m = Macro()
    with m.program("discus"):
        m.command("finished")
    assert m.render().splitlines() == ["discus", "finished"]
    with pytest.raises(MacroError):
        Macro().command("finished")  # at the suite level it ends the run
