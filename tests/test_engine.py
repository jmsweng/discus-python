"""Tests for the shared engine (discuspdf.engine): launching, error detection, fork choice.

The engine is a fake Python script; no DISCUS_SUITE is run."""
import json
import subprocess
import sys
from pathlib import Path

import pytest

from discuspdf import Diffuse3DRunner, PDFRunner
from discuspdf.mc import MCRunner
from discuspdf.engine import (
    DEFAULT_LAUNCHER,
    Engine,
    EngineChoiceError,
    choose_engine,
    current_engine,
    fork_reasons,
)

REPO = Path(__file__).resolve().parents[1]
COMMITTED_MACRO = REPO / "examples" / "run_sample" / "run.mac"


def make_script(tmp_path: Path, body: str, name: str = "fake_discus.py") -> Path:
    script = tmp_path / name
    script.write_text(body, encoding="utf-8")
    return script


def printing_engine(tmp_path: Path, text: str, *, exit_code: int = 0) -> Engine:
    """A fake engine that prints ``text`` (taken from its environment) and exits with a code."""
    script = make_script(tmp_path, (
        "import os, sys\n"
        "sys.stdout.write(os.environ.get('FAKE_OUT', '') + '\\n')\n"
        f"sys.exit({exit_code})\n"))
    return Engine(launcher=(sys.executable, str(script)), env={"FAKE_OUT": text})


# --- validation and description ----------------------------------------------

def test_defaults_are_the_stock_engine():
    e = Engine()
    assert e.launcher == DEFAULT_LAUNCHER
    assert e.label == "stock"
    assert e.version == "6.23.01"
    assert e.timeout_seconds == 3600


def test_launcher_is_normalised_to_a_tuple_of_strings():
    e = Engine(launcher=["python3", Path("/x/y.py")])
    assert e.launcher == ("python3", "/x/y.py")


def test_command_appends_the_macro_name():
    assert Engine(("a", "-m")).command() == ["a", "-m", "run.mac"]
    assert Engine(("a",)).command("other.mac") == ["a", "other.mac"]


@pytest.mark.parametrize("kwargs, message", [
    ({"launcher": ()}, "launcher must not be empty"),
    ({"launcher": ("x",), "label": "beta"}, "label must be one of"),
    ({"launcher": ("x",), "timeout_seconds": 0}, "timeout_seconds must be positive"),
    ({"launcher": ("x",), "timeout_seconds": -5}, "timeout_seconds must be positive"),
])
def test_invalid_engines_are_refused(kwargs, message):
    with pytest.raises(ValueError, match=message):
        Engine(**kwargs)


def test_describe_is_json_serialisable_and_names_the_engine():
    e = Engine(("discus_suite", "-noautorun", "-macro"), "fork")
    d = e.describe()
    assert d == {"label": "fork", "version": "6.23.01", "launcher": ["discus_suite", "-noautorun", "-macro"]}
    assert json.loads(json.dumps(d)) == d


# --- from_env ------------------------------------------------------------------

def test_from_env_stock_reads_discuspdf_engine(monkeypatch):
    monkeypatch.setenv("DISCUSPDF_ENGINE", "/opt/d/discus_suite  -noautorun -macro")
    e = Engine.from_env()
    assert e.label == "stock"
    assert e.launcher == ("/opt/d/discus_suite", "-noautorun", "-macro")


def test_from_env_stock_defaults_when_unset(monkeypatch):
    monkeypatch.delenv("DISCUSPDF_ENGINE", raising=False)
    assert Engine.from_env().launcher == DEFAULT_LAUNCHER


def test_from_env_fork_reads_discuspdf_engine_fork(monkeypatch):
    monkeypatch.setenv("DISCUSPDF_ENGINE_FORK", "/opt/fork/bin/discus_suite -noautorun -macro")
    e = Engine.from_env("fork", timeout_seconds=120)
    assert e.label == "fork"
    assert e.launcher == ("/opt/fork/bin/discus_suite", "-noautorun", "-macro")
    assert e.timeout_seconds == 120


def test_from_env_fork_without_variable_is_a_lookup_error(monkeypatch):
    monkeypatch.delenv("DISCUSPDF_ENGINE_FORK", raising=False)
    with pytest.raises(LookupError, match="DISCUSPDF_ENGINE_FORK"):
        Engine.from_env("fork")


# --- Engine.run -------------------------------------------------------------------

def test_run_environment_and_stdin_and_cwd(tmp_path, monkeypatch):
    monkeypatch.setenv("DISPLAY", ":42")
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    script = make_script(tmp_path, (
        "import os, sys, json\n"
        "json.dump({'DISPLAY': os.environ.get('DISPLAY'), 'PGPLOT': os.environ.get('PGPLOT'),\n"
        "           'stdin': sys.stdin.read(), 'cwd': os.getcwd(), 'argv': sys.argv[1:]},\n"
        "          open('probe.json', 'w'))\n"))
    result = Engine(launcher=(sys.executable, str(script))).run(run_dir)
    probe = json.loads((run_dir / "probe.json").read_text())
    assert probe["DISPLAY"] is None
    assert probe["PGPLOT"] == "none"
    assert probe["stdin"] == "exit\n"
    assert Path(probe["cwd"]).resolve() == run_dir.resolve()
    assert probe["argv"] == ["run.mac"]
    assert result.returncode == 0 and result.ok


def test_run_log_is_stdout_separator_and_stderr(tmp_path):
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    script = make_script(tmp_path, "import sys\nsys.stdout.write('hello\\n')\nsys.stderr.write('warn\\n')\n")
    result = Engine(launcher=(sys.executable, str(script))).run(run_dir)
    expected = "hello\n" + "\n--- STDERR ---\n" + "warn\n"
    assert result.log == expected
    assert (run_dir / "run.log").read_text(encoding="utf-8") == expected


def test_log_name_none_writes_no_log(tmp_path):
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    engine = printing_engine(tmp_path, "fine")
    engine.run(run_dir, log_name=None)
    assert not (run_dir / "run.log").exists()


def test_returncode_is_propagated(tmp_path):
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    result = printing_engine(tmp_path, "all good", exit_code=7).run(run_dir)
    assert result.returncode == 7
    assert result.ok is False
    assert result.engine_error is None


@pytest.mark.parametrize("text, found", [
    (" ***FOUR*** something went wrong", "***FOUR***"),
    ("Unknown command 'frobnicate'", "Unknown command"),
    ("*****************", None),
    ("xdpyinfo: unable to open display", None),
    ("discus: finished, 3 atoms written", None),
])
def test_engine_error_detection(tmp_path, text, found):
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    result = printing_engine(tmp_path, text).run(run_dir)
    assert result.engine_error == found
    assert result.ok is (found is None)


def test_startup_banner_and_display_noise_are_not_errors(tmp_path):
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    banner = "\n".join(["*" * 17, "DISCUS_SUITE 6.23.01", "*" * 17,
                        "xdpyinfo: unable to open display", "pgxwin_server: cannot connect"])
    result = printing_engine(tmp_path, banner).run(run_dir)
    assert result.engine_error is None


def test_strip_ansi_finds_an_error_hidden_by_colour_codes(tmp_path):
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    coloured = "\x1b[31m***\x1b[1mFOUR\x1b[0m*** bad\x1b[0m"
    engine = printing_engine(tmp_path, coloured)
    assert engine.run(run_dir).engine_error is None
    stripped = engine.run(run_dir, strip_ansi=True)
    assert stripped.engine_error == "***FOUR***"


def test_error_pattern_none_gives_no_engine_error(tmp_path):
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    result = printing_engine(tmp_path, " ***FOUR*** something").run(run_dir, error_pattern=None)
    assert result.engine_error is None


def test_missing_launcher_raises_file_not_found(tmp_path):
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    engine = Engine(launcher=(str(tmp_path / "no_such_discus_suite"),))
    with pytest.raises(FileNotFoundError):
        engine.run(run_dir)


def test_sleeping_engine_raises_timeout_expired(tmp_path):
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    script = make_script(tmp_path, "import time\ntime.sleep(30)\n")
    engine = Engine(launcher=(sys.executable, str(script)))
    with pytest.raises(subprocess.TimeoutExpired):
        engine.run(run_dir, timeout_seconds=1)


# --- runners take an Engine and expose its launcher / timeout ------------------------

def test_runners_accept_an_engine(tmp_path):
    engine = Engine(("my_discus", "-macro"), timeout_seconds=99)
    runners = [
        PDFRunner(runs_directory=tmp_path / "a", engine=engine),
        Diffuse3DRunner(runs_directory=tmp_path / "b", engine=engine),
        MCRunner(runs_directory=tmp_path / "c", engine=engine),
    ]
    for runner in runners:
        assert runner.launcher == engine.launcher
        assert runner.timeout_seconds == engine.timeout_seconds
        assert current_engine(runner) == engine


@pytest.mark.parametrize("make", [
    lambda p: PDFRunner(runs_directory=p, engine="discus_suite"),
    lambda p: Diffuse3DRunner(runs_directory=p, engine=("discus_suite",)),
    lambda p: MCRunner(runs_directory=p, engine={"launcher": "x"}),
])
def test_runners_refuse_a_non_engine(tmp_path, make):
    with pytest.raises(TypeError, match="must be an Engine"):
        make(tmp_path)


def test_current_engine_honours_later_assignments(tmp_path):
    runner = PDFRunner(runs_directory=tmp_path / "runs")
    runner.launcher = ["/opt/alt/discus_suite", "-macro"]
    runner.timeout_seconds = 42
    engine = current_engine(runner)
    assert engine.launcher == ("/opt/alt/discus_suite", "-macro")
    assert engine.timeout_seconds == 42
    assert engine.label == "stock"


def test_current_engine_for_diffuse3d_and_unchanged_runner(tmp_path):
    runner = Diffuse3DRunner(runs_directory=tmp_path / "runs")
    assert current_engine(runner) is runner.engine
    runner.timeout_seconds = 7
    assert current_engine(runner).timeout_seconds == 7


# --- which macros need the fork ----------------------------------------------------

def test_committed_macro_needs_no_fork():
    assert fork_reasons(COMMITTED_MACRO.read_text()) == []


@pytest.mark.parametrize("text", [
    "refine\nrun\nexit\n",
])
def test_refine_run_needs_fork(text):
    reasons = fork_reasons(text)
    assert len(reasons) == 1 and "line 2" in reasons[0] and "REFINE 'run'" in reasons[0]


def test_refine_points_two_needs_fork():
    reasons = fork_reasons("refine\nnewparam a,points:2\nexit\n")
    assert len(reasons) == 1 and "points:2" in reasons[0]


def test_rmc_lots_needs_fork_but_lots_off_does_not():
    on = fork_reasons("discus\nrmc\nset lots,on\nexit\nexit\n")
    off = fork_reasons("discus\nrmc\nset lots,off\nexit\nexit\n")
    assert len(on) == 1 and "RMC with lots" in on[0]
    assert off == []


def test_nufft_needs_fork():
    reasons = fork_reasons("discus\nfourier\nset technique:nufft\nexit\nexit\n")
    assert len(reasons) == 1 and "NUFFT" in reasons[0]


def test_two_fourier_runs_after_temp_use_need_fork():
    two = fork_reasons("discus\nfourier\ntemp use\nrun\nrun\nexit\nexit\n")
    assert len(two) == 1 and "Debye-Waller" in two[0]


def test_one_fourier_run_after_temp_use_is_fine():
    assert fork_reasons("discus\nfourier\ntemp use\nrun\nexit\nexit\n") == []


def test_two_fourier_runs_after_temp_ignore_are_fine():
    assert fork_reasons("discus\nfourier\ntemp ignore\nrun\nrun\nexit\nexit\n") == []


# --- choose_engine ------------------------------------------------------------------

def test_choose_engine_returns_stock_for_committed_macro():
    stock, fork = Engine(), Engine(("fork_discus",), "fork")
    assert choose_engine(COMMITTED_MACRO.read_text(), stock, fork) is stock
    assert choose_engine(COMMITTED_MACRO.read_text(), stock) is stock


def test_choose_engine_returns_fork_when_needed():
    stock, fork = Engine(), Engine(("fork_discus",), "fork")
    assert choose_engine("refine\nrun\nexit\n", stock, fork) is fork


def test_choose_engine_raises_when_fork_needed_but_missing():
    with pytest.raises(EngineChoiceError, match="needs the fork engine"):
        choose_engine("refine\nrun\nexit\n", Engine())


def test_choose_engine_raises_when_given_fork_is_labelled_stock():
    with pytest.raises(EngineChoiceError, match="labelled 'stock'"):
        choose_engine("refine\nrun\nexit\n", Engine(), Engine(("x",), "stock"))


def test_run_rejects_non_positive_timeout(tmp_path):
    engine = Engine((sys.executable, "-c", "pass"))
    with pytest.raises(ValueError):
        engine.run(tmp_path, timeout_seconds=0)


def test_choose_engine_checks_labels_even_for_stock_macros():
    from discuspdf.engine import choose_engine
    stock = Engine(("x",), "stock")
    with pytest.raises(EngineChoiceError):
        choose_engine("discus\nexit\n", stock, Engine(("y",), "stock"))
    with pytest.raises(EngineChoiceError):
        choose_engine("discus\nexit\n", Engine(("y",), "fork"))
