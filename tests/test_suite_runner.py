"""Tests for SuiteRunner (discuspdf.suite_runner): prepare / execute with a fake engine.

No DISCUS_SUITE is run and discus_suite is never called; the fake engines are Python scripts."""
import hashlib
import json
import sys
from pathlib import Path

import pytest

from discuspdf.engine import Engine, EngineChoiceError
from discuspdf.macro import Macro
from discuspdf.suite_runner import SuiteRunner


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def plain_macro(output: str = "result.csv") -> Macro:
    m = Macro("suite runner test")
    with m.program("discus"):
        with m.menu("read"):
            m.command("cell", "model.cll", 1, 1, 1)
        m.command("save", output)
    return m


def refine_macro() -> Macro:
    m = Macro("needs the fork")
    with m.program("refine"):
        m.command("run")
    return m


def fake_engine(tmp_path: Path, *, outputs=None, stdout="", exit_code=0, label="stock",
                name="fake_discus.py", marker: Path | None = None) -> Engine:
    """An engine that writes ``outputs`` into its run directory, prints ``stdout``, touches
    ``marker`` (so tests can see whether it ran) and exits with ``exit_code``."""
    body = ("import sys\n"
            "from pathlib import Path\n"
            f"for name, text in {dict(outputs or {})!r}.items():\n"
            "    Path(name).write_text(text)\n"
            f"sys.stdout.write({stdout!r})\n"
            + (f"Path({str(marker)!r}).write_text('ran')\n" if marker else "")
            + f"sys.exit({exit_code})\n")
    script = tmp_path / name
    script.write_text(body, encoding="utf-8")
    return Engine(launcher=(sys.executable, str(script)), label=label)


def make_runner(tmp_path: Path, engine: Engine, fork: Engine | None = None) -> SuiteRunner:
    return SuiteRunner(runs_directory=tmp_path / "runs", engine=engine, fork=fork)


def model_file(tmp_path: Path) -> Path:
    path = tmp_path / "model.cll"
    path.write_text("cell data\n", encoding="utf-8")
    return path


def marker_engine(tmp_path: Path, **kwargs):
    marker = tmp_path / "engine_ran.flag"
    return fake_engine(tmp_path, marker=marker, **kwargs), marker


# --- prepare ------------------------------------------------------------------------

def test_prepare_writes_macro_copies_inputs_and_records_hashes(tmp_path):
    engine, marker = marker_engine(tmp_path)
    runner = make_runner(tmp_path, engine)
    src = model_file(tmp_path)
    macro = plain_macro()
    run = runner.prepare(macro, {"model.cll": src}, outputs=["result.csv"])

    text = macro.render().encode("utf-8")
    assert run.macro.read_bytes() == text
    assert b"\r" not in run.macro.read_bytes()
    assert (run.directory / "model.cll").read_text() == "cell data\n"

    record = run.record
    assert record["engine"]["label"] == "stock"
    assert record["fork_reasons"] == []
    assert record["inputs"] == ["model.cll"]
    assert record["outputs"] == ["result.csv"]
    assert record["sha256"]["run.mac"] == sha(text)
    assert record["sha256"]["model.cll"] == sha(src.read_bytes())

    assert run.status() == {"status": "prepared"}
    assert not marker.exists(), "prepare must not start the engine"


def test_prepare_accepts_a_list_of_inputs_keeping_their_names(tmp_path):
    engine, _ = marker_engine(tmp_path)
    runner = make_runner(tmp_path, engine)
    src = model_file(tmp_path)
    run = runner.prepare(plain_macro(), [src], outputs=["result.csv"])
    assert run.record["inputs"] == ["model.cll"]
    assert (run.directory / "model.cll").is_file()


def test_prepare_accepts_macro_text(tmp_path):
    engine, _ = marker_engine(tmp_path)
    runner = make_runner(tmp_path, engine)
    run = runner.prepare("discus\nread\ncell model.cll,1,1,1\nexit\nexit\n", {"model.cll": model_file(tmp_path)})
    assert run.macro.read_text().endswith("exit\n")


def test_unknown_command_in_macro_text_is_refused(tmp_path):
    engine, marker = marker_engine(tmp_path)
    runner = make_runner(tmp_path, engine)
    with pytest.raises(ValueError, match="does not know"):
        runner.prepare("discus\nno_such_command 1,2\nexit\n")
    assert not marker.exists()


def test_missing_input_file_is_refused(tmp_path):
    runner = make_runner(tmp_path, fake_engine(tmp_path))
    with pytest.raises(FileNotFoundError):
        runner.prepare(plain_macro(), {"model.cll": tmp_path / "absent.cll"}, outputs=["result.csv"])


@pytest.mark.parametrize("inputs, outputs", [
    ({"sub/model.cll": "src"}, []),          # path separator in input name
    ({"my model.cll": "src"}, []),           # space in input name
    ({}, ["out dir.csv"]),                   # space in output name
    ({}, ["a/b.csv"]),                       # path separator in output name
    ({}, ["..", ]),                          # dot-dot
    ({"run.mac": "src"}, []),                # reserved input name
    ({"status.json": "src"}, []),            # reserved input name
    ({}, ["run.log"]),                       # reserved output name
    ({}, ["run.json"]),                      # reserved output name
    ({}, ["result.csv", "result.csv"]),      # duplicate outputs
])
def test_bad_file_names_are_refused(tmp_path, inputs, outputs):
    runner = make_runner(tmp_path, fake_engine(tmp_path))
    src = model_file(tmp_path)
    staged = {name: (src if value == "src" else value) for name, value in inputs.items()}
    with pytest.raises(ValueError):
        runner.prepare(plain_macro(), staged, outputs=outputs)


# --- engine choice ------------------------------------------------------------------------

def test_refine_macro_without_fork_is_refused(tmp_path):
    runner = make_runner(tmp_path, fake_engine(tmp_path))
    with pytest.raises(EngineChoiceError, match="needs the fork"):
        runner.prepare(refine_macro())


def test_refine_macro_with_fork_is_recorded_as_fork(tmp_path):
    fork = fake_engine(tmp_path, label="fork", name="fork_discus.py")
    runner = make_runner(tmp_path, fake_engine(tmp_path), fork=fork)
    run = runner.prepare(refine_macro())
    record = run.record
    assert record["engine"]["label"] == "fork"
    assert record["engine"]["launcher"] == list(fork.launcher)
    assert record["fork_reasons"] and "REFINE" in record["fork_reasons"][0]


def test_engine_stock_forces_stock_even_for_a_fork_macro(tmp_path):
    fork = fake_engine(tmp_path, label="fork", name="fork_discus.py")
    runner = make_runner(tmp_path, fake_engine(tmp_path), fork=fork)
    run = runner.prepare(refine_macro(), engine="stock")
    assert run.record["engine"]["label"] == "stock"


def test_engine_fork_without_a_fork_is_refused(tmp_path):
    runner = make_runner(tmp_path, fake_engine(tmp_path))
    with pytest.raises(EngineChoiceError):
        runner.prepare(plain_macro(), engine="fork")


def test_unknown_engine_choice_is_refused(tmp_path):
    runner = make_runner(tmp_path, fake_engine(tmp_path))
    with pytest.raises(ValueError, match="'stock', 'fork' or None"):
        runner.prepare(plain_macro(), engine="beta")


# --- execute ----------------------------------------------------------------------------

def test_execute_completes_when_all_outputs_are_written(tmp_path):
    engine = fake_engine(tmp_path, outputs={"result.csv": "r,G\n0,1\n"}, stdout="done\n")
    runner = make_runner(tmp_path, engine)
    run = runner.prepare(plain_macro(), {"model.cll": model_file(tmp_path)}, outputs=["result.csv"])
    result = runner.execute(run)
    assert result is run
    status = run.status()
    assert status["status"] == "complete"
    assert status["returncode"] == 0
    assert status["outputs"] == {"result.csv": sha(b"r,G\n0,1\n")}
    assert "done" in run.log.read_text()
    assert run.output("result.csv").read_text() == "r,G\n0,1\n"


def test_run_is_prepare_then_execute(tmp_path):
    engine = fake_engine(tmp_path, outputs={"result.csv": "x\n"})
    runner = make_runner(tmp_path, engine)
    run = runner.run(plain_macro(), outputs=["result.csv"])
    assert run.status()["status"] == "complete"


def test_nonzero_exit_fails(tmp_path):
    engine = fake_engine(tmp_path, outputs={"result.csv": "x\n"}, exit_code=3)
    runner = make_runner(tmp_path, engine)
    run = runner.prepare(plain_macro(), outputs=["result.csv"])
    with pytest.raises(RuntimeError, match="DISCUS_SUITE run failed"):
        runner.execute(run)
    status = run.status()
    assert status["status"] == "failed"
    assert status["reason"] == "nonzero exit"
    assert status["returncode"] == 3


def test_engine_error_in_log_fails_even_with_exit_zero(tmp_path):
    engine = fake_engine(tmp_path, outputs={"result.csv": "x\n"}, stdout=" ***FOUR*** boom\n")
    runner = make_runner(tmp_path, engine)
    run = runner.prepare(plain_macro(), outputs=["result.csv"])
    with pytest.raises(RuntimeError):
        runner.execute(run)
    status = run.status()
    assert status["status"] == "failed"
    assert "engine error in log" in status["reason"]
    assert status["engine_error"] == "***FOUR***"


def test_missing_output_fails(tmp_path):
    engine = fake_engine(tmp_path, outputs={})
    runner = make_runner(tmp_path, engine)
    run = runner.prepare(plain_macro(), outputs=["result.csv"])
    with pytest.raises(RuntimeError):
        runner.execute(run)
    status = run.status()
    assert status["status"] == "failed"
    assert status["reason"] == "missing or empty output: result.csv"


def test_empty_output_fails(tmp_path):
    engine = fake_engine(tmp_path, outputs={"result.csv": ""})
    runner = make_runner(tmp_path, engine)
    run = runner.prepare(plain_macro(), outputs=["result.csv"])
    with pytest.raises(RuntimeError):
        runner.execute(run)
    assert run.status()["status"] == "failed"
    assert "result.csv" in run.status()["reason"]


def test_modified_macro_is_refused_before_the_engine_starts(tmp_path):
    engine, marker = marker_engine(tmp_path, outputs={"result.csv": "x\n"})
    runner = make_runner(tmp_path, engine)
    run = runner.prepare(plain_macro(), outputs=["result.csv"])
    run.macro.write_text(run.macro.read_text() + "# edited\n")
    with pytest.raises(RuntimeError):
        runner.execute(run)
    status = run.status()
    assert status["status"] == "failed"
    assert status["engine_started"] is False
    assert "run.mac" in status["reason"]
    assert not marker.exists()


def test_modified_input_is_refused_before_the_engine_starts(tmp_path):
    engine, marker = marker_engine(tmp_path, outputs={"result.csv": "x\n"})
    runner = make_runner(tmp_path, engine)
    run = runner.prepare(plain_macro(), {"model.cll": model_file(tmp_path)}, outputs=["result.csv"])
    (run.directory / "model.cll").write_text("different cell\n")
    with pytest.raises(RuntimeError):
        runner.execute(run)
    status = run.status()
    assert status["status"] == "failed"
    assert status["engine_started"] is False
    assert "model.cll" in status["reason"]
    assert not marker.exists()


def test_existing_output_is_refused_before_the_engine_starts(tmp_path):
    engine, marker = marker_engine(tmp_path, outputs={"result.csv": "x\n"})
    runner = make_runner(tmp_path, engine)
    run = runner.prepare(plain_macro(), outputs=["result.csv"])
    (run.directory / "result.csv").write_text("stale\n")
    with pytest.raises(RuntimeError):
        runner.execute(run)
    status = run.status()
    assert status["status"] == "failed"
    assert status["engine_started"] is False
    assert "exists already" in status["reason"]
    assert not marker.exists()


def test_output_lookup_only_knows_declared_names(tmp_path):
    engine = fake_engine(tmp_path, outputs={"result.csv": "x\n"})
    runner = make_runner(tmp_path, engine)
    run = runner.run(plain_macro(), outputs=["result.csv"])
    assert run.output("result.csv") == run.directory / "result.csv"
    with pytest.raises(KeyError, match="not a declared output"):
        run.output("run.log")
    with pytest.raises(KeyError):
        run.output("unknown.csv")


# --- fixes after the first test pass --------------------------------------------------

def test_timeout_is_recorded_as_started_and_stopped(tmp_path):
    script = tmp_path / "sleepy.py"
    script.write_text("import time\ntime.sleep(30)\n", encoding="utf-8")
    runner = make_runner(tmp_path, Engine((sys.executable, str(script)), timeout_seconds=1))
    run = runner.prepare(plain_macro(), {"model.cll": model_file(tmp_path)}, ["result.csv"])
    with pytest.raises(RuntimeError):
        runner.execute(run)
    status = run.status()
    assert status["status"] == "failed" and status["engine_started"] is True
    assert "timed out" in status["reason"]


def test_forced_stock_is_flagged(tmp_path):
    engine = fake_engine(tmp_path, outputs={"out.txt": "x"})
    runner = make_runner(tmp_path, engine)
    run = runner.run(refine_macro(), (), ["out.txt"], engine="stock")
    assert run.record["stock_despite_fork_reasons"] is True
    assert run.status()["warnings"] and "fork" in run.status()["warnings"][0]
    plain = make_runner(tmp_path, fake_engine(tmp_path, outputs={"a.stru": "x"}, name="f2.py")).run(
        plain_macro("a.stru"), {"model.cll": model_file(tmp_path)}, ["a.stru"])
    assert plain.record["stock_despite_fork_reasons"] is False
    assert plain.status()["warnings"] == []


def test_tagged_warning_is_not_an_error(tmp_path):
    engine = fake_engine(tmp_path, outputs={"a.stru": "x"}, stdout=" ***MAC *** Warning in macro line 3\n")
    run = make_runner(tmp_path, engine).run(plain_macro("a.stru"), {"model.cll": model_file(tmp_path)}, ["a.stru"])
    assert run.status()["status"] == "complete"
    assert any("Warning in macro" in w for w in run.status()["warnings"])


def test_tagged_error_still_fails(tmp_path):
    engine = fake_engine(tmp_path, outputs={"a.stru": "x"}, stdout=" ***MAC *** Error in macro line 3\n")
    runner = make_runner(tmp_path, engine)
    run = runner.prepare(plain_macro("a.stru"), {"model.cll": model_file(tmp_path)}, ["a.stru"])
    with pytest.raises(RuntimeError):
        runner.execute(run)


def test_finished_run_is_never_re_executed(tmp_path):
    engine = fake_engine(tmp_path, outputs={"a.stru": "x"})
    runner = make_runner(tmp_path, engine)
    run = runner.run(plain_macro("a.stru"), {"model.cll": model_file(tmp_path)}, ["a.stru"])
    (run.directory / "a.stru").unlink()
    with pytest.raises(RuntimeError, match="not 'prepared'"):
        runner.execute(run)
    assert run.status()["status"] == "complete"


def test_duplicate_input_basenames_are_refused(tmp_path):
    a = tmp_path / "a" / "model.cll"; b = tmp_path / "b" / "model.cll"
    for p in (a, b):
        p.parent.mkdir(); p.write_text("x")
    with pytest.raises(ValueError, match="share a file name"):
        make_runner(tmp_path, fake_engine(tmp_path)).prepare(plain_macro(), [a, b], [])
