"""Tests for discuspdf.discus_chem: the chem analysis macro, the log parser and run_chem.

The macro lines are compared with the committed probe macros (benchmarks/discus_6_23_01/c_probes/
inputs/*.mac), the parser is checked on the real DISCUS 6.23.01 logs of the probe run and against
predictions.json. No DISCUS_SUITE is run: run_chem uses a Python script that prints a real log.
"""
from __future__ import annotations

import json
import re
import shutil
import sys
from pathlib import Path

import pytest

from discuspdf.discus_chem import (ChemAnalysis, ChemDisplacement, ChemIncomplete, ChemOccupational,
                                   CorrelationVector, DistanceNeighbourhood, Neighbourhood,
                                   VectorNeighbourhood, parse_chem_log, run_chem)
from discuspdf.discus_ops import StructureJob
from discuspdf.engine import Engine
from discuspdf.fourier_macro import StructureSource
from discuspdf.suite_catalogue import load_catalogue
from discuspdf.suite_runner import SuiteRunner

PROBE = Path(__file__).resolve().parents[1] / "benchmarks" / "discus_6_23_01" / "c_probes"
INPUTS = PROBE / "inputs"
EVIDENCE = sorted((PROBE / "evidence").iterdir())[-1]
PREDICTIONS = json.loads((PROBE / "predictions.json").read_text(encoding="ascii"))
V = CorrelationVector


def _analysis_occ() -> ChemAnalysis:
    return ChemAnalysis(
        vectors=(V(1, 1, 2, (0, 0, 0)), V(2, 1, 2, (-1, -1, -1)), V(3, 1, 1, (1, 0, 0)),
                 V(4, 2, 2, (0, 1, 0)), V(5, 2, 1, (1, 0, 1))),
        neighbourhoods=tuple(Neighbourhood(VectorNeighbourhood(i)) for i in range(1, 6))
        + (Neighbourhood(DistanceNeighbourhood((0.5, 0.5, 0.5), 0.1)),
           Neighbourhood(DistanceNeighbourhood((1.0, 0.0, 0.0), 0.1))),
        calculations=(ChemOccupational("Au", "Cu"), ChemOccupational("Cu", "Au")),
    )


def _analysis_disp_sc() -> ChemAnalysis:
    dirs = [((1, 0, 0), (1, 0, 0)), ((1, 0, 0), (1, 0, 0)), ((1, 0, 0), (1, 0, 0)),
            ((1, 0, 0), (1, 0, 0)), ((1, 0, 0), (0, 1, 0)), ((0, 1, 0), (0, 1, 0))]
    return ChemAnalysis(
        vectors=(V(1, 1, 1, (1, 0, 0)), V(2, 1, 1, (2, 0, 0)), V(3, 1, 1, (3, 0, 0)),
                 V(4, 1, 1, (0, 1, 0)), V(5, 1, 1, (1, 0, 0)), V(6, 1, 1, (1, 0, 0))),
        neighbourhoods=tuple(Neighbourhood(VectorNeighbourhood(k), d) for k, d in enumerate(dirs, start=1)),
        calculations=(ChemDisplacement("Cu", "Cu"),),
    )


def _analysis_mixed_bcc() -> ChemAnalysis:
    return ChemAnalysis(
        vectors=(V(1, 1, 2, (0, 0, 0)), V(2, 1, 1, (0, 0, 1)), V(3, 1, 2, (0, 0, 0))),
        neighbourhoods=(Neighbourhood(VectorNeighbourhood(1), ((0, 0, 1), (0, 0, 1))),
                        Neighbourhood(VectorNeighbourhood(2), ((0, 0, 1), (0, 0, 1))),
                        Neighbourhood(VectorNeighbourhood(3), ((1, 1, 0), (1, 1, 0)))),
        calculations=(ChemDisplacement("Au", "Cu"), ChemDisplacement("Cu", "Cu"),
                      ChemDisplacement("all", "all")),
    )


ANALYSES = {"occ_bcc": _analysis_occ, "disp_sc": _analysis_disp_sc, "mixed_bcc": _analysis_mixed_bcc}


def _job(name: str, analysis: ChemAnalysis | None = None) -> StructureJob:
    analysis = analysis or ANALYSES[name]()
    return StructureJob(StructureSource.stru(f"{name}.stru"), [analysis], readback="none")


def _chem_block(lines: list[str]) -> list[str]:
    """From the 'chem' line to the 'exit' that leaves it."""
    start = lines.index("chem")
    end = lines.index("exit", start)
    return lines[start:end + 1]


def _norm(line: str) -> str:
    """Spaces around commas and the spelling of numbers (1 vs 1.0) do not change a DISCUS line."""
    parts = []
    for tok in re.split(r"\s*,\s*", line.strip()):
        try:
            parts.append(repr(float(tok)))
        except ValueError:
            parts.append(tok)
    return ",".join(parts)


def _probe_lines(name: str) -> list[str]:
    return (INPUTS / f"{name}.mac").read_text().splitlines()


# --- macro ------------------------------------------------------------------------------------

@pytest.mark.parametrize("name", sorted(ANALYSES))
def test_chem_block_matches_the_probe_macro_line_by_line(name):
    rendered = _job(name).macro().render().splitlines()
    mine = [_norm(x) for x in _chem_block(rendered)]
    probe = [_norm(x) for x in _chem_block(_probe_lines(name))]
    assert mine == probe


@pytest.mark.parametrize("name", sorted(ANALYSES))
def test_job_macro_is_accepted_by_the_catalogue(name):
    text = _job(name).macro().render()
    assert load_catalogue().check(text) == []


def test_job_macro_reads_the_crystal_then_chem_then_saves():
    lines = _job("occ_bcc").macro().render().splitlines()
    assert lines[1:4] == ["discus", "read", "structure occ_bcc.stru"]
    assert lines[4] == "chem"
    assert lines[-3:] == ["exit", "save result.stru", "exit"]


def test_occupational_and_displacement_lines_use_the_probe_commands():
    lines = _job("disp_sc").macro().render().splitlines()
    assert "corr disp,Cu,Cu" in lines
    assert "set neig,dir,1,0,0,1,0,0,number:1" in lines
    assert "show" in lines
    assert sum(1 for x in lines if x.startswith("eval res[")) == 7  # res[0..6] for 6 neighbourhoods


@pytest.mark.parametrize("make, message", [
    (lambda: ChemOccupational("Au", "Au"), "two different"),
    (lambda: ChemOccupational("all", "Cu"), "only allowed for displacement"),
    (lambda: DistanceNeighbourhood((0.5, 0.5, 0.5), 0.0), "sigma must be > 0"),
    (lambda: Neighbourhood(DistanceNeighbourhood((1, 0, 0), 0.1), ((1, 0, 0), (1, 0, 0))), "only written for a VectorNeighbourhood"),
    (lambda: ChemAnalysis(vectors=(V(1, 1, 2, (0, 0, 0)),), neighbourhoods=(Neighbourhood(VectorNeighbourhood(2)),),
                          calculations=(ChemOccupational("Au", "Cu"),)), "not defined"),
    (lambda: ChemAnalysis(vectors=(V(1, 1, 2, (0, 0, 0)), V(1, 2, 1, (0, 0, 0))), neighbourhoods=(Neighbourhood(VectorNeighbourhood(1)),),
                          calculations=(ChemOccupational("Au", "Cu"),)), "unique"),
    (lambda: ChemAnalysis(vectors=(V(1, 1, 2, (0, 0, 0)),), neighbourhoods=(Neighbourhood(VectorNeighbourhood(1)),),
                          calculations=()), "at least one calculation"),
    (lambda: ChemAnalysis(vectors=(V(1, 1, 2, (0, 0, 0)),), neighbourhoods=(),
                          calculations=(ChemOccupational("Au", "Cu"),)), "at least one neighbourhood"),
])
def test_invalid_definitions_are_refused(make, message):
    with pytest.raises((ValueError, TypeError), match=message):
        make()


# --- parser on the committed evidence ----------------------------------------------------------

def _log(name: str) -> str:
    return (EVIDENCE / f"{name}.log").read_text(errors="replace")


def test_parse_occ_bcc_matches_predictions_and_printed_table():
    results = parse_chem_log(_log("occ_bcc"))
    preds = PREDICTIONS["occ_bcc"]["neighborhoods"]
    assert len(results) == 14
    for block, key in ((0, "corr_AuCu"), (1, "corr_CuAu")):
        rows = [r for r in results if r.block == block]
        assert [r.neighbour for r in rows] == list(range(1, 8))
        for r, p in zip(rows, preds):
            assert r.kind == "occ"
            assert (r.species_a, r.species_b) == (("Au", "Cu") if block == 0 else ("Cu", "Au"))
            assert r.value_exact
            assert abs(r.value - p[key]) <= 1e-9
            assert abs(r.table_value - p[key]) <= 0.5e-4 + 1e-12  # printed with 4 decimals
            assert r.pairs == p.get("n_pairs_discus", p.get("n_pairs_ours"))
            aa, ab, bb = r.occupation_percent
            assert abs(aa + ab + bb - 100.0) <= 0.05


def test_parse_occ_bcc_definitions_come_from_the_show_table():
    rows = [r for r in parse_chem_log(_log("occ_bcc")) if r.block == 0]
    assert rows[0].definition == "vec 1: site 1 -> site 2, neig 0 0 0"
    assert rows[1].definition == "vec 2: site 1 -> site 2, neig -1 -1 -1"
    assert rows[5].definition == "dis 0.50 0.50 0.50 fsig 0.10"
    assert rows[6].definition == "dis 1.00 0.00 0.00 fsig 0.10"


def test_parse_disp_sc_matches_predictions():
    results = parse_chem_log(_log("disp_sc"))
    preds = PREDICTIONS["disp_sc"]["neighborhoods"]
    assert len(results) == 6
    for r, p in zip(results, preds):
        assert r.kind == "disp" and (r.species_a, r.species_b) == ("Cu", "Cu") and r.block == 0
        assert abs(r.value - p["Cu_Cu"]["discus_convention"]) <= 1e-9
        assert abs(r.table_value - p["Cu_Cu"]["discus_convention"]) <= 0.5e-4 + 1e-12
        assert r.pairs == p["Cu_Cu"]["n_pairs"]
        assert r.occupation_percent is None
    assert results[4].direction_a == (1.0, 0.0, 0.0) and results[4].direction_b == (0.0, 1.0, 0.0)
    assert results[5].definition == "vec 6: site 1 -> site 1, neig 1 0 0"


@pytest.mark.parametrize("block, key", [(0, "Au_Cu"), (1, "Cu_Cu"), (2, "all_all")])
def test_parse_mixed_bcc_matches_predictions(block, key):
    results = [r for r in parse_chem_log(_log("mixed_bcc")) if r.block == block]
    preds = PREDICTIONS["mixed_bcc"]["neighborhoods"]
    assert [r.neighbour for r in results] == [1, 2, 3]
    for r, p in zip(results, preds):
        assert abs(r.value - p[key]["discus_convention"]) <= 1e-9
        assert r.pairs == p[key]["n_pairs"]


def test_all_committed_logs_give_all_29_correlations():
    total = sum(len(parse_chem_log(_log(name))) for name in ("occ_bcc", "disp_sc", "mixed_bcc"))
    assert total == 29


def test_parse_ignores_ansi_colour_codes():
    plain = parse_chem_log(_log("disp_sc"))
    coloured = parse_chem_log("\x1b[31m" + _log("disp_sc").replace("\n", "\x1b[0m\n"))
    assert [r.value for r in coloured] == [r.value for r in plain]


def test_truncated_log_keeps_only_the_complete_blocks_it_has():
    text = _log("occ_bcc")
    results = parse_chem_log(text[:text.index("corr occ, Cu, Au")])
    assert {r.block for r in results} == {0} and len(results) == 7


# --- run_chem with a fake engine -------------------------------------------------------------

def fake_engine(tmp_path: Path, *, stdout: str, exit_code: int, outputs=None) -> Engine:
    """An engine that writes ``outputs`` into its run directory, prints ``stdout`` and exits."""
    body = ("import sys\n"
            "from pathlib import Path\n"
            f"for name, text in {dict(outputs or {})!r}.items():\n"
            "    Path(name).write_bytes(text.encode())\n"
            f"sys.stdout.write({stdout!r})\n"
            f"sys.exit({exit_code})\n")
    script = tmp_path / "fake_discus.py"
    script.write_text(body, encoding="utf-8")
    return Engine(launcher=(sys.executable, str(script)), label="stock")


def _runner(tmp_path: Path, engine: Engine) -> SuiteRunner:
    return SuiteRunner(runs_directory=tmp_path / "runs", engine=engine)


def _crystal(tmp_path: Path, name: str) -> dict:
    dst = tmp_path / f"{name}.stru"
    shutil.copyfile(INPUTS / f"{name}.stru", dst)
    return {f"{name}.stru": dst}


OCC_LOG = _log("occ_bcc")
DISP_LOG = _log("disp_sc")
RESULT_TEXT = "stru placeholder written by the fake engine\n"


def test_run_chem_returns_the_results_when_engine_aborts_after_printing(tmp_path):
    """occ_bcc on the probe run: exit 134 after every result was printed."""
    from discuspdf.discus_chem import CRASH_RETURN_CODE
    engine = fake_engine(tmp_path, stdout=OCC_LOG, exit_code=CRASH_RETURN_CODE,
                         outputs={"result.stru": RESULT_TEXT})
    outcome = run_chem(_job("occ_bcc"), _runner(tmp_path, engine), _crystal(tmp_path, "occ_bcc"))
    assert outcome.engine_crashed_at_exit is True
    assert len(outcome.results) == 14
    assert [r.value for r in outcome.results] == [r.value for r in parse_chem_log(OCC_LOG)]
    status = outcome.run.status()
    assert status["status"] == "failed" and status["returncode"] == 134
    assert outcome.run.output("result.stru").read_text() == RESULT_TEXT
    assert outcome.run.macro.read_text().count("corr occ,") == 2


def test_run_chem_raises_when_the_engine_aborts_without_results(tmp_path):
    engine = fake_engine(tmp_path, stdout="free(): invalid pointer\n", exit_code=134,
                         outputs={"result.stru": RESULT_TEXT})
    with pytest.raises(RuntimeError, match="DISCUS_SUITE run failed"):
        run_chem(_job("occ_bcc"), _runner(tmp_path, engine), _crystal(tmp_path, "occ_bcc"))


def test_run_chem_raises_on_abort_when_results_are_truncated(tmp_path):
    from discuspdf.discus_chem import CRASH_RETURN_CODE
    text = OCC_LOG[:OCC_LOG.index("corr occ, Cu, Au")]  # no signature, second block missing
    engine = fake_engine(tmp_path, stdout=text, exit_code=CRASH_RETURN_CODE,
                         outputs={"result.stru": RESULT_TEXT})
    with pytest.raises(RuntimeError) as info:
        run_chem(_job("occ_bcc"), _runner(tmp_path, engine), _crystal(tmp_path, "occ_bcc"))
    assert not isinstance(info.value, ChemIncomplete)


def test_run_chem_raises_on_abort_without_the_signature(tmp_path):
    from discuspdf.discus_chem import CRASH_RETURN_CODE
    text = OCC_LOG.replace("free(): invalid pointer", "some other abort")
    engine = fake_engine(tmp_path, stdout=text, exit_code=CRASH_RETURN_CODE,
                         outputs={"result.stru": RESULT_TEXT})
    with pytest.raises(RuntimeError, match="DISCUS_SUITE run failed"):
        run_chem(_job("occ_bcc"), _runner(tmp_path, engine), _crystal(tmp_path, "occ_bcc"))


def test_run_chem_raises_on_abort_when_the_saved_crystal_is_missing(tmp_path):
    from discuspdf.discus_chem import CRASH_RETURN_CODE
    engine = fake_engine(tmp_path, stdout=OCC_LOG, exit_code=CRASH_RETURN_CODE)
    with pytest.raises(RuntimeError, match="DISCUS_SUITE run failed"):
        run_chem(_job("occ_bcc"), _runner(tmp_path, engine), _crystal(tmp_path, "occ_bcc"))


def test_run_chem_clean_exit_returns_results_without_crash_flag(tmp_path):
    engine = fake_engine(tmp_path, stdout=DISP_LOG, exit_code=0, outputs={"result.stru": RESULT_TEXT})
    outcome = run_chem(_job("disp_sc"), _runner(tmp_path, engine), _crystal(tmp_path, "disp_sc"))
    assert outcome.engine_crashed_at_exit is False
    assert len(outcome.results) == 6
    assert outcome.run.status()["status"] == "complete"


def test_run_chem_clean_exit_with_missing_values_is_incomplete(tmp_path):
    text = DISP_LOG[:DISP_LOG.index("Value of res[5]")]  # rows 5 and 6 have no eval value
    engine = fake_engine(tmp_path, stdout=text, exit_code=0, outputs={"result.stru": RESULT_TEXT})
    with pytest.raises(ChemIncomplete, match="missing eval value"):
        run_chem(_job("disp_sc"), _runner(tmp_path, engine), _crystal(tmp_path, "disp_sc"))


def test_run_chem_needs_the_crystal_file(tmp_path):
    engine = fake_engine(tmp_path, stdout="", exit_code=0)
    with pytest.raises(ValueError, match="inputs must contain"):
        run_chem(_job("disp_sc"), _runner(tmp_path, engine), {})


def test_run_chem_needs_exactly_one_chem_step(tmp_path):
    engine = fake_engine(tmp_path, stdout="", exit_code=0)
    job = StructureJob(StructureSource.stru("disp_sc.stru"), [], readback="none")
    with pytest.raises(ValueError, match="exactly one ChemAnalysis"):
        run_chem(job, _runner(tmp_path, engine), _crystal(tmp_path, "disp_sc"))
