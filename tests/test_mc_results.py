"""MMC log parsing and result loading on the committed DISCUS 6.23.01 probe evidence (no engine needed)."""
from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pytest

from discuspdf.discus_io import read_stru
from discuspdf.fourier_macro import Seed
from discuspdf.mc import MCConfig, MCRunner, Neighborhood, OccupationalTarget, Swap, parse_mmc_log

ROOT = Path(__file__).resolve().parents[1]
KIT = ROOT / "benchmarks/discus_6_23_01/d_mc_probes"
EV = KIT / "evidence/run1"
SIX = Neighborhood.from_translations([(1, 0, 0), (0, 1, 0), (0, 0, 1)])
R1 = MCConfig((SIX,), (OccupationalTarget(0, "Cu", "Au", correlation=-0.3),), (Swap(),),
              Seed((11, 22, 33, 44, 55, 66, 77, 88)), moves_per_atom=10, feedback_per_atom=1)


def _log(name):
    return parse_mmc_log((EV / f"{name}.log").read_text(errors="replace"))


def test_parse_occupational_run():
    L = _log("r1_occ_corr")
    assert len(L["initial"]) == 1 and L["initial"][0]["achieved"] == -0.003
    assert [f["gen"] for f in L["feedback"]] == list(range(4096, 40961, 4096))
    assert L["feedback"][0]["accepted"] == (793, 362, 286) and L["feedback"][0]["fraction"] == 0.1
    fin = L["final"]
    assert fin["gen"] == 40960 and fin["tries"] == 40961 and fin["accepted"] == (4955, 2521, 3987)
    assert fin["rows"] == [{"neighborhood": 1, "kind": "occupational", "a": "AU", "b": "CU", "target": -0.3,
                            "achieved": -0.306, "sigma": None, "diff": 0.006, "rel_diff": 0.021, "pairs": 24576}]
    assert L["energies"] == [{"energy": "Occupation correlation", "n_down": 4955, "mean_down": -1.437552,
                              "n_zero": 2521, "n_up": 3987, "mean_up": 1.486729}]
    assert L["errors"] == []


@pytest.mark.parametrize("name,kind,target,achieved,sigma,pairs", [
    ("r4_shift_spring", "spring", 3.0, 3.004, 0.096, 20736),
    ("r5_shift_lennard", "lennard_jones", 3.0, 3.004, 0.094, 20736),
    ("r6_swdisp_cd", "displacement", 0.3, 0.052, None, 10368),
])
def test_parse_other_energies(name, kind, target, achieved, sigma, pairs):
    row = _log(name)["final"]["rows"][0]
    assert (row["kind"], row["target"], row["achieved"], row["sigma"], row["pairs"]) == \
        (kind, target, achieved, sigma, pairs)


def test_parse_all_probe_logs_have_complete_tables():
    for log in sorted(EV.glob("*.log")):
        L = parse_mmc_log(log.read_text(errors="replace"))
        assert L["final"] is not None and L["final"]["rows"], log.name
        assert all(len(f["rows"]) == 1 for f in L["feedback"]), log.name
        assert L["errors"] == [], log.name


def test_parse_flags_engine_errors():
    L = parse_mmc_log(" ***RMC *** RMC did not find a valid pair\n")
    assert L["errors"] and L["final"] is None


def _staged_r1(tmp_path):
    inp = read_stru(EV / "work" / "cuau.stru").structure
    runner = MCRunner(runs_directory=tmp_path, launcher=(sys.executable, "-c", "pass"))
    run = runner.prepare(R1, inp)
    shutil.copy(EV / "r1_occ_corr.log", run.directory / "run.log")
    shutil.copy(EV / "work" / "r1_occ_corr_out.stru", run.directory / "mc_out.stru")
    return runner, run, inp


def test_load_probe_run_as_result(tmp_path):
    runner, run, inp = _staged_r1(tmp_path)
    res = runner.load(run, returncode=0)
    assert res.status == "complete"
    assert json.loads((run.directory / "status.json").read_text())["status"] == "complete"
    st = res.structure
    assert tuple(st.supercell) == (16, 16, 16) and st.n_atoms == 4096
    assert sorted(set(st.species)) == ["Au", "Cu"]                      # input spelling restored
    assert {s: int(np.sum(st.species == s)) for s in ("Au", "Cu")} == {"Au": 2048, "Cu": 2048}
    rep = res.target_report(tol=0.02)
    assert rep[0]["reached"] is True and rep[0]["achieved"] == -0.306
    assert res.target_report(tol=0.001)[0]["reached"] is False
    assert res.acceptance().shape == (10, 5) and res.history(0)[-1].tolist() == [40960, -0.3, -0.306]


def test_load_refuses_incomplete_or_failed_runs(tmp_path):
    runner, run, _ = _staged_r1(tmp_path)
    text = (run.directory / "run.log").read_text()
    (run.directory / "run.log").write_text(text.replace("Gen:      40960 try:      40961", "Gen:      40000 try:      40001"))
    with pytest.raises(RuntimeError, match="40000 moves done"):
        runner.load(run, returncode=0)
    assert json.loads((run.directory / "status.json").read_text())["status"] == "failed"
    (run.directory / "run.log").write_text(text + "\n ***RMC *** RMC did not find a valid pair\n")
    with pytest.raises(RuntimeError, match="engine error"):
        runner.load(run, returncode=0)
    (run.directory / "run.log").write_text(text)
    (run.directory / "mc_out.stru").unlink()
    with pytest.raises(RuntimeError, match="not written"):
        runner.load(run, returncode=0)


def test_execute_with_fake_engine(tmp_path):
    """execute() runs the launcher in the run directory, keeps the log and loads the output."""
    script = tmp_path / "fake_discus.py"
    script.write_text(
        "import shutil, sys\n"
        f"shutil.copy({str(EV / 'work' / 'r1_occ_corr_out.stru')!r}, 'mc_out.stru')\n"
        f"sys.stdout.write(open({str(EV / 'r1_occ_corr.log')!r}, errors='replace').read())\n")
    runner = MCRunner(runs_directory=tmp_path / "runs", launcher=(sys.executable, str(script)))
    run = runner.prepare(R1, read_stru(EV / "work" / "cuau.stru").structure)
    res = runner.execute(run)
    assert res.status == "complete" and (run.directory / "run.log").is_file()
    with pytest.raises(FileExistsError):
        runner.execute(run)


def test_prepare_stages_inputs_and_detects_tampering(tmp_path):
    runner, run, _ = _staged_r1(tmp_path)
    cfg = json.loads((run.directory / "mc_config.json").read_text())
    assert cfg["counts"] == [40960, 4096] and "heat-bath" in " ".join(cfg["cautions"])
    assert "set cyc, 40960" in (run.directory / "run.mac").read_text()
    (run.directory / "run.mac").write_text("tampered")
    (run.directory / "mc_out.stru").unlink()
    with pytest.raises(RuntimeError, match="modified"):
        runner.execute(run)


def _load_probe(tmp_path, name, cfg, stru):
    runner = MCRunner(runs_directory=tmp_path, launcher=(sys.executable, "-c", "pass"))
    run = runner.prepare(cfg, read_stru(EV / "work" / stru).structure)
    shutil.copy(EV / f"{name}.log", run.directory / "run.log")
    shutil.copy(EV / "work" / f"{name}_out.stru", run.directory / "mc_out.stru")
    return runner.load(run, returncode=0)


def test_independent_measurement_matches_discus_printed_values(tmp_path):
    from discuspdf.mc import DisplacementTarget, Shift, SpringTarget, SwapDisplacement
    seed = Seed((11, 22, 33, 44, 55, 66, 77, 88))
    occ = _load_probe(tmp_path, "r1_occ_corr", R1, "cuau.stru").measure()[0]
    assert occ["pairs"] == 24576 and abs(occ["value"] - (-0.306)) <= 5e-4
    spring = MCConfig((SIX,), (SpringTarget(0, "Cu", "Cu", 3.0, 1.0),), (Shift({"Cu": (0.01, 0.01, 0.01)}),), seed)
    m = _load_probe(tmp_path, "r4_shift_spring", spring, "cu_disp.stru").measure()[0]
    assert abs(m["value"] - 3.004) <= 5e-4 and abs(m["std"] - 0.096) <= 5e-4   # DISCUS prints mean, sigma
    cd = MCConfig((SIX,), (DisplacementTarget(0, "Cu", "Cu", (1, 0, 0), correlation=0.3),), (SwapDisplacement(),), seed)
    m = _load_probe(tmp_path, "r6_swdisp_cd", cd, "cu_disp.stru").measure()[0]
    assert abs(m["value"] - 0.052) <= 5e-4


def test_replica_summary_checks_inputs(tmp_path):
    from discuspdf.mc import replica_summary
    a = _load_probe(tmp_path, "r1_occ_corr", R1, "cuau.stru")
    with pytest.raises(ValueError, match="at least two"):
        replica_summary([a])
    with pytest.raises(ValueError, match="different seeds"):
        replica_summary([a, a])
