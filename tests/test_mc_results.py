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


# ------------------------------------------------------------------------------ review follow-ups
def _row(nb, kind, a, b, target, achieved, pairs=100, sigma=None):
    return {"neighborhood": nb, "kind": kind, "a": a, "b": b, "target": target, "achieved": achieved,
            "sigma": sigma, "diff": target - achieved, "rel_diff": 0.0, "pairs": pairs}


def test_rows_are_matched_not_taken_by_position():
    from discuspdf.mc import MCResult, SpringTarget, match_row
    nb2 = Neighborhood.from_translations([(1, 1, 0)])
    cfg = MCConfig((SIX, nb2), (SpringTarget(0, "Cu", "Cu", 3.0, 1.0), OccupationalTarget(1, "Cu", "Au", correlation=-0.3),
                                OccupationalTarget(0, "Cu", "Au", energy=0.2)), (Swap(),),
                   Seed((1, 2, 3, 4, 5, 6, 7, 8)))
    # DISCUS order: neighbourhood 1 (occupancy before Hooke), then neighbourhood 2
    rows = [_row(1, "occupational", "AU", "CU", 0.0, -0.11), _row(1, "spring", "CU", "CU", 3.0, 3.004, sigma=0.09),
            _row(2, "occupational", "AU", "CU", -0.3, -0.29)]
    log = {"initial": [], "feedback": [{"gen": 10, "tries": 10, "accepted": (1, 1, 1), "fraction": 1.0, "rows": rows}],
           "final": {"gen": 10, "tries": 11, "accepted": (1, 1, 1), "fraction": None, "rows": rows}, "energies": [], "errors": []}
    res = MCResult(Path("."), cfg, None, log, "complete", 10)
    rep = res.target_report(tol=0.02)
    assert [(r["kind"], r["neighborhood"], r["mode"], r["reached"]) for r in rep] == \
        [("spring", 1, "ENER", None), ("occupational", 2, "CORR", True), ("occupational", 1, "ENER", None)]
    assert res.history(1).tolist() == [[10, -0.3, -0.29]]
    with pytest.raises(ValueError, match="no printed row"):
        match_row(rows[:2], cfg.targets[1])


def test_measure_follows_discus_definitions_for_two_species_two_sites():
    """Brute force of DISCUS's formulas (mmc_basic.f90) on a 2-site Cu/Au crystal with displacements."""
    from discuspdf.lattice import Lattice
    from discuspdf.mc import DisplacementTarget, MCResult, SpringTarget, Vector
    from discuspdf.neighbors import label_sites_by_slot
    from discuspdf.structure import Structure
    from discuspdf.structure_ops import repeat
    rng = np.random.default_rng(1)
    cell = Structure(Lattice.from_parameters(3, 3.2, 3.4, 90, 90, 90), ["Cu", "Cu"], [[0, 0, 0], [0.5, 0.5, 0.5]],
                     model="configuration")
    n = np.array([5, 4, 4])
    st = label_sites_by_slot(repeat(cell, tuple(n)).structure)
    N = st.n_atoms
    sp = np.where(rng.random(N) < 0.35, "Au", "Cu").astype(object)
    u = rng.normal(0, 0.03, (N, 3))
    u[sp == "Au"] *= 3
    u += np.where(st.site_label == "site2", 1, 0)[:, None] * np.array([0.02, 0, 0])
    st = st.replace(species=sp, frac=st.frac + u / n, min_separation=None)
    nb = Neighborhood((Vector(1, 2, (0, 0, 0)), Vector(1, 1, (1, 0, 0))))
    nb2 = Neighborhood((Vector(1, 2, (0, 0, 0)), Vector(1, 1, (1, 0, 0))))
    targets = (OccupationalTarget(0, "Cu", "Au", energy=1.0), DisplacementTarget(0, "Cu", "Au", (1, 0, 0), energy=1.0),
               DisplacementTarget(1, "Au", "Au", (1, 1, 0), energy=1.0), SpringTarget(0, "Cu", "Au", 2.7, 1.0))
    cfg = MCConfig((nb, nb2), targets, (Swap(),), Seed((1, 2, 3, 4, 5, 6, 7, 8)))
    meas = MCResult(Path("."), cfg, st, {}, "complete", N).measure()
    # ---- brute force
    ci = st.cell_index
    site = np.array([int(s[4:]) for s in st.site_label])
    idx = {(tuple(c), s): r for r, (c, s) in enumerate(zip(map(tuple, ci), site))}
    pairs = []
    for r in range(N):
        for v in nb.vectors:
            if site[r] == v.origin_site:
                c = ci[r] + np.array(v.translation)
                off = np.floor_divide(c, n)
                pairs.append((r, idx[(tuple(c - off * n), v.target_site)], off))
    xf = st.frac * n
    ave = {s: (xf[site == s] - ci[site == s]).mean(0) for s in (1, 2)}
    G = np.diag([3, 3.2, 3.4])
    U = np.array([xf[r] - ci[r] - ave[site[r]] for r in range(N)]) @ G

    def occ(a, b):
        p = np.zeros((3, 3))
        g = lambda s: 0 if s == a else (1 if s == b else 2)   # noqa: E731
        for i, j, _ in pairs:
            p[g(sp[i]), g(sp[j])] += 1
        R = p.sum()
        return 0.5 * ((1 - p[1, 0] / p[1].sum() / (p[0].sum() / R)) + (1 - p[0, 1] / p[0].sum() / (p[1].sum() / R)))

    def cd(a, b, d):
        d = np.array(d, float) @ G
        d /= np.linalg.norm(d)
        s = np.zeros(3)
        for i, j, _ in pairs:
            if (sp[i], sp[j]) == (a, b) or (a != b and (sp[i], sp[j]) == (b, a)):
                x, y = U[i] @ d, U[j] @ d
                s += [x * y, x * x, y * y]
        return s[0] / np.sqrt(s[1] * s[2])

    ds = [np.linalg.norm((xf[j] + off * n - xf[i]) @ G) for i, j, off in pairs if {sp[i], sp[j]} == {"Cu", "Au"}]
    assert abs(meas[0]["value"] - occ("Cu", "Au")) < 1e-12
    assert abs(meas[1]["value"] - cd("Cu", "Au", (1, 0, 0))) < 1e-12
    assert abs(meas[2]["value"] - cd("Au", "Au", (1, 1, 0))) < 1e-12
    assert abs(meas[3]["value"] - np.mean(ds)) < 1e-12 and abs(meas[3]["std"] - np.std(ds)) < 1e-12
    assert meas[3]["pairs"] == len(ds)


def test_load_failure_modes(tmp_path):
    runner, run, _ = _staged_r1(tmp_path)
    text = (run.directory / "run.log").read_text()
    with pytest.raises(RuntimeError, match="exit code unknown"):
        runner.load(run)                                        # nothing recorded yet
    with pytest.raises(RuntimeError, match="exit code 3"):
        runner.load(run, returncode=3)
    with pytest.raises(RuntimeError, match="exit code 3"):
        runner.load(run)                                        # the recorded exit code is used
    (run.directory / "run.log").write_text(text.replace("Gen:      40960 try:      40961", "Gen:      40959 try:      40960"))
    with pytest.raises(RuntimeError, match="40959 moves done"):
        runner.load(run, returncode=0)
    (run.directory / "run.log").write_text(text + "\n ***MAC *** Warning in macro line 3\n")
    assert runner.load(run, returncode=0).status == "complete"   # warnings are not errors
    out = run.directory / "mc_out.stru"
    good = out.read_text()
    out.write_text(good.replace("\nAU ", "\nNI ", 1))
    with pytest.raises(RuntimeError, match="species"):
        runner.load(run, returncode=0)
    out.write_text("\n".join(good.splitlines()[:-1]) + "\n")
    with pytest.raises(RuntimeError, match="output not usable"):
        runner.load(run, returncode=0)
    assert json.loads((run.directory / "status.json").read_text())["status"] == "failed"
    out.write_text(good)
    assert runner.load(run, returncode=0).status == "complete"
    assert runner.load(run).status == "complete"                # recorded exit code 0 is reused


def test_parser_tolerates_overflowed_pair_count():
    line = "   1   Occupancy AU      CU           -0.300    -0.306               0.006     0.021   ********"
    L = parse_mmc_log(" Gen:  10 try:  10 acc: (g/n/b):  1 /  2 /  3  MC moves 1.0000\n" + line + "\n")
    assert L["feedback"][0]["rows"][0]["pairs"] is None and L["feedback"][0]["rows"][0]["achieved"] == -0.306
