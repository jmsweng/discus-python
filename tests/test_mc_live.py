"""Live DISCUS 6.23.01 tests for MCConfig + MCRunner (opt-in: DISCUS_LIVE=1).

Each test generates a macro with MCConfig, runs it with MCRunner and compares the result with the committed
probe evidence (benchmarks/discus_6_23_01/d_mc_probes/evidence/run1), which was produced from hand-written
macros with the same seed, vectors (same order), moves and targets. The generated macro adds `set vec, reset`
and `set output, feed:on, final:on`; these consume no random numbers, so the final structure must be identical.

Environment: DISCUS_LIVE=1; DISCUS_LIVE_WORK (no spaces; default C:\\tmp\\mclive on Windows, /tmp/mclive);
DISCUS_WSL_DISTRO (default Ubuntu-22.04).
"""
from __future__ import annotations

import json
import os
import platform
from pathlib import Path

import numpy as np
import pytest

from discuspdf.discus_io import read_stru
from discuspdf.fourier_macro import Seed
from discuspdf.mc import (DisplacementTarget, LennardJonesTarget, MCConfig, MCRunner, Neighborhood, OccupationalTarget, Shift,
                          SpringTarget, Swap, SwapDisplacement, Vector, parse_mmc_log)

pytestmark = pytest.mark.discus_live
if os.environ.get("DISCUS_LIVE") != "1":
    pytest.skip("live DISCUS test; set DISCUS_LIVE=1 to enable", allow_module_level=True)

ROOT = Path(__file__).resolve().parents[1]
EV = ROOT / "benchmarks/discus_6_23_01/d_mc_probes/evidence/run1"
WORK = Path(os.environ.get("DISCUS_LIVE_WORK", r"C:\tmp\mclive" if platform.system() == "Windows" else "/tmp/mclive"))
SEED = Seed((11, 22, 33, 44, 55, 66, 77, 88))
# same vector order as the probe macros: +x, -x, +y, -y, +z, -z
SIX = Neighborhood(tuple(Vector(1, 1, t) for t in [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]))
TITLES = {"cuau.stru": "D MMC probe Cu0.5Au0.5 sc 16^3", "cu_disp.stru": "D MMC probe Cu sc 12^3 displaced"}


def launcher():
    if platform.system() == "Windows":
        return ("wsl.exe", "-d", os.environ.get("DISCUS_WSL_DISTRO", "Ubuntu-22.04"), "--",
                "discus_suite", "-noautorun", "-macro")
    return ("discus_suite", "-noautorun", "-macro")


CASES = {
    "r1_occ_corr": ("cuau.stru", MCConfig((SIX,), (OccupationalTarget(0, "Cu", "Au", correlation=-0.3),), (Swap(),),
                                          SEED, moves_per_atom=10, feedback_per_atom=1)),
    "r3b_occ_ener_cold": ("cuau.stru", MCConfig((SIX,), (OccupationalTarget(0, "Cu", "Au", energy=0.5),), (Swap(),),
                                                SEED, moves_per_atom=10, feedback_per_atom=1)),
    "r4_shift_spring": ("cu_disp.stru", MCConfig((SIX,), (SpringTarget(0, "Cu", "Cu", 3.0, 1.0),),
                                                 (Shift({"Cu": (0.01, 0.01, 0.01)}),), SEED, 10, 1)),
    "r6_swdisp_cd": ("cu_disp.stru", MCConfig((SIX,), (DisplacementTarget(0, "Cu", "Cu", (1, 0, 0), correlation=0.3),),
                                              (SwapDisplacement(),), SEED, 10, 1)),
}


@pytest.mark.parametrize("name", sorted(CASES))
def test_generated_macro_reproduces_probe(name):
    if " " in str(WORK):
        pytest.fail(f"DISCUS_LIVE_WORK must not contain spaces: {WORK}")
    stru, cfg = CASES[name]
    structure = read_stru(EV / "work" / stru).structure
    runner = MCRunner(WORK / name, launcher=launcher(), timeout_seconds=900)
    run = runner.prepare(cfg, structure, title=TITLES[stru])
    res = runner.execute(run)
    probe = read_stru(EV / "work" / f"{name}_out.stru", supercell=tuple(structure.supercell)).structure
    assert [s.upper() for s in res.structure.species] == list(probe.species)
    np.testing.assert_array_equal(res.structure.frac, probe.frac)
    probe_log = parse_mmc_log((EV / f"{name}.log").read_text(errors="replace"))
    assert res.log["final"] == probe_log["final"]
    assert [f["rows"] for f in res.log["feedback"]] == [f["rows"] for f in probe_log["feedback"]]
    (WORK / f"{name}_result.json").write_text(json.dumps({"directory": str(run.directory), "measure": res.measure(),
                                                          "final": res.log["final"]}, indent=1) + "\n")


# ---------------------------------------------------------------------------------- combinations not in the probe
def _b2(n=8, seed=5):
    from discuspdf.lattice import Lattice
    from discuspdf.structure import Structure
    from discuspdf.structure_ops import repeat
    cell = Structure(Lattice.from_parameters(3.0, 3.0, 3.0, 90, 90, 90), ["Cu", "Cu"], [[0, 0, 0], [0.5, 0.5, 0.5]],
                     model="configuration", b_iso=[0.5, 0.5])
    st = repeat(cell, (n, n, n)).structure
    sp = np.array(st.species, dtype=object)
    sp[np.random.default_rng(seed).permutation(st.n_atoms)[: st.n_atoms // 2]] = "Au"
    return st.replace(species=sp)


EIGHT = [(0, 0, 0), (-1, 0, 0), (0, -1, 0), (0, 0, -1), (-1, -1, 0), (-1, 0, -1), (0, -1, -1), (-1, -1, -1)]
COMBOS = {
    # two neighbourhoods on a two-site cell, a feedback and a fixed-energy occupational target, local swaps
    "c1_b2_two_targets": (_b2, "D MC combo B2", MCConfig(
        (Neighborhood(tuple(Vector(1, 2, t) for t in EIGHT)),
         Neighborhood((Vector(1, 1, (1, 0, 0)), Vector(1, 1, (0, 1, 0)), Vector(1, 1, (0, 0, 1)),
                       Vector(2, 2, (1, 0, 0)), Vector(2, 2, (0, 1, 0)), Vector(2, 2, (0, 0, 1))))),
        (OccupationalTarget(0, "Cu", "Au", correlation=-0.4), OccupationalTarget(1, "Au", "Cu", energy=-0.3)),
        (Swap(1.0, "local"),), SEED, moves_per_atom=8, feedback_per_atom=1)),
    # three energies on two neighbourhoods, two move types, kT != 1
    "c2_disp_three_energies": (lambda: read_stru(EV / "work" / "cu_disp.stru").structure, "D MC combo disp",
                               MCConfig((SIX, Neighborhood.from_translations([(1, 1, 0), (1, -1, 0), (1, 0, 1),
                                                                              (1, 0, -1), (0, 1, 1), (0, 1, -1)])),
                                        (SpringTarget(0, "Cu", "Cu", 3.0, 50.0),
                                         DisplacementTarget(0, "Cu", "Cu", (1, 0, 0), energy=0.5),
                                         LennardJonesTarget(1, "Cu", "Cu", 4.2426, 0.5)),
                                        (Shift({"Cu": (0.005, 0.005, 0.005)}, 0.6), SwapDisplacement(0.4, "slocal")),
                                        SEED, moves_per_atom=8, feedback_per_atom=2, temperature=0.5)),
}


@pytest.mark.parametrize("name", sorted(COMBOS))
def test_generated_combination_runs_and_measures(name):
    make, title, cfg = COMBOS[name]
    structure = make()
    runner = MCRunner(WORK / name, launcher=launcher(), timeout_seconds=900)
    run = runner.prepare(cfg, structure, title=title)
    res = runner.execute(run)                          # raises on any failure, incl. a missing target row
    counts_in = {s: int(np.sum(structure.species == s)) for s in set(structure.species)}
    counts_out = {s: int(np.sum(res.structure.species == s)) for s in set(res.structure.species)}
    assert counts_in == counts_out
    report = res.target_report()
    measured = res.measure()
    for rep, mea in zip(report, measured):
        assert abs(rep["achieved"] - mea["value"]) <= 0.0015, (rep, mea)
        if mea.get("std") is not None and rep["sigma"] is not None:
            assert abs(rep["sigma"] - mea["std"]) <= 0.0015, (rep, mea)
    (WORK / f"{name}_result.json").write_text(json.dumps({"directory": str(run.directory), "report": report,
                                                          "measure": measured}, indent=1, default=str) + "\n")
