"""MCConfig rendering, validation and persistence (Package D, Feature 15). No engine needed."""
from __future__ import annotations

import json

import numpy as np
import pytest

from discuspdf.fourier_macro import Seed
from discuspdf.lattice import Lattice
from discuspdf.mc import (BuckinghamTarget, Convergence, DisplacementTarget, LennardJonesTarget, MCConfig,
                          Neighborhood, OccupationalTarget, RepulsiveTarget, Shift, SpringTarget, Swap,
                          SwapDisplacement, Vector, buckingham_minimum, match_row, parse_mmc_log)
from discuspdf.structure import Structure
from discuspdf.structure_ops import repeat

SEED = Seed((11, 22, 33, 44, 55, 66, 77, 88))
NN = Neighborhood.from_translations([(1, 0, 0), (0, 1, 0), (0, 0, 1)])


def _cfg(**kw):
    base = dict(neighborhoods=(NN,), targets=(OccupationalTarget(0, "Cu", "Au", correlation=-0.3),),
                moves=(Swap(),), seed=SEED, moves_per_atom=10, feedback_per_atom=1)
    base.update(kw)
    return MCConfig(**base)


def _alloy(n=4):
    cell = Structure(Lattice.from_parameters(3, 3, 3, 90, 90, 90), ["Cu"], [[0, 0, 0]], model="configuration")
    st = repeat(cell, (n, n, n)).structure
    sp = np.array(["Cu", "Au"] * (st.n_atoms // 2), dtype=object)
    return st.replace(species=sp)


def test_symmetric_neighborhood_adds_reverse_vectors_once():
    nb = Neighborhood((Vector(1, 1, (1, 0, 0)), Vector(1, 1, (-1, 0, 0)), Vector(1, 2, (0, 0, 0))))
    assert len(nb.vectors) == 4
    assert Vector(2, 1, (0, 0, 0)) in nb.vectors
    assert len(Neighborhood((Vector(1, 1, (1, 0, 0)),), symmetric=False).vectors) == 1


def test_render_matches_probe_macro_lines():
    """The generated mmc block uses exactly the command forms of the engine probe (d_mc_probes/r1)."""
    text = _cfg().render("mc_in.stru", 4096)
    lines = text.splitlines()
    for expected in ["discus", "seed 11,22,33,44,55,66,77,88", "set parallel, useomp:off", "  stru mc_in.stru",
                     "set vec, reset", "set output, feed:on, final:on",
                     "mmc", "set vec, 1, 1, 1, 1, 0, 0", "set vec, 4, 1, 1, -1, 0, 0",
                     "set neig, vec, 1, 2, 3, 4, 5, 6, number:1", "set mode, 1.0, swchem, all",
                     "set target, 1, corr, Cu, Au, -0.3, 0.0, CORR", "set temp, 1.0", "set cyc, 40960",
                     "set feed, 4096", "run", "save mc_out.stru"]:
        assert expected in lines, expected
    assert lines.index("set neig, rese") > lines.index("set vec, 6, 1, 1, 0, 0, -1")
    assert lines.index("set target, 1, corr, Cu, Au, -0.3, 0.0, CORR") > lines.index("set neig, vec, 1, 2, 3, 4, 5, 6, number:1")
    assert text.endswith("exit\nexit\n") and "\r" not in text


def test_all_target_and_move_lines():
    nb2 = Neighborhood((Vector(1, 2, (0, 0, 0)),))
    cfg = MCConfig((NN, nb2),
                   (OccupationalTarget(0, "Cu", "Au", energy=0.5),
                    DisplacementTarget(0, "Au", "Cu", (1, 0, 0), correlation=0.3),
                    SpringTarget(1, "Cu", "Au", 2.6, 1.0),
                    LennardJonesTarget(1, "Au", "Au", 2.9, 1.0)),
                   (Swap(0.5, "site"), Shift({"Cu": (0.01, 0.01, 0.02)}, 0.3), SwapDisplacement(0.2)), SEED)
    lines = cfg.render("in.stru", 100).splitlines()
    for expected in ["set target, 1, corr, Cu, Au, 0.0, 0.5, ENER",
                     "set target, 1, cd, Au, Cu, 0.3, 0.0, CORR, dir:[1.0,0.0,0.0]",
                     "set target, 2, spring, Cu, Au, 2.6, 1.0",
                     "set target, 2, lennard, Au, Au, 2.9, 1.0, 12.0, 6.0",
                     "set mode, 0.5, swchem, site", "set mode, 0.3, shift, all", "set move, Cu, 0.01, 0.01, 0.02",
                     "set mode, 0.2, swdisp, all", "set neig, vec, 7, 8, number:2"]:
        assert expected in lines, expected


def test_roundtrip_dict():
    cfg = MCConfig((NN,), (SpringTarget(0, "Cu", "Cu", 3.0, 1.0),), (Shift({"Cu": (0.01, 0.01, 0.01)}),), SEED,
                   moves_per_atom=5, feedback_per_atom=0.5, temperature=0.8)
    assert MCConfig.from_dict(cfg.to_dict()) == cfg
    with pytest.raises(ValueError, match="version 1"):
        MCConfig.from_dict({**cfg.to_dict(), "version": 2})


def test_counts():
    assert _cfg().counts(4096) == (40960, 4096)
    assert _cfg(moves_per_atom=0.5, feedback_per_atom=0.5).counts(3) == (2, 2)


@pytest.mark.parametrize("kw,match", [
    (dict(seed=None), "seed"),
    (dict(neighborhoods=(NN, NN), targets=(OccupationalTarget(0, "Cu", "Au", correlation=-0.3),
                                           OccupationalTarget(1, "Cu", "Ni", correlation=0.2))), "at most one"),
    (dict(targets=(OccupationalTarget(0, "Cu", "Au", energy=1.0),
                   OccupationalTarget(0, "Cu", "Ni", energy=0.5))), "only the last occupational"),
    (dict(targets=(DisplacementTarget(0, "Cu", "Cu", (1, 0, 0), energy=1.0),
                   DisplacementTarget(0, "Au", "Au", (0, 1, 0), energy=1.0))), "one displacement direction"),
    (dict(targets=(OccupationalTarget(1, "Cu", "Au", energy=1.0),)), "out of range"),
    (dict(moves=(Swap(), Swap(0.5))), "once"),
    (dict(feedback_per_atom=20), "exceed"),
    (dict(temperature=0), "temperature"),
    (dict(output_file="a b.stru"), "output_file"),
])
def test_config_refusals(kw, match):
    with pytest.raises(ValueError, match=match):
        _cfg(**kw)


def test_target_refusals():
    with pytest.raises(ValueError, match="exactly one"):
        OccupationalTarget(0, "Cu", "Au")
    with pytest.raises(ValueError, match="exactly one"):
        OccupationalTarget(0, "Cu", "Au", correlation=0.1, energy=1.0)
    with pytest.raises(ValueError, match="different"):
        OccupationalTarget(0, "Cu", "CU", energy=1.0)
    with pytest.raises(ValueError, match="atom name"):
        OccupationalTarget(0, "Copper", "Au", energy=1.0)
    with pytest.raises(ValueError, match=r"\[-1, 1\]"):
        OccupationalTarget(0, "Cu", "Au", correlation=1.5)
    with pytest.raises(ValueError, match="n_rep > n_att"):
        LennardJonesTarget(0, "Cu", "Cu", 3.0, 1.0, 6, 12)
    with pytest.raises(ValueError, match="d0 = 0"):
        SpringTarget(0, "Cu", "Cu", 0.0, 1.0)
    with pytest.raises(ValueError, match="sigma"):
        Shift({"Cu": (0, 0, 0)})
    with pytest.raises(ValueError, match="pairing"):
        Swap(1.0, "neighbours")
    with pytest.raises(ValueError, match="itself"):
        Neighborhood((Vector(1, 1, (0, 0, 0)),))


def test_check_structure():
    st = _alloy()
    _cfg().check_structure(st)
    with pytest.raises(ValueError, match="not in the structure"):
        _cfg(targets=(OccupationalTarget(0, "Cu", "Ni", energy=1.0),)).check_structure(st)
    with pytest.raises(ValueError, match="site 2"):
        _cfg(neighborhoods=(Neighborhood((Vector(1, 2, (0, 0, 0)),)),)).check_structure(st)
    flat = repeat(Structure(Lattice.from_parameters(3, 3, 3, 90, 90, 90), ["Cu"], [[0, 0, 0]],
                            model="configuration"), (4, 4, 1)).structure
    with pytest.raises(ValueError, match="one cell"):
        _cfg(targets=(SpringTarget(0, "Cu", "Cu", 3.0, 1.0),)).check_structure(flat)


def test_cautions_name_known_engine_behaviour():
    c = " ".join(_cfg(targets=(OccupationalTarget(0, "Cu", "Au", correlation=0.05),)).cautions())
    assert "heat-bath" in c and "0.005" in c
    assert "J = 0" in " ".join(_cfg(targets=(OccupationalTarget(0, "Cu", "Au", correlation=0.0),)).cautions())


def test_numbers_are_written_without_exponents():
    cfg = _cfg(targets=(SpringTarget(0, "Cu", "Cu", 3.0, 1e-5),), moves=(Shift({"Cu": (1e-4, 0.01, 0.02)}),),
               temperature=2e-3)
    text = cfg.render("in.stru", 10)
    assert "set target, 1, spring, Cu, Cu, 3.0, 0.00001" in text
    assert "set move, Cu, 0.0001, 0.01, 0.02" in text and "set temp, 0.002" in text
    assert "e-" not in text and "e+" not in text


# ------------------------------------------------------------------------------------- phase 2b extensions
def test_repulsive_line_and_defaults():
    cfg = MCConfig((NN,), (RepulsiveTarget(0, "Cu", "Au", -0.5, 1.2, 0.8, 6),), (Swap(),), SEED)
    lines = cfg.render("in.stru", 100).splitlines()
    assert "set target, 1, repulsive, Cu, Au, 0.5, 1.2, 0.8, 6.0" in lines     # |depth| is written
    assert RepulsiveTarget(0, "Cu", "Cu", 1.0).line(1) == "set target, 1, repulsive, Cu, Cu, 1.0, 1.0, 0.0, 1.0"
    assert RepulsiveTarget(0, "Cu", "Cu", 1.0).mode == "ENER"
    assert any("averaged in with -|depth_min|" in c for c in cfg.cautions())


@pytest.mark.parametrize("kw,match", [
    (dict(depth=0.0), "non-zero"), (dict(depth=float("nan")), "finite"),
    (dict(b_scale=0.0), "b_scale"), (dict(c=-0.1), "c must"), (dict(c=float("inf")), "finite"),
    (dict(m=0.0), "m must"),
])
def test_repulsive_refusals(kw, match):
    args = dict(depth=1.0)
    args.update(kw)
    with pytest.raises(ValueError, match=match):
        RepulsiveTarget(0, "Cu", "Au", **args)


def test_buckingham_line_and_hand_checked_minimum():
    bt = BuckinghamTarget(0, "Cu", "Au", A=1000.0, rho=0.3, B=10.0)
    assert bt.line(1) == "set target, 1, bucking, Cu, Au, 0.0, 1000.0, 0.3, 10.0"
    assert bt.mode == "ENER"
    d, e, rmin, atmin = buckingham_minimum(1000.0, 0.3, 10.0)
    # analytic minimum: root of V'(x) = -(A/rho)·exp(-x/rho) + 6B/x^7 in (3, 6), by bisection (no scipy)
    dv = lambda x: -(1000.0 / 0.3) * np.exp(-x / 0.3) + 60.0 / x ** 7   # noqa: E731
    lo, hi = 3.0, 6.0
    assert dv(lo) < 0 < dv(hi)
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if dv(mid) < 0 else (lo, mid)
    assert abs(d - 0.5 * (lo + hi)) < 2e-3
    assert e < 0 and abs(e - (1000.0 * np.exp(-d / 0.3) - 10.0 / d ** 6)) < 1e-9
    # rmin: local maximum below the minimum (root of V' in (0.7, 1.0)); DISCUS walks in steps of d/100
    lo, hi = 0.7, 1.0
    assert dv(lo) > 0 > dv(hi)
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if dv(mid) > 0 else (lo, mid)
    assert abs(rmin - 0.5 * (lo + hi)) < d / 100
    assert atmin > 0 and bt.minimum == (d, e, rmin, atmin)


@pytest.mark.parametrize("kw,match", [
    (dict(A=0.0), "A > 0"), (dict(rho=-1.0), "rho > 0"), (dict(B=-1.0), "B >= 0"),
    (dict(A=float("inf")), "finite"),
])
def test_buckingham_refusals(kw, match):
    args = dict(A=1000.0, rho=0.3, B=10.0)
    args.update(kw)
    with pytest.raises(ValueError, match=match):
        BuckinghamTarget(0, "Cu", "Au", **args)


def test_buckingham_without_minimum_inside_range_is_refused():
    with pytest.raises(ValueError, match="attractive"):          # V = exp(-d): no well
        BuckinghamTarget(0, "Cu", "Au", A=1.0, rho=1.0, B=0.0)
    with pytest.raises(ValueError, match="attractive"):          # tiny B: the minimum lies beyond 10 Å
        BuckinghamTarget(0, "Cu", "Au", A=1000.0, rho=0.3, B=1e-6)


def test_directional_shift_lines_and_mixed_species():
    sh = Shift({"Cu": 0.05, "Au": (0.01, 0.01, 0.02)}, 0.3, directions={"Cu": (1, 0, 0)})
    assert sh.lines() == ["set mode, 0.3, shift, all", "set move, Cu, 1.0, 0.0, 0.0, 0.05",
                          "set move, Au, 0.01, 0.01, 0.02"]
    assert Shift({"Cu": 0.05}, directions={"Cu": (0, 2, -1)}).directions == {"Cu": (0.0, 2.0, -1.0)}   # not normalised


@pytest.mark.parametrize("kw,match", [
    (dict(sigma={"Cu": (0.01, 0.01, 0.01)}, directions={"Cu": (1, 0, 0)}), "one positive number"),
    (dict(sigma={"Cu": 0.0}, directions={"Cu": (1, 0, 0)}), "must be positive"),
    (dict(sigma={"Cu": 0.05}, directions={"Cu": (0, 0, 0)}), "nonzero"),
    (dict(sigma={"Cu": (0.01, 0.01, 0.01)}, directions={"Au": (1, 0, 0)}), "no sigma"),
    (dict(sigma={"Cu": 0.05}), "three numbers"),
    (dict(sigma={"Cu": 0.05}, directions={"Cu": (float("nan"), 0, 1)}), "finite"),
    (dict(sigma={"Cu": 0.05}, directions={"Cu": (1, 0)}), r"\[u, v, w\]"),
])
def test_directional_shift_refusals(kw, match):
    with pytest.raises(ValueError, match=match):
        Shift(**kw)


def test_finish_line_sits_after_feed_and_never_writes_feed_keyword():
    cfg = MCConfig((NN,), (OccupationalTarget(0, "Cu", "Au", correlation=-0.3),), (Swap(),), SEED,
                   finish=Convergence(diff=0.01, aver=0.002))
    lines = cfg.render("in.stru", 4096).splitlines()
    i = lines.index("set feed, 4096")
    assert lines[i + 1] == "set finish, stop:converge, diff:0.01, rdiff:99999, change:99999, aver:0.002"
    assert lines[i + 2] == "set output, feed:on, final:on"
    assert not any("feed:" in x for x in lines if x.startswith("set finish"))
    assert "set finish" not in MCConfig((NN,), cfg.targets, (Swap(),), SEED).render("in.stru", 4096)


def test_convergence_refusals_and_restrictions():
    with pytest.raises(ValueError, match="at least one"):
        Convergence()
    for bad in (0.0, -1.0, float("nan")):
        with pytest.raises(ValueError):
            Convergence(diff=bad)
    cfg = dict(neighborhoods=(NN,), targets=(OccupationalTarget(0, "Cu", "Au", correlation=-0.3),),
               moves=(Swap(),), seed=SEED, finish=Convergence(diff=0.01))
    for targets in ((RepulsiveTarget(0, "Cu", "Au", 1.0),),
                    (OccupationalTarget(0, "Cu", "Au", energy=0.5),),
                    (OccupationalTarget(0, "Cu", "Au", correlation=-0.3), SpringTarget(0, "Cu", "Au", 2.5, 1.0))):
        with pytest.raises(ValueError, match="one occupational target"):
            MCConfig(**{**cfg, "targets": targets})
    MCConfig(**{**cfg, "finish": Convergence(aver=0.01)})          # one target: allowed
    with pytest.raises(ValueError, match="finish must be"):
        MCConfig(**{**cfg, "finish": {"diff": 0.01}})


def test_finish_caution_and_roundtrip():
    cfg = MCConfig((NN,), (OccupationalTarget(0, "Cu", "Au", correlation=-0.3),), (Swap(),), SEED,
                   finish=Convergence(diff=0.01))
    assert ("set finish stops at the earliest after 3 feedbacks; convergence is tested only at feedbacks"
            in cfg.cautions())
    assert MCConfig.from_dict(json.loads(json.dumps(cfg.to_dict()))) == cfg
    assert not any("set finish" in c for c in _cfg().cautions())


def test_roundtrip_new_fields_and_old_dictionaries():
    cfg = _cfg(moves=(Shift({"Cu": 0.05, "Au": (0.01, 0.01, 0.01)}, directions={"Cu": (1, 0, 0)}),),
               finish=Convergence(diff=0.01, aver=0.002))
    d = json.loads(json.dumps(cfg.to_dict()))
    assert MCConfig.from_dict(d) == cfg
    assert d["finish"] == {"diff": 0.01, "rdiff": None, "change": None, "aver": 0.002}
    old = dict(d)
    del old["finish"]
    old["moves"] = [{"type": "Shift", "sigma": {"Cu": [0.01, 0.01, 0.01]}, "probability": 1.0}]
    back = MCConfig.from_dict(old)
    assert back.finish is None and back.moves[0].directions is None


def test_roundtrip_repulsive_and_buckingham():
    cfg = MCConfig((NN,), (RepulsiveTarget(0, "Cu", "Au", 0.5, 1.2, 0.8, 6.0),
                           BuckinghamTarget(0, "Au", "Au", A=1000.0, rho=0.3, B=10.0)), (Swap(),), SEED)
    text = cfg.render("in.stru", 100)
    assert "set target, 1, repulsive, Cu, Au, 0.5, 1.2, 0.8, 6.0" in text
    assert "set target, 1, bucking, Au, Au, 0.0, 1000.0, 0.3, 10.0" in text
    assert MCConfig.from_dict(json.loads(json.dumps(cfg.to_dict()))) == cfg


def _dist_row(label, ic, a, b, nums, n):
    """Format 3900 (``'Repuls.  '``) and format 3100 (title ``a9``): label in 9 columns, then two ``a5`` atoms."""
    return (f" {ic:3d}   {label:<9}{a:>5}   {b:>5}" + " " * 8
            + "".join(f"{x:7.3f}   " for x in nums) + f"{n:8d}")


def _buck_row(ic, a, b, ach, sig, n):
    """Format 2100 (mmc_basic.f90:1590): no label, two ``a9`` atoms, achieved and sigma, count."""
    return f" {ic:3d}   {a:>9}   {b:>9}     {ach:7.3f}   {sig:7.3f}   {n:8d}"


def test_parse_repulsive_and_buckingham_rows():
    nums = (100.0, 2.812, 0.095, -97.188, -0.972)    # the target column is DISCUS's gate 100, not a target
    text = "\n".join([
        " --- Initial multiple energy configuration ---",          # mmc_correlation_write: no Buckingham row
        _dist_row("Repulsive", 1, "CU", "AU", nums, 40),
        " ---------------------------------------------",
        " Gen:       1728 try:       1728 acc: (g/n/b):      610 /        5 /      521  MC moves 0.1000",
        _dist_row("Repuls.", 1, "CU", "AU", nums, 40),
        _buck_row(1, "CU", "AU", 3.010, 0.070, 20736),
        " --- Final multiple energy configuration ---",
        "",
        " Gen:       3456 try:       3456 acc: (g/n/b):     1226 /       11 /     1056  MC moves 0.2000",
        _dist_row("Repuls.", 1, "CU", "AU", nums, 40),
        _buck_row(1, "CU", "AU", 3.011, 0.069, 20736),
        "",
    ])
    L = parse_mmc_log(text)
    assert [r["kind"] for r in L["initial"]] == ["repulsive"]
    assert [r["kind"] for r in L["feedback"][0]["rows"]] == ["repulsive", "buckingham"]
    rep = match_row(L["final"]["rows"], RepulsiveTarget(0, "Cu", "Au", depth=1.0))
    assert rep["target"] is None and rep["achieved"] == 2.812 and rep["pairs"] == 40
    buck = match_row(L["final"]["rows"], BuckinghamTarget(0, "Cu", "Au", A=1000.0, rho=0.3, B=10.0))
    assert buck["achieved"] == 3.011 and buck["sigma"] == 0.069 and buck["pairs"] == 20736


def test_shift_refuses_case_duplicate_species():
    with pytest.raises(ValueError, match="twice"):
        Shift({"Cu": 0.05, "CU": (0.01, 0.01, 0.01)}, directions={"Cu": (1, 0, 0)})
