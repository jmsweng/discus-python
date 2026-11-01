"""MCConfig rendering, validation and persistence (Package D, Feature 15). No engine needed."""
from __future__ import annotations

import numpy as np
import pytest

from discuspdf.fourier_macro import Seed
from discuspdf.lattice import Lattice
from discuspdf.mc import (DisplacementTarget, LennardJonesTarget, MCConfig, Neighborhood, OccupationalTarget,
                          Shift, SpringTarget, Swap, SwapDisplacement, Vector)
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
                    DisplacementTarget(0, "Cu", "Cu", (1, 0, 0), correlation=0.3),
                    SpringTarget(1, "Cu", "Au", 2.6, 1.0),
                    LennardJonesTarget(1, "Au", "Au", 2.9, 1.0)),
                   (Swap(0.5, "site"), Shift({"Cu": (0.01, 0.01, 0.02)}, 0.3), SwapDisplacement(0.2)), SEED)
    lines = cfg.render("in.stru", 100).splitlines()
    for expected in ["set target, 1, corr, Cu, Au, 0.0, 0.5, ENER",
                     "set target, 1, cd, Cu, Cu, 0.3, 0.0, CORR, dir:[1.0,0.0,0.0]",
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
    (dict(targets=(OccupationalTarget(0, "Cu", "Au", correlation=-0.3),
                   OccupationalTarget(0, "Cu", "Ni", correlation=0.2))), "at most one"),
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
