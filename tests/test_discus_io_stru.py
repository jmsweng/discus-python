"""Tests for the DISCUS ``.stru`` reader/writer (Package A, Feature A3).

Uses only the fixtures directory and ``tmp_path``; never writes into the
repository tree.
"""
from __future__ import annotations

import numpy as np
import pytest

from discuspdf import (
    DiscusFormatError,
    DiscusStru,
    Lattice,
    Structure,
    read_stru,
    parse_stru,
    format_stru,
    write_stru,
)

FIX_222 = "tests/fixtures/discus_6_23_01/gaas_222.stru"
FIX_111 = "tests/fixtures/discus_6_23_01/gaas_111.stru"
FIX_P1 = "tests/fixtures/discus_6_23_01/p1_2x2x2_a4.stru"
FIX_C1 = "benchmarks/discus_6_23_01/gate1_c1_dyadic/inputs/C1_frozen.stru"


# -- 1. fixture fields -------------------------------------------------------
def test_read_fixture_gaas_222():
    s = read_stru(FIX_222)
    assert s.title == "GaAs"
    assert "F -4 3 m" in s.spacegroup
    assert s.ncell is None
    assert s.structure.n_atoms == 64
    assert tuple(s.structure.supercell) == (1, 1, 1)
    a, b, c, alpha, beta, gamma = s.structure.lattice.parameters
    assert a == pytest.approx(5.6537)
    assert b == pytest.approx(5.6537)
    assert c == pytest.approx(5.6537)
    assert (alpha, beta, gamma) == (90.0, 90.0, 90.0)
    from collections import Counter
    assert dict(Counter(s.structure.species)) == {"GA": 32, "AS": 32}
    assert s.structure.model == "configuration"
    assert s.structure.pbc == (True, True, True)


# -- 2. round-trip ----------------------------------------------------------
def test_roundtrip_gaas_222(tmp_path):
    s = read_stru(FIX_222)
    out = tmp_path / "rt.stru"
    write_stru(out, s)
    s2 = read_stru(out)
    assert s.structure.isclose(s2.structure, rtol=1e-8, atol=1e-6)
    assert s.title == s2.title


# -- 3. ncell fixture -------------------------------------------------------
def test_read_fixture_p1_2x2x2():
    s = read_stru(FIX_P1)
    assert s.title == "p1_2x2x2_a4"
    assert s.ncell == (2, 2, 2, 1)
    assert s.structure.n_atoms == 8
    assert tuple(s.structure.supercell) == (2, 2, 2)
    a, b, c, alpha, beta, gamma = s.structure.lattice.parameters
    assert a == pytest.approx(8.0)
    assert b == pytest.approx(8.0)
    assert c == pytest.approx(8.0)
    from collections import Counter
    assert dict(Counter(s.structure.species)) == {"CU": 4, "AG": 4}


# -- 4. C1 no-ncell with supercell ------------------------------------------
def test_read_c1_no_ncell_with_supercell():
    s = read_stru(FIX_C1, supercell=(8, 8, 8))
    assert s.structure.n_atoms == 4096
    from collections import Counter
    assert dict(Counter(s.structure.species)) == {"GA": 1723, "AS": 2373}


# -- 5. ncell/supercell conflict -------------------------------------------
def test_ncell_supercell_conflict():
    text = open(FIX_P1, encoding="ascii").read()
    with pytest.raises(DiscusFormatError, match="disagree"):
        parse_stru(text, supercell=(3, 3, 3))


# -- 6. ncell count mismatch ------------------------------------------------
def test_ncell_count_mismatch():
    text = (
        "title t\n"
        "spcgr P1\n"
        "cell 4.0 4.0 4.0 90.0 90.0 90.0\n"
        "ncell 2, 2, 2, 4\n"
        "format numbers,XYZBPMMOS\n"
        "atoms x,y,z,Biso,Property,MoleNo,MoleAt,Occ,St,Sh,Sk,Sl\n"
        "CU 0.000000, 0.000000, 0.000000, 0.000000, 1, 0, 0, 1.000000, _, 0, 0, 0\n"
        "CU 1.000000, 0.000000, 0.000000, 0.000000, 1, 0, 0, 1.000000, _, 0, 0, 0\n"
    )
    with pytest.raises(DiscusFormatError, match="atom count"):
        parse_stru(text)


# -- 7. occupancy != 1 ------------------------------------------------------
def test_occupancy_not_one():
    # In the saved .stru format the Occ column is the 8th value (index 8 in
    # tokens, after the species name).  Set Occ = 0.5 on one atom.
    text = (
        "title t\n"
        "spcgr P1\n"
        "cell 4.0 4.0 4.0 90.0 90.0 90.0\n"
        "ncell 1, 1, 1, 2\n"
        "format numbers,XYZBPMMOS\n"
        "atoms x,y,z,Biso,Property,MoleNo,MoleAt,Occ,St,Sh,Sk,Sl\n"
        "CU 0.000000, 0.000000, 0.000000, 0.000000, 1, 0, 0, 0.500000, _, 0, 0, 0\n"
        "CU 1.000000, 0.000000, 0.000000, 0.000000, 1, 0, 0, 1.000000, _, 0, 0, 0\n"
    )
    with pytest.raises(DiscusFormatError, match="occupancy"):
        parse_stru(text)


# -- 8. C1 order check ------------------------------------------------------
def test_c1_order_check():
    s = read_stru(FIX_C1, supercell=(8, 8, 8))
    assert s.structure.n_atoms == 4096


# -- 9. allow_loss ----------------------------------------------------------
def test_allow_loss_drops_nonconstant_column():
    text = (
        "title t\n"
        "spcgr P1\n"
        "cell 4.0 4.0 4.0 90.0 90.0 90.0\n"
        "ncell 1, 1, 1, 2\n"
        "format numbers,XYZBPMMOS\n"
        "atoms x,y,z,Biso,Property,MoleNo,MoleAt,Occ,St,Sh,Sk,Sl\n"
        "CU 0.000000, 0.000000, 0.000000, 0.000000, 1, 0, 0, 1.000000, _, 0, 0, 0\n"
        "CU 1.000000, 0.000000, 0.000000, 0.000000, 2, 0, 0, 1.000000, _, 0, 0, 0\n"
    )
    with pytest.raises(DiscusFormatError, match="Property"):
        parse_stru(text)
    s = parse_stru(text, allow_loss=True)
    assert "Property" in s.dropped


# -- 10. reject .cll --------------------------------------------------------
def test_rejects_cll():
    text = open("tests/fixtures/discus_6_23_01/gaas.cll", encoding="ascii").read()
    with pytest.raises(DiscusFormatError, match="read_cll"):
        parse_stru(text)


# -- 11. write semantics ----------------------------------------------------
def test_write_refuses_overwrite(tmp_path):
    s = read_stru(FIX_P1)
    out = tmp_path / "x.stru"
    write_stru(out, s)
    with pytest.raises(FileExistsError):
        write_stru(out, s)


def test_overwrite_replaces_atomically(tmp_path):
    s = read_stru(FIX_P1)
    out = tmp_path / "x.stru"
    write_stru(out, s)
    write_stru(out, s, overwrite=True)
    assert read_stru(out).title == "p1_2x2x2_a4"


def test_failed_write_leaves_target_and_no_temp(tmp_path, monkeypatch):
    s = read_stru(FIX_P1)
    out = tmp_path / "x.stru"
    real_replace = __import__("os").replace

    def boom(*a, **k):
        raise OSError("simulated replace failure")

    monkeypatch.setattr(__import__("os"), "replace", boom)
    with pytest.raises(OSError):
        write_stru(out, s)
    assert not out.exists()
    leftovers = [p.name for p in tmp_path.iterdir() if p.name.endswith(".tmp")]
    assert leftovers == []
    monkeypatch.setattr(__import__("os"), "replace", real_replace)


# -- 12. reject average model ----------------------------------------------
def test_write_rejects_average_model():
    lat = Lattice.from_parameters(4.0, 4.0, 4.0, 90.0, 90.0, 90.0)
    s = Structure(
        lat,
        np.array(["Cu"], dtype=object),
        np.array([[0.0, 0.0, 0.0]]),
        model="average",
        min_separation=None,
    )
    with pytest.raises(ValueError, match="realised crystal"):
        format_stru(s)


# -- 13. reject out-of-order rows ------------------------------------------
def test_write_rejects_out_of_order_rows():
    lat = Lattice.from_parameters(8.0, 8.0, 8.0, 90.0, 90.0, 90.0)
    species = np.array(["CU", "AG", "CU", "AG", "CU", "AG", "CU", "AG"], dtype=object)
    frac = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.5, 0.0, 0.0],
            [0.0, 0.5, 0.0],
            [0.5, 0.5, 0.0],
            [0.0, 0.0, 0.5],
            [0.5, 0.0, 0.5],
            [0.0, 0.5, 0.5],
            [0.5, 0.5, 0.5],
        ]
    )
    # cell_index in DISCUS order: x fastest
    cell_index = np.array(
        [
            [0, 0, 0],
            [1, 0, 0],
            [0, 1, 0],
            [1, 1, 0],
            [0, 0, 1],
            [1, 0, 1],
            [0, 1, 1],
            [1, 1, 1],
        ],
        dtype=np.int64,
    )
    # Swap cell_index rows 0 and 1 to break the DISCUS cell-major order.
    cell_index_bad = cell_index.copy()
    cell_index_bad[0], cell_index_bad[1] = cell_index_bad[1], cell_index_bad[0]
    s_bad = Structure(
        lat,
        species,
        frac,
        model="configuration",
        cell_index=cell_index_bad,
        supercell=(2, 2, 2),
        min_separation=None,
    )
    with pytest.raises(DiscusFormatError, match="not in DISCUS order"):
        format_stru(s_bad)
