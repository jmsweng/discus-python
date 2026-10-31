"""Tests for the DISCUS ``.cll`` reader/writer (Package A, Feature A3).

Uses only the fixtures directory and ``tmp_path``; never writes into the
repository tree.
"""
from __future__ import annotations

import numpy as np
import pytest

from discuspdf import (
    DiscusFormatError,
    DiscusCell,
    Lattice,
    Structure,
    read_cll,
    parse_cll,
    format_cll,
    write_cll,
)

FIX = "tests/fixtures/discus_6_23_01/gaas.cll"
FIX_STRU = "tests/fixtures/discus_6_23_01/gaas_111.stru"


def _gaas() -> DiscusCell:
    return read_cll(FIX)


# -- 1. fixture fields -------------------------------------------------------
def test_read_fixture_gaas():
    c = _gaas()
    assert c.title == "GaAs"
    assert c.spacegroup == "F-43m"
    assert c.structure.n_atoms == 2
    a, b, cc, alpha, beta, gamma = c.structure.lattice.parameters
    assert a == pytest.approx(5.6537)
    assert b == pytest.approx(5.6537)
    assert cc == pytest.approx(5.6537)
    assert (alpha, beta, gamma) == (90.0, 90.0, 90.0)
    assert list(c.structure.species) == ["GA", "AS"]
    assert np.allclose(c.structure.frac, [[0.0, 0.0, 0.0], [0.25, 0.25, 0.25]])
    assert np.allclose(c.structure.b_iso, [1.0, 1.0])
    assert c.structure.model == "average"
    assert tuple(c.structure.supercell) == (1, 1, 1)
    assert c.structure.pbc == (True, True, True)
    assert c.defaults_applied == ()


# -- 2. round-trip ----------------------------------------------------------
def test_roundtrip_fixture(tmp_path):
    c = _gaas()
    out = tmp_path / "rt.cll"
    write_cll(out, c)
    c2 = read_cll(out)
    assert c.structure.isclose(c2.structure, rtol=1e-8, atol=1e-6)
    assert c.title == c2.title
    assert c.spacegroup == c2.spacegroup


# -- 3. Structure -> .cll -> Structure --------------------------------------
def test_structure_to_cll_to_structure():
    lat = Lattice.from_parameters(5.0, 6.0, 7.0, 80.0, 85.0, 95.0)
    species = np.array(["Ga", "Ga", "As", "As", "O"], dtype=object)
    frac = np.array(
        [
            [0.1, 0.2, 0.3],
            [0.4, 0.5, 0.6],
            [0.7, 0.8, 0.9],
            [0.05, 0.15, 0.25],
            [0.33, 0.44, 0.55],
        ]
    )
    s = Structure(lat, species, frac, model="average", min_separation=None)
    text = format_cll(s)
    c = parse_cll(text)
    a1, b1, c1, al1, be1, g1 = c.structure.lattice.parameters
    a0, b0, c0, al0, be0, g0 = s.lattice.parameters
    assert a1 == pytest.approx(a0, rel=1e-9)
    assert b1 == pytest.approx(b0, rel=1e-9)
    assert c1 == pytest.approx(c0, rel=1e-9)
    assert al1 == pytest.approx(al0, rel=1e-9)
    assert be1 == pytest.approx(be0, rel=1e-9)
    assert g1 == pytest.approx(g0, rel=1e-9)
    assert np.allclose(c.structure.frac, s.frac, atol=0.5e-6)
    assert list(c.structure.species) == list(species)
    # order identical
    assert np.array_equal(c.structure.frac, s.frac, equal_nan=True) or np.allclose(
        c.structure.frac, s.frac, atol=0.5e-6
    )


# -- 4. b_iso None -> explicit 0.0 ------------------------------------------
def test_b_iso_none_becomes_explicit_zero(tmp_path):
    lat = Lattice.from_parameters(4.0, 4.0, 4.0, 90.0, 90.0, 90.0)
    s = Structure(
        lat,
        np.array(["Cu", "Ag"], dtype=object),
        np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]),
        model="average",
        b_iso=None,
        min_separation=None,
    )
    text = format_cll(s)
    # the written text contains an explicit B column (0.000000)
    assert "0.000000" in text
    out = tmp_path / "b0.cll"
    write_cll(out, s)
    c = read_cll(out)
    assert np.allclose(c.structure.b_iso, [0.0, 0.0])


# -- 5. missing B column -> DISCUS default ----------------------------------
def test_missing_b_column_uses_discus_default():
    text = (
        "title t\n"
        "spcgr P1\n"
        "cell 4.0 4.0 4.0 90.0 90.0 90.0\n"
        "atoms\n"
        "Cu 0.0 0.0 0.0\n"
        "Ag 0.5 0.5 0.5\n"
    )
    c = parse_cll(text)
    assert np.allclose(c.structure.b_iso, [1.0, 1.0])
    assert c.defaults_applied == ("b_iso",)


# -- 6. rejection rules -----------------------------------------------------
@pytest.mark.parametrize(
    "text, frag",
    [
        (
            "title t\nspcgr P1\ncell 4 4 4 90 90 90\nbogus 1 2 3\natoms\nCu 0 0 0 1\n",
            "unknown directive",
        ),
        (
            "title t\nspcgr P1\natoms\nCu 0 0 0 1\n",
            "missing required 'cell'",
        ),
        (
            "title t\nspcgr P1\ncell 4 4 4 90 90\natoms\nCu 0 0 0 1\n",
            "exactly six numbers",
        ),
        (
            "title t\nspcgr P1\ncell 4 4 4 90 90 90\natoms\n",
            "no atom lines",
        ),
        (
            "title t\nspcgr P1\ncell 4 4 4 90 90 90\natoms\nCu 0 0\n",
            "name + x,y,z",
        ),
        (
            "title t\nspcgr P1\ncell 4 4 4 90 90 90\natoms\nCu nan 0 0 1\n",
            "non-finite",
        ),
        (
            "title t\nspcgr P1\ncell 4 4 4 90 90 90\natoms\nCu 0 0 0 -1\n",
            "negative B",
        ),
        (
            # 5 numeric values (name + x,y,z,B,occ) is a .stru-only marker.
            "title t\nspcgr P1\ncell 4 4 4 90 90 90\natoms\nCu 0 0 0 1 1.5\n",
            "extended .stru columns",
        ),
        (
            "title t\ntitle u\nspcgr P1\ncell 4 4 4 90 90 90\natoms\nCu 0 0 0 1\n",
            "repeated 'title'",
        ),
    ],
)
def test_rejects(text, frag):
    with pytest.raises(DiscusFormatError) as exc:
        parse_cll(text)
    msg = str(exc.value)
    assert frag in msg
    # message names the line number and the offending text
    assert ":" in msg  # source:line prefix


# -- 7. reject a saved structure -------------------------------------------
def test_rejects_saved_structure():
    text = open(FIX_STRU, encoding="ascii").read()
    with pytest.raises(DiscusFormatError) as exc:
        parse_cll(text)
    assert "read_stru" in str(exc.value)


# -- 8. write semantics -----------------------------------------------------
def test_write_refuses_overwrite(tmp_path):
    out = tmp_path / "x.cll"
    write_cll(out, _gaas())
    with pytest.raises(FileExistsError):
        write_cll(out, _gaas())


def test_overwrite_replaces_atomically(tmp_path):
    out = tmp_path / "x.cll"
    write_cll(out, _gaas())
    write_cll(out, _gaas(), overwrite=True)
    assert read_cll(out).title == "GaAs"


def test_failed_write_leaves_target_and_no_temp(tmp_path, monkeypatch):
    out = tmp_path / "x.cll"
    real_replace = __import__("os").replace

    def boom(*a, **k):
        raise OSError("simulated replace failure")

    monkeypatch.setattr(__import__("os"), "replace", boom)
    with pytest.raises(OSError):
        write_cll(out, _gaas())
    assert not out.exists()
    # no stray temp file in the folder
    leftovers = [p.name for p in tmp_path.iterdir() if p.name.endswith(".tmp")]
    assert leftovers == []
    monkeypatch.setattr(__import__("os"), "replace", real_replace)


# -- 9. no CR in written bytes ---------------------------------------------
def test_written_bytes_have_no_cr(tmp_path):
    out = tmp_path / "x.cll"
    write_cll(out, _gaas())
    assert b"\r" not in out.read_bytes()


# -- 10. reject supercell and mixed pbc ------------------------------------
def test_write_rejects_supercell_and_mixed_pbc():
    lat = Lattice.from_parameters(4.0, 4.0, 4.0, 90.0, 90.0, 90.0)
    s = Structure(
        lat,
        np.array(["Cu"], dtype=object),
        np.array([[0.0, 0.0, 0.0]]),
        model="configuration",
        supercell=(2, 2, 2),
        min_separation=None,
    )
    with pytest.raises(ValueError, match="unit cell"):
        format_cll(s)
    s2 = Structure(
        lat,
        np.array(["Cu"], dtype=object),
        np.array([[0.0, 0.0, 0.0]]),
        model="configuration",
        pbc=(True, False, True),
        min_separation=None,
    )
    with pytest.raises(ValueError, match="pbc"):
        format_cll(s2)


# -- 11. CRLF and comment handling ------------------------------------------
def test_crlf_and_comment_handling():
    text = (
        "title t\r\n"
        "spcgr P1\r\n"
        "cell 4.0 4.0 4.0 90.0 90.0 90.0\r\n"
        "# a comment\r\n"
        "atoms\r\n"
        "Cu 0.0 0.0 0.0 1.0\r\n"
        "! another comment\r\n"
        "Ag 0.5 0.5 0.5 1.0\r\n"
    )
    c = parse_cll(text)
    assert c.structure.n_atoms == 2
    assert list(c.structure.species) == ["Cu", "Ag"]


# -- 12. species name policy ------------------------------------------------
def test_species_names_policy():
    # 4-char names with charges round-trip verbatim
    text = (
        "title t\n"
        "spcgr P1\n"
        "cell 4.0 4.0 4.0 90.0 90.0 90.0\n"
        "atoms\n"
        "O2- 0.0 0.0 0.0 1.0\n"
        "Fe3+ 0.5 0.5 0.5 1.0\n"
    )
    c = parse_cll(text)
    assert list(c.structure.species) == ["O2-", "Fe3+"]
    # over-length (5 chars) rejected
    bad = (
        "title t\n"
        "spcgr P1\n"
        "cell 4.0 4.0 4.0 90.0 90.0 90.0\n"
        "atoms\n"
        "ABCDE 0.0 0.0 0.0 1.0\n"
    )
    with pytest.raises(DiscusFormatError, match="4 characters"):
        parse_cll(bad)
