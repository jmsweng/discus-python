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
    s = read_stru(FIX_222, supercell="infer")
    assert s.title == "GaAs"
    assert "F -4 3 m" in s.spacegroup
    assert s.ncell is None
    assert s.structure.n_atoms == 64
    assert tuple(s.structure.supercell) == (2, 2, 2)
    a, b, c, alpha, beta, gamma = s.structure.lattice.parameters
    assert a == pytest.approx(11.3074)
    assert b == pytest.approx(11.3074)
    assert c == pytest.approx(11.3074)
    assert (alpha, beta, gamma) == (90.0, 90.0, 90.0)
    from collections import Counter
    assert dict(Counter(s.structure.species)) == {"GA": 32, "AS": 32}
    assert s.structure.model == "configuration"
    assert s.structure.pbc == (True, True, True)
    # the fixture has Biso 1.000000 on every atom: read back as 1.0, no
    # default applied
    assert np.allclose(s.structure.b_iso, 1.0)
    assert "b_iso" not in s.defaults_applied


# -- 2. round-trip ----------------------------------------------------------
def test_roundtrip_gaas_222(tmp_path):
    s = read_stru(FIX_222, supercell="infer")
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
    # the fixture has Biso 0.000000 on every atom: the reader maps that to the
    # DISCUS default (1.0) and records the applied default
    assert np.allclose(s.structure.b_iso, 1.0)
    assert "b_iso" in s.defaults_applied


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
        "CU 0.500000, 0.000000, 0.000000, 0.000000, 2, 0, 0, 1.000000, _, 0, 0, 0\n"
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


# -- A3-6a: ncell-less supercell policy -------------------------------------
def test_gaas_222_is_zinc_blende():
    """gaas_222.stru read with supercell='infer' is a proper 2x2x2 zinc-blende
    crystal: 64 atoms, 32 GA + 32 AS, 4 neighbours at 2.4480 A per atom."""
    import numpy as np
    from collections import Counter
    s = read_stru(FIX_222, supercell="infer").structure
    assert tuple(s.supercell) == (2, 2, 2)
    assert s.n_atoms == 64
    assert dict(Counter(s.species)) == {"GA": 32, "AS": 32}
    # minimum distance = 5.6537 * sqrt(3) / 4 = 2.4480 A
    D = s.distance_matrix()
    np.fill_diagonal(D, 99.0)
    assert D.min() == pytest.approx(5.6537 * 3 ** 0.5 / 4, abs=1e-4)
    # every atom has exactly 4 neighbours within 2.6 A
    counts = (D < 2.6).sum(axis=1)
    assert set(counts.tolist()) == {4}
    # cell_index contains each of the 8 cells exactly 8 times
    cells = [tuple(v) for v in s.cell_index]
    from collections import Counter as C
    assert dict(C(cells)) == {
        (ix, iy, iz): 8
        for iz in range(2)
        for iy in range(2)
        for ix in range(2)
    }


def test_no_ncell_policy():
    """The ncell-less supercell policy (A3-6a decision 4)."""
    import hashlib
    import numpy as np
    # default read of gaas_222 (spans 2 cells) raises with supercell=
    with pytest.raises(DiscusFormatError) as exc:
        read_stru(FIX_222)
    assert "supercell=" in str(exc.value)
    # "infer" gives (2,2,2)
    assert tuple(read_stru(FIX_222, supercell="infer").structure.supercell) == (2, 2, 2)
    # explicit (2,2,2) agrees
    s_inf = read_stru(FIX_222, supercell="infer").structure
    s_exp = read_stru(FIX_222, supercell=(2, 2, 2)).structure
    assert np.allclose(s_inf.frac, s_exp.frac, atol=1e-9)
    # wrong tuple (4,1,1): 64 % 4 == 0 so the divisibility check passes, but
    # the cell-major order check fails (the file is a 2x2x2 layout).
    with pytest.raises((DiscusFormatError, ValueError)):
        read_stru(FIX_222, supercell=(4, 1, 1))
    # gaas_111 (0..0.75) reads with the default as (1,1,1)
    s111 = read_stru(FIX_111).structure
    assert tuple(s111.supercell) == (1, 1, 1)
    # a hand-made text with lo = -0.3, span 3, supercell="infer" raises
    # (n = round(-2 * -0.3) = 1; hi = 2.7 is not < 1/2)
    text = (
        "title t\n"
        "spcgr P1\n"
        "cell 4.0 4.0 4.0 90.0 90.0 90.0\n"
        "format numbers,XYZBPMMOS\n"
        "atoms x,y,z,Biso,Property,MoleNo,MoleAt,Occ,St,Sh,Sk,Sl\n"
        "CU -0.300000, 0.000000, 0.000000, 0.000000, 1, 0, 0, 1.000000, _, 0, 0, 0\n"
        "CU 2.700000, 0.000000, 0.000000, 0.000000, 1, 0, 0, 1.000000, _, 0, 0, 0\n"
    )
    with pytest.raises(DiscusFormatError):
        parse_stru(text, supercell="infer")
    # C1_frozen: default raises; "infer" and (8,8,8) agree; GA 1723 / AS 2373.
    with pytest.raises(DiscusFormatError):
        read_stru(FIX_C1)
    c1_inf = read_stru(FIX_C1, supercell="infer").structure
    c1_exp = read_stru(FIX_C1, supercell=(8, 8, 8)).structure
    assert tuple(c1_inf.supercell) == (8, 8, 8)
    assert np.allclose(c1_inf.frac, c1_exp.frac, atol=1e-9)
    from collections import Counter
    assert dict(Counter(c1_inf.species)) == {"GA": 1723, "AS": 2373}
    # the C1 file's bytes are unchanged (compare before and after, not with a
    # constant: a CRLF checkout changes the bytes).
    before = hashlib.sha256(open(FIX_C1, "rb").read()).hexdigest()
    read_stru(FIX_C1, supercell="infer")
    after = hashlib.sha256(open(FIX_C1, "rb").read()).hexdigest()
    assert before == after


def test_overlap_check_is_active():
    """The 0.1 A overlap check runs for .stru and .cll (min_separation is not
    disabled by the parsers)."""
    # two .stru atoms 0.05 A apart (x 0 and 0.0125 in a 4 A cell)
    text = (
        "title t\n"
        "spcgr P1\n"
        "cell 4.0 4.0 4.0 90.0 90.0 90.0\n"
        "ncell 1, 1, 1, 2\n"
        "format numbers,XYZBPMMOS\n"
        "atoms x,y,z,Biso,Property,MoleNo,MoleAt,Occ,St,Sh,Sk,Sl\n"
        "CU 0.000000, 0.000000, 0.000000, 0.000000, 1, 0, 0, 1.000000, _, 0, 0, 0\n"
        "CU 0.012500, 0.000000, 0.000000, 0.000000, 1, 0, 0, 1.000000, _, 0, 0, 0\n"
    )
    with pytest.raises((DiscusFormatError, ValueError)) as exc:
        parse_stru(text)
    assert "apart" in str(exc.value)
    # two identical listed .cll positions
    cll = (
        "title t\n"
        "spcgr P1\n"
        "cell 4.0 4.0 4.0 90.0 90.0 90.0\n"
        "atoms\n"
        "CU 0.0 0.0 0.0 1.0\n"
        "CU 0.0 0.0 0.0 1.0\n"
    )
    with pytest.raises((DiscusFormatError, ValueError)) as exc:
        from discuspdf import parse_cll
        parse_cll(cll)
    assert "apart" in str(exc.value)


def test_kind_confusion():
    """parse_cll / parse_stru reject the other kind by name."""
    import numpy as np
    stru_text = open(FIX_111, encoding="ascii").read()
    with pytest.raises(DiscusFormatError) as exc:
        from discuspdf import parse_cll
        parse_cll(stru_text)
    assert "read_stru" in str(exc.value)
    cll_text = open("tests/fixtures/discus_6_23_01/gaas.cll", encoding="ascii").read()
    with pytest.raises(DiscusFormatError) as exc:
        parse_stru(cll_text)
    assert "read_cll" in str(exc.value)
    # the p1_2x2x2_a4 .stru is rejected by parse_cll
    p1_text = open(FIX_P1, encoding="ascii").read()
    with pytest.raises(DiscusFormatError) as exc:
        from discuspdf import parse_cll
        parse_cll(p1_text)
    assert "read_stru" in str(exc.value)


def test_cartesian_positions_preserved_triclinic():
    """A triclinic 3x2x2 .stru round-trips its Cartesian positions to 1e-9."""
    import numpy as np
    a, b, c, alpha, beta, gamma = 4.1, 5.3, 6.2, 80.0, 85.0, 95.0
    unit = Lattice.from_parameters(a, b, c, alpha, beta, gamma)
    n = (3, 2, 2)
    m = 2
    N = n[0] * n[1] * n[2] * m
    slot_pos = np.array([[0.1, 0.1, 0.1], [0.6, 0.5, 0.4]])
    # build the .stru text by hand (6 decimals), cell-major, x fastest
    lines = [
        f"title t",
        f"spcgr P1",
        f"cell    {a:.6f},  {b:.6f},  {c:.6f}, {alpha:.6f}, {beta:.6f}, {gamma:.6f}",
        f"ncell       {n[0]},       {n[1]},       {n[2]},       {m}",
        "format numbers,XYZBPMMOS",
        "atoms      x,              y,              z,             Biso,    "
        "Property,  MoleNo,  MoleAt,   Occ,     St,  Sh,  Sk,  Sl",
    ]
    for r in range(N):
        cell_num = r // m
        slot = r % m
        iz = cell_num // (n[0] * n[1])
        rem = cell_num % (n[0] * n[1])
        iy = rem // n[0]
        ix = rem % n[0]
        x = slot_pos[slot][0] + ix
        y = slot_pos[slot][1] + iy
        z = slot_pos[slot][2] + iz
        lines.append(
            f"CU{' ' * max(1, 4 - 2)}     {x:>10.6f},{y:>15.6f},{z:>15.6f},"
            f"{'0.000000':>14},{'1':>8},{'0':>8},{'0':>8},"
            f"{'1.000000':>11},{'_':>2},{'0':>4},{'0':>4},{'0':>4}"
        )
    text = "\n".join(lines) + "\n"
    s = parse_stru(text).structure
    # Cartesian positions equal unit_lattice.frac_to_cart(x_file) to 1e-9
    n_arr = np.array(n, dtype=float)
    x_file = s.frac * n_arr
    cart = s.cartesian
    ref = unit.frac_to_cart(x_file)
    assert np.allclose(cart, ref, atol=1e-9)
    # write_stru -> read_stru reproduces them to 1e-9
    import tempfile, os
    with tempfile.TemporaryDirectory() as tmp:
        p = os.path.join(tmp, "rt.stru")
        write_stru(p, s)
        s2 = read_stru(p).structure
        assert np.allclose(s2.cartesian, ref, atol=1e-9)
        assert np.allclose(s2.frac, s.frac, atol=1e-9)


def test_disordered_roundtrip():
    """4x4x4 simple cubic, m=1, random species + small displacements: the
    round trip keeps species per row, cell_index and Cartesian positions."""
    import numpy as np
    rng = np.random.default_rng(7)
    lat = Lattice.from_parameters(4.0, 4.0, 4.0, 90.0, 90.0, 90.0)
    n = (4, 4, 4)
    N = 64
    species = np.array(
        ["CU" if rng.random() < 0.5 else "AG" for _ in range(N)], dtype=object
    )
    frac = np.zeros((N, 3))
    cell_index = np.zeros((N, 3), dtype=np.int64)
    for r in range(N):
        iz = r // 16
        rem = r % 16
        iy = rem // 4
        ix = rem % 4
        cell_index[r] = (ix, iy, iz)
        frac[r] = [ix / 4, iy / 4, iz / 4]
    # random displacements <= 0.05 A (in fractional units)
    disp = rng.uniform(-0.05, 0.05, size=(N, 3)) / 4.0
    frac = frac + disp
    s = Structure(
        lat, species, frac,
        model="configuration",
        b_iso=None,
        cell_index=cell_index,
        supercell=n,
        pbc=(True, True, True),
    )
    import tempfile, os
    with tempfile.TemporaryDirectory() as tmp:
        p = os.path.join(tmp, "rt.stru")
        write_stru(p, s)
        s2 = read_stru(p).structure
    assert np.array_equal(s2.species, s.species)
    assert np.array_equal(s2.cell_index, s.cell_index)
    assert np.allclose(s2.cartesian, s.cartesian, atol=1e-6)


def test_writer_header_golden():
    """The .stru header is exactly the DISCUS layout for a 2-atom structure."""
    lat = Lattice.from_parameters(4.0, 4.0, 4.0, 90.0, 90.0, 90.0)
    s = Structure(
        lat,
        np.array(["CU", "AG"], dtype=object),
        np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]),
        model="configuration",
        b_iso=np.array([1.0, 0.5]),
        cell_index=np.zeros((2, 3), dtype=np.int64),
        supercell=(1, 1, 1),
        pbc=(True, True, True),
    )
    text = format_stru(s, title="t", spacegroup="P1")
    lines = text.splitlines()
    assert lines[0] == "title t"
    assert lines[1] == "spcgr P1"
    assert lines[2] == "cell    4.000000,  4.000000,  4.000000, 90.000000, 90.000000, 90.000000"
    assert lines[3] == "ncell       1,       1,       1,       2"
    assert lines[4] == "format numbers,XYZBPMMOS"
    assert (
        lines[5]
        == "atoms      x,              y,              z,             Biso,    "
        "Property,  MoleNo,  MoleAt,   Occ,     St,  Sh,  Sk,  Sl"
    )
    assert b"\r" not in text.encode("ascii")


def test_atom_ids_not_stored():
    """atom_id is not stored in a .stru; re-reading gives 0..N-1, order kept."""
    lat = Lattice.from_parameters(4.0, 4.0, 4.0, 90.0, 90.0, 90.0)
    s = Structure(
        lat,
        np.array(["CU", "AG"], dtype=object),
        np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]),
        model="configuration",
        atom_id=np.array([10, 3], dtype=np.int64),
        cell_index=np.zeros((2, 3), dtype=np.int64),
        supercell=(1, 1, 1),
        pbc=(True, True, True),
    )
    import tempfile, os
    with tempfile.TemporaryDirectory() as tmp:
        p = os.path.join(tmp, "rt.stru")
        write_stru(p, s)
        s2 = read_stru(p).structure
    assert list(s2.atom_id) == [0, 1]
    # row order kept
    assert list(s2.species) == ["CU", "AG"]
    assert np.allclose(s2.frac, s.frac, atol=1e-9)


@pytest.mark.parametrize(
    "text, frag",
    [
        # unknown directive
        (
            "title t\nspcgr P1\ncell 4.0 4.0 4.0 90.0 90.0 90.0\nbogus 1\n"
            "format numbers,XYZBPMMOS\n"
            "atoms x,y,z,Biso,Property,MoleNo,MoleAt,Occ,St,Sh,Sk,Sl\n"
            "CU 0.000000, 0.000000, 0.000000, 0.000000, 1, 0, 0, 1.000000, _, 0, 0, 0\n",
            "unknown directive",
        ),
        # bad format line (a non-numeric token where a number is expected)
        (
            "title t\nspcgr P1\ncell 4.0 4.0 4.0 90.0 90.0 90.0\n"
            "format 12345,XYZBPMMOS\n"
            "atoms x,y,z,Biso,Property,MoleNo,MoleAt,Occ,St,Sh,Sk,Sl\n"
            "CU 0.000000, 0.000000, 0.000000, 0.000000, 1, 0, 0, 1.000000, _, 0, 0, 0\n",
            "format",
        ),
        # nan coordinate
        (
            "title t\nspcgr P1\ncell 4.0 4.0 4.0 90.0 90.0 90.0\n"
            "ncell 1, 1, 1, 1\nformat numbers,XYZBPMMOS\n"
            "atoms x,y,z,Biso,Property,MoleNo,MoleAt,Occ,St,Sh,Sk,Sl\n"
            "CU nan, 0.000000, 0.000000, 0.000000, 1, 0, 0, 1.000000, _, 0, 0, 0\n",
            "non-finite",
        ),
        # ncell with three numbers
        (
            "title t\nspcgr P1\ncell 4.0 4.0 4.0 90.0 90.0 90.0\n"
            "ncell 1, 1, 1\nformat numbers,XYZBPMMOS\n"
            "atoms x,y,z,Biso,Property,MoleNo,MoleAt,Occ,St,Sh,Sk,Sl\n"
            "CU 0.000000, 0.000000, 0.000000, 0.000000, 1, 0, 0, 1.000000, _, 0, 0, 0\n",
            "four integers",
        ),
        # ncell with a 0
        (
            "title t\nspcgr P1\ncell 4.0 4.0 4.0 90.0 90.0 90.0\n"
            "ncell 0, 1, 1, 1\nformat numbers,XYZBPMMOS\n"
            "atoms x,y,z,Biso,Property,MoleNo,MoleAt,Occ,St,Sh,Sk,Sl\n"
            "CU 0.000000, 0.000000, 0.000000, 0.000000, 1, 0, 0, 1.000000, _, 0, 0, 0\n",
            "non-positive",
        ),
        # wrong column count
        (
            "title t\nspcgr P1\ncell 4.0 4.0 4.0 90.0 90.0 90.0\n"
            "ncell 1, 1, 1, 1\nformat numbers,XYZBPMMOS\n"
            "atoms x,y,z,Biso,Property,MoleNo,MoleAt,Occ,St,Sh,Sk,Sl\n"
            "CU 0.000000, 0.000000, 0.000000, 1, 0, 0, 1.000000, _, 0, 0, 0\n",
            "12 fields",
        ),
        # duplicate atoms line
        (
            "title t\nspcgr P1\ncell 4.0 4.0 4.0 90.0 90.0 90.0\n"
            "ncell 1, 1, 1, 1\nformat numbers,XYZBPMMOS\n"
            "atoms x,y,z,Biso,Property,MoleNo,MoleAt,Occ,St,Sh,Sk,Sl\n"
            "atoms x,y,z,Biso,Property,MoleNo,MoleAt,Occ,St,Sh,Sk,Sl\n"
            "CU 0.000000, 0.000000, 0.000000, 0.000000, 1, 0, 0, 1.000000, _, 0, 0, 0\n",
            "repeated 'atoms'",
        ),
        # species of 5 characters
        (
            "title t\nspcgr P1\ncell 4.0 4.0 4.0 90.0 90.0 90.0\n"
            "ncell 1, 1, 1, 1\nformat numbers,XYZBPMMOS\n"
            "atoms x,y,z,Biso,Property,MoleNo,MoleAt,Occ,St,Sh,Sk,Sl\n"
            "ABCDE 0.000000, 0.000000, 0.000000, 0.000000, 1, 0, 0, 1.000000, _, 0, 0, 0\n",
            "4 characters",
        ),
    ],
)
def test_parse_rejects(text, frag):
    with pytest.raises(DiscusFormatError) as exc:
        parse_stru(text)
    msg = str(exc.value)
    assert frag in msg
    # the message names the line number
    assert ":" in msg


def test_order_check_runs_for_m1():
    """The order check runs for m=1 too: an ncell 2,2,2,1 file with two rows
    swapped (cells out of order) must be rejected."""
    # 8 atoms, m=1, ncell 2,2,2,1.  Build a valid cell-major file, then swap
    # two rows so the cells are out of order.
    n = (2, 2, 2)
    m = 1
    N = 8
    # valid: row r in cell r; slot positions distinct per cell
    lines = [
        "title t",
        "spcgr P1",
        "cell 4.000000,  4.000000,  4.000000, 90.000000, 90.000000, 90.000000",
        f"ncell       {n[0]},       {n[1]},       {n[2]},       {m}",
        "format numbers,XYZBPMMOS",
        "atoms      x,              y,              z,             Biso,    "
        "Property,  MoleNo,  MoleAt,   Occ,     St,  Sh,  Sk,  Sl",
    ]
    rows = []
    for r in range(N):
        iz = r // 4
        rem = r % 4
        iy = rem // 2
        ix = rem % 2
        # slot 0 position in cell (ix,iy,iz): (ix*0.1, iy*0.1, iz*0.1)
        rows.append((ix * 0.1, iy * 0.1, iz * 0.1))
    # swap rows 0 and 1 (cells (0,0,0) and (1,0,0) out of order)
    rows[0], rows[1] = rows[1], rows[0]
    for (x, y, z) in rows:
        lines.append(
            f"CU{' ' * 2}     {x:>10.6f},{y:>15.6f},{z:>15.6f},"
            f"{'0.000000':>14},{'1':>8},{'0':>8},{'0':>8},"
            f"{'1.000000':>11},{'_':>2},{'0':>4},{'0':>4},{'0':>4}"
        )
    text = "\n".join(lines) + "\n"
    with pytest.raises(DiscusFormatError, match="order"):
        parse_stru(text)


def test_large_roundtrip_speed():
    """GaAs 16x16x16 (32768 atoms, 8 per cell): write+read under 15 s."""
    import time
    import numpy as np
    lat = Lattice.from_parameters(5.6537, 5.6537, 5.6537, 90.0, 90.0, 90.0)
    n = (16, 16, 16)
    m = 8
    N = n[0] * n[1] * n[2] * m
    unit_frac = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.5, 0.5, 0.0],
            [0.0, 0.5, 0.5],
            [0.5, 0.0, 0.5],
            [0.25, 0.25, 0.25],
            [0.75, 0.75, 0.25],
            [0.25, 0.75, 0.75],
            [0.75, 0.25, 0.75],
        ]
    )
    unit_species = np.array(
        ["GA", "GA", "GA", "GA", "AS", "AS", "AS", "AS"], dtype=object
    )
    # cell-major, x fastest (vectorised)
    cell_num = np.arange(N) // m
    ix = (cell_num % n[0]).astype(np.int64)
    iy = ((cell_num // n[0]) % n[1]).astype(np.int64)
    iz = (cell_num // (n[0] * n[1])).astype(np.int64)
    cell_index = np.stack([ix, iy, iz], axis=1)
    slot = np.arange(N) % m
    species = unit_species[slot]
    frac = (unit_frac[slot] + cell_index) / np.array(n, dtype=float)
    s = Structure(
        lat, species, frac,
        model="configuration",
        b_iso=np.full(N, 0.35),
        cell_index=cell_index,
        supercell=n,
        pbc=(True, True, True),
    )
    import tempfile, os
    with tempfile.TemporaryDirectory() as tmp:
        p = os.path.join(tmp, "big.stru")
        t0 = time.perf_counter()
        write_stru(p, s)
        s2 = read_stru(p).structure
        dt = time.perf_counter() - t0
    assert dt < 15.0, f"round trip took {dt:.2f} s (limit 15 s)"
    assert np.allclose(s2.frac, s.frac, atol=2e-6)
    assert np.array_equal(s2.species, s.species)


def test_species_policy():
    """Ga, O2-, Fe3+ are written verbatim and read back verbatim; a 5-char
    name raises on write."""
    lat = Lattice.from_parameters(4.0, 4.0, 4.0, 90.0, 90.0, 90.0)
    species = np.array(["Ga", "O2-", "Fe3+"], dtype=object)
    frac = np.array([[0.0, 0.0, 0.0], [0.25, 0.25, 0.25], [0.5, 0.5, 0.5]])
    s = Structure(
        lat, species, frac,
        model="configuration",
        b_iso=np.ones(3),
        cell_index=np.zeros((3, 3), dtype=np.int64),
        supercell=(1, 1, 1),
        pbc=(True, True, True),
    )
    import tempfile, os
    with tempfile.TemporaryDirectory() as tmp:
        p = os.path.join(tmp, "sp.stru")
        write_stru(p, s)
        # written verbatim
        text = open(p, encoding="ascii").read()
        assert "Ga " in text or "Ga," in text
        assert "O2-" in text
        assert "Fe3+" in text
        s2 = read_stru(p).structure
        assert list(s2.species) == ["Ga", "O2-", "Fe3+"]
    # 5-char name raises on write
    bad = Structure(
        lat,
        np.array(["ABCDE"], dtype=object),
        np.array([[0.0, 0.0, 0.0]]),
        model="configuration",
        b_iso=np.ones(1),
        cell_index=np.zeros((1, 3), dtype=np.int64),
        supercell=(1, 1, 1),
        pbc=(True, True, True),
    )
    with pytest.raises(ValueError):
        format_stru(bad)


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


# -- A3-9: B-factor policy ---------------------------------------------------
def _stru_two_atom(b_iso) -> Structure:
    """Two-atom P1 structure (model="configuration") for the B-policy tests."""
    lat = Lattice.from_parameters(4.0, 4.0, 4.0, 90.0, 90.0, 90.0)
    return Structure(
        lat,
        np.array(["Cu", "Ag"], dtype=object),
        np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]),
        model="configuration",
        b_iso=b_iso,
        cell_index=np.zeros((2, 3), dtype=np.int64),
        supercell=(1, 1, 1),
        pbc=(True, True, True),
    )


# -- 1. b_iso None -> explicit 1.0 (DISCUS default) --------------------------
def test_b_iso_none_becomes_explicit_default(tmp_path):
    s = _stru_two_atom(None)
    text = format_stru(s)
    # the written text contains an explicit Biso column (1.000000, not 0.000000)
    assert "1.000000" in text
    # the Biso value is the 4th field on each atom line
    atom_lines = [l for l in text.splitlines() if l.split(",")[0].split()[0] in ("Cu", "Ag")]
    assert len(atom_lines) == 2
    assert all(l.split(",")[3].strip() == "1.000000" for l in atom_lines)
    out = tmp_path / "b1.stru"
    write_stru(out, s)
    s2 = read_stru(out)
    assert np.allclose(s2.structure.b_iso, [1.0, 1.0])
    # B was written explicitly, so no default is applied on read-back
    assert "b_iso" not in s2.defaults_applied


# -- 2. writer rejects B below the DISCUS representable floor ---------------
@pytest.mark.parametrize("bad_b", [0.0, 1e-6, 5e-4])
def test_writer_rejects_unrepresentable_b(tmp_path, bad_b):
    s = _stru_two_atom(np.array([bad_b, 0.35]))
    with pytest.raises(DiscusFormatError) as exc:
        format_stru(s)
    msg = str(exc.value)
    assert "DISCUS" in msg
    assert "0.001" in msg
    assert "2" in msg  # atom count
    # write_stru raises too, and leaves neither the target nor a temp file
    out = tmp_path / "bad.stru"
    with pytest.raises(DiscusFormatError):
        write_stru(out, s)
    assert not out.exists()
    leftovers = [p.name for p in tmp_path.iterdir() if p.name.endswith(".tmp")]
    assert leftovers == []


# -- 2b. the writer's B check rejects a negative B directly ------------------
def test_validate_b_for_write_rejects_negative():
    from discuspdf.discus_io import _validate_b_for_write
    with pytest.raises(DiscusFormatError) as exc:
        _validate_b_for_write(np.array([-0.1, 1.0]), 2)
    msg = str(exc.value)
    assert "DISCUS" in msg
    assert "0.001" in msg


# -- 2c. Structure itself still rejects a negative B -------------------------
def test_structure_still_rejects_negative_b():
    lat = Lattice.from_parameters(4.0, 4.0, 4.0, 90.0, 90.0, 90.0)
    with pytest.raises(ValueError, match="b_iso"):
        Structure(
            lat,
            np.array(["Cu"], dtype=object),
            np.array([[0.0, 0.0, 0.0]]),
            model="configuration",
            b_iso=np.array([-1.0]),
            min_separation=None,
        )


# -- 3. small-but-safe B round-trips -----------------------------------------
@pytest.mark.parametrize("good_b", [0.001, 0.35, 2.0])
def test_small_but_safe_b_roundtrips(tmp_path, good_b):
    s = _stru_two_atom(np.array([good_b, 1.0]))
    out = tmp_path / "b_safe.stru"
    write_stru(out, s)
    s2 = read_stru(out)
    assert np.allclose(s2.structure.b_iso, [good_b, 1.0], atol=1e-6)
    assert "b_iso" not in s2.defaults_applied


# -- 4. reader maps a 0 B to the DISCUS default ------------------------------
def test_reader_maps_zero_b_to_discus_default():
    # the tracked fixture has Biso 0.000000 on all 8 atoms
    s = read_stru(FIX_P1)
    assert np.allclose(s.structure.b_iso, 1.0)
    assert "b_iso" in s.defaults_applied
    # a hand-made .stru with explicit 0.0 Biso behaves the same
    text = (
        "title t\n"
        "spcgr P1\n"
        "cell 4.0 4.0 4.0 90.0 90.0 90.0\n"
        "ncell 1, 1, 1, 2\n"
        "format numbers,XYZBPMMOS\n"
        "atoms x,y,z,Biso,Property,MoleNo,MoleAt,Occ,St,Sh,Sk,Sl\n"
        "CU 0.000000, 0.000000, 0.000000, 0.000000, 1, 0, 0, 1.000000, _, 0, 0, 0\n"
        "CU 0.500000, 0.000000, 0.000000, 0.000000, 1, 0, 0, 1.000000, _, 0, 0, 0\n"
    )
    s2 = parse_stru(text)
    assert np.allclose(s2.structure.b_iso, [1.0, 1.0])
    assert "b_iso" in s2.defaults_applied


# -- 5. reader rejects a tiny (0, 0.001) B -----------------------------------
def test_reader_rejects_tiny_b():
    text = (
        "title t\n"
        "spcgr P1\n"
        "cell 4.0 4.0 4.0 90.0 90.0 90.0\n"
        "ncell 1, 1, 1, 2\n"
        "format numbers,XYZBPMMOS\n"
        "atoms x,y,z,Biso,Property,MoleNo,MoleAt,Occ,St,Sh,Sk,Sl\n"
        "CU 0.000000, 0.000000, 0.000000, 0.000100, 1, 0, 0, 1.000000, _, 0, 0, 0\n"
        "CU 0.500000, 0.000000, 0.000000, 1.000000, 1, 0, 0, 1.000000, _, 0, 0, 0\n"
    )
    with pytest.raises(DiscusFormatError) as exc:
        parse_stru(text)
    msg = str(exc.value)
    assert "0.001" in msg
    # the message names a line number
    assert ":" in msg
