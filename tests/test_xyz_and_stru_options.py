"""Extended XYZ writer and the .stru options for molecule configurations (check_order, min_separation)."""
import numpy as np
import pytest

from discuspdf import DiscusFormatError, Lattice, Structure, parse_stru, read_stru, write_stru
from discuspdf.structure_ops import repeat
from discuspdf.xyz_io import format_extxyz, write_extxyz


def _crystal(n=(2, 1, 1)):
    lat = Lattice.from_parameters(4.0, 5.0, 6.0, 90, 100, 90)
    cell = Structure(lat, ["C", "O"], [[0.0, 0.0, 0.0], [0.25, 0.0, 0.0]], model="configuration",
                     site_label=["a", "b"])
    return repeat(cell, n).structure


def test_extxyz_format():
    st = _crystal()
    text = format_extxyz(st)
    lines = text.splitlines()
    assert lines[0] == "4"
    lat = " ".join(f"{x:.6f}" for x in st.lattice.matrix.ravel())
    assert lines[1] == f'Lattice="{lat}" Properties=species:S:1:pos:R:3 pbc="T T T"'
    x, y, z = st.cartesian[1]
    assert lines[3] == f"O {x:.6f} {y:.6f} {z:.6f}"
    assert text.endswith("\n") and len(lines) == 6


def test_extxyz_mask_and_write(tmp_path):
    st = _crystal()
    text = format_extxyz(st, mask=[True, False, True, False])
    assert text.splitlines()[0] == "2" and all(l.startswith("C ") for l in text.splitlines()[2:])
    p = write_extxyz(tmp_path / "a.xyz", st)
    assert p.read_text() == format_extxyz(st)
    with pytest.raises(FileExistsError):
        write_extxyz(tmp_path / "a.xyz", st)
    with pytest.raises(ValueError, match="mask"):
        format_extxyz(st, mask=[True])


def test_extxyz_reads_back_with_ase(tmp_path):
    ase_io = pytest.importorskip("ase.io")
    st = _crystal()
    write_extxyz(tmp_path / "a.xyz", st)
    atoms = ase_io.read(tmp_path / "a.xyz")
    assert np.allclose(atoms.cell[:], st.lattice.matrix, atol=1e-6)
    assert np.allclose(atoms.positions, st.cartesian, atol=1e-6)


def _molecule_like():
    """Three cells; the second slot's atom lies 0.7 cell away from its position in the other cell (as atoms of a
    molecule do in another orientation), and a VOID sits exactly on a real atom (3 cells, 6 rows)."""
    st = _crystal((3, 1, 1))
    frac = st.frac.copy()
    frac[3, 0] += 0.7 / 3          # supercell fractions: 0.7 of a unit cell along a
    species = np.array(st.species, dtype=object)
    species[2] = "VOID"
    frac[2] = frac[1]              # VOID on top of the O of cell 0
    return Structure(st.lattice, species, frac, model="configuration", site_label=st.site_label,
                     cell_index=st.cell_index, supercell=st.supercell, b_iso=np.ones(6), min_separation=None)


def test_write_stru_refuses_molecule_layout_by_default(tmp_path):
    with pytest.raises((DiscusFormatError, ValueError)):
        write_stru(tmp_path / "m.stru", _molecule_like())
    assert not (tmp_path / "m.stru").exists()      # nothing left behind


def test_write_stru_molecule_options_round_trip(tmp_path):
    st = _molecule_like()
    p = write_stru(tmp_path / "m.stru", st, check_order=False, min_separation=None)
    back = read_stru(p, check_order=False, min_separation=None).structure
    assert list(back.species) == list(st.species)
    assert np.allclose(back.frac, st.frac, atol=1e-6)
    assert np.array_equal(back.cell_index, st.cell_index)
    with pytest.raises((DiscusFormatError, ValueError)):
        read_stru(p)                                 # default checks still apply when reading


def test_parse_stru_options_only_relax_what_they_name(tmp_path):
    st = _molecule_like()
    text = write_stru(tmp_path / "m.stru", st, check_order=False, min_separation=None).read_text()
    with pytest.raises(ValueError):                  # order relaxed, overlap check still on
        parse_stru(text, check_order=False)
    with pytest.raises(DiscusFormatError):           # overlap relaxed, order check still on
        parse_stru(text, min_separation=None)
