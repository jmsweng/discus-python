"""DISCUS atom-name table and name guessing (discuspdf.atom_names)."""
import pytest

from discuspdf.atom_names import (DISCUS_NAMES, AtomName, atomic_number, discus_atom_name, discus_name_rows,
                                  element_of, is_isotope_label, mass_number)
from discuspdf._discus_names import DISCUS_NAME_TABLE


def test_table_shape_and_known_values():
    assert len(DISCUS_NAME_TABLE) == 225 and len(DISCUS_NAMES) == 225
    assert DISCUS_NAME_TABLE[0] == ("VOID", "VOID", 0.0)
    # values used by earlier engine probes (D reference probe, B_DISCUS) and the builder
    assert DISCUS_NAMES["CU"] == ("CU", 7.718) and DISCUS_NAMES["AU"] == ("AU", 7.63)
    assert DISCUS_NAMES["GA"] == ("GA", 7.288) and DISCUS_NAMES["AS"] == ("AS", 6.58)
    assert DISCUS_NAMES["NI"] == ("NI", 10.3) and DISCUS_NAMES["H"] == ("H", -3.739)
    assert DISCUS_NAMES["D"] == ("D", 6.671) and DISCUS_NAMES["D1-"] == ("D", 6.671)
    assert DISCUS_NAMES["AL3+"] == ("AL", 4.1491)   # as in the source (Si's value), not "fixed"
    assert all(n == n.upper() and len(n) <= 4 for n, _, _ in DISCUS_NAME_TABLE)


def test_scattering_rows_drop_void_names():
    rows = discus_name_rows()
    assert len(rows) == 215
    assert not {"VOID", "POSI", "XAXI", "CENT", "E1-"} & {r[0] for r in rows}
    assert len(discus_name_rows(scattering_only=False)) == 225


@pytest.mark.parametrize("name,expected", [
    ("Ni", AtomName("NI", True, "NI", "Ni", 10.3)),
    ("Ni62", AtomName("NI62", False, "NI", "Ni", 10.3)),
    ("C1", AtomName("C1", False, "C", "C", 6.646)),
    ("FE3+", AtomName("FE3+", True, "FE3+", "Fe", DISCUS_NAMES["FE3+"][1])),
    ("O2-", AtomName("O2-", True, "O2-", "O", 5.803)),
    ("D", AtomName("D", True, "D", "D", 6.671)),
    ("VOID", AtomName("VOID", True, "VOID", "Void", 0.0)),
    ("OH", AtomName("OH", False, "O", "O", 5.803)),            # first two letters fail, first letter wins
    ("FM", AtomName("FM", False, "F", "F", 5.654)),            # DISCUS has no Fm entry: guessed as F
])
def test_discus_atom_name(name, expected):
    assert discus_atom_name(name) == expected


def test_sign_truncation_before_two_letters():
    # "NA1+" is a table name; "NA3+" is not, and is read via "NA" (no "NA3+" entry), like DISCUS
    assert discus_atom_name("NA3+").read_as == "NA"
    assert discus_atom_name("O1-").known and discus_atom_name("O3-").read_as == "O"


@pytest.mark.parametrize("name", ["", None, "XX", "Q", "E", "ES", "  "])
def test_unreadable(name):
    assert discus_atom_name(name) is None and element_of(name) is None
    with pytest.raises(ValueError):
        atomic_number(name)


@pytest.mark.parametrize("name,el,z", [("D", "H", 1), ("Ni62", "Ni", 28), ("cu", "Cu", 29), ("SI.", "Si", 14)])
def test_element_and_z(name, el, z):
    assert element_of(name) == el and atomic_number(name) == z


def test_void_is_not_an_element():
    assert element_of("VOID") is None and element_of("void") is None
    assert element_of("CENT") is None   # special name, not cerium


@pytest.mark.parametrize("name,a", [("D", 2), ("D1-", 2), ("Ni62", 62), ("O18", 18), ("Li7", 7), ("U238", 238),
                                    ("C1", None), ("Fe2", None), ("H2", None), ("Ni", None), ("FE3+", None),
                                    ("C13", 13), ("C40", None)])
def test_mass_number(name, a):
    assert mass_number(name) == a


def test_isotope_label():
    assert is_isotope_label("Ni62") and is_isotope_label("C13")
    assert not is_isotope_label("D") and not is_isotope_label("Ni") and not is_isotope_label("C1")


def test_long_names_truncated_and_spaces_refused():
    assert discus_atom_name("Ni62x").name == "NI62"      # DISCUS keeps 4 characters
    assert discus_atom_name("N i") is None
