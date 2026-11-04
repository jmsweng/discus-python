"""Unit tests for the shared validators (discuspdf._validate) and the shared DISCUS atom-name check."""
from __future__ import annotations

import pytest

from discuspdf import _validate as v
from discuspdf.atom_names import check_discus_name


def test_numbers():
    assert v.real("x", 2) == 2.0 and v.positive("x", 0, allow_zero=True) == 0.0 and v.integer("n", 3) == 3
    assert v.probability("p", 1) == 1.0 and v.vec3("v", (1, 2, 3)) == (1.0, 2.0, 3.0) and v.int3("t", [1, -1, 0]) == (1, -1, 0)
    for f, args in ((v.real, (True,)), (v.real, (float("nan"),)), (v.positive, (0,)), (v.integer, (0,)),
                    (v.integer, (1.0,)), (v.probability, (1.5,)), (v.vec3, ((1, 2),)), (v.vec3, (5,)),
                    (v.int3, ((1, 2, 3.5),))):
        with pytest.raises(ValueError):
            f("x", *args)


def test_plain_file_name():
    assert v.plain_file_name("f", "a_b-1.nipl") == "a_b-1.nipl"
    for bad in ("a b", "a,b", "../x", ".", "", 3):
        with pytest.raises(ValueError):
            v.plain_file_name("f", bad)


def test_check_discus_name():
    assert check_discus_name("a", " Cu ") == "Cu" and check_discus_name("a", "O2-") == "O2-"
    assert check_discus_name("a", 3, allow_number=True) == 3 and check_discus_name("a", "all", allow_all=True) == "all"
    for bad, kw in (("all", {}), ("COPPER", {}), ("C u", {}), ("Cu,", {}), (3, {}), (-1, {"allow_number": True}),
                    ("", {}), (True, {"allow_number": True})):
        with pytest.raises(ValueError):
            check_discus_name("a", bad, **kw)
