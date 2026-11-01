"""Tests for structure_ops results through the DISCUS .stru / .cll writers.

Covers the format_stru / parse_stru round trip for structures produced by
``repeat``, ``substitute`` and ``displace`` (Package A, Feature A4), plus the
unit-cell ``format_cll`` path and an opt-in live DISCUS round trip.

The live test is opt-in: set ``DISCUS_LIVE=1`` to enable it.  Evidence is
written under ``tests/fixtures/structure_ops/live_evidence/<timestamp>/`` when
``DISCUS_LIVE_EVIDENCE`` points at that directory.
"""
from __future__ import annotations

import hashlib
import os
import re
import time
from collections import Counter
from pathlib import Path

import numpy as np
import pytest

from discuspdf import (
    DiscusFormatError,
    displace,
    format_cll,
    format_stru,
    parse_stru,
    repeat,
    sample_configuration,
    substitute,
    write_stru,
)
from discuspdf.cif_io import read_cif

from tests.discus_live import evidence_log, evidence_note, run_discus

FIX_GAAS = "tests/fixtures/interop/gaas_cod_9008845.cif"

# Skip the live test unless DISCUS_LIVE=1.
skip_live = pytest.mark.skipif(
    os.environ.get("DISCUS_LIVE") != "1",
    reason="live DISCUS test; set DISCUS_LIVE=1 to enable",
)


# -- helpers ---------------------------------------------------------------
def _gaas_cell() -> object:
    """GaAs 8-atom cell from the CIF fixture (average model)."""
    cif = read_cif(FIX_GAAS)
    return cif.structure


def _gaas_configuration() -> object:
    """GaAs 8-atom cell as a configuration (fixed seed)."""
    s = _gaas_cell()
    return sample_configuration(s, seed=42).structure


def _gaas_repeated() -> object:
    """GaAs 8-atom configuration cell repeated 4x4x4 (512 atoms)."""
    return repeat(_gaas_configuration(), (4, 4, 4)).structure


def _periodic_key(frac: np.ndarray) -> tuple:
    """Periodic key of a fractional position: (frac mod 1) rounded to 1e-6.

    DISCUS saves supercell coordinates centred (−n/2 …, see
    docs/DISCUS_FORMATS.md), so positions are compared modulo 1 instead of
    row by row.
    """
    return tuple(round(v % 1.0, 6) for v in frac)


def _periodic_positions(structure: object, species: str) -> list[tuple]:
    """Sorted list of periodic positions of ``species`` in ``structure``."""
    mask = structure.species == species
    keys = [_periodic_key(structure.frac[i]) for i in np.where(mask)[0]]
    return sorted(keys)


# -- 1. repeat round trip ---------------------------------------------------
def test_format_stru_roundtrip_repeated():
    """repeat 4x4x4 -> format_stru -> parse_stru: species, frac, occupancy and
    the b_iso None -> 1.0 policy all survive."""
    s = _gaas_repeated()
    text = format_stru(s)
    back = parse_stru(text).structure

    assert back.n_atoms == s.n_atoms
    assert np.array_equal(back.species, s.species)
    np.testing.assert_allclose(back.frac, s.frac, atol=1e-6)
    np.testing.assert_array_equal(back.occupancy, s.occupancy)
    # b_iso None -> explicit 1.0 on write, 1.0 on read-back (no default applied).
    assert s.b_iso is None
    assert np.allclose(back.b_iso, 1.0)


# -- 2. substitute round trip -------------------------------------------------
def test_format_stru_roundtrip_substituted():
    """repeat 4x4x4, substitute 25% Ga -> In (seed 42) -> format_stru ->
    parse_stru: species counts and frac survive."""
    s = _gaas_repeated()
    sub = substitute(s, old="Ga", new="In", seed=42, fraction=0.25).structure

    text = format_stru(sub)
    back = parse_stru(text).structure

    assert back.n_atoms == sub.n_atoms
    assert dict(Counter(back.species)) == dict(Counter(sub.species))
    np.testing.assert_allclose(back.frac, sub.frac, atol=1e-6)


# -- 3/4. displace round trip -------------------------------------------------
@pytest.mark.parametrize("sigma", [0.05, 0.4])
def test_format_stru_roundtrip_displaced(tmp_path, sigma):
    """repeat 4x4x4, displace with sigma (seed 42) -> format_stru -> parse_stru.

    The writer (format_stru) checks row order from cell_index only and never
    rejects a displacement (probe P5b); up to sigma = 0.4 A the reader
    (parse_stru) accepts the rows because every atom stays within 0.5 unit cell
    of its slot's median on each axis.
    """
    s = _gaas_repeated()
    disp = displace(s, seed=42, sigma=sigma).structure

    text = format_stru(disp)
    back = parse_stru(text).structure
    assert np.array_equal(back.species, disp.species)
    assert np.array_equal(back.cell_index, disp.cell_index)
    np.testing.assert_allclose(back.frac, disp.frac, atol=1e-6)


def test_parse_stru_rejects_large_displacement():
    """repeat 4x4x4, displace with sigma = 1.0 A (seed 42): format_stru
    succeeds (the writer checks cell_index only, probe P5b) but parse_stru
    raises DiscusFormatError because atoms leave the 0.5-unit-cell median
    window."""
    s = _gaas_repeated()
    disp = displace(s, seed=42, sigma=1.0).structure

    text = format_stru(disp)  # the writer never rejects displacements

    with pytest.raises(DiscusFormatError, match="atom order does not match"):
        parse_stru(text)


# -- 5. unit-cell .cll --------------------------------------------------------
def test_format_cll_unit_cell_parent():
    """format_cll works for the 1x1x1 GaAs parent and contains the lattice
    parameters."""
    s = _gaas_configuration()
    assert tuple(s.supercell) == (1, 1, 1)

    text = format_cll(s)

    a, b, c, alpha, beta, gamma = s.lattice.parameters
    assert f"{a:.6f}" in text
    assert f"{b:.6f}" in text
    assert f"{c:.6f}" in text
    assert f"{alpha:.6f}" in text
    assert f"{beta:.6f}" in text
    assert f"{gamma:.6f}" in text


# -- 6. live DISCUS round trip (opt-in) ----------------------------------------
@skip_live
@pytest.mark.discus_live
def test_live_discus_stru_roundtrip(tmp_path):
    """DISCUS ``read stru`` -> ``eval n[1]`` -> ``save back_<name>`` on the
    Python-written 4x4x4 GaAs and the substituted (25% Ga -> In) structure.

    The round trip is verified against DISCUS's own ``save`` output (not our
    written file): the saved .stru is read back with ``parse_stru`` and its
    atom count, species counts and (periodic) coordinates are compared with
    the structure we wrote.

    Evidence is saved under
    ``tests/fixtures/structure_ops/live_evidence/<timestamp>/`` via
    ``evidence_log``/``evidence_note`` (requires ``DISCUS_LIVE_EVIDENCE`` to
    point there); a ``MANIFEST.sha256`` of the written .stru files is written
    into the evidence folder.
    """
    rep = _gaas_repeated()
    sub = substitute(rep, old="Ga", new="In", seed=42, fraction=0.25).structure

    # Write both .stru files into a temp workdir.
    workdir = str(tmp_path / "live")
    os.makedirs(workdir, exist_ok=True)
    rep_name = "gaas_rep4.stru"
    sub_name = "gaas_sub4.stru"
    write_stru(os.path.join(workdir, rep_name), rep)
    write_stru(os.path.join(workdir, sub_name), sub)

    # Evidence folder: tests/fixtures/structure_ops/live_evidence/<timestamp>/
    # (yyyyMMddTHHmmZ, e.g. 20261101T0111Z — no seconds).
    ts = time.strftime("%Y%m%dT%H%MZ", time.gmtime())
    evidence_dir = os.path.join(
        "tests", "fixtures", "structure_ops", "live_evidence", ts
    )
    os.environ["DISCUS_LIVE_EVIDENCE"] = evidence_dir

    for name, s in (
        (rep_name, rep),
        (sub_name, sub),
    ):
        macro = (
            "discus\n"
            "read\n"
            f"  stru {name}\n"
            "eval n[1]\n"
            f"save back_{name}\n"
            "exit\n"
            "exit\n"
        )
        result = run_discus(macro, workdir)

        # DISCUS atom count (eval n[1] prints the total number of atoms).
        values = re.findall(r"Value of n\[1\] =\s+([\d.]+)", result.stdout)
        assert values, f"no n[1] value in DISCUS output: {result.stdout[:500]!r}"
        n_discus = int(float(values[0]))
        assert n_discus == s.n_atoms, (
            f"DISCUS atom count {n_discus} != Python {s.n_atoms}"
        )

        # Read DISCUS's own saved file.  DISCUS `save` drops the ncell line
        # and writes centred coordinates (−n/2 …, docs/DISCUS_FORMATS.md),
        # so pass the supercell explicitly.
        back_path = os.path.join(workdir, f"back_{name}")
        assert os.path.exists(back_path), (
            f"DISCUS did not write {back_path}: {result.stdout[:500]!r}"
        )
        back_text = Path(back_path).read_text(encoding="ascii")
        back = parse_stru(back_text, supercell=(4, 4, 4)).structure

        # Atom count of DISCUS's saved file.
        assert back.n_atoms == s.n_atoms, (
            f"saved-file atom count {back.n_atoms} != Python {s.n_atoms}"
        )

        # Species counts: DISCUS upper-cases species names on save.
        back_counts = {k.upper(): v for k, v in back.counts().items()}
        py_counts = {k.upper(): v for k, v in s.counts().items()}
        assert back_counts == py_counts, (
            f"DISCUS species counts {back_counts} != Python {py_counts}"
        )

        # Coordinates: compare sorted periodic positions per species
        # ((frac mod 1) rounded to 1e-6) — the saved file is centred, so
        # positions are compared modulo 1 instead of row by row.
        # DISCUS upper-cases species names on save, so iterate the original
        # (mixed-case) names from the Python structure and look them up in
        # the saved structure by their upper-cased form.
        for sp in s.counts():
            ours = _periodic_positions(s, sp)
            theirs = _periodic_positions(back, sp.upper())
            assert len(ours) == len(theirs), (
                f"species {sp}: {len(ours)} positions vs {len(theirs)}"
            )
            for o, t in zip(ours, theirs):
                np.testing.assert_allclose(o, t, atol=1e-6)

        # Evidence: macro, input file hash, DISCUS output (head of the
        # saved file and the n[1] value).
        evidence_log(
            f"test_live_discus_stru_roundtrip_{name}",
            macro,
            {name: Path(workdir, name).read_bytes()},
            result,
        )
        evidence_note(
            f"test_live_discus_stru_roundtrip_{name}",
            f"--- DISCUS saved {f'back_{name}'} ---\n"
            f"n[1] = {n_discus}\n"
            f"{back_text[:2000]}",
        )

    # MANIFEST.sha256: sha256 of each written .stru file in the evidence folder.
    manifest_lines = []
    for name in (rep_name, sub_name):
        data = Path(workdir, name).read_bytes()
        sha = hashlib.sha256(data).hexdigest()
        manifest_lines.append(f"{sha}  {name}")
    manifest_path = Path(evidence_dir) / "MANIFEST.sha256"
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text("\n".join(manifest_lines) + "\n", encoding="ascii")
