"""Tests for structure_ops.sample_configuration (Package A, Feature A4).

All use the average mixed fixture (12 rows, 8 unique sites: 4 As
single-member sites with occupancy 1.0, 4 Ga/In sites with two members each:
Ga at 0.7 and In at 0.3, total 1.0).
"""
import math
import subprocess
import sys
import warnings

import numpy as np
import pytest

from discuspdf import Lattice, Structure, repeat, substitute, displace, sample_configuration
from discuspdf.cif_io import read_cif
from discuspdf.discus_io import format_stru, parse_stru


# ------------------------------ helpers ------------------------------
def _mixed_cell() -> Structure:
    """Ga0.7In0.3As average cell from the CIF fixture (12 rows, 8 sites)."""
    cif = read_cif("tests/fixtures/interop/ga07in03as_mixed.cif")
    s = cif.structure
    assert s.model == "average"
    assert s.n_atoms == 12
    assert len(s.sites()) == 8
    return s


def _vacancy_cell() -> Structure:
    """Average structure with one site of total occupancy 0.6 (Ga 0.4 + In 0.2)
    plus one fully occupied site (As 1.0)."""
    a = 5.8
    return Structure(
        Lattice.from_parameters(a, a, a, 90.0, 90.0, 90.0),
        np.array(["Ga", "In", "As"], dtype=object),
        np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.25, 0.25, 0.25]]),
        model="average",
        occupancy=np.array([0.4, 0.2, 1.0]),
        site_label=np.array(["Ga1", "Ga1", "As1"], dtype=object),
        cell_index=np.zeros((3, 3), dtype=np.int64),
        supercell=(1, 1, 1),
    )


def _near_full_cell() -> Structure:
    """Average structure with one site of total occupancy 1 - 1e-12
    (Ga 0.5 - 5e-13 + In 0.5 - 5e-13)."""
    a = 5.8
    return Structure(
        Lattice.from_parameters(a, a, a, 90.0, 90.0, 90.0),
        np.array(["Ga", "In"], dtype=object),
        np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]),
        model="average",
        occupancy=np.array([0.5 - 5e-13, 0.5 - 5e-13]),
        site_label=np.array(["Ga1", "Ga1"], dtype=object),
        cell_index=np.zeros((2, 3), dtype=np.int64),
        supercell=(1, 1, 1),
    )


# ------------------------------ oracle ------------------------------
def test_sample_configuration_oracle():
    """Independent numpy re-implementation of the documented algorithm:
    sites in order of first appearance (min member row), one
    default_rng(seed).random(n_sites) draw, cumulative occupancies with the
    last clamped to 1.0 when the total is >= 1 - 1e-9, first member with
    u < c_j (vacancy if none).  Must match sample_configuration for seeds 0..49.
    """
    s = _mixed_cell()
    sites = sorted(s.sites(), key=lambda d: min(d["members"]))
    assert len(sites) == 8

    for seed in range(50):
        u = np.random.default_rng(seed).random(len(sites))
        oracle_rows = []
        for k, site in enumerate(sites):
            members = site["members"]
            cum = np.cumsum(s.occupancy[members])
            if cum[-1] >= 1.0 - 1e-9:
                cum = cum.copy()
                cum[-1] = 1.0
            j = int(np.searchsorted(cum, u[k], side="right"))
            if j < len(members):
                oracle_rows.append(members[j])
        result = sample_configuration(s, seed=seed)
        np.testing.assert_array_equal(np.sort(result.structure.atom_id),
                                     np.sort(s.atom_id[oracle_rows]))


# ------------------------------ structure ------------------------------
def test_sample_configuration_structure():
    """Structure of the result for a single seed."""
    s = _mixed_cell()
    result = sample_configuration(s, seed=42)
    new_s = result.structure

    # model = configuration, all occupancies 1.0.
    assert new_s.model == "configuration"
    np.testing.assert_array_equal(new_s.occupancy, np.ones(new_s.n_atoms))

    # atom_ids are a subset of the input ids, in input order.
    assert set(new_s.atom_id).issubset(set(s.atom_id))
    keep_rows = [i for i in range(s.n_atoms) if s.atom_id[i] in set(new_s.atom_id)]
    np.testing.assert_array_equal(new_s.atom_id, s.atom_id[keep_rows])

    # site_label / cell_index / frac of kept rows equal the input's.
    np.testing.assert_array_equal(new_s.site_label, s.site_label[keep_rows])
    np.testing.assert_array_equal(new_s.cell_index, s.cell_index[keep_rows])
    np.testing.assert_array_equal(new_s.frac, s.frac[keep_rows])

    # Exactly one member survives per non-vacant site.
    sites = s.sites()
    for site in sites:
        kept = [m for m in site["members"] if m in keep_rows]
        assert len(kept) <= 1

    # Fully occupied single-member sites (As) always kept.
    as_rows = [i for i in range(s.n_atoms) if s.species[i] == "As"]
    assert len(as_rows) == 4
    for i in as_rows:
        assert i in keep_rows

    # ids = removed ids (sorted).
    removed = [i for i in range(s.n_atoms) if i not in keep_rows]
    np.testing.assert_array_equal(result.ids, np.sort(s.atom_id[removed]))

    # n_vacancies == 0 when every site totals 1.0.
    assert all(site["total_occupancy"] >= 1.0 - 1e-9 for site in sites)
    assert result.info["n_vacancies"] == 0


# ------------------------------ statistics ------------------------------
def test_sample_configuration_statistics():
    """Realised Ga fraction among cation sites over 20 seeds.

    The 6x6x6 repeat gives 8 * 216 = 1728 cation sites, each independently
    realised as Ga with probability 0.7.  The per-seed Ga fraction is a
    binomial mean, so its standard error over the seeds is the binomial
    standard error sigma = sqrt(p * (1 - p) / n_cation_sites) with
    p = 0.7 and n_cation_sites = 1728; the mean over seeds must be within
    5 sigma of 0.7.
    """
    s = repeat(_mixed_cell(), (6, 6, 6)).structure
    n_cation_sites = len(s.sites())
    assert n_cation_sites == 1728

    p = 0.7
    sigma = math.sqrt(p * (1 - p) / n_cation_sites)

    fractions = []
    for seed in range(20):
        new_s = sample_configuration(s, seed=seed).structure
        n_cation = new_s.n_atoms - int(np.sum(new_s.species == "As"))
        n_ga = int(np.sum(new_s.species == "Ga"))
        fractions.append(n_ga / n_cation)
    mean_fraction = float(np.mean(fractions))
    print(f"  Ga fraction: mean={mean_fraction:.6f} (expected {p}, 5-sigma={5 * sigma:.6f})")

    assert abs(mean_fraction - p) < 5 * sigma, (
        f"mean Ga fraction {mean_fraction:.6f} not within 5 sigma of {p}"
    )

    # The sample is never identical for two different seeds.
    for i in range(20):
        for j in range(i + 1, 20):
            assert sample_configuration(s, seed=i).structure != \
                sample_configuration(s, seed=j).structure, (
                f"seeds {i} and {j} produced identical structures"
            )


# ------------------------------ vacancy ------------------------------
def test_sample_configuration_vacancy():
    """A site of total occupancy 0.6 (Ga 0.4 + In 0.2) is vacant with
    probability 0.4.  Over 1000 seeds the vacancy fraction must be within
    5 sigma of 0.4, where sigma = sqrt(0.4 * 0.6 / 1000) is the binomial
    standard error of the vacancy indicator.
    """
    s = _vacancy_cell()
    n_seeds = 1000
    p_vacant = 0.4
    sigma = math.sqrt(p_vacant * (1 - p_vacant) / n_seeds)

    n_vacant = 0
    for seed in range(n_seeds):
        result = sample_configuration(s, seed=seed)
        new_s = result.structure
        # Vacant site: no row kept at the site position.
        site_rows = [i for i in range(s.n_atoms) if s.species[i] in ("Ga", "In")]
        kept_site_rows = [i for i in site_rows if s.atom_id[i] in set(new_s.atom_id)]
        if not kept_site_rows:
            n_vacant += 1
        # info["n_vacancies"] is consistent (0 or 1 for this 2-site structure).
        assert result.info["n_vacancies"] in (0, 1)
        assert result.info["n_vacancies"] == (0 if kept_site_rows else 1)

    vacancy_fraction = n_vacant / n_seeds
    print(f"  vacancy fraction: {vacancy_fraction:.4f} (expected {p_vacant}, 5-sigma={5 * sigma:.4f})")
    assert abs(vacancy_fraction - p_vacant) < 5 * sigma, (
        f"vacancy fraction {vacancy_fraction:.4f} not within 5 sigma of {p_vacant}"
    )


# ------------------------------ tolerance ------------------------------
def test_sample_configuration_tolerance():
    """A site with total occupancy 1 - 1e-12 is treated as 1 (never vacant
    over 2000 seeds)."""
    s = _near_full_cell()
    assert abs(s.sites()[0]["total_occupancy"] - 1.0) < 1e-11

    for seed in range(2000):
        result = sample_configuration(s, seed=seed)
        assert result.info["n_vacancies"] == 0
        assert result.structure.n_atoms == 1


# ------------------------------ errors ------------------------------
def test_sample_configuration_errors():
    """Configuration input -> ValueError; all-vacant result -> ValueError."""
    s = _mixed_cell()
    # Configuration input.
    config_s = sample_configuration(s, seed=42).structure
    with pytest.raises(ValueError, match="model='average'"):
        sample_configuration(config_s, seed=42)

    # Empty result: a one-site structure with occupancy 1e-9 is always vacant
    # (the draw is never < 1e-9), so the result is empty and raises.
    a = 5.8
    tiny = Structure(
        Lattice.from_parameters(a, a, a, 90.0, 90.0, 90.0),
        np.array(["Ga"], dtype=object),
        np.array([[0.0, 0.0, 0.0]]),
        model="average",
        occupancy=np.array([1e-9]),
        site_label=np.array(["Ga1"], dtype=object),
        cell_index=np.zeros((1, 3), dtype=np.int64),
        supercell=(1, 1, 1),
    )
    with pytest.raises(ValueError, match="empty"):
        sample_configuration(tiny, seed=0)


# ------------------------------ round trip ------------------------------
def test_sample_configuration_round_trip():
    """sample_configuration(read_cif(mixed)).structure round-trips through
    format_stru / parse_stru (species, frac, occupancy) and converts to ASE
    without a partial-occupancy warning.
    """
    s = _mixed_cell()
    new_s = sample_configuration(s, seed=42).structure

    text = format_stru(new_s)
    back = parse_stru(text).structure
    assert back.n_atoms == new_s.n_atoms
    np.testing.assert_array_equal(back.species, new_s.species)
    np.testing.assert_allclose(back.frac, new_s.frac, atol=1e-9)
    np.testing.assert_array_equal(back.occupancy, new_s.occupancy)

    # A configuration structure has no partial occupancy, so to_ase must not
    # warn.
    from discuspdf import to_ase
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        to_ase(new_s)


# ------------------------------ methods ------------------------------
def test_structure_methods_equal_functions():
    """The four Structure methods equal the functions."""
    s = _mixed_cell()

    assert s.repeat((2, 1, 1)) == repeat(s, (2, 1, 1)).structure

    config_s = s.sample_configuration(seed=42)
    assert config_s.substitute(old="Ga", new="In", seed=42, count=1) == \
        substitute(config_s, old="Ga", new="In", seed=42, count=1).structure
    assert config_s.displace(seed=42, sigma=0.05) == \
        displace(config_s, seed=42, sigma=0.05).structure

    assert s.sample_configuration(seed=42) == sample_configuration(s, seed=42).structure


# ------------------------------ imports ------------------------------
def test_import_structure_does_not_import_structure_ops():
    """import discuspdf.structure alone must not import structure_ops
    (checked in a subprocess so the parent test's imports do not interfere).
    """
    code = (
        "import sys, discuspdf.structure\n"
        "print('structure_ops' in sys.modules)\n"
    )
    result = subprocess.run([sys.executable, "-c", code],
                           capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "False"
