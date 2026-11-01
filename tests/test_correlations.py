"""Tests for discuspdf.correlations (Package C, Features 13 and 14).

Oracles:
* ordered alloys with textbook Warren–Cowley values (B2 on bcc, L1_2 Cu3Au on fcc);
* random occupation: values within their null-model standard errors, and the
  spread over many independent configurations equals the predicted standard error;
* an explicit Python loop over the pairs for counts, probabilities and covariances;
* constructed displacement waves u(R) = A cos(2π q·R), whose correlation for the
  neighbor vector t is cos(2π q·t) (and sin(2π q·t) for the cross term).
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from discuspdf import Lattice, Structure, neighbor_list, sample_configuration
from discuspdf.correlations import (
    concentrations,
    displacement_correlation,
    displacement_correlation_table,
    displacements,
    node_species,
    occupational_correlation,
    occupational_correlation_table,
)


# ------------------------------ helpers ------------------------------
def _cell(kind: str, species, a: float = 4.0) -> Structure:
    base = {
        "sc": [[0, 0, 0]],
        "bcc": [[0, 0, 0], [0.5, 0.5, 0.5]],
        "fcc": [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]],
    }[kind]
    return Structure(Lattice.from_parameters(a, a, a, 90, 90, 90), np.array(species, object), base,
                     model="configuration", site_label=[f"S{k + 1}" for k in range(len(base))])


def _with_species(structure: Structure, species) -> Structure:
    return structure.replace(species=np.asarray(species, dtype=object))


def _random_binary(n=10, c=0.3, seed=0, kind="sc"):
    ref = _cell(kind, ["A"] * {"sc": 1, "bcc": 2, "fcc": 4}[kind]).repeat((n, n, n))
    rng = np.random.default_rng(seed)
    spec = np.where(rng.random(ref.n_atoms) < c, "B", "A")
    return ref, _with_species(ref, spec)


def _loop_counts(nl, spec):
    counts = {}
    for i, j in zip(nl.i, nl.j):
        key = (spec[i], spec[j])
        counts[key] = counts.get(key, 0) + 1
    return counts


# ------------------------------ occupational: ordered alloys ------------------------------
def test_b2_order_on_bcc():
    s = _cell("bcc", ["Cs", "Cl"]).repeat((4, 4, 4))
    nl = neighbor_list(s, 4.1)
    first = occupational_correlation(nl.select(shell=1))
    assert first.value("warren_cowley", "Cs", "Cl") == pytest.approx(-1.0)
    assert first.value("warren_cowley", "Cs", "Cs") == pytest.approx(-1.0)
    assert first.value("pearson", "Cs", "Cs") == pytest.approx(-1.0)
    assert first.discus("Cs", "Cl") == pytest.approx(-1.0)
    second = occupational_correlation(nl.select(shell=2))
    assert second.value("warren_cowley", "Cs", "Cl") == pytest.approx(1.0)
    assert second.value("warren_cowley", "Cs", "Cs") == pytest.approx(1.0)
    assert second.discus("Cs", "Cl") == pytest.approx(1.0)
    assert first.n_pairs == 128 * 8 and first.n_unique_pairs == 128 * 4 and first.n_reversed == first.n_pairs


def test_l12_cu3au():
    s = _cell("fcc", ["Au", "Cu", "Cu", "Cu"]).repeat((3, 3, 3))
    nl = neighbor_list(s, 4.1)
    first = occupational_correlation(nl.select(shell=1))
    assert first.value("warren_cowley", "Au", "Cu") == pytest.approx(-1.0 / 3.0)
    assert first.value("warren_cowley", "Cu", "Au") == pytest.approx(-1.0 / 3.0)
    assert first.discus("Au", "Cu") == pytest.approx(-1.0 / 3.0)
    assert first.value("pearson", "Au", "Au") == pytest.approx(-1.0 / 3.0)
    np.testing.assert_allclose(first.origin_concentration, [0.25, 0.75])
    second = occupational_correlation(nl.select(shell=2))
    assert second.value("warren_cowley", "Au", "Au") == pytest.approx(1.0)
    assert second.value("conditional", "Au", "Au") == pytest.approx(1.0)


def test_sublattice_resolved_pure_origin_is_undefined_for_pearson_only():
    s = _cell("bcc", ["Cs", "Cl"]).repeat((3, 3, 3))
    nl = neighbor_list(s, 3.5).select(origin_site="S1")       # all origins are Cs
    res = occupational_correlation(nl)
    assert math.isnan(res.value("pearson", "Cs", "Cl"))
    assert "pure end member" in res.undefined[("Cs", "Cl")] or "absent" in res.undefined[("Cs", "Cl")]
    # Warren–Cowley is defined: P(Cl|Cs) = 1, c_t(Cl) = 1  ->  (1 - 1)/(0 - 1) = 0
    assert res.value("warren_cowley", "Cs", "Cl") == pytest.approx(0.0)
    assert math.isnan(res.value("warren_cowley", "Cl", "Cs"))
    assert "absent at the origins" in res.undefined[("Cl", "Cs")]


def test_pure_end_member_and_no_pairs():
    s = _cell("sc", ["A"]).repeat((3, 3, 3))
    res = occupational_correlation(neighbor_list(s, 4.1))
    assert res.species == ("A",)
    assert math.isnan(res.value("pearson", "A", "A")) and math.isnan(res.value("warren_cowley", "A", "A"))
    assert "pure end member" in res.undefined[("A", "A")]
    empty = occupational_correlation(neighbor_list(s, 4.1).select(distance=100.0))
    assert empty.n_pairs == 0 and empty.undefined[("A", "A")] == "no pairs"
    _, cfg = _random_binary(4, 0.5, 1)
    res = occupational_correlation(neighbor_list(cfg, 4.1), species=("A", "C"))
    assert math.isnan(res.discus("A", "C"))


# ------------------------------ occupational: definitions against a loop ------------------------------
def test_counts_and_formulas_match_an_explicit_loop_ternary():
    ref = _cell("fcc", ["A"] * 4).repeat((3, 3, 3))
    rng = np.random.default_rng(4)
    spec = rng.choice(np.array(["A", "B", "C"]), size=ref.n_atoms, p=[0.5, 0.3, 0.2])
    cfg = _with_species(ref, spec)
    nl = neighbor_list(ref, 2.9).select(origin_site="S1")       # asymmetric selection
    res = occupational_correlation(nl, cfg)
    loop = _loop_counts(nl, spec[nl.node_row])
    n = nl.n_pairs
    for a in "ABC":
        for b in "ABC":
            nab = loop.get((a, b), 0)
            assert res.value("counts", a, b) == nab
            co = sum(loop.get((a, x), 0) for x in "ABC") / n
            ct = sum(loop.get((x, b), 0) for x in "ABC") / n
            p_cond = nab / (co * n)
            assert res.value("conditional", a, b) == pytest.approx(p_cond)
            assert res.value("warren_cowley", a, b) == pytest.approx((p_cond - ct) / ((a == b) - ct))
            assert res.value("pearson", a, b) == pytest.approx(
                (nab / n - co * ct) / math.sqrt(co * (1 - co) * ct * (1 - ct)))
    # identities: rows of P(b|a) sum to 1; Σ_b α_ab (δ_ab − c_b) = 0
    np.testing.assert_allclose(res.conditional.sum(axis=1), 1.0)
    delta = np.eye(3) - res.target_concentration[None, :]
    np.testing.assert_allclose((res.warren_cowley * delta).sum(axis=1), 0.0, atol=1e-12)


def test_discus_convention_on_an_asymmetric_list_matches_its_formula():
    """θ pools both ends: (2 N_aa + N_ab + N_ba) / 2N; P_aa = N_aa / N, pairs restricted to {a, b}."""
    ref = _cell("bcc", ["A", "A"]).repeat((5, 5, 5))
    rng = np.random.default_rng(21)
    on_s1 = ref.site_label == "S1"
    spec = np.where(rng.random(ref.n_atoms) < np.where(on_s1, 0.7, 0.2), "B", "A")
    spec[rng.random(ref.n_atoms) < 0.1] = "C"
    cfg = _with_species(ref, spec)
    nl = neighbor_list(ref, 3.5).select(origin_site="S1", target_site="S2")
    loop = _loop_counts(nl, spec[nl.node_row])
    naa, nab, nba, nbb = (loop.get(k, 0) for k in [("B", "B"), ("B", "A"), ("A", "B"), ("A", "A")])
    n = naa + nab + nba + nbb
    theta = (2 * naa + nab + nba) / (2 * n)
    expected = (naa / n - theta ** 2) / (theta * (1 - theta))
    res = occupational_correlation(nl, cfg)
    assert res.discus("B", "A") == pytest.approx(expected)
    assert nab != nba                                   # really asymmetric


def test_discus_convention_equals_pearson_for_symmetric_lists():
    _, cfg = _random_binary(6, 0.35, 2)
    nl = neighbor_list(cfg, 4.1)
    sel = nl.select(origin_site="S1", target_site="S1", translation=(1, 0, 0), symmetric=True)
    res = occupational_correlation(sel)
    assert res.n_reversed == res.n_pairs
    assert res.discus("B", "A") == pytest.approx(res.value("pearson", "B", "B"))
    assert res.discus("A", "B") == pytest.approx(-res.value("pearson", "A", "B"))
    assert res.value("warren_cowley", "A", "B") == pytest.approx(res.value("pearson", "A", "A"))


# ------------------------------ occupational: random disorder ------------------------------
def test_random_disorder_is_zero_within_standard_error():
    _, cfg = _random_binary(20, 0.3, 11)
    nl = neighbor_list(cfg, 7.0)
    for shell in (1, 2, 3):
        res = occupational_correlation(nl.select(shell=shell))
        z = res.value("pearson", "A", "A") / res.value("stderr_pearson", "A", "A")
        assert abs(z) < 4, (shell, z)
        zab = res.value("warren_cowley", "A", "B") / res.value("stderr_warren_cowley", "A", "B")
        assert abs(zab) < 4


@pytest.mark.parametrize("symmetric", [True, False])
def test_null_standard_error_matches_spread_over_configurations(symmetric):
    """The predicted null SE equals the observed standard deviation over 300 random crystals (±15 %)."""
    ref = _cell("sc", ["A"]).repeat((8, 8, 8))
    nl = neighbor_list(ref, 4.1).select(translation=(1, 0, 0), symmetric=symmetric)
    rng = np.random.default_rng(123)
    vals, ses_ab, vals_ab = [], [], []
    for _ in range(300):
        spec = rng.choice(np.array(["A", "B", "C"]), size=ref.n_atoms, p=[0.6, 0.3, 0.1])
        res = occupational_correlation(nl, _with_species(ref, spec))
        vals.append(res.value("pearson", "A", "A"))
        vals_ab.append(res.value("pearson", "A", "B"))
        ses_ab.append(res.value("stderr_pearson", "A", "B"))
    se_aa = 1 / math.sqrt(res.n_unique_pairs)
    assert res.value("stderr_pearson", "A", "A") == pytest.approx(se_aa, rel=0.05)
    assert np.std(vals) == pytest.approx(se_aa, rel=0.15)
    assert np.std(vals_ab) == pytest.approx(np.mean(ses_ab), rel=0.15)
    assert abs(np.mean(vals)) < 4 * se_aa / math.sqrt(300) + 3 / ref.n_atoms


# ------------------------------ vacancies, concentrations, tables ------------------------------
def test_vacancies_from_sample_configuration():
    lat = Lattice.from_parameters(3, 3, 3, 90, 90, 90)
    avg = Structure(lat, np.array(["Ni", "Fe"], object), [[0, 0, 0], [0, 0, 0]], model="average",
                    occupancy=[0.6, 0.3], site_label=["M", "M"]).repeat((8, 8, 8))
    cfg = sample_configuration(avg, seed=5).structure
    nl = neighbor_list(avg, 3.1)
    spec = node_species(nl, cfg)
    n_va = int(np.sum(spec == "Va"))
    assert n_va == avg.n_atoms // 2 - cfg.n_atoms and n_va > 0
    res = occupational_correlation(nl, cfg)
    assert res.species == ("Fe", "Ni", "Va")
    assert res.counts.sum() == nl.n_pairs
    conc = concentrations(nl, cfg)
    assert conc["global"]["Va"] == pytest.approx(n_va / nl.n_nodes)
    assert conc["by_site"]["M"] == conc["global"]
    with pytest.raises(ValueError, match="vacancy=None"):
        occupational_correlation(nl, cfg, vacancy=None)
    with pytest.raises(ValueError, match="also a species"):
        occupational_correlation(nl, cfg, vacancy="Ni")
    with pytest.raises(ValueError, match="average-model"):
        occupational_correlation(nl)


def test_site_resolved_concentrations():
    s = _cell("fcc", ["Au", "Cu", "Cu", "Cu"]).repeat((2, 2, 2))
    conc = concentrations(neighbor_list(s, 3.0))
    assert conc["global"] == {"Au": 0.25, "Cu": 0.75}
    assert conc["by_site"]["S1"] == {"Au": 1.0, "Cu": 0.0}
    assert conc["n_by_site"]["S2"] == 8


def test_occupational_table_by_shell_and_vector():
    s = _cell("bcc", ["Cs", "Cl"]).repeat((3, 3, 3))
    nl = neighbor_list(s, 4.1)
    rows = occupational_correlation_table(nl, by="shell")
    first = [r for r in rows if r["shell"] == 1 and r["a"] == "Cs" and r["b"] == "Cl"][0]
    assert first["warren_cowley"] == pytest.approx(-1.0) and first["distance"] == pytest.approx(2 * math.sqrt(3))
    vec = occupational_correlation_table(nl, by="vector")
    groups = {(r["origin_site"], r["target_site"], r["translation"]) for r in vec}
    assert ("S1", "S2", (0, 0, 0)) in groups and ("S2", "S1", (1, 1, 1)) in groups
    assert all(r["n_pairs"] == 27 for r in vec)
    with pytest.raises(ValueError, match="by must be"):
        occupational_correlation_table(nl, by="angle")


def test_configuration_must_match_the_neighbor_list():
    ref = _cell("sc", ["A"]).repeat((3, 3, 3))
    nl = neighbor_list(ref, 3.1)
    with pytest.raises(ValueError, match="supercell"):
        occupational_correlation(nl, _cell("sc", ["A"]).repeat((3, 3, 2)))
    other = Structure(Lattice.from_parameters(12.1, 12, 12, 90, 90, 90), ref.species, ref.frac,
                      model="configuration", site_label=ref.site_label, cell_index=ref.cell_index,
                      supercell=(3, 3, 3))
    with pytest.raises(ValueError, match="lattices"):
        occupational_correlation(nl, other)


# ------------------------------ displacements ------------------------------
def _wave_configuration(reps=(16, 4, 4), a=3.0, amp=0.05, h=2, transverse=False, lattice=None):
    lat = lattice or Lattice.from_parameters(a, a, a, 90, 90, 90)
    ref = Structure(lat, np.array(["A"], object), [[0, 0, 0]], model="configuration",
                    site_label=["S1"]).repeat(reps)
    phase = 2 * np.pi * h * ref.cell_index[:, 0] / reps[0]
    u = np.zeros((ref.n_atoms, 3))
    u[:, 0] = amp * np.cos(phase)
    if transverse:
        u[:, 1] = amp * np.sin(phase)
    cfg = ref.with_cartesian(ref.cartesian + u)
    return ref, cfg, u


@pytest.mark.parametrize("k", [1, 2, 3, 4])
def test_longitudinal_wave_gives_cosine_correlation(k):
    ref, cfg, _ = _wave_configuration()
    nl = neighbor_list(ref, 3.0 * 4 + 0.1)
    sel = nl.select(translation=(k, 0, 0), symmetric=True)
    for centered in (True, False):
        res = displacement_correlation(sel, cfg, direction=(1, 0, 0), centered=centered)
        assert res.correlation == pytest.approx(math.cos(2 * math.pi * 2 * k / 16), abs=1e-12)
    # projections perpendicular to the wave have zero variance
    perp = displacement_correlation(sel, cfg, direction=(0, 1, 0))
    assert math.isnan(perp.correlation) and "zero variance" in perp.undefined


@pytest.mark.parametrize("k", [1, 3])
def test_cross_correlation_of_a_circular_wave_gives_sine(k):
    ref, cfg, _ = _wave_configuration(transverse=True)
    nl = neighbor_list(ref, 3.0 * 3 + 0.1)
    forward = nl.select(translation=(k, 0, 0))                    # one direction only
    res = displacement_correlation(forward, cfg, direction=(1, 0, 0), target_direction=(0, 1, 0))
    assert res.correlation == pytest.approx(math.sin(2 * math.pi * 2 * k / 16), abs=1e-12)
    # both directions together cancel the odd (sine) part
    both = displacement_correlation(nl.select(translation=(k, 0, 0), symmetric=True), cfg,
                                    direction=(1, 0, 0), target_direction=(0, 1, 0))
    assert both.correlation == pytest.approx(0.0, abs=1e-12)
    np.testing.assert_allclose(res.correlation_matrix[0, 1], res.correlation, atol=1e-12)


def test_wave_along_a_fractional_direction_in_a_triclinic_cell():
    lat = Lattice.from_parameters(3.0, 3.5, 4.0, 80, 95, 110)
    ref, cfg_unused, _ = _wave_configuration(reps=(12, 3, 3), lattice=lat)
    # displacement along the a axis (fractional [100]) with a wave along a
    phase = 2 * np.pi * 3 * ref.cell_index[:, 0] / 12
    a_hat = lat.matrix[0] / np.linalg.norm(lat.matrix[0])
    cfg = ref.with_cartesian(ref.cartesian + 0.04 * np.cos(phase)[:, None] * a_hat)
    nl = neighbor_list(ref, 2 * 3.0 + 0.1)
    sel = nl.select(translation=(2, 0, 0), symmetric=True)
    res = displacement_correlation(sel, cfg, direction=(1, 0, 0), direction_frame="fractional")
    assert res.correlation == pytest.approx(math.cos(2 * math.pi * 3 * 2 / 12), abs=1e-12)
    np.testing.assert_allclose(res.origin_direction, a_hat, atol=1e-12)


def test_displacement_covariance_matches_an_explicit_loop_and_wrapping_does_not_matter():
    ref = _cell("fcc", ["A"] * 4).repeat((3, 3, 3))
    rng = np.random.default_rng(8)
    u = rng.normal(scale=0.05, size=(ref.n_atoms, 3)) + np.array([0.02, 0.0, -0.01])
    cfg = ref.with_cartesian(ref.cartesian + u)
    wrapped = cfg.wrapped()
    nl = neighbor_list(ref, 2.9).select(origin_site="S2")
    for c in (cfg, wrapped):
        res = displacement_correlation(nl, c, direction=(1, 1, 0), target_direction=(0, 0, 1))
        uo, ut = u[nl.node_row[nl.i]], u[nl.node_row[nl.j]]
        mo, mt = uo.mean(axis=0), ut.mean(axis=0)
        cov = np.zeros((3, 3))
        for p in range(nl.n_pairs):
            cov += np.outer(uo[p] - mo, ut[p] - mt)
        np.testing.assert_allclose(res.covariance, cov / nl.n_pairs, atol=1e-14)
        d1 = np.array([1, 1, 0]) / math.sqrt(2)
        x, y = (uo - mo) @ d1, (ut - mt)[:, 2]
        assert res.correlation == pytest.approx(np.mean(x * y) / math.sqrt(np.mean(x * x) * np.mean(y * y)))
    np.testing.assert_allclose(displacements(nl, cfg), displacements(nl, wrapped), atol=1e-12)


def test_independent_displacements_are_uncorrelated_within_standard_error():
    ref = _cell("sc", ["A"]).repeat((8, 8, 8))
    nl = neighbor_list(ref, 4.1).select(translation=(0, 1, 0), symmetric=True)
    rng = np.random.default_rng(77)
    vals = []
    for _ in range(300):
        cfg = ref.with_cartesian(ref.cartesian + rng.normal(scale=0.05, size=(ref.n_atoms, 3)))
        res = displacement_correlation(nl, cfg, direction=(0, 1, 0))
        vals.append(res.correlation)
    assert res.stderr_null == pytest.approx(1 / math.sqrt(res.n_unique_pairs), rel=0.05)
    assert np.std(vals) == pytest.approx(res.stderr_null, rel=0.15)


def test_uncentered_discus_convention_and_site_average_reference():
    ref, cfg, u = _wave_configuration()
    shifted = cfg.with_cartesian(cfg.cartesian + np.array([0.03, 0.0, 0.0]))   # uniform offset
    sel = neighbor_list(ref, 3.1).select(translation=(1, 0, 0), symmetric=True)
    centred = displacement_correlation(sel, shifted, direction=(1, 0, 0))
    raw = displacement_correlation(sel, shifted, direction=(1, 0, 0), centered=False)
    site_avg = displacement_correlation(sel, shifted, direction=(1, 0, 0), centered=False, reference="site_average")
    expected = math.cos(2 * math.pi * 2 / 16)
    assert centred.correlation == pytest.approx(expected, abs=1e-12)
    assert raw.correlation > expected + 0.1                     # offset makes the uncentred value biased
    assert site_avg.correlation == pytest.approx(expected, abs=1e-12)
    np.testing.assert_allclose(displacements(sel, shifted, reference="site_average").mean(axis=0), 0, atol=1e-15)


def test_vacancies_and_species_filters_in_displacements():
    lat = Lattice.from_parameters(3, 3, 3, 90, 90, 90)
    avg = Structure(lat, np.array(["Ni", "Fe"], object), [[0, 0, 0], [0, 0, 0]], model="average",
                    occupancy=[0.5, 0.4], site_label=["M", "M"]).repeat((6, 6, 6))
    cfg = sample_configuration(avg, seed=2).structure
    rng = np.random.default_rng(3)
    cfg = cfg.with_cartesian(cfg.cartesian + rng.normal(scale=0.03, size=(cfg.n_atoms, 3)))
    nl = neighbor_list(avg, 3.1)
    res = displacement_correlation(nl, cfg, direction=(1, 0, 0))
    spec = node_species(nl, cfg)
    occupied = (spec[nl.i] != "Va") & (spec[nl.j] != "Va")
    assert res.n_pairs == int(occupied.sum()) and res.n_skipped == nl.n_pairs - res.n_pairs
    ni = displacement_correlation(nl, cfg, direction=(1, 0, 0), origin_species="Ni", target_species="Ni")
    assert ni.n_pairs == int(np.sum((spec[nl.i] == "Ni") & (spec[nl.j] == "Ni")))
    assert np.isnan(displacements(nl, cfg)[spec == "Va"]).all()
    none = displacement_correlation(nl, cfg)
    assert math.isnan(none.correlation) and none.undefined == "no direction given"
    with pytest.raises(ValueError, match="target_direction needs direction"):
        displacement_correlation(nl, cfg, target_direction=(1, 0, 0))


def test_displacement_table():
    ref, cfg, _ = _wave_configuration()
    nl = neighbor_list(ref, 6.1)
    rows = displacement_correlation_table(nl, cfg, direction=(1, 0, 0), by="vector")
    by_t = {r["translation"]: r for r in rows}
    for k in (1, 2):
        assert by_t[(k, 0, 0)]["correlation"] == pytest.approx(math.cos(2 * math.pi * 2 * k / 16), abs=1e-12)
        assert by_t[(0, k, 0)]["correlation"] == pytest.approx(1.0, abs=1e-12)
    shells = displacement_correlation_table(nl, cfg, direction=(1, 0, 0), by="shell")
    assert shells[0]["shell"] == 1 and shells[0]["n_pairs"] == ref.n_atoms * 6
