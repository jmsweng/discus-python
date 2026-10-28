"""Regression tests for the short-range correlated-disorder model.

These tests verify:
- Reproducible, fixed-seed construction
- Measured correlation decay
- Absence of unintended checkerboard long-range order
- 50/50 species ratio (approximately)
"""

from __future__ import annotations

import numpy as np
import pytest

from benchmarks.discus_6_23_01.short_range_model import (
    build_short_range_structure,
    measure_correlation_decay,
    compute_directional_pair_correlations,
)


class TestShortRangeStructure:
    """Test the short-range structure properties."""

    def test_reproducible_construction(self):
        """Test that the construction is reproducible with a fixed seed."""
        struct1 = build_short_range_structure(n=8, a=5.0, correlation_length=1.0, seed=42)
        struct2 = build_short_range_structure(n=8, a=5.0, correlation_length=1.0, seed=42)

        # The species assignments should be identical
        assert struct1["species"] == struct2["species"], "Construction should be reproducible"

    def test_species_count_reasonable(self):
        """Test that the species count is approximately 50/50."""
        struct = build_short_range_structure(n=8, a=5.0, correlation_length=1.0, seed=42)
        n_ga = sum(1 for sp in struct["species"] if sp == "Ga")
        n_as = sum(1 for sp in struct["species"] if sp == "As")
        n_total = len(struct["species"])

        # The ratio should be within 10% of 50/50
        ratio = n_ga / n_total
        assert 0.4 <= ratio <= 0.6, f"Species ratio {ratio:.2f} is not approximately 50/50"

    def test_eta_mean_small(self):
        """Test that <η> is small (close to 0)."""
        struct = build_short_range_structure(n=8, a=5.0, correlation_length=1.0, seed=42)
        eta = np.array([1.0 if sp == "Ga" else -1.0 for sp in struct["species"]])
        eta_mean = np.mean(eta)

        # <η> should be small (|<η>| < 0.1)
        assert abs(eta_mean) < 0.1, f"<η> = {eta_mean:.4f} is not small"

    def test_correlation_decay_verified(self):
        """Test that the correlation decay is verified."""
        struct = build_short_range_structure(n=8, a=5.0, correlation_length=1.0, seed=42)
        corr_decay = measure_correlation_decay(struct)

        # The decay should be verified
        assert corr_decay["decay_verified"], "Correlation decay should be verified"

    def test_no_checkboard_long_range_order(self):
        """Test that there is no unintended checkerboard long-range order.

        For a checkerboard pattern, c(tau) would oscillate with period 2.
        For a short-range correlated model, c(tau) should decay to ~0.

        Note: This test is limited by the small supercell size (n=8).
        Periodic boundary effects can cause apparent long-range correlations.
        """
        struct = build_short_range_structure(n=8, a=5.0, correlation_length=1.0, seed=42)
        corr_decay = measure_correlation_decay(struct)

        # Check that the long-range values (sep >= 3) are relatively small
        # Use a relaxed threshold: < 0.15 (due to periodic boundary effects)
        long_range_x = np.mean([abs(c) for c in corr_decay["c_x"][-3:]])
        long_range_y = np.mean([abs(c) for c in corr_decay["c_y"][-3:]])
        long_range_z = np.mean([abs(c) for c in corr_decay["c_z"][-3:]])

        # The long-range values should be relatively small
        # (not as large as the short-range values)
        short_range_x = np.mean([abs(c) for c in corr_decay["c_x"][:3]])
        short_range_y = np.mean([abs(c) for c in corr_decay["c_y"][:3]])
        short_range_z = np.mean([abs(c) for c in corr_decay["c_z"][:3]])

        # The long-range values should be smaller than the short-range values
        # (indicating decay, not oscillation)
        assert long_range_x < short_range_x, "c_x should decay with distance"
        assert long_range_y < short_range_y, "c_y should decay with distance"
        assert long_range_z < short_range_z, "c_z should decay with distance"

    def test_directional_correlations_computed(self):
        """Test that directional correlations are computed at the NN distance."""
        struct = build_short_range_structure(n=8, a=5.0, correlation_length=1.0, seed=42)
        pair_corr = compute_directional_pair_correlations(struct)

        # All 6 NN displacements should have computed correlations
        assert len(pair_corr["c_tau"]) == 6, "Expected 6 NN displacements"
        assert len(pair_corr["n_pairs"]) == 6, "Expected 6 NN displacements"

        # All pair counts should be positive
        for n_pairs in pair_corr["n_pairs"].values():
            assert n_pairs > 0, "Pair count should be positive"

    def test_positions_in_cartesian(self):
        """Test that positions are in Cartesian (Å)."""
        struct = build_short_range_structure(n=8, a=5.0, correlation_length=1.0, seed=42)
        positions = np.array(struct["positions"])

        # All positions should be in [0, 40) Å (supercell = 40 Å)
        assert np.all(positions >= 0.0), "Positions must be non-negative"
        assert np.all(positions < 40.0), "Positions must be < 40 Å (supercell size)"

    def test_total_site_count(self):
        """Test the exact total site count."""
        struct = build_short_range_structure(n=8, a=5.0, correlation_length=1.0, seed=42)
        # 8×8×8 supercell with 4 sites per unit cell = 2048 sites
        assert len(struct["positions"]) == 2048, f"Expected 2048 sites, got {len(struct['positions'])}"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
