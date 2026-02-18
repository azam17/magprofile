"""Tests for mag.net_correlation module."""

import numpy as np
import pytest

from mag.io import AbundanceTable
from mag.net_correlation import (
    NetworkResult,
    ThresholdSensitivityResult,
    build_network,
    proportionality,
    threshold_sensitivity,
)
from mag.tests.fixtures import (
    generate_edge_case_single_mag,
    generate_synthetic_abundance_table,
)


class TestProportionality:
    def test_output_shape(self):
        table, _, _ = generate_synthetic_abundance_table()
        result = proportionality(table)
        assert result.shape == (table.n_mags, table.n_mags)

    def test_symmetric(self):
        table, _, _ = generate_synthetic_abundance_table()
        result = proportionality(table)
        # NaN == NaN would be False, so use assert_array_equal which handles NaN
        np.testing.assert_array_equal(result, result.T)

    def test_diagonal_zero(self):
        table, _, _ = generate_synthetic_abundance_table()
        result = proportionality(table)
        diag = np.diag(result)
        # Valid MAGs have 0 on diagonal, masked MAGs have NaN
        for val in diag:
            assert val == 0.0 or np.isnan(val)

    def test_proportional_mags_low_phi(self):
        """MAGs with proportional abundances should have lower phi."""
        rng = np.random.default_rng(123)
        n_samples = 50
        # M0: base signal
        base = rng.poisson(100, n_samples).astype(np.float64)
        base = np.maximum(base, 1.0)
        # M1: proportional to M0 (scaled by 3x with small noise)
        m1 = base * 3.0 + rng.normal(0, 2, n_samples)
        m1 = np.maximum(m1, 1.0)
        # M2: independent random pattern
        m2 = rng.poisson(200, n_samples).astype(np.float64)
        m2 = np.maximum(m2, 1.0)

        table = AbundanceTable(
            mag_ids=["M0", "M1", "M2"],
            sample_ids=[f"s{i}" for i in range(n_samples)],
            abundances=np.vstack([base, m1, m2]),
        )
        result = proportionality(table)
        phi_01 = result[0, 1]  # proportional pair
        phi_02 = result[0, 2]  # unrelated pair
        assert phi_01 < phi_02, (
            f"Proportional MAGs should have lower phi: "
            f"phi(M0,M1)={phi_01:.4f} vs phi(M0,M2)={phi_02:.4f}"
        )

    def test_min_prevalence_filter(self):
        """MAGs present in too few samples should be masked with NaN."""
        # 3 MAGs, 10 samples.  M2 is present in only 2/10 = 20% of samples.
        abundances = np.array([
            [10, 20, 30, 40, 50, 60, 70, 80, 90, 100],  # M0: 100% prevalence
            [5, 15, 25, 35, 45, 55, 65, 75, 85, 95],     # M1: 100% prevalence
            [0, 0, 0, 0, 0, 0, 0, 0, 10, 20],            # M2: 20% prevalence
        ], dtype=np.float64)
        table = AbundanceTable(
            mag_ids=["M0", "M1", "M2"],
            sample_ids=[f"s{i}" for i in range(10)],
            abundances=abundances,
        )
        result = proportionality(table, min_prevalence=0.5)
        # M0-M1 should be a valid number
        assert not np.isnan(result[0, 1])
        # M2 row/col should be NaN (except possibly diagonal)
        assert np.isnan(result[0, 2])
        assert np.isnan(result[1, 2])
        assert np.isnan(result[2, 0])
        assert np.isnan(result[2, 1])


class TestBuildNetwork:
    def test_returns_network_result(self):
        table, _, _ = generate_synthetic_abundance_table()
        phi = proportionality(table)
        net = build_network(phi, table.mag_ids)
        assert isinstance(net, NetworkResult)
        assert net.adjacency.shape == (table.n_mags, table.n_mags)
        assert isinstance(net.edges, list)
        assert isinstance(net.threshold, float)
        assert net.mag_ids == table.mag_ids

    def test_edges_below_threshold(self):
        table, _, _ = generate_synthetic_abundance_table()
        phi = proportionality(table)
        net = build_network(phi, table.mag_ids, threshold_percentile=10)
        for mag_a, mag_b, phi_val in net.edges:
            assert phi_val <= net.threshold, (
                f"Edge ({mag_a}, {mag_b}) has phi={phi_val:.4f} "
                f"> threshold={net.threshold:.4f}"
            )

    def test_no_self_edges(self):
        table, _, _ = generate_synthetic_abundance_table()
        phi = proportionality(table)
        net = build_network(phi, table.mag_ids)
        for mag_a, mag_b, _ in net.edges:
            assert mag_a != mag_b, f"Self-edge found: ({mag_a}, {mag_b})"
        # Also check adjacency diagonal
        np.testing.assert_array_equal(np.diag(net.adjacency), np.zeros(table.n_mags))

    def test_edge_completeness(self):
        """Every valid pair at or below threshold must appear in edges."""
        table, _, _ = generate_synthetic_abundance_table()
        phi = proportionality(table)
        net = build_network(phi, table.mag_ids, threshold_percentile=10)
        n = len(table.mag_ids)
        edge_set = {(a, b) for a, b, _ in net.edges}
        for i in range(n):
            for j in range(i + 1, n):
                v = phi[i, j]
                pair = (table.mag_ids[i], table.mag_ids[j])
                if not np.isnan(v) and v <= net.threshold:
                    assert pair in edge_set, (
                        f"Pair {pair} with phi={v:.6f} <= threshold "
                        f"{net.threshold:.6f} missing from edges"
                    )
                else:
                    assert pair not in edge_set

    def test_adjacency_symmetric(self):
        """Adjacency matrix must be symmetric."""
        table, _, _ = generate_synthetic_abundance_table()
        phi = proportionality(table)
        net = build_network(phi, table.mag_ids)
        np.testing.assert_array_equal(net.adjacency, net.adjacency.T)


class TestProportionalityExactValues:
    """Verify vectorized phi matches hand-computed var(log(x_i/x_j), ddof=1)."""

    def test_exact_phi_3x4(self):
        """3 MAGs x 4 samples — compare vectorized vs loop-based phi."""
        abundances = np.array([
            [10.0, 20.0, 30.0, 40.0],
            [20.0, 40.0, 60.0, 80.0],
            [5.0, 100.0, 3.0, 50.0],
        ])
        table = AbundanceTable(
            mag_ids=["A", "B", "C"],
            sample_ids=["s0", "s1", "s2", "s3"],
            abundances=abundances,
        )
        result = proportionality(table, pseudocount=1.0)

        # Hand-compute expected phi: var(log((x_i+1)/(x_j+1)), ddof=1)
        log_x = np.log(abundances + 1.0)
        for i in range(3):
            for j in range(3):
                if i == j:
                    assert result[i, j] == 0.0
                else:
                    expected = float(np.var(log_x[i] - log_x[j], ddof=1))
                    np.testing.assert_allclose(
                        result[i, j], expected, rtol=1e-10,
                        err_msg=f"phi[{i},{j}] mismatch",
                    )

    def test_single_mag(self):
        """Single MAG should produce 1x1 matrix with phi=0."""
        table = generate_edge_case_single_mag()
        result = proportionality(table)
        assert result.shape == (1, 1)
        assert result[0, 0] == 0.0

    def test_phi_nonnegative(self):
        """Phi (variance of log-ratios) must be >= 0 for all valid entries."""
        table, _, _ = generate_synthetic_abundance_table()
        result = proportionality(table)
        valid = ~np.isnan(result)
        assert np.all(result[valid] >= -1e-15), (
            f"Negative phi found: min={result[valid].min()}"
        )


class TestThresholdSensitivity:
    def test_sensitivity_returns_result(self):
        table, _, _ = generate_synthetic_abundance_table()
        phi = proportionality(table)
        result = threshold_sensitivity(phi, table.mag_ids)
        assert isinstance(result, ThresholdSensitivityResult)
        assert len(result.percentiles) == 20  # default 1..20
        assert result.thresholds.shape == result.percentiles.shape
        assert result.n_edges.shape == result.percentiles.shape
        assert result.modularities.shape == result.percentiles.shape
        assert result.n_modules.shape == result.percentiles.shape

    def test_sensitivity_monotonic_edges(self):
        """More permissive threshold (higher percentile) -> more or equal edges."""
        table, _, _ = generate_synthetic_abundance_table()
        phi = proportionality(table)
        result = threshold_sensitivity(phi, table.mag_ids)
        for i in range(1, len(result.n_edges)):
            assert result.n_edges[i] >= result.n_edges[i - 1], (
                f"Edge count decreased at percentile {result.percentiles[i]}: "
                f"{result.n_edges[i]} < {result.n_edges[i - 1]}"
            )

    def test_sensitivity_custom_percentiles(self):
        table, _, _ = generate_synthetic_abundance_table()
        phi = proportionality(table)
        percs = np.array([5.0, 10.0, 15.0])
        result = threshold_sensitivity(phi, table.mag_ids, percentiles=percs)
        assert len(result.percentiles) == 3
        np.testing.assert_array_equal(result.percentiles, percs)
