"""Tests for mag.diversity module."""

import numpy as np
import pytest

from mag.diversity import (
    RarefactionResult,
    compute_alpha_diversity,
    pielou_evenness,
    rarefaction_curves,
    richness,
    shannon_index,
    simpson_index,
)
from mag.io import AbundanceTable
from mag.tests.fixtures import (
    generate_edge_case_all_zeros,
    generate_edge_case_single_mag,
    generate_synthetic_abundance_table,
)


class TestShannonIndex:
    def test_uniform_distribution(self):
        abundances = np.array([25.0, 25.0, 25.0, 25.0])
        assert abs(shannon_index(abundances) - np.log(4)) < 1e-10

    def test_single_species(self):
        abundances = np.array([100.0, 0.0, 0.0])
        assert shannon_index(abundances) == 0.0

    def test_zero_vector(self):
        assert shannon_index(np.zeros(5)) == 0.0


class TestSimpsonIndex:
    def test_uniform_distribution(self):
        abundances = np.array([25.0, 25.0, 25.0, 25.0])
        assert abs(simpson_index(abundances) - 0.75) < 1e-10

    def test_single_species(self):
        abundances = np.array([100.0, 0.0, 0.0])
        assert simpson_index(abundances) == 0.0

    def test_zero_vector(self):
        assert simpson_index(np.zeros(5)) == 0.0


class TestRichness:
    def test_basic(self):
        assert richness(np.array([1.0, 0.0, 5.0, 0.0])) == 2

    def test_all_zeros(self):
        assert richness(np.zeros(5)) == 0


class TestPielouEvenness:
    def test_uniform(self):
        abundances = np.array([25.0, 25.0, 25.0, 25.0])
        assert abs(pielou_evenness(abundances) - 1.0) < 1e-10

    def test_single_species(self):
        assert pielou_evenness(np.array([100.0])) == 0.0

    def test_zero_vector(self):
        assert pielou_evenness(np.zeros(5)) == 0.0


class TestComputeAlphaDiversity:
    def test_synthetic_data(self):
        table, _, _ = generate_synthetic_abundance_table()
        result = compute_alpha_diversity(table)
        assert len(result.sample_ids) == table.n_samples
        assert result.shannon.shape == (table.n_samples,)
        assert np.all(result.shannon >= 0)
        assert np.all(result.simpson >= 0)
        assert np.all(result.simpson <= 1)

    def test_all_zeros(self):
        table = generate_edge_case_all_zeros()
        result = compute_alpha_diversity(table)
        np.testing.assert_array_equal(result.shannon, np.zeros(3))

    def test_single_mag(self):
        table = generate_edge_case_single_mag()
        result = compute_alpha_diversity(table)
        assert np.all(result.richness == 1)
        np.testing.assert_array_equal(result.shannon, np.zeros(3))


class TestRarefactionCurves:
    @pytest.fixture
    def abundance(self):
        rng = np.random.default_rng(42)
        return AbundanceTable(
            mag_ids=[f"M{i}" for i in range(10)],
            sample_ids=["s0", "s1", "s2"],
            abundances=rng.poisson(50, (10, 3)).astype(float),
        )

    def test_output_shape(self, abundance):
        result = rarefaction_curves(abundance, n_depths=10, n_iterations=5, seed=42)
        assert isinstance(result, RarefactionResult)
        n_d = len(result.depths)
        assert result.richness.shape == (n_d, 3)
        assert result.shannon.shape == (n_d, 3)
        assert result.richness_sd.shape == (n_d, 3)
        assert result.shannon_sd.shape == (n_d, 3)

    def test_richness_monotonic(self, abundance):
        result = rarefaction_curves(abundance, n_depths=15, n_iterations=20, seed=42)
        for j in range(abundance.n_samples):
            r = result.richness[:, j]
            valid = ~np.isnan(r)
            r_valid = r[valid]
            for i in range(1, len(r_valid)):
                assert r_valid[i] >= r_valid[i - 1] - 0.5, (
                    f"Richness decreased at sample {j}, depth index {i}"
                )

    def test_deterministic(self, abundance):
        r1 = rarefaction_curves(abundance, n_depths=5, n_iterations=5, seed=99)
        r2 = rarefaction_curves(abundance, n_depths=5, n_iterations=5, seed=99)
        np.testing.assert_array_equal(r1.richness, r2.richness)
        np.testing.assert_array_equal(r1.shannon, r2.shannon)

    def test_max_depth_matches_full(self, abundance):
        result = rarefaction_curves(abundance, n_depths=20, n_iterations=30, seed=42)
        for j in range(abundance.n_samples):
            full_richness = richness(abundance.abundances[:, j])
            # At max depth, mean richness should be close to true richness
            last_valid = result.richness[:, j]
            valid = ~np.isnan(last_valid)
            if valid.any():
                at_max = last_valid[valid][-1]
                assert abs(at_max - full_richness) <= 2.0, (
                    f"Sample {j}: rarefied richness {at_max} vs true {full_richness}"
                )
