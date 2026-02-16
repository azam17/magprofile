"""Tests for mag.diversity module."""

import numpy as np
import pytest

from mag.diversity import (
    compute_alpha_diversity,
    pielou_evenness,
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
