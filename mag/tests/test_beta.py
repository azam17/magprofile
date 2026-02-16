"""Tests for mag.beta module."""

import numpy as np
import pytest

from mag.beta import bray_curtis, jaccard
from mag.io import AbundanceTable
from mag.tests.fixtures import generate_synthetic_abundance_table


class TestBrayCurtis:
    def test_identical_samples(self):
        table = AbundanceTable(
            mag_ids=["a", "b"],
            sample_ids=["s1", "s2"],
            abundances=np.array([[10.0, 10.0], [20.0, 20.0]]),
        )
        result = bray_curtis(table)
        assert result.distance_matrix[0, 1] == pytest.approx(0.0)

    def test_different_samples(self):
        table = AbundanceTable(
            mag_ids=["a", "b"],
            sample_ids=["s1", "s2"],
            abundances=np.array([[100.0, 0.0], [0.0, 100.0]]),
        )
        result = bray_curtis(table)
        assert result.distance_matrix[0, 1] == pytest.approx(1.0)

    def test_symmetry(self):
        table, _, _ = generate_synthetic_abundance_table()
        result = bray_curtis(table)
        np.testing.assert_array_almost_equal(
            result.distance_matrix, result.distance_matrix.T
        )

    def test_diagonal_zero(self):
        table, _, _ = generate_synthetic_abundance_table()
        result = bray_curtis(table)
        np.testing.assert_array_almost_equal(
            np.diag(result.distance_matrix), np.zeros(table.n_samples)
        )


class TestJaccard:
    def test_identical_presence(self):
        table = AbundanceTable(
            mag_ids=["a", "b"],
            sample_ids=["s1", "s2"],
            abundances=np.array([[10.0, 5.0], [20.0, 15.0]]),
        )
        result = jaccard(table)
        assert result.distance_matrix[0, 1] == pytest.approx(0.0)

    def test_no_overlap(self):
        table = AbundanceTable(
            mag_ids=["a", "b"],
            sample_ids=["s1", "s2"],
            abundances=np.array([[100.0, 0.0], [0.0, 100.0]]),
        )
        result = jaccard(table)
        assert result.distance_matrix[0, 1] == pytest.approx(1.0)

    def test_symmetry(self):
        table, _, _ = generate_synthetic_abundance_table()
        result = jaccard(table)
        np.testing.assert_array_almost_equal(
            result.distance_matrix, result.distance_matrix.T
        )
