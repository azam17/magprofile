"""Tests for mag.ordination module."""

import numpy as np
import pytest

from mag.beta import bray_curtis
from mag.ordination import nmds, pcoa
from mag.tests.fixtures import generate_synthetic_abundance_table


class TestPCoA:
    def test_basic(self):
        table, _, _ = generate_synthetic_abundance_table()
        bc = bray_curtis(table)
        result = pcoa(bc)
        assert result.coordinates.shape == (table.n_samples, 2)
        assert result.explained_variance is not None
        assert len(result.explained_variance) == 2
        assert np.all(result.explained_variance >= 0)
        assert result.method == "PCoA"

    def test_explained_variance_sums_to_one_or_less(self):
        table, _, _ = generate_synthetic_abundance_table()
        bc = bray_curtis(table)
        result = pcoa(bc, n_axes=5)
        assert result.explained_variance.sum() <= 1.0 + 1e-10

    def test_two_identical_samples(self):
        from mag.io import AbundanceTable

        table = AbundanceTable(
            mag_ids=["a", "b"],
            sample_ids=["s1", "s2", "s3"],
            abundances=np.array([[10.0, 10.0, 100.0], [20.0, 20.0, 5.0]]),
        )
        bc = bray_curtis(table)
        result = pcoa(bc)
        np.testing.assert_array_almost_equal(
            result.coordinates[0], result.coordinates[1]
        )


class TestNMDS:
    def test_basic(self):
        table, _, _ = generate_synthetic_abundance_table()
        bc = bray_curtis(table)
        result = nmds(bc)
        assert result.coordinates.shape == (table.n_samples, 2)
        assert result.stress is not None
        assert result.method == "NMDS"
