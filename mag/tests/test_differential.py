"""Tests for mag.differential module."""

import numpy as np
import pytest

from mag.differential import _bh_fdr, clr_transform, cohens_d, differential_abundance
from mag.io import AbundanceTable, SampleMetadata
from mag.tests.fixtures import generate_synthetic_abundance_table


class TestCLRTransform:
    def test_basic(self):
        x = np.array([[10.0, 20.0], [30.0, 40.0]])
        clr = clr_transform(x, pseudocount=0)
        np.testing.assert_allclose(clr.mean(axis=0), [0, 0], atol=1e-10)

    def test_with_zeros(self):
        x = np.array([[0.0, 100.0], [100.0, 0.0]])
        clr = clr_transform(x)
        assert np.all(np.isfinite(clr))


class TestCohensD:
    def test_identical_groups(self):
        g = np.array([1.0, 2.0, 3.0])
        assert cohens_d(g, g) == pytest.approx(0.0)

    def test_large_effect(self):
        g1 = np.array([10.0, 11.0, 12.0])
        g2 = np.array([0.0, 1.0, 2.0])
        d = cohens_d(g1, g2)
        assert d > 5.0


class TestBHFDR:
    def test_basic(self):
        pvals = np.array([0.01, 0.04, 0.03, 0.5])
        qvals = _bh_fdr(pvals)
        assert np.all(qvals >= pvals)
        assert np.all(qvals <= 1.0)

    def test_monotonicity(self):
        pvals = np.array([0.001, 0.01, 0.05, 0.1, 0.5])
        qvals = _bh_fdr(pvals)
        sorted_q = qvals[np.argsort(pvals)]
        assert np.all(np.diff(sorted_q) >= -1e-10)


class TestDifferentialAbundance:
    def test_synthetic_data(self):
        table, meta, _ = generate_synthetic_abundance_table(seed=42)
        result = differential_abundance(table, meta, "compartment", "topsoil", "bulk")
        assert len(result.mag_ids) == table.n_mags
        assert result.log_fold_changes.shape == (table.n_mags,)
        assert np.all(result.q_values >= 0)
        assert np.all(result.q_values <= 1)

    def test_empty_group_raises(self):
        table, meta, _ = generate_synthetic_abundance_table()
        with pytest.raises(ValueError):
            differential_abundance(table, meta, "compartment", "topsoil", "nonexistent")
