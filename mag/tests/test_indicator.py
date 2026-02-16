"""Tests for mag.indicator module."""

import numpy as np
import pytest

from mag.indicator import indicator_species
from mag.tests.fixtures import generate_synthetic_abundance_table


class TestIndicatorSpecies:
    def test_basic(self):
        table, meta, _ = generate_synthetic_abundance_table(seed=42)
        result = indicator_species(table, meta, "compartment", n_permutations=99, seed=42)
        assert len(result.mag_ids) == table.n_mags
        assert np.all(result.indval_scores >= 0)
        assert np.all(result.indval_scores <= 1)
        assert len(result.best_group) == table.n_mags
        assert np.all(result.p_values >= 0)
        assert np.all(result.p_values <= 1)

    def test_specialists_have_high_indval(self):
        table, meta, _ = generate_synthetic_abundance_table(seed=42)
        result = indicator_species(table, meta, "compartment", n_permutations=99, seed=42)
        topsoil_scores = result.indval_scores[:6]
        assert np.mean(topsoil_scores) > 0.3

    def test_specificity_fidelity_range(self):
        table, meta, _ = generate_synthetic_abundance_table(seed=42)
        result = indicator_species(table, meta, "compartment", n_permutations=99, seed=42)
        assert np.all(result.specificity >= 0)
        assert np.all(result.specificity <= 1)
        assert np.all(result.fidelity >= 0)
        assert np.all(result.fidelity <= 1)
