"""Tests for mag.compartment_specificity module."""

import numpy as np
import pytest

from mag.compartment_specificity import compartment_specificity
from mag.io import AbundanceTable, SampleMetadata
from mag.tests.fixtures import generate_synthetic_abundance_table


class TestCompartmentSpecificity:
    def test_basic(self):
        table, meta, _ = generate_synthetic_abundance_table(seed=42)
        result = compartment_specificity(table, meta, "compartment")
        assert len(result.mag_ids) == table.n_mags
        assert np.all(result.scores >= 0)
        assert np.all(result.scores <= 1.0 + 1e-10)
        assert len(result.dominant_compartment) == table.n_mags

    def test_perfect_specialist(self):
        table = AbundanceTable(
            mag_ids=["a"],
            sample_ids=["s1", "s2", "s3"],
            abundances=np.array([[100.0, 0.0, 0.0]]),
        )
        meta = SampleMetadata(
            records={
                "s1": {"grp": "A"},
                "s2": {"grp": "B"},
                "s3": {"grp": "C"},
            }
        )
        result = compartment_specificity(table, meta, "grp")
        assert result.scores[0] == pytest.approx(1.0)
        assert result.dominant_compartment[0] == "A"

    def test_perfect_generalist(self):
        table = AbundanceTable(
            mag_ids=["a"],
            sample_ids=["s1", "s2", "s3"],
            abundances=np.array([[100.0, 100.0, 100.0]]),
        )
        meta = SampleMetadata(
            records={
                "s1": {"grp": "A"},
                "s2": {"grp": "B"},
                "s3": {"grp": "C"},
            }
        )
        result = compartment_specificity(table, meta, "grp")
        assert result.scores[0] == pytest.approx(0.0, abs=1e-10)

    def test_specialists_vs_generalists(self):
        table, meta, _ = generate_synthetic_abundance_table(seed=42)
        result = compartment_specificity(table, meta, "compartment")
        specialist_scores = result.scores[:6].mean()
        generalist_scores = result.scores[14:].mean()
        assert specialist_scores > generalist_scores
