"""Tests for mag.core_microbiome module."""

import numpy as np
import pytest

from mag.core_microbiome import core_microbiome
from mag.io import AbundanceTable, SampleMetadata
from mag.tests.fixtures import generate_synthetic_abundance_table


class TestCoreMicrobiome:
    def test_basic_thresholds(self):
        table, meta, _ = generate_synthetic_abundance_table(seed=42)
        result = core_microbiome(table, meta, "compartment")
        assert result.thresholds == [0.5, 0.75, 0.9, 1.0]
        for t in result.thresholds:
            assert t in result.core_mags_per_threshold
            assert t in result.shared_across_groups
            for g, mags in result.core_mags_per_threshold[t].items():
                assert isinstance(mags, list)
            assert isinstance(result.shared_across_groups[t], list)

    def test_stricter_threshold_fewer_core(self):
        table, meta, _ = generate_synthetic_abundance_table(seed=42)
        result = core_microbiome(table, meta, "compartment")
        # Stricter thresholds should yield <= core MAGs per group
        for g in ["topsoil", "bulk", "rhizosphere"]:
            n_50 = len(result.core_mags_per_threshold[0.5][g])
            n_100 = len(result.core_mags_per_threshold[1.0][g])
            assert n_100 <= n_50

    def test_shared_subset_of_each_group(self):
        table, meta, _ = generate_synthetic_abundance_table(seed=42)
        result = core_microbiome(table, meta, "compartment")
        for t in result.thresholds:
            shared = set(result.shared_across_groups[t])
            for g, mags in result.core_mags_per_threshold[t].items():
                assert shared.issubset(set(mags))

    def test_custom_thresholds(self):
        table, meta, _ = generate_synthetic_abundance_table(seed=42)
        result = core_microbiome(table, meta, "compartment", thresholds=[0.25, 0.5])
        assert result.thresholds == [0.25, 0.5]
        assert 0.25 in result.core_mags_per_threshold

    def test_all_present_everywhere(self):
        """When all MAGs are present in all samples, all should be core at 100%."""
        rng = np.random.default_rng(99)
        table = AbundanceTable(
            mag_ids=["M0", "M1", "M2"],
            sample_ids=["s0", "s1", "s2", "s3"],
            abundances=rng.uniform(1, 100, (3, 4)),
        )
        meta = SampleMetadata(records={
            "s0": {"g": "A"}, "s1": {"g": "A"},
            "s2": {"g": "B"}, "s3": {"g": "B"},
        })
        result = core_microbiome(table, meta, "g", thresholds=[1.0])
        assert set(result.shared_across_groups[1.0]) == {"M0", "M1", "M2"}

    def test_group_specific_mags(self):
        """MAGs present only in one group should not be shared."""
        table = AbundanceTable(
            mag_ids=["M0", "M1"],
            sample_ids=["s0", "s1", "s2", "s3"],
            abundances=np.array([
                [100.0, 100.0, 0.0, 0.0],  # Only in group A
                [0.0, 0.0, 100.0, 100.0],  # Only in group B
            ]),
        )
        meta = SampleMetadata(records={
            "s0": {"g": "A"}, "s1": {"g": "A"},
            "s2": {"g": "B"}, "s3": {"g": "B"},
        })
        result = core_microbiome(table, meta, "g", thresholds=[0.5, 1.0])
        # M0 core in A, M1 core in B, neither shared
        assert "M0" in result.core_mags_per_threshold[1.0]["A"]
        assert "M1" in result.core_mags_per_threshold[1.0]["B"]
        assert result.shared_across_groups[1.0] == []
        assert result.shared_across_groups[0.5] == []
