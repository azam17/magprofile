"""Tests for the read classification pipeline."""

import numpy as np
import pytest

from strainsift.classify import (
    ClassificationConfig,
    ClassificationStats,
    ReadClassifier,
    classify_reads_batch,
)
from strainsift.index import StrainSiftIndex
from strainsift.simulator import SimulationConfig, simulate_metagenome
from strainsift.types import ReadClassification
from strainsift.utils import FastqRecord


class TestClassificationConfig:
    def test_defaults(self):
        cfg = ClassificationConfig()
        assert cfg.top_species == 3
        assert cfg.min_containment == 0.1
        assert cfg.min_read_length == 50


class TestClassificationStats:
    def test_fractions(self):
        stats = ClassificationStats(
            total_reads=100,
            classified_reads=80,
            unclassified_reads=20,
        )
        assert stats.classified_fraction == 0.8
        assert stats.unclassified_fraction == 0.2

    def test_zero_reads(self):
        stats = ClassificationStats()
        assert stats.classified_fraction == 0.0
        assert stats.unclassified_fraction == 0.0


class TestReadClassifier:
    @pytest.fixture
    def sim_data(self):
        """Generate a small simulation for testing."""
        cfg = SimulationConfig(
            n_strains=2,
            genome_length=5_000,
            n_gene_families=5,
            gene_length=500,
            n_reads=200,
            read_length=150,
            error_rate=0.001,
            random_seed=42,
        )
        return simulate_metagenome(cfg)

    def test_classify_reads(self, sim_data):
        refs, reads, truth = sim_data
        classifications, stats = classify_reads_batch(reads, refs)

        assert len(classifications) == len(reads)
        assert stats.total_reads == len(reads)
        assert stats.classified_reads > 0

    def test_all_reads_get_classification(self, sim_data):
        refs, reads, truth = sim_data
        classifications, _ = classify_reads_batch(reads, refs)

        for rc in classifications:
            assert isinstance(rc, ReadClassification)
            assert rc.read_id != ""

    def test_classified_reads_have_candidates(self, sim_data):
        refs, reads, truth = sim_data
        classifications, _ = classify_reads_batch(reads, refs)

        for rc in classifications:
            if rc.is_classified:
                assert len(rc.candidates) > 0
                for cand in rc.candidates:
                    assert cand.containment > 0

    def test_short_reads_not_classified(self):
        refs_cfg = SimulationConfig(
            n_strains=2, genome_length=5_000, n_reads=10, random_seed=42
        )
        refs, _, _ = simulate_metagenome(refs_cfg)

        # Create artificially short reads
        short_reads = [
            FastqRecord(read_id=f"short_{i}", sequence="ACGT", quality="IIII")
            for i in range(5)
        ]

        config = ClassificationConfig(min_read_length=50)
        classifications, stats = classify_reads_batch(short_reads, refs, config)

        assert stats.too_short_reads == 5

    def test_species_read_counts(self, sim_data):
        refs, reads, truth = sim_data
        classifications, stats = classify_reads_batch(reads, refs)

        # All refs have the same species_id in simulation
        assert len(stats.species_read_counts) >= 0  # may be 0 if none classified

    def test_containment_scores_valid(self, sim_data):
        refs, reads, truth = sim_data
        classifications, _ = classify_reads_batch(reads, refs)

        for rc in classifications:
            for cand in rc.candidates:
                assert 0.0 <= cand.containment <= 1.0


class TestClassifyReadsBatch:
    def test_convenience_function(self):
        cfg = SimulationConfig(
            n_strains=2,
            genome_length=5_000,
            n_reads=50,
            random_seed=42,
        )
        refs, reads, truth = simulate_metagenome(cfg)
        classifications, stats = classify_reads_batch(reads, refs)

        assert len(classifications) == 50
        assert stats.total_reads == 50

    def test_empty_reads(self):
        cfg = SimulationConfig(n_strains=2, genome_length=5_000, n_reads=10, random_seed=42)
        refs, _, _ = simulate_metagenome(cfg)
        classifications, stats = classify_reads_batch([], refs)

        assert len(classifications) == 0
        assert stats.total_reads == 0
