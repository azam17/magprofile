"""Tests for the metagenome simulator."""

import numpy as np
import pytest

from strainsift.simulator import SimulationConfig, simulate_metagenome
from strainsift.types import StrainReference, SimulationGroundTruth
from strainsift.utils import FastqRecord


class TestSimulationConfig:
    def test_default_config(self):
        cfg = SimulationConfig()
        assert cfg.n_strains == 3
        assert cfg.genome_length == 50_000
        assert cfg.n_gene_families == 10
        assert cfg.read_length == 150
        assert cfg.error_rate == 0.005
        assert cfg.random_seed == 42


class TestSimulateMetagenome:
    def test_basic_simulation(self):
        cfg = SimulationConfig(n_strains=2, n_reads=100, random_seed=42)
        refs, reads, truth = simulate_metagenome(cfg)

        assert len(refs) == 2
        assert len(reads) == 100
        assert isinstance(truth, SimulationGroundTruth)

    def test_returns_correct_types(self):
        cfg = SimulationConfig(n_strains=2, n_reads=50, random_seed=42)
        refs, reads, truth = simulate_metagenome(cfg)

        assert all(isinstance(r, StrainReference) for r in refs)
        assert all(isinstance(r, FastqRecord) for r in reads)
        assert isinstance(truth.true_pi, np.ndarray)
        assert isinstance(truth.true_phi, np.ndarray)

    def test_deterministic_output(self):
        cfg = SimulationConfig(n_strains=2, n_reads=50, random_seed=42)
        _, reads1, truth1 = simulate_metagenome(cfg)
        _, reads2, truth2 = simulate_metagenome(cfg)

        assert reads1[0].sequence == reads2[0].sequence
        np.testing.assert_array_equal(truth1.true_pi, truth2.true_pi)

    def test_different_seeds_differ(self):
        cfg1 = SimulationConfig(n_strains=2, n_reads=50, random_seed=42)
        cfg2 = SimulationConfig(n_strains=2, n_reads=50, random_seed=99)
        _, reads1, _ = simulate_metagenome(cfg1)
        _, reads2, _ = simulate_metagenome(cfg2)

        # Very unlikely to be identical
        assert reads1[0].sequence != reads2[0].sequence

    def test_ground_truth_shapes(self):
        cfg = SimulationConfig(
            n_strains=3, n_gene_families=5, n_reads=100, random_seed=42
        )
        refs, reads, truth = simulate_metagenome(cfg)

        assert truth.true_pi.shape == (3,)
        assert np.isclose(truth.true_pi.sum(), 1.0)
        assert len(truth.strain_ids) == 3
        # phi shape: (n_strains, n_placed_genes) — n_placed may be <= n_gene_families
        assert truth.true_phi.shape[0] == 3
        assert truth.true_phi.shape[1] <= 5

    def test_read_assignments(self):
        cfg = SimulationConfig(n_strains=2, n_reads=100, random_seed=42)
        refs, reads, truth = simulate_metagenome(cfg)

        assert len(truth.read_assignments) == 100
        for read in reads:
            assert read.read_id in truth.read_assignments

    def test_read_quality_valid(self):
        cfg = SimulationConfig(n_strains=2, n_reads=50, random_seed=42)
        _, reads, _ = simulate_metagenome(cfg)

        for read in reads:
            assert len(read.quality) == len(read.sequence)
            # Quality scores should be valid Phred+33
            for q_char in read.quality:
                phred = ord(q_char) - 33
                assert 0 <= phred <= 60

    def test_abundances_sum_to_one(self):
        cfg = SimulationConfig(n_strains=3, n_reads=50, random_seed=42)
        _, _, truth = simulate_metagenome(cfg)

        assert np.isclose(truth.true_pi.sum(), 1.0)

    def test_custom_abundances(self):
        cfg = SimulationConfig(
            n_strains=3,
            abundances=np.array([0.5, 0.3, 0.2]),
            n_reads=50,
            random_seed=42,
        )
        _, _, truth = simulate_metagenome(cfg)

        np.testing.assert_allclose(truth.true_pi, [0.5, 0.3, 0.2])

    def test_open_universe(self):
        cfg = SimulationConfig(
            n_strains=4,
            open_universe_fraction=0.5,
            n_reads=100,
            random_seed=42,
        )
        refs, reads, truth = simulate_metagenome(cfg)

        # 2 of 4 strains should be hidden from DB
        assert len(refs) == 2
        assert len(truth.strains_not_in_db) == 2
        assert len(truth.strains_in_db) == 2
        assert len(truth.strain_ids) == 4  # truth knows about all strains

    def test_mobile_genes(self):
        cfg = SimulationConfig(
            n_strains=2,
            n_gene_families=10,
            mobile_gene_fraction=0.3,
            n_reads=50,
            random_seed=42,
        )
        _, _, truth = simulate_metagenome(cfg)

        assert len(truth.mobile_genes) >= 1  # at least some mobile genes

    def test_genome_length_respected(self):
        cfg = SimulationConfig(
            n_strains=2,
            genome_length=10_000,
            n_reads=50,
            random_seed=42,
        )
        refs, _, _ = simulate_metagenome(cfg)

        for ref in refs:
            assert ref.genome_length == 10_000

    def test_zero_error_rate(self):
        cfg = SimulationConfig(
            n_strains=2,
            n_reads=50,
            error_rate=0.0,
            random_seed=42,
        )
        refs, reads, truth = simulate_metagenome(cfg)

        # With 0 error rate, reads should be exact substrings of genomes
        ref_genomes = {r.strain_id: r.genome for r in refs}
        for read in reads[:10]:
            strain = truth.read_assignments[read.read_id]
            if strain in ref_genomes:
                genome = ref_genomes[strain]
                # Check that the read is a substring (may wrap around)
                assert read.sequence in genome or read.sequence in genome + genome

    def test_gene_annotations_present(self):
        cfg = SimulationConfig(
            n_strains=2,
            n_gene_families=5,
            n_reads=50,
            random_seed=42,
        )
        refs, _, _ = simulate_metagenome(cfg)

        for ref in refs:
            assert len(ref.gene_annotations) > 0
            for ann in ref.gene_annotations:
                assert ann.start >= 0
                assert ann.end > ann.start
                assert ann.end <= ref.genome_length

    def test_serialization(self, tmp_path):
        cfg = SimulationConfig(n_strains=2, n_reads=50, random_seed=42)
        _, _, truth = simulate_metagenome(cfg)

        path = tmp_path / "truth.json"
        truth.save(path)
        assert path.exists()
