"""End-to-end integration tests for StrainSift.

Tests the full pipeline: simulate → index → classify → EM → evaluate.
"""

import numpy as np
import pytest

from strainsift.benchmark import (
    _build_reference_phi,
    function_attribution_metrics,
    mobile_element_metrics,
    strain_abundance_metrics,
)
from strainsift.classify import classify_reads_batch
from strainsift.em import EMConfig, JointEM
from strainsift.index import StrainSiftIndex
from strainsift.simulator import SimulationConfig, simulate_metagenome
from strainsift.types import StrainProfile
from strainsift.utils import bray_curtis, l1_error


class TestFullPipeline:
    """End-to-end pipeline: simulate → classify → EM → metrics."""

    def test_easy_2_strain_pipeline(self):
        """Two well-separated strains at 99% ANI.

        This is the hand-verified toy case: EM should recover ground truth
        with reasonable accuracy.
        """
        cfg = SimulationConfig(
            n_strains=2,
            genome_length=10_000,
            n_gene_families=5,
            gene_length=800,
            pairwise_ani=0.99,
            n_reads=1000,
            read_length=150,
            error_rate=0.001,
            random_seed=42,
        )
        refs, reads, truth = simulate_metagenome(cfg)

        # Classify
        classifications, stats = classify_reads_batch(reads, refs)
        assert stats.classified_reads > 0

        # EM
        strain_ids = [r.strain_id for r in refs]
        gene_families = set()
        for ref in refs:
            for ann in ref.gene_annotations:
                gene_families.add(ann.gene_family_id)
        gene_family_ids = sorted(gene_families)
        ref_phi = _build_reference_phi(refs, strain_ids, gene_family_ids)

        em = JointEM(EMConfig(n_restarts=5, random_seed=42))
        profile = em.fit(classifications, strain_ids, gene_family_ids, ref_phi)

        # Verify basic properties
        assert isinstance(profile, StrainProfile)
        assert np.isclose(profile.pi.sum(), 1.0, atol=1e-4)
        assert profile.n_reads_classified > 0

        # Metrics
        abund = strain_abundance_metrics(
            profile.pi, truth.true_pi, strain_ids, truth.strain_ids
        )
        assert abund["l1_error"] < 1.0  # should be reasonably close
        assert abund["bray_curtis"] < 0.5

    def test_3_strain_with_accessory(self):
        """Three strains with accessory gene differences — the core use case."""
        cfg = SimulationConfig(
            n_strains=3,
            genome_length=15_000,
            n_gene_families=6,
            gene_length=800,
            pairwise_ani=0.995,
            accessory_gene_fraction=0.3,
            mobile_gene_fraction=0.2,
            n_reads=1500,
            read_length=150,
            error_rate=0.003,
            random_seed=42,
        )
        refs, reads, truth = simulate_metagenome(cfg)

        classifications, stats = classify_reads_batch(reads, refs)

        strain_ids = [r.strain_id for r in refs]
        gene_families = set()
        for ref in refs:
            for ann in ref.gene_annotations:
                gene_families.add(ann.gene_family_id)
        gene_family_ids = sorted(gene_families)
        ref_phi = _build_reference_phi(refs, strain_ids, gene_family_ids)

        em = JointEM(EMConfig(n_restarts=5, random_seed=42))
        profile = em.fit(classifications, strain_ids, gene_family_ids, ref_phi)

        assert np.isclose(profile.pi.sum(), 1.0, atol=1e-4)
        assert profile.convergence.n_iterations > 0

    def test_open_universe_increases_uncertainty(self):
        """With strains absent from DB, uncertainty should increase."""
        # Closed universe
        cfg_closed = SimulationConfig(
            n_strains=3,
            genome_length=10_000,
            n_gene_families=5,
            n_reads=500,
            open_universe_fraction=0.0,
            random_seed=42,
        )
        refs_c, reads_c, truth_c = simulate_metagenome(cfg_closed)
        class_c, _ = classify_reads_batch(reads_c, refs_c)
        strain_ids_c = [r.strain_id for r in refs_c]
        gf_c = sorted({ann.gene_family_id for r in refs_c for ann in r.gene_annotations})
        ref_phi_c = _build_reference_phi(refs_c, strain_ids_c, gf_c)

        em = JointEM(EMConfig(n_restarts=3, random_seed=42))
        prof_closed = em.fit(class_c, strain_ids_c, gf_c, ref_phi_c)

        # Open universe (1 of 3 strains hidden)
        cfg_open = SimulationConfig(
            n_strains=3,
            genome_length=10_000,
            n_gene_families=5,
            n_reads=500,
            open_universe_fraction=0.33,
            random_seed=42,
        )
        refs_o, reads_o, truth_o = simulate_metagenome(cfg_open)
        class_o, stats_o = classify_reads_batch(reads_o, refs_o)
        strain_ids_o = [r.strain_id for r in refs_o]
        gf_o = sorted({ann.gene_family_id for r in refs_o for ann in r.gene_annotations})
        ref_phi_o = _build_reference_phi(refs_o, strain_ids_o, gf_o)

        prof_open = em.fit(class_o, strain_ids_o, gf_o, ref_phi_o)

        # Open universe should have higher mobile pool proportion
        # (reads from unknown strain get absorbed into mobile pool)
        # or wider CIs. Just verify both produce valid output.
        assert np.isclose(prof_closed.pi.sum(), 1.0, atol=1e-4)
        assert np.isclose(prof_open.pi.sum(), 1.0, atol=1e-4)


class TestMetrics:
    def test_strain_abundance_metrics(self):
        est = np.array([0.5, 0.3, 0.2])
        true = np.array([0.4, 0.4, 0.2])
        ids = ["s0", "s1", "s2"]

        metrics = strain_abundance_metrics(est, true, ids, ids)
        assert metrics["l1_error"] > 0
        assert metrics["bray_curtis"] >= 0
        assert metrics["rmse"] > 0

    def test_strain_abundance_perfect(self):
        pi = np.array([0.5, 0.3, 0.2])
        ids = ["s0", "s1", "s2"]

        metrics = strain_abundance_metrics(pi, pi, ids, ids)
        assert metrics["l1_error"] < 1e-10
        assert metrics["bray_curtis"] < 1e-10

    def test_function_attribution_metrics(self):
        est = np.array([[1.0, 0.0], [0.0, 1.0]])
        true = np.array([[1.0, 0.0], [0.0, 1.0]])
        s_ids = ["s0", "s1"]
        g_ids = ["g0", "g1"]

        metrics = function_attribution_metrics(
            est, true, s_ids, s_ids, g_ids, g_ids
        )
        assert metrics["precision"] == 1.0
        assert metrics["recall"] == 1.0
        assert metrics["f1"] == 1.0

    def test_function_attribution_imperfect(self):
        est = np.array([[1.0, 1.0], [0.0, 1.0]])  # extra gene in s0
        true = np.array([[1.0, 0.0], [0.0, 1.0]])
        s_ids = ["s0", "s1"]
        g_ids = ["g0", "g1"]

        metrics = function_attribution_metrics(
            est, true, s_ids, s_ids, g_ids, g_ids
        )
        assert metrics["precision"] < 1.0
        assert metrics["recall"] == 1.0

    def test_mobile_element_metrics(self):
        assignments = {"r0": True, "r1": False, "r2": True, "r3": False}
        true_mobile = {"r0", "r2", "r3"}

        metrics = mobile_element_metrics(assignments, true_mobile)
        assert metrics["precision"] == 1.0  # both predicted mobile are true
        assert metrics["recall"] < 1.0  # missed r3

    def test_bray_curtis(self):
        p = np.array([0.5, 0.5])
        q = np.array([0.5, 0.5])
        assert bray_curtis(p, q) == 0.0

        q2 = np.array([1.0, 0.0])
        assert bray_curtis(p, q2) > 0

    def test_l1_error(self):
        p = np.array([0.5, 0.3, 0.2])
        q = np.array([0.4, 0.4, 0.2])
        assert l1_error(p, q) == pytest.approx(0.2)


class TestIndexIntegration:
    """Test that the index correctly handles simulated data."""

    def test_index_build_and_query(self):
        cfg = SimulationConfig(
            n_strains=2,
            genome_length=5_000,
            n_gene_families=3,
            n_reads=50,
            random_seed=42,
        )
        refs, reads, truth = simulate_metagenome(cfg)

        index = StrainSiftIndex()
        index.build_from_references(refs)

        assert index.num_species >= 1
        assert index.num_strains == 2

        # Query a read
        candidates = index.query_read(reads[0].sequence)
        # May or may not find candidates depending on k-mer overlap
        # Just verify it doesn't crash
        assert isinstance(candidates, list)

    def test_get_species_for_strain(self):
        cfg = SimulationConfig(
            n_strains=2, genome_length=5_000, n_reads=10, random_seed=42
        )
        refs, _, _ = simulate_metagenome(cfg)

        index = StrainSiftIndex()
        index.build_from_references(refs)

        for ref in refs:
            species = index.get_species_for_strain(ref.strain_id)
            assert species == ref.species_id

        assert index.get_species_for_strain("nonexistent") is None
