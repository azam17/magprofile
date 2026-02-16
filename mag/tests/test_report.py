"""Tests for mag.report module (integration test)."""

import pytest

from mag.report import generate_report
from mag.tests.fixtures import generate_synthetic_abundance_table


class TestReport:
    def test_full_pipeline(self, tmp_path):
        table, meta, tax = generate_synthetic_abundance_table(seed=42)
        out = tmp_path / "report_output"
        generate_report(table, tax, meta, "compartment", out, n_permutations=9)

        # Original outputs
        assert (out / "alpha_diversity.csv").exists()
        assert (out / "beta_diversity_bray_curtis.csv").exists()
        assert (out / "beta_diversity_jaccard.csv").exists()
        assert (out / "ordination_coordinates.csv").exists()
        assert (out / "permanova_results.txt").exists()
        assert (out / "anosim_results.txt").exists()
        assert (out / "indicator_species.csv").exists()
        assert (out / "compartment_specificity.csv").exists()
        assert (out / "summary_report.txt").exists()
        assert (out / "alpha_diversity.pdf").exists()
        assert (out / "ordination_pcoa.pdf").exists()
        assert (out / "taxonomy_bars.pdf").exists()
        assert (out / "venn_diagram.pdf").exists()
        assert (out / "indicator_species.pdf").exists()

        # New outputs (v0.2)
        assert (out / "permdisp_results.txt").exists()
        assert (out / "pairwise_permanova.csv").exists()
        assert (out / "core_microbiome.csv").exists()
        assert (out / "report.html").exists()
        assert (out / "taxonomy_bars_grouped.pdf").exists()

        # Verify HTML report is non-trivial
        html = (out / "report.html").read_text()
        assert "MAG Community Profiling Report" in html
        assert "Executive" in html or "executive" in html
        assert "PERMANOVA" in html
        assert "Core Microbiome" in html

    def test_without_taxonomy(self, tmp_path):
        table, meta, _ = generate_synthetic_abundance_table(seed=42)
        out = tmp_path / "report_no_tax"
        generate_report(table, None, meta, "compartment", out, n_permutations=9)
        assert (out / "summary_report.txt").exists()
        assert not (out / "taxonomy_bars.pdf").exists()
        assert not (out / "taxonomy_bars_grouped.pdf").exists()
        # HTML report still generated
        assert (out / "report.html").exists()

    def test_prevalence_filtering(self, tmp_path):
        table, meta, tax = generate_synthetic_abundance_table(seed=42)
        out = tmp_path / "report_filtered"
        generate_report(table, tax, meta, "compartment", out, n_permutations=9, min_prevalence=0.5)
        assert (out / "summary_report.txt").exists()
        summary = (out / "summary_report.txt").read_text()
        assert "Prevalence" in summary or "prevalence" in summary

    def test_no_filtering(self, tmp_path):
        table, meta, tax = generate_synthetic_abundance_table(seed=42)
        out = tmp_path / "report_nofilter"
        generate_report(table, tax, meta, "compartment", out, n_permutations=9, min_prevalence=0.0)
        assert (out / "summary_report.txt").exists()
