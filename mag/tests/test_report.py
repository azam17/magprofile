"""Tests for mag.report module (integration test)."""

import csv

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

        # v0.2.1 outputs — phylum composition and rank-level differential
        assert (out / "phylum_composition.csv").exists()
        # 3 groups → 3 pairwise comparisons
        group_names = sorted(meta.get_groups("compartment").keys())
        for i, g1 in enumerate(group_names):
            for g2 in group_names[i + 1 :]:
                assert (out / f"phylum_differential_{g1}_vs_{g2}.csv").exists()

        # Verify HTML report is non-trivial
        html = (out / "report.html").read_text()
        assert "MAG Community Profiling Report" in html
        assert "Executive" in html or "executive" in html
        assert "PERMANOVA" in html
        assert "Core Microbiome" in html
        # v0.2.1: taxonomy content in HTML
        assert "phylum" in html.lower()

    def test_taxonomy_columns_in_csvs(self, tmp_path):
        """Verify taxonomy columns are present in CSVs when taxonomy is provided."""
        table, meta, tax = generate_synthetic_abundance_table(seed=42)
        out = tmp_path / "report_tax_cols"
        generate_report(table, tax, meta, "compartment", out, n_permutations=9)

        # Differential CSV should have phylum/class/genus columns
        group_names = sorted(meta.get_groups("compartment").keys())
        diff_path = out / f"differential_{group_names[0]}_vs_{group_names[1]}.csv"
        with open(diff_path, newline="") as f:
            reader = csv.reader(f)
            header = next(reader)
        assert "phylum" in header
        assert "class" in header
        assert "genus" in header

        # Indicator CSV should have taxonomy columns
        with open(out / "indicator_species.csv", newline="") as f:
            reader = csv.reader(f)
            header = next(reader)
        assert "phylum" in header
        assert "class" in header
        assert "genus" in header

        # Specificity CSV should have taxonomy columns
        with open(out / "compartment_specificity.csv", newline="") as f:
            reader = csv.reader(f)
            header = next(reader)
        assert "phylum" in header
        assert "class" in header
        assert "genus" in header

        # Core microbiome CSV should have core_mag_phyla column
        with open(out / "core_microbiome.csv", newline="") as f:
            reader = csv.reader(f)
            header = next(reader)
        assert "core_mag_phyla" in header

        # Phylum composition CSV header check
        with open(out / "phylum_composition.csv", newline="") as f:
            reader = csv.reader(f)
            header = next(reader)
        assert header[0] == "phylum"
        # Should have mean/sd columns for each group
        for g in group_names:
            assert f"{g}_mean" in header
            assert f"{g}_sd" in header

    def test_taxonomy_absent_when_none(self, tmp_path):
        """Verify taxonomy columns are absent when taxonomy=None."""
        table, meta, _ = generate_synthetic_abundance_table(seed=42)
        out = tmp_path / "report_no_tax"
        generate_report(table, None, meta, "compartment", out, n_permutations=9)

        # Differential CSV should NOT have phylum column
        group_names = sorted(meta.get_groups("compartment").keys())
        diff_path = out / f"differential_{group_names[0]}_vs_{group_names[1]}.csv"
        with open(diff_path, newline="") as f:
            reader = csv.reader(f)
            header = next(reader)
        assert "phylum" not in header

        # Indicator CSV should NOT have phylum column
        with open(out / "indicator_species.csv", newline="") as f:
            reader = csv.reader(f)
            header = next(reader)
        assert "phylum" not in header

        # No phylum composition or rank-level diff files
        assert not (out / "phylum_composition.csv").exists()
        for i, g1 in enumerate(group_names):
            for g2 in group_names[i + 1 :]:
                assert not (out / f"phylum_differential_{g1}_vs_{g2}.csv").exists()

        # HTML report still generated
        assert (out / "report.html").exists()
        assert not (out / "taxonomy_bars.pdf").exists()

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

    def test_html_taxonomy_tables(self, tmp_path):
        """Verify HTML report contains taxonomy columns in tables."""
        table, meta, tax = generate_synthetic_abundance_table(seed=42)
        out = tmp_path / "report_html_tax"
        generate_report(table, tax, meta, "compartment", out, n_permutations=9)

        html = (out / "report.html").read_text()
        # Taxonomy columns should appear in table headers
        assert "Phylum" in html
        assert "Class" in html
        assert "Genus" in html
        # Phylum composition section
        assert "Phylum Composition by Group" in html
        # Phylum-level differential section
        assert "Phylum-Level Differential Abundance" in html
