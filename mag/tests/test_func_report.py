# mag/tests/test_func_report.py
"""Tests for functional report generation."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from mag.func_io import DRAMAnnotation, FunctionalTable, build_functional_table
from mag.func_report import generate_func_report
from mag.tests.fixtures import generate_synthetic_abundance_table


@pytest.fixture
def synth():
    return generate_synthetic_abundance_table()


@pytest.fixture
def annotations(synth) -> list[DRAMAnnotation]:
    """Generate synthetic DRAM annotations for test MAGs."""
    table = synth[0]
    rng = np.random.default_rng(77)
    annots: list[DRAMAnnotation] = []
    ko_pool = ["K02588", "K02586", "K02591", "K06137", "K00111", "K01505", "K99999"]
    cazy_pool = ["GH5", "GH13", "GT2", "AA3", "PL1", "CE1", None, None]

    for mag_id in table.mag_ids:
        n_genes = rng.integers(5, 15)
        for g in range(n_genes):
            annots.append(
                DRAMAnnotation(
                    gene_id=f"{mag_id}_{g:05d}",
                    mag_id=mag_id,
                    ko_id=rng.choice(ko_pool) if rng.random() > 0.3 else None,
                    cazy_id=rng.choice(cazy_pool),
                )
            )
    return annots


class TestGenerateFuncReport:
    def test_creates_output_files(self, synth, annotations, tmp_path) -> None:
        table, meta, tax = synth
        generate_func_report(
            abundance=table,
            annotations=annotations,
            taxonomy=tax,
            metadata=meta,
            grouping_var="compartment",
            output_dir=tmp_path,
        )
        assert (tmp_path / "pathway_abundance.csv").exists()
        assert (tmp_path / "pgpr_traits.csv").exists()
        assert (tmp_path / "cazyme_summary.csv").exists()
        assert (tmp_path / "functional_redundancy.csv").exists()
        assert (tmp_path / "func_report.html").exists()

    def test_differential_csvs_created(self, synth, annotations, tmp_path) -> None:
        table, meta, tax = synth
        generate_func_report(
            abundance=table,
            annotations=annotations,
            taxonomy=tax,
            metadata=meta,
            grouping_var="compartment",
            output_dir=tmp_path,
        )
        # Should have pairwise pathway differential CSVs
        diff_files = list(tmp_path.glob("pathway_differential_*.csv"))
        assert len(diff_files) >= 1

    def test_plots_created(self, synth, annotations, tmp_path) -> None:
        table, meta, tax = synth
        generate_func_report(
            abundance=table,
            annotations=annotations,
            taxonomy=tax,
            metadata=meta,
            grouping_var="compartment",
            output_dir=tmp_path,
        )
        assert (tmp_path / "pgpr_traits.pdf").exists()
        assert (tmp_path / "cazyme_bars.pdf").exists()

    def test_no_taxonomy_ok(self, synth, annotations, tmp_path) -> None:
        table, meta, _ = synth
        generate_func_report(
            abundance=table,
            annotations=annotations,
            taxonomy=None,
            metadata=meta,
            grouping_var="compartment",
            output_dir=tmp_path,
        )
        assert (tmp_path / "func_report.html").exists()
