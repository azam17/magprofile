# mag/tests/test_net_report.py
"""Tests for network report generation."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from mag.net_report import generate_net_report
from mag.tests.fixtures import generate_synthetic_abundance_table


@pytest.fixture
def synth():
    return generate_synthetic_abundance_table()


class TestGenerateNetReport:
    def test_creates_output_files(self, synth, tmp_path) -> None:
        table, meta, tax = synth
        generate_net_report(
            abundance=table,
            taxonomy=tax,
            metadata=meta,
            grouping_var="compartment",
            output_dir=tmp_path,
        )
        assert (tmp_path / "network_edges.csv").exists()
        assert (tmp_path / "network_topology.csv").exists()
        assert (tmp_path / "keystone_taxa.csv").exists()
        assert (tmp_path / "net_report.html").exists()

    def test_per_group_networks(self, synth, tmp_path) -> None:
        table, meta, tax = synth
        generate_net_report(
            abundance=table,
            taxonomy=tax,
            metadata=meta,
            grouping_var="compartment",
            output_dir=tmp_path,
        )
        # Should have per-group network edge files
        group_files = list(tmp_path.glob("network_edges_*.csv"))
        assert len(group_files) >= 1

    def test_differential_network_csv(self, synth, tmp_path) -> None:
        table, meta, tax = synth
        generate_net_report(
            abundance=table,
            taxonomy=tax,
            metadata=meta,
            grouping_var="compartment",
            output_dir=tmp_path,
        )
        diff_files = list(tmp_path.glob("differential_network_*.csv"))
        assert len(diff_files) >= 1

    def test_no_taxonomy_ok(self, synth, tmp_path) -> None:
        table, meta, _ = synth
        generate_net_report(
            abundance=table,
            taxonomy=None,
            metadata=meta,
            grouping_var="compartment",
            output_dir=tmp_path,
        )
        assert (tmp_path / "net_report.html").exists()

    def test_plots_created(self, synth, tmp_path) -> None:
        table, meta, tax = synth
        generate_net_report(
            abundance=table,
            taxonomy=tax,
            metadata=meta,
            grouping_var="compartment",
            output_dir=tmp_path,
        )
        assert (tmp_path / "network_graph.pdf").exists()
        assert (tmp_path / "degree_distribution.pdf").exists()
