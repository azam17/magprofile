"""Tests for mag.io module."""

from pathlib import Path

import numpy as np
import pytest

from mag.io import (
    AbundanceTable,
    SampleMetadata,
    TaxonomyRecord,
    TaxonomyTable,
    load_abundance_table,
    load_metadata,
    load_taxonomy,
)


class TestAbundanceTable:
    def test_shape_validation(self):
        with pytest.raises(ValueError):
            AbundanceTable(
                mag_ids=["a", "b"],
                sample_ids=["s1"],
                abundances=np.zeros((3, 1)),
            )

    def test_normalize(self):
        table = AbundanceTable(
            mag_ids=["a", "b"],
            sample_ids=["s1", "s2"],
            abundances=np.array([[100.0, 200.0], [300.0, 800.0]]),
        )
        norm = table.normalize()
        np.testing.assert_allclose(norm.abundances.sum(axis=0), [1.0, 1.0])

    def test_normalize_zero_column(self):
        table = AbundanceTable(
            mag_ids=["a"],
            sample_ids=["s1", "s2"],
            abundances=np.array([[0.0, 100.0]]),
        )
        norm = table.normalize()
        assert norm.abundances[0, 0] == 0.0

    def test_filter_low_abundance(self):
        table = AbundanceTable(
            mag_ids=["a", "b", "c"],
            sample_ids=["s1", "s2"],
            abundances=np.array([[1.0, 0.0], [100.0, 200.0], [0.0, 0.0]]),
        )
        filtered = table.filter_low_abundance(min_total=5.0)
        assert filtered.mag_ids == ["b"]

    def test_filter_prevalence(self):
        table = AbundanceTable(
            mag_ids=["a", "b"],
            sample_ids=["s1", "s2", "s3"],
            abundances=np.array([[1.0, 0.0, 0.0], [1.0, 1.0, 1.0]]),
        )
        filtered = table.filter_prevalence(min_prevalence=0.5)
        assert filtered.mag_ids == ["b"]

    def test_subset_samples(self):
        table = AbundanceTable(
            mag_ids=["a"],
            sample_ids=["s1", "s2", "s3"],
            abundances=np.array([[10.0, 20.0, 30.0]]),
        )
        sub = table.subset_samples(["s1", "s3"])
        assert sub.sample_ids == ["s1", "s3"]
        np.testing.assert_array_equal(sub.abundances, [[10.0, 30.0]])


class TestLoadAbundanceTable:
    def test_load_tsv(self, tmp_path):
        p = tmp_path / "abundance.tsv"
        p.write_text("MAG_ID\ts1\ts2\nMAG_001\t100\t200\nMAG_002\t50\t75\n")
        table = load_abundance_table(p)
        assert table.mag_ids == ["MAG_001", "MAG_002"]
        assert table.sample_ids == ["s1", "s2"]
        assert table.abundances.shape == (2, 2)
        assert table.abundances[0, 0] == 100.0


class TestLoadTaxonomy:
    def test_load_taxonomy(self, tmp_path):
        p = tmp_path / "taxonomy.tsv"
        p.write_text(
            "MAG_ID\tdomain\tphylum\tclass\torder\tfamily\tgenus\tspecies\n"
            "MAG_001\tBacteria\tProteo\tGamma\tO1\tF1\tG1\tsp1\n"
        )
        tax = load_taxonomy(p)
        rec = tax.get("MAG_001")
        assert rec is not None
        assert rec.phylum == "Proteo"
        assert rec.rank("class") == "Gamma"


class TestLoadMetadata:
    def test_load_metadata(self, tmp_path):
        p = tmp_path / "metadata.tsv"
        p.write_text("sample_id\tcompartment\treplicate\ns1\ttopsoil\t1\ns2\tbulk\t1\n")
        meta = load_metadata(p)
        groups = meta.get_groups("compartment")
        assert "topsoil" in groups
        assert groups["topsoil"] == ["s1"]


class TestTaxonomyAggregation:
    def test_aggregate_at_rank(self):
        table = AbundanceTable(
            mag_ids=["a", "b", "c"],
            sample_ids=["s1"],
            abundances=np.array([[10.0], [20.0], [30.0]]),
        )
        tax = TaxonomyTable(
            records={
                "a": TaxonomyRecord(mag_id="a", phylum="P1"),
                "b": TaxonomyRecord(mag_id="b", phylum="P1"),
                "c": TaxonomyRecord(mag_id="c", phylum="P2"),
            }
        )
        agg = tax.aggregate_at_rank(table, "phylum")
        assert agg["P1"][0] == 30.0
        assert agg["P2"][0] == 30.0
