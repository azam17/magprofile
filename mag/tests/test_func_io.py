"""Tests for mag.func_io — DRAM parsing and FunctionalTable."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from mag.func_io import (
    PGPR_MARKERS,
    DRAMAnnotation,
    FunctionalTable,
    build_functional_table,
    load_dram_annotations,
)
from mag.io import AbundanceTable

# ---------------------------------------------------------------------------
# Shared fixture: minimal DRAM annotations.tsv with 5 genes across 2 MAGs
# ---------------------------------------------------------------------------

_HEADER = "\t".join(
    ["", "gene_id", "kegg_id", "kegg_hit", "cazy_hits", "pfam_hits"]
)

# MAG_001 has 3 genes, MAG_002 has 2 genes.
# Mix of KO IDs, CAZy IDs, and one PGPR-relevant KO (nifH → K02588).
_ROWS = [
    "MAG_001_00001\tMAG_001_00001\tK02588\tnifH; nitrogenase\t\tPF00142",
    "MAG_001_00002\tMAG_001_00002\tK00001\talcohol dehydrogenase\tGH13\t",
    "MAG_001_00003\tMAG_001_00003\t\t\tGH5\tPF00150",
    "MAG_002_00001\tMAG_002_00001\tK00002\tsome enzyme\tGH13\tPF00560",
    "MAG_002_00002\tMAG_002_00002\tK02586\tnifD; nitrogenase\t\t",
]


@pytest.fixture
def dram_tsv(tmp_path: Path) -> Path:
    """Write a minimal DRAM annotations.tsv and return its path."""
    tsv = tmp_path / "annotations.tsv"
    tsv.write_text(_HEADER + "\n" + "\n".join(_ROWS) + "\n")
    return tsv


@pytest.fixture
def annotations(dram_tsv: Path) -> list[DRAMAnnotation]:
    return load_dram_annotations(dram_tsv)


# ===== TestLoadDRAMAnnotations ==============================================


class TestLoadDRAMAnnotations:
    def test_parse_basic(self, annotations: list[DRAMAnnotation]) -> None:
        """Five annotation rows should produce 5 DRAMAnnotation objects."""
        assert len(annotations) == 5

    def test_cazy_parsed(self, annotations: list[DRAMAnnotation]) -> None:
        """Three of the five genes have CAZy annotations."""
        cazy_count = sum(1 for a in annotations if a.cazy_id is not None)
        assert cazy_count == 3

    def test_mag_id_from_gene_id(self, annotations: list[DRAMAnnotation]) -> None:
        """MAG IDs should be derived by stripping the trailing _NNNNN suffix."""
        mag_ids = sorted(set(a.mag_id for a in annotations))
        assert mag_ids == ["MAG_001", "MAG_002"]


# ===== TestBuildFunctionalTable =============================================


class TestBuildFunctionalTable:
    def test_ko_table(self, annotations: list[DRAMAnnotation]) -> None:
        ft = build_functional_table(annotations, "ko")
        assert ft.function_type == "ko"
        # 4 unique KOs in the fixture: K02588, K00001, K00002, K02586
        assert len(ft.function_ids) == 4
        assert len(ft.mag_ids) == 2
        # All values are binary (0 or 1)
        assert set(np.unique(ft.values)).issubset({0.0, 1.0})

    def test_cazy_table(self, annotations: list[DRAMAnnotation]) -> None:
        ft = build_functional_table(annotations, "cazy")
        assert ft.function_type == "cazy"
        # 2 unique CAZy families: GH13 (count 2), GH5 (count 1)
        assert len(ft.function_ids) == 2
        # GH13 appears in MAG_001 and MAG_002 → counts should reflect that
        gh13_idx = ft.function_ids.index("GH13")
        assert ft.values[gh13_idx].sum() == 2.0

    def test_pgpr_table(self, annotations: list[DRAMAnnotation]) -> None:
        ft = build_functional_table(annotations, "pgpr")
        assert ft.function_type == "pgpr"
        # Rows = all 13 PGPR traits; columns = 2 MAGs
        assert len(ft.function_ids) == len(PGPR_MARKERS)
        assert len(ft.mag_ids) == 2
        # nifH (K02588) present in MAG_001, nifD (K02586) present in MAG_002
        nifh_row = ft.function_ids.index("nifH")
        nifd_row = ft.function_ids.index("nifD")
        mag1_col = ft.mag_ids.index("MAG_001")
        mag2_col = ft.mag_ids.index("MAG_002")
        assert ft.values[nifh_row, mag1_col] == 1.0
        assert ft.values[nifd_row, mag2_col] == 1.0
        # Other traits should be 0
        assert ft.values[ft.function_ids.index("pqqC")].sum() == 0.0

    def test_unknown_type_raises(self, annotations: list[DRAMAnnotation]) -> None:
        with pytest.raises(ValueError, match="Unknown function_type"):
            build_functional_table(annotations, "nonexistent")


# ===== TestFunctionalTableToSampleAbundance =================================


class TestFunctionalTableToSampleAbundance:
    def test_matrix_multiplication(self) -> None:
        """F @ A should produce correct function-by-sample matrix."""
        ft = FunctionalTable(
            function_ids=["KO1", "KO2"],
            mag_ids=["M1", "M2"],
            values=np.array([[1.0, 0.0], [0.0, 1.0]]),
            function_type="ko",
        )
        abund = AbundanceTable(
            mag_ids=["M1", "M2"],
            sample_ids=["S1", "S2"],
            abundances=np.array([[10.0, 20.0], [30.0, 40.0]]),
        )
        result = ft.to_sample_abundance(abund)
        expected = np.array([[10.0, 20.0], [30.0, 40.0]])
        np.testing.assert_array_almost_equal(result.abundances, expected)
        assert result.mag_ids == ["KO1", "KO2"]
        assert result.sample_ids == ["S1", "S2"]

    def test_shared_mag_subset(self) -> None:
        """Only MAGs present in both tables should be used."""
        ft = FunctionalTable(
            function_ids=["KO1"],
            mag_ids=["M1", "M3"],  # M3 not in abundance table
            values=np.array([[1.0, 1.0]]),
            function_type="ko",
        )
        abund = AbundanceTable(
            mag_ids=["M1", "M2"],  # M2 not in functional table
            sample_ids=["S1"],
            abundances=np.array([[5.0], [99.0]]),
        )
        result = ft.to_sample_abundance(abund)
        # Only M1 overlaps: 1.0 * 5.0 = 5.0
        np.testing.assert_array_almost_equal(result.abundances, np.array([[5.0]]))


# ===== TestPGPRMarkers ======================================================


class TestPGPRMarkers:
    def test_markers_exist(self) -> None:
        """All 13 PGPR markers should be defined."""
        assert len(PGPR_MARKERS) == 13
        assert "nifH" in PGPR_MARKERS
        assert PGPR_MARKERS["nifH"] == "K02588"
        assert "phzF" in PGPR_MARKERS
        assert PGPR_MARKERS["phzF"] == "K18000"
