"""Tests for mag.func_profile — functional profiling analysis."""

from __future__ import annotations

import numpy as np
import pytest

from mag.func_io import FunctionalTable
from mag.func_profile import (
    cazyme_summary,
    differential_pathway,
    functional_redundancy,
    pathway_abundance,
)
from mag.io import AbundanceTable, SampleMetadata
from mag.tests.fixtures import generate_synthetic_abundance_table


# ---------------------------------------------------------------------------
# Helpers: build a small FunctionalTable compatible with synthetic fixtures
# ---------------------------------------------------------------------------


def _make_func_table(mag_ids: list[str], n_functions: int = 5) -> FunctionalTable:
    """Create a binary function table with deterministic pattern."""
    rng = np.random.default_rng(99)
    values = (rng.random((n_functions, len(mag_ids))) > 0.5).astype(np.float64)
    function_ids = [f"KO_{i:04d}" for i in range(n_functions)]
    return FunctionalTable(
        function_ids=function_ids,
        mag_ids=list(mag_ids),
        values=values,
        function_type="ko",
    )


# ===== TestPathwayAbundance =================================================


class TestPathwayAbundance:
    def test_output_shape(self) -> None:
        """Result should have n_functions rows and n_samples columns."""
        abundance, metadata, _ = generate_synthetic_abundance_table()
        ft = _make_func_table(abundance.mag_ids, n_functions=5)
        result = pathway_abundance(ft, abundance)

        assert result.n_mags == 5  # n_functions
        assert result.n_samples == abundance.n_samples
        assert result.mag_ids == ft.function_ids
        assert result.sample_ids == abundance.sample_ids

    def test_zero_function_zero_abundance(self) -> None:
        """If a function is absent from all MAGs, its pathway abundance is 0."""
        abundance, _, _ = generate_synthetic_abundance_table()
        n_functions = 3
        # All-zero function row for the first function
        values = np.zeros((n_functions, len(abundance.mag_ids)), dtype=np.float64)
        values[1, 0] = 1.0  # only the second function is present in MAG_000
        values[2, :] = 1.0  # third function present in all MAGs

        ft = FunctionalTable(
            function_ids=["F_zero", "F_one", "F_all"],
            mag_ids=list(abundance.mag_ids),
            values=values,
            function_type="ko",
        )
        result = pathway_abundance(ft, abundance)

        # First function row should be all zeros
        np.testing.assert_array_equal(result.abundances[0, :], 0.0)
        # Second function row should equal MAG_000's abundances
        mag0_idx = abundance.mag_ids.index("MAG_000")
        np.testing.assert_array_almost_equal(
            result.abundances[1, :], abundance.abundances[mag0_idx, :]
        )


# ===== TestDifferentialPathway ==============================================


class TestDifferentialPathway:
    def test_returns_result(self) -> None:
        """Result should have correct shape and q-values in [0, 1]."""
        abundance, metadata, _ = generate_synthetic_abundance_table()
        ft = _make_func_table(abundance.mag_ids, n_functions=8)

        result = differential_pathway(
            ft,
            abundance,
            metadata,
            grouping_var="compartment",
            group1="topsoil",
            group2="bulk",
        )

        # Number of results should match number of functions
        assert len(result.mag_ids) == 8
        assert len(result.log_fold_changes) == 8
        assert len(result.p_values) == 8
        assert len(result.q_values) == 8
        assert len(result.effect_sizes) == 8

        # q-values must be in [0, 1]
        assert np.all(result.q_values >= 0.0)
        assert np.all(result.q_values <= 1.0)

        # p-values must be in [0, 1]
        assert np.all(result.p_values >= 0.0)
        assert np.all(result.p_values <= 1.0)


# ===== TestFunctionalRedundancy =============================================


class TestFunctionalRedundancy:
    def test_returns_per_group(self) -> None:
        """All groups should be present, with valid n_carriers and shannon."""
        abundance, metadata, _ = generate_synthetic_abundance_table()
        ft = _make_func_table(abundance.mag_ids, n_functions=4)

        result = functional_redundancy(
            ft, abundance, metadata, grouping_var="compartment"
        )

        # All three groups from the fixture should be present
        expected_groups = {"topsoil", "bulk", "rhizosphere"}
        assert set(result.keys()) == expected_groups

        for group_name, func_dict in result.items():
            assert len(func_dict) == 4  # n_functions
            for func_id, metrics in func_dict.items():
                assert "n_carriers" in metrics
                assert "shannon" in metrics
                assert metrics["n_carriers"] >= 0
                assert metrics["shannon"] >= 0.0


# ===== TestCAZymeSummary ====================================================


class TestCAZymeSummary:
    def test_aggregates_by_class(self) -> None:
        """CAZyme families should be aggregated by their class prefix."""
        # Build a functional table with known CAZy families
        function_ids = ["GH5", "GH13", "GT2", "AA3", "PL1", "CE1"]
        mag_ids = ["M1", "M2", "M3"]
        # M1 has GH5(2), GH13(1), GT2(3)
        # M2 has GH13(1), AA3(2), PL1(1)
        # M3 has GH5(1), CE1(4)
        values = np.array(
            [
                [2.0, 0.0, 1.0],  # GH5
                [1.0, 1.0, 0.0],  # GH13
                [3.0, 0.0, 0.0],  # GT2
                [0.0, 2.0, 0.0],  # AA3
                [0.0, 1.0, 0.0],  # PL1
                [0.0, 0.0, 4.0],  # CE1
            ],
            dtype=np.float64,
        )
        ft = FunctionalTable(
            function_ids=function_ids,
            mag_ids=mag_ids,
            values=values,
            function_type="cazy",
        )

        result = cazyme_summary(ft)

        # Should have entries for GH, GT, AA, PL, CE
        assert set(result.keys()) == {"GH", "GT", "AA", "PL", "CE"}

        # GH: GH5 + GH13 → total = 2+0+1+1+1+0 = 5, 2 families, MAGs M1,M2,M3
        gh = result["GH"]
        assert gh["total_genes"] == 5.0
        assert gh["n_families"] == 2.0
        assert gh["n_mags"] == 3.0  # M1 has both, M2 has GH13, M3 has GH5

        # GT: GT2 → total = 3, 1 family, 1 MAG (M1)
        gt = result["GT"]
        assert gt["total_genes"] == 3.0
        assert gt["n_families"] == 1.0
        assert gt["n_mags"] == 1.0

        # AA: AA3 → total = 2, 1 family, 1 MAG (M2)
        aa = result["AA"]
        assert aa["total_genes"] == 2.0
        assert aa["n_families"] == 1.0
        assert aa["n_mags"] == 1.0

        # PL: PL1 → total = 1, 1 family, 1 MAG (M2)
        pl = result["PL"]
        assert pl["total_genes"] == 1.0
        assert pl["n_families"] == 1.0
        assert pl["n_mags"] == 1.0

        # CE: CE1 → total = 4, 1 family, 1 MAG (M3)
        ce = result["CE"]
        assert ce["total_genes"] == 4.0
        assert ce["n_families"] == 1.0
        assert ce["n_mags"] == 1.0
