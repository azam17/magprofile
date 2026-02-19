"""Tests for mag.func_profile — functional profiling analysis."""

from __future__ import annotations

import numpy as np
import pytest

from mag.func_io import FunctionalTable
from mag.func_profile import (
    cazyme_summary,
    differential_pathway,
    functional_redundancy,
    keystone_pgpr_cross_reference,
    pathway_abundance,
    pgpr_enrichment,
    redundancy_comparison,
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

    def test_exact_carrier_count(self) -> None:
        """n_carriers must equal the number of MAGs with value > 0 for that function."""
        # 3 MAGs, 2 functions
        ft = FunctionalTable(
            function_ids=["F1", "F2"],
            mag_ids=["A", "B", "C"],
            values=np.array([
                [1.0, 0.0, 1.0],  # F1: carried by A, C
                [0.0, 0.0, 0.0],  # F2: carried by nobody
            ]),
            function_type="ko",
        )
        abund = AbundanceTable(
            mag_ids=["A", "B", "C"],
            sample_ids=["s1", "s2"],
            abundances=np.array([
                [10.0, 20.0],
                [30.0, 40.0],
                [50.0, 60.0],
            ]),
        )
        meta = SampleMetadata(records={
            "s1": {"group": "g1"},
            "s2": {"group": "g1"},
        })

        result = functional_redundancy(ft, abund, meta, grouping_var="group")
        assert result["g1"]["F1"]["n_carriers"] == 2.0
        assert result["g1"]["F2"]["n_carriers"] == 0.0
        assert result["g1"]["F2"]["shannon"] == 0.0

    def test_exact_shannon(self) -> None:
        """Shannon diversity must match hand computation."""
        # 3 MAGs carry F1; mean abundances across the group:
        # A: mean([10,20])=15, B: mean([30,40])=35, C: mean([50,60])=55
        # total = 105, props = [15/105, 35/105, 55/105]
        # shannon = -sum(p * ln(p))
        ft = FunctionalTable(
            function_ids=["F1"],
            mag_ids=["A", "B", "C"],
            values=np.array([[1.0, 1.0, 1.0]]),
            function_type="ko",
        )
        abund = AbundanceTable(
            mag_ids=["A", "B", "C"],
            sample_ids=["s1", "s2"],
            abundances=np.array([
                [10.0, 20.0],
                [30.0, 40.0],
                [50.0, 60.0],
            ]),
        )
        meta = SampleMetadata(records={
            "s1": {"group": "g1"},
            "s2": {"group": "g1"},
        })

        result = functional_redundancy(ft, abund, meta, grouping_var="group")

        # Hand-compute Shannon
        means = np.array([15.0, 35.0, 55.0])
        props = means / means.sum()
        expected_shannon = float(-np.sum(props * np.log(props)))

        np.testing.assert_allclose(
            result["g1"]["F1"]["shannon"], expected_shannon, rtol=1e-10,
        )

    def test_single_carrier_shannon_zero(self) -> None:
        """Shannon must be 0.0 when only one MAG carries the function."""
        ft = FunctionalTable(
            function_ids=["F1"],
            mag_ids=["A", "B"],
            values=np.array([[1.0, 0.0]]),
            function_type="ko",
        )
        abund = AbundanceTable(
            mag_ids=["A", "B"],
            sample_ids=["s1"],
            abundances=np.array([[100.0], [200.0]]),
        )
        meta = SampleMetadata(records={"s1": {"group": "g1"}})

        result = functional_redundancy(ft, abund, meta, grouping_var="group")
        assert result["g1"]["F1"]["n_carriers"] == 1.0
        assert result["g1"]["F1"]["shannon"] == 0.0


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


# ===== TestPGPREnrichment ====================================================


def _make_pgpr_table() -> FunctionalTable:
    """PGPR table with 3 traits, 6 MAGs."""
    return FunctionalTable(
        function_ids=["nifH", "pqqC", "gcd"],
        mag_ids=["M1", "M2", "M3", "M4", "M5", "M6"],
        values=np.array([
            [1, 1, 0, 0, 0, 0],  # nifH: M1, M2
            [0, 0, 1, 1, 0, 0],  # pqqC: M3, M4
            [1, 0, 0, 0, 1, 1],  # gcd:  M1, M5, M6
        ], dtype=np.float64),
        function_type="pgpr",
    )


class TestPGPREnrichment:
    def test_returns_result(self) -> None:
        pgpr = _make_pgpr_table()
        # M1-M3 dominant in group A, M4-M6 dominant in group B
        abund = AbundanceTable(
            mag_ids=["M1", "M2", "M3", "M4", "M5", "M6"],
            sample_ids=["sA1", "sA2", "sB1", "sB2"],
            abundances=np.array([
                [100, 90, 1, 2],    # M1: A-dominant
                [80, 70, 5, 3],     # M2: A-dominant
                [60, 50, 10, 8],    # M3: A-dominant
                [2, 3, 100, 90],    # M4: B-dominant
                [1, 2, 80, 70],     # M5: B-dominant
                [5, 3, 60, 50],     # M6: B-dominant
            ], dtype=np.float64),
        )
        meta = SampleMetadata(records={
            "sA1": {"comp": "A"}, "sA2": {"comp": "A"},
            "sB1": {"comp": "B"}, "sB2": {"comp": "B"},
        })

        result = pgpr_enrichment(pgpr, abund, meta, "comp")
        assert set(result.trait_names) == {"nifH", "pqqC", "gcd"}
        assert set(result.group_names) == {"A", "B"}
        for trait in result.trait_names:
            assert 0 <= result.p_values[trait] <= 1
            assert 0 <= result.q_values[trait] <= 1

    def test_counts_correct(self) -> None:
        pgpr = _make_pgpr_table()
        # All MAGs dominant in one group
        abund = AbundanceTable(
            mag_ids=["M1", "M2", "M3", "M4", "M5", "M6"],
            sample_ids=["s1", "s2"],
            abundances=np.ones((6, 2), dtype=np.float64) * 10,
        )
        meta = SampleMetadata(records={
            "s1": {"g": "X"}, "s2": {"g": "X"},
        })
        result = pgpr_enrichment(pgpr, abund, meta, "g")
        # All MAGs assigned to group X
        assert result.counts["nifH"]["X"] == 2  # M1, M2
        assert result.counts["pqqC"]["X"] == 2  # M3, M4


# ===== TestKeystonePGPRCrossReference ========================================


class TestKeystonePGPRCrossReference:
    def test_cross_reference(self) -> None:
        pgpr = _make_pgpr_table()
        # M1 and M3 are keystones
        mag_ids = ["M1", "M2", "M3", "M4", "M5", "M6"]
        mask = np.array([True, False, True, False, False, False])

        result = keystone_pgpr_cross_reference(mag_ids, mask, pgpr)
        assert result.n_keystones_total == 2
        # M1 carries nifH + gcd, M3 carries pqqC
        assert result.n_keystones_with_pgpr == 2
        assert "nifH" in result.keystone_traits["M1"]
        assert "gcd" in result.keystone_traits["M1"]
        assert "pqqC" in result.keystone_traits["M3"]
        assert 0 <= result.enrichment_p_value <= 1

    def test_no_overlap(self) -> None:
        """Keystones carry no PGPR traits."""
        pgpr = FunctionalTable(
            function_ids=["nifH"],
            mag_ids=["A", "B", "C"],
            values=np.array([[0, 1, 0]], dtype=np.float64),
            function_type="pgpr",
        )
        # A and C are keystones, but only B has nifH
        result = keystone_pgpr_cross_reference(
            ["A", "B", "C"],
            np.array([True, False, True]),
            pgpr,
        )
        assert result.n_keystones_with_pgpr == 0
        assert result.n_non_keystones_with_pgpr == 1


# ===== TestRedundancyComparison ==============================================


class TestRedundancyComparison:
    def test_returns_result(self) -> None:
        redundancy = {
            "bulk": {"F1": {"n_carriers": 5, "shannon": 1.2},
                     "F2": {"n_carriers": 3, "shannon": 0.8}},
            "endo": {"F1": {"n_carriers": 2, "shannon": 0.5},
                     "F2": {"n_carriers": 1, "shannon": 0.0}},
        }
        result = redundancy_comparison(redundancy)
        assert set(result.group_names) == {"bulk", "endo"}
        assert result.mean_shannon["bulk"] > result.mean_shannon["endo"]
        assert 0 <= result.p_value <= 1
        assert ("bulk", "endo") in result.pairwise

    def test_single_group(self) -> None:
        redundancy = {
            "only": {"F1": {"n_carriers": 3, "shannon": 1.0}},
        }
        result = redundancy_comparison(redundancy)
        assert result.p_value == 1.0
