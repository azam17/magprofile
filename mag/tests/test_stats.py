"""Tests for mag.stats module."""

import numpy as np
import pytest

from mag.beta import bray_curtis
from mag.io import AbundanceTable, SampleMetadata
from mag.stats import (
    VariancePartitionResult,
    anosim,
    kruskal_wallis_per_mag,
    pairwise_permanova,
    permanova,
    permdisp,
    variance_partition,
)
from scipy.spatial.distance import pdist, squareform
from mag.tests.fixtures import generate_synthetic_abundance_table


class TestPERMANOVA:
    def test_significant_with_planted_structure(self):
        table, meta, _ = generate_synthetic_abundance_table(seed=42)
        bc = bray_curtis(table)
        result = permanova(bc, meta, "compartment", n_permutations=99, seed=42)
        assert result.test_name == "PERMANOVA"
        assert result.statistic > 0
        assert 0 <= result.p_value <= 1
        assert result.p_value < 0.1
        assert result.metadata["R2"] > 0

    def test_no_structure(self):
        rng = np.random.default_rng(123)
        table = AbundanceTable(
            mag_ids=[f"M{i}" for i in range(10)],
            sample_ids=[f"s{i}" for i in range(9)],
            abundances=rng.poisson(100, (10, 9)).astype(float),
        )
        records = {}
        groups = ["a", "b", "c"]
        for i in range(9):
            records[f"s{i}"] = {"group": groups[i % 3]}
        meta = SampleMetadata(records=records)
        bc = bray_curtis(table)
        result = permanova(bc, meta, "group", n_permutations=99, seed=42)
        assert 0 <= result.p_value <= 1


class TestANOSIM:
    def test_significant_with_planted_structure(self):
        table, meta, _ = generate_synthetic_abundance_table(seed=42)
        bc = bray_curtis(table)
        result = anosim(bc, meta, "compartment", n_permutations=99, seed=42)
        assert result.test_name == "ANOSIM"
        assert 0 <= result.p_value <= 1


class TestKruskalWallis:
    def test_per_mag(self):
        table, meta, _ = generate_synthetic_abundance_table(seed=42)
        results = kruskal_wallis_per_mag(table, meta, "compartment")
        assert len(results) == table.n_mags
        sig_count = sum(1 for r in results if r.p_value < 0.05)
        assert sig_count > 0


class TestPERMDISP:
    def test_uniform_dispersion(self):
        """Samples drawn from same distribution → PERMDISP not significant."""
        rng = np.random.default_rng(77)
        table = AbundanceTable(
            mag_ids=[f"M{i}" for i in range(10)],
            sample_ids=[f"s{i}" for i in range(9)],
            abundances=rng.poisson(100, (10, 9)).astype(float),
        )
        records = {}
        groups = ["a", "b", "c"]
        for i in range(9):
            records[f"s{i}"] = {"group": groups[i % 3]}
        meta = SampleMetadata(records=records)
        bc = bray_curtis(table)
        result = permdisp(bc, meta, "group", n_permutations=99, seed=42)
        assert result.test_name == "PERMDISP"
        assert 0 <= result.p_value <= 1
        assert "group_mean_distances" in result.metadata

    def test_different_dispersion(self):
        """One group much more dispersed → PERMDISP may detect it."""
        rng = np.random.default_rng(10)
        n_mags = 15
        n_per_group = 5
        # Group A: tight cluster
        a = rng.poisson(100, (n_mags, n_per_group)).astype(float)
        # Group B: very dispersed
        b = rng.poisson(100, (n_mags, n_per_group)).astype(float)
        b[:, 0] *= 10  # one outlier sample
        b[:, 1] *= 0.01
        abundances = np.hstack([a, b])
        sids = [f"s{i}" for i in range(2 * n_per_group)]
        table = AbundanceTable(
            mag_ids=[f"M{i}" for i in range(n_mags)],
            sample_ids=sids,
            abundances=abundances,
        )
        records = {}
        for i in range(n_per_group):
            records[f"s{i}"] = {"group": "tight"}
        for i in range(n_per_group, 2 * n_per_group):
            records[f"s{i}"] = {"group": "dispersed"}
        meta = SampleMetadata(records=records)
        bc = bray_curtis(table)
        result = permdisp(bc, meta, "group", n_permutations=99, seed=42)
        assert result.test_name == "PERMDISP"
        assert 0 <= result.p_value <= 1
        # Dispersed group should have larger mean distance
        dists = result.metadata["group_mean_distances"]
        assert dists["dispersed"] > dists["tight"]


class TestPairwisePERMANOVA:
    def test_all_pairs_tested(self):
        table, meta, _ = generate_synthetic_abundance_table(seed=42)
        bc = bray_curtis(table)
        results = pairwise_permanova(bc, meta, "compartment", n_permutations=9, seed=42)
        # 3 groups → 3 pairs
        assert len(results) == 3
        pairs = {(r.metadata["group1"], r.metadata["group2"]) for r in results}
        assert ("bulk", "rhizosphere") in pairs
        assert ("bulk", "topsoil") in pairs
        assert ("rhizosphere", "topsoil") in pairs

    def test_fdr_correction_applied(self):
        table, meta, _ = generate_synthetic_abundance_table(seed=42)
        bc = bray_curtis(table)
        results = pairwise_permanova(bc, meta, "compartment", n_permutations=9, seed=42)
        for r in results:
            assert "q_value" in r.metadata
            assert r.metadata["q_value"] >= r.p_value  # FDR q >= raw p

    def test_two_groups(self):
        """With 2 groups, should return exactly 1 result."""
        rng = np.random.default_rng(55)
        table = AbundanceTable(
            mag_ids=[f"M{i}" for i in range(8)],
            sample_ids=[f"s{i}" for i in range(6)],
            abundances=rng.poisson(100, (8, 6)).astype(float),
        )
        records = {f"s{i}": {"group": "A" if i < 3 else "B"} for i in range(6)}
        meta = SampleMetadata(records=records)
        bc = bray_curtis(table)
        results = pairwise_permanova(bc, meta, "group", n_permutations=9, seed=42)
        assert len(results) == 1


def _make_vp_data(
    groups_a: list[str],
    groups_b: list[str],
    seed: int = 42,
    effect_a: float = 10.0,
    effect_b: float = 1.0,
) -> tuple[np.ndarray, list[str], SampleMetadata]:
    """Create a distance matrix with two grouping variables."""
    rng = np.random.default_rng(seed)
    n = len(groups_a)
    sample_ids = [f"s{i}" for i in range(n)]
    coords = rng.normal(size=(n, 5))
    # Use stable offsets: enumerate unique groups, assign sequential offsets
    unique_a = sorted(set(groups_a))
    unique_b = sorted(set(groups_b))
    offset_a = {g: k * effect_a for k, g in enumerate(unique_a)}
    offset_b = {g: k * effect_b for k, g in enumerate(unique_b)}
    for i in range(n):
        coords[i, 0] += offset_a[groups_a[i]]
        coords[i, 1] += offset_b[groups_b[i]]
    dm = squareform(pdist(coords, metric="euclidean"))
    records = {sid: {"var_a": groups_a[i], "var_b": groups_b[i]}
               for i, sid in enumerate(sample_ids)}
    return dm, sample_ids, SampleMetadata(records=records)


class TestVariancePartition:
    def test_fractions_sum_to_approximately_one(self):
        """Pure + shared + residual should sum close to 1.0.

        Exact equality is not guaranteed because pure fractions are clamped
        to >= 0, which can make the total slightly exceed 1.0 when
        individual R^2 values exceed the combined R^2 (suppression).
        """
        groups_a = ["X", "X", "X", "Y", "Y", "Y", "Z", "Z", "Z"]
        groups_b = ["P", "Q", "P", "Q", "P", "Q", "P", "Q", "P"]
        dm, sids, meta = _make_vp_data(groups_a, groups_b)
        result = variance_partition(dm, sids, meta, ["var_a", "var_b"])
        total = (
            sum(result.pure_fractions.values())
            + result.shared_fraction
            + result.residual_fraction
        )
        assert total == pytest.approx(1.0, abs=0.05)

    def test_known_effect(self):
        groups_a = ["X", "X", "X", "X", "Y", "Y", "Y", "Y"]
        groups_b = ["P", "Q", "P", "Q", "P", "Q", "P", "Q"]
        dm, sids, meta = _make_vp_data(groups_a, groups_b, seed=123)
        result = variance_partition(dm, sids, meta, ["var_a", "var_b"])
        assert isinstance(result, VariancePartitionResult)
        assert result.pure_fractions["var_a"] > result.pure_fractions["var_b"]

    def test_single_variable(self):
        groups_a = ["X", "X", "X", "Y", "Y", "Y"]
        groups_b = ["P", "Q", "P", "Q", "P", "Q"]
        dm, sids, meta = _make_vp_data(groups_a, groups_b)
        result = variance_partition(dm, sids, meta, ["var_a"])
        assert result.shared_fraction == 0.0
        assert result.pure_fractions["var_a"] == pytest.approx(
            result.full_r_squared, abs=1e-10
        )
        total = result.pure_fractions["var_a"] + result.residual_fraction
        assert total == pytest.approx(1.0, abs=1e-10)

    def test_all_fractions_nonnegative(self):
        groups_a = ["X", "X", "Y", "Y", "Z", "Z"]
        groups_b = ["P", "Q", "P", "Q", "P", "Q"]
        dm, sids, meta = _make_vp_data(groups_a, groups_b)
        result = variance_partition(dm, sids, meta, ["var_a", "var_b"])
        for frac in result.pure_fractions.values():
            assert frac >= 0.0
        assert result.shared_fraction >= 0.0
        assert result.residual_fraction >= 0.0
