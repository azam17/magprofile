"""Tests for mag.net_topology module."""

import numpy as np
import pytest

from mag.io import AbundanceTable, TaxonomyTable, TaxonomyRecord
from mag.net_correlation import NetworkResult
from mag.net_topology import (
    HubBridgeResult,
    KeystoneTaxaResult,
    ModuleCompositionResult,
    NetworkTopology,
    NullModelResult,
    compute_topology,
    differential_network,
    hub_bridge_classification,
    identify_keystones,
    module_composition,
    network_null_model,
)


@pytest.fixture
def star_network() -> NetworkResult:
    """Star topology: M0 connected to M1,M2,M3,M4 plus M1-M2 edge.

    M0 is the central hub with degree 4.
    M1 and M2 each have degree 2 (connected to M0 and each other).
    M3 and M4 each have degree 1 (connected only to M0).
    """
    mag_ids = ["M0", "M1", "M2", "M3", "M4"]
    n = len(mag_ids)
    adjacency = np.zeros((n, n), dtype=np.float64)
    edges = [
        ("M0", "M1", 0.1),
        ("M0", "M2", 0.2),
        ("M0", "M3", 0.15),
        ("M0", "M4", 0.18),
        ("M1", "M2", 0.25),
    ]
    idx = {m: i for i, m in enumerate(mag_ids)}
    for a, b, _ in edges:
        adjacency[idx[a], idx[b]] = 1.0
        adjacency[idx[b], idx[a]] = 1.0

    return NetworkResult(
        mag_ids=mag_ids,
        adjacency=adjacency,
        edges=edges,
        threshold=0.3,
    )


@pytest.fixture
def star_abundance() -> AbundanceTable:
    """Abundance table where M0 is rare (low abundance) and others are common."""
    return AbundanceTable(
        mag_ids=["M0", "M1", "M2", "M3", "M4"],
        sample_ids=["s1", "s2", "s3"],
        abundances=np.array([
            [1.0, 2.0, 1.0],     # M0: low abundance (mean=1.33)
            [100.0, 120.0, 110.0],  # M1: high
            [90.0, 95.0, 100.0],    # M2: high
            [80.0, 85.0, 90.0],     # M3: high
            [70.0, 75.0, 80.0],     # M4: high
        ]),
    )


@pytest.fixture
def empty_network() -> NetworkResult:
    """Network with nodes but no edges."""
    mag_ids = ["M0", "M1", "M2"]
    n = len(mag_ids)
    return NetworkResult(
        mag_ids=mag_ids,
        adjacency=np.zeros((n, n), dtype=np.float64),
        edges=[],
        threshold=0.0,
    )


class TestComputeTopology:
    def test_returns_topology(self, star_network):
        topo = compute_topology(star_network)
        assert isinstance(topo, NetworkTopology)
        assert topo.mag_ids == star_network.mag_ids
        assert len(topo.degree) == 5
        assert len(topo.betweenness) == 5
        assert len(topo.closeness) == 5
        assert len(topo.hub_scores) == 5
        assert len(topo.module_assignments) == 5
        assert isinstance(topo.modularity, float)

    def test_degree_correct(self, star_network):
        topo = compute_topology(star_network)
        idx = {m: i for i, m in enumerate(topo.mag_ids)}
        # Degree centrality: degree / (n-1) where n=5
        # M0 has 4 connections -> 4/4 = 1.0
        # M3 has 1 connection -> 1/4 = 0.25
        assert topo.degree[idx["M0"]] == pytest.approx(1.0)
        assert topo.degree[idx["M3"]] == pytest.approx(0.25)

    def test_betweenness_hub(self, star_network):
        """M0 (center of star) should have the highest betweenness."""
        topo = compute_topology(star_network)
        idx = {m: i for i, m in enumerate(topo.mag_ids)}
        m0_betweenness = topo.betweenness[idx["M0"]]
        for mag_id in ["M1", "M2", "M3", "M4"]:
            assert m0_betweenness >= topo.betweenness[idx[mag_id]], (
                f"M0 betweenness ({m0_betweenness:.4f}) should be >= "
                f"{mag_id} betweenness ({topo.betweenness[idx[mag_id]]:.4f})"
            )

    def test_modularity_computed(self, star_network):
        topo = compute_topology(star_network)
        # Modularity should be a finite number
        assert np.isfinite(topo.modularity)
        # Module assignments should be non-negative integers
        assert topo.module_assignments.dtype == int
        assert (topo.module_assignments >= 0).all()

    def test_empty_network(self, empty_network):
        topo = compute_topology(empty_network)
        assert isinstance(topo, NetworkTopology)
        np.testing.assert_array_equal(topo.degree, np.zeros(3))
        np.testing.assert_array_equal(topo.betweenness, np.zeros(3))
        np.testing.assert_array_equal(topo.closeness, np.zeros(3))
        np.testing.assert_array_equal(topo.hub_scores, np.zeros(3))
        assert topo.modularity == 0.0
        np.testing.assert_array_equal(topo.module_assignments, np.zeros(3, dtype=int))


class TestIdentifyKeystones:
    def test_hub_is_keystone(self, star_network, star_abundance):
        """M0 has highest centrality and lowest abundance -> keystone."""
        topo = compute_topology(star_network)
        result = identify_keystones(
            topo,
            star_abundance,
            betweenness_threshold=0.5,
            abundance_threshold=0.5,
        )
        assert isinstance(result, KeystoneTaxaResult)
        idx = {m: i for i, m in enumerate(result.mag_ids)}

        # M0 should be a keystone (high betweenness, low abundance)
        assert result.is_keystone[idx["M0"]], "M0 should be keystone"
        # M0 should have the highest keystone score
        assert result.keystone_scores[idx["M0"]] == result.keystone_scores.max()

        # Metrics dict should have expected keys
        assert "betweenness" in result.metrics
        assert "mean_abundance" in result.metrics
        assert "normalized_betweenness" in result.metrics
        assert "normalized_abundance" in result.metrics

    def test_empty_network_no_keystones(self, empty_network):
        abundance = AbundanceTable(
            mag_ids=["M0", "M1", "M2"],
            sample_ids=["s1", "s2"],
            abundances=np.array([[10.0, 20.0], [30.0, 40.0], [50.0, 60.0]]),
        )
        topo = compute_topology(empty_network)
        result = identify_keystones(topo, abundance)
        # No edges -> no betweenness -> no keystones
        assert not result.is_keystone.any()
        np.testing.assert_array_equal(result.keystone_scores, np.zeros(3))


class TestDifferentialNetwork:
    def test_finds_gained_lost_edges(self):
        """net1 has M0-M1 and M0-M2; net2 has M0-M1 and M1-M2.

        Conserved: M0-M1
        Lost:      M0-M2 (in net1, not in net2)
        Gained:    M1-M2 (in net2, not in net1)
        """
        mag_ids = ["M0", "M1", "M2"]
        n = len(mag_ids)

        adj1 = np.zeros((n, n), dtype=np.float64)
        adj1[0, 1] = adj1[1, 0] = 1.0
        adj1[0, 2] = adj1[2, 0] = 1.0
        net1 = NetworkResult(
            mag_ids=mag_ids,
            adjacency=adj1,
            edges=[("M0", "M1", 0.1), ("M0", "M2", 0.2)],
            threshold=0.3,
        )

        adj2 = np.zeros((n, n), dtype=np.float64)
        adj2[0, 1] = adj2[1, 0] = 1.0
        adj2[1, 2] = adj2[2, 1] = 1.0
        net2 = NetworkResult(
            mag_ids=mag_ids,
            adjacency=adj2,
            edges=[("M0", "M1", 0.1), ("M1", "M2", 0.15)],
            threshold=0.3,
        )

        result = differential_network(net1, net2)
        assert ("M0", "M1") in result["conserved"]
        assert ("M0", "M2") in result["lost"]
        assert ("M1", "M2") in result["gained"]
        assert len(result["conserved"]) == 1
        assert len(result["lost"]) == 1
        assert len(result["gained"]) == 1


class TestNetworkNullModel:
    def test_null_model_returns_result(self, star_network):
        result = network_null_model(star_network, n_iterations=50, seed=42)
        assert isinstance(result, NullModelResult)
        assert isinstance(result.observed_modularity, float)
        assert result.null_modularities.shape == (50,)
        assert 0.0 <= result.p_value <= 1.0
        assert np.isfinite(result.z_score)

    def test_null_model_empty_network(self, empty_network):
        result = network_null_model(empty_network, n_iterations=10)
        assert result.observed_modularity == 0.0
        assert result.z_score == 0.0
        assert result.p_value == 1.0

    def test_null_model_z_score_sign(self):
        """A well-structured clustered network should have positive z-score."""
        # Two clusters of 8 nodes each, densely connected within,
        # with only 2 bridge edges between clusters.
        cluster_a = [f"A{i}" for i in range(8)]
        cluster_b = [f"B{i}" for i in range(8)]
        mag_ids = cluster_a + cluster_b
        n = len(mag_ids)
        adjacency = np.zeros((n, n), dtype=np.float64)
        edges = []
        idx = {m: i for i, m in enumerate(mag_ids)}
        # Within-cluster edges (all pairs within each cluster)
        for cluster in [cluster_a, cluster_b]:
            for i in range(len(cluster)):
                for j in range(i + 1, len(cluster)):
                    edges.append((cluster[i], cluster[j], 0.1))
                    adjacency[idx[cluster[i]], idx[cluster[j]]] = 1.0
                    adjacency[idx[cluster[j]], idx[cluster[i]]] = 1.0
        # 2 bridge edges
        edges.append(("A0", "B0", 0.5))
        adjacency[idx["A0"], idx["B0"]] = 1.0
        adjacency[idx["B0"], idx["A0"]] = 1.0
        edges.append(("A1", "B1", 0.5))
        adjacency[idx["A1"], idx["B1"]] = 1.0
        adjacency[idx["B1"], idx["A1"]] = 1.0

        net = NetworkResult(
            mag_ids=mag_ids, adjacency=adjacency, edges=edges, threshold=0.6,
        )
        result = network_null_model(net, n_iterations=200, seed=42)
        assert result.z_score > 0, (
            f"Clustered network should have positive z-score, got {result.z_score}"
        )

    def test_null_model_deterministic(self, star_network):
        r1 = network_null_model(star_network, n_iterations=50, seed=99)
        r2 = network_null_model(star_network, n_iterations=50, seed=99)
        np.testing.assert_array_equal(r1.null_modularities, r2.null_modularities)
        assert r1.z_score == r2.z_score
        assert r1.p_value == r2.p_value


class TestModuleComposition:
    def test_returns_composition(self, star_network):
        taxonomy = TaxonomyTable(records={
            "M0": TaxonomyRecord(mag_id="M0", phylum="Proteobacteria"),
            "M1": TaxonomyRecord(mag_id="M1", phylum="Proteobacteria"),
            "M2": TaxonomyRecord(mag_id="M2", phylum="Actinobacteria"),
            "M3": TaxonomyRecord(mag_id="M3", phylum="Firmicutes"),
            "M4": TaxonomyRecord(mag_id="M4", phylum="Firmicutes"),
        })
        topo = compute_topology(star_network)
        result = module_composition(topo, taxonomy)
        assert isinstance(result, ModuleCompositionResult)
        assert len(result.module_ids) > 0
        # Total MAGs across modules should equal 5
        assert sum(result.n_mags) == 5
        # Dominant fractions in [0, 1]
        for frac in result.dominant_fraction:
            assert 0.0 <= frac <= 1.0

    def test_empty_network(self, empty_network):
        taxonomy = TaxonomyTable(records={
            "M0": TaxonomyRecord(mag_id="M0", phylum="P1"),
            "M1": TaxonomyRecord(mag_id="M1", phylum="P2"),
            "M2": TaxonomyRecord(mag_id="M2", phylum="P3"),
        })
        topo = compute_topology(empty_network)
        result = module_composition(topo, taxonomy)
        # All 3 nodes in module 0
        assert sum(result.n_mags) == 3


class TestHubBridgeClassification:
    def test_returns_result(self, star_network):
        topo = compute_topology(star_network)
        result = hub_bridge_classification(topo, star_network)
        assert isinstance(result, HubBridgeResult)
        assert len(result.mag_ids) == 5
        assert len(result.roles) == 5
        assert result.within_module_degree_z.shape == (5,)
        assert result.participation_coefficient.shape == (5,)
        # All roles should be valid
        for role in result.roles:
            assert role in {"hub", "bridge", "peripheral", "kinless"}

    def test_participation_in_range(self, star_network):
        topo = compute_topology(star_network)
        result = hub_bridge_classification(topo, star_network)
        # P must be in [0, 1]
        assert (result.participation_coefficient >= 0).all()
        assert (result.participation_coefficient <= 1).all()

    def test_clustered_network_has_bridges(self):
        """Two clusters connected by a single bridge node should yield bridge role."""
        cluster_a = [f"A{i}" for i in range(5)]
        cluster_b = [f"B{i}" for i in range(5)]
        bridge = ["X0"]
        mag_ids = cluster_a + bridge + cluster_b
        n = len(mag_ids)
        adjacency = np.zeros((n, n), dtype=np.float64)
        edges = []
        idx = {m: i for i, m in enumerate(mag_ids)}

        # Dense within-cluster
        for cluster in [cluster_a, cluster_b]:
            for i in range(len(cluster)):
                for j in range(i + 1, len(cluster)):
                    edges.append((cluster[i], cluster[j], 0.1))
                    adjacency[idx[cluster[i]], idx[cluster[j]]] = 1.0
                    adjacency[idx[cluster[j]], idx[cluster[i]]] = 1.0
        # X0 connects to A0 and B0
        for target in ["A0", "B0"]:
            edges.append(("X0", target, 0.2))
            adjacency[idx["X0"], idx[target]] = 1.0
            adjacency[idx[target], idx["X0"]] = 1.0

        net = NetworkResult(mag_ids=mag_ids, adjacency=adjacency, edges=edges, threshold=0.3)
        topo = compute_topology(net)
        result = hub_bridge_classification(topo, net)
        idx_x0 = result.mag_ids.index("X0")
        # X0 should have high participation (connects two modules)
        assert result.participation_coefficient[idx_x0] > 0.3

    def test_empty_network(self, empty_network):
        topo = compute_topology(empty_network)
        result = hub_bridge_classification(topo, empty_network)
        # All nodes peripheral with no edges
        assert all(r == "peripheral" for r in result.roles)
        np.testing.assert_array_equal(result.participation_coefficient, np.zeros(3))
        np.testing.assert_array_equal(result.within_module_degree_z, np.zeros(3))
