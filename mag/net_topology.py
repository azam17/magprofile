"""Network topology metrics and keystone taxa identification.

Computes graph-theoretic properties (degree, betweenness, closeness,
modularity, hub scores) from co-occurrence networks and identifies
keystone taxa as high-centrality, low-abundance MAGs.
"""

from __future__ import annotations

from dataclasses import dataclass

import networkx as nx
import numpy as np

from .io import AbundanceTable
from .net_correlation import NetworkResult


@dataclass(frozen=True)
class NullModelResult:
    """Result of network null model analysis (degree-preserving rewiring)."""

    observed_modularity: float
    null_modularities: np.ndarray  # shape (n_iterations,)
    z_score: float
    p_value: float  # fraction of null >= observed


@dataclass(frozen=True)
class NetworkTopology:
    """Graph-theoretic properties of a co-occurrence network."""

    mag_ids: list[str]
    degree: np.ndarray           # shape (n_mags,)
    betweenness: np.ndarray      # shape (n_mags,)
    closeness: np.ndarray        # shape (n_mags,)
    modularity: float
    module_assignments: np.ndarray  # shape (n_mags,) int
    hub_scores: np.ndarray       # shape (n_mags,)


@dataclass(frozen=True)
class KeystoneTaxaResult:
    """Keystone taxa identification results."""

    mag_ids: list[str]
    keystone_scores: np.ndarray   # shape (n_mags,)
    is_keystone: np.ndarray       # shape (n_mags,) bool
    metrics: dict[str, np.ndarray]


def compute_topology(network: NetworkResult) -> NetworkTopology:
    """Compute topology metrics for a co-occurrence network.

    Builds a networkx Graph from the edge list and computes degree,
    betweenness centrality, closeness centrality, greedy modularity
    communities, and HITS hub scores.

    Parameters
    ----------
    network : NetworkResult
        Network from :func:`~mag.net_correlation.build_network`.

    Returns
    -------
    NetworkTopology
        Frozen dataclass with all topology metrics.
    """
    mag_ids = list(network.mag_ids)
    n = len(mag_ids)

    G = nx.Graph()
    G.add_nodes_from(mag_ids)
    for mag_i, mag_j, weight in network.edges:
        G.add_edge(mag_i, mag_j, weight=weight)

    # Handle empty network (no edges)
    if G.number_of_edges() == 0:
        return NetworkTopology(
            mag_ids=mag_ids,
            degree=np.zeros(n, dtype=np.float64),
            betweenness=np.zeros(n, dtype=np.float64),
            closeness=np.zeros(n, dtype=np.float64),
            modularity=0.0,
            module_assignments=np.zeros(n, dtype=int),
            hub_scores=np.zeros(n, dtype=np.float64),
        )

    # Compute all centrality dicts
    degree_dict = nx.degree_centrality(G)
    betweenness_dict = nx.betweenness_centrality(G)
    closeness_dict = nx.closeness_centrality(G)

    # Single pass to build all arrays
    degree = np.empty(n, dtype=np.float64)
    betweenness = np.empty(n, dtype=np.float64)
    closeness = np.empty(n, dtype=np.float64)
    for i, m in enumerate(mag_ids):
        degree[i] = degree_dict.get(m, 0.0)
        betweenness[i] = betweenness_dict.get(m, 0.0)
        closeness[i] = closeness_dict.get(m, 0.0)

    # Community detection (greedy modularity)
    communities = list(nx.community.greedy_modularity_communities(G))
    modularity = nx.community.modularity(G, communities)
    module_assignments = np.zeros(n, dtype=int)
    mag_to_idx = {m: i for i, m in enumerate(mag_ids)}
    for comm_idx, community in enumerate(communities):
        for node in community:
            module_assignments[mag_to_idx[node]] = comm_idx

    # HITS hub scores (fall back to degree centrality on convergence failure)
    try:
        hubs, _authorities = nx.hits(G, max_iter=1000)
        hub_scores = np.array([hubs[m] for m in mag_ids])
    except nx.PowerIterationFailedConvergence:
        hub_scores = degree.copy()

    return NetworkTopology(
        mag_ids=mag_ids,
        degree=degree,
        betweenness=betweenness,
        closeness=closeness,
        modularity=modularity,
        module_assignments=module_assignments,
        hub_scores=hub_scores,
    )


def identify_keystones(
    topology: NetworkTopology,
    abundance: AbundanceTable,
    betweenness_threshold: float = 0.5,
    abundance_threshold: float = 0.5,
) -> KeystoneTaxaResult:
    """Identify keystone taxa from topology and abundance data.

    A keystone taxon has high betweenness centrality (important network
    connector) but low mean abundance (rare yet structurally important).

    keystone_score = normalized_betweenness * (1 - normalized_abundance)

    A MAG is classified as keystone if its betweenness is at or above the
    ``betweenness_threshold`` percentile AND its mean abundance is at or
    below the ``abundance_threshold`` percentile.

    Parameters
    ----------
    topology : NetworkTopology
        Topology metrics from :func:`compute_topology`.
    abundance : AbundanceTable
        MAG-by-sample abundance matrix.
    betweenness_threshold : float
        Percentile (0-1) for betweenness centrality to qualify as high.
    abundance_threshold : float
        Percentile (0-1) for abundance to qualify as low.

    Returns
    -------
    KeystoneTaxaResult
        Keystone scores and boolean classification for each MAG.
    """
    mag_ids = topology.mag_ids
    n = len(mag_ids)

    # Align abundance to topology MAG ordering
    abund_idx = {m: i for i, m in enumerate(abundance.mag_ids)}
    mean_abundance = np.zeros(n, dtype=np.float64)
    for i, mag_id in enumerate(mag_ids):
        if mag_id in abund_idx:
            mean_abundance[i] = abundance.abundances[abund_idx[mag_id]].mean()

    betweenness = topology.betweenness

    # Normalize to [0, 1] range
    bmax = betweenness.max()
    norm_betweenness = betweenness / bmax if bmax > 0 else np.zeros(n)

    amax = mean_abundance.max()
    norm_abundance = mean_abundance / amax if amax > 0 else np.zeros(n)

    # Keystone score
    keystone_scores = norm_betweenness * (1.0 - norm_abundance)

    # Boolean classification using percentile thresholds
    if bmax > 0:
        b_cutoff = np.percentile(betweenness, betweenness_threshold * 100)
        high_betweenness = betweenness >= b_cutoff
    else:
        high_betweenness = np.zeros(n, dtype=bool)

    if amax > 0:
        a_cutoff = np.percentile(mean_abundance, abundance_threshold * 100)
        low_abundance = mean_abundance <= a_cutoff
    else:
        low_abundance = np.zeros(n, dtype=bool)

    is_keystone = high_betweenness & low_abundance

    metrics: dict[str, np.ndarray] = {
        "betweenness": betweenness,
        "mean_abundance": mean_abundance,
        "normalized_betweenness": norm_betweenness,
        "normalized_abundance": norm_abundance,
    }

    return KeystoneTaxaResult(
        mag_ids=mag_ids,
        keystone_scores=keystone_scores,
        is_keystone=is_keystone,
        metrics=metrics,
    )


def differential_network(
    net1: NetworkResult,
    net2: NetworkResult,
) -> dict[str, list[tuple[str, str]]]:
    """Compare two networks to find conserved, gained, and lost edges.

    Edges are normalized as sorted tuples (alphabetical order) for
    consistent comparison.

    Parameters
    ----------
    net1 : NetworkResult
        First (reference/baseline) network.
    net2 : NetworkResult
        Second (comparison) network.

    Returns
    -------
    dict with keys "conserved", "gained", "lost"
        Each value is a list of (mag_a, mag_b) tuples.
        - conserved: edges present in both networks
        - gained: edges in net2 but not in net1
        - lost: edges in net1 but not in net2
    """
    def _edge_set(network: NetworkResult) -> set[tuple[str, str]]:
        edges: set[tuple[str, str]] = set()
        for mag_i, mag_j, _weight in network.edges:
            edge = tuple(sorted((mag_i, mag_j)))
            edges.add(edge)  # type: ignore[arg-type]
        return edges

    edges1 = _edge_set(net1)
    edges2 = _edge_set(net2)

    conserved = sorted(edges1 & edges2)
    gained = sorted(edges2 - edges1)
    lost = sorted(edges1 - edges2)

    return {
        "conserved": conserved,
        "gained": gained,
        "lost": lost,
    }


def network_null_model(
    network: NetworkResult,
    n_iterations: int = 1000,
    seed: int = 42,
) -> NullModelResult:
    """Test whether observed modularity exceeds a degree-preserving null model.

    Uses double-edge-swap rewiring to generate random graphs that preserve
    the degree sequence, then compares observed modularity against the null
    distribution.

    Parameters
    ----------
    network : NetworkResult
        Network from :func:`~mag.net_correlation.build_network`.
    n_iterations : int
        Number of null-model randomisations.
    seed : int
        Random seed for reproducibility.

    Returns
    -------
    NullModelResult
        Observed modularity, null distribution, z-score, and p-value.
    """
    G = nx.Graph()
    G.add_nodes_from(network.mag_ids)
    for mag_i, mag_j, weight in network.edges:
        G.add_edge(mag_i, mag_j, weight=weight)

    n_edges = G.number_of_edges()

    # Empty network: nothing to test
    if n_edges == 0:
        return NullModelResult(
            observed_modularity=0.0,
            null_modularities=np.zeros(n_iterations),
            z_score=0.0,
            p_value=1.0,
        )

    # Observed modularity
    communities = list(nx.community.greedy_modularity_communities(G))
    observed_mod = nx.community.modularity(G, communities)

    # Null distribution via degree-preserving rewiring
    rng = np.random.default_rng(seed)
    null_mods = np.empty(n_iterations)
    for i in range(n_iterations):
        G_copy = G.copy()
        try:
            nx.double_edge_swap(
                G_copy,
                nswap=n_edges,
                max_tries=n_edges * 10,
                seed=int(rng.integers(0, 2**31)),
            )
        except nx.NetworkXAlgorithmError:
            pass  # keep partially-rewired graph
        comms = list(nx.community.greedy_modularity_communities(G_copy))
        null_mods[i] = nx.community.modularity(G_copy, comms)

    # Statistics
    null_mean = null_mods.mean()
    null_std = null_mods.std()
    z_score = (observed_mod - null_mean) / null_std if null_std > 0 else 0.0
    p_value = (np.sum(null_mods >= observed_mod) + 1) / (n_iterations + 1)

    return NullModelResult(
        observed_modularity=observed_mod,
        null_modularities=null_mods,
        z_score=z_score,
        p_value=p_value,
    )
