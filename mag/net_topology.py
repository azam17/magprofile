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

    # Degree centrality
    degree_dict = nx.degree_centrality(G)
    degree = np.array([degree_dict[m] for m in mag_ids])

    # Betweenness centrality
    betweenness_dict = nx.betweenness_centrality(G)
    betweenness = np.array([betweenness_dict[m] for m in mag_ids])

    # Closeness centrality
    closeness_dict = nx.closeness_centrality(G)
    closeness = np.array([closeness_dict[m] for m in mag_ids])

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
