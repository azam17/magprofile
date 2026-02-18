"""Compositionally-aware co-occurrence network construction.

Implements the phi proportionality metric (Lovell et al., 2015) which is
appropriate for compositional data (e.g. relative abundances) unlike
Pearson/Spearman correlation.  Lower phi values indicate stronger
proportionality between two MAGs across samples.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .io import AbundanceTable


@dataclass(frozen=True)
class ThresholdSensitivityResult:
    """Sensitivity of network properties across phi thresholds."""

    percentiles: np.ndarray        # tested percentile values
    thresholds: np.ndarray         # phi threshold at each percentile
    n_edges: np.ndarray            # edge count at each percentile
    modularities: np.ndarray       # modularity at each percentile
    n_modules: np.ndarray          # module count at each percentile


@dataclass(frozen=True)
class NetworkResult:
    """Result of network construction from a proportionality matrix."""

    mag_ids: list[str]
    adjacency: np.ndarray  # shape (n_mags, n_mags)
    edges: list[tuple[str, str, float]]  # (mag_i, mag_j, phi) filtered by threshold
    threshold: float


def proportionality(
    table: AbundanceTable,
    pseudocount: float = 1.0,
    min_prevalence: float = 0.0,
) -> np.ndarray:
    """Compute the phi proportionality metric between all MAG pairs.

    phi(x, y) = var(log(x / y), ddof=1)

    Lower phi indicates stronger proportionality (co-occurrence).

    Parameters
    ----------
    table : AbundanceTable
        MAG-by-sample abundance matrix.
    pseudocount : float
        Added to all counts before log-transform to avoid log(0).
    min_prevalence : float
        Minimum fraction of samples in which a MAG must be present
        (abundance > 0) to be included.  MAGs below this threshold
        receive NaN in all their matrix entries.

    Returns
    -------
    np.ndarray
        Symmetric n_mags x n_mags matrix of phi values.  Diagonal is 0.
        Masked MAGs have NaN.
    """
    n_mags = table.n_mags
    abundances = table.abundances + pseudocount
    log_abundances = np.log(abundances)

    # Prevalence mask: fraction of samples where raw abundance > 0
    prevalence = (table.abundances > 0).sum(axis=1) / table.n_samples
    valid = prevalence >= min_prevalence

    # Vectorized: var(log_xi - log_xj) = var(log_xi) + var(log_xj) - 2*cov(log_xi, log_xj)
    cov_matrix = np.cov(log_abundances)  # (n_mags, n_mags), ddof=1
    if cov_matrix.ndim == 0:
        # Single MAG edge case
        cov_matrix = cov_matrix.reshape(1, 1)
    var_diag = np.diag(cov_matrix)  # (n_mags,)
    phi_matrix = var_diag[:, None] + var_diag[None, :] - 2 * cov_matrix
    np.fill_diagonal(phi_matrix, 0.0)

    # Mask low-prevalence MAGs with NaN
    phi_matrix[~valid, :] = np.nan
    phi_matrix[:, ~valid] = np.nan

    return phi_matrix


def build_network(
    corr_matrix: np.ndarray,
    mag_ids: list[str],
    threshold_percentile: float = 5,
) -> NetworkResult:
    """Build a co-occurrence network from a proportionality matrix.

    Edges are retained where phi <= threshold (low phi = strong
    association).  The threshold is computed at the given percentile
    of all valid (non-NaN, non-diagonal) phi values.

    Parameters
    ----------
    corr_matrix : np.ndarray
        Symmetric n_mags x n_mags phi matrix (from :func:`proportionality`).
    mag_ids : list[str]
        MAG identifiers matching matrix rows/columns.
    threshold_percentile : float
        Percentile (0-100) of phi values to use as edge threshold.

    Returns
    -------
    NetworkResult
        Network with adjacency matrix and edge list.
    """
    n = len(mag_ids)

    # Extract upper-triangle values with numpy indexing
    upper_i, upper_j = np.triu_indices(n, k=1)
    upper_vals = corr_matrix[upper_i, upper_j]
    valid_mask = ~np.isnan(upper_vals)
    valid_values = upper_vals[valid_mask]

    if len(valid_values) == 0:
        threshold = 0.0
    else:
        threshold = float(np.percentile(valid_values, threshold_percentile))

    # Build adjacency and edge list
    adjacency = np.zeros((n, n), dtype=np.float64)
    edge_mask = valid_mask & (upper_vals <= threshold)
    edge_indices = np.where(edge_mask)[0]

    # Set adjacency matrix (symmetric)
    ei = upper_i[edge_indices]
    ej = upper_j[edge_indices]
    adjacency[ei, ej] = 1.0
    adjacency[ej, ei] = 1.0

    edges: list[tuple[str, str, float]] = [
        (mag_ids[upper_i[k]], mag_ids[upper_j[k]], float(upper_vals[k]))
        for k in edge_indices
    ]

    return NetworkResult(
        mag_ids=list(mag_ids),
        adjacency=adjacency,
        edges=edges,
        threshold=threshold,
    )


def threshold_sensitivity(
    corr_matrix: np.ndarray,
    mag_ids: list[str],
    percentiles: np.ndarray | None = None,
) -> ThresholdSensitivityResult:
    """Assess network robustness across a range of phi thresholds.

    For each percentile, builds a network and computes topology metrics
    using :func:`build_network` and :func:`~mag.net_topology.compute_topology`.

    Parameters
    ----------
    corr_matrix : np.ndarray
        Symmetric phi proportionality matrix.
    mag_ids : list[str]
        MAG identifiers matching matrix rows/columns.
    percentiles : np.ndarray or None
        Percentile values to test (0–100). Defaults to 1–20.

    Returns
    -------
    ThresholdSensitivityResult
        Network properties at each tested percentile.
    """
    from .net_topology import compute_topology

    if percentiles is None:
        percentiles = np.arange(1, 21, dtype=np.float64)

    n_p = len(percentiles)
    thresholds = np.empty(n_p)
    n_edges = np.empty(n_p, dtype=int)
    modularities = np.empty(n_p)
    n_modules = np.empty(n_p, dtype=int)

    for i, p in enumerate(percentiles):
        net = build_network(corr_matrix, mag_ids, threshold_percentile=float(p))
        topo = compute_topology(net)
        thresholds[i] = net.threshold
        n_edges[i] = len(net.edges)
        modularities[i] = topo.modularity
        n_modules[i] = len(set(topo.module_assignments))

    return ThresholdSensitivityResult(
        percentiles=percentiles,
        thresholds=thresholds,
        n_edges=n_edges,
        modularities=modularities,
        n_modules=n_modules,
    )
