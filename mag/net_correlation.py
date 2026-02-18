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

    phi_matrix = np.full((n_mags, n_mags), np.nan)

    for i in range(n_mags):
        if not valid[i]:
            continue
        for j in range(i, n_mags):
            if not valid[j]:
                continue
            if i == j:
                phi_matrix[i, j] = 0.0
            else:
                log_ratio = log_abundances[i] - log_abundances[j]
                phi_val = np.var(log_ratio, ddof=1)
                phi_matrix[i, j] = phi_val
                phi_matrix[j, i] = phi_val

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

    # Collect all valid upper-triangle phi values (non-NaN, non-diagonal)
    valid_values: list[float] = []
    for i in range(n):
        for j in range(i + 1, n):
            val = corr_matrix[i, j]
            if not np.isnan(val):
                valid_values.append(float(val))

    if not valid_values:
        threshold = 0.0
    else:
        threshold = float(np.percentile(valid_values, threshold_percentile))

    # Build adjacency and edge list
    adjacency = np.zeros((n, n), dtype=np.float64)
    edges: list[tuple[str, str, float]] = []

    for i in range(n):
        for j in range(i + 1, n):
            val = corr_matrix[i, j]
            if not np.isnan(val) and val <= threshold:
                adjacency[i, j] = 1.0
                adjacency[j, i] = 1.0
                edges.append((mag_ids[i], mag_ids[j], float(val)))

    return NetworkResult(
        mag_ids=list(mag_ids),
        adjacency=adjacency,
        edges=edges,
        threshold=threshold,
    )
