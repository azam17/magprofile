"""Ordination methods (PCoA, NMDS) for MAG community profiling."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .beta import BetaDiversityResult


@dataclass
class OrdinationResult:
    """Ordination coordinates and diagnostics."""

    sample_ids: list[str]
    coordinates: np.ndarray  # shape (n_samples, n_axes)
    explained_variance: np.ndarray | None  # per axis (PCoA only)
    stress: float | None  # NMDS only
    method: str


def pcoa(beta: BetaDiversityResult, n_axes: int = 2) -> OrdinationResult:
    """Principal Coordinates Analysis via classical MDS.

    Double-centers the squared distance matrix, eigendecomposes,
    and returns top-k axes with explained variance.
    """
    dm = beta.distance_matrix
    n = dm.shape[0]

    # Double-center the squared distance matrix
    d2 = dm**2
    row_mean = d2.mean(axis=1, keepdims=True)
    col_mean = d2.mean(axis=0, keepdims=True)
    grand_mean = d2.mean()
    B = -0.5 * (d2 - row_mean - col_mean + grand_mean)

    # Eigendecompose
    eigenvalues, eigenvectors = np.linalg.eigh(B)

    # Sort descending
    idx = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]

    # Keep positive eigenvalues only for coordinate calculation
    n_axes = min(n_axes, n - 1)
    pos = eigenvalues[:n_axes].clip(min=0)
    coords = eigenvectors[:, :n_axes] * np.sqrt(pos)[np.newaxis, :]

    # Explained variance as proportion of sum of positive eigenvalues
    total_pos = eigenvalues[eigenvalues > 0].sum()
    if total_pos > 0:
        explained = pos / total_pos
    else:
        explained = np.zeros(n_axes)

    return OrdinationResult(
        sample_ids=list(beta.sample_ids),
        coordinates=coords,
        explained_variance=explained,
        stress=None,
        method="PCoA",
    )


def nmds(
    beta: BetaDiversityResult,
    n_axes: int = 2,
    random_state: int = 42,
) -> OrdinationResult:
    """Non-metric Multidimensional Scaling via sklearn."""
    from sklearn.manifold import MDS

    mds = MDS(
        n_components=n_axes,
        metric=False,
        dissimilarity="precomputed",
        random_state=random_state,
        normalized_stress="auto",
    )
    coords = mds.fit_transform(beta.distance_matrix)
    return OrdinationResult(
        sample_ids=list(beta.sample_ids),
        coordinates=coords,
        explained_variance=None,
        stress=float(mds.stress_),
        method="NMDS",
    )
