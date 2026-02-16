"""Beta diversity metrics for MAG community profiling."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.spatial.distance import pdist, squareform

from .io import AbundanceTable


@dataclass
class BetaDiversityResult:
    """Beta diversity dissimilarity matrices."""

    sample_ids: list[str]
    distance_matrix: np.ndarray  # shape (n_samples, n_samples)
    metric: str


def bray_curtis(table: AbundanceTable) -> BetaDiversityResult:
    """Compute Bray-Curtis dissimilarity between all sample pairs."""
    # pdist expects rows=observations, so transpose: (n_samples, n_mags)
    mat = table.abundances.T
    dists = pdist(mat, metric="braycurtis")
    dm = squareform(dists)
    return BetaDiversityResult(
        sample_ids=list(table.sample_ids),
        distance_matrix=dm,
        metric="bray_curtis",
    )


def jaccard(table: AbundanceTable) -> BetaDiversityResult:
    """Compute Jaccard dissimilarity on presence/absence."""
    pa = (table.abundances > 0).astype(float).T  # (n_samples, n_mags)
    dists = pdist(pa, metric="jaccard")
    dm = squareform(dists)
    return BetaDiversityResult(
        sample_ids=list(table.sample_ids),
        distance_matrix=dm,
        metric="jaccard",
    )
