"""Alpha diversity metrics for MAG community profiling."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .io import AbundanceTable


@dataclass
class AlphaDiversityResult:
    """Alpha diversity metrics per sample."""

    sample_ids: list[str]
    shannon: np.ndarray
    simpson: np.ndarray
    richness: np.ndarray
    evenness: np.ndarray


def shannon_index(abundances: np.ndarray) -> float:
    """Shannon diversity H = -sum(pi * ln(pi)) for a single sample vector."""
    total = abundances.sum()
    if total == 0:
        return 0.0
    p = abundances / total
    p = p[p > 0]
    return -float(np.sum(p * np.log(p)))


def simpson_index(abundances: np.ndarray) -> float:
    """Simpson diversity D = 1 - sum(pi^2) for a single sample vector."""
    total = abundances.sum()
    if total == 0:
        return 0.0
    p = abundances / total
    return 1.0 - float(np.sum(p**2))


def richness(abundances: np.ndarray) -> int:
    """Number of MAGs with abundance > 0."""
    return int(np.sum(abundances > 0))


def pielou_evenness(abundances: np.ndarray) -> float:
    """Pielou's evenness J = H / ln(S)."""
    s = int(np.sum(abundances > 0))
    if s <= 1:
        return 0.0
    h = shannon_index(abundances)
    return h / np.log(s)


def compute_alpha_diversity(table: AbundanceTable) -> AlphaDiversityResult:
    """Compute alpha diversity for all samples in the abundance table."""
    n = table.n_samples
    h = np.zeros(n)
    d = np.zeros(n)
    r = np.zeros(n)
    j = np.zeros(n)
    for i in range(n):
        col = table.abundances[:, i]
        h[i] = shannon_index(col)
        d[i] = simpson_index(col)
        r[i] = richness(col)
        j[i] = pielou_evenness(col)
    return AlphaDiversityResult(
        sample_ids=list(table.sample_ids),
        shannon=h,
        simpson=d,
        richness=r,
        evenness=j,
    )
