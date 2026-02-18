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


@dataclass(frozen=True)
class RarefactionResult:
    """Rarefaction curve data for richness and Shannon diversity."""

    depths: np.ndarray                # subsampling depths
    sample_ids: list[str]
    richness: np.ndarray              # shape (n_depths, n_samples) mean richness
    shannon: np.ndarray               # shape (n_depths, n_samples) mean Shannon
    richness_sd: np.ndarray           # shape (n_depths, n_samples) std
    shannon_sd: np.ndarray            # shape (n_depths, n_samples) std


def rarefaction_curves(
    abundance: AbundanceTable,
    n_depths: int = 20,
    n_iterations: int = 10,
    seed: int = 42,
) -> RarefactionResult:
    """Compute rarefaction curves by multinomial subsampling.

    For each subsampling depth, repeatedly draws counts from a multinomial
    distribution parameterised by the observed proportions, then computes
    richness and Shannon diversity.

    Parameters
    ----------
    abundance : AbundanceTable
        MAG-by-sample abundance matrix (raw counts).
    n_depths : int
        Number of evenly spaced subsampling depths.
    n_iterations : int
        Replicates per depth for variance estimation.
    seed : int
        Random seed for reproducibility.

    Returns
    -------
    RarefactionResult
        Mean and SD of richness and Shannon at each depth.
    """
    rng = np.random.default_rng(seed)
    n_samples = abundance.n_samples
    col_sums = abundance.abundances.sum(axis=0)

    # Depth range: 10% of smallest nonzero column sum to smallest column sum
    nonzero_sums = col_sums[col_sums > 0]
    if len(nonzero_sums) == 0:
        depths = np.array([0])
        return RarefactionResult(
            depths=depths,
            sample_ids=list(abundance.sample_ids),
            richness=np.zeros((1, n_samples)),
            shannon=np.zeros((1, n_samples)),
            richness_sd=np.zeros((1, n_samples)),
            shannon_sd=np.zeros((1, n_samples)),
        )

    min_sum = nonzero_sums.min()
    min_depth = max(1, int(min_sum * 0.1))
    max_depth = int(min_sum)
    if max_depth <= min_depth:
        max_depth = min_depth + 1
    depths = np.linspace(min_depth, max_depth, n_depths).astype(int)
    depths = np.unique(depths)  # remove duplicates from rounding
    n_d = len(depths)

    rich_mean = np.zeros((n_d, n_samples))
    shan_mean = np.zeros((n_d, n_samples))
    rich_sd = np.zeros((n_d, n_samples))
    shan_sd = np.zeros((n_d, n_samples))

    for j in range(n_samples):
        counts = abundance.abundances[:, j]
        total = counts.sum()
        if total == 0:
            continue
        props = counts / total

        for di, depth in enumerate(depths):
            if depth > total:
                rich_mean[di, j] = np.nan
                shan_mean[di, j] = np.nan
                rich_sd[di, j] = np.nan
                shan_sd[di, j] = np.nan
                continue

            iter_rich = np.empty(n_iterations)
            iter_shan = np.empty(n_iterations)
            for it in range(n_iterations):
                subsample = rng.multinomial(depth, props)
                iter_rich[it] = np.count_nonzero(subsample)
                p = subsample[subsample > 0] / subsample.sum()
                iter_shan[it] = -float(np.sum(p * np.log(p)))

            rich_mean[di, j] = iter_rich.mean()
            shan_mean[di, j] = iter_shan.mean()
            rich_sd[di, j] = iter_rich.std()
            shan_sd[di, j] = iter_shan.std()

    return RarefactionResult(
        depths=depths,
        sample_ids=list(abundance.sample_ids),
        richness=rich_mean,
        shannon=shan_mean,
        richness_sd=rich_sd,
        shannon_sd=shan_sd,
    )


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
