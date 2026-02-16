"""Differential abundance analysis with CLR transform."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import stats as sp_stats

from .io import AbundanceTable, SampleMetadata


@dataclass
class DifferentialAbundanceResult:
    """Results from differential abundance testing."""

    mag_ids: list[str]
    log_fold_changes: np.ndarray
    p_values: np.ndarray
    q_values: np.ndarray
    effect_sizes: np.ndarray  # Cohen's d


def clr_transform(abundances: np.ndarray, pseudocount: float = 1.0) -> np.ndarray:
    """Centered log-ratio transform.

    Args:
        abundances: shape (n_mags, n_samples)
        pseudocount: added to handle zeros

    Returns:
        CLR-transformed array of same shape.
    """
    x = abundances + pseudocount
    log_x = np.log(x)
    geo_mean = log_x.mean(axis=0, keepdims=True)
    return log_x - geo_mean


def cohens_d(group1: np.ndarray, group2: np.ndarray) -> float:
    """Cohen's d effect size (pooled standard deviation)."""
    n1, n2 = len(group1), len(group2)
    if n1 < 2 or n2 < 2:
        return 0.0
    var1 = np.var(group1, ddof=1)
    var2 = np.var(group2, ddof=1)
    pooled_std = np.sqrt(((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2))
    if pooled_std == 0:
        return 0.0
    return float((np.mean(group1) - np.mean(group2)) / pooled_std)


def differential_abundance(
    table: AbundanceTable,
    metadata: SampleMetadata,
    grouping_var: str,
    group1: str,
    group2: str,
) -> DifferentialAbundanceResult:
    """Differential abundance via CLR + Welch's t-test + BH-FDR.

    Compares group1 vs group2 for each MAG.
    """
    groups = metadata.get_groups(grouping_var)
    g1_samples = groups.get(group1, [])
    g2_samples = groups.get(group2, [])

    if not g1_samples or not g2_samples:
        raise ValueError(f"One or both groups empty: {group1}={len(g1_samples)}, {group2}={len(g2_samples)}")

    # Get sample indices
    sid_to_idx = {s: i for i, s in enumerate(table.sample_ids)}
    g1_idx = [sid_to_idx[s] for s in g1_samples if s in sid_to_idx]
    g2_idx = [sid_to_idx[s] for s in g2_samples if s in sid_to_idx]

    # CLR transform
    clr = clr_transform(table.abundances)

    n_mags = table.n_mags
    lfc = np.zeros(n_mags)
    pvals = np.zeros(n_mags)
    effects = np.zeros(n_mags)

    for i in range(n_mags):
        vals1 = clr[i, g1_idx]
        vals2 = clr[i, g2_idx]
        lfc[i] = np.mean(vals1) - np.mean(vals2)
        effects[i] = cohens_d(vals1, vals2)
        if len(vals1) >= 2 and len(vals2) >= 2:
            _, p = sp_stats.ttest_ind(vals1, vals2, equal_var=False)
            pvals[i] = p
        else:
            pvals[i] = 1.0

    # Benjamini-Hochberg FDR correction
    qvals = _bh_fdr(pvals)

    return DifferentialAbundanceResult(
        mag_ids=list(table.mag_ids),
        log_fold_changes=lfc,
        p_values=pvals,
        q_values=qvals,
        effect_sizes=effects,
    )


def _bh_fdr(p_values: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR correction."""
    n = len(p_values)
    if n == 0:
        return np.array([])
    sorted_idx = np.argsort(p_values)
    sorted_p = p_values[sorted_idx]
    q = np.zeros(n)
    for i in range(n):
        rank = i + 1
        q[i] = sorted_p[i] * n / rank
    # Enforce monotonicity from the end
    for i in range(n - 2, -1, -1):
        q[i] = min(q[i], q[i + 1])
    q = np.clip(q, 0, 1)
    # Unsort
    result = np.zeros(n)
    result[sorted_idx] = q
    return result
