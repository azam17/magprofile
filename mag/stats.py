"""Statistical tests for MAG community profiling."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy import stats as sp_stats

from .beta import BetaDiversityResult
from .io import AbundanceTable, SampleMetadata
from .ordination import OrdinationResult, pcoa


@dataclass
class StatisticalTestResult:
    """Result from a statistical test."""

    test_name: str
    statistic: float
    p_value: float
    metadata: dict[str, object] = field(default_factory=dict)


def permanova(
    beta: BetaDiversityResult,
    sample_metadata: SampleMetadata,
    grouping_var: str,
    n_permutations: int = 999,
    seed: int = 42,
) -> StatisticalTestResult:
    """PERMANOVA: partition distance matrix variance by group.

    F = (SS_between / (g-1)) / (SS_within / (n-g))
    Permutation test by shuffling group labels.
    """
    dm = beta.distance_matrix
    n = dm.shape[0]
    groups = sample_metadata.get_groups(grouping_var)

    # Map sample_id -> group label
    sample_to_group: dict[str, str] = {}
    for grp, sids in groups.items():
        for sid in sids:
            sample_to_group[sid] = grp

    labels = np.array([sample_to_group[s] for s in beta.sample_ids])
    unique_groups = np.unique(labels)
    g = len(unique_groups)

    def compute_f(lab: np.ndarray) -> float:
        # Total sum of squares (normalized by n)
        ss_total = 0.0
        for i in range(n):
            for j in range(i + 1, n):
                ss_total += dm[i, j] ** 2
        ss_total /= n

        # Within-group sum of squares (each group normalized by group size)
        ss_within = 0.0
        for grp in unique_groups:
            idx = np.where(lab == grp)[0]
            ng = len(idx)
            if ng > 0:
                grp_ss = sum(
                    dm[idx[a], idx[b]] ** 2
                    for a in range(len(idx))
                    for b in range(a + 1, len(idx))
                )
                ss_within += grp_ss / ng

        ss_between = ss_total - ss_within
        df_between = g - 1
        df_within = n - g

        if df_within == 0 or ss_within == 0:
            return 0.0
        return (ss_between / df_between) / (ss_within / df_within)

    observed_f = compute_f(labels)

    # Permutation test
    rng = np.random.default_rng(seed)
    count = 0
    for _ in range(n_permutations):
        perm_labels = rng.permutation(labels)
        if compute_f(perm_labels) >= observed_f:
            count += 1

    p_value = (count + 1) / (n_permutations + 1)

    # R² = SS_between / SS_total
    ss_total = sum(dm[i, j] ** 2 for i in range(n) for j in range(i + 1, n)) / n
    ss_within_norm = 0.0
    for grp in unique_groups:
        idx = np.where(labels == grp)[0]
        ng = len(idx)
        if ng > 0:
            grp_ss = sum(
                dm[idx[a], idx[b]] ** 2
                for a in range(len(idx))
                for b in range(a + 1, len(idx))
            )
            ss_within_norm += grp_ss / ng
    r_squared = (ss_total - ss_within_norm) / ss_total if ss_total > 0 else 0.0

    return StatisticalTestResult(
        test_name="PERMANOVA",
        statistic=observed_f,
        p_value=p_value,
        metadata={"R2": r_squared, "n_permutations": n_permutations},
    )


def anosim(
    beta: BetaDiversityResult,
    sample_metadata: SampleMetadata,
    grouping_var: str,
    n_permutations: int = 999,
    seed: int = 42,
) -> StatisticalTestResult:
    """ANOSIM: rank-based R statistic with permutation test."""
    dm = beta.distance_matrix
    n = dm.shape[0]
    groups = sample_metadata.get_groups(grouping_var)

    sample_to_group: dict[str, str] = {}
    for grp, sids in groups.items():
        for sid in sids:
            sample_to_group[sid] = grp

    labels = np.array([sample_to_group[s] for s in beta.sample_ids])

    # Get upper triangle distances and rank them
    tri_idx = np.triu_indices(n, k=1)
    dists = dm[tri_idx]
    ranks = sp_stats.rankdata(dists)

    def compute_r(lab: np.ndarray) -> float:
        # Separate within-group and between-group ranks
        within_ranks = []
        between_ranks = []
        k = 0
        for i in range(n):
            for j in range(i + 1, n):
                if lab[i] == lab[j]:
                    within_ranks.append(ranks[k])
                else:
                    between_ranks.append(ranks[k])
                k += 1
        if not within_ranks or not between_ranks:
            return 0.0
        r_b = np.mean(between_ranks)
        r_w = np.mean(within_ranks)
        m = len(dists)
        return (r_b - r_w) / (m / 2)

    observed_r = compute_r(labels)

    rng = np.random.default_rng(seed)
    count = 0
    for _ in range(n_permutations):
        perm_labels = rng.permutation(labels)
        if compute_r(perm_labels) >= observed_r:
            count += 1

    p_value = (count + 1) / (n_permutations + 1)

    return StatisticalTestResult(
        test_name="ANOSIM",
        statistic=observed_r,
        p_value=p_value,
        metadata={"n_permutations": n_permutations},
    )


def kruskal_wallis_per_mag(
    table: AbundanceTable,
    sample_metadata: SampleMetadata,
    grouping_var: str,
) -> list[StatisticalTestResult]:
    """Kruskal-Wallis test per MAG across groups."""
    groups = sample_metadata.get_groups(grouping_var)
    sample_to_group: dict[str, str] = {}
    for grp, sids in groups.items():
        for sid in sids:
            sample_to_group[sid] = grp

    # Build group index arrays
    group_names = sorted(groups.keys())
    group_indices: dict[str, list[int]] = {g: [] for g in group_names}
    for j, sid in enumerate(table.sample_ids):
        grp = sample_to_group.get(sid)
        if grp is not None:
            group_indices[grp].append(j)

    results = []
    for i, mag_id in enumerate(table.mag_ids):
        group_values = [table.abundances[i, group_indices[g]] for g in group_names]
        # Need at least 2 groups with data
        non_empty = [gv for gv in group_values if len(gv) > 0]
        if len(non_empty) < 2:
            results.append(
                StatisticalTestResult(
                    test_name="Kruskal-Wallis",
                    statistic=0.0,
                    p_value=1.0,
                    metadata={"mag_id": mag_id},
                )
            )
            continue
        stat, pval = sp_stats.kruskal(*non_empty)
        results.append(
            StatisticalTestResult(
                test_name="Kruskal-Wallis",
                statistic=float(stat),
                p_value=float(pval),
                metadata={"mag_id": mag_id},
            )
        )
    return results


def permdisp(
    beta: BetaDiversityResult,
    metadata: SampleMetadata,
    grouping_var: str,
    n_permutations: int = 999,
    seed: int = 42,
) -> StatisticalTestResult:
    """PERMDISP: test for homogeneity of multivariate dispersions.

    Computes group centroids in PCoA space, then tests whether
    distances from samples to their group centroid differ across groups.
    Uses ANOVA F-test on centroid distances with permutation p-value.
    """
    # Get PCoA coordinates
    ord_result = pcoa(beta)
    coords = ord_result.coordinates  # (n_samples, n_axes)

    groups = metadata.get_groups(grouping_var)
    sample_to_group: dict[str, str] = {}
    for grp, sids in groups.items():
        for sid in sids:
            sample_to_group[sid] = grp

    labels = np.array([sample_to_group[s] for s in beta.sample_ids])
    unique_groups = sorted(set(labels))

    # Compute group centroids
    centroids: dict[str, np.ndarray] = {}
    for g in unique_groups:
        mask = labels == g
        centroids[g] = coords[mask].mean(axis=0)

    # Distance of each sample to its group centroid
    distances = np.zeros(len(labels))
    for i, (s, g) in enumerate(zip(beta.sample_ids, labels)):
        distances[i] = np.linalg.norm(coords[i] - centroids[g])

    # Per-group mean distances
    group_mean_distances: dict[str, float] = {}
    for g in unique_groups:
        mask = labels == g
        group_mean_distances[g] = float(distances[mask].mean())

    # ANOVA F-test on distances
    def compute_f(lab: np.ndarray) -> float:
        group_vals = [distances[lab == g] for g in unique_groups]
        group_vals = [gv for gv in group_vals if len(gv) > 0]
        if len(group_vals) < 2:
            return 0.0
        grand_mean = distances.mean()
        n = len(distances)
        k = len(group_vals)
        ss_between = sum(
            len(gv) * (gv.mean() - grand_mean) ** 2 for gv in group_vals
        )
        ss_within = sum(np.sum((gv - gv.mean()) ** 2) for gv in group_vals)
        if ss_within == 0 or (n - k) == 0:
            return 0.0
        return (ss_between / (k - 1)) / (ss_within / (n - k))

    observed_f = compute_f(labels)

    # Permutation test
    rng = np.random.default_rng(seed)
    count = 0
    for _ in range(n_permutations):
        perm_labels = rng.permutation(labels)
        if compute_f(perm_labels) >= observed_f:
            count += 1

    p_value = (count + 1) / (n_permutations + 1)

    return StatisticalTestResult(
        test_name="PERMDISP",
        statistic=observed_f,
        p_value=p_value,
        metadata={
            "group_mean_distances": group_mean_distances,
            "n_permutations": n_permutations,
        },
    )


def pairwise_permanova(
    beta: BetaDiversityResult,
    metadata: SampleMetadata,
    grouping_var: str,
    n_permutations: int = 999,
    seed: int = 42,
) -> list[StatisticalTestResult]:
    """Pairwise PERMANOVA between all pairs of groups with BH-FDR correction.

    Subsets the distance matrix for each pair and runs PERMANOVA.
    Applies Benjamini-Hochberg correction across all pairwise p-values.
    """
    groups = metadata.get_groups(grouping_var)
    group_names = sorted(groups.keys())
    sid_to_idx = {s: i for i, s in enumerate(beta.sample_ids)}

    pairs: list[tuple[str, str]] = []
    raw_results: list[StatisticalTestResult] = []

    for i, g1 in enumerate(group_names):
        for g2 in group_names[i + 1 :]:
            pairs.append((g1, g2))

            # Subset samples for this pair
            g1_sids = [s for s in groups[g1] if s in sid_to_idx]
            g2_sids = [s for s in groups[g2] if s in sid_to_idx]
            all_sids = g1_sids + g2_sids
            indices = [sid_to_idx[s] for s in all_sids]

            # Subset distance matrix
            sub_dm = beta.distance_matrix[np.ix_(indices, indices)]
            sub_beta = BetaDiversityResult(
                sample_ids=all_sids,
                distance_matrix=sub_dm,
                metric=beta.metric,
            )

            # Build metadata for just these two groups
            sub_records: dict[str, dict[str, str]] = {}
            for s in g1_sids:
                sub_records[s] = {grouping_var: g1}
            for s in g2_sids:
                sub_records[s] = {grouping_var: g2}
            sub_meta = SampleMetadata(records=sub_records)

            result = permanova(
                sub_beta, sub_meta, grouping_var,
                n_permutations=n_permutations, seed=seed,
            )
            result.metadata["group1"] = g1
            result.metadata["group2"] = g2
            raw_results.append(result)

    # BH-FDR correction across all pairwise p-values
    if raw_results:
        p_values = np.array([r.p_value for r in raw_results])
        q_values = _bh_fdr(p_values)
        for r, q in zip(raw_results, q_values):
            r.metadata["q_value"] = float(q)

    return raw_results


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
    for i in range(n - 2, -1, -1):
        q[i] = min(q[i], q[i + 1])
    q = np.clip(q, 0, 1)
    result = np.zeros(n)
    result[sorted_idx] = q
    return result
