"""Functional profiling analysis: pathway abundance, differential pathways,
functional redundancy, CAZyme summaries, PGPR trait enrichment, and
keystone-PGPR cross-referencing."""

from __future__ import annotations

import re
from dataclasses import dataclass

import numpy as np
from scipy import stats as sp_stats

from .differential import DifferentialAbundanceResult, differential_abundance
from .func_io import FunctionalTable
from .io import AbundanceTable, SampleMetadata, TaxonomyTable


def pathway_abundance(
    func_table: FunctionalTable,
    abundance: AbundanceTable,
) -> AbundanceTable:
    """Project functional profiles into sample space.

    Computes the function-by-sample abundance matrix by multiplying the
    function-by-MAG matrix (*func_table*) with the MAG-by-sample
    *abundance* matrix.

    Parameters
    ----------
    func_table
        Function-by-MAG presence/count matrix.
    abundance
        MAG-by-sample abundance table.

    Returns
    -------
    AbundanceTable
        Function-by-sample abundance table (``mag_ids`` field holds
        function IDs).
    """
    return func_table.to_sample_abundance(abundance)


def differential_pathway(
    func_table: FunctionalTable,
    abundance: AbundanceTable,
    metadata: SampleMetadata,
    grouping_var: str,
    group1: str,
    group2: str,
) -> DifferentialAbundanceResult:
    """Differential pathway abundance between two sample groups.

    First projects the functional table into sample space via
    :func:`pathway_abundance`, then runs CLR + Welch's t-test +
    Benjamini-Hochberg FDR correction through
    :func:`~mag.differential.differential_abundance`.

    Parameters
    ----------
    func_table
        Function-by-MAG presence/count matrix.
    abundance
        MAG-by-sample abundance table.
    metadata
        Sample metadata with grouping information.
    grouping_var
        Metadata variable name that defines the groups.
    group1, group2
        The two group labels to compare.

    Returns
    -------
    DifferentialAbundanceResult
        Per-function log-fold changes, p-values, q-values, and effect sizes.
    """
    pathway_abund = pathway_abundance(func_table, abundance)
    return differential_abundance(
        pathway_abund, metadata, grouping_var, group1, group2
    )


def functional_redundancy(
    func_table: FunctionalTable,
    abundance: AbundanceTable,
    metadata: SampleMetadata,
    grouping_var: str,
) -> dict[str, dict[str, dict[str, float]]]:
    """Compute functional redundancy per group per function.

    For each sample group and each function, counts how many MAGs carry
    that function (carrier count) and computes the Shannon diversity of
    mean carrier abundances within the group.

    Parameters
    ----------
    func_table
        Function-by-MAG matrix (values > 0 indicate presence/counts).
    abundance
        MAG-by-sample abundance table.
    metadata
        Sample metadata.
    grouping_var
        Metadata variable that defines sample groups.

    Returns
    -------
    dict
        ``{group -> {func_id -> {"n_carriers": int, "shannon": float}}}``.
        Shannon diversity is 0.0 when there are zero or one carriers.
    """
    groups = metadata.get_groups(grouping_var)

    # Build index lookups for the abundance table
    abund_sid_to_idx = {s: i for i, s in enumerate(abundance.sample_ids)}
    abund_mid_to_idx = {m: i for i, m in enumerate(abundance.mag_ids)}

    # Pre-compute carrier mask: (n_functions, n_func_mags) boolean
    carrier_mask = func_table.values > 0

    # Pre-compute aligned MAG indices: func_table MAGs → abundance rows
    # For MAGs not in abundance table, use -1 as sentinel
    func_to_abund = np.array([
        abund_mid_to_idx.get(m, -1) for m in func_table.mag_ids
    ], dtype=np.intp)
    has_abund = func_to_abund >= 0

    result: dict[str, dict[str, dict[str, float]]] = {}

    for group_name, sample_ids in groups.items():
        group_result: dict[str, dict[str, float]] = {}

        # Indices of this group's samples in the abundance table
        sample_indices = np.array([
            abund_sid_to_idx[s] for s in sample_ids if s in abund_sid_to_idx
        ], dtype=np.intp)

        if len(sample_indices) > 0:
            # Pre-compute mean abundance per func_table MAG for this group
            # Only for MAGs that exist in the abundance table
            valid_abund_rows = func_to_abund[has_abund]
            mean_by_mag = np.zeros(len(func_table.mag_ids), dtype=np.float64)
            mean_by_mag[has_abund] = abundance.abundances[
                np.ix_(valid_abund_rows, sample_indices)
            ].mean(axis=1)

        for fi, func_id in enumerate(func_table.function_ids):
            carriers = carrier_mask[fi]  # boolean mask on func_table MAGs
            n_carriers = int(carriers.sum())

            if n_carriers == 0 or len(sample_indices) == 0:
                shannon = 0.0
            else:
                carrier_abund = mean_by_mag[carriers]
                total = carrier_abund.sum()
                if total > 0 and n_carriers > 1:
                    proportions = carrier_abund / total
                    proportions = proportions[proportions > 0]
                    shannon = float(-np.sum(proportions * np.log(proportions)))
                else:
                    shannon = 0.0

            group_result[func_id] = {
                "n_carriers": float(n_carriers),
                "shannon": shannon,
            }

        result[group_name] = group_result

    return result


# Regex to extract the CAZy class prefix (uppercase letters before the number)
_CAZY_CLASS_RE = re.compile(r"^([A-Z]+)")

# Recognised CAZy class prefixes
_CAZY_CLASSES = {"GH", "GT", "PL", "CE", "AA", "CBM"}


def cazyme_summary(
    func_table: FunctionalTable,
) -> dict[str, dict[str, float]]:
    """Aggregate CAZyme families by class prefix.

    Groups function IDs that start with a recognised CAZy class prefix
    (GH, GT, PL, CE, AA, CBM) and reports total gene counts, number of
    distinct families, and number of MAGs carrying at least one family
    in the class.

    Parameters
    ----------
    func_table
        A :class:`FunctionalTable` (typically of ``function_type="cazy"``).

    Returns
    -------
    dict
        ``{class -> {"total_genes": float, "n_families": int, "n_mags": int}}``.
    """
    # Map: class prefix -> list of row indices
    class_to_rows: dict[str, list[int]] = {}

    for fi, func_id in enumerate(func_table.function_ids):
        m = _CAZY_CLASS_RE.match(func_id)
        if m:
            prefix = m.group(1)
            if prefix in _CAZY_CLASSES:
                class_to_rows.setdefault(prefix, []).append(fi)

    result: dict[str, dict[str, float]] = {}

    for cls, row_indices in class_to_rows.items():
        # Total genes: sum of all values across those rows
        total_genes = float(func_table.values[row_indices, :].sum())

        # Number of distinct families in this class
        n_families = len(row_indices)

        # Number of MAGs that carry at least one family in this class
        # A MAG carries a family if value > 0 for any row in the class
        class_matrix = func_table.values[row_indices, :]  # (n_class_funcs, n_mags)
        mag_has_any = (class_matrix > 0).any(axis=0)  # (n_mags,)
        n_mags = int(mag_has_any.sum())

        result[cls] = {
            "total_genes": total_genes,
            "n_families": float(n_families),
            "n_mags": float(n_mags),
        }

    return result


# ---------------------------------------------------------------------------
# BH-FDR correction (local copy to avoid circular import)
# ---------------------------------------------------------------------------


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


# ---------------------------------------------------------------------------
# PGPR trait distribution by compartment
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class PGPREnrichmentResult:
    """Result of PGPR trait enrichment test across groups."""

    trait_names: list[str]
    group_names: list[str]
    # counts[trait][group] = number of MAGs carrying that trait in the group
    counts: dict[str, dict[str, int]]
    p_values: dict[str, float]
    q_values: dict[str, float]


def pgpr_enrichment(
    pgpr_table: FunctionalTable,
    abundance: AbundanceTable,
    metadata: SampleMetadata,
    grouping_var: str,
) -> PGPREnrichmentResult:
    """Test whether PGPR trait prevalence differs across sample groups.

    For each trait, builds a contingency table of carrier counts per group
    (a MAG is assigned to the group where it has highest mean abundance)
    and runs a chi-squared test (or Fisher's exact for 2x2).

    Parameters
    ----------
    pgpr_table
        Binary PGPR trait-by-MAG matrix (function_type="pgpr").
    abundance
        MAG-by-sample abundance table.
    metadata
        Sample metadata with grouping information.
    grouping_var
        Metadata variable that defines sample groups.

    Returns
    -------
    PGPREnrichmentResult
        Per-trait p-values and carrier counts per group.
    """
    groups = metadata.get_groups(grouping_var)
    group_names = sorted(groups.keys())
    abund_sid_to_idx = {s: i for i, s in enumerate(abundance.sample_ids)}
    abund_mid_to_idx = {m: i for i, m in enumerate(abundance.mag_ids)}

    # Assign each MAG in the PGPR table to a group (highest mean abundance)
    mag_group: dict[str, str] = {}
    for mag_id in pgpr_table.mag_ids:
        if mag_id not in abund_mid_to_idx:
            continue
        mi = abund_mid_to_idx[mag_id]
        best_group = ""
        best_mean = -1.0
        for g in group_names:
            sidx = [abund_sid_to_idx[s] for s in groups[g] if s in abund_sid_to_idx]
            if sidx:
                mean_a = float(abundance.abundances[mi, sidx].mean())
                if mean_a > best_mean:
                    best_mean = mean_a
                    best_group = g
        if best_group:
            mag_group[mag_id] = best_group

    # Build counts and run tests
    counts: dict[str, dict[str, int]] = {}
    raw_p: dict[str, float] = {}

    for ti, trait in enumerate(pgpr_table.function_ids):
        counts[trait] = {g: 0 for g in group_names}
        for mi, mag_id in enumerate(pgpr_table.mag_ids):
            if pgpr_table.values[ti, mi] > 0 and mag_id in mag_group:
                counts[trait][mag_group[mag_id]] += 1

        # Contingency: [carriers, non-carriers] x groups
        n_per_group: dict[str, int] = {}
        for g in group_names:
            n_per_group[g] = sum(1 for m, mg in mag_group.items() if mg == g)
        observed = np.array(
            [[counts[trait][g] for g in group_names],
             [n_per_group[g] - counts[trait][g] for g in group_names]]
        )
        if observed.sum() == 0 or observed.shape[1] < 2:
            raw_p[trait] = 1.0
        elif observed.shape[1] == 2:
            _, p = sp_stats.fisher_exact(observed)
            raw_p[trait] = float(p)
        else:
            try:
                _, p, _, _ = sp_stats.chi2_contingency(observed)
                raw_p[trait] = float(p)
            except ValueError:
                # Expected frequencies have zeros (sparse data)
                raw_p[trait] = 1.0

    # BH-FDR correction
    trait_list = list(raw_p.keys())
    p_arr = np.array([raw_p[t] for t in trait_list])
    q_arr = _bh_fdr(p_arr)
    q_values = {t: float(q_arr[i]) for i, t in enumerate(trait_list)}

    return PGPREnrichmentResult(
        trait_names=list(pgpr_table.function_ids),
        group_names=group_names,
        counts=counts,
        p_values=raw_p,
        q_values=q_values,
    )


# ---------------------------------------------------------------------------
# Keystone × PGPR cross-reference
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class KeystonePGPRResult:
    """Cross-reference of keystone taxa and PGPR traits."""

    keystone_mag_ids: list[str]
    # keystone_traits[mag_id] = list of trait names carried by this keystone
    keystone_traits: dict[str, list[str]]
    enrichment_p_value: float  # Fisher's exact: are keystones enriched for PGPR?
    n_keystones_with_pgpr: int
    n_keystones_total: int
    n_non_keystones_with_pgpr: int
    n_non_keystones_total: int


def keystone_pgpr_cross_reference(
    keystone_mag_ids: list[str],
    keystone_mask: np.ndarray,
    pgpr_table: FunctionalTable,
) -> KeystonePGPRResult:
    """Cross-reference keystone taxa with PGPR trait carriers.

    Tests whether keystone MAGs are enriched for PGPR traits compared to
    non-keystone MAGs using Fisher's exact test.

    Parameters
    ----------
    keystone_mag_ids
        MAG IDs from topology analysis (full list).
    keystone_mask
        Boolean array of same length indicating keystone status.
    pgpr_table
        Binary PGPR trait-by-MAG matrix.

    Returns
    -------
    KeystonePGPRResult
    """
    # Determine which MAGs carry any PGPR trait
    has_any_pgpr: set[str] = set()
    mag_traits: dict[str, list[str]] = {}
    for mi, mag_id in enumerate(pgpr_table.mag_ids):
        traits = []
        for ti, trait in enumerate(pgpr_table.function_ids):
            if pgpr_table.values[ti, mi] > 0:
                traits.append(trait)
        if traits:
            has_any_pgpr.add(mag_id)
            mag_traits[mag_id] = traits

    # Count keystones with/without PGPR
    ks_with = 0
    ks_total = 0
    non_ks_with = 0
    non_ks_total = 0
    keystone_trait_map: dict[str, list[str]] = {}

    for i, mag_id in enumerate(keystone_mag_ids):
        is_ks = bool(keystone_mask[i])
        has_pgpr = mag_id in has_any_pgpr
        if is_ks:
            ks_total += 1
            if has_pgpr:
                ks_with += 1
            keystone_trait_map[mag_id] = mag_traits.get(mag_id, [])
        else:
            non_ks_total += 1
            if has_pgpr:
                non_ks_with += 1

    # Fisher's exact test: keystones vs non-keystones enriched for PGPR
    table = np.array([
        [ks_with, ks_total - ks_with],
        [non_ks_with, non_ks_total - non_ks_with],
    ])
    if table.sum() == 0:
        p_value = 1.0
    else:
        _, p_value = sp_stats.fisher_exact(table, alternative="greater")

    keystone_ids = [
        keystone_mag_ids[i] for i in range(len(keystone_mag_ids)) if keystone_mask[i]
    ]

    return KeystonePGPRResult(
        keystone_mag_ids=keystone_ids,
        keystone_traits=keystone_trait_map,
        enrichment_p_value=float(p_value),
        n_keystones_with_pgpr=ks_with,
        n_keystones_total=ks_total,
        n_non_keystones_with_pgpr=non_ks_with,
        n_non_keystones_total=non_ks_total,
    )


# ---------------------------------------------------------------------------
# Functional redundancy comparison across compartments
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class RedundancyComparisonResult:
    """Comparison of functional redundancy across groups."""

    group_names: list[str]
    mean_shannon: dict[str, float]  # mean Shannon per group
    sd_shannon: dict[str, float]
    statistic: float  # Kruskal-Wallis H
    p_value: float
    pairwise: dict[tuple[str, str], float]  # pair -> Mann-Whitney p-value


def redundancy_comparison(
    redundancy: dict[str, dict[str, dict[str, float]]],
) -> RedundancyComparisonResult:
    """Compare functional redundancy (Shannon diversity) across groups.

    Takes the output of :func:`functional_redundancy` and tests whether
    mean redundancy differs across groups using Kruskal-Wallis, with
    pairwise Mann-Whitney U follow-up.

    Parameters
    ----------
    redundancy
        Output of :func:`functional_redundancy`:
        ``{group -> {func_id -> {"n_carriers": float, "shannon": float}}}``.

    Returns
    -------
    RedundancyComparisonResult
    """
    group_names = sorted(redundancy.keys())

    # Collect Shannon values per group (one value per function)
    group_shannons: dict[str, np.ndarray] = {}
    for g in group_names:
        vals = [v["shannon"] for v in redundancy[g].values()]
        group_shannons[g] = np.array(vals)

    # Descriptive stats
    mean_shannon = {g: float(group_shannons[g].mean()) for g in group_names}
    sd_shannon = {g: float(group_shannons[g].std()) for g in group_names}

    # Kruskal-Wallis
    arrays = [group_shannons[g] for g in group_names]
    if len(arrays) >= 2 and all(len(a) > 0 for a in arrays):
        try:
            h_stat, kw_p = sp_stats.kruskal(*arrays)
        except ValueError:
            h_stat, kw_p = 0.0, 1.0
    else:
        h_stat, kw_p = 0.0, 1.0

    # Pairwise Mann-Whitney U
    pairwise: dict[tuple[str, str], float] = {}
    for i, g1 in enumerate(group_names):
        for g2 in group_names[i + 1:]:
            try:
                _, p = sp_stats.mannwhitneyu(
                    group_shannons[g1], group_shannons[g2],
                    alternative="two-sided",
                )
                pairwise[(g1, g2)] = float(p)
            except ValueError:
                pairwise[(g1, g2)] = 1.0

    return RedundancyComparisonResult(
        group_names=group_names,
        mean_shannon=mean_shannon,
        sd_shannon=sd_shannon,
        statistic=float(h_stat),
        p_value=float(kw_p),
        pairwise=pairwise,
    )
