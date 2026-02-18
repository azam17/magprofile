"""Functional profiling analysis: pathway abundance, differential pathways,
functional redundancy, and CAZyme summaries."""

from __future__ import annotations

import re

import numpy as np

from .differential import DifferentialAbundanceResult, differential_abundance
from .func_io import FunctionalTable
from .io import AbundanceTable, SampleMetadata


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

    # Map func_table MAG columns to abundance rows
    func_mag_to_idx = {m: i for i, m in enumerate(func_table.mag_ids)}

    result: dict[str, dict[str, dict[str, float]]] = {}

    for group_name, sample_ids in groups.items():
        group_result: dict[str, dict[str, float]] = {}

        # Indices of this group's samples in the abundance table
        sample_indices = [
            abund_sid_to_idx[s] for s in sample_ids if s in abund_sid_to_idx
        ]

        for fi, func_id in enumerate(func_table.function_ids):
            # Find carrier MAGs: those with value > 0 for this function
            carrier_mag_ids = [
                func_table.mag_ids[mi]
                for mi in range(len(func_table.mag_ids))
                if func_table.values[fi, mi] > 0
            ]

            n_carriers = len(carrier_mag_ids)

            # Compute Shannon diversity of mean abundances of carriers
            if n_carriers == 0 or len(sample_indices) == 0:
                shannon = 0.0
            else:
                # Mean abundance across group samples for each carrier MAG
                mean_abundances = []
                for mag_id in carrier_mag_ids:
                    if mag_id in abund_mid_to_idx:
                        mag_row = abund_mid_to_idx[mag_id]
                        vals = abundance.abundances[mag_row, sample_indices]
                        mean_abundances.append(float(np.mean(vals)))
                    else:
                        mean_abundances.append(0.0)

                total = sum(mean_abundances)
                if total > 0 and n_carriers > 1:
                    proportions = np.array(mean_abundances) / total
                    # Filter out zero proportions to avoid log(0)
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
