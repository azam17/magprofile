"""Core microbiome analysis for MAG community profiling."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .io import AbundanceTable, SampleMetadata


@dataclass
class CoreMicrobiomeResult:
    """Core microbiome analysis results."""

    thresholds: list[float]
    core_mags_per_threshold: dict[float, dict[str, list[str]]]
    # threshold -> group -> list of core MAG IDs
    shared_across_groups: dict[float, list[str]]
    # threshold -> MAG IDs core in ALL groups


def core_microbiome(
    table: AbundanceTable,
    metadata: SampleMetadata,
    grouping_var: str,
    thresholds: list[float] | None = None,
) -> CoreMicrobiomeResult:
    """Identify core MAGs at various prevalence thresholds.

    For each threshold, a MAG is "core" in a group if it is present
    (abundance > 0) in >= threshold fraction of samples within that group.
    """
    if thresholds is None:
        thresholds = [0.5, 0.75, 0.9, 1.0]

    groups = metadata.get_groups(grouping_var)
    sid_to_idx = {s: i for i, s in enumerate(table.sample_ids)}

    core_per_threshold: dict[float, dict[str, list[str]]] = {}
    shared: dict[float, list[str]] = {}

    for t in thresholds:
        core_per_threshold[t] = {}
        for g, sids in sorted(groups.items()):
            indices = [sid_to_idx[s] for s in sids if s in sid_to_idx]
            if not indices:
                core_per_threshold[t][g] = []
                continue
            n_samples = len(indices)
            presence = (table.abundances[:, indices] > 0).sum(axis=1)
            prevalence = presence / n_samples
            core_mags = [
                m for m, prev in zip(table.mag_ids, prevalence) if prev >= t
            ]
            core_per_threshold[t][g] = core_mags

        # MAGs core in ALL groups at this threshold
        group_sets = [
            set(mags) for mags in core_per_threshold[t].values() if mags is not None
        ]
        if group_sets:
            shared[t] = sorted(set.intersection(*group_sets))
        else:
            shared[t] = []

    return CoreMicrobiomeResult(
        thresholds=thresholds,
        core_mags_per_threshold=core_per_threshold,
        shared_across_groups=shared,
    )
