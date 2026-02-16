"""Information-theoretic compartment specificity score."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .io import AbundanceTable, SampleMetadata


@dataclass
class CompartmentSpecificityResult:
    """Compartment specificity for each MAG."""

    mag_ids: list[str]
    scores: np.ndarray  # 0 = generalist, 1 = specialist
    dominant_compartment: list[str]
    preferences: np.ndarray  # shape (n_mags, n_compartments), normalized mean abundance


def compartment_specificity(
    table: AbundanceTable,
    metadata: SampleMetadata,
    grouping_var: str,
) -> CompartmentSpecificityResult:
    """Compute information-theoretic specificity for each MAG.

    Specificity = 1 - H(P) / H_max
    P = normalized mean abundance per compartment
    H_max = ln(C), C = number of compartments
    """
    groups = metadata.get_groups(grouping_var)
    group_names = sorted(groups.keys())
    n_groups = len(group_names)

    sid_to_idx = {s: i for i, s in enumerate(table.sample_ids)}

    # Compute mean abundance per compartment for each MAG
    prefs = np.zeros((table.n_mags, n_groups))
    for gi, g in enumerate(group_names):
        indices = [sid_to_idx[s] for s in groups[g] if s in sid_to_idx]
        if indices:
            prefs[:, gi] = table.abundances[:, indices].mean(axis=1)

    # Normalize to probability distribution per MAG
    row_sums = prefs.sum(axis=1, keepdims=True)
    row_sums = np.where(row_sums == 0, 1.0, row_sums)
    P = prefs / row_sums

    # Shannon entropy per MAG
    H_max = np.log(n_groups) if n_groups > 1 else 1.0
    H = np.zeros(table.n_mags)
    for i in range(table.n_mags):
        p = P[i]
        p = p[p > 0]
        if len(p) > 0:
            H[i] = -np.sum(p * np.log(p))

    scores = 1.0 - H / H_max if H_max > 0 else np.ones(table.n_mags)

    # Dominant compartment
    dominant_idx = P.argmax(axis=1)
    dominant = [group_names[i] for i in dominant_idx]

    return CompartmentSpecificityResult(
        mag_ids=list(table.mag_ids),
        scores=scores,
        dominant_compartment=dominant,
        preferences=P,
    )
