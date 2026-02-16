"""Indicator species analysis (IndVal) for MAG community profiling."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .io import AbundanceTable, SampleMetadata


@dataclass
class IndicatorSpeciesResult:
    """IndVal indicator species results."""

    mag_ids: list[str]
    indval_scores: np.ndarray
    best_group: list[str]
    specificity: np.ndarray
    fidelity: np.ndarray
    p_values: np.ndarray


def indicator_species(
    table: AbundanceTable,
    metadata: SampleMetadata,
    grouping_var: str,
    n_permutations: int = 999,
    seed: int = 42,
) -> IndicatorSpeciesResult:
    """IndVal = Specificity x Fidelity with permutation test.

    Specificity = mean_abundance_in_group / sum(mean_abundance_all_groups)
    Fidelity = n_samples_present_in_group / n_samples_in_group
    """
    groups = metadata.get_groups(grouping_var)
    group_names = sorted(groups.keys())
    n_groups = len(group_names)

    # Build sample -> group index mapping
    sid_to_idx = {s: i for i, s in enumerate(table.sample_ids)}
    group_sample_indices: list[list[int]] = []
    for g in group_names:
        indices = [sid_to_idx[s] for s in groups[g] if s in sid_to_idx]
        group_sample_indices.append(indices)

    def compute_indval(
        abundances: np.ndarray, sample_groups: list[list[int]]
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str]]:
        """Compute IndVal for each MAG. Returns (indval, spec, fidel, best_group)."""
        n_mags = abundances.shape[0]
        spec = np.zeros((n_mags, n_groups))
        fidel = np.zeros((n_mags, n_groups))

        for gi, indices in enumerate(sample_groups):
            if not indices:
                continue
            grp_abund = abundances[:, indices]
            spec[:, gi] = grp_abund.mean(axis=1)
            fidel[:, gi] = (grp_abund > 0).sum(axis=1) / len(indices)

        # Normalize specificity
        spec_sum = spec.sum(axis=1, keepdims=True)
        spec_sum = np.where(spec_sum == 0, 1.0, spec_sum)
        spec_norm = spec / spec_sum

        indval = spec_norm * fidel
        best_idx = indval.argmax(axis=1)
        best_indval = indval[np.arange(n_mags), best_idx]
        best_spec = spec_norm[np.arange(n_mags), best_idx]
        best_fidel = fidel[np.arange(n_mags), best_idx]
        best_grp = [group_names[i] for i in best_idx]

        return best_indval, best_spec, best_fidel, best_grp

    observed_indval, obs_spec, obs_fidel, obs_best = compute_indval(
        table.abundances, group_sample_indices
    )

    # Permutation test
    rng = np.random.default_rng(seed)
    all_indices = list(range(table.n_samples))
    counts = np.zeros(table.n_mags)

    for _ in range(n_permutations):
        perm = rng.permutation(all_indices)
        # Reassign samples to groups
        perm_groups: list[list[int]] = []
        offset = 0
        for gi in range(n_groups):
            size = len(group_sample_indices[gi])
            perm_groups.append(list(perm[offset : offset + size]))
            offset += size
        perm_indval, _, _, _ = compute_indval(table.abundances, perm_groups)
        counts += (perm_indval >= observed_indval).astype(float)

    p_values = (counts + 1) / (n_permutations + 1)

    return IndicatorSpeciesResult(
        mag_ids=list(table.mag_ids),
        indval_scores=observed_indval,
        best_group=obs_best,
        specificity=obs_spec,
        fidelity=obs_fidel,
        p_values=p_values,
    )
