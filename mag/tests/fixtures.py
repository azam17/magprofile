"""Synthetic data generation for MAG profiling tests."""

from __future__ import annotations

import numpy as np

from mag.io import AbundanceTable, SampleMetadata, TaxonomyRecord, TaxonomyTable


def generate_synthetic_abundance_table(
    n_mags: int = 20,
    n_samples: int = 12,
    n_groups: int = 3,
    seed: int = 42,
) -> tuple[AbundanceTable, SampleMetadata, TaxonomyTable]:
    """Generate synthetic data with planted group structure.

    Creates n_groups compartments with n_samples/n_groups samples each.
    Some MAGs are group-specific (high in one group, low elsewhere),
    others are generalists (similar across groups).
    """
    rng = np.random.default_rng(seed)
    samples_per_group = n_samples // n_groups
    group_names = ["topsoil", "bulk", "rhizosphere"][:n_groups]

    mag_ids = [f"MAG_{i:03d}" for i in range(n_mags)]
    sample_ids = []
    metadata_records: dict[str, dict[str, str]] = {}

    for gi, g in enumerate(group_names):
        for rep in range(samples_per_group):
            sid = f"{g}_rep{rep+1}"
            sample_ids.append(sid)
            metadata_records[sid] = {
                "compartment": g,
                "replicate": str(rep + 1),
                "site": "A",
            }

    abundances = np.zeros((n_mags, n_samples), dtype=np.float64)

    # First third: specialists for group 0
    # Second third: specialists for group 1
    # Last third: generalists (if n_groups >= 2)
    mags_per_class = n_mags // 3

    for i in range(n_mags):
        if i < mags_per_class:
            # Specialist for group 0
            target_group = 0
        elif i < 2 * mags_per_class:
            # Specialist for group 1
            target_group = 1
        else:
            # Generalist or specialist for remaining groups
            target_group = -1  # generalist

        for gi in range(n_groups):
            start = gi * samples_per_group
            end = start + samples_per_group
            if target_group == gi:
                abundances[i, start:end] = rng.poisson(500, samples_per_group)
            elif target_group == -1:
                abundances[i, start:end] = rng.poisson(200, samples_per_group)
            else:
                abundances[i, start:end] = rng.poisson(10, samples_per_group)

    # Taxonomy
    phyla = ["Proteobacteria", "Actinobacteriota", "Firmicutes", "Bacteroidota", "Acidobacteriota"]
    taxonomy_records: dict[str, TaxonomyRecord] = {}
    for i, mag_id in enumerate(mag_ids):
        taxonomy_records[mag_id] = TaxonomyRecord(
            mag_id=mag_id,
            domain="Bacteria",
            phylum=phyla[i % len(phyla)],
            class_=f"Class_{i % 3}",
            order=f"Order_{i % 4}",
            family=f"Family_{i % 5}",
            genus=f"Genus_{i % 7}",
            species=f"sp_{i:03d}",
        )

    return (
        AbundanceTable(mag_ids=mag_ids, sample_ids=sample_ids, abundances=abundances),
        SampleMetadata(records=metadata_records),
        TaxonomyTable(records=taxonomy_records),
    )


def generate_edge_case_all_zeros(n_mags: int = 5, n_samples: int = 3) -> AbundanceTable:
    """Abundance table with all zeros."""
    return AbundanceTable(
        mag_ids=[f"MAG_{i}" for i in range(n_mags)],
        sample_ids=[f"s_{i}" for i in range(n_samples)],
        abundances=np.zeros((n_mags, n_samples)),
    )


def generate_edge_case_single_mag() -> AbundanceTable:
    """Single MAG, multiple samples."""
    return AbundanceTable(
        mag_ids=["MAG_0"],
        sample_ids=["s_0", "s_1", "s_2"],
        abundances=np.array([[100.0, 200.0, 300.0]]),
    )


def generate_edge_case_single_sample() -> AbundanceTable:
    """Multiple MAGs, single sample."""
    return AbundanceTable(
        mag_ids=["MAG_0", "MAG_1", "MAG_2"],
        sample_ids=["s_0"],
        abundances=np.array([[100.0], [200.0], [300.0]]),
    )
