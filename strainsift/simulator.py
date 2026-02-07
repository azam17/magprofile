"""Synthetic metagenome simulator for StrainSift benchmarking.

Generates synthetic metagenomes with known ground truth, including:
- Random strain reference genomes with controlled pairwise ANI
- Gene family annotations (core, accessory, mobile)
- Illumina-like reads with position-dependent error profiles
- Open-universe mode where some strains are absent from the reference DB

All randomness flows through numpy's Generator API for full reproducibility.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from strainsift.types import (
    GeneFamilyAnnotation,
    SimulationGroundTruth,
    StrainReference,
)
from strainsift.utils import FastqRecord

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_BASES = np.array(list("ACGT"))
_BASE_TO_IDX = {"A": 0, "C": 1, "G": 2, "T": 3}

# Illumina error model: substitutions dominate, with rare indels.
_ERROR_TYPE_PROBS = np.array([0.90, 0.05, 0.05])  # sub, ins, del


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------


@dataclass
class SimulationConfig:
    """Configuration for synthetic metagenome generation.

    Attributes:
        n_strains: Number of strains to simulate.
        genome_length: Length of each strain genome in base pairs.
        n_gene_families: Number of distinct gene families to place.
        gene_length: Length in bp of each gene family instance.
        mobile_gene_fraction: Fraction of gene families annotated as mobile.
        accessory_gene_fraction: Fraction of gene families that are accessory
            (present in a random subset of strains rather than all).
        pairwise_ani: Target average nucleotide identity between strain pairs.
        abundances: Explicit abundance vector (length n_strains). If None,
            abundances are drawn from a symmetric Dirichlet(1) distribution.
        n_reads: Total number of reads to simulate.
        read_length: Length of each sequenced read in base pairs.
        error_rate: Mean per-base sequencing error rate.
        open_universe_fraction: Fraction of strains excluded from the
            returned reference database (simulates novel strains).
        random_seed: Seed for the numpy random number generator.
    """

    n_strains: int = 3
    genome_length: int = 50_000
    n_gene_families: int = 10
    gene_length: int = 1000
    mobile_gene_fraction: float = 0.2
    accessory_gene_fraction: float = 0.3
    pairwise_ani: float = 0.995
    abundances: np.ndarray | None = None
    n_reads: int = 10_000
    read_length: int = 150
    error_rate: float = 0.005
    open_universe_fraction: float = 0.0
    random_seed: int = 42


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _generate_template_genome(length: int, rng: np.random.Generator) -> str:
    """Generate a random DNA template genome with uniform base composition.

    Args:
        length: Number of bases.
        rng: Numpy random generator.

    Returns:
        A string of uppercase ACGT characters.
    """
    indices = rng.integers(0, 4, size=length)
    return "".join(_BASES[indices])


def _mutate_genome(template: str, ani: float, rng: np.random.Generator) -> str:
    """Create a mutant copy of *template* at a target ANI.

    Each position is independently mutated with probability ``1 - ani``.
    Mutations are always to a *different* base (uniform among the three
    alternatives), which matches the typical SNP model.

    Args:
        template: The source genome string.
        ani: Target average nucleotide identity (e.g. 0.995).
        rng: Numpy random generator.

    Returns:
        Mutated genome string.
    """
    mutation_rate = 1.0 - ani
    seq = list(template)
    n = len(seq)
    # Decide which positions to mutate.
    mutate_mask = rng.random(n) < mutation_rate
    n_mutations = int(mutate_mask.sum())
    if n_mutations == 0:
        return template

    for i in np.where(mutate_mask)[0]:
        original_idx = _BASE_TO_IDX[seq[i]]
        # Pick one of the 3 alternative bases uniformly.
        alt_offset = rng.integers(1, 4)
        new_idx = (original_idx + alt_offset) % 4
        seq[i] = _BASES[new_idx]
    return "".join(seq)


def _place_gene_families(
    n_gene_families: int,
    gene_length: int,
    genome_length: int,
    rng: np.random.Generator,
) -> list[tuple[int, int]]:
    """Choose non-overlapping positions for gene families on a genome.

    Genes are placed greedily at random positions. If the genome is too
    short to accommodate all genes, as many as possible are placed.

    Args:
        n_gene_families: Number of genes to place.
        gene_length: Length of each gene.
        genome_length: Total genome length.
        rng: Numpy random generator.

    Returns:
        Sorted list of (start, end) tuples.
    """
    if gene_length > genome_length:
        return []

    # Available space between genes (minimum 1 bp spacer).
    max_genes = genome_length // (gene_length + 1)
    actual_n = min(n_gene_families, max_genes)
    if actual_n == 0:
        return []

    # Evenly space candidate starts, then jitter within the available window.
    spacing = genome_length // actual_n
    positions: list[tuple[int, int]] = []
    for i in range(actual_n):
        window_start = i * spacing
        window_end = window_start + spacing - gene_length
        if window_end <= window_start:
            window_end = window_start + 1
        start = int(rng.integers(window_start, max(window_start + 1, window_end)))
        # Clamp to genome boundary.
        start = min(start, genome_length - gene_length)
        positions.append((start, start + gene_length))

    positions.sort()
    return positions


def _assign_gene_family_categories(
    n_gene_families: int,
    mobile_gene_fraction: float,
    accessory_gene_fraction: float,
    rng: np.random.Generator,
) -> tuple[list[bool], list[bool]]:
    """Classify each gene family as core/accessory and mobile/non-mobile.

    Args:
        n_gene_families: Total number of gene families.
        mobile_gene_fraction: Fraction of gene families that are mobile.
        accessory_gene_fraction: Fraction of gene families that are accessory
            (not shared by all strains).
        rng: Numpy random generator.

    Returns:
        Tuple of (is_mobile, is_accessory) boolean lists, each of length
        *n_gene_families*.
    """
    n_mobile = max(1, int(round(n_gene_families * mobile_gene_fraction)))
    n_accessory = max(0, int(round(n_gene_families * accessory_gene_fraction)))

    indices = rng.permutation(n_gene_families)
    mobile_set = set(indices[:n_mobile])
    # Accessory genes are drawn from the remaining pool *plus* some overlap
    # with mobile genes is allowed (mobile genes are often accessory).
    accessory_set = set(indices[:n_accessory])

    is_mobile = [i in mobile_set for i in range(n_gene_families)]
    is_accessory = [i in accessory_set for i in range(n_gene_families)]
    return is_mobile, is_accessory


def _build_strain_gene_presence(
    n_strains: int,
    n_gene_families: int,
    is_accessory: list[bool],
    rng: np.random.Generator,
) -> np.ndarray:
    """Build a boolean matrix (n_strains x n_gene_families) of gene presence.

    Core genes are present in every strain.  Accessory genes are present
    independently with probability 0.5 in each strain, but at least one
    strain always carries each gene.

    Args:
        n_strains: Number of strains.
        n_gene_families: Number of gene families.
        is_accessory: Per-gene boolean indicating accessory status.
        rng: Numpy random generator.

    Returns:
        Boolean ndarray of shape (n_strains, n_gene_families).
    """
    presence = np.ones((n_strains, n_gene_families), dtype=bool)
    for g in range(n_gene_families):
        if is_accessory[g]:
            # Each strain carries this gene independently with p=0.5.
            col = rng.random(n_strains) < 0.5
            # Ensure at least one strain carries the gene.
            if not col.any():
                col[rng.integers(n_strains)] = True
            presence[:, g] = col
    return presence


def _illumina_error_rate_at_position(
    pos: int,
    read_length: int,
    mean_error_rate: float,
) -> float:
    """Compute position-dependent Illumina error rate.

    Error rate increases roughly linearly from the 5' to the 3' end.
    At position 0 the rate is ``mean_error_rate * 0.5`` and at the last
    position it is ``mean_error_rate * 1.5``.

    Args:
        pos: 0-based position within the read.
        read_length: Total read length.
        mean_error_rate: Average per-base error rate across the read.

    Returns:
        Error rate for the given position.
    """
    if read_length <= 1:
        return mean_error_rate
    fraction = pos / (read_length - 1)
    # Linear ramp: 0.5x at start to 1.5x at end; integrates to 1x mean.
    scale = 0.5 + fraction
    return mean_error_rate * scale


def _apply_errors(
    sequence: str,
    error_rate: float,
    read_length: int,
    rng: np.random.Generator,
) -> tuple[str, str]:
    """Apply Illumina-like errors to a read and generate quality scores.

    Error model:
        - 90% substitution, 5% insertion, 5% deletion
        - Position-dependent rate (increases toward 3' end)
        - Quality scores anti-correlate with the actual error probability
          at each position, with stochastic noise to mimic real data.

    The output sequence length may differ slightly from *read_length* due
    to indels, but the quality string is always the same length as the
    returned sequence.

    Args:
        sequence: Error-free read sequence.
        error_rate: Mean per-base error rate.
        read_length: Nominal read length (for position-dependent model).
        rng: Numpy random generator.

    Returns:
        Tuple of (sequence_with_errors, quality_string).
    """
    result_bases: list[str] = []
    error_flags: list[bool] = []
    position_rates: list[float] = []

    i = 0  # index into the original (error-free) sequence
    while i < len(sequence):
        pos_rate = _illumina_error_rate_at_position(i, read_length, error_rate)
        position_rates.append(pos_rate)

        if rng.random() < pos_rate:
            # Error occurs -- choose error type.
            etype = rng.choice(3, p=_ERROR_TYPE_PROBS)

            if etype == 0:
                # Substitution.
                original_idx = _BASE_TO_IDX.get(sequence[i], 0)
                alt_offset = rng.integers(1, 4)
                new_idx = (original_idx + alt_offset) % 4
                result_bases.append(_BASES[new_idx])
                error_flags.append(True)
                i += 1

            elif etype == 1:
                # Insertion: insert a random base *before* consuming i.
                inserted_base = _BASES[rng.integers(0, 4)]
                result_bases.append(inserted_base)
                error_flags.append(True)
                position_rates.append(pos_rate)
                # Also emit the correct base at position i.
                result_bases.append(sequence[i])
                error_flags.append(False)
                i += 1

            else:
                # Deletion: skip the base at position i.
                error_flags.append(True)
                i += 1
                # We flagged an error but did not emit a base. To keep
                # error_flags aligned with result_bases, pop the flag and
                # just record nothing.  The position_rate was already
                # appended, so pop that too.
                error_flags.pop()
                position_rates.pop()
        else:
            result_bases.append(sequence[i])
            error_flags.append(False)
            i += 1

    seq_out = "".join(result_bases)

    # --- Generate quality scores ---
    # True Phred for each emitted base: Q = -10 log10(p_error).
    # We compute a "nominal" Q from the position-dependent rate, then add
    # Gaussian noise, and intentionally give *lower* Q where an error
    # actually occurred.
    quals: list[int] = []
    for j in range(len(result_bases)):
        # Use the position rate corresponding to this emitted base.
        pr = position_rates[j] if j < len(position_rates) else error_rate
        pr = max(pr, 1e-9)  # avoid log(0)
        nominal_q = -10.0 * np.log10(pr)
        # Add stochastic noise.
        noise = rng.normal(0, 2.0)
        q = nominal_q + noise
        # If an actual error occurred here, depress quality further.
        if error_flags[j]:
            q -= rng.uniform(5, 15)
        # Clamp to valid Phred range [2, 41] (Illumina 1.8+ encoding).
        q = max(2, min(41, int(round(q))))
        quals.append(q)

    quality_str = "".join(chr(q + 33) for q in quals)
    return seq_out, quality_str


def _generate_quality_for_perfect_read(
    read_length: int,
    error_rate: float,
    rng: np.random.Generator,
) -> str:
    """Generate a quality string for a read with no actual errors.

    Used as a fast path when error_rate is zero.

    Args:
        read_length: Length of the read.
        error_rate: Nominal error rate (used for Q computation).
        rng: Numpy random generator.

    Returns:
        Phred+33 encoded quality string.
    """
    quals: list[int] = []
    for pos in range(read_length):
        pr = _illumina_error_rate_at_position(pos, read_length, error_rate)
        pr = max(pr, 1e-9)
        nominal_q = -10.0 * np.log10(pr)
        noise = rng.normal(0, 2.0)
        q = max(2, min(41, int(round(nominal_q + noise))))
        quals.append(q)
    return "".join(chr(q + 33) for q in quals)


# ---------------------------------------------------------------------------
# Main simulation entry point
# ---------------------------------------------------------------------------


def simulate_metagenome(
    config: SimulationConfig,
) -> tuple[list[StrainReference], list[FastqRecord], SimulationGroundTruth]:
    """Generate a complete synthetic metagenome with ground truth.

    The simulation proceeds in five stages:

    1. **Template & Strain Genomes** -- A random template genome is created,
       then mutated *n_strains* times at the desired ANI to produce strain
       genomes.
    2. **Gene Family Placement** -- Non-overlapping gene regions are placed
       on the template; each gene family is classified as core, accessory,
       or mobile.  A presence/absence matrix records which genes each
       strain carries.
    3. **Abundance Assignment** -- If no explicit abundance vector is
       provided, one is drawn from a symmetric Dirichlet(1) distribution.
    4. **Read Sampling** -- Reads are drawn proportionally to strain
       abundances, each from a uniformly random start position on the
       source strain's genome.  Illumina-like errors are applied.
    5. **Open-Universe Partitioning** -- A fraction of strains is withheld
       from the returned reference database.

    Args:
        config: A :class:`SimulationConfig` controlling all parameters.

    Returns:
        A tuple of ``(references, reads, truth)`` where:

        - **references** is a list of :class:`StrainReference` objects for
          strains *in the database* (excludes open-universe strains).
        - **reads** is a list of :class:`FastqRecord` objects.
        - **truth** is a :class:`SimulationGroundTruth` recording the
          full ground truth for benchmarking.
    """
    rng = np.random.default_rng(config.random_seed)

    # ------------------------------------------------------------------
    # 1. Generate template and strain genomes
    # ------------------------------------------------------------------
    template = _generate_template_genome(config.genome_length, rng)
    strain_ids = [f"strain_{i:03d}" for i in range(config.n_strains)]
    species_id = "simulated_species"

    strain_genomes: list[str] = []
    for _ in range(config.n_strains):
        mutant = _mutate_genome(template, config.pairwise_ani, rng)
        strain_genomes.append(mutant)

    # ------------------------------------------------------------------
    # 2. Gene family placement and classification
    # ------------------------------------------------------------------
    gene_family_ids = [f"gene_{j:03d}" for j in range(config.n_gene_families)]

    # Place genes on the template coordinate system.
    gene_positions = _place_gene_families(
        config.n_gene_families,
        config.gene_length,
        config.genome_length,
        rng,
    )
    # If fewer positions were placed than families requested, truncate.
    n_placed = len(gene_positions)
    gene_family_ids = gene_family_ids[:n_placed]

    is_mobile, is_accessory = _assign_gene_family_categories(
        n_placed,
        config.mobile_gene_fraction,
        config.accessory_gene_fraction,
        rng,
    )
    mobile_gene_set = {
        gene_family_ids[g] for g in range(n_placed) if is_mobile[g]
    }

    # Strain-gene presence matrix.
    presence = _build_strain_gene_presence(
        config.n_strains, n_placed, is_accessory, rng
    )

    # Build StrainReference objects with annotations.
    all_references: list[StrainReference] = []
    for s in range(config.n_strains):
        annotations: list[GeneFamilyAnnotation] = []
        for g in range(n_placed):
            if presence[s, g]:
                start, end = gene_positions[g]
                annotations.append(
                    GeneFamilyAnnotation(
                        gene_family_id=gene_family_ids[g],
                        start=start,
                        end=end,
                        is_mobile=is_mobile[g],
                    )
                )
        ref = StrainReference(
            strain_id=strain_ids[s],
            species_id=species_id,
            genome=strain_genomes[s],
            gene_annotations=annotations,
            metadata={"pairwise_ani_target": config.pairwise_ani},
        )
        all_references.append(ref)

    # ------------------------------------------------------------------
    # 3. Abundance assignment
    # ------------------------------------------------------------------
    if config.abundances is not None:
        abundances = np.array(config.abundances, dtype=np.float64)
        if len(abundances) != config.n_strains:
            raise ValueError(
                f"abundances length ({len(abundances)}) != n_strains "
                f"({config.n_strains})"
            )
        abundances = abundances / abundances.sum()
    else:
        abundances = rng.dirichlet(np.ones(config.n_strains))

    # ------------------------------------------------------------------
    # 4. Read sampling
    # ------------------------------------------------------------------
    reads: list[FastqRecord] = []
    read_assignments: dict[str, str] = {}
    read_gene_families: dict[str, str | None] = {}

    # Pre-compute cumulative distribution for strain selection.
    cum_abundances = np.cumsum(abundances)

    for r_idx in range(config.n_reads):
        # Pick a source strain.
        u = rng.random()
        strain_idx = int(np.searchsorted(cum_abundances, u))
        strain_idx = min(strain_idx, config.n_strains - 1)
        genome = strain_genomes[strain_idx]
        sid = strain_ids[strain_idx]

        # Pick a random start position.  Wrap-around for circular genomes.
        start = rng.integers(0, config.genome_length)
        end = start + config.read_length
        if end <= config.genome_length:
            fragment = genome[start:end]
        else:
            # Wrap around.
            fragment = genome[start:] + genome[: end - config.genome_length]

        # Apply sequencing errors.
        if config.error_rate > 0:
            seq_out, qual_out = _apply_errors(
                fragment, config.error_rate, config.read_length, rng
            )
        else:
            seq_out = fragment
            qual_out = _generate_quality_for_perfect_read(
                config.read_length, config.error_rate, rng
            )

        read_id = f"read_{r_idx:07d}"
        reads.append(
            FastqRecord(read_id=read_id, sequence=seq_out, quality=qual_out)
        )
        read_assignments[read_id] = sid

        # Determine which gene family (if any) the read overlaps.
        gene_hit: str | None = None
        for g in range(n_placed):
            if not presence[strain_idx, g]:
                continue
            g_start, g_end = gene_positions[g]
            # Check overlap between [start, end) and [g_start, g_end).
            if start < g_end and end > g_start:
                gene_hit = gene_family_ids[g]
                break  # First overlapping gene wins.
        read_gene_families[read_id] = gene_hit

    # ------------------------------------------------------------------
    # 5. Open-universe partitioning
    # ------------------------------------------------------------------
    n_hidden = max(0, int(round(config.n_strains * config.open_universe_fraction)))
    # Choose which strains to hide (remove from DB).
    if n_hidden > 0:
        hidden_indices = set(
            rng.choice(config.n_strains, size=n_hidden, replace=False)
        )
    else:
        hidden_indices = set()

    strains_in_db = {
        strain_ids[i]
        for i in range(config.n_strains)
        if i not in hidden_indices
    }
    strains_not_in_db = {
        strain_ids[i]
        for i in range(config.n_strains)
        if i in hidden_indices
    }

    db_references = [
        ref
        for i, ref in enumerate(all_references)
        if i not in hidden_indices
    ]

    # ------------------------------------------------------------------
    # Build the true phi matrix (n_strains x n_gene_families).
    # phi[s, g] = 1.0 if strain s carries gene g, else 0.0.
    # ------------------------------------------------------------------
    true_phi = presence.astype(np.float64)

    truth = SimulationGroundTruth(
        true_pi=abundances,
        true_phi=true_phi,
        strain_ids=strain_ids,
        gene_family_ids=gene_family_ids,
        read_assignments=read_assignments,
        read_gene_families=read_gene_families,
        mobile_genes=mobile_gene_set,
        strains_in_db=strains_in_db,
        strains_not_in_db=strains_not_in_db,
    )

    return db_references, reads, truth
