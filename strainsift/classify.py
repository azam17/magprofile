"""Hierarchical read classification for StrainSift.

Two-stage classification pipeline:
1. Coarse screening (k=21, FracMinHash) -> candidate species
2. Fine classification (k=31, exact k-mer containment) -> candidate strains + gene families

The output feeds into the joint EM inference engine.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from strainsift.index import StrainSiftIndex
from strainsift.types import ReadCandidate, ReadClassification, StrainReference
from strainsift.utils import FastqRecord, parse_fastq

logger = logging.getLogger(__name__)


@dataclass
class ClassificationConfig:
    """Configuration for read classification.

    Attributes:
        top_species: Number of top species candidates from coarse screening.
        min_containment: Minimum k-mer containment to retain a strain candidate.
        min_read_length: Minimum read length to attempt classification.
        coarse_k: k-mer size for coarse screening.
        fine_k: k-mer size for fine strain-level classification.
        scale: FracMinHash scale parameter.
    """

    top_species: int = 3
    min_containment: float = 0.1
    min_read_length: int = 50
    coarse_k: int = 21
    fine_k: int = 31
    scale: float = 1 / 1000


@dataclass
class ClassificationStats:
    """Summary statistics from a classification run.

    Attributes:
        total_reads: Total reads processed.
        classified_reads: Reads with at least one strain candidate.
        unclassified_reads: Reads with no strain candidates.
        too_short_reads: Reads shorter than min_read_length.
        mean_candidates_per_read: Average number of strain candidates per classified read.
        species_read_counts: Dict mapping species_id -> number of reads assigned.
    """

    total_reads: int = 0
    classified_reads: int = 0
    unclassified_reads: int = 0
    too_short_reads: int = 0
    mean_candidates_per_read: float = 0.0
    species_read_counts: dict[str, int] = field(default_factory=dict)

    @property
    def classified_fraction(self) -> float:
        if self.total_reads == 0:
            return 0.0
        return self.classified_reads / self.total_reads

    @property
    def unclassified_fraction(self) -> float:
        if self.total_reads == 0:
            return 0.0
        return self.unclassified_reads / self.total_reads


class ReadClassifier:
    """Hierarchical coarse→fine read classification engine.

    Classifies reads against a StrainSiftIndex using a two-stage approach:
    1. Coarse: FracMinHash sketch at k=21 identifies candidate species
    2. Fine: Exact k-mer containment at k=31 identifies candidate strains
       and their associated gene families

    Usage:
        index = StrainSiftIndex()
        index.build_from_references(references)
        classifier = ReadClassifier(index)
        results, stats = classifier.classify_reads(reads)
    """

    def __init__(
        self,
        index: StrainSiftIndex,
        config: ClassificationConfig | None = None,
    ):
        self.index = index
        self.config = config or ClassificationConfig()

    def classify_read(self, read: FastqRecord) -> ReadClassification:
        """Classify a single read using hierarchical coarse→fine approach.

        Args:
            read: A FastqRecord to classify.

        Returns:
            ReadClassification with candidate strain assignments.
        """
        result = ReadClassification(read_id=read.read_id)

        # Skip reads that are too short for meaningful k-mer analysis
        if len(read.sequence) < self.config.min_read_length:
            return result

        # Query the hierarchical index
        candidates = self.index.query_read(
            read_seq=read.sequence,
            top_species=self.config.top_species,
            min_containment=self.config.min_containment,
        )

        if candidates:
            result.candidates = candidates
            result.is_classified = True
            # Extract unique species from candidates
            species_seen = set()
            for c in candidates:
                # Species is encoded in strain_id or available from index
                species_id = self.index.get_species_for_strain(c.strain_id)
                if species_id:
                    species_seen.add(species_id)
            result.species_candidates = list(species_seen)

        return result

    def classify_reads(
        self, reads: list[FastqRecord]
    ) -> tuple[list[ReadClassification], ClassificationStats]:
        """Classify a batch of reads.

        Args:
            reads: List of FastqRecord objects.

        Returns:
            Tuple of (classifications, stats).
        """
        stats = ClassificationStats()
        classifications = []
        total_candidates = 0

        for read in reads:
            stats.total_reads += 1

            if len(read.sequence) < self.config.min_read_length:
                stats.too_short_reads += 1
                stats.unclassified_reads += 1
                classifications.append(
                    ReadClassification(read_id=read.read_id)
                )
                continue

            result = self.classify_read(read)
            classifications.append(result)

            if result.is_classified:
                stats.classified_reads += 1
                total_candidates += len(result.candidates)
                for sp in result.species_candidates:
                    stats.species_read_counts[sp] = (
                        stats.species_read_counts.get(sp, 0) + 1
                    )
            else:
                stats.unclassified_reads += 1

        if stats.classified_reads > 0:
            stats.mean_candidates_per_read = (
                total_candidates / stats.classified_reads
            )

        logger.info(
            "Classification complete: %d/%d reads classified (%.1f%%)",
            stats.classified_reads,
            stats.total_reads,
            stats.classified_fraction * 100,
        )

        return classifications, stats

    def classify_fastq(
        self, fastq_path: Path
    ) -> tuple[list[ReadClassification], ClassificationStats]:
        """Classify reads from a FASTQ file.

        Args:
            fastq_path: Path to FASTQ file (plain or gzipped).

        Returns:
            Tuple of (classifications, stats).
        """
        reads = list(parse_fastq(fastq_path))
        return self.classify_reads(reads)


def _auto_scale(references: list[StrainReference], read_length: int = 150) -> float:
    """Choose a FracMinHash scale that works for the given genome sizes.

    For large genomes (>1Mb), scale=1/1000 is appropriate.
    For smaller genomes (simulated data), we need a higher scale so that
    short reads produce at least a few sketch hashes.

    Target: a read of `read_length` with k=21 should produce >=5 hashes.
    A read has (read_length - k + 1) k-mers. We want:
        scale * n_kmers >= 5  =>  scale >= 5 / n_kmers
    """
    n_kmers = max(read_length - 21 + 1, 1)
    min_scale = 5.0 / n_kmers  # ~0.038 for 150bp reads

    if not references:
        return min_scale

    avg_genome = sum(r.genome_length for r in references) / len(references)

    # For real genomes (>1Mb), use standard 1/1000
    if avg_genome >= 1_000_000:
        return max(1 / 1000, min_scale)

    # For smaller genomes, scale up proportionally
    # At 50kb, use ~1/50; at 500kb, use ~1/500
    genome_scale = 1.0 / max(avg_genome / 50, 1)
    return max(genome_scale, min_scale)


def classify_reads_batch(
    reads: list[FastqRecord],
    references: list[StrainReference],
    config: ClassificationConfig | None = None,
) -> tuple[list[ReadClassification], ClassificationStats]:
    """Convenience function: build index and classify reads in one call.

    Useful for small datasets and testing. For large datasets, build the
    index once and reuse it across samples.

    Automatically selects an appropriate FracMinHash scale if the config
    uses the default scale and the genomes are small (e.g., simulated data).

    Args:
        reads: List of reads to classify.
        references: List of strain references for the database.
        config: Classification configuration.

    Returns:
        Tuple of (classifications, stats).
    """
    config = config or ClassificationConfig()

    # Auto-tune scale for small genomes
    scale = config.scale
    if references and config.scale == 1 / 1000:
        avg_genome = sum(r.genome_length for r in references) / len(references)
        if avg_genome < 1_000_000:
            scale = _auto_scale(references)
            logger.info("Auto-tuned FracMinHash scale to %.4f for %d bp genomes",
                        scale, int(avg_genome))

    index = StrainSiftIndex(
        coarse_k=config.coarse_k,
        fine_k=config.fine_k,
        scale=scale,
    )
    index.build_from_references(references)
    classifier = ReadClassifier(index, config)
    return classifier.classify_reads(reads)
