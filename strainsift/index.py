"""Dual-annotated inverted index for hierarchical classification.

StrainSift uses a two-tier indexing strategy:

1. **Coarse index** (k=21, FracMinHash) -- species-level screening.
   A small in-memory inverted index that maps each retained hash value to
   the set of species IDs whose reference genomes contain that k-mer.
   Querying a read against this index produces a ranked list of candidate
   species via FracMinHash containment.

2. **Fine index** (k=31, exact k-mers) -- strain-level classification.
   A larger inverted index loaded on demand per candidate species.  Each
   k-mer hash maps to a list of ``(strain_id, gene_family_id)`` tuples,
   enabling simultaneous strain and functional attribution.

The combined ``StrainSiftIndex`` orchestrates both tiers: a read is first
screened against the coarse index to identify plausible species, then the
fine index for those species is consulted for strain-level resolution.
"""

from __future__ import annotations

from collections import defaultdict
from typing import Any

import numpy as np

from strainsift.sketch import FracMinHash
from strainsift.types import (
    GeneFamilyAnnotation,
    ReadCandidate,
    StrainReference,
)
from strainsift.utils import (
    _MAX_HASH,
    canonical_kmer,
    extract_kmers,
    hash_kmers,
    murmurhash3_64,
)


# ======================================================================
# Coarse index -- species-level FracMinHash screening
# ======================================================================


class CoarseIndex:
    """Species-level inverted index using FracMinHash at k=21.

    For each species represented in the reference database, a FracMinHash
    sketch is built.  The inverted index maps individual hash values to the
    set of species that contain them, enabling fast containment queries
    without iterating over every species sketch.

    Attributes:
        k: k-mer size (default 21).
        scale: FracMinHash scale factor (default 1/1000).
        species_sketches: mapping from species_id to its merged FracMinHash.
        hash_to_species: inverted index from hash value to set of species IDs.
        strain_to_species: mapping from strain_id to species_id.
    """

    def __init__(self, k: int = 21, scale: float = 1 / 1000) -> None:
        self.k: int = k
        self.scale: float = scale
        # Per-species merged sketch (union of all strains in that species).
        self.species_sketches: dict[str, FracMinHash] = {}
        # Inverted index: hash -> {species_id, ...}
        self.hash_to_species: dict[int, set[str]] = defaultdict(set)
        # Book-keeping
        self.strain_to_species: dict[str, str] = {}

    # ------------------------------------------------------------------
    # Building
    # ------------------------------------------------------------------

    def add_strain(self, strain: StrainReference) -> None:
        """Add a strain's genome to the coarse index.

        The strain's k-mers are sketched with FracMinHash and merged into
        the species-level sketch.  Every retained hash is recorded in the
        inverted index.

        Parameters:
            strain: A :class:`StrainReference` with populated ``genome``
                and ``species_id`` fields.
        """
        species_id = strain.species_id
        self.strain_to_species[strain.strain_id] = species_id

        # Build a sketch for this strain's genome.
        strain_sketch = FracMinHash(k=self.k, scale=self.scale)
        strain_sketch.add_sequence(strain.genome)

        # Merge into species-level sketch.
        if species_id not in self.species_sketches:
            self.species_sketches[species_id] = strain_sketch.copy()
        else:
            self.species_sketches[species_id].merge(strain_sketch)

        # Update the inverted index with any new hashes.
        for h in strain_sketch._hash_set:
            self.hash_to_species[h].add(species_id)

    # ------------------------------------------------------------------
    # Querying
    # ------------------------------------------------------------------

    def query_read(
        self, read_seq: str, top_n: int = 3
    ) -> list[tuple[str, float]]:
        """Query a read against the coarse index.

        The read is sketched at the same ``k`` and ``scale``.  For every
        species whose sketch shares at least one hash with the read sketch,
        the FracMinHash containment of the read in the species is computed.

        Parameters:
            read_seq: DNA sequence of the read.
            top_n: Maximum number of candidate species to return.

        Returns:
            List of ``(species_id, containment_score)`` tuples sorted by
            score in descending order, limited to *top_n* entries.  Only
            species with a positive containment score are included.
        """
        # Build sketch for the read.
        read_sketch = FracMinHash(k=self.k, scale=self.scale)
        read_sketch.add_sequence(read_seq)

        if read_sketch.num_hashes == 0:
            # Read too short or all k-mers filtered -- fall back to exact
            # hash lookup without the FracMinHash threshold so we still
            # get a signal from short reads.
            return self._query_read_exact(read_seq, top_n)

        # Identify candidate species via inverted index.
        candidate_species: set[str] = set()
        for h in read_sketch._hash_set:
            if h in self.hash_to_species:
                candidate_species.update(self.hash_to_species[h])

        if not candidate_species:
            return []

        # Compute containment for each candidate.
        scores: list[tuple[str, float]] = []
        for sp_id in candidate_species:
            sp_sketch = self.species_sketches[sp_id]
            c = read_sketch.containment(sp_sketch)
            if c > 0.0:
                scores.append((sp_id, c))

        # Sort descending by score, take top_n.
        scores.sort(key=lambda x: x[1], reverse=True)
        return scores[:top_n]

    def _query_read_exact(
        self, read_seq: str, top_n: int
    ) -> list[tuple[str, float]]:
        """Fallback for reads whose FracMinHash sketch is empty.

        Computes exact k-mer containment against every species sketch by
        hashing all read k-mers (without the FracMinHash threshold) and
        checking membership.
        """
        kmers = extract_kmers(read_seq, self.k, canonical=True)
        if not kmers:
            return []
        read_hashes = set(hash_kmers(kmers))
        total = len(read_hashes)

        # Count shared hashes per species.
        species_hits: dict[str, int] = defaultdict(int)
        for h in read_hashes:
            if h in self.hash_to_species:
                for sp_id in self.hash_to_species[h]:
                    species_hits[sp_id] += 1

        scores = [
            (sp_id, hits / total)
            for sp_id, hits in species_hits.items()
            if hits > 0
        ]
        scores.sort(key=lambda x: x[1], reverse=True)
        return scores[:top_n]

    # ------------------------------------------------------------------
    # Introspection
    # ------------------------------------------------------------------

    @property
    def num_species(self) -> int:
        """Number of distinct species in the index."""
        return len(self.species_sketches)

    @property
    def num_hashes(self) -> int:
        """Total number of distinct hashes in the inverted index."""
        return len(self.hash_to_species)

    def __repr__(self) -> str:
        return (
            f"CoarseIndex(k={self.k}, scale={self.scale}, "
            f"species={self.num_species}, hashes={self.num_hashes})"
        )


# ======================================================================
# Fine index -- strain-level exact k-mer classification
# ======================================================================


class FineIndex:
    """Strain-level exact k-mer index at k=31.

    Each k-mer hash maps to a list of ``(strain_id, gene_family_id)``
    tuples.  The ``gene_family_id`` is ``None`` for intergenic k-mers,
    enabling reads to be attributed to both a strain and (optionally) a
    gene family in a single lookup pass.

    The index is designed to be loaded on demand for a single species at a
    time, so memory is proportional to one species' pan-genome rather than
    the entire database.

    Attributes:
        k: k-mer size (default 31).
        hash_to_annotations: inverted index from hash to list of
            (strain_id, gene_family_id | None) tuples.
        strains: set of strain IDs present in this fine index.
    """

    def __init__(self, k: int = 31) -> None:
        self.k: int = k
        # hash -> list[(strain_id, gene_family_id | None)]
        self.hash_to_annotations: dict[
            int, list[tuple[str, str | None]]
        ] = defaultdict(list)
        self.strains: set[str] = set()
        # Per-strain total distinct k-mer count (for diagnostics).
        self._strain_kmer_counts: dict[str, int] = {}

    # ------------------------------------------------------------------
    # Building
    # ------------------------------------------------------------------

    def add_strain(self, strain: StrainReference) -> None:
        """Add a strain's genome with gene annotations to the fine index.

        For every k-mer position in the genome, the k-mer is hashed and
        recorded together with the strain ID and the gene family that
        covers that genomic position (if any).

        Parameters:
            strain: A :class:`StrainReference` with ``genome`` and
                optionally ``gene_annotations`` populated.
        """
        self.strains.add(strain.strain_id)
        genome = strain.genome.upper()
        k = self.k

        # Pre-compute a position-to-gene lookup structure for efficiency.
        # We build a sorted interval list and do a simple scan, which is
        # adequate for typical bacterial genomes (< 20 k genes).
        gene_intervals = _build_gene_intervals(strain.gene_annotations)

        kmer_count = 0
        gene_idx = 0  # pointer into sorted gene_intervals for sweep

        for pos in range(len(genome) - k + 1):
            kmer_str = genome[pos: pos + k]
            if "N" in kmer_str:
                continue

            canon = canonical_kmer(kmer_str)
            h = murmurhash3_64(canon.encode())

            # Determine which gene family (if any) covers this position.
            gene_family_id = _gene_at_position(
                pos, gene_intervals, gene_idx
            )

            self.hash_to_annotations[h].append(
                (strain.strain_id, gene_family_id)
            )
            kmer_count += 1

        self._strain_kmer_counts[strain.strain_id] = kmer_count

    # ------------------------------------------------------------------
    # Querying
    # ------------------------------------------------------------------

    def query_read(
        self, read_seq: str, min_containment: float = 0.1
    ) -> list[ReadCandidate]:
        """Query a read against the fine index.

        All canonical k-mers of the read are hashed and looked up.  For
        each strain, we compute:

        * **containment** -- fraction of the read's distinct k-mer hashes
          found in the strain.
        * **gene family** -- the most frequently annotated gene family
          among the matching k-mers (``None`` if most hits are intergenic).

        Parameters:
            read_seq: DNA sequence of the read.
            min_containment: Minimum containment to include a strain in the
                result list (default 0.1).

        Returns:
            List of :class:`ReadCandidate` objects for every strain whose
            containment exceeds *min_containment*, sorted by containment
            descending.
        """
        kmers = extract_kmers(read_seq, self.k, canonical=True)
        if not kmers:
            return []

        hashes = hash_kmers(kmers)
        unique_hashes = set(hashes)
        total_unique = len(unique_hashes)

        # Accumulate per-strain hit counts and gene family tallies.
        # strain_id -> number of read k-mers found
        strain_hits: dict[str, int] = defaultdict(int)
        # strain_id -> {gene_family_id | None -> count}
        strain_gene_votes: dict[str, dict[str | None, int]] = defaultdict(
            lambda: defaultdict(int)
        )

        for h in unique_hashes:
            annotations = self.hash_to_annotations.get(h)
            if annotations is None:
                continue
            # Record one hit per strain (deduplicated per hash per strain).
            seen_strains_for_hash: set[str] = set()
            for strain_id, gf_id in annotations:
                if strain_id not in seen_strains_for_hash:
                    strain_hits[strain_id] += 1
                    seen_strains_for_hash.add(strain_id)
                strain_gene_votes[strain_id][gf_id] += 1

        # Build result list.
        candidates: list[ReadCandidate] = []
        for strain_id, hits in strain_hits.items():
            containment = hits / total_unique
            if containment < min_containment:
                continue

            # Determine the dominant gene family annotation.
            gene_votes = strain_gene_votes[strain_id]
            best_gene: str | None = max(gene_votes, key=gene_votes.get)  # type: ignore[arg-type]
            # If the dominant annotation is None (intergenic), keep it None.
            # If there is a tie between None and a gene, prefer the gene.
            if best_gene is None and len(gene_votes) > 1:
                # Check if any actual gene family ties with intergenic.
                none_count = gene_votes[None]
                gene_only = {
                    gf: c for gf, c in gene_votes.items() if gf is not None
                }
                if gene_only:
                    top_gene = max(gene_only, key=gene_only.get)  # type: ignore[arg-type]
                    if gene_only[top_gene] >= none_count:
                        best_gene = top_gene

            candidates.append(
                ReadCandidate(
                    strain_id=strain_id,
                    gene_family_id=best_gene,
                    containment=containment,
                    k=self.k,
                )
            )

        # Sort by containment descending.
        candidates.sort(key=lambda c: c.containment, reverse=True)
        return candidates

    # ------------------------------------------------------------------
    # Introspection
    # ------------------------------------------------------------------

    @property
    def num_strains(self) -> int:
        """Number of distinct strains in this fine index."""
        return len(self.strains)

    @property
    def num_kmers(self) -> int:
        """Number of distinct k-mer hashes in the inverted index."""
        return len(self.hash_to_annotations)

    def __repr__(self) -> str:
        return (
            f"FineIndex(k={self.k}, strains={self.num_strains}, "
            f"kmers={self.num_kmers})"
        )


# ======================================================================
# Combined hierarchical index
# ======================================================================


class StrainSiftIndex:
    """Combined coarse + fine index for hierarchical classification.

    Manages species-level screening via :class:`CoarseIndex` and on-demand
    strain-level resolution via per-species :class:`FineIndex` instances.

    Workflow:
        1. :meth:`build_from_references` ingests reference genomes and
           populates both index tiers.
        2. :meth:`query_read` performs hierarchical classification: first
           the coarse index identifies candidate species, then the fine
           index for each candidate species resolves strain and gene
           family.

    Attributes:
        coarse: The species-level :class:`CoarseIndex`.
        fine_indices: Mapping from species_id to :class:`FineIndex`.
        coarse_k: k-mer size for coarse screening (default 21).
        fine_k: k-mer size for fine classification (default 31).
        scale: FracMinHash scale (default 1/1000).
    """

    def __init__(
        self,
        coarse_k: int = 21,
        fine_k: int = 31,
        scale: float = 1 / 1000,
    ) -> None:
        self.coarse_k: int = coarse_k
        self.fine_k: int = fine_k
        self.scale: float = scale
        self.coarse: CoarseIndex = CoarseIndex(k=coarse_k, scale=scale)
        # Per-species fine indices, loaded/built on demand.
        self.fine_indices: dict[str, FineIndex] = {}
        # Track which strains belong to which species for on-demand loading.
        self._species_strains: dict[str, list[StrainReference]] = defaultdict(
            list
        )

    # ------------------------------------------------------------------
    # Building
    # ------------------------------------------------------------------

    def build_from_references(
        self, references: list[StrainReference]
    ) -> None:
        """Build both coarse and fine indices from reference genomes.

        Parameters:
            references: List of :class:`StrainReference` objects.  Each
                reference must have ``species_id``, ``strain_id``, and
                ``genome`` populated.
        """
        for strain in references:
            # --- Coarse tier ---
            self.coarse.add_strain(strain)

            # --- Fine tier ---
            species_id = strain.species_id
            if species_id not in self.fine_indices:
                self.fine_indices[species_id] = FineIndex(k=self.fine_k)
            self.fine_indices[species_id].add_strain(strain)

            # Keep a reference for potential future on-demand rebuilding.
            self._species_strains[species_id].append(strain)

    def add_strain(self, strain: StrainReference) -> None:
        """Add a single strain to both index tiers.

        Convenience method for incremental index building.
        """
        self.coarse.add_strain(strain)
        species_id = strain.species_id
        if species_id not in self.fine_indices:
            self.fine_indices[species_id] = FineIndex(k=self.fine_k)
        self.fine_indices[species_id].add_strain(strain)
        self._species_strains[species_id].append(strain)

    # ------------------------------------------------------------------
    # Querying
    # ------------------------------------------------------------------

    def query_read(
        self,
        read_seq: str,
        top_species: int = 3,
        min_containment: float = 0.1,
    ) -> list[ReadCandidate]:
        """Hierarchical query: coarse screening then fine classification.

        1. The read is screened against the coarse index to identify the
           *top_species* most likely species (by FracMinHash containment).
        2. For each candidate species, the read is queried against the
           corresponding fine index to obtain strain-level candidates.
        3. All strain-level candidates exceeding *min_containment* are
           returned, sorted by containment descending.

        Parameters:
            read_seq: DNA sequence of the read.
            top_species: Number of candidate species to carry forward from
                coarse screening (default 3).
            min_containment: Minimum strain-level containment to report
                (default 0.1).

        Returns:
            List of :class:`ReadCandidate` objects sorted by containment
            descending.  May be empty if no species passes coarse screening
            or no strain exceeds the containment threshold.
        """
        # --- Stage 1: coarse species screening ---
        species_candidates = self.coarse.query_read(
            read_seq, top_n=top_species
        )
        if not species_candidates:
            return []

        # --- Stage 2: fine strain classification ---
        all_candidates: list[ReadCandidate] = []
        for species_id, _coarse_score in species_candidates:
            fine_idx = self.fine_indices.get(species_id)
            if fine_idx is None:
                continue
            strain_candidates = fine_idx.query_read(
                read_seq, min_containment=min_containment
            )
            all_candidates.extend(strain_candidates)

        # Global sort by containment, best first.
        all_candidates.sort(key=lambda c: c.containment, reverse=True)
        return all_candidates

    # ------------------------------------------------------------------
    # On-demand fine index management
    # ------------------------------------------------------------------

    def load_fine_index(self, species_id: str) -> FineIndex | None:
        """Return the fine index for a species, building it if needed.

        Returns ``None`` if no strains have been registered for the
        given species.
        """
        if species_id in self.fine_indices:
            return self.fine_indices[species_id]

        strains = self._species_strains.get(species_id)
        if not strains:
            return None

        fine = FineIndex(k=self.fine_k)
        for s in strains:
            fine.add_strain(s)
        self.fine_indices[species_id] = fine
        return fine

    def get_species_for_strain(self, strain_id: str) -> str | None:
        """Return the species ID for a given strain ID.

        Searches across all registered species groups. Returns ``None``
        if the strain is not found.
        """
        for species_id, strains in self._species_strains.items():
            for s in strains:
                if s.strain_id == strain_id:
                    return species_id
        return None

    def unload_fine_index(self, species_id: str) -> None:
        """Remove a species' fine index from memory.

        The index can be rebuilt later via :meth:`load_fine_index` provided
        the strain references are still stored.
        """
        self.fine_indices.pop(species_id, None)

    # ------------------------------------------------------------------
    # Introspection
    # ------------------------------------------------------------------

    @property
    def num_species(self) -> int:
        """Number of species in the coarse index."""
        return self.coarse.num_species

    @property
    def num_strains(self) -> int:
        """Total number of strains across all fine indices."""
        return sum(fi.num_strains for fi in self.fine_indices.values())

    def __repr__(self) -> str:
        fine_loaded = len(self.fine_indices)
        return (
            f"StrainSiftIndex(coarse_k={self.coarse_k}, "
            f"fine_k={self.fine_k}, "
            f"species={self.num_species}, "
            f"fine_indices_loaded={fine_loaded})"
        )


# ======================================================================
# Internal helpers
# ======================================================================


def _build_gene_intervals(
    annotations: list[GeneFamilyAnnotation],
) -> list[tuple[int, int, str]]:
    """Build a sorted list of (start, end, gene_family_id) intervals.

    Sorted by start position for efficient sweep-line lookup.
    """
    intervals = [
        (ann.start, ann.end, ann.gene_family_id) for ann in annotations
    ]
    intervals.sort(key=lambda x: x[0])
    return intervals


def _gene_at_position(
    pos: int,
    intervals: list[tuple[int, int, str]],
    hint_idx: int = 0,
) -> str | None:
    """Return the gene family covering *pos*, or ``None`` if intergenic.

    Uses a simple linear scan starting from *hint_idx* for cache-friendly
    sequential access patterns (which is the common case when iterating
    k-mer positions left-to-right along a genome).

    Parameters:
        pos: 0-based genomic position.
        intervals: sorted list of (start, end, gene_family_id).
        hint_idx: starting index into *intervals* for the scan.

    Returns:
        Gene family ID string, or ``None`` if no interval covers *pos*.
    """
    n = len(intervals)
    if n == 0:
        return None

    # Start from hint; back up if needed.
    idx = min(hint_idx, n - 1)

    # Scan forward until we pass pos or find a covering interval.
    for i in range(n):
        start, end, gf_id = intervals[i]
        if start > pos:
            # Past the position -- no gene covers it.
            return None
        if start <= pos < end:
            return gf_id

    return None
