"""FracMinHash sketch implementation for coarse species-level screening.

FracMinHash retains all hash values below a threshold determined by
``max_hash * scale``, rather than keeping a fixed number of minimum hashes
(as traditional MinHash does).  This allows meaningful set-algebraic
operations (union, intersection, containment) between sketches of different
sizes, which is critical for comparing metagenomic reads of varying length
against large reference genomes.

Default parameters:
    k = 21          (species-level resolution)
    scale = 1/1000  (retain ~0.1 % of all k-mer hashes)

References:
    Irber et al., "Lightweight compositional analysis of metagenomes with
    FracMinHash and minimum metagenome covers", bioRxiv (2022).
"""

from __future__ import annotations

import numpy as np

from strainsift.utils import (
    _MAX_HASH,
    extract_kmers,
    hash_kmers,
)


class FracMinHash:
    """FracMinHash sketch of a DNA sequence.

    A FracMinHash sketch keeps every k-mer hash whose value falls below a
    deterministic threshold ``max_hash * scale``.  Because the threshold is
    global and independent of the input, two sketches built with the same
    ``k`` and ``scale`` are directly comparable via set intersection.

    Attributes:
        k: k-mer size used during sketching (default 21).
        scale: fraction of the 64-bit hash space to retain (default 1/1000).
        max_hash_threshold: the integer cutoff -- hashes strictly below this
            value are kept.
        hashes: sorted numpy array (uint64) of retained hash values.
    """

    # ------------------------------------------------------------------
    # Construction
    # ------------------------------------------------------------------

    def __init__(self, k: int = 21, scale: float = 1 / 1000) -> None:
        if k < 1:
            raise ValueError(f"k must be >= 1, got {k}")
        if not 0 < scale <= 1:
            raise ValueError(f"scale must be in (0, 1], got {scale}")

        self.k: int = k
        self.scale: float = scale
        # Integer threshold: keep h when h < max_hash_threshold.
        # Using int() floors the result, which is the conservative choice.
        self.max_hash_threshold: int = int(_MAX_HASH * scale)
        # Internal storage -- a Python set for fast O(1) insertion during
        # building; converted to a sorted numpy array on demand.
        self._hash_set: set[int] = set()
        self._sorted_cache: np.ndarray | None = None

    # ------------------------------------------------------------------
    # Adding sequences
    # ------------------------------------------------------------------

    def add_sequence(self, seq: str) -> None:
        """Add all k-mers from a DNA sequence to the sketch.

        Each k-mer is canonicalised (lexicographically smaller of forward and
        reverse complement), hashed, and retained if the hash value is below
        ``self.max_hash_threshold``.

        Parameters:
            seq: DNA string (ACGT).  Non-ACGT characters cause the
                overlapping k-mers to be silently skipped (handled inside
                :func:`extract_kmers`).
        """
        kmers = extract_kmers(seq, self.k, canonical=True)
        if not kmers:
            return
        hashes = hash_kmers(kmers)
        threshold = self.max_hash_threshold
        new_hashes = {h for h in hashes if h < threshold}
        if new_hashes:
            self._hash_set.update(new_hashes)
            self._sorted_cache = None  # invalidate cache

    def add_hash(self, h: int) -> None:
        """Add a single pre-computed hash value (if below threshold).

        This is useful when building a sketch from an external hash stream.
        """
        if h < self.max_hash_threshold:
            self._hash_set.add(h)
            self._sorted_cache = None

    # ------------------------------------------------------------------
    # Sorted hash array (lazy)
    # ------------------------------------------------------------------

    @property
    def hashes(self) -> np.ndarray:
        """Sorted numpy array (dtype uint64) of retained hash values."""
        if self._sorted_cache is None:
            if self._hash_set:
                self._sorted_cache = np.array(
                    sorted(self._hash_set), dtype=np.uint64
                )
            else:
                self._sorted_cache = np.empty(0, dtype=np.uint64)
        return self._sorted_cache

    @property
    def num_hashes(self) -> int:
        """Number of hash values currently in the sketch."""
        return len(self._hash_set)

    # ------------------------------------------------------------------
    # Set operations
    # ------------------------------------------------------------------

    def intersection_size(self, other: FracMinHash) -> int:
        """Return ``|self ∩ other|`` (number of shared hashes).

        Both sketches must use the same ``k`` and ``scale``.
        """
        self._check_compatible(other)
        return len(self._hash_set & other._hash_set)

    def union_size(self, other: FracMinHash) -> int:
        """Return ``|self ∪ other|`` (number of distinct hashes in union)."""
        self._check_compatible(other)
        return len(self._hash_set | other._hash_set)

    # ------------------------------------------------------------------
    # Similarity / containment
    # ------------------------------------------------------------------

    def containment(self, other: FracMinHash) -> float:
        """Estimate containment of *self* in *other*.

        .. math::

            C(A, B) = \\frac{|A \\cap B|}{|A|}

        A containment of 1.0 means every hash in *self* also appears in
        *other* (i.e. all k-mers of the query are found in the reference).

        Returns 0.0 when *self* is empty.
        """
        if self.num_hashes == 0:
            return 0.0
        self._check_compatible(other)
        return self.intersection_size(other) / self.num_hashes

    def jaccard(self, other: FracMinHash) -> float:
        """Estimate Jaccard similarity between two sketches.

        .. math::

            J(A, B) = \\frac{|A \\cap B|}{|A \\cup B|}

        Returns 0.0 when both sketches are empty.
        """
        u = self.union_size(other)
        if u == 0:
            return 0.0
        return self.intersection_size(other) / u

    # ------------------------------------------------------------------
    # Merge
    # ------------------------------------------------------------------

    def merge(self, other: FracMinHash) -> None:
        """Merge another sketch into this one (in-place union).

        After merging, this sketch contains the union of both hash sets.
        """
        self._check_compatible(other)
        if other._hash_set:
            self._hash_set.update(other._hash_set)
            self._sorted_cache = None

    def copy(self) -> FracMinHash:
        """Return an independent deep copy of this sketch."""
        new = FracMinHash(k=self.k, scale=self.scale)
        new._hash_set = set(self._hash_set)
        return new

    # ------------------------------------------------------------------
    # Internals
    # ------------------------------------------------------------------

    def _check_compatible(self, other: FracMinHash) -> None:
        """Raise if two sketches have incompatible parameters."""
        if self.k != other.k:
            raise ValueError(
                f"k-mer size mismatch: {self.k} vs {other.k}"
            )
        if self.scale != other.scale:
            raise ValueError(
                f"scale mismatch: {self.scale} vs {other.scale}"
            )

    # ------------------------------------------------------------------
    # Dunder helpers
    # ------------------------------------------------------------------

    def __len__(self) -> int:
        return self.num_hashes

    def __bool__(self) -> bool:
        return self.num_hashes > 0

    def __repr__(self) -> str:
        return (
            f"FracMinHash(k={self.k}, scale={self.scale}, "
            f"num_hashes={self.num_hashes})"
        )

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, FracMinHash):
            return NotImplemented
        return (
            self.k == other.k
            and self.scale == other.scale
            and self._hash_set == other._hash_set
        )
