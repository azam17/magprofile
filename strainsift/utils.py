"""Utility functions for StrainSift.

k-mer extraction, hash functions, FASTQ I/O, and helper routines.
"""

from __future__ import annotations

import gzip
import hashlib
import struct
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path

import numpy as np

# MurmurHash3-like hash using built-in hashlib for portability.
# For production, mmh3 would be used; here we use a deterministic hash
# that gives consistent results without external C dependencies.

_HASH_SEED = 42
_MAX_HASH = (1 << 64) - 1


def murmurhash3_64(key: bytes, seed: int = _HASH_SEED) -> int:
    """Deterministic 64-bit hash of a byte string.

    Uses a fast xor-shift mixing of MD5 for portability.
    For production, replace with mmh3.hash64().
    """
    h = hashlib.md5(seed.to_bytes(8, "little") + key).digest()
    return struct.unpack("<Q", h[:8])[0]


def canonical_kmer(kmer: str) -> str:
    """Return the lexicographically smaller of a k-mer and its reverse complement."""
    comp = str.maketrans("ACGT", "TGCA")
    rc = kmer[::-1].translate(comp)
    return min(kmer, rc)


def extract_kmers(seq: str, k: int, canonical: bool = True) -> list[str]:
    """Extract all k-mers from a DNA sequence.

    Args:
        seq: DNA sequence (ACGT only, uppercase).
        k: k-mer size.
        canonical: If True, return canonical (min of forward/revcomp) k-mers.

    Returns:
        List of k-mer strings.
    """
    seq = seq.upper()
    kmers = []
    for i in range(len(seq) - k + 1):
        kmer = seq[i : i + k]
        if "N" in kmer:
            continue
        if canonical:
            kmer = canonical_kmer(kmer)
        kmers.append(kmer)
    return kmers


def hash_kmers(kmers: list[str]) -> list[int]:
    """Hash a list of k-mer strings to 64-bit integers."""
    return [murmurhash3_64(k.encode()) for k in kmers]


def kmer_containment(query_hashes: set[int], ref_hashes: set[int]) -> float:
    """Compute k-mer containment of query in reference.

    containment = |query ∩ ref| / |query|

    Returns 0.0 if query is empty.
    """
    if not query_hashes:
        return 0.0
    return len(query_hashes & ref_hashes) / len(query_hashes)


@dataclass
class FastqRecord:
    """A single FASTQ read."""

    read_id: str
    sequence: str
    quality: str

    @property
    def mean_quality(self) -> float:
        """Mean Phred quality score."""
        if not self.quality:
            return 0.0
        return np.mean([ord(c) - 33 for c in self.quality])


def parse_fastq(path: Path) -> Iterator[FastqRecord]:
    """Parse a FASTQ file (plain or gzipped).

    Yields FastqRecord objects.
    """
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        while True:
            header = fh.readline().strip()
            if not header:
                break
            seq = fh.readline().strip()
            fh.readline()  # + line
            qual = fh.readline().strip()
            read_id = header[1:].split()[0]  # strip @ and take first field
            yield FastqRecord(read_id=read_id, sequence=seq, quality=qual)


def write_fastq(records: list[FastqRecord], path: Path) -> None:
    """Write FASTQ records to a file (plain or gzipped)."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "wt") as fh:
        for rec in records:
            fh.write(f"@{rec.read_id}\n")
            fh.write(f"{rec.sequence}\n")
            fh.write("+\n")
            fh.write(f"{rec.quality}\n")


def bray_curtis(p: np.ndarray, q: np.ndarray) -> float:
    """Bray-Curtis dissimilarity between two abundance vectors."""
    num = np.sum(np.abs(p - q))
    denom = np.sum(p) + np.sum(q)
    if denom == 0:
        return 0.0
    return float(num / denom)


def l1_error(p: np.ndarray, q: np.ndarray) -> float:
    """L1 (Manhattan) distance between two vectors."""
    return float(np.sum(np.abs(p - q)))
