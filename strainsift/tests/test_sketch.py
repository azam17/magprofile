"""Tests for the FracMinHash sketch implementation."""

import numpy as np
import pytest

from strainsift.sketch import FracMinHash
from strainsift.utils import extract_kmers, hash_kmers


class TestFracMinHash:
    def test_construction(self):
        sketch = FracMinHash(k=21, scale=1 / 1000)
        assert sketch.k == 21
        assert sketch.scale == 1 / 1000
        assert sketch.num_hashes == 0

    def test_invalid_k(self):
        with pytest.raises(ValueError):
            FracMinHash(k=0)

    def test_invalid_scale(self):
        with pytest.raises(ValueError):
            FracMinHash(k=21, scale=0)
        with pytest.raises(ValueError):
            FracMinHash(k=21, scale=1.5)

    def test_add_sequence(self):
        sketch = FracMinHash(k=5, scale=1 / 10)  # high scale for testing
        seq = "ACGTACGTACGTACGTACGT"  # 20 bp
        sketch.add_sequence(seq)
        assert sketch.num_hashes > 0

    def test_deterministic(self):
        s1 = FracMinHash(k=5, scale=1 / 10)
        s2 = FracMinHash(k=5, scale=1 / 10)
        seq = "ACGTACGTACGTACGTACGT"
        s1.add_sequence(seq)
        s2.add_sequence(seq)
        assert s1 == s2

    def test_hashes_sorted(self):
        sketch = FracMinHash(k=5, scale=1 / 10)
        sketch.add_sequence("ACGTACGTACGTACGTACGT")
        hashes = sketch.hashes
        assert np.all(hashes[:-1] <= hashes[1:])

    def test_containment_self(self):
        sketch = FracMinHash(k=5, scale=1 / 10)
        sketch.add_sequence("ACGTACGTACGTACGTACGT")
        if sketch.num_hashes > 0:
            assert sketch.containment(sketch) == 1.0

    def test_containment_empty(self):
        s1 = FracMinHash(k=5, scale=1 / 10)
        s2 = FracMinHash(k=5, scale=1 / 10)
        s2.add_sequence("ACGTACGTACGTACGTACGT")
        assert s1.containment(s2) == 0.0

    def test_containment_subset(self):
        # A short sequence should be contained in a longer one that includes it
        s1 = FracMinHash(k=5, scale=0.5)  # very high scale
        s2 = FracMinHash(k=5, scale=0.5)
        short = "ACGTACGTACGT"
        long = "ACGTACGTACGTACGTACGT"
        s1.add_sequence(short)
        s2.add_sequence(long)
        if s1.num_hashes > 0:
            c = s1.containment(s2)
            assert c >= 0.5  # short is mostly contained in long

    def test_jaccard_identical(self):
        s1 = FracMinHash(k=5, scale=0.5)
        seq = "ACGTACGTACGTACGTACGT"
        s1.add_sequence(seq)
        if s1.num_hashes > 0:
            assert s1.jaccard(s1) == 1.0

    def test_jaccard_disjoint(self):
        s1 = FracMinHash(k=5, scale=0.5)
        s2 = FracMinHash(k=5, scale=0.5)
        s1.add_sequence("AAAAAAAAAAA")
        s2.add_sequence("CCCCCCCCCCC")
        j = s1.jaccard(s2)
        # Should have very low Jaccard for sequences with no shared k-mers
        assert j < 0.5

    def test_merge(self):
        s1 = FracMinHash(k=5, scale=0.5)
        s2 = FracMinHash(k=5, scale=0.5)
        s1.add_sequence("ACGTACGT")
        s2.add_sequence("TGCATGCA")
        n1 = s1.num_hashes
        s1.merge(s2)
        assert s1.num_hashes >= n1

    def test_copy(self):
        s1 = FracMinHash(k=5, scale=0.5)
        s1.add_sequence("ACGTACGT")
        s2 = s1.copy()
        assert s1 == s2
        s2.add_sequence("TGCATGCA")
        assert s1 != s2  # copy is independent

    def test_incompatible_k(self):
        s1 = FracMinHash(k=5, scale=0.5)
        s2 = FracMinHash(k=7, scale=0.5)
        s1.add_sequence("ACGTACGT")
        s2.add_sequence("ACGTACGTACGT")
        with pytest.raises(ValueError, match="k-mer size mismatch"):
            s1.containment(s2)

    def test_incompatible_scale(self):
        s1 = FracMinHash(k=5, scale=0.5)
        s2 = FracMinHash(k=5, scale=0.1)
        s1.add_sequence("ACGTACGT")
        s2.add_sequence("ACGTACGT")
        with pytest.raises(ValueError, match="scale mismatch"):
            s1.containment(s2)

    def test_len(self):
        sketch = FracMinHash(k=5, scale=0.5)
        assert len(sketch) == 0
        sketch.add_sequence("ACGTACGTACGT")
        assert len(sketch) == sketch.num_hashes

    def test_bool(self):
        sketch = FracMinHash(k=5, scale=0.5)
        assert not sketch
        sketch.add_sequence("ACGTACGTACGT")
        if sketch.num_hashes > 0:
            assert sketch

    def test_threshold_filtering(self):
        # With very low scale, few hashes should be retained
        sketch_low = FracMinHash(k=5, scale=1 / 100000)
        # With high scale, many hashes should be retained
        sketch_high = FracMinHash(k=5, scale=0.9)
        seq = "ACGTACGTACGTACGTACGT"
        sketch_low.add_sequence(seq)
        sketch_high.add_sequence(seq)
        assert sketch_high.num_hashes >= sketch_low.num_hashes
