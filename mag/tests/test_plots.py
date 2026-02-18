"""Tests for mag.plots module."""

from pathlib import Path

import numpy as np
import pytest

from mag.diversity import rarefaction_curves
from mag.io import AbundanceTable
from mag.plots import plot_rarefaction


class TestPlotRarefaction:
    def test_creates_file(self, tmp_path):
        rng = np.random.default_rng(42)
        table = AbundanceTable(
            mag_ids=[f"M{i}" for i in range(8)],
            sample_ids=["s0", "s1", "s2"],
            abundances=rng.poisson(50, (8, 3)).astype(float),
        )
        result = rarefaction_curves(table, n_depths=5, n_iterations=3, seed=42)
        out = tmp_path / "rarefaction.pdf"
        plot_rarefaction(result, metric="richness", output=out)
        assert out.exists()
        assert out.stat().st_size > 0

    def test_shannon_metric(self, tmp_path):
        rng = np.random.default_rng(42)
        table = AbundanceTable(
            mag_ids=[f"M{i}" for i in range(8)],
            sample_ids=["s0", "s1", "s2"],
            abundances=rng.poisson(50, (8, 3)).astype(float),
        )
        result = rarefaction_curves(table, n_depths=5, n_iterations=3, seed=42)
        out = tmp_path / "rarefaction_shannon.pdf"
        plot_rarefaction(result, metric="shannon", output=out)
        assert out.exists()
        assert out.stat().st_size > 0
