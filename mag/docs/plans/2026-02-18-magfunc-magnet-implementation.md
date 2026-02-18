# magfunc + magnet Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add functional profiling (magfunc) and network analysis (magnet) modules to the existing mag/ package, with DRAM-based input, CAZyme profiling, proportionality networks, and full report generation.

**Architecture:** Both tools are modules within mag/, sharing existing data types (AbundanceTable, TaxonomyTable, SampleMetadata) and utilities (CLR transform, BH-FDR). Each has its own I/O, analysis, and report modules. CLI adds `func-report` and `net-report` subcommands.

**Tech Stack:** Python 3.10+, numpy, scipy, matplotlib, seaborn, click, networkx (new dep for magnet)

---

## Task 1: FunctionalTable data type and DRAM parsing

**Files:**
- Create: `mag/func_io.py`
- Test: `mag/tests/test_func_io.py`

**Step 1: Write the failing tests**

```python
# mag/tests/test_func_io.py
"""Tests for DRAM parsing and FunctionalTable."""

from __future__ import annotations

import textwrap
from pathlib import Path

import numpy as np
import pytest

from mag.func_io import (
    DRAMAnnotation,
    FunctionalTable,
    build_functional_table,
    load_dram_annotations,
    PGPR_MARKERS,
)
from mag.io import AbundanceTable


@pytest.fixture
def dram_tsv(tmp_path: Path) -> Path:
    """Create a minimal DRAM annotations.tsv."""
    content = textwrap.dedent("""\
        \tgene_id\tgene_description\tmodule\theader\tsubheader\tkegg_id\tkegg_hit\tpfam_hits\tcazy_hits
        MAG_001_00001\tMAG_001_00001\tnifH protein\tM00175\tNitrogen metabolism\tNitrogen fixation\tK02588\tnitrogenase iron protein\tPF00142\t
        MAG_001_00002\tMAG_001_00002\tGH5 cellulase\t\tCarbohydrate metabolism\t\t\t\t\tGH5
        MAG_002_00001\tMAG_002_00001\tnifD protein\tM00175\tNitrogen metabolism\tNitrogen fixation\tK02586\tnitrogenase molybdenum-iron protein\t\t
        MAG_002_00002\tMAG_002_00002\tGT2 glycosyl\t\t\t\t\t\t\tGT2
        MAG_002_00003\tMAG_002_00003\tAA3 oxidase\t\t\t\tK00111\tsome enzyme\t\tAA3
    """)
    p = tmp_path / "annotations.tsv"
    p.write_text(content)
    return p


class TestLoadDRAMAnnotations:
    def test_parse_basic(self, dram_tsv: Path) -> None:
        annots = load_dram_annotations(dram_tsv)
        assert len(annots) == 5
        assert annots[0].mag_id == "MAG_001"
        assert annots[0].ko_id == "K02588"
        assert annots[0].cazy_id is None

    def test_cazy_parsed(self, dram_tsv: Path) -> None:
        annots = load_dram_annotations(dram_tsv)
        cazy_annots = [a for a in annots if a.cazy_id]
        assert len(cazy_annots) == 3
        assert {a.cazy_id for a in cazy_annots} == {"GH5", "GT2", "AA3"}

    def test_mag_id_from_gene_id(self, dram_tsv: Path) -> None:
        annots = load_dram_annotations(dram_tsv)
        mag_ids = {a.mag_id for a in annots}
        assert mag_ids == {"MAG_001", "MAG_002"}


class TestBuildFunctionalTable:
    def test_ko_table(self, dram_tsv: Path) -> None:
        annots = load_dram_annotations(dram_tsv)
        ft = build_functional_table(annots, "ko")
        assert ft.function_type == "ko"
        assert ft.values.shape[0] > 0  # some KO IDs
        assert ft.values.shape[1] == 2  # 2 MAGs
        # K02588 should be present in MAG_001
        ko_idx = ft.function_ids.index("K02588")
        mag_idx = ft.mag_ids.index("MAG_001")
        assert ft.values[ko_idx, mag_idx] == 1

    def test_cazy_table(self, dram_tsv: Path) -> None:
        annots = load_dram_annotations(dram_tsv)
        ft = build_functional_table(annots, "cazy")
        assert ft.function_type == "cazy"
        assert "GH5" in ft.function_ids
        assert "GT2" in ft.function_ids
        assert "AA3" in ft.function_ids

    def test_pgpr_table(self, dram_tsv: Path) -> None:
        annots = load_dram_annotations(dram_tsv)
        ft = build_functional_table(annots, "pgpr")
        assert ft.function_type == "pgpr"
        # nifH (K02588) should be in PGPR markers
        assert "nifH" in ft.function_ids

    def test_unknown_type_raises(self, dram_tsv: Path) -> None:
        annots = load_dram_annotations(dram_tsv)
        with pytest.raises(ValueError, match="Unknown function_type"):
            build_functional_table(annots, "invalid")


class TestFunctionalTableToSampleAbundance:
    def test_matrix_multiplication(self) -> None:
        ft = FunctionalTable(
            function_ids=["K1", "K2"],
            mag_ids=["M1", "M2"],
            values=np.array([[1, 0], [0, 1]], dtype=np.float64),
            function_type="ko",
        )
        abundance = AbundanceTable(
            mag_ids=["M1", "M2"],
            sample_ids=["s1", "s2"],
            abundances=np.array([[100, 200], [300, 400]], dtype=np.float64),
        )
        result = ft.to_sample_abundance(abundance)
        assert result.mag_ids == ["K1", "K2"]
        assert result.sample_ids == ["s1", "s2"]
        np.testing.assert_array_equal(result.abundances[0], [100, 200])  # K1 = M1
        np.testing.assert_array_equal(result.abundances[1], [300, 400])  # K2 = M2

    def test_shared_mag_subset(self) -> None:
        """FunctionalTable and AbundanceTable may have different MAG sets."""
        ft = FunctionalTable(
            function_ids=["K1"],
            mag_ids=["M1", "M3"],  # M3 not in abundance
            values=np.array([[1, 1]], dtype=np.float64),
            function_type="ko",
        )
        abundance = AbundanceTable(
            mag_ids=["M1", "M2"],  # M2 not in func table
            sample_ids=["s1"],
            abundances=np.array([[100], [200]], dtype=np.float64),
        )
        result = ft.to_sample_abundance(abundance)
        # Only M1 is shared
        np.testing.assert_array_equal(result.abundances[0], [100])


class TestPGPRMarkers:
    def test_markers_exist(self) -> None:
        assert "nifH" in PGPR_MARKERS
        assert "pqqC" in PGPR_MARKERS
        assert "acdS" in PGPR_MARKERS
        for name, ko_ids in PGPR_MARKERS.items():
            assert isinstance(ko_ids, list)
            assert all(k.startswith("K") for k in ko_ids)
```

**Step 2: Run tests to verify they fail**

Run: `.venv/bin/python -m pytest mag/tests/test_func_io.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'mag.func_io'`

**Step 3: Implement func_io.py**

```python
# mag/func_io.py
"""DRAM annotation parsing and functional data types."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

from .io import AbundanceTable


# Curated PGPR marker genes (trait name -> list of KEGG Orthology IDs)
PGPR_MARKERS: dict[str, list[str]] = {
    "nifH": ["K02588"],
    "nifD": ["K02586"],
    "nifK": ["K02591"],
    "pqqC": ["K06137"],
    "gcd": ["K00111"],
    "acdS": ["K01505"],
    "ipdC": ["K04103"],
    "acsA": ["K14187"],
    "budB": ["K01652"],
    "entA": ["K02361"],
    "entB": ["K01252"],
    "entC": ["K02362"],
    "phzF": ["K18000"],
}


@dataclass
class DRAMAnnotation:
    """Single gene annotation from DRAM."""

    gene_id: str
    mag_id: str
    ko_id: str | None = None
    kegg_hit: str | None = None
    cazy_id: str | None = None
    pfam_id: str | None = None


@dataclass(frozen=True)
class FunctionalTable:
    """Function-by-MAG matrix.

    ``values`` has shape (n_functions, n_mags).  Entries are binary (0/1)
    for KO/PGPR tables, integer counts for CAZyme tables, or completeness
    scores (0.0–1.0) for module tables.
    """

    function_ids: list[str]
    mag_ids: list[str]
    values: np.ndarray  # shape (n_functions, n_mags)
    function_type: str  # "ko" | "module" | "cazy" | "pgpr"

    def to_sample_abundance(self, abundance: AbundanceTable) -> AbundanceTable:
        """Weight functions by MAG abundances to get function x sample matrix.

        Only MAGs present in both this table and the AbundanceTable are used.
        Returns an AbundanceTable with function_ids as mag_ids.
        """
        # Find shared MAGs and align indices
        shared = [m for m in self.mag_ids if m in set(abundance.mag_ids)]
        if not shared:
            return AbundanceTable(
                mag_ids=list(self.function_ids),
                sample_ids=list(abundance.sample_ids),
                abundances=np.zeros(
                    (len(self.function_ids), abundance.n_samples), dtype=np.float64
                ),
            )

        func_idx = [self.mag_ids.index(m) for m in shared]
        abund_idx = [abundance.mag_ids.index(m) for m in shared]

        func_sub = self.values[:, func_idx]  # (n_func, n_shared)
        abund_sub = abundance.abundances[abund_idx, :]  # (n_shared, n_samples)

        result = func_sub @ abund_sub  # (n_func, n_samples)

        return AbundanceTable(
            mag_ids=list(self.function_ids),
            sample_ids=list(abundance.sample_ids),
            abundances=result.astype(np.float64),
        )


def load_dram_annotations(path: str | Path) -> list[DRAMAnnotation]:
    """Parse DRAM annotations.tsv into a list of DRAMAnnotation.

    MAG ID is derived from gene_id by dropping the last ``_NNNNN`` suffix.
    """
    path = Path(path)
    annotations: list[DRAMAnnotation] = []

    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            gene_id = row.get("gene_id", row.get("", "")).strip()
            if not gene_id:
                # First column might be unnamed (index column)
                gene_id = list(row.values())[0].strip()
            if not gene_id:
                continue

            # MAG ID = gene_id minus last _NNNNN segment
            parts = gene_id.rsplit("_", 1)
            mag_id = parts[0] if len(parts) > 1 else gene_id

            ko_id = row.get("kegg_id", "").strip() or None
            kegg_hit = row.get("kegg_hit", "").strip() or None
            cazy_raw = row.get("cazy_hits", row.get("cazy_id", "")).strip() or None
            pfam_raw = row.get("pfam_hits", row.get("pfam_id", "")).strip() or None

            annotations.append(
                DRAMAnnotation(
                    gene_id=gene_id,
                    mag_id=mag_id,
                    ko_id=ko_id,
                    kegg_hit=kegg_hit,
                    cazy_id=cazy_raw,
                    pfam_id=pfam_raw,
                )
            )

    return annotations


def build_functional_table(
    annotations: Sequence[DRAMAnnotation],
    function_type: str,
) -> FunctionalTable:
    """Build a function x MAG matrix from DRAM annotations.

    Args:
        annotations: Parsed DRAM annotations.
        function_type: One of "ko", "cazy", or "pgpr".

    Returns:
        FunctionalTable with binary (ko/pgpr) or count (cazy) values.
    """
    if function_type == "ko":
        return _build_ko_table(annotations)
    elif function_type == "cazy":
        return _build_cazy_table(annotations)
    elif function_type == "pgpr":
        return _build_pgpr_table(annotations)
    else:
        raise ValueError(f"Unknown function_type: {function_type!r}")


def _build_ko_table(annotations: Sequence[DRAMAnnotation]) -> FunctionalTable:
    """Binary KO presence/absence per MAG."""
    ko_mags: dict[str, set[str]] = {}
    all_mags: set[str] = set()
    for a in annotations:
        all_mags.add(a.mag_id)
        if a.ko_id:
            ko_mags.setdefault(a.ko_id, set()).add(a.mag_id)

    ko_ids = sorted(ko_mags.keys())
    mag_ids = sorted(all_mags)
    mag_idx = {m: i for i, m in enumerate(mag_ids)}

    values = np.zeros((len(ko_ids), len(mag_ids)), dtype=np.float64)
    for ki, ko in enumerate(ko_ids):
        for m in ko_mags[ko]:
            values[ki, mag_idx[m]] = 1

    return FunctionalTable(
        function_ids=ko_ids, mag_ids=mag_ids, values=values, function_type="ko"
    )


def _build_cazy_table(annotations: Sequence[DRAMAnnotation]) -> FunctionalTable:
    """CAZyme family counts per MAG."""
    cazy_mags: dict[str, dict[str, int]] = {}
    all_mags: set[str] = set()
    for a in annotations:
        all_mags.add(a.mag_id)
        if a.cazy_id:
            cazy_mags.setdefault(a.cazy_id, {})
            cazy_mags[a.cazy_id][a.mag_id] = (
                cazy_mags[a.cazy_id].get(a.mag_id, 0) + 1
            )

    cazy_ids = sorted(cazy_mags.keys())
    mag_ids = sorted(all_mags)
    mag_idx = {m: i for i, m in enumerate(mag_ids)}

    values = np.zeros((len(cazy_ids), len(mag_ids)), dtype=np.float64)
    for ci, cazy in enumerate(cazy_ids):
        for m, count in cazy_mags[cazy].items():
            values[ci, mag_idx[m]] = count

    return FunctionalTable(
        function_ids=cazy_ids, mag_ids=mag_ids, values=values, function_type="cazy"
    )


def _build_pgpr_table(annotations: Sequence[DRAMAnnotation]) -> FunctionalTable:
    """PGPR trait binary presence per MAG based on curated KO markers."""
    # Invert PGPR_MARKERS: KO -> trait name
    ko_to_trait: dict[str, str] = {}
    for trait, ko_list in PGPR_MARKERS.items():
        for ko in ko_list:
            ko_to_trait[ko] = trait

    trait_mags: dict[str, set[str]] = {t: set() for t in PGPR_MARKERS}
    all_mags: set[str] = set()
    for a in annotations:
        all_mags.add(a.mag_id)
        if a.ko_id and a.ko_id in ko_to_trait:
            trait_mags[ko_to_trait[a.ko_id]].add(a.mag_id)

    trait_ids = sorted(trait_mags.keys())
    mag_ids = sorted(all_mags)
    mag_idx = {m: i for i, m in enumerate(mag_ids)}

    values = np.zeros((len(trait_ids), len(mag_ids)), dtype=np.float64)
    for ti, trait in enumerate(trait_ids):
        for m in trait_mags[trait]:
            values[ti, mag_idx[m]] = 1

    return FunctionalTable(
        function_ids=trait_ids, mag_ids=mag_ids, values=values, function_type="pgpr"
    )
```

**Step 4: Run tests to verify they pass**

Run: `.venv/bin/python -m pytest mag/tests/test_func_io.py -v`
Expected: All PASS

**Step 5: Commit**

```bash
git add mag/func_io.py mag/tests/test_func_io.py
git commit -m "feat: add FunctionalTable data type and DRAM parsing (magfunc)"
```

---

## Task 2: Functional profiling analysis module

**Files:**
- Create: `mag/func_profile.py`
- Test: `mag/tests/test_func_profile.py`

**Step 1: Write the failing tests**

```python
# mag/tests/test_func_profile.py
"""Tests for functional profiling analysis."""

from __future__ import annotations

import numpy as np
import pytest

from mag.func_io import FunctionalTable
from mag.func_profile import (
    pathway_abundance,
    differential_pathway,
    functional_redundancy,
    cazyme_summary,
)
from mag.io import AbundanceTable, SampleMetadata
from mag.tests.fixtures import generate_synthetic_abundance_table


@pytest.fixture
def synth():
    return generate_synthetic_abundance_table()


@pytest.fixture
def func_table(synth) -> FunctionalTable:
    """Create a FunctionalTable aligned with synthetic data."""
    table = synth[0]
    rng = np.random.default_rng(99)
    n_functions = 10
    values = rng.integers(0, 2, size=(n_functions, table.n_mags)).astype(np.float64)
    return FunctionalTable(
        function_ids=[f"K{i:05d}" for i in range(n_functions)],
        mag_ids=list(table.mag_ids),
        values=values,
        function_type="ko",
    )


class TestPathwayAbundance:
    def test_output_shape(self, func_table, synth) -> None:
        table, meta, _ = synth
        result = pathway_abundance(func_table, table)
        assert result.n_mags == len(func_table.function_ids)
        assert result.n_samples == table.n_samples

    def test_zero_function_zero_abundance(self, synth) -> None:
        table = synth[0]
        ft = FunctionalTable(
            function_ids=["K1"],
            mag_ids=list(table.mag_ids),
            values=np.zeros((1, table.n_mags)),
            function_type="ko",
        )
        result = pathway_abundance(ft, table)
        assert result.abundances.sum() == 0.0


class TestDifferentialPathway:
    def test_returns_result(self, func_table, synth) -> None:
        table, meta, _ = synth
        result = differential_pathway(
            func_table, table, meta, "compartment", "topsoil", "bulk"
        )
        assert len(result.mag_ids) == len(func_table.function_ids)
        assert len(result.q_values) == len(func_table.function_ids)
        assert all(0 <= q <= 1 for q in result.q_values)


class TestFunctionalRedundancy:
    def test_returns_per_group(self, func_table, synth) -> None:
        table, meta, _ = synth
        result = functional_redundancy(func_table, table, meta, "compartment")
        groups = meta.get_groups("compartment")
        for g in groups:
            assert g in result
        for g, funcs in result.items():
            for func_id, info in funcs.items():
                assert "n_carriers" in info
                assert "shannon" in info
                assert info["n_carriers"] >= 0


class TestCAZymeSummary:
    def test_aggregates_by_class(self) -> None:
        ft = FunctionalTable(
            function_ids=["GH5", "GH13", "GT2", "AA3", "PL1", "CE1"],
            mag_ids=["M1", "M2"],
            values=np.array([
                [2, 0], [1, 3], [0, 1], [1, 1], [0, 2], [1, 0],
            ], dtype=np.float64),
            function_type="cazy",
        )
        result = cazyme_summary(ft)
        assert "GH" in result
        assert "GT" in result
        assert "AA" in result
        assert "PL" in result
        assert "CE" in result
        # GH total across both MAGs: 2+0+1+3 = 6
        assert result["GH"]["total_genes"] == 6
```

**Step 2: Run tests to verify they fail**

Run: `.venv/bin/python -m pytest mag/tests/test_func_profile.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'mag.func_profile'`

**Step 3: Implement func_profile.py**

```python
# mag/func_profile.py
"""Functional profiling analysis for MAG communities."""

from __future__ import annotations

import numpy as np

from .differential import DifferentialAbundanceResult, differential_abundance
from .func_io import FunctionalTable
from .io import AbundanceTable, SampleMetadata


def pathway_abundance(
    func_table: FunctionalTable,
    abundance: AbundanceTable,
) -> AbundanceTable:
    """Compute function x sample abundance matrix.

    Multiplies function presence (or completeness) by MAG abundances.
    """
    return func_table.to_sample_abundance(abundance)


def differential_pathway(
    func_table: FunctionalTable,
    abundance: AbundanceTable,
    metadata: SampleMetadata,
    grouping_var: str,
    group1: str,
    group2: str,
) -> DifferentialAbundanceResult:
    """Differential pathway abundance between two groups.

    Computes function x sample abundance, then runs standard CLR +
    Welch's t-test + BH-FDR pipeline.
    """
    func_abundance = pathway_abundance(func_table, abundance)
    return differential_abundance(func_abundance, metadata, grouping_var, group1, group2)


def functional_redundancy(
    func_table: FunctionalTable,
    abundance: AbundanceTable,
    metadata: SampleMetadata,
    grouping_var: str,
) -> dict[str, dict[str, dict[str, float]]]:
    """Compute functional redundancy per function per group.

    For each function and group, counts how many MAGs carry it and
    computes the Shannon diversity of their abundance distribution
    (higher = more evenly distributed = more redundant).

    Returns:
        dict[group -> dict[function_id -> {"n_carriers": int, "shannon": float}]]
    """
    groups = metadata.get_groups(grouping_var)
    result: dict[str, dict[str, dict[str, float]]] = {}

    for group_name, sample_ids in groups.items():
        sub = abundance.subset_samples(sample_ids)
        group_result: dict[str, dict[str, float]] = {}

        for fi, func_id in enumerate(func_table.function_ids):
            # MAGs carrying this function
            carrier_mask = func_table.values[fi] > 0
            # Align to abundance table
            shared = [
                (mi_f, sub.mag_ids.index(m))
                for mi_f, m in enumerate(func_table.mag_ids)
                if carrier_mask[mi_f] and m in set(sub.mag_ids)
            ]

            n_carriers = len(shared)
            if n_carriers == 0:
                group_result[func_id] = {"n_carriers": 0, "shannon": 0.0}
                continue

            # Mean abundance per carrier across group samples
            carrier_abundances = np.array(
                [sub.abundances[idx_a].mean() for _, idx_a in shared]
            )
            total = carrier_abundances.sum()
            if total == 0:
                group_result[func_id] = {
                    "n_carriers": n_carriers,
                    "shannon": 0.0,
                }
                continue

            props = carrier_abundances / total
            props = props[props > 0]
            shannon = float(-np.sum(props * np.log(props)))

            group_result[func_id] = {
                "n_carriers": n_carriers,
                "shannon": shannon,
            }

        result[group_name] = group_result

    return result


def cazyme_summary(
    func_table: FunctionalTable,
) -> dict[str, dict[str, float]]:
    """Aggregate CAZyme families by class (GH, GT, PL, CE, AA, CBM).

    Returns dict[class -> {"total_genes": int, "n_families": int, "n_mags": int}].
    """
    class_stats: dict[str, dict[str, float]] = {}

    for fi, fam_id in enumerate(func_table.function_ids):
        # Extract class prefix (e.g. "GH" from "GH5", "GT" from "GT2")
        prefix = ""
        for ch in fam_id:
            if ch.isalpha():
                prefix += ch
            else:
                break
        if not prefix:
            prefix = "Other"

        row = func_table.values[fi]
        total_genes = float(row.sum())
        n_mags_with = int((row > 0).sum())

        if prefix not in class_stats:
            class_stats[prefix] = {"total_genes": 0, "n_families": 0, "n_mags": 0}

        class_stats[prefix]["total_genes"] += total_genes
        class_stats[prefix]["n_families"] += 1
        class_stats[prefix]["n_mags"] = max(
            class_stats[prefix]["n_mags"], n_mags_with
        )

    return class_stats
```

**Step 4: Run tests to verify they pass**

Run: `.venv/bin/python -m pytest mag/tests/test_func_profile.py -v`
Expected: All PASS

**Step 5: Commit**

```bash
git add mag/func_profile.py mag/tests/test_func_profile.py
git commit -m "feat: add functional profiling analysis module (magfunc)"
```

---

## Task 3: Functional report generation (CSV + plots + HTML)

**Files:**
- Create: `mag/func_report.py`
- Modify: `mag/plots.py` — add `plot_pathway_heatmap`, `plot_pgpr_matrix`, `plot_cazyme_bars`
- Modify: `mag/cli.py` — add `func-report` command
- Test: `mag/tests/test_func_report.py`

**Step 1: Write the failing tests**

```python
# mag/tests/test_func_report.py
"""Tests for functional report generation."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from mag.func_io import DRAMAnnotation, FunctionalTable, build_functional_table
from mag.func_report import generate_func_report
from mag.tests.fixtures import generate_synthetic_abundance_table


@pytest.fixture
def synth():
    return generate_synthetic_abundance_table()


@pytest.fixture
def annotations(synth) -> list[DRAMAnnotation]:
    """Generate synthetic DRAM annotations for test MAGs."""
    table = synth[0]
    rng = np.random.default_rng(77)
    annots: list[DRAMAnnotation] = []
    ko_pool = ["K02588", "K02586", "K02591", "K06137", "K00111", "K01505", "K99999"]
    cazy_pool = ["GH5", "GH13", "GT2", "AA3", "PL1", "CE1", None, None]

    for mag_id in table.mag_ids:
        n_genes = rng.integers(5, 15)
        for g in range(n_genes):
            annots.append(
                DRAMAnnotation(
                    gene_id=f"{mag_id}_{g:05d}",
                    mag_id=mag_id,
                    ko_id=rng.choice(ko_pool) if rng.random() > 0.3 else None,
                    cazy_id=rng.choice(cazy_pool),
                )
            )
    return annots


class TestGenerateFuncReport:
    def test_creates_output_files(self, synth, annotations, tmp_path) -> None:
        table, meta, tax = synth
        generate_func_report(
            abundance=table,
            annotations=annotations,
            taxonomy=tax,
            metadata=meta,
            grouping_var="compartment",
            output_dir=tmp_path,
        )
        assert (tmp_path / "pathway_abundance.csv").exists()
        assert (tmp_path / "pgpr_traits.csv").exists()
        assert (tmp_path / "cazyme_summary.csv").exists()
        assert (tmp_path / "functional_redundancy.csv").exists()
        assert (tmp_path / "func_report.html").exists()

    def test_differential_csvs_created(self, synth, annotations, tmp_path) -> None:
        table, meta, tax = synth
        generate_func_report(
            abundance=table,
            annotations=annotations,
            taxonomy=tax,
            metadata=meta,
            grouping_var="compartment",
            output_dir=tmp_path,
        )
        # Should have pairwise pathway differential CSVs
        diff_files = list(tmp_path.glob("pathway_differential_*.csv"))
        assert len(diff_files) >= 1

    def test_plots_created(self, synth, annotations, tmp_path) -> None:
        table, meta, tax = synth
        generate_func_report(
            abundance=table,
            annotations=annotations,
            taxonomy=tax,
            metadata=meta,
            grouping_var="compartment",
            output_dir=tmp_path,
        )
        assert (tmp_path / "pgpr_traits.pdf").exists()
        assert (tmp_path / "cazyme_bars.pdf").exists()

    def test_no_taxonomy_ok(self, synth, annotations, tmp_path) -> None:
        table, meta, _ = synth
        generate_func_report(
            abundance=table,
            annotations=annotations,
            taxonomy=None,
            metadata=meta,
            grouping_var="compartment",
            output_dir=tmp_path,
        )
        assert (tmp_path / "func_report.html").exists()
```

**Step 2: Run tests to verify they fail**

Run: `.venv/bin/python -m pytest mag/tests/test_func_report.py -v`
Expected: FAIL

**Step 3: Implement func_report.py and add plots**

Create `mag/func_report.py`:

```python
# mag/func_report.py
"""Report generation for functional profiling."""

from __future__ import annotations

import base64
import csv
import logging
from itertools import combinations
from pathlib import Path

import numpy as np

from .differential import differential_abundance
from .func_io import DRAMAnnotation, FunctionalTable, build_functional_table
from .func_profile import (
    cazyme_summary,
    differential_pathway,
    functional_redundancy,
    pathway_abundance,
)
from .io import AbundanceTable, SampleMetadata, TaxonomyTable

logger = logging.getLogger(__name__)


def generate_func_report(
    abundance: AbundanceTable,
    annotations: list[DRAMAnnotation],
    taxonomy: TaxonomyTable | None,
    metadata: SampleMetadata,
    grouping_var: str,
    output_dir: str | Path,
) -> None:
    """Run functional profiling pipeline and write results."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    # Build tables
    ko_table = build_functional_table(annotations, "ko")
    cazy_table = build_functional_table(annotations, "cazy")
    pgpr_table = build_functional_table(annotations, "pgpr")

    # Pathway abundance
    pa = pathway_abundance(ko_table, abundance)
    _write_abundance_csv(pa, out / "pathway_abundance.csv")

    # PGPR traits
    _write_functional_table_csv(pgpr_table, out / "pgpr_traits.csv")

    # CAZyme summary
    cazy_stats = cazyme_summary(cazy_table)
    _write_cazyme_csv(cazy_stats, out / "cazyme_summary.csv")

    # Functional redundancy
    redundancy = functional_redundancy(ko_table, abundance, metadata, grouping_var)
    _write_redundancy_csv(redundancy, out / "functional_redundancy.csv")

    # Differential pathway abundance for all group pairs
    groups = metadata.get_groups(grouping_var)
    group_names = sorted(groups.keys())
    diff_results = {}

    for g1, g2 in combinations(group_names, 2):
        try:
            result = differential_pathway(
                ko_table, abundance, metadata, grouping_var, g1, g2,
            )
            diff_results[(g1, g2)] = result
            _write_diff_csv(result, g1, g2, out / f"pathway_differential_{g1}_vs_{g2}.csv")
        except Exception as e:
            logger.warning("Pathway differential %s vs %s failed: %s", g1, g2, e)

    # Plots
    from . import plots

    plots.plot_pgpr_matrix(pgpr_table, taxonomy, output=out / "pgpr_traits.pdf")
    plots.plot_cazyme_bars(cazy_table, metadata, grouping_var, abundance, output=out / "cazyme_bars.pdf")

    # HTML report
    _write_func_html(
        out,
        ko_table=ko_table,
        cazy_table=cazy_table,
        pgpr_table=pgpr_table,
        cazy_stats=cazy_stats,
        redundancy=redundancy,
        diff_results=diff_results,
        taxonomy=taxonomy,
    )


def _write_abundance_csv(table: AbundanceTable, path: Path) -> None:
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["function_id"] + table.sample_ids)
        for i, fid in enumerate(table.mag_ids):
            w.writerow([fid] + [f"{v:.4f}" for v in table.abundances[i]])


def _write_functional_table_csv(ft: FunctionalTable, path: Path) -> None:
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["function_id"] + ft.mag_ids)
        for i, fid in enumerate(ft.function_ids):
            w.writerow([fid] + [f"{v:g}" for v in ft.values[i]])


def _write_cazyme_csv(stats: dict, path: Path) -> None:
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["cazy_class", "total_genes", "n_families", "n_mags"])
        for cls in sorted(stats.keys()):
            s = stats[cls]
            w.writerow([cls, int(s["total_genes"]), int(s["n_families"]), int(s["n_mags"])])


def _write_redundancy_csv(redundancy: dict, path: Path) -> None:
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["group", "function_id", "n_carriers", "shannon"])
        for group, funcs in sorted(redundancy.items()):
            for func_id, info in sorted(funcs.items()):
                w.writerow([group, func_id, int(info["n_carriers"]), f"{info['shannon']:.4f}"])


def _write_diff_csv(result, g1: str, g2: str, path: Path) -> None:
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["function_id", "CLR_diff", "p_value", "q_value", "cohens_d"])
        for i, fid in enumerate(result.mag_ids):
            w.writerow([
                fid,
                f"{result.log_fold_changes[i]:.4f}",
                f"{result.p_values[i]:.6f}",
                f"{result.q_values[i]:.6f}",
                f"{result.effect_sizes[i]:.4f}",
            ])


def _embed_figure(path: Path) -> str:
    if not path.exists():
        return ""
    data = path.read_bytes()
    b64 = base64.b64encode(data).decode("ascii")
    return f'<img src="data:image/png;base64,{b64}" style="max-width:100%">'


def _write_func_html(
    out: Path,
    ko_table: FunctionalTable,
    cazy_table: FunctionalTable,
    pgpr_table: FunctionalTable,
    cazy_stats: dict,
    redundancy: dict,
    diff_results: dict,
    taxonomy: TaxonomyTable | None,
) -> None:
    """Generate self-contained HTML report for functional profiling."""
    html_parts = [
        "<!DOCTYPE html><html><head>",
        "<meta charset='utf-8'>",
        "<title>magfunc — Functional Profiling Report</title>",
        "<style>",
        "body{font-family:sans-serif;max-width:1200px;margin:0 auto;padding:20px;line-height:1.6}",
        "h1{border-bottom:3px solid #2c3e50;padding-bottom:10px}",
        "h2{color:#2c3e50;cursor:pointer;border-bottom:1px solid #ddd;padding:8px 0}",
        "table{border-collapse:collapse;width:100%;margin:10px 0}",
        "th,td{border:1px solid #ddd;padding:6px 10px;text-align:left}",
        "th{background:#f5f5f5;cursor:pointer}",
        "tr:nth-child(even){background:#fafafa}",
        ".section{margin-bottom:30px}",
        ".sig{color:#e74c3c;font-weight:bold}",
        "</style></head><body>",
        "<h1>magfunc — Functional Profiling Report</h1>",
    ]

    # Summary
    html_parts.append("<div class='section'>")
    html_parts.append("<h2>Summary</h2>")
    html_parts.append(f"<p>KO functions detected: <strong>{len(ko_table.function_ids)}</strong></p>")
    html_parts.append(f"<p>CAZyme families detected: <strong>{len(cazy_table.function_ids)}</strong></p>")
    html_parts.append(f"<p>PGPR traits screened: <strong>{len(pgpr_table.function_ids)}</strong></p>")
    html_parts.append(f"<p>MAGs annotated: <strong>{len(ko_table.mag_ids)}</strong></p>")
    html_parts.append("</div>")

    # PGPR traits
    html_parts.append("<div class='section'>")
    html_parts.append("<h2>Plant Growth-Promoting Traits (PGPR)</h2>")
    pgpr_fig = _embed_figure(out / "pgpr_traits.pdf")
    if pgpr_fig:
        html_parts.append(pgpr_fig)
    n_with_trait = int((pgpr_table.values.sum(axis=0) > 0).sum())
    html_parts.append(f"<p>{n_with_trait} of {len(pgpr_table.mag_ids)} MAGs carry at least one PGPR trait.</p>")
    html_parts.append("</div>")

    # CAZymes
    html_parts.append("<div class='section'>")
    html_parts.append("<h2>Carbohydrate-Active Enzymes (CAZymes)</h2>")
    cazy_fig = _embed_figure(out / "cazyme_bars.pdf")
    if cazy_fig:
        html_parts.append(cazy_fig)
    html_parts.append("<table><tr><th>Class</th><th>Families</th><th>Total Genes</th><th>MAGs</th></tr>")
    for cls in sorted(cazy_stats.keys()):
        s = cazy_stats[cls]
        html_parts.append(f"<tr><td>{cls}</td><td>{int(s['n_families'])}</td><td>{int(s['total_genes'])}</td><td>{int(s['n_mags'])}</td></tr>")
    html_parts.append("</table></div>")

    # Differential pathways
    if diff_results:
        html_parts.append("<div class='section'>")
        html_parts.append("<h2>Differential Pathway Abundance</h2>")
        for (g1, g2), result in sorted(diff_results.items()):
            n_sig = int((result.q_values < 0.05).sum())
            html_parts.append(f"<h3>{g1} vs {g2}: {n_sig} significant pathways (FDR &lt; 0.05)</h3>")
            if n_sig > 0:
                html_parts.append("<table><tr><th>Function</th><th>CLR diff</th><th>q-value</th><th>Cohen's d</th></tr>")
                sig_idx = np.where(result.q_values < 0.05)[0]
                sorted_idx = sig_idx[np.argsort(np.abs(result.effect_sizes[sig_idx]))[::-1]]
                for i in sorted_idx[:20]:
                    q_class = ' class="sig"' if result.q_values[i] < 0.05 else ''
                    html_parts.append(
                        f"<tr><td>{result.mag_ids[i]}</td>"
                        f"<td>{result.log_fold_changes[i]:.3f}</td>"
                        f"<td{q_class}>{result.q_values[i]:.4f}</td>"
                        f"<td>{result.effect_sizes[i]:.3f}</td></tr>"
                    )
                html_parts.append("</table>")
        html_parts.append("</div>")

    html_parts.append("</body></html>")
    (out / "func_report.html").write_text("\n".join(html_parts))
```

Add plot functions to `mag/plots.py` (append to end of file):

```python
def plot_pgpr_matrix(
    pgpr_table: FunctionalTable,
    taxonomy: TaxonomyTable | None = None,
    top_n: int = 30,
    output: str | Path = "pgpr_traits.pdf",
) -> None:
    """Binary heatmap of PGPR traits per MAG."""
    from .func_io import FunctionalTable as FT

    _setup_style()
    # Select MAGs with at least one trait
    has_trait = pgpr_table.values.sum(axis=0) > 0
    mag_indices = np.where(has_trait)[0][:top_n]

    if len(mag_indices) == 0:
        fig, ax = plt.subplots(figsize=(8, 3))
        ax.text(0.5, 0.5, "No MAGs with PGPR traits detected", ha="center", va="center", transform=ax.transAxes)
        fig.savefig(str(output), bbox_inches="tight")
        plt.close(fig)
        return

    subset = pgpr_table.values[:, mag_indices]
    mag_labels = [pgpr_table.mag_ids[i] for i in mag_indices]
    if taxonomy:
        mag_labels = [
            f"{m} ({taxonomy.get(m).rank('phylum') if taxonomy.get(m) else ''})"
            for m in mag_labels
        ]

    fig, ax = plt.subplots(figsize=(max(8, len(mag_labels) * 0.4), max(4, len(pgpr_table.function_ids) * 0.5)))
    ax.imshow(subset, cmap="YlGn", aspect="auto", vmin=0, vmax=1)
    ax.set_xticks(range(len(mag_labels)))
    ax.set_xticklabels(mag_labels, rotation=90, fontsize=7)
    ax.set_yticks(range(len(pgpr_table.function_ids)))
    ax.set_yticklabels(pgpr_table.function_ids, fontsize=9)
    ax.set_title("PGPR Trait Presence")
    fig.savefig(str(output), bbox_inches="tight")
    plt.close(fig)


def plot_cazyme_bars(
    cazy_table: FunctionalTable,
    metadata: SampleMetadata,
    grouping_var: str,
    abundance: AbundanceTable,
    output: str | Path = "cazyme_bars.pdf",
) -> None:
    """Stacked bar chart of CAZyme classes per group."""
    from .func_io import FunctionalTable as FT
    from .func_profile import pathway_abundance

    _setup_style()
    groups = metadata.get_groups(grouping_var)
    group_names = sorted(groups.keys())

    # Get CAZyme abundance per sample
    cazy_abund = pathway_abundance(cazy_table, abundance)

    # Aggregate by class prefix per group
    class_totals: dict[str, dict[str, float]] = {}
    for fi, fam_id in enumerate(cazy_abund.mag_ids):
        prefix = ""
        for ch in fam_id:
            if ch.isalpha():
                prefix += ch
            else:
                break
        if not prefix:
            prefix = "Other"
        if prefix not in class_totals:
            class_totals[prefix] = {g: 0.0 for g in group_names}
        for g in group_names:
            g_samples = groups[g]
            g_idx = [cazy_abund.sample_ids.index(s) for s in g_samples if s in set(cazy_abund.sample_ids)]
            if g_idx:
                class_totals[prefix][g] += float(cazy_abund.abundances[fi, g_idx].mean())

    classes = sorted(class_totals.keys())
    fig, ax = plt.subplots(figsize=(max(6, len(group_names) * 1.5), 5))
    x = np.arange(len(group_names))
    bottom = np.zeros(len(group_names))
    colors = sns.color_palette("Set2", len(classes))

    for ci, cls in enumerate(classes):
        vals = np.array([class_totals[cls][g] for g in group_names])
        ax.bar(x, vals, bottom=bottom, label=cls, color=colors[ci % len(colors)])
        bottom += vals

    ax.set_xticks(x)
    ax.set_xticklabels(group_names, rotation=45, ha="right")
    ax.set_ylabel("Mean weighted abundance")
    ax.set_title("CAZyme Class Abundance by Group")
    ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
    fig.savefig(str(output), bbox_inches="tight")
    plt.close(fig)
```

Add `func-report` command to `mag/cli.py`:

```python
@main.command("func-report")
@click.option("--abundance", "-a", required=True, type=click.Path(exists=True), help="Abundance table TSV")
@click.option("--dram-annotations", "-d", required=True, type=click.Path(exists=True), help="DRAM annotations.tsv")
@click.option("--taxonomy", "-t", default=None, type=click.Path(exists=True), help="Taxonomy TSV (optional)")
@click.option("--metadata", "-m", required=True, type=click.Path(exists=True), help="Sample metadata TSV")
@click.option("--group", "-g", default="compartment", help="Metadata variable for grouping")
@click.option("--output", "-o", default="func_results", help="Output directory")
def func_report(abundance: str, dram_annotations: str, taxonomy: str | None, metadata: str, group: str, output: str) -> None:
    """Run functional profiling pipeline (magfunc)."""
    from .func_io import load_dram_annotations
    from .func_report import generate_func_report

    table = load_abundance_table(abundance)
    annots = load_dram_annotations(dram_annotations)
    tax = load_taxonomy(taxonomy) if taxonomy else None
    meta = load_metadata(metadata)

    generate_func_report(table, annots, tax, meta, group, output)
    click.echo(f"Functional profiling report written to {output}/")
```

**Step 4: Run tests to verify they pass**

Run: `.venv/bin/python -m pytest mag/tests/test_func_report.py -v`
Expected: All PASS

**Step 5: Commit**

```bash
git add mag/func_report.py mag/plots.py mag/cli.py mag/tests/test_func_report.py
git commit -m "feat: add functional report generation with PGPR, CAZyme, and HTML (magfunc)"
```

---

## Task 4: Network correlation and construction

**Files:**
- Create: `mag/net_correlation.py`
- Test: `mag/tests/test_net_correlation.py`

**Step 1: Write the failing tests**

```python
# mag/tests/test_net_correlation.py
"""Tests for network correlation and construction."""

from __future__ import annotations

import numpy as np
import pytest

from mag.net_correlation import (
    NetworkResult,
    proportionality,
    build_network,
)
from mag.io import AbundanceTable
from mag.tests.fixtures import generate_synthetic_abundance_table


@pytest.fixture
def synth():
    return generate_synthetic_abundance_table()


class TestProportionality:
    def test_output_shape(self, synth) -> None:
        table = synth[0]
        corr = proportionality(table)
        assert corr.shape == (table.n_mags, table.n_mags)

    def test_symmetric(self, synth) -> None:
        table = synth[0]
        corr = proportionality(table)
        np.testing.assert_array_almost_equal(corr, corr.T)

    def test_diagonal_zero(self, synth) -> None:
        table = synth[0]
        corr = proportionality(table)
        np.testing.assert_array_almost_equal(np.diag(corr), np.zeros(table.n_mags))

    def test_proportional_mags_low_phi(self) -> None:
        """Two MAGs with proportional abundances should have low phi."""
        abundances = np.array([
            [100, 200, 300, 400, 500],
            [50, 100, 150, 200, 250],   # proportional to MAG 0
            [500, 10, 300, 20, 400],     # unrelated
        ], dtype=np.float64)
        table = AbundanceTable(
            mag_ids=["M0", "M1", "M2"],
            sample_ids=[f"s{i}" for i in range(5)],
            abundances=abundances,
        )
        corr = proportionality(table)
        # phi(M0, M1) should be near 0 (proportional)
        # phi(M0, M2) should be larger
        assert corr[0, 1] < corr[0, 2]

    def test_min_prevalence_filter(self) -> None:
        """MAGs absent in most samples should be excluded."""
        abundances = np.array([
            [100, 200, 300, 400],
            [0, 0, 0, 50],  # only in 1/4 = 25%
            [50, 100, 150, 200],
        ], dtype=np.float64)
        table = AbundanceTable(
            mag_ids=["M0", "M1", "M2"],
            sample_ids=["s0", "s1", "s2", "s3"],
            abundances=abundances,
        )
        corr = proportionality(table, min_prevalence=0.5)
        # M1 excluded -> row/col should be NaN or zero
        assert np.isnan(corr[1, 0]) or corr[1, 0] == 0


class TestBuildNetwork:
    def test_returns_network_result(self, synth) -> None:
        table = synth[0]
        corr = proportionality(table)
        net = build_network(corr, table.mag_ids, threshold_percentile=10)
        assert isinstance(net, NetworkResult)
        assert len(net.mag_ids) == table.n_mags

    def test_edges_below_threshold(self, synth) -> None:
        table = synth[0]
        corr = proportionality(table)
        net = build_network(corr, table.mag_ids, threshold_percentile=10)
        if net.edges:
            max_weight = max(w for _, _, w in net.edges)
            assert max_weight <= net.threshold

    def test_no_self_edges(self, synth) -> None:
        table = synth[0]
        corr = proportionality(table)
        net = build_network(corr, table.mag_ids, threshold_percentile=50)
        for m1, m2, _ in net.edges:
            assert m1 != m2
```

**Step 2: Run tests to verify they fail**

Run: `.venv/bin/python -m pytest mag/tests/test_net_correlation.py -v`
Expected: FAIL

**Step 3: Implement net_correlation.py**

```python
# mag/net_correlation.py
"""Compositionally-aware co-occurrence network construction."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .io import AbundanceTable


@dataclass(frozen=True)
class NetworkResult:
    """Co-occurrence network from proportionality analysis."""

    mag_ids: list[str]
    adjacency: np.ndarray  # (n_mags, n_mags), phi values (lower = more proportional)
    edges: list[tuple[str, str, float]]  # (mag1, mag2, phi) filtered by threshold
    threshold: float


def proportionality(
    table: AbundanceTable,
    pseudocount: float = 1.0,
    min_prevalence: float = 0.0,
) -> np.ndarray:
    """Compute phi proportionality metric between all MAG pairs.

    phi(x, y) = var(log(x/y)) — lower values indicate more proportional
    (co-occurring) pairs.  This is compositionally aware unlike Pearson/Spearman.

    MAGs present in fewer than ``min_prevalence`` fraction of samples are
    masked (set to NaN in the output matrix).

    Returns:
        Symmetric matrix of shape (n_mags, n_mags) with phi values.
        Diagonal is 0.  Masked MAGs have NaN.
    """
    x = table.abundances + pseudocount
    log_x = np.log(x)
    n_mags = table.n_mags

    # Prevalence filter
    prevalence = (table.abundances > 0).sum(axis=1) / table.n_samples
    valid = prevalence >= min_prevalence

    phi = np.full((n_mags, n_mags), np.nan)

    for i in range(n_mags):
        if not valid[i]:
            continue
        phi[i, i] = 0.0
        for j in range(i + 1, n_mags):
            if not valid[j]:
                continue
            log_ratio = log_x[i] - log_x[j]
            val = float(np.var(log_ratio, ddof=1))
            phi[i, j] = val
            phi[j, i] = val

    return phi


def build_network(
    corr_matrix: np.ndarray,
    mag_ids: list[str],
    threshold_percentile: float = 5,
) -> NetworkResult:
    """Build a network by thresholding the proportionality matrix.

    Edges are retained where phi <= threshold (low phi = strong association).
    The threshold is set at the given percentile of all valid (non-NaN,
    non-diagonal) phi values.

    Args:
        corr_matrix: Symmetric phi matrix from proportionality().
        mag_ids: MAG identifiers matching matrix rows/columns.
        threshold_percentile: Percentile of phi values to use as cutoff.

    Returns:
        NetworkResult with filtered edges.
    """
    n = len(mag_ids)

    # Collect all valid off-diagonal values
    valid_values = []
    for i in range(n):
        for j in range(i + 1, n):
            v = corr_matrix[i, j]
            if not np.isnan(v):
                valid_values.append(v)

    if not valid_values:
        return NetworkResult(
            mag_ids=list(mag_ids),
            adjacency=corr_matrix,
            edges=[],
            threshold=0.0,
        )

    threshold = float(np.percentile(valid_values, threshold_percentile))

    edges: list[tuple[str, str, float]] = []
    for i in range(n):
        for j in range(i + 1, n):
            v = corr_matrix[i, j]
            if not np.isnan(v) and v <= threshold:
                edges.append((mag_ids[i], mag_ids[j], float(v)))

    return NetworkResult(
        mag_ids=list(mag_ids),
        adjacency=corr_matrix,
        edges=edges,
        threshold=threshold,
    )
```

**Step 4: Run tests to verify they pass**

Run: `.venv/bin/python -m pytest mag/tests/test_net_correlation.py -v`
Expected: All PASS

**Step 5: Commit**

```bash
git add mag/net_correlation.py mag/tests/test_net_correlation.py
git commit -m "feat: add proportionality correlation and network construction (magnet)"
```

---

## Task 5: Network topology and keystone taxa

**Files:**
- Create: `mag/net_topology.py`
- Test: `mag/tests/test_net_topology.py`

**Step 1: Write the failing tests**

```python
# mag/tests/test_net_topology.py
"""Tests for network topology and keystone taxa identification."""

from __future__ import annotations

import numpy as np
import pytest

from mag.net_correlation import NetworkResult
from mag.net_topology import (
    NetworkTopology,
    KeystoneTaxaResult,
    compute_topology,
    identify_keystones,
    differential_network,
)
from mag.io import AbundanceTable


@pytest.fixture
def simple_network() -> NetworkResult:
    """A small network with known topology."""
    # Star topology: M0 connected to all others
    mag_ids = ["M0", "M1", "M2", "M3", "M4"]
    adj = np.full((5, 5), 1.0)  # high phi (no edge)
    np.fill_diagonal(adj, 0.0)
    # Strong edges from M0 to all
    adj[0, 1] = adj[1, 0] = 0.1
    adj[0, 2] = adj[2, 0] = 0.1
    adj[0, 3] = adj[3, 0] = 0.1
    adj[0, 4] = adj[4, 0] = 0.1
    # One edge between M1-M2
    adj[1, 2] = adj[2, 1] = 0.1

    edges = [
        ("M0", "M1", 0.1), ("M0", "M2", 0.1),
        ("M0", "M3", 0.1), ("M0", "M4", 0.1),
        ("M1", "M2", 0.1),
    ]
    return NetworkResult(mag_ids=mag_ids, adjacency=adj, edges=edges, threshold=0.5)


class TestComputeTopology:
    def test_returns_topology(self, simple_network) -> None:
        topo = compute_topology(simple_network)
        assert isinstance(topo, NetworkTopology)
        assert len(topo.mag_ids) == 5

    def test_degree_correct(self, simple_network) -> None:
        topo = compute_topology(simple_network)
        # M0 has degree 4 (connected to M1,M2,M3,M4)
        m0_idx = topo.mag_ids.index("M0")
        assert topo.degree[m0_idx] == 4
        # M3 has degree 1 (only connected to M0)
        m3_idx = topo.mag_ids.index("M3")
        assert topo.degree[m3_idx] == 1

    def test_betweenness_hub(self, simple_network) -> None:
        topo = compute_topology(simple_network)
        m0_idx = topo.mag_ids.index("M0")
        # M0 is the hub -> highest betweenness
        assert topo.betweenness[m0_idx] == max(topo.betweenness)

    def test_modularity_computed(self, simple_network) -> None:
        topo = compute_topology(simple_network)
        assert isinstance(topo.modularity, float)

    def test_empty_network(self) -> None:
        net = NetworkResult(mag_ids=["M0", "M1"], adjacency=np.zeros((2, 2)), edges=[], threshold=0.5)
        topo = compute_topology(net)
        assert all(d == 0 for d in topo.degree)


class TestIdentifyKeystones:
    def test_hub_is_keystone(self, simple_network) -> None:
        topo = compute_topology(simple_network)
        abundance = AbundanceTable(
            mag_ids=["M0", "M1", "M2", "M3", "M4"],
            sample_ids=["s0"],
            abundances=np.array([[10], [100], [100], [100], [100]], dtype=np.float64),
        )
        result = identify_keystones(topo, abundance)
        assert isinstance(result, KeystoneTaxaResult)
        # M0 has high centrality + low abundance -> keystone
        m0_idx = result.mag_ids.index("M0")
        assert result.is_keystone[m0_idx]

    def test_empty_network_no_keystones(self) -> None:
        net = NetworkResult(mag_ids=["M0"], adjacency=np.zeros((1, 1)), edges=[], threshold=0.5)
        topo = compute_topology(net)
        abundance = AbundanceTable(
            mag_ids=["M0"], sample_ids=["s0"],
            abundances=np.array([[100]], dtype=np.float64),
        )
        result = identify_keystones(topo, abundance)
        assert not any(result.is_keystone)


class TestDifferentialNetwork:
    def test_finds_gained_lost_edges(self) -> None:
        net1 = NetworkResult(
            mag_ids=["M0", "M1", "M2"],
            adjacency=np.zeros((3, 3)),
            edges=[("M0", "M1", 0.1), ("M0", "M2", 0.2)],
            threshold=0.5,
        )
        net2 = NetworkResult(
            mag_ids=["M0", "M1", "M2"],
            adjacency=np.zeros((3, 3)),
            edges=[("M0", "M1", 0.1), ("M1", "M2", 0.15)],
            threshold=0.5,
        )
        diff = differential_network(net1, net2)
        assert ("M0", "M2") in diff["lost"]  # in net1 but not net2
        assert ("M1", "M2") in diff["gained"]  # in net2 but not net1
        assert ("M0", "M1") in diff["conserved"]
```

**Step 2: Run tests to verify they fail**

Run: `.venv/bin/python -m pytest mag/tests/test_net_topology.py -v`
Expected: FAIL

**Step 3: Implement net_topology.py**

```python
# mag/net_topology.py
"""Network topology metrics and keystone taxa identification."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .io import AbundanceTable
from .net_correlation import NetworkResult


@dataclass(frozen=True)
class NetworkTopology:
    """Network topology metrics per MAG."""

    mag_ids: list[str]
    degree: np.ndarray
    betweenness: np.ndarray
    closeness: np.ndarray
    modularity: float
    module_assignments: np.ndarray  # integer labels
    hub_scores: np.ndarray


@dataclass(frozen=True)
class KeystoneTaxaResult:
    """Identified keystone taxa."""

    mag_ids: list[str]
    keystone_scores: np.ndarray
    is_keystone: np.ndarray  # boolean
    metrics: dict[str, np.ndarray]


def compute_topology(network: NetworkResult) -> NetworkTopology:
    """Compute network topology metrics using networkx.

    Metrics: degree, betweenness centrality, closeness centrality,
    Louvain modularity, hub scores.
    """
    import networkx as nx

    G = nx.Graph()
    G.add_nodes_from(network.mag_ids)
    for m1, m2, w in network.edges:
        G.add_edge(m1, m2, weight=w)

    n = len(network.mag_ids)
    mag_idx = {m: i for i, m in enumerate(network.mag_ids)}

    # Degree
    degree = np.zeros(n)
    for m, d in G.degree():
        degree[mag_idx[m]] = d

    # Betweenness centrality
    betweenness = np.zeros(n)
    if G.number_of_edges() > 0:
        bc = nx.betweenness_centrality(G)
        for m, v in bc.items():
            betweenness[mag_idx[m]] = v

    # Closeness centrality
    closeness = np.zeros(n)
    if G.number_of_edges() > 0:
        cc = nx.closeness_centrality(G)
        for m, v in cc.items():
            closeness[mag_idx[m]] = v

    # Modularity (Louvain / greedy modularity)
    modularity = 0.0
    module_assignments = np.zeros(n, dtype=int)
    if G.number_of_edges() > 0:
        try:
            communities = nx.community.greedy_modularity_communities(G)
            for ci, comm in enumerate(communities):
                for m in comm:
                    module_assignments[mag_idx[m]] = ci
            modularity = float(nx.community.modularity(G, communities))
        except Exception:
            pass

    # Hub scores (HITS algorithm)
    hub_scores = np.zeros(n)
    if G.number_of_edges() > 0:
        try:
            hubs, _ = nx.hits(G, max_iter=100)
            for m, v in hubs.items():
                hub_scores[mag_idx[m]] = v
        except nx.PowerIterationFailedConvergence:
            # Fall back to degree centrality as hub proxy
            for m, d in G.degree():
                hub_scores[mag_idx[m]] = d / max(1, n - 1)

    return NetworkTopology(
        mag_ids=list(network.mag_ids),
        degree=degree,
        betweenness=betweenness,
        closeness=closeness,
        modularity=modularity,
        module_assignments=module_assignments,
        hub_scores=hub_scores,
    )


def identify_keystones(
    topology: NetworkTopology,
    abundance: AbundanceTable,
    betweenness_threshold: float = 0.5,
    abundance_threshold: float = 0.5,
) -> KeystoneTaxaResult:
    """Identify keystone taxa: high centrality + low relative abundance.

    Keystone score = normalized_betweenness * (1 - normalized_abundance).
    A MAG is keystone if its score is above the median of non-zero scores
    and its betweenness is in the top percentile and abundance in the bottom.
    """
    n = len(topology.mag_ids)

    # Align abundance to topology MAGs
    mean_abund = np.zeros(n)
    abund_idx = {m: i for i, m in enumerate(abundance.mag_ids)}
    for i, m in enumerate(topology.mag_ids):
        if m in abund_idx:
            mean_abund[i] = abundance.abundances[abund_idx[m]].mean()

    # Normalize to [0, 1]
    max_between = topology.betweenness.max()
    norm_between = topology.betweenness / max_between if max_between > 0 else np.zeros(n)

    total_abund = mean_abund.sum()
    norm_abund = mean_abund / total_abund if total_abund > 0 else np.zeros(n)

    # Keystone score
    keystone_scores = norm_between * (1 - norm_abund)

    # Classify: keystone if high betweenness and low abundance
    is_keystone = np.zeros(n, dtype=bool)
    if max_between > 0 and total_abund > 0:
        high_between = norm_between >= np.percentile(norm_between[norm_between > 0], 100 * betweenness_threshold) if (norm_between > 0).any() else np.zeros(n, dtype=bool)
        low_abund = norm_abund <= np.percentile(norm_abund[norm_abund > 0], 100 * abundance_threshold) if (norm_abund > 0).any() else np.zeros(n, dtype=bool)
        is_keystone = high_between & low_abund

    return KeystoneTaxaResult(
        mag_ids=list(topology.mag_ids),
        keystone_scores=keystone_scores,
        is_keystone=is_keystone,
        metrics={
            "betweenness": topology.betweenness,
            "closeness": topology.closeness,
            "hub_score": topology.hub_scores,
            "mean_abundance": mean_abund,
        },
    )


def differential_network(
    net1: NetworkResult,
    net2: NetworkResult,
) -> dict[str, list[tuple[str, str]]]:
    """Compare two networks: find conserved, gained, and lost edges.

    Returns dict with keys "conserved", "gained", "lost".
    """
    def _edge_set(edges: list[tuple[str, str, float]]) -> set[tuple[str, str]]:
        s = set()
        for m1, m2, _ in edges:
            key = tuple(sorted([m1, m2]))
            s.add(key)
        return s

    edges1 = _edge_set(net1.edges)
    edges2 = _edge_set(net2.edges)

    return {
        "conserved": sorted(edges1 & edges2),
        "lost": sorted(edges1 - edges2),
        "gained": sorted(edges2 - edges1),
    }
```

**Step 4: Run tests to verify they pass**

Run: `.venv/bin/python -m pytest mag/tests/test_net_topology.py -v`
Expected: All PASS

**Step 5: Commit**

```bash
git add mag/net_topology.py mag/tests/test_net_topology.py
git commit -m "feat: add network topology metrics and keystone taxa (magnet)"
```

---

## Task 6: Network report generation (CSV + plots + HTML)

**Files:**
- Create: `mag/net_report.py`
- Modify: `mag/plots.py` — add `plot_network_graph`, `plot_degree_distribution`
- Modify: `mag/cli.py` — add `net-report` command
- Test: `mag/tests/test_net_report.py`

**Step 1: Write the failing tests**

```python
# mag/tests/test_net_report.py
"""Tests for network report generation."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from mag.net_report import generate_net_report
from mag.tests.fixtures import generate_synthetic_abundance_table


@pytest.fixture
def synth():
    return generate_synthetic_abundance_table()


class TestGenerateNetReport:
    def test_creates_output_files(self, synth, tmp_path) -> None:
        table, meta, tax = synth
        generate_net_report(
            abundance=table,
            taxonomy=tax,
            metadata=meta,
            grouping_var="compartment",
            output_dir=tmp_path,
        )
        assert (tmp_path / "network_edges.csv").exists()
        assert (tmp_path / "network_topology.csv").exists()
        assert (tmp_path / "keystone_taxa.csv").exists()
        assert (tmp_path / "net_report.html").exists()

    def test_per_group_networks(self, synth, tmp_path) -> None:
        table, meta, tax = synth
        generate_net_report(
            abundance=table,
            taxonomy=tax,
            metadata=meta,
            grouping_var="compartment",
            output_dir=tmp_path,
        )
        # Should have per-group network edge files
        group_files = list(tmp_path.glob("network_edges_*.csv"))
        assert len(group_files) >= 1

    def test_differential_network_csv(self, synth, tmp_path) -> None:
        table, meta, tax = synth
        generate_net_report(
            abundance=table,
            taxonomy=tax,
            metadata=meta,
            grouping_var="compartment",
            output_dir=tmp_path,
        )
        diff_files = list(tmp_path.glob("differential_network_*.csv"))
        assert len(diff_files) >= 1

    def test_no_taxonomy_ok(self, synth, tmp_path) -> None:
        table, meta, _ = synth
        generate_net_report(
            abundance=table,
            taxonomy=None,
            metadata=meta,
            grouping_var="compartment",
            output_dir=tmp_path,
        )
        assert (tmp_path / "net_report.html").exists()

    def test_plots_created(self, synth, tmp_path) -> None:
        table, meta, tax = synth
        generate_net_report(
            abundance=table,
            taxonomy=tax,
            metadata=meta,
            grouping_var="compartment",
            output_dir=tmp_path,
        )
        assert (tmp_path / "network_graph.pdf").exists()
        assert (tmp_path / "degree_distribution.pdf").exists()
```

**Step 2: Run tests to verify they fail**

Run: `.venv/bin/python -m pytest mag/tests/test_net_report.py -v`
Expected: FAIL

**Step 3: Implement net_report.py and add plots**

Create `mag/net_report.py`:

```python
# mag/net_report.py
"""Report generation for network analysis."""

from __future__ import annotations

import base64
import csv
import logging
from itertools import combinations
from pathlib import Path

import numpy as np

from .io import AbundanceTable, SampleMetadata, TaxonomyTable
from .net_correlation import NetworkResult, build_network, proportionality
from .net_topology import (
    KeystoneTaxaResult,
    NetworkTopology,
    compute_topology,
    differential_network,
    identify_keystones,
)

logger = logging.getLogger(__name__)


def generate_net_report(
    abundance: AbundanceTable,
    taxonomy: TaxonomyTable | None,
    metadata: SampleMetadata,
    grouping_var: str,
    output_dir: str | Path,
    threshold_percentile: float = 5,
    min_prevalence: float = 0.5,
) -> None:
    """Run network analysis pipeline and write results."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    # Global network
    corr = proportionality(abundance, min_prevalence=min_prevalence)
    global_net = build_network(corr, abundance.mag_ids, threshold_percentile)
    _write_edges_csv(global_net, out / "network_edges.csv", taxonomy)

    global_topo = compute_topology(global_net)
    _write_topology_csv(global_topo, out / "network_topology.csv", taxonomy)

    keystones = identify_keystones(global_topo, abundance)
    _write_keystones_csv(keystones, out / "keystone_taxa.csv", taxonomy)

    # Per-group networks
    groups = metadata.get_groups(grouping_var)
    group_names = sorted(groups.keys())
    group_nets: dict[str, NetworkResult] = {}

    for g in group_names:
        try:
            sub = abundance.subset_samples(groups[g])
            g_corr = proportionality(sub, min_prevalence=0.0)
            g_net = build_network(g_corr, sub.mag_ids, threshold_percentile)
            group_nets[g] = g_net
            _write_edges_csv(g_net, out / f"network_edges_{g}.csv", taxonomy)
        except Exception as e:
            logger.warning("Network for group %s failed: %s", g, e)

    # Differential networks
    diff_results: dict[tuple[str, str], dict] = {}
    for g1, g2 in combinations(group_names, 2):
        if g1 in group_nets and g2 in group_nets:
            diff = differential_network(group_nets[g1], group_nets[g2])
            diff_results[(g1, g2)] = diff
            _write_diff_network_csv(diff, g1, g2, out / f"differential_network_{g1}_vs_{g2}.csv")

    # Plots
    from . import plots

    plots.plot_network_graph(global_net, taxonomy, output=out / "network_graph.pdf")
    plots.plot_degree_distribution(global_topo, output=out / "degree_distribution.pdf")

    # HTML report
    _write_net_html(
        out,
        global_net=global_net,
        global_topo=global_topo,
        keystones=keystones,
        group_nets=group_nets,
        diff_results=diff_results,
        taxonomy=taxonomy,
    )


def _write_edges_csv(net: NetworkResult, path: Path, taxonomy: TaxonomyTable | None) -> None:
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        header = ["MAG_1", "MAG_2", "phi"]
        if taxonomy:
            header.extend(["phylum_1", "phylum_2"])
        w.writerow(header)
        for m1, m2, phi in sorted(net.edges, key=lambda e: e[2]):
            row = [m1, m2, f"{phi:.6f}"]
            if taxonomy:
                r1 = taxonomy.get(m1)
                r2 = taxonomy.get(m2)
                row.extend([r1.phylum if r1 else "", r2.phylum if r2 else ""])
            w.writerow(row)


def _write_topology_csv(topo: NetworkTopology, path: Path, taxonomy: TaxonomyTable | None) -> None:
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        header = ["MAG_ID"]
        if taxonomy:
            header.extend(["phylum", "class", "genus"])
        header.extend(["degree", "betweenness", "closeness", "hub_score", "module"])
        w.writerow(header)
        for i, m in enumerate(topo.mag_ids):
            row = [m]
            if taxonomy:
                rec = taxonomy.get(m)
                row.extend([
                    rec.phylum if rec else "",
                    rec.rank("class") if rec else "",
                    rec.genus if rec else "",
                ])
            row.extend([
                int(topo.degree[i]),
                f"{topo.betweenness[i]:.6f}",
                f"{topo.closeness[i]:.6f}",
                f"{topo.hub_scores[i]:.6f}",
                int(topo.module_assignments[i]),
            ])
            w.writerow(row)


def _write_keystones_csv(ks: KeystoneTaxaResult, path: Path, taxonomy: TaxonomyTable | None) -> None:
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        header = ["MAG_ID"]
        if taxonomy:
            header.extend(["phylum", "class", "genus"])
        header.extend(["keystone_score", "is_keystone", "betweenness", "mean_abundance"])
        w.writerow(header)
        # Sort by keystone score descending
        order = np.argsort(ks.keystone_scores)[::-1]
        for i in order:
            m = ks.mag_ids[i]
            row = [m]
            if taxonomy:
                rec = taxonomy.get(m)
                row.extend([
                    rec.phylum if rec else "",
                    rec.rank("class") if rec else "",
                    rec.genus if rec else "",
                ])
            row.extend([
                f"{ks.keystone_scores[i]:.6f}",
                "yes" if ks.is_keystone[i] else "no",
                f"{ks.metrics['betweenness'][i]:.6f}",
                f"{ks.metrics['mean_abundance'][i]:.2f}",
            ])
            w.writerow(row)


def _write_diff_network_csv(diff: dict, g1: str, g2: str, path: Path) -> None:
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["MAG_1", "MAG_2", "status"])
        for m1, m2 in diff.get("conserved", []):
            w.writerow([m1, m2, "conserved"])
        for m1, m2 in diff.get("gained", []):
            w.writerow([m1, m2, f"gained_in_{g2}"])
        for m1, m2 in diff.get("lost", []):
            w.writerow([m1, m2, f"lost_in_{g2}"])


def _embed_figure(path: Path) -> str:
    if not path.exists():
        return ""
    data = path.read_bytes()
    b64 = base64.b64encode(data).decode("ascii")
    return f'<img src="data:image/png;base64,{b64}" style="max-width:100%">'


def _write_net_html(
    out: Path,
    global_net: NetworkResult,
    global_topo: NetworkTopology,
    keystones: KeystoneTaxaResult,
    group_nets: dict[str, NetworkResult],
    diff_results: dict[tuple[str, str], dict],
    taxonomy: TaxonomyTable | None,
) -> None:
    """Generate self-contained HTML report for network analysis."""
    html_parts = [
        "<!DOCTYPE html><html><head>",
        "<meta charset='utf-8'>",
        "<title>magnet — Network Analysis Report</title>",
        "<style>",
        "body{font-family:sans-serif;max-width:1200px;margin:0 auto;padding:20px;line-height:1.6}",
        "h1{border-bottom:3px solid #2c3e50;padding-bottom:10px}",
        "h2{color:#2c3e50;cursor:pointer;border-bottom:1px solid #ddd;padding:8px 0}",
        "table{border-collapse:collapse;width:100%;margin:10px 0}",
        "th,td{border:1px solid #ddd;padding:6px 10px;text-align:left}",
        "th{background:#f5f5f5}",
        "tr:nth-child(even){background:#fafafa}",
        ".section{margin-bottom:30px}",
        ".keystone{color:#e74c3c;font-weight:bold}",
        "</style></head><body>",
        "<h1>magnet — Network Analysis Report</h1>",
    ]

    # Summary
    html_parts.append("<div class='section'>")
    html_parts.append("<h2>Network Summary</h2>")
    html_parts.append(f"<p>MAGs in network: <strong>{len(global_net.mag_ids)}</strong></p>")
    html_parts.append(f"<p>Edges (phi &le; {global_net.threshold:.4f}): <strong>{len(global_net.edges)}</strong></p>")
    html_parts.append(f"<p>Modularity: <strong>{global_topo.modularity:.4f}</strong></p>")
    n_modules = len(set(global_topo.module_assignments))
    html_parts.append(f"<p>Modules detected: <strong>{n_modules}</strong></p>")
    n_ks = int(keystones.is_keystone.sum())
    html_parts.append(f"<p>Keystone taxa: <strong>{n_ks}</strong></p>")
    html_parts.append("</div>")

    # Network graph
    html_parts.append("<div class='section'>")
    html_parts.append("<h2>Network Graph</h2>")
    fig = _embed_figure(out / "network_graph.pdf")
    if fig:
        html_parts.append(fig)
    html_parts.append("</div>")

    # Degree distribution
    html_parts.append("<div class='section'>")
    html_parts.append("<h2>Degree Distribution</h2>")
    fig = _embed_figure(out / "degree_distribution.pdf")
    if fig:
        html_parts.append(fig)
    html_parts.append("</div>")

    # Keystone taxa table
    if n_ks > 0:
        html_parts.append("<div class='section'>")
        html_parts.append("<h2>Keystone Taxa</h2>")
        html_parts.append("<table><tr><th>MAG</th>")
        if taxonomy:
            html_parts.append("<th>Phylum</th><th>Genus</th>")
        html_parts.append("<th>Keystone Score</th><th>Betweenness</th><th>Mean Abundance</th></tr>")
        order = np.argsort(keystones.keystone_scores)[::-1]
        for i in order:
            if not keystones.is_keystone[i]:
                continue
            m = keystones.mag_ids[i]
            html_parts.append(f"<tr><td class='keystone'>{m}</td>")
            if taxonomy:
                rec = taxonomy.get(m)
                html_parts.append(f"<td>{rec.phylum if rec else ''}</td>")
                html_parts.append(f"<td>{rec.genus if rec else ''}</td>")
            html_parts.append(
                f"<td>{keystones.keystone_scores[i]:.4f}</td>"
                f"<td>{keystones.metrics['betweenness'][i]:.4f}</td>"
                f"<td>{keystones.metrics['mean_abundance'][i]:.1f}</td></tr>"
            )
        html_parts.append("</table></div>")

    # Per-group comparison
    if group_nets:
        html_parts.append("<div class='section'>")
        html_parts.append("<h2>Per-Group Networks</h2>")
        html_parts.append("<table><tr><th>Group</th><th>Edges</th></tr>")
        for g in sorted(group_nets.keys()):
            html_parts.append(f"<tr><td>{g}</td><td>{len(group_nets[g].edges)}</td></tr>")
        html_parts.append("</table></div>")

    # Differential networks
    if diff_results:
        html_parts.append("<div class='section'>")
        html_parts.append("<h2>Differential Networks</h2>")
        for (g1, g2), diff in sorted(diff_results.items()):
            html_parts.append(f"<h3>{g1} vs {g2}</h3>")
            html_parts.append(f"<p>Conserved edges: {len(diff['conserved'])}, "
                            f"Gained in {g2}: {len(diff['gained'])}, "
                            f"Lost in {g2}: {len(diff['lost'])}</p>")
        html_parts.append("</div>")

    html_parts.append("</body></html>")
    (out / "net_report.html").write_text("\n".join(html_parts))
```

Add plot functions to `mag/plots.py` (append to end):

```python
def plot_network_graph(
    network: NetworkResult,
    taxonomy: TaxonomyTable | None = None,
    output: str | Path = "network_graph.pdf",
) -> None:
    """Force-directed network graph with nodes colored by phylum."""
    import networkx as nx
    from .net_correlation import NetworkResult as NR

    _setup_style()
    G = nx.Graph()
    G.add_nodes_from(network.mag_ids)
    for m1, m2, w in network.edges:
        G.add_edge(m1, m2, weight=w)

    if G.number_of_edges() == 0:
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.text(0.5, 0.5, "No edges in network", ha="center", va="center", transform=ax.transAxes)
        fig.savefig(str(output), bbox_inches="tight")
        plt.close(fig)
        return

    fig, ax = plt.subplots(figsize=(10, 8))
    pos = nx.spring_layout(G, seed=42, k=2/np.sqrt(len(network.mag_ids)))

    # Color by phylum
    if taxonomy:
        phyla = set()
        for m in network.mag_ids:
            rec = taxonomy.get(m)
            phyla.add(rec.phylum if rec else "Unknown")
        phyla_list = sorted(phyla)
        phyla_colors = {p: PALETTE[i % len(PALETTE)] for i, p in enumerate(phyla_list)}

        node_colors = []
        for m in G.nodes():
            rec = taxonomy.get(m)
            p = rec.phylum if rec else "Unknown"
            node_colors.append(phyla_colors[p])
    else:
        node_colors = [PALETTE[0]] * len(G.nodes())

    # Node size by degree
    degrees = dict(G.degree())
    node_sizes = [max(20, degrees.get(m, 0) * 15) for m in G.nodes()]

    nx.draw_networkx_edges(G, pos, alpha=0.2, ax=ax)
    nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_sizes, alpha=0.8, ax=ax)
    ax.set_title(f"Co-occurrence Network ({len(network.edges)} edges)")
    ax.axis("off")

    if taxonomy:
        # Legend for top phyla
        from matplotlib.lines import Line2D
        phyla_in_net = set()
        for m in G.nodes():
            rec = taxonomy.get(m)
            phyla_in_net.add(rec.phylum if rec else "Unknown")
        legend_elements = [
            Line2D([0], [0], marker="o", color="w", markerfacecolor=phyla_colors[p], markersize=8, label=p)
            for p in sorted(phyla_in_net)[:10]
        ]
        ax.legend(handles=legend_elements, loc="upper left", fontsize=7)

    fig.savefig(str(output), bbox_inches="tight")
    plt.close(fig)


def plot_degree_distribution(
    topology: NetworkTopology,
    output: str | Path = "degree_distribution.pdf",
) -> None:
    """Histogram of node degree distribution."""
    from .net_topology import NetworkTopology as NT

    _setup_style()
    fig, ax = plt.subplots(figsize=(7, 5))
    degrees = topology.degree[topology.degree > 0]
    if len(degrees) == 0:
        ax.text(0.5, 0.5, "No connected nodes", ha="center", va="center", transform=ax.transAxes)
    else:
        ax.hist(degrees, bins=min(30, int(degrees.max())), color=PALETTE[0], edgecolor="white")
        ax.set_xlabel("Degree")
        ax.set_ylabel("Number of MAGs")
        ax.set_title("Degree Distribution")
        ax.axvline(np.mean(degrees), color="red", linestyle="--", label=f"Mean={np.mean(degrees):.1f}")
        ax.legend()
    fig.savefig(str(output), bbox_inches="tight")
    plt.close(fig)
```

Add `net-report` command to `mag/cli.py`:

```python
@main.command("net-report")
@click.option("--abundance", "-a", required=True, type=click.Path(exists=True), help="Abundance table TSV")
@click.option("--taxonomy", "-t", default=None, type=click.Path(exists=True), help="Taxonomy TSV (optional)")
@click.option("--metadata", "-m", required=True, type=click.Path(exists=True), help="Sample metadata TSV")
@click.option("--group", "-g", default="compartment", help="Metadata variable for grouping")
@click.option("--threshold", default=5.0, type=float, help="Phi percentile threshold for edges (default 5)")
@click.option("--min-prevalence", default=0.5, type=float, help="Min prevalence to include MAG in network (default 0.5)")
@click.option("--output", "-o", default="net_results", help="Output directory")
def net_report(abundance: str, taxonomy: str | None, metadata: str, group: str, threshold: float, min_prevalence: float, output: str) -> None:
    """Run network analysis pipeline (magnet)."""
    from .net_report import generate_net_report

    table = load_abundance_table(abundance)
    tax = load_taxonomy(taxonomy) if taxonomy else None
    meta = load_metadata(metadata)

    generate_net_report(table, tax, meta, group, output, threshold_percentile=threshold, min_prevalence=min_prevalence)
    click.echo(f"Network analysis report written to {output}/")
```

**Step 4: Run tests to verify they pass**

Run: `.venv/bin/python -m pytest mag/tests/test_net_report.py -v`
Expected: All PASS

**Step 5: Commit**

```bash
git add mag/net_report.py mag/plots.py mag/cli.py mag/tests/test_net_report.py
git commit -m "feat: add network report generation with graph, topology, keystones (magnet)"
```

---

## Task 7: Install networkx dependency and run full test suite

**Step 1: Install networkx**

Run: `.venv/bin/pip install networkx`

**Step 2: Run full test suite**

Run: `.venv/bin/python -m pytest mag/tests/ -v`
Expected: All tests PASS (existing 76 + new ~35 = ~111 total)

**Step 3: Final commit if any fixes needed**

```bash
git add -A mag/
git commit -m "chore: install networkx and fix any test issues"
```

---

## Verification Checklist

After all tasks complete:

1. `.venv/bin/python -m pytest mag/tests/ -v` — all tests pass
2. Check all new files exist:
   - `mag/func_io.py`, `mag/func_profile.py`, `mag/func_report.py`
   - `mag/net_correlation.py`, `mag/net_topology.py`, `mag/net_report.py`
   - 6 new test files
3. CLI commands work: `python -m mag.cli func-report --help` and `net-report --help`
4. Push to magprofile remote
