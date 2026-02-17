"""Data loading and validation for MAG community profiling."""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path
from typing import Sequence

import numpy as np


TAXONOMY_RANKS = ("domain", "phylum", "class", "order", "family", "genus", "species")


@dataclass
class TaxonomyRecord:
    """Taxonomy annotation for a single MAG."""

    mag_id: str
    domain: str = ""
    phylum: str = ""
    class_: str = ""
    order: str = ""
    family: str = ""
    genus: str = ""
    species: str = ""

    def rank(self, level: str) -> str:
        if level == "class":
            return self.class_
        return getattr(self, level)


@dataclass
class AbundanceTable:
    """MAG-by-sample abundance matrix."""

    mag_ids: list[str]
    sample_ids: list[str]
    abundances: np.ndarray  # shape (n_mags, n_samples)

    def __post_init__(self) -> None:
        n_mags, n_samples = self.abundances.shape
        if n_mags != len(self.mag_ids):
            raise ValueError(
                f"Row count {n_mags} != len(mag_ids) {len(self.mag_ids)}"
            )
        if n_samples != len(self.sample_ids):
            raise ValueError(
                f"Col count {n_samples} != len(sample_ids) {len(self.sample_ids)}"
            )

    @property
    def n_mags(self) -> int:
        return len(self.mag_ids)

    @property
    def n_samples(self) -> int:
        return len(self.sample_ids)

    def normalize(self) -> AbundanceTable:
        """Return relative-abundance table (columns sum to 1)."""
        col_sums = self.abundances.sum(axis=0, keepdims=True)
        col_sums = np.where(col_sums == 0, 1.0, col_sums)
        return AbundanceTable(
            mag_ids=list(self.mag_ids),
            sample_ids=list(self.sample_ids),
            abundances=self.abundances / col_sums,
        )

    def filter_low_abundance(self, min_total: float = 0.0) -> AbundanceTable:
        """Remove MAGs with total abundance <= min_total."""
        totals = self.abundances.sum(axis=1)
        mask = totals > min_total
        return AbundanceTable(
            mag_ids=[m for m, keep in zip(self.mag_ids, mask) if keep],
            sample_ids=list(self.sample_ids),
            abundances=self.abundances[mask],
        )

    def filter_prevalence(self, min_prevalence: float = 0.0) -> AbundanceTable:
        """Remove MAGs present in fewer than min_prevalence fraction of samples."""
        prevalence = (self.abundances > 0).sum(axis=1) / self.n_samples
        mask = prevalence >= min_prevalence
        return AbundanceTable(
            mag_ids=[m for m, keep in zip(self.mag_ids, mask) if keep],
            sample_ids=list(self.sample_ids),
            abundances=self.abundances[mask],
        )

    def subset_samples(self, sample_ids: Sequence[str]) -> AbundanceTable:
        """Return table with only the specified samples."""
        idx_map = {s: i for i, s in enumerate(self.sample_ids)}
        indices = [idx_map[s] for s in sample_ids]
        return AbundanceTable(
            mag_ids=list(self.mag_ids),
            sample_ids=list(sample_ids),
            abundances=self.abundances[:, indices],
        )


@dataclass
class TaxonomyTable:
    """Collection of taxonomy records keyed by MAG ID."""

    records: dict[str, TaxonomyRecord] = field(default_factory=dict)

    def get(self, mag_id: str) -> TaxonomyRecord | None:
        return self.records.get(mag_id)

    def aggregate_at_rank(
        self, abundance: AbundanceTable, rank: str
    ) -> dict[str, np.ndarray]:
        """Aggregate abundances at a given taxonomic rank.

        Returns dict mapping taxon name -> array of shape (n_samples,).
        """
        result: dict[str, np.ndarray] = {}
        for i, mag_id in enumerate(abundance.mag_ids):
            rec = self.records.get(mag_id)
            taxon = rec.rank(rank) if rec else ""
            if not taxon:
                taxon = "Unclassified"
            if taxon not in result:
                result[taxon] = np.zeros(abundance.n_samples)
            result[taxon] += abundance.abundances[i]
        return result

    def lookup_columns(
        self,
        mag_ids: Sequence[str],
        ranks: Sequence[str] = ("phylum", "class", "genus"),
    ) -> list[dict[str, str]]:
        """Return taxonomy columns for a list of MAG IDs.

        Each entry is a dict mapping rank name -> taxon string.
        Missing MAGs get empty strings for all ranks.
        """
        rows: list[dict[str, str]] = []
        for mag_id in mag_ids:
            rec = self.records.get(mag_id)
            if rec:
                rows.append({r: rec.rank(r) for r in ranks})
            else:
                rows.append({r: "" for r in ranks})
        return rows

    def aggregate_to_abundance_table(
        self, abundance: AbundanceTable, rank: str
    ) -> AbundanceTable:
        """Aggregate abundances at a rank and return as an AbundanceTable.

        The returned table has taxon names as ``mag_ids``.
        """
        agg = self.aggregate_at_rank(abundance, rank)
        taxon_names = sorted(agg.keys())
        matrix = np.array([agg[t] for t in taxon_names])
        return AbundanceTable(
            mag_ids=taxon_names,
            sample_ids=list(abundance.sample_ids),
            abundances=matrix,
        )


@dataclass
class SampleMetadata:
    """Sample metadata keyed by sample ID."""

    records: dict[str, dict[str, str]] = field(default_factory=dict)

    def get_groups(self, variable: str) -> dict[str, list[str]]:
        """Group sample IDs by a metadata variable.

        Returns dict mapping group_value -> list of sample_ids.
        """
        groups: dict[str, list[str]] = {}
        for sample_id, meta in self.records.items():
            val = meta.get(variable, "")
            groups.setdefault(val, []).append(sample_id)
        return groups

    @property
    def sample_ids(self) -> list[str]:
        return list(self.records.keys())


def load_abundance_table(path: str | Path) -> AbundanceTable:
    """Load a tab-separated abundance table.

    First column = MAG_ID, remaining columns = sample counts.
    """
    path = Path(path)
    with open(path, newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        sample_ids = header[1:]
        mag_ids: list[str] = []
        rows: list[list[float]] = []
        for row in reader:
            if not row or not row[0].strip():
                continue
            mag_ids.append(row[0].strip())
            rows.append([float(x) for x in row[1:]])
    abundances = np.array(rows, dtype=np.float64)
    return AbundanceTable(mag_ids=mag_ids, sample_ids=sample_ids, abundances=abundances)


def load_taxonomy(path: str | Path) -> TaxonomyTable:
    """Load a tab-separated taxonomy table (GTDB-Tk format)."""
    path = Path(path)
    records: dict[str, TaxonomyRecord] = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            mag_id = row.get("MAG_ID", row.get("mag_id", "")).strip()
            if not mag_id:
                continue
            records[mag_id] = TaxonomyRecord(
                mag_id=mag_id,
                domain=row.get("domain", ""),
                phylum=row.get("phylum", ""),
                class_=row.get("class", ""),
                order=row.get("order", ""),
                family=row.get("family", ""),
                genus=row.get("genus", ""),
                species=row.get("species", ""),
            )
    return TaxonomyTable(records=records)


def load_metadata(path: str | Path) -> SampleMetadata:
    """Load a tab-separated sample metadata table.

    First column = sample_id, remaining columns = covariates.
    """
    path = Path(path)
    records: dict[str, dict[str, str]] = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            sample_id = row.get("sample_id", "").strip()
            if not sample_id:
                continue
            records[sample_id] = {k: v for k, v in row.items() if k != "sample_id"}
    return SampleMetadata(records=records)
