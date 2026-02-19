"""Functional annotation I/O: DRAM parsing and FunctionalTable data type."""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np

from .io import AbundanceTable


# PGPR (Plant Growth-Promoting Rhizobacteria) marker genes → KEGG Orthology IDs
PGPR_MARKERS: dict[str, str] = {
    "nifH": "K02588",
    "nifD": "K02586",
    "nifK": "K02591",
    "pqqC": "K06137",
    "gcd": "K00111",
    "acdS": "K01505",
    "ipdC": "K04103",
    "acsA": "K14187",
    "budB": "K01652",
    "entA": "K02361",
    "entB": "K01252",
    "entC": "K02362",
    "phzF": "K18000",
}


@dataclass
class DRAMAnnotation:
    """A single gene annotation row from DRAM output."""

    gene_id: str
    mag_id: str
    ko_id: Optional[str] = None
    kegg_hit: Optional[str] = None
    cazy_id: Optional[str] = None
    pfam_id: Optional[str] = None


@dataclass(frozen=True)
class FunctionalTable:
    """Function-by-MAG presence/count matrix.

    Parameters
    ----------
    function_ids : list[str]
        Row labels (e.g. KO IDs, CAZy families, or PGPR trait names).
    mag_ids : list[str]
        Column labels (MAG identifiers).
    values : np.ndarray
        Matrix of shape ``(len(function_ids), len(mag_ids))``.
    function_type : str
        One of ``"ko"``, ``"module"``, ``"cazy"``, ``"pgpr"``.
    """

    function_ids: list[str]
    mag_ids: list[str]
    values: np.ndarray
    function_type: str

    def to_sample_abundance(self, abundance: AbundanceTable) -> AbundanceTable:
        """Project functional profiles into sample space.

        Computes ``F @ A`` where *F* is the function-by-MAG matrix and *A* is
        the MAG-by-sample abundance matrix, aligning on the intersection of
        MAG IDs present in both tables.

        Returns an :class:`AbundanceTable` whose ``mag_ids`` are
        ``function_ids`` (overloading the field for row labels).
        """
        # Build index maps for shared MAGs
        func_mag_idx = {m: i for i, m in enumerate(self.mag_ids)}
        abund_mag_idx = {m: i for i, m in enumerate(abundance.mag_ids)}

        shared_mags = [m for m in self.mag_ids if m in abund_mag_idx]
        if not shared_mags:
            # No overlap → zero matrix
            return AbundanceTable(
                mag_ids=list(self.function_ids),
                sample_ids=list(abundance.sample_ids),
                abundances=np.zeros(
                    (len(self.function_ids), abundance.n_samples), dtype=np.float64
                ),
            )

        func_cols = [func_mag_idx[m] for m in shared_mags]
        abund_rows = [abund_mag_idx[m] for m in shared_mags]

        F_sub = self.values[:, func_cols]  # (n_functions, n_shared)
        A_sub = abundance.abundances[abund_rows, :]  # (n_shared, n_samples)

        result = F_sub @ A_sub  # (n_functions, n_samples)

        return AbundanceTable(
            mag_ids=list(self.function_ids),
            sample_ids=list(abundance.sample_ids),
            abundances=result.astype(np.float64),
        )


# ---------------------------------------------------------------------------
# DRAM annotation loader
# ---------------------------------------------------------------------------

_GENE_ID_SUFFIX_RE = re.compile(r"_\d+$")
_CONTIG_GENE_SUFFIX_RE = re.compile(r"_k\d+_\d+$")


def _normalize_mag_id(mag_id: str) -> str:
    """Normalize MAG identifiers across tools (dot/underscore variants)."""
    return mag_id.strip().replace(".", "_")


def _mag_id_from_gene_id(gene_id: str) -> str:
    """Derive MAG ID from gene IDs emitted by DRAM/prodigal.

    Handles both patterns:
    - ``MAG_12345``
    - ``MAG_k141_1063138``
    """
    mag = _CONTIG_GENE_SUFFIX_RE.sub("", gene_id)
    mag = _GENE_ID_SUFFIX_RE.sub("", mag)
    return _normalize_mag_id(mag)


def load_dram_annotations(path: str | Path) -> list[DRAMAnnotation]:
    """Parse a DRAM ``annotations.tsv`` file.

    The first column is typically an unnamed index that duplicates the gene ID.
    We detect this by checking whether the first header field is empty or looks
    like an unnamed pandas index column.

    Returns a list of :class:`DRAMAnnotation` objects.
    """
    path = Path(path)
    annotations: list[DRAMAnnotation] = []

    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fieldnames = reader.fieldnames or []

        # Determine which column holds the gene ID.  DRAM often writes an
        # unnamed first column (header is "" or "Unnamed: 0") that contains
        # the gene ID, with a subsequent "gene_id" column that duplicates it.
        has_unnamed_index = (
            len(fieldnames) > 0 and fieldnames[0].strip() in ("", "Unnamed: 0")
        )

        for row in reader:
            if has_unnamed_index:
                gene_id = row.get(fieldnames[0], "").strip()
            else:
                gene_id = row.get("gene_id", row.get("", "")).strip()

            if not gene_id:
                continue

            # Prefer explicit MAG/bin identifier from DRAM when present.
            fasta_mag = row.get("fasta", "").strip()
            if fasta_mag:
                mag_id = _normalize_mag_id(fasta_mag)
            else:
                mag_id = _mag_id_from_gene_id(gene_id)

            # DRAM commonly uses `ko_id`; keep `kegg_id` fallback for compatibility.
            ko_id = (
                row.get("ko_id", "").strip()
                or row.get("kegg_id", "").strip()
                or None
            )
            kegg_hit = row.get("kegg_hit", "").strip() or None
            cazy_id = row.get("cazy_hits", "").strip() or None
            pfam_id = row.get("pfam_hits", "").strip() or None

            annotations.append(
                DRAMAnnotation(
                    gene_id=gene_id,
                    mag_id=mag_id,
                    ko_id=ko_id,
                    kegg_hit=kegg_hit,
                    cazy_id=cazy_id,
                    pfam_id=pfam_id,
                )
            )

    return annotations


# ---------------------------------------------------------------------------
# Functional table builders
# ---------------------------------------------------------------------------


def _build_ko_table(annotations: list[DRAMAnnotation]) -> FunctionalTable:
    """Binary KO presence matrix (1 if any gene in MAG has this KO)."""
    mag_set: dict[str, int] = {}
    ko_set: dict[str, int] = {}

    for ann in annotations:
        if ann.ko_id:
            mag_set.setdefault(ann.mag_id, len(mag_set))
            ko_set.setdefault(ann.ko_id, len(ko_set))

    mag_ids = sorted(mag_set, key=mag_set.get)  # type: ignore[arg-type]
    ko_ids = sorted(ko_set, key=ko_set.get)  # type: ignore[arg-type]

    mag_idx = {m: i for i, m in enumerate(mag_ids)}
    ko_idx = {k: i for i, k in enumerate(ko_ids)}

    values = np.zeros((len(ko_ids), len(mag_ids)), dtype=np.float64)
    for ann in annotations:
        if ann.ko_id and ann.ko_id in ko_idx:
            values[ko_idx[ann.ko_id], mag_idx[ann.mag_id]] = 1.0

    return FunctionalTable(
        function_ids=ko_ids,
        mag_ids=mag_ids,
        values=values,
        function_type="ko",
    )


def _build_cazy_table(annotations: list[DRAMAnnotation]) -> FunctionalTable:
    """Integer CAZy family count matrix."""
    mag_set: dict[str, int] = {}
    cazy_set: dict[str, int] = {}

    for ann in annotations:
        if ann.cazy_id:
            mag_set.setdefault(ann.mag_id, len(mag_set))
            cazy_set.setdefault(ann.cazy_id, len(cazy_set))

    mag_ids = sorted(mag_set, key=mag_set.get)  # type: ignore[arg-type]
    cazy_ids = sorted(cazy_set, key=cazy_set.get)  # type: ignore[arg-type]

    mag_idx = {m: i for i, m in enumerate(mag_ids)}
    cazy_idx = {c: i for i, c in enumerate(cazy_ids)}

    values = np.zeros((len(cazy_ids), len(mag_ids)), dtype=np.float64)
    for ann in annotations:
        if ann.cazy_id and ann.cazy_id in cazy_idx:
            values[cazy_idx[ann.cazy_id], mag_idx[ann.mag_id]] += 1.0

    return FunctionalTable(
        function_ids=cazy_ids,
        mag_ids=mag_ids,
        values=values,
        function_type="cazy",
    )


def _build_pgpr_table(annotations: list[DRAMAnnotation]) -> FunctionalTable:
    """Binary PGPR trait presence matrix based on :data:`PGPR_MARKERS`."""
    # Reverse lookup: KO → trait name(s)
    ko_to_trait: dict[str, str] = {ko: trait for trait, ko in PGPR_MARKERS.items()}

    mag_set: dict[str, int] = {}
    for ann in annotations:
        mag_set.setdefault(ann.mag_id, len(mag_set))

    mag_ids = sorted(mag_set, key=mag_set.get)  # type: ignore[arg-type]
    mag_idx = {m: i for i, m in enumerate(mag_ids)}

    trait_names = list(PGPR_MARKERS.keys())
    trait_idx = {t: i for i, t in enumerate(trait_names)}

    values = np.zeros((len(trait_names), len(mag_ids)), dtype=np.float64)
    for ann in annotations:
        if ann.ko_id and ann.ko_id in ko_to_trait:
            trait = ko_to_trait[ann.ko_id]
            values[trait_idx[trait], mag_idx[ann.mag_id]] = 1.0

    return FunctionalTable(
        function_ids=trait_names,
        mag_ids=mag_ids,
        values=values,
        function_type="pgpr",
    )


def build_all_functional_tables(
    annotations: list[DRAMAnnotation],
) -> dict[str, FunctionalTable]:
    """Build KO, CAZy, and PGPR tables in a single pass over annotations.

    Returns a dict with keys ``"ko"``, ``"cazy"``, ``"pgpr"``.
    """
    ko_to_trait: dict[str, str] = {ko: trait for trait, ko in PGPR_MARKERS.items()}

    # Collect data in single pass
    ko_mags: dict[str, set[str]] = {}
    cazy_counts: dict[str, dict[str, float]] = {}
    pgpr_mags: dict[str, set[str]] = {trait: set() for trait in PGPR_MARKERS}
    all_mags: dict[str, int] = {}

    for ann in annotations:
        all_mags.setdefault(ann.mag_id, len(all_mags))
        if ann.ko_id:
            ko_mags.setdefault(ann.ko_id, set()).add(ann.mag_id)
            if ann.ko_id in ko_to_trait:
                pgpr_mags[ko_to_trait[ann.ko_id]].add(ann.mag_id)
        if ann.cazy_id:
            cazy_counts.setdefault(ann.cazy_id, {})
            cazy_counts[ann.cazy_id][ann.mag_id] = (
                cazy_counts[ann.cazy_id].get(ann.mag_id, 0.0) + 1.0
            )

    # Build KO table
    ko_mag_set: dict[str, int] = {}
    for mags in ko_mags.values():
        for m in mags:
            ko_mag_set.setdefault(m, len(ko_mag_set))
    ko_mag_ids = sorted(ko_mag_set, key=ko_mag_set.get)  # type: ignore[arg-type]
    ko_ids = sorted(ko_mags, key=lambda k: next(
        i for i, ann in enumerate(annotations) if ann.ko_id == k
    ))
    # Simpler: preserve insertion order
    ko_ids = list(ko_mags.keys())
    ko_mag_idx = {m: i for i, m in enumerate(ko_mag_ids)}
    ko_func_idx = {k: i for i, k in enumerate(ko_ids)}
    ko_values = np.zeros((len(ko_ids), len(ko_mag_ids)), dtype=np.float64)
    for ko, mags in ko_mags.items():
        fi = ko_func_idx[ko]
        for m in mags:
            ko_values[fi, ko_mag_idx[m]] = 1.0
    ko_table = FunctionalTable(
        function_ids=ko_ids, mag_ids=ko_mag_ids,
        values=ko_values, function_type="ko",
    )

    # Build CAZy table
    cazy_mag_set: dict[str, int] = {}
    for mag_counts in cazy_counts.values():
        for m in mag_counts:
            cazy_mag_set.setdefault(m, len(cazy_mag_set))
    cazy_mag_ids = sorted(cazy_mag_set, key=cazy_mag_set.get)  # type: ignore[arg-type]
    cazy_ids = list(cazy_counts.keys())
    cazy_mag_idx = {m: i for i, m in enumerate(cazy_mag_ids)}
    cazy_func_idx = {c: i for i, c in enumerate(cazy_ids)}
    cazy_values = np.zeros((len(cazy_ids), len(cazy_mag_ids)), dtype=np.float64)
    for cazy, mag_counts in cazy_counts.items():
        fi = cazy_func_idx[cazy]
        for m, count in mag_counts.items():
            cazy_values[fi, cazy_mag_idx[m]] = count
    cazy_table = FunctionalTable(
        function_ids=cazy_ids, mag_ids=cazy_mag_ids,
        values=cazy_values, function_type="cazy",
    )

    # Build PGPR table
    all_mag_ids = sorted(all_mags, key=all_mags.get)  # type: ignore[arg-type]
    pgpr_mag_idx = {m: i for i, m in enumerate(all_mag_ids)}
    trait_names = list(PGPR_MARKERS.keys())
    pgpr_values = np.zeros((len(trait_names), len(all_mag_ids)), dtype=np.float64)
    for ti, trait in enumerate(trait_names):
        for m in pgpr_mags[trait]:
            pgpr_values[ti, pgpr_mag_idx[m]] = 1.0
    pgpr_table = FunctionalTable(
        function_ids=trait_names, mag_ids=all_mag_ids,
        values=pgpr_values, function_type="pgpr",
    )

    return {"ko": ko_table, "cazy": cazy_table, "pgpr": pgpr_table}


def build_functional_table(
    annotations: list[DRAMAnnotation],
    function_type: str,
) -> FunctionalTable:
    """Build a :class:`FunctionalTable` from DRAM annotations.

    Parameters
    ----------
    annotations
        Parsed DRAM annotation records.
    function_type
        ``"ko"`` for binary KO presence, ``"cazy"`` for CAZy family counts,
        or ``"pgpr"`` for PGPR marker trait presence.

    Raises
    ------
    ValueError
        If *function_type* is not one of the recognised types.
    """
    builders = {
        "ko": _build_ko_table,
        "cazy": _build_cazy_table,
        "pgpr": _build_pgpr_table,
    }
    if function_type not in builders:
        raise ValueError(
            f"Unknown function_type {function_type!r}; "
            f"expected one of {sorted(builders)}"
        )
    return builders[function_type](annotations)
