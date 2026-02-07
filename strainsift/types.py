"""Core data types for StrainSift.

Dataclasses for strain profiles, read assignments, classification results,
and EM outputs. All probabilistic outputs carry uncertainty estimates.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np


@dataclass
class GeneFamilyAnnotation:
    """A gene family annotation on a reference genome.

    Attributes:
        gene_family_id: Unique identifier for the gene family.
        start: Start position on the genome (0-based).
        end: End position on the genome (exclusive).
        is_mobile: Whether this gene is annotated as mobile (plasmid/phage/ICE).
    """

    gene_family_id: str
    start: int
    end: int
    is_mobile: bool = False


@dataclass
class StrainReference:
    """Reference genome for a single strain.

    Attributes:
        strain_id: Unique identifier for the strain.
        species_id: Species this strain belongs to.
        genome: DNA sequence (uppercase ACGT).
        gene_annotations: List of gene family annotations on this genome.
        metadata: Additional metadata (e.g., source, ANI to type strain).
    """

    strain_id: str
    species_id: str
    genome: str
    gene_annotations: list[GeneFamilyAnnotation] = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)

    @property
    def genome_length(self) -> int:
        return len(self.genome)

    def gene_at_position(self, pos: int) -> str | None:
        """Return gene family ID covering position, or None if intergenic."""
        for ann in self.gene_annotations:
            if ann.start <= pos < ann.end:
                return ann.gene_family_id
        return None


@dataclass
class ReadCandidate:
    """A candidate (strain, gene_family) assignment for a single read.

    Attributes:
        strain_id: Candidate strain identifier.
        gene_family_id: Gene family the read maps to (None if intergenic).
        containment: k-mer containment score (fraction of read k-mers found
            in the strain reference). Range [0, 1].
        k: The k-mer size used for containment calculation.
    """

    strain_id: str
    gene_family_id: str | None
    containment: float
    k: int = 31


@dataclass
class ReadClassification:
    """Classification result for a single read.

    Attributes:
        read_id: Read identifier from FASTQ header.
        candidates: List of candidate (strain, gene_family) pairs with scores.
        species_candidates: Species identified at coarse screening stage.
        is_classified: Whether any candidate exceeds minimum threshold.
    """

    read_id: str
    candidates: list[ReadCandidate] = field(default_factory=list)
    species_candidates: list[str] = field(default_factory=list)
    is_classified: bool = False


@dataclass
class ReadAssignment:
    """Posterior strain assignment for a single read (EM output).

    Attributes:
        read_id: Read identifier.
        posteriors: Dict mapping strain_id -> posterior probability.
        gene_family_id: Most likely gene family (None if intergenic).
        is_mobile: Whether assigned to the mobile element pool.
    """

    read_id: str
    posteriors: dict[str, float] = field(default_factory=dict)
    gene_family_id: str | None = None
    is_mobile: bool = False

    @property
    def map_strain(self) -> str:
        """Maximum a posteriori strain assignment."""
        return max(self.posteriors, key=self.posteriors.get)

    @property
    def map_probability(self) -> float:
        """Probability of MAP strain assignment."""
        return max(self.posteriors.values())


@dataclass
class ConvergenceDiagnostics:
    """Diagnostics for EM convergence.

    Attributes:
        log_likelihoods: Log-likelihood trace per iteration.
        n_restarts: Number of random restarts performed.
        restart_log_likelihoods: Final log-likelihood for each restart.
        best_restart: Index of the restart with highest log-likelihood.
        converged: Whether the best run converged.
        n_iterations: Number of iterations in the best run.
        multimodal_warning: True if restarts disagree substantially.
    """

    log_likelihoods: list[float] = field(default_factory=list)
    n_restarts: int = 1
    restart_log_likelihoods: list[float] = field(default_factory=list)
    best_restart: int = 0
    converged: bool = False
    n_iterations: int = 0
    multimodal_warning: bool = False


@dataclass
class StrainProfile:
    """Complete StrainSift output for one species in one sample.

    Attributes:
        species_id: Species identifier.
        strain_ids: Ordered list of strain identifiers (length K).
        pi: Strain abundance vector (length K+1, last is mobile pool).
        pi_ci_lower: Lower 95% CI for pi.
        pi_ci_upper: Upper 95% CI for pi.
        phi: Strain-function matrix (K x G), P(strain k carries gene g).
        phi_ci_lower: Lower 95% CI for phi.
        phi_ci_upper: Upper 95% CI for phi.
        gene_family_ids: Ordered list of gene family IDs (length G).
        mobile_pool_proportion: Fraction of reads in mobile pool (= pi[-1]).
        fisher_information: Fisher information matrix (K x K) at convergence.
        identifiability_scores: Per-strain-pair identifiability (from Fisher info).
        convergence: Convergence diagnostics.
        bic: Bayesian Information Criterion for this K.
        n_reads_classified: Number of reads contributing to inference.
        n_reads_total: Total reads assigned to this species at coarse stage.
        read_assignments: Per-read posterior assignments (optional, large).
    """

    species_id: str
    strain_ids: list[str]
    pi: np.ndarray  # (K+1,)
    pi_ci_lower: np.ndarray  # (K+1,)
    pi_ci_upper: np.ndarray  # (K+1,)
    phi: np.ndarray  # (K, G)
    phi_ci_lower: np.ndarray  # (K, G)
    phi_ci_upper: np.ndarray  # (K, G)
    gene_family_ids: list[str]
    mobile_pool_proportion: float
    fisher_information: np.ndarray  # (K, K)
    identifiability_scores: dict[tuple[str, str], float] = field(
        default_factory=dict
    )
    convergence: ConvergenceDiagnostics = field(
        default_factory=ConvergenceDiagnostics
    )
    bic: float = 0.0
    n_reads_classified: int = 0
    n_reads_total: int = 0
    read_assignments: list[ReadAssignment] | None = None

    def to_dict(self) -> dict[str, Any]:
        """Serialize to JSON-compatible dict (numpy arrays -> lists)."""
        return {
            "species_id": self.species_id,
            "strain_ids": self.strain_ids,
            "pi": self.pi.tolist(),
            "pi_ci_lower": self.pi_ci_lower.tolist(),
            "pi_ci_upper": self.pi_ci_upper.tolist(),
            "phi": self.phi.tolist(),
            "gene_family_ids": self.gene_family_ids,
            "mobile_pool_proportion": self.mobile_pool_proportion,
            "identifiability_scores": {
                f"{k[0]}__{k[1]}": v
                for k, v in self.identifiability_scores.items()
            },
            "bic": self.bic,
            "n_reads_classified": self.n_reads_classified,
            "n_reads_total": self.n_reads_total,
        }

    def save(self, path: Path) -> None:
        """Save profile to JSON file."""
        path.write_text(json.dumps(self.to_dict(), indent=2))


@dataclass
class SimulationGroundTruth:
    """Ground truth from the simulator for benchmarking.

    Attributes:
        true_pi: True strain abundance vector.
        true_phi: True strain-function matrix (K x G).
        strain_ids: Strain identifiers.
        gene_family_ids: Gene family identifiers.
        read_assignments: True strain assignment per read (read_id -> strain_id).
        read_gene_families: True gene family per read (read_id -> gene_family_id | None).
        mobile_genes: Set of gene family IDs that are mobile.
        strains_in_db: Set of strain IDs present in the reference database.
        strains_not_in_db: Set of strain IDs NOT in the reference database.
    """

    true_pi: np.ndarray
    true_phi: np.ndarray
    strain_ids: list[str]
    gene_family_ids: list[str]
    read_assignments: dict[str, str] = field(default_factory=dict)
    read_gene_families: dict[str, str | None] = field(default_factory=dict)
    mobile_genes: set[str] = field(default_factory=set)
    strains_in_db: set[str] = field(default_factory=set)
    strains_not_in_db: set[str] = field(default_factory=set)

    def to_dict(self) -> dict[str, Any]:
        """Serialize to JSON-compatible dict."""
        return {
            "true_pi": self.true_pi.tolist(),
            "true_phi": self.true_phi.tolist(),
            "strain_ids": self.strain_ids,
            "gene_family_ids": self.gene_family_ids,
            "read_assignments": self.read_assignments,
            "read_gene_families": {
                k: v for k, v in self.read_gene_families.items()
            },
            "mobile_genes": list(self.mobile_genes),
            "strains_in_db": list(self.strains_in_db),
            "strains_not_in_db": list(self.strains_not_in_db),
        }

    def save(self, path: Path) -> None:
        """Save ground truth to JSON file."""
        path.write_text(json.dumps(self.to_dict(), indent=2))


@dataclass
class BenchmarkResult:
    """Results from a benchmarking comparison.

    Attributes:
        method: Name of the method being evaluated.
        scenario: Description of the benchmark scenario.
        strain_l1_error: L1 error in strain abundance estimation.
        strain_bray_curtis: Bray-Curtis dissimilarity for strain abundances.
        function_precision: Precision for function-to-strain attribution.
        function_recall: Recall for function-to-strain attribution.
        function_f1: F1 score for function-to-strain attribution.
        mobile_precision: Precision for mobile element detection.
        mobile_recall: Recall for mobile element detection.
        unclassified_fraction: Fraction of reads not classified.
        runtime_seconds: Wall-clock runtime.
        peak_memory_mb: Peak memory usage in MB.
        metadata: Additional result metadata.
    """

    method: str
    scenario: str
    strain_l1_error: float = 0.0
    strain_bray_curtis: float = 0.0
    function_precision: float = 0.0
    function_recall: float = 0.0
    function_f1: float = 0.0
    mobile_precision: float = 0.0
    mobile_recall: float = 0.0
    unclassified_fraction: float = 0.0
    runtime_seconds: float = 0.0
    peak_memory_mb: float = 0.0
    metadata: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        """Serialize to dict."""
        return {
            "method": self.method,
            "scenario": self.scenario,
            "strain_l1_error": self.strain_l1_error,
            "strain_bray_curtis": self.strain_bray_curtis,
            "function_precision": self.function_precision,
            "function_recall": self.function_recall,
            "function_f1": self.function_f1,
            "mobile_precision": self.mobile_precision,
            "mobile_recall": self.mobile_recall,
            "unclassified_fraction": self.unclassified_fraction,
            "runtime_seconds": self.runtime_seconds,
            "peak_memory_mb": self.peak_memory_mb,
            "metadata": self.metadata,
        }
