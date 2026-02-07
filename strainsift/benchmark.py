"""Evaluation framework and metrics for StrainSift benchmarking.

Implements all metrics and comparison workflows needed for:
- Ablation study (joint vs strain-only vs function-only)
- Identifiability phase diagram
- Head-to-head comparisons
- Mock community validation
"""

from __future__ import annotations

import json
import logging
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

from strainsift.classify import (
    ClassificationConfig,
    ReadClassifier,
    classify_reads_batch,
)
from strainsift.em import EMConfig, JointEM
from strainsift.index import StrainSiftIndex
from strainsift.simulator import SimulationConfig, simulate_metagenome
from strainsift.types import (
    BenchmarkResult,
    ReadClassification,
    SimulationGroundTruth,
    StrainProfile,
    StrainReference,
)
from strainsift.utils import FastqRecord, bray_curtis, l1_error

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------


def strain_abundance_metrics(
    estimated_pi: np.ndarray,
    true_pi: np.ndarray,
    strain_ids: list[str],
    true_strain_ids: list[str],
) -> dict[str, float]:
    """Compute strain abundance estimation metrics.

    Aligns estimated and true strain vectors by strain ID before comparison.

    Returns dict with: l1_error, bray_curtis, rmse, pearson_r, max_abs_error.
    """
    # Build aligned vectors
    all_ids = sorted(set(strain_ids) | set(true_strain_ids))
    est_map = dict(zip(strain_ids, estimated_pi[: len(strain_ids)]))
    true_map = dict(zip(true_strain_ids, true_pi[: len(true_strain_ids)]))

    est_aligned = np.array([est_map.get(sid, 0.0) for sid in all_ids])
    true_aligned = np.array([true_map.get(sid, 0.0) for sid in all_ids])

    # Normalize in case of mobile pool inclusion
    if est_aligned.sum() > 0:
        est_aligned = est_aligned / est_aligned.sum()
    if true_aligned.sum() > 0:
        true_aligned = true_aligned / true_aligned.sum()

    metrics = {
        "l1_error": l1_error(est_aligned, true_aligned),
        "bray_curtis": bray_curtis(est_aligned, true_aligned),
        "rmse": float(np.sqrt(np.mean((est_aligned - true_aligned) ** 2))),
        "max_abs_error": float(np.max(np.abs(est_aligned - true_aligned))),
    }

    # Pearson correlation (requires >= 2 non-zero entries)
    if len(all_ids) >= 2 and np.std(est_aligned) > 0 and np.std(true_aligned) > 0:
        metrics["pearson_r"] = float(np.corrcoef(est_aligned, true_aligned)[0, 1])
    else:
        metrics["pearson_r"] = 0.0

    return metrics


def function_attribution_metrics(
    estimated_phi: np.ndarray,
    true_phi: np.ndarray,
    strain_ids: list[str],
    true_strain_ids: list[str],
    gene_family_ids: list[str],
    true_gene_family_ids: list[str],
    threshold: float = 0.5,
) -> dict[str, float]:
    """Compute function-to-strain attribution metrics.

    Binarizes phi at threshold and computes precision/recall/F1 on the
    resulting strain×gene presence/absence matrix.

    Returns dict with: precision, recall, f1, accuracy.
    """
    # Align strains and genes
    common_strains = sorted(set(strain_ids) & set(true_strain_ids))
    common_genes = sorted(set(gene_family_ids) & set(true_gene_family_ids))

    if not common_strains or not common_genes:
        return {"precision": 0.0, "recall": 0.0, "f1": 0.0, "accuracy": 0.0}

    # Index maps
    est_strain_idx = {s: i for i, s in enumerate(strain_ids)}
    true_strain_idx = {s: i for i, s in enumerate(true_strain_ids)}
    est_gene_idx = {g: i for i, g in enumerate(gene_family_ids)}
    true_gene_idx = {g: i for i, g in enumerate(true_gene_family_ids)}

    # Extract aligned submatrices
    est_rows = [est_strain_idx[s] for s in common_strains]
    est_cols = [est_gene_idx[g] for g in common_genes]
    true_rows = [true_strain_idx[s] for s in common_strains]
    true_cols = [true_gene_idx[g] for g in common_genes]

    est_sub = estimated_phi[np.ix_(est_rows, est_cols)]
    true_sub = true_phi[np.ix_(true_rows, true_cols)]

    # Binarize
    est_binary = (est_sub >= threshold).astype(int)
    true_binary = (true_sub >= threshold).astype(int)

    tp = np.sum((est_binary == 1) & (true_binary == 1))
    fp = np.sum((est_binary == 1) & (true_binary == 0))
    fn = np.sum((est_binary == 0) & (true_binary == 1))
    tn = np.sum((est_binary == 0) & (true_binary == 0))

    precision = float(tp / (tp + fp)) if (tp + fp) > 0 else 0.0
    recall = float(tp / (tp + fn)) if (tp + fn) > 0 else 0.0
    f1 = (
        2 * precision * recall / (precision + recall)
        if (precision + recall) > 0
        else 0.0
    )
    accuracy = float((tp + tn) / (tp + fp + fn + tn)) if (tp + fp + fn + tn) > 0 else 0.0

    return {
        "precision": precision,
        "recall": recall,
        "f1": f1,
        "accuracy": accuracy,
    }


def mobile_element_metrics(
    read_assignments: dict[str, bool],
    true_mobile_reads: set[str],
) -> dict[str, float]:
    """Compute mobile element detection precision/recall.

    Args:
        read_assignments: Dict mapping read_id -> is_mobile (from EM).
        true_mobile_reads: Set of read IDs truly from mobile elements.

    Returns dict with: precision, recall, f1.
    """
    predicted_mobile = {rid for rid, is_mob in read_assignments.items() if is_mob}

    tp = len(predicted_mobile & true_mobile_reads)
    fp = len(predicted_mobile - true_mobile_reads)
    fn = len(true_mobile_reads - predicted_mobile)

    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1 = (
        2 * precision * recall / (precision + recall)
        if (precision + recall) > 0
        else 0.0
    )

    return {"precision": precision, "recall": recall, "f1": f1}


# ---------------------------------------------------------------------------
# Ablation study
# ---------------------------------------------------------------------------


@dataclass
class AblationResult:
    """Result of a single ablation experiment run.

    Attributes:
        condition: One of 'joint', 'strain_only', 'function_only'.
        scenario: Description of the test scenario.
        abundance_metrics: Strain abundance metrics.
        attribution_metrics: Function attribution metrics.
        mobile_metrics: Mobile element metrics.
        profile: The full StrainProfile from EM.
        config: EMConfig used.
    """

    condition: str
    scenario: str
    abundance_metrics: dict[str, float] = field(default_factory=dict)
    attribution_metrics: dict[str, float] = field(default_factory=dict)
    mobile_metrics: dict[str, float] = field(default_factory=dict)
    profile: StrainProfile | None = None
    config: EMConfig | None = None


def run_ablation_experiment(
    sim_config: SimulationConfig,
    em_config: EMConfig | None = None,
    scenario_name: str = "default",
) -> list[AblationResult]:
    """Run the ablation study: joint vs strain-only vs function-only.

    This is THE critical experiment for StrainSift.

    Three conditions:
    1. Joint: full model (fix_phi=False, fix_pi=False)
    2. Strain-only: phi fixed from reference (fix_phi=True)
    3. Function-only: pi uniform, no strain feedback (fix_pi=True)

    Args:
        sim_config: Configuration for synthetic metagenome.
        em_config: Base EM configuration (ablation flags will be overridden).
        scenario_name: Name for this scenario.

    Returns:
        List of 3 AblationResult objects.
    """
    em_config = em_config or EMConfig()

    # Generate synthetic data
    references, reads, truth = simulate_metagenome(sim_config)

    # Build index and classify reads
    classifications, stats = classify_reads_batch(reads, references)

    # Get gene families and strain IDs from references
    strain_ids = [r.strain_id for r in references]
    gene_families = set()
    for ref in references:
        for ann in ref.gene_annotations:
            gene_families.add(ann.gene_family_id)
    gene_family_ids = sorted(gene_families)

    # Build reference phi from annotations
    reference_phi = _build_reference_phi(references, strain_ids, gene_family_ids)

    # Identify true mobile reads
    true_mobile_reads = {
        rid
        for rid, gf in truth.read_gene_families.items()
        if gf in truth.mobile_genes
    }

    results = []

    conditions = [
        ("joint", False, False),
        ("strain_only", True, False),
        ("function_only", False, True),
    ]

    for condition_name, fix_phi, fix_pi in conditions:
        cfg = EMConfig(
            max_iterations=em_config.max_iterations,
            convergence_threshold=em_config.convergence_threshold,
            n_restarts=em_config.n_restarts,
            alpha=em_config.alpha,
            beta_a=em_config.beta_a,
            beta_b=em_config.beta_b,
            min_containment=em_config.min_containment,
            mobile_pool=em_config.mobile_pool,
            fix_phi=fix_phi,
            fix_pi=fix_pi,
            random_seed=em_config.random_seed,
        )

        em = JointEM(config=cfg)
        profile = em.fit(
            read_classifications=classifications,
            strain_ids=strain_ids,
            gene_family_ids=gene_family_ids,
            reference_phi=reference_phi,
        )

        # Compute metrics
        abund_metrics = strain_abundance_metrics(
            profile.pi,
            truth.true_pi,
            profile.strain_ids,
            truth.strain_ids,
        )

        attr_metrics = function_attribution_metrics(
            profile.phi,
            truth.true_phi,
            profile.strain_ids,
            truth.strain_ids,
            profile.gene_family_ids,
            truth.gene_family_ids,
        )

        # Mobile metrics from read assignments
        mob_assignments = {}
        if profile.read_assignments:
            mob_assignments = {
                ra.read_id: ra.is_mobile for ra in profile.read_assignments
            }
        mob_metrics = mobile_element_metrics(mob_assignments, true_mobile_reads)

        result = AblationResult(
            condition=condition_name,
            scenario=scenario_name,
            abundance_metrics=abund_metrics,
            attribution_metrics=attr_metrics,
            mobile_metrics=mob_metrics,
            profile=profile,
            config=cfg,
        )
        results.append(result)

        logger.info(
            "Ablation [%s/%s]: L1=%.4f, F1=%.4f, mobile_F1=%.4f",
            scenario_name,
            condition_name,
            abund_metrics["l1_error"],
            attr_metrics["f1"],
            mob_metrics["f1"],
        )

    return results


# ---------------------------------------------------------------------------
# Identifiability phase diagram
# ---------------------------------------------------------------------------


def run_identifiability_sweep(
    ani_values: list[float] | None = None,
    coverage_values: list[float] | None = None,
    accessory_divergence_values: list[float] | None = None,
    n_replicates: int = 10,
    base_config: SimulationConfig | None = None,
    em_config: EMConfig | None = None,
) -> dict[str, Any]:
    """Sweep ANI × coverage × accessory divergence for identifiability phase diagram.

    Returns a nested dict structure suitable for heatmap plotting:
    {
        "ani_values": [...],
        "coverage_values": [...],
        "accessory_values": [...],
        "rmse": 3D array [ani][cov][acc],
        "fisher_min_eigenvalue": 3D array,
        "merge_fraction": 3D array (fraction of runs where Fisher merge was triggered),
    }
    """
    if ani_values is None:
        ani_values = [0.990, 0.995, 0.998, 0.999, 0.9995, 0.9999]
    if coverage_values is None:
        coverage_values = [0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0]
    if accessory_divergence_values is None:
        accessory_divergence_values = [0.0, 0.05, 0.10, 0.15, 0.20]

    base_config = base_config or SimulationConfig(n_strains=2, genome_length=50_000)
    em_config = em_config or EMConfig(n_restarts=5)

    n_ani = len(ani_values)
    n_cov = len(coverage_values)
    n_acc = len(accessory_divergence_values)

    rmse_grid = np.zeros((n_ani, n_cov, n_acc))
    fisher_grid = np.zeros((n_ani, n_cov, n_acc))
    merge_grid = np.zeros((n_ani, n_cov, n_acc))

    for i, ani in enumerate(ani_values):
        for j, cov in enumerate(coverage_values):
            for m, acc in enumerate(accessory_divergence_values):
                rmse_vals = []
                fisher_vals = []
                merge_count = 0

                for rep in range(n_replicates):
                    # Compute n_reads from coverage
                    read_len = base_config.read_length
                    genome_len = base_config.genome_length
                    n_reads = int(cov * genome_len / read_len)

                    cfg = SimulationConfig(
                        n_strains=base_config.n_strains,
                        genome_length=genome_len,
                        n_gene_families=base_config.n_gene_families,
                        gene_length=base_config.gene_length,
                        mobile_gene_fraction=base_config.mobile_gene_fraction,
                        accessory_gene_fraction=acc,
                        pairwise_ani=ani,
                        n_reads=max(n_reads, 10),
                        read_length=read_len,
                        error_rate=base_config.error_rate,
                        random_seed=base_config.random_seed + rep,
                    )

                    try:
                        references, reads, truth = simulate_metagenome(cfg)
                        classifications, _ = classify_reads_batch(reads, references)

                        strain_ids = [r.strain_id for r in references]
                        gene_families = set()
                        for ref in references:
                            for ann in ref.gene_annotations:
                                gene_families.add(ann.gene_family_id)
                        gene_family_ids = sorted(gene_families)
                        reference_phi = _build_reference_phi(
                            references, strain_ids, gene_family_ids
                        )

                        em = JointEM(config=em_config)
                        profile = em.fit(
                            classifications, strain_ids, gene_family_ids, reference_phi
                        )

                        # Compute RMSE
                        abund = strain_abundance_metrics(
                            profile.pi, truth.true_pi, strain_ids, truth.strain_ids
                        )
                        rmse_vals.append(abund["rmse"])

                        # Fisher info minimum eigenvalue
                        eigvals = np.linalg.eigvalsh(profile.fisher_information)
                        fisher_vals.append(float(np.min(eigvals)))

                        # Check if merge was triggered
                        if len(profile.strain_ids) < len(strain_ids):
                            merge_count += 1

                    except Exception as e:
                        logger.warning(
                            "Replicate failed (ANI=%.4f, cov=%.1f, acc=%.2f, rep=%d): %s",
                            ani, cov, acc, rep, e,
                        )
                        continue

                rmse_grid[i, j, m] = np.mean(rmse_vals) if rmse_vals else np.nan
                fisher_grid[i, j, m] = np.mean(fisher_vals) if fisher_vals else np.nan
                merge_grid[i, j, m] = merge_count / max(len(rmse_vals), 1)

    return {
        "ani_values": ani_values,
        "coverage_values": coverage_values,
        "accessory_values": accessory_divergence_values,
        "rmse": rmse_grid.tolist(),
        "fisher_min_eigenvalue": fisher_grid.tolist(),
        "merge_fraction": merge_grid.tolist(),
    }


# ---------------------------------------------------------------------------
# Full benchmark runner
# ---------------------------------------------------------------------------


def run_full_benchmark(
    scenarios: list[dict[str, Any]] | None = None,
    output_dir: Path | None = None,
) -> list[BenchmarkResult]:
    """Run the full StrainSift benchmark suite.

    Scenarios define simulation parameters and expected outcomes.
    If None, uses the default scenario set.

    Returns list of BenchmarkResult objects.
    """
    if scenarios is None:
        scenarios = _default_scenarios()

    results = []

    for scenario in scenarios:
        name = scenario.get("name", "unnamed")
        sim_config = SimulationConfig(**scenario.get("sim_params", {}))
        em_config = EMConfig(**scenario.get("em_params", {}))

        logger.info("Running benchmark scenario: %s", name)
        t0 = time.time()

        # Simulate
        references, reads, truth = simulate_metagenome(sim_config)

        # Classify
        classifications, class_stats = classify_reads_batch(reads, references)

        # EM inference
        strain_ids = [r.strain_id for r in references]
        gene_families = set()
        for ref in references:
            for ann in ref.gene_annotations:
                gene_families.add(ann.gene_family_id)
        gene_family_ids = sorted(gene_families)
        reference_phi = _build_reference_phi(references, strain_ids, gene_family_ids)

        em = JointEM(config=em_config)
        profile = em.fit(classifications, strain_ids, gene_family_ids, reference_phi)

        runtime = time.time() - t0

        # Metrics
        abund = strain_abundance_metrics(
            profile.pi, truth.true_pi, strain_ids, truth.strain_ids
        )
        attr = function_attribution_metrics(
            profile.phi, truth.true_phi,
            strain_ids, truth.strain_ids,
            gene_family_ids, truth.gene_family_ids,
        )

        true_mobile_reads = {
            rid for rid, gf in truth.read_gene_families.items()
            if gf in truth.mobile_genes
        }
        mob_assignments = {}
        if profile.read_assignments:
            mob_assignments = {ra.read_id: ra.is_mobile for ra in profile.read_assignments}
        mob = mobile_element_metrics(mob_assignments, true_mobile_reads)

        result = BenchmarkResult(
            method="strainsift_joint",
            scenario=name,
            strain_l1_error=abund["l1_error"],
            strain_bray_curtis=abund["bray_curtis"],
            function_precision=attr["precision"],
            function_recall=attr["recall"],
            function_f1=attr["f1"],
            mobile_precision=mob["precision"],
            mobile_recall=mob["recall"],
            unclassified_fraction=class_stats.unclassified_fraction,
            runtime_seconds=runtime,
            metadata={
                "abundance_metrics": abund,
                "attribution_metrics": attr,
                "convergence": {
                    "converged": profile.convergence.converged,
                    "n_iterations": profile.convergence.n_iterations,
                    "multimodal": profile.convergence.multimodal_warning,
                },
            },
        )
        results.append(result)

        logger.info(
            "Scenario %s: L1=%.4f, BC=%.4f, F1=%.4f, runtime=%.1fs",
            name, abund["l1_error"], abund["bray_curtis"], attr["f1"], runtime,
        )

    if output_dir:
        output_dir.mkdir(parents=True, exist_ok=True)
        out_path = output_dir / "benchmark_results.json"
        out_path.write_text(
            json.dumps([r.to_dict() for r in results], indent=2)
        )
        logger.info("Results saved to %s", out_path)

    return results


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _build_reference_phi(
    references: list[StrainReference],
    strain_ids: list[str],
    gene_family_ids: list[str],
) -> np.ndarray:
    """Build reference phi matrix from strain annotations.

    phi[k][g] = 1.0 if strain k's reference genome carries gene family g, else 0.0.
    """
    K = len(strain_ids)
    G = len(gene_family_ids)
    phi = np.zeros((K, G))

    strain_idx = {sid: i for i, sid in enumerate(strain_ids)}
    gene_idx = {gid: i for i, gid in enumerate(gene_family_ids)}

    for ref in references:
        k = strain_idx.get(ref.strain_id)
        if k is None:
            continue
        for ann in ref.gene_annotations:
            g = gene_idx.get(ann.gene_family_id)
            if g is not None:
                phi[k, g] = 1.0

    return phi


def _default_scenarios() -> list[dict[str, Any]]:
    """Default benchmark scenarios covering key evaluation axes."""
    return [
        {
            "name": "easy_3_strains",
            "sim_params": {
                "n_strains": 3,
                "pairwise_ani": 0.99,
                "n_reads": 10_000,
                "accessory_gene_fraction": 0.3,
                "random_seed": 42,
            },
            "em_params": {"n_restarts": 5, "random_seed": 42},
        },
        {
            "name": "hard_high_ani",
            "sim_params": {
                "n_strains": 3,
                "pairwise_ani": 0.999,
                "n_reads": 10_000,
                "accessory_gene_fraction": 0.3,
                "random_seed": 43,
            },
            "em_params": {"n_restarts": 10, "random_seed": 43},
        },
        {
            "name": "high_ani_with_accessory",
            "sim_params": {
                "n_strains": 3,
                "pairwise_ani": 0.995,
                "n_reads": 20_000,
                "accessory_gene_fraction": 0.4,
                "mobile_gene_fraction": 0.2,
                "random_seed": 44,
            },
            "em_params": {"n_restarts": 10, "random_seed": 44},
        },
        {
            "name": "low_coverage",
            "sim_params": {
                "n_strains": 2,
                "pairwise_ani": 0.995,
                "n_reads": 1_000,
                "accessory_gene_fraction": 0.2,
                "random_seed": 45,
            },
            "em_params": {"n_restarts": 10, "random_seed": 45},
        },
        {
            "name": "open_universe",
            "sim_params": {
                "n_strains": 4,
                "pairwise_ani": 0.995,
                "n_reads": 10_000,
                "accessory_gene_fraction": 0.3,
                "open_universe_fraction": 0.25,
                "random_seed": 46,
            },
            "em_params": {"n_restarts": 10, "random_seed": 46},
        },
    ]
