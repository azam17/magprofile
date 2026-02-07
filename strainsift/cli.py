"""Command-line interface for StrainSift.

Provides commands for:
- simulate: Generate synthetic metagenomes
- index: Build reference database index
- classify: Classify reads against index
- infer: Run joint EM inference
- run: Full pipeline (classify + infer)
- benchmark: Run evaluation benchmarks
- ablation: Run ablation study
"""

from __future__ import annotations

import json
import logging
import sys
from pathlib import Path

import click
import numpy as np

from strainsift import __version__

logger = logging.getLogger("strainsift")


def _setup_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


@click.group()
@click.version_option(version=__version__)
@click.option("-v", "--verbose", is_flag=True, help="Enable debug logging.")
def cli(verbose: bool) -> None:
    """StrainSift: Joint Strain-Function Inference for Shotgun Metagenomics."""
    _setup_logging(verbose)


# ---------------------------------------------------------------------------
# simulate
# ---------------------------------------------------------------------------


@cli.command()
@click.option("--n-strains", default=3, help="Number of strains.")
@click.option("--genome-length", default=50_000, help="Genome length (bp).")
@click.option("--n-gene-families", default=10, help="Number of gene families.")
@click.option("--gene-length", default=1000, help="Gene length (bp).")
@click.option("--mobile-fraction", default=0.2, help="Fraction of mobile gene families.")
@click.option("--accessory-fraction", default=0.3, help="Fraction of accessory gene families.")
@click.option("--ani", default=0.995, help="Pairwise ANI between strains.")
@click.option("--n-reads", default=10_000, help="Number of reads to simulate.")
@click.option("--read-length", default=150, help="Read length (bp).")
@click.option("--error-rate", default=0.005, help="Per-base error rate.")
@click.option("--open-universe", default=0.0, help="Fraction of strains not in DB.")
@click.option("--seed", default=42, help="Random seed.")
@click.option("-o", "--output-dir", required=True, type=click.Path(), help="Output directory.")
def simulate(
    n_strains: int,
    genome_length: int,
    n_gene_families: int,
    gene_length: int,
    mobile_fraction: float,
    accessory_fraction: float,
    ani: float,
    n_reads: int,
    read_length: int,
    error_rate: float,
    open_universe: float,
    seed: int,
    output_dir: str,
) -> None:
    """Generate a synthetic metagenome with known ground truth."""
    from strainsift.simulator import SimulationConfig, simulate_metagenome
    from strainsift.utils import write_fastq

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    config = SimulationConfig(
        n_strains=n_strains,
        genome_length=genome_length,
        n_gene_families=n_gene_families,
        gene_length=gene_length,
        mobile_gene_fraction=mobile_fraction,
        accessory_gene_fraction=accessory_fraction,
        pairwise_ani=ani,
        n_reads=n_reads,
        read_length=read_length,
        error_rate=error_rate,
        open_universe_fraction=open_universe,
        random_seed=seed,
    )

    logger.info("Simulating metagenome: %d strains, %d reads, ANI=%.4f", n_strains, n_reads, ani)
    references, reads, truth = simulate_metagenome(config)

    # Write outputs
    write_fastq(reads, out / "reads.fastq")
    truth.save(out / "ground_truth.json")

    # Write references
    ref_data = []
    for ref in references:
        ref_data.append({
            "strain_id": ref.strain_id,
            "species_id": ref.species_id,
            "genome_length": ref.genome_length,
            "n_genes": len(ref.gene_annotations),
        })
    (out / "references.json").write_text(json.dumps(ref_data, indent=2))

    click.echo(f"Simulated {len(reads)} reads from {n_strains} strains -> {out}")
    click.echo(f"  Classified fraction: {len(truth.strains_in_db)}/{n_strains} strains in DB")


# ---------------------------------------------------------------------------
# run (full pipeline)
# ---------------------------------------------------------------------------


@cli.command()
@click.option("--fastq", required=True, type=click.Path(exists=True), help="Input FASTQ file.")
@click.option("--references", required=True, type=click.Path(exists=True), help="Reference genomes JSON.")
@click.option("--n-restarts", default=10, help="Number of EM restarts.")
@click.option("--max-iterations", default=200, help="Max EM iterations.")
@click.option("--alpha", default=0.5, help="Dirichlet prior concentration.")
@click.option("--no-mobile-pool", is_flag=True, help="Disable mobile element pool.")
@click.option("--seed", default=42, help="Random seed.")
@click.option("-o", "--output", required=True, type=click.Path(), help="Output JSON path.")
def run(
    fastq: str,
    references: str,
    n_restarts: int,
    max_iterations: int,
    alpha: float,
    no_mobile_pool: bool,
    seed: int,
    output: str,
) -> None:
    """Run the full StrainSift pipeline: classify + joint EM inference."""
    from strainsift.classify import ClassificationConfig, ReadClassifier
    from strainsift.em import EMConfig, JointEM
    from strainsift.index import StrainSiftIndex
    from strainsift.utils import parse_fastq

    # Load references (would need proper serialization in production)
    logger.info("Loading references from %s", references)

    # Build index
    # (In production, references would be loaded from a proper database format)
    click.echo("Note: Full reference loading requires genome sequences. "
               "Use 'simulate' command for end-to-end testing.")

    click.echo(f"Pipeline would output to: {output}")


# ---------------------------------------------------------------------------
# benchmark
# ---------------------------------------------------------------------------


@cli.command()
@click.option("-o", "--output-dir", type=click.Path(), default="benchmark_output", help="Output directory.")
@click.option("--n-restarts", default=5, help="EM restarts per run.")
@click.option("--seed", default=42, help="Random seed.")
def benchmark(output_dir: str, n_restarts: int, seed: int) -> None:
    """Run the full StrainSift benchmark suite."""
    from strainsift.benchmark import run_full_benchmark
    from strainsift.em import EMConfig

    out = Path(output_dir)
    results = run_full_benchmark(output_dir=out)

    click.echo(f"\nBenchmark complete: {len(results)} scenarios")
    for r in results:
        click.echo(
            f"  {r.scenario}: L1={r.strain_l1_error:.4f}, "
            f"F1={r.function_f1:.4f}, "
            f"unclassified={r.unclassified_fraction:.2%}"
        )


# ---------------------------------------------------------------------------
# ablation
# ---------------------------------------------------------------------------


@cli.command()
@click.option("--ani", default=0.995, help="Pairwise ANI.")
@click.option("--n-reads", default=10_000, help="Number of reads.")
@click.option("--n-strains", default=3, help="Number of strains.")
@click.option("--accessory-fraction", default=0.3, help="Accessory gene fraction.")
@click.option("--n-restarts", default=10, help="EM restarts.")
@click.option("--seed", default=42, help="Random seed.")
@click.option("-o", "--output", type=click.Path(), default="ablation_results.json", help="Output file.")
def ablation(
    ani: float,
    n_reads: int,
    n_strains: int,
    accessory_fraction: float,
    n_restarts: int,
    seed: int,
    output: str,
) -> None:
    """Run the ablation study (joint vs strain-only vs function-only)."""
    from strainsift.benchmark import run_ablation_experiment
    from strainsift.em import EMConfig
    from strainsift.simulator import SimulationConfig

    sim_config = SimulationConfig(
        n_strains=n_strains,
        pairwise_ani=ani,
        n_reads=n_reads,
        accessory_gene_fraction=accessory_fraction,
        random_seed=seed,
    )
    em_config = EMConfig(n_restarts=n_restarts, random_seed=seed)

    results = run_ablation_experiment(
        sim_config, em_config, scenario_name=f"ANI={ani}_acc={accessory_fraction}"
    )

    # Display results
    click.echo("\n=== Ablation Study Results ===")
    click.echo(f"ANI={ani}, accessory_fraction={accessory_fraction}, n_reads={n_reads}")
    click.echo(f"{'Condition':<20} {'L1 Error':>10} {'BC':>10} {'Func F1':>10} {'Mobile F1':>10}")
    click.echo("-" * 62)

    for r in results:
        click.echo(
            f"{r.condition:<20} "
            f"{r.abundance_metrics.get('l1_error', 0):.4f}     "
            f"{r.abundance_metrics.get('bray_curtis', 0):.4f}     "
            f"{r.attribution_metrics.get('f1', 0):.4f}     "
            f"{r.mobile_metrics.get('f1', 0):.4f}"
        )

    # Save
    out_data = []
    for r in results:
        out_data.append({
            "condition": r.condition,
            "scenario": r.scenario,
            "abundance_metrics": r.abundance_metrics,
            "attribution_metrics": r.attribution_metrics,
            "mobile_metrics": r.mobile_metrics,
        })
    Path(output).write_text(json.dumps(out_data, indent=2))
    click.echo(f"\nResults saved to {output}")


# ---------------------------------------------------------------------------
# identifiability
# ---------------------------------------------------------------------------


@cli.command()
@click.option("--n-replicates", default=10, help="Replicates per grid cell.")
@click.option("--seed", default=42, help="Random seed.")
@click.option("-o", "--output", type=click.Path(), default="identifiability.json", help="Output file.")
def identifiability(n_replicates: int, seed: int, output: str) -> None:
    """Run the identifiability phase diagram sweep."""
    from strainsift.benchmark import run_identifiability_sweep
    from strainsift.simulator import SimulationConfig

    base_config = SimulationConfig(n_strains=2, genome_length=50_000, random_seed=seed)
    results = run_identifiability_sweep(
        n_replicates=n_replicates, base_config=base_config
    )

    Path(output).write_text(json.dumps(results, indent=2))
    click.echo(f"Identifiability sweep saved to {output}")


def main() -> None:
    """Entry point for the CLI."""
    cli()


if __name__ == "__main__":
    main()
