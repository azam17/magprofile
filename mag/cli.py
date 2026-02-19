"""Click CLI for magprofile — genome-resolved community profiling."""

from __future__ import annotations

from pathlib import Path

import click

from .io import load_abundance_table, load_metadata, load_taxonomy


@click.group()
def main() -> None:
    """magprofile — Genome-Resolved Community Profiling Tool."""


@main.command()
@click.option("--abundance", "-a", required=True, type=click.Path(exists=True), help="Abundance table TSV")
@click.option("--metadata", "-m", required=True, type=click.Path(exists=True), help="Sample metadata TSV")
@click.option("--group", "-g", default="compartment", help="Metadata variable for grouping")
@click.option("--output", "-o", default="results", help="Output directory")
def diversity(abundance: str, metadata: str, group: str, output: str) -> None:
    """Compute alpha diversity metrics."""
    from .diversity import compute_alpha_diversity
    from .plots import plot_alpha_diversity
    from .report import _write_alpha_csv

    table = load_abundance_table(abundance)
    meta = load_metadata(metadata)
    out = Path(output)
    out.mkdir(parents=True, exist_ok=True)

    result = compute_alpha_diversity(table)
    _write_alpha_csv(result, out / "alpha_diversity.csv")
    plot_alpha_diversity(result, meta, group, out / "alpha_diversity.pdf")
    click.echo(f"Alpha diversity results written to {out}/")


@main.command()
@click.option("--abundance", "-a", required=True, type=click.Path(exists=True), help="Abundance table TSV")
@click.option("--metadata", "-m", required=True, type=click.Path(exists=True), help="Sample metadata TSV")
@click.option("--group", "-g", default="compartment", help="Metadata variable for grouping")
@click.option("--method", type=click.Choice(["pcoa", "nmds"]), default="pcoa", help="Ordination method")
@click.option("--output", "-o", default="results", help="Output directory")
def ordination(abundance: str, metadata: str, group: str, method: str, output: str) -> None:
    """Run ordination analysis (PCoA or NMDS)."""
    from .beta import bray_curtis
    from .ordination import nmds as run_nmds
    from .ordination import pcoa as run_pcoa
    from .plots import plot_ordination
    from .report import _write_ordination

    table = load_abundance_table(abundance)
    meta = load_metadata(metadata)
    out = Path(output)
    out.mkdir(parents=True, exist_ok=True)

    bc = bray_curtis(table)
    if method == "pcoa":
        result = run_pcoa(bc)
    else:
        result = run_nmds(bc)

    _write_ordination(result, out / "ordination_coordinates.csv")
    plot_ordination(result, meta, group, out / f"ordination_{method}.pdf")
    click.echo(f"Ordination results written to {out}/")


@main.command()
@click.option("--abundance", "-a", required=True, type=click.Path(exists=True), help="Abundance table TSV")
@click.option("--metadata", "-m", required=True, type=click.Path(exists=True), help="Sample metadata TSV")
@click.option("--group", "-g", default="compartment", help="Metadata variable for grouping")
@click.option("--group1", required=True, help="First group for comparison")
@click.option("--group2", required=True, help="Second group for comparison")
@click.option("--output", "-o", default="results", help="Output directory")
def differential(abundance: str, metadata: str, group: str, group1: str, group2: str, output: str) -> None:
    """Differential abundance analysis between two groups."""
    from .differential import differential_abundance
    from .plots import plot_volcano
    from .report import _write_differential

    table = load_abundance_table(abundance)
    meta = load_metadata(metadata)
    out = Path(output)
    out.mkdir(parents=True, exist_ok=True)

    result = differential_abundance(table, meta, group, group1, group2)
    _write_differential(result, group1, group2, out / f"differential_{group1}_vs_{group2}.csv")
    plot_volcano(result, group1, group2, output=out / f"volcano_{group1}_vs_{group2}.pdf")

    n_sig = int((result.q_values < 0.05).sum())
    click.echo(f"Found {n_sig} significant MAGs (FDR < 0.05)")
    click.echo(f"Results written to {out}/")


@main.command()
@click.option("--abundance", "-a", required=True, type=click.Path(exists=True), help="Abundance table TSV")
@click.option("--metadata", "-m", required=True, type=click.Path(exists=True), help="Sample metadata TSV")
@click.option("--group", "-g", default="compartment", help="Metadata variable for grouping")
@click.option("--output", "-o", default="results", help="Output directory")
def indicator(abundance: str, metadata: str, group: str, output: str) -> None:
    """Indicator species analysis."""
    from .indicator import indicator_species
    from .plots import plot_indicator_species
    from .report import _write_indicator

    table = load_abundance_table(abundance)
    meta = load_metadata(metadata)
    out = Path(output)
    out.mkdir(parents=True, exist_ok=True)

    result = indicator_species(table, meta, group)
    _write_indicator(result, out / "indicator_species.csv")
    plot_indicator_species(result, output=out / "indicator_species.pdf")
    click.echo(f"Indicator species results written to {out}/")


@main.command()
@click.option("--abundance", "-a", required=True, type=click.Path(exists=True), help="Abundance table TSV")
@click.option("--taxonomy", "-t", default=None, type=click.Path(exists=True), help="Taxonomy TSV (optional)")
@click.option("--metadata", "-m", required=True, type=click.Path(exists=True), help="Sample metadata TSV")
@click.option("--group", "-g", default="compartment", help="Metadata variable for grouping")
@click.option("--permutations", "-p", default=999, help="Number of permutations for tests")
@click.option("--min-prevalence", default=0.1, type=float, help="Min fraction of samples a MAG must be present in (default 0.1)")
@click.option("--output", "-o", default="results", help="Output directory")
def report(abundance: str, taxonomy: str | None, metadata: str, group: str, permutations: int, min_prevalence: float, output: str) -> None:
    """Run full analysis pipeline and generate report."""
    from .report import generate_report

    table = load_abundance_table(abundance)
    tax = load_taxonomy(taxonomy) if taxonomy else None
    meta = load_metadata(metadata)

    generate_report(table, tax, meta, group, output, n_permutations=permutations, min_prevalence=min_prevalence)
    click.echo(f"Full report written to {output}/")


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


@main.command("net-report")
@click.option("--abundance", "-a", required=True, type=click.Path(exists=True), help="Abundance table TSV")
@click.option("--taxonomy", "-t", default=None, type=click.Path(exists=True), help="Taxonomy TSV (optional)")
@click.option("--metadata", "-m", required=True, type=click.Path(exists=True), help="Sample metadata TSV")
@click.option("--group", "-g", default="compartment", help="Metadata variable for grouping")
@click.option("--threshold", default=5.0, type=float, help="Phi percentile threshold (used for threshold-mode=group/global calibration)")
@click.option("--threshold-mode", type=click.Choice(["global", "group", "fixed"]), default="global", show_default=True, help="How to set group thresholds")
@click.option("--phi-threshold", default=None, type=float, help="Absolute phi threshold (required for threshold-mode=fixed; optional override for global)")
@click.option("--min-prevalence", default=0.5, type=float, help="Min prevalence for global network (default 0.5)")
@click.option("--group-min-prevalence", default=None, type=float, help="Min prevalence for group networks (default: same as --min-prevalence)")
@click.option("--output", "-o", default="net_results", help="Output directory")
def net_report(abundance: str, taxonomy: str | None, metadata: str, group: str, threshold: float, threshold_mode: str, phi_threshold: float | None, min_prevalence: float, group_min_prevalence: float | None, output: str) -> None:
    """Run network analysis pipeline (magnet)."""
    from .net_report import generate_net_report

    table = load_abundance_table(abundance)
    tax = load_taxonomy(taxonomy) if taxonomy else None
    meta = load_metadata(metadata)

    generate_net_report(
        table,
        tax,
        meta,
        group,
        output,
        threshold_percentile=threshold,
        min_prevalence=min_prevalence,
        threshold_mode=threshold_mode,
        phi_threshold=phi_threshold,
        group_min_prevalence=group_min_prevalence,
    )
    click.echo(f"Network analysis report written to {output}/")
