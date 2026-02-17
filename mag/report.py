"""Full analysis report generator for MAG community profiling."""

from __future__ import annotations

import csv
import logging
from pathlib import Path

import numpy as np

from . import beta as beta_mod
from . import compartment_specificity as cs_mod
from . import core_microbiome as cm_mod
from . import differential as diff_mod
from . import diversity as div_mod
from . import indicator as ind_mod
from . import ordination as ord_mod
from . import plots
from . import stats as stats_mod
from .html_report import generate_html_report
from .io import AbundanceTable, SampleMetadata, TaxonomyTable

logger = logging.getLogger(__name__)


def generate_report(
    abundance: AbundanceTable,
    taxonomy: TaxonomyTable | None,
    metadata: SampleMetadata,
    grouping_var: str,
    output_dir: str | Path,
    n_permutations: int = 999,
    min_prevalence: float = 0.1,
) -> None:
    """Run all analyses and write results to output directory."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    # Prevalence filtering
    n_before = abundance.n_mags
    if min_prevalence > 0:
        abundance = abundance.filter_prevalence(min_prevalence)
    n_after = abundance.n_mags
    logger.info(
        "Prevalence filtering (>= %.1f%%): %d -> %d MAGs (%d removed)",
        min_prevalence * 100, n_before, n_after, n_before - n_after,
    )

    # Alpha diversity
    alpha = div_mod.compute_alpha_diversity(abundance)
    _write_alpha_csv(alpha, out / "alpha_diversity.csv")
    plots.plot_alpha_diversity(alpha, metadata, grouping_var, out / "alpha_diversity.pdf")

    # Beta diversity
    bc = beta_mod.bray_curtis(abundance)
    jc = beta_mod.jaccard(abundance)
    _write_distance_matrix(bc, out / "beta_diversity_bray_curtis.csv")
    _write_distance_matrix(jc, out / "beta_diversity_jaccard.csv")

    # Ordination
    pcoa_result = ord_mod.pcoa(bc)
    _write_ordination(pcoa_result, out / "ordination_coordinates.csv")
    plots.plot_ordination(pcoa_result, metadata, grouping_var, out / "ordination_pcoa.pdf")

    try:
        nmds_result = ord_mod.nmds(bc)
        plots.plot_ordination(nmds_result, metadata, grouping_var, out / "ordination_nmds.pdf")
    except Exception:
        nmds_result = None

    # Statistical tests
    perm = stats_mod.permanova(bc, metadata, grouping_var, n_permutations=n_permutations)
    anos = stats_mod.anosim(bc, metadata, grouping_var, n_permutations=n_permutations)
    _write_stat_result(perm, out / "permanova_results.txt")
    _write_stat_result(anos, out / "anosim_results.txt")

    # PERMDISP
    pd_result = stats_mod.permdisp(bc, metadata, grouping_var, n_permutations=n_permutations)
    _write_stat_result(pd_result, out / "permdisp_results.txt")

    # Pairwise PERMANOVA
    pw_results = stats_mod.pairwise_permanova(
        bc, metadata, grouping_var, n_permutations=n_permutations
    )
    _write_pairwise_permanova(pw_results, out / "pairwise_permanova.csv")

    # Differential abundance — all pairwise comparisons
    groups = metadata.get_groups(grouping_var)
    group_names = sorted(groups.keys())
    diff_results = []
    for i, g1 in enumerate(group_names):
        for g2 in group_names[i + 1 :]:
            try:
                diff = diff_mod.differential_abundance(
                    abundance, metadata, grouping_var, g1, g2
                )
                diff_results.append((g1, g2, diff))
                _write_differential(
                    diff, g1, g2,
                    out / f"differential_{g1}_vs_{g2}.csv",
                    taxonomy=taxonomy,
                )
                plots.plot_volcano(diff, g1, g2, output=out / f"volcano_{g1}_vs_{g2}.pdf")
            except ValueError:
                pass

    # Heatmap of significant MAGs (from first comparison if available)
    if diff_results:
        g1, g2, diff = diff_results[0]
        sig_mags = [
            m for m, q in zip(diff.mag_ids, diff.q_values) if q < 0.05
        ]
        if sig_mags:
            plots.plot_heatmap(
                abundance, metadata, grouping_var, sig_mags,
                taxonomy=taxonomy, output=out / "heatmap_significant.pdf",
            )
        else:
            plots.plot_heatmap(
                abundance, metadata, grouping_var,
                taxonomy=taxonomy, output=out / "heatmap_all.pdf",
            )

    # Indicator species
    ind = ind_mod.indicator_species(
        abundance, metadata, grouping_var, n_permutations=n_permutations
    )
    _write_indicator(ind, out / "indicator_species.csv", taxonomy=taxonomy)
    plots.plot_indicator_species(
        ind, taxonomy=taxonomy, output=out / "indicator_species.pdf",
    )

    # Core microbiome
    core = cm_mod.core_microbiome(abundance, metadata, grouping_var)
    _write_core_microbiome(core, out / "core_microbiome.csv", taxonomy=taxonomy)

    # Compartment specificity
    cs = cs_mod.compartment_specificity(abundance, metadata, grouping_var)
    _write_specificity(cs, out / "compartment_specificity.csv", taxonomy=taxonomy)

    # Rank-level (phylum) differential abundance
    rank_diff_results: list[tuple[str, str, diff_mod.DifferentialAbundanceResult]] = []
    if taxonomy:
        for i, g1 in enumerate(group_names):
            for g2 in group_names[i + 1 :]:
                try:
                    rdiff = diff_mod.rank_level_differential_abundance(
                        abundance, taxonomy, metadata, grouping_var, g1, g2,
                        rank="phylum",
                    )
                    rank_diff_results.append((g1, g2, rdiff))
                    _write_differential(
                        rdiff, g1, g2,
                        out / f"phylum_differential_{g1}_vs_{g2}.csv",
                    )
                except ValueError:
                    pass

    # Phylum composition
    phylum_composition: dict[str, dict[str, tuple[float, float]]] | None = None
    if taxonomy:
        phylum_composition = _write_phylum_composition(
            abundance, taxonomy, metadata, grouping_var,
            out / "phylum_composition.csv",
        )

    # Taxonomy plots
    if taxonomy:
        plots.plot_taxonomy_bars(
            abundance, taxonomy, metadata, grouping_var, output=out / "taxonomy_bars.pdf"
        )
        plots.plot_taxonomy_bars(
            abundance, taxonomy, metadata, grouping_var,
            average_by_group=True, output=out / "taxonomy_bars_grouped.pdf",
        )

    # Venn diagram
    plots.plot_venn(abundance, metadata, grouping_var, output=out / "venn_diagram.pdf")

    # Summary report
    _write_summary(
        alpha, perm, anos, pd_result, pw_results, diff_results, ind, core, cs,
        group_names, n_before, n_after, min_prevalence, out / "summary_report.txt",
    )

    # HTML report
    report_data = {
        "output_dir": str(out),
        "grouping_var": grouping_var,
        "n_mags_before": n_before,
        "n_mags_after": n_after,
        "min_prevalence": min_prevalence,
        "alpha": alpha,
        "metadata": metadata,
        "permanova": perm,
        "anosim": anos,
        "permdisp": pd_result,
        "pairwise_permanova": pw_results,
        "differential": diff_results,
        "indicator": ind,
        "core_microbiome": core,
        "specificity": cs,
        "group_names": group_names,
        "taxonomy": taxonomy,
        "rank_diff_results": rank_diff_results,
        "phylum_composition": phylum_composition,
    }
    generate_html_report(str(out), report_data)


def _write_alpha_csv(result, path: Path) -> None:
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["sample_id", "shannon", "simpson", "richness", "evenness"])
        for i, sid in enumerate(result.sample_ids):
            w.writerow([
                sid,
                f"{result.shannon[i]:.4f}",
                f"{result.simpson[i]:.4f}",
                f"{result.richness[i]:.0f}",
                f"{result.evenness[i]:.4f}",
            ])


def _write_distance_matrix(result, path: Path) -> None:
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([""] + result.sample_ids)
        for i, sid in enumerate(result.sample_ids):
            w.writerow([sid] + [f"{result.distance_matrix[i, j]:.6f}" for j in range(len(result.sample_ids))])


def _write_ordination(result, path: Path) -> None:
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        n_axes = result.coordinates.shape[1]
        header = ["sample_id"] + [f"Axis{k+1}" for k in range(n_axes)]
        w.writerow(header)
        for i, sid in enumerate(result.sample_ids):
            w.writerow([sid] + [f"{result.coordinates[i, k]:.6f}" for k in range(n_axes)])


def _write_stat_result(result, path: Path) -> None:
    with open(path, "w") as f:
        f.write(f"Test: {result.test_name}\n")
        f.write(f"Statistic: {result.statistic:.4f}\n")
        f.write(f"P-value: {result.p_value:.4f}\n")
        for k, v in result.metadata.items():
            f.write(f"{k}: {v}\n")


def _write_pairwise_permanova(results: list, path: Path) -> None:
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["group1", "group2", "F_statistic", "p_value", "q_value", "R2"])
        for r in results:
            w.writerow([
                r.metadata.get("group1", ""),
                r.metadata.get("group2", ""),
                f"{r.statistic:.4f}",
                f"{r.p_value:.4f}",
                f"{r.metadata.get('q_value', r.p_value):.4f}",
                f"{r.metadata.get('R2', 0):.4f}",
            ])


def _write_differential(
    result, g1: str, g2: str, path: Path,
    taxonomy: TaxonomyTable | None = None,
) -> None:
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        if taxonomy:
            w.writerow(["MAG_ID", "phylum", "class", "genus",
                         "CLR_diff", "p_value", "q_value", "cohens_d"])
            tax_cols = taxonomy.lookup_columns(result.mag_ids)
        else:
            w.writerow(["MAG_ID", "CLR_diff", "p_value", "q_value", "cohens_d"])
            tax_cols = None
        for i, mag_id in enumerate(result.mag_ids):
            row = [mag_id]
            if tax_cols:
                row.extend([tax_cols[i]["phylum"], tax_cols[i]["class"], tax_cols[i]["genus"]])
            row.extend([
                f"{result.log_fold_changes[i]:.4f}",
                f"{result.p_values[i]:.6f}",
                f"{result.q_values[i]:.6f}",
                f"{result.effect_sizes[i]:.4f}",
            ])
            w.writerow(row)


def _write_indicator(
    result, path: Path,
    taxonomy: TaxonomyTable | None = None,
) -> None:
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        if taxonomy:
            w.writerow(["MAG_ID", "phylum", "class", "genus",
                         "best_group", "indval", "specificity", "fidelity", "p_value"])
            tax_cols = taxonomy.lookup_columns(result.mag_ids)
        else:
            w.writerow(["MAG_ID", "best_group", "indval", "specificity", "fidelity", "p_value"])
            tax_cols = None
        for i, mag_id in enumerate(result.mag_ids):
            row = [mag_id]
            if tax_cols:
                row.extend([tax_cols[i]["phylum"], tax_cols[i]["class"], tax_cols[i]["genus"]])
            row.extend([
                result.best_group[i],
                f"{result.indval_scores[i]:.4f}",
                f"{result.specificity[i]:.4f}",
                f"{result.fidelity[i]:.4f}",
                f"{result.p_values[i]:.4f}",
            ])
            w.writerow(row)


def _write_core_microbiome(
    result, path: Path,
    taxonomy: TaxonomyTable | None = None,
) -> None:
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        if taxonomy:
            w.writerow(["threshold", "group", "n_core_mags", "core_mag_ids", "core_mag_phyla"])
        else:
            w.writerow(["threshold", "group", "n_core_mags", "core_mag_ids"])
        for t in result.thresholds:
            for g, mags in sorted(result.core_mags_per_threshold[t].items()):
                row = [f"{t:.2f}", g, len(mags), ";".join(mags)]
                if taxonomy:
                    phyla = []
                    for m in mags:
                        rec = taxonomy.get(m)
                        phyla.append(rec.phylum if rec and rec.phylum else "Unclassified")
                    row.append(";".join(phyla))
                w.writerow(row)
            # Shared across all groups
            shared = result.shared_across_groups[t]
            row = [f"{t:.2f}", "_ALL_GROUPS_", len(shared), ";".join(shared)]
            if taxonomy:
                phyla = []
                for m in shared:
                    rec = taxonomy.get(m)
                    phyla.append(rec.phylum if rec and rec.phylum else "Unclassified")
                row.append(";".join(phyla))
            w.writerow(row)


def _write_specificity(
    result, path: Path,
    taxonomy: TaxonomyTable | None = None,
) -> None:
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        if taxonomy:
            w.writerow(["MAG_ID", "phylum", "class", "genus",
                         "specificity_score", "dominant_compartment"])
            tax_cols = taxonomy.lookup_columns(result.mag_ids)
        else:
            w.writerow(["MAG_ID", "specificity_score", "dominant_compartment"])
            tax_cols = None
        for i, mag_id in enumerate(result.mag_ids):
            row = [mag_id]
            if tax_cols:
                row.extend([tax_cols[i]["phylum"], tax_cols[i]["class"], tax_cols[i]["genus"]])
            row.extend([
                f"{result.scores[i]:.4f}",
                result.dominant_compartment[i],
            ])
            w.writerow(row)


def _write_phylum_composition(
    abundance: AbundanceTable,
    taxonomy: TaxonomyTable,
    metadata: SampleMetadata,
    grouping_var: str,
    path: Path,
) -> dict[str, dict[str, tuple[float, float]]]:
    """Write mean +/- SD relative abundance per phylum per group.

    Returns dict mapping phylum -> group -> (mean, sd).
    """
    rel = abundance.normalize()
    agg = taxonomy.aggregate_at_rank(rel, "phylum")
    groups = metadata.get_groups(grouping_var)
    group_names = sorted(groups.keys())
    sid_to_idx = {s: i for i, s in enumerate(abundance.sample_ids)}

    phylum_names = sorted(agg.keys())
    composition: dict[str, dict[str, tuple[float, float]]] = {}

    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        header = ["phylum"]
        for g in group_names:
            header.extend([f"{g}_mean", f"{g}_sd"])
        w.writerow(header)

        for phylum in phylum_names:
            composition[phylum] = {}
            row: list[str] = [phylum]
            for g in group_names:
                indices = [sid_to_idx[s] for s in groups[g] if s in sid_to_idx]
                if indices:
                    vals = agg[phylum][indices]
                    m, sd = float(vals.mean()), float(vals.std())
                else:
                    m, sd = 0.0, 0.0
                composition[phylum][g] = (m, sd)
                row.extend([f"{m:.6f}", f"{sd:.6f}"])
            w.writerow(row)

    return composition


def _taxonomy_breakdown(
    mag_ids: list[str],
    taxonomy: TaxonomyTable,
    rank: str = "phylum",
) -> dict[str, int]:
    """Count MAGs by taxon at *rank*."""
    counts: dict[str, int] = {}
    for mag_id in mag_ids:
        rec = taxonomy.get(mag_id)
        taxon = rec.rank(rank) if rec else ""
        if not taxon:
            taxon = "Unclassified"
        counts[taxon] = counts.get(taxon, 0) + 1
    return counts


def _write_summary(
    alpha, perm, anos, pd_result, pw_results, diff_results, ind, core, cs,
    group_names, n_before, n_after, min_prevalence, path: Path,
) -> None:
    with open(path, "w") as f:
        f.write("=" * 60 + "\n")
        f.write("MAG Community Profiling Summary Report\n")
        f.write("=" * 60 + "\n\n")

        f.write("PREVALENCE FILTERING\n")
        f.write("-" * 40 + "\n")
        f.write(f"Minimum prevalence: {min_prevalence*100:.0f}%\n")
        f.write(f"MAGs before filtering: {n_before}\n")
        f.write(f"MAGs after filtering: {n_after}\n")
        f.write(f"MAGs removed: {n_before - n_after}\n\n")

        f.write("ALPHA DIVERSITY\n")
        f.write("-" * 40 + "\n")
        f.write(f"Samples analyzed: {len(alpha.sample_ids)}\n")
        f.write(f"Shannon (mean +/- SD): {alpha.shannon.mean():.3f} +/- {alpha.shannon.std():.3f}\n")
        f.write(f"Simpson (mean +/- SD): {alpha.simpson.mean():.3f} +/- {alpha.simpson.std():.3f}\n")
        f.write(f"Richness (mean +/- SD): {alpha.richness.mean():.1f} +/- {alpha.richness.std():.1f}\n\n")

        f.write("BETA DIVERSITY & STATISTICAL TESTS\n")
        f.write("-" * 40 + "\n")
        f.write(f"PERMANOVA: F={perm.statistic:.3f}, p={perm.p_value:.4f}, R2={perm.metadata['R2']:.3f}\n")
        f.write(f"ANOSIM: R={anos.statistic:.3f}, p={anos.p_value:.4f}\n")
        f.write(f"PERMDISP: F={pd_result.statistic:.3f}, p={pd_result.p_value:.4f}\n\n")

        if pw_results:
            f.write("PAIRWISE PERMANOVA\n")
            f.write("-" * 40 + "\n")
            for r in pw_results:
                g1 = r.metadata.get("group1", "")
                g2 = r.metadata.get("group2", "")
                q = r.metadata.get("q_value", r.p_value)
                f.write(f"{g1} vs {g2}: F={r.statistic:.3f}, p={r.p_value:.4f}, q={q:.4f}\n")
            f.write("\n")

        if diff_results:
            f.write("DIFFERENTIAL ABUNDANCE\n")
            f.write("-" * 40 + "\n")
            for g1, g2, diff in diff_results:
                n_sig = int(np.sum(diff.q_values < 0.05))
                f.write(f"{g1} vs {g2}: {n_sig} significant MAGs (FDR < 0.05)\n")
            f.write("\n")

        f.write("INDICATOR SPECIES\n")
        f.write("-" * 40 + "\n")
        for g in group_names:
            mask = [bg == g for bg in ind.best_group]
            sig_mask = [m and p < 0.05 for m, p in zip(mask, ind.p_values)]
            n_sig = sum(sig_mask)
            f.write(f"{g}: {n_sig} significant indicators\n")
        f.write("\n")

        f.write("CORE MICROBIOME\n")
        f.write("-" * 40 + "\n")
        for t in core.thresholds:
            n_shared = len(core.shared_across_groups[t])
            f.write(f"Threshold {t*100:.0f}%: {n_shared} MAGs core across all groups\n")
        f.write("\n")

        f.write("COMPARTMENT SPECIFICITY\n")
        f.write("-" * 40 + "\n")
        specialists = np.sum(cs.scores > 0.5)
        generalists = np.sum(cs.scores <= 0.5)
        f.write(f"Specialists (score > 0.5): {specialists}\n")
        f.write(f"Generalists (score <= 0.5): {generalists}\n")
