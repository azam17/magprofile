"""Publication-quality visualizations for MAG community profiling."""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import numpy as np
import seaborn as sns
from scipy import stats as sp_stats

from .compartment_specificity import CompartmentSpecificityResult
from .differential import DifferentialAbundanceResult
from .diversity import AlphaDiversityResult
from .indicator import IndicatorSpeciesResult
from .io import AbundanceTable, SampleMetadata, TaxonomyTable
from .ordination import OrdinationResult

# Consistent style
PALETTE = sns.color_palette("Set2")
DPI = 300


def _setup_style() -> None:
    sns.set_theme(style="whitegrid", palette="Set2")
    plt.rcParams.update({"figure.dpi": DPI, "savefig.dpi": DPI, "font.size": 10})


def _significance_stars(p: float) -> str:
    """Convert p-value to significance stars."""
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def _bh_fdr(p_values: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR correction."""
    n = len(p_values)
    if n == 0:
        return np.array([])
    sorted_idx = np.argsort(p_values)
    sorted_p = p_values[sorted_idx]
    q = np.zeros(n)
    for i in range(n):
        rank = i + 1
        q[i] = sorted_p[i] * n / rank
    for i in range(n - 2, -1, -1):
        q[i] = min(q[i], q[i + 1])
    q = np.clip(q, 0, 1)
    result = np.zeros(n)
    result[sorted_idx] = q
    return result


def plot_alpha_diversity(
    result: AlphaDiversityResult,
    metadata: SampleMetadata,
    grouping_var: str,
    output: str | Path,
) -> None:
    """2x2 boxplot grid: Shannon, Simpson, Richness, Evenness by compartment.

    Includes Kruskal-Wallis test and pairwise Mann-Whitney U with FDR
    correction, displayed as significance brackets.
    """
    _setup_style()

    groups = metadata.get_groups(grouping_var)
    sample_to_group = {}
    for g, sids in groups.items():
        for s in sids:
            sample_to_group[s] = g

    labels = [sample_to_group.get(s, "Unknown") for s in result.sample_ids]
    group_names = sorted(set(labels))
    metrics = {
        "Shannon": result.shannon,
        "Simpson": result.simpson,
        "Richness": result.richness,
        "Evenness": result.evenness,
    }

    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    for ax, (name, values) in zip(axes.flat, metrics.items()):
        data = {"Group": labels, name: values}
        sns.boxplot(x="Group", y=name, data=data, ax=ax, palette="Set2")
        ax.set_title(name)
        ax.set_xlabel("")

        # Statistical annotations
        if len(group_names) >= 2:
            group_values = {}
            for g in group_names:
                mask = [l == g for l in labels]
                group_values[g] = values[np.array(mask)]

            # Kruskal-Wallis overall test
            kw_groups = [group_values[g] for g in group_names if len(group_values[g]) > 0]
            if len(kw_groups) >= 2:
                try:
                    _, kw_p = sp_stats.kruskal(*kw_groups)
                except ValueError:
                    kw_p = 1.0

                if kw_p < 0.05 and len(group_names) >= 2:
                    # Pairwise Mann-Whitney U with FDR
                    pairs = []
                    raw_p = []
                    for i, g1 in enumerate(group_names):
                        for g2 in group_names[i + 1:]:
                            pairs.append((g1, g2))
                            try:
                                _, p = sp_stats.mannwhitneyu(
                                    group_values[g1], group_values[g2],
                                    alternative="two-sided",
                                )
                            except ValueError:
                                p = 1.0
                            raw_p.append(p)

                    q_values = _bh_fdr(np.array(raw_p))

                    # Draw significance brackets for significant pairs
                    y_max = max(values)
                    y_range = max(values) - min(values)
                    if y_range == 0:
                        y_range = 1.0
                    bracket_offset = y_range * 0.08
                    y_pos = y_max + bracket_offset

                    group_positions = {g: i for i, g in enumerate(group_names)}
                    for (g1, g2), q in zip(pairs, q_values):
                        stars = _significance_stars(q)
                        if stars == "ns":
                            continue
                        x1 = group_positions[g1]
                        x2 = group_positions[g2]
                        ax.plot(
                            [x1, x1, x2, x2],
                            [y_pos, y_pos + bracket_offset * 0.3,
                             y_pos + bracket_offset * 0.3, y_pos],
                            color="black", linewidth=0.8,
                        )
                        ax.text(
                            (x1 + x2) / 2, y_pos + bracket_offset * 0.35,
                            stars, ha="center", va="bottom", fontsize=9,
                        )
                        y_pos += bracket_offset * 1.2

    fig.suptitle("Alpha Diversity by Compartment", fontsize=14, y=1.02)
    fig.tight_layout()
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def plot_ordination(
    result: OrdinationResult,
    metadata: SampleMetadata,
    grouping_var: str,
    output: str | Path,
) -> None:
    """Scatter plot of ordination coordinates with 95% confidence ellipses."""
    _setup_style()

    groups = metadata.get_groups(grouping_var)
    sample_to_group = {}
    for g, sids in groups.items():
        for s in sids:
            sample_to_group[s] = g

    fig, ax = plt.subplots(figsize=(8, 6))
    group_names = sorted(set(sample_to_group.values()))
    colors = {g: PALETTE[i % len(PALETTE)] for i, g in enumerate(group_names)}

    for g in group_names:
        mask = [sample_to_group.get(s) == g for s in result.sample_ids]
        coords = result.coordinates[mask]
        ax.scatter(coords[:, 0], coords[:, 1], label=g, color=colors[g], s=60, alpha=0.8)

        # 95% confidence ellipse (need >= 3 samples)
        if coords.shape[0] >= 3:
            mean = coords[:, :2].mean(axis=0)
            cov = np.cov(coords[:, 0], coords[:, 1])
            # Eigendecomposition for ellipse parameters
            eigenvalues, eigenvectors = np.linalg.eigh(cov)
            # Sort descending
            order = eigenvalues.argsort()[::-1]
            eigenvalues = eigenvalues[order]
            eigenvectors = eigenvectors[:, order]
            # Angle of major axis
            angle = np.degrees(np.arctan2(eigenvectors[1, 0], eigenvectors[0, 0]))
            # 95% CI: chi-squared with 2 df at 0.05 = 5.991
            chi2_val = 5.991
            width = 2 * np.sqrt(chi2_val * max(eigenvalues[0], 0))
            height = 2 * np.sqrt(chi2_val * max(eigenvalues[1], 0))
            ellipse = Ellipse(
                xy=mean, width=width, height=height, angle=angle,
                facecolor=colors[g], alpha=0.15, edgecolor=colors[g],
                linewidth=1.5, linestyle="--",
            )
            ax.add_patch(ellipse)

    xlabel = "Axis 1"
    ylabel = "Axis 2"
    if result.explained_variance is not None:
        xlabel += f" ({result.explained_variance[0]*100:.1f}%)"
        ylabel += f" ({result.explained_variance[1]*100:.1f}%)"

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(f"{result.method} Ordination")
    if result.stress is not None:
        ax.text(
            0.02, 0.02, f"Stress: {result.stress:.4f}",
            transform=ax.transAxes, fontsize=9,
        )
    ax.legend()
    fig.tight_layout()
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def plot_taxonomy_bars(
    table: AbundanceTable,
    taxonomy: TaxonomyTable,
    metadata: SampleMetadata,
    grouping_var: str,
    rank: str = "phylum",
    top_n: int = 10,
    average_by_group: bool = False,
    output: str | Path = "taxonomy_bars.pdf",
) -> None:
    """Stacked bar chart of relative abundance at given rank, top-N + Other.

    When average_by_group=True, averages relative abundances within each
    group before plotting (produces cleaner bars for many-sample datasets).
    """
    _setup_style()

    rel = table.normalize()
    agg = taxonomy.aggregate_at_rank(rel, rank)

    # Sort taxa by total abundance
    taxon_totals = {t: v.sum() for t, v in agg.items()}
    sorted_taxa = sorted(taxon_totals, key=taxon_totals.get, reverse=True)
    top_taxa = sorted_taxa[:top_n]

    if average_by_group:
        # Average within groups
        groups = metadata.get_groups(grouping_var)
        group_names = sorted(groups.keys())
        sid_to_idx = {s: i for i, s in enumerate(table.sample_ids)}

        sample_data = np.zeros((len(top_taxa) + 1, len(group_names)))
        for i, t in enumerate(top_taxa):
            for gi, g in enumerate(group_names):
                indices = [sid_to_idx[s] for s in groups[g] if s in sid_to_idx]
                if indices:
                    sample_data[i, gi] = agg[t][indices].mean()
        for t in sorted_taxa[top_n:]:
            for gi, g in enumerate(group_names):
                indices = [sid_to_idx[s] for s in groups[g] if s in sid_to_idx]
                if indices:
                    sample_data[-1, gi] += agg[t][indices].mean()

        x_labels = group_names
    else:
        # Per-sample bars
        sample_data = np.zeros((len(top_taxa) + 1, table.n_samples))
        for i, t in enumerate(top_taxa):
            sample_data[i] = agg[t]
        for t in sorted_taxa[top_n:]:
            sample_data[-1] += agg[t]

        # Sort samples by group
        groups = metadata.get_groups(grouping_var)
        ordered_samples: list[str] = []
        for g in sorted(groups.keys()):
            ordered_samples.extend(sorted(groups[g]))
        order_idx = [table.sample_ids.index(s) for s in ordered_samples if s in table.sample_ids]
        sample_data = sample_data[:, order_idx]
        x_labels = [table.sample_ids[i] for i in order_idx]

    labels = list(top_taxa) + ["Other"]

    fig, ax = plt.subplots(figsize=(max(8, len(x_labels) * 0.5), 6))
    x = np.arange(len(x_labels))
    bottom = np.zeros(len(x_labels))
    bar_colors = sns.color_palette("tab20", len(labels))

    for i, (label, color) in enumerate(zip(labels, bar_colors)):
        ax.bar(x, sample_data[i], bottom=bottom, label=label, color=color, width=0.8)
        bottom += sample_data[i]

    ax.set_xticks(x)
    rotation = 0 if average_by_group else 90
    fontsize = 10 if average_by_group else 7
    ax.set_xticklabels(x_labels, rotation=rotation, fontsize=fontsize)
    ax.set_ylabel("Relative Abundance")
    title_suffix = " (Group Averages)" if average_by_group else ""
    ax.set_title(f"Taxonomy at {rank.capitalize()} Level{title_suffix}")
    ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=8)
    fig.tight_layout()
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def plot_volcano(
    result: DifferentialAbundanceResult,
    group1: str,
    group2: str,
    fdr_threshold: float = 0.05,
    lfc_threshold: float = 1.0,
    output: str | Path = "volcano.pdf",
) -> None:
    """Volcano plot: CLR difference vs -log10(q-value)."""
    _setup_style()

    lfc = result.log_fold_changes
    neg_log_q = -np.log10(np.clip(result.q_values, 1e-300, 1.0))

    sig = (result.q_values < fdr_threshold) & (np.abs(lfc) > lfc_threshold)
    up = sig & (lfc > 0)
    down = sig & (lfc < 0)
    ns = ~sig

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(lfc[ns], neg_log_q[ns], c="grey", alpha=0.4, s=20, label="NS")
    ax.scatter(lfc[up], neg_log_q[up], c="firebrick", alpha=0.7, s=30, label=f"Up in {group1}")
    ax.scatter(lfc[down], neg_log_q[down], c="steelblue", alpha=0.7, s=30, label=f"Up in {group2}")

    ax.axhline(-np.log10(fdr_threshold), ls="--", c="grey", lw=0.8)
    ax.axvline(lfc_threshold, ls="--", c="grey", lw=0.8)
    ax.axvline(-lfc_threshold, ls="--", c="grey", lw=0.8)

    ax.set_xlabel("CLR Difference")
    ax.set_ylabel("-log10(q-value)")
    ax.set_title(f"Differential Abundance: {group1} vs {group2}")
    ax.legend()
    fig.tight_layout()
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def plot_heatmap(
    table: AbundanceTable,
    metadata: SampleMetadata,
    grouping_var: str,
    mag_subset: list[str] | None = None,
    taxonomy: TaxonomyTable | None = None,
    output: str | Path = "heatmap.pdf",
) -> None:
    """Heatmap of CLR-transformed abundances, samples sorted by compartment.

    When taxonomy is provided, appends phylum/class to MAG IDs on y-axis.
    """
    from .differential import clr_transform

    _setup_style()

    clr = clr_transform(table.abundances)

    # Subset MAGs if specified
    if mag_subset:
        idx = [table.mag_ids.index(m) for m in mag_subset if m in table.mag_ids]
        clr = clr[idx]
        row_mag_ids = [table.mag_ids[i] for i in idx]
    else:
        row_mag_ids = list(table.mag_ids)

    # Build row labels with taxonomy annotation
    if taxonomy:
        row_labels = []
        for m in row_mag_ids:
            rec = taxonomy.get(m)
            if rec and (rec.phylum or rec.class_):
                tax_str = rec.phylum or ""
                if rec.class_:
                    tax_str += f"/{rec.class_}"
                row_labels.append(f"{m}  ({tax_str})")
            else:
                row_labels.append(m)
    else:
        row_labels = row_mag_ids

    # Sort samples by group
    groups = metadata.get_groups(grouping_var)
    ordered_samples: list[str] = []
    for g in sorted(groups.keys()):
        ordered_samples.extend(sorted(groups[g]))
    order_idx = [table.sample_ids.index(s) for s in ordered_samples if s in table.sample_ids]
    clr = clr[:, order_idx]
    col_labels = [table.sample_ids[i] for i in order_idx]

    fig, ax = plt.subplots(figsize=(max(8, len(col_labels) * 0.4), max(6, len(row_labels) * 0.3)))
    sns.heatmap(
        clr,
        xticklabels=col_labels,
        yticklabels=row_labels,
        cmap="RdBu_r",
        center=0,
        ax=ax,
    )
    # Use smaller font for taxonomy-annotated labels
    if taxonomy:
        ax.set_yticklabels(ax.get_yticklabels(), fontsize=7)
    ax.set_title("CLR-Transformed Abundance")
    ax.set_xlabel("Sample")
    ax.set_ylabel("MAG")
    fig.tight_layout()
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def plot_venn(
    table: AbundanceTable,
    metadata: SampleMetadata,
    grouping_var: str,
    output: str | Path = "venn.pdf",
) -> None:
    """Venn diagram of shared/unique MAGs across compartments."""
    from matplotlib_venn import venn2, venn3

    _setup_style()

    groups = metadata.get_groups(grouping_var)
    sid_to_idx = {s: i for i, s in enumerate(table.sample_ids)}

    compartment_mags: dict[str, set[str]] = {}
    for g, sids in groups.items():
        indices = [sid_to_idx[s] for s in sids if s in sid_to_idx]
        if indices:
            present = table.abundances[:, indices].sum(axis=1) > 0
            compartment_mags[g] = {m for m, p in zip(table.mag_ids, present) if p}
        else:
            compartment_mags[g] = set()

    comp_names = sorted(compartment_mags.keys())
    sets = [compartment_mags[c] for c in comp_names]

    fig, ax = plt.subplots(figsize=(8, 6))
    if len(sets) == 2:
        venn2(sets, set_labels=comp_names, ax=ax)
    elif len(sets) == 3:
        venn3(sets, set_labels=comp_names, ax=ax)
    else:
        ax.text(0.5, 0.5, f"Venn requires 2-3 groups, got {len(sets)}", ha="center", va="center")

    ax.set_title("Shared and Unique MAGs by Compartment")
    fig.tight_layout()
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def plot_indicator_species(
    result: IndicatorSpeciesResult,
    top_n: int = 10,
    taxonomy: TaxonomyTable | None = None,
    output: str | Path = "indicator_species.pdf",
) -> None:
    """Barplot of top indicator MAGs per compartment with IndVal scores."""
    _setup_style()

    # Get top indicators per group
    group_names = sorted(set(result.best_group))
    plot_data: list[tuple[str, str, float]] = []  # (mag_id, group, indval)

    for g in group_names:
        mask = [bg == g for bg in result.best_group]
        indices = [i for i, m in enumerate(mask) if m]
        # Sort by indval within group
        sorted_idx = sorted(indices, key=lambda i: result.indval_scores[i], reverse=True)
        for i in sorted_idx[:top_n]:
            plot_data.append((result.mag_ids[i], g, float(result.indval_scores[i])))

    if not plot_data:
        return

    # Build labels with optional taxonomy annotation
    mag_labels = []
    for d in plot_data:
        label = d[0]
        if taxonomy:
            rec = taxonomy.get(d[0])
            if rec and rec.phylum:
                label = f"{d[0]} ({rec.phylum})"
        mag_labels.append(label)
    grp_labels = [d[1] for d in plot_data]
    scores = [d[2] for d in plot_data]

    fig, ax = plt.subplots(figsize=(10, max(4, len(plot_data) * 0.3)))
    colors_map = {g: PALETTE[i % len(PALETTE)] for i, g in enumerate(group_names)}
    bar_colors = [colors_map[g] for g in grp_labels]

    y = np.arange(len(mag_labels))
    ax.barh(y, scores, color=bar_colors)
    ax.set_yticks(y)
    ax.set_yticklabels(mag_labels, fontsize=8)
    ax.set_xlabel("IndVal Score")
    ax.set_title("Indicator Species by Compartment")
    ax.invert_yaxis()

    # Legend
    from matplotlib.patches import Patch
    handles = [Patch(color=colors_map[g], label=g) for g in group_names]
    ax.legend(handles=handles, loc="lower right")

    fig.tight_layout()
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def plot_pgpr_matrix(
    pgpr_table,
    taxonomy=None,
    top_n: int = 30,
    output: str | Path = "pgpr_traits.pdf",
) -> None:
    """Binary heatmap of PGPR traits per MAG."""
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
        new_labels = []
        for m in mag_labels:
            rec = taxonomy.get(m)
            if rec and rec.phylum:
                new_labels.append(f"{m} ({rec.phylum})")
            else:
                new_labels.append(m)
        mag_labels = new_labels

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
    cazy_table,
    metadata,
    grouping_var: str,
    abundance,
    output: str | Path = "cazyme_bars.pdf",
) -> None:
    """Stacked bar chart of CAZyme classes per group."""
    from .func_profile import pathway_abundance as _pathway_abundance

    _setup_style()
    groups = metadata.get_groups(grouping_var)
    group_names = sorted(groups.keys())

    # Get CAZyme abundance per sample
    cazy_abund = _pathway_abundance(cazy_table, abundance)

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


def plot_network_graph(
    network,
    taxonomy=None,
    output: str | Path = "network_graph.pdf",
) -> None:
    """Force-directed network graph with nodes colored by phylum."""
    import networkx as nx

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
    topology,
    output: str | Path = "degree_distribution.pdf",
) -> None:
    """Histogram of node degree distribution."""
    _setup_style()
    fig, ax = plt.subplots(figsize=(7, 5))
    degrees = topology.degree[topology.degree > 0]
    if len(degrees) == 0:
        ax.text(0.5, 0.5, "No connected nodes", ha="center", va="center", transform=ax.transAxes)
    else:
        n_bins = max(1, min(30, int(np.ceil(degrees.max() * 10))))
        ax.hist(degrees, bins=n_bins, color=PALETTE[0], edgecolor="white")
        ax.set_xlabel("Degree")
        ax.set_ylabel("Number of MAGs")
        ax.set_title("Degree Distribution")
        ax.axvline(np.mean(degrees), color="red", linestyle="--", label=f"Mean={np.mean(degrees):.1f}")
        ax.legend()
    fig.savefig(str(output), bbox_inches="tight")
    plt.close(fig)
