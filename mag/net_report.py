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
from .net_correlation import (
    NetworkResult,
    ThresholdSensitivityResult,
    build_network,
    proportionality,
    threshold_sensitivity,
)
from .net_topology import (
    HubBridgeResult,
    KeystoneTaxaResult,
    ModuleCompositionResult,
    NullModelResult,
    NetworkTopology,
    compute_topology,
    differential_network,
    hub_bridge_classification,
    identify_keystones,
    module_composition,
    network_null_model,
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
    threshold_mode: str = "global",
    phi_threshold: float | None = None,
    group_min_prevalence: float | None = None,
) -> None:
    """Run network analysis pipeline and write results."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    threshold_mode = str(threshold_mode).lower()
    if threshold_mode not in {"global", "group", "fixed"}:
        raise ValueError("threshold_mode must be one of: global, group, fixed")

    # Global network
    corr = proportionality(abundance, min_prevalence=min_prevalence)
    global_net = build_network(corr, abundance.mag_ids, threshold_percentile)
    if phi_threshold is not None:
        global_threshold = float(phi_threshold)
    else:
        global_threshold = float(global_net.threshold)
    _write_edges_csv(global_net, out / "network_edges.csv", taxonomy)

    global_topo = compute_topology(global_net)
    _write_topology_csv(global_topo, out / "network_topology.csv", taxonomy)

    keystones = identify_keystones(global_topo, abundance)
    _write_keystones_csv(keystones, out / "keystone_taxa.csv", taxonomy)

    # Null model test for modularity
    null_model = network_null_model(global_net)
    _write_null_model_csv(null_model, out / "null_model.csv")

    # Threshold sensitivity analysis
    sensitivity = threshold_sensitivity(corr, abundance.mag_ids)
    _write_sensitivity_csv(sensitivity, out / "threshold_sensitivity.csv")

    # Module taxonomic composition
    if taxonomy:
        mod_comp = module_composition(global_topo, taxonomy)
        _write_module_composition_csv(mod_comp, out / "module_composition.csv")

    # Hub vs bridge classification (Guimera-Amaral z-P)
    hub_bridge = hub_bridge_classification(global_topo, global_net)
    _write_hub_bridge_csv(hub_bridge, out / "hub_bridge_classification.csv", taxonomy)

    # Keystone × PGPR cross-reference (if pgpr_table provided via kwargs)
    # This is wired from generate_combined_report if both func + net are run

    # Per-group networks
    groups = metadata.get_groups(grouping_var)
    group_names = sorted(groups.keys())
    group_nets: dict[str, NetworkResult] = {}

    eff_group_min_prev = min_prevalence if group_min_prevalence is None else float(group_min_prevalence)

    for g in group_names:
        try:
            sub = abundance.subset_samples(groups[g])
            g_corr = proportionality(sub, min_prevalence=eff_group_min_prev)
            if threshold_mode == "group":
                g_net = build_network(g_corr, sub.mag_ids, threshold_percentile)
            elif threshold_mode == "fixed":
                if phi_threshold is None:
                    raise ValueError("phi_threshold is required when threshold_mode='fixed'")
                g_net = _build_network_with_absolute_threshold(g_corr, sub.mag_ids, float(phi_threshold))
            else:  # global
                g_net = _build_network_with_absolute_threshold(g_corr, sub.mag_ids, global_threshold)
            group_nets[g] = g_net
            _write_edges_csv(g_net, out / f"network_edges_{g}.csv", taxonomy)
        except Exception as e:
            logger.warning("Network for group %s failed: %s", g, e)

    _write_network_summary_csv(
        out / "network_summary.csv",
        global_net,
        group_nets,
        threshold_mode=threshold_mode,
        threshold_percentile=threshold_percentile,
        min_prevalence=min_prevalence,
        group_min_prevalence=eff_group_min_prev,
    )

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


def _build_network_with_absolute_threshold(
    corr_matrix: np.ndarray,
    mag_ids: list[str],
    threshold: float,
) -> NetworkResult:
    """Build network using a fixed absolute phi threshold."""
    n = len(mag_ids)
    upper_i, upper_j = np.triu_indices(n, k=1)
    upper_vals = corr_matrix[upper_i, upper_j]
    valid_mask = ~np.isnan(upper_vals)
    edge_mask = valid_mask & (upper_vals <= threshold)
    edge_indices = np.where(edge_mask)[0]

    adjacency = np.zeros((n, n), dtype=np.float64)
    ei = upper_i[edge_indices]
    ej = upper_j[edge_indices]
    adjacency[ei, ej] = 1.0
    adjacency[ej, ei] = 1.0

    edges: list[tuple[str, str, float]] = [
        (mag_ids[upper_i[k]], mag_ids[upper_j[k]], float(upper_vals[k]))
        for k in edge_indices
    ]

    return NetworkResult(
        mag_ids=list(mag_ids),
        adjacency=adjacency,
        edges=edges,
        threshold=float(threshold),
    )


def _write_network_summary_csv(
    path: Path,
    global_net: NetworkResult,
    group_nets: dict[str, NetworkResult],
    *,
    threshold_mode: str,
    threshold_percentile: float,
    min_prevalence: float,
    group_min_prevalence: float,
) -> None:
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "scope",
            "group",
            "n_mags",
            "n_possible_edges",
            "n_edges",
            "density",
            "threshold_mode",
            "threshold_percentile",
            "threshold_phi",
            "min_prevalence",
        ])

        def _row(scope: str, group: str, net: NetworkResult, prev: float) -> list:
            n = len(net.mag_ids)
            n_possible = n * (n - 1) // 2
            density = (len(net.edges) / n_possible) if n_possible else 0.0
            return [
                scope,
                group,
                n,
                n_possible,
                len(net.edges),
                f"{density:.6f}",
                threshold_mode,
                f"{threshold_percentile:.2f}",
                f"{net.threshold:.6f}",
                f"{prev:.3f}",
            ]

        w.writerow(_row("global", "all", global_net, min_prevalence))
        for g in sorted(group_nets.keys()):
            w.writerow(_row("group", g, group_nets[g], group_min_prevalence))


def _write_null_model_csv(result: NullModelResult, path: Path) -> None:
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["metric", "value"])
        w.writerow(["observed_modularity", f"{result.observed_modularity:.6f}"])
        w.writerow(["null_mean", f"{result.null_modularities.mean():.6f}"])
        w.writerow(["null_sd", f"{result.null_modularities.std():.6f}"])
        w.writerow(["z_score", f"{result.z_score:.4f}"])
        w.writerow(["p_value", f"{result.p_value:.6f}"])


def _write_sensitivity_csv(result: ThresholdSensitivityResult, path: Path) -> None:
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["percentile", "threshold", "n_edges", "modularity", "n_modules"])
        for i in range(len(result.percentiles)):
            w.writerow([
                f"{result.percentiles[i]:.1f}",
                f"{result.thresholds[i]:.6f}",
                int(result.n_edges[i]),
                f"{result.modularities[i]:.6f}",
                int(result.n_modules[i]),
            ])


def _write_edges_csv(net: NetworkResult, path: Path, taxonomy: TaxonomyTable | None) -> None:
    # Build phylum lookup once
    phylum_map: dict[str, str] = {}
    if taxonomy:
        for m in net.mag_ids:
            rec = taxonomy.get(m)
            phylum_map[m] = rec.phylum if rec else ""

    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        header = ["MAG_1", "MAG_2", "phi"]
        if taxonomy:
            header.extend(["phylum_1", "phylum_2"])
        w.writerow(header)
        for m1, m2, phi in sorted(net.edges, key=lambda e: e[2]):
            row = [m1, m2, f"{phi:.6f}"]
            if taxonomy:
                row.extend([phylum_map.get(m1, ""), phylum_map.get(m2, "")])
            w.writerow(row)


def _write_topology_csv(topo: NetworkTopology, path: Path, taxonomy: TaxonomyTable | None) -> None:
    # Build taxonomy lookup once
    tax_map: dict[str, tuple[str, str, str]] = {}
    if taxonomy:
        for m in topo.mag_ids:
            rec = taxonomy.get(m)
            tax_map[m] = (
                rec.phylum if rec else "",
                rec.rank("class") if rec else "",
                rec.genus if rec else "",
            )

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
                row.extend(tax_map[m])
            row.extend([
                int(topo.degree[i]),
                f"{topo.betweenness[i]:.6f}",
                f"{topo.closeness[i]:.6f}",
                f"{topo.hub_scores[i]:.6f}",
                int(topo.module_assignments[i]),
            ])
            w.writerow(row)


def _write_keystones_csv(ks: KeystoneTaxaResult, path: Path, taxonomy: TaxonomyTable | None) -> None:
    # Build taxonomy lookup once
    tax_map: dict[str, tuple[str, str, str]] = {}
    if taxonomy:
        for m in ks.mag_ids:
            rec = taxonomy.get(m)
            tax_map[m] = (
                rec.phylum if rec else "",
                rec.rank("class") if rec else "",
                rec.genus if rec else "",
            )

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
                row.extend(tax_map[m])
            row.extend([
                f"{ks.keystone_scores[i]:.6f}",
                "yes" if ks.is_keystone[i] else "no",
                f"{ks.metrics['betweenness'][i]:.6f}",
                f"{ks.metrics['mean_abundance'][i]:.2f}",
            ])
            w.writerow(row)


def _write_module_composition_csv(result: ModuleCompositionResult, path: Path) -> None:
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["module_id", "n_mags", "dominant_phylum", "dominant_fraction", "phylum_breakdown"])
        for i, mod_id in enumerate(result.module_ids):
            breakdown = "; ".join(
                f"{p}:{c}" for p, c in sorted(
                    result.phylum_counts[mod_id].items(), key=lambda x: -x[1]
                )
            )
            w.writerow([
                mod_id,
                result.n_mags[i],
                result.dominant_phylum[i],
                f"{result.dominant_fraction[i]:.3f}",
                breakdown,
            ])


def _write_hub_bridge_csv(
    result: HubBridgeResult, path: Path, taxonomy: TaxonomyTable | None,
) -> None:
    tax_map: dict[str, tuple[str, str]] = {}
    if taxonomy:
        for m in result.mag_ids:
            rec = taxonomy.get(m)
            tax_map[m] = (rec.phylum if rec else "", rec.genus if rec else "")

    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        header = ["MAG_ID"]
        if taxonomy:
            header.extend(["phylum", "genus"])
        header.extend(["z_score", "participation_coeff", "role"])
        w.writerow(header)
        for i, m in enumerate(result.mag_ids):
            row = [m]
            if taxonomy:
                row.extend(tax_map.get(m, ("", "")))
            row.extend([
                f"{result.within_module_degree_z[i]:.4f}",
                f"{result.participation_coefficient[i]:.4f}",
                result.roles[i],
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
        "<title>magnet -- Network Analysis Report</title>",
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
        "<h1>magnet -- Network Analysis Report</h1>",
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
        # Build taxonomy lookup once for keystones
        ks_tax: dict[str, tuple[str, str]] = {}
        if taxonomy:
            for m in keystones.mag_ids:
                rec = taxonomy.get(m)
                ks_tax[m] = (rec.phylum if rec else "", rec.genus if rec else "")

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
                phylum, genus = ks_tax[m]
                html_parts.append(f"<td>{phylum}</td>")
                html_parts.append(f"<td>{genus}</td>")
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
