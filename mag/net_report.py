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
