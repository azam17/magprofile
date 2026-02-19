# mag/func_report.py
"""Report generation for functional profiling."""

from __future__ import annotations

import base64
import csv
import logging
from itertools import combinations
from pathlib import Path

import numpy as np

from .differential import differential_abundance
from .func_io import DRAMAnnotation, FunctionalTable, build_all_functional_tables, build_functional_table
from .func_profile import (
    cazyme_summary,
    differential_pathway,
    functional_redundancy,
    pathway_abundance,
    pgpr_enrichment,
    redundancy_comparison,
)
from .io import AbundanceTable, SampleMetadata, TaxonomyTable

logger = logging.getLogger(__name__)


def generate_func_report(
    abundance: AbundanceTable,
    annotations: list[DRAMAnnotation],
    taxonomy: TaxonomyTable | None,
    metadata: SampleMetadata,
    grouping_var: str,
    output_dir: str | Path,
) -> None:
    """Run functional profiling pipeline and write results."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    # Build all tables in single pass over annotations
    tables = build_all_functional_tables(annotations)
    ko_table = tables["ko"]
    cazy_table = tables["cazy"]
    pgpr_table = tables["pgpr"]

    # Pathway abundance
    pa = pathway_abundance(ko_table, abundance)
    _write_abundance_csv(pa, out / "pathway_abundance.csv")

    # PGPR traits
    _write_functional_table_csv(pgpr_table, out / "pgpr_traits.csv")

    # CAZyme summary
    cazy_stats = cazyme_summary(cazy_table)
    _write_cazyme_csv(cazy_stats, out / "cazyme_summary.csv")

    # Functional redundancy
    redundancy = functional_redundancy(ko_table, abundance, metadata, grouping_var)
    _write_redundancy_csv(redundancy, out / "functional_redundancy.csv")

    # Redundancy comparison across groups
    red_cmp = redundancy_comparison(redundancy)
    _write_redundancy_comparison_csv(red_cmp, out / "redundancy_comparison.csv")

    # PGPR enrichment by group
    pgpr_enrich = pgpr_enrichment(pgpr_table, abundance, metadata, grouping_var)
    _write_pgpr_enrichment_csv(pgpr_enrich, out / "pgpr_enrichment.csv")

    # Group setup for pairwise comparisons
    groups = metadata.get_groups(grouping_var)
    group_names = sorted(groups.keys())

    # CAZyme differential (same as pathway differential but for CAZy table)
    cazy_diff_results = {}
    for g1, g2 in combinations(group_names, 2):
        try:
            cazy_diff = differential_pathway(
                cazy_table, abundance, metadata, grouping_var, g1, g2,
            )
            cazy_diff_results[(g1, g2)] = cazy_diff
            _write_diff_csv(cazy_diff, g1, g2, out / f"cazyme_differential_{g1}_vs_{g2}.csv")
        except Exception as e:
            logger.warning("CAZyme differential %s vs %s failed: %s", g1, g2, e)

    # Differential pathway abundance for all group pairs
    diff_results = {}

    for g1, g2 in combinations(group_names, 2):
        try:
            result = differential_pathway(
                ko_table, abundance, metadata, grouping_var, g1, g2,
            )
            diff_results[(g1, g2)] = result
            _write_diff_csv(result, g1, g2, out / f"pathway_differential_{g1}_vs_{g2}.csv")
        except Exception as e:
            logger.warning("Pathway differential %s vs %s failed: %s", g1, g2, e)

    # Plots
    from . import plots

    plots.plot_pgpr_matrix(pgpr_table, taxonomy, output=out / "pgpr_traits.pdf")
    plots.plot_cazyme_bars(cazy_table, metadata, grouping_var, abundance, output=out / "cazyme_bars.pdf")

    # HTML report
    _write_func_html(
        out,
        ko_table=ko_table,
        cazy_table=cazy_table,
        pgpr_table=pgpr_table,
        cazy_stats=cazy_stats,
        redundancy=redundancy,
        diff_results=diff_results,
        taxonomy=taxonomy,
    )


def _write_abundance_csv(table: AbundanceTable, path: Path) -> None:
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["function_id"] + table.sample_ids)
        for i, fid in enumerate(table.mag_ids):
            w.writerow([fid] + [f"{v:.4f}" for v in table.abundances[i]])


def _write_functional_table_csv(ft: FunctionalTable, path: Path) -> None:
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["function_id"] + ft.mag_ids)
        for i, fid in enumerate(ft.function_ids):
            w.writerow([fid] + [f"{v:g}" for v in ft.values[i]])


def _write_cazyme_csv(stats: dict, path: Path) -> None:
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["cazy_class", "total_genes", "n_families", "n_mags"])
        for cls in sorted(stats.keys()):
            s = stats[cls]
            w.writerow([cls, int(s["total_genes"]), int(s["n_families"]), int(s["n_mags"])])


def _write_redundancy_csv(redundancy: dict, path: Path) -> None:
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["group", "function_id", "n_carriers", "shannon"])
        for group, funcs in sorted(redundancy.items()):
            for func_id, info in sorted(funcs.items()):
                w.writerow([group, func_id, int(info["n_carriers"]), f"{info['shannon']:.4f}"])


def _write_diff_csv(result, g1: str, g2: str, path: Path) -> None:
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["function_id", "CLR_diff", "p_value", "q_value", "cohens_d"])
        for i, fid in enumerate(result.mag_ids):
            w.writerow([
                fid,
                f"{result.log_fold_changes[i]:.4f}",
                f"{result.p_values[i]:.6f}",
                f"{result.q_values[i]:.6f}",
                f"{result.effect_sizes[i]:.4f}",
            ])


def _write_redundancy_comparison_csv(
    result, path: Path,
) -> None:
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["group", "mean_shannon", "sd_shannon"])
        for g in result.group_names:
            w.writerow([g, f"{result.mean_shannon[g]:.4f}", f"{result.sd_shannon[g]:.4f}"])
        w.writerow([])
        w.writerow(["test", "statistic", "p_value"])
        w.writerow(["Kruskal-Wallis", f"{result.statistic:.4f}", f"{result.p_value:.6f}"])
        for (g1, g2), p in sorted(result.pairwise.items()):
            w.writerow([f"Mann-Whitney_{g1}_vs_{g2}", "", f"{p:.6f}"])


def _write_pgpr_enrichment_csv(result, path: Path) -> None:
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        header = ["trait"] + result.group_names + ["p_value", "q_value"]
        w.writerow(header)
        for trait in result.trait_names:
            row = [trait]
            row.extend([str(result.counts[trait][g]) for g in result.group_names])
            row.append(f"{result.p_values[trait]:.6f}")
            row.append(f"{result.q_values[trait]:.6f}")
            w.writerow(row)


def _embed_figure(path: Path) -> str:
    if not path.exists():
        return ""
    data = path.read_bytes()
    b64 = base64.b64encode(data).decode("ascii")
    return f'<img src="data:image/png;base64,{b64}" style="max-width:100%">'


def _write_func_html(
    out: Path,
    ko_table: FunctionalTable,
    cazy_table: FunctionalTable,
    pgpr_table: FunctionalTable,
    cazy_stats: dict,
    redundancy: dict,
    diff_results: dict,
    taxonomy: TaxonomyTable | None,
) -> None:
    """Generate self-contained HTML report for functional profiling."""
    html_parts = [
        "<!DOCTYPE html><html><head>",
        "<meta charset='utf-8'>",
        "<title>magfunc — Functional Profiling Report</title>",
        "<style>",
        "body{font-family:sans-serif;max-width:1200px;margin:0 auto;padding:20px;line-height:1.6}",
        "h1{border-bottom:3px solid #2c3e50;padding-bottom:10px}",
        "h2{color:#2c3e50;cursor:pointer;border-bottom:1px solid #ddd;padding:8px 0}",
        "table{border-collapse:collapse;width:100%;margin:10px 0}",
        "th,td{border:1px solid #ddd;padding:6px 10px;text-align:left}",
        "th{background:#f5f5f5;cursor:pointer}",
        "tr:nth-child(even){background:#fafafa}",
        ".section{margin-bottom:30px}",
        ".sig{color:#e74c3c;font-weight:bold}",
        "</style></head><body>",
        "<h1>magfunc — Functional Profiling Report</h1>",
    ]

    # Summary
    html_parts.append("<div class='section'>")
    html_parts.append("<h2>Summary</h2>")
    html_parts.append(f"<p>KO functions detected: <strong>{len(ko_table.function_ids)}</strong></p>")
    html_parts.append(f"<p>CAZyme families detected: <strong>{len(cazy_table.function_ids)}</strong></p>")
    html_parts.append(f"<p>PGPR traits screened: <strong>{len(pgpr_table.function_ids)}</strong></p>")
    html_parts.append(f"<p>MAGs annotated: <strong>{len(ko_table.mag_ids)}</strong></p>")
    html_parts.append("</div>")

    # PGPR traits
    html_parts.append("<div class='section'>")
    html_parts.append("<h2>Plant Growth-Promoting Traits (PGPR)</h2>")
    pgpr_fig = _embed_figure(out / "pgpr_traits.pdf")
    if pgpr_fig:
        html_parts.append(pgpr_fig)
    n_with_trait = int((pgpr_table.values.sum(axis=0) > 0).sum())
    html_parts.append(f"<p>{n_with_trait} of {len(pgpr_table.mag_ids)} MAGs carry at least one PGPR trait.</p>")
    html_parts.append("</div>")

    # CAZymes
    html_parts.append("<div class='section'>")
    html_parts.append("<h2>Carbohydrate-Active Enzymes (CAZymes)</h2>")
    cazy_fig = _embed_figure(out / "cazyme_bars.pdf")
    if cazy_fig:
        html_parts.append(cazy_fig)
    html_parts.append("<table><tr><th>Class</th><th>Families</th><th>Total Genes</th><th>MAGs</th></tr>")
    for cls in sorted(cazy_stats.keys()):
        s = cazy_stats[cls]
        html_parts.append(f"<tr><td>{cls}</td><td>{int(s['n_families'])}</td><td>{int(s['total_genes'])}</td><td>{int(s['n_mags'])}</td></tr>")
    html_parts.append("</table></div>")

    # Differential pathways
    if diff_results:
        html_parts.append("<div class='section'>")
        html_parts.append("<h2>Differential Pathway Abundance</h2>")
        for (g1, g2), result in sorted(diff_results.items()):
            n_sig = int((result.q_values < 0.05).sum())
            html_parts.append(f"<h3>{g1} vs {g2}: {n_sig} significant pathways (FDR &lt; 0.05)</h3>")
            if n_sig > 0:
                html_parts.append("<table><tr><th>Function</th><th>CLR diff</th><th>q-value</th><th>Cohen's d</th></tr>")
                sig_idx = np.where(result.q_values < 0.05)[0]
                sorted_idx = sig_idx[np.argsort(np.abs(result.effect_sizes[sig_idx]))[::-1]]
                for i in sorted_idx[:20]:
                    q_class = ' class="sig"' if result.q_values[i] < 0.05 else ''
                    html_parts.append(
                        f"<tr><td>{result.mag_ids[i]}</td>"
                        f"<td>{result.log_fold_changes[i]:.3f}</td>"
                        f"<td{q_class}>{result.q_values[i]:.4f}</td>"
                        f"<td>{result.effect_sizes[i]:.3f}</td></tr>"
                    )
                html_parts.append("</table>")
        html_parts.append("</div>")

    html_parts.append("</body></html>")
    (out / "func_report.html").write_text("\n".join(html_parts))
