"""Self-contained HTML report generator for MAG community profiling."""

from __future__ import annotations

import base64
from pathlib import Path

import numpy as np


def _embed_figure(path: str | Path) -> str:
    """Embed a figure as base64 PNG in an img tag. Returns empty string if missing."""
    p = Path(path)
    if not p.exists():
        return ""
    data = p.read_bytes()
    b64 = base64.b64encode(data).decode("ascii")
    suffix = p.suffix.lstrip(".").lower()
    mime = "image/png" if suffix == "png" else "image/jpeg" if suffix in ("jpg", "jpeg") else "application/pdf"
    # For PDFs, skip embedding (not renderable in img)
    if mime == "application/pdf":
        # Try to find a png equivalent
        png_path = p.with_suffix(".png")
        if png_path.exists():
            data = png_path.read_bytes()
            b64 = base64.b64encode(data).decode("ascii")
            mime = "image/png"
        else:
            return f'<p class="muted">Figure saved as PDF: {p.name}</p>'
    return f'<img src="data:{mime};base64,{b64}" alt="{p.stem}" style="max-width:100%;height:auto;">'


def _sig_color(p: float) -> str:
    """Return CSS color class based on p-value significance."""
    if p < 0.001:
        return "sig-high"
    if p < 0.01:
        return "sig-med"
    if p < 0.05:
        return "sig-low"
    return "sig-ns"


def _sig_text(p: float) -> str:
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


CSS = """
:root {
    --primary: #2c3e50;
    --accent: #3498db;
    --bg: #fdfdfd;
    --card: #fff;
    --border: #e0e0e0;
    --text: #333;
    --muted: #777;
}
* { box-sizing: border-box; margin: 0; padding: 0; }
body {
    font-family: 'Segoe UI', system-ui, -apple-system, sans-serif;
    color: var(--text); background: var(--bg);
    line-height: 1.6; max-width: 1100px; margin: 0 auto; padding: 24px;
}
h1 { color: var(--primary); border-bottom: 3px solid var(--accent); padding-bottom: 8px; margin-bottom: 24px; }
h2 { color: var(--primary); margin-top: 32px; margin-bottom: 12px; cursor: pointer; }
h2::before { content: '\\25BC '; font-size: 0.7em; color: var(--muted); }
h2.collapsed::before { content: '\\25B6 '; }
.section { background: var(--card); border: 1px solid var(--border); border-radius: 8px; padding: 20px; margin-bottom: 20px; }
.section-content { overflow: hidden; transition: max-height 0.3s ease; }
.section-content.collapsed { max-height: 0; padding: 0; overflow: hidden; }
.executive { background: #eaf2f8; border-left: 4px solid var(--accent); padding: 16px 20px; margin-bottom: 24px; border-radius: 4px; }
.executive p { font-size: 1.05em; }
table { width: 100%; border-collapse: collapse; margin: 12px 0; font-size: 0.92em; }
th { background: var(--primary); color: white; padding: 8px 12px; text-align: left; cursor: pointer; user-select: none; }
th:hover { background: #34495e; }
td { padding: 8px 12px; border-bottom: 1px solid var(--border); }
tr:hover { background: #f5f5f5; }
.sig-high { background: #c0392b; color: white; padding: 2px 8px; border-radius: 3px; font-weight: bold; }
.sig-med { background: #e67e22; color: white; padding: 2px 8px; border-radius: 3px; }
.sig-low { background: #f1c40f; color: #333; padding: 2px 8px; border-radius: 3px; }
.sig-ns { background: #95a5a6; color: white; padding: 2px 8px; border-radius: 3px; }
.interpretation { background: #fef9e7; border-left: 3px solid #f39c12; padding: 12px 16px; margin: 12px 0; border-radius: 4px; font-style: italic; }
.warning { background: #fdedec; border-left: 3px solid #e74c3c; padding: 12px 16px; margin: 12px 0; border-radius: 4px; }
.figure { text-align: center; margin: 16px 0; }
.figure img { border: 1px solid var(--border); border-radius: 4px; }
.muted { color: var(--muted); font-size: 0.9em; }
.methods { background: #f9f9f9; font-size: 0.9em; line-height: 1.5; }
@media print { .section { break-inside: avoid; } h2::before { content: ''; } }
@media (max-width: 768px) { body { padding: 12px; } table { font-size: 0.8em; } }
"""

JS = """
document.addEventListener('DOMContentLoaded', function() {
    // Collapsible sections
    document.querySelectorAll('h2[data-toggle]').forEach(function(h2) {
        h2.addEventListener('click', function() {
            var target = document.getElementById(h2.getAttribute('data-toggle'));
            if (target) {
                target.classList.toggle('collapsed');
                h2.classList.toggle('collapsed');
            }
        });
    });
    // Sortable tables
    document.querySelectorAll('th[data-sort]').forEach(function(th) {
        th.addEventListener('click', function() {
            var table = th.closest('table');
            var idx = Array.from(th.parentNode.children).indexOf(th);
            var tbody = table.querySelector('tbody');
            var rows = Array.from(tbody.querySelectorAll('tr'));
            var dir = th.getAttribute('data-dir') === 'asc' ? 'desc' : 'asc';
            th.setAttribute('data-dir', dir);
            rows.sort(function(a, b) {
                var av = a.children[idx].textContent.trim();
                var bv = b.children[idx].textContent.trim();
                var an = parseFloat(av), bn = parseFloat(bv);
                if (!isNaN(an) && !isNaN(bn)) { return dir === 'asc' ? an - bn : bn - an; }
                return dir === 'asc' ? av.localeCompare(bv) : bv.localeCompare(av);
            });
            rows.forEach(function(r) { tbody.appendChild(r); });
        });
    });
});
"""


def generate_html_report(output_dir: str, report_data: dict) -> str:
    """Generate a self-contained HTML report with embedded figures and narratives.

    Returns the path to the generated HTML file.
    """
    out = Path(output_dir)
    sections = []

    # Extract data
    grouping_var = report_data["grouping_var"]
    group_names = report_data["group_names"]
    alpha = report_data["alpha"]
    metadata = report_data["metadata"]
    perm = report_data["permanova"]
    anos = report_data["anosim"]
    pd_result = report_data["permdisp"]
    pw_results = report_data["pairwise_permanova"]
    diff_results = report_data["differential"]
    ind = report_data["indicator"]
    core = report_data["core_microbiome"]
    cs = report_data["specificity"]
    n_before = report_data["n_mags_before"]
    n_after = report_data["n_mags_after"]
    min_prev = report_data["min_prevalence"]

    # --- Executive Summary ---
    exec_lines = []
    exec_lines.append(
        f"This report analyzes <strong>{n_after} MAGs</strong> across "
        f"<strong>{len(alpha.sample_ids)} samples</strong> in "
        f"<strong>{len(group_names)} {grouping_var} groups</strong> "
        f"({', '.join(group_names)})."
    )
    if n_before != n_after:
        exec_lines.append(
            f" Prevalence filtering (>={min_prev*100:.0f}%) removed "
            f"{n_before - n_after} of {n_before} MAGs."
        )
    if perm.p_value < 0.05:
        exec_lines.append(
            f" Community composition differs significantly across {grouping_var} "
            f"(PERMANOVA p={perm.p_value:.4f}, R<sup>2</sup>={perm.metadata['R2']:.3f})."
        )
    else:
        exec_lines.append(
            f" No significant compositional differences detected "
            f"(PERMANOVA p={perm.p_value:.4f})."
        )
    if pd_result.p_value < 0.05:
        exec_lines.append(
            " Note: significant dispersion differences detected (PERMDISP "
            f"p={pd_result.p_value:.4f}); PERMANOVA result may partly reflect "
            "dispersion rather than location effects."
        )
    n_total_sig = sum(
        int(np.sum(d.q_values < 0.05)) for _, _, d in diff_results
    )
    if n_total_sig > 0:
        exec_lines.append(
            f" A total of <strong>{n_total_sig} differentially abundant MAGs</strong> "
            "were identified across all pairwise comparisons (FDR < 0.05)."
        )
    n_core_90 = len(core.shared_across_groups.get(0.9, []))
    if n_core_90 > 0:
        exec_lines.append(
            f" <strong>{n_core_90} MAGs</strong> form the universal core "
            "(present in 90%+ of all samples across all groups)."
        )

    sections.append(
        '<div class="executive"><p>' + " ".join(exec_lines) + "</p></div>"
    )

    # --- 1. Alpha Diversity ---
    sid = "alpha"
    s = f'<h2 data-toggle="{sid}">Alpha Diversity</h2>\n<div class="section" id="{sid}"><div class="section-content">'

    # Figure
    fig_path = out / "alpha_diversity.pdf"
    s += f'<div class="figure">{_embed_figure(fig_path)}</div>'

    # Table: per-group summary
    s += "<table><thead><tr>"
    s += '<th data-sort="0">Group</th><th data-sort="1">N</th>'
    s += '<th data-sort="2">Shannon (mean)</th><th data-sort="3">Simpson (mean)</th>'
    s += '<th data-sort="4">Richness (mean)</th><th data-sort="5">Evenness (mean)</th>'
    s += "</tr></thead><tbody>"
    groups = metadata.get_groups(grouping_var)
    sample_to_group = {}
    for g, sids in groups.items():
        for sid_val in sids:
            sample_to_group[sid_val] = g
    for g in group_names:
        mask = np.array([sample_to_group.get(sid_val) == g for sid_val in alpha.sample_ids])
        n = int(mask.sum())
        s += f"<tr><td>{g}</td><td>{n}</td>"
        s += f"<td>{alpha.shannon[mask].mean():.3f}</td>"
        s += f"<td>{alpha.simpson[mask].mean():.3f}</td>"
        s += f"<td>{alpha.richness[mask].mean():.1f}</td>"
        s += f"<td>{alpha.evenness[mask].mean():.3f}</td></tr>"
    s += "</tbody></table>"

    # Interpretation
    group_medians = {g: np.median(alpha.shannon[np.array([sample_to_group.get(sid_val) == g for sid_val in alpha.sample_ids])]) for g in group_names}
    highest = max(group_medians, key=group_medians.get)
    lowest = min(group_medians, key=group_medians.get)
    s += f'<div class="interpretation">Highest median Shannon diversity: <strong>{highest}</strong> ({group_medians[highest]:.3f}); '
    s += f'lowest: <strong>{lowest}</strong> ({group_medians[lowest]:.3f}).</div>'

    s += "</div></div>"
    sections.append(s)

    # --- 2. Beta Diversity & Ordination ---
    sid = "beta"
    s = f'<h2 data-toggle="{sid}">Beta Diversity &amp; Ordination</h2>\n<div class="section" id="{sid}"><div class="section-content">'

    fig_path = out / "ordination_pcoa.pdf"
    s += f'<div class="figure">{_embed_figure(fig_path)}</div>'

    # Stats table
    s += "<table><thead><tr>"
    s += '<th data-sort="0">Test</th><th data-sort="1">Statistic</th>'
    s += '<th data-sort="2">p-value</th><th data-sort="3">Significance</th>'
    s += '<th>Details</th></tr></thead><tbody>'

    for test, label in [(perm, "PERMANOVA"), (anos, "ANOSIM"), (pd_result, "PERMDISP")]:
        sig_cls = _sig_color(test.p_value)
        extra = ""
        if "R2" in test.metadata:
            extra = f"R<sup>2</sup>={test.metadata['R2']:.3f}"
        s += f"<tr><td>{label}</td><td>{test.statistic:.4f}</td>"
        s += f'<td>{test.p_value:.4f}</td><td><span class="{sig_cls}">{_sig_text(test.p_value)}</span></td>'
        s += f"<td>{extra}</td></tr>"
    s += "</tbody></table>"

    # Interpretations
    if perm.p_value < 0.05:
        s += f'<div class="interpretation">Community composition differs significantly across {grouping_var} '
        s += f'(PERMANOVA F={perm.statistic:.3f}, R<sup>2</sup>={perm.metadata["R2"]:.3f}, p={perm.p_value:.4f}). '
        s += f'The grouping variable explains {perm.metadata["R2"]*100:.1f}% of community variation.</div>'
    if pd_result.p_value < 0.05:
        s += '<div class="warning"><strong>Warning:</strong> Significant differences in group dispersions '
        s += f'(PERMDISP F={pd_result.statistic:.3f}, p={pd_result.p_value:.4f}). '
        s += "The PERMANOVA result may partly reflect dispersion differences rather than true location shifts.</div>"

    s += "</div></div>"
    sections.append(s)

    # --- 3. Pairwise Comparisons ---
    if pw_results:
        sid = "pairwise"
        s = f'<h2 data-toggle="{sid}">Pairwise PERMANOVA</h2>\n<div class="section" id="{sid}"><div class="section-content">'
        s += "<table><thead><tr>"
        s += '<th data-sort="0">Group 1</th><th data-sort="1">Group 2</th>'
        s += '<th data-sort="2">F</th><th data-sort="3">p-value</th>'
        s += '<th data-sort="4">q-value (FDR)</th><th data-sort="5">R&sup2;</th>'
        s += '<th>Significance</th></tr></thead><tbody>'
        for r in pw_results:
            g1 = r.metadata.get("group1", "")
            g2 = r.metadata.get("group2", "")
            q = r.metadata.get("q_value", r.p_value)
            sig_cls = _sig_color(q)
            r2 = r.metadata.get("R2", 0)
            s += f"<tr><td>{g1}</td><td>{g2}</td><td>{r.statistic:.4f}</td>"
            s += f"<td>{r.p_value:.4f}</td><td>{q:.4f}</td><td>{r2:.3f}</td>"
            s += f'<td><span class="{sig_cls}">{_sig_text(q)}</span></td></tr>'
        s += "</tbody></table>"

        sig_pairs = [r for r in pw_results if r.metadata.get("q_value", 1) < 0.05]
        if sig_pairs:
            pair_strs = [f"{r.metadata['group1']} vs {r.metadata['group2']}" for r in sig_pairs]
            s += f'<div class="interpretation">Significant pairwise differences (FDR < 0.05): {", ".join(pair_strs)}.</div>'
        else:
            s += '<div class="interpretation">No pairwise comparisons reach significance after FDR correction.</div>'

        s += "</div></div>"
        sections.append(s)

    # --- 4. Differential Abundance ---
    if diff_results:
        sid = "differential"
        s = f'<h2 data-toggle="{sid}">Differential Abundance</h2>\n<div class="section" id="{sid}"><div class="section-content">'

        for g1, g2, diff in diff_results:
            n_sig = int(np.sum(diff.q_values < 0.05))
            n_up = int(np.sum((diff.q_values < 0.05) & (diff.log_fold_changes > 0)))
            n_down = int(np.sum((diff.q_values < 0.05) & (diff.log_fold_changes < 0)))

            fig_path = out / f"volcano_{g1}_vs_{g2}.pdf"
            s += f"<h3>{g1} vs {g2}</h3>"
            s += f'<div class="figure">{_embed_figure(fig_path)}</div>'
            s += f'<div class="interpretation">{n_sig} MAGs differentially abundant: '
            s += f"{n_up} enriched in {g1}, {n_down} enriched in {g2}.</div>"

            # Top significant MAGs table
            if n_sig > 0:
                sig_idx = np.where(diff.q_values < 0.05)[0]
                top_idx = sig_idx[np.argsort(np.abs(diff.effect_sizes[sig_idx]))[::-1][:10]]
                s += "<table><thead><tr>"
                s += '<th data-sort="0">MAG</th><th data-sort="1">CLR Diff</th>'
                s += '<th data-sort="2">q-value</th><th data-sort="3">Cohen\'s d</th>'
                s += "</tr></thead><tbody>"
                for i in top_idx:
                    s += f"<tr><td>{diff.mag_ids[i]}</td>"
                    s += f"<td>{diff.log_fold_changes[i]:.3f}</td>"
                    s += f"<td>{diff.q_values[i]:.4f}</td>"
                    s += f"<td>{diff.effect_sizes[i]:.3f}</td></tr>"
                s += "</tbody></table>"

        # Heatmap
        heatmap_path = out / "heatmap_significant.pdf"
        if not heatmap_path.exists():
            heatmap_path = out / "heatmap_all.pdf"
        if heatmap_path.exists():
            s += f'<div class="figure">{_embed_figure(heatmap_path)}</div>'

        s += "</div></div>"
        sections.append(s)

    # --- 5. Indicator Species ---
    sid = "indicator"
    s = f'<h2 data-toggle="{sid}">Indicator Species</h2>\n<div class="section" id="{sid}"><div class="section-content">'

    fig_path = out / "indicator_species.pdf"
    s += f'<div class="figure">{_embed_figure(fig_path)}</div>'

    # Top 10 per group table
    s += "<table><thead><tr>"
    s += '<th data-sort="0">MAG</th><th data-sort="1">Group</th>'
    s += '<th data-sort="2">IndVal</th><th data-sort="3">p-value</th>'
    s += '<th>Significance</th></tr></thead><tbody>'

    for g in group_names:
        mask = [bg == g for bg in ind.best_group]
        indices = [i for i, m in enumerate(mask) if m]
        sorted_idx = sorted(indices, key=lambda i: ind.indval_scores[i], reverse=True)[:10]
        for i in sorted_idx:
            sig_cls = _sig_color(ind.p_values[i])
            s += f"<tr><td>{ind.mag_ids[i]}</td><td>{g}</td>"
            s += f"<td>{ind.indval_scores[i]:.3f}</td>"
            s += f'<td>{ind.p_values[i]:.4f}</td>'
            s += f'<td><span class="{sig_cls}">{_sig_text(ind.p_values[i])}</span></td></tr>'
    s += "</tbody></table>"
    s += "</div></div>"
    sections.append(s)

    # --- 6. Core Microbiome ---
    sid = "core"
    s = f'<h2 data-toggle="{sid}">Core Microbiome</h2>\n<div class="section" id="{sid}"><div class="section-content">'

    s += "<table><thead><tr>"
    s += '<th data-sort="0">Threshold</th>'
    for g in group_names:
        s += f'<th data-sort="{group_names.index(g)+1}">{g}</th>'
    s += '<th data-sort="last">Shared (all groups)</th></tr></thead><tbody>'

    for t in core.thresholds:
        s += f"<tr><td>{t*100:.0f}%</td>"
        for g in group_names:
            n = len(core.core_mags_per_threshold[t].get(g, []))
            s += f"<td>{n}</td>"
        n_shared = len(core.shared_across_groups[t])
        s += f"<td><strong>{n_shared}</strong></td></tr>"
    s += "</tbody></table>"

    n_core_90 = len(core.shared_across_groups.get(0.9, []))
    n_core_50 = len(core.shared_across_groups.get(0.5, []))
    s += f'<div class="interpretation">{n_core_90} MAGs form the strict core (present in >=90% of samples in all groups). '
    s += f'{n_core_50} MAGs form the relaxed core (>=50% threshold).</div>'
    s += "</div></div>"
    sections.append(s)

    # --- 7. Compartment Specificity ---
    sid = "specificity"
    s = f'<h2 data-toggle="{sid}">Compartment Specificity</h2>\n<div class="section" id="{sid}"><div class="section-content">'

    specialists = int(np.sum(cs.scores > 0.5))
    generalists = int(np.sum(cs.scores <= 0.5))
    s += f"<p>Specialists (score > 0.5): <strong>{specialists}</strong> | "
    s += f"Generalists (score <= 0.5): <strong>{generalists}</strong></p>"

    # Top specialists table
    top_spec_idx = np.argsort(cs.scores)[::-1][:15]
    s += "<table><thead><tr>"
    s += '<th data-sort="0">MAG</th><th data-sort="1">Specificity Score</th>'
    s += '<th data-sort="2">Dominant Compartment</th></tr></thead><tbody>'
    for i in top_spec_idx:
        s += f"<tr><td>{cs.mag_ids[i]}</td>"
        s += f"<td>{cs.scores[i]:.3f}</td>"
        s += f"<td>{cs.dominant_compartment[i]}</td></tr>"
    s += "</tbody></table>"
    s += "</div></div>"
    sections.append(s)

    # --- 8. Taxonomy ---
    tax_fig = out / "taxonomy_bars.pdf"
    tax_grp_fig = out / "taxonomy_bars_grouped.pdf"
    if tax_fig.exists() or tax_grp_fig.exists():
        sid = "taxonomy"
        s = f'<h2 data-toggle="{sid}">Taxonomy</h2>\n<div class="section" id="{sid}"><div class="section-content">'
        if tax_grp_fig.exists():
            s += f'<div class="figure">{_embed_figure(tax_grp_fig)}</div>'
        if tax_fig.exists():
            s += f'<div class="figure">{_embed_figure(tax_fig)}</div>'
        s += "</div></div>"
        sections.append(s)

    # --- 9. Methods ---
    sid = "methods"
    s = f'<h2 data-toggle="{sid}">Methods</h2>\n<div class="section methods" id="{sid}"><div class="section-content">'
    s += "<p>MAG community profiling was performed using magprofile. "
    if min_prev > 0:
        s += f"MAGs present in fewer than {min_prev*100:.0f}% of samples were removed prior to analysis "
        s += f"(retaining {n_after} of {n_before} MAGs). "
    s += "Alpha diversity was quantified using Shannon entropy, Simpson's diversity index, "
    s += "observed richness, and Pielou's evenness. "
    s += "Beta diversity was computed using Bray-Curtis dissimilarity and visualized "
    s += "via Principal Coordinates Analysis (PCoA) with 95% confidence ellipses. "
    s += "Differences in community composition were assessed using PERMANOVA and ANOSIM "
    s += f"with {perm.metadata.get('n_permutations', 999)} permutations. "
    s += "Homogeneity of multivariate dispersions was tested using PERMDISP. "
    s += "Pairwise group comparisons were performed with PERMANOVA and p-values were "
    s += "adjusted using Benjamini-Hochberg FDR correction. "
    s += "Differential abundance was assessed using centered log-ratio (CLR) transformation "
    s += "followed by Welch's t-test with BH-FDR correction. "
    s += "Indicator species were identified using the IndVal index with permutation tests. "
    s += "Core microbiome was defined as MAGs present in at least 50-100% of samples "
    s += "within each group. "
    s += "Compartment specificity was quantified using an information-theoretic approach "
    s += "(1 - H/H_max).</p>"
    s += "</div></div>"
    sections.append(s)

    # --- Assemble HTML ---
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>MAG Community Profiling Report</title>
<style>{CSS}</style>
</head>
<body>
<h1>MAG Community Profiling Report</h1>
{"".join(sections)}
<script>{JS}</script>
</body>
</html>"""

    report_path = out / "report.html"
    report_path.write_text(html, encoding="utf-8")
    return str(report_path)
