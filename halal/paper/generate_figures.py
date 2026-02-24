#!/usr/bin/env python3
"""
SpeciesID Paper Figure Generation
Generates Figures 1–5 and Supplementary Figure S1 for speciesid_paper.md
Run from the halal/ directory: python3 paper/generate_figures.py
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
from pathlib import Path

# ── palette ──────────────────────────────────────────────────────────────────
C = {
    'green':  '#4A7C59',
    'rust':   '#C1694F',
    'gold':   '#D4A843',
    'blue':   '#5B8FA8',
    'cream':  '#F5F0E8',
    'brown':  '#8B6914',
    'sage':   '#87A878',
    'clay':   '#B8860B',
    'dark':   '#2B2B2B',
    'mid':    '#555555',
    'light':  '#999999',
}

# ── global rcParams ───────────────────────────────────────────────────────────
matplotlib.rcParams.update({
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 10,
    'axes.labelsize': 11,
    'axes.titlesize': 12,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.facecolor': C['cream'],
    'axes.facecolor': '#FAFAF5',
    'axes.edgecolor': '#666666',
    'axes.grid': True,
    'grid.alpha': 0.3,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'axes.linewidth': 0.8,
})

OUT  = Path('paper/figures')
DATA = Path('benchmark_results')
OUT.mkdir(parents=True, exist_ok=True)


# ─────────────────────────────────────────────────────────────────────────────
# Figure 1: SpeciesID Pipeline Architecture
# ─────────────────────────────────────────────────────────────────────────────
def fig1_pipeline():
    fig, ax = plt.subplots(figsize=(12, 4.5))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 4.5)
    ax.axis('off')
    fig.patch.set_facecolor(C['cream'])

    # ── helper to draw a rounded-rect box ────────────────────────────────────
    def box(cx, cy, w, h, color, text, subtext=None, fontsize=10):
        fancy = mpatches.FancyBboxPatch(
            (cx - w/2, cy - h/2), w, h,
            boxstyle="round,pad=0.08",
            facecolor=color, edgecolor='white',
            linewidth=1.5, zorder=3,
        )
        ax.add_patch(fancy)
        ax.text(cx, cy + (0.15 if subtext else 0), text,
                ha='center', va='center', fontsize=fontsize,
                fontweight='bold', color='white', zorder=4, wrap=True)
        if subtext:
            ax.text(cx, cy - 0.28, subtext,
                    ha='center', va='center', fontsize=8,
                    color='white', alpha=0.85, zorder=4)

    def arrow(x0, y0, x1, y1, label=None):
        ax.annotate('', xy=(x1, y1), xytext=(x0, y0),
                    arrowprops=dict(arrowstyle='->', color=C['mid'],
                                   lw=1.5, mutation_scale=16),
                    zorder=2)
        if label:
            mx, my = (x0+x1)/2, (y0+y1)/2 + 0.18
            ax.text(mx, my, label, ha='center', va='bottom',
                    fontsize=8, color=C['mid'], style='italic')

    # ── pipeline boxes (y = 2.25 = vertical centre) ──────────────────────────
    y = 2.25
    # 1. Input
    box(0.9, y, 1.5, 1.3, C['brown'],  'FASTQ', 'Amplicon reads')
    # 2. Stage 1 – FracMinHash
    box(3.0, y, 1.9, 1.3, C['green'],  'Stage 1', 'FracMinHash\n(k=21, s=1/1000)')
    # 3. Stage 2 – k-mer classify
    box(5.4, y, 1.9, 1.3, '#3D6B50',   'Stage 2', 'Exact k-mer\ncontainment (k=31)')
    # 4. EM core
    box(7.9, y, 2.0, 1.3, C['rust'],   'EM Algorithm', 'π, φ, λ, bias\n(Dirichlet prior)')
    # 5. Output
    box(10.7, y, 1.9, 1.3, C['blue'],  'Output', 'Weight fractions\n+ LRT / BIC')

    # ── arrows ────────────────────────────────────────────────────────────────
    arrow(1.65, y, 2.05, y, 'raw reads')
    arrow(3.95, y, 4.45, y, 'candidates')
    arrow(6.35, y, 6.9,  y, 'read table')
    arrow(8.9,  y, 9.75, y, 'estimates')

    # ── annotation bubbles (below main row) ──────────────────────────────────
    def note(cx, cy, text, color):
        ax.text(cx, cy, text, ha='center', va='center', fontsize=7.5,
                color=color, style='italic',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                          edgecolor=color, alpha=0.7, linewidth=0.8))

    note(3.0,  0.85, '~10× speed-up\nvs exact k-mer', C['green'])
    note(5.4,  0.85, 'Species + strain\ndiscrimination',  C['green'])
    note(7.9,  0.85, 'Bias correction:\nmtDNA CN, PCR, λ', C['rust'])
    note(10.7, 0.85, 'BIC model selection\n+ LRT presence test', C['blue'])

    # ── vertical connectors ───────────────────────────────────────────────────
    for cx in [3.0, 5.4, 7.9, 10.7]:
        ax.annotate('', xy=(cx, 1.18), xytext=(cx, y - 0.65),
                    arrowprops=dict(arrowstyle='->', color=C['light'],
                                   lw=0.8, mutation_scale=10), zorder=1)

    # ── title ─────────────────────────────────────────────────────────────────
    ax.set_title('Figure 1.  SpeciesID two-stage classification and bias-aware EM pipeline',
                 fontsize=11, loc='left', pad=14, color=C['dark'])

    fig.tight_layout()
    out = OUT / 'fig1_pipeline.png'
    fig.savefig(out, bbox_inches='tight', facecolor=C['cream'])
    plt.close()
    print(f'  Saved {out}')


# ─────────────────────────────────────────────────────────────────────────────
# Figure 2: Binary Mixture Quantification Accuracy (3-panel scatter)
# ─────────────────────────────────────────────────────────────────────────────
def fig2_binary_accuracy():
    pairs = [
        ('binary_beef_pork.tsv',    'Sus_scrofa',     'Beef–Pork',    'Pork'),
        ('binary_beef_horse.tsv',   'Equus_caballus', 'Beef–Horse',   'Horse'),
        ('binary_chicken_pork.tsv', 'Sus_scrofa',     'Chicken–Pork', 'Pork'),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(12, 4.2), sharey=False)
    fig.patch.set_facecolor(C['cream'])

    for ax, (fname, target_sp, pair_label, sp_short) in zip(axes, pairs):
        df = pd.read_csv(DATA / fname, sep='\t')
        df_t = df[df['species'] == target_sp].copy()

        grp = (df_t.groupby('experiment')
                   .agg(true_w=('true_w', 'first'),
                        est_mean=('est_w', 'mean'),
                        est_sd=('est_w', 'std'))
                   .reset_index()
                   .sort_values('true_w'))

        x = grp['true_w'].values * 100
        y = grp['est_mean'].values * 100
        e = grp['est_sd'].values * 100

        # identity line
        mx = max(x.max(), y.max()) * 1.12
        ax.plot([0, mx], [0, mx], '--', color=C['mid'], lw=1.2, alpha=0.6,
                zorder=1, label='Identity (y=x)')

        # ±10pp tolerance band
        ax.fill_between([0, mx], [-10, mx-10], [10, mx+10],
                        color=C['sage'], alpha=0.10, zorder=0)

        # data
        ax.errorbar(x, y, yerr=e,
                    fmt='o', color=C['rust'], ecolor=C['clay'],
                    elinewidth=1.4, capsize=3.5, capthick=1.2,
                    markersize=6.5, zorder=3,
                    label=f'{sp_short} fraction')

        # per-point annotation at low concentrations
        for xi, yi, ei in zip(x, y, e):
            if xi <= 10:
                ax.annotate(f'{yi:.1f}±{ei:.1f}',
                            xy=(xi, yi), xytext=(xi + 1.5, yi + 1.5),
                            fontsize=7, color=C['clay'],
                            arrowprops=dict(arrowstyle='-', color=C['light'],
                                           lw=0.6))

        # r² annotation
        r2 = np.corrcoef(x, y)[0, 1] ** 2
        mae = np.mean(np.abs(y - x))
        ax.text(0.05, 0.93,
                f'$R^2$ = {r2:.4f}\nMAE = {mae:.2f} pp',
                transform=ax.transAxes, fontsize=8.5,
                color=C['dark'],
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                          edgecolor=C['light'], alpha=0.8))

        ax.set_xlabel('True weight fraction (%)', fontsize=10)
        ax.set_ylabel('Estimated weight fraction (%)', fontsize=10)
        ax.set_title(pair_label, fontsize=11, fontweight='bold')
        ax.legend(fontsize=8, loc='lower right')
        ax.set_xlim(-0.5, mx)
        ax.set_ylim(-0.5, mx)

    fig.suptitle('Figure 2.  Quantification accuracy on simulated binary meat mixtures '
                 '(mean ± SD, n=3 seeds)',
                 fontsize=10.5, y=1.02)
    fig.tight_layout()
    out = OUT / 'fig2_binary_accuracy.png'
    fig.savefig(out, bbox_inches='tight', facecolor=C['cream'])
    plt.close()
    print(f'  Saved {out}')


# ─────────────────────────────────────────────────────────────────────────────
# Figure 3: mtDNA CN-driven error heteroscedasticity across species pairs
# ─────────────────────────────────────────────────────────────────────────────
def fig3_ablation():
    """
    Shows that species pairs with high mtDNA CN differential (beef-pork)
    have systematically larger quantification errors than pairs with low
    CN differential (chicken-pork), even when the adultulterant is the
    same species (Sus scrofa). This demonstrates the practical impact of
    CN bias and the need for CN correction.
    """
    # Three pairs, each measuring Sus_scrofa as the adultulterant
    pairs = [
        ('binary_beef_pork.tsv',    'Sus_scrofa', 'Beef–Pork',    C['rust'],  'high CN differential\n(cattle ~10× pig)'),
        ('binary_beef_horse.tsv',   'Equus_caballus', 'Beef–Horse', C['gold'], 'moderate CN differential'),
        ('binary_chicken_pork.tsv', 'Sus_scrofa', 'Chicken–Pork', C['green'], 'low CN differential\n(chicken ≈ pig)'),
    ]

    all_grps = []
    for fname, target_sp, label, color, cn_note in pairs:
        df = pd.read_csv(DATA / fname, sep='\t')
        df = df[df['species'] == target_sp].copy()
        grp = (df.groupby('experiment')
                 .agg(true_w=('true_w', 'first'),
                      mae_mean=('abs_error', 'mean'),
                      mae_sd=('abs_error', 'std'))
                 .reset_index()
                 .sort_values('true_w'))
        all_grps.append((grp, label, color, cn_note))

    fracs  = all_grps[0][0]['true_w'].values * 100
    labels = [f'{int(v)}%' if v == int(v) else f'{v:.0f}%' for v in fracs]
    x      = np.arange(len(fracs))
    w      = 0.25

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(13, 4.8),
                                   gridspec_kw={'width_ratios': [2, 1]})
    fig.patch.set_facecolor(C['cream'])

    # ── panel A: grouped bars ─────────────────────────────────────────────────
    offsets = [-w, 0, w]
    for (grp, lbl, col, _), offset in zip(all_grps, offsets):
        ax.bar(x + offset,
               grp['mae_mean'].values * 100,
               w,
               yerr=grp['mae_sd'].values * 100,
               error_kw=dict(elinewidth=1.1, capsize=2.5),
               color=col, alpha=0.88, label=lbl)

    ax.axhline(1.0, color=C['brown'], lw=1.5, ls='--', alpha=0.8,
               label='EU 1% labelling threshold')

    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_xlabel('True adultulterant fraction')
    ax.set_ylabel('Mean absolute error (percentage points)')
    ax.set_title('(A) MAE by species pair and mixture proportion', fontsize=10.5)
    ax.legend(fontsize=8.5, loc='upper left')

    # ── panel B: scatter — MAE at 50% vs CN differential rank ────────────────
    # pairs order: [beef-pork(0), beef-horse(1), chicken-pork(2)]
    # reorder to ascending CN rank: chicken-pork(low), beef-horse(mod), beef-pork(high)
    rank_idx  = [2, 1, 0]   # index into all_grps for low→high CN order
    cn_ranks  = [1, 2, 3]   # 1=low, 3=high CN differential
    cn_labels = ['Chicken–Pork\n(low CN diff)', 'Beef–Horse\n(mod.)',
                 'Beef–Pork\n(high CN diff)']
    mae_50    = [all_grps[i][0]['mae_mean'].values[-1] * 100 for i in rank_idx]
    sd_50     = [all_grps[i][0]['mae_sd'].values[-1]   * 100 for i in rank_idx]
    colors_50 = [C['green'], C['gold'], C['rust']]

    for rank, mae, sd, lbl, col in zip(cn_ranks, mae_50, sd_50,
                                        cn_labels, colors_50):
        ax2.errorbar(rank, mae, yerr=sd, fmt='o', color=col,
                     ecolor=col, elinewidth=1.5, capsize=4,
                     markersize=10, zorder=3)
        ax2.text(rank, mae + sd + 0.5, lbl, ha='center', va='bottom',
                 fontsize=7.5, color=col)

    ax2.set_xticks([1, 2, 3])
    ax2.set_xticklabels(['Low', 'Mod.', 'High'])
    ax2.set_xlabel('mtDNA CN differential (rank)')
    ax2.set_ylabel('MAE at 50% adultulterant (%)')
    ax2.set_title('(B) Error at 50% vs CN differential', fontsize=10.5)
    ax2.set_xlim(0.5, 3.5)
    ax2.axhline(1.0, color=C['brown'], lw=1.2, ls='--', alpha=0.6)

    # fold annotation: mae_50 now ordered [chicken-pork, beef-horse, beef-pork]
    fold = mae_50[2] / mae_50[0]  # beef-pork / chicken-pork
    ax2.annotate(f'{fold:.0f}× higher error\nin beef–pork',
                 xy=(3, mae_50[2]), xytext=(2.3, mae_50[2] - 4),
                 fontsize=8, color=C['rust'],
                 arrowprops=dict(arrowstyle='->', color=C['rust'], lw=0.9))

    fig.suptitle('Figure 3.  Mitochondrial CN differential drives quantification '
                 'error heteroscedasticity across species pairs',
                 fontsize=10.5, y=1.02)
    fig.tight_layout()
    out = OUT / 'fig3_cn_heteroscedasticity.png'
    fig.savefig(out, bbox_inches='tight', facecolor=C['cream'])
    plt.close()
    print(f'  Saved {out}')


# ─────────────────────────────────────────────────────────────────────────────
# Figure 4: Trace Detection across Sequencing Depths
# ─────────────────────────────────────────────────────────────────────────────
def fig4_trace_detection():
    df = pd.read_csv(DATA / 'trace_detection.tsv', sep='\t')

    # parse depth and trace level from experiment name, e.g. trace_0.005_500rpm
    def parse_exp(s):
        parts = s.split('_')
        trace = float(parts[1])
        depth = int(parts[2].replace('rpm', ''))
        return trace, depth

    df[['trace_level', 'depth']] = pd.DataFrame(
        df['experiment'].apply(parse_exp).tolist(), index=df.index)

    # keep only Sus_scrofa (the trace species)
    df = df[df['species'] == 'Sus_scrofa'].copy()

    depths  = sorted(df['depth'].unique())
    traces  = sorted(df['trace_level'].unique())  # [0.005, 0.01]
    colors  = [C['rust'], C['blue']]
    markers = ['o', 's']

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), sharey=False)
    fig.patch.set_facecolor(C['cream'])

    for ax, (trace_lv, color, marker) in zip(axes, zip(traces, colors, markers)):
        sub = df[df['trace_level'] == trace_lv]
        grp = (sub.groupby('depth')
                  .agg(est_mean=('est_w', 'mean'),
                       est_sd=('est_w', 'std'))
                  .reset_index())

        depth_vals = grp['depth'].values
        y_mean = grp['est_mean'].values * 100
        y_sd   = grp['est_sd'].values * 100

        # true level dashed
        ax.axhline(trace_lv * 100, color=C['mid'], ls='--', lw=1.3, alpha=0.7,
                   label=f'True = {trace_lv*100:.1f}%')

        # EU 1% threshold
        ax.axhline(1.0, color=C['gold'], ls=':', lw=1.2, label='EU 1%')

        ax.errorbar(depth_vals, y_mean, yerr=y_sd,
                    fmt=marker + '-', color=color,
                    ecolor=color, elinewidth=1.3, capsize=3.5,
                    markersize=7, linewidth=1.8,
                    label='Estimated (mean ± SD)')

        ax.set_xscale('log')
        ax.set_xticks(depths)
        ax.set_xticklabels([f'{int(d):,}' for d in depths], rotation=30)
        ax.set_xlabel('Reads per marker (log scale)')
        ax.set_ylabel('Estimated trace fraction (%)')
        ax.set_title(f'Trace level: {trace_lv*100:.1f}% ({trace_lv*1000:.0f} rpm)',
                     fontsize=10.5)
        ax.legend(fontsize=8)

        # annotation at highest depth
        best_y = y_mean[-1]
        true_y = trace_lv * 100
        ax.annotate(f'Δ = {abs(best_y - true_y):.2f} pp\n@ {depth_vals[-1]:,} reads',
                    xy=(depth_vals[-1], best_y),
                    xytext=(depth_vals[-1] * 0.5, best_y + 0.12),
                    fontsize=7.5, color=color,
                    arrowprops=dict(arrowstyle='->', color=color, lw=0.8))

    fig.suptitle('Figure 4.  Trace adultulterant detection across sequencing depths',
                 fontsize=10.5, y=1.02)
    fig.tight_layout()
    out = OUT / 'fig4_trace_detection.png'
    fig.savefig(out, bbox_inches='tight', facecolor=C['cream'])
    plt.close()
    print(f'  Saved {out}')


# ─────────────────────────────────────────────────────────────────────────────
# Figure 5: Real-Data Validation — spike-in and certified reference materials
# ─────────────────────────────────────────────────────────────────────────────
def fig5_real_data():
    det  = pd.read_csv(DATA / 'real_data_details.tsv', sep='\t')
    full = pd.read_csv(DATA / 'real_data.tsv',         sep='\t')

    # ── panel A: classified read % by sample (coloured by category) ──────────
    cat_order   = ['spike', 'spike_outofdb', 'spike_unknown', 'lgc', 'equiden', 'gemisch']
    cat_labels  = ['Spike\n(in-DB)', 'Spike\n(OOD)', 'Spike\n(unknown)',
                   'LGC\nReference', 'Equine\nAdulteration', 'Mixed\n(Gemisch)']
    cat_colors  = [C['green'], C['clay'], C['sage'], C['blue'], C['rust'], C['gold']]

    verdict_marker = {'PASS': 'o', 'FAIL': 'X', 'INCONCLUSIVE': 's'}
    verdict_color  = {'PASS': C['green'], 'FAIL': C['rust'], 'INCONCLUSIVE': C['gold']}

    fig = plt.figure(figsize=(13, 5.5))
    fig.patch.set_facecolor(C['cream'])
    gs  = gridspec.GridSpec(1, 2, width_ratios=[1.6, 1], wspace=0.35)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])

    # -- strip chart: classified % per sample, coloured by category ───────────
    x_pos = []
    for i, cat in enumerate(cat_order):
        sub = det[det['category'] == cat].reset_index(drop=True)
        for j, row in sub.iterrows():
            jitter = (j - len(sub)/2) * 0.06
            verdict = row['verdict']
            mark    = verdict_marker.get(verdict, 'o')
            col     = cat_colors[i]
            ax1.scatter(i + jitter, row['classified_pct'],
                        marker=mark, s=55, color=col,
                        edgecolors=verdict_color[verdict],
                        linewidths=1.5, zorder=3)
            x_pos.append((i + jitter, row['classified_pct']))

    ax1.set_xticks(range(len(cat_order)))
    ax1.set_xticklabels(cat_labels, fontsize=8.5)
    ax1.set_ylabel('Classified reads (%)')
    ax1.set_title('(A) Classification rate by sample type', fontsize=10.5)
    ax1.set_ylim(0, 108)
    ax1.axhline(80, color=C['mid'], ls=':', lw=1, alpha=0.5)
    ax1.text(len(cat_order) - 0.1, 81, '80% guide', fontsize=7.5,
             color=C['mid'], ha='right')

    # legend for verdict symbols
    for v, m in verdict_marker.items():
        ax1.scatter([], [], marker=m, s=45,
                    color=C['mid'], edgecolors=verdict_color[v],
                    linewidths=1.5, label=v)
    ax1.legend(title='Verdict', fontsize=8, title_fontsize=8.5)

    # -- panel B: n_species_detected distribution ─────────────────────────────
    det['cat_label'] = det['category'].map(dict(zip(cat_order, cat_labels)))
    grp = (det.groupby('category')['n_species_detected']
              .agg(['mean', 'std', 'min', 'max'])
              .reindex(cat_order))

    y = np.arange(len(cat_order))
    ax2.barh(y, grp['mean'].values,
             xerr=grp['std'].values,
             color=cat_colors, alpha=0.85, height=0.55,
             error_kw=dict(elinewidth=1.2, capsize=3))

    ax2.set_yticks(y)
    ax2.set_yticklabels(cat_labels, fontsize=8.5)
    ax2.set_xlabel('Detected species (mean ± SD)')
    ax2.set_title('(B) Species detected per sample', fontsize=10.5)
    ax2.axvline(1, color=C['mid'], ls=':', lw=1, alpha=0.5)

    fig.suptitle('Figure 5.  SpeciesID performance on independent real-world datasets',
                 fontsize=10.5, y=1.02)
    fig.tight_layout()
    out = OUT / 'fig5_real_data.png'
    fig.savefig(out, bbox_inches='tight', facecolor=C['cream'])
    plt.close()
    print(f'  Saved {out}')


# ─────────────────────────────────────────────────────────────────────────────
# Supplementary Figure S1: Runtime Scaling
# ─────────────────────────────────────────────────────────────────────────────
def figS1_runtime():
    df = pd.read_csv(DATA / 'performance.tsv', sep='\t')
    grp = (df.groupby('n_reads_per_marker')['wall_seconds']
             .agg(['mean', 'std'])
             .reset_index())

    x = grp['n_reads_per_marker'].values
    y = grp['mean'].values
    e = grp['std'].values

    # linear fit in log-log space
    log_x = np.log10(x)
    log_y = np.log10(y)
    slope, intercept = np.polyfit(log_x, log_y, 1)
    x_fit  = np.logspace(log_x.min() - 0.1, log_x.max() + 0.3, 100)
    y_fit  = 10 ** (intercept + slope * np.log10(x_fit))

    fig, ax = plt.subplots(figsize=(6, 4))
    fig.patch.set_facecolor(C['cream'])

    ax.plot(x_fit, y_fit * 1000, '--', color=C['mid'], lw=1.2, alpha=0.7,
            label=f'Power-law fit (slope = {slope:.2f})')
    ax.errorbar(x, y * 1000, yerr=e * 1000,
                fmt='o', color=C['green'], ecolor=C['brown'],
                elinewidth=1.4, capsize=4, markersize=7.5,
                zorder=3, label='Observed (mean ± SD, n=3)')

    # annotate each point
    for xi, yi, ei in zip(x, y * 1000, e * 1000):
        ax.annotate(f'{yi:.0f} ms',
                    xy=(xi, yi), xytext=(xi * 1.12, yi * 1.18),
                    fontsize=8, color=C['dark'])

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Reads per marker')
    ax.set_ylabel('Wall time (ms)')
    ax.set_title('Supplementary Figure S1.  SpeciesID runtime scaling',
                 fontsize=11, loc='left')
    ax.legend(fontsize=8.5)

    ax.text(0.05, 0.05,
            f'Sub-second for ≤1,000 reads/marker\n'
            f'~{y[x == 1000][0]*1000:.0f} ms at 1,000 reads',
            transform=ax.transAxes, fontsize=8.5, color=C['dark'],
            bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                      edgecolor=C['light'], alpha=0.8))

    fig.tight_layout()
    out = OUT / 'figS1_runtime.png'
    fig.savefig(out, bbox_inches='tight', facecolor=C['cream'])
    plt.close()
    print(f'  Saved {out}')


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────
if __name__ == '__main__':
    print('Generating SpeciesID paper figures...')
    fig1_pipeline()
    fig2_binary_accuracy()
    fig3_ablation()
    fig4_trace_detection()
    fig5_real_data()
    figS1_runtime()
    print(f'\nDone. All figures saved to {OUT}/')
