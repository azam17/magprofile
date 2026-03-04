#!/usr/bin/env python3
"""
Figure 2 — SpeciesID Methodology Detail (Food Chemistry)
Detailed workflow figure showing:
  Panel A: Two-stage classification with concrete example
  Panel B: Bias correction concept (raw reads vs corrected weights)
  Panel C: Output and regulatory verdict

Run from halal/: python3 paper/generate_fig2_methodology.py
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from pathlib import Path

OUT = Path('paper/figures')
OUT.mkdir(parents=True, exist_ok=True)

PAL = {
    'text':     '#1A1A1A',
    'grey':     '#666666',
    'ltgrey':   '#E0E0E0',
    'midgrey':  '#AAAAAA',
    'input':    '#6B4C3B',
    'stage1':   '#4A7C59',
    'stage2':   '#3A6B5A',
    'em':       '#B85C3E',
    'output':   '#4A7896',
    'beef':     '#C1694F',
    'pork':     '#7B68AE',
    'chicken':  '#D4A843',
    'raw':      '#CC6666',
    'corrected':'#4A7C59',
    'arrow':    '#888888',
    'panel_bg': '#F9F9F9',
    'panel_bdr':'#DDDDDD',
    'haram_bg': '#FFEBEE',
    'haram':    '#C62828',
}

matplotlib.rcParams.update({
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'figure.facecolor': 'white',
    'axes.facecolor': 'white',
})


def box(ax, xy, w, h, fc, ec=None, lw=1.5, alpha=1.0, zorder=2):
    if ec is None:
        ec = fc
    p = mpatches.FancyBboxPatch(xy, w, h, boxstyle="round,pad=0.06",
                                 facecolor=fc, edgecolor=ec,
                                 linewidth=lw, alpha=alpha, zorder=zorder)
    ax.add_patch(p)


def arrow(ax, x0, y0, x1, y1, c='#888888', lw=1.5):
    ax.annotate('', xy=(x1, y1), xytext=(x0, y0),
                arrowprops=dict(arrowstyle='->', color=c, lw=lw,
                                mutation_scale=14))


def fig2():
    fig = plt.figure(figsize=(17, 22))
    fig.patch.set_facecolor('white')

    gs = fig.add_gridspec(3, 1, height_ratios=[1.0, 0.52, 0.52],
                          hspace=0.06, left=0.03, right=0.97, top=0.975, bottom=0.015)

    # ═══════════════════════════════════════════════════════════════════
    # PANEL A: Two-stage classification
    # ═══════════════════════════════════════════════════════════════════
    ax = fig.add_subplot(gs[0])
    ax.set_xlim(0, 17)
    ax.set_ylim(0, 11)
    ax.axis('off')

    ax.text(0.15, 10.7, 'A', fontsize=20, fontweight='bold', color=PAL['text'])
    ax.text(0.7, 10.7, 'Two-stage read classification', fontsize=14,
            fontweight='bold', color=PAL['text'])

    # ── Row 1: Input flow ────────────────────────────────────────────
    # Meat sample
    box(ax, (0.3, 8.2), 2.6, 2.0, '#F5F0E8', ec=PAL['input'], lw=2)
    ax.text(1.6, 9.7, 'Meat Sample', fontsize=11, fontweight='bold',
            ha='center', color=PAL['input'])
    ax.text(1.6, 9.2, '50% beef + 50% pork', fontsize=9,
            ha='center', color=PAL['grey'])
    ax.text(1.6, 8.8, '(true composition)', fontsize=8,
            ha='center', color=PAL['midgrey'], style='italic')

    arrow(ax, 3.0, 9.2, 3.7, 9.2, PAL['arrow'])

    # Sequencing
    box(ax, (3.8, 8.2), 3.0, 2.0, '#F0F4F8', ec=PAL['output'], lw=2)
    ax.text(5.3, 9.7, 'DNA Sequencing', fontsize=11, fontweight='bold',
            ha='center', color=PAL['output'])
    ax.text(5.3, 9.2, '1,500 short DNA reads', fontsize=9,
            ha='center', color=PAL['grey'])
    ax.text(5.3, 8.8, '(each ~150 bases long)', fontsize=8,
            ha='center', color=PAL['midgrey'], style='italic')

    arrow(ax, 6.9, 9.2, 7.6, 9.2, PAL['arrow'])

    # Reads
    box(ax, (7.7, 8.2), 3.8, 2.0, '#FFF8F0', ec=PAL['em'], lw=2)
    ax.text(9.6, 9.7, 'Each Read Processed Individually',
            fontsize=10, fontweight='bold', ha='center', color=PAL['em'])
    # Read bars
    rcolors = [PAL['beef'], PAL['pork'], PAL['beef'], PAL['pork'],
               PAL['beef'], PAL['pork'], PAL['beef'], PAL['pork']]
    for i, c in enumerate(rcolors):
        x = 8.0 + i * 0.42
        box(ax, (x, 8.7), 0.32, 0.55, c, alpha=0.7, lw=0.8)
        ax.text(x + 0.16, 8.98, f'R{i+1}', fontsize=6, ha='center',
                va='center', color='white', fontweight='bold')
    ax.text(9.6, 8.45, '... x 1,500 reads', fontsize=8,
            ha='center', color=PAL['midgrey'], style='italic')

    arrow(ax, 11.6, 9.2, 12.3, 9.2, PAL['arrow'])

    # Reference DB
    box(ax, (12.4, 8.2), 4.0, 2.0, '#F0F8F0', ec=PAL['stage1'], lw=2)
    ax.text(14.4, 9.8, 'Reference Database', fontsize=11, fontweight='bold',
            ha='center', color=PAL['stage1'])
    spp = ['Cattle  (2000 CN)', 'Pig  (1800 CN)', 'Horse  (1500 CN)',
           'Chicken  (1000 CN)', 'Sheep  (1700 CN)', 'Goat  (1600 CN)']
    for i, sp in enumerate(spp):
        col = i % 2
        row = i // 2
        ax.text(12.7 + col * 2.0, 9.3 - row * 0.3, sp, fontsize=7.5,
                color=PAL['grey'])
    ax.text(14.4, 8.35, '19 species x 3 markers', fontsize=8,
            ha='center', color=PAL['midgrey'], style='italic')

    # ── Row 2: STAGE 1 ──────────────────────────────────────────────
    box(ax, (0.3, 4.6), 7.6, 3.1, PAL['stage1'], alpha=0.06,
        ec=PAL['stage1'], lw=2)
    ax.text(0.55, 7.4, 'STAGE 1', fontsize=12, fontweight='bold',
            color=PAL['stage1'])
    ax.text(2.3, 7.4, 'Rapid Screening  (FracMinHash, k = 21)',
            fontsize=11, color=PAL['stage1'])

    # Step boxes
    # 1) k-mers
    box(ax, (0.5, 5.4), 2.2, 1.6, 'white', ec=PAL['stage1'], lw=1.2)
    ax.text(1.6, 6.75, '1. Cut read into', fontsize=9, fontweight='bold',
            ha='center', color=PAL['stage1'])
    ax.text(1.6, 6.4, '21-base fragments', fontsize=9, fontweight='bold',
            ha='center', color=PAL['stage1'])
    ax.text(1.6, 5.95, 'ACGTACGTACGTACGTACGTG', fontsize=6, ha='center',
            color=PAL['grey'], family='monospace')
    ax.text(1.6, 5.7, 'CGTACGTACGTACGTACGTGA', fontsize=6, ha='center',
            color=PAL['midgrey'], family='monospace')
    ax.text(1.6, 5.45, '...130 fragments total', fontsize=7, ha='center',
            color=PAL['midgrey'], style='italic')

    arrow(ax, 2.8, 6.2, 3.2, 6.2, PAL['stage1'])

    # 2) Hash & sketch
    box(ax, (3.3, 5.4), 2.2, 1.6, 'white', ec=PAL['stage1'], lw=1.2)
    ax.text(4.4, 6.75, '2. Hash & keep', fontsize=9, fontweight='bold',
            ha='center', color=PAL['stage1'])
    ax.text(4.4, 6.4, 'smallest values', fontsize=9, fontweight='bold',
            ha='center', color=PAL['stage1'])
    ax.text(4.4, 5.95, '130 fragments', fontsize=8, ha='center',
            color=PAL['grey'])
    ax.text(4.4, 5.65, '  --> 3-5 "fingerprints"', fontsize=8, ha='center',
            color=PAL['stage1'], fontweight='bold')
    ax.text(4.4, 5.4, '(compressed signature)', fontsize=7.5, ha='center',
            color=PAL['midgrey'], style='italic')

    arrow(ax, 5.6, 6.2, 6.0, 6.2, PAL['stage1'])

    # 3) Compare
    box(ax, (6.1, 5.4), 1.6, 1.6, 'white', ec=PAL['stage1'], lw=1.2)
    ax.text(6.9, 6.75, '3. Compare', fontsize=9, fontweight='bold',
            ha='center', color=PAL['stage1'])
    ax.text(6.9, 6.4, 'vs 19 species', fontsize=9, fontweight='bold',
            ha='center', color=PAL['stage1'])
    # Mini species list with match indicator
    matches = [('Cattle', True), ('Pig', True), ('Horse', False),
               ('Chicken', False), ('...x15 more', False)]
    for i, (name, hit) in enumerate(matches):
        y = 6.0 - i * 0.14
        color = PAL['stage1'] if hit else PAL['ltgrey']
        marker = 'MATCH' if hit else '--'
        ax.text(6.2, y, name, fontsize=6, color=PAL['grey'] if hit else PAL['midgrey'])
        ax.text(7.5, y, marker, fontsize=5.5, ha='right',
                color=PAL['stage1'] if hit else PAL['ltgrey'],
                fontweight='bold' if hit else 'normal')

    # Stage 1 result
    box(ax, (0.5, 4.75), 7.2, 0.5, '#E8F5E9', ec=PAL['stage1'], lw=1)
    ax.text(4.1, 5.0,
            'Result:  17/19 species eliminated  -->  2 candidates remain (Cattle + Pig)',
            fontsize=9.5, ha='center', va='center', color=PAL['stage1'],
            fontweight='bold')

    # ── Row 2: STAGE 2 ──────────────────────────────────────────────
    box(ax, (8.3, 4.6), 8.2, 3.1, PAL['stage2'], alpha=0.06,
        ec=PAL['stage2'], lw=2)
    ax.text(8.55, 7.4, 'STAGE 2', fontsize=12, fontweight='bold',
            color=PAL['stage2'])
    ax.text(10.3, 7.4, 'Precise Scoring  (Exact k-mers, k = 31)',
            fontsize=11, color=PAL['stage2'])

    # 4) Full comparison
    box(ax, (8.5, 5.4), 3.2, 1.6, 'white', ec=PAL['stage2'], lw=1.2)
    ax.text(10.1, 6.75, '4. Compare full read', fontsize=9, fontweight='bold',
            ha='center', color=PAL['stage2'])
    ax.text(10.1, 6.4, 'vs each candidate', fontsize=9, fontweight='bold',
            ha='center', color=PAL['stage2'])
    ax.text(10.1, 5.95, 'Count matching 31-base', fontsize=8,
            ha='center', color=PAL['grey'])
    ax.text(10.1, 5.7, 'subsequences in reference', fontsize=8,
            ha='center', color=PAL['grey'])
    ax.text(10.1, 5.4, '(exact DNA letter match)', fontsize=7.5,
            ha='center', color=PAL['midgrey'], style='italic')

    arrow(ax, 11.8, 6.2, 12.3, 6.2, PAL['stage2'])

    # 5) Containment scores
    box(ax, (12.4, 5.1), 3.8, 1.9, 'white', ec=PAL['stage2'], lw=1.2)
    ax.text(14.3, 6.75, '5. Containment Scores', fontsize=9.5, fontweight='bold',
            ha='center', color=PAL['stage2'])
    ax.text(14.3, 6.35, '(fraction of read matching reference)',
            fontsize=7.5, ha='center', color=PAL['midgrey'], style='italic')

    # Cattle bar
    bw = 0.85 * 2.8
    ax.barh(5.9, bw, height=0.35, left=12.65,
            color=PAL['beef'], alpha=0.85, zorder=3)
    ax.text(12.55, 5.9, 'Cattle', fontsize=8, ha='right', va='center',
            color=PAL['beef'], fontweight='bold')
    ax.text(12.65 + bw + 0.1, 5.9, '0.85', fontsize=9,
            ha='left', va='center', color=PAL['beef'], fontweight='bold')

    # Pig bar
    bw2 = 0.12 * 2.8
    ax.barh(5.4, bw2, height=0.35, left=12.65,
            color=PAL['pork'], alpha=0.85, zorder=3)
    ax.text(12.55, 5.4, 'Pig', fontsize=8, ha='right', va='center',
            color=PAL['pork'], fontweight='bold')
    ax.text(12.65 + bw2 + 0.1, 5.4, '0.12', fontsize=9,
            ha='left', va='center', color=PAL['pork'], fontweight='bold')

    # Scale
    ax.plot([12.65, 15.45], [5.15, 5.15], color=PAL['ltgrey'], lw=0.8)
    for v in [0, 0.5, 1.0]:
        x = 12.65 + v * 2.8
        ax.plot([x, x], [5.1, 5.15], color=PAL['ltgrey'], lw=0.8)
        ax.text(x, 5.02, f'{v:.0%}', fontsize=6.5, ha='center',
                color=PAL['midgrey'])

    # Stage 2 result
    box(ax, (8.5, 4.75), 7.8, 0.5, '#E0F2F1', ec=PAL['stage2'], lw=1)
    ax.text(12.4, 5.0,
            'Result:  This read is most likely from Cattle (score 0.85 >> 0.12)',
            fontsize=9.5, ha='center', va='center', color=PAL['stage2'],
            fontweight='bold')

    # ── Summary bar ──────────────────────────────────────────────────
    box(ax, (0.3, 3.5), 16.2, 0.8, '#F5F5F5', ec=PAL['ltgrey'], lw=1)
    steps = [
        ('1,500 reads', 1.8, PAL['input']),
        ('Stage 1: filter\n19 --> 2 species/read', 5.2, PAL['stage1']),
        ('Stage 2: score each\nread vs candidates', 9.2, PAL['stage2']),
        ('1,500 reads x 2 scores\npassed to EM algorithm', 13.5, PAL['em']),
    ]
    for label, x, color in steps:
        ax.text(x, 3.9, label, fontsize=8.5, ha='center', va='center',
                color=color, fontweight='bold', linespacing=1.2)
    for x in [3.4, 7.2, 11.3]:
        ax.text(x, 3.9, '-->', fontsize=14, ha='center', va='center',
                color=PAL['arrow'], fontweight='bold')

    ax.text(8.4, 3.15,
            'Note: Each of 1,500 reads is processed independently. '
            'Stage 1 eliminates >95% of species per read; '
            'Stage 2 produces precise similarity scores.',
            fontsize=8, ha='center', va='center', color=PAL['midgrey'],
            style='italic')

    # ═══════════════════════════════════════════════════════════════════
    # PANEL B: Bias correction
    # ═══════════════════════════════════════════════════════════════════
    ax2 = fig.add_subplot(gs[1])
    ax2.set_xlim(0, 17)
    ax2.set_ylim(0, 6.5)
    ax2.axis('off')

    ax2.text(0.15, 6.3, 'B', fontsize=20, fontweight='bold', color=PAL['text'])
    ax2.text(0.7, 6.3, 'Why bias correction is needed: raw read counts do not equal true weight',
             fontsize=14, fontweight='bold', color=PAL['text'])

    # ── Left: The problem ────────────────────────────────────────────
    box(ax2, (0.3, 0.2), 5.5, 5.5, PAL['panel_bg'], ec=PAL['panel_bdr'], lw=1.5)
    ax2.text(3.05, 5.4, 'The Problem', fontsize=12, fontweight='bold',
             ha='center', color=PAL['raw'])

    # True composition
    ax2.text(0.6, 4.8, 'True weight (by scale):', fontsize=10,
             color=PAL['text'], fontweight='bold')
    ax2.barh(4.3, 2.3, height=0.45, left=0.6,
             color=PAL['beef'], alpha=0.85)
    ax2.barh(4.3, 2.3, height=0.45, left=2.9,
             color=PAL['pork'], alpha=0.85)
    ax2.text(1.75, 4.3, 'Beef 50%', fontsize=9, ha='center', va='center',
             color='white', fontweight='bold')
    ax2.text(4.05, 4.3, 'Pork 50%', fontsize=9, ha='center', va='center',
             color='white', fontweight='bold')

    # Bias explanation
    ax2.text(0.6, 3.65, 'But cattle cells contain 2,000 mtDNA copies', fontsize=9,
             color=PAL['grey'])
    ax2.text(0.6, 3.3, 'while pig cells contain 1,800 copies', fontsize=9,
             color=PAL['grey'])
    ax2.text(0.6, 2.85, 'Cattle generates ~10% more reads', fontsize=10,
             color=PAL['raw'], fontweight='bold')
    ax2.text(0.6, 2.5, 'per gram of tissue', fontsize=10,
             color=PAL['raw'], fontweight='bold')

    # Biased reads
    ax2.text(0.6, 1.9, 'Raw sequencing read counts:', fontsize=10,
             color=PAL['text'], fontweight='bold')
    ax2.barh(1.4, 2.65, height=0.45, left=0.6,
             color=PAL['beef'], alpha=0.85)
    ax2.barh(1.4, 1.95, height=0.45, left=3.25,
             color=PAL['pork'], alpha=0.85)
    ax2.text(1.93, 1.4, 'Beef 56%', fontsize=9, ha='center', va='center',
             color='white', fontweight='bold')
    ax2.text(4.23, 1.4, 'Pork 44%', fontsize=9, ha='center', va='center',
             color='white', fontweight='bold')

    ax2.text(3.05, 0.75, 'Without correction:', fontsize=10,
             ha='center', color=PAL['raw'], fontweight='bold')
    ax2.text(3.05, 0.4, 'WRONG ANSWER', fontsize=11,
             ha='center', color=PAL['raw'], fontweight='bold')

    # ── Middle: EM ───────────────────────────────────────────────────
    box(ax2, (6.1, 1.3), 2.4, 3.6, PAL['em'], alpha=0.1,
        ec=PAL['em'], lw=2)
    ax2.text(7.3, 4.55, 'EM Algorithm', fontsize=11, fontweight='bold',
             ha='center', color=PAL['em'])
    ax2.text(7.3, 4.1, 'Corrects for:', fontsize=9.5, ha='center',
             color=PAL['em'])

    corrections = [
        ('1. mtDNA copy number\n    (species-specific)', 3.5),
        ('2. PCR amplification\n    efficiency (2-10x)', 2.65),
        ('3. DNA degradation\n    (processed food)', 1.8),
    ]
    for label, y in corrections:
        ax2.text(7.3, y, label, fontsize=8.5, ha='center',
                 color=PAL['text'], linespacing=1.2)

    arrow(ax2, 5.85, 2.9, 6.15, 2.9, PAL['em'], lw=2)
    arrow(ax2, 8.55, 2.9, 8.85, 2.9, PAL['em'], lw=2)

    # ── Right: Solution ──────────────────────────────────────────────
    box(ax2, (9.0, 0.2), 7.5, 5.5, PAL['panel_bg'], ec=PAL['panel_bdr'], lw=1.5)
    ax2.text(12.75, 5.4, 'The Solution: Iterative Correction',
             fontsize=12, fontweight='bold', ha='center', color=PAL['corrected'])

    iterations = [
        ('Start',   56, 44, 'Raw reads (biased)'),
        ('Iter 1',  53, 47, 'First correction'),
        ('Iter 5',  51, 49, 'Converging...'),
        ('Iter 20', 50, 50, 'Converged!'),
    ]

    for i, (label, beef, pork, note) in enumerate(iterations):
        y = 4.5 - i * 0.85
        ax2.text(9.3, y, label, fontsize=9, fontweight='bold',
                 color=PAL['grey'], va='center')
        bx = 10.4
        bw_total = 5.5
        w1 = beef / 100 * bw_total
        w2 = pork / 100 * bw_total
        ax2.barh(y, w1, height=0.4, left=bx,
                 color=PAL['beef'], alpha=0.85, zorder=3)
        ax2.barh(y, w2, height=0.4, left=bx + w1,
                 color=PAL['pork'], alpha=0.85, zorder=3)
        ax2.text(bx + w1/2, y, f'{beef}%', fontsize=8,
                 ha='center', va='center', color='white', fontweight='bold')
        ax2.text(bx + w1 + w2/2, y, f'{pork}%', fontsize=8,
                 ha='center', va='center', color='white', fontweight='bold')
        # Arrows between iterations
        if i < len(iterations) - 1:
            ax2.annotate('', xy=(10.1, y - 0.4),
                         xytext=(10.1, y - 0.25),
                         arrowprops=dict(arrowstyle='->', color=PAL['midgrey'],
                                         lw=0.8, mutation_scale=8))

    # Final result
    box(ax2, (9.3, 0.5), 6.9, 0.75, '#E8F5E9', ec=PAL['corrected'], lw=1.5)
    ax2.text(12.75, 0.88,
             'Corrected: 50% Beef + 50% Pork  (matches true weight)',
             fontsize=10.5, ha='center', va='center', color=PAL['corrected'],
             fontweight='bold')

    # Legend
    leg = [mpatches.Patch(facecolor=PAL['beef'], alpha=0.85,
                          label='Cattle (Bos taurus)'),
           mpatches.Patch(facecolor=PAL['pork'], alpha=0.85,
                          label='Pig (Sus scrofa)')]
    ax2.legend(handles=leg, loc='upper left', bbox_to_anchor=(0.3, 0.18),
               ncol=2, fontsize=9, frameon=False)

    # ═══════════════════════════════════════════════════════════════════
    # PANEL C: Output and verdict
    # ═══════════════════════════════════════════════════════════════════
    ax3 = fig.add_subplot(gs[2])
    ax3.set_xlim(0, 17)
    ax3.set_ylim(0, 6.5)
    ax3.axis('off')

    ax3.text(0.15, 6.3, 'C', fontsize=20, fontweight='bold', color=PAL['text'])
    ax3.text(0.7, 6.3, 'Output: species verdict with statistical confidence',
             fontsize=14, fontweight='bold', color=PAL['text'])

    # ── Left: EM output ──────────────────────────────────────────────
    box(ax3, (0.3, 0.2), 5.5, 5.5, PAL['panel_bg'], ec=PAL['panel_bdr'], lw=1.5)
    ax3.text(3.05, 5.4, 'EM Output (3-species example)',
             fontsize=11, fontweight='bold', ha='center', color=PAL['em'])

    species_out = [
        ('Cattle', 62.3, PAL['beef'], 'Halal', PAL['corrected']),
        ('Pig', 35.1, PAL['pork'], 'Haram', PAL['haram']),
        ('Chicken', 2.6, PAL['chicken'], 'Halal', PAL['corrected']),
    ]
    for i, (name, pct, color, status, sc) in enumerate(species_out):
        y = 4.5 - i * 1.3
        ax3.text(0.6, y + 0.3, name, fontsize=10, fontweight='bold',
                 color=PAL['text'])
        ax3.text(0.6, y - 0.05, f'{pct:.1f}%', fontsize=13, fontweight='bold',
                 color=color)

        # Bar
        bw = pct / 100 * 3.2
        ax3.barh(y - 0.45, bw, height=0.28, left=1.8,
                 color=color, alpha=0.8, zorder=3)

        # Status badge
        box(ax3, (3.8, y - 0.2), 1.4, 0.45, sc, alpha=0.12, ec=sc, lw=1.2)
        ax3.text(4.5, y + 0.03, status, fontsize=9, fontweight='bold',
                 ha='center', va='center', color=sc)

    # LRT
    ax3.text(3.05, 0.65, 'All species confirmed present', fontsize=9,
             ha='center', color=PAL['grey'], fontweight='bold')
    ax3.text(3.05, 0.35, '(likelihood ratio test: p < 0.001)',
             fontsize=8.5, ha='center', color=PAL['midgrey'], style='italic')

    # ── Arrow ────────────────────────────────────────────────────────
    arrow(ax3, 6.0, 3.0, 6.8, 3.0, PAL['output'], lw=2.5)

    # ── Right: Regulatory report ─────────────────────────────────────
    box(ax3, (7.0, 0.2), 9.5, 5.5, PAL['panel_bg'], ec=PAL['panel_bdr'], lw=1.5)
    ax3.text(11.75, 5.4, 'Regulatory Report',
             fontsize=12, fontweight='bold', ha='center', color=PAL['output'])

    # Haram alert
    box(ax3, (7.4, 4.0), 4.5, 1.0, PAL['haram_bg'], ec=PAL['haram'], lw=2.5)
    ax3.text(9.65, 4.65, 'HARAM SPECIES DETECTED', fontsize=12,
             fontweight='bold', ha='center', va='center', color=PAL['haram'])
    ax3.text(9.65, 4.25, 'Sus scrofa (pig): 35.1% +/- 2.3 pp',
             fontsize=9.5, ha='center', va='center', color=PAL['haram'])

    # Results table
    headers = ['Species', 'Weight', 'CI (95%)', 'LRT p', 'Status']
    rows = [
        ['Cattle',  '62.3%', '+/-3.1 pp', '<0.001', 'Halal'],
        ['Pig',     '35.1%', '+/-2.3 pp', '<0.001', 'Haram'],
        ['Chicken', '2.6%',  '+/-0.8 pp', '<0.001', 'Halal'],
    ]

    for j, h in enumerate(headers):
        ax3.text(7.6 + j * 1.6, 3.6, h, fontsize=8.5, fontweight='bold',
                 color=PAL['text'])
    ax3.plot([7.5, 15.6], [3.42, 3.42], color=PAL['ltgrey'], lw=1)

    for i, row in enumerate(rows):
        y = 3.1 - i * 0.45
        for j, val in enumerate(row):
            c = PAL['haram'] if i == 1 else PAL['grey']
            fw = 'bold' if i == 1 else 'normal'
            ax3.text(7.6 + j * 1.6, y, val, fontsize=8.5,
                     fontweight=fw, color=c)

    # Threshold note
    box(ax3, (12.2, 3.7), 4.0, 1.4, '#F0F4F8', ec=PAL['output'], lw=1.2)
    ax3.text(14.2, 4.75, 'EU Threshold', fontsize=9, fontweight='bold',
             ha='center', color=PAL['output'])
    ax3.text(14.2, 4.35, '1% (w/w)', fontsize=10, fontweight='bold',
             ha='center', color=PAL['output'])
    ax3.text(14.2, 3.95, 'Pig at 35.1% exceeds\nthreshold by 34x',
             fontsize=8.5, ha='center', color=PAL['grey'], linespacing=1.3)

    # Footer
    ax3.text(11.75, 0.65, 'Processing time: 0.64 s on commodity laptop',
             fontsize=9, ha='center', color=PAL['grey'], style='italic')
    ax3.text(11.75, 0.3, 'Reads classified: 1,247 / 1,500 (83%)',
             fontsize=9, ha='center', color=PAL['midgrey'], style='italic')

    # ── Save ─────────────────────────────────────────────────────────
    out = OUT / 'fig2_methodology_detail.png'
    fig.savefig(out, bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.close()
    print(f'Saved {out}  ({out.stat().st_size / 1024:.0f} KB)')


if __name__ == '__main__':
    fig2()
