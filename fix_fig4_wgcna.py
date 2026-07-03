"""Rebuild Fig 4 WGCNA module-trait heatmaps (Leaf + Adv-Root).
The original R figures rendered empty (tiles/labels missing). Regenerated here in
matplotlib from rnaseq/wgcna_*/module_trait_correlations.tsv, matching the study's style.
"""
import os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['DejaVu Sans'],
                     'pdf.fonttype': 42, 'svg.fonttype': 'none'})
HERE = os.path.dirname(os.path.abspath(__file__))
F = os.path.join(HERE, 'figures')
CMAP = LinearSegmentedColormap.from_list('bwy', ['#3B6DC4', '#FFFFFF', '#E8E23A'])
NORM = TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)
TRAITS = [('cor_flight', 'padj_flight', 'Flight'),
          ('cor_light', 'padj_light', 'Light'),
          ('cor_dysbiosis_16S', 'padj_dysbiosis_16S', 'Dysbiosis\n(16S)'),
          ('cor_dysbiosis_ITS', 'padj_dysbiosis_ITS', 'Dysbiosis\n(ITS)')]

def stars(p):
    return '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else ''

def build(tissue, tsv, out):
    d = pd.read_csv(tsv, sep='\t')
    d = d.sort_values('cor_flight', ascending=False).reset_index(drop=True)
    C = d[[t[0] for t in TRAITS]].values.astype(float)
    P = d[[t[1] for t in TRAITS]].values.astype(float)
    rows = [f"{c}  (n={n})" for c, n in zip(d['module_color'], d['n_genes'])]
    nr, nc = C.shape
    fig, ax = plt.subplots(figsize=(5.6, max(3.5, 0.42 * nr + 1.2)))
    im = ax.imshow(C, cmap=CMAP, norm=NORM, aspect='auto')
    for i in range(nr):
        for j in range(nc):
            val = C[i, j]
            txt = f"{val:.2f}\n{stars(P[i, j])}".rstrip()
            ax.text(j, i, txt, ha='center', va='center', fontsize=6.5,
                    color='black' if abs(val) < 0.6 else 'white' if val < 0 else 'black')
    ax.set_xticks(range(nc)); ax.set_xticklabels([t[2] for t in TRAITS], fontsize=8.5)
    ax.set_yticks(range(nr)); ax.set_yticklabels(rows, fontsize=7.5)
    ax.set_xticks(np.arange(-.5, nc, 1), minor=True)
    ax.set_yticks(np.arange(-.5, nr, 1), minor=True)
    ax.grid(which='minor', color='white', linewidth=1.2); ax.tick_params(which='minor', length=0)
    ax.set_ylabel('WGCNA Module', fontsize=10)
    ax.set_title(f'{tissue}: Module–Trait Correlations', fontsize=12, loc='left', pad=8)
    cb = fig.colorbar(im, ax=ax, shrink=0.5, aspect=16, pad=0.03, ticks=[-1, -0.5, 0, 0.5, 1])
    cb.set_label('Spearman rho', fontsize=9); cb.ax.tick_params(labelsize=8)
    fig.text(0.01, 0.01, '* padj<0.05   ** <0.01   *** <0.001', fontsize=6.5, color='#555')
    plt.tight_layout()
    for ext in ('png', 'svg', 'pdf'):
        fig.savefig(os.path.join(F, f'{out}.{ext}'), bbox_inches='tight', dpi=300)
    plt.close(fig)
    print(f'{out}: {nr} modules x {nc} traits')

build('Leaf', os.path.join(HERE, 'rnaseq/wgcna_Leaf/module_trait_correlations.tsv'), 'fig4_module_traits_Leaf')
build('Adv-Root', os.path.join(HERE, 'rnaseq/wgcna_Adv-Root/module_trait_correlations.tsv'), 'fig4_module_traits_AdvRoot')
