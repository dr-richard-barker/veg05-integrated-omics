"""Derive human-readable biology labels for each WGCNA module from per-module GO
enrichment, write results/module_annotations.tsv, and regenerate Fig 4 with labelled axes.
Curated short labels are based on the top GO terms (Description) per module.
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

# Curated GO-based biology labels (from top enriched terms per module; see module_annotations.tsv)
LABELS = {
    'Leaf': {
        'black': 'Photosynthesis & ribosome', 'red': 'RNA processing',
        'turquoise': 'Chromatin / nucleosome', 'green': 'Protein ubiquitination',
        'pink': 'Cell wall / glucan', 'purple': 'Fe–S cluster / biosynthesis',
        'greenyellow': 'Chitinase / defense', 'magenta': 'Transcription / ethylene',
        'tan': 'Oxidative stress (peroxidase)', 'blue': 'Phosphate starvation / transport',
        'brown': 'Lignin / phenolic metabolism',
    },
    'Adv-Root': {
        'cyan': 'Photosynthesis (light harvesting)', 'blue': 'RNA processing',
        'turquoise': 'Oxidative stress (peroxidase)', 'salmon': 'Cell cycle / mitosis',
        'lightcyan': 'Protein phosphorylation', 'purple': 'Cellulose biosynthesis',
        'magenta': 'Terpenoid biosynthesis', 'pink': 'Transcription / ubiquitination',
        'midnightblue': 'Membrane transport', 'brown': 'Lipid / GTPase signalling',
        'greenyellow': 'Microtubule / nucleus', 'green': 'Cell wall biogenesis',
        'black': 'Nucleus',
    },
}
GO_FILES = {'Leaf': 'rnaseq/go_enrichment/go_enrichment_Leaf.tsv',
            'Adv-Root': 'rnaseq/go_enrichment/go_enrichment_AdvRoot.tsv'}
MT_FILES = {'Leaf': 'rnaseq/wgcna_Leaf/module_trait_correlations.tsv',
            'Adv-Root': 'rnaseq/wgcna_Adv-Root/module_trait_correlations.tsv'}
TRAITS = [('cor_flight', 'padj_flight', 'Flight'), ('cor_light', 'padj_light', 'Light'),
          ('cor_dysbiosis_16S', 'padj_dysbiosis_16S', 'Dysbiosis\n(16S)'),
          ('cor_dysbiosis_ITS', 'padj_dysbiosis_ITS', 'Dysbiosis\n(ITS)')]
CMAP = LinearSegmentedColormap.from_list('bwy', ['#3B6DC4', '#FFFFFF', '#E8E23A'])
NORM = TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)

def stars(p):
    return '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else ''

# ---- 1. annotations table (label + top-3 GO terms) ----
ann_rows = []
top_terms = {}
for tissue, gf in GO_FILES.items():
    go = pd.read_csv(os.path.join(HERE, gf), sep='\t'); go['Description'] = go['Description'].astype(str)
    for m, g in go.groupby('module'):
        g = g.sort_values('p.adjust')
        terms = [t for t in g['Description'].head(3).tolist() if t != 'nan']
        top_terms[(tissue, m)] = '; '.join(terms)
        ann_rows.append({'tissue': tissue, 'module': m,
                         'biology_label': LABELS[tissue].get(m, ''),
                         'n_GO_terms': len(g),
                         'top_GO_terms': '; '.join(terms),
                         'top_padj': f"{g['p.adjust'].min():.2e}"})
ann = pd.DataFrame(ann_rows)
os.makedirs(os.path.join(HERE, 'results'), exist_ok=True)
ann.to_csv(os.path.join(HERE, 'results', 'module_annotations.tsv'), sep='\t', index=False)
print('wrote module_annotations.tsv:', len(ann), 'modules with GO enrichment')

# ---- 2. regenerate Fig 4 with labelled y-axis ----
def build(tissue, out):
    d = pd.read_csv(os.path.join(HERE, MT_FILES[tissue]), sep='\t')
    d = d.sort_values('cor_flight', ascending=False).reset_index(drop=True)
    C = d[[t[0] for t in TRAITS]].values.astype(float)
    P = d[[t[1] for t in TRAITS]].values.astype(float)
    rows = []
    for c, n in zip(d['module_color'], d['n_genes']):
        lab = LABELS[tissue].get(c, '(no GO enrichment)')
        rows.append(f"{c} (n={n})\n{lab}")
    nr, nc = C.shape
    fig, ax = plt.subplots(figsize=(6.6, max(3.8, 0.52 * nr + 1.2)))
    im = ax.imshow(C, cmap=CMAP, norm=NORM, aspect='auto')
    for i in range(nr):
        for j in range(nc):
            v = C[i, j]
            ax.text(j, i, f"{v:.2f}\n{stars(P[i, j])}".rstrip(), ha='center', va='center',
                    fontsize=6.5, color='black' if abs(v) < 0.6 else ('white' if v < 0 else 'black'))
    ax.set_xticks(range(nc)); ax.set_xticklabels([t[2] for t in TRAITS], fontsize=8.5)
    ax.set_yticks(range(nr)); ax.set_yticklabels(rows, fontsize=6.8)
    ax.set_xticks(np.arange(-.5, nc, 1), minor=True); ax.set_yticks(np.arange(-.5, nr, 1), minor=True)
    ax.grid(which='minor', color='white', linewidth=1.2); ax.tick_params(which='minor', length=0)
    ax.set_ylabel('WGCNA module  (colour + dominant GO biology)', fontsize=9.5)
    ax.set_title(f'{tissue}: Module–Trait Correlations', fontsize=12, loc='left', pad=8)
    cb = fig.colorbar(im, ax=ax, shrink=0.5, aspect=16, pad=0.03, ticks=[-1, -.5, 0, .5, 1])
    cb.set_label('Spearman rho', fontsize=9); cb.ax.tick_params(labelsize=8)
    fig.text(0.01, 0.005, '* padj<0.05   ** <0.01   *** <0.001', fontsize=6.5, color='#555')
    plt.tight_layout()
    for ext in ('png', 'svg', 'pdf'):
        fig.savefig(os.path.join(F, f'{out}.{ext}'), bbox_inches='tight', dpi=300)
    plt.close(fig); print(f'{out}: {nr} modules relabelled')

build('Leaf', 'fig4_module_traits_Leaf')
build('Adv-Root', 'fig4_module_traits_AdvRoot')
