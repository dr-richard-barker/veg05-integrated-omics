"""Rebuild Fig 6 (module-taxon) as a clean bubble grid with GO-based module biology labels.
Shows the significant (padj<0.05) module-taxon correlations for Leaf 16S (the 28 of 29 hits);
the single Adv-Root ITS hit (blue x Ralstonia) is noted in the caption/footnote.
"""
import os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['DejaVu Sans'],
                     'pdf.fonttype': 42, 'svg.fonttype': 'none'})
HERE = os.path.dirname(os.path.abspath(__file__)); F = os.path.join(HERE, 'figures')
LAB = {'turquoise': 'Chromatin/nucleosome', 'blue': 'Phosphate starvation/transport',
       'red': 'RNA processing', 'brown': 'Lignin/phenolics', 'purple': 'Fe–S/biosynthesis',
       'black': 'Photosynthesis/ribosome', 'yellow': '(no GO enrichment)'}
SHORT = {'Burkholderia-Caballeronia-Paraburkholderia': 'Burkholderia\n(complex)',
         'Methylobacterium-Methylorubrum': 'Methylobacterium'}
CMAP = LinearSegmentedColormap.from_list('bwy', ['#3B6DC4', '#FFFFFF', '#E8E23A'])
NORM = TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)

d = pd.read_csv(os.path.join(HERE, 'integration/networks/module_taxon_correlations.tsv'), sep='\t')
sig = d[(d['padj'] < 0.05) & (d['tissue'] == 'Leaf') & (d['amplicon'] == '16S')].copy()
# one point per (module, genus): keep most significant
sig = sig.sort_values('padj').drop_duplicates(['module', 'Genus'])
sig['glab'] = sig['Genus'].map(lambda g: SHORT.get(g, g))
mods = sig.groupby('module')['rho'].apply(lambda s: s.abs().max()).sort_values(ascending=False).index.tolist()
gens = sorted(sig['glab'].unique())
mi = {m: i for i, m in enumerate(mods)}; gi = {g: i for i, g in enumerate(gens)}

fig, ax = plt.subplots(figsize=(8.4, 0.62 * len(mods) + 2.2))
sc = ax.scatter([gi[g] for g in sig['glab']], [mi[m] for m in sig['module']],
                c=sig['rho'], cmap=CMAP, norm=NORM, s=-np.log10(sig['padj']) * 90,
                edgecolor='k', linewidth=0.5, zorder=3)
ax.set_xticks(range(len(gens))); ax.set_xticklabels(gens, rotation=40, ha='right', fontsize=8)
ax.set_yticks(range(len(mods)))
ax.set_yticklabels([f"{m}\n{LAB.get(m, '')}" for m in mods], fontsize=8)
ax.set_ylim(-0.6, len(mods) - 0.4); ax.set_xlim(-0.6, len(gens) - 0.4)
ax.grid(True, color='#eee', zorder=0)
ax.set_xlabel('Bacterial genus (16S)', fontsize=10)
ax.set_ylabel('WGCNA module  (colour + GO biology)', fontsize=10)
ax.set_title('Leaf module–taxon correlations (significant, BH padj<0.05)', fontsize=12, loc='left', pad=8)
cb = fig.colorbar(sc, ax=ax, shrink=0.55, aspect=16, pad=0.02, ticks=[-1, -.5, 0, .5, 1])
cb.set_label('Spearman rho', fontsize=9)
# size legend
for p, lab in [(0.05, 'padj .05'), (0.005, 'padj .005')]:
    ax.scatter([], [], s=-np.log10(p) * 90, c='#bbb', edgecolor='k', linewidth=0.5, label=lab)
ax.legend(title='dot size', loc='upper left', bbox_to_anchor=(1.02, 1), fontsize=7.5, title_fontsize=8, frameon=False)
fig.text(0.01, 0.005, 'One additional significant hit not shown: Adv-Root ITS blue × Ralstonia (rho=-0.98).',
         fontsize=6.5, color='#555')
plt.tight_layout()
for ext in ('png', 'svg', 'pdf'):
    fig.savefig(os.path.join(F, f'fig6_module_taxon_network.{ext}'), bbox_inches='tight', dpi=300)
plt.close(fig)
print('fig6 rebuilt:', len(sig), 'Leaf-16S pairs,', len(mods), 'modules x', len(gens), 'genera')
