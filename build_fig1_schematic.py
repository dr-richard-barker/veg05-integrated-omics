"""Figure 1 — VEG-05 multi-omics study-design / overview schematic."""
import os
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['DejaVu Sans'],
                     'pdf.fonttype': 42, 'svg.fonttype': 'none'})
HERE = os.path.dirname(os.path.abspath(__file__)); F = os.path.join(HERE, 'figures')
C_EXP, C_RNA, C_MIC, C_INT, C_OUT, C_DK = '#4C6E91', '#C58A2E', '#4A8C6F', '#6E5A8C', '#B0472F', '#222222'

fig, ax = plt.subplots(figsize=(11, 5.6)); ax.set_xlim(0, 116); ax.set_ylim(0, 56); ax.axis('off')

def box(x, y, w, h, fc, ec, title, body, tcol='white', fs=8.4):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.4,rounding_size=1.6',
                                facecolor=fc, edgecolor=ec, linewidth=1.4))
    ax.text(x + w/2, y + h - 3.0, title, ha='center', va='top', fontsize=fs+1.1, fontweight='bold', color=tcol)
    ax.text(x + w/2, y + h - 8.0, body, ha='center', va='top', fontsize=fs, color=tcol, linespacing=1.4)

def arrow(x1, y1, x2, y2, col=C_DK, lw=2.0):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle='-|>', mutation_scale=17, lw=lw, color=col))

ax.text(58, 54, 'ISS VEG-05 multi-omics integration: host transcriptome × microbiome',
        ha='center', fontsize=13, fontweight='bold', color=C_DK)

# experiment
box(1.5, 20, 26, 26, C_EXP, '#33506e', 'VEG-05 experiment',
    'Tomato cv. Red Robin\naboard the ISS\n\nFlight vs Ground\nRed vs Blue LED\n\nLeaf +\nAdventitious root')
# two data streams
box(33, 33, 27, 13, C_RNA, '#8f6216', 'Host transcriptome', 'RNA-seq · OSD-767\nLeaf & Adv-Root')
box(33, 15, 27, 13, C_MIC, '#316049', 'Microbiome', '16S + ITS amplicons · OSD-766\n348 bacterial / 77 fungal ASVs')
arrow(27.5, 39, 32.5, 40); arrow(27.5, 27, 32.5, 21)
# integration
box(65, 20, 27, 26, C_INT, '#4a3d61', 'Integration',
    'DESeq2 (DEGs)\nWGCNA co-expression\nDysbiosis index\nMOFA+ (3 views)\nModule–taxon networks\nGO annotation')
arrow(60.5, 39, 64.5, 36); arrow(60.5, 21, 64.5, 26)
# outcome
box(97, 20, 17.5, 26, C_OUT, '#7d2f1f', 'Findings',
    'Spaceflight\ndysbiosis\n\nLight-\ndependent\nresponse\n\nMicrobe-\ndriven host\ndefense', fs=8.2)
arrow(92.5, 33, 96.5, 33)

for ext in ('png', 'svg', 'pdf'):
    fig.savefig(os.path.join(F, f'fig1_study_design.{ext}'), bbox_inches='tight', dpi=300)
plt.close(fig)
print('fig1_study_design written')
