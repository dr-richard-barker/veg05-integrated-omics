"""Build Space_Tomato_Multi-Omics_Research.png, the summary infographic shown at
the top of README.md. Every number is taken from the manuscript and from
rnaseq/RERUN_2026-10-14.md and microbiome/RERUN_2026-10-15.md.

This replaces a NotebookLM-generated infographic whose claims (87% of root
genes up, 29 module-taxon links, adaptation to microgravity) the analyses do
not support.
"""
import os
import textwrap
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Ellipse, FancyBboxPatch, Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'Space_Tomato_Multi-Omics_Research.png')

W, H = 16.0, 8.93                      # inches; 2752 x 1536 px at 172 dpi
SURFACE, PANEL, RULE = '#fcfcfb', '#f3f2ee', '#dddcd6'
INK, INK2, MUTED = '#0b0b0b', '#52514e', '#7a7974'
BLUE, RED = '#2a78d6', '#e34948'       # blue-rich and red-rich light (validated pair)
LEAF, SOIL, ROOT = '#008300', '#8a7a66', '#b8a58c'
plt.rcParams.update({'font.family': 'DejaVu Sans'})

fig = plt.figure(figsize=(W, H), dpi=172)
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, W)
ax.set_ylim(0, H)
ax.axis('off')
fig.patch.set_facecolor(SURFACE)


def panel(x, y, w, h, fill=PANEL):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0,rounding_size=0.12',
                                fc=fill, ec=RULE, lw=1.2))


BOXES = []   # (text artist, (x0, y0, x1, y1) it must stay inside), checked after drawing


def chars_for(width_in, size):
    # characters per line for DejaVu Sans at this size (average glyph about 0.6 em)
    return max(10, int(width_in * 72 / (size * 0.6)))


def wrapped(x, y, text, width_in, size, color=INK2, weight='normal', spacing=1.35):
    lines = textwrap.wrap(text, chars_for(width_in, size))
    return ax.text(x, y, '\n'.join(lines), fontsize=size, color=color, fontweight=weight,
                   va='top', ha='left', linespacing=spacing)


def stat_tile(x, y, w, h, label, value, body):
    panel(x, y, w, h)
    inner = (x + 0.1, y + 0.05, x + w - 0.1, y + h - 0.05)
    BOXES.append((ax.text(x + 0.22, y + h - 0.22, label, fontsize=11, color=MUTED, va='top',
                          fontweight='bold'), inner))
    BOXES.append((ax.text(x + 0.22, y + h - 0.52, value, fontsize=30, color=INK, va='top',
                          fontweight='bold'), inner))
    BOXES.append((wrapped(x + 0.22, y + h - 1.12, body, w - 0.44, 10.5), inner))


# --- Title ---------------------------------------------------------------------
ax.text(0.4, H - 0.35, 'The Space Tomato: multi-omics of ISS VEG-05', fontsize=30,
        fontweight='bold', color=INK, va='top')
BOXES.append((wrapped(0.4, H - 0.98, 'Tomato cv. Red Robin grown on the ISS under red-rich or blue-rich '
                      'LED light and compared with ground controls: host transcriptome (OSD-767) and leaf and '
                      'root microbiome (OSD-766)', W - 0.8, 12.5), (0.3, 0, W - 0.3, H)))

# --- Light comparison: hero tiles with bars on one shared scale -----------------
TOP = H - 1.75
BAR_MAX = 4.0


def light_tile(x, colour, value, n, light):
    panel(x, TOP - 2.05, 4.9, 2.05)
    ax.add_patch(Rectangle((x, TOP - 2.05), 0.1, 2.05, fc=colour, ec='none'))
    ax.text(x + 0.35, TOP - 0.22, f'Leaf, {light} light', fontsize=11, color=MUTED,
            va='top', fontweight='bold')
    ax.text(x + 0.35, TOP - 0.55, value, fontsize=40, color=INK, va='top', fontweight='bold')
    ax.text(x + 0.35, TOP - 1.33, 'genes changed in flight vs ground', fontsize=11.5,
            color=INK2, va='top')
    ax.add_patch(FancyBboxPatch((x + 0.35, TOP - 1.85), BAR_MAX * n / 4716, 0.22,
                                boxstyle='round,pad=0,rounding_size=0.05', fc=colour, ec='none'))


light_tile(0.4, BLUE, '4,716', 4716, 'blue-rich')
light_tile(W - 0.4 - 4.9, RED, '523', 523, 'red-rich')
ax.text(W / 2, TOP - 2.22, 'Bars share one scale: about nine times more leaf genes respond to '
        'spaceflight under blue-rich light (separate models within each light)',
        fontsize=10.5, color=INK2, ha='center', va='top')

# --- Growth chamber with a tomato plant (decorative) ---------------------------
cx, cy0, cw, ch = W / 2, TOP - 2.05, 4.2, 2.05
left = cx - cw / 2
ax.add_patch(FancyBboxPatch((left, cy0), cw / 2, ch, boxstyle='round,pad=0,rounding_size=0.12',
                            fc=BLUE, alpha=0.10, ec='none'))
ax.add_patch(FancyBboxPatch((cx, cy0), cw / 2, ch, boxstyle='round,pad=0,rounding_size=0.12',
                            fc=RED, alpha=0.10, ec='none'))
ax.add_patch(FancyBboxPatch((left, cy0), cw, ch, boxstyle='round,pad=0,rounding_size=0.12',
                            fc='none', ec=RULE, lw=1.2))
ax.add_patch(Rectangle((left + 0.3, cy0 + ch - 0.16), cw / 2 - 0.45, 0.07, fc=BLUE, ec='none'))
ax.add_patch(Rectangle((cx + 0.15, cy0 + ch - 0.16), cw / 2 - 0.45, 0.07, fc=RED, ec='none'))
ax.text(left + 0.3, cy0 + ch - 0.25, 'blue-rich', fontsize=9.5, color=INK2, va='top')
ax.text(cx + cw / 2 - 0.3, cy0 + ch - 0.25, 'red-rich', fontsize=9.5, color=INK2, va='top', ha='right')
soil = cy0 + 0.45
ax.add_patch(Rectangle((left + 0.25, cy0 + 0.12), cw - 0.5, soil - cy0 - 0.12, fc=SOIL, alpha=0.35, ec='none'))
for dx, dy in [(-0.35, -0.25), (-0.12, -0.3), (0.1, -0.28), (0.32, -0.22), (0.0, -0.32)]:
    ax.plot([cx, cx + dx], [soil, soil + dy], color=ROOT, lw=1.6, solid_capstyle='round')
ax.plot([cx, cx], [soil, cy0 + ch - 0.45], color=LEAF, lw=3, solid_capstyle='round')
for i, (y, side) in enumerate([(0.75, -1), (0.95, 1), (1.15, -1), (1.32, 1), (1.47, -1)]):
    ax.add_patch(Ellipse((cx + side * 0.32, cy0 + y), 0.6, 0.2, angle=side * 25, fc=LEAF, alpha=0.75, ec='none'))
for x, y in [(cx - 0.22, cy0 + 1.0), (cx + 0.25, cy0 + 1.22), (cx - 0.18, cy0 + 1.36)]:
    ax.add_patch(Circle((x, y), 0.07, fc=RED, ec='none'))

# --- Second row: main effects and interaction -----------------------------------
ROW2 = TOP - 2.55
stat_tile(0.4, ROW2 - 1.6, 4.9, 1.6, 'Flight main effect (averaged over light)', '2,683 · 1,501',
          'genes changed in leaves · in adventitious roots, where 83% are higher in flight')
stat_tile(W / 2 - 2.1, ROW2 - 1.6, 4.2, 1.6, 'Flight × light interaction', '3,189 · 6',
          'genes in leaves · in roots: light shapes the leaf response far more')
stat_tile(W - 0.4 - 4.9, ROW2 - 1.6, 4.9, 1.6, 'MOFA+ shared factor (18 leaf samples)', '44%',
          'of transcriptome variance; tracks flight (ρ=−0.68) and light, which are partly confounded here')

# --- Third row: microbiome and host-microbe links -------------------------------
ROW3 = ROW2 - 1.8
tw, gap = (W - 0.8 - 3 * 0.25) / 4, 0.25
tiles = [
    ('Microbiome composition', '28% · 36%',
     'of leaf bacterial and fungal community variation explained by spaceflight (PERMANOVA, p = 0.001 and 0.007)'),
    ('Leaf dysbiosis, flight vs ground', '2.1× · 3.7×',
     'distance from ground communities for bacteria and fungi (ground = 1)'),
    ('Root module and fungi', 'r = −0.85',
     'a 169-gene root module tracks fungal dysbiosis (12 plants) but not flight (r=−0.38); dysbiosis is itself higher in flight'),
    ('Module–taxon links', '31',
     'significant correlations; the flight-associated leaf module rises with Methylobacterium-Methylorubrum, Burkholderia and Azospirillum'),
]
for i, (label, value, body) in enumerate(tiles):
    stat_tile(0.4 + i * (tw + gap), ROW3 - 2.1, tw, 2.1, label, value, body)

# --- Footer --------------------------------------------------------------------
ax.text(0.4, 0.42, 'Correlational study with small samples (21 leaf and 15 root RNA-seq samples). '
        'Data: NASA OSDR OSD-766 and OSD-767. Code: github.com/dr-richard-barker/veg05-integrated-omics '
        '(built by build_infographic.py).', fontsize=10, color=MUTED, va='center')

# Every text block must sit inside its tile (or the page): measure, don't guess
fig.canvas.draw()
to_data = ax.transData.inverted()
bad = []
for artist, (x0, y0, x1, y1) in BOXES:
    (bx0, by0), (bx1, by1) = to_data.transform(artist.get_window_extent())
    if bx0 < x0 - 1e-6 or by0 < y0 - 1e-6 or bx1 > x1 + 1e-6 or by1 > y1 + 1e-6:
        bad.append(artist.get_text()[:50])
if bad:
    raise SystemExit('text outside its box: ' + ' | '.join(bad))
fig.savefig(OUT, dpi=172, facecolor=SURFACE)
print(f'wrote {os.path.basename(OUT)}')
