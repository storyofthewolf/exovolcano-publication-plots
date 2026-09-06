"""
plot_feature_summary.py -- one heatmap: signed change in feature amplitude.

Each cell is the change in a band's peak-to-trough transit depth, in ppm,
against that case's own pre-eruption spectrum:

    amp(t) = max - min of transit depth within the band at epoch t
    value  = amp(t*) - amp(0)    at whichever t* gives the larger excursion

Positive means the eruption DEEPENED the feature; negative that the sulfate
aerosol's raised continuum FLATTENED it. Both are real outcomes, and the
diverging scale keeps them visibly distinct kinds of event.

Why amplitude and not a difference spectrum: a raised continuum lifts every
point in a band, so depth(t) - depth(0) returns its LARGEST values exactly
where the aerosol is erasing the feature. Amplitude asks how far a band's core
sits below its own wings, which is what detectability follows. This is the
continuum-absorber behaviour in Fauchez et al. (2019), Fig. 9.

Why ppm and not a ratio: a ratio normalises each band to itself, so losing
half of a 4 ppm feature and half of a 45 ppm feature look identical. ppm is
the unit an instrument measures.

Bands and windows come from exovolcano-spectra's feature_table_clean.py, which
measures species footprints from the spectra rather than assuming line
positions. Five bands and not eight: gas overlap makes the rest
unattributable, the same reason Fauchez et al. report a handful of named lines
rather than every band.

    python plot_feature_summary.py
"""

import os
import argparse

import numpy as np
import pandas as pd
import yaml
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 8,
    'axes.labelsize': 9,
    'xtick.labelsize': 8.5,
    'ytick.labelsize': 8,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})

here = os.path.dirname(os.path.abspath(__file__))
ap = argparse.ArgumentParser()
ap.add_argument('--config', default='config_feature_summary.yaml')
args = ap.parse_args()
with open(os.path.join(here, args.config)) as f:
    cfg = yaml.safe_load(f)

df = pd.read_csv(cfg['table_csv'])
ATMS = cfg['atmospheres']
RUNGS = cfg['rungs']
BANDS = cfg['bands']
VCOL = cfg['value_column']
SCOL = cfg.get('secondary_column')


def case_name(ladder, rung):
    return ladder if rung == 1 else f"{ladder}_{rung}x"


# Rows blocked by atmosphere, each holding its six eruptions, so reading down a
# block is one planet's full response. The Tambora|Hunga seam inside each block
# is ruled lightly: the two ladders scale different quantities and are not
# comparable rung for rung.
order, block_edges, block_starts = [], [], []
for atm in ATMS:
    block_starts.append((len(order), cfg['atm_labels'][atm]))
    for ladder in ('tambora', 'hunga'):
        for rung in RUNGS:
            order.append((case_name(ladder, rung), atm))
    block_edges.append(len(order))

M = np.full((len(order), len(BANDS)), np.nan)
S = np.full((len(order), len(BANDS)), np.nan)
for i, (case, atm) in enumerate(order):
    s = df[(df.case == case) & (df.atm == atm)]
    if s.empty:
        continue
    for j, b in enumerate(BANDS):
        M[i, j] = s[f"{b['band']} {VCOL}"].values[0]
        if SCOL:
            S[i, j] = s[f"{b['band']} {SCOL}"].values[0]

v = cfg['vlim']
fig, ax = plt.subplots(figsize=cfg['figsize'])
im = ax.imshow(M, cmap=cfg['cmap'], vmin=-v, vmax=v,
               aspect='auto', interpolation='nearest')

if cfg.get('annotate', True):
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            if np.isnan(M[i, j]):
                continue
            # Text colour from the cell's own luminance, so it stays legible
            # at both saturated ends of a diverging map.
            r, g, b_, _ = im.cmap(im.norm(M[i, j]))
            shade = ('white' if (0.2126 * r + 0.7152 * g + 0.0722 * b_) < 0.5
                     else 'black')
            val = M[i, j]
            txt = f"{val:+.0f}" if abs(val) >= 1 else f"{val:+.1f}"
            if SCOL and not np.isnan(S[i, j]):
                # Two numbers per cell: the amplitude change sets the colour,
                # the difference-spectrum peak sits under it in parentheses so
                # the figure can be reconciled against A5-A10 by eye.
                sv = S[i, j]
                stxt = f"({sv:+.0f})" if abs(sv) >= 1 else f"({sv:+.1f})"
                ax.text(j, i - 0.16, txt, ha='center', va='center',
                        fontsize=7, color=shade, zorder=3)
                # No alpha: the PostScript backend drops transparency, so an
                # alpha'd label renders full-strength in the EPS and faded in
                # the PDF, and the manuscript would disagree with the
                # reference figure built here.
                ax.text(j, i + 0.20, stxt, ha='center', va='center',
                        fontsize=5.8, color=shade, zorder=3)
            else:
                ax.text(j, i, txt, ha='center', va='center', fontsize=7,
                        color=shade, zorder=3)

ax.set_xticks(range(len(BANDS)))
ax.set_xticklabels([b['label'] for b in BANDS])
ax.set_yticks(range(len(order)))
ax.set_yticklabels([c for c, _ in order], fontsize=7.5)

for k in block_edges[:-1]:
    ax.axhline(k - 0.5, color='black', lw=1.6)
for start, _ in block_starts:
    ax.axhline(start + len(RUNGS) - 0.5, color='0.5', lw=0.6)
for start, lab in block_starts:
    ax.text(-0.235, start + (2 * len(RUNGS) - 1) / 2.0, lab,
            transform=ax.get_yaxis_transform(), rotation=90,
            ha='center', va='center', fontsize=9.5, fontweight='semibold')

cb = fig.colorbar(im, ax=ax, pad=0.025, fraction=0.048)
cb.set_label('change in feature amplitude [ppm]')
# Name the two directions on the bar itself, so the sign convention does not
# have to be carried in from the caption.
cb.ax.text(0.5, 0.985, 'deepened', transform=cb.ax.transAxes, ha='center',
           va='top', fontsize=7, rotation=90)
cb.ax.text(0.5, 0.015, 'flattened', transform=cb.ax.transAxes, ha='center',
           va='bottom', fontsize=7, rotation=90)

_sub = ('change in band amplitude, peak-to-trough'
        if not SCOL else
        'top: change in band amplitude   (bottom): '
        + cfg.get('secondary_label', SCOL))
ax.set_title('Spectral feature response to the eruption\n' + _sub
             + '  [ppm, R = 250]', fontsize=9.5, pad=10)

fig.subplots_adjust(left=0.16, right=0.99, top=0.93, bottom=0.05)
stem = os.path.join(here, cfg['outfile'])
for ext in ('pdf', 'eps'):
    fig.savefig(f'{stem}.{ext}', bbox_inches='tight')
    print(f'wrote {cfg["outfile"]}.{ext}')
