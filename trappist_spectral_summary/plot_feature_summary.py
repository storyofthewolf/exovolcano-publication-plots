"""
plot_feature_summary.py -- the spectral summary on the day-0-referenced metric.

Supersedes plot_amplitude_summary.py, whose windowed continuum fits were
audited on 2026-09-06 and found wrong in six places
(exovolcano-spectra/notes/BAND_WINDOW_AUDIT.md). Nothing is fitted here: each
case's day-0 spectrum is the pre-eruption state of that same atmosphere and is
identical across all six of its eruptions, so it IS the control, and
differencing against it removes that atmosphere's own bands exactly. See the
config header.

TWO PANELS, because the two views answer different questions:

  LEFT   signal     day-0-referenced change, ppm. What the eruption ADDED.
  RIGHT  amplitude  peak-to-trough at the most muted epoch, as a fraction of
                    pre-eruption. What the aerosol FLATTENED.

CO2 4.3 is why both are needed. Saturated on every atmosphere, so the eruption
cannot deepen it and its signal is ~0 -- yet its amplitude falls to 0.12 of
pre-eruption on hab1 tambora_100x. Either panel alone calls that band inert.

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
from matplotlib.patches import Patch

plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 8,
    'axes.labelsize': 9,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'legend.fontsize': 7.5,
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
ALAB = cfg['atm_labels']
RUNGS = cfg['rungs']
BANDS = cfg['bands']


def case_name(ladder, rung):
    return ladder if rung == 1 else f"{ladder}_{rung}x"


# Rows: blocked by atmosphere, each holding its six eruptions. Reading down a
# block is one planet's full response; the Tambora|Hunga seam inside each block
# is ruled lightly, since the two ladders scale different quantities.
order, block_edges, block_starts = [], [], []
for atm in ATMS:
    block_starts.append((len(order), ALAB[atm]))
    for ladder in ('tambora', 'hunga'):
        for rung in RUNGS:
            order.append((case_name(ladder, rung), atm))
    block_edges.append(len(order))

nb = len(BANDS)
SIG = np.full((len(order), nb), np.nan)
AMP = np.full((len(order), nb), np.nan)
CAU = np.zeros((len(order), nb), dtype=bool)
for i, (case, atm) in enumerate(order):
    s = df[(df.case == case) & (df.atm == atm)]
    if s.empty:
        continue
    for j, b in enumerate(BANDS):
        SIG[i, j] = s[f"{b['band']} sig"].values[0]
        AMP[i, j] = s[f"{b['band']} amp_ratio"].values[0]
        CAU[i, j] = bool(s[f"{b['band']} caution"].values[0])

fig, axes = plt.subplots(1, 2, figsize=cfg['figsize'], sharey=True)

panels = [
    (axes[0], SIG, cfg['signal_cmap'], 0.0, cfg['signal_vmax'],
     cfg['signal_label'], 'signal'),
    (axes[1], AMP, cfg['amplitude_cmap'], 0.0, 1.0,
     cfg['amplitude_label'], 'amplitude'),
]

for ax, M, cmap, vmin, vmax, label, kind in panels:
    im = ax.imshow(M, cmap=cmap, vmin=vmin, vmax=vmax,
                   aspect='auto', interpolation='nearest')

    # Hatch the cells where another absorber overlaps this band on this
    # atmosphere. The value is still measured and still drawn; the hatch says
    # only that it cannot be attributed to one species.
    for i, j in np.argwhere(CAU):
        ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False,
                                   hatch=cfg['caution_hatch'], lw=0.0,
                                   edgecolor='0.25', alpha=0.55, zorder=4))

    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            if np.isnan(M[i, j]):
                continue
            r, g, b_, _ = im.cmap(im.norm(M[i, j]))
            shade = ('white' if (0.2126 * r + 0.7152 * g + 0.0722 * b_) < 0.5
                     else 'black')
            if kind == 'signal':
                v = M[i, j]
                # Parenthesised negatives, the accounting convention, so the
                # sign survives on a sequential colour scale.
                txt = (f"({abs(v):.0f})" if v < -0.5
                       else (f"{v:.0f}" if abs(v) >= 1 else f"{v:.1f}"))
            else:
                txt = f"{M[i, j]:.2f}"
            ax.text(j, i, txt, ha='center', va='center', fontsize=6.4,
                    color=shade, zorder=5)

    ax.set_xticks(range(nb))
    ax.set_xticklabels([b['label'] for b in BANDS], rotation=35, ha='right')
    for k in block_edges[:-1]:
        ax.axhline(k - 0.5, color='black', lw=1.6)
    for start, _ in block_starts:
        ax.axhline(start + len(RUNGS) - 0.5, color='0.45', lw=0.7)

    # Colorbar padded well clear of the rotated band labels; at the default
    # pad it sat on top of them.
    cb = fig.colorbar(im, ax=ax, pad=0.16, fraction=0.045,
                      orientation='horizontal', location='bottom')
    cb.set_label(label, fontsize=8)

axes[0].set_yticks(range(len(order)))
axes[0].set_yticklabels([c for c, _ in order], fontsize=7)
# Atmosphere name down the left of each block. Pushed further out than the
# case labels and given the short name only: the full "aquaplanet, CO2" gloss
# collided with its neighbours at this row height.
for start, lab in block_starts:
    axes[0].text(-0.245, start + (2 * len(RUNGS) - 1) / 2.0,
                 lab.split(':')[0], transform=axes[0].get_yaxis_transform(),
                 rotation=90, ha='center', va='center', fontsize=9,
                 fontweight='semibold')

axes[0].set_title('What the eruption ADDED\n'
                  'day-0 referenced signal', fontsize=9, pad=8)
axes[1].set_title('What the aerosol FLATTENED\n'
                  'amplitude relative to pre-eruption', fontsize=9, pad=8)

fig.legend(handles=[Patch(facecolor='white', edgecolor='0.25',
                          hatch=cfg['caution_hatch'],
                          label='another absorber overlaps this band on this '
                                'atmosphere')],
           loc='lower center', frameon=False, bbox_to_anchor=(0.5, -0.06))

fig.subplots_adjust(wspace=0.05, bottom=0.20, top=0.90, left=0.12)
stem = os.path.join(here, cfg['outfile'])
for ext in ('pdf', 'eps'):
    fig.savefig(f'{stem}.{ext}', bbox_inches='tight')
    print(f'wrote {cfg["outfile"]}.{ext}')
