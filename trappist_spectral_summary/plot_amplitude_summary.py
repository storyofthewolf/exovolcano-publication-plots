"""
plot_amplitude_summary.py -- figures A and B rebuilt on the AMPLITUDE metric.

The earlier pair (plot_trappist_spectral_summary.py) plots
max |depth(lam,t) - depth(lam,0)| within a band. That is the right question
for a gas that adds absorption and the wrong one whenever aerosol is present:
sulfate raises the continuum bodily, so every point in an H2O band moves up
and the difference metric returns its LARGEST value exactly when the aerosol
is flattening the band into the continuum. On hab1 tambora_100x at 1.4 um the
difference peaks at +64 ppm -- the biggest H2O entry in that table -- while
the band's peak-to-trough amplitude falls to 4.9% of pre-eruption. The feature
is all but erased and the old number called it the strongest water signal in
the suite.

These two plot feature AMPLITUDE instead, as a ratio to the same band's
pre-eruption amplitude, which is what detectability actually follows. See
Fauchez et al. (2019), arXiv:1911.08596, Fig. 9.

    python plot_amplitude_summary.py           # both
    python plot_amplitude_summary.py --only a

A RATIO IS MULTIPLICATIVE ABOUT 1. Both figures work in log space for that
reason -- see the config header. A linear axis would squeeze every degree of
suppression into [0,1] and hand all of [1,inf) to enhancement, making a factor
of 2 down look smaller than a factor of 2 up.
"""

import os
import argparse

import numpy as np
import pandas as pd
import yaml
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import LogNorm

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
ap.add_argument('--config', default='config_amplitude_summary.yaml')
ap.add_argument('--only', default=None, choices=('a', 'b'))
args = ap.parse_args()
with open(os.path.join(here, args.config)) as f:
    cfg = yaml.safe_load(f)

df = pd.read_csv(cfg['table_csv'])
ATMS = cfg['atmospheres']
ACOL, ALAB, AMRK = cfg['atm_colors'], cfg['atm_labels'], cfg['atm_markers']
RUNGS = cfg['rungs']


def case_name(ladder, rung):
    return ladder if rung == 1 else f"{ladder}_{rung}x"


def save(fig, stem):
    path = os.path.join(here, stem)
    for ext in ('pdf', 'eps'):
        fig.savefig(f'{path}.{ext}', bbox_inches='tight')
        print(f'  wrote {stem}.{ext}')
    plt.close(fig)


# ---------------------------------------------------------------------------
# A -- ladder on the amplitude ratio
# ---------------------------------------------------------------------------
def figure_a():
    print('A: amplitude ladder')
    rows, cols = cfg['ladder_rows'], cfg['ladder_cols']
    fig, axes = plt.subplots(len(rows), len(cols),
                             figsize=cfg['ladder_figsize'], sharex=True,
                             sharey=cfg.get('ladder_logy', True))
    x = np.arange(len(RUNGS))

    for r, rowspec in enumerate(rows):
        col_key = f"{rowspec['band']} {rowspec['series']}"
        for c, colspec in enumerate(cols):
            ax = axes[r, c]
            # Unity is the meaningful reference: at 1.0 the eruption left the
            # feature's amplitude untouched. Everything below is a feature
            # being veiled, everything above one being deepened.
            ax.axhline(1.0, color='0.45', lw=0.9, ls='--', zorder=2)

            for atm in ATMS:
                ys = []
                for rung in RUNGS:
                    s = df[(df.case == case_name(colspec['ladder'], rung))
                           & (df.atm == atm)]
                    ys.append(s[col_key].values[0] if not s.empty else np.nan)
                ys = np.array(ys, dtype=float)
                ax.plot(x, ys, marker=AMRK[atm], color=ACOL[atm], lw=1.3,
                        ms=4.5, zorder=3,
                        label=ALAB[atm] if (r == 0 and c == 0) else None)

                # A band that was flat before the eruption has no amplitude to
                # take a ratio against (ben1's 6.3 um on a dry N2 planet), so
                # amplitude_table_24.py leaves it NaN. Say so on the panel
                # rather than letting a gap read as missing data.
                if np.all(np.isnan(ys)):
                    ax.text(0.5, 0.06,
                            f"{atm}: no pre-eruption feature",
                            transform=ax.transAxes, ha='center', va='bottom',
                            fontsize=6.5, color=ACOL[atm], style='italic')

            if cfg.get('ladder_logy', True):
                ax.set_yscale('log')
                ax.set_ylim(*cfg['ladder_ylim'])
            ax.set_xticks(x)
            ax.set_xticklabels(cfg['rung_labels'])
            ax.set_xlim(-0.25, len(RUNGS) - 0.75)
            ax.grid(axis='y', color='0.9', lw=0.5, zorder=0)
            ax.set_axisbelow(True)

            if r == 0:
                ax.set_title(colspec['label'], fontsize=8.5, pad=7)
            if c == 0:
                ax.set_ylabel(rowspec['label'] + "\n"
                              r"amplitude / pre-eruption")
            else:
                ax.tick_params(labelleft=False)
            if r == len(rows) - 1:
                ax.set_xlabel('magnitude rung')

    # Shade the suppressed half once per panel, so "below the line is worse"
    # is legible without reading the axis.
    for axrow in axes:
        for ax in axrow:
            ax.axhspan(cfg['ladder_ylim'][0], 1.0, color='0.93', zorder=0)

    axes[0, 0].legend(loc='lower left', frameon=False, handlelength=1.8)
    fig.text(0.5, 0.955,
             'Feature amplitude relative to pre-eruption '
             '(below 1 = muted by aerosol continuum)',
             ha='center', fontsize=9)
    fig.subplots_adjust(hspace=0.16, wspace=0.06, top=0.925)
    save(fig, cfg['outfile_a'])


# ---------------------------------------------------------------------------
# B -- heatmap on the amplitude ratio
# ---------------------------------------------------------------------------
def figure_b():
    print('B: amplitude heatmap')
    bands = cfg['heatmap_bands']
    series = cfg.get('heatmap_series', 'minratio')

    order, block_edges, block_labels = [], [], []
    if cfg.get('heatmap_group_by', 'atmosphere') == 'atmosphere':
        for atm in ATMS:
            block_labels.append((len(order), ALAB[atm]))
            for ladder in ('tambora', 'hunga'):
                for rung in RUNGS:
                    order.append((case_name(ladder, rung), atm))
            block_edges.append(len(order))
    else:
        for ladder in ('tambora', 'hunga'):
            for rung in RUNGS:
                block_labels.append((len(order), case_name(ladder, rung)))
                for atm in ATMS:
                    order.append((case_name(ladder, rung), atm))
                block_edges.append(len(order))

    M = np.full((len(order), len(bands)), np.nan)
    for i, (case, atm) in enumerate(order):
        s = df[(df.case == case) & (df.atm == atm)]
        if s.empty:
            continue
        for j, b in enumerate(bands):
            M[i, j] = s[f"{b['band']} {series}"].values[0]

    hr = cfg['heatmap_log_halfrange']
    fig, ax = plt.subplots(figsize=cfg['heatmap_figsize'])
    # LogNorm centred on 1: a factor-of-N suppression and a factor-of-N
    # enhancement then sit equally far either side of the neutral colour.
    im = ax.imshow(M, cmap=cfg['heatmap_cmap'],
                   norm=LogNorm(vmin=1.0 / hr, vmax=hr),
                   aspect='auto', interpolation='nearest')

    # Undefined cells (pre-eruption band below the amplitude floor) get an
    # explicit flat colour, so an absent ratio never reads as a small one.
    undef = cfg.get('heatmap_undefined_color')
    if undef:
        for i, j in np.argwhere(np.isnan(M)):
            ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1,
                                       facecolor=undef, edgecolor='none',
                                       zorder=2))
            ax.text(j, i, 'n/a', ha='center', va='center',
                    fontsize=5.8, color='0.45', style='italic', zorder=3)

    ax.set_xticks(range(len(bands)))
    ax.set_xticklabels([b['label'] for b in bands], rotation=35, ha='right')
    ax.set_yticks(range(len(order)))
    if cfg.get('heatmap_group_by', 'atmosphere') == 'atmosphere':
        ax.set_yticklabels([c for c, _ in order], fontsize=7)
    else:
        ax.set_yticklabels([f"{c}  {a}" for c, a in order], fontsize=7)

    if cfg.get('heatmap_annotate'):
        for i in range(M.shape[0]):
            for j in range(M.shape[1]):
                if np.isnan(M[i, j]):
                    continue
                r, g, b_, _ = im.cmap(im.norm(M[i, j]))
                shade = ('white'
                         if (0.2126 * r + 0.7152 * g + 0.0722 * b_) < 0.5
                         else 'black')
                # Two decimals below 1 (0.05 and 0.09 are different stories),
                # one above, where the interesting range is 1.5-20.
                v = M[i, j]
                txt = f"{v:.2f}" if v < 1 else f"{v:.1f}"
                ax.text(j, i, txt, ha='center', va='center',
                        fontsize=6.2, color=shade, zorder=3)

    for k in block_edges[:-1]:
        ax.axhline(k - 0.5, color='black', lw=1.6)
    if cfg.get('heatmap_group_by', 'atmosphere') == 'atmosphere':
        for start, _ in block_labels:
            ax.axhline(start + len(RUNGS) - 0.5, color='0.45', lw=0.7)
        for start, label in block_labels:
            ax.text(-0.155, start + (2 * len(RUNGS) - 1) / 2.0, label,
                    transform=ax.get_yaxis_transform(), rotation=90,
                    ha='center', va='center', fontsize=8)

    # Ticks derived from the configured half-range, so changing it in the
    # YAML cannot leave the labels describing a scale that is no longer drawn.
    _t, _v = [1.0], 1.0
    while _v * 2 <= hr + 1e-9:
        _v *= 2
        _t = [1.0 / _v] + _t + [_v]
    cb = fig.colorbar(im, ax=ax, pad=0.02, fraction=0.046, ticks=_t)
    cb.ax.set_yticklabels([('1' if abs(t - 1) < 1e-9 else
                            (f'1/{int(round(1 / t))}' if t < 1
                             else f'{int(round(t))}')) for t in _t])
    cb.set_label('feature amplitude / pre-eruption amplitude')
    ax.set_title('Feature amplitude at its most-muted epoch\n'
                 '(R = 250; below 1 = suppressed by aerosol continuum)',
                 fontsize=9, pad=9)
    save(fig, cfg['outfile_b'])


if args.only in (None, 'a'):
    figure_a()
if args.only in (None, 'b'):
    figure_b()
