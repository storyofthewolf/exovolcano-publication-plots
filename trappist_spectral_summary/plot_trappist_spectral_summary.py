"""
plot_trappist_spectral_summary.py -- three candidate summary figures for the
24-case spectral table, drawn from one data source so they can be compared.

The table (exovolcano-spectra/data/band_table_24.csv) holds, per case and
band, the largest absolute departure of the R=250-binned transit spectrum from
that case's own day-0 spectrum, signed at the peak, with the epoch it occurred
on. Three axes are crossed: eruption template (Tambora / Hunga) x magnitude
rung (1x / 10x / 100x) x atmosphere (ben1 / ben2 / hab1 / hab2).

  A  fig_spectral_ladder   peak ppm vs rung; 3 bands x 2 eruptions, one line
                           per atmosphere. Reads the SCALING -- ben1's SO2
                           climbs while its aerosol stays flat, hab1's does
                           the reverse.
  B  fig_spectral_heatmap  24 cases x 7 bands, colour = ppm, diverging about
                           zero. Reads EVERYTHING, at the cost of reducing
                           each number to a colour.
  C  fig_spectral_timing   peak ppm vs day-of-peak. Reads the TIMING regimes,
                           which neither the table nor A nor B shows.

    python plot_trappist_spectral_summary.py            # all three
    python plot_trappist_spectral_summary.py --only a

EVERY FIGURE HERE INHERITS ONE CAVEAT. Spectra exist only at the 11 prescribed
epochs, so each "peak" is a peak over those epochs -- a lower bound on the true
maximum, not the maximum. Figure C marks the cases where that matters most by
drawing endpoint peaks hollow. See the config header for the measured
breakdown of how many maxima are broad, sharp, or at an endpoint.
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
ap.add_argument('--config', default='config_trappist_spectral_summary.yaml')
ap.add_argument('--only', default=None, choices=('a', 'b', 'c'),
                help="draw just one of the three")
args = ap.parse_args()
with open(os.path.join(here, args.config)) as f:
    cfg = yaml.safe_load(f)

df = pd.read_csv(cfg['table_csv'])
ATMS = cfg['atmospheres']
ACOL = cfg['atm_colors']
ALAB = cfg['atm_labels']
AMRK = cfg['atm_markers']
RUNGS = cfg['rungs']


def case_name(ladder, rung):
    """'tambora' + 10 -> 'tambora_10x'; the 1x rung carries no suffix."""
    return ladder if rung == 1 else f"{ladder}_{rung}x"


def save(fig, stem):
    path = os.path.join(here, stem)
    for ext in ('pdf', 'eps'):
        fig.savefig(f'{path}.{ext}', bbox_inches='tight')
        print(f'  wrote {stem}.{ext}')
    plt.close(fig)


# ---------------------------------------------------------------------------
# A -- ladder plot
# ---------------------------------------------------------------------------
def figure_a():
    print('A: ladder (peak ppm vs magnitude rung)')
    rows, cols = cfg['ladder_rows'], cfg['ladder_cols']
    fig, axes = plt.subplots(len(rows), len(cols),
                             figsize=cfg['ladder_figsize'],
                             sharex=True)
    x = np.arange(len(RUNGS))

    for r, rowspec in enumerate(rows):
        band = rowspec['band']
        # One y-scale across the two eruption columns of a row, so the reader
        # can see that (e.g.) the Hunga ladder barely moves the SO2 band that
        # the Tambora ladder drives hard. Not shared across rows: the aerosol
        # continuum and SO2 differ by an order of magnitude on dry cases.
        rowmax = 0.0
        for colspec in cols:
            for rung in RUNGS:
                sub = df[df.case == case_name(colspec['ladder'], rung)]
                if not sub.empty:
                    rowmax = max(rowmax, np.abs(sub[band].values).max())

        for c, colspec in enumerate(cols):
            ax = axes[r, c]
            for atm in ATMS:
                ys = []
                for rung in RUNGS:
                    s = df[(df.case == case_name(colspec['ladder'], rung))
                           & (df.atm == atm)]
                    ys.append(s[band].values[0] if not s.empty else np.nan)
                ax.plot(x, ys, marker=AMRK[atm], color=ACOL[atm],
                        lw=1.3, ms=4.5, label=ALAB[atm] if (r == 0 and c == 0)
                        else None, zorder=3)

            ax.axhline(0.0, color='0.7', lw=0.6, zorder=1)
            ax.set_xticks(x)
            ax.set_xticklabels(cfg['rung_labels'])
            ax.set_xlim(-0.25, len(RUNGS) - 0.75)
            ax.set_ylim(min(-0.06 * rowmax, 1.15 * min(
                0.0, *[df[(df.case == case_name(cs['ladder'], rg))][band].min()
                       for cs in cols for rg in RUNGS])),
                1.18 * rowmax)
            ax.grid(axis='y', color='0.9', lw=0.5, zorder=0)
            ax.set_axisbelow(True)

            if r == 0:
                ax.set_title(colspec['label'], fontsize=8.5, pad=7)
            if c == 0:
                ax.set_ylabel(rowspec['label'] + "\npeak $\\Delta$ depth [ppm]")
            else:
                ax.tick_params(labelleft=False)
            if r == len(rows) - 1:
                ax.set_xlabel('magnitude rung')

    axes[0, 0].legend(loc='upper left', frameon=False, handlelength=1.8)
    fig.subplots_adjust(hspace=0.14, wspace=0.06)
    save(fig, cfg['outfile_a'])


# ---------------------------------------------------------------------------
# B -- heatmap
# ---------------------------------------------------------------------------
def figure_b():
    print('B: heatmap (24 cases x 7 bands)')
    bands = cfg['heatmap_bands']
    # Row grouping. "atmosphere" blocks the four atmospheres, each holding its
    # six eruptions -- so a block is one planet and reading DOWN it is that
    # planet's full response to every eruption, which is how the suite is
    # actually argued. "eruption" was the original, blocking by eruption and
    # rung with the four atmospheres inside; it answers the transposed
    # question (how do the planets differ at fixed forcing) and is kept
    # because that is the comparison the per-eruption figures make.
    grouping = cfg.get('heatmap_group_by', 'atmosphere')
    order, block_edges, block_labels = [], [], []
    if grouping == 'atmosphere':
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
            M[i, j] = s[b['band']].values[0]

    v = cfg['heatmap_vlim']
    fig, ax = plt.subplots(figsize=cfg['heatmap_figsize'])
    # Sequential from 0, so the full colour range covers the positive values
    # that are 157 of the 168 entries. Negatives (11 entries, none beyond
    # -2.3 ppm) would otherwise clamp to the bottom colour and be
    # indistinguishable from a true zero, so they are overpainted below.
    im = ax.imshow(M, cmap=cfg['heatmap_cmap'], vmin=0.0, vmax=v,
                   aspect='auto', interpolation='nearest')

    negc = cfg.get('heatmap_negative_color')
    if negc:
        neg = np.argwhere(M < 0)
        for i, j in neg:
            ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1,
                                       facecolor=negc, edgecolor='none',
                                       zorder=2))

    ax.set_xticks(range(len(bands)))
    ax.set_xticklabels([b['label'] for b in bands], rotation=35, ha='right')
    ax.set_yticks(range(len(order)))
    # Within an atmosphere block the atmosphere is constant, so repeating it on
    # every row is noise -- the block label at the left carries it instead.
    if grouping == 'atmosphere':
        ax.set_yticklabels([c for c, _ in order], fontsize=7)
    else:
        ax.set_yticklabels([f"{c}  {a}" for c, a in order], fontsize=7)

    if cfg.get('heatmap_annotate'):
        for i in range(M.shape[0]):
            for j in range(M.shape[1]):
                if np.isnan(M[i, j]):
                    continue
                # Text colour follows the CELL'S OWN luminance rather than a
                # fixed value threshold. A threshold only works for a colormap
                # that darkens monotonically with value; on viridis, which is
                # darkest at ZERO, it put black text on the darkest cells in
                # the map. Relative luminance (Rec. 709) asks the question
                # that actually matters -- is this background dark? -- and is
                # correct for any colormap.
                if M[i, j] < 0 and negc:
                    # Negative cells are overpainted in a flat colour, so ask
                    # that colour, not the colormap it never receives.
                    r, g, b_ = matplotlib.colors.to_rgb(negc)
                else:
                    r, g, b_, _ = im.cmap(im.norm(M[i, j]))
                shade = ('white'
                         if (0.2126 * r + 0.7152 * g + 0.0722 * b_) < 0.55
                         else 'black')
                # Parenthesised, the accounting convention for a negative, so
                # the sign survives even in a greyscale print where the
                # negative-cell colour is just another grey.
                txt = (f"({abs(M[i, j]):.0f})" if M[i, j] < 0
                       else f"{M[i, j]:.0f}")
                ax.text(j, i, txt, ha='center', va='center',
                        fontsize=6.2, color=shade, zorder=3)

    # Rule between blocks so the eye does not read across a boundary. Under
    # atmosphere grouping the heavy rules separate the four planets and a
    # lighter rule marks the Tambora|Hunga seam inside each, since the two
    # ladders scale different quantities and are not rung-for-rung comparable.
    for k in block_edges[:-1]:
        ax.axhline(k - 0.5, color='black', lw=1.6)
    if grouping == 'atmosphere':
        for start, _ in block_labels:
            ax.axhline(start + len(RUNGS) - 0.5, color='0.45', lw=0.7)
        # Name each block once, outside the axes on the left.
        for start, label in block_labels:
            ax.text(-0.145, start + (2 * len(RUNGS) - 1) / 2.0, label,
                    transform=ax.get_yaxis_transform(),
                    rotation=90, ha='center', va='center', fontsize=8)

    cb = fig.colorbar(im, ax=ax, pad=0.02, fraction=0.046)
    cb.set_label('peak $\\Delta$ transit depth from day 0 [ppm]')
    # Plain text, not TeX: matplotlib's default renderer has no \, spacing
    # macro and printed it literally as "R\,=\,250".
    ax.set_title('Peak spectral response by case and band\n'
                 '(R = 250; peak over sampled epochs)',
                 fontsize=9, pad=9)
    save(fig, cfg['outfile_b'])


# ---------------------------------------------------------------------------
# C -- peak value against time of peak
# ---------------------------------------------------------------------------
def figure_c():
    print('C: timing (peak ppm vs day of peak)')
    bands = cfg['scatter_bands']
    ends = set(cfg['scatter_endpoint_days'])
    smap = {int(k): v for k, v in cfg['scatter_size_by_rung'].items()}

    fig, axes = plt.subplots(1, len(bands), figsize=cfg['scatter_figsize'],
                             sharey=False)
    for k, b in enumerate(bands):
        ax = axes[k]
        for ladder in ('tambora', 'hunga'):
            for rung in RUNGS:
                case = case_name(ladder, rung)
                for atm in ATMS:
                    s = df[(df.case == case) & (df.atm == atm)]
                    if s.empty:
                        continue
                    day = float(s[b['band'] + ' day'].values[0])
                    val = float(s[b['band']].values[0])
                    if day <= 0:
                        continue
                    endpoint = int(day) in ends
                    # Hollow marker = the peak sits at an endpoint of the
                    # sampled epoch set, so the true maximum may lie outside
                    # the sampled range. This is the sampling caveat made
                    # visible rather than left to the caption.
                    ax.scatter(day, val,
                               s=smap.get(rung, 40),
                               marker=AMRK[atm],
                               facecolor='none' if endpoint else ACOL[atm],
                               edgecolor=ACOL[atm],
                               linewidth=1.1 if endpoint else 0.5,
                               zorder=3)
        ax.set_xscale('log')
        ax.axhline(0.0, color='0.7', lw=0.6, zorder=1)
        ax.grid(color='0.92', lw=0.5, zorder=0)
        ax.set_axisbelow(True)
        ax.set_xlabel('day of peak')
        ax.set_title(b['label'], fontsize=8.5, pad=6)
        if k == 0:
            ax.set_ylabel('peak $\\Delta$ transit depth [ppm]')

    handles = [Line2D([], [], marker=AMRK[a], color=ACOL[a], ls='none',
                      ms=5.5, label=ALAB[a]) for a in ATMS]
    handles += [
        Line2D([], [], marker='o', color='0.35', ls='none', ms=3.4,
               label='1$\\times$'),
        Line2D([], [], marker='o', color='0.35', ls='none', ms=5.2,
               label='10$\\times$'),
        Line2D([], [], marker='o', color='0.35', ls='none', ms=7.2,
               label='100$\\times$'),
        Line2D([], [], marker='o', mfc='none', mec='0.35', color='0.35',
               ls='none', ms=5.5, label='peak at sampling endpoint'),
    ]
    fig.legend(handles=handles, loc='lower center', ncol=4, frameon=False,
               bbox_to_anchor=(0.5, -0.16))
    fig.subplots_adjust(wspace=0.26)
    save(fig, cfg['outfile_c'])


which = args.only
if which in (None, 'a'):
    figure_a()
if which in (None, 'b'):
    figure_b()
if which in (None, 'c'):
    figure_c()
