"""
plot_trappist_spectra_balances.py -- the four balances that govern the volcanic
transmission signal on TRAPPIST-1 e, in one 2 x 4 figure.

The per-case spectral time series carry far more information than a manuscript
figure can spend space on. This script reduces the sixteen of them to the four
contrasts that actually make the argument, one per column, with the GCM burdens
that cause them directly underneath:

  col 1  DRY  + sulfur-rich (ben1 tambora_100x)
         Water is the limiting reactant. Aerosol production is capped at what
         the eruption's own ~1 Tg of H2O allows, and the excess SO2 -- 6.0e12 kg
         of it -- simply sits in the gas phase for the whole six years, where
         the spectrum sees it as growing nu3 and nu1 bands.

  col 2  DRY  + water-rich (ben1 hunga_100x)
         The same dry planet, the opposite eruption. Now sulfur is limiting and
         is consumed to zero, so little aerosol forms; but the injected water
         has nothing to remove it and persists in the stratosphere for years,
         showing up in the H2O bands rather than the SO2 ones.

  col 3  WET  + sulfur-rich (hab1 tambora_100x)
         Surface water lifts the oxidant limit entirely. SO2 converts nearly
         completely, aerosol grows to ~3.1e12 kg -- some 580x the dry cap -- and
         the spectrum shows both gas bands and a strong scattering continuum.

  col 4  WET, water-rich aloft (hab2 tambora_100x)
         Conversion is likewise unlimited, but hab2's large ambient water aloft
         rains the aerosol out within months. Every band decays to zero by
         ~day 300: a bright, brief event rather than an enduring one.

Columns 3 and 4 are deliberately the same eruption on two different wet
atmospheres. A literal 2x2 of the case matrix would put hab2 hunga_100x in the
corner -- the weakest signal in the suite -- and lose the rainout contrast,
which is the fourth point being made.

TWO REPOSITORIES FEED THIS FIGURE. Spectra come from exovolcano-spectra's
computed .npz files, burdens from remote_analysis' scalar CSVs. The .npz files
are gitignored there (regenerable, multi-GB); if they are missing this script
says which command rebuilds them instead of drawing an incomplete figure.

Configuration is in config_trappist_spectra_balances.yaml; see that file's
header for the R = 250 binning convention and the two-ladder trap.
"""

import os
import sys
import argparse
import yaml
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.patheffects as pe

plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 8,
    'axes.labelsize': 9,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'legend.fontsize': 7,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})

here = os.path.dirname(os.path.abspath(__file__))
parser = argparse.ArgumentParser()
parser.add_argument('--config', default='config_trappist_spectra_balances.yaml')
args = parser.parse_args()
with open(os.path.join(here, args.config)) as f:
    cfg = yaml.safe_load(f)

COLS = cfg['columns']
EPOCHS = cfg['epochs_days']
ECOL = cfg['epoch_colors']


def bin_to_resolution(wl, depth, resolution):
    """Bin spectra to a constant resolving power by AVERAGING within each bin.

    Averaging, not sampling: a transit depth is an area, and the mean over a bin
    is the area a spectrograph of that resolving power would record. Sampling
    the nearest point instead lands on arbitrary line cores and produces a
    spectrum that depends on grid alignment rather than on physics.
    """
    edges = [wl.min()]
    while edges[-1] < wl.max():
        edges.append(edges[-1] * (1.0 + 1.0 / resolution))
    edges = np.array(edges)
    centres, binned = [], []
    for b in range(len(edges) - 1):
        sel = (wl >= edges[b]) & (wl < edges[b + 1])
        if not sel.any():
            continue
        centres.append(wl[sel].mean())
        binned.append(depth[:, sel].mean(axis=1))
    return np.array(centres), np.array(binned).T


def load_spectra(col):
    path = os.path.join(cfg['spectra_dir'], col['key'] + cfg['spectra_suffix'])
    if not os.path.exists(path):
        sys.exit(
            f"MISSING computed spectra: {path}\n"
            f"These are gitignored in exovolcano-spectra (regenerable). Rebuild with:\n"
            f"  cd ~/Desktop/projects/volcanos/exovolcano-spectra && "
            f"scripts/run_all_suites.sh\n"
            f"Refusing to draw a partial figure."
        )
    z = np.load(path, allow_pickle=True)
    wl = z['wavelength_um']
    days = z['days']
    depth = z['depth'] * 1e6                      # fraction -> ppm
    wlb, depb = bin_to_resolution(wl, depth, cfg['resolution'])
    return wlb, days, depb


def load_burden(col, var):
    """One scalar burden time series, in Tg. Returns (days, Tg) or (None, None)."""
    path = os.path.join(cfg['burden_root'], f"{col['atm']}_suite1",
                        f"exovolc_{col['atm']}_{col['case']}",
                        cfg['subpath_scalar'], f"{var}.csv")
    if not os.path.exists(path):
        print(f"  WARNING: no {var} for {col['key']}")
        return None, None
    df = pd.read_csv(path, comment='#')
    return df.iloc[:, 0].values, df.iloc[:, 1].values * cfg['kg_to_tg']


ncol = len(COLS)
fig, axes = plt.subplots(2, ncol, figsize=cfg['figsize'],
                         gridspec_kw={'height_ratios': [1.0, 0.85]})

wl_lo, wl_hi = cfg['wl_range_um']

# ---------------------------------------------------------------------------
# TOP ROW -- difference spectra
# ---------------------------------------------------------------------------
print('top row: transmission spectra (difference from day 0)')
spec_lims = []
for i, col in enumerate(COLS):
    ax = axes[0, i]
    wlb, days, depb = load_spectra(col)
    d0 = depb[0]

    m = (wlb >= wl_lo) & (wlb <= wl_hi)

    # Band spans first, so the spectra draw over them.
    for band in cfg['bands']:
        ax.axvspan(band['lo'], band['hi'], color=band['color'], alpha=0.13, lw=0)
        if i == 0 and band.get('label'):
            # Two rows, so the adjacent 5.5-8.9 um spans do not collide.
            y = 0.972 - 0.062 * band.get('row', 0)
            ax.text(0.5 * (band['lo'] + band['hi']), y, band['label'],
                    transform=ax.get_xaxis_transform(), ha='center', va='top',
                    fontsize=6.5, color=band['color'], zorder=5,
                    path_effects=[pe.withStroke(linewidth=1.8,
                                                foreground='white')])

    ax.axhline(0.0, color='0.55', lw=0.6, zorder=1)

    peak = 0.0
    for e, day in enumerate(EPOCHS):
        j = int(np.argmin(np.abs(days - day)))
        # Snapping to the nearest computed epoch is silent by nature, and a
        # stale .npz would otherwise plot day 500 under a "day 2000" label.
        if days[j] != day:
            print(f"    NOTE: epoch {day} d not computed; using {days[j]} d. "
                  f"Available: {list(days)}")
        diff = depb[j] - d0
        ax.plot(wlb[m], diff[m], color=ECOL[e], lw=0.9,
                label=f"day {days[j]}", zorder=3)
        peak = max(peak, np.abs(diff[m]).max())
    spec_lims.append(peak)

    ax.set_xscale('log')
    ax.set_xlim(wl_lo, wl_hi)
    ax.set_xticks([1, 2, 3, 5, 7, 10])
    ax.set_xticklabels(['1', '2', '3', '5', '7', '10'])
    ax.set_title(f"{col['balance']}\n"
                 + r"$\mathrm{\mathsf{" + col['case'].replace('_', r'\_') + r"}}$"
                 + f"  ({col['atm']})",
                 fontsize=8.5, pad=6)
    if i == 0:
        ax.set_ylabel(cfg['spectra_ylabel'])
        ax.legend(loc='upper left', frameon=False, handlelength=1.4,
                  bbox_to_anchor=(0.0, 0.86))
    ax.set_xlabel(r'wavelength [$\mu$m]')
    print(f"  {col['key']:30s} peak |diff| {peak:6.1f} ppm")

# One shared y-scale across the top row, so panel heights mean the same thing.
ytop = 1.12 * max(spec_lims)
for i in range(ncol):
    axes[0, i].set_ylim(-0.18 * ytop, ytop)
    if i > 0:
        axes[0, i].set_yticklabels([])

# ---------------------------------------------------------------------------
# BOTTOM ROW -- GCM burdens
# ---------------------------------------------------------------------------
print('bottom row: GCM burdens')
for i, col in enumerate(COLS):
    ax = axes[1, i]
    for spec in cfg['burden_vars']:
        d, y = load_burden(col, spec['var'])
        if d is None:
            continue
        # A log axis cannot show zero. On the dry atmospheres Q and the aerosol
        # are exactly zero before the eruption delivers them; masking rather
        # than clipping leaves an honest gap instead of drawing a floor value
        # that was never in the model.
        y = np.where(y > 0, y, np.nan)
        ax.plot(d, y, color=spec['color'], ls=spec['ls'], lw=1.1,
                label=spec['label'])
        finite = y[np.isfinite(y)]
        if finite.size:
            print(f"  {col['key']:30s} {spec['var']:9s} "
                  f"max {finite.max():.3e} Tg  end {y[-1]:.3e} Tg")

    if cfg.get('mark_epochs'):
        for e, day in enumerate(EPOCHS):
            ax.axvline(day, color=ECOL[e], lw=0.6, alpha=0.55, zorder=0)

    ax.set_yscale('log')
    ax.set_xlim(*cfg['burden_xlim'])
    ax.set_ylim(*cfg['burden_ylim'])
    ax.set_xlabel(cfg['burden_xlabel'])
    ax.text(0.035, 0.055, col['note'], transform=ax.transAxes,
            ha='left', va='bottom', fontsize=7.5, style='italic', color='0.25')
    if i == 0:
        ax.set_ylabel(cfg['burden_ylabel'])
        ax.legend(loc='upper right', frameon=False, handlelength=1.6)
    else:
        ax.set_yticklabels([])

fig.tight_layout(rect=[0, 0, 1, 0.99])
stem = os.path.join(here, cfg['outfile_stem'])
for ext in ('pdf', 'eps'):
    fig.savefig(f'{stem}.{ext}', bbox_inches='tight')
    print(f'wrote {stem}.{ext}')
