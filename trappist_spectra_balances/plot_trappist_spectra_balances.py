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
import matplotlib.colors as mcolors

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


def _blend(color, frac, bg='white'):
    """Opaque stand-in for `color` at opacity `frac` over `bg`.

    Used instead of alpha= so that the EPS and PDF backends agree; see the
    epoch-line call site for why transparency cannot be used here.
    """
    c = np.array(mcolors.to_rgb(color))
    b = np.array(mcolors.to_rgb(bg))
    return tuple(frac * c + (1.0 - frac) * b)


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
# Two different row gaps are needed, and a single hspace cannot give them:
# rows 1 and 2 share the wavelength abscissa (row 1 hides its tick labels and
# row 2 carries the one x-label for both), so they close up tight and read as
# one block; row 3 is a different abscissa in days and has to clear row 2's
# x-label as well as its own ticks. An outer 2-row gridspec sets the wide gap,
# and a nested 2-row spec inside its first cell sets the tight one.
_hr = cfg.get('height_ratios', [1, 1, 0.85])
fig = plt.figure(figsize=cfg['figsize'])
_outer = fig.add_gridspec(2, 1,
                          height_ratios=[_hr[0] + _hr[1], _hr[2]],
                          hspace=cfg.get('hspace_rows23', 0.17))
_top = _outer[0].subgridspec(2, ncol, height_ratios=_hr[:2],
                             hspace=cfg.get('hspace_rows12', 0.05),
                             wspace=cfg.get('wspace', 0.08))
_bot = _outer[1].subgridspec(1, ncol, wspace=cfg.get('wspace', 0.08))
axes = np.empty((3, ncol), dtype=object)
for _r in range(2):
    for _c in range(ncol):
        axes[_r, _c] = fig.add_subplot(_top[_r, _c])
for _c in range(ncol):
    axes[2, _c] = fig.add_subplot(_bot[0, _c])

wl_lo, wl_hi = cfg['wl_range_um']

def draw_bands(ax, label_them, y_base=0.972):
    """Shade the diagnostic bands; label them only on the panel that asks."""
    for band in cfg['bands']:
        # The band colours are already very light, so they are laid down at
        # full opacity rather than at the old alpha=0.13 over a saturated
        # colour: alpha on a pale fill washes it out to nothing. The label
        # takes its own darker colour, since a 6.5 pt glyph in the fill colour
        # would be illegible.
        ax.axvspan(band['lo'], band['hi'], color=band['color'], zorder=0, lw=0)
        if label_them and band.get('label'):
            # Two rows, so the adjacent 5.5-8.9 um spans do not collide.
            y = y_base - 0.062 * band.get('row', 0)
            # Centre across this band and any UNLABELLED spans that follow it,
            # which are its continuation: SO2 is two bands (7.0-7.7, 8.4-8.9)
            # and centring on the first alone threw the label to the left.
            # Grouping by colour instead would be wrong, since H2O shares its
            # blue with the unrelated 5.5-7.0 span.
            _members = [band]
            _bands = cfg['bands']
            for _b in _bands[_bands.index(band) + 1:]:
                if _b.get('label'):
                    break
                if _b['color'] == band['color']:
                    _members.append(_b)
            _xc = 0.5 * (min(b['lo'] for b in _members)
                         + max(b['hi'] for b in _members))
            ax.text(_xc, y, band['label'],
                    transform=ax.get_xaxis_transform(), ha='center', va='top',
                    fontsize=cfg.get('band_label_size', 8.0),
                    color=cfg.get('band_label_color', 'black'), zorder=5,
                    path_effects=[pe.withStroke(linewidth=2.2,
                                                foreground='white')])


def style_spectral_axis(ax):
    ax.set_xscale('log')
    ax.set_xlim(wl_lo, wl_hi)
    ax.set_xticks([1, 2, 3, 5, 7, 10])
    ax.set_xticklabels(['1', '2', '3', '5', '7', '10'])


# ---------------------------------------------------------------------------
# ROW 1 -- RAW transit depth
#
# Cached per column so the difference row below does not re-read and re-bin the
# same .npz; binning 47k points x 11 epochs is the slow step here.
# ---------------------------------------------------------------------------
print('row 1: raw transmission spectra')
SPEC = {}
raw_lo, raw_hi = np.inf, -np.inf
for i, col in enumerate(COLS):
    ax = axes[0, i]
    wlb, days, depb = load_spectra(col)
    SPEC[i] = (wlb, days, depb)
    m = (wlb >= wl_lo) & (wlb <= wl_hi)

    # Band labels ride on the top-left panel (author request 2026-09-04).
    # The y-limit block below adds headroom so they clear both the spectra
    # and this panel's epoch legend.
    draw_bands(ax, label_them=(i == 0), y_base=0.985)

    geo = None
    if cfg.get('show_geometric'):
        z = np.load(os.path.join(cfg['spectra_dir'],
                                 col['key'] + cfg['spectra_suffix']),
                    allow_pickle=True)
        geo = float(z['geometric_ppm'])
        ax.axhline(geo, color='0.45', lw=0.7, ls=':', zorder=2)

    for e, day in enumerate(EPOCHS):
        j = int(np.argmin(np.abs(days - day)))
        base = (day == 0)
        ax.plot(wlb[m], depb[j][m], color=ECOL[e],
                lw=1.0 if base else 0.9,
                ls='--' if base else '-',
                zorder=2 if base else 3,
                label=f"day {days[j]}" + (" (pre-eruption)" if base else ""))
        raw_lo = min(raw_lo, depb[j][m].min())
        raw_hi = max(raw_hi, depb[j][m].max())
    if geo is not None:
        raw_lo = min(raw_lo, geo)

    style_spectral_axis(ax)
    ax.set_title(f"{col['balance']}\n"
                 + r"$\mathrm{\mathsf{" + col['case'].replace('_', r'\_') + r"}}$"
                 + f"  ({col['atm']})",
                 fontsize=8.5, pad=6)
    if i == 0:
        ax.set_ylabel(cfg['raw_ylabel'])
        ax.legend(loc='upper left', frameon=False, handlelength=1.4,
                  bbox_to_anchor=(0.0, 0.86))

# One shared scale across the raw row: the four sit on the same geometric floor,
# and per-column limits would conceal that hab2's quiet continuum already stands
# well above ben1's before any eruption.
pad = 0.06 * (raw_hi - raw_lo)
# Extra headroom at the top only, so the band labels on the leftmost panel clear
# both the spectra and the epoch legend. Author sanctioned ~10 ppm for this.
head = float(cfg.get('raw_headroom_ppm', 10.0))
for i in range(ncol):
    axes[0, i].set_ylim(raw_lo - pad, raw_hi + pad + head)
    axes[0, i].tick_params(labelbottom=False)
    if i > 0:
        axes[0, i].set_yticklabels([])
if cfg.get('show_geometric') and cfg.get('geometric_label'):
    # Label is optional and off by default: the dotted reference line reads
    # clearly on its own and the caption already says what it is, while the
    # text had to sit outside the axes (inside, it lands on a curve, since the
    # line runs along the bottom where every case's continuum also runs).
    axes[0, ncol - 1].annotate(
        cfg['geometric_label'],
        xy=(1.005, geo), xycoords=('axes fraction', 'data'),
        fontsize=6.0, color='0.4', ha='left', va='center')
print(f"  shared raw range {raw_lo:.1f}-{raw_hi:.1f} ppm")

# ---------------------------------------------------------------------------
# ROW 2 -- difference spectra
# ---------------------------------------------------------------------------
print('row 2: transmission spectra (difference from day 0)')
spec_lims = []
for i, col in enumerate(COLS):
    ax = axes[1, i]
    wlb, days, depb = SPEC[i]
    d0 = depb[0]

    m = (wlb >= wl_lo) & (wlb <= wl_hi)
    # Labelled once, on row 2: its top margin is clear, whereas row 1's
    # curves run right to the axis top.
    draw_bands(ax, label_them=False)

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
        base = (day == 0)
        # Day 0 IS the reference here, so its difference is identically zero;
        # it coincides with the zero line and is drawn only to make that
        # explicit, never as a curve carrying information.
        ax.plot(wlb[m], diff[m], color=ECOL[e],
                lw=1.0 if base else 0.9,
                ls='--' if base else '-',
                zorder=2 if base else 3)
        if not base:
            peak = max(peak, np.abs(diff[m]).max())
    spec_lims.append(peak)

    style_spectral_axis(ax)
    if i == 0:
        ax.set_ylabel(cfg['spectra_ylabel'])
    ax.set_xlabel(r'wavelength [$\mu$m]')
    print(f"  {col['key']:30s} peak |diff| {peak:6.1f} ppm")


# One shared y-scale across the top row, so panel heights mean the same thing.
ytop = 1.26 * max(spec_lims)
for i in range(ncol):
    axes[1, i].set_ylim(-0.18 * ytop, ytop)
    if i > 0:
        axes[1, i].set_yticklabels([])

# ---------------------------------------------------------------------------
# BOTTOM ROW -- GCM burdens
# ---------------------------------------------------------------------------
print('row 3: GCM burdens')
for i, col in enumerate(COLS):
    ax = axes[2, i]
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
        # Pre-blend toward white instead of using alpha. The PostScript backend
        # drops transparency entirely ("partially transparent artists will be
        # rendered opaque"), so an alpha'd line renders full-strength in the EPS
        # and washed-out in the PDF, and the manuscript, which includes the EPS,
        # then disagrees with the reference PDF built here. Blending gives both
        # backends the same colour.
        for e, day in enumerate(EPOCHS):
            ax.axvline(day, color=_blend(ECOL[e], 0.55), lw=0.6, zorder=0)

    ax.set_yscale('log')
    ax.set_xlim(*cfg['burden_xlim'])
    ax.set_ylim(*cfg['burden_ylim'])
    ax.set_xlabel(cfg['burden_xlabel'])
    ax.text(0.035, 0.055, col['note'], transform=ax.transAxes,
            ha='left', va='bottom', fontsize=7.5, style='italic', color='0.25')
    if i == 0:
        ax.set_ylabel(cfg['burden_ylabel'])
        # Left-aligned per author request. Placed below the water line
        # rather than at the top: on the hunga_wet variant (Fig. A8) the
        # H2O burden runs flat at ~1e6-1e7 across the whole panel and an
        # upper legend collides with it, while the panel below ~1e4 is
        # empty in every case here.
        ax.legend(loc=cfg.get('burden_legend_loc', 'lower left'),
                  frameon=False, handlelength=1.6,
                  bbox_to_anchor=tuple(cfg['burden_legend_anchor'])
                  if cfg.get('burden_legend_anchor') else None)
    else:
        ax.set_yticklabels([])

# No tight_layout here: it recomputes every gap and would discard the two
# distinct row spacings set on the gridspecs above. Margins are set explicitly
# instead, and bbox_inches='tight' at save time trims whatever is left over.
fig.subplots_adjust(left=0.055, right=0.995, top=0.945, bottom=0.075)
stem = os.path.join(here, cfg['outfile_stem'])
for ext in ('pdf', 'eps'):
    fig.savefig(f'{stem}.{ext}', bbox_inches='tight')
    print(f'wrote {stem}.{ext}')
