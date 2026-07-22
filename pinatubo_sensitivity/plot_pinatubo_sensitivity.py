"""
plot_pinatubo_sensitivity.py - Sensitivity of global-mean AOD to the four
Pinatubo sweep parameters: SO2 -> H2SO4 conversion timescale K (the
prescribed e-folding time for SO2 oxidation), aerosol effective radius Reff,
injected SO2 mass, and injection pressure.

2x2 panels, one sweep per panel, each drawn on a shared y-axis so peak
magnitudes are directly comparable across panels. The fiducial case is
highlighted with a heavier line in every panel it appears in.

Reads:  base_dir/<case>/data/aod/aod_550nm_band.csv

Configuration (panels, cases, paths, output) is in
config_pinatubo_sensitivity.yaml.
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import argparse
import yaml
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.dates as mdates

plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 8,
    'axes.labelsize': 9,
    'axes.titlesize': 9,
    'xtick.labelsize': 7.5,
    'ytick.labelsize': 8,
    'legend.fontsize': 6.5,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})

here = os.path.dirname(os.path.abspath(__file__))

parser = argparse.ArgumentParser()
parser.add_argument('--config', default='config_pinatubo_sensitivity.yaml',
                    help='YAML config file (relative to script directory)')
args = parser.parse_args()

cfg = yaml.safe_load(open(os.path.join(here, args.config)))

base_dir = cfg['base_dir']
fiducial = cfg['fiducial']
t0 = pd.Timestamp(cfg['eruption_date'])
xlim = (pd.Timestamp(cfg['xlim'][0]), pd.Timestamp(cfg['xlim'][1]))

# Sequential, colorblind-safe color sets per panel (low -> high parameter
# value), independent of which member is the fiducial. Chosen from the
# Okabe-Ito / print-safe families used elsewhere in the repo.
panel_colors = {
    'k':    ['#56B4E9', '#0072B2', '#D55E00', '#7A0000'],   # 4 members
    'reff': ['#009E73', '#0072B2', '#D55E00'],               # 3 members
    'so2':  ['#009E73', '#0072B2', '#D55E00'],
    'pinj': ['#009E73', '#0072B2', '#D55E00'],
}

# ---------------------------------------------------------------------------
# Load data and print diagnostics
# ---------------------------------------------------------------------------
def load_ts(case):
    p = os.path.join(base_dir, case, 'data', cfg['subpath_timeseries'])
    if not os.path.exists(p):
        raise SystemExit(f"Missing time series CSV: {p}")
    df = pd.read_csv(p, header=0)
    days = df.iloc[:, 0].values.astype(float)
    aod = df.iloc[:, 1].values.astype(float)
    return days, aod

print("Peak AOD diagnostics by sweep:")
loaded = {}   # panel key -> list of (case, label, days, aod)
for panel in cfg['panels']:
    print(f"  -- {panel['title']} --")
    entries = []
    for c in panel['cases']:
        days, aod = load_ts(c['case'])
        pk_idx = np.nanargmax(aod)
        print(f"     {c['label']:14s} ({c['case']}): peak AOD = {aod[pk_idx]:.4f} "
              f"at day {days[pk_idx]:.0f}")
        entries.append((c['case'], c['label'], days, aod))
    loaded[panel['key']] = entries

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(7.1, 5.6), sharex=True, sharey=True)
axes = axes.ravel()

for ax, panel in zip(axes, cfg['panels']):
    colors = panel_colors[panel['key']]
    entries = loaded[panel['key']]
    for (case, label, days, aod), color in zip(entries, colors):
        dates = t0 + pd.to_timedelta(days, unit='D')
        is_fid = (case == fiducial)
        lw = 2.2 if is_fid else 1.3
        zorder = 5 if is_fid else 3
        ax.plot(dates, aod, lw=lw, color=color, zorder=zorder, label=label)

    ax.legend(fontsize=6.5, frameon=False, loc='upper right',
              bbox_to_anchor=(1.0, 0.88), handlelength=1.6,
              labelspacing=0.25, borderaxespad=0.2)
    ax.text(0.03, 0.96, panel['title'], transform=ax.transAxes, fontsize=8,
            fontweight='bold', va='top', ha='left')
    ax.set_xlim(*xlim)
    ax.set_ylim(0, 0.30)
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.xaxis.set_major_locator(mdates.YearLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    ax.xaxis.set_minor_locator(mdates.MonthLocator())

for ax in (axes[0], axes[2]):
    ax.set_ylabel('AOD at 550 nm')
for ax in (axes[2], axes[3]):
    ax.set_xlabel('Year')

fig.align_ylabels([axes[0], axes[2]])
fig.subplots_adjust(hspace=0.15, wspace=0.12)

# The shared x tick at the left/right column boundary (1994 on the left
# column's right edge, 1991 on the right column's left edge) collides
# because both columns' major tick labels are drawn at full width. Blank
# out the rightmost major tick label on the left-column bottom axis so
# only the right column's "1991" shows at that boundary.
fig.canvas.draw()
xticklabels = axes[2].xaxis.get_majorticklabels()
if xticklabels:
    xticklabels[-1].set_visible(False)

for ext in ('pdf', 'eps'):
    outfile = os.path.join(here, f"{cfg['outfile_stem']}.{ext}")
    plt.savefig(outfile, dpi=300, bbox_inches='tight')
    print(f"Saved: {outfile}")
plt.close()
