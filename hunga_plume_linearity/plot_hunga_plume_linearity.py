"""
plot_hunga_plume_linearity.py - Does the plume-tracked AOD scale linearly with
injected SO2, as the global-mean AOD does exactly?

Three panels:
(a) Plume-core AOD (the maximum over the longitude-resolved grid) versus days
    since injection, for the 0.5 Tg fiducial and the 1.0/1.5/2.0 Tg Sellitto
    SO2 axis. Shows the plume rise, peak near day 10, and decay.
(b) The same three Sellitto cases normalized by the 1.0 Tg case, against the
    dashed lines that exact proportionality would give (1.5x and 2.0x). The
    departure below those lines with plume age IS the result.
(c) Global-mean AOD for the same cases, normalized identically -- flat on the
    proportionality lines at every time, which is the contrast that makes
    panel (b) meaningful rather than a numerical artifact.

Reads:  <base_dir>/<case>/data/lonlat/aod_day<DDDD.DD>.csv   (lat x lon grid)
        <base_dir>/<case>/data/aod/aod_550nm_band.csv        (global mean)

The lon-lat CSVs come from exovolcano-analysis' lonlat_aod.py, which is
separate from the run_time_series.py pipeline.

This figure is a MODEL-BEHAVIOR diagnostic, not an observational validation:
plume-tracked AOD is not comparable to global- or zonal-mean SAOD. See the
config YAML header before adding any reference line.

Configuration is in config_hunga_plume_linearity.yaml.
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

plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 8,
    'axes.labelsize': 9,
    'axes.titlesize': 9,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'legend.fontsize': 7,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})

here = os.path.dirname(os.path.abspath(__file__))

parser = argparse.ArgumentParser()
parser.add_argument('--config', default='config_hunga_plume_linearity.yaml',
                    help='YAML config file (relative to script directory)')
args = parser.parse_args()

cfg = yaml.safe_load(open(os.path.join(here, args.config)))
base_dirs = cfg['base_dirs']
days = [float(d) for d in cfg['days']]


def lonlat_path(entry, day):
    """Path to one lon-lat snapshot CSV. Day is formatted as lonlat_aod.py writes it."""
    base = base_dirs[entry['phase']]
    return os.path.join(base, entry['case'], 'data', cfg['subpath_lonlat'],
                        f'aod_day{day:07.2f}.csv')


def plume_peak(entry, day):
    """Maximum AOD over the (lat, lon) grid -- the plume core."""
    p = lonlat_path(entry, day)
    if not os.path.exists(p):
        raise SystemExit(f"Missing lon-lat CSV: {p}")
    d = pd.read_csv(p, index_col=0)
    return float(np.nanmax(d.values))


def global_mean_at(entry, day):
    """Global-mean AOD interpolated to the same day, from the scalar time series."""
    base = base_dirs[entry['phase']]
    p = os.path.join(base, entry['case'], 'data', cfg['subpath_timeseries'])
    if not os.path.exists(p):
        raise SystemExit(f"Missing global-mean CSV: {p}")
    d = pd.read_csv(p)
    t, v = d.iloc[:, 0].values, d.iloc[:, 1].values
    return float(np.interp(day, t, v))


axis = cfg['so2_axis']
fid = cfg['fiducial']

peaks = {e['case']: np.array([plume_peak(e, d) for d in days]) for e in axis}
peaks[fid['case']] = np.array([plume_peak(fid, d) for d in days])
gmeans = {e['case']: np.array([global_mean_at(e, d) for d in days]) for e in axis}

ref = axis[0]           # the 1.0 Tg case both panels normalize against
ref_peak = peaks[ref['case']]
ref_gmean = gmeans[ref['case']]

print('Plume-core AOD (max over lon-lat grid):')
hdr = '  ' + 'case'.ljust(24) + ''.join(f'{d:>9.0f}' for d in days)
print(hdr)
for e in [fid] + axis:
    print('  ' + e['label'].ljust(24)
          + ''.join(f'{v:>9.4f}' for v in peaks[e['case']]))
print('\nRatio to 1.0 Tg  (plume core | global mean):')
for e in axis[1:]:
    rp = peaks[e['case']] / ref_peak
    rg = gmeans[e['case']] / ref_gmean
    print('  ' + e['label'].ljust(10)
          + 'plume ' + ''.join(f'{v:>8.3f}' for v in rp))
    print('  ' + ''.ljust(10)
          + 'gmean ' + ''.join(f'{v:>8.3f}' for v in rg))

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
fig, (ax_a, ax_b, ax_c) = plt.subplots(
    1, 3, figsize=(7.3, 2.5),
    gridspec_kw={'wspace': 0.34},
)

colors = ['C0', 'C1', 'C3']

# --- (a) Plume-core AOD vs time ---
ax_a.plot(days, peaks[fid['case']], marker='s', ms=3.2, lw=1.1,
          color='0.45', ls='--', label=fid['label'])
for e, c in zip(axis, colors):
    ax_a.plot(days, peaks[e['case']], marker='o', ms=3.2, lw=1.3, color=c,
              label=e['label'])

ax_a.set_xscale('log')
ax_a.set_xlabel('Days since injection')
ax_a.set_ylabel('Plume-core AOD at 550 nm')
ax_a.set_ylim(0, 1.0)
ax_a.yaxis.set_minor_locator(ticker.AutoMinorLocator())
ax_a.legend(fontsize=6, frameon=False, loc='upper right', handlelength=1.8,
            labelspacing=0.25, borderaxespad=0.3)
ax_a.text(0.04, 0.95, '(a)', transform=ax_a.transAxes, fontsize=9,
          fontweight='bold', va='top', ha='left')

# --- (b) Plume core, normalized to the 1.0 Tg case ---
# Dashed guides mark exact proportionality; the curves falling below them with
# plume age is the finding.
for e, c, expect in zip(axis[1:], colors[1:], [1.5, 2.0]):
    ax_b.axhline(expect, color=c, lw=0.8, ls=':', zorder=1)
    ax_b.plot(days, peaks[e['case']] / ref_peak, marker='o', ms=3.2, lw=1.3,
              color=c, zorder=3, label=f"{e['label']} / 1.0 Tg")

ax_b.set_xscale('log')
ax_b.set_xlabel('Days since injection')
ax_b.set_ylabel('Plume-core AOD ratio')
ax_b.set_ylim(0.9, 2.3)
ax_b.yaxis.set_minor_locator(ticker.AutoMinorLocator())
ax_b.legend(fontsize=6, frameon=False, loc='lower left', handlelength=1.8,
            labelspacing=0.25, borderaxespad=0.3)
ax_b.text(0.04, 0.95, '(b)', transform=ax_b.transAxes, fontsize=9,
          fontweight='bold', va='top', ha='left')

# --- (c) Global mean, normalized identically ---
for e, c, expect in zip(axis[1:], colors[1:], [1.5, 2.0]):
    ax_c.axhline(expect, color=c, lw=0.8, ls=':', zorder=1)
    ax_c.plot(days, gmeans[e['case']] / ref_gmean, marker='o', ms=3.2, lw=1.3,
              color=c, zorder=3, label=f"{e['label']} / 1.0 Tg")

ax_c.set_xscale('log')
ax_c.set_xlabel('Days since injection')
ax_c.set_ylabel('Global-mean AOD ratio')
ax_c.set_ylim(0.9, 2.3)
ax_c.yaxis.set_minor_locator(ticker.AutoMinorLocator())
ax_c.legend(fontsize=6, frameon=False, loc='lower left', handlelength=1.8,
            labelspacing=0.25, borderaxespad=0.3)
ax_c.text(0.04, 0.95, '(c)', transform=ax_c.transAxes, fontsize=9,
          fontweight='bold', va='top', ha='left')

for ax in (ax_a, ax_b, ax_c):
    ax.set_xlim(0.9, 130)
    ax.set_xticks([1, 10, 100])
    ax.set_xticklabels(['1', '10', '100'])

for ext in ('pdf', 'eps'):
    outfile = os.path.join(here, f"{cfg['outfile_stem']}.{ext}")
    plt.savefig(outfile, dpi=300, bbox_inches='tight')
    print(f"Saved: {outfile}")
plt.close()
