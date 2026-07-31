"""
plot_hunga_schoeberl_compare.py - First-year 550 nm AOD for the Hunga Tonga
ensemble, on the global-mean basis Schoeberl et al. (2023) report.

Two panels:
(a) Global-mean 550 nm AOD versus time over the first 365 days, with the full
    28-case ensemble drawn faintly, the fiducial and the Raikoke-equivalent
    calibration control highlighted, and the Schoeberl et al. (2023)
    global-mean peak of 0.018 (March-April 2022) marked.
(b) Zonal-mean 550 nm AOD as a latitude-time contour over the same 365 days
    for the fiducial case, showing the tropical/SH confinement that makes the
    global mean an ungenerous comparison during the first months.

NOT a replication of any Schoeberl figure -- see the config YAML header. The
only quantity taken from that paper is the verified 0.018 global-mean peak.

Reads:  <base_dir>/<case>/data/aod/aod_550nm_band.csv        (global mean)
        <base_dir>/<case>/data/aod/aod_zonal_550nm_band.csv  (zonal mean)

Configuration is in config_hunga_schoeberl_compare.yaml.
"""

import sys
import os
import glob
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
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'legend.fontsize': 7,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})

here = os.path.dirname(os.path.abspath(__file__))

parser = argparse.ArgumentParser()
parser.add_argument('--config', default='config_hunga_schoeberl_compare.yaml',
                    help='YAML config file (relative to script directory)')
args = parser.parse_args()

cfg = yaml.safe_load(open(os.path.join(here, args.config)))
base_dirs = cfg['base_dirs']
n_days = float(cfg['n_days'])
t0 = pd.Timestamp(cfg['eruption_date'])

# ---------------------------------------------------------------------------
# Load the ensemble (both phases), truncated to the first year
# ---------------------------------------------------------------------------
series = {}     # case -> (days, aod), restricted to days <= n_days
for phase, base in base_dirs.items():
    for d in sorted(glob.glob(os.path.join(base, cfg['case_glob']))):
        if not os.path.isdir(d):
            continue
        case = os.path.basename(d)
        p = os.path.join(d, 'data', cfg['subpath_timeseries'])
        if not os.path.exists(p):
            continue
        df = pd.read_csv(p)
        t, v = df.iloc[:, 0].values.astype(float), df.iloc[:, 1].values.astype(float)
        m = t <= n_days
        # Later phases must not silently overwrite an identically named phase-1
        # case; the two suites are disjoint by construction, so warn if not.
        if case in series:
            print(f"  WARNING: {case} appears in more than one phase; keeping the first")
            continue
        series[case] = (t[m], v[m])

print(f"ensemble cases loaded: {len(series)}")


def get(entry):
    case = entry['case']
    if case not in series:
        raise SystemExit(f"Missing case in ensemble: {case}")
    return series[case]


fid_days, fid_aod = get(cfg['fiducial'])
cal_days, cal_aod = get(cfg['calibration'])

# --- Zonal field for panel (b) ---
zbase = base_dirs[cfg['zonal_phase']]
zpath = os.path.join(zbase, cfg['zonal_case'], 'data', cfg['subpath_zonal'])
if not os.path.exists(zpath):
    raise SystemExit(f"Missing zonal CSV: {zpath}")

zon = pd.read_csv(zpath, index_col=0)
zon_days = zon.index.values.astype(float)
zmask = zon_days <= n_days
zon_days = zon_days[zmask]
lats = zon.columns.values.astype(float)
zon_aod = zon.values.astype(float)[zmask]
zon_aod[zon_days < 0, :] = 0.0

fid_dates = t0 + pd.to_timedelta(fid_days, unit='D')
cal_dates = t0 + pd.to_timedelta(cal_days, unit='D')
zon_dates = t0 + pd.to_timedelta(zon_days, unit='D')

ref = cfg['reference']
xlim = (t0, t0 + pd.Timedelta(days=n_days))

# ---------------------------------------------------------------------------
# Diagnostics
# ---------------------------------------------------------------------------
i_f, i_c = int(np.argmax(fid_aod)), int(np.argmax(cal_aod))
print(f"fiducial       peak {fid_aod[i_f]:.5f} on {fid_dates[i_f]:%Y-%m-%d} "
      f"({fid_aod[i_f]/ref['value']:.2f}x the {ref['value']} reference)")
print(f"raikoke-equiv  peak {cal_aod[i_c]:.5f} on {cal_dates[i_c]:%Y-%m-%d} "
      f"({cal_aod[i_c]/0.016:.2f}x the 0.016 Raikoke value)")
iz = np.unravel_index(np.argmax(zon_aod), zon_aod.shape)
print(f"zonal peak     {zon_aod[iz]:.5f} at lat {lats[iz[1]]:+.0f} "
      f"on {zon_dates[iz[0]]:%Y-%m-%d}")

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
fig, (ax_a, ax_b) = plt.subplots(
    1, 2, figsize=(7.1, 3.0),
    gridspec_kw={'width_ratios': [1, 1.15], 'wspace': 0.30},
)

# --- (a) Global-mean AOD, ensemble + highlights + reference ---
n_other = 0
for case, (t, v) in series.items():
    if case in (cfg['fiducial']['case'], cfg['calibration']['case']):
        continue
    ax_a.plot(t0 + pd.to_timedelta(t, unit='D'), v, lw=0.5, color='0.78', zorder=1)
    n_other += 1
ax_a.plot([], [], lw=1.2, color='0.78', label=f'Ensemble (n={n_other})')

# The reference window is drawn as a solid pale band rather than an alpha
# fill: EPS has no real transparency.
w0, w1 = pd.Timestamp(ref['window'][0]), pd.Timestamp(ref['window'][1])
ax_a.axvspan(w0, w1, color='#DDE8F0', zorder=0, lw=0)
ax_a.axhline(ref['value'], color='C0', lw=1.1, ls='--', zorder=3,
             label=ref['label'])

ax_a.plot(cal_dates, cal_aod, lw=1.5, color='C2', zorder=4,
          label=cfg['calibration']['label'])
ax_a.plot(fid_dates, fid_aod, lw=1.6, color='C3', zorder=5,
          label=cfg['fiducial']['label'])

ax_a.set_xlabel('Date')
ax_a.set_ylabel('Global-mean AOD at 550 nm')
ax_a.set_xlim(*xlim)
ax_a.set_ylim(0, 0.024)
ax_a.yaxis.set_minor_locator(ticker.AutoMinorLocator())
ax_a.text(0.03, 0.95, '(a)', transform=ax_a.transAxes, fontsize=9,
          fontweight='bold', va='top', ha='left')
ax_a.legend(fontsize=6, frameon=False, loc='upper center',
            bbox_to_anchor=(0.5, -0.24), ncol=1, handlelength=2.0,
            labelspacing=0.3, borderaxespad=0.0)

# --- (b) Zonal-mean AOD contour ---
# Round the top of the scale up to a clean value so the colorbar carries
# readable ticks rather than the raw data maximum's decimals.
vmax = float(np.ceil(np.nanmax(zon_aod) * 100.0) / 100.0)
levels = np.linspace(0, vmax, 21)
# vmax is rounded UP from the data maximum, so nothing is clipped and the
# colorbar needs no extend arrow.
cf = ax_b.contourf(zon_dates, lats, zon_aod.T, levels=levels,
                   cmap='magma_r', antialiased=False)
cf.set_edgecolor('face')
cf.set_linewidth(0)

cbar = fig.colorbar(cf, ax=ax_b, pad=0.02, aspect=30)
cbar.set_ticks(np.linspace(0, vmax, 5))
cbar.set_label('Zonal-mean AOD at 550 nm')
cbar.ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())

# Vent latitude, for reading the SH confinement against.
ax_b.axhline(-20.5, color='cyan', lw=0.8, ls=':', zorder=3)

ax_b.set_xlabel('Date')
ax_b.set_ylabel('Latitude (°)')
ax_b.set_xlim(*xlim)
ax_b.set_ylim(-90, 90)
ax_b.set_yticks(np.arange(-90, 91, 30))
ax_b.yaxis.set_minor_locator(ticker.MultipleLocator(10))
ax_b.text(0.03, 0.95, '(b)', transform=ax_b.transAxes, fontsize=9,
          fontweight='bold', va='top', ha='left', color='white')

for ax in (ax_a, ax_b):
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
    ax.xaxis.set_minor_locator(mdates.MonthLocator())

for ext in ('pdf', 'eps'):
    outfile = os.path.join(here, f"{cfg['outfile_stem']}.{ext}")
    plt.savefig(outfile, dpi=300, bbox_inches='tight')
    print(f"Saved: {outfile}")
plt.close()
