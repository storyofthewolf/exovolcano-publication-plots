"""
plot_tambora_validation.py - Figure 1: headline Tambora VolMIP validation figure.

Two panels:
(a) Global-mean AOD at 550 nm vs time for the fiducial case, with the other
    15 sweep cases drawn faintly behind it as a sensitivity envelope. A
    shaded band marks the VolMIP-Tambora inter-model peak SAOD range
    (Clyne et al. 2021), with a dashed line at the VolMIP ensemble mean.
(b) Latitude-time Hovmoller contour of zonal-mean AOD at 550 nm for the
    fiducial case.

Reads (current on-disk layout):
    base_dir/<eruption>/<case_name>/data/aod/aod_550nm_band.csv
    base_dir/<eruption>/<case_name>/data/aod/aod_zonal_550nm_band.csv

Configuration (cases, paths, VolMIP band, output) is in
config_tambora_validation.yaml.
"""

import os
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
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'legend.fontsize': 7,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})

here = os.path.dirname(os.path.abspath(__file__))

parser = argparse.ArgumentParser()
parser.add_argument('--config', default='config_tambora_validation.yaml',
                    help='YAML config file (relative to script directory)')
args = parser.parse_args()

cfg = yaml.safe_load(open(os.path.join(here, args.config)))

base_dir = cfg['base_dir']
eruption = cfg['eruption']
case = cfg['case']

def case_dir(c):
    return os.path.join(base_dir, eruption, c, 'data')

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
t0 = pd.Timestamp(cfg['eruption_date'])
xlim = (pd.Timestamp(cfg['xlim'][0]), pd.Timestamp(cfg['xlim'][1]))

# --- Panel (a): global-mean AOD time series, all 16 cases ---
ts_series = {}   # case_name -> pd.Series indexed by days_since_start
for c in cfg['cases']:
    p = os.path.join(case_dir(c), cfg['subpath_timeseries'])
    if not os.path.exists(p):
        print(f"  (skip, no time series) {c}: {p}")
        continue
    df = pd.read_csv(p, header=0)
    ts_series[c] = pd.Series(df.iloc[:, 1].values, index=df.iloc[:, 0].values)

if case not in ts_series:
    raise SystemExit(f"Fiducial case {case} has no time series CSV")

ts_frame = pd.DataFrame(ts_series).sort_index()
ts_days  = ts_frame.index.values.astype(float)
ts_dates = t0 + pd.to_timedelta(ts_days, unit='D')
ts_fid   = ts_frame[case].values

# --- Panel (b): zonal-mean AOD contour, fiducial only ---
zonal_path = os.path.join(case_dir(case), cfg['subpath_zonal'])
if not os.path.exists(zonal_path):
    raise SystemExit(f"Missing zonal CSV: {zonal_path}")

zon = pd.read_csv(zonal_path, header=0, index_col=0)
zon_days  = zon.index.values.astype(float)
zon_dates = t0 + pd.to_timedelta(zon_days, unit='D')
lats      = zon.columns.values.astype(float)
zon_aod   = zon.values.astype(float)          # shape (n_times, n_lats)

# Force all pre-eruption data (days_since_start < 0) strictly to zero.
zon_aod[zon_days < 0, :] = 0.0

# Pad the zonal grid with 0.0 out to the axis limits so the contour fills the
# pre-eruption (and any post-data) span with zeros instead of leaving it blank.
pad_dates = [d for d in (xlim[0], xlim[1])
             if d < zon_dates.min() or d > zon_dates.max()]
if pad_dates:
    pad_dates = np.array(pad_dates, dtype='datetime64[ns]')
    pad_rows  = np.zeros((len(pad_dates), zon_aod.shape[1]))
    all_dates = np.concatenate([zon_dates.values, pad_dates])
    all_aod   = np.concatenate([zon_aod, pad_rows], axis=0)
    order     = np.argsort(all_dates)
    zon_dates = pd.DatetimeIndex(all_dates[order])
    zon_aod   = all_aod[order]

# ---------------------------------------------------------------------------
# Diagnostics
# ---------------------------------------------------------------------------
fid_peak_idx = np.nanargmax(ts_fid)
fid_peak_aod = ts_fid[fid_peak_idx]
fid_peak_day = ts_days[fid_peak_idx]

zon_peak_flat = np.nanargmax(zon_aod)
zi, zj = np.unravel_index(zon_peak_flat, zon_aod.shape)
zon_peak_aod = zon_aod[zi, zj]
zon_peak_day = zon_days[zi]
zon_peak_lat = lats[zj]

env_peaks = ts_frame.max(axis=0)

print(f"fiducial case = {case}")
print(f"  n cases in envelope     = {ts_frame.shape[1]}")
print(f"  fiducial peak AOD       = {fid_peak_aod:.4f} at day {fid_peak_day:.0f}")
print(f"  zonal peak AOD          = {zon_peak_aod:.4f} at day {zon_peak_day:.0f}, "
      f"lat {zon_peak_lat:.0f}")
print(f"  16-case envelope peak AOD: min = {env_peaks.min():.4f} "
      f"({env_peaks.idxmin()}), max = {env_peaks.max():.4f} ({env_peaks.idxmax()})")

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
fig, (ax_ts, ax_zon) = plt.subplots(
    1, 2, figsize=(7.1, 3.0),
    gridspec_kw={'width_ratios': [1, 1.15], 'wspace': 0.30},
)

# --- Panel (a): global-mean AOD, fiducial + envelope + VolMIP band ---
# VolMIP-Tambora inter-model peak SAOD range (Clyne et al. 2021).
# Use a solid pale fill (not alpha) so the EPS backend renders it without a
# transparency warning/fallback.
ax_ts.axhspan(cfg['volmip_saod_min'], cfg['volmip_saod_max'],
              color='#d6e6f5', zorder=0,
              label='VolMIP-Tambora peak SAOD range')
ax_ts.axhline(cfg['volmip_saod_mean'], color='C0', lw=1.0, ls='--', zorder=2,
              label='VolMIP ensemble mean')

# Faint grey line per non-fiducial case for the sensitivity envelope.
for c in ts_frame.columns:
    if c == case:
        continue
    ax_ts.plot(ts_dates, ts_frame[c].values, lw=0.6, color='0.75', zorder=1)
# One representative grey line for the legend entry (proxy handle).
ax_ts.plot([], [], lw=0.6, color='0.75', label='Sensitivity sweeps (15 cases)')

ax_ts.plot(ts_dates, ts_fid, lw=1.8, color='C3', zorder=4,
           label='Fiducial')

ax_ts.legend(fontsize=6.5, frameon=False, loc='upper right')
ax_ts.set_xlabel('Year')
ax_ts.set_ylabel('AOD at 550 nm')
ax_ts.set_xlim(*xlim)
ax_ts.set_ylim(bottom=0)
ax_ts.yaxis.set_minor_locator(ticker.AutoMinorLocator())
ax_ts.text(0.02, 0.96, '(a)', transform=ax_ts.transAxes,
           fontsize=9, fontweight='bold', va='top', ha='left')

# --- Panel (b): zonal AOD contour, fiducial only ---
levels = np.linspace(0, 1.4, 29)

cf = ax_zon.contourf(
    zon_dates, lats, zon_aod.T,
    levels=levels,
    cmap='magma_r',
    extend='max',
)
cf.set_edgecolor('face')
cf.set_linewidth(0)

cbar = fig.colorbar(cf, ax=ax_zon, pad=0.02, aspect=30)
cbar.set_ticks(np.linspace(0, 1.4, 8))
cbar.set_label('AOD at 550 nm')
cbar.ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())

ax_zon.set_xlabel('Year')
ax_zon.set_ylabel('Latitude (°)')
ax_zon.set_xlim(*xlim)
ax_zon.set_yticks(np.arange(-90, 91, 30))
ax_zon.yaxis.set_minor_locator(ticker.MultipleLocator(10))
ax_zon.text(0.02, 0.96, '(b)', transform=ax_zon.transAxes,
            fontsize=9, fontweight='bold', va='top', ha='left', color='black')

# Calendar ticks: major = each year (labelled), minor = each month.
for ax in (ax_ts, ax_zon):
    ax.xaxis.set_major_locator(mdates.YearLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    ax.xaxis.set_minor_locator(mdates.MonthLocator())

outfile_pdf = os.path.join(here, cfg['outfile_base'] + '.pdf')
outfile_eps = os.path.join(here, cfg['outfile_base'] + '.eps')
plt.savefig(outfile_pdf, dpi=300, bbox_inches='tight')
plt.savefig(outfile_eps, dpi=300, bbox_inches='tight')
plt.close()
print(f"Saved: {outfile_pdf}")
print(f"Saved: {outfile_eps}")
