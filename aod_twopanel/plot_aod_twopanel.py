"""
plot_aod_twopanel.py - Two-panel AOD figure for the Tambora case ensemble.

Left panel:  global-mean AOD at 550 nm vs time for every case matching
             `case_glob`, drawn faintly, with the ensemble mean and the
             fiducial case (`case`) highlighted.
Right panel: latitude vs time filled contour of zonal-mean AOD at 550 nm for
             the fiducial case.

Reads:  base_dir/<case>/data/<case>/aod/aod_0p550um_mie.csv        (time series)
        base_dir/<case>/data/<case>/aod/aod_zonal_0p550um_mie.csv  (zonal)

Configuration (case, paths, output) is in config_aod_twopanel.yaml.
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

here = os.path.dirname(os.path.abspath(__file__))

parser = argparse.ArgumentParser()
parser.add_argument('--config', default='config_aod_twopanel.yaml',
                    help='YAML config file (relative to script directory)')
args = parser.parse_args()

cfg = yaml.safe_load(open(os.path.join(here, args.config)))

base_dir = cfg['base_dir']
case = cfg['case']                    # fiducial case
data_dir = os.path.join(base_dir, case, 'data', case)
zonal_path = os.path.join(data_dir, cfg['subpath_zonal'])

if not os.path.exists(zonal_path):
    raise SystemExit(f"Missing zonal CSV: {zonal_path}")

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
# t=0 corresponds to the eruption date; days_since_start map to real calendar dates.
t0 = pd.Timestamp(cfg['eruption_date'])

# --- Ensemble of global-mean time series (left panel) ---
# Every case dir matching case_glob that has the time-series CSV.
case_dirs = sorted(
    d for d in glob.glob(os.path.join(base_dir, cfg['case_glob']))
    if os.path.isdir(d)
)
ts_series = {}   # case_name -> pd.Series indexed by days_since_start
for d in case_dirs:
    c = os.path.basename(d)
    p = os.path.join(base_dir, c, 'data', c, cfg['subpath_timeseries'])
    if not os.path.exists(p):
        print(f"  (skip, no time series) {c}")
        continue
    df = pd.read_csv(p, header=0)
    ts_series[c] = pd.Series(df.iloc[:, 1].values, index=df.iloc[:, 0].values)

if case not in ts_series:
    raise SystemExit(f"Fiducial case {case} has no time series CSV")

# Align every case on the union of day indices and take the ensemble mean.
ts_frame = pd.DataFrame(ts_series).sort_index()
ts_days  = ts_frame.index.values.astype(float)
ts_dates = t0 + pd.to_timedelta(ts_days, unit='D')
ts_mean  = ts_frame.mean(axis=1).values
ts_fid   = ts_frame[case].values

zon = pd.read_csv(zonal_path, header=0, index_col=0)
zon_days  = zon.index.values.astype(float)
zon_dates = t0 + pd.to_timedelta(zon_days, unit='D')
lats      = zon.columns.values.astype(float)
zon_aod   = zon.values.astype(float)          # shape (n_times, n_lats)

# Force all pre-eruption data (days_since_start < 0) strictly to zero.
zon_aod[zon_days < 0, :] = 0.0

# Calendar x-axis range and tick configuration (shared by both panels).
xlim = (pd.Timestamp(cfg['xlim'][0]), pd.Timestamp(cfg['xlim'][1]))

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

print(f"fiducial case = {case}")
print(f"  ensemble cases       = {ts_frame.shape[1]}")
print(f"  ensemble peak AOD    = {ts_frame.max().max():.4f}")
print(f"  fiducial peak AOD    = {np.nanmax(ts_fid):.4f}")
print(f"  zonal peak AOD       = {zon_aod.max():.4f}")

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
fig, (ax_ts, ax_zon) = plt.subplots(
    1, 2, figsize=(10, 3.8),
    gridspec_kw={'width_ratios': [1, 1.15], 'wspace': 0.28},
)

# --- Left: global-mean AOD time series ensemble ---
# Faint grey line per case for the ensemble spread.
for c in ts_frame.columns:
    if c == case:
        continue
    ax_ts.plot(ts_dates, ts_frame[c].values, lw=0.6, color='0.75', zorder=1)
# Highlight the ensemble mean and the fiducial case.
ax_ts.plot(ts_dates, ts_mean, lw=2.0, color='k', zorder=3,
           label='Ensemble mean')
ax_ts.plot(ts_dates, ts_fid, lw=2.0, color='C3', zorder=4,
           label='Fiducial (k35d_r0.5)')
ax_ts.legend(fontsize=7, frameon=False, loc='upper right')
ax_ts.set_xlabel('Year')
ax_ts.set_ylabel('AOD at 550 nm')
ax_ts.set_xlim(*xlim)
ax_ts.set_ylim(bottom=0)
ax_ts.yaxis.set_minor_locator(ticker.AutoMinorLocator())
ax_ts.set_title('(a) Global-mean AOD', fontsize=9, loc='left')

# --- Right: zonal AOD contour ---
# Linear color scale from 0 to 1.4. Shading only, no contour lines.
levels = np.linspace(0, 1.4, 29)

cf = ax_zon.contourf(
    zon_dates, lats, zon_aod.T,
    levels=levels,
    cmap='magma_r',
    extend='max',
)
# Remove the thin antialiased band edges so no faint lines show between levels.
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
ax_zon.set_title('(b) Zonal-mean AOD', fontsize=9, loc='left')

# Calendar ticks: major = each year (labelled), minor = each month.
for ax in (ax_ts, ax_zon):
    ax.xaxis.set_major_locator(mdates.YearLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    ax.xaxis.set_minor_locator(mdates.MonthLocator())

outfile = os.path.join(here, cfg['outfile'])
plt.savefig(outfile, dpi=300, bbox_inches='tight')
plt.close()
print(f"Saved: {outfile}")
