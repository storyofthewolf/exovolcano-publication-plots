"""
plot_aod_twopanel.py - Two-panel AOD figure for a single case.

Left panel:  global-mean AOD at 550 nm vs time.
Right panel: latitude vs time filled contour of zonal-mean AOD at 550 nm.

Reads:  base_dir/<case>/data/<case>/aod/aod_0p550um_mie.csv        (time series)
        base_dir/<case>/data/<case>/aod/aod_zonal_0p550um_mie.csv  (zonal)

Configuration (case, paths, output) is in config_aod_twopanel.yaml.
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

here = os.path.dirname(os.path.abspath(__file__))

parser = argparse.ArgumentParser()
parser.add_argument('--config', default='config_aod_twopanel.yaml',
                    help='YAML config file (relative to script directory)')
args = parser.parse_args()

cfg = yaml.safe_load(open(os.path.join(here, args.config)))

case = cfg['case']
data_dir = os.path.join(cfg['base_dir'], case, 'data', case)
ts_path = os.path.join(data_dir, cfg['subpath_timeseries'])
zonal_path = os.path.join(data_dir, cfg['subpath_zonal'])

for p in (ts_path, zonal_path):
    if not os.path.exists(p):
        raise SystemExit(f"Missing CSV: {p}")

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
# t=0 corresponds to the eruption date; days_since_start map to real calendar dates.
t0 = pd.Timestamp(cfg['eruption_date'])

ts = pd.read_csv(ts_path, header=0)
ts_dates = t0 + pd.to_timedelta(ts.iloc[:, 0].values, unit='D')
ts_aod   = ts.iloc[:, 1].values

zon = pd.read_csv(zonal_path, header=0, index_col=0)
zon_dates = t0 + pd.to_timedelta(zon.index.values.astype(float), unit='D')
lats      = zon.columns.values.astype(float)
zon_aod   = zon.values.astype(float)          # shape (n_times, n_lats)

# Calendar x-axis range and tick configuration (shared by both panels).
xlim = (pd.Timestamp(cfg['xlim'][0]), pd.Timestamp(cfg['xlim'][1]))

print(f"case = {case}")
print(f"  time-series peak AOD = {ts_aod.max():.4f}")
print(f"  zonal peak AOD       = {zon_aod.max():.4f}")

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
fig, (ax_ts, ax_zon) = plt.subplots(
    1, 2, figsize=(10, 3.8),
    gridspec_kw={'width_ratios': [1, 1.15], 'wspace': 0.28},
)

# --- Left: global-mean AOD time series ---
ax_ts.plot(ts_dates, ts_aod, lw=1.5, color='C3')
ax_ts.set_xlabel('Year')
ax_ts.set_ylabel('AOD at 550 nm')
ax_ts.set_xlim(*xlim)
ax_ts.set_ylim(bottom=0)
ax_ts.yaxis.set_minor_locator(ticker.AutoMinorLocator())
ax_ts.set_title('(a) Global-mean AOD', fontsize=9, loc='left')

# --- Right: zonal AOD contour ---
levels = np.concatenate([
    np.linspace(0, 0.1, 6),
    np.linspace(0.2, 1.0, 5),
])

cf = ax_zon.contourf(
    zon_dates, lats, zon_aod.T,
    levels=levels,
    cmap='magma_r',
    extend='max',
)
ax_zon.contour(
    zon_dates, lats, zon_aod.T,
    levels=levels,
    colors='k',
    linewidths=0.3,
    alpha=0.4,
)

cbar = fig.colorbar(cf, ax=ax_zon, pad=0.02, aspect=30)
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
