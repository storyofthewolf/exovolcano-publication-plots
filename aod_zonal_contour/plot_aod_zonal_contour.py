"""
plot_aod_zonal_contour.py - Ensemble-mean zonal AOD contour plot.

Reads:  remote_analysis/<case>/data/<case>/aod/aod_zonal_550nm_band.csv
Plots:  latitude vs time filled contour of AOD at 550 nm, averaged across all
        listed cases.

Configuration (cases, paths, output) is in config_aod_zonal_contour.yaml.
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
from pub_data import find_csvs_list

here = os.path.dirname(os.path.abspath(__file__))

parser = argparse.ArgumentParser()
parser.add_argument('--config', default='config_aod_zonal_contour.yaml',
                    help='YAML config file (relative to script directory)')
args = parser.parse_args()

cfg = yaml.safe_load(open(os.path.join(here, args.config)))

cases = find_csvs_list(cfg['base_dir'], cfg['cases'], cfg['subpath'])
if not cases:
    raise SystemExit("No matching CSV files found.")

# Accumulate arrays for averaging
stack = []
for case_name, csv_path in cases.items():
    short_name = case_name.removeprefix('exovolc_tambora').lstrip('_')
    df = pd.read_csv(csv_path, header=0, index_col=0)
    stack.append(df.values.astype(float))
    print(f"  loaded {short_name}")

# Build coordinate arrays from the last-loaded df (all files share the same grid)
days = df.index.values.astype(float)
lats = df.columns.values.astype(float)
years = days / cfg['days_per_year']

# Ensemble mean
mean_aod = np.mean(stack, axis=0)   # shape (n_times, n_lats)
print(f"\nAveraged {len(stack)} cases. Peak ensemble-mean AOD = {mean_aod.max():.4f}")

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
levels = np.concatenate([
    np.linspace(0, 0.1, 6),
    np.linspace(0.2, 1.0, 5),
])

fig, ax = plt.subplots(figsize=(4, 3.5))

cf = ax.contourf(
    years, lats, mean_aod.T,
    levels=levels,
    cmap='magma_r',
    extend='max',
)
ax.contour(
    years, lats, mean_aod.T,
    levels=levels,
    colors='k',
    linewidths=0.3,
    alpha=0.4,
)

cbar = fig.colorbar(cf, ax=ax, pad=0.02, aspect=30)
cbar.set_label('AOD at 550 nm')
cbar.ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())

ax.set_xlabel('Time (years)')
ax.set_ylabel('Latitude (°)')
ax.set_yticks(np.arange(-90, 91, 30))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(10))
ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
ax.set_title(f'Ensemble-mean zonal AOD (n={len(stack)})', fontsize=9)

plt.tight_layout()

outfile = os.path.join(here, cfg['outfile'])
plt.savefig(outfile, dpi=300, bbox_inches='tight')
plt.close()
print(f"Saved: {outfile}")
