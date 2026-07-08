"""
plot_aod_timeseries.py - AOD vs time for selected cases.

Reads:  remote_analysis/<case>/data/<case>/aod/aod_550nm_band.csv
Plots:  all cases on one axes, x = time in years, y = AOD at 550 nm.

Configuration (cases, paths, output) is in config_aod_timeseries.yaml.
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import argparse
import yaml
import pandas as pd
import matplotlib.pyplot as plt
from pub_data import find_csvs_list

here = os.path.dirname(os.path.abspath(__file__))

parser = argparse.ArgumentParser()
parser.add_argument('--config', default='config_aod_timeseries.yaml',
                    help='YAML config file (relative to script directory)')
args = parser.parse_args()

cfg = yaml.safe_load(open(os.path.join(here, args.config)))

cases = find_csvs_list(cfg['base_dir'], cfg['cases'], cfg['subpath'])
if not cases:
    raise SystemExit("No matching CSV files found.")

fig, ax = plt.subplots(figsize=(10, 4))

for case_name, csv_path in cases.items():
    df   = pd.read_csv(csv_path, header=0)
    days = df.iloc[:, 0].values
    aod  = df.iloc[:, 1].values
    short_name = case_name.removeprefix('exovolc_tambora').lstrip('_')
    print(f"{short_name} peak AOD = {aod.max():.4f}")
    ax.plot(days / cfg['days_per_year'], aod, lw=1.5, label=short_name)

ax.set_xlabel('Time (years)')
ax.set_ylabel('AOD at 550 nm')
ax.legend(frameon=False, fontsize=7, ncols=2, labelspacing=0.3)
plt.tight_layout()

outfile = os.path.join(here, cfg['outfile'])
plt.savefig(outfile, dpi=300, bbox_inches='tight')
plt.close()
print(f"Saved: {outfile}")
