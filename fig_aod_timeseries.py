"""
fig_aod_timeseries.py - AOD vs time for all matching cases.

Reads:  remote_analysis/exovolc_tambora_*/data/exovolc_tambora_*/aod/aod_550nm_band.csv
Plots:  all cases on one axes, x = time in years, y = AOD at 550 nm.

Edit the constants at the top to target different cases or variables.
"""

import pandas as pd
import matplotlib.pyplot as plt
from pub_data import find_csvs

# ---------------------------------------------------------------------------
# Configure here
# ---------------------------------------------------------------------------
BASE_DIR   = '/Users/wolfe/Desktop/projects/volcanos/remote_analysis'
CASE_GLOB  = 'exovolc_tambora_*'
SUBPATH    = 'aod/aod_550nm_band.csv'
DAYS_PER_YEAR = 365.0
OUTFILE    = 'fig_aod_timeseries.pdf'
EXCLUDE = {'k40d'}
# ---------------------------------------------------------------------------

cases = find_csvs(BASE_DIR, CASE_GLOB, SUBPATH)
if not cases:
    raise SystemExit("No matching CSV files found.")

fig, ax = plt.subplots(figsize=(7, 4))


for case_name, csv_path in cases.items():
    df   = pd.read_csv(csv_path, header=0)
    days = df.iloc[:, 0].values
    aod  = df.iloc[:, 1].values
    short_name = case_name.removeprefix('exovolc_tambora').lstrip('_')
    if short_name in EXCLUDE:
        continue
    print(f"{short_name} peak AOD ={aod.max():.4f}")
    ax.plot(days / DAYS_PER_YEAR, aod, lw=1.5, label=short_name)

ax.set_xlabel('Time (years)')
ax.set_ylabel('AOD at 550 nm')
ax.legend(frameon=False, fontsize=7, ncols=2, labelspacing=0.3)
plt.tight_layout()
plt.savefig(OUTFILE, dpi=300, bbox_inches='tight')
plt.close()
print(f"Saved: {OUTFILE}")
