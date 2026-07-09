"""
plot_sulfur_burden.py - Global sulfur burden (Tg S) for the Tambora ensemble.

Plots the H2SO4 aerosol (VOLCHZ) sulfur burden in Tg of sulfur atoms only, to
match the burden figure in Clyne et al. Every case matching `case_glob` is
drawn as a faint line; the fiducial case (`case`) is overplotted in bold.

For the fiducial only, the SO2 and H2SO4 gas reservoirs are also read and
reported to stdout (converted to Tg S) as a sulfur-conservation sanity check;
only the VOLCHZ aerosol curve is plotted. Mass fractions:
    SO2       gas                  (MW 64  -> S fraction 32/64)
    H2SO4     gas                  (MW 98  -> S fraction 32/98)
    VOLCHZMD  H2SO4 aerosol        (MW 98  -> S fraction 32/98)

Each scalar CSV holds a global integrated mass time series in kg. (The CSV
column headers list the wrong units — treat the values as kg.) Only the mass of
the sulfur atoms is counted, per the mass-fraction factors in the YAML.

Configuration (case, paths, output, S fractions) is in config_sulfur_burden.yaml.
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
parser.add_argument('--config', default='config_sulfur_burden.yaml',
                    help='YAML config file (relative to script directory)')
args = parser.parse_args()

cfg = yaml.safe_load(open(os.path.join(here, args.config)))

base_dir = cfg['base_dir']
case = cfg['case']                    # fiducial case

KG_PER_TG = 1.0e9   # 1 Tg = 1e9 kg


def load_kg(case_name, subpath):
    """Load a scalar CSV; return (days, mass_kg). Values are kg despite header."""
    p = os.path.join(base_dir, case_name, 'data', case_name, subpath)
    if not os.path.exists(p):
        raise SystemExit(f"Missing CSV: {p}")
    df = pd.read_csv(p, header=0)
    return df.iloc[:, 0].values.astype(float), df.iloc[:, 1].values.astype(float)


def volc_burden_tgs(case_name):
    """VOLCHZ (H2SO4 aerosol) sulfur burden in Tg S for one case."""
    days, volc_kg = load_kg(case_name, cfg['subpath_volc'])
    return days, volc_kg * cfg['s_frac_volc'] / KG_PER_TG

# ---------------------------------------------------------------------------
# Load data and convert to Tg S
# ---------------------------------------------------------------------------
# t=0 corresponds to the eruption date; days_since_start map to calendar dates.
t0   = pd.Timestamp(cfg['eruption_date'])
xlim = (pd.Timestamp(cfg['xlim'][0]), pd.Timestamp(cfg['xlim'][1]))

# Ensemble VOLCHZ burden: every case dir matching case_glob with the CSV.
case_dirs = sorted(
    d for d in glob.glob(os.path.join(base_dir, cfg['case_glob']))
    if os.path.isdir(d)
)
volc_series = {}   # case_name -> pd.Series of Tg S indexed by days
for d in case_dirs:
    c = os.path.basename(d)
    p = os.path.join(base_dir, c, 'data', c, cfg['subpath_volc'])
    if not os.path.exists(p):
        print(f"  (skip, no VOLCHZMD) {c}")
        continue
    days, burden = volc_burden_tgs(c)
    volc_series[c] = pd.Series(burden, index=days)

if case not in volc_series:
    raise SystemExit(f"Fiducial case {case} has no VOLCHZMD CSV")

# Fiducial reservoirs (for the conservation sanity check printed to stdout).
days,  so2_kg   = load_kg(case, cfg['subpath_so2'])
_,     h2so4_kg = load_kg(case, cfg['subpath_h2so4'])
_,     volc_kg  = load_kg(case, cfg['subpath_volc'])
so2_s   = so2_kg   * cfg['s_frac_so2']   / KG_PER_TG
h2so4_s = h2so4_kg * cfg['s_frac_h2so4'] / KG_PER_TG
volc_s  = volc_kg  * cfg['s_frac_volc']  / KG_PER_TG
total_s = so2_s + h2so4_s + volc_s

print(f"fiducial case = {case}")
print(f"  ensemble cases      = {len(volc_series)}")
print(f"  peak SO2   S burden = {so2_s.max():.3f} Tg S")
print(f"  peak H2SO4 S burden = {h2so4_s.max():.4f} Tg S")
print(f"  peak VOLC  S burden = {volc_s.max():.3f} Tg S")
print(f"  peak TOTAL S burden = {total_s.max():.3f} Tg S")

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(6, 3.8))

# Faint grey line per case for the ensemble spread.
for c, s in volc_series.items():
    if c == case:
        continue
    ax.plot(t0 + pd.to_timedelta(s.index.values, unit='D'), s.values,
            lw=0.6, color='0.75', zorder=1)
# Fiducial VOLCHZ sulfur burden in bold, to match Clyne et al.
fid = volc_series[case]
ax.plot(t0 + pd.to_timedelta(fid.index.values, unit='D'), fid.values,
        lw=1.8, color='k', zorder=3, label='Fiducial (k35d_r0.5)')

ax.set_xlabel('Year')
ax.set_ylabel('Global sulfur burden (Tg S)')
ax.set_xlim(*xlim)
ax.set_ylim(bottom=0)
ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
ax.legend(fontsize=8, frameon=False, loc='upper right')

# Calendar ticks: major = each year (labelled), minor = each month.
ax.xaxis.set_major_locator(mdates.YearLocator())
ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
ax.xaxis.set_minor_locator(mdates.MonthLocator())

outfile = os.path.join(here, cfg['outfile'])
plt.savefig(outfile, dpi=300, bbox_inches='tight')
plt.close()
print(f"Saved: {outfile}")
