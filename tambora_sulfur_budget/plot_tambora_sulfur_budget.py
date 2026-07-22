"""
plot_tambora_sulfur_budget.py - Figure 2: sulfur/aerosol mass budget for the
fiducial Tambora case.

Top panel:    SO2 gas, H2SO4 gas, and condensed sulfate aerosol (VOLCHZMD)
              burdens vs time, each converted to Tg of sulfur (Tg S) so the
              three reservoirs sum meaningfully. The total (sum of all three)
              is overplotted in black, with the injected 30 Tg S total
              annotated as a horizontal reference line.
Bottom panel: global-mean AOD at 550 nm vs time for the same case, sharing
              the x-axis, so the optical signature lines up in time with the
              mass budget above it.

Reads (current on-disk layout):
    base_dir/<eruption>/<case_name>/data/scalar/SO2.csv
    base_dir/<eruption>/<case_name>/data/scalar/H2SO4.csv
    base_dir/<eruption>/<case_name>/data/scalar/VOLCHZMD.csv
    base_dir/<eruption>/<case_name>/data/aod/aod_550nm_band.csv

Each scalar CSV holds a global integrated mass time series in kg (the column
header units are mislabeled as e.g. "SO2 [kg]", which is correct -- despite
similar issues in other figures' source data, these are kg as labeled).

Configuration (case, paths, molar masses, output) is in
config_tambora_sulfur_budget.yaml.
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
parser.add_argument('--config', default='config_tambora_sulfur_budget.yaml',
                    help='YAML config file (relative to script directory)')
args = parser.parse_args()

cfg = yaml.safe_load(open(os.path.join(here, args.config)))

base_dir = cfg['base_dir']
eruption = cfg['eruption']
case = cfg['case']
data_dir = os.path.join(base_dir, eruption, case, 'data')

KG_PER_TG = 1.0e9   # 1 Tg = 1e9 kg

s_frac_so2   = cfg['mw_s'] / cfg['mw_so2']
s_frac_h2so4 = cfg['mw_s'] / cfg['mw_h2so4']
s_frac_volc  = cfg['mw_s'] / cfg['mw_volc']


def load_kg(subpath):
    p = os.path.join(data_dir, subpath)
    if not os.path.exists(p):
        raise SystemExit(f"Missing CSV: {p}")
    df = pd.read_csv(p, header=0)
    return df.iloc[:, 0].values.astype(float), df.iloc[:, 1].values.astype(float)

# ---------------------------------------------------------------------------
# Load data and convert to Tg S
# ---------------------------------------------------------------------------
t0   = pd.Timestamp(cfg['eruption_date'])
xlim = (pd.Timestamp(cfg['xlim'][0]), pd.Timestamp(cfg['xlim'][1]))

days, so2_kg   = load_kg(cfg['subpath_so2'])
_,    h2so4_kg = load_kg(cfg['subpath_h2so4'])
_,    volc_kg  = load_kg(cfg['subpath_volc'])
_,    aod      = load_kg(cfg['subpath_aod'])   # not a mass, but same 2-col shape

so2_s   = so2_kg   * s_frac_so2   / KG_PER_TG
h2so4_s = h2so4_kg * s_frac_h2so4 / KG_PER_TG
volc_s  = volc_kg  * s_frac_volc  / KG_PER_TG
total_s = so2_s + h2so4_s + volc_s

dates = t0 + pd.to_timedelta(days, unit='D')

# ---------------------------------------------------------------------------
# Diagnostics
# ---------------------------------------------------------------------------
so2_peak_idx = np.nanargmax(so2_kg)
volc_peak_idx = np.nanargmax(volc_kg)
tot_peak_idx = np.nanargmax(total_s)

print(f"fiducial case = {case}")
print(f"  peak SO2 burden      = {so2_kg.max()/KG_PER_TG:.3f} Tg "
      f"({so2_s.max():.3f} Tg S) at day {days[so2_peak_idx]:.0f}")
print(f"  peak H2SO4 S burden  = {h2so4_s.max():.4f} Tg S")
print(f"  peak VOLCHZMD burden = {volc_kg.max()/KG_PER_TG:.3f} Tg "
      f"({volc_s.max():.3f} Tg S) at day {days[volc_peak_idx]:.0f}")
print(f"  peak TOTAL S burden  = {total_s.max():.3f} Tg S at day {days[tot_peak_idx]:.0f}")
print(f"  injected total S     = {cfg['injected_total_s_tg']:.1f} Tg S")

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
# Linear y-axis chosen for the mass-budget panel: the three reservoirs and
# their sum are all within a factor of a few of one another over the active
# period, and a linear scale makes the total-S conservation/exchange between
# SO2 -> H2SO4 -> aerosol visually intuitive (areas trade off against each
# other). A log axis would compress that trade-off and obscure the intended
# reading -- that the black total curve tracks the injected 30 TgS line.
fig, (ax_s, ax_aod) = plt.subplots(
    2, 1, figsize=(3.5, 5.2), sharex=True,
    gridspec_kw={'height_ratios': [1.3, 1], 'hspace': 0.08},
)

color_so2   = 'C0'
color_h2so4 = 'C1'
color_volc  = 'C2'

ax_s.plot(dates, so2_s,   lw=1.4, color=color_so2,   label=r'SO$_2$ (gas)')
ax_s.plot(dates, h2so4_s, lw=1.4, color=color_h2so4, label=r'H$_2$SO$_4$ (gas)')
ax_s.plot(dates, volc_s,  lw=1.4, color=color_volc,  label='Sulfate aerosol')
ax_s.plot(dates, total_s, lw=1.8, color='k', label='Total S', zorder=5)

ax_s.axhline(cfg['injected_total_s_tg'], color='0.4', lw=0.9, ls=':', zorder=1)
# Label sits right-aligned at the end of the dotted line, in the clean space
# above it (after the curves have decayed well below 30 TgS by ~1818-1819),
# clear of the (a) tag, the legend, and the total-S curve.
ax_s.text(0.97, cfg['injected_total_s_tg'],
          f"Injected: {cfg['injected_total_s_tg']:.0f} Tg S",
          transform=ax_s.get_yaxis_transform(),
          fontsize=6.5, color='0.3', va='bottom', ha='right')

ax_s.set_ylabel('Sulfur burden (Tg S)')
ax_s.set_ylim(0, 33)
ax_s.yaxis.set_minor_locator(ticker.AutoMinorLocator())
ax_s.legend(fontsize=6.5, frameon=False, loc='center right')
ax_s.text(0.02, 0.96, '(a)', transform=ax_s.transAxes,
          fontsize=9, fontweight='bold', va='top', ha='left')

ax_aod.plot(dates, aod, lw=1.6, color='C3')
ax_aod.set_ylabel('AOD at 550 nm')
ax_aod.set_xlabel('Year')
ax_aod.set_ylim(bottom=0)
ax_aod.yaxis.set_minor_locator(ticker.AutoMinorLocator())
ax_aod.text(0.02, 0.94, '(b)', transform=ax_aod.transAxes,
            fontsize=9, fontweight='bold', va='top', ha='left')

ax_aod.set_xlim(*xlim)
ax_aod.xaxis.set_major_locator(mdates.YearLocator())
ax_aod.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
ax_aod.xaxis.set_minor_locator(mdates.MonthLocator())

outfile_pdf = os.path.join(here, cfg['outfile_base'] + '.pdf')
outfile_eps = os.path.join(here, cfg['outfile_base'] + '.eps')
plt.savefig(outfile_pdf, dpi=300, bbox_inches='tight')
plt.savefig(outfile_eps, dpi=300, bbox_inches='tight')
plt.close()
print(f"Saved: {outfile_pdf}")
print(f"Saved: {outfile_eps}")
