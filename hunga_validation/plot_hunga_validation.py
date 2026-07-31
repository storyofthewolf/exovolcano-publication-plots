"""
plot_hunga_validation.py - Headline observational validation figure for the
ExoCAM Hunga Tonga-Hunga Ha'apai (2022) case.

Two panels:
(a) Global-mean AOD at 550 nm vs time for the fiducial case, with the other
    19 sweep members drawn faintly behind it to show the sensitivity
    envelope, plus a horizontal reference line (with a shaded +/- band) for
    the Khaykin et al. (2022) DISPERSED (not plume-tracked) 525 nm SAOD
    anomaly, 60S-60N monthly mean for March 2022.
(b) Latitude-time Hovmoller contour of zonal-mean AOD at 550 nm for the
    fiducial case, showing the aerosol confined toward the vent's own (SH)
    hemisphere (vent at 20.5 S).

Note: the model has no H2O-dependent sulfate pathway, so the retired
"water injection roughly doubles AOD" target is NOT reproducible and is not
referenced anywhere in this figure.

Reads:  base_dir/<case>/data/aod/aod_550nm_band.csv        (time series)
        base_dir/<case>/data/aod/aod_zonal_550nm_band.csv  (zonal)

Configuration (case, paths, references, output) is in
config_hunga_validation.yaml.
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
parser.add_argument('--config', default='config_hunga_validation.yaml',
                    help='YAML config file (relative to script directory)')
args = parser.parse_args()

cfg = yaml.safe_load(open(os.path.join(here, args.config)))

# The fiducial moved into the phase-2 suite, so cases now live under either
# of two roots; resolve each case against the first root that has it.
base_dirs = cfg['base_dirs']
case = cfg['case']


def resolve(case_name, subpath):
    """Absolute path to one case's CSV, searching each base_dir in order."""
    for base in base_dirs:
        p = os.path.join(base, case_name, 'data', subpath)
        if os.path.exists(p):
            return p
    return None


zonal_path = resolve(case, cfg['subpath_zonal'])
if zonal_path is None:
    raise SystemExit(f"Missing zonal CSV for {case} under {base_dirs}")

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
# t=0 corresponds to the eruption date; days_since_start map to real calendar dates.
t0 = pd.Timestamp(cfg['eruption_date'])

# --- Ensemble of global-mean time series (panel a) ---
case_dirs = sorted(
    d
    for base in base_dirs
    for d in glob.glob(os.path.join(base, cfg['case_glob']))
    if os.path.isdir(d)
)
ts_series = {}   # case_name -> pd.Series indexed by days_since_start
for d in case_dirs:
    c = os.path.basename(d)
    if c in ts_series:      # the two suites are disjoint; never double-count
        continue
    p = os.path.join(d, 'data', cfg['subpath_timeseries'])
    if not os.path.exists(p):
        print(f"  (skip, no time series) {c}")
        continue
    df = pd.read_csv(p, header=0)
    ts_series[c] = pd.Series(df.iloc[:, 1].values, index=df.iloc[:, 0].values)

if case not in ts_series:
    raise SystemExit(f"Fiducial case {case} has no time series CSV")

ts_frame = pd.DataFrame(ts_series).sort_index()
ts_days  = ts_frame.index.values.astype(float)
ts_dates = t0 + pd.to_timedelta(ts_days, unit='D')
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

# ---------------------------------------------------------------------------
# Diagnostics
# ---------------------------------------------------------------------------
fid_peak_idx = np.nanargmax(ts_fid)
fid_peak_val = ts_fid[fid_peak_idx]
fid_peak_day = ts_days[fid_peak_idx]

zon_peak_flat = np.nanargmax(zon_aod)
zon_peak_it, zon_peak_ilat = np.unravel_index(zon_peak_flat, zon_aod.shape)
zon_peak_val = zon_aod[zon_peak_it, zon_peak_ilat]
zon_peak_day = zon_days[zon_peak_it]
zon_peak_lat = lats[zon_peak_ilat]

print(f"fiducial case = {case}")
print(f"  ensemble cases (incl. fiducial) = {ts_frame.shape[1]}")
print(f"  fiducial peak global AOD = {fid_peak_val:.5f} at day {fid_peak_day:.0f}")
print(f"  zonal peak AOD           = {zon_peak_val:.5f} at day {zon_peak_day:.0f}, lat {zon_peak_lat:+.0f}")

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
fig, (ax_ts, ax_zon) = plt.subplots(
    1, 2, figsize=(7.1, 3.0),
    gridspec_kw={'width_ratios': [1, 1.15], 'wspace': 0.30},
)

# --- (a) Global-mean AOD time series, fiducial + sweep envelope ---
n_other = 0
for c in ts_frame.columns:
    if c == case:
        continue
    ax_ts.plot(ts_dates, ts_frame[c].values, lw=0.5, color='0.75', zorder=1)
    n_other += 1
# Dummy handle for the grey envelope, so it appears once in the legend.
ax_ts.plot([], [], lw=1.2, color='0.75', label=f'Sweep members (n={n_other})')

ax_ts.plot(ts_dates, ts_fid, lw=1.6, color='C3', zorder=4,
           label=cfg.get('case_label', 'Fiducial'))

# Observational reference: dispersed (not plume-tracked) SAOD anomaly, drawn
# as a solid line plus a solid pale band for +/- the quoted uncertainty.
# (EPS has no real transparency, so the band is a solid pale fill rather
# than an alpha-blended axhspan.)
for ref in cfg['references']:
    err = ref.get('error', 0.0)
    if err:
        ax_ts.axhspan(ref['value'] - err, ref['value'] + err,
                      color='#F5CBA7', zorder=2, lw=0)
    ax_ts.axhline(ref['value'], color=ref['color'], lw=1.1, ls='--', zorder=3,
                  label=f"{ref['label']}")

ax_ts.set_xlabel('Year')
ax_ts.set_ylabel('AOD at 550 nm')
ax_ts.set_xlim(*xlim)
# Derived from the ensemble rather than hardcoded: the fiducial changed in
# 2026-07 and a fixed ceiling silently clipped the new one.
ts_top = float(np.ceil(np.nanmax(ts_frame.values) * 50.0) / 50.0)
ax_ts.set_ylim(0, ts_top)
ax_ts.yaxis.set_minor_locator(ticker.AutoMinorLocator())
# Panel label in the upper-left; the eruption-driven curves rise from zero
# near the eruption date and stay well below this figure's y-max, so the
# upper-left corner is clear here (unlike the Pinatubo panel).
ax_ts.text(0.03, 0.93, '(a)', transform=ax_ts.transAxes, fontsize=9,
           fontweight='bold', va='top', ha='left')

# Legend placed below the panel (outside the axes) so it cannot obscure the
# y-axis tick labels or overrun the right-hand axis spine.
ax_ts.legend(fontsize=6, frameon=False, loc='upper center',
             bbox_to_anchor=(0.5, -0.22), ncol=1, handlelength=2.0,
             labelspacing=0.3, borderaxespad=0.0)

# --- (b) Zonal AOD contour, fiducial ---
zon_top = float(np.ceil(np.nanmax(zon_aod) * 50.0) / 50.0)
levels = np.linspace(0, zon_top, 25)

cf = ax_zon.contourf(
    zon_dates, lats, zon_aod.T,
    levels=levels,
    cmap='magma_r',
    extend='max',
    antialiased=False,
)
cf.set_edgecolor('face')
cf.set_linewidth(0)

cbar = fig.colorbar(cf, ax=ax_zon, pad=0.02, aspect=30)
cbar.set_ticks(np.linspace(0, zon_top, 6))
cbar.set_label('AOD at 550 nm')
cbar.ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())

ax_zon.set_xlabel('Year')
ax_zon.set_ylabel('Latitude (°)')
ax_zon.set_xlim(*xlim)
ax_zon.set_yticks(np.arange(-90, 91, 30))
ax_zon.yaxis.set_minor_locator(ticker.MultipleLocator(10))
ax_zon.text(0.03, 0.93, '(b)', transform=ax_zon.transAxes, fontsize=9,
           fontweight='bold', va='top', ha='left', color='white')

# Calendar ticks: major = each year (labelled), minor = each month.
for ax in (ax_ts, ax_zon):
    ax.xaxis.set_major_locator(mdates.YearLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    ax.xaxis.set_minor_locator(mdates.MonthLocator())

for ext in ('pdf', 'eps'):
    outfile = os.path.join(here, f"{cfg['outfile_stem']}.{ext}")
    plt.savefig(outfile, dpi=300, bbox_inches='tight')
    print(f"Saved: {outfile}")
plt.close()
