"""
plot_hunga_water_persistence.py - Figure: Hunga Tonga stratospheric water
persistence validation.

This model has NO H2O-dependent sulfate/aerosol coupling (confirmed: AOD is
identical across the entire water-injection sweep). This figure therefore
validates ONLY the model's ability to reproduce the observationally defining
feature of Hunga Tonga: an anomalous, long-lived stratospheric water
enhancement that scales with injected H2O mass and persists for years
(Millan et al. 2022; Khaykin et al. 2022). No aerosol effect of the water is
claimed or implied anywhere in this figure.

Two panels:
(a) Time series of the injection-layer (~25-40 km) area-mean water vapor
    mixing ratio (ppmv) vs calendar time, one line per water-mass sweep case
    (0, 50, 100, 146 Tg H2O), colorblind-safe palette, fiducial (146 Tg) drawn
    heaviest. Spans the full 6-year run to show both the scaling with
    injected mass and the multi-year decay/persistence.
(b) Height-time (Hovmoller) contour of the water vapor mixing ratio ANOMALY
    (relative to each level's day-0 background) for the fiducial case only,
    over the first ~2 years, showing the vertical structure and slow
    decay/descent of the water-anomaly layer. This is the closest reachable
    analog to the literature's rising-water discussion (Schoeberl et al.
    2022); we do not claim to reproduce their specific rising-water /
    descending-aerosol vertical separation, only our own model's water
    anomaly structure, shown honestly.

Reads (current on-disk layout):
    base_dir/<eruption>/<case_name>/data/profiles/Q.csv

Profile CSV format: first two lines are `# pressure_Pa: ...` and
`# altitude_m: ...` comment lines carrying the 51-level vertical grid,
followed by a `days_since_start,Q_lev0..Q_lev50` header. Level 0 is the model
top (~1.08 Pa, ~75 km); level 50 is the surface.

Configuration (cases, paths, altitude band, output) is in
config_hunga_water_persistence.yaml.
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
parser.add_argument('--config', default='config_hunga_water_persistence.yaml',
                    help='YAML config file (relative to script directory)')
args = parser.parse_args()

cfg = yaml.safe_load(open(os.path.join(here, args.config)))

base_dir = cfg['base_dir']
eruption = cfg['eruption']
fiducial = cfg['fiducial']

PPMV_PER_KGKG = 1.0e6   # kg/kg mixing ratio -> ppmv


def case_dir(c):
    return os.path.join(base_dir, eruption, c, 'data')


def load_q_profile(case):
    """Read a profiles/Q.csv file. Returns (days, alt_m, Q[ntime, nlev])."""
    p = os.path.join(case_dir(case), cfg['subpath_q_profile'])
    if not os.path.exists(p):
        raise SystemExit(f"Missing CSV: {p}")
    with open(p) as f:
        f.readline()                 # pressure_Pa comment, unused here
        alt_line = f.readline()
    alt_m = np.array([float(x) for x in
                       alt_line.split(':', 1)[1].strip().strip(',').split(',')])
    df = pd.read_csv(p, skiprows=2)
    days = df.iloc[:, 0].values.astype(float)
    q = df.iloc[:, 1:].values.astype(float)   # shape (ntime, nlev)
    return days, alt_m, q


# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
t0 = pd.Timestamp(cfg['eruption_date'])
xlim_ts  = (pd.Timestamp(cfg['xlim_ts'][0]),  pd.Timestamp(cfg['xlim_ts'][1]))
xlim_hov = (pd.Timestamp(cfg['xlim_hov'][0]), pd.Timestamp(cfg['xlim_hov'][1]))

band_lo = cfg['band_alt_min_m']
band_hi = cfg['band_alt_max_m']

# --- Panel (a): injection-layer band-mean Q time series, all sweep cases ---
band_series = {}   # case -> pd.Series indexed by days_since_start (ppmv)
alt_ref = None
band_idx_ref = None
for entry in cfg['cases']:
    c = entry['case']
    days, alt_m, q = load_q_profile(c)
    band_mask = (alt_m >= band_lo) & (alt_m <= band_hi)
    band_mean_ppmv = q[:, band_mask].mean(axis=1) * PPMV_PER_KGKG
    band_series[c] = pd.Series(band_mean_ppmv, index=days)
    if c == fiducial:
        alt_ref = alt_m
        band_idx_ref = np.where(band_mask)[0]

ts_frame = pd.DataFrame(band_series).sort_index()
ts_days  = ts_frame.index.values.astype(float)
ts_dates = t0 + pd.to_timedelta(ts_days, unit='D')

# --- Panel (b): fiducial full-profile Q anomaly, Hovmoller ---
fid_days, fid_alt_m, fid_q = load_q_profile(fiducial)
fid_q_ppmv = fid_q * PPMV_PER_KGKG
fid_bg_ppmv = fid_q_ppmv[0, :]                 # per-level day-0 background
fid_anom = fid_q_ppmv - fid_bg_ppmv[None, :]   # anomaly relative to each level
fid_dates = t0 + pd.to_timedelta(fid_days, unit='D')
fid_alt_km = fid_alt_m / 1000.0

# ---------------------------------------------------------------------------
# Diagnostics
# ---------------------------------------------------------------------------
print(f"Injection-layer band: {band_lo/1000:.0f}-{band_hi/1000:.0f} km "
      f"(levels {band_idx_ref.min()}-{band_idx_ref.max()} of the fiducial grid)")
print()
print("Panel (a) sweep diagnostics (band-mean water vapor mixing ratio):")
for entry in cfg['cases']:
    c = entry['case']
    s = ts_frame[c]
    bg = s.iloc[0]
    peak_idx = np.nanargmax(s.values)
    peak = s.values[peak_idx]
    peak_day = s.index.values[peak_idx]
    print(f"  {c:28s} day-0 = {bg:.4f} ppmv, peak = {peak:.4f} ppmv "
          f"at day {peak_day:.0f}, ratio to background = {peak/bg:.2f}")

# Fiducial vertical-structure diagnostics: peak anomaly level/altitude, and
# e-folding decay time of the injection-layer band-mean anomaly. Restrict the
# level search to the stratosphere (> 15 km) -- the raw full-column max is
# dominated by tropospheric water-cycle noise near the surface (natural
# mixing ratios there are orders of magnitude larger than the injection
# signal), not the volcanic anomaly.
strat_mask = fid_alt_m > 15000.0
flat_peak = np.nanargmax(fid_anom[:, strat_mask])
strat_alt_km = fid_alt_km[strat_mask]
pi, pj = np.unravel_index(flat_peak, fid_anom[:, strat_mask].shape)
peak_anom_val = fid_anom[:, strat_mask][pi, pj]
peak_anom_day = fid_days[pi]
peak_anom_alt_km = strat_alt_km[pj]

band_anom_fid = (fid_q_ppmv[:, band_idx_ref].mean(axis=1)
                  - fid_bg_ppmv[band_idx_ref].mean())
band_peak_idx = np.nanargmax(band_anom_fid)
band_peak_val = band_anom_fid[band_peak_idx]
band_peak_day = fid_days[band_peak_idx]
efold_thresh = band_peak_val / np.e
post_peak = np.where(
    (fid_days >= band_peak_day) & (band_anom_fid <= efold_thresh)
)[0]
if len(post_peak) > 0:
    efold_day = fid_days[post_peak[0]]
    efold_elapsed = efold_day - band_peak_day
    efold_str = f"day {efold_day:.0f} ({efold_elapsed:.0f} d after peak)"
else:
    efold_str = "not reached within run (>2190 d after peak)"

print()
print(f"Fiducial ({fiducial}) vertical-structure diagnostics:")
print(f"  peak stratospheric (>15 km) grid-level anomaly = {peak_anom_val:.4f} "
      f"ppmv at day {peak_anom_day:.0f}, altitude {peak_anom_alt_km:.1f} km")
print(f"  injection-layer band-mean peak anomaly = {band_peak_val:.4f} ppmv "
      f"at day {band_peak_day:.0f}")
print(f"  decays to 1/e of peak ({efold_thresh:.4f} ppmv) at {efold_str}")

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
fig, (ax_ts, ax_hov) = plt.subplots(
    1, 2, figsize=(7.1, 3.0),
    gridspec_kw={'width_ratios': [1, 1.15], 'wspace': 0.32},
)

# --- Panel (a): band-mean water vapor mixing ratio, sweep cases ---
# Colorblind-safe (Okabe-Ito-derived) palette, light -> heavy H2O mass.
sweep_colors = ['#999999', '#0072B2', '#E69F00', '#D55E00']

for entry, color in zip(cfg['cases'], sweep_colors):
    c = entry['case']
    is_fid = (c == fiducial)
    ax_ts.plot(ts_dates, ts_frame[c].values,
               lw=2.0 if is_fid else 1.2,
               color=color,
               zorder=5 if is_fid else 3,
               label=entry['label'])

ax_ts.set_xlabel('Year')
ax_ts.set_ylabel(r'H$_2$O mixing ratio, 25$-$40 km (ppmv)')
ax_ts.set_xlim(*xlim_ts)
ax_ts.set_ylim(bottom=0)
ax_ts.yaxis.set_minor_locator(ticker.AutoMinorLocator())
ax_ts.legend(fontsize=6.5, frameon=False, loc='upper right')
ax_ts.text(0.02, 0.96, '(a)', transform=ax_ts.transAxes,
           fontsize=9, fontweight='bold', va='top', ha='left')

# --- Panel (b): fiducial Q anomaly Hovmoller ---
# Set the color scale from the plotted stratospheric altitude/time window
# only. The raw full-column anomaly is dominated by tropospheric water-cycle
# noise near the surface (orders of magnitude larger than the volcanic
# signal), which would otherwise wash out the stratospheric structure this
# panel is meant to show.
hov_alt_mask = (fid_alt_km >= 15) & (fid_alt_km <= 50)
hov_time_mask = (fid_dates >= xlim_hov[0]) & (fid_dates <= xlim_hov[1])
vmax = np.nanmax(fid_anom[np.ix_(hov_time_mask, hov_alt_mask)])
levels = np.linspace(0, vmax, 25)

cf = ax_hov.contourf(
    fid_dates, fid_alt_km, fid_anom.T,
    levels=levels,
    cmap='magma_r',
    extend='max',
)
cf.set_edgecolor('face')
cf.set_linewidth(0)

cbar = fig.colorbar(cf, ax=ax_hov, pad=0.02, aspect=30)
cbar.set_label(r'H$_2$O anomaly (ppmv)')
cbar.ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())

# Mark the injection-layer band for visual cross-reference with panel (a).
ax_hov.axhline(band_lo / 1000.0, color='0.9', lw=0.7, ls=':', zorder=4)
ax_hov.axhline(band_hi / 1000.0, color='0.9', lw=0.7, ls=':', zorder=4)

ax_hov.set_xlabel('Year')
ax_hov.set_ylabel('Altitude (km)')
ax_hov.set_xlim(*xlim_hov)
ax_hov.set_ylim(15, 50)
ax_hov.yaxis.set_minor_locator(ticker.AutoMinorLocator())
ax_hov.text(0.02, 0.96, '(b)', transform=ax_hov.transAxes,
            fontsize=9, fontweight='bold', va='top', ha='left', color='black')

# Calendar ticks: major = each year (labelled), minor = each month.
ax_ts.xaxis.set_major_locator(mdates.YearLocator())
ax_ts.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
ax_ts.xaxis.set_minor_locator(mdates.MonthLocator())

ax_hov.xaxis.set_major_locator(mdates.MonthLocator(bymonth=[1, 7]))
ax_hov.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
ax_hov.xaxis.set_minor_locator(mdates.MonthLocator())
for label in ax_hov.get_xticklabels():
    label.set_rotation(30)
    label.set_ha('right')

outfile_pdf = os.path.join(here, cfg['outfile_base'] + '.pdf')
outfile_eps = os.path.join(here, cfg['outfile_base'] + '.eps')
plt.savefig(outfile_pdf, dpi=300, bbox_inches='tight')
plt.savefig(outfile_eps, dpi=300, bbox_inches='tight')
plt.close()
print()
print(f"Saved: {outfile_pdf}")
print(f"Saved: {outfile_eps}")
