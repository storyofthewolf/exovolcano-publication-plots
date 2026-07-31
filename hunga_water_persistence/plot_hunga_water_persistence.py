"""
plot_hunga_water_persistence.py - Figure: Hunga Tonga stratospheric water
persistence validation.

This model has NO H2O-dependent sulfate/aerosol coupling (confirmed: AOD is
identical across the entire water-injection sweep). This figure therefore
validates ONLY the model's ability to reproduce the observationally defining
feature of Hunga Tonga: an anomalous, long-lived stratospheric water
enhancement that scales with injected H2O mass and persists for years
(Millan et al. 2022; Khaykin et al. 2022; Zhou et al. 2026). No aerosol
effect of the water is claimed or implied anywhere in this figure.

Two panels:
(a) Time series of the global stratospheric (p < strat_p_max_pa, default
    68 hPa) excess H2O mass [Tg], differenced against the dry-injection
    control case to remove the background seasonal cycle, one line per
    water-mass sweep case (0, 50, 100, 146 Tg H2O), colorblind-safe palette,
    fiducial (146 Tg) drawn heaviest. This is the same quantity and pressure
    definition Zhou et al. (2026, their Fig. 1/4) use for the observed HEW
    (Hunga excess water) mass, so their Fig. 4a exponential decay fit
    (MLS/TOMCAT e-folding times, stated in their text) is overlaid directly
    as a reference curve, not a digitized pixel reading. Air mass per level
    is reconstructed from the profile CSV's pressure_Pa grid by midpoint
    bisection; cross-checked against the exact CAM hybrid-sigma interfaces
    for exovolc_hunga_fid and found accurate to ~0.6% in the summed
    stratospheric dp -- adequate for this comparison.
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
# Renamed from 'fiducial' in 2026-07: this is the 146 Tg member of the water
# sweep, which is no longer the manuscript's fiducial configuration. The
# adopted fiducial is overlaid separately (fiducial_overlay) and lives in a
# different batch directory.
reference_case = cfg['reference_case']
fid_overlay = cfg.get('fiducial_overlay')

PPMV_PER_KGKG = 1.0e6   # kg/kg mixing ratio -> ppmv
KG_PER_TG = 1.0e9

G_CONST = cfg['g_const']
R_EARTH_M = cfg['r_earth_m']
AREA_EARTH = 4.0 * np.pi * R_EARTH_M**2
STRAT_P_MAX_PA = cfg['strat_p_max_pa']


def case_dir(c, batch=None):
    """Data directory for a case. `batch` overrides the default eruption
    directory, which the phase-2 overlay case needs."""
    return os.path.join(base_dir, batch or eruption, c, 'data')


def load_q_profile(case, batch=None):
    """Read a profiles/Q.csv file. Returns (days, pres_Pa, alt_m, Q[ntime, nlev])."""
    p = os.path.join(case_dir(case, batch), cfg['subpath_q_profile'])
    if not os.path.exists(p):
        raise SystemExit(f"Missing CSV: {p}")
    with open(p) as f:
        pres_line = f.readline()
        alt_line = f.readline()
    pres_pa = np.array([float(x) for x in
                         pres_line.split(':', 1)[1].strip().strip(',').split(',')])
    alt_m = np.array([float(x) for x in
                       alt_line.split(':', 1)[1].strip().strip(',').split(',')])
    df = pd.read_csv(p, skiprows=2)
    days = df.iloc[:, 0].values.astype(float)
    q = df.iloc[:, 1:].values.astype(float)   # shape (ntime, nlev)
    return days, pres_pa, alt_m, q


def layer_dp(pres_pa):
    """Interface pressures by midpoint bisection (level 0 = top, last = surface),
    then per-level thickness dp [Pa]. Cross-checked against exact CAM
    hybrid-sigma interfaces for exovolc_hunga_fid: ~0.6% accurate in the
    summed stratospheric dp (see module docstring)."""
    nlev = len(pres_pa)
    interfaces = np.zeros(nlev + 1)
    interfaces[1:-1] = 0.5 * (pres_pa[:-1] + pres_pa[1:])
    interfaces[0] = 0.0
    interfaces[-1] = pres_pa[-1] + (pres_pa[-1] - interfaces[-2])
    return np.diff(interfaces)


def strat_h2o_mass_tg_batch(case, batch=None):
    """Global stratospheric (p < STRAT_P_MAX_PA) H2O mass [Tg] vs days_since_start.
    `batch` selects a directory other than the default eruption one."""
    days, pres_pa, alt_m, q = load_q_profile(case, batch)
    dp = layer_dp(pres_pa)
    dm_level = dp * AREA_EARTH / G_CONST   # kg of air per level, global
    strat_mask = pres_pa < STRAT_P_MAX_PA
    mass_kg = (q[:, strat_mask] * dm_level[None, strat_mask]).sum(axis=1)
    return days, mass_kg / KG_PER_TG


def strat_h2o_mass_tg(case):
    return strat_h2o_mass_tg_batch(case, None)


# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
t0 = pd.Timestamp(cfg['eruption_date'])
xlim_ts  = (pd.Timestamp(cfg['xlim_ts'][0]),  pd.Timestamp(cfg['xlim_ts'][1]))
xlim_hov = (pd.Timestamp(cfg['xlim_hov'][0]), pd.Timestamp(cfg['xlim_hov'][1]))

band_lo = cfg['band_alt_min_m']
band_hi = cfg['band_alt_max_m']

# --- Panel (a): global stratospheric excess H2O mass [Tg], all sweep cases ---
# Differenced against the dry-injection control to remove the background
# seasonal cycle, matching Zhou et al. (2026)'s HEW definition.
dry_days, dry_mass_tg = strat_h2o_mass_tg(cfg['dry_control'])
dry_series = pd.Series(dry_mass_tg, index=dry_days)

mass_series = {}   # case -> pd.Series indexed by days_since_start (Tg, excess)
for entry in cfg['cases']:
    c = entry['case']
    days, mass_tg = strat_h2o_mass_tg(c)
    excess_tg = mass_tg - dry_series.reindex(days).values
    mass_series[c] = pd.Series(excess_tg, index=days)

ts_frame = pd.DataFrame(mass_series).sort_index()
ts_days  = ts_frame.index.values.astype(float)
ts_dates = t0 + pd.to_timedelta(ts_days, unit='D')

# --- Zhou et al. (2026) Fig. 4a reference curve: exponential decay from the
# stated onset date and e-folding times, peak Tg from their Fig. 1 plateau. ---
zcfg = cfg['zhou2026']
zhou_onset = pd.Timestamp(zcfg['decay_onset_date'])
zhou_dates = pd.date_range(zhou_onset, xlim_ts[1], freq='30D')
zhou_years_since_onset = (zhou_dates - zhou_onset).days / 365.25
zhou_mls_tg = zcfg['peak_tg'] * np.exp(-zhou_years_since_onset / zcfg['efold_mls_yr'])
zhou_tomcat_tg = zcfg['peak_tg'] * np.exp(-zhou_years_since_onset / zcfg['efold_tomcat_yr'])

# --- Panel (b): fiducial full-profile Q anomaly, Hovmoller ---
fid_days, fid_pres_pa, fid_alt_m, fid_q = load_q_profile(reference_case)
fid_q_ppmv = fid_q * PPMV_PER_KGKG
fid_bg_ppmv = fid_q_ppmv[0, :]                 # per-level day-0 background
fid_anom = fid_q_ppmv - fid_bg_ppmv[None, :]   # anomaly relative to each level
fid_dates = t0 + pd.to_timedelta(fid_days, unit='D')
fid_alt_km = fid_alt_m / 1000.0

# ---------------------------------------------------------------------------
# Diagnostics
# ---------------------------------------------------------------------------
print(f"Global stratospheric H2O mass: p < {STRAT_P_MAX_PA/100:.0f} hPa, "
      f"dry-control-differenced (control = {cfg['dry_control']})")
print()
print("Panel (a) sweep diagnostics (global stratospheric excess H2O mass):")
for entry in cfg['cases']:
    c = entry['case']
    s = ts_frame[c]
    peak_idx = np.nanargmax(s.values)
    peak = s.values[peak_idx]
    peak_day = s.index.values[peak_idx]
    print(f"  {c:28s} peak excess = {peak:7.1f} Tg at day {peak_day:.0f}")

# Fiducial e-folding, computed from the same decay-onset convention Zhou et
# al. (2026) use (their May-2023 onset), not from the model's own peak, so
# the two e-folding times are apples-to-apples.
fid_excess = ts_frame[reference_case].values
fid_peak_idx = np.nanargmax(fid_excess)
fid_peak_day = ts_days[fid_peak_idx]
fid_peak_val = fid_excess[fid_peak_idx]

onset_day = (zhou_onset - t0).days
onset_idx = np.searchsorted(ts_days, onset_day)
onset_val = fid_excess[onset_idx]
efold_thresh = onset_val / np.e
post_onset = np.where((ts_days >= onset_day) & (fid_excess <= efold_thresh))[0]
if len(post_onset) > 0:
    efold_day = ts_days[post_onset[0]]
    efold_yr = (efold_day - onset_day) / 365.25
    efold_str = f"{efold_yr:.2f} yr after onset (day {efold_day:.0f})"
else:
    efold_str = "not reached within run"

print()
print(f"Reference case ({reference_case}) mass diagnostics:")
print(f"  peak excess = {fid_peak_val:.1f} Tg at day {fid_peak_day:.0f} "
      f"({fid_peak_day/365.25:.2f} yr)")
print(f"  excess at Zhou et al. decay-onset date ({zhou_onset.date()}) = "
      f"{onset_val:.1f} Tg")
print(f"  e-folds (from onset) to {efold_thresh:.1f} Tg at {efold_str}")
print(f"  cf. Zhou et al. (2026): e-folding {zcfg['efold_mls_yr']} yr (MLS), "
      f"{zcfg['efold_tomcat_yr']} yr (TOMCAT), from the same onset convention")

# Vertical-structure diagnostics for panel (b), unchanged in method.
strat_mask_b = fid_alt_m > 15000.0
flat_peak = np.nanargmax(fid_anom[:, strat_mask_b])
strat_alt_km = fid_alt_km[strat_mask_b]
pi, pj = np.unravel_index(flat_peak, fid_anom[:, strat_mask_b].shape)
peak_anom_val = fid_anom[:, strat_mask_b][pi, pj]
peak_anom_day = fid_days[pi]
peak_anom_alt_km = strat_alt_km[pj]
print()
print(f"Reference case ({reference_case}) vertical-structure diagnostics (panel b):")
print(f"  peak stratospheric (>15 km) grid-level anomaly = {peak_anom_val:.4f} "
      f"ppmv at day {peak_anom_day:.0f}, altitude {peak_anom_alt_km:.1f} km")

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
fig, (ax_ts, ax_hov) = plt.subplots(
    1, 2, figsize=(7.1, 3.0),
    gridspec_kw={'width_ratios': [1, 1.15], 'wspace': 0.32},
)

# --- Panel (a): global stratospheric excess H2O mass, sweep cases ---
# Colorblind-safe (Okabe-Ito-derived) palette, light -> heavy H2O mass.
sweep_colors = ['#999999', '#0072B2', '#E69F00', '#D55E00']

for entry, color in zip(cfg['cases'], sweep_colors):
    c = entry['case']
    is_fid = (c == reference_case)
    ax_ts.plot(ts_dates, ts_frame[c].values,
               lw=2.0 if is_fid else 1.2,
               color=color,
               zorder=5 if is_fid else 3,
               label=entry['label'])

# The manuscript fiducial, overlaid rather than substituted into the sweep:
# it differs from the 146 Tg member in SO2 mass and Reff as well as being the
# adopted configuration, so putting it IN the sweep would break the
# one-parameter control that makes the sweep's linearity meaningful. Drawn
# thin over the 146 Tg line it nearly coincides with.
if fid_overlay:
    ov_days, ov_mass = strat_h2o_mass_tg_batch(
        fid_overlay['case'], fid_overlay.get('batch'))
    ov_excess = ov_mass - dry_series.reindex(ov_days).values
    ov_dates = t0 + pd.to_timedelta(ov_days, unit='D')
    ax_ts.plot(ov_dates, ov_excess, lw=1.0, color='#009E73', ls='-',
               zorder=6, label=fid_overlay['label'])
    ref_peak = np.nanmax(ts_frame[reference_case].values)
    print(f"\nFiducial overlay ({fid_overlay['case']}):")
    print(f"  peak excess          {np.nanmax(ov_excess):.1f} Tg "
          f"(146 Tg sweep member: {ref_peak:.1f} Tg)")

# Zhou et al. (2026) Fig. 4a reference curves (their stated exponential fit,
# not digitized pixel data) -- MLS observations and TOMCAT model, both scaled
# to the same peak Tg as a shape/timescale comparison, not an amplitude claim.
ax_ts.plot(zhou_dates, zhou_mls_tg, color='black', lw=1.1, ls='--',
           zorder=4, label='Zhou et al. (2026), MLS')
ax_ts.plot(zhou_dates, zhou_tomcat_tg, color='black', lw=1.1, ls=':',
           zorder=4, label='Zhou et al. (2026), TOMCAT')

ax_ts.set_xlabel('Year')
ax_ts.set_ylabel(r'Stratospheric excess H$_2$O mass (Tg)')
ax_ts.set_xlim(*xlim_ts)
ax_ts.set_ylim(bottom=0)
ax_ts.yaxis.set_minor_locator(ticker.AutoMinorLocator())
ax_ts.legend(fontsize=6, frameon=False, loc='upper right')
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
# Rounded up to a clean value so the colorbar carries readable ticks rather
# than the raw data maximum's decimals.
vmax = float(np.ceil(np.nanmax(fid_anom[np.ix_(hov_time_mask, hov_alt_mask)]) * 2.0) / 2.0)
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
