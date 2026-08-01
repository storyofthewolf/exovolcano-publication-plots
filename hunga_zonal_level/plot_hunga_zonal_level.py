"""
plot_hunga_zonal_level.py - Zonal-mean latitude-time fields at fixed altitude
for the Hunga Tonga case, on the basis of Schoeberl et al. (2024) Figure 3.

Two panels:
(a) Sulfate aerosol mass density at 20 km, latitude vs time, for the fiducial.
    Analogous in structure to their Figure 3a (OMPS-LP aerosol extinction
    coefficient at 20 km) but NOT on the same magnitude axis -- see the
    config header.
(b) Water vapor volume mixing ratio anomaly at 25 km, latitude vs time,
    differenced against the dry-injection control. Directly comparable in
    quantity and units to their Figure 3b (MLS water vapor at 25 km), once
    CAM's mass mixing ratio is converted to ppmv.

Reads: <base_dir>/<case>/data/zonal_level/{Q_25km.csv,VOLCHZMD_20km.csv}
produced by exovolcano-analysis' zonal_level_timeseries.py.

Configuration is in config_hunga_zonal_level.yaml.
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
import matplotlib.colors as mcolors

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
parser.add_argument('--config', default='config_hunga_zonal_level.yaml')
args = parser.parse_args()

cfg = yaml.safe_load(open(os.path.join(here, args.config)))
base_dirs = cfg['base_dirs']
t0 = pd.Timestamp(cfg['eruption_date'])
xlim = (pd.Timestamp(cfg['xlim'][0]), pd.Timestamp(cfg['xlim'][1]))


def load(entry, fname):
    """Read one zonal_level CSV. Returns (days, lats, values[time, lat], header)."""
    p = os.path.join(base_dirs[entry['phase']], entry['case'], cfg['subpath'], fname)
    if not os.path.exists(p):
        raise SystemExit(f"Missing CSV: {p}")
    meta = {}
    with open(p) as f:
        for line in f:
            if not line.startswith('#'):
                break
            k, _, v = line[1:].partition(':')
            meta[k.strip()] = v.strip()
    df = pd.read_csv(p, comment='#', index_col=0)
    lats = np.array([float(c) for c in df.columns])
    return df.index.values.astype(float), lats, df.values.astype(float), meta


fid = cfg['fiducial']
dry = cfg['dry_control']

# --- Panel (a): aerosol at 20 km ---
a_days, a_lats, a_vals, a_meta = load(fid, cfg['aerosol_file'])

# --- Panel (b): water vapor anomaly at 25 km, dry-control differenced ---
w_days, w_lats, w_vals, w_meta = load(fid, cfg['water_file'])
d_days, d_lats, d_vals, _ = load(dry, cfg['water_file'])
n = min(len(w_days), len(d_days))
# CAM Q is a mass mixing ratio; MLS reports ppmv. Convert before differencing
# so the anomaly is in the observed units.
w_anom = (w_vals[:n] - d_vals[:n]) * cfg['mmr_to_vmr'] * 1.0e6
w_days = w_days[:n]

a_dates = t0 + pd.to_timedelta(a_days, unit='D')
w_dates = t0 + pd.to_timedelta(w_days, unit='D')

print(f"aerosol level : {a_meta.get('level_altitude_km')} km, "
      f"{a_meta.get('level_pressure_Pa')} Pa")
print(f"water level   : {w_meta.get('level_altitude_km')} km, "
      f"{w_meta.get('level_pressure_Pa')} Pa")
ia = np.unravel_index(np.nanargmax(a_vals), a_vals.shape)
iw = np.unravel_index(np.nanargmax(w_anom), w_anom.shape)
print(f"aerosol peak  : {a_vals[ia]:.4g} g/cm3 at lat {a_lats[ia[1]]:+.0f}, "
      f"{a_dates[ia[0]]:%Y-%m-%d}")
print(f"water peak    : {w_anom[iw]:.3f} ppmv at lat {w_lats[iw[1]]:+.0f}, "
      f"{w_dates[iw[0]]:%Y-%m-%d}")

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
fig, (ax_a, ax_b) = plt.subplots(
    2, 1, figsize=(5.2, 5.0),
    gridspec_kw={'hspace': 0.32},
)

# --- (a) aerosol ---
amax = float(np.nanmax(a_vals))
levels_a = np.linspace(0, amax, 21)
cf_a = ax_a.contourf(a_dates, a_lats, a_vals.T, levels=levels_a,
                     cmap='magma_r', antialiased=False)
cf_a.set_edgecolor('face')
cb_a = fig.colorbar(cf_a, ax=ax_a, pad=0.02, aspect=22)
cb_a.set_label(r'Sulfate mass density (g cm$^{-3}$)')
cb_a.formatter.set_powerlimits((0, 0))
cb_a.update_ticks()
ax_a.set_ylabel('Latitude (°)')
ax_a.set_title(f"Aerosol at {float(a_meta.get('level_altitude_km', 20)):.1f} km",
               fontsize=8)
# Dark text: magma_r is pale at low values, so white is unreadable here.
ax_a.text(0.02, 0.94, '(a)', transform=ax_a.transAxes, fontsize=9,
          fontweight='bold', va='top', ha='left', color='black')

# --- (b) water vapor anomaly ---
wmax = float(np.ceil(np.nanmax(w_anom) * 2.0) / 2.0)
levels_b = np.linspace(0, wmax, 21)
# Sequential map anchored on white so the near-zero background reads as empty,
# matching panel (a)'s convention rather than a dark-background viridis.
cf_b = ax_b.contourf(w_dates, w_lats, np.clip(w_anom, 0, None).T,
                     levels=levels_b, cmap='YlGnBu', antialiased=False)
cf_b.set_edgecolor('face')
cb_b = fig.colorbar(cf_b, ax=ax_b, pad=0.02, aspect=22)
cb_b.set_label(r'$\Delta$H$_2$O (ppmv)')
ax_b.set_ylabel('Latitude (°)')
ax_b.set_xlabel('Date')
ax_b.set_title(f"Water vapor anomaly at "
               f"{float(w_meta.get('level_altitude_km', 25)):.1f} km", fontsize=8)
ax_b.text(0.02, 0.94, '(b)', transform=ax_b.transAxes, fontsize=9,
          fontweight='bold', va='top', ha='left', color='black')

for ax in (ax_a, ax_b):
    # Vent latitude, so the SH confinement can be read against it.
    ax.axhline(-20.5, color='0.35', lw=0.7, ls=':', zorder=3)
    ax.set_xlim(*xlim)
    ax.set_ylim(-90, 90)
    ax.set_yticks(np.arange(-90, 91, 30))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(10))
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=3))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b\n%Y'))
    ax.xaxis.set_minor_locator(mdates.MonthLocator())

for ext in ('pdf', 'eps'):
    outfile = os.path.join(here, f"{cfg['outfile_stem']}.{ext}")
    plt.savefig(outfile, dpi=300, bbox_inches='tight')
    print(f"Saved: {outfile}")
plt.close()
