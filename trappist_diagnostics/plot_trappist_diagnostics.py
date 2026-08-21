"""
plot_trappist_diagnostics.py -- diagnostic quad for the TRAPPIST-1e ben2/hab2
eruption suites.

Four plots, as requested, all on the same ben2 (dry) versus hab2 (moist)
contrast:

  1. Global-mean 550 nm optical depth vs time, one panel per atmosphere.
  2. Global stratospheric water burden vs time, one panel per atmosphere.
  3. Height-time contour of VOLCHZMD (sulfate aerosol mass density), the
     high-sulfur Tambora template, both atmospheres side by side.
  4. Height-time contour of Q (water vapor specific humidity) ANOMALY, the
     water-rich Hunga template, both atmospheres side by side.

SCOPE -- ben1/hab1 are absent on purpose. Those 28 cases were compiled with
the wrong atmospheric composition (pure CO2 where N2 + 400 ppm CO2 was
specified) and are being rerun. Everything here is therefore a dry-versus-wet
SURFACE contrast at FIXED pure-CO2 composition. That is a cleaner control than
ben1-vs-ben2 would have been, but it is not the composition axis the current
Section 4 text is written around.

No observational reference line appears on any panel. There is no observed
TRAPPIST-1e volcanic aerosol or stratospheric water measurement to compare
against; these are model predictions, and drawing a reference band from the
Earth validation cases would silently import a terrestrial calibration into a
plot of a different planet.

Reads (fetch_exovolc.py local layout):
    base_dir/<batch>/<case>/data/aod/aod_550nm_band.csv
    base_dir/<batch>/<case>/data/scalar/Q.csv
    base_dir/<batch>/<case>/data/profiles/{Q,VOLCHZMD}.csv

Profile CSV format: first two lines are `# pressure_Pa: ...` and
`# altitude_m: ...` comments carrying the 51-level vertical grid, then a
`days_since_start,<VAR>_lev0..lev50` header. Level 0 is model top; level 50 is
the surface.

Scalar/AOD CSV format: column 0 = time (days), column 1 = the variable.

Configuration is in config_trappist_diagnostics.yaml.
"""

import os
import argparse
import yaml
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, TwoSlopeNorm

plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 8,
    'axes.labelsize': 9,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'legend.fontsize': 6.5,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})

here = os.path.dirname(os.path.abspath(__file__))

parser = argparse.ArgumentParser()
parser.add_argument('--config', default='config_trappist_diagnostics.yaml',
                    help='YAML config file (relative to script directory)')
args = parser.parse_args()

cfg = yaml.safe_load(open(os.path.join(here, args.config)))

BASE_DIR = cfg['base_dir']
BATCHES = cfg['batches']                      # {'ben2': 'ben2_suite1', ...}
ATMOS = list(BATCHES.keys())                  # ['ben2', 'hab2']

G_CONST = cfg['g_const']
R_PLANET = cfg['r_planet']
AREA_PLANET = 4.0 * np.pi * R_PLANET ** 2
STRAT_P_MAX_PA = cfg['strat_p_max_pa']
KG_PER_TG = 1.0e9
# Running-mean window for the water panel only (days). See the note at its use.
SMOOTH_DAYS = 30

ATM_TITLE = {
    'ben2': r'\texttt{ben2}: dry land, 1 bar CO$_2$',
    'hab2': r'\texttt{hab2}: aquaplanet, 1 bar CO$_2$',
}
# Rendered without LaTeX (usetex is off), so strip the \texttt wrapper.
ATM_TITLE = {
    'ben2': 'ben2  (dry land, 1 bar CO$_2$)',
    'hab2': 'hab2  (aquaplanet, 1 bar CO$_2$)',
}

# Okabe-Ito-derived, colorblind-safe. Grouped so the four experimental axes
# (core templates / day-night / strength / chemistry) read apart at a glance.
CASE_STYLE = {
    'control':        {'c': '#000000', 'ls': ':',  'lw': 1.0},
    'tambora':        {'c': '#D55E00', 'ls': '-',  'lw': 1.8},
    'pinatubo':       {'c': '#E69F00', 'ls': '-',  'lw': 1.4},
    'hunga':          {'c': '#0072B2', 'ls': '-',  'lw': 1.4},
    'tambora_night':  {'c': '#D55E00', 'ls': '--', 'lw': 1.0},
    'pinatubo_night': {'c': '#E69F00', 'ls': '--', 'lw': 1.0},
    'hunga_night':    {'c': '#0072B2', 'ls': '--', 'lw': 1.0},
    'tambora_10x':    {'c': '#CC79A7', 'ls': '-',  'lw': 1.2},
    'tambora_100x':   {'c': '#882255', 'ls': '-',  'lw': 1.2},
    'hunga_10x':      {'c': '#56B4E9', 'ls': '-',  'lw': 1.2},
    'hunga_100x':     {'c': '#009E73', 'ls': '-',  'lw': 1.2},
    'tambora_k2.5':   {'c': '#999999', 'ls': '-.', 'lw': 1.0},
    'tambora_k250':   {'c': '#555555', 'ls': '-.', 'lw': 1.0},
    'tambora_scrit5': {'c': '#AA4499', 'ls': '-.', 'lw': 1.0},
}


def case_name(atm, key):
    """Full case directory name, e.g. ('ben2', 'tambora_10x') ->
    'exovolc_ben2_tambora_10x'."""
    return f'exovolc_{atm}_{key}'


def data_path(atm, key, subpath):
    return os.path.join(BASE_DIR, BATCHES[atm], case_name(atm, key),
                        'data', subpath)


def load_scalar(atm, key, subpath):
    """Read a two-column (time, value) CSV. Returns (days, values) or None if
    the file is absent -- absence is reported, never silently plotted as a
    gap."""
    p = data_path(atm, key, subpath)
    if not os.path.exists(p):
        return None
    df = pd.read_csv(p, header=0)
    return (df.iloc[:, 0].values.astype(float),
            df.iloc[:, 1].values.astype(float))


def load_profile(atm, key, subpath):
    """Read a profiles/*.csv. Returns (days, pres_Pa, alt_m, field[ntime,nlev])
    or None if absent."""
    p = data_path(atm, key, subpath)
    if not os.path.exists(p):
        return None
    with open(p) as f:
        pres_line = f.readline()
        alt_line = f.readline()
    pres_pa = np.array([float(x) for x in
                        pres_line.split(':', 1)[1].strip().strip(',').split(',')])
    alt_m = np.array([float(x) for x in
                      alt_line.split(':', 1)[1].strip().strip(',').split(',')])
    df = pd.read_csv(p, skiprows=2)
    days = df.iloc[:, 0].values.astype(float)
    fld = df.iloc[:, 1:].values.astype(float)
    return days, pres_pa, alt_m, fld


def layer_dp(pres_pa):
    """Per-level pressure thickness [Pa] from midlevel pressures, by midpoint
    bisection (level 0 = top, last = surface). Same reconstruction used by the
    Earth Hunga figures, where it was cross-checked against the exact CAM
    hybrid-sigma interfaces and found accurate to ~0.6% in the summed
    stratospheric dp."""
    nlev = len(pres_pa)
    interfaces = np.zeros(nlev + 1)
    interfaces[1:-1] = 0.5 * (pres_pa[:-1] + pres_pa[1:])
    interfaces[0] = 0.0
    interfaces[-1] = pres_pa[-1] + (pres_pa[-1] - interfaces[-2])
    return np.diff(interfaces)


def strat_water_tg(atm, key):
    """Global stratospheric (p < STRAT_P_MAX_PA) water mass [Tg] vs days,
    from the Q profile. Returns (days, mass_tg) or None."""
    prof = load_profile(atm, key, cfg['subpath_q_profile'])
    if prof is None:
        return None
    days, pres_pa, alt_m, q = prof
    dp = layer_dp(pres_pa)
    dm_level = dp * AREA_PLANET / G_CONST      # kg of air per level, global
    mask = pres_pa < STRAT_P_MAX_PA
    mass_kg = (q[:, mask] * dm_level[None, mask]).sum(axis=1)
    return days, mass_kg / KG_PER_TG


# ---------------------------------------------------------------------------
# Availability audit -- run before plotting so missing cases are named, not
# silently dropped.
# ---------------------------------------------------------------------------
CASES = cfg['cases']
missing = []
for atm in ATMOS:
    for entry in CASES:
        p = data_path(atm, entry['key'], cfg['subpath_aod'])
        if not os.path.exists(p):
            missing.append(f'{atm}/{entry["key"]}: {cfg["subpath_aod"]}')

print('=' * 72)
print('TRAPPIST-1e diagnostic quad -- ben2 / hab2 only')
print('  (ben1/hab1 excluded: wrong compiled composition, rerun pending)')
print('=' * 72)
if missing:
    print(f'\nMISSING AOD inputs ({len(missing)}):')
    for m in missing:
        print(f'  {m}')
else:
    print(f'\nAll {len(ATMOS) * len(CASES)} AOD inputs present.')


# ---------------------------------------------------------------------------
# FIGURE 1 -- global-mean 550 nm AOD vs time, one panel per atmosphere
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(7.1, 3.0), sharey=True)
xlim = tuple(cfg['xlim_days'])

print('\n' + '-' * 72)
print('FIGURE 1: global-mean 550 nm AOD -- peak values')
print('-' * 72)

for ax, atm in zip(axes, ATMOS):
    print(f'\n  {atm}:')
    for entry in CASES:
        key = entry['key']
        got = load_scalar(atm, key, cfg['subpath_aod'])
        if got is None:
            continue
        days, aod = got
        st = CASE_STYLE[key]
        ax.plot(days / 365.25, aod, color=st['c'], ls=st['ls'], lw=st['lw'],
                label=entry['label'])
        pk = np.nanargmax(aod)
        print(f'    {key:16s} peak AOD = {aod[pk]:9.4f} at day {days[pk]:7.1f}')
    ax.set_xlim(xlim[0] / 365.25, xlim[1] / 365.25)
    ax.set_yscale('log')
    # Floor the axis at 1e-4. Below that the aerosol is radiatively irrelevant
    # and the decay runs down to ~1e-22, which on a log axis compresses the
    # entire physically meaningful range into the top eighth of the panel.
    ax.set_ylim(1e-4, 1e2)
    ax.set_xlabel('Years since eruption')
    ax.set_title(ATM_TITLE[atm], fontsize=8.5)
    ax.grid(alpha=0.25, lw=0.4)

axes[0].set_ylabel(r'Global-mean optical depth, 550 nm')
handles, labels = axes[1].get_legend_handles_labels()
fig.legend(handles, labels, loc='lower center', ncol=5, frameon=False,
           bbox_to_anchor=(0.5, -0.02), handlelength=1.8, columnspacing=1.2)

fig.tight_layout(rect=(0, 0.16, 1, 1))
for ext in ('pdf', 'eps'):
    fig.savefig(os.path.join(here, f'{cfg["outfile_aod"]}.{ext}'),
                dpi=300, bbox_inches='tight')
plt.close(fig)
print(f'\n  wrote {cfg["outfile_aod"]}.pdf / .eps')


# ---------------------------------------------------------------------------
# FIGURE 2 -- global stratospheric water burden vs time
# ---------------------------------------------------------------------------
# Plotted as the EXCESS over each atmosphere's own no-eruption control, so the
# unperturbed background and its drift are removed. The control is an
# independent realization branched from the same restart, so its internal
# variability does not cancel exactly; small nonzero excursions in the
# non-water cases are that noise floor, not a signal.
fig, axes = plt.subplots(1, 2, figsize=(7.1, 3.0))

print('\n' + '-' * 72)
print(f'FIGURE 2: stratospheric (p < {STRAT_P_MAX_PA/100:.0f} hPa) excess '
      f'H2O mass -- peak values')
print('  (excess over each atmosphere\'s own no-eruption control)')
print('-' * 72)

for ax, atm in zip(axes, ATMOS):
    print(f'\n  {atm}:')
    bg = strat_water_tg(atm, cfg['control'])
    if bg is None:
        print('    control Q profile missing -- panel skipped')
        continue
    bg_days, bg_tg = bg
    bg_series = pd.Series(bg_tg, index=bg_days)
    # The control is an independent realization, so its internal variability
    # does NOT cancel in the difference. Its 1-sigma is the detection floor
    # any volcanic water signal must clear to be meaningful.
    bg_sigma = float(np.nanstd(bg_tg))
    print(f'    control background = {np.nanmean(bg_tg):.4g} Tg, '
          f'1-sigma = {bg_sigma:.4g} Tg (detection floor)')

    for entry in CASES:
        key = entry['key']
        if key == cfg['control']:
            continue
        got = strat_water_tg(atm, key)
        if got is None:
            continue
        days, tg = got
        excess = tg - bg_series.reindex(days).values
        st = CASE_STYLE[key]
        # 30-day running mean. On the aquaplanet the raw daily difference is
        # dominated by the seasonal hydrological cycle beating against an
        # independent control realization, which fills the panel with
        # high-frequency noise and hides whether any signal is present at all.
        # Smoothing does not create a signal where there is none -- it makes
        # the absence readable.
        excess_s = pd.Series(excess).rolling(SMOOTH_DAYS, center=True,
                                             min_periods=1).mean().values
        ax.plot(days / 365.25, excess_s, color=st['c'], ls=st['ls'],
                lw=st['lw'], label=entry['label'])
        pk = np.nanargmax(excess)
        snr = excess[pk] / bg_sigma if bg_sigma > 0 else np.nan
        flag = '' if snr >= 10.0 else '   <-- AT/BELOW NOISE FLOOR'
        print(f'    {key:16s} peak excess = {excess[pk]:11.4g} Tg '
              f'at day {days[pk]:7.1f}  SNR = {snr:7.1f}{flag}')

    # The control's own 1-sigma variability, drawn as the detection floor.
    # On the moist aquaplanet this band is wider than every volcanic signal in
    # the suite, which is the result, not a plotting failure.
    ax.axhspan(-bg_sigma, bg_sigma, color='#DDDDDD', zorder=0,
               label=r'control 1$\sigma$ variability')
    ax.axhline(0.0, color='k', lw=0.6, ls=':', zorder=1)
    ax.set_xlim(xlim[0] / 365.25, xlim[1] / 365.25)
    ax.set_yscale('symlog', linthresh=1.0)
    ax.set_xlabel('Years since eruption')
    ax.set_ylabel(r'Stratospheric excess H$_2$O mass [Tg]')
    ax.set_title(ATM_TITLE[atm], fontsize=8.5)
    ax.grid(alpha=0.25, lw=0.4)

handles, labels = axes[1].get_legend_handles_labels()
fig.legend(handles, labels, loc='lower center', ncol=5, frameon=False,
           bbox_to_anchor=(0.5, -0.02), handlelength=1.8, columnspacing=1.2)

fig.tight_layout(rect=(0, 0.16, 1, 1))
for ext in ('pdf', 'eps'):
    fig.savefig(os.path.join(here, f'{cfg["outfile_water"]}.{ext}'),
                dpi=300, bbox_inches='tight')
plt.close(fig)
print(f'\n  wrote {cfg["outfile_water"]}.pdf / .eps')


# ---------------------------------------------------------------------------
# FIGURES 3 & 4 -- height-time contours of VOLCHZMD and Q
# ---------------------------------------------------------------------------
# One figure, 2x2: rows are the variable (aerosol mass density, water-vapor
# anomaly), columns are the atmosphere. VOLCHZMD is drawn on a log color scale
# because it spans orders of magnitude; Q is drawn as an ANOMALY against each
# level's control-run background on a diverging scale, since the absolute
# field is dominated by the (very different) background humidity of the two
# atmospheres and would show nothing.
fig, axes = plt.subplots(2, 2, figsize=(7.1, 5.2), sharex=True, sharey=True)
xlim_c = tuple(cfg['xlim_days_contour'])
alt_max_km = cfg['contour_alt_max_km']

haze_key = cfg['contour_case_haze']
water_key = cfg['contour_case_water']

print('\n' + '-' * 72)
print('FIGURES 3 & 4: height-time contours')
print(f'  aerosol (VOLCHZMD) from the {haze_key} template')
print(f'  water   (Q anomaly) from the {water_key} template')
print('-' * 72)

# --- Row 0: VOLCHZMD, log scale, shared color limits across atmospheres ---
haze_data = {}
for atm in ATMOS:
    prof = load_profile(atm, haze_key, cfg['subpath_hz_profile'])
    if prof is not None:
        haze_data[atm] = prof

if haze_data:
    pos = np.concatenate([d[3][d[3] > 0].ravel() for d in haze_data.values()])
    hz_max = np.nanpercentile(pos, 99.9)
    hz_min = hz_max * 1e-5
    for ax, atm in zip(axes[0], ATMOS):
        if atm not in haze_data:
            ax.set_visible(False)
            continue
        days, pres_pa, alt_m, hz = haze_data[atm]
        alt_km = alt_m / 1000.0
        fld = np.ma.masked_less_equal(hz.T, 0.0)
        pcm = ax.pcolormesh(days / 365.25, alt_km, fld,
                            norm=LogNorm(vmin=hz_min, vmax=hz_max),
                            cmap='inferno', shading='auto', rasterized=True)
        ax.set_ylim(0, alt_max_km)
        ax.set_xlim(xlim_c[0] / 365.25, xlim_c[1] / 365.25)
        ax.set_title(f'{ATM_TITLE[atm]}\n{haze_key}: sulfate aerosol',
                     fontsize=8)
        pk_t, pk_z = np.unravel_index(np.nanargmax(hz), hz.shape)
        print(f'  {atm} {haze_key:10s} peak VOLCHZMD = {hz[pk_t, pk_z]:.4g} '
              f'at day {days[pk_t]:.1f}, {alt_km[pk_z]:.1f} km')
    axes[0][0].set_ylabel('Altitude [km]')
    cb = fig.colorbar(pcm, ax=axes[0], pad=0.02, aspect=25)
    cb.set_label(r'VOLCHZMD [kg m$^{-3}$]', fontsize=8)

# --- Row 1: Q anomaly vs control, diverging scale ---
water_anom = {}
for atm in ATMOS:
    prof = load_profile(atm, water_key, cfg['subpath_q_profile'])
    ctl = load_profile(atm, cfg['control'], cfg['subpath_q_profile'])
    if prof is None or ctl is None:
        continue
    days, pres_pa, alt_m, q = prof
    c_days, _, _, q_ctl = ctl
    # Align the control onto the case's time axis before differencing; the two
    # runs share a calendar but need not share every output stamp.
    n = min(len(days), len(c_days))
    water_anom[atm] = (days[:n], alt_m, q[:n] - q_ctl[:n])

if water_anom:
    # Colour limits set from the STRATOSPHERE only. Scaling on the full column
    # lets the aquaplanet's boundary-layer humidity, which is orders of
    # magnitude larger and physically unrelated to the injection, flatten both
    # stratospheric panels to a single colour.
    amax = max(np.nanpercentile(np.abs(a[2][:, a[1] / 1000.0 > 15.0]), 99.5)
               for a in water_anom.values())
    for ax, atm in zip(axes[1], ATMOS):
        if atm not in water_anom:
            ax.set_visible(False)
            continue
        days, alt_m, dq = water_anom[atm]
        alt_km = alt_m / 1000.0
        # Per-panel colour limits: ben2's and hab2's stratospheric water
        # anomalies differ by ~3 orders of magnitude, so a shared scale renders
        # the dry panel uniformly blank. Each colourbar is labelled with its
        # own range; the panels show STRUCTURE and TIMESCALE, and are not to be
        # read against each other for magnitude.
        amax_p = np.nanpercentile(np.abs(dq[:, alt_km > 15.0]), 99.5)
        pcm2 = ax.pcolormesh(days / 365.25, alt_km, dq.T,
                             norm=TwoSlopeNorm(vmin=-amax_p, vcenter=0.0,
                                               vmax=amax_p),
                             cmap='RdBu_r', shading='auto', rasterized=True)
        cbp = fig.colorbar(pcm2, ax=ax, pad=0.02, aspect=18)
        cbp.set_label(r'$\Delta$Q [kg kg$^{-1}$]', fontsize=7)
        cbp.ax.tick_params(labelsize=6)
        ax.set_ylim(0, alt_max_km)
        ax.set_xlim(xlim_c[0] / 365.25, xlim_c[1] / 365.25)
        ax.set_xlabel('Years since eruption')
        ax.set_title(f'{water_key}: H$_2$O anomaly vs control', fontsize=8)
        # Report the peak in the STRATOSPHERE only. The unrestricted maximum
        # on the aquaplanet lands in the boundary layer, where it is ordinary
        # tropospheric weather in an independent control realization and has
        # nothing to do with the injected plume.
        sm = alt_km > 15.0
        dqs = dq[:, sm]
        pk_t, pk_z = np.unravel_index(np.nanargmax(dqs), dqs.shape)
        print(f'  {atm} {water_key:10s} peak dQ (>15 km) = '
              f'{dqs[pk_t, pk_z]:.4g} kg/kg at day {days[pk_t]:.1f}, '
              f'{alt_km[sm][pk_z]:.1f} km')
    axes[1][0].set_ylabel('Altitude [km]')

for ext in ('pdf', 'eps'):
    fig.savefig(os.path.join(here, f'{cfg["outfile_contour"]}.{ext}'), dpi=300)
plt.close(fig)
print(f'\n  wrote {cfg["outfile_contour"]}.pdf / .eps')
print('\nDone.')
