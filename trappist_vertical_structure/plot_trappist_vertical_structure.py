"""
plot_trappist_vertical_structure.py -- height-time structure of water vapor and
sulfate aerosol for the TRAPPIST-1 e suites.

Writes ONE figure per atmosphere (ben1, ben2, hab1, hab2), each a 2 x N grid of
filled contours:

    row 1   water vapor      Q         [kg/kg]
    row 2   sulfate aerosol  VOLCHZMD  [kg/m3]

with time (years since eruption) on x and altitude (km) on y. Both rows plot the
RAW field in native units on a logarithmic scale. Nothing is differenced against
a control: the panels show what each atmosphere actually holds.

The profiles are GLOBAL MEANS -- exovolcano-analysis already took the
area-weighted horizontal average when it wrote profiles/*.csv -- so these show
the global column evolving, not a resolved plume.

Scales are shared across the panels of a row WITHIN a figure, so cases are
directly comparable there, and are NOT shared between figures, because the four
atmospheres differ in background water by many orders of magnitude and a common
scale would render the dry pair blank. The one exception is the water row on the
dry atmospheres, where the four cases themselves span ~13 orders of magnitude
and no shared scale can render them together; there each panel is scaled to its
own peak and annotated with it.

Configuration is in config_trappist_vertical_structure.yaml; see that file's
header for scope constraints, the ben/hab-vs-1/2 axis definition, and why the
dry Tambora water panels legitimately empty out after ~day 20.
"""

import os
import argparse
import yaml
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

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
parser.add_argument('--config', default='config_trappist_vertical_structure.yaml')
args = parser.parse_args()
cfg = yaml.safe_load(open(os.path.join(here, args.config)))

BASE_DIR = cfg['base_dir']
BATCHES = cfg['batches']
ATMOS = list(BATCHES.keys())
ATM_LABEL = cfg['atm_label']
CASES = cfg['cases']
ALT_MAX_CFG = cfg['contour_alt_max_km']
XLIM = cfg['xlim_years']
# Atmospheres whose water row needs a per-panel scale (see module docstring).
PER_PANEL_Q = set(cfg.get('per_panel_water_atm', []))


def dpath(atm, key, sub):
    return os.path.join(BASE_DIR, BATCHES[atm], f'exovolc_{atm}_{key}',
                        'data', sub)


def load_profile(atm, key, sub):
    """Return (days, pressure_Pa, altitude_m, field[ntime, nlev]) or None.

    Profile CSVs carry the 51-level vertical coordinate in two leading comment
    lines, then a days_since_start + one-column-per-level table.
    """
    p = dpath(atm, key, sub)
    if not os.path.exists(p):
        return None
    with open(p) as f:
        pl, al = f.readline(), f.readline()
    pres = np.array([float(x) for x in
                     pl.split(':', 1)[1].strip().strip(',').split(',')])
    alt = np.array([float(x) for x in
                    al.split(':', 1)[1].strip().strip(',').split(',')])
    d = pd.read_csv(p, skiprows=2)
    return (d.iloc[:, 0].values.astype(float), pres, alt,
            d.iloc[:, 1:].values.astype(float))


def pos_scale(fields, decades):
    """Shared (vmin, vmax) for a row: vmax is the row's own 99.9th percentile of
    positive values, vmin that many decades below. None if the row is everywhere
    non-positive (e.g. a control's aerosol, or a dry atmosphere's water)."""
    vals = [f[f > 0].ravel() for f in fields if f is not None]
    pos = np.concatenate(vals) if vals else np.array([])
    if pos.size == 0:
        return None
    vmax = float(np.nanpercentile(pos, 99.9))
    return vmax * 10.0 ** (-decades), vmax


# name, subpath, cmap, colorbar label, decades
ROWS = [
    ('q',  cfg['subpath_q_profile'],  cfg['cmap_q_pos'],
     r'Water vapor [kg kg$^{-1}$]', cfg['q_decades']),
    ('hz', cfg['subpath_hz_profile'], cfg['cmap_hz'],
     r'Sulfate aerosol [kg m$^{-3}$]', cfg['hz_decades']),
]

print('=' * 74)
print('TRAPPIST-1 e vertical structure -- water vapor and sulfate aerosol')
print('  global-mean profiles, raw fields in native units')
print('=' * 74)

for atm in ATMOS:
    print(f'\n--- {atm}: {ATM_LABEL[atm]} '.ljust(70, '-'))

    data = {}
    for rkey, sub, _, _, _ in ROWS:
        for c in CASES:
            data[(rkey, c['key'])] = load_profile(atm, c['key'], sub)

    ncase = len(CASES)
    fig, axes = plt.subplots(len(ROWS), ncase,
                             figsize=(2.0 * ncase + 1.0, 4.6),
                             sharex=True, sharey=True, squeeze=False)

    # Altitude ceiling: fit this atmosphere's own model top unless overridden.
    tops = [d[2].max() / 1000.0 for d in data.values() if d is not None]
    alt_max = (min(tops) if (ALT_MAX_CFG in (None, 'auto') and tops)
               else float(ALT_MAX_CFG))

    for irow, (rkey, sub, cmap, cblabel, decades) in enumerate(ROWS):
        fields = [data[(rkey, c['key'])][3]
                  if data[(rkey, c['key'])] is not None else None
                  for c in CASES]

        # The dry atmospheres' water spans ~13 orders of magnitude across the
        # four cases: the Hunga template injects 146 Tg H2O and holds it, while
        # Tambora injects 1 Tg against a ~17 Tg stoichiometric demand and is
        # drawn down to ~1e-21 kg/kg within weeks. No shared scale renders both,
        # so each panel is normalized to its own peak, which is printed in it.
        # Colour is then not comparable between those panels; the printed peaks
        # are. The aerosol row always keeps a shared scale, where the
        # cross-case comparison is meaningful.
        per_panel = (rkey == 'q' and atm in PER_PANEL_Q)
        dec = cfg['q_decades_dry'] if per_panel else decades
        scale = pos_scale(fields, dec)

        pcm = None
        for icol, c in enumerate(CASES):
            ax = axes[irow][icol]
            pr = data[(rkey, c['key'])]
            if pr is None:
                ax.text(0.5, 0.5, 'no data', ha='center', va='center',
                        transform=ax.transAxes, fontsize=7, color='0.4')
                continue
            days, pres, alt, f = pr
            akm = alt / 1000.0
            yrs = days / 365.25

            if scale is None or not np.any(f > 0):
                # Identically zero: a control's aerosol, or a dry atmosphere's
                # water. Say so rather than drawing an empty box.
                ax.set_facecolor('0.93')
                ax.text(0.5, 0.5, 'identically zero', ha='center', va='center',
                        transform=ax.transAxes, fontsize=6.5, color='0.35')
            else:
                vmin, vmax = scale
                if per_panel:
                    pmax = float(np.nanpercentile(f[f > 0], 99.9))
                    vmin, vmax = pmax * 10.0 ** (-dec), pmax
                    ax.text(0.97, 0.04, f'max {np.nanmax(f):.0e}',
                            transform=ax.transAxes, ha='right', va='bottom',
                            fontsize=5.6, color='0.15',
                            bbox=dict(boxstyle='round,pad=0.18', fc='white',
                                      ec='none', alpha=0.75))
                pcm = ax.pcolormesh(yrs, akm,
                                    np.ma.masked_less_equal(f.T, 0.0),
                                    norm=LogNorm(vmin=vmin, vmax=vmax),
                                    cmap=cmap, shading='auto', rasterized=True)
                it, iz = np.unravel_index(int(np.nanargmax(f)), f.shape)
                print(f'  {rkey:2s} {c["key"]:13s} peak {f[it, iz]:.3e} '
                      f'at day {days[it]:6.0f}, {akm[iz]:5.1f} km')

            ax.set_ylim(0, alt_max)
            ax.set_xlim(*XLIM)
            if irow == 0:
                ax.set_title(c['label'], fontsize=7.5)
            if irow == len(ROWS) - 1:
                ax.set_xlabel('Years since eruption')
            if icol == 0:
                ax.set_ylabel('Altitude [km]')

        if pcm is not None:
            cb = fig.colorbar(pcm, ax=axes[irow].tolist(), pad=0.015,
                              aspect=16, fraction=0.035)
            label = cblabel + (', per-panel scale' if per_panel else '')
            cb.set_label(label, fontsize=7)
            cb.ax.tick_params(labelsize=6)

    fig.suptitle(ATM_LABEL[atm], fontsize=9.5, y=0.98)

    stem = f'{cfg["outfile_stem"]}_{atm}'
    for ext in ('pdf', 'eps'):
        fig.savefig(os.path.join(here, f'{stem}.{ext}'),
                    dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'  wrote {stem}.pdf / .eps')

print('\nDone.')
