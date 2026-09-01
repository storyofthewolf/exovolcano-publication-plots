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

Colour scales are FIXED PER SUITE, set in the config: within one figure every
panel shares one scale per variable, so the five cases are directly comparable
and colour reads as value; between figures the scales differ, the four
atmospheres spanning far too much dynamic range to share one legibly. Values
below the floor are clamped to the lowest colour rather than masked, so no panel
shows white holes where the field is merely small; a field that is identically
zero is drawn as a labelled grey panel instead, which is a different statement
from "small".

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


# name, subpath, cmap, colorbar label; limits come from cfg['scales'][atm]
ROWS = [
    ('q',  cfg['subpath_q_profile'],  cfg['cmap_q'],
     r'Water vapor [kg kg$^{-1}$]'),
    ('hz', cfg['subpath_hz_profile'], cfg['cmap_hz'],
     r'Sulfate aerosol [kg m$^{-3}$]'),
]

print('=' * 74)
print('TRAPPIST-1 e vertical structure -- water vapor and sulfate aerosol')
print('  global-mean profiles, raw fields in native units')
print('=' * 74)

for atm in ATMOS:
    print(f'\n--- {atm}: {ATM_LABEL[atm]} '.ljust(70, '-'))

    data = {}
    for rkey, sub, _, _ in ROWS:
        for c in CASES:
            data[(rkey, c['key'])] = load_profile(atm, c['key'], sub)

    ncase = len(CASES)
    fig, axes = plt.subplots(len(ROWS), ncase,
                             figsize=(1.85 * ncase + 1.2, 4.6),
                             sharex=True, sharey=True, squeeze=False)

    # Altitude ceiling: fit this atmosphere's own model top unless overridden.
    tops = [d[2].max() / 1000.0 for d in data.values() if d is not None]
    alt_max = (min(tops) if (ALT_MAX_CFG in (None, 'auto') and tops)
               else float(ALT_MAX_CFG))

    for irow, (rkey, sub, cmap, cblabel) in enumerate(ROWS):
        vmin, vmax = (float(v) for v in cfg['scales'][atm][rkey])

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

            if not np.any(f > 0):
                # Identically zero everywhere -- a control's aerosol, or a dry
                # control's water. That is a different statement from "small",
                # so it gets a labelled panel rather than a floor-coloured one.
                ax.set_facecolor('0.93')
                ax.text(0.5, 0.5, 'identically zero', ha='center', va='center',
                        transform=ax.transAxes, fontsize=6.5, color='0.35')
            else:
                # Clamp into range rather than masking: a value below the floor
                # is small, not missing, and masking it punched white holes
                # through the aerosol panels.
                plot_f = np.clip(f.T, vmin, vmax)
                pcm = ax.pcolormesh(yrs, akm, plot_f,
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
            cb.set_label(cblabel, fontsize=7)
            cb.ax.tick_params(labelsize=6)

    fig.suptitle(ATM_LABEL[atm], fontsize=9.5, y=0.98)

    stem = f'{cfg["outfile_stem"]}_{atm}'
    for ext in ('pdf', 'eps'):
        fig.savefig(os.path.join(here, f'{stem}.{ext}'),
                    dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'  wrote {stem}.pdf / .eps')

print('\nDone.')
