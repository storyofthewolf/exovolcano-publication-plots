"""
plot_trappist_vertical_structure.py -- height-time structure of water vapor and
sulfate aerosol for the TRAPPIST-1 e suites.

Writes ONE figure per atmosphere (ben1, ben2, hab1, hab2), each a 2 x N grid of
filled contours:

    row 1   water vapor ANOMALY   Q - Q_control   [kg/kg]
    row 2   sulfate aerosol       VOLCHZMD        [kg/m3]

with time (years since eruption) on x and altitude (km) on y. The profiles are
GLOBAL MEANS -- exovolcano-analysis already took the area-weighted horizontal
average when it wrote profiles/*.csv -- so these show the global column
evolving, not a resolved plume.

The water row is differenced against each atmosphere's own no-eruption control
because on the moist atmospheres the ambient column exceeds the injected plume
by five or more orders of magnitude, and a raw-Q panel of an eruption is
visually identical to its control. It is drawn on a symmetric-log diverging
scale, the anomaly being signed. On the dry atmospheres the control holds
Q == 0 identically, so there the anomaly and the raw field coincide exactly.

Colour scales are shared across the panels of a row WITHIN a figure, so cases
are directly comparable there, and are NOT shared between figures, because the
four atmospheres differ in background water by many orders of magnitude and a
common scale would render the dry pair blank.

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
from matplotlib.colors import LogNorm, SymLogNorm

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
CONTROL_KEY = 'control'


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
    """Shared (vmin, vmax) for a strictly positive row: vmax is the row's own
    99.9th percentile of positive values, vmin that many decades below.
    None if the row is everywhere non-positive (e.g. a control's aerosol)."""
    vals = [f[f > 0].ravel() for f in fields if f is not None]
    pos = np.concatenate(vals) if vals else np.array([])
    if pos.size == 0:
        return None
    vmax = float(np.nanpercentile(pos, 99.9))
    return vmax * 10.0 ** (-decades), vmax


def sym_scale(fields, decades, noise=None):
    """Shared (linthresh, vmax) for a signed row.

    vmax is the 99.9th percentile of |anomaly|. linthresh -- where the
    symmetric-log scale goes linear into zero -- is set to the control run own
    internal variability when that is known, so that anomalies indistinguishable
    from noise render as neutral colour instead of as vivid red/blue mottling.
    Without that floor the moist atmospheres show nothing but their own
    variability, which is larger than the injected plume. Falls back to
    `decades` below vmax when no noise estimate is available (the dry
    atmospheres, whose controls are identically zero and carry no noise).
    None if every field is identically zero."""
    vals = [np.abs(f[f != 0]).ravel() for f in fields if f is not None]
    mag = np.concatenate(vals) if vals else np.array([])
    if mag.size == 0:
        return None
    vmax = float(np.nanpercentile(mag, 99.9))
    lin = noise if (noise is not None and noise > 0) else vmax * 10.0 ** (-decades)
    return min(lin, vmax * 0.5), vmax


# name, subpath, cmap, colorbar label, decades, signed?
ROWS = [
    ('q',  cfg['subpath_q_profile'],  cfg['cmap_q'],
     r'$\Delta$ water vapor [kg kg$^{-1}$]', cfg['q_decades'], True),
    ('hz', cfg['subpath_hz_profile'], cfg['cmap_hz'],
     r'Sulfate aerosol [kg m$^{-3}$]', cfg['hz_decades'], False),
]

print('=' * 74)
print('TRAPPIST-1 e vertical structure -- water vapor and sulfate aerosol')
print('  global-mean profiles; one figure per atmosphere')
print('  water row is an anomaly against each atmosphere own control')
print('=' * 74)

for atm in ATMOS:
    print(f'\n--- {atm}: {ATM_LABEL[atm]} '.ljust(70, '-'))

    # Load every panel up front so each row can be scaled against itself.
    data = {}
    for rkey, sub, _, _, _, signed in ROWS:
        ctrl = load_profile(atm, CONTROL_KEY, sub)
        for c in CASES:
            pr = load_profile(atm, c['key'], sub)
            if pr is not None and signed:
                if c['key'] == CONTROL_KEY:
                    # Zero by construction; kept so the column stays aligned.
                    pr = (pr[0], pr[1], pr[2], np.zeros_like(pr[3]))
                elif ctrl is not None:
                    n = min(len(pr[0]), len(ctrl[0]))
                    pr = (pr[0][:n], pr[1], pr[2],
                          pr[3][:n] - ctrl[3][:n])
            data[(rkey, c['key'])] = pr

    ncase = len(CASES)
    fig, axes = plt.subplots(len(ROWS), ncase,
                             figsize=(2.0 * ncase + 1.0, 4.6),
                             sharex=True, sharey=True, squeeze=False)

    # Altitude ceiling: fit this atmosphere's own model top unless overridden.
    tops = [d[2].max() / 1000.0 for d in data.values() if d is not None]
    alt_max = (min(tops) if (ALT_MAX_CFG in (None, 'auto') and tops)
               else float(ALT_MAX_CFG))

    for irow, (rkey, sub, cmap, cblabel, decades, signed) in enumerate(ROWS):
        fields = [data[(rkey, c['key'])][3]
                  if data[(rkey, c['key'])] is not None else None
                  for c in CASES]

        # A signed row only needs the diverging symmetric-log treatment where
        # the anomaly actually goes negative -- i.e. on the moist atmospheres,
        # where it is a perturbation on a large background. On the dry
        # atmospheres the control is Q == 0, so the anomaly IS the raw field:
        # strictly positive, spanning ~13 orders of magnitude across cases, and
        # far better served by a plain log scale. Forcing symlog there wastes
        # half the colourmap on negatives that do not exist and flattens the
        # Hunga plume against the vanishing Tambora one.
        has_neg = any(f is not None and np.any(f < 0) for f in fields)
        use_sym = signed and has_neg

        if use_sym:
            # Native kg/kg on a symmetric-log scale. The anomaly is signed and
            # spans a wide range with height -- on hab1 the stratospheric plume
            # is ~1.6e-5 kg/kg while boundary-layer variability reaches ~1e-3,
            # some 60x larger -- so a linear scale would let the wet
            # troposphere swamp the thin stratosphere. `q_linthresh` sets where
            # the scale turns linear into zero; below it the anomaly is small
            # enough to be indistinguishable from the control's own
            # variability at most levels, and renders near-white.
            scale = (float(cfg['q_linthresh_by_atm'][atm]),
                     float(cfg['q_vmax_by_atm'][atm]))
        else:
            # The dry water row spans far more decades than the moist one.
            dec = (cfg['q_decades_dry'] if (signed and not has_neg)
                   else decades)
            scale = pos_scale(fields, dec)
            # On the dry atmospheres the four water panels differ by up to
            # ~13 orders of magnitude (Hunga injects 146 Tg H2O and keeps it;
            # Tambora injects 1 Tg against a ~17 Tg demand and is stripped to
            # ~1e-21 within weeks). No shared scale can render both, so each
            # panel gets its own and is annotated with its peak. Cross-panel
            # comparison of colour is meaningless here and is not intended;
            # compare the printed peaks instead. The aerosol row keeps a
            # shared scale, where the comparison IS meaningful.
            per_panel = signed and not has_neg
        if use_sym:
            print(f'  {rkey:2s} anomaly in native kg/kg, symlog: linear within '
                  f'+/-{cfg["q_linthresh_by_atm"][atm]:.0e}, log to '
                  f'+/-{cfg["q_vmax_by_atm"][atm]:.0e}')

        if not (signed and not has_neg):
            per_panel = False

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

            blank = signed and c['key'] == CONTROL_KEY
            if scale is None or blank or not np.any(f != 0):
                # A dry control's water, any control's aerosol, or the water
                # row's control panel, which is zero by construction. Say so
                # rather than drawing an empty box.
                ax.set_facecolor('0.93')
                msg = ('zero by construction' if blank else 'identically zero')
                ax.text(0.5, 0.5, msg, ha='center', va='center',
                        transform=ax.transAxes, fontsize=6.5, color='0.35')
            else:
                if use_sym:
                    lin, vmax = scale
                    norm = SymLogNorm(linthresh=lin, vmin=-vmax, vmax=vmax,
                                      base=10)
                    plot_f = f.T
                else:
                    vmin, vmax = scale
                    if per_panel:
                        pmax = float(np.nanpercentile(f[f > 0], 99.9))
                        vmin, vmax = pmax * 10.0 ** (-decades), pmax
                    norm = LogNorm(vmin=vmin, vmax=vmax)
                    plot_f = np.ma.masked_less_equal(f.T, 0.0)
                row_cmap = cmap if use_sym or not signed else cfg['cmap_q_pos']
                pcm = ax.pcolormesh(yrs, akm, plot_f, norm=norm, cmap=row_cmap,
                                    shading='auto', rasterized=True)
                if per_panel:
                    ax.text(0.97, 0.04, f'max {np.nanmax(f):.0e}',
                            transform=ax.transAxes, ha='right', va='bottom',
                            fontsize=5.6, color='0.15',
                            bbox=dict(boxstyle='round,pad=0.18', fc='white',
                                      ec='none', alpha=0.75))
                it, iz = np.unravel_index(int(np.nanargmax(np.abs(f))), f.shape)
                print(f'  {rkey:2s} {c["key"]:13s} peak {f[it, iz]:+.3e} '
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
            # On the dry atmospheres the "anomaly" is the raw field, so label
            # it as such rather than implying a difference was taken.
            if use_sym or not signed:
                label = cblabel
            elif per_panel:
                label = r'Water vapor [kg kg$^{-1}$], per-panel scale'
            else:
                label = r'Water vapor [kg kg$^{-1}$]'
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
