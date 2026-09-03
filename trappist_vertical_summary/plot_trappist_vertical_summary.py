"""
plot_trappist_vertical_summary.py -- one 2 x 2 summary of where volcanic
aerosol sits and how long it stays, across all four TRAPPIST-1 e atmospheres.

This is the main-text replacement for the four per-atmosphere height-time
figures, which move to an appendix. It answers the question those four are
read for -- how does the removal regime differ between atmospheres? -- on a
single shared colour scale, which the per-suite scaling of the appendix
figures cannot do.

Each panel is the sulfate aerosol mass density of one atmosphere's
tambora_100x case, annotated with the e-folding time of its own
column-integrated burden. Those four numbers are the removal-regime result in
its most compact form: ~776 d on ben1 down to ~38 d on hab2, a factor of 20
set by surface water alone.

Configuration is in config_trappist_vertical_summary.yaml; see that file's
header for why the scale is shared here but per-suite in the appendix figures,
and why aerosol rather than water.
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
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})

here = os.path.dirname(os.path.abspath(__file__))
ap = argparse.ArgumentParser()
ap.add_argument('--config', default='config_trappist_vertical_summary.yaml')
args = ap.parse_args()
with open(os.path.join(here, args.config)) as f:
    cfg = yaml.safe_load(f)


FIELD = cfg.get('field', 'aerosol')
SUBPATH = (cfg['subpath_q_profile'] if FIELD == 'water'
           else cfg['subpath_hz_profile'])
LIMITS = (cfg['q_limits'] if FIELD == 'water' else cfg['hz_limits'])
CBLABEL = (r'water vapour [kg kg$^{-1}$]' if FIELD == 'water'
           else r'sulfate aerosol [kg m$^{-3}$]')
CMAP = (cfg.get('cmap_q', 'viridis') if FIELD == 'water' else cfg['cmap_hz'])

#: Printed on each panel, so the figure names its own eruption rather than
#: relying on the caption.
CASE_LABEL = {
    'tambora': 'Tambora 1$\\times$', 'tambora_10x': 'Tambora 10$\\times$',
    'tambora_100x': 'Tambora 100$\\times$',
    'hunga': 'Hunga 1$\\times$', 'hunga_10x': 'Hunga 10$\\times$',
    'hunga_100x': 'Hunga 100$\\times$',
    'pinatubo': 'Pinatubo', 'control': 'control (no eruption)',
}


def load(atm, case):
    """Global-mean height-time field: (days, altitude_km, field)."""
    path = os.path.join(cfg['base_dir'], f"{atm}_suite1",
                        f"exovolc_{atm}_{case}", SUBPATH)
    alt = None
    with open(path) as fh:
        for line in fh:
            if line.startswith('# altitude_m:'):
                alt = np.array([float(x) for x in
                                line.split(':', 1)[1].split(',') if x.strip()])
            if not line.startswith('#'):
                break
    df = pd.read_csv(path, comment='#')
    return df.iloc[:, 0].values, alt / 1000.0, df.iloc[:, 1:].values


def efold_days(days, field):
    """e-folding time of the column-integrated burden, from its own peak.

    Returns NaN if the burden never falls to 1/e within the record, which is
    the honest answer for a case still near peak at the end of the run.
    """
    col = field.sum(axis=1)
    ipk = int(np.argmax(col))
    post, dpost = col[ipk:], days[ipk:]
    below = np.where(post < col[ipk] / np.e)[0]
    return (dpost[below[0]] - days[ipk]) if below.size else np.nan


vmin, vmax = (float(v) for v in LIMITS)
panels = cfg['panels']
fig, axes = plt.subplots(2, 2, figsize=cfg['figsize'], sharex=True, sharey=True)

tops = []
for spec in panels:
    _, alt, _ = load(spec['atm'], cfg['case'])
    tops.append(alt.max())
alt_max = min(tops) if cfg['contour_alt_max_km'] in (None, 'auto') else \
    float(cfg['contour_alt_max_km'])

print(f"{FIELD} summary, case {cfg['case']}")
mesh = None
for k, spec in enumerate(panels):
    ax = axes[k // 2, k % 2]
    days, alt, f = load(spec['atm'], cfg['case'])
    years = days / 365.0

    # Clip rather than mask, so a panel never shows white holes where the
    # field is merely small; below-floor reads as the lowest colour.
    mesh = ax.pcolormesh(years, alt, np.clip(f.T, vmin, vmax),
                         norm=LogNorm(vmin=vmin, vmax=vmax),
                         cmap=CMAP, shading='auto')

    ax.set_ylim(0, alt_max)
    ax.set_xlim(*cfg['xlim_years'])
    ax.set_title(spec['label'].replace('\\texttt{', '').replace('}', ''),
                 fontsize=8.5, pad=4)

    # The eruption, named on every panel.
    ax.text(0.035, 0.93, CASE_LABEL.get(cfg['case'], cfg['case']),
            transform=ax.transAxes, ha='left', va='top', fontsize=7.5,
            color='white',
            bbox=dict(boxstyle='round,pad=0.25', fc='black', ec='none',
                      alpha=0.45))

    # An e-folding time is a statement about REMOVAL of injected material. On
    # the moist atmospheres the water field is dominated by the ambient
    # reservoir, so the same number would describe the background climate
    # rather than the eruption; it is reported for aerosol only.
    if cfg.get('annotate_efold') and FIELD != 'water':
        e = efold_days(days, f)
        txt = (f"$e$-fold {e:.0f} d" if np.isfinite(e)
               else "no $e$-fold in 6 yr")
        ax.text(0.965, 0.93, txt, transform=ax.transAxes, ha='right',
                va='top', fontsize=7.5, color='white',
                bbox=dict(boxstyle='round,pad=0.25', fc='black', ec='none',
                          alpha=0.45))
        print(f"  {spec['atm']:5s} e-folding {e:7.1f} d   peak {f.max():.3e}")
    else:
        print(f"  {spec['atm']:5s} peak {f.max():.3e}  "
              f"min-positive {f[f > 0].min():.3e}" if (f > 0).any()
              else f"  {spec['atm']:5s} identically zero")

    if k // 2 == 1:
        ax.set_xlabel('years since eruption')
    if k % 2 == 0:
        ax.set_ylabel('altitude [km]')

fig.subplots_adjust(right=0.86, hspace=0.24, wspace=0.08)
cax = fig.add_axes([0.885, 0.13, 0.022, 0.74])
cb = fig.colorbar(mesh, cax=cax)
cb.set_label(CBLABEL, fontsize=8.5)

stem = os.path.join(here, cfg['outfile_stem'])
for ext in ('pdf', 'eps'):
    fig.savefig(f'{stem}.{ext}', bbox_inches='tight')
    print(f'wrote {stem}.{ext}')
