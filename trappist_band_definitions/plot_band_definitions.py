"""
plot_band_definitions.py -- show where the fixed band wavelengths sit.

A diagnostic, not a manuscript figure. Its job is to let the band definitions
in exovolcano-spectra/scripts/band_depth_fixed.py be checked by eye against
the spectra they are applied to: every core and shoulder is drawn on the
spectrum of the case it was located from, and on the aerosol-rich case where
it is most at risk of being contaminated.

Each panel shows one band. Solid vertical line = the core wavelength; dashed
= the two shoulders; the straight segment joining the shoulder points is the
baseline the depth is measured against, so the vertical gap between it and the
spectrum at the core IS the reported band depth, drawn as an arrow.

Two epochs per panel: day 0 (grey) and the epoch of largest response
(coloured), so a reader can see both the pre-eruption state and what the
eruption did to it.

    python plot_band_definitions.py
"""

import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, '/Users/wolfe/Desktop/projects/volcanos/exovolcano-spectra')
from scripts.band_depth_fixed import (BANDS, REFERENCE, bin_to_resolution,
                                      band_depth, at, ROOT, SUFFIX)

plt.rcParams.update({
    'font.family': 'serif', 'font.size': 8, 'axes.labelsize': 9,
    'xtick.labelsize': 8, 'ytick.labelsize': 8, 'legend.fontsize': 7,
    'pdf.fonttype': 42, 'ps.fonttype': 42,
})

here = os.path.dirname(os.path.abspath(__file__))

# Each band is shown twice: on its own reference case, and on the case where
# the aerosol veil is heaviest, which is where a badly placed shoulder does
# the most damage.
STRESS = "hab1/tambora_100x"


def load(case_key):
    atm, case = case_key.split('/')
    z = np.load(os.path.join(ROOT, f"exovolc_{atm}_{case}" + SUFFIX),
                allow_pickle=True)
    wl, dep = bin_to_resolution(z['wavelength_um'], z['depth'] * 1e6)
    return z['days'], wl, dep


names = list(BANDS)
fig, axes = plt.subplots(len(names), 2, figsize=(11.2, 2.35 * len(names)))

for r, name in enumerate(names):
    lb, lc, lr = BANDS[name]
    ref_case, ref_day, *_ = REFERENCE[name]
    for c, case_key in enumerate([ref_case, STRESS]):
        ax = axes[r, c]
        days, wl, dep = load(case_key)
        # Epoch of largest |change| in this band's depth.
        d0 = band_depth(wl, dep[0], lb, lc, lr)
        best = max(range(1, len(days)),
                   key=lambda j: abs(band_depth(wl, dep[j], lb, lc, lr) - d0))

        pad = 0.35 * (lr - lb)
        m = (wl >= lb - pad) & (wl <= lr + pad)
        ax.plot(wl[m], dep[0][m], color='0.6', lw=0.9, label='day 0')
        ax.plot(wl[m], dep[best][m], color='#b2182b', lw=1.0,
                label=f'day {int(days[best])}')

        for j, col in [(0, '0.6'), (best, '#b2182b')]:
            yb, yr = at(wl, dep[j], lb), at(wl, dep[j], lr)
            # The baseline: a straight line between the two shoulder points.
            ax.plot([lb, lr], [yb, yr], color=col, lw=0.8, ls=':')
            frac = (lc - lb) / (lr - lb)
            base = yb + (yr - yb) * frac
            top = at(wl, dep[j], lc)
            ax.annotate('', xy=(lc, top), xytext=(lc, base),
                        arrowprops=dict(arrowstyle='<->', color=col, lw=0.9))

        ax.axvline(lc, color='0.25', lw=0.8)
        for lam in (lb, lr):
            ax.axvline(lam, color='0.45', lw=0.8, ls='--')
        ax.set_xlim(lb - pad, lr + pad)
        ax.set_title(f"{name}   {case_key}"
                     + ("   [reference]" if c == 0 else "   [aerosol-rich]"),
                     fontsize=8.5)
        if c == 0:
            ax.set_ylabel('transit depth [ppm]')
            ax.legend(frameon=False, loc='upper left')
        if r == len(names) - 1:
            ax.set_xlabel(r'wavelength [$\mu$m]')
        dd = band_depth(wl, dep[best], lb, lc, lr) - d0
        ax.text(0.985, 0.05, f"$D_0$={d0:.1f}   $\\Delta D$={dd:+.1f} ppm",
                transform=ax.transAxes, ha='right', va='bottom', fontsize=7.5)

fig.suptitle('Fixed band wavelengths: core (solid), shoulders (dashed), '
             'baseline (dotted), measured depth (arrow)', fontsize=9.5,
             y=0.997)
fig.tight_layout(rect=[0, 0, 1, 0.985])
for ext in ('pdf', 'eps'):
    fig.savefig(os.path.join(here, f'fig_band_definitions.{ext}'),
                bbox_inches='tight')
    print(f'wrote fig_band_definitions.{ext}')
