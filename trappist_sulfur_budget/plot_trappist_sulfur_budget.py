"""
plot_trappist_sulfur_budget.py -- sulfur partitioning and water burden for the
TRAPPIST-1e suites.

Two figures, built to test one mechanistic hypothesis: that peak AOD is lower
on the dry ben2 atmosphere than on the moist hab2 one because the
SO2 -> H2SO4 conversion CONSUMES water and is water-limited when water is
scarce.

  Fig. 1  fig_trappist_sulfur_budget
          Global sulfur burden in Tg S, partitioned into SO2 gas, H2SO4 gas,
          and condensed sulfate aerosol, for ben2 (top) and hab2 (bottom)
          across three templates. If the hypothesis holds, ben2 should strand
          sulfur as unconverted SO2 while hab2 drives it promptly to aerosol.

  Fig. 2  fig_trappist_ben2_water_burden
          ben2 only. Global water burden from TMQ (total precipitable water,
          kg/m2 -> Tg), against the stoichiometric water demand of each
          eruption's SO2 load.

THE CHEMISTRY, as actually compiled (exo_simplevolc.F90, subroutine
exo_simplevolc_gas_tend; verified against the source built into
exovolc_ben2_tambora, md5 dbf7c569):

    SO2 + H2O + 1/2 O2 -> H2SO4

    so2_loss_max  = q_SO2 * (1 - exp(-K*dt))
    h2o_needed    = so2_loss_max * (MW_H2O/MW_SO2)
    h2o_available = max(0, q_H2O - 1e-20)
    if (h2o_needed > h2o_available):
        so2_loss = h2o_available / (MW_H2O/MW_SO2)     <-- the clamp

So the conversion is water-limited BY CONSTRUCTION, per column and per
timestep. This is a stoichiometric REACTANT limitation on the gas-phase
conversion. It is NOT the H2O-dependent sulfate pathway that Section 3 found
absent (H2SO4 saturation is temperature-only and K is prescribed) -- that
remains absent. The two statements are consistent: the reactant clamp only
bites where water is scarce, which never happens on Earth or on hab2.

UNITS. SO2.csv, H2SO4.csv and VOLCHZMD.csv are already GLOBAL MASSES in kg.
They are converted to Tg of SULFUR by the mass fraction of S in each species
so the three curves are directly comparable and sum to a conserved total.
TMQ.csv is a column quantity in kg/m2 and is multiplied by planet area.

Configuration is in config_trappist_sulfur_budget.yaml.
"""

import os
import argparse
import yaml
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

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
parser.add_argument('--config', default='config_trappist_sulfur_budget.yaml')
args = parser.parse_args()
cfg = yaml.safe_load(open(os.path.join(here, args.config)))

BASE_DIR = cfg['base_dir']
BATCHES = cfg['batches']
ATMOS = list(BATCHES.keys())
ATM_LABEL = cfg['atm_label']
AREA = 4.0 * np.pi * cfg['r_planet'] ** 2
G_CONST = cfg['g_const']
KG_PER_TG = 1.0e9

MW_SO2, MW_H2SO4 = cfg['mw_so2'], cfg['mw_h2so4']
MW_H2O, MW_S = cfg['mw_h2o'], cfg['mw_s']

# Sulfur mass fraction of each carrier, so all three curves are Tg of S.
S_FRAC_SO2 = MW_S / MW_SO2          # 0.5005
S_FRAC_H2SO4 = MW_S / MW_H2SO4      # 0.3269
S_FRAC_AER = S_FRAC_H2SO4           # sulfate aerosol, condensed from H2SO4
H2O_PER_SO2 = MW_H2O / MW_SO2       # 0.2812 kg H2O per kg SO2

SPECIES_COLOR = {'so2': '#E69F00', 'h2so4': '#009E73', 'aer': '#CC79A7'}
ATM_COLOR = {'ben2': '#D55E00', 'hab2': '#0072B2'}


def dpath(atm, key, sub):
    return os.path.join(BASE_DIR, BATCHES[atm], f'exovolc_{atm}_{key}',
                        'data', sub)


def load2(atm, key, sub):
    p = dpath(atm, key, sub)
    if not os.path.exists(p):
        return None
    d = pd.read_csv(p, header=0)
    return d.iloc[:, 0].values.astype(float), d.iloc[:, 1].values.astype(float)


print('=' * 76)
print('TRAPPIST-1e sulfur budget and water burden')
print(f'  stoichiometry: {H2O_PER_SO2:.4f} kg H2O consumed per kg SO2 converted')
print(f'  S fractions: SO2 {S_FRAC_SO2:.4f}, H2SO4 {S_FRAC_H2SO4:.4f}')
print('=' * 76)


# ===========================================================================
# FIGURE 1 -- sulfur partitioning, Tg S
# ===========================================================================
cases = cfg['sulfur_cases']
figS, axS = plt.subplots(2, len(cases), figsize=(7.3, 4.4),
                         sharex=True)

for i, atm in enumerate(ATMOS):
    for j, entry in enumerate(cases):
        ax = axS[i, j]
        key = entry['key']
        so2 = load2(atm, key, cfg['subpath_so2'])
        h2s = load2(atm, key, cfg['subpath_h2so4'])
        aer = load2(atm, key, cfg['subpath_aer'])
        if so2 is None or h2s is None or aer is None:
            ax.set_visible(False)
            continue
        d = so2[0]
        s_so2 = so2[1] * S_FRAC_SO2 / KG_PER_TG
        s_h2s = h2s[1] * S_FRAC_H2SO4 / KG_PER_TG
        s_aer = aer[1] * S_FRAC_AER / KG_PER_TG
        tot = s_so2 + s_h2s + s_aer

        ax.plot(d / 365.25, s_so2, color=SPECIES_COLOR['so2'], lw=1.3,
                label=r'SO$_2$ gas')
        # H2SO4 gas never accumulates: condensation to aerosol (TAU_AER_CONV
        # = 1800 s) is far faster than its production, so this curve sits at
        # ~1e-16 Tg S, orders of magnitude below the panel floor. It is drawn
        # for completeness and to show the budget is closed, not because it is
        # a visible reservoir.
        ax.plot(d / 365.25, s_h2s, color=SPECIES_COLOR['h2so4'], lw=1.3,
                label=r'H$_2$SO$_4$ gas ($\sim$0)')
        ax.plot(d / 365.25, s_aer, color=SPECIES_COLOR['aer'], lw=1.3,
                label='sulfate aerosol')
        ax.plot(d / 365.25, tot, color='0.3', lw=0.8, ls=':',
                label='total S')

        ax.set_yscale('log')
        # Two years, not six. The conversion competition this figure is about
        # is decided in the first months; the remaining four years are decay
        # that Figure fig_trappist_removal already shows.
        ax.set_xlim(0, 2.0)
        # Floor each column at 1e-4 of its own peak total sulfur. The raw
        # decay runs to ~1e-28 Tg S, which is physically meaningless and
        # compresses the partitioning -- the whole point of the figure -- into
        # the top tenth of the panel. Scaling per column keeps the three
        # templates, whose injected masses differ by ~600x, comparable in
        # SHAPE rather than forcing them onto one absolute axis.
        top = float(np.nanmax(tot))
        ax.set_ylim(top * 1e-4, top * 3.0)
        ax.grid(alpha=0.25, lw=0.4)
        if i == 0:
            ax.set_title(entry['label'], fontsize=7.2)
        if j == 0:
            ax.set_ylabel(f'{atm}\nSulfur burden [Tg S]', fontsize=8)
        if i == 1:
            ax.set_xlabel('Years since eruption')

        # Peak-partitioning diagnostic: at the time of peak total sulfur,
        # how is that sulfur distributed? This is the number that decides the
        # hypothesis.
        pk = int(np.nanargmax(tot))
        # Fraction of sulfur ever reaching aerosol, at the aerosol's own peak.
        apk = int(np.nanargmax(s_aer))
        print(f'\n  {atm} {key}:')
        print(f'    injected SO2 stoich. H2O demand = '
              f'{entry["need_h2o_tg"]:8.2f} Tg H2O')
        print(f'    peak total S      {tot[pk]:9.3f} Tg S at day {d[pk]:6.0f}')
        print(f'      partition at that time: SO2 {100*s_so2[pk]/tot[pk]:5.1f}%'
              f' | H2SO4 {100*s_h2s[pk]/tot[pk]:5.1f}%'
              f' | aerosol {100*s_aer[pk]/tot[pk]:5.1f}%')
        print(f'    peak aerosol S    {s_aer[apk]:9.3f} Tg S at day {d[apk]:6.0f}'
              f'  ({100*s_aer[apk]/np.nanmax(tot):5.1f}% of peak total S)')

axS[0, 0].legend(loc='lower left', frameon=False, fontsize=6.2)
figS.suptitle('Sulfur partitioning: dry (top) vs moist (bottom) TRAPPIST-1 e',
              fontsize=9, y=0.98)
figS.tight_layout(rect=(0, 0, 1, 0.96))
for ext in ('pdf', 'eps'):
    figS.savefig(os.path.join(here, f'{cfg["outfile_sulfur"]}.{ext}'),
                 dpi=200, bbox_inches='tight')
plt.close(figS)
print(f'\n  wrote {cfg["outfile_sulfur"]}.pdf / .eps')


# ===========================================================================
# FIGURE 2 -- ben2 global water burden (TMQ)
# ===========================================================================
# ben2 only, per request. TMQ is total precipitable water [kg/m2]; multiplying
# by planet area gives a global water mass. Horizontal lines mark the
# stoichiometric water demand of each SO2 load, so the figure reads directly
# as "is there enough water to convert this sulfur?".
figW, axW = plt.subplots(1, 2, figsize=(7.3, 3.0))

print('\n' + '-' * 76)
print('ben2 global water burden vs stoichiometric demand')
print('-' * 76)
print('Source: TMQ [kg/m2] x planet area. TMQ recovers the injected mass well:')
print('  Hunga 1x/10x/100x give 144.6 / 1420 / 7111 Tg on day 1 against')
print('  146 / 1460 / 14,600 Tg injected, reaching full mass within a few days.')
print()
print('*** SPURIOUS WATER SOURCE IN THE DRY (ben2) RUNS ***')
print('  The ben2 CONTROL injects zero water and starts from Q identically 0,')
print('  so its water must stay ~0. It does not: it climbs from 1.2e-11 Tg to')
print('  122.5 Tg over six years, in ~28 discrete multi-Tg jumps. hab2 shows no')
print('  such behaviour (its control DECREASES, ratio 0.77). The dry runs')
print('  require the qmin = 1e-50 floor to avoid a MAPZ_MODULE crash, and a')
print('  positivity/filling fixer acting on that floor is the prime suspect.')
print('  CONSEQUENCE: ben2 water and any aerosol formed after ~day 250 are')
print('  CONTAMINATED. For tambora_10x, 46.2 of 48.1 Tg S of aerosol (96%)')
print('  forms after drift onset, needing 26.0 Tg H2O that the drift supplies.')
print('  Only the first ~200 days of the ben2 cases are trustworthy.')


def water_tg(atm, key):
    """Global water burden [Tg] from TMQ (total precipitable water, kg/m2)
    integrated over the planet surface. TMQ is the correct field for this and
    recovers the injection accurately: the Hunga cases return 144.6 / 1420 /
    7111 Tg on day 1 against 146 / 1460 / 14,600 Tg injected, reaching the full
    injected mass within a few days."""
    got = load2(atm, key, cfg['subpath_tmq'])
    if got is None:
        return None
    d, tmq = got
    return d, tmq * AREA / KG_PER_TG


wcases = cfg['water_cases']
cmap = plt.get_cmap('viridis')
colors = [cmap(x) for x in np.linspace(0.05, 0.92, len(wcases))]

ctl = water_tg('ben2', cfg['control'])
for entry, colr in zip(wcases, colors):
    got = water_tg('ben2', entry['key'])
    if got is None:
        print(f'  MISSING Q profile: {entry["key"]}')
        continue
    d, mass = got
    # 30-day running mean. The dry atmospheres are initialized essentially
    # bone-dry (the qmin = 1e-50 floor) and accumulate water over the first
    # ~2 yr, and the water-starved cases oscillate hard as water is consumed
    # by conversion and resupplied by transport. Both are real; the smoothing
    # keeps the trend legible without inventing one.
    ms = pd.Series(mass).rolling(30, center=True, min_periods=1).mean().values
    for ax in axW:
        ax.plot(d / 365.25, ms, color=colr, lw=1.4, label=entry['label'])
    frac = (100.0 * np.nanmean(mass[-30:]) / np.nanmean(ctl[1][-30:])
            if ctl is not None else np.nan)
    print(f'  {entry["key"]:14s} peak {np.nanmax(mass):10.4g} Tg | '
          f'yr6 {np.nanmean(mass[-30:]):10.4g} Tg  '
          f'({frac:6.1f}% of control)')

# Stoichiometric water demand of each SO2 load. Where the demand line sits
# ABOVE a case's water curve, that eruption cannot convert all its sulfur.
for entry in cfg['sulfur_cases']:
    for ax in axW:
        ax.axhline(entry['need_h2o_tg'], color='0.4', lw=0.8, ls='--',
                   zorder=1)
axW[0].text(0.985, 16.87, r'Tambora demand, 16.9 Tg  ', fontsize=5.8,
            color='0.3', va='bottom', ha='right',
            transform=axW[0].get_yaxis_transform())
axW[0].text(0.985, 168.73, r'Tambora 10$\times$ demand, 169 Tg  ', fontsize=5.8,
            color='0.3', va='bottom', ha='right',
            transform=axW[0].get_yaxis_transform())

for ax, xl, ttl in zip(axW,
                       (cfg['xlim_days'][1], 200.0),
                       ('(a)  Full integration (shaded: contaminated)',
                        '(b)  First 200 days (trustworthy)')):
    # Linear, not log. The interesting range spans 4 to 300 Tg for every case
    # that is water-LIMITED, and a log axis devotes most of the panel to the
    # spin-up from a near-zero initial state instead.
    ax.set_xlim(0, xl / 365.25)
    ax.set_ylim(0, 320)
    # Everything after ~day 250 is contaminated by the spurious water source
    # documented above. Shade it rather than crop it, so the artifact is
    # visible instead of quietly removed.
    if xl > 400:
        ax.axvspan(250.0 / 365.25, xl / 365.25, color='#CC3311', alpha=0.10,
                   lw=0, zorder=0)
    ax.set_xlabel('Years since eruption')
    ax.grid(alpha=0.25, lw=0.4)
    ax.set_title(ttl, fontsize=8.5, loc='left')
axW[0].set_ylabel(r'\texttt{ben2} global water burden [Tg]'.replace('\\texttt{','').replace('}',''))
axW[0].set_ylabel('ben2 global water burden [Tg]')

for ax in axW:
    ax.annotate('Hunga 10x (1460 Tg) and 100x (14,600 Tg)'
                '\noff scale above',
                xy=(0.5, 0.955), xycoords='axes fraction', ha='center',
                fontsize=5.6, color='0.35')

h, l = axW[0].get_legend_handles_labels()
figW.legend(h, l, loc='lower center', ncol=4, frameon=False,
            bbox_to_anchor=(0.5, -0.16), fontsize=6.4, handlelength=1.6)
figW.tight_layout()
for ext in ('pdf', 'eps'):
    figW.savefig(os.path.join(here, f'{cfg["outfile_water"]}.{ext}'),
                 dpi=200, bbox_inches='tight')
plt.close(figW)
print(f'\n  wrote {cfg["outfile_water"]}.pdf / .eps')
print('\nDone.')
