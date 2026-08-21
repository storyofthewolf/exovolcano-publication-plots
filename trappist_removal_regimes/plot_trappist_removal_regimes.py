"""
plot_trappist_removal_regimes.py -- the four TRAPPIST-1e publication figures.

Selected from the 28-case ben2/hab2 sweep as the figures that carry the
science: what a large eruption does to a rocky M-dwarf planet, and how that
answer splits on whether the planet holds surface water.

  Fig. A  fig_trappist_removal      Headline. (a) global-mean AOD, both
                                    atmospheres, all templates; (b) peak AOD
                                    against injected SO2 across the magnitude
                                    ladder; (c) fraction of peak AOD surviving
                                    at year 6.
  Fig. B  fig_trappist_water_plume  (a) stratospheric excess H2O mass against
                                    each atmosphere's own control-variability
                                    floor; (b,c) height-time H2O anomaly.
  Fig. C  fig_trappist_sensitivity  What the eruption parameters control:
                                    (a) conversion timescale K, (b) day vs
                                    night injection, (c) nucleation barrier.
  Fig. D  fig_trappist_vertical     Height-time sulfate aerosol structure,
                                    which is what a transmission spectrum
                                    records.

SCOPE. ben2/hab2 only; ben1/hab1 were compiled with the wrong composition and
are being rerun. Both atmospheres here are 1 bar CO2 and differ in SURFACE
WATER alone. No observational reference line appears anywhere: there is no
measured TRAPPIST-1e volcanic aerosol or stratospheric water abundance to
compare against.

The 100x rung injects over 48 h rather than 24 h and is marked as such; the
injected mass is as labeled but the rate is halved, so it is not a clean
member of the magnitude ladder.

Configuration is in config_trappist_removal_regimes.yaml.
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
    'legend.fontsize': 7,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})

here = os.path.dirname(os.path.abspath(__file__))
parser = argparse.ArgumentParser()
parser.add_argument('--config', default='config_trappist_removal_regimes.yaml')
args = parser.parse_args()
cfg = yaml.safe_load(open(os.path.join(here, args.config)))

BASE_DIR = cfg['base_dir']
BATCHES = cfg['batches']
ATMOS = list(BATCHES.keys())
ATM_LABEL = cfg['atm_label']
G_CONST = cfg['g_const']
R_PLANET = cfg['r_planet']
AREA_PLANET = 4.0 * np.pi * R_PLANET ** 2
STRAT_P_MAX_PA = cfg['strat_p_max_pa']
KG_PER_TG = 1.0e9
AOD_FLOOR = cfg['aod_floor']

# One colour per atmosphere, held across all four figures so the dry/wet
# contrast is the same visual axis everywhere. Okabe-Ito.
ATM_COLOR = {'ben2': '#D55E00', 'hab2': '#0072B2'}
TPL_STYLE = {'tambora': '-', 'pinatubo': '--', 'hunga': ':'}


def dpath(atm, key, sub):
    return os.path.join(BASE_DIR, BATCHES[atm], f'exovolc_{atm}_{key}',
                        'data', sub)


def load2(atm, key, sub):
    p = dpath(atm, key, sub)
    if not os.path.exists(p):
        return None
    d = pd.read_csv(p, header=0)
    return d.iloc[:, 0].values.astype(float), d.iloc[:, 1].values.astype(float)


def load_profile(atm, key, sub):
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


def layer_dp(pres):
    n = len(pres)
    i = np.zeros(n + 1)
    i[1:-1] = 0.5 * (pres[:-1] + pres[1:])
    i[0] = 0.0
    i[-1] = pres[-1] + (pres[-1] - i[-2])
    return np.diff(i)


def strat_water_tg(atm, key):
    pr = load_profile(atm, key, cfg['subpath_q_profile'])
    if pr is None:
        return None
    days, pres, alt, q = pr
    w = layer_dp(pres) * AREA_PLANET / G_CONST
    m = pres < STRAT_P_MAX_PA
    return days, (q[:, m] * w[m]).sum(axis=1) / KG_PER_TG


def aod(atm, key):
    return load2(atm, key, cfg['subpath_aod'])


print('=' * 74)
print('TRAPPIST-1e publication figures -- ben2 / hab2')
print('  ben1/hab1 excluded: compiled composition wrong, rerun pending')
print('  Both atmospheres are 1 bar CO2; they differ in SURFACE WATER alone.')
print('=' * 74)


# ===========================================================================
# FIGURE A -- removal regimes (the headline)
# ===========================================================================
figA = plt.figure(figsize=(7.3, 4.6))
gs = figA.add_gridspec(2, 2, height_ratios=[1.25, 1.0], hspace=0.42,
                       wspace=0.28)
axA = figA.add_subplot(gs[0, :])
axB = figA.add_subplot(gs[1, 0])
axC = figA.add_subplot(gs[1, 1])

print('\n--- Fig. A(a): global-mean AOD, three templates -------------------')
for atm in ATMOS:
    for tpl in ('tambora', 'pinatubo', 'hunga'):
        got = aod(atm, tpl)
        if got is None:
            continue
        d, v = got
        axA.plot(d / 365.25, v, color=ATM_COLOR[atm], ls=TPL_STYLE[tpl],
                 lw=1.5 if tpl == 'tambora' else 1.0,
                 label=f'{atm} {tpl}')
        pk = int(np.nanargmax(v))
        tail = float(np.nanmean(v[-30:]))
        print(f'  {atm} {tpl:9s} peak {v[pk]:8.4f} at day {d[pk]:6.0f} | '
              f'yr6 {tail:9.3g} | retained {100*tail/v[pk]:5.1f}%')

axA.set_yscale('log')
axA.set_ylim(AOD_FLOOR, 1e1)
axA.set_xlim(0, cfg['xlim_days'][1] / 365.25)
axA.set_xlabel('Years since eruption')
axA.set_ylabel('Global-mean AOD, 550 nm')
axA.grid(alpha=0.25, lw=0.4)
axA.set_title('(a)  Aerosol lifetime splits on surface water', fontsize=8.5,
              loc='left')

# Legend: atmosphere by colour, template by dash.
from matplotlib.lines import Line2D
leg = [Line2D([], [], color=ATM_COLOR[a], lw=1.6, label=ATM_LABEL[a])
       for a in ATMOS]
leg += [Line2D([], [], color='0.35', ls=TPL_STYLE[t], lw=1.2,
               label=t.capitalize()) for t in ('tambora', 'pinatubo', 'hunga')]
axA.legend(handles=leg, loc='upper right', ncol=2, frameon=False,
           handlelength=2.2, fontsize=6.8)

# --- (b) peak AOD against injected SO2, magnitude ladder --------------------
print('\n--- Fig. A(b): magnitude ladder (Tambora, SO2) --------------------')
for atm in ATMOS:
    xs, ys, marks = [], [], []
    for step in cfg['ladder_tambora']:
        got = aod(atm, step['key'])
        if got is None:
            continue
        xs.append(step['so2_tg'])
        ys.append(float(np.nanmax(got[1])))
        marks.append(step['duration_h'])
    axB.plot(xs, ys, 'o-', color=ATM_COLOR[atm], ms=4, lw=1.3,
             label=ATM_LABEL[atm].split(':')[0])
    # Ring the 48 h rung, which is not a clean member of the ladder.
    for x, y, h in zip(xs, ys, marks):
        if h != 24:
            axB.plot([x], [y], 'o', ms=9, mfc='none',
                     mec=ATM_COLOR[atm], mew=1.0)
    for i in range(1, len(xs)):
        print(f'  {atm}: {xs[i-1]:7.0f} -> {xs[i]:7.0f} Tg SO2 : '
              f'AOD {ys[i-1]:8.4f} -> {ys[i]:8.4f}  '
              f'(x{ys[i]/ys[i-1]:.2f} for x{xs[i]/xs[i-1]:.0f} mass)')

# Linear-scaling guide, anchored on each atmosphere's 1x point.
for atm in ATMOS:
    got = aod(atm, 'tambora')
    if got is None:
        continue
    y0 = float(np.nanmax(got[1]))
    xg = np.array([60.0, 6000.0])
    axB.plot(xg, y0 * xg / 60.0, ls=(0, (1, 2)), lw=0.8, color=ATM_COLOR[atm],
             zorder=1)

axB.set_xscale('log')
axB.set_yscale('log')
axB.set_xlabel(r'Injected SO$_2$ [Tg]')
axB.set_ylabel('Peak global-mean AOD')
axB.grid(alpha=0.25, lw=0.4)
axB.set_title('(b)  Magnitude scaling', fontsize=8.5, loc='left')
axB.legend(loc='upper left', frameon=False)

# --- (c) fraction of peak surviving at year 6 ------------------------------
print('\n--- Fig. A(c): retained fraction at year 6 ------------------------')
keys = ['tambora', 'pinatubo', 'hunga', 'tambora_10x', 'tambora_100x']
klab = ['Tamb.', 'Pin.', 'Hunga', r'Tamb. 10$\times$', r'Tamb. 100$\times$']
w = 0.36
xpos = np.arange(len(keys))
for j, atm in enumerate(ATMOS):
    fr = []
    for k in keys:
        got = aod(atm, k)
        if got is None:
            fr.append(np.nan)
            continue
        v = got[1]
        fr.append(100.0 * float(np.nanmean(v[-30:])) / float(np.nanmax(v)))
    axC.bar(xpos + (j - 0.5) * w, fr, w, color=ATM_COLOR[atm],
            label=ATM_LABEL[atm].split(':')[0])
    for k, f in zip(keys, fr):
        print(f'  {atm} {k:14s} retained {f:6.2f}% of peak at yr 6')

axC.set_xticks(xpos)
axC.set_xticklabels(klab, fontsize=6.2, rotation=20, ha='right')
axC.set_ylabel('Peak AOD retained\nat year 6 [%]')
axC.grid(alpha=0.25, lw=0.4, axis='y')
axC.set_title('(c)  Aerosol survival', fontsize=8.5, loc='left')
axC.legend(loc='upper left', frameon=False)

for ext in ('pdf', 'eps'):
    figA.savefig(os.path.join(here, f'{cfg["outfile_removal"]}.{ext}'),
                 dpi=300, bbox_inches='tight')
plt.close(figA)
print(f'\n  wrote {cfg["outfile_removal"]}.pdf / .eps')


# ===========================================================================
# FIGURE B -- stratospheric water plume
# ===========================================================================
figB = plt.figure(figsize=(7.3, 4.4))
gsB = figB.add_gridspec(2, 2, height_ratios=[1.0, 1.0], hspace=0.45,
                        wspace=0.26)
axW = figB.add_subplot(gsB[0, :])
axQ = [figB.add_subplot(gsB[1, 0]), figB.add_subplot(gsB[1, 1])]

print('\n--- Fig. B(a): stratospheric excess H2O vs control noise floor ----')
lad_lw = {1: 1.0, 10: 1.4, 100: 1.8}
for atm in ATMOS:
    bg = strat_water_tg(atm, cfg['control'])
    if bg is None:
        continue
    bgd, bgv = bg
    sigma = float(np.nanstd(bgv))
    bgs = pd.Series(bgv, index=bgd)
    print(f'  {atm}: control mean {np.nanmean(bgv):10.4g} Tg, '
          f'1-sigma {sigma:10.4g} Tg  <- detection floor')
    # Every atmosphere gets its noise band. On hab2 that band IS the result:
    # the seasonal hydrological cycle, beating against an independent control
    # realization, is larger than any injection in the suite.
    axW.axhspan(-sigma, sigma, color=ATM_COLOR[atm], alpha=0.16, lw=0,
                zorder=1)
    for step in cfg['ladder_hunga']:
        got = strat_water_tg(atm, step['key'])
        if got is None:
            continue
        d, v = got
        ex = v - bgs.reindex(d).values
        snr = float(np.nanmax(ex)) / sigma
        # Draw the curve only where it means something. On hab2 the raw daily
        # difference is seasonal noise spanning +/-1e4 Tg, and plotting it
        # buries the dry-atmosphere signal under an oscillation that carries
        # no volcanic information. The shaded band already states hab2's
        # result -- that nothing clears the floor -- so the noise itself is
        # summarized rather than drawn.
        if snr >= 10.0:
            axW.plot(d / 365.25, ex, color=ATM_COLOR[atm],
                     lw=lad_lw[step['mult']], alpha=0.95, zorder=3)
        print(f'    {step["key"]:12s} injected {step["h2o_tg"]:8.0f} Tg -> '
              f'peak excess {np.nanmax(ex):10.4g} Tg, SNR {snr:8.1f}'
              f'{"" if snr >= 10 else "   <-- BELOW FLOOR, curve not drawn"}')

axW.set_yscale('symlog', linthresh=1.0)
axW.set_xlim(0, cfg['xlim_days'][1] / 365.25)
axW.set_xlabel('Years since eruption')
axW.set_ylabel(r'Stratospheric excess' '\n' r'H$_2$O mass [Tg]')
axW.grid(alpha=0.25, lw=0.4)
axW.set_title(r'(a)  Hunga-template water, 1$\times$/10$\times$/100$\times$ '
              r'(line weight); shaded band = that atmosphere'
              r"'s control 1$\sigma$", fontsize=7.5, loc='left')
legW = [Line2D([], [], color=ATM_COLOR[a], lw=1.6, label=ATM_LABEL[a])
        for a in ATMOS]
axW.legend(handles=legW, loc='upper right', frameon=False)

print('\n--- Fig. B(b,c): height-time H2O anomaly -------------------------')
wkey = cfg['water_case']
for ax, atm in zip(axQ, ATMOS):
    pr = load_profile(atm, wkey, cfg['subpath_q_profile'])
    ct = load_profile(atm, cfg['control'], cfg['subpath_q_profile'])
    if pr is None or ct is None:
        ax.set_visible(False)
        continue
    d, pres, alt, q = pr
    cd, _, _, qc = ct
    n = min(len(d), len(cd))
    dq = q[:n] - qc[:n]
    akm = alt / 1000.0
    # Colour limits from the stratosphere only. On the aquaplanet the
    # boundary layer carries an anomaly orders of magnitude larger, driven by
    # ordinary weather in an independent control realization and unrelated to
    # the injection; scaling on the full column erases both plumes.
    amax = float(np.nanpercentile(np.abs(dq[:, akm > 15.0]), 99.5))
    pcm = ax.pcolormesh(d[:n] / 365.25, akm, dq.T,
                        norm=TwoSlopeNorm(vmin=-amax, vcenter=0.0, vmax=amax),
                        cmap='RdBu_r', shading='auto', rasterized=True)
    ax.set_ylim(0, cfg['contour_alt_max_km'])
    ax.set_xlim(0, cfg['xlim_days_short'][1] / 365.25)
    ax.set_xlabel('Years since eruption')
    ax.set_title(f'({"bc"[ATMOS.index(atm)]})  {ATM_LABEL[atm]}',
                 fontsize=8, loc='left')
    cb = figB.colorbar(pcm, ax=ax, pad=0.02, aspect=16)
    cb.set_label(r'$\Delta$Q [kg kg$^{-1}$]', fontsize=7)
    cb.ax.tick_params(labelsize=6)
    sm = akm > 15.0
    dqs = dq[:, sm]
    it, iz = np.unravel_index(int(np.nanargmax(dqs)), dqs.shape)
    print(f'  {atm} {wkey}: peak stratospheric dQ {dqs[it, iz]:.4g} kg/kg '
          f'at day {d[it]:.0f}, {akm[sm][iz]:.1f} km')
axQ[0].set_ylabel('Altitude [km]')

for ext in ('pdf', 'eps'):
    figB.savefig(os.path.join(here, f'{cfg["outfile_water"]}.{ext}'),
                 dpi=150, bbox_inches='tight')
plt.close(figB)
print(f'\n  wrote {cfg["outfile_water"]}.pdf / .eps')


# ===========================================================================
# FIGURE C -- what the eruption parameters control
# ===========================================================================
figC, axC3 = plt.subplots(1, 3, figsize=(7.3, 2.5))

print('\n--- Fig. C(a): conversion timescale K ----------------------------')
for atm in ATMOS:
    ks, ps = [], []
    for step in cfg['k_ladder']:
        got = aod(atm, step['key'])
        if got is None:
            continue
        ks.append(step['k_days'])
        ps.append(float(np.nanmax(got[1])))
    axC3[0].plot(ks, ps, 'o-', color=ATM_COLOR[atm], ms=4, lw=1.3,
                 label=ATM_LABEL[atm].split(':')[0])
    print(f'  {atm}: ' + ', '.join(f'K={k:g}d -> {p:.4f}'
                                   for k, p in zip(ks, ps)))
axC3[0].set_xscale('log')
axC3[0].set_yscale('log')
axC3[0].set_xlabel(r'Conversion timescale $K$ [d]')
axC3[0].set_ylabel('Peak global-mean AOD')
axC3[0].grid(alpha=0.25, lw=0.4)
axC3[0].set_title('(a)  Chemistry', fontsize=8.5, loc='left')
axC3[0].legend(loc='lower left', frameon=False)

print('\n--- Fig. C(b): dayside vs nightside injection --------------------')
pairs = cfg['daynight_pairs']
xp = np.arange(len(pairs))
wb = 0.36
for j, atm in enumerate(ATMOS):
    rat = []
    for pr_ in pairs:
        dgot, ngot = aod(atm, pr_['day']), aod(atm, pr_['night'])
        if dgot is None or ngot is None:
            rat.append(np.nan)
            continue
        r = float(np.nanmax(ngot[1])) / float(np.nanmax(dgot[1]))
        rat.append(r)
        print(f'  {atm} {pr_["label"]:9s} night/day peak AOD = {r:.3f}')
    axC3[1].bar(xp + (j - 0.5) * wb, rat, wb, color=ATM_COLOR[atm])
axC3[1].axhline(1.0, color='k', lw=0.8, ls=':')
axC3[1].set_xticks(xp)
axC3[1].set_xticklabels([p['label'] for p in pairs], fontsize=7)
axC3[1].set_ylabel('Peak AOD ratio,\nnight / day')
axC3[1].set_ylim(0, 1.6)
axC3[1].grid(alpha=0.25, lw=0.4, axis='y')
axC3[1].set_title('(b)  Injection hemisphere', fontsize=8.5, loc='left')

print('\n--- Fig. C(c): nucleation barrier S_crit -------------------------')
sp = cfg['scrit_pair']
for j, atm in enumerate(ATMOS):
    b, p_ = aod(atm, sp['base']), aod(atm, sp['probe'])
    if b is None or p_ is None:
        continue
    d0, v0 = b
    d1, v1 = p_
    axC3[2].plot(d0 / 365.25, v0, color=ATM_COLOR[atm], lw=1.4)
    axC3[2].plot(d1 / 365.25, v1, color=ATM_COLOR[atm], lw=1.0, ls='--')
    print(f'  {atm}: baseline peak {np.nanmax(v0):.4f}, '
          f'S_crit=5 peak {np.nanmax(v1):.4f}, '
          f'ratio {np.nanmax(v1)/np.nanmax(v0):.4f}')
axC3[2].set_yscale('log')
axC3[2].set_ylim(AOD_FLOOR, 1e1)
axC3[2].set_xlim(0, cfg['xlim_days_short'][1] / 365.25)
axC3[2].set_xlabel('Years since eruption')
axC3[2].set_ylabel('Global-mean AOD')
axC3[2].grid(alpha=0.25, lw=0.4)
axC3[2].set_title('(c)  Nucleation barrier', fontsize=8.5, loc='left')
axC3[2].legend(handles=[
    Line2D([], [], color='0.35', lw=1.4, label='baseline'),
    Line2D([], [], color='0.35', lw=1.0, ls='--', label=sp['label'])],
    loc='lower left', frameon=False)

figC.tight_layout()
for ext in ('pdf', 'eps'):
    figC.savefig(os.path.join(here, f'{cfg["outfile_sensitivity"]}.{ext}'),
                 dpi=300, bbox_inches='tight')
plt.close(figC)
print(f'\n  wrote {cfg["outfile_sensitivity"]}.pdf / .eps')


# ===========================================================================
# FIGURE D -- vertical aerosol structure
# ===========================================================================
figD, axD = plt.subplots(1, 2, figsize=(7.3, 2.7), sharey=True)
vkey = cfg['vertical_case']

print('\n--- Fig. D: height-time sulfate aerosol --------------------------')
hz = {}
for atm in ATMOS:
    pr = load_profile(atm, vkey, cfg['subpath_hz_profile'])
    if pr is not None:
        hz[atm] = pr

if hz:
    pos = np.concatenate([d[3][d[3] > 0].ravel() for d in hz.values()])
    vmax = float(np.nanpercentile(pos, 99.9))
    vmin = vmax * 1e-5
    for ax, atm in zip(axD, ATMOS):
        if atm not in hz:
            ax.set_visible(False)
            continue
        d, pres, alt, f = hz[atm]
        akm = alt / 1000.0
        pcm = ax.pcolormesh(d / 365.25, akm,
                            np.ma.masked_less_equal(f.T, 0.0),
                            norm=LogNorm(vmin=vmin, vmax=vmax),
                            cmap='inferno', shading='auto', rasterized=True)
        ax.set_ylim(0, cfg['contour_alt_max_km'])
        ax.set_xlim(0, cfg['xlim_days'][1] / 365.25)
        ax.set_xlabel('Years since eruption')
        ax.set_title(f'({"ab"[ATMOS.index(atm)]})  {ATM_LABEL[atm]}',
                     fontsize=8, loc='left')
        it, iz = np.unravel_index(int(np.nanargmax(f)), f.shape)
        print(f'  {atm} {vkey}: peak {f[it, iz]:.4g} kg/m3 at '
              f'day {d[it]:.0f}, {akm[iz]:.1f} km')
    axD[0].set_ylabel('Altitude [km]')
    cb = figD.colorbar(pcm, ax=axD, pad=0.02, aspect=22)
    cb.set_label(r'Sulfate aerosol [kg m$^{-3}$]', fontsize=7.5)
    cb.ax.tick_params(labelsize=6.5)

for ext in ('pdf', 'eps'):
    figD.savefig(os.path.join(here, f'{cfg["outfile_vertical"]}.{ext}'),
                 dpi=150, bbox_inches='tight')
plt.close(figD)
print(f'\n  wrote {cfg["outfile_vertical"]}.pdf / .eps')
print('\nDone.')
