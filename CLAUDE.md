# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Purpose

Bespoke publication-quality plotting scripts for:
> Wolf, E.T. et al. (2026) "Simulating volcanic eruptions in diverse planetary context with ExoCAM: Descriptions and baselines" *(In preparation)*

These are **not** general-purpose routines — each figure script is tailored for a specific plot.

## Data pipeline

CSV input files are produced by [exovolcano-analysis](https://github.com/storyofthewolf/exovolcano-analysis) from ExoCAM `.h1.` time-series output.

**Current on-disk layout** — note the eruption directory, and that the case name is NOT repeated:

```
remote_analysis/
  <eruption>/            # tambora | pinatubo | hunga | hunga_phase2
    <case_name>/         # e.g. exovolc_tambora_fid
      data/
        aod/aod_550nm_band.csv
        aod/aod_zonal_550nm_band.csv
        scalar/{SO2,H2SO4,VOLCHZMD,Q,TS,TMQ,TGCLDLWP}.csv
        profiles/{T,Q,SO2,H2SO4,VOLCHZMD}.csv
        zonal/<VAR>_day<NNNN.NN>.csv
        zonal_level/<VAR>_<NN>km.csv     # e.g. Q_25km.csv, VOLCHZMD_20km.csv
        lonlat/aod_day<NNNN.NN>.csv      # full (lat, lon) AOD field
```

`base_dir` in each figure's YAML points to the local copy of `remote_analysis` (currently `/Users/wolfe/Desktop/projects/volcanos/remote_analysis`).

The `zonal_level/` and `lonlat/` diagnostics come from separate scripts in
exovolcano-analysis (`zonal_level_timeseries.py` and `lonlat_aod.py`), not from
`run_time_series.py`, and exist only for the cases they were run on. `lonlat`
snapshots exist on a fixed day set (1, 4, 10, 30, 107).

### The Hunga suite spans two roots

The Hunga cases live under **two** sibling directories, `hunga` (phase 1) and
`hunga_phase2`, and `hunga_phase2` sits at the same level as the eruption
directories rather than inside `hunga`. The phase-2 extension added the SO$_2$
linearity axis, the Raikoke-equivalent calibration control, and the current
manuscript fiducial `exovolc_hunga_r0.4_so2_1.0tg`.

Hunga figure configs therefore carry a `base_dirs` mapping (`phase1`/`phase2`)
instead of a single `base_dir`, and each case entry names the phase it belongs
to. A case named without its phase will not resolve.

### KNOWN TRAP: `pub_data` encodes the OLD layout

`pub_data.find_csvs()` and `find_csvs_list()` both resolve paths as
`base_dir/<case>/data/<case>/<subpath>` — the *superseded* convention, with no
eruption directory and the case name repeated. They no longer resolve against
the layout above. The CSV basenames also changed (`aod_550nm_band.csv`, not
`aod_0p550um_mie.csv`).

Consequently four figure folders are **stale**: `aod_timeseries` and
`aod_zonal_contour` still import `find_csvs_list`; `aod_twopanel` hardcodes the
old path and the old CSV names; and `sulfur_burden` uses the same superseded
`data/<case>/` convention and an obsolete fiducial
(`exovolc_tambora_k35d_r0.5`). All four will report every case as missing until
repointed. Treat them as style references, not working scripts.

`aod_twopanel` is superseded by `tambora_validation`, and `sulfur_burden` by
`tambora_sulfur_budget`; prefer the replacements over repointing the originals.

The nine figure folders added in 2026-07 and later (the four Tambora/Pinatubo
figures, plus the five Hunga figures) build paths directly and are the pattern
to copy — see the directory listing below.

## Running a figure script

From inside the figure folder:
```bash
cd tambora_validation
python plot_tambora_validation.py
```

Or from the project root:
```bash
python tambora_validation/plot_tambora_validation.py
```

Every script also accepts `--config <file>` (resolved relative to the script's
own directory) to run an alternate YAML.

Output is written into the figure's own folder. Current figures set an
`outfile_stem` / `outfile_base` in the YAML and save both `.pdf` and `.eps` from
it; the stale folders instead set a single `outfile` carrying its own extension.
Change that key to version outputs.

## Config split: YAML vs Python

**Edit the YAML** (`config_<name>.yaml`) to change: which cases are plotted,
`base_dir`/`base_dirs`, `subpath`, the output stem, `eruption_date`, `xlim`,
reference values and their labels, and unit-conversion constants.

Legend labels for highlighted cases belong in the YAML next to the case they
name, not hardcoded in the script — a previous fiducial's label survived a
change of fiducial that way and briefly mislabelled a figure.

**Edit the Python** (`plot_<name>.py`) to change: plot type, axes structure, data transforms, legend style, figure size. These are fixed once a figure is established.

## Directory structure

One folder per figure, each holding a `plot_<name>.py`, a `config_<name>.yaml`,
and its gitignored `.pdf`/`.eps` output.

```
pub_data.py                     ← legacy path helper; see the trap above

tambora_validation/             ← current pattern (2026-07)
tambora_sulfur_budget/
pinatubo_validation/
pinatubo_sensitivity/

hunga_validation/               ← current pattern (2026-07/08); two-root configs
hunga_schoeberl_compare/
hunga_zonal_level/
hunga_water_persistence/
hunga_plume_linearity/

trappist_removal_regimes/       ← TRAPPIST-1 e publication figures (2026-08)
trappist_sulfur_budget/         ← TRAPPIST-1 e mechanism figures (2026-08)
trappist_diagnostics/           ← exploratory quad; NOT a manuscript source

aod_timeseries/                 ← stale, old path convention
aod_zonal_contour/              ← stale
aod_twopanel/                   ← stale; superseded by tambora_validation
sulfur_burden/                  ← stale; superseded by tambora_sulfur_budget
```

What each current figure is:

| Folder | Figure |
| --- | --- |
| `tambora_validation` | Headline VolMIP validation: global-mean AOD with the 16-case envelope over the Clyne et al. (2021) peak-SAOD band, plus a latitude–time zonal Hovmöller. |
| `tambora_sulfur_budget` | SO₂ / H₂SO₄ / condensed-sulfate burdens in Tg S against the injected 30 Tg S line, over a shared-x global-mean AOD panel. |
| `pinatubo_validation` | Same two-panel structure, 18-case envelope, with Mills et al. (2016) and GloSSAC reference lines. |
| `pinatubo_sensitivity` | 2×2 shared-y sweeps: conversion timescale K, effective radius, SO₂ mass, injection pressure. |
| `hunga_validation` | Headline Hunga validation, 28-case envelope across both phases. Plots **no** observational reference line — see the scope notes below. |
| `hunga_schoeberl_compare` | First-year global-mean AOD against Schoeberl et al. (2023)'s 0.018 global-mean peak, plus a zonal contour. |
| `hunga_zonal_level` | Latitude–time fields at fixed altitude on the basis of Schoeberl et al. (2024) Fig. 3: sulfate at 20 km, water-vapor anomaly at 25 km. |
| `hunga_water_persistence` | Stratospheric excess H₂O mass (p < 68 hPa, differenced against the no-eruption control `exovolc_hunga_control`) across the water sweep with the Zhou et al. (2026) decay fit, plus a height–time anomaly Hovmöller. |
| `hunga_plume_linearity` | Model-behavior diagnostic: plume-core AOD does **not** scale linearly with injected SO₂, though the global mean does exactly. |
| `trappist_removal_regimes` | The four TRAPPIST-1 e publication figures: removal regimes (AOD lifetime + magnitude scaling + year-6 retention), stratospheric water plume, eruption-parameter sensitivities, and height–time aerosol structure. |
| `trappist_sulfur_budget` | Mechanism figures for the water-limited conversion: sulfur partitioned into SO₂ / H₂SO₄ gas / aerosol in Tg S, and the ben2 water burden against each eruption's stoichiometric water demand. |
| `trappist_diagnostics` | Exploratory quad run first against the fresh ben2/hab2 output. Superseded by the two folders above; keep as a scratch reference, do not source manuscript figures from it. |

### Scope constraints encoded in the Hunga configs

Several Hunga YAML headers carry retracted or easily-misapplied reference
values, and the constraints matter more than they look:

- **No global-mean observational target** belongs on `hunga_validation`. The
  often-quoted Khaykin et al. (2022) 0.047 ± 0.011 dispersed SAOD does not
  appear as stated text in that paper; `references:` is empty on purpose.
- **No observational line at all** belongs on `hunga_plume_linearity` —
  in-plume AOD is not comparable to zonal- or global-mean SAOD.
- The Schoeberl figures are **not** replications. The only quantity taken from
  Schoeberl et al. (2023) is the verified 0.018 global-mean peak.
- `hunga_zonal_level`'s aerosol panel plots sulfate mass density, not
  extinction: comparable in *shape* to Schoeberl Fig. 3a, not in magnitude. Its
  water panel needs the kg/kg → ppmv factor 1.6076 (CAM's Q is a mass mixing
  ratio; MLS reports volume) — omitting it understates the model by ~38%.
- The scheme has **no H₂O-dependent sulfate pathway**, so AOD is identical
  across the entire water sweep. `hunga_water_persistence` validates the water
  anomaly alone and claims no aerosol effect of the water.
- `hunga_water_persistence` differences against `exovolc_hunga_control`, a true
  no-eruption run, **not** against the dry-injection case
  `exovolc_hunga_h2o_none` — the latter still injects the full SO₂ mass, so it
  is not a clean background. The control is an independent realization, so its
  internal variability no longer cancels in the difference: the excess carries
  a noise floor of roughly ±10 Tg (~2% of the ~620 Tg background reservoir,
  which itself swings ~24 Tg seasonally). Consequences: the 0 Tg H₂O sweep
  member is no longer identically zero but wanders to ~9 Tg — that is the
  noise floor, **not** an SO₂-driven moistening signal, and should not be
  quoted as one. The late-time tail wobbles for the same reason. The peak
  (154 Tg) and the decay fit are well above this floor and unaffected.

Note that `hunga_water_persistence` deliberately keeps `exovolc_hunga_fid` as
its reference case rather than the manuscript fiducial: there it is the 146 Tg
member of a sweep that varies water at *fixed* SO₂ and Reff, and substituting a
case differing in three parameters would destroy that control. The manuscript
fiducial is drawn as a separate overlay instead.

### Scope constraints encoded in the TRAPPIST-1 configs

- **ben2 and hab2 only.** The ben1/hab1 suites were compiled with pure CO₂ where
  N₂ + 400 ppm CO₂ was specified, so those 28 cases are invalid and are being
  rerun. Everything here is a dry-vs-wet **surface** contrast at **fixed**
  composition. Do not add ben1/hab1 cases to these configs until they rerun.
- **No observational reference line belongs on any TRAPPIST-1 panel.** There is
  no measured TRAPPIST-1 e volcanic aerosol or stratospheric water abundance.
  Importing a band from the Earth validation cases would carry a terrestrial
  calibration onto a different planet.
- **The 100× rung is not a clean ladder member.** Those cases inject over 48 h
  rather than 24 h; the injected mass is as labeled but the rate is halved.
  `trappist_removal_regimes` rings that point in panel (b); keep the marking.
- **Planet constants come from the compiled `exoplanet_mod.F90`**, not from the
  stale `exovolc_ben.yaml` / `exovolc_hab.yaml` in exovolcano-analysis, which
  have `r_air` crossed between the suites.
- **ben2 output after ~day 250 is contaminated** by a suspected spurious water
  source (the dry-atmosphere controls gain ~122 Tg of water from an initial
  `Q ≡ 0`). Because the conversion consumes water, that drift feeds late aerosol
  production. `trappist_sulfur_budget` shades the affected window rather than
  cropping it. Treat any ben2 quantity past day ~250 as unusable until the
  underlying bug is resolved.
- **Use `TMQ` × planet area for water burden**, not a `profiles/Q.csv` column
  integral. TMQ recovers the injection accurately (144.6 Tg on day 1 against
  146 Tg injected for Hunga).
- **Convert sulfur species to Tg of S** via the S mass fraction (SO₂ 0.5005,
  H₂SO₄ and sulfate aerosol 0.3269) so the three curves are comparable and the
  budget closes. H₂SO₄ gas always sits at ~1e-16 Tg S: condensation
  (`TAU_AER_CONV` = 1800 s) far outpaces its production, so it never
  accumulates. That flat near-zero line is correct, not a missing field.

Each `plot_*.py` inserts the project root into `sys.path` so `pub_data` is importable regardless of invocation directory.

CSV column convention: column 0 = time (days), column 1 = the variable of interest. Read with `pd.read_csv(..., header=0)` and index via `.iloc`.

Profile CSVs are the exception: the first two lines are `# pressure_Pa: ...` and
`# altitude_m: ...` comments carrying the 51-level vertical coordinate, followed
by a `days_since_start,<VAR>_lev0..lev50` header. Read with `skiprows=2` and
parse the two comment lines separately for the vertical grid.

Zonal AOD CSVs carry one column per latitude, and the header names *are* the
latitude values (46 lats, -90..90 step 4).

## Figure output formats

The manuscript is AASTeX and includes figures as `.eps`, so each script saves
**both** `.pdf` and `.eps`. Note that matplotlib's PostScript backend does not
support real transparency: use solid pale fills rather than alpha-blended
`axhspan`/`fill_between`, or the EPS will warn and rasterize.

Rendered `.eps` files are gitignored here and committed to the manuscript repo
(`ExoVolcano_ExoCAM_Part1`) instead, so the two can drift. After regenerating a
figure, re-copy the `.eps` across.

## Adding a new figure

1. Create a new folder `<name>/` in the project root.
2. Copy `config_*.yaml` and `plot_*.py` from one of the nine current figures
   (not from the stale ones), update the config filename in the `yaml.safe_load`
   call, and adjust the plotting logic. For a Hunga figure, copy from a Hunga
   folder so you inherit the two-root `base_dirs` handling.
3. Build data paths directly as `base_dir/<eruption>/<case>/data/<subpath>`.
   Do not use `pub_data.find_csvs_list` unless you first repoint it.
4. The `sys.path` block and `here`-relative file paths can be copied verbatim.
5. Record in the config header what the figure does and does not claim, and
   why any reference line is legitimate. Several Hunga headers exist to stop a
   retracted value being reintroduced; keep that convention.

## Dependencies

`pandas`, `matplotlib`, `pyyaml` (standard scientific Python stack — no special install steps documented).
