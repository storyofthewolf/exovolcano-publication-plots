# DEVELOPER_NOTES.md

Implementation reference for exovolcano-publication-plots.

## Data layout (authoritative)

```
remote_analysis/<eruption>/<case_name>/data/<subpath>
```

`<eruption>` is `tambora`, `pinatubo`, or `hunga`. The case name is **not**
repeated inside `data/`.

Per-case files:

| Subpath | Contents |
|---|---|
| `aod/aod_550nm_band.csv` | `days_since_start, AOD_550nm_band` — global mean |
| `aod/aod_zonal_550nm_band.csv` | `days_since_start`, then one column per latitude; header names are the latitude values (46 lats, -90..90 step 4) |
| `scalar/SO2.csv` | `days_since_start, SO2 [kg]` — global burden |
| `scalar/H2SO4.csv` | H2SO4 gas burden [kg] |
| `scalar/VOLCHZMD.csv` | condensed sulfate aerosol burden [kg] |
| `scalar/Q.csv` | total water mass [kg] |
| `scalar/{TS,TMQ,TGCLDLWP}.csv` | area-mean surface temperature, precipitable water, cloud LWP |
| `profiles/{T,Q,SO2,H2SO4,VOLCHZMD}.csv` | area-mean vertical profiles, 51 levels |
| `zonal/<VAR>_day<NNNN.NN>.csv` | zonal-mean snapshots at days 0, 1, 4, 10, 50, 100 |

All runs are 2190 days (6 years) at daily cadence.

### Profile CSV format

First two lines are comments carrying the vertical coordinate, then the header:

```
# pressure_Pa: <51 comma-separated values>
# altitude_m: <51 comma-separated values>
days_since_start,<VAR>_lev0,<VAR>_lev1,...,<VAR>_lev50
```

Read with `pd.read_csv(path, skiprows=2)` and parse the two comment lines
manually for the pressure/altitude grid. Level 0 is the model top
(~1.08 Pa); level 50 is the surface (~97577 Pa).

## pub_data.py — LEGACY, does not resolve current paths

**`find_csvs(base_dir, case_glob, subpath)`** and
**`find_csvs_list(base_dir, cases, subpath)`** both resolve paths as:

```
base_dir / case_name / data / case_name / subpath
```

This is the **superseded** convention: no eruption directory, case name
repeated. Against the current layout every case reports as missing (printed to
stdout, not raised). `aod_timeseries/` and `aod_zonal_contour/` still import
`find_csvs_list` and are therefore broken; `aod_twopanel/` hardcodes the same
old convention plus the old CSV basenames (`aod_0p550um_mie.csv`).

Either repoint these helpers or build paths directly. The four 2026-07 figure
scripts build paths directly.

## Current figure scripts

Each folder holds `plot_<name>.py`, `config_<name>.yaml`, and gitignored
`.pdf`/`.eps` output. Every script accepts `--config <file>` (relative to the
script's own directory) and saves both formats at 300 dpi.

### tambora_validation/

Two panels: (a) global-mean 550 nm AOD, fiducial heavy over the 15-case
envelope, with the VolMIP peak-SAOD band and ensemble mean overlaid;
(b) latitude-time Hovmoller of zonal-mean AOD.

Config keys: `base_dir`, `eruption`, `case`, `cases` (list),
`subpath_timeseries`, `subpath_zonal`, `outfile_base`, `eruption_date`, `xlim`,
`volmip_saod_min` / `_max` / `_mean`.

### tambora_sulfur_budget/

Two stacked panels sharing the x-axis: (a) SO2, H2SO4 gas, and condensed
aerosol burdens plus their total, all converted to Tg of **sulfur** so the
reservoirs are comparable; (b) global-mean AOD.

Config keys: `base_dir`, `eruption`, `case`, `subpath_so2`, `subpath_h2so4`,
`subpath_volc`, `subpath_aod`, `outfile_base`, `eruption_date`, `xlim`,
`mw_s`, `mw_so2`, `mw_h2so4`, `mw_volc`, `injected_total_s_tg`.

Sulfur conversion: SO2 kg x 32.06/64.06; H2SO4 gas and VOLCHZMD kg x
32.06/98.08. VOLCHZMD is condensed sulfate carried as H2SO4 mass. The H2SO4
gas reservoir is a short-lived intermediate and is ~0 at plot scale; its legend
entry is kept deliberately so a reader sees it was checked, not omitted.

### pinatubo_validation/

Same two-panel structure as `tambora_validation`, with horizontal reference
lines instead of a band.

Config keys: `base_dir` (**includes the eruption directory here**), `case`,
`case_glob`, `subpath_timeseries`, `subpath_zonal`, `outfile_stem`,
`eruption_date`, `xlim`, and a `references` list of
`{label, value, kind, color}` where `kind` is `model` or `obs`.

### pinatubo_sensitivity/

2x2 one-at-a-time sweeps, shared y-axis, fiducial drawn heavier in every panel.

Config keys: `base_dir` (includes eruption), `subpath_timeseries`,
`outfile_stem`, `eruption_date`, `xlim`, `fiducial`, and `panels` — a list of
`{key, title, cases: [{case, label}, ...]}`.

### Config-key inconsistency

The Tambora scripts use `base_dir` plus a separate `eruption` key and name the
output `outfile_base`; the Pinatubo scripts fold the eruption into `base_dir`
and use `outfile_stem`. Both work. Do not assume uniformity when copying a
config between them.

## Style conventions

- Serif fonts; ~7.1 in wide for two-column-spanning figures, ~3.5 in single
  column; 8-9 pt labels and ticks.
- Panel labels `(a)`, `(b)` as in-axes text annotations, not `set_title`.
- No figure titles: captions live in the LaTeX.
- Colorblind-safe Okabe-Ito-derived palette for lines; `magma_r` for the
  Hovmoller contours.
- `bbox_inches='tight'` on save.
- **EPS has no transparency.** Use solid pale fills rather than alpha-blended
  `axhspan`/`fill_between`, or the PostScript backend warns and rasterizes.

## Verification habit

Each script prints diagnostics to stdout on every run (peak AOD and its day,
zonal peak and its latitude, peak burdens, per-sweep peaks). These were checked
against independently computed values before the figures were accepted:

| Quantity | Tambora fid | Pinatubo fid |
|---|---|---|
| peak global-mean AOD | 0.5708 @ day 246 | 0.1695 @ day 186 |
| peak zonal-mean AOD | 1.2186 @ lat -42 | 0.5397 @ lat +90 |
| peak SO2 burden | 59.15 Tg | 17.74 Tg |
| peak aerosol (VOLCHZMD) | 91.21 Tg | 27.08 Tg |

Keep this practice: print the numbers, check them against an independent
calculation, and do not silently adjust a script until a discrepancy is
understood.

## Dependencies

`pandas`, `matplotlib`, `pyyaml`, `numpy`.
