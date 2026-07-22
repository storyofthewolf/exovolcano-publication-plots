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
  <eruption>/            # tambora | pinatubo | hunga
    <case_name>/         # e.g. exovolc_tambora_fid
      data/
        aod/aod_550nm_band.csv
        aod/aod_zonal_550nm_band.csv
        scalar/{SO2,H2SO4,VOLCHZMD,Q,TS,TMQ,TGCLDLWP}.csv
        profiles/{T,Q,SO2,H2SO4,VOLCHZMD}.csv
        zonal/<VAR>_day<NNNN.NN>.csv
```

`base_dir` in each figure's YAML points to the local copy of `remote_analysis` (currently `/Users/wolfe/Desktop/projects/volcanos/remote_analysis`).

### KNOWN TRAP: `pub_data` encodes the OLD layout

`pub_data.find_csvs()` and `find_csvs_list()` both resolve paths as
`base_dir/<case>/data/<case>/<subpath>` — the *superseded* convention, with no
eruption directory and the case name repeated. They no longer resolve against
the layout above. The CSV basenames also changed (`aod_550nm_band.csv`, not
`aod_0p550um_mie.csv`).

Consequently the three older figure folders (`aod_timeseries`,
`aod_zonal_contour`, `aod_twopanel`) are **stale**: the first two still import
`find_csvs_list`, and `aod_twopanel` hardcodes the old path and old CSV names.
They will report every case as missing until repointed. Treat them as style
references, not working scripts.

The four figure folders added in 2026-07 (`tambora_validation`,
`tambora_sulfur_budget`, `pinatubo_validation`, `pinatubo_sensitivity`) build
paths directly and are the pattern to copy.

## Running a figure script

From inside the figure folder:
```bash
cd aod_timeseries
python plot_aod_timeseries.py
```

Or from the project root:
```bash
python aod_timeseries/plot_aod_timeseries.py
```

Output PDF is written into the figure's own folder. Change `outfile` in the YAML to version outputs (e.g. `fig_aod_timeseries_v2.pdf`).

## Config split: YAML vs Python

**Edit the YAML** (`config_<name>.yaml`) to change: which cases are plotted, `base_dir`, `subpath`, `outfile`, `days_per_year`.

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
aod_timeseries/                 ← stale, old path convention
aod_zonal_contour/              ← stale
aod_twopanel/                   ← stale
sulfur_burden/
```

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
2. Copy `config_*.yaml` and `plot_*.py` from one of the four current figures
   (not from the stale ones), update the config filename in the `yaml.safe_load`
   call, and adjust the plotting logic.
3. Build data paths directly as `base_dir/<eruption>/<case>/data/<subpath>`.
   Do not use `pub_data.find_csvs_list` unless you first repoint it.
4. The `sys.path` block and `here`-relative file paths can be copied verbatim.

## Dependencies

`pandas`, `matplotlib`, `pyyaml` (standard scientific Python stack — no special install steps documented).
