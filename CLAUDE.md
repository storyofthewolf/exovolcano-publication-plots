# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Purpose

Bespoke publication-quality plotting scripts for:
> Wolf, E.T. et al. (2026) "Simulating volcanic eruptions in diverse planetary context with ExoCAM: Descriptions and baselines" *(In preparation)*

These are **not** general-purpose routines — each figure script is tailored for a specific plot.

## Data pipeline

CSV input files are produced by [exovolcano-analysis](https://github.com/storyofthewolf/exovolcano-analysis) from ExoCAM `.h1.` time-series output. The expected on-disk layout is:

```
remote_analysis/
  <case_name>/
    data/
      <case_name>/
        aod/aod_550nm_band.csv
        ...
```

`base_dir` in each figure's YAML points to the local copy of `remote_analysis` (currently `/Users/wolfe/Desktop/projects/volcanos/remote_analysis`).

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

```
pub_data.py
aod_timeseries/
  plot_aod_timeseries.py
  config_aod_timeseries.yaml
  fig_aod_timeseries.pdf        ← gitignored output
aod_zonal_contour/
  plot_aod_zonal_contour.py
  config_aod_zonal_contour.yaml
  fig_aod_zonal_contour.pdf
```

Each `plot_*.py` inserts the project root into `sys.path` so `pub_data` is importable regardless of invocation directory.

`find_csvs_list` resolves paths as `base_dir/case_name/data/case_name/subpath`. Missing CSVs are printed to stdout and omitted (not raised).

CSV column convention: column 0 = time (days), column 1 = the variable of interest. Read with `pd.read_csv(..., header=0)` and index via `.iloc`.

## Adding a new figure

1. Create a new folder `<name>/` in the project root.
2. Copy `config_*.yaml` from an existing figure; update `base_dir`, `subpath`, `outfile`, and the `cases` list.
3. Copy `plot_*.py` from a similar figure; update the config filename in the `yaml.safe_load` call and adjust plotting logic.
4. The `sys.path` block and `here`-relative file paths can be copied verbatim.

## Dependencies

`pandas`, `matplotlib`, `pyyaml` (standard scientific Python stack — no special install steps documented).
