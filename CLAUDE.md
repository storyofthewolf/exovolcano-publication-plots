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

`BASE_DIR` in each figure script points to the local copy of `remote_analysis` (currently hardcoded to `/Users/wolfe/Desktop/projects/volcanos/remote_analysis`).

## Running a figure script

```bash
python fig_aod_timeseries.py
```

Each script is self-contained: edit the constants block at the top (BASE_DIR, CASE_GLOB, SUBPATH, EXCLUDE, etc.), then run it directly. Output is a PDF saved to the working directory.

## Code structure

| File | Role |
|------|------|
| `pub_data.py` | Shared utility — `find_csvs(base_dir, case_glob, subpath)` returns `{case_name: csv_path}` |
| `fig_aod_timeseries.py` | Figure: AOD at 550 nm vs time for all `exovolc_tambora_*` cases |

New figure scripts should import `find_csvs` from `pub_data` and follow the same constants-block pattern.

## Dependencies

`pandas`, `matplotlib` (standard scientific Python stack — no special install steps documented).
