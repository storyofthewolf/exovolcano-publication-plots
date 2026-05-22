# exovolcano-publication-plots

Publication-quality plotting scripts for:
> Wolf, E.T. et al. (2026) "Simulating volcanic eruptions in diverse planetary context with ExoCAM: Descriptions and baselines" *(In preparation)*

## Overview

This repo contains **bespoke plotting scripts** tailored for specific publication figures. Each script is self-contained with a constants block at the top for configuration.

### Data source

CSV inputs are produced by [exovolcano-analysis](https://github.com/storyofthewolf/exovolcano-analysis) from ExoCAM `.h1.` time-series output files. See that repo for how to generate the CSV files used here.

### Quick start

```bash
# Edit BASE_DIR, CASE_GLOB, etc. at the top of the script
python fig_aod_timeseries.py

# Output: fig_aod_timeseries.pdf
```

## Scripts

- **fig_aod_timeseries.py**: AOD at 550 nm vs time for all `exovolc_tambora_*` cases

## Dependencies

- `pandas`
- `matplotlib`

See [CLAUDE.md](CLAUDE.md) for architecture and [DEVELOPER_NOTES.md](DEVELOPER_NOTES.md) for implementation details.
