# exovolcano-publication-plots

Publication-quality plotting scripts for:
> Wolf, E.T. et al. (2026) "Simulating volcanic eruptions in diverse planetary context with ExoCAM: Descriptions and baselines" *(In preparation)*

## Overview

This repo contains **bespoke plotting scripts** tailored for specific publication figures. Each figure lives in its own folder as a `plot_<name>.py` / `config_<name>.yaml` pair: data selection and paths go in the YAML, plot structure in the Python.

### Data source

CSV inputs are produced by [exovolcano-analysis](https://github.com/storyofthewolf/exovolcano-analysis) from ExoCAM `.h1.` time-series output files. See that repo for how to generate them. The expected layout is:

```
remote_analysis/<eruption>/<case_name>/data/<subpath>
```

### Quick start

```bash
cd tambora_validation
python plot_tambora_validation.py

# Outputs fig_tambora_validation.pdf and .eps into the same folder
```

Point `base_dir` in the YAML at your local copy of `remote_analysis`. Each script also accepts `--config <file>` to run a variant.

## Figures

Current (2026-07), built against the current data layout:

- **tambora_validation/** — global-mean 550 nm AOD against the VolMIP-Tambora peak-SAOD band, plus a latitude-time Hovmoller
- **tambora_sulfur_budget/** — SO2 / H2SO4 / sulfate-aerosol partitioning in Tg S, with AOD on a shared time axis
- **pinatubo_validation/** — global-mean AOD against WACCM and GloSSAC reference values, plus a Hovmoller
- **pinatubo_sensitivity/** — 2x2 one-at-a-time parameter sweeps about the fiducial

TRAPPIST-1 e (2026-08), built against the `ben2_suite1` / `hab2_suite1` output:

- **trappist_removal_regimes/** — the four publication figures: aerosol removal regimes, stratospheric water plume, eruption-parameter sensitivities, height-time aerosol structure
- **trappist_sulfur_budget/** — sulfur partitioning in Tg S, and the ben2 water burden against each eruption's stoichiometric water demand
- **trappist_diagnostics/** — exploratory quad, superseded by the two above; not a source of manuscript figures

Older folders (`aod_timeseries/`, `aod_zonal_contour/`, `aod_twopanel/`, `sulfur_burden/`) predate a change in the `remote_analysis` directory layout and no longer resolve their input paths. Keep them as style references; see DEVELOPER_NOTES.md before reusing them.

## Output formats

Scripts save both `.pdf` and `.eps` at 300 dpi. The manuscript is AASTeX and includes the `.eps` versions, which are committed to the manuscript repo rather than here (outputs are gitignored in this repo). After regenerating a figure, re-copy the `.eps` across.

## Dependencies

- `pandas`
- `matplotlib`
- `pyyaml`
- `numpy`

See [CLAUDE.md](CLAUDE.md) for architecture and [DEVELOPER_NOTES.md](DEVELOPER_NOTES.md) for implementation details.
