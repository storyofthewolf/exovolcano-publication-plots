# DEVELOPER_NOTES.md

Implementation reference for exovolcano-publication-plots.

## pub_data.py

**Function: `find_csvs(base_dir, case_glob, subpath)`**

Discovers CSV files across experiment case directories.

**Parameters:**
- `base_dir` (str): Root directory (e.g., `'remote_analysis'`)
- `case_glob` (str): Glob pattern for case directories (e.g., `'exovolc_tambora_*'`)
- `subpath` (str): Path from case data dir to target CSV (e.g., `'aod/aod_550nm_band.csv'`)

**Returns:**
- `dict`: `{case_name: absolute_csv_path}`, sorted by case_name. Only includes cases where the CSV exists.

**Directory convention:**
```
base_dir / case_name / data / case_name / subpath
```

Missing CSVs are printed to stdout (not raised as errors).

## fig_aod_timeseries.py

Reads AOD 550nm time series from all matching Tambora cases and plots on one axes.

**Configuration constants (edit to customize):**
- `BASE_DIR`: Root of remote_analysis directory
- `CASE_GLOB`: Glob pattern for case dirs to include
- `SUBPATH`: Relative path to CSV within each case
- `DAYS_PER_YEAR`: Conversion factor for x-axis (default 365.0)
- `OUTFILE`: Output PDF filename
- `EXCLUDE`: Set of case short names to skip

**Data format (CSV):**
- Column 0: Time (days)
- Column 1: AOD at 550 nm

**Output:**
- High-resolution PDF (300 dpi) with publication-quality styling
- Legend: two-column layout, frameless, 7pt font
- Peak AOD values printed to stdout for each case

## Adding a new figure script

1. Create `fig_<name>.py` in the root directory
2. Import `find_csvs` from `pub_data`
3. Add constants block at the top matching the pattern in `fig_aod_timeseries.py`
4. Call `find_csvs(BASE_DIR, CASE_GLOB, SUBPATH)` to get `{case_name: csv_path}` dict
5. Iterate over cases, read CSV with `pd.read_csv()`, plot as needed
6. Save output to PDF with `plt.savefig(OUTFILE, dpi=300, bbox_inches='tight')`
