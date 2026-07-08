"""
pub_data.py - CSV discovery for publication plotting.

find_csvs()      — glob-based discovery (legacy/exploratory use)
find_csvs_list() — explicit case list (preferred for figure scripts)

Usage example:
    from pub_data import find_csvs_list

    cases = find_csvs_list(
        base_dir  = '/path/to/remote_analysis',
        cases     = ['exovolc_tambora_k10d', 'exovolc_tambora_k20d'],
        subpath   = 'aod/aod_550nm_band.csv',
    )
    # cases = {'exovolc_tambora_k10d': '/path/to/.../aod_550nm_band.csv', ...}
"""

import glob
import os


def find_csvs(base_dir, case_glob, subpath):
    """
    Discover CSV files across multiple experiment cases.

    Parameters
    ----------
    base_dir  : str   root directory, e.g. 'remote_analysis'
    case_glob : str   glob pattern for case directories, e.g. 'exovolc_tambora_*'
    subpath   : str   path from the case data dir to the target CSV,
                      e.g. 'aod/aod_550nm_band.csv'

    Returns
    -------
    dict  {case_name: absolute_csv_path}  sorted by case_name.
          Only cases where the CSV actually exists are included.

    Directory convention assumed:
        base_dir / case_name / data / case_name / subpath
    """
    pattern = os.path.join(base_dir, case_glob)
    case_dirs = sorted(glob.glob(pattern))

    result = {}
    for case_dir in case_dirs:
        case_name = os.path.basename(case_dir)
        csv_path  = os.path.join(case_dir, 'data', case_name, subpath)
        if os.path.exists(csv_path):
            result[case_name] = csv_path
        else:
            print(f"  missing: {csv_path}")

    return result


def find_csvs_list(base_dir, cases, subpath):
    """
    Resolve CSV paths for an explicit list of case names.

    Parameters
    ----------
    base_dir : str   root directory, e.g. '/path/to/remote_analysis'
    cases    : list  ordered list of case names, e.g. ['exovolc_tambora_k10d', ...]
    subpath  : str   path from the case data dir to the target CSV

    Returns
    -------
    dict  {case_name: absolute_csv_path}  in the order given by `cases`.
          Missing CSVs are printed to stdout and omitted from the result.
    """
    result = {}
    for case_name in cases:
        csv_path = os.path.join(base_dir, case_name, 'data', case_name, subpath)
        if os.path.exists(csv_path):
            result[case_name] = csv_path
        else:
            print(f"  missing: {csv_path}")
    return result
