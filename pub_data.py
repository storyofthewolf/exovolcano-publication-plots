"""
pub_data.py - CSV discovery for publication plotting.

Single function: find_csvs() scans a remote_analysis directory tree
and returns a mapping of case_name -> csv_path for any matching files.

Usage example:
    from pub_data import find_csvs

    cases = find_csvs(
        base_dir  = 'remote_analysis',
        case_glob = 'exovolc_tambora_*',
        subpath   = 'aod/aod_550nm_band.csv',
    )
    # cases = {'exovolc_tambora_k50d': '/path/to/.../aod_550nm_band.csv', ...}
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
