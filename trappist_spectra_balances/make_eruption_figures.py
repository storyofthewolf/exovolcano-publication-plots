#!/usr/bin/env python3
"""Generate the six per-eruption balances figures: one eruption, four atmospheres.

This transposes the manuscript figure's organisation. There, the four columns
are four DIFFERENT contrasts chosen to make one argument, and the case name
varies down the row. Here each figure holds ONE eruption fixed and puts the
four TRAPPIST-1 e atmospheres side by side:

    columns:   ben1   ben2   hab1   hab2
    figures:   hunga, hunga_10x, hunga_100x,
               tambora, tambora_10x, tambora_100x

so that reading across a row isolates the atmosphere's effect at fixed
forcing, and reading across the six figures isolates the magnitude ladder at
fixed atmosphere. The three rows are unchanged from the template: raw
transmission spectra, difference-from-day-0 spectra, and the GCM burdens that
drive them.

Everything except `columns` and `outfile_stem` is inherited from
config_trappist_spectra_balances.yaml, exactly as make_variants.py does it, so
the R = 250 binning, the epoch set, and the shared axis limits cannot drift
away from the figure these are compared against.

THE TWO LADDERS SCALE DIFFERENT QUANTITIES. Tambora's rungs scale SO2 at
roughly fixed H2O; Hunga's scale H2O at roughly fixed SO2 (~1 Tg). A
"hunga_100x" therefore carries the SAME sulfur as "hunga" and differs only in
water. The rung is labelled on each panel, but the two ladders are not
comparable rung-for-rung and should not be read as if they were.

    python make_eruption_figures.py            # write configs and render
    python make_eruption_figures.py --list     # show what would be made
    python make_eruption_figures.py --only tambora_10x

DATA AVAILABILITY: the 10x rungs have no computed spectra. See the module
constant SPECTRA_MISSING below and the note this script prints.
"""

from __future__ import annotations

import argparse
import copy
import os
import subprocess
import sys

import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.join(HERE, "config_trappist_spectra_balances.yaml")

#: The four TRAPPIST-1 e atmospheres, in column order. Dry pair first, moist
#: pair second, and within each pair N2 before CO2 -- so moving right crosses
#: the surface-water boundary once, at the ben2 | hab1 seam, rather than
#: alternating. The composition contrast then reads vertically between the
#: adjacent columns of each pair.
ATMOSPHERES = ["ben1", "ben2", "hab1", "hab2"]

ATM_LABEL = {
    "ben1": "ben1: dry land, N$_2$",
    "ben2": "ben2: dry land, CO$_2$",
    "hab1": "hab1: aquaplanet, N$_2$",
    "hab2": "hab2: aquaplanet, CO$_2$",
}

#: Short human descriptor per atmosphere, for the panel title's first line.
ATM_SHORT = {
    "ben1": "dry land, N$_2$",
    "ben2": "dry land, CO$_2$",
    "hab1": "aquaplanet, N$_2$",
    "hab2": "aquaplanet, CO$_2$",
}

#: The six eruptions, in figure order: the water-rich ladder then the
#: sulfur-rich one, each ascending. (case, eruption label, rung label).
ERUPTIONS = [
    ("hunga",        "Hunga-like (water-rich)",  r"1$\times$"),
    ("hunga_10x",    "Hunga-like (water-rich)",  r"10$\times$"),
    ("hunga_100x",   "Hunga-like (water-rich)",  r"100$\times$"),
    ("tambora",      "Tambora-like (sulfur-rich)", r"1$\times$"),
    ("tambora_10x",  "Tambora-like (sulfur-rich)", r"10$\times$"),
    ("tambora_100x", "Tambora-like (sulfur-rich)", r"100$\times$"),
]

#: Rungs with NO computed transmission spectra on disk (2026-09-04). The
#: exovolcano-spectra sweep (scripts/run_all_suites.sh) covers only the 1x and
#: 100x rungs of each ladder -- 16 cases -- and the 10x raw h1 extracts are not
#: under data/raw/ either, so this is not a rerun-the-sweep away: the epochs
#: must first be pulled from Discover. The burden row draws from
#: remote_analysis and IS complete for all 24 cases.
SPECTRA_MISSING = {"hunga_10x", "tambora_10x"}


def make_columns(case, rung):
    """One column per atmosphere, at fixed eruption."""
    cols = []
    for atm in ATMOSPHERES:
        wet = atm.startswith("hab")
        cols.append({
            "key": f"exovolc_{atm}_{case}",
            "atm": atm,
            "case": case,
            # `balance` is retained because the template reads it for the
            # default title; the explicit `title` below supersedes it here.
            "balance": f"{'wet' if wet else 'dry'} ({rung})",
            "atm_label": ATM_LABEL[atm],
            # Title leads with the atmosphere, since that is what varies across
            # this figure's columns. The case name is constant here and lives
            # in the suptitle instead of being repeated four times.
            "title": f"{atm}\n{ATM_SHORT[atm]}",
            # No per-column note: the template's italic notes are claims about
            # a specific case's behaviour, and asserting one for 24 cases
            # would be inventing results. Left blank for the author to fill.
            "note": "",
        })
    return cols


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--list", action="store_true",
                    help="show the figures without writing or rendering")
    ap.add_argument("--only", default=None,
                    help="render just this eruption by case name")
    args = ap.parse_args()

    with open(BASE) as f:
        base = yaml.safe_load(f)

    rc = 0
    skipped = []
    for case, label, rung in ERUPTIONS:
        if args.only and case != args.only:
            continue
        print(f"\n=== {case}: {label}, {rung} rung")
        print(f"    columns: {', '.join(ATMOSPHERES)}")

        if case in SPECTRA_MISSING:
            # Refuse rather than draw two empty rows over a valid burden row.
            # The template already exits on a missing .npz; catching it here
            # keeps the message about the SET of figures rather than the first
            # column that happens to fail.
            print("    SKIPPED: no computed transmission spectra for this rung.")
            print("             Rows 1-2 have no data; the burden row alone is")
            print("             not this figure. See SPECTRA_MISSING in this file.")
            skipped.append(case)
            continue

        if args.list:
            continue

        cfg = copy.deepcopy(base)
        cfg["columns"] = make_columns(case, rung)
        cfg["outfile_stem"] = f"fig_eruption_{case}"
        # Names the eruption once, since the columns no longer carry it.
        # Literal en-dashes, not TeX "---": matplotlib does not apply TeX
        # ligatures, so a triple hyphen renders as three hyphens.
        cfg["suptitle"] = f"{label}  \u2013  {rung} rung  \u2013  TRAPPIST-1 e"

        # Burden legend moved UP for these figures. The base config's "center
        # left" was chosen against the manuscript figure's case mix, where the
        # middle of the leftmost panel is empty. Here the leftmost column is
        # always a DRY atmosphere, whose total water sits near 1e2 Tg -- right
        # where that legend lands, so it collided with the H2O curve on the
        # Hunga figures. The top-left corner is clear in all four: the dry
        # columns hold nothing above ~1e4 Tg.
        cfg["burden_legend_loc"] = "upper left"
        cfg["burden_legend_anchor"] = [0.01, 0.98]

        path = os.path.join(HERE, f"config_eruption_{case}.yaml")
        with open(path, "w") as f:
            f.write(f"# GENERATED by make_eruption_figures.py -- do not edit.\n"
                    f"# {label}, {rung} rung, across ben1/ben2/hab1/hab2.\n"
                    f"# Everything but `columns`, `outfile_stem` and `suptitle`\n"
                    f"# is inherited from config_trappist_spectra_balances.yaml;\n"
                    f"# change it there.\n")
            yaml.safe_dump(cfg, f, sort_keys=False, allow_unicode=True,
                           default_flow_style=False)

        r = subprocess.run(
            [sys.executable,
             os.path.join(HERE, "plot_trappist_spectra_balances.py"),
             "--config", os.path.basename(path)],
            cwd=HERE, capture_output=True, text=True,
        )
        for line in r.stdout.splitlines():
            if ("NOTE: epoch" in line or "peak |diff|" in line
                    or "shared raw" in line or "WARNING" in line):
                print("   ", line.strip())
        if r.returncode != 0:
            print(f"    FAILED:\n{r.stderr[-1500:]}")
            rc = 1
        else:
            print(f"    wrote fig_eruption_{case}.pdf / .eps")

    if skipped:
        print("\n" + "=" * 66)
        print("NOT RENDERED (no computed spectra): " + ", ".join(skipped))
        print("To add them, from exovolcano-spectra with the ssh master up:")
        print('  CASES="$(for a in ben1 ben2 hab1 hab2; do '
              'for e in hunga_10x tambora_10x; do '
              'echo exovolc_${a}_${e}; done; done)" scripts/pull_epochs.sh')
        print("  # then run scripts/run_timeseries.py per case (see")
        print("  # scripts/run_all_suites.sh for the exact invocation)")
        print("=" * 66)

    return rc


if __name__ == "__main__":
    raise SystemExit(main())
