"""
Varve-counted age-depth model plot.

Reads varve counts for one or more cores from an Excel workbook and draws
all cores on one age-depth figure: age (yr CE) on a top x-axis, depth (cm)
increasing downward, one stepped line per core.

Accepted workbook layouts (detected automatically)
--------------------------------------------------
1. One sheet per core (sheet name = core ID), with a depth column and a
   year column, e.g.

       Depth (cm) | Year (CE)
       0.00       | 2023
       0.42       | 2022
       ...

2. One sheet holding several cores side by side, as depth/year column
   pairs whose headers start with the core ID, e.g.

       GUAC-22A-1G-1 Depth | GUAC-22A-1G-1 Year | GUAC-23A-1G-1 Depth | ...

3. Varve thickness only (no year column): a column of varve thicknesses
   in mm, one row per varve from the top down. Depth is the cumulative
   thickness and the year is counted back from --top-year.

Column headers are matched case-insensitively: "depth"; "year", "age",
"ce" or "ad"; "thick" (thickness in mm). Rows with blanks are skipped.

Usage
-----
    python plot_age_depth.py varve_counts.xlsx
    python plot_age_depth.py varve_counts.xlsx -o age_depth.png \
        --xlim 1700 2050 --ylim 90 0 --style steps
"""

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Colours matched to the Excel chart (Office theme, in series order).
CORE_COLORS = [
    "#156082",  # dark teal-blue
    "#E97132",  # orange
    "#196B24",  # dark green
    "#0F9ED5",  # light blue
    "#A02B93",  # purple
    "#4EA72E",  # light green
    "#0E2841",  # navy
    "#8B3A0F",  # brown
    "#5E5E5E",  # grey (extra cores)
    "#C9A100",  # mustard
]

DEPTH_RE = re.compile(r"depth", re.I)
YEAR_RE = re.compile(r"\b(year|age|ce|ad)\b|year|age", re.I)
THICK_RE = re.compile(r"thick", re.I)


def _find_cols(columns, pattern):
    return [c for c in columns if pattern.search(str(c))]


def _core_name(header, fallback):
    """Strip the depth/year words and units from a header to get the core ID."""
    name = re.sub(r"\(.*?\)|depth|year|age|thickness|\bce\b|\bad\b|\bcm\b|\bmm\b",
                  "", str(header), flags=re.I)
    name = name.strip(" _-:")
    return name or fallback


def _clean(depth, year):
    df = pd.DataFrame({"depth": pd.to_numeric(depth, errors="coerce"),
                       "year": pd.to_numeric(year, errors="coerce")}).dropna()
    return df.sort_values(["depth", "year"], ascending=[True, False]).reset_index(drop=True)


def read_cores(xlsx_path, top_year=None):
    """Return {core_id: DataFrame(depth, year)} from every sheet of the workbook."""
    sheets = pd.read_excel(xlsx_path, sheet_name=None)
    cores = {}
    for sheet_name, df in sheets.items():
        df = df.dropna(how="all").dropna(axis=1, how="all")
        if df.empty:
            continue
        depth_cols = _find_cols(df.columns, DEPTH_RE)
        year_cols = [c for c in _find_cols(df.columns, YEAR_RE) if c not in depth_cols]
        thick_cols = _find_cols(df.columns, THICK_RE)

        if depth_cols and year_cols:
            if len(depth_cols) != len(year_cols):
                raise ValueError(f"Sheet '{sheet_name}': {len(depth_cols)} depth columns "
                                 f"but {len(year_cols)} year columns.")
            for d, y in zip(depth_cols, year_cols):
                name = sheet_name if len(depth_cols) == 1 else _core_name(d, sheet_name)
                cores[name] = _clean(df[d], df[y])
        elif thick_cols:
            if top_year is None:
                raise ValueError(f"Sheet '{sheet_name}' has only varve thickness; "
                                 "pass --top-year (year of the topmost varve).")
            for t in thick_cols:
                name = sheet_name if len(thick_cols) == 1 else _core_name(t, sheet_name)
                thick_cm = pd.to_numeric(df[t], errors="coerce").dropna().to_numpy() / 10.0
                depth = np.concatenate([[0.0], np.cumsum(thick_cm)])
                year = top_year - np.arange(len(depth))
                cores[name] = _clean(depth, year)
        else:
            print(f"Skipping sheet '{sheet_name}': no depth/year or thickness columns found.")
    if not cores:
        raise ValueError("No core data found in the workbook.")
    return cores


def _nice_limits(lo, hi, step):
    return np.floor(lo / step) * step, np.ceil(hi / step) * step


def plot_age_depth(cores, out_path, xlim=None, ylim=None, style="steps",
                   title=None, xstep=50, ystep=10):
    fig, ax = plt.subplots(figsize=(11, 8))

    for i, (name, df) in enumerate(cores.items()):
        ax.plot(df["year"], df["depth"],
                color=CORE_COLORS[i % len(CORE_COLORS)], lw=1.6, label=name,
                drawstyle="steps-post" if style == "steps" else "default")

    all_years = pd.concat([d["year"] for d in cores.values()])
    all_depths = pd.concat([d["depth"] for d in cores.values()])
    if xlim is None:
        xlim = _nice_limits(all_years.min(), all_years.max(), xstep)
    if ylim is None:
        ylim = (_nice_limits(0, all_depths.max(), ystep)[1], 0)
    ax.set_xlim(xlim)
    ax.set_ylim(max(ylim), min(ylim))  # depth increases downward

    # Age axis on top, like the Excel chart.
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position("top")
    ax.set_xticks(np.arange(min(xlim), max(xlim) + xstep, xstep))
    ax.set_yticks(np.arange(min(ylim), max(ylim) + ystep, ystep))
    ax.set_xlabel("Age (yr CE)", fontsize=12, labelpad=8)
    ax.set_ylabel("Depth (cm)", fontsize=12)
    ax.tick_params(labelsize=10, direction="out")
    for side in ("right", "bottom"):
        ax.spines[side].set_visible(False)

    ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.6), frameon=False,
              fontsize=11, handlelength=1.5)
    if title:
        ax.set_title(title, fontsize=13, pad=40)

    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    print(f"Saved {out_path}")
    return fig, ax


def main():
    p = argparse.ArgumentParser(description="Plot a varve-counted age-depth model from Excel.")
    p.add_argument("excel", help="Varve counting workbook (.xlsx)")
    p.add_argument("-o", "--out", default="age_depth_model.png",
                   help="Output figure (.png, .pdf or .svg). Default: age_depth_model.png")
    p.add_argument("--xlim", nargs=2, type=float, metavar=("MIN_YEAR", "MAX_YEAR"),
                   help="Age axis limits, e.g. 1700 2050")
    p.add_argument("--ylim", nargs=2, type=float, metavar=("MAX_DEPTH", "MIN_DEPTH"),
                   help="Depth axis limits in cm, e.g. 90 0")
    p.add_argument("--style", choices=["steps", "line"], default="steps",
                   help="'steps' draws the staircase look of the Excel chart (default)")
    p.add_argument("--top-year", type=int,
                   help="Year of the topmost varve (only needed for thickness-only sheets)")
    p.add_argument("--title", help="Optional figure title")
    args = p.parse_args()

    cores = read_cores(args.excel, top_year=args.top_year)
    for name, df in cores.items():
        print(f"  {name}: {len(df)} points, {df['depth'].max():.1f} cm, "
              f"{df['year'].min():.0f}-{df['year'].max():.0f} CE")
    plot_age_depth(cores, Path(args.out), xlim=args.xlim, ylim=args.ylim,
                   style=args.style, title=args.title)


if __name__ == "__main__":
    main()
