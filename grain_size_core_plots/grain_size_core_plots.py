# -*- coding: utf-8 -*-
"""
Grain-size down-core plots at TRUE (real cm) scale — Google Colab script
=======================================================================

Plots pre-calculated grain-size statistics (D10, D50, D90, Mean, Sorting)
for two sediment cores. Every parameter is saved as its OWN figure, and the
depth axis is drawn at real scale (1 cm of core = 1 cm on the printed page
by default), so each plot can be placed directly next to the core photo
(e.g. in Illustrator / Inkscape / CorelDRAW / PowerPoint) without rescaling.

Outputs per core:
  * full-core plot for each parameter
  * one plot per event-depth interval for each parameter (exact top/bottom)
  * SVG (editable vectors), PDF and 600-dpi PNG for each
  * a zip file of everything, downloaded automatically in Colab

HOW TO USE IN COLAB
-------------------
1. Paste this whole file into a Colab cell (or upload it and `%run` it).
2. Edit the "USER SETTINGS" block below: file names, column names,
   core top/bottom depths and your event depths.
3. Run. If the data files are not found, Colab asks you to upload them.
   If you just want to test, set USE_DEMO_DATA = True.

INPUT FILE FORMAT (CSV or Excel, one file per core)
---------------------------------------------------
    Depth_cm, D10, D50, D90, Mean, Sorting
    0.5,      4.1, 22.3, 88.0, 5.6, 1.82
    1.5,      ...
Column names can be anything — map them in COLUMN_MAP below.
"""

import os
import shutil

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, AutoMinorLocator

# =============================================================================
# USER SETTINGS — edit this block only
# =============================================================================

USE_DEMO_DATA = False          # True = synthetic data, just to test the layout

CORES = {
    # Core name : settings
    "Core_1": {
        "file": "core1_grainsize.csv",   # .csv, .xlsx or .xls
        "sheet": 0,                      # Excel sheet (ignored for CSV)
        "top_cm": 0.0,                   # depth of core top (cm)
        "bottom_cm": 100.0,              # depth of core base (cm)
        # Event layers: (top_cm, bottom_cm, "label"). Put your exact depths here.
        "events": [
            (12.0, 18.5, "E1"),
            (41.0, 47.0, "E2"),
            (73.5, 80.0, "E3"),
        ],
    },
    "Core_2": {
        "file": "core2_grainsize.csv",
        "sheet": 0,
        "top_cm": 0.0,
        "bottom_cm": 120.0,
        "events": [
            (20.0, 26.0, "E1"),
            (58.0, 66.5, "E2"),
        ],
    },
}

# Map the standard names (left) to the column headers in YOUR file (right).
COLUMN_MAP = {
    "depth":   "Depth_cm",
    "D10":     "D10",
    "D50":     "D50",
    "D90":     "D90",
    "Mean":    "Mean",
    "Sorting": "Sorting",
}

# Axis label + units for each parameter (edit units to match your data:
# µm or φ). log_x = True is usual for D-values in µm.
PARAMS = {
    "D10":     {"label": "D10 (µm)",       "log_x": False},
    "D50":     {"label": "D50 (µm)",       "log_x": False},
    "D90":     {"label": "D90 (µm)",       "log_x": False},
    "Mean":    {"label": "Mean (µm)",      "log_x": False},
    "Sorting": {"label": "Sorting (σ, φ)", "log_x": False},
}

# Fixed x-limits per parameter so both cores share the same scale and can be
# compared side by side. None = automatic from the data of BOTH cores.
X_LIMITS = {
    "D10": None, "D50": None, "D90": None, "Mean": None, "Sorting": None,
}

# ---- True-scale settings ----------------------------------------------------
# PRINT_CM_PER_CORE_CM = 1.0 -> 1 cm of core = 1 cm on the figure (1:1).
# Use e.g. 0.5 if your core photo is printed at 1:2.
PRINT_CM_PER_CORE_CM = 1.0
PLOT_WIDTH_CM = 3.0            # width of the data panel itself (cm)

# Margins around the data panel (cm). They do NOT change the depth scale.
MARGIN_LEFT_CM = 1.6
MARGIN_RIGHT_CM = 0.3
MARGIN_TOP_CM = 1.3
MARGIN_BOTTOM_CM = 0.3

DEPTH_MAJOR_TICK_CM = 10       # full-core depth tick spacing
DEPTH_MINOR_TICK_CM = 1
EVENT_MAJOR_TICK_CM = 1        # event-plot depth tick spacing
EVENT_MINOR_TICK_CM = 0.5

MAKE_FULL_CORE_PLOTS = True
MAKE_EVENT_PLOTS = True
SHADE_EVENTS_ON_FULL_CORE = True
SHOW_EVENT_LABELS = True

# Style
LINE_COLOR = "#1f3b73"
EVENT_FILL = "#d9a441"
EVENT_ALPHA = 0.22
LINE_WIDTH = 1.0
MARKER_SIZE = 2.5
FONT_SIZE = 7
FONT_FAMILY = "Arial"          # falls back to DejaVu Sans if not installed

OUTPUT_DIR = "grain_size_plots"
FORMATS = ("svg", "pdf", "png")
PNG_DPI = 600

# =============================================================================
# END OF USER SETTINGS
# =============================================================================

CM = 1 / 2.54  # inches per cm

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": [FONT_FAMILY, "Liberation Sans", "DejaVu Sans"],
    "font.size": FONT_SIZE,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.minor.width": 0.4,
    "ytick.minor.width": 0.4,
    "svg.fonttype": "none",       # keep text editable in SVG
    "pdf.fonttype": 42,           # keep text editable in PDF
})


def in_colab():
    try:
        import google.colab  # noqa: F401
        return True
    except ImportError:
        return False


# -----------------------------------------------------------------------------
# Data loading
# -----------------------------------------------------------------------------
def make_demo_data(top, bottom, events, seed):
    rng = np.random.default_rng(seed)
    depth = np.arange(top + 0.5, bottom, 1.0)
    d50 = 25 + 5 * np.sin(depth / 9) + rng.normal(0, 2, depth.size)
    for et, eb, _ in events:  # coarser event beds
        inside = (depth >= et) & (depth <= eb)
        d50[inside] += 60 * np.exp(-((depth[inside] - et) / (eb - et + 1e-9)) * 2)
    d10 = d50 * 0.18 + rng.normal(0, 0.5, depth.size)
    d90 = d50 * 3.5 + rng.normal(0, 5, depth.size)
    mean = d50 * 1.1
    sorting = 1.4 + 0.012 * d50 + rng.normal(0, 0.08, depth.size)
    return pd.DataFrame({"depth": depth, "D10": d10, "D50": d50, "D90": d90,
                         "Mean": mean, "Sorting": sorting})


def ensure_files_present(paths):
    missing = [p for p in paths if not os.path.exists(p)]
    if missing and in_colab():
        from google.colab import files
        print("Upload these data files:", missing)
        uploaded = files.upload()
        for name in uploaded:
            print("  uploaded:", name)
    still_missing = [p for p in paths if not os.path.exists(p)]
    if still_missing:
        raise FileNotFoundError(
            f"Data file(s) not found: {still_missing}. Upload them or fix the "
            "'file' entries in CORES (or set USE_DEMO_DATA = True).")


def load_core(cfg):
    path = cfg["file"]
    if path.lower().endswith((".xlsx", ".xls")):
        raw = pd.read_excel(path, sheet_name=cfg.get("sheet", 0))
    else:
        raw = pd.read_csv(path, sep=None, engine="python")  # auto , ; or tab
    raw.columns = [str(c).strip() for c in raw.columns]

    missing = [v for v in COLUMN_MAP.values() if v not in raw.columns]
    if missing:
        raise KeyError(f"{path}: columns {missing} not found. "
                       f"Available columns: {list(raw.columns)}")

    df = raw[list(COLUMN_MAP.values())].copy()
    df.columns = list(COLUMN_MAP.keys())
    df = df.apply(pd.to_numeric, errors="coerce")
    df = df.dropna(subset=["depth"]).sort_values("depth").reset_index(drop=True)
    return df


# -----------------------------------------------------------------------------
# Plotting
# -----------------------------------------------------------------------------
def shared_xlims(data):
    """Same x-range for a parameter in both cores (unless set in X_LIMITS)."""
    lims = {}
    for p, opts in PARAMS.items():
        if X_LIMITS.get(p) is not None:
            lims[p] = X_LIMITS[p]
            continue
        vals = pd.concat([df[p] for df in data.values()]).dropna()
        if opts["log_x"]:
            vals = vals[vals > 0]
            lo, hi = vals.min() / 1.2, vals.max() * 1.2
        else:
            pad = 0.05 * (vals.max() - vals.min() or 1)
            lo, hi = vals.min() - pad, vals.max() + pad
        lims[p] = (lo, hi)
    return lims


def true_scale_figure(top_cm, bottom_cm):
    """Figure whose data panel height is exactly (bottom-top) * scale in cm."""
    panel_h_cm = (bottom_cm - top_cm) * PRINT_CM_PER_CORE_CM
    fig_w_cm = MARGIN_LEFT_CM + PLOT_WIDTH_CM + MARGIN_RIGHT_CM
    fig_h_cm = MARGIN_TOP_CM + panel_h_cm + MARGIN_BOTTOM_CM
    fig = plt.figure(figsize=(fig_w_cm * CM, fig_h_cm * CM))
    ax = fig.add_axes([
        MARGIN_LEFT_CM / fig_w_cm,
        MARGIN_BOTTOM_CM / fig_h_cm,
        PLOT_WIDTH_CM / fig_w_cm,
        panel_h_cm / fig_h_cm,
    ])
    return fig, ax, panel_h_cm


def plot_parameter(df, param, top_cm, bottom_cm, xlim, title,
                   events=(), major=10, minor=1):
    opts = PARAMS[param]
    fig, ax, panel_h_cm = true_scale_figure(top_cm, bottom_cm)

    # Event shading (clipped to the plotted window)
    for et, eb, lab in events:
        a, b = max(et, top_cm), min(eb, bottom_cm)
        if a >= b:
            continue
        ax.axhspan(a, b, color=EVENT_FILL, alpha=EVENT_ALPHA, lw=0, zorder=0)
        if SHOW_EVENT_LABELS and lab:
            ax.text(1.0, (a + b) / 2, f" {lab}", transform=ax.get_yaxis_transform(),
                    ha="left", va="center", fontsize=FONT_SIZE - 1, clip_on=False)

    # Include the points just outside the window so the line runs to the edge
    sel = df[(df["depth"] >= top_cm) & (df["depth"] <= bottom_cm)]
    before = df[df["depth"] < top_cm].tail(1)
    after = df[df["depth"] > bottom_cm].head(1)
    line_df = pd.concat([before, sel, after]).dropna(subset=[param])

    ax.plot(line_df[param], line_df["depth"], color=LINE_COLOR, lw=LINE_WIDTH,
            solid_joinstyle="round", zorder=2)
    ax.plot(sel[param], sel["depth"], "o", ms=MARKER_SIZE, color=LINE_COLOR,
            mec="white", mew=0.3, zorder=3)

    # Depth axis: exact window, depth increasing downward
    ax.set_ylim(bottom_cm, top_cm)
    ax.yaxis.set_major_locator(MultipleLocator(major))
    ax.yaxis.set_minor_locator(MultipleLocator(minor))
    ax.set_ylabel("Depth (cm)")

    # Value axis on top (standard for down-core logs)
    if opts["log_x"]:
        ax.set_xscale("log")
    ax.set_xlim(*xlim)
    if not opts["log_x"]:
        ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.xaxis.set_label_position("top")
    ax.xaxis.tick_top()
    ax.set_xlabel(opts["label"])
    ax.tick_params(axis="both", which="both", direction="out")
    ax.grid(axis="x", which="major", lw=0.3, color="#c8c8c8", zorder=1)
    ax.set_axisbelow(True)

    fig.suptitle(title, fontsize=FONT_SIZE + 1, y=1 - 0.15 / (fig.get_figheight() / CM))
    return fig, panel_h_cm


def save(fig, out_base):
    os.makedirs(os.path.dirname(out_base), exist_ok=True)
    for fmt in FORMATS:
        # NO bbox_inches="tight": it would crop the canvas and break true scale
        fig.savefig(f"{out_base}.{fmt}", dpi=PNG_DPI if fmt == "png" else None)
    plt.close(fig)


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
def main():
    # Load
    data = {}
    if USE_DEMO_DATA:
        for i, (name, cfg) in enumerate(CORES.items()):
            data[name] = make_demo_data(cfg["top_cm"], cfg["bottom_cm"],
                                        cfg["events"], seed=i)
    else:
        ensure_files_present([cfg["file"] for cfg in CORES.values()])
        for name, cfg in CORES.items():
            data[name] = load_core(cfg)

    for name, df in data.items():
        print(f"{name}: {len(df)} samples, depth {df['depth'].min():g}–"
              f"{df['depth'].max():g} cm")

    xlims = shared_xlims(data)
    if os.path.isdir(OUTPUT_DIR):
        shutil.rmtree(OUTPUT_DIR)

    summary = []
    for name, cfg in CORES.items():
        df = data[name]
        top, bottom = cfg["top_cm"], cfg["bottom_cm"]

        if MAKE_FULL_CORE_PLOTS:
            for p in PARAMS:
                fig, h = plot_parameter(
                    df, p, top, bottom, xlims[p], f"{name}  {p}",
                    events=cfg["events"] if SHADE_EVENTS_ON_FULL_CORE else (),
                    major=DEPTH_MAJOR_TICK_CM, minor=DEPTH_MINOR_TICK_CM)
                save(fig, os.path.join(OUTPUT_DIR, name, "full_core",
                                       f"{name}_{p}_full_{top:g}-{bottom:g}cm"))
                summary.append((name, "full core", p, top, bottom, h))

        if MAKE_EVENT_PLOTS:
            for et, eb, lab in cfg["events"]:
                n = ((df["depth"] >= et) & (df["depth"] <= eb)).sum()
                if n == 0:
                    print(f"  WARNING {name} {lab}: no samples between {et} and {eb} cm")
                for p in PARAMS:
                    fig, h = plot_parameter(
                        df, p, et, eb, xlims[p], f"{name} {lab}  {p}",
                        events=[(et, eb, "")],
                        major=EVENT_MAJOR_TICK_CM, minor=EVENT_MINOR_TICK_CM)
                    save(fig, os.path.join(OUTPUT_DIR, name, f"event_{lab}",
                                           f"{name}_{lab}_{p}_{et:g}-{eb:g}cm"))
                    summary.append((name, lab, p, et, eb, h))

    summary = pd.DataFrame(summary, columns=["core", "interval", "parameter",
                                             "top_cm", "bottom_cm",
                                             "panel_height_cm"])
    summary.to_csv(os.path.join(OUTPUT_DIR, "plot_scale_summary.csv"), index=False)
    print(f"\nSaved {len(summary)} plots x {len(FORMATS)} formats to '{OUTPUT_DIR}/'")
    print(f"Scale: 1 cm core = {PRINT_CM_PER_CORE_CM:g} cm on figure. "
          f"Data panel starts {MARGIN_TOP_CM:g} cm below the top edge of each file "
          "and fills exactly top→bottom depth: align that edge with the core photo.")

    zip_path = shutil.make_archive(OUTPUT_DIR, "zip", OUTPUT_DIR)
    if in_colab():
        from google.colab import files
        files.download(zip_path)
    return summary


if __name__ == "__main__":
    main()
