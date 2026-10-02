"""
Varve-counted age-depth model figure.

Reads varve counts for one or more cores from an Excel workbook and draws a
figure with linked scales:

  (a) age-depth model: age (yr CE) vs depth (cm), one line per core. Every
      event layer is drawn as a thick bar on its core's curve. Each event that
      is correlated between cores (from an "Events" sheet) is shaded as a
      column hanging from the age axis: as wide as the event's age range across
      the cores, with its lower edge running through each core's event base, so
      it matches both the age and depth scales. Columns are labelled E1, E2, ..., and a table beside the plot gives each event's
      age and depth range across the cores
  (b) thickness through time: one strip per core, on the same age axis as (a),
      with each event's age range shaded in the same colour

Event layers are found from the counts themselves: when the same year is
listed on consecutive rows with increasing depth, the depth between them is
an instantaneous deposit (e.g. a turbidite) inside that year. All other depth
between consecutive years is background varve thickness.

Accepted workbook layouts (detected automatically)
--------------------------------------------------
1. Several cores side by side on one sheet, with the core ID in a title row
   above each block of columns:

       GUAC-22A-1G-1 |                 |                  | | GUAC-23A-1G-1 | ...
       Year          | Depth (in core) | Depth (adjusted) | | Year          | ...
       2023          | 14.80           | 0.00             | | 2023          | ...

   When a core has several depth columns, the one containing --depth-col
   (default "adjusted") is used.

2. One sheet per core (sheet name = core ID) with a year and a depth column.

3. Varve thickness only (no year column): one row per varve from the top, in
   mm. Depth is the cumulative thickness and years count back from --top-year.

Header cells are matched case-insensitively: "year"/"age"/"CE"/"AD",
"depth", "thick". If a core appears on more than one sheet (e.g. an extra
"Events" sheet), only the first is used for the counts. Depths are in cm.

Correlated events (optional): a sheet whose name contains "event", in the same
layout, with two rows per event (its top and base in every core). See
read_events().

Usage
-----
    python plot_age_depth.py Lachua_Varve_counts.xlsx -o age_depth_model.pdf
    python plot_age_depth.py counts.xlsx --xlim 1700 2050 --ylim 90 0 --event-mm 10
"""

import argparse
import io
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Polygon
from matplotlib.ticker import FuncFormatter

PLAIN = FuncFormatter(lambda v, _: f"{v:g}")  # 0.1, 1, 10 instead of 10^-1

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

# Soft background shades for correlated events (cycled; neighbours always differ).
EVENT_COLORS = ["#F4A261", "#2A9D8F", "#E9C46A", "#6D8EAD", "#E76F51", "#8AB17D",
                "#B5838D", "#9C89B8"]

YEAR_RE = re.compile(r"year|\bage\b|\bce\b|\bad\b", re.I)
DEPTH_RE = re.compile(r"depth", re.I)
THICK_RE = re.compile(r"thick", re.I)


# ----------------------------------------------------------------------------
# Reading
# ----------------------------------------------------------------------------

def _header_row(raw, max_rows=10):
    """Index of the first row holding both a year header and a depth/thickness header."""
    for i in range(min(max_rows, len(raw))):
        cells = [str(v) for v in raw.iloc[i] if pd.notna(v)]
        if any(YEAR_RE.search(c) for c in cells) and \
                any(DEPTH_RE.search(c) or THICK_RE.search(c) for c in cells):
            return i
    return None


def _column_groups(raw, hdr):
    """Split columns into per-core blocks: [(core_name or None, [col, ...]), ...]."""
    names = raw.iloc[hdr - 1].ffill() if hdr > 0 else pd.Series([None] * raw.shape[1])
    groups, current = [], None
    for c in range(raw.shape[1]):
        label = raw.iat[hdr, c]
        if pd.isna(label):
            current = None  # a blank header cell ends a block
            continue
        name = names.iat[c] if pd.notna(names.iat[c]) else None
        if YEAR_RE.search(str(label)) or current is None or current[0] != name:
            current = (name, [])
            groups.append(current)
        current[1].append(c)
    return groups


def _series(year, depth):
    """Numeric (year, depth) rows in file order; rows with blanks dropped."""
    df = pd.DataFrame({"year": pd.to_numeric(pd.Series(year), errors="coerce").to_numpy(),
                       "depth": pd.to_numeric(pd.Series(depth), errors="coerce").to_numpy()})
    return df.dropna().reset_index(drop=True)


def read_cores(source, top_year=None, depth_col="adjusted", sheet=None):
    """Return {core_id: DataFrame(year, depth)} from an Excel workbook.

    `source` is a file path, bytes, or a file-like object.
    """
    if isinstance(source, (bytes, bytearray)):
        source = io.BytesIO(source)
    book = pd.ExcelFile(source, engine="openpyxl")
    sheets = [sheet] if sheet else book.sheet_names
    cores = {}
    for sheet_name in sheets:
        raw = book.parse(sheet_name, header=None)
        raw = raw.dropna(how="all").dropna(axis=1, how="all").reset_index(drop=True)
        raw.columns = range(raw.shape[1])
        hdr = _header_row(raw) if not raw.empty else None
        if hdr is None:
            print(f"Skipping sheet '{sheet_name}': no Year/Depth header row found.")
            continue
        body = raw.iloc[hdr + 1:]
        groups = [g for g in _column_groups(raw, hdr)
                  if any(YEAR_RE.search(str(raw.iat[hdr, c])) or
                         THICK_RE.search(str(raw.iat[hdr, c])) for c in g[1])]
        for k, (name, cols) in enumerate(groups):
            name = str(name).strip() if name else (sheet_name if len(groups) == 1
                                                   else f"{sheet_name} {k + 1}")
            if name in cores:
                continue  # first occurrence wins (e.g. skip a repeated "Events" sheet)
            labels = {c: str(raw.iat[hdr, c]) for c in cols}
            year_c = [c for c in cols if YEAR_RE.search(labels[c])]
            depth_c = [c for c in cols if DEPTH_RE.search(labels[c])]
            thick_c = [c for c in cols if THICK_RE.search(labels[c])]
            if year_c and depth_c:
                preferred = [c for c in depth_c if depth_col and depth_col.lower() in labels[c].lower()]
                d = (preferred or depth_c)[0]
                df = _series(body[year_c[0]], body[d])
            elif thick_c:
                if top_year is None:
                    raise ValueError(f"'{name}' has only varve thickness; "
                                     "pass top_year (year of the topmost varve).")
                t_cm = pd.to_numeric(body[thick_c[0]], errors="coerce").dropna().to_numpy() / 10
                depth = np.concatenate([[0.0], np.cumsum(t_cm)])
                df = _series(top_year - np.arange(len(depth)), depth)
            else:
                continue
            if len(df) >= 2:
                cores[name] = df
    if not cores:
        raise ValueError("No core data found. Check that each core has a 'Year' and a "
                         "'Depth' column header.")
    return cores


def split_varves_events(df):
    """Split one core's (year, depth) rows into background varves and event layers.

    Returns two DataFrames with columns year, top, base, thickness (mm).
    Consecutive rows with the same year and increasing depth are an event layer;
    the depth between consecutive different years is that varve's thickness.
    """
    year, depth = df["year"].to_numpy(), df["depth"].to_numpy()
    rows_v, rows_e = [], []
    for i in range(1, len(df)):
        inc = depth[i] - depth[i - 1]
        if inc <= 0:
            continue
        rec = (year[i - 1], depth[i - 1], depth[i], inc * 10.0)
        (rows_e if year[i] == year[i - 1] else rows_v).append(rec)
    cols = ["year", "top", "base", "thickness"]
    return pd.DataFrame(rows_v, columns=cols), pd.DataFrame(rows_e, columns=cols)


def read_events(source, depth_col="adjusted", sheet=None):
    """Read correlated event layers from an "Events" sheet, if the workbook has one.

    The sheet uses the same layout as the counts (core ID row, then Year/Depth
    headers). Rows come in pairs, one pair per event: the top row and the base
    row of that event in every core, e.g. rows 1-2 = E1, rows 3-4 = E2. A core
    with no layer for an event has blank cells. Reading stops at the first
    fully blank row.

    Returns [{"name": "E1", "cores": {core_id: (year, top, base)}}, ...], or []
    when there is no such sheet.
    """
    if isinstance(source, (bytes, bytearray)):
        source = io.BytesIO(source)
    if hasattr(source, "seek"):
        source.seek(0)
    book = pd.ExcelFile(source, engine="openpyxl")
    names = [sheet] if sheet else [n for n in book.sheet_names if "event" in n.lower()]
    for sheet_name in names:
        raw = book.parse(sheet_name, header=None)
        raw = raw.dropna(axis=1, how="all")
        raw.columns = range(raw.shape[1])
        hdr = _header_row(raw)
        if hdr is None:
            continue
        series = {}
        for name, cols in _column_groups(raw, hdr):
            labels = {c: str(raw.iat[hdr, c]) for c in cols}
            year_c = [c for c in cols if YEAR_RE.search(labels[c])]
            depth_c = [c for c in cols if DEPTH_RE.search(labels[c])]
            if name is None or not (year_c and depth_c):
                continue
            preferred = [c for c in depth_c if depth_col and depth_col.lower() in labels[c].lower()]
            series[str(name).strip()] = (
                pd.to_numeric(raw[year_c[0]].iloc[hdr + 1:], errors="coerce").to_numpy(),
                pd.to_numeric(raw[(preferred or depth_c)[0]].iloc[hdr + 1:], errors="coerce").to_numpy())
        if not series:
            continue
        years = np.array([y for y, _ in series.values()])
        blank = np.isnan(years).all(axis=0)
        n_rows = int(np.argmax(blank)) if blank.any() else years.shape[1]
        events = []
        for k in range(n_rows // 2):
            per_core = {}
            for core, (y, d) in series.items():
                y1, y2, d1, d2 = y[2 * k], y[2 * k + 1], d[2 * k], d[2 * k + 1]
                if not np.isnan([y1, d1, d2]).any() and d2 > d1:
                    per_core[core] = (y1, d1, d2)
            if per_core:
                events.append({"name": f"E{k + 1}", "cores": per_core})
        if events:
            print(f"Read {len(events)} correlated events from sheet '{sheet_name}'.")
            return events
    return []


def event_envelope(event, top=None):
    """Polygon (list of (age, depth)) enclosing one correlated event in every core.

    Each core's layer is the rectangle [year - 0.5, year + 0.5] x [top, base]
    (a layer fills its whole varve year); the envelope is the convex hull of
    those corners, so it touches every core's event top and base exactly.

    With `top` (a depth, e.g. 0), the polygon is extended straight up to that
    depth: a column as wide as the event's age range across the cores, hanging
    from the age axis, whose lower edge runs through the cores' event bases.
    """
    pts = {(y + dx, d) for y, t, b in event["cores"].values()
           for dx in (-0.5, 0.5) for d in (t, b)}
    if top is not None:
        xs = [x for x, _ in pts]
        pts |= {(min(xs), top), (max(xs), top)}
    pts = sorted(pts)

    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    lower, upper = [], []
    for p in pts:  # Andrew's monotone chain
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)
    for p in reversed(pts):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)
    return lower[:-1] + upper[:-1]


# ----------------------------------------------------------------------------
# Plotting
# ----------------------------------------------------------------------------

def _nice_limits(lo, hi, step):
    return np.floor(lo / step) * step, np.ceil(hi / step) * step


def _auto_limits(cores, xlim, ylim, xstep, ystep):
    years = pd.concat([d["year"] for d in cores.values()])
    depths = pd.concat([d["depth"] for d in cores.values()])
    if xlim is None:
        xlim = _nice_limits(years.min(), years.max(), xstep)
    if ylim is None:
        ylim = (_nice_limits(0, depths.max(), ystep)[1], 0)
    return xlim, ylim


def _style_age_depth(ax, xlim, ylim, xstep, ystep):
    ax.set_xlim(xlim)
    ax.set_ylim(max(ylim), min(ylim))  # depth increases downward
    ax.set_xticks(np.arange(min(xlim), max(xlim) + xstep, xstep))
    ax.set_yticks(np.arange(min(ylim), max(ylim) + ystep, ystep))
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position("top")
    ax.set_xlabel("Age (yr CE)", fontsize=12, labelpad=8)
    ax.set_ylabel("Depth (cm)", fontsize=12)
    ax.grid(color="0.9", lw=0.6)
    ax.set_axisbelow(True)


def plot_age_depth(cores, out_path, xlim=None, ylim=None, style="line",
                   title=None, xstep=50, ystep=10):
    """Single panel: age-depth curves only (like the original Excel chart)."""
    xlim, ylim = _auto_limits(cores, xlim, ylim, xstep, ystep)
    fig, ax = plt.subplots(figsize=(11, 8))
    for i, (name, df) in enumerate(cores.items()):
        ax.plot(df["year"], df["depth"], color=CORE_COLORS[i % len(CORE_COLORS)], lw=1.6,
                label=name, drawstyle="steps-post" if style == "steps" else "default")
    _style_age_depth(ax, xlim, ylim, xstep, ystep)
    ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.6), frameon=False, fontsize=11)
    if title:
        ax.set_title(title, fontsize=13, pad=40)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    print(f"Saved {out_path}")
    return fig


def _event_table(ax, events, cores, depth_label="Depth\n(cm)"):
    """Table of each correlated event's age and depth range across the cores."""
    rows = []
    for ev in events:
        y = [v[0] for v in ev["cores"].values()]
        tops = [v[1] for v in ev["cores"].values()]
        bases = [v[2] for v in ev["cores"].values()]
        age = f"{min(y):.0f}" if min(y) == max(y) else f"{min(y):.0f}–{max(y):.0f}"
        rows.append([ev["name"], age, f"{min(tops):.1f}–{max(bases):.1f}",
                     f"{len(ev['cores'])}/{len(cores)}"])
    ax.axis("off")
    table = ax.table(cellText=rows, colLabels=["Event", "Age\n(yr CE)", depth_label, "Cores"],
                     loc="upper center", cellLoc="center", colWidths=[0.2, 0.32, 0.32, 0.18])
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1, 0.85)
    for (r, c), cell in table.get_celld().items():
        cell.set_edgecolor("0.8")
        cell.set_linewidth(0.5)
        if r == 0:
            cell.set_text_props(fontweight="bold")
            cell.set_facecolor("0.93")
            cell.set_height(cell.get_height() * 1.8)
        elif c == 0:  # event name in that event's shade, as in the plot
            cell.set_facecolor(EVENT_COLORS[(r - 1) % len(EVENT_COLORS)])
            cell.set_alpha(0.45)
            cell.set_text_props(fontweight="bold")


def plot_age_depth_panel(cores, out_path, events=None, xlim=None, ylim=None, style="line",
                         title=None, xstep=50, ystep=10, event_mm=None):
    """(a) age-depth curves with event layers and shaded event columns hanging
    from the age axis, plus a table of event ages and depths; (b) thickness through time,
    one strip per core, with each event's age range shaded."""
    names = list(cores)
    events = events or []
    colors = {n: CORE_COLORS[i % len(CORE_COLORS)] for i, n in enumerate(names)}
    parts = {n: split_varves_events(df) for n, df in cores.items()}
    xlim, ylim = _auto_limits(cores, xlim, ylim, xstep, ystep)
    all_t = pd.concat([p["thickness"] for v_e in parts.values() for p in v_e])
    tlim = (10 ** np.floor(np.log10(all_t.min())), 10 ** np.ceil(np.log10(all_t.max())))

    n = len(names)
    fig = plt.figure(figsize=(12, 8.5 + 0.45 * n))
    outer = fig.add_gridspec(2, 2, width_ratios=[3.2, 1] if events else [1, 0.001],
                             height_ratios=[3, 0.16 * n + 0.4], wspace=0.04, hspace=0.08)
    ax = fig.add_subplot(outer[0, 0])
    strips = outer[1, 0].subgridspec(n, 1, hspace=0)
    axs = [fig.add_subplot(strips[i], sharex=ax) for i in range(n)]

    # (a) age-depth model
    for name, df in cores.items():
        ax.plot(df["year"], df["depth"], color=colors[name], lw=1.5, label=name, zorder=3,
                drawstyle="steps-post" if style == "steps" else "default")
        # every event layer in the counts, drawn as a thick bar over the curve
        evt = parts[name][1]
        ax.vlines(evt["year"], evt["top"], evt["base"], color=colors[name], lw=4.5,
                  alpha=0.45, zorder=2)
    _style_age_depth(ax, xlim, ylim, xstep, ystep)

    # correlated events: a shaded column hanging from the age axis for each event.
    # Its width is the event's age range across the cores (each layer spanning
    # its varve year, +/- 0.5 yr) and its lower edge runs through the cores'
    # event bases, so it matches both the age and the depth scale.
    placed = []
    depth_top = min(ylim)
    for k, ev in enumerate(events):
        color = EVENT_COLORS[k % len(EVENT_COLORS)]
        hull = event_envelope(ev, top=depth_top)
        ax.add_patch(Polygon(hull, closed=True, facecolor=color, alpha=0.2,
                             edgecolor="none", zorder=1))
        ax.add_patch(Polygon(hull, closed=True, facecolor="none", edgecolor=color,
                             alpha=0.9, lw=0.9, zorder=1))
        # coloured tick on the age axis marking the event's age range
        x_lo, x_hi = min(x for x, _ in hull), max(x for x, _ in hull)
        ax.plot([x_lo, x_hi], [depth_top, depth_top], color=color, lw=5,
                solid_capstyle="butt", clip_on=False, zorder=6)
        # label just right of the column's deepest corner, nudged if crowded
        x0, y0 = max(hull, key=lambda p: (p[1], p[0]))
        lx, ly = x0 + 2, y0
        while any(abs(lx - px) < 13 and abs(ly - py) < 1.6 for px, py in placed):
            ly += 1.6
        placed.append((lx, ly))
        if ly != y0:  # leader line from a nudged label back to its column
            ax.plot([x0, lx], [y0, ly], color="0.4", lw=0.5, zorder=4)
        ax.text(lx, ly, ev["name"], fontsize=8, fontweight="bold", color="0.15",
                va="center", zorder=5, clip_on=True,
                bbox=dict(facecolor=color, edgecolor="none", alpha=0.45, pad=0.6))

    tax = None
    if events:
        tax = fig.add_subplot(outer[0, 1])
        _event_table(tax, events, cores)

    # (b) thickness through time, one strip per core
    for a, name in zip(axs, names):
        var, evt = parts[name]
        for k, ev in enumerate(events):  # event age range, same shade as panel a
            yrs = [v[0] for v in ev["cores"].values()]
            a.axvspan(min(yrs) - 0.5, max(yrs) + 0.5, color=EVENT_COLORS[k % len(EVENT_COLORS)],
                      alpha=0.28, lw=0, zorder=0)
        a.fill_between(var["year"], tlim[0], var["thickness"], step="mid",
                       color=colors[name], alpha=0.25, lw=0)
        a.plot(var["year"], var["thickness"], color=colors[name], lw=0.7, drawstyle="steps-mid")
        if len(evt):
            a.vlines(evt["year"], tlim[0], evt["thickness"], color=colors[name], lw=1.2)
            a.plot(evt["year"], evt["thickness"], "o", ms=3, color=colors[name])
        a.set_yscale("log")
        a.set_ylim(tlim)
        a.set_yticks([v for v in (1, 10, 100) if tlim[0] < v < tlim[1]])
        a.yaxis.set_major_formatter(PLAIN)
        a.minorticks_off()
        a.tick_params(axis="y", labelsize=7, length=2, pad=1)
        a.text(1.005, 0.5, name, transform=a.transAxes, va="center", fontsize=9,
               color=colors[name])
        a.grid(axis="x", color="0.9", lw=0.6)
        a.set_axisbelow(True)
        if event_mm:
            a.axhline(event_mm, color="0.35", ls="--", lw=0.7)
        if a is not axs[-1]:
            a.tick_params(labelbottom=False, bottom=False)
    axs[-1].set_xlabel("Age (yr CE)", fontsize=12)
    axs[n // 2].set_ylabel("Thickness\n(mm, log)", fontsize=10,
                           rotation=0, ha="right", va="center", labelpad=10)

    for a, label in ((ax, "a"), (axs[0], "b")):
        a.text(0.012, 0.985 if a is ax else 0.92, label, transform=a.transAxes,
               fontsize=13, fontweight="bold", va="top", zorder=6,
               bbox=dict(facecolor="white", edgecolor="none", pad=1.5))

    handles = ax.get_legend_handles_labels()[0] + [
        Line2D([], [], color="0.5", lw=4.5, alpha=0.45, label="Event layer")]
    if events:
        handles.append(Patch(facecolor=EVENT_COLORS[0], edgecolor=EVENT_COLORS[0], alpha=0.4,
                             label="Event age range\n(from age axis to\nevent base)"))
    if tax is not None:  # below the event table, clear of the curves and labels
        tax.legend(handles=handles, loc="lower center", frameon=True, edgecolor="0.8",
                   fontsize=8, ncol=1, title="Core", title_fontsize=8.5)
    else:
        ax.legend(handles=handles, loc="lower right", frameon=True, framealpha=0.9,
                  edgecolor="0.8", fontsize=9, ncol=2, title="Core", title_fontsize=9)
    if title:
        fig.suptitle(title, fontsize=14)

    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    print(f"Saved {out_path}")
    return fig


def summarize(cores):
    for name, df in cores.items():
        var, evt = split_varves_events(df)
        print(f"  {name}: {df['year'].min():.0f}-{df['year'].max():.0f} CE, "
              f"{df['depth'].max():.1f} cm, {len(var)} varves "
              f"(median {var['thickness'].median():.1f} mm), {len(evt)} event layers "
              f"(max {evt['thickness'].max() if len(evt) else 0:.0f} mm)")


def main():
    p = argparse.ArgumentParser(description="Plot a varve-counted age-depth model from Excel.")
    p.add_argument("excel", help="Varve counting workbook (.xlsx)")
    p.add_argument("-o", "--out", default="age_depth_model.png",
                   help="Output figure (.png, .pdf or .svg). Default: age_depth_model.png")
    p.add_argument("--xlim", nargs=2, type=float, metavar=("MIN_YEAR", "MAX_YEAR"),
                   help="Age axis limits, e.g. 1700 2050")
    p.add_argument("--ylim", nargs=2, type=float, metavar=("MAX_DEPTH", "MIN_DEPTH"),
                   help="Depth axis limits in cm, e.g. 90 0")
    p.add_argument("--depth-col", default="adjusted",
                   help="If a core has several depth columns, use the one whose header "
                        "contains this text (default: adjusted)")
    p.add_argument("--sheet", help="Read only this sheet (default: all sheets)")
    p.add_argument("--events-sheet",
                   help="Sheet with correlated events (default: a sheet named like 'Events')")
    p.add_argument("--layout", choices=["panel", "simple"], default="panel",
                   help="'panel' (default): age-depth with events + thickness strips; "
                        "'simple': age-depth curves only")
    p.add_argument("--style", choices=["line", "steps"], default="line",
                   help="Age-depth line style (default: line, as in Excel)")
    p.add_argument("--event-mm", type=float,
                   help="Draw a dashed thickness threshold (mm) on the thickness strips")
    p.add_argument("--top-year", type=int,
                   help="Year of the topmost varve (only for thickness-only sheets)")
    p.add_argument("--title", help="Optional figure title")
    args = p.parse_args()

    cores = read_cores(args.excel, top_year=args.top_year, depth_col=args.depth_col,
                       sheet=args.sheet)
    summarize(cores)
    events = read_events(args.excel, depth_col=args.depth_col, sheet=args.events_sheet)
    kw = dict(xlim=args.xlim, ylim=args.ylim, style=args.style, title=args.title)
    if args.layout == "panel":
        plot_age_depth_panel(cores, Path(args.out), events=events, event_mm=args.event_mm, **kw)
    else:
        plot_age_depth(cores, Path(args.out), **kw)


if __name__ == "__main__":
    main()
