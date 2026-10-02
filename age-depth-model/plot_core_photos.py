"""
Core photos with event layers, at true scale.

Places one photo per core side by side on a shared depth axis (cm). Each photo
keeps its real proportions: its width in cm is worked out from the same
cm-per-pixel scale as its length, and the plot uses equal x and y scales, so
1 cm across = 1 cm down. Event layers are drawn on each photo at their depth,
with a coloured bar beside the core, and the same event is linked between
neighbouring cores by a shaded band. Event colours are the ones used in
plot_age_depth.py, so E1, E2, ... match the age-depth figure.

Inputs
------
1. Core photos (JPG, PNG, TIFF). Crop each photo to the core itself. A photo is
   matched to its core when the file name contains the core ID (e.g.
   "GUAC-22A-1G-1.jpg" or "guac_22a_1g_1_photo.tif").
2. Photo depths: the depth (cm) at the top and bottom edge of each photo, as a
   table with columns Core, Top (cm), Bottom (cm) and optionally File,
   Px per cm (to work out the bottom from the photo's scale instead) and
   Top side (top/bottom/left/right: which edge of the photo is the core top;
   default top for upright photos, left for sideways ones).
3. Event depths, either
   - the "Events" sheet layout used by plot_age_depth.py (core ID row, then
     Year / Depth headers, two rows per event: top then base), or
   - a long table with columns Core, Event, Top (cm), Base (cm).
   Use depths measured on the same reference as the photos (by default the
   "Depth (in core)" column).

Usage
-----
    python plot_core_photos.py --photos photos/*.jpg --photo-depths photo_depths.csv \
        --events Lachua_Varve_counts.xlsx -o core_photos_events.pdf
"""

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Polygon, Rectangle
from PIL import Image

from plot_age_depth import EVENT_COLORS, read_events

Image.MAX_IMAGE_PIXELS = None  # core scans can be very large
MAX_PIXELS_LONG_SIDE = 4000  # photos are downsampled to this for plotting


def _key(text):
    """Core ID / file name reduced to lowercase letters and digits for matching."""
    return re.sub(r"[^a-z0-9]", "", str(text).lower())


def event_color(name, k):
    """Same colour as plot_age_depth.py: E1 -> first colour, E2 -> second, ..."""
    m = re.search(r"(\d+)", str(name))
    i = int(m.group(1)) - 1 if m else k
    return EVENT_COLORS[i % len(EVENT_COLORS)]


# ----------------------------------------------------------------------------
# Reading
# ----------------------------------------------------------------------------

def _find(columns, *words):
    for c in columns:
        if all(w in str(c).lower() for w in words):
            return c
    return None


def _read_table(source):
    if isinstance(source, pd.DataFrame):
        return source
    if str(source).lower().endswith((".csv", ".txt")):
        return pd.read_csv(source)
    return pd.read_excel(source, engine="openpyxl")


def read_photo_depths(source):
    """{core: {"top": cm, "bottom": cm or None, "px_per_cm": float or None,
    "file": name or None, "top_side": str or None}} from a CSV/Excel table or DataFrame."""
    df = _read_table(source)
    core_c = _find(df.columns, "core")
    top_c = next((c for c in df.columns
                  if "top" in str(c).lower() and "side" not in str(c).lower()), None)
    bot_c = _find(df.columns, "bottom") or _find(df.columns, "base")
    ppc_c = _find(df.columns, "px") or _find(df.columns, "pixel")
    file_c = _find(df.columns, "file") or _find(df.columns, "photo")
    side_c = _find(df.columns, "side")
    if core_c is None or top_c is None or (bot_c is None and ppc_c is None):
        raise ValueError("Photo depth table needs columns Core, Top (cm) and Bottom (cm) "
                         "(or Px per cm).")
    out = {}
    for _, r in df.iterrows():
        if pd.isna(r[core_c]):
            continue
        num = lambda c: float(r[c]) if c is not None and pd.notna(r[c]) else None
        out[str(r[core_c]).strip()] = {
            "top": num(top_c) or 0.0, "bottom": num(bot_c), "px_per_cm": num(ppc_c),
            "file": str(r[file_c]).strip() if file_c is not None and pd.notna(r[file_c]) else None,
            "top_side": str(r[side_c]).strip().lower() if side_c is not None and pd.notna(r[side_c]) else None,
        }
    return out


def read_event_depths(source, depth_col="in core"):
    """Events as [{"name": "E1", "cores": {core: (year or None, top, base)}}, ...].

    Accepts a long table (Core, Event, Top, Base) or the Events-sheet layout
    read by plot_age_depth.read_events.
    """
    if not isinstance(source, pd.DataFrame):
        try:
            events = read_events(source, depth_col=depth_col)
        except Exception:
            events = []
        if events:
            return events
    df = _read_table(source)
    core_c, ev_c = _find(df.columns, "core"), _find(df.columns, "event")
    top_c, base_c = _find(df.columns, "top"), _find(df.columns, "base") or _find(df.columns, "bottom")
    if None in (core_c, ev_c, top_c, base_c):
        raise ValueError("Event file needs an 'Events' sheet (Year/Depth layout) or columns "
                         "Core, Event, Top (cm), Base (cm).")
    events = {}
    for _, r in df.dropna(subset=[core_c, ev_c, top_c, base_c]).iterrows():
        name = str(r[ev_c]).strip()
        name = f"E{int(float(name))}" if re.fullmatch(r"\d+(\.0)?", name) else name
        events.setdefault(name, {})[str(r[core_c]).strip()] = (None, float(r[top_c]), float(r[base_c]))
    return [{"name": n, "cores": c} for n, c in events.items()]


def load_photo(path, top_side=None):
    """Open a photo and turn it so the core top is at the top edge."""
    img = Image.open(path)
    img = img.convert("RGB")
    w, h = img.size
    side = top_side or ("left" if w > h else "top")
    rotate = {"top": 0, "left": -90, "right": 90, "bottom": 180}[side]
    if rotate:
        img = img.rotate(rotate, expand=True)
    scale = MAX_PIXELS_LONG_SIDE / max(img.size)
    if scale < 1:
        img = img.resize((round(img.size[0] * scale), round(img.size[1] * scale)), Image.LANCZOS)
    return img, (h if side in ("top", "bottom") else w)  # original length in pixels


def match_photos(files, cores):
    """{core: file} by finding each core ID inside the file names."""
    out = {}
    for core in cores:
        hits = [f for f in files if _key(core) in _key(Path(f).stem)]
        if hits:
            out[core] = sorted(hits, key=lambda f: len(Path(f).stem))[0]
    return out


# ----------------------------------------------------------------------------
# Plotting
# ----------------------------------------------------------------------------

def plot_core_photos(photos, photo_depths, events, out_path, cores=None, ylim=None,
                     gap_cm=4.0, default_width_cm=7.0, ystep=5, fill_alpha=0.22,
                     link_alpha=0.35, height_in=12, title=None):
    """Draw the core photos side by side at true scale with the event layers.

    photos: {core: image file path}; photo_depths: from read_photo_depths();
    events: from read_event_depths(); cores: plotting order (default: the
    photo-depth table order, then any other cores that have events).
    """
    if cores is None:
        cores = list(photo_depths)
        cores += [c for ev in events for c in ev["cores"] if c not in cores]
        cores = list(dict.fromkeys(cores))

    # lay the cores out left to right, each as wide as its photo is in cm
    columns = []
    x = 0.0
    for core in cores:
        spec = photo_depths.get(core, {})
        img = None
        path = photos.get(core) or spec.get("file")
        top = spec.get("top", 0.0) or 0.0
        bottom = spec.get("bottom")
        width = default_width_cm
        if path and Path(path).exists():
            img, length_px = load_photo(path, spec.get("top_side"))
            if bottom is None and spec.get("px_per_cm"):
                bottom = top + length_px / spec["px_per_cm"]
            if bottom is None:
                raise ValueError(f"{core}: give the photo's Bottom (cm) or Px per cm.")
            width = (bottom - top) * img.size[0] / img.size[1]  # same cm per pixel both ways
        else:
            if path:
                print(f"{core}: photo '{path}' not found; drawing a blank column.")
            depths = [d for ev in events if core in ev["cores"] for d in ev["cores"][core][1:]]
            bottom = bottom or (max(depths) if depths else 10.0)
        columns.append({"core": core, "img": img, "x0": x, "x1": x + width,
                        "top": top, "bottom": bottom})
        x += width + gap_cm
    total_w = x - gap_cm

    d_top = min(c["top"] for c in columns)
    d_bot = max(c["bottom"] for c in columns)
    if ylim is None:
        ylim = (np.ceil(d_bot / ystep) * ystep, np.floor(d_top / ystep) * ystep)
    depth_span = max(ylim) - min(ylim)
    fig_w = max(6.0, height_in * (total_w + 6) / depth_span + 2.5)
    fig, ax = plt.subplots(figsize=(fig_w, height_in))

    for col in columns:
        if col["img"] is not None:
            ax.imshow(col["img"], extent=(col["x0"], col["x1"], col["bottom"], col["top"]),
                      interpolation="lanczos", zorder=1)
        else:
            ax.add_patch(Rectangle((col["x0"], col["top"]), col["x1"] - col["x0"],
                                   col["bottom"] - col["top"], facecolor="0.92",
                                   edgecolor="0.6", lw=0.6, ls="--", zorder=1))
        ax.add_patch(Rectangle((col["x0"], col["top"]), col["x1"] - col["x0"],
                               col["bottom"] - col["top"], facecolor="none",
                               edgecolor="0.25", lw=0.8, zorder=4))
        ax.text((col["x0"] + col["x1"]) / 2, min(ylim) - 0.012 * depth_span, col["core"],
                ha="left", va="bottom", rotation=40, rotation_mode="anchor", fontsize=9)

    by_core = {c["core"]: c for c in columns}
    bar_w = min(0.9, gap_cm * 0.22)
    placed = []
    for k, ev in enumerate(events):
        color = event_color(ev["name"], k)
        present = [c for c in cores if c in ev["cores"]]
        for core in present:
            col, (_, t, b) = by_core[core], ev["cores"][core]
            ax.add_patch(Rectangle((col["x0"], t), col["x1"] - col["x0"], b - t,
                                   facecolor=color, alpha=fill_alpha, edgecolor="none", zorder=2))
            for d in (t, b):  # crisp top and base lines across the photo
                ax.plot([col["x0"], col["x1"]], [d, d], color=color, lw=0.9, zorder=3)
            ax.add_patch(Rectangle((col["x1"], t), bar_w, b - t, facecolor=color,
                                   edgecolor="none", zorder=3))
        # shaded band linking the event to the next core that has it; a band that
        # skips cores (event missing there) runs behind them, never over a photo
        for c1, c2 in zip(present, present[1:]):
            a, b_ = by_core[c1], by_core[c2]
            (_, t1, s1), (_, t2, s2) = ev["cores"][c1], ev["cores"][c2]
            adjacent = cores.index(c2) - cores.index(c1) == 1
            ax.add_patch(Polygon([(a["x1"] + bar_w, t1), (b_["x0"], t2), (b_["x0"], s2),
                                  (a["x1"] + bar_w, s1)], closed=True, facecolor=color,
                                 alpha=link_alpha if adjacent else link_alpha * 0.6,
                                 edgecolor=color, lw=0.5, zorder=2 if adjacent else 0.5))
        # label right of the last core with this event, nudged if crowded
        if present:
            col, (_, t, b) = by_core[present[-1]], ev["cores"][present[-1]]
            lx, ly = col["x1"] + bar_w + 0.4, (t + b) / 2
            while any(abs(lx - px) < 3 and abs(ly - py) < 0.018 * depth_span for px, py in placed):
                ly += 0.018 * depth_span
            placed.append((lx, ly))
            if ly != (t + b) / 2:
                ax.plot([col["x1"] + bar_w, lx], [(t + b) / 2, ly], color="0.4", lw=0.5, zorder=4)
            ax.text(lx, ly, ev["name"], fontsize=7.5, fontweight="bold", va="center",
                    color="0.15", zorder=5, clip_on=False,
                    bbox=dict(facecolor=color, edgecolor="none", alpha=0.6, pad=0.6))

    ax.set_xlim(-1.5, total_w + 5)
    ax.set_ylim(max(ylim), min(ylim))
    ax.set_aspect("equal")  # 1 cm across = 1 cm down: photos at true proportions
    ax.set_yticks(np.arange(min(ylim), max(ylim) + ystep, ystep))
    ax.set_ylabel("Depth (cm)", fontsize=12)
    ax.set_xticks([])
    for side in ("top", "right", "bottom"):
        ax.spines[side].set_visible(False)
    ax.grid(axis="y", color="0.88", lw=0.5)
    ax.set_axisbelow(True)
    if title:
        ax.set_title(title, fontsize=13, pad=60)

    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    print(f"Saved {out_path}")
    return fig


def main():
    p = argparse.ArgumentParser(description="Core photos with event layers at true scale.")
    p.add_argument("--photos", nargs="+", default=[], help="Core photo files")
    p.add_argument("--photo-depths", required=True,
                   help="CSV/Excel with Core, Top (cm), Bottom (cm) [, File, Px per cm, Top side]")
    p.add_argument("--events", required=True,
                   help="Event depths: workbook with an 'Events' sheet, or CSV/Excel with "
                        "Core, Event, Top (cm), Base (cm)")
    p.add_argument("--depth-col", default="in core",
                   help="Depth column to use from an Events sheet (default: in core)")
    p.add_argument("--ylim", nargs=2, type=float, metavar=("MAX_DEPTH", "MIN_DEPTH"))
    p.add_argument("--cores", nargs="+", help="Cores to plot, in order (default: table order)")
    p.add_argument("--title")
    p.add_argument("-o", "--out", default="core_photos_events.png")
    args = p.parse_args()

    depths = read_photo_depths(args.photo_depths)
    events = read_event_depths(args.events, depth_col=args.depth_col)
    cores = args.cores or list(depths)
    photos = match_photos(args.photos, cores)
    for spec_core, spec in depths.items():  # a File column overrides name matching
        if spec.get("file"):
            photos[spec_core] = spec["file"]
    print(f"{len(photos)} photos matched; {len(events)} events.")
    plot_core_photos(photos, depths, events, Path(args.out), cores=cores, ylim=args.ylim,
                     title=args.title)


if __name__ == "__main__":
    main()
