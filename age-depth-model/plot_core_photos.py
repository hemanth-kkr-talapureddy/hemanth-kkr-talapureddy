"""
Core photos with event layers, at the photos' own true scale.

Each photo is used exactly as provided (not cropped): the whole image,
including any ruler and background, is placed on one shared depth axis (cm)
using a linear depth scale for that photo,

    depth (cm) = depth at the photo's top edge + pixel position / pixels per cm

where the pixel position is counted along the core from the photo's core-top
edge. The same pixels-per-cm sets the photo's width, and the plot uses equal x
and y scales, so the photo keeps its real proportions (1 cm across = 1 cm
down). Event layers are then drawn at their depths on that same scale, with a
coloured bar beside the core, and the same event is linked between
neighbouring cores. Event colours are those of plot_age_depth.py, so E1, E2,
... match the age-depth figure.

Use plot_scale_check() first: it draws the code's depth ticks over each full
photo so you can see that they sit on the photo's own ruler.

Depth scale of each photo (first available wins)
------------------------------------------------
1. Two ruler readings: two depths read off the ruler in the photo
   (Ref1 cm, Ref2 cm) and the pixel positions where they are (Ref1 px,
   Ref2 px). Pixel positions are the coordinates your image viewer shows
   along the core: the y coordinate for an upright photo, the x coordinate for
   a sideways one.
2. One ruler reading plus Px per cm (or the photo's stored DPI).
3. Px per cm (or stored DPI) plus Top (cm), the depth at the photo's top edge
   (default 0).
4. Top (cm) and Bottom (cm): the depths at the photo's top and bottom edges.
5. Last resort, APPROXIMATE: the photo is stretched from 0 to the deepest
   counted depth of that core in the varve workbook.

Inputs
------
- Core photos (JPG, PNG, TIFF) as provided. A photo is matched to its core
  when the file name contains the core ID ("GUAC-22A-1G-1.jpg",
  "guac_22a_1g_1.tif", ...). Sideways photos are turned upright; set Top side
  (top/bottom/left/right) for the edge where the core top is (default: top for
  upright photos, left for sideways ones).
- Photo depth table (optional, CSV/Excel): Core and any of Ref1 px, Ref1 cm,
  Ref2 px, Ref2 cm, Px per cm, Top (cm), Bottom (cm), Top side, File.
- Event depths: the "Events" sheet of the varve workbook (two rows per event:
  top then base; the "Depth (in core)" column by default), or a table with
  Core, Event, Top (cm), Base (cm). Use the same depth reference as the photo
  rulers (e.g. both measured from the top of the core liner).

Usage
-----
    python plot_core_photos.py --photos photos/*.jpg --photo-depths photo_depths.csv \
        --events Lachua_Varve_counts.xlsx -o core_photos_events.pdf \
        --check scale_check.png
"""

import argparse
import contextlib
import io
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Polygon, Rectangle
from PIL import Image

from plot_age_depth import EVENT_COLORS, read_cores, read_events

Image.MAX_IMAGE_PIXELS = None  # core scans can be very large
MAX_PIXELS_LONG_SIDE = 6000  # only for drawing; the depth scale uses the full-size photo
SPEC_KEYS = ("ref1_px", "ref1_cm", "ref2_px", "ref2_cm", "px_per_cm", "top", "bottom",
             "top_side", "file")


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

def _read_table(source):
    if isinstance(source, pd.DataFrame):
        return source
    if str(source).lower().endswith((".csv", ".txt")):
        return pd.read_csv(source)
    return pd.read_excel(source, engine="openpyxl")


def _col(columns, test):
    return next((c for c in columns if test(re.sub(r"[^a-z0-9]", "", str(c).lower()))), None)


def read_photo_depths(source):
    """{core: spec} from a CSV/Excel table or DataFrame (see module docstring)."""
    df = _read_table(source)
    cols = {
        "core": _col(df.columns, lambda s: s.startswith("core")),
        "ref1_px": _col(df.columns, lambda s: s.startswith("ref1") and "px" in s),
        "ref1_cm": _col(df.columns, lambda s: s.startswith("ref1") and "cm" in s),
        "ref2_px": _col(df.columns, lambda s: s.startswith("ref2") and "px" in s),
        "ref2_cm": _col(df.columns, lambda s: s.startswith("ref2") and "cm" in s),
        "px_per_cm": _col(df.columns, lambda s: "percm" in s or s.startswith("pxper")),
        "top": _col(df.columns, lambda s: s.startswith("top") and "side" not in s),
        "bottom": _col(df.columns, lambda s: s.startswith("bottom") or s.startswith("base")),
        "top_side": _col(df.columns, lambda s: "side" in s),
        "file": _col(df.columns, lambda s: s.startswith("file") or s.startswith("photo")),
    }
    if cols["core"] is None:
        raise ValueError("The photo depth table needs a 'Core' column.")
    out = {}
    for _, r in df.iterrows():
        if pd.isna(r[cols["core"]]):
            continue
        spec = {}
        for k in SPEC_KEYS:
            c = cols[k]
            v = r[c] if c is not None else None
            if v is None or (not isinstance(v, str) and pd.isna(v)) or str(v).strip() == "":
                spec[k] = None
            elif k in ("top_side", "file"):
                spec[k] = str(v).strip().lower() if k == "top_side" else str(v).strip()
            else:
                spec[k] = float(v)
        out[str(r[cols["core"]]).strip()] = spec
    return out


def read_event_depths(source, depth_col="in core"):
    """Events as [{"name": "E1", "cores": {core: (year or None, top, base)}}, ...].

    Accepts the Events-sheet layout read by plot_age_depth.read_events, or a
    long table with columns Core, Event, Top (cm), Base (cm).
    """
    if not isinstance(source, pd.DataFrame):
        try:
            events = read_events(source, depth_col=depth_col)
        except Exception:
            events = []
        if events:
            return events
    df = _read_table(source)
    core_c = _col(df.columns, lambda s: s.startswith("core"))
    ev_c = _col(df.columns, lambda s: s.startswith("event"))
    top_c = _col(df.columns, lambda s: s.startswith("top"))
    base_c = _col(df.columns, lambda s: s.startswith("base") or s.startswith("bottom"))
    if None in (core_c, ev_c, top_c, base_c):
        raise ValueError("Event file needs an 'Events' sheet (Year/Depth layout) or columns "
                         "Core, Event, Top (cm), Base (cm).")
    events = {}
    for _, r in df.dropna(subset=[core_c, ev_c, top_c, base_c]).iterrows():
        name = str(r[ev_c]).strip()
        name = f"E{int(float(name))}" if re.fullmatch(r"\d+(\.0)?", name) else name
        events.setdefault(name, {})[str(r[core_c]).strip()] = (None, float(r[top_c]), float(r[base_c]))
    return [{"name": n, "cores": c} for n, c in events.items()]


def counts_bottoms(source, depth_col="in core"):
    """{core: deepest counted depth} from the varve workbook, or {} if unavailable."""
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            cores = read_cores(source, depth_col=depth_col)
        return {c: float(df["depth"].max()) for c, df in cores.items()}
    except Exception:
        return {}


def match_photos(files, cores):
    """{core: file} by finding each core ID inside the file names."""
    out = {}
    for core in cores:
        hits = [f for f in files if _key(core) in _key(Path(f).stem)]
        if hits:
            out[core] = sorted(hits, key=lambda f: len(Path(f).stem))[0]
    return out


# ----------------------------------------------------------------------------
# Photo geometry and depth scale
# ----------------------------------------------------------------------------

def photo_geometry(path, top_side=None):
    """Size of the photo as provided, and which edge is the core top.

    Returns (width_px, height_px, side, length_px, across_px), where length is
    measured along the core and across perpendicular to it.
    """
    w, h = Image.open(path).size
    side = (top_side or ("left" if w > h else "top")).lower()
    if side not in ("top", "bottom", "left", "right"):
        raise ValueError(f"Top side must be top, bottom, left or right (got {top_side!r}).")
    length, across = (h, w) if side in ("top", "bottom") else (w, h)
    return w, h, side, length, across


def along_core(px, side, w, h):
    """Viewer coordinate (y for upright, x for sideways) -> pixels from the core-top edge."""
    return {"top": px, "bottom": h - px, "left": px, "right": w - px}[side]


def photo_px_per_cm(path):
    """Pixels per cm from the photo's stored resolution (DPI), or None.

    72 or 96 DPI is a screen default that says nothing about the real scale,
    so values of 100 DPI or less are ignored.
    """
    try:
        dpi = Image.open(path).info.get("dpi")
    except Exception:
        return None
    if not dpi:
        return None
    d = float(dpi[1] if len(dpi) > 1 else dpi[0])
    return d / 2.54 if d > 100 else None


def photo_scale(path, spec, counts_bottom=None):
    """Linear depth scale of one photo: (depth at top edge, px per cm, source).

    `spec` holds any of ref1_px/ref1_cm/ref2_px/ref2_cm, px_per_cm, top,
    bottom, top_side (see module docstring for the order used).
    """
    w, h, side, length, _ = photo_geometry(path, spec.get("top_side"))
    g = lambda k: spec.get(k)
    refs = [(along_core(g(f"ref{i}_px"), side, w, h), g(f"ref{i}_cm")) for i in (1, 2)
            if g(f"ref{i}_px") is not None and g(f"ref{i}_cm") is not None]
    ppc_given, ppc_dpi = g("px_per_cm"), photo_px_per_cm(path)
    if len(refs) == 2:
        (p1, d1), (p2, d2) = refs
        if d1 == d2 or p1 == p2:
            raise ValueError(f"{Path(path).name}: the two ruler readings must differ.")
        ppc = (p2 - p1) / (d2 - d1)
        if ppc <= 0:
            raise ValueError(f"{Path(path).name}: ruler readings give depth decreasing "
                             "down the core; check Top side or the pixel values.")
        return d1 - p1 / ppc, ppc, "ruler (2 readings)"
    if len(refs) == 1 and (ppc_given or ppc_dpi):
        (p1, d1), ppc = refs[0], ppc_given or ppc_dpi
        return d1 - p1 / ppc, ppc, "ruler (1 reading) + " + ("px per cm" if ppc_given else "photo DPI")
    top = g("top") if g("top") is not None else 0.0
    if ppc_given:
        return top, ppc_given, "px per cm + top depth"
    if g("bottom") is not None:
        if g("bottom") <= top:
            raise ValueError(f"{Path(path).name}: Bottom (cm) must be deeper than Top (cm).")
        return top, length / (g("bottom") - top), "top + bottom depths"
    if ppc_dpi:
        return top, ppc_dpi, "photo DPI + top depth"
    if counts_bottom:
        return 0.0, length / counts_bottom, "deepest counted depth (APPROXIMATE - check)"
    raise ValueError(f"{Path(path).name}: no depth scale. Give two ruler readings "
                     "(Ref1/Ref2 px and cm), or Px per cm, or Top and Bottom (cm).")


def resolve_photo_depths(cores, photos, photo_depths=None, counts_bottom=None):
    """Depth scale for every core's photo.

    Returns ({core: {"file", "top_side", "top", "px_per_cm", "bottom"}},
    summary DataFrame showing where each scale came from).
    """
    photo_depths, counts_bottom = photo_depths or {}, counts_bottom or {}
    scales, rows = {}, []
    for core in cores:
        spec = {k: None for k in SPEC_KEYS}
        spec.update(photo_depths.get(core, {}))
        path = photos.get(core) or spec.get("file")
        if not path or not Path(path).exists():
            rows.append({"Core": core, "Photo": "-", "Px per cm": None, "Top edge (cm)": None,
                         "Bottom edge (cm)": None, "Depth scale from": "no photo (blank column)"})
            continue
        top, ppc, source = photo_scale(path, spec, counts_bottom.get(core))
        _, _, side, length, _ = photo_geometry(path, spec.get("top_side"))
        scales[core] = {"file": path, "top_side": side, "top": top, "px_per_cm": ppc,
                        "bottom": top + length / ppc}
        rows.append({"Core": core, "Photo": Path(path).name, "Px per cm": round(ppc, 3),
                     "Top edge (cm)": round(top, 2), "Bottom edge (cm)": round(top + length / ppc, 2),
                     "Depth scale from": source})
    return scales, pd.DataFrame(rows)


def load_upright(path, side):
    """The whole photo (nothing cropped), turned so the core top is at the top edge."""
    img = Image.open(path).convert("RGB")
    rotate = {"top": 0, "left": -90, "right": 90, "bottom": 180}[side]
    if rotate:
        img = img.rotate(rotate, expand=True)
    scale = MAX_PIXELS_LONG_SIDE / max(img.size)
    if scale < 1:  # smaller copy for drawing only; it still spans the same cm extent
        img = img.resize((round(img.size[0] * scale), round(img.size[1] * scale)), Image.LANCZOS)
    return img


def show_photo_pixels(path, start_px=None, end_px=None, out_path=None, height_in=10):
    """Show a photo as provided with its pixel coordinates, to read ruler marks.

    The axes are pixel coordinates exactly as an image viewer reports them
    (x across, y down). start_px/end_px zoom to part of the photo along its
    long side, e.g. around the ruler's 0 and 80 cm marks.
    """
    img = Image.open(path)
    w, h = img.size
    sideways = w > h
    fig, ax = plt.subplots(figsize=(height_in, height_in * 0.35) if sideways
                           else (height_in * 0.35, height_in))
    ax.imshow(img.convert("RGB"), interpolation="nearest")
    lo, hi = start_px or 0, end_px or (w if sideways else h)
    if sideways:
        ax.set_xlim(lo, hi)
    else:
        ax.set_ylim(hi, lo)
    ax.minorticks_on()
    ax.grid(which="major", color="cyan", lw=0.5, alpha=0.7)
    ax.grid(which="minor", color="cyan", lw=0.3, alpha=0.3)
    ax.set_title(f"{Path(path).name}  ({w} x {h} px)  - read ruler marks as "
                 f"{'x' if sideways else 'y'} pixel coordinates", fontsize=9)
    if out_path:
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
    return fig


# ----------------------------------------------------------------------------
# Plotting
# ----------------------------------------------------------------------------

def _layout(cores, scales, events, gap_cm, default_width_cm):
    """x position and depth extent (cm) of each core's column."""
    columns, x = [], 0.0
    for core in cores:
        s = scales.get(core)
        if s:
            _, _, _, length, across = photo_geometry(s["file"], s["top_side"])
            width = across / s["px_per_cm"]  # same scale across as along the core
            top, bottom = s["top"], s["bottom"]
        else:
            depths = [d for ev in events if core in ev["cores"] for d in ev["cores"][core][1:]]
            width, top, bottom = default_width_cm, 0.0, (max(depths) if depths else 10.0)
        columns.append({"core": core, "x0": x, "x1": x + width, "top": top, "bottom": bottom,
                        "scale": s})
        x += width + gap_cm
    return columns, x - gap_cm


def _draw_photo(ax, col):
    s = col["scale"]
    if s:
        ax.imshow(load_upright(s["file"], s["top_side"]),
                  extent=(col["x0"], col["x1"], col["bottom"], col["top"]),
                  interpolation="lanczos", zorder=1)
    else:
        ax.add_patch(Rectangle((col["x0"], col["top"]), col["x1"] - col["x0"],
                               col["bottom"] - col["top"], facecolor="0.92", edgecolor="0.6",
                               lw=0.6, ls="--", zorder=1))
    ax.add_patch(Rectangle((col["x0"], col["top"]), col["x1"] - col["x0"],
                           col["bottom"] - col["top"], facecolor="none", edgecolor="0.25",
                           lw=0.8, zorder=4))


def _depth_axis(ax, ylim, ystep):
    ax.set_ylim(max(ylim), min(ylim))
    ax.set_aspect("equal")  # 1 cm across = 1 cm down
    ax.set_yticks(np.arange(min(ylim), max(ylim) + ystep / 2, ystep))
    ax.set_ylabel("Depth (cm)", fontsize=12)
    ax.set_xticks([])
    for side in ("top", "right", "bottom"):
        ax.spines[side].set_visible(False)


def _auto_ylim(columns, ystep):
    lo = min(c["top"] for c in columns)
    hi = max(c["bottom"] for c in columns)
    return np.ceil(hi / ystep) * ystep, np.floor(lo / ystep) * ystep


def plot_core_photos(scales, events, out_path, cores=None, ylim=None, gap_cm=4.0,
                     default_width_cm=7.0, ystep=5, fill_alpha=0.22, link_alpha=0.35,
                     height_in=12, title=None):
    """Full photos side by side at true scale with the event layers.

    scales: from resolve_photo_depths(); events: from read_event_depths().
    """
    cores = cores or list(dict.fromkeys(list(scales) + [c for e in events for c in e["cores"]]))
    columns, total_w = _layout(cores, scales, events, gap_cm, default_width_cm)
    ylim = ylim or _auto_ylim(columns, ystep)
    span = max(ylim) - min(ylim)
    fig_w = max(6.0, height_in * (total_w + 6) / span + 2.5)
    fig, ax = plt.subplots(figsize=(fig_w, height_in))

    for col in columns:
        _draw_photo(ax, col)
        ax.text((col["x0"] + col["x1"]) / 2, min(ylim) - 0.012 * span, col["core"],
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
            for d in (t, b):  # event top and base lines across the photo
                ax.plot([col["x0"], col["x1"]], [d, d], color=color, lw=0.9, zorder=3)
            ax.add_patch(Rectangle((col["x1"], t), bar_w, b - t, facecolor=color,
                                   edgecolor="none", zorder=3))
        # band linking the event to the next core that has it; a band that skips
        # cores (event missing there) runs behind them, never over a photo
        for c1, c2 in zip(present, present[1:]):
            a, b_ = by_core[c1], by_core[c2]
            (_, t1, s1), (_, t2, s2) = ev["cores"][c1], ev["cores"][c2]
            adjacent = cores.index(c2) - cores.index(c1) == 1
            ax.add_patch(Polygon([(a["x1"] + bar_w, t1), (b_["x0"], t2), (b_["x0"], s2),
                                  (a["x1"] + bar_w, s1)], closed=True, facecolor=color,
                                 alpha=link_alpha if adjacent else link_alpha * 0.6,
                                 edgecolor=color, lw=0.5, zorder=2 if adjacent else 0.5))
        if present:  # label right of the last core with this event, nudged if crowded
            col, (_, t, b) = by_core[present[-1]], ev["cores"][present[-1]]
            lx, ly = col["x1"] + bar_w + 0.4, (t + b) / 2
            while any(abs(lx - px) < 3 and abs(ly - py) < 0.018 * span for px, py in placed):
                ly += 0.018 * span
            placed.append((lx, ly))
            if ly != (t + b) / 2:
                ax.plot([col["x1"] + bar_w, lx], [(t + b) / 2, ly], color="0.4", lw=0.5, zorder=4)
            ax.text(lx, ly, ev["name"], fontsize=7.5, fontweight="bold", va="center",
                    color="0.15", zorder=5, clip_on=False,
                    bbox=dict(facecolor=color, edgecolor="none", alpha=0.6, pad=0.6))

    ax.set_xlim(-1.5, total_w + 5)
    _depth_axis(ax, ylim, ystep)
    ax.grid(axis="y", color="0.88", lw=0.5)
    ax.set_axisbelow(True)
    if title:
        ax.set_title(title, fontsize=13, pad=60)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    print(f"Saved {out_path}")
    return fig


def plot_scale_check(scales, events, out_path, cores=None, spec=None, height_in=14):
    """Check figure: each full photo with the code's depth ticks drawn over it.

    Red ticks every 1 cm (longer and labelled every 5 cm) are drawn from the
    depth scale the code uses; they should sit on the photo's own ruler marks.
    Green crosses mark the ruler readings you entered, and coloured lines the
    event tops and bases.
    """
    cores = [c for c in (cores or list(scales)) if c in scales]
    spec = spec or {}
    columns, total_w = _layout(cores, scales, events, 3.0, 7.0)
    ylim = _auto_ylim(columns, 5)
    span = max(ylim) - min(ylim)
    fig, ax = plt.subplots(figsize=(max(6.0, height_in * (total_w + 4) / span + 2), height_in))
    for col in columns:
        _draw_photo(ax, col)
        s, x0, x1 = col["scale"], col["x0"], col["x1"]
        for d in np.arange(np.ceil(col["top"]), np.floor(col["bottom"]) + 1):
            major = d % 5 == 0
            ax.plot([x1 - (1.2 if major else 0.6), x1], [d, d], color="red",
                    lw=0.9 if major else 0.5, zorder=6)
            if major:  # faint line across the whole photo, through the photo's own ruler
                ax.plot([x0, x1], [d, d], color="red", lw=0.4, alpha=0.55, zorder=6)
                ax.text(x1 + 0.15, d, f"{d:g}", color="red", fontsize=6.5, va="center", zorder=6)
        sp = spec.get(col["core"], {}) or {}
        for i in (1, 2):
            if sp.get(f"ref{i}_cm") is not None and sp.get(f"ref{i}_px") is not None:
                ax.plot((x0 + x1) / 2, sp[f"ref{i}_cm"], "+", color="lime", ms=12, mew=1.8, zorder=7)
        for k, ev in enumerate(events):
            if col["core"] in ev["cores"]:
                _, t, b = ev["cores"][col["core"]]
                for d in (t, b):
                    ax.plot([x0, x1], [d, d], color=event_color(ev["name"], k), lw=0.8, zorder=5)
        ax.text((x0 + x1) / 2, min(ylim) - 0.012 * span,
                f"{col['core']}\n{s['px_per_cm']:.2f} px/cm", ha="center", va="bottom", fontsize=8)
    ax.set_xlim(-1, total_w + 2)
    _depth_axis(ax, ylim, 5)
    ax.set_title("Scale check: red ticks = depth scale used by the code; they should sit on "
                 "the photo's ruler", fontsize=10, pad=34)
    fig.savefig(out_path, dpi=250, bbox_inches="tight")
    print(f"Saved {out_path}")
    return fig


def main():
    p = argparse.ArgumentParser(description="Core photos with event layers at true scale.")
    p.add_argument("--photos", nargs="+", default=[], help="Core photo files (as provided)")
    p.add_argument("--photo-depths",
                   help="CSV/Excel: Core and Ref1 px, Ref1 cm, Ref2 px, Ref2 cm "
                        "(or Px per cm / Top (cm) / Bottom (cm)), optional Top side, File")
    p.add_argument("--events", required=True,
                   help="Varve workbook with an 'Events' sheet, or CSV/Excel with "
                        "Core, Event, Top (cm), Base (cm)")
    p.add_argument("--depth-col", default="in core",
                   help="Depth column to use from an Events sheet (default: in core)")
    p.add_argument("--ylim", nargs=2, type=float, metavar=("MAX_DEPTH", "MIN_DEPTH"))
    p.add_argument("--cores", nargs="+", help="Cores to plot, in order")
    p.add_argument("--check", help="Also save the scale-check figure to this file")
    p.add_argument("--title")
    p.add_argument("-o", "--out", default="core_photos_events.png")
    args = p.parse_args()

    spec = read_photo_depths(args.photo_depths) if args.photo_depths else {}
    events = read_event_depths(args.events, depth_col=args.depth_col)
    cores = args.cores or list(spec) or list(dict.fromkeys(c for e in events for c in e["cores"]))
    photos = match_photos(args.photos, cores)
    photos.update({c: s["file"] for c, s in spec.items() if s.get("file")})
    scales, summary = resolve_photo_depths(cores, photos, spec,
                                           counts_bottoms(args.events, args.depth_col))
    print(summary.to_string(index=False))
    if args.check:
        plot_scale_check(scales, events, args.check, cores=cores, spec=spec)
    plot_core_photos(scales, events, Path(args.out), cores=cores, ylim=args.ylim,
                     title=args.title)


if __name__ == "__main__":
    main()
