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
There is deliberately no fallback: a photo without one of these is refused,
because guessing its scale (e.g. stretching it to the counted depths) puts the
sediment, and so every event, at the wrong depth.

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
import io
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Polygon, Rectangle
from PIL import Image

from plot_age_depth import EVENT_COLORS, _event_table, read_cores, read_events

Image.MAX_IMAGE_PIXELS = None  # core scans can be very large
CHECK_PIXELS_LONG_SIDE = 6000  # the scale-check figure uses a smaller copy (faster)
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


def _oriented_gray(path, side):
    """Grey photo turned so the core top is on the LEFT (core axis along x)."""
    img = Image.open(path).convert("L")
    rotate = {"left": 0, "top": 90, "right": 180, "bottom": -90}[side]
    if rotate:
        img = img.rotate(rotate, expand=True)
    return np.asarray(img, float)


def _ruler_band(a, bottom=True):
    """Rows of the ruler strip: beyond the darkest row near that photo edge."""
    h = a.shape[0]
    if not bottom:
        a = a[::-1]
    band = a[int(h * 0.7):]
    dark = int(np.argmin(np.median(band, axis=1)))
    return band[dark + 2:]


def _block_edges(rows):
    """Sub-pixel x positions of white/grey block changes along the ruler."""
    p = np.median(rows, axis=0)
    k = max(3, len(p) // 400)
    ps = np.convolve(p, np.ones(k) / k, mode="same")
    lo, hi = np.percentile(ps, 10), np.percentile(ps, 90)
    thr = (lo + hi) / 2
    idx = np.where(np.diff((ps > thr).astype(int)) != 0)[0]
    idx = idx[(idx > k) & (idx < len(ps) - k - 1)]
    edges = idx + (thr - ps[idx]) / (ps[idx + 1] - ps[idx])
    rising = ps[idx + 1] > ps[idx]
    return edges, rising, hi - lo


def _tick_px_per_cm(rows, guess):
    """Independent px/cm from the repeating 1 cm tick pattern (Fourier peak)."""
    p = rows.mean(axis=0)
    win = max(3, int(guess * 2))
    p = p - np.convolve(p, np.ones(win) / win, mode="same")
    n = len(p) * 16
    F = np.abs(np.fft.rfft(p * np.hanning(len(p)), n=n))
    f = np.fft.rfftfreq(n)
    band = (f > 1 / (guess * 1.3)) & (f < 1 / (guess * 0.7))
    return 1 / f[band][np.argmax(F[band])] if band.any() else None


def _fit_ruler(a, block_cm=10.0):
    """Fit the ruler in an image oriented with the core top on the LEFT.

    Returns the best fit over both long sides of the photo, as a dict with
    x0 (pixel of ruler 0), ppc (px per cm), the edges used and a direction
    score: on these rulers 0-10 cm is grey, 10-20 white, 20-30 grey, ..., so
    going down-core grey->white changes fall at 10, 30, 50 ... cm and
    white->grey at 20, 40 ... cm. The score is the fraction of edges that
    follow this pattern (1 = ruler reads the right way, ~0 = reversed).
    """
    best = None
    for bottom in (True, False):  # ruler on either long side of the core
        rows = _ruler_band(a, bottom)
        if rows.shape[0] < 3:
            continue
        edges, rising, contrast = _block_edges(rows)
        if len(edges) < 4:
            continue
        gap = np.median(np.diff(edges))
        keep = np.array([any(abs(abs(edges[j] - edges[i]) - gap) < 0.12 * gap
                             for j in (i - 1, i + 1) if 0 <= j < len(edges))
                         for i in range(len(edges))])
        e, rs = edges[keep], rising[keep]
        if len(e) < 4:
            continue
        n = np.round(e / gap)
        for _ in range(5):  # fit, drop the worst misfit > 1.5 px, refit
            A = np.column_stack([np.ones(len(e)), n * block_cm, np.where(rs, 1.0, -1.0)])
            coef, *_ = np.linalg.lstsq(A, e, rcond=None)
            resid = e - A @ coef
            worst = int(np.argmax(np.abs(resid)))
            if abs(resid[worst]) <= 1.5 or len(e) <= 4:
                break
            e, rs, n = np.delete(e, worst), np.delete(rs, worst), np.delete(n, worst)
        direction = float(np.mean(np.where(rs, n % 2 == 1, n % 2 == 0)))
        score = contrast * len(e)
        if best is None or score > best["score"]:
            best = {"score": score, "x0": coef[0], "ppc": coef[1], "edges": e, "n": n,
                    "resid": resid, "rows": rows, "direction": direction}
    return best


def detect_ruler_scale(path, top_side=None, block_cm=10.0, ruler_start_cm=0.0):
    """Read the depth scale - and which end is the core top - from the ruler.

    The ruler alternates grey and white blocks every `block_cm` (10 cm). The
    block edges along the ruler strip (beside the core) are fitted with a
    straight line (depth -> pixel) after removing edges that do not fit (e.g. a
    differently printed 100 cm mark), and the scale is checked against the 1 cm
    tick spacing measured independently. If `top_side` is not given, both ends
    of the photo are tried and the one where the ruler reads the right way
    (grey 0-10 cm, white 10-20 cm, ...) is taken as the core top. The ruler is
    assumed to start at `ruler_start_cm` at the core-top edge of the photo.

    Returns (depth at the photo's core-top edge, px per cm, report dict with
    "side" = the core-top edge), or (None, None, report) if no ruler is found.
    """
    w, h = Image.open(path).size
    sides = [top_side.lower()] if top_side else (["left", "right"] if w > h else ["top", "bottom"])
    fits = {}
    for side in sides:
        fit = _fit_ruler(_oriented_gray(path, side), block_cm)
        if fit:
            fits[side] = fit
    if not fits:
        return None, None, {"ok": False, "why": "no ruler found"}
    # core top = the end where the ruler reads the right way and its 0 mark sits at
    # the photo edge (on these photos the ruler starts exactly at the core-top edge)
    edge_cm = {k: abs(v["x0"]) / v["ppc"] for k, v in fits.items()}
    side = max(fits, key=lambda k: (round(fits[k]["direction"], 1), -edge_cm[k]))
    f = fits[side]
    tick = _tick_px_per_cm(f["rows"], f["ppc"])
    report = {
        "ok": True, "side": side, "px_per_cm": f["ppc"], "block_edges_used": len(f["edges"]),
        "direction_score": f["direction"],
        "ruler0_from_edge_cm": edge_cm[side],
        "other_end": None if len(fits) < 2 else {
            "direction_score": next(v["direction"] for k, v in fits.items() if k != side),
            "ruler0_from_edge_cm": next(edge_cm[k] for k in fits if k != side)},
        "marks_cm": f"{ruler_start_cm + f['n'].min() * block_cm:g}-"
                    f"{ruler_start_cm + f['n'].max() * block_cm:g}",
        "max_misfit_mm": float(np.abs(f["resid"]).max() / f["ppc"] * 10),
        "tick_px_per_cm": tick,
        "scale_agreement_pct": None if tick is None else float(abs(tick - f["ppc"]) / f["ppc"] * 100),
        "ruler0_offset_px": float(f["x0"]),
    }
    # x0 = pixel position (from the core-top edge) of ruler mark ruler_start_cm
    return ruler_start_cm - f["x0"] / f["ppc"], f["ppc"], report


def photo_scale(path, spec):
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
    t_auto, ppc_auto, rep = detect_ruler_scale(path, spec.get("top_side"))
    if ppc_auto:
        spec["top_side"] = rep["side"]  # core-top end found from the ruler
        return t_auto, ppc_auto, f"ruler (read automatically; core top = {rep['side']} edge)"
    return None, None, "NO DEPTH SCALE - add two ruler readings"


def resolve_photo_depths(cores, photos, photo_depths=None):
    """Depth scale for every core's photo.

    Returns ({core: {"file", "top_side", "top", "px_per_cm", "bottom"}},
    summary DataFrame showing where each scale came from).
    """
    photo_depths = photo_depths or {}
    scales, rows, missing = {}, [], []
    for core in cores:
        spec = {k: None for k in SPEC_KEYS}
        spec.update(photo_depths.get(core, {}))
        path = photos.get(core) or spec.get("file")
        if not path or not Path(path).exists():
            rows.append({"Core": core, "Photo": "-", "Px per cm": None, "Top edge (cm)": None,
                         "Bottom edge (cm)": None, "Depth scale from": "no photo (blank column)"})
            continue
        top, ppc, source = photo_scale(path, spec)
        _, _, side, length, _ = photo_geometry(path, spec.get("top_side"))
        if ppc is None:
            missing.append(core)
            rows.append({"Core": core, "Photo": Path(path).name, "Px per cm": None,
                         "Top edge (cm)": None, "Bottom edge (cm)": None,
                         "Depth scale from": source})
            continue
        scales[core] = {"file": path, "top_side": side, "top": top, "px_per_cm": ppc,
                        "bottom": top + length / ppc}
        rows.append({"Core": core, "Photo": Path(path).name, "Px per cm": round(ppc, 3),
                     "Top edge (cm)": round(top, 2), "Bottom edge (cm)": round(top + length / ppc, 2),
                     "Depth scale from": source})
    summary = pd.DataFrame(rows)
    if missing:
        raise ValueError(
            "No depth scale for: " + ", ".join(missing) + ".\nEnter two ruler readings "
            "(depth in cm and pixel position of two marks on the photo's ruler) for each, "
            "or Px per cm with Top (cm), or Top and Bottom (cm).\n" + summary.to_string(index=False))
    return scales, summary


def load_upright(path, side, max_px=None):
    """The whole photo (nothing cropped), turned so the core top is at the top edge.

    Turning is an exact pixel transpose (no resampling). The photo is kept at its
    original size unless `max_px` (long side) asks for a smaller copy.
    """
    img = Image.open(path).convert("RGB")
    turn = {"top": None, "left": Image.Transpose.ROTATE_270, "right": Image.Transpose.ROTATE_90,
            "bottom": Image.Transpose.ROTATE_180}[side]
    if turn is not None:
        img = img.transpose(turn)
    if max_px and max(img.size) > max_px:
        scale = max_px / max(img.size)
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


def _draw_photo(ax, col, max_px=None, interpolation="none"):
    s = col["scale"]
    if s:
        # interpolation="none": PDF/SVG embed the photo's own pixels, not a resampled copy
        ax.imshow(load_upright(s["file"], s["top_side"], max_px),
                  extent=(col["x0"], col["x1"], col["bottom"], col["top"]),
                  interpolation=interpolation, zorder=1)
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
                     default_width_cm=7.0, ystep=5, fill_alpha=0.18, link_alpha=0.35,
                     height_in=12, title=None, table=True, png_dpi="native",
                     max_png_megapixels=400):
    """Full photos side by side at true scale with the event layers.

    scales: from resolve_photo_depths(); events: from read_event_depths().

    Resolution: in a PDF or SVG every photo is embedded at its ORIGINAL pixel
    size (nothing resampled or cropped), and the event drawing stays vector -
    best for editing in a design canvas. For a PNG, png_dpi="native" picks the
    dpi that keeps the photos' own pixels per cm, limited to max_png_megapixels
    (all 8 full-size cores side by side can be several hundred megapixels); a
    number sets the dpi directly.
    """
    vector = str(out_path).lower().endswith((".pdf", ".svg", ".eps"))
    cores = cores or list(dict.fromkeys(list(scales) + [c for e in events for c in e["cores"]]))
    columns, total_w = _layout(cores, scales, events, gap_cm, default_width_cm)
    ylim = ylim or _auto_ylim(columns, ystep)
    span = max(ylim) - min(ylim)
    fig_w = max(6.0, height_in * (total_w + 6) / span + 2.5)
    fig, ax = plt.subplots(figsize=(fig_w, height_in))

    for col in columns:
        _draw_photo(ax, col, interpolation="none" if vector else "antialiased")
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
                ax.plot([col["x0"], col["x1"]], [d, d], color=color, lw=1.4, zorder=3)
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
                    bbox=dict(facecolor=color, edgecolor="none", alpha=0.45, pad=0.6))

    ax.set_xlim(-1.5, total_w + 5)
    _depth_axis(ax, ylim, ystep)
    ax.grid(axis="y", color="0.88", lw=0.5)
    ax.set_axisbelow(True)
    if table and events:  # same event table and colours as the age-depth figure
        fig.canvas.draw()
        pos = ax.get_position()
        tax = fig.add_axes([pos.x1 + 0.01, pos.y0 + 0.35 * pos.height, 0.22, 0.65 * pos.height])
        all_cores = list(dict.fromkeys(c for e in events for c in e["cores"]))  # whole workbook
        _event_table(tax, events, all_cores, depth_label="Depth in\ncore (cm)")
    if title:
        ax.set_title(title, fontsize=13, pad=60)
    dpi = 300
    if not vector:
        fig.canvas.draw()
        in_per_cm = ax.get_window_extent().height / fig.dpi / span
        ppc_max = max((c["scale"]["px_per_cm"] for c in columns if c["scale"]), default=0)
        cap = np.sqrt(max_png_megapixels * 1e6 / (fig.get_figwidth() * fig.get_figheight()))
        dpi = (min(ppc_max / in_per_cm, cap) if png_dpi == "native" else float(png_dpi)) if ppc_max else 300
        dpi = max(dpi, 150)
        kept = dpi * in_per_cm
        print(f"PNG at {dpi:.0f} dpi = {kept:.1f} px per cm of core "
              f"(photos have up to {ppc_max:.1f} px/cm"
              + (")" if kept >= ppc_max - 0.5 else "; use the PDF/SVG for full resolution)"))
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
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
        _draw_photo(ax, col, max_px=CHECK_PIXELS_LONG_SIDE, interpolation="antialiased")
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
    scales, summary = resolve_photo_depths(cores, photos, spec)
    print(summary.to_string(index=False))
    if args.check:
        plot_scale_check(scales, events, args.check, cores=cores, spec=spec)
    plot_core_photos(scales, events, Path(args.out), cores=cores, ylim=args.ylim,
                     title=args.title)


if __name__ == "__main__":
    main()
