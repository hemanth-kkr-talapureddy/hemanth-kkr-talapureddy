# -*- coding: utf-8 -*-
"""
GUAC-24A and GUAC-29A event panels: core photo on a true depth scale + grain size in phi.
Works in Google Colab or plain Python.

HOW TO RUN IN GOOGLE COLAB
  1. Check the file names in the SETTINGS block below (they must match your files).
  2. Run all cells. When asked, upload the files listed in SETTINGS.
  3. Four figures are shown, saved to /content/output/ and downloaded as a zip.
  Only the data in these files is plotted: event depths and ages come from the Excel file,
  grain size from the two CSV files, depth scale from the ruler on each core photo.

Every figure (2 cores x E2/E3) shares:
  * the same 9 cm depth window  -> identical cm-per-inch vertical scale
  * the same x-axis limits and the same panel widths (set in inches)
  * 24 pt Liberation Sans for every label, tick and legend entry
  * a height of 575 pt (7.99 in), same as the cartoon PDFs

Panels: core photo | mean grain size (phi): continuous line through the samples, with a
        mean +/- 1 sorting envelope | class fractions (%), continuous
Missing samples inside an event are filled by linear interpolation between neighbours.
Before plotting, run_checks() tests the depth calibration and the data; any problem stops the script.
"""
import glob
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image

Image.MAX_IMAGE_PIXELS = None
IN_COLAB = "google.colab" in sys.modules
if not IN_COLAB:
    import matplotlib
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib import font_manager

WORK = Path("/content") if IN_COLAB else Path(__file__).resolve().parent
OUT = WORK / "output"

# ================================================================== SETTINGS: your file names
FILES = {
    "csv_29A": "GUAC_29A_GSD_Removedless1um_FolkWard_Stats.csv",
    "csv_24A": "GUAC_24A_GSD_Removedless1um_FolkWard_Stats.csv",
    "img_29A": "GUAC-29A-1G-1-W.jpg",          # core photo with the ruler
    "img_24A": "GUAC-24A-1G-1-W.jpg",          # core photo with the ruler (or GUAC-24A_realscale.pdf)
    "events":  "Lachua_Graine_Siza_Samples_N_Events_with_Ages.xlsx",   # event depths + ages
}
EVENTS_TO_PLOT = ["E2", "E3"]
# Colab renames a re-uploaded file "name (1).csv"; such copies and "name_1.csv" are accepted too.


# ------------------------------------------------------------------ font (24 pt Liberation Sans)
def ensure_font(name="Liberation Sans"):
    if any(f.name == name for f in font_manager.fontManager.ttflist):
        return name
    if IN_COLAB:  # Colab: install the font, then register it
        subprocess.run(["apt-get", "-qq", "install", "-y", "fonts-liberation"],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        for f in glob.glob("/usr/share/fonts/truetype/liberation/LiberationSans-*.ttf"):
            font_manager.fontManager.addfont(f)
        if any(f.name == name for f in font_manager.fontManager.ttflist):
            return name
    print(f"WARNING: {name} not found, using DejaVu Sans (size is still 24 pt)")
    return "DejaVu Sans"


FONT = 24
plt.rcParams.update({
    "font.family": ensure_font(),
    "font.size": FONT, "axes.titlesize": FONT, "axes.labelsize": FONT,
    "xtick.labelsize": FONT, "ytick.labelsize": FONT, "legend.fontsize": FONT,
    "axes.linewidth": 1.5, "xtick.major.width": 1.5, "ytick.major.width": 1.5,
    "xtick.major.size": 7, "ytick.major.size": 7, "xtick.minor.size": 4, "ytick.minor.size": 4,
    "svg.fonttype": "none", "pdf.fonttype": 42,
})
INK = "#1F2D3A"


# ------------------------------------------------------------------ input files
def _resolve(name):
    """Path of the file called `name` (or Colab's 'name (1).ext' / 'name_1.ext' copy), else None."""
    stem, ext = Path(name).stem, Path(name).suffix
    for root in (WORK, WORK / "data"):
        exact = root / name
        if exact.exists():
            return exact
        copies = sorted(glob.glob(str(root / f"{glob.escape(stem)} (*){ext}"))
                        + glob.glob(str(root / f"{glob.escape(stem)}_[0-9]{ext}")))
        if copies:
            return Path(copies[-1])
    return None


def find_inputs():
    """Paths of the files named in FILES; in Colab, ask for an upload of any that are missing."""
    found = {k: _resolve(v) for k, v in FILES.items()}
    if None in found.values() and IN_COLAB:
        from google.colab import files
        print("Please upload:\n  " + "\n  ".join(FILES[k] for k, v in found.items() if v is None))
        files.upload()
        found = {k: _resolve(v) for k, v in FILES.items()}
    missing = [FILES[k] for k, v in found.items() if v is None]
    if missing:
        raise FileNotFoundError(f"file(s) not found: {missing}. Upload them, or change the names "
                                f"in FILES (SETTINGS block) to your file names.")
    for k, v in found.items():
        print(f"{k:8s}: {v.name}")
    if any(found[k].suffix.lower() == ".pdf" for k in ("img_29A", "img_24A")):
        ensure_pymupdf()
    return found


EVENT_COLS = ["event", "age_CE", "age_unc"]


def ensure_pymupdf():
    try:
        import pymupdf  # noqa: F401
    except ImportError:
        subprocess.run([sys.executable, "-m", "pip", "install", "-q", "pymupdf"], check=True)


def image_from_pdf(path):
    """Largest image on page 1, oriented exactly as the PDF displays it."""
    import pymupdf
    with pymupdf.open(path) as d:
        page = d[0]
        info = max(page.get_image_info(xrefs=True),
                   key=lambda i: (i["bbox"][2] - i["bbox"][0]) * (i["bbox"][3] - i["bbox"][1]))
        raw = d.extract_image(info["xref"])["image"]
    import io
    arr = np.asarray(Image.open(io.BytesIO(raw)).convert("RGB"))
    a, b, c, dd, _, _ = info["transform"]
    assert abs(b) < 1e-6 and abs(c) < 1e-6, "rotated image in PDF is not supported"
    if dd < 0:      # PDF stores the image bottom-up and flips it for display
        arr = arr[::-1]
    if a < 0:
        arr = arr[:, ::-1]
    return Image.fromarray(np.ascontiguousarray(arr))


# ------------------------------------------------------------------ depth calibration
# Works on any core photo with the grey/white ruler along one side (e.g. the original
# GUAC-24A-1G-1-W.jpg / GUAC-29A-1G-1-W.jpg at full resolution, or smaller copies):
#   1. the long side of the photo is the depth axis; the ruler strip is found automatically
#   2. the 10 cm grey/white ruler blocks give a first scale; the ruler's 0 cm must sit at one
#      end of the photo, which also tells which end is the top
#   3. every cm and half-cm tick line on the ruler is then used for the final straight-line fit
# GUAC-24A_realscale.pdf (image cropped at 10.5 cm) is calibrated from its ticks directly.
REF24_PX_PER_CM, REF24_TOP_CM = 104.209, 10.512   # expected values for the realscale PDF image


def _runs(mask):
    """(start, end) of each run of True values."""
    m = np.r_[False, np.asarray(mask, bool), False]
    d = np.flatnonzero(np.diff(m.astype(np.int8)))
    return list(zip(d[::2], d[1::2]))


def _tick_centres(prof, dark_thr):
    """Centre (px) of each dark tick line in a ruler brightness profile."""
    out = []
    for i, j in _runs(prof < dark_thr):
        y = np.arange(max(i - 2, 0), min(j + 2, len(prof)))
        w = np.clip(dark_thr + 40 - prof[y], 0, None)
        if w.sum() > 0:
            out.append((w * (y + 0.5)).sum() / w.sum())
    return np.array(out)


def _ruler_fit(arr):
    """Fit depth = (px - b) / a along the rows of arr (rows = depth). Raises ValueError."""
    rgb = arr.astype(np.int16)
    g = rgb.mean(2).astype(np.float32)
    chroma = (rgb.max(2) - rgb.min(2)).astype(np.float32)
    H, W = g.shape
    # ruler strip: neutral (grey/white) columns that are white along a large part of the length
    wf = ((g > 200) & (chroma < 30)).mean(0)
    nf = (chroma < 30).mean(0)
    runs = _runs((wf > 0.2) & (wf < 0.85) & (nf > 0.6))
    if not runs:
        raise ValueError("no ruler strip found")
    r0, r1 = max(runs, key=lambda r: r[1] - r[0])
    pad = int(0.08 * (r1 - r0))
    prof = np.median(g[:, r0 + pad:max(r1 - pad, r0 + pad + 1)], axis=1)
    # 10 cm blocks: white vs grey, thin tick lines removed with a running majority filter
    thr = 0.5 * (np.percentile(prof, 95) + np.percentile(prof, 25))
    k = max(3, H // 400) | 1
    white = np.convolve((prof > thr).astype(float), np.ones(k) / k, "same") > 0.5
    rl = _runs(white) + _runs(~white)
    long_runs = [e - s for s, e in rl if e - s > H / 40]
    if len(long_runs) < 4:
        raise ValueError("ruler blocks not found")
    L10 = float(np.median(long_runs))
    edges = []                                            # (px, +1 grey->white / -1 white->grey)
    for e in np.flatnonzero(np.diff(white.astype(np.int8))) + 1:
        before = e - np.flatnonzero(white[:e] != white[e - 1])[-1] - 1 if (white[:e] != white[e - 1]).any() else e
        after_idx = np.flatnonzero(white[e:] != white[e])
        after = after_idx[0] if len(after_idx) else H - e
        if before >= 0.3 * L10 and after >= 0.3 * L10:
            edges.append((float(e), 1 if white[e] else -1))
    ups = [e for e, t in edges if t > 0]
    if len(edges) < 4 or not ups:
        raise ValueError("ruler block edges not found")
    e10 = ups[0]                                          # first grey->white edge = 10 cm
    cm_e = np.array([10 + 10 * round((e - e10) / L10) for e, _ in edges])
    px_e = np.array([e for e, _ in edges])
    parity_ok = all((round(c / 10) % 2 == 1) == (t > 0) for c, (_, t) in zip(cm_e, edges))
    a, b = np.polyfit(cm_e, px_e, 1)
    res_blocks = np.abs((px_e - b) / a - cm_e)
    # refine on the cm / half-cm tick lines
    grey = np.median(prof[(prof < thr) & (prof > 0.5 * np.percentile(prof, 25))]) \
        if ((prof < thr) & (prof > 0.5 * np.percentile(prof, 25))).any() else thr / 2
    c = _tick_centres(prof, min(70.0, 0.6 * grey))
    method, res = f"{len(edges)} 10 cm ruler-block edges", res_blocks
    if len(c) >= 20:
        for tol in (0.08, 0.04):
            cc = (c - b) / a
            near = np.round(cc * 2) / 2
            keep = np.abs(cc - near) < tol
            if keep.sum() < 20 or np.ptp(near[keep]) < 30:
                break
            a, b = np.polyfit(near[keep], c[keep], 1)
        else:
            cc = (c - b) / a
            near = np.round(cc * 2) / 2
            keep = np.abs(cc - near) < 0.04
            method, res = f"{int(keep.sum())} ruler ticks (cm + half cm)", np.abs(cc - near)[keep]
    block_err_mm = 10 * np.abs((px_e - b) / a - cm_e)     # independent check of the final fit
    return dict(px_per_cm=a, px_at_0=b, fit_res_mm=10 * np.asarray(res), method=method,
                block_err_mm=block_err_mm,
                ruler=(r0, r1), parity_ok=parity_ok, edge_depth_cm=-b / a)


def _sediment_band(arr, cal):
    """Across-core pixel range of the sediment (coloured, not ruler, not black liner gaps)."""
    H = arr.shape[0]
    y0 = int(np.clip(cm_to_px(cal, 15), 0, H - 1)); y1 = int(np.clip(cm_to_px(cal, 85), y0 + 1, H))
    sl = arr[y0:y1:max(1, (y1 - y0) // 2000)].astype(np.int16)
    g = np.median(sl.mean(2), axis=0)
    chroma = np.median(sl.max(2) - sl.min(2), axis=0)
    r0, r1 = cal["ruler"]
    chroma[max(r0 - 2, 0):r1 + 2] = 0                       # never the ruler
    ok = (chroma >= 0.8 * np.percentile(chroma, 95)) & (g > 35) & (g < 220)
    runs = _runs(ok)
    if not runs:
        raise ValueError("sediment strip not found")
    s0, s1 = max(runs, key=lambda r: r[1] - r[0])
    return float(s0), float(s1)


def calibrate_photo(img, core):
    """Calibrate an original core photo (any size); returns the photo turned so depth runs down."""
    rgb = np.asarray(img.convert("RGB"))
    if rgb.shape[1] > rgb.shape[0]:
        rgb = rgb.transpose(1, 0, 2)                      # depth along rows
    tried = []
    for flip in (False, True):
        arr = np.ascontiguousarray(rgb[::-1] if flip else rgb)
        try:
            cal = _ruler_fit(arr)
        except ValueError as e:
            tried.append(str(e)); continue
        tried.append(f"photo end = {cal['edge_depth_cm']:.1f} cm")
        if cal["parity_ok"] and -2.0 <= cal["edge_depth_cm"] <= 1.5:   # ruler 0 cm at this end
            cal.update(core=core, arr=arr, source="photo", core_band=_sediment_band(arr, cal))
            return cal
    raise ValueError(f"GUAC-{core}: could not read the ruler on this photo ({'; '.join(tried)}). "
                     "Use the original core photo with the ruler starting at 0 cm"
                     + (", or GUAC-24A_realscale.pdf." if core == "24A" else "."))


def calibrate_24A_pdf(img):
    """GUAC-24A_realscale.pdf image: starts at 10.5 cm, ruler on the left; fit on all ticks."""
    arr = np.asarray(img.convert("RGB"))
    g = arr.mean(2)
    prof = np.median(g[:, int(10 / 819 * img.width):int(125 / 819 * img.width)], axis=1)
    c = _tick_centres(prof, 50)
    assert len(c) > 50, "24A ruler ticks not found - is this the GUAC-24A_realscale.pdf image?"
    scale = img.height / 8679
    ref = c[np.argmin(np.abs(c - REF24_PX_PER_CM * scale * (20 - REF24_TOP_CM)))]   # 20 cm tick
    cm = 20 + np.round((c - ref) / (REF24_PX_PER_CM * scale / 2)) * 0.5
    a, b = np.polyfit(cm, c, 1)
    assert abs(a / scale / REF24_PX_PER_CM - 1) < 0.005 and abs(-b / a - REF24_TOP_CM) < 0.05
    return dict(core="24A", arr=arr, source="pdf", px_per_cm=a, px_at_0=b,
                fit_res_mm=10 * np.abs((c - b) / a - cm), method=f"{len(c)} ruler ticks (cm + half cm)",
                core_band=(160 / 819 * img.width, 800 / 819 * img.width))


def cm_to_px(cal, cm):
    return cal["px_per_cm"] * np.asarray(cm, float) + cal["px_at_0"]


def px_to_cm(cal, px):
    return (np.asarray(px, float) - cal["px_at_0"]) / cal["px_per_cm"]


# ------------------------------------------------------------------ layout (inches; same for every figure)
FIG_H = 575 / 72
DEPTH_SPAN_CM = 9.0
BOTTOM, PLOT_H = 1.55, 5.45
CORE_W_CM = 3.0
CORE_W_IN = CORE_W_CM * PLOT_H / DEPTH_SPAN_CM      # true aspect ratio for the photo
L_CORE = 1.75
PANEL_W = 3.6
GAP = 0.45
L_PHI = L_CORE + CORE_W_IN + GAP
L_FRAC = L_PHI + PANEL_W + GAP
L_LEG = L_FRAC + PANEL_W + 0.3
FIG_W = L_LEG + 6.2
XLABEL_Y = -0.115
MARKER = 9
SHOW_SORTING = True   # False -> mean line only, without the sorting envelope

PHI_LIM = (9.0, 2.0)           # inverted: coarser to the right
PHI_TICKS = [9, 8, 7, 6, 5, 4, 3, 2]

CLASSES = [  # Folk-Ward classes, fine -> coarse (coarse sand, v. coarse sand, granule merged)
    ("Clay", ["Clay(>8)"], "> 8 φ", "#5B4636"),
    ("V. fine silt", ["V.fine silt(7-8)"], "7–8 φ", "#7E6650"),
    ("Fine silt", ["Fine silt(6-7)"], "6–7 φ", "#A08468"),
    ("Medium silt", ["Medium silt(5-6)"], "5–6 φ", "#C2A584"),
    ("Coarse silt", ["Coarse silt(4-5)"], "4–5 φ", "#DCC6A2"),
    ("V. fine sand", ["V.fine sand(3-4)"], "3–4 φ", "#F1E2B8"),
    ("Fine sand", ["Fine sand(2-3)"], "2–3 φ", "#F2C66D"),
    ("Medium sand", ["Medium sand(1-2)"], "1–2 φ", "#E39A3B"),
    ("Coarse sand +", ["Coarse sand(0-1)", "V.coarse sand(-1to0)", "Granule(<-1)"],
     "< 1 φ", "#A8521C"),
]

EVENT_STYLE = {"E2": ("Flood-triggered", "#E39A3B"), "E3": ("Earthquake-triggered", "#3E7CA6")}


def window_top(ev, cal):
    """Top (cm) of the DEPTH_SPAN_CM window: event centred, kept inside the photo."""
    if ev.base - ev.top > DEPTH_SPAN_CM:
        raise ValueError(f"event {ev.top}-{ev.base} cm is longer than DEPTH_SPAN_CM = {DEPTH_SPAN_CM}")
    photo_top = px_to_cm(cal, 0)
    photo_base = px_to_cm(cal, cal["arr"].shape[0])
    d0 = round((ev.top + ev.base) / 2 - DEPTH_SPAN_CM / 2, 1)
    d0 = max(d0, np.ceil(photo_top * 10) / 10)
    d0 = min(d0, np.floor((photo_base - DEPTH_SPAN_CM) * 10) / 10)
    return float(d0)
EVENT_TOL_CM = 0.15  # sample centre may sit this far outside the event top/base (slice width)
POINT_HALF_CM = 0.1  # 24A labels are single depths of 2 mm slices: drawn as depth +/- 1 mm


# ------------------------------------------------------------------ data
def parse_interval(label):
    nums = [float(n) for n in re.findall(r"\d+(?:\.\d+)?", str(label))]
    if len(nums) >= 2 and 0 < nums[1] - nums[0] <= 1.0:         # "15.5 - 16cm" (29A)
        return nums[0], nums[1]
    return nums[0] - POINT_HALF_CM, nums[0] + POINT_HALF_CM     # "13.0cm" (24A)


def load_gsd(path):
    d = pd.read_csv(path)
    d[["top", "base"]] = [parse_interval(s) for s in d["Sample"]]
    num = d.columns.drop("Sample")
    dup = d[d.duplicated(["top", "base"], keep=False)]
    if len(dup):  # replicate measurements at one depth -> mean of the replicates
        print(f"  {Path(path).name}: replicates averaged at {sorted(set(dup.Sample))}")
        d = d.groupby(["top", "base"], as_index=False).agg(
            {**{c: "mean" for c in num if c not in ("top", "base")}, "Sample": "first"})
        d["Mean_um"] = 1000 * 2.0 ** -d["Mean_phi"]
    d = d.sort_values("top").reset_index(drop=True)
    for name, cols, *_ in CLASSES:
        d[name] = d[cols].sum(axis=1)
    return d


def load_events(path, core):
    """Event depths in `core` and ages, from the sheet of the Excel file that has these columns."""
    need = EVENT_COLS + [f"{core}_top_cm", f"{core}_base_cm"]
    sheets = pd.read_excel(path, sheet_name=None)
    for name, x in sheets.items():
        x.columns = [str(c).strip() for c in x.columns]
        if set(need) <= set(x.columns):
            ev = pd.DataFrame({
                "event": x["event"].astype(str).str.strip(), "age": x["age_CE"], "unc": x["age_unc"],
                "top": pd.to_numeric(x[f"{core}_top_cm"], errors="coerce"),
                "base": pd.to_numeric(x[f"{core}_base_cm"], errors="coerce"),
            })
            return ev.dropna(subset=["event"]).set_index("event")
    found = {n: list(x.columns)[:12] for n, x in sheets.items()}
    raise KeyError(f"{Path(path).name}: no sheet has the columns {need}. Sheets/columns found: {found}")


def profile(df, values):
    """Continuous depth profile through the sample centres.

    Values are held flat from the top of the first sample to its centre and from the
    centre of the last sample to its base; between centres they are joined by straight
    lines, so depths with no sample (e.g. GUAC-24A 14.4-14.8 cm) are filled by linear
    interpolation between the neighbouring samples.
    """
    v = np.asarray(values, float)
    y = np.r_[df.top.iloc[0], (df.top + df.base).values / 2, df.base.iloc[-1]]
    return y, np.r_[v[0], v, v[-1]]


def sample_gaps(df):
    """(from_cm, to_cm) depth ranges between samples that were not measured."""
    return [(b, t) for b, t in zip(df.base.values[:-1], df.top.values[1:]) if t - b > 0.05]


# ------------------------------------------------------------------ plotting
def axes_at(fig, left, width):
    return fig.add_axes([left / FIG_W, BOTTOM / FIG_H, width / FIG_W, PLOT_H / FIG_H])


def core_crop(cal, d0, d1):
    """Central CORE_W_CM-wide strip between d0 and d1 cm, depth increasing downward."""
    arr = cal["arr"]                                  # photo turned so depth runs down the rows
    p0, p1 = int(np.floor(cm_to_px(cal, d0))), int(np.ceil(cm_to_px(cal, d1)))
    assert p0 >= 0 and p1 <= arr.shape[0], f"{cal['core']}: window {d0}-{d1} cm runs off the photo"
    lo, hi = cal["core_band"]
    mid, half = (lo + hi) / 2, CORE_W_CM * cal["px_per_cm"] / 2
    c0, c1 = int(round(mid - half)), int(round(mid + half))
    assert c0 >= lo - 1 and c1 <= hi + 1, f"{cal['core']}: sediment strip narrower than {CORE_W_CM} cm"
    strip = arr[p0:p1, c0:c1]
    extent = [0, (c1 - c0) / cal["px_per_cm"], float(px_to_cm(cal, p1)), float(px_to_cm(cal, p0))]
    return strip, extent


def depth_axis(ax, d0, d1):
    ax.set_ylim(d1, d0)
    ax.yaxis.set_major_locator(mticker.MultipleLocator(1))
    ax.yaxis.set_minor_locator(mticker.MultipleLocator(0.5))


def event_lines(ax, ev, color):
    for z in (ev.top, ev.base):
        ax.axhline(z, color=color, lw=2.5, ls="--", zorder=5, clip_on=False)


def tint(color, alpha=0.22):
    """Opaque colour that looks like `color` at `alpha` over white (no PDF transparency)."""
    from matplotlib.colors import to_rgb
    return tuple(1 - alpha * (1 - c) for c in to_rgb(color))


def unclip(ax):
    """Editable-PDF friendly: no clipping masks on the grid / tick lines of this axes."""
    for line in ax.get_xgridlines() + ax.get_ygridlines():
        line.set_clip_on(False)


def plot_event(core, code, gsd, events, cal):
    title, color = EVENT_STYLE[code]
    ev = events.loc[code]
    d0 = window_top(ev, cal)
    d1 = d0 + DEPTH_SPAN_CM
    mid = (gsd.top + gsd.base) / 2           # only samples belonging to this event
    sub = gsd[(mid >= ev.top - EVENT_TOL_CM) & (mid <= ev.base + EVENT_TOL_CM)]
    for g0, g1 in sample_gaps(sub):
        print(f"  {core} {code}: no samples {g0:.1f}-{g1:.1f} cm -> filled by interpolation")
    fig = plt.figure(figsize=(FIG_W, FIG_H))

    # core photo
    ax = axes_at(fig, L_CORE, CORE_W_IN)
    strip, ext = core_crop(cal, d0, d1)
    ax.imshow(strip, extent=ext, aspect="equal", interpolation="lanczos")
    ax.set_xlim(0, CORE_W_CM); depth_axis(ax, d0, d1)
    ax.set_xticks([])
    ax.set_ylabel("Depth (cm)")
    ax.set_xlabel(f"GUAC-{core}"); ax.xaxis.set_label_coords(0.5, XLABEL_Y)
    event_lines(ax, ev, color)

    # mean (phi): continuous line through the samples (style of the 24A reference figure)
    # Everything is drawn inside the axes limits without clipping masks, so each object stays
    # a plain editable path in the PDF (Canvas X Draw, Illustrator, Inkscape).
    a = axes_at(fig, L_PHI, PANEL_W)
    lo_x, hi_x = min(PHI_LIM), max(PHI_LIM)
    centre = (sub.top + sub.base) / 2
    if SHOW_SORTING:
        y, lo = profile(sub, sub.Mean_phi + sub.Sorting_phi)
        _, hi = profile(sub, sub.Mean_phi - sub.Sorting_phi)
        a.fill_betweenx(y, np.clip(lo, lo_x, hi_x), np.clip(hi, lo_x, hi_x),
                        color=tint(color), lw=0, clip_on=False, zorder=1)
    x = np.clip(sub.Mean_phi, lo_x, hi_x)
    a.plot(x, centre, "-", color=color, lw=3, zorder=6, clip_on=False)
    a.plot(x, centre, "o", ms=MARKER, mfc=color, mec=INK, mew=1.5, zorder=7, clip_on=False)
    a.set_xlim(*PHI_LIM); depth_axis(a, d0, d1)
    a.set_xticks(PHI_TICKS)
    a.tick_params(labelleft=False)
    a.set_axisbelow(False)
    a.grid(axis="x", color="#C8CDD3", lw=1, zorder=2)
    a.set_xlabel("Mean (φ)"); a.xaxis.set_label_coords(0.5, XLABEL_Y)
    event_lines(a, ev, color)
    unclip(a)

    # stacked class fractions, continuous (gaps between samples filled by interpolation)
    a = axes_at(fig, L_FRAC, PANEL_W)
    left = np.zeros(len(sub))
    for name, _, _, c in CLASSES:
        y, v0 = profile(sub, left)
        _, v1 = profile(sub, left + sub[name].values)
        a.fill_betweenx(y, np.clip(v0, 0, 100), np.clip(v1, 0, 100), color=c, lw=0, clip_on=False)
        left = left + sub[name].values
    a.set_xlim(0, 100); depth_axis(a, d0, d1)
    a.set_xticks([0, 50, 100]); a.set_xticks([25, 75], minor=True)
    a.tick_params(labelleft=False)
    a.set_xlabel("Fraction (%)"); a.xaxis.set_label_coords(0.5, XLABEL_Y)
    event_lines(a, ev, color)

    # legend: classes coarse -> fine, then envelope and event lines
    handles = [plt.Rectangle((0, 0), 1, 1, fc=c) for *_, c in CLASSES][::-1]
    labels = [f"{n}  ({p})" for n, _, p, _ in CLASSES][::-1]
    handles.append(plt.Line2D([], [], color=color, lw=3, marker="o", ms=MARKER,
                              mfc=color, mec=INK, mew=1.5))
    labels.append("Mean grain size")
    if SHOW_SORTING:
        handles.append(plt.Rectangle((0, 0), 1, 1, fc=tint(color)))
        labels.append("Mean ± 1σ sorting")
    handles.append(plt.Line2D([], [], color=color, lw=2.5, ls="--"))
    labels.append(f"{code} top / base")
    fig.legend(handles, labels, loc="upper left", frameon=False,
               bbox_to_anchor=(L_LEG / FIG_W, (BOTTOM + PLOT_H) / FIG_H),
               handlelength=1.2, handleheight=0.9, labelspacing=0.3, borderaxespad=0)

    fig.text(L_CORE / FIG_W, (FIG_H - 0.25) / FIG_H,
             f"{code} · {title} ({ev.age:.1f} ± {ev.unc:.1f} CE)",
             color=color, fontweight="bold", va="top", ha="left")

    OUT.mkdir(exist_ok=True)
    stem = OUT / f"GUAC{core}_{code}_core_grain_size_phi"
    for ext_ in ("png", "pdf", "svg"):      # dpi=300 also sets the core-photo resolution in PDF/SVG
        fig.savefig(f"{stem}.{ext_}", dpi=300)
    if IN_COLAB:
        plt.show()
    plt.close(fig)
    return stem, sub


# ------------------------------------------------------------------ checks
MAX_DEPTH_ERR_MM = 1.0   # run stops if any ruler mark is further than this from the fit


def run_checks(core, gsd, events, cal):
    """Stop with an AssertionError if the calibration or the data are not right."""
    # 1. depth calibration: fit residuals + independent check on the 10 cm block edges
    err_mm = np.r_[cal["fit_res_mm"], cal.get("block_err_mm", [])]
    what = cal["method"]
    assert err_mm.max() < MAX_DEPTH_ERR_MM, \
        f"{core}: ruler error too large: {err_mm.max():.2f} mm (limit {MAX_DEPTH_ERR_MM} mm)"
    # 2. crop maps back to the requested window (within 1 pixel)
    for code in EVENTS_TO_PLOT:
        assert code in events.index, f"{core}: event {code} not in the events file"
        ev = events.loc[code]
        assert np.isfinite(ev.top) and np.isfinite(ev.base), f"{core} {code}: no depths in the events file"
        assert ev.top < ev.base, f"{core} {code}: top {ev.top} >= base {ev.base}"
        d0 = window_top(ev, cal)
        _, ext = core_crop(cal, d0, d0 + DEPTH_SPAN_CM)
        px_cm = 1 / cal["px_per_cm"]
        assert abs(ext[3] - d0) < px_cm and abs(ext[2] - d0 - DEPTH_SPAN_CM) < px_cm
    # 3. grain-size data
    assert np.allclose(gsd[[c[0] for c in CLASSES]].sum(axis=1), 100, atol=0.01)
    assert np.allclose(1000 * 2.0 ** -gsd.Mean_phi, gsd.Mean_um, rtol=1e-6)
    # 5. continuous profiles: depth increasing, no holes, fractions still sum to 100 %
    mid = (gsd.top + gsd.base) / 2
    for code in EVENTS_TO_PLOT:
        ev = events.loc[code]
        sub = gsd[(mid >= ev.top - EVENT_TOL_CM) & (mid <= ev.base + EVENT_TOL_CM)]
        assert len(sub), f"{core} {code}: no grain-size samples between {ev.top} and {ev.base} cm"
        y, _ = profile(sub, sub.Mean_phi)
        assert np.all(np.diff(y) >= 0) and np.isfinite(y).all()
        grid = np.linspace(y[0], y[-1], 500)
        tot = sum(np.interp(grid, *profile(sub, sub[c[0]].values)) for c in CLASSES)
        assert np.allclose(tot, 100, atol=0.01), f"{core} {code}: interpolated total != 100 %"

    max_mm = err_mm.max()
    print(f"GUAC-{core}: ALL CHECKS PASSED")
    print(f"  core image      : {cal['source']}, {cal['arr'].shape[1]} x {cal['arr'].shape[0]} px")
    print(f"  depth scale     : {cal['px_per_cm']:.3f} px/cm (1 px = {10 / cal['px_per_cm']:.2f} mm),"
          f" fitted on {what}")
    if len(cal.get("block_err_mm", [])):
        print(f"  check vs blocks : max {cal['block_err_mm'].max():.2f} mm on "
              f"{len(cal['block_err_mm'])} independent 10 cm block edges")
    print(f"  depth error     : max {max_mm:.2f} mm, sd {err_mm.std():.2f} mm")
    print(f"  depth accuracy  : {100 * (1 - max_mm / 10 / DEPTH_SPAN_CM):.1f} % of the 9 cm window"
          f" (worst case)")
    print(f"  grain-size data : {len(gsd)} samples, 100 % as in the CSV (no smoothing; "
          f"fractions sum to 100 %)")


def main():
    paths = find_inputs()
    for core in ("29A", "24A"):
        print()
        gsd, events = load_gsd(paths[f"csv_{core}"]), load_events(paths["events"], core)
        img_path = paths[f"img_{core}"]
        if img_path.suffix.lower() == ".pdf":
            cal = calibrate_24A_pdf(image_from_pdf(img_path))
        else:
            cal = calibrate_photo(Image.open(img_path), core)
        run_checks(core, gsd, events, cal)
        for code in EVENTS_TO_PLOT:
            stem, sub = plot_event(core, code, gsd, events, cal)
            d0 = window_top(events.loc[code], cal)
            print(f"  {code}: {d0:.1f}-{d0 + DEPTH_SPAN_CM:.1f} cm, {len(sub)} samples, mean "
                  f"{sub.Mean_phi.max():.2f}-{sub.Mean_phi.min():.2f} φ -> {stem.name}.png/.pdf/.svg")
    if IN_COLAB:
        import shutil
        from google.colab import files
        shutil.make_archive(str(WORK / "GUAC_event_figures"), "zip", OUT)
        files.download(str(WORK / "GUAC_event_figures.zip"))


if __name__ == "__main__":
    main()
