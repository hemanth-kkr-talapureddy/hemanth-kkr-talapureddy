# -*- coding: utf-8 -*-
"""
GUAC-24A and GUAC-29A event panels: core photo on a true depth scale + grain size in phi.
Works in Google Colab or plain Python.

HOW TO RUN IN GOOGLE COLAB
  1. Upload this .py file to Colab (Files pane, left), or paste it all into one cell.
  2. Run:   %run GUAC_event_grain_size_colab.py      (or just run the pasted cell)
  3. When asked, upload these 5 files (names are matched loosely, see find_inputs):
       - GUAC_29A_GSD_..._FolkWard_Stats.csv        (29A grain-size stats)
       - GUAC_24A_GSD_..._FolkWard_Stats.csv        (24A grain-size stats)
       - Lachua_..._Events_with_Ages.xlsx           (events, depths in both cores, ages)
       - the original 29A core photo with the ruler (2000 x 156 px JPG)
       - GUAC-24A_realscale.pdf                     (24A core image is taken from this PDF)
  4. Four figures are shown, saved to /content/output/ and downloaded as a zip.

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
def find_inputs():
    """Return a dict of input paths; in Colab ask for an upload if any are missing."""
    def pick(patterns, test=None):
        for root in (WORK, WORK / "data"):
            for p in patterns:
                for h in sorted(glob.glob(str(root / p))):
                    if test is None or test(Path(h)):
                        return Path(h)
        return None

    def is_29a_photo(p):
        with Image.open(p) as im:
            return im.width > 5 * im.height          # horizontal core photo

    def is_24a_pdf(p):
        import pymupdf
        with pymupdf.open(p) as d:
            return any(h > 5000 for _, _, _, h, *_ in d[0].get_images(full=True))

    def locate():
        return {
            "csv_29A": pick(["*29A*GSD*.csv", "*29A*.csv"]),
            "csv_24A": pick(["*24A*GSD*.csv", "*24A*.csv"]),
            "xlsx": pick(["*Events*Ages*.xlsx", "*.xlsx"]),
            "img_29A": pick(["*29A*.jp*g", "*.jp*g", "*.png"], is_29a_photo),
            "img_24A": pick(["*24A*realscale*.pdf", "*24A*.pdf", "*.pdf"], is_24a_pdf),
        }

    ensure_pymupdf()
    found = locate()
    if None in found.values() and IN_COLAB:
        from google.colab import files
        print("Upload: 29A CSV, 24A CSV, events Excel, 29A core photo (JPG), GUAC-24A_realscale.pdf")
        files.upload()
        found = locate()
    missing = [k for k, v in found.items() if v is None]
    if missing:
        raise FileNotFoundError(f"missing input file(s): {missing} (put them next to this script)")
    for k, v in found.items():
        print(f"{k:8s}: {v.name}")
    return found


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
# GUAC-29A: horizontal 2000 x 156 px photo. The grey/white ruler blocks change every 10 cm;
# x = pixel edge of each change at 10, 20 ... 100 cm (measured once, re-checked every run).
RULER29_CM = np.array([10, 20, 30, 40, 50, 60, 70, 80, 90, 100], float)
RULER29_PX = np.array([200, 398, 598, 796, 995, 1193, 1393, 1591, 1792, 1989], float)
REF29_SIZE = (2000, 156)

# GUAC-24A: vertical 819 x 8679 px image from GUAC-24A_realscale.pdf, ruler on the left.
# Calibrated on every run from all cm and half-cm ruler ticks; expected values:
REF24_PX_PER_CM, REF24_TOP_CM = 104.209, 10.512


def calibrate_29A(img):
    sx, sy = img.width / REF29_SIZE[0], img.height / REF29_SIZE[1]
    px = RULER29_PX * sx
    a, b = np.polyfit(RULER29_CM, px, 1)
    return dict(core="29A", horizontal=True, px_per_cm=a, px_at_0=b,
                fit_res_mm=10 * np.abs((px - b) / a - RULER29_CM),
                core_band=(14 * sy, 124 * sy),          # image rows with sediment
                ruler_line=int(round(150 * sy)), sx=sx)


def ruler_ticks_24A(img):
    """Centre (px) of every dark ruler tick line along the 24A ruler strip."""
    g = np.asarray(img.convert("L"), float)
    prof = np.median(g[:, int(10 / 819 * img.width):int(125 / 819 * img.width)], axis=1)
    dark, centres, i = prof < 50, [], 0
    while i < len(prof):
        if dark[i]:
            j = i
            while j < len(prof) and dark[j]:
                j += 1
            if j - i >= 3:
                y = np.arange(max(i - 2, 0), min(j + 2, len(prof)))
                w = np.clip(120 - prof[y], 0, None)
                centres.append((w * (y + 0.5)).sum() / w.sum())
            i = j
        else:
            i += 1
    return np.array(centres)


def calibrate_24A(img):
    c = ruler_ticks_24A(img)
    assert len(c) > 50, "24A ruler ticks not found - is this the GUAC-24A_realscale.pdf image?"
    guess = lambda cm: REF24_PX_PER_CM * (cm - REF24_TOP_CM) * img.height / 8679
    ref = c[np.argmin(np.abs(c - guess(20)))]                    # 20 cm tick (white->grey block)
    step = REF24_PX_PER_CM * img.height / 8679 / 2
    cm = 20 + np.round((c - ref) / step) * 0.5
    a, b = np.polyfit(cm, c, 1)
    return dict(core="24A", horizontal=False, px_per_cm=a, px_at_0=b,
                fit_res_mm=10 * np.abs((c - b) / a - cm), n_ticks=len(c),
                core_band=(160 / 819 * img.width, 800 / 819 * img.width))  # columns with sediment


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
WINDOWS = {  # top of the 9 cm window (cm) per core and event; each event sits inside its window
    ("29A", "E2"): 12.6, ("29A", "E3"): 21.9,
    ("24A", "E2"): 10.6, ("24A", "E3"): 18.3,
}
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
    x = pd.read_excel(path, header=0)
    ev = pd.DataFrame({
        "event": x["event"], "age": x["age_CE"], "unc": x["age_unc"],
        "top": pd.to_numeric(x[f"{core}_top_cm"], errors="coerce"),
        "base": pd.to_numeric(x[f"{core}_base_cm"], errors="coerce"),
    })
    return ev.set_index("event")


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


def core_crop(img, cal, d0, d1):
    """Central CORE_W_CM-wide strip between d0 and d1 cm, depth increasing downward."""
    p0, p1 = int(np.floor(cm_to_px(cal, d0))), int(np.ceil(cm_to_px(cal, d1)))
    assert p0 >= 0 and p1 <= (img.width if cal["horizontal"] else img.height), \
        f"{cal['core']}: window {d0}-{d1} cm runs off the photo"
    lo, hi = cal["core_band"]
    mid, half = (lo + hi) / 2, CORE_W_CM * cal["px_per_cm"] / 2
    c0, c1 = int(round(mid - half)), int(round(mid + half))
    arr = np.asarray(img)
    strip = np.transpose(arr[c0:c1, p0:p1], (1, 0, 2)) if cal["horizontal"] else arr[p0:p1, c0:c1]
    extent = [0, (c1 - c0) / cal["px_per_cm"], float(px_to_cm(cal, p1)), float(px_to_cm(cal, p0))]
    return strip, extent


def depth_axis(ax, d0, d1):
    ax.set_ylim(d1, d0)
    ax.yaxis.set_major_locator(mticker.MultipleLocator(1))
    ax.yaxis.set_minor_locator(mticker.MultipleLocator(0.5))


def event_lines(ax, ev, color):
    for z in (ev.top, ev.base):
        ax.axhline(z, color=color, lw=2.5, ls="--", zorder=5)


def plot_event(core, code, gsd, events, img, cal):
    title, color = EVENT_STYLE[code]
    ev = events.loc[code]
    d0 = WINDOWS[(core, code)]
    d1 = d0 + DEPTH_SPAN_CM
    mid = (gsd.top + gsd.base) / 2           # only samples belonging to this event
    sub = gsd[(mid >= ev.top - EVENT_TOL_CM) & (mid <= ev.base + EVENT_TOL_CM)]
    for g0, g1 in sample_gaps(sub):
        print(f"  {core} {code}: no samples {g0:.1f}-{g1:.1f} cm -> filled by interpolation")
    fig = plt.figure(figsize=(FIG_W, FIG_H))

    # core photo
    ax = axes_at(fig, L_CORE, CORE_W_IN)
    strip, ext = core_crop(img, cal, d0, d1)
    ax.imshow(strip, extent=ext, aspect="equal", interpolation="lanczos")
    ax.set_xlim(0, CORE_W_CM); depth_axis(ax, d0, d1)
    ax.set_xticks([])
    ax.set_ylabel("Depth (cm)")
    ax.set_xlabel(f"GUAC-{core}"); ax.xaxis.set_label_coords(0.5, XLABEL_Y)
    event_lines(ax, ev, color)

    # mean (phi): continuous line through the samples (style of the 24A reference figure)
    a = axes_at(fig, L_PHI, PANEL_W)
    centre = (sub.top + sub.base) / 2
    if SHOW_SORTING:
        y, lo = profile(sub, sub.Mean_phi + sub.Sorting_phi)
        _, hi = profile(sub, sub.Mean_phi - sub.Sorting_phi)
        a.fill_betweenx(y, lo, hi, color=color, alpha=0.22, lw=0)
    a.plot(sub.Mean_phi, centre, "-", color=color, lw=3, zorder=6)
    a.plot(sub.Mean_phi, centre, "o", ms=MARKER, mfc=color, mec=INK, mew=1.5, zorder=7)
    a.set_xlim(*PHI_LIM); depth_axis(a, d0, d1)
    a.set_xticks(PHI_TICKS)
    a.tick_params(labelleft=False)
    a.grid(axis="x", color="#C8CDD3", lw=1)
    a.set_xlabel("Mean (φ)"); a.xaxis.set_label_coords(0.5, XLABEL_Y)
    event_lines(a, ev, color)

    # stacked class fractions, continuous (gaps between samples filled by interpolation)
    a = axes_at(fig, L_FRAC, PANEL_W)
    left = np.zeros(len(sub))
    for name, _, _, c in CLASSES:
        y, v0 = profile(sub, left)
        _, v1 = profile(sub, left + sub[name].values)
        a.fill_betweenx(y, v0, v1, color=c, lw=0)
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
        handles.append(plt.Rectangle((0, 0), 1, 1, fc=color, alpha=0.22))
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
    for ext_ in ("png", "pdf", "svg"):
        fig.savefig(f"{stem}.{ext_}", dpi=300 if ext_ == "png" else None)
    if IN_COLAB:
        plt.show()
    plt.close(fig)
    return stem, sub


# ------------------------------------------------------------------ checks
EXPECTED = {  # from the events spreadsheet; the run stops if the file disagrees
    "29A": {"n": 86, "E2": (15.1, 19.1), "E3": (22.6, 30.15)},
    "24A": {"n": 231, "E2": (13.0, 17.2), "E3": (19.6, 26.0)},
}


def check_depth_29A(img, cal):
    """Find the 10 cm block edges in the photo again; return their errors (mm)."""
    row = np.asarray(img.convert("L"), float)[cal["ruler_line"]]
    d = np.abs(np.diff(np.convolve(row, np.ones(5) / 5, "same")))
    win, err = max(6, int(6 * cal["sx"])), []
    for cm in RULER29_CM[:-1]:
        i = np.arange(int(cm_to_px(cal, cm)) - win, int(cm_to_px(cal, cm)) + win)
        w = np.where(d[i] > 20, d[i], 0)
        assert w.sum() > 0, (f"29A: no ruler edge near {cm} cm - is this the original GUAC-29A "
                             f"photo with the ruler along the bottom?")
        err.append(10 * abs(float(px_to_cm(cal, (w * (i + 1)).sum() / w.sum())) - cm))
    return np.array(err)


def run_checks(core, gsd, events, img, cal):
    """Stop with an AssertionError if the calibration or the data are not right."""
    exp = EXPECTED[core]
    # 1. depth calibration
    if core == "29A":
        err_mm = np.r_[cal["fit_res_mm"], check_depth_29A(img, cal)]
        what = "10 cm ruler blocks"
    else:
        err_mm = cal["fit_res_mm"]
        what = f"{cal['n_ticks']} ruler ticks (cm + half cm)"
        scale = img.height / 8679
        assert abs(cal["px_per_cm"] / scale / REF24_PX_PER_CM - 1) < 0.005, cal["px_per_cm"]
        assert abs(-cal["px_at_0"] / cal["px_per_cm"] - REF24_TOP_CM) < 0.05
    assert err_mm.max() < 1.0, f"{core}: ruler error too large (mm): {err_mm.max():.2f}"
    # 2. crop maps back to the requested window (within 1 pixel)
    for code in ("E2", "E3"):
        d0 = WINDOWS[(core, code)]
        _, ext = core_crop(img, cal, d0, d0 + DEPTH_SPAN_CM)
        px_cm = 1 / cal["px_per_cm"]
        assert abs(ext[3] - d0) < px_cm and abs(ext[2] - d0 - DEPTH_SPAN_CM) < px_cm
    # 3. grain-size data
    assert len(gsd) == exp["n"], f"{core}: expected {exp['n']} samples, got {len(gsd)}"
    assert np.allclose(gsd[[c[0] for c in CLASSES]].sum(axis=1), 100, atol=0.01)
    assert np.allclose(1000 * 2.0 ** -gsd.Mean_phi, gsd.Mean_um, rtol=1e-6)
    # 4. events
    for code in ("E2", "E3"):
        ev = events.loc[code]
        assert (ev.top, ev.base) == exp[code], f"{core} {code}: {ev.top}-{ev.base} cm"
        d0 = WINDOWS[(core, code)]
        assert d0 <= ev.top and ev.base <= d0 + DEPTH_SPAN_CM, f"{core} {code} outside window"
    # 5. continuous profiles: depth increasing, no holes, fractions still sum to 100 %
    mid = (gsd.top + gsd.base) / 2
    for code in ("E2", "E3"):
        ev = events.loc[code]
        sub = gsd[(mid >= ev.top - EVENT_TOL_CM) & (mid <= ev.base + EVENT_TOL_CM)]
        y, _ = profile(sub, sub.Mean_phi)
        assert np.all(np.diff(y) >= 0) and np.isfinite(y).all()
        grid = np.linspace(y[0], y[-1], 500)
        tot = sum(np.interp(grid, *profile(sub, sub[c[0]].values)) for c in CLASSES)
        assert np.allclose(tot, 100, atol=0.01), f"{core} {code}: interpolated total != 100 %"
    assert events.loc["E2", "age"] == 2000.1 and events.loc["E3", "age"] == 1977.9

    max_mm = err_mm.max()
    print(f"GUAC-{core}: ALL CHECKS PASSED")
    print(f"  depth scale     : {cal['px_per_cm']:.3f} px/cm (1 px = {10 / cal['px_per_cm']:.2f} mm),"
          f" checked on {what}")
    print(f"  depth error     : max {max_mm:.2f} mm, sd {err_mm.std():.2f} mm")
    print(f"  depth accuracy  : {100 * (1 - max_mm / 10 / DEPTH_SPAN_CM):.1f} % of the 9 cm window"
          f" (worst case)")
    print("  grain-size data : 100 % (plotted directly from the CSV, no smoothing)")


def main():
    paths = find_inputs()
    cores = {
        "29A": (paths["csv_29A"], Image.open(paths["img_29A"]).convert("RGB"), calibrate_29A),
        "24A": (paths["csv_24A"], image_from_pdf(paths["img_24A"]), calibrate_24A),
    }
    for core, (csv, img, calibrate) in cores.items():
        print()
        gsd, events = load_gsd(csv), load_events(paths["xlsx"], core)
        cal = calibrate(img)
        run_checks(core, gsd, events, img, cal)
        for code in ("E2", "E3"):
            stem, sub = plot_event(core, code, gsd, events, img, cal)
            d0 = WINDOWS[(core, code)]
            print(f"  {code}: {d0:.1f}-{d0 + DEPTH_SPAN_CM:.1f} cm, {len(sub)} samples, mean "
                  f"{sub.Mean_phi.max():.2f}-{sub.Mean_phi.min():.2f} φ -> {stem.name}.png/.pdf/.svg")
    if IN_COLAB:
        import shutil
        from google.colab import files
        shutil.make_archive(str(WORK / "GUAC_event_figures"), "zip", OUT)
        files.download(str(WORK / "GUAC_event_figures.zip"))


if __name__ == "__main__":
    main()
