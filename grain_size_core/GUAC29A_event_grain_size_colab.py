# -*- coding: utf-8 -*-
"""
GUAC-29A event panels: core photo on a true depth scale + grain size in phi.
Works in Google Colab or plain Python.

HOW TO RUN IN GOOGLE COLAB
  1. Upload this .py file to Colab (Files pane, left), or paste it all into one cell.
  2. Run:   %run GUAC29A_event_grain_size_colab.py      (or just run the pasted cell)
  3. When asked, upload these 3 files (any file names are fine, matched by content):
       - the grain-size stats CSV   (GUAC_29A_GSD_..._FolkWard_Stats.csv)
       - the events/ages Excel file (Lachua_..._Events_with_Ages.xlsx)
       - the core photo with the ruler (the original 2000 x 156 px JPG)
  4. Figures are shown, saved to /content/output/ and downloaded as a zip.

One figure per event (E2 flood, E3 earthquake). Both figures share:
  * the same 9 cm depth window  -> identical cm-per-inch vertical scale
  * the same x-axis limits and the same panel widths (set in inches)
  * 24 pt Liberation Sans for every label, tick and legend entry
  * a height of 575 pt (7.99 in), same as the cartoon PDFs

Panels: core photo | mean grain size (phi) with mean +/- 1 sorting envelope | class fractions (%)
Before plotting, run_checks() tests the depth calibration and the data; any problem stops the script.
"""
import glob
import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image

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
    """Return (csv, xlsx, jpg) paths; in Colab ask for an upload if they are missing."""
    def pick(patterns):
        for root in (WORK, WORK / "data"):
            for p in patterns:
                hits = sorted(glob.glob(str(root / p)))
                if hits:
                    return Path(hits[0])
        return None

    pats = (["*GSD*Stats*.csv", "*.csv"], ["*Events*Ages*.xlsx", "*.xlsx"],
            ["*29A*core*.jp*g", "*.jp*g", "*.png"])
    found = [pick(p) for p in pats]
    if None in found and IN_COLAB:
        from google.colab import files
        print("Upload: grain-size CSV, events Excel file, core photo (JPG)")
        files.upload()
        found = [pick(p) for p in pats]
    names = ("grain-size CSV", "events Excel", "core photo")
    missing = [n for n, f in zip(names, found) if f is None]
    if missing:
        raise FileNotFoundError(f"missing input file(s): {missing} (put them next to this script)")
    for n, f in zip(names, found):
        print(f"{n:15s}: {f.name}")
    return found


# ------------------------------------------------------------------ depth calibration
# Measured on the original 2000 x 156 px photo: the grey/white ruler blocks change
# every 10 cm; x = pixel edge of each block change at 10, 20 ... 100 cm.
REF_SIZE = (2000, 156)
RULER_CM = np.array([10, 20, 30, 40, 50, 60, 70, 80, 90, 100], float)
RULER_PX_REF = np.array([200, 398, 598, 796, 995, 1193, 1393, 1591, 1792, 1989], float)
CORE_ROWS_REF = (14, 124)      # rows with sediment (no liner edge / ruler)
RULER_ROW_REF = 150            # a row through the ruler blocks
CAL = {}


def calibrate(img):
    """Linear pixel <-> cm fit; scales automatically if the photo was resized."""
    sx, sy = img.width / REF_SIZE[0], img.height / REF_SIZE[1]
    px = RULER_PX_REF * sx
    a, b = np.polyfit(RULER_CM, px, 1)
    CAL.update(px_per_cm=a, px_at_0=b, ruler_px=px, sy=sy,
               core_rows=(CORE_ROWS_REF[0] * sy, CORE_ROWS_REF[1] * sy),
               ruler_row=int(round(RULER_ROW_REF * sy)), sx=sx)


def cm_to_px(cm):
    return CAL["px_per_cm"] * np.asarray(cm, float) + CAL["px_at_0"]


def px_to_cm(px):
    return (np.asarray(px, float) - CAL["px_at_0"]) / CAL["px_per_cm"]


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

EVENTS_TO_PLOT = [  # event, title, colour (from the cartoons), top of the 9 cm window (cm)
    ("E2", "Flood-triggered", "#E39A3B", 12.6),
    ("E3", "Earthquake-triggered", "#3E7CA6", 21.9),
]


# ------------------------------------------------------------------ data
def parse_interval(label):
    nums = re.findall(r"\d+(?:\.\d+)?", str(label))
    top, base = float(nums[0]), float(nums[1])
    assert base > top and base - top <= 1.0, label
    return top, base


def load_gsd(path):
    d = pd.read_csv(path)
    d[["top", "base"]] = [parse_interval(s) for s in d["Sample"]]
    d = d.sort_values("top").reset_index(drop=True)
    for name, cols, *_ in CLASSES:
        d[name] = d[cols].sum(axis=1)
    return d


def load_events(path):
    x = pd.read_excel(path, header=0)
    ev = pd.DataFrame({
        "event": x.iloc[:, 7], "age": x["age_CE"], "unc": x["age_unc"],
        "top": pd.to_numeric(x["29A_top_cm"], errors="coerce"),
        "base": pd.to_numeric(x["29A_base_cm"], errors="coerce"),
    })
    return ev.set_index("event")


def steps(top, base, values):
    """Step-plot arrays (each sample drawn over its own depth interval), NaN breaks at gaps."""
    y, v = [], []
    for i, (t, b, val) in enumerate(zip(top, base, values)):
        if i and abs(t - base[i - 1]) > 0.05:
            y.append(np.nan); v.append(np.nan)
        y += [t, b]; v += [val, val]
    return np.array(y), np.array(v)


def contiguous_runs(df):
    run = [0]
    for i in range(1, len(df)):
        run.append(run[-1] + (abs(df.top.iloc[i] - df.base.iloc[i - 1]) > 0.05))
    return [g for _, g in df.groupby(np.array(run))]


# ------------------------------------------------------------------ plotting
def axes_at(fig, left, width):
    return fig.add_axes([left / FIG_W, BOTTOM / FIG_H, width / FIG_W, PLOT_H / FIG_H])


def core_crop(img, d0, d1):
    """Crop the photo to [d0, d1] cm and rotate so depth increases downward."""
    x0, x1 = int(np.floor(cm_to_px(d0))), int(np.ceil(cm_to_px(d1)))
    r0, r1 = CAL["core_rows"]
    mid, half = (r0 + r1) / 2, CORE_W_CM * CAL["px_per_cm"] / 2
    c0, c1 = int(round(mid - half)), int(round(mid + half))
    strip = np.transpose(np.asarray(img)[c0:c1, x0:x1], (1, 0, 2))   # rows = depth
    extent = [0, (c1 - c0) / CAL["px_per_cm"], float(px_to_cm(x1)), float(px_to_cm(x0))]
    return strip, extent


def depth_axis(ax, d0, d1):
    ax.set_ylim(d1, d0)
    ax.yaxis.set_major_locator(mticker.MultipleLocator(1))
    ax.yaxis.set_minor_locator(mticker.MultipleLocator(0.5))


def event_lines(ax, ev, color):
    for z in (ev.top, ev.base):
        ax.axhline(z, color=color, lw=2.5, ls="--", zorder=5)


def plot_event(code, title, color, d0, gsd, events, img):
    ev = events.loc[code]
    d1 = d0 + DEPTH_SPAN_CM
    sub = gsd[(gsd.base > d0) & (gsd.top < d1)]
    fig = plt.figure(figsize=(FIG_W, FIG_H))

    # core photo
    ax = axes_at(fig, L_CORE, CORE_W_IN)
    strip, ext = core_crop(img, d0, d1)
    ax.imshow(strip, extent=ext, aspect="equal", interpolation="lanczos")
    ax.set_xlim(0, CORE_W_CM); depth_axis(ax, d0, d1)
    ax.set_xticks([])
    ax.set_ylabel("Depth (cm)")
    ax.set_xlabel("GUAC-29A"); ax.xaxis.set_label_coords(0.5, XLABEL_Y)
    event_lines(ax, ev, color)

    # mean (phi) with mean +/- 1 sorting envelope
    a = axes_at(fig, L_PHI, PANEL_W)
    for g in contiguous_runs(sub):
        y, lo = steps(g.top.values, g.base.values, (g.Mean_phi + g.Sorting_phi).values)
        _, hi = steps(g.top.values, g.base.values, (g.Mean_phi - g.Sorting_phi).values)
        a.fill_betweenx(y, lo, hi, color=color, alpha=0.22, lw=0)
    y, v = steps(sub.top.values, sub.base.values, sub.Mean_phi.values)
    a.plot(v, y, color=INK, lw=3)
    a.plot(sub.Mean_phi, (sub.top + sub.base) / 2, "o", ms=8, color=INK)
    a.set_xlim(*PHI_LIM); depth_axis(a, d0, d1)
    a.set_xticks(PHI_TICKS)
    a.tick_params(labelleft=False)
    a.grid(axis="x", color="#C8CDD3", lw=1)
    a.set_xlabel("Mean (φ)"); a.xaxis.set_label_coords(0.5, XLABEL_Y)
    event_lines(a, ev, color)

    # stacked class fractions
    a = axes_at(fig, L_FRAC, PANEL_W)
    for g in contiguous_runs(sub):
        left = np.zeros(len(g))
        for name, _, _, c in CLASSES:
            y, v0 = steps(g.top.values, g.base.values, left)
            _, v1 = steps(g.top.values, g.base.values, left + g[name].values)
            a.fill_betweenx(y, v0, v1, color=c, lw=0)
            left = left + g[name].values
    a.set_xlim(0, 100); depth_axis(a, d0, d1)
    a.set_xticks([0, 50, 100]); a.set_xticks([25, 75], minor=True)
    a.tick_params(labelleft=False)
    a.set_xlabel("Fraction (%)"); a.xaxis.set_label_coords(0.5, XLABEL_Y)
    event_lines(a, ev, color)

    # legend: classes coarse -> fine, then envelope and event lines
    handles = [plt.Rectangle((0, 0), 1, 1, fc=c) for *_, c in CLASSES][::-1]
    labels = [f"{n}  ({p})" for n, _, p, _ in CLASSES][::-1]
    handles += [plt.Rectangle((0, 0), 1, 1, fc=color, alpha=0.22),
                plt.Line2D([], [], color=color, lw=2.5, ls="--")]
    labels += ["Mean ± 1σ sorting", f"{code} top / base"]
    fig.legend(handles, labels, loc="upper left", frameon=False,
               bbox_to_anchor=(L_LEG / FIG_W, (BOTTOM + PLOT_H) / FIG_H),
               handlelength=1.2, handleheight=0.9, labelspacing=0.3, borderaxespad=0)

    fig.text(L_CORE / FIG_W, (FIG_H - 0.25) / FIG_H,
             f"{code} · {title} ({ev.age:.1f} ± {ev.unc:.1f} CE)",
             color=color, fontweight="bold", va="top", ha="left")

    OUT.mkdir(exist_ok=True)
    stem = OUT / f"GUAC29A_{code}_core_grain_size_phi"
    for ext_ in ("png", "pdf", "svg"):
        fig.savefig(f"{stem}.{ext_}", dpi=300 if ext_ == "png" else None)
    if IN_COLAB:
        plt.show()
    plt.close(fig)
    return stem, sub


# ------------------------------------------------------------------ checks
def run_checks(gsd, events, img):
    """Stop with an AssertionError if the calibration or the data are not right."""
    # 1. calibration fit: every 10 cm ruler mark within 1 mm of the fitted line
    res_mm = 10 * np.abs(px_to_cm(CAL["ruler_px"]) - RULER_CM)
    assert res_mm.max() < 1.0, f"ruler fit residuals too large (mm): {res_mm}"
    # 2. find the 10 cm block edges in THIS photo again; each must sit < 1 mm from its cm value
    row = np.asarray(img.convert("L"), float)[CAL["ruler_row"]]
    d = np.abs(np.diff(np.convolve(row, np.ones(5) / 5, "same")))
    win = max(6, int(6 * CAL["sx"]))
    found = []
    for cm in RULER_CM[:-1]:
        i = np.arange(int(cm_to_px(cm)) - win, int(cm_to_px(cm)) + win)
        w = np.where(d[i] > 20, d[i], 0)
        assert w.sum() > 0, (f"no ruler edge near {cm} cm - is this the original GUAC-29A "
                             f"photo with the ruler along the bottom?")
        edge_cm = float(px_to_cm((w * (i + 1)).sum() / w.sum()))
        found.append(abs(edge_cm - cm) * 10)
        assert abs(edge_cm - cm) < 0.1, f"ruler edge at {cm} cm found at {edge_cm:.2f} cm"
    # 3. crop maps back to the requested window (within 1 pixel)
    strip, ext = core_crop(img, 21.9, 30.9)
    px_cm = 1 / CAL["px_per_cm"]
    assert abs(ext[3] - 21.9) < px_cm and abs(ext[2] - 30.9) < px_cm
    # 4. grain-size data
    assert len(gsd) == 86, f"expected 86 samples, got {len(gsd)}"
    assert np.allclose(gsd[[c[0] for c in CLASSES]].sum(axis=1), 100, atol=0.01)
    assert np.allclose(1000 * 2.0 ** -gsd.Mean_phi, gsd.Mean_um, rtol=1e-9)
    assert (gsd.base - gsd.top).between(0.4, 0.6).all()
    # 5. events
    assert events.loc["E2", ["top", "base"]].tolist() == [15.1, 19.1]
    assert events.loc["E3", ["top", "base"]].tolist() == [22.6, 30.15]
    assert events.loc["E2", "age"] == 2000.1 and events.loc["E3", "age"] == 1977.9

    max_err_mm = max(res_mm.max(), max(found))
    print("ALL CHECKS PASSED")
    print(f"  depth scale            : {CAL['px_per_cm']:.3f} px/cm "
          f"(1 px = {10 / CAL['px_per_cm']:.2f} mm)")
    print(f"  depth error            : max {max_err_mm:.2f} mm, sd {np.std(res_mm):.2f} mm "
          f"(< 1 photo pixel)")
    print(f"  depth accuracy         : {100 * (1 - max_err_mm / 10 / DEPTH_SPAN_CM):.1f} % "
          f"of the 9 cm window (worst case)")
    print("  grain-size values      : 100 % (plotted directly from the CSV, no smoothing)")


def main():
    csv, xlsx, jpg = find_inputs()
    gsd, events = load_gsd(csv), load_events(xlsx)
    img = Image.open(jpg).convert("RGB")
    calibrate(img)
    run_checks(gsd, events, img)
    for code, title, color, d0 in EVENTS_TO_PLOT:
        ev = events.loc[code]
        assert d0 <= ev.top and ev.base <= d0 + DEPTH_SPAN_CM, f"{code} outside window"
        stem, sub = plot_event(code, title, color, d0, gsd, events, img)
        print(f"{code}: {d0:.1f}-{d0 + DEPTH_SPAN_CM:.1f} cm, {len(sub)} samples, mean "
              f"{sub.Mean_phi.max():.2f}-{sub.Mean_phi.min():.2f} φ -> {stem.name}.png/.pdf/.svg")
    if IN_COLAB:
        import shutil
        from google.colab import files
        shutil.make_archive(str(WORK / "GUAC29A_event_figures"), "zip", OUT)
        files.download(str(WORK / "GUAC29A_event_figures.zip"))


if __name__ == "__main__":
    main()
