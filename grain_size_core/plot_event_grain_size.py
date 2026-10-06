"""
Core crop (true depth scale) + grain-size panels for event layers in GUAC-29A.

One figure per event (E2 flood, E3 earthquake). Both figures share:
  * the same depth window LENGTH  -> identical cm-per-inch vertical scale
  * the same x-axis limits on every panel
  * the same panel widths / positions (set in inches, not fractions)
  * 24 pt Liberation Sans for all text and tick labels
  * a height of 575 pt (7.99 in), matching the 24 pt cartoon PDFs

Panels (left -> right):
  core photo | mean grain size (um, log) | mean grain size (phi) | class fractions (%) | legend
The shaded envelope on the two mean panels is mean +/- 1 sorting (Folk & Ward, phi),
i.e. the width of the grain-size spread, shown in both um and phi.

Run:  python plot_event_grain_size.py
"""
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
from PIL import Image

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
OUT = HERE / "output"

GSD_CSV = DATA / "GUAC_29A_GSD_Removedless1um_FolkWard_Stats.csv"
EVENTS_XLSX = DATA / "Lachua_Grain_Size_Samples_N_Events_with_Ages.xlsx"
CORE_JPG = DATA / "GUAC_29A_core.jpg"

# ---------------------------------------------------------------- calibration
# Ruler on the core photo: the grey/white blocks change every 10 cm. Edge
# positions (pixel boundary, x) measured at 10, 20 ... 100 cm. Linear fit:
# max residual 1.5 px (0.75 mm), sd 0.34 mm.
RULER_CM = np.array([10, 20, 30, 40, 50, 60, 70, 80, 90, 100], float)
RULER_PX = np.array([200, 398, 598, 796, 995, 1193, 1393, 1591, 1792, 1989], float)
PX_PER_CM, PX_AT_0CM = np.polyfit(RULER_CM, RULER_PX, 1)
CORE_ROWS = (14, 124)          # image rows covering the sediment (excludes liner edge / ruler)


def cm_to_px(cm):
    return PX_PER_CM * np.asarray(cm, float) + PX_AT_0CM


def px_to_cm(px):
    return (np.asarray(px, float) - PX_AT_0CM) / PX_PER_CM


# ---------------------------------------------------------------- figure style
FONT = 24
plt.rcParams.update({
    "font.family": "Liberation Sans",
    "font.size": FONT, "axes.titlesize": FONT, "axes.labelsize": FONT,
    "xtick.labelsize": FONT, "ytick.labelsize": FONT, "legend.fontsize": FONT,
    "axes.linewidth": 1.5, "xtick.major.width": 1.5, "ytick.major.width": 1.5,
    "xtick.major.size": 7, "ytick.major.size": 7, "xtick.minor.size": 4, "ytick.minor.size": 4,
    "svg.fonttype": "none", "pdf.fonttype": 42,
})
INK = "#1F2D3A"

FIG_H = 575 / 72               # same height as the cartoon PDFs
DEPTH_SPAN_CM = 9.0            # identical vertical scale for every event figure
# panel layout in inches: (left, width)
BOTTOM, PLOT_H = 1.55, 5.45    # bottom margin holds x labels; top holds title
CORE_W_CM = 3.0                # width of the core strip shown (true aspect)
CM_PER_IN = DEPTH_SPAN_CM / PLOT_H
CORE_W_IN = CORE_W_CM / CM_PER_IN
L_CORE = 1.75
L_UM = L_CORE + CORE_W_IN + 0.45
PANEL_W = 3.6
L_PHI = L_UM + PANEL_W + 0.45
L_FRAC = L_PHI + PANEL_W + 0.45
L_LEG = L_FRAC + PANEL_W + 0.3
FIG_W = L_LEG + 6.2

# shared x scales (um and phi are the same mapping: phi = -log2(um / 1000))
UM_LIM = (2.0, 250.0)
UM_TICKS = [4, 16, 63, 250]           # = 8, 6, 4, 2 phi
UM_MINOR = [2, 8, 31, 125]             # = 9, 7, 5, 3 phi
PHI_LIM = (9.0, 2.0)           # inverted so coarser plots to the right in both panels
PHI_TICKS = [8, 6, 4, 2]
PHI_MINOR = [9, 7, 5, 3]
XLABEL_Y = -0.115                      # same x-label position (axes fraction) on every panel

# Folk-Ward classes merged so the legend fits; fixed order fine -> coarse
CLASSES = [
    ("Clay", ["Clay(>8)"], "> 8 φ", "< 3.9 µm", "#5B4636"),
    ("V. fine silt", ["V.fine silt(7-8)"], "7–8 φ", "3.9–7.8 µm", "#7E6650"),
    ("Fine silt", ["Fine silt(6-7)"], "6–7 φ", "7.8–16 µm", "#A08468"),
    ("Medium silt", ["Medium silt(5-6)"], "5–6 φ", "16–31 µm", "#C2A584"),
    ("Coarse silt", ["Coarse silt(4-5)"], "4–5 φ", "31–63 µm", "#DCC6A2"),
    ("V. fine sand", ["V.fine sand(3-4)"], "3–4 φ", "63–125 µm", "#F1E2B8"),
    ("Fine sand", ["Fine sand(2-3)"], "2–3 φ", "125–250 µm", "#F2C66D"),
    ("Medium sand", ["Medium sand(1-2)"], "1–2 φ", "250–500 µm", "#E39A3B"),
    ("Coarse sand +", ["Coarse sand(0-1)", "V.coarse sand(-1to0)", "Granule(<-1)"],
     "< 1 φ", "> 500 µm", "#A8521C"),
]

EVENTS_TO_PLOT = [
    # event, figure title, colour (from the cartoons), window top (cm)
    ("E2", "Flood-triggered", "#E39A3B", 12.6),
    ("E3", "Earthquake-triggered", "#3E7CA6", 21.9),
]


# ---------------------------------------------------------------- data
def parse_interval(label):
    nums = re.findall(r"\d+(?:\.\d+)?", str(label))
    top, base = float(nums[0]), float(nums[1])
    assert base > top and base - top <= 1.0, label
    return top, base


def load_gsd():
    d = pd.read_csv(GSD_CSV)
    d[["top", "base"]] = [parse_interval(s) for s in d["Sample"]]
    d = d.sort_values("top").reset_index(drop=True)
    for name, cols, *_ in CLASSES:
        d[name] = d[cols].sum(axis=1)
    return d


def load_events():
    x = pd.read_excel(EVENTS_XLSX, header=0)
    ev = pd.DataFrame({
        "event": x.iloc[:, 7], "age": x["age_CE"], "unc": x["age_unc"],
        "top": pd.to_numeric(x["29A_top_cm"], errors="coerce"),
        "base": pd.to_numeric(x["29A_base_cm"], errors="coerce"),
    })
    return ev.set_index("event")


def steps(top, base, values):
    """Depth/value arrays for a step plot, broken (NaN) where samples are not contiguous."""
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


# ---------------------------------------------------------------- plotting
def axes_at(fig, left, width):
    return fig.add_axes([left / FIG_W, BOTTOM / FIG_H, width / FIG_W, PLOT_H / FIG_H])


def core_crop(img, d0, d1):
    """Crop the photo to [d0, d1] cm and rotate so depth increases downward."""
    x0, x1 = int(np.floor(cm_to_px(d0))), int(np.ceil(cm_to_px(d1)))
    r0, r1 = CORE_ROWS
    mid, half = (r0 + r1) / 2, CORE_W_CM * PX_PER_CM / 2
    c0, c1 = int(round(mid - half)), int(round(mid + half))
    crop = np.asarray(img)[c0:c1, x0:x1]
    strip = np.transpose(crop, (1, 0, 2))          # rows = depth
    extent = [0, (c1 - c0) / PX_PER_CM, float(px_to_cm(x1)), float(px_to_cm(x0))]
    return strip, extent


def event_lines(ax, ev, color):
    for z in (ev.top, ev.base):
        ax.axhline(z, color=color, lw=2.5, ls="--", zorder=5)


def plot_event(code, title, color, d0, gsd, events, img, suffix=""):
    ev = events.loc[code]
    d1 = d0 + DEPTH_SPAN_CM
    sub = gsd[(gsd.base > d0) & (gsd.top < d1)]
    fig = plt.figure(figsize=(FIG_W, FIG_H))

    # --- core photo
    ax = axes_at(fig, L_CORE, CORE_W_IN)
    strip, ext = core_crop(img, d0, d1)
    ax.imshow(strip, extent=ext, aspect="equal", interpolation="lanczos")
    ax.set_xlim(0, CORE_W_CM); ax.set_ylim(d1, d0)
    ax.set_xticks([])
    ax.set_ylabel("Depth (cm)")
    ax.yaxis.set_major_locator(mticker.MultipleLocator(1))
    ax.yaxis.set_minor_locator(mticker.MultipleLocator(0.5))
    event_lines(ax, ev, color)
    ax.set_xlabel("GUAC-29A"); ax.xaxis.set_label_coords(0.5, XLABEL_Y)
    core_ax = ax

    # --- mean (um) with +/- sorting envelope
    def mean_panel(left, xform, xlim, ticks, minor, xlabel, log):
        a = axes_at(fig, left, PANEL_W)
        for g in contiguous_runs(sub):
            lo = xform(g.Mean_phi + g.Sorting_phi)
            hi = xform(g.Mean_phi - g.Sorting_phi)
            y, vlo = steps(g.top.values, g.base.values, lo.values)
            _, vhi = steps(g.top.values, g.base.values, hi.values)
            a.fill_betweenx(y, vlo, vhi, color=color, alpha=0.22, lw=0, step=None)
        y, v = steps(sub.top.values, sub.base.values, xform(sub.Mean_phi).values)
        a.plot(v, y, color=INK, lw=3)
        a.plot(xform(sub.Mean_phi), (sub.top + sub.base) / 2, "o", ms=8, color=INK)
        if log:
            a.set_xscale("log", base=2)
        a.set_xlim(*xlim); a.set_ylim(d1, d0)
        a.set_xticks(ticks); a.set_xticklabels([str(t) for t in ticks])
        a.set_xticks(minor, minor=True); a.set_xticklabels([], minor=True)
        a.yaxis.set_major_locator(core_ax.yaxis.get_major_locator())
        a.yaxis.set_minor_locator(core_ax.yaxis.get_minor_locator())
        a.tick_params(labelleft=False)
        a.grid(axis="x", which="both", color="#C8CDD3", lw=1)
        a.set_xlabel(xlabel); a.xaxis.set_label_coords(0.5, XLABEL_Y)
        event_lines(a, ev, color)
        return a

    mean_panel(L_UM, lambda p: 1000 * 2.0 ** (-p), UM_LIM, UM_TICKS, UM_MINOR, "Mean (µm)", True)
    mean_panel(L_PHI, lambda p: p, PHI_LIM, PHI_TICKS, PHI_MINOR, "Mean (φ)", False)

    # --- stacked class fractions
    a = axes_at(fig, L_FRAC, PANEL_W)
    for g in contiguous_runs(sub):
        left = np.zeros(len(g))
        for name, _, _, _, c in CLASSES:
            y, v0 = steps(g.top.values, g.base.values, left)
            _, v1 = steps(g.top.values, g.base.values, left + g[name].values)
            a.fill_betweenx(y, v0, v1, color=c, lw=0)
            left = left + g[name].values
    a.set_xlim(0, 100); a.set_ylim(d1, d0)
    a.set_xticks([0, 50, 100]); a.set_xticks([25, 75], minor=True)
    a.yaxis.set_major_locator(core_ax.yaxis.get_major_locator())
    a.yaxis.set_minor_locator(core_ax.yaxis.get_minor_locator())
    a.tick_params(labelleft=False)
    a.set_xlabel("Fraction (%)"); a.xaxis.set_label_coords(0.5, XLABEL_Y)
    event_lines(a, ev, color)

    # --- legend (fractions, coarse on top) + envelope key
    handles = [plt.Rectangle((0, 0), 1, 1, fc=c) for *_, c in CLASSES][::-1]
    labels = [f"{n}  {p} · {u}" for n, _, p, u, _ in CLASSES][::-1]
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
    stem = OUT / f"GUAC29A_{code}_core_grain_size{suffix}"
    for ext_ in ("png", "pdf", "svg"):
        fig.savefig(f"{stem}.{ext_}", dpi=300 if ext_ == "png" else None)
    plt.close(fig)
    return stem, sub


# ---------------------------------------------------------------- checks
def run_checks(gsd, events, img):
    """Self-tests; raise AssertionError if anything is off."""
    # 1. ruler calibration: residuals < 1 mm
    res_mm = 10 * np.abs(px_to_cm(RULER_PX) - RULER_CM)
    assert res_mm.max() < 1.0, res_mm
    # 2. re-detect the 10 cm block edges in the image and check they land on 10 cm multiples
    row = np.asarray(img.convert("L"), float)[150]
    s = np.convolve(row, np.ones(5) / 5, "same")
    d = np.abs(np.diff(s))
    for cm in RULER_CM[:-1]:
        i = np.arange(int(cm_to_px(cm)) - 6, int(cm_to_px(cm)) + 6)
        w = np.where(d[i] > 20, d[i], 0)          # block edge, |d|-weighted centroid
        assert w.sum() > 0, f"no ruler edge near {cm} cm"
        edge = (w * (i + 1)).sum() / w.sum()
        assert abs(px_to_cm(edge) - cm) < 0.1, (cm, px_to_cm(edge))
    # 3. the crop maps back exactly: first/last rows of the strip sit at the window limits
    strip, ext = core_crop(img, 21.9, 30.9)
    assert abs(ext[3] - 21.9) < 1 / PX_PER_CM and abs(ext[2] - 30.9) < 1 / PX_PER_CM
    assert strip.shape[0] == int(np.ceil(cm_to_px(30.9))) - int(np.floor(cm_to_px(21.9)))
    # 4. grain-size data
    assert len(gsd) == 86
    tot = gsd[[c[0] for c in CLASSES]].sum(axis=1)
    assert np.allclose(tot, 100, atol=0.01), tot.describe()
    assert np.allclose(1000 * 2.0 ** -gsd.Mean_phi, gsd.Mean_um, rtol=1e-9)
    assert (gsd.base - gsd.top).between(0.4, 0.6).all()
    # 5. events
    assert events.loc["E2", ["top", "base"]].tolist() == [15.1, 19.1]
    assert events.loc["E3", ["top", "base"]].tolist() == [22.6, 30.15]
    assert events.loc["E2", "age"] == 2000.1 and events.loc["E3", "age"] == 1977.9
    print("all checks passed: ruler residual max %.2f mm, px/cm %.3f" % (res_mm.max(), PX_PER_CM))


def main():
    gsd, events = load_gsd(), load_events()
    img = Image.open(CORE_JPG).convert("RGB")
    run_checks(gsd, events, img)
    for code, title, color, d0 in EVENTS_TO_PLOT:
        ev = events.loc[code]
        assert d0 <= ev.top and ev.base <= d0 + DEPTH_SPAN_CM, f"{code} outside window"
        stem, sub = plot_event(code, title, color, d0, gsd, events, img)
        print(f"{code}: window {d0:.1f}-{d0 + DEPTH_SPAN_CM:.1f} cm, {len(sub)} samples "
              f"({sub.top.min()}-{sub.base.max()} cm), mean {sub.Mean_um.min():.1f}-"
              f"{sub.Mean_um.max():.1f} µm -> {stem.name}.png/.pdf/.svg")


if __name__ == "__main__":
    main()
