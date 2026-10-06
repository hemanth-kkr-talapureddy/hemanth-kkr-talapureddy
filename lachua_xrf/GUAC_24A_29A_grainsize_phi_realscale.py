# =============================================================================
# GUAC-24A + GUAC-29A - core photo | D10-D50-D90 (phi) | Sorting (phi)
# at REAL depth scale (1 cm of core = 1 cm on the page), both cores on ONE
# common depth frame so they line up with each other.
# Paste this whole file into ONE Colab cell and run it.
#
# Makes:
#   GUAC-24A_29A_realscale.pdf   24A | common depth scale | 29A   (one page)
#   GUAC-24A_realscale.pdf       24A alone  (same page height and depth frame)
#   GUAC-29A_realscale.pdf       29A alone  (same page height and depth frame)
#   + a zip with all PDFs
#
# Both cores use the SAME x-axis limits, the SAME panel widths and the SAME
# font size (FONT = 24 pt for all text and scales). Every event of each core
# is drawn as a coloured band, also the ones without grain-size samples
# (their label gets a *). Event bands = varve-count depths
# (Lachua_Varve_counts.xlsx); grain-size points = the samples in each event's
# GSD interval (Lachua_Graine_Siza_Samples_N_Events_with_Ages.xlsx).
# =============================================================================
import os, re, glob, zipfile
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')                       # draw to files only (robust in Colab)
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.ticker import MultipleLocator, MaxNLocator
from matplotlib.lines import Line2D

# ----------------------------------------------------------------------------
# SETTINGS
# ----------------------------------------------------------------------------
FOLDER = '/content'                 # where your files are
UPLOAD = True                       # True = show the upload box first
# Two event files (words in their names):
#  VARVE_FILE = event depths from the varve counting (sheet 'Events'): the
#               CORRECT event depths -> used for the coloured bands, the bars
#               next to the photo, the labels and the contact lines
#  GSD_EVENT_FILE = event depths of the grain-size samples (wide table
#               event | in_24A | 24A_top_cm | 24A_base_cm | ...): decides
#               WHICH grain-size points belong to which event
VARVE_FILE, VARVE_SHEET = 'varve counts', 'Events'
GSD_EVENT_FILE = 'samples events'

# one entry per core, in the order they are drawn (left -> right)
#   stats  = words in the grain-size stats csv name
#   photo  = words in the core photo name (jpg/png/tif)
#   window = depth window plotted, cm (photo cropped to it)
#   photo_top_cm / photo_bot_cm = depth at the image's top / bottom edge,
#            used ONLY if the ruler cannot be read automatically
#   keep   = part of the photo width kept (0 = left edge, 1 = right edge)
CORES = {
 'GUAC-24A': dict(stats='stats 24A', photo='24A', window=(10.5, 93.8),
                  photo_top_cm=-0.09, photo_bot_cm=105.65, keep=(0.0, 1.0)),
 'GUAC-29A': dict(stats='stats 29A', photo='29A', window=(12.5, 90.8),
                  photo_top_cm=0.04, photo_bot_cm=100.55, keep=(0.0, 1.0)),
}

# event depths (cm) - used ONLY if the event files are not found.
# From Lachua_Varve_counts.xlsx, sheet 'Events' (29A E19 is 87.75-87.75 there,
# zero thickness, so its GSD interval is used instead)
VARVE_BACKUP = {
 'GUAC-24A': {1: (11.75, 11.9), 2: (12.9, 17.4), 3: (19.35, 26.2), 4: (27.5, 27.85),
              5: (28.35, 30.8), 6: (30.9, 33.0), 7: (33.0, 48.0), 8: (51.7, 53.0),
             10: (55.7, 55.85), 11: (56.5, 62.0), 14: (63.5, 63.75), 15: (66.25, 66.5),
             16: (66.7, 67.1), 17: (68.8, 69.75), 18: (70.15, 70.55), 19: (71.3, 72.7),
             20: (75.2, 79.8)},
 'GUAC-29A': {1: (13.4, 14.1), 2: (15.1, 19.1), 3: (22.6, 30.15), 4: (31.3, 33.2),
              5: (33.55, 36.8), 6: (36.8, 41.0), 7: (41.0, 49.1), 8: (53.25, 55.8),
              9: (56.8, 58.1), 10: (60.2, 60.6), 11: (61.2, 64.1), 12: (64.65, 66.6),
             13: (67.15, 67.5), 14: (67.9, 68.3), 15: (70.8, 71.15), 16: (71.25, 71.5),
             17: (73.15, 75.4), 18: (76.65, 80.0), 19: (87.75, 87.75)},
}
# From Lachua_Graine_Siza_Samples_N_Events_with_Ages.xlsx
GSD_BACKUP = {
 'GUAC-24A': {1: (11.5, 12.4), 2: (13.0, 17.2), 3: (19.6, 26.0), 4: (27.4, 28.0),
              5: (28.6, 30.6), 6: (30.6, 33.0), 7: (33.0, 47.8), 8: (51.6, 52.8),
             10: (55.6, 56.2), 11: (56.6, 61.6), 14: (63.4, 63.8), 15: (66.2, 66.5),
             16: (66.5, 66.9), 17: (68.6, 69.6), 18: (70.2, 70.6), 19: (71.4, 72.4),
             20: (75.2, 79.6)},                                  # E9 E12 E13 absent
 'GUAC-29A': {1: (13.4, 14.2), 2: (15.1, 19.1), 3: (22.6, 30.15), 4: (31.3, 33.2),
              5: (33.55, 36.8), 6: (36.8, 41.0), 7: (41.0, 49.1), 8: (53.25, 55.8),
              9: (57.0, 58.0), 10: (60.0, 60.6), 11: (61.5, 64.0), 12: (65.0, 66.5),
             13: (67.2, 67.6), 14: (67.9, 68.3), 15: (70.8, 71.15), 16: (71.25, 71.5),
             17: (73.5, 75.5), 18: (77.0, 80.0), 19: (82.5, 87.5)},  # E20 absent
}

# photo handling
PHOTO_ROTATE = 'auto'     # 'auto' = image wider than tall -> core lying down
                          # (top on the LEFT) -> turned 90 deg clockwise.
                          # Or force 'cw', 'ccw' or None.
AUTO_CALIBRATE = True     # read the image-edge depths from the ruler's
                          # 10-cm white/grey blocks in the photo itself
PHOTO_PX_PER_CM = 100     # photo reduced to this before drawing (fast)

# event colours from your event table (8 colours, repeating every 8 events)
BASE8 = ['#f6a362', '#2a9e8f', '#eac46b', '#6d8fae',
         '#e86f52', '#8ab27d', '#b5838f', '#9c8ab9']
BAND_ALPHA, LINE_DARKEN = 0.45, 0.35  # bands as in the table; curves darker

# grain size
CLAY_FINE_UM = 1.0        # clay class runs 3.9 -> 1 µm ("<1 µm removed")
PHI_COARSE_RIGHT = True   # phi axis reversed: COARSER grains plot to the RIGHT

# page layout (cm) - identical for both cores
FONT = 24                 # pt, ALL text and scale numbers
SCALE = 1.0               # 1.0 = true size (1 cm core = 1 cm page)
D_W_CM, SORT_W_CM = 8.0, 5.0     # width of the D10-D50-D90 and Sorting boxes
GAP_CM = 0.8                     # gap between boxes
LABEL_W_CM = 3.6                 # room for the E## labels right of Sorting
LEFT_CM = 3.4                    # room for each core's depth numbers + title
MID_W_CM = 4.2                   # room for the common depth scale between cores
TOP_CM, BOTTOM_CM = 7.0, 2.6     # room for titles/legend above, note below
CONTACT_TOL = 0.15               # events closer than this (cm) = touching
OUT = '/content/realscale_24A_29A'

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['font.size'] = FONT
for k in ('axes.labelsize', 'axes.titlesize', 'xtick.labelsize', 'ytick.labelsize',
          'legend.fontsize', 'figure.titlesize'):
    matplotlib.rcParams[k] = FONT
matplotlib.rcParams['axes.linewidth'] = 1.5
for a in ('x', 'y'):
    matplotlib.rcParams[f'{a}tick.major.size'] = 8
    matplotlib.rcParams[f'{a}tick.minor.size'] = 4
    matplotlib.rcParams[f'{a}tick.major.width'] = 1.5
    matplotlib.rcParams[f'{a}tick.minor.width'] = 1.0

# ----------------------------------------------------------------------------
# 0. upload + file finder
# ----------------------------------------------------------------------------
if UPLOAD:
    try:
        from google.colab import files
        print("Upload: both stats csv files, BOTH event workbooks (varve counts + GSD samples) and both core photos.")
        files.upload()
    except Exception as ex:
        print("(no upload box:", ex, ")")
print("Files in", FOLDER, ":", sorted(os.path.basename(f) for f in glob.glob(FOLDER + '/*')
                                      if os.path.isfile(f)))

def norm(x):
    return re.sub(r'[^a-z0-9]', '', str(x).lower())

def find(spec, exts):
    words = [norm(w) for w in re.split(r'[\s_\-]+', spec) if norm(w)]
    hits = [f for f in glob.glob(FOLDER + '/*')
            if f.lower().endswith(exts) and all(w in norm(os.path.basename(f)) for w in words)
            and 'realscale' not in os.path.basename(f).lower()
            and not os.path.basename(f).startswith(('GS_STATS_', 'EVENT_TABLE_'))]
    return max(hits, key=os.path.getmtime) if hits else None

# ----------------------------------------------------------------------------
# 1. events
# ----------------------------------------------------------------------------
def read_gsd_events():
    """GSD event workbook, wide layout: event | in_24A | 24A_top_cm |
    24A_base_cm | ... | event | in_29A | 29A_top_cm | 29A_base_cm | ..."""
    p = find(GSD_EVENT_FILE, ('.xlsx', '.xls', '.csv'))
    if p is None:
        print(f"   !!! no file with '{GSD_EVENT_FILE}' in its name - using GSD_BACKUP")
        return GSD_BACKUP
    t = pd.read_excel(p) if p.lower().endswith(('.xlsx', '.xls')) else pd.read_csv(p)
    cols = list(t.columns)
    out = {}
    for core in CORES:
        tag = norm(core.split('-')[-1])                       # '24a'
        top = next((c for c in cols if tag in norm(c) and 'top' in norm(c)), None)
        bas = next((c for c in cols if tag in norm(c) and ('base' in norm(c) or 'bot' in norm(c))), None)
        if top is None or bas is None:
            print(f"   {core}: no top/base columns in {os.path.basename(p)} - using GSD_BACKUP")
            out[core] = GSD_BACKUP[core]; continue
        # event-name column = the nearest 'event' column left of the top column
        evc = [c for c in cols[:cols.index(top)] if norm(c).startswith('event')][-1]
        ev = {}
        for _, r in t.iterrows():
            m = re.search(r'(\d+)', str(r[evc]))
            a, b = pd.to_numeric(r[top], errors='coerce'), pd.to_numeric(r[bas], errors='coerce')
            if m and np.isfinite(a) and np.isfinite(b):
                ev[int(m.group(1))] = (float(min(a, b)), float(max(a, b)))
        out[core] = ev
    print(f"   GSD sample intervals from {os.path.basename(p)}")
    return out

def read_varve_events():
    """Varve-count workbook, sheet 'Events': row 1 = core names
    (GUAC-24A-1G-1 ...), row 2 = Year | Depth (in core) | Depth (adjusted),
    then TWO rows per event (top, base) for E1, E2, E3 ... in order; an
    empty pair = event absent in that core."""
    p = find(VARVE_FILE, ('.xlsx', '.xls'))
    if p is None:
        print(f"   !!! no file with '{VARVE_FILE}' in its name - using VARVE_BACKUP")
        return VARVE_BACKUP
    xl = pd.ExcelFile(p)
    sheet = next((n for n in xl.sheet_names if norm(n) == norm(VARVE_SHEET)), None)
    if sheet is None:
        print(f"   !!! no sheet '{VARVE_SHEET}' in {os.path.basename(p)} - using VARVE_BACKUP")
        return VARVE_BACKUP
    t = pd.read_excel(p, sheet_name=sheet, header=None)
    # the event rows end at the first fully empty row block / the 'Lat' table
    stop = next((i for i in range(2, len(t)) if norm(t.iat[i, 0]) == 'lat'), len(t))
    dcols = {}
    for core in CORES:
        tag = norm(core.split('-')[-1])
        j = next((j for j in range(t.shape[1]) if tag in norm(t.iat[0, j])), None)
        dc = None if j is None else next((k for k in range(j, min(j + 3, t.shape[1]))
                                          if 'incore' in norm(t.iat[1, k])), None)
        dcols[core] = dc
    num = t.iloc[2:stop].apply(pd.to_numeric, errors='coerce')
    used = [c for c in dcols.values() if c is not None]
    last = max(i for i in num.index if num.loc[i, used].notna().any())
    out = {}
    for core, dc in dcols.items():
        if dc is None:
            print(f"   {core}: not found in {os.path.basename(p)} - using VARVE_BACKUP")
            out[core] = VARVE_BACKUP[core]; continue
        ev = {}
        for k, i in enumerate(range(2, last + 1, 2)):
            a, b = num.at[i, dc], num.at[i + 1, dc] if i + 1 in num.index else np.nan
            if np.isfinite(a) and np.isfinite(b):
                ev[k + 1] = (float(min(a, b)), float(max(a, b)))
        out[core] = ev
    print(f"   event depths from {os.path.basename(p)}, sheet '{sheet}'")
    return out

print("\n1. events")
GSD_EV, VARVE_EV = read_gsd_events(), read_varve_events()
EVENTS = {}          # depths used for bands / labels / photo bars
for core in CORES:
    EVENTS[core] = {}
    print(f"   {core}:   event   varve depth (band)   GSD interval (points)")
    for e in sorted(set(VARVE_EV.get(core, {})) | set(GSD_EV.get(core, {}))):
        v, g = VARVE_EV.get(core, {}).get(e), GSD_EV.get(core, {}).get(e)
        note = ''
        if v is not None and v[1] - v[0] < 0.05:              # zero thickness
            note = '  <- varve depth has no thickness: band uses the GSD interval'
            v = None
        if v is None and g is not None and not note:
            note = '  <- not in the varve file: band uses the GSD interval'
        if v is not None and g is not None and abs((v[0] + v[1]) - (g[0] + g[1])) / 2 > 2:
            note = '  <- CHECK: varve and GSD depths differ by more than 2 cm'
        EVENTS[core][e] = v if v is not None else g
        f_ = lambda x: f"{x[0]:6.2f}-{x[1]:6.2f}" if x else '     absent  '
        print(f"            E{e:<3d}   {f_(VARVE_EV.get(core, {}).get(e))}        "
              f"{f_(g)}{note}")

COLOR = {e: BASE8[(e - 1) % 8] for e in range(1, 41)}
DARK = {e: tuple(v * (1 - LINE_DARKEN) for v in mcolors.to_rgb(c)) for e, c in COLOR.items()}

# ----------------------------------------------------------------------------
# 2. grain-size stats of each core -> D10 / D50 / D90 / Sorting in phi
# ----------------------------------------------------------------------------
def label_depth(s):                   # '13.0cm', '15.5 - 16cm', '16.5-17.0cm-G2'
    s = str(s)
    m = re.findall(r'(\d+\.?\d*)\s*-\s*(\d+\.?\d*)', s)
    if m: a, b = float(m[-1][0]), float(m[-1][1]); return (a + b) / 2
    m = re.findall(r'(\d+\.?\d*)\s*cm', s, re.I) or re.findall(r'\d+\.\d+', s)
    return float(m[-1]) if m else np.nan

CLASSES = [('granule', -2, -1), ('vcoarsesand', -1, 0), ('coarsesand', 0, 1),
           ('mediumsand', 1, 2), ('finesand', 2, 3), ('vfinesand', 3, 4),
           ('coarsesilt', 4, 5), ('mediumsilt', 5, 6), ('finesilt', 6, 7),
           ('vfinesilt', 7, 8), ('clay', 8, -np.log2(CLAY_FINE_UM / 1000))]
EDGES = [CLASSES[0][1]] + [hi for _, _, hi in CLASSES]

def phi_at(row, pct_coarser):
    """phi where pct_coarser % of the sample is coarser (cumulative curve over
    the 1-phi Udden-Wentworth classes, straight line inside each class)."""
    if row.sum() <= 0: return np.nan
    cum = np.concatenate([[0], np.cumsum(row) / row.sum() * 100])
    return np.interp(pct_coarser, cum, EDGES)

def read_stats(core, spec):
    path = find(spec, ('.csv', '.xlsx', '.xls'))
    if path is None:
        raise SystemExit(f"STOP: no file with all of the words '{spec}' in its name in "
                         f"{FOLDER} (stats file of {core}).")
    raw = pd.read_excel(path) if path.lower().endswith(('.xlsx', '.xls')) \
          else pd.read_csv(path, encoding='utf-8-sig')
    label_col = next((c for c in raw.columns if norm(c) in ('sample', 'samplename', 'sampleid')),
                     raw.columns[0])
    depth = raw[label_col].map(label_depth).astype(float)
    col = lambda *k: next((c for c in raw.columns if norm(c).startswith(k)), None)
    c_sort = col('sortingphi', 'sorting')
    if c_sort is None:
        raise SystemExit(f"STOP: no Sorting column in {os.path.basename(path)}")
    cls = []
    for key, _, _ in CLASSES:
        c = next((c for c in raw.columns if norm(c).startswith(key)), None)
        if c is None:
            raise SystemExit(f"STOP: class column '{key}' not in {os.path.basename(path)}")
        cls.append(c)
    P = raw[cls].apply(pd.to_numeric, errors='coerce').fillna(0).to_numpy(float)
    pct = {p: np.array([phi_at(r, p) for r in P]) for p in (10, 16, 50, 84, 90)}
    d = pd.DataFrame({'depth': depth.values,
                      'D10': pct[90],      # phi: 10 % of the sample is FINER
                      'D50': pct[50],
                      'D90': pct[10],      # phi: 90 % of the sample is finer
                      'Sorting': pd.to_numeric(raw[c_sort], errors='coerce').values})
    d = d.dropna(subset=['depth']).groupby('depth', as_index=False).mean()
    print(f"   {core}: {os.path.basename(path)} - {len(raw)} rows, "
          f"{d.depth.min():.2f}-{d.depth.max():.2f} cm")
    c_mean = col('meanphi')
    if c_mean:                          # check of the method against your file
        mz = (pct[16] + pct[50] + pct[84]) / 3
        diff = np.nanmax(np.abs(mz - pd.to_numeric(raw[c_mean], errors='coerce').values))
        print(f"      check: Folk & Ward mean from the classes vs '{c_mean}': max "
              f"difference {diff:.3f} phi {'(OK)' if diff < 0.1 else '(LARGE!)'}")
    return d

print("\n2. grain-size stats")
C = {}                                  # everything per core
for core, s in CORES.items():
    t0, b0 = s['window']
    data = read_stats(core, s['stats'])
    ev = sorted([(e, t, b) for e, (t, b) in EVENTS[core].items() if b > t0 and t < b0],
                key=lambda x: x[1])
    contacts = [(ev[i][2] + ev[i + 1][1]) / 2 for i in range(len(ev) - 1)
                if abs(ev[i + 1][1] - ev[i][2]) <= CONTACT_TOL]
    # grain-size points of each event: the samples inside its GSD interval
    # (from the GSD event file); the varve interval if it has none there
    gi = {e: GSD_EV.get(core, {}).get(e, (t, b)) for e, t, b in ev}
    seg = {e: data[(data.depth >= gi[e][0] - 1e-6) & (data.depth <= gi[e][1] + 1e-6)]
           for e, t, b in ev}
    C[core] = dict(window=(t0, b0), data=data, ev=ev, contacts=contacts, seg=seg)
    print(f"      window {t0}-{b0} cm, {len(ev)} events: " +
          ", ".join(f"E{e} {len(seg[e])}" for e, _, _ in ev) + "  (samples)")
    print(f"      no grain-size samples: " +
          (", ".join(f"E{e}" for e, _, _ in ev if len(seg[e]) == 0) or "none"))

# SAME x-limits for both cores (from all plotted samples of both cores)
def nice(lo, hi, step):
    return np.floor(lo / step) * step, np.ceil(hi / step) * step
dv = np.concatenate([c['seg'][e][v].dropna().values for c in C.values()
                     for e, _, _ in c['ev'] for v in ('D10', 'D50', 'D90')])
sv = np.concatenate([c['seg'][e]['Sorting'].dropna().values for c in C.values()
                     for e, _, _ in c['ev']])
D_LIM = nice(dv.min() - 0.1, dv.max() + 0.1, 0.5)
S_LIM = nice(sv.min() - 0.05, sv.max() + 0.05, 0.25)
print(f"   common x-limits: D10-D50-D90 {D_LIM[0]}-{D_LIM[1]} phi, "
      f"Sorting {S_LIM[0]}-{S_LIM[1]} phi")

# ----------------------------------------------------------------------------
# 3. photos
# ----------------------------------------------------------------------------
def ruler_calibration(gray):
    """Depth (cm) at the top and bottom image edge, read from the ruler.
    The ruler alternates white / grey every 10 cm. Every column in the left
    and right 12 % of the image is tried; the one with the most evenly spaced
    white/grey changes wins. Returns (top_cm, bottom_cm) or None."""
    h, w = gray.shape
    best, best_cv = None, None
    win = max(9, int(0.006 * h) | 1)       # ~0.6 cm: smooths over the ruler digits
    for x in list(range(0, max(1, int(0.12 * w)))) + list(range(int(0.88 * w), w)):
        prof = gray[:, x]
        lo_, hi_ = np.percentile(prof, [15, 85])
        if hi_ - lo_ < 40: continue                    # no white/grey contrast here
        thr = (lo_ + hi_) / 2                           # adapts to photo brightness
        on = np.convolve((prof > thr).astype(float), np.ones(win) / win, 'same') > 0.5
        edges = np.where(np.diff(on.astype(int)) != 0)[0] + 1
        if len(edges) < 5: continue
        step = np.median(np.diff(edges))
        if not 0.06 * h < step < 0.14 * h: continue    # ~10 cm of a ~1-m photo
        good = edges[np.abs(np.diff(np.concatenate([[edges[0] - step], edges])) - step) < 0.12 * step]
        if len(good) < 5: continue
        cv = np.std(np.diff(good)) / np.mean(np.diff(good))
        if best is None or (len(good), -cv) > (len(best), -best_cv):
            best, best_cv = good, cv
    if best is None or best_cv > 0.08: return None
    step = np.median(np.diff(best))
    n0 = int(round(best[0] / step))                       # first change = n0 x 10 cm
    depths = 10.0 * (n0 + np.round((best - best[0]) / step))
    a, b = np.polyfit(best, depths, 1)                    # depth = a * row + b
    if not 0.5 < a * step / 10 < 2: return None
    return b, a * h + b

def load_photo(core, s):
    p = find(s['photo'], ('.jpg', '.jpeg', '.png', '.tif', '.tiff'))
    if p is None:
        raise FileNotFoundError(f"no photo with '{s['photo']}' in its name in {FOLDER}")
    from PIL import Image
    Image.MAX_IMAGE_PIXELS = None
    img = Image.open(p)
    w0, h0 = img.size
    rot = PHOTO_ROTATE
    if rot == 'auto':
        rot = 'cw' if w0 > h0 else None
    print(f"   {core} photo {os.path.basename(p)}: {w0} x {h0} px "
          f"({'lying down -> turned 90 deg clockwise' if rot == 'cw' else 'turned 90 deg anticlockwise' if rot == 'ccw' else 'upright'})")
    long_px = max(w0, h0)                       # shrink big scans first
    k = min(1.0, PHOTO_PX_PER_CM * 110 / long_px)
    if k < 0.95:
        img.draft('RGB', (int(w0 * k), int(h0 * k)))
        img = img.convert('RGB')
        k2 = PHOTO_PX_PER_CM * 110 / max(img.size)
        if k2 < 0.95:
            img = img.resize((round(img.size[0] * k2), round(img.size[1] * k2)), Image.LANCZOS)
    img = img.convert('RGB')
    if rot == 'cw':    img = img.transpose(Image.ROTATE_270)
    elif rot == 'ccw': img = img.transpose(Image.ROTATE_90)
    w, h = img.size
    top_cm, bot_cm = s['photo_top_cm'], s['photo_bot_cm']
    if AUTO_CALIBRATE:
        cal = ruler_calibration(np.asarray(img.convert('L')).astype(float))
        if cal:
            top_cm, bot_cm = cal
            print(f"      ruler found: image edges at {top_cm:.2f} cm (top) and {bot_cm:.2f} cm (bottom)")
        else:
            print(f"      !!! ruler not found - using photo_top_cm / photo_bot_cm = "
                  f"{top_cm} / {bot_cm} cm. Check the photo against the depth axis!")
    per_cm = h / (bot_cm - top_cm)
    x0, x1 = int(round(s['keep'][0] * w)), int(round(s['keep'][1] * w))
    width_cm = (x1 - x0) / per_cm * SCALE
    if width_cm > 25:
        raise ValueError(f"the photo would be {width_cm:.0f} cm wide - its orientation "
                         f"is wrong. Set PHOTO_ROTATE to 'cw', 'ccw' or None.")
    t0, b0 = s['window']
    r0 = int(round((t0 - top_cm) * per_cm))
    r1 = int(round((b0 - top_cm) * per_cm))
    crop = np.asarray(img.crop((x0, max(r0, 0), x1, min(r1, h))))
    if r0 < 0 or r1 > h:
        print(f"      note: the photo covers only {top_cm:.1f}-{bot_cm:.1f} cm")
    print(f"      {per_cm:.1f} px/cm, true width {width_cm:.2f} cm")
    return dict(img=crop, top=top_cm + max(r0, 0) / per_cm,
                bot=top_cm + min(r1, h) / per_cm, width_cm=width_cm)

print("\n3. photos")
for core, s in CORES.items():
    try:
        C[core]['photo'] = load_photo(core, s)
    except Exception as ex:
        C[core]['photo'] = None
        print(f"   !!! {core} photo NOT used: {type(ex).__name__}: {ex}")
# the same photo-box width for both cores (the wider photo decides; the
# narrower one is centred in it, never stretched)
PHOTO_W_CM = max([c['photo']['width_cm'] for c in C.values() if c['photo']] or [6.0])

# ----------------------------------------------------------------------------
# 4. drawing
# ----------------------------------------------------------------------------
CM = 1 / 2.54
FRAME = (min(s['window'][0] for s in CORES.values()),      # common depth frame
         max(s['window'][1] for s in CORES.values()))
FRAME_H = (FRAME[1] - FRAME[0]) * SCALE
LABEL_GAP = 0.95 / SCALE            # min. distance between E## labels (cm of core)
STYLE = {'D10': ('o', ':'), 'D50': ('s', '-'), 'D90': ('^', '--')}
CORE_W = LEFT_CM + PHOTO_W_CM + 2 * GAP_CM + D_W_CM + SORT_W_CM + GAP_CM + LABEL_W_CM
os.makedirs(OUT, exist_ok=True)

def depth_ticks(ax):
    ax.yaxis.set_major_locator(MultipleLocator(5))
    ax.yaxis.set_minor_locator(MultipleLocator(1))

def draw_photo(ax, c):
    ph = c['photo']
    t0, b0 = c['window']
    if ph is None:
        ax.text(0.5, (t0 + b0) / 2, 'photo\nnot found', ha='center', va='center',
                rotation=90, color='0.5')
    else:
        pad = (PHOTO_W_CM - ph['width_cm']) / PHOTO_W_CM / 2   # centre a narrower photo
        ax.imshow(ph['img'], extent=[pad, 1 - pad, ph['bot'], ph['top']], aspect='auto',
                  interpolation='none')
    ax.set_xlim(0, 1); ax.set_xticks([])
    xb = 1 + 0.25 / PHOTO_W_CM                    # event bars just right of the photo
    for e, t, b in c['ev']:
        ax.plot([xb, xb], [max(t, t0), min(b, b0)], color=COLOR[e], lw=7,
                solid_capstyle='butt', clip_on=False, zorder=4)
    for y in c['contacts']:
        ax.plot([xb - 0.15 / PHOTO_W_CM, xb + 0.15 / PHOTO_W_CM], [y, y],
                color='k', lw=1.8, clip_on=False, zorder=5)

def draw_data(ax, c, vars_):
    t0, b0 = c['window']
    for e, t, b in c['ev']:
        ax.axhspan(max(t, t0), min(b, b0), color=COLOR[e], alpha=BAND_ALPHA, lw=0, zorder=0)
        s = c['seg'][e]
        for v in vars_:
            sv_ = s[s[v].notna()]
            if sv_.empty: continue
            mk, ls = STYLE.get(v, ('o', '-'))
            ax.plot(sv_[v], sv_.depth, color=DARK[e], lw=2.0, ls=ls, marker=mk, ms=7,
                    markeredgewidth=0, zorder=2)
    for y in c['contacts']:
        ax.axhline(y, color='0.25', lw=1.5, ls=(0, (4, 2)), zorder=3)
    if vars_[0].startswith('D'):
        ax.set_xlim(D_LIM[::-1] if PHI_COARSE_RIGHT else D_LIM)
        ax.xaxis.set_major_locator(MaxNLocator(5, steps=[1, 2, 5, 10]))
        ax.xaxis.set_minor_locator(MultipleLocator(0.5))
        ax.set_xlabel('D10, D50, D90 (φ)')
    else:
        ax.set_xlim(S_LIM)
        ax.xaxis.set_major_locator(MaxNLocator(3, steps=[1, 2, 2.5, 5, 10]))
        ax.xaxis.set_minor_locator(MultipleLocator(0.25))
        ax.set_xlabel('Sorting (φ)')
    ax.xaxis.tick_top(); ax.xaxis.set_label_position('top')
    ax.tick_params(axis='x', which='both', top=True, bottom=False)
    for sp in ('right', 'bottom'): ax.spines[sp].set_visible(False)

def spread(targets, gap, lo, hi):
    """Move label positions apart so they are at least `gap` apart, staying
    as close as possible to their targets and inside lo..hi."""
    y = list(targets)
    for _ in range(200):
        moved = False
        for i in range(1, len(y)):
            if y[i] - y[i - 1] < gap - 1e-9:
                mid = (y[i] + y[i - 1]) / 2
                y[i - 1], y[i] = mid - gap / 2, mid + gap / 2
                moved = True
        y = [min(max(v, lo + gap / 2), hi - gap / 2) for v in y]
        if not moved: break
    return y

def event_labels(fig, ax, c, x_right_cm, W):
    """E## labels right of the Sorting box, pulled apart where events are thin
    and close together, with a thin leader line to their band."""
    t0, b0 = c['window']
    tg = [(max(t, t0) + min(b, b0)) / 2 for e, t, b in c['ev']]
    ys = spread(tg, LABEL_GAP, FRAME[0], FRAME[1])
    w_ax = SORT_W_CM
    for (e, t, b), y0, y in zip(c['ev'], tg, ys):
        x_a, x_b, x_t = 1.0, 1 + 0.45 / w_ax, 1 + 0.65 / w_ax
        ax.plot([x_a, x_b], [y0, y], color=DARK[e], lw=1.5, clip_on=False,
                transform=ax.get_yaxis_transform())
        star = '*' if len(c['seg'][e]) == 0 else ''
        ax.text(x_t, y, f'E{e:02d}{star}', color=DARK[e], fontweight='bold',
                va='center', ha='left', transform=ax.get_yaxis_transform(), clip_on=False)

def place_core(fig, W, H, x, core, own_axis_label=True):
    """Draw one core starting at x (cm from the left page edge). The boxes of
    the core cover its own depth window, placed on the common depth frame."""
    c = C[core]
    t0, b0 = c['window']
    y_bot = BOTTOM_CM + (FRAME[1] - b0) * SCALE          # cm from page bottom
    box_h = (b0 - t0) * SCALE
    axes, x = [], x + LEFT_CM
    for kind, w in (('photo', PHOTO_W_CM), ('D', D_W_CM), ('S', SORT_W_CM)):
        ax = fig.add_axes([x / W, y_bot / H, w / W, box_h / H])
        if kind == 'photo': draw_photo(ax, c)
        else: draw_data(ax, c, ['D10', 'D50', 'D90'] if kind == 'D' else ['Sorting'])
        ax.set_ylim(b0, t0); depth_ticks(ax)
        if kind == 'photo':
            ax.set_ylabel(f'{core} depth (cm)' if own_axis_label else 'Depth (cm)')
        else:
            ax.tick_params(axis='y', which='both', labelleft=False)
        axes.append(ax); x += w + GAP_CM
    event_labels(fig, axes[-1], c, x, W)
    # core name and legend at the SAME page height for both cores
    # (just above the top of the common depth frame)
    frame_top = BOTTOM_CM + FRAME_H
    fig.text(axes[0].get_position().x0, (frame_top + 4.6) / H, core,
             fontweight='bold', va='bottom', ha='left')
    d_pos = axes[1].get_position()
    fig.legend(handles=[Line2D([], [], color='0.25', lw=2.0, ms=7, marker=STYLE[v][0],
                               ls=STYLE[v][1], label=v) for v in STYLE],
               loc='lower center', ncol=3, frameon=False, handlelength=1.4,
               handletextpad=0.3, columnspacing=0.6, borderaxespad=0,
               bbox_to_anchor=((d_pos.x0 + d_pos.x1) / 2, (frame_top + 3.1) / H))
    return axes

def common_scale(fig, W, H, x):
    """The extra depth scale: a black/white ruler over the whole common depth
    frame, with numbers every 5 cm. Both cores are drawn on this same frame,
    so any depth reads straight across from one core to the other."""
    w = 0.45
    ax = fig.add_axes([(x + MID_W_CM - 1.2) / W, BOTTOM_CM / H, w / W, FRAME_H / H])
    ax.set_ylim(FRAME[1], FRAME[0]); ax.set_xlim(0, 1); ax.set_xticks([])
    d = np.floor(FRAME[0] / 10) * 10
    while d < FRAME[1]:                                     # 10-cm blocks
        ax.axhspan(max(d, FRAME[0]), min(d + 10, FRAME[1]),
                   color='k' if int(d // 10) % 2 == 0 else 'w', lw=0)
        d += 10
    depth_ticks(ax)
    ax.tick_params(axis='y', which='both', left=True, right=True)
    ax.set_ylabel('Common depth (cm)')
    return ax

def legend_note(fig, W, H):
    fig.text(0.5 / W, 0.9 / H, "* event without grain-size samples (band only)   "
             "- - - contact between two events   "
             f"φ axis: {'coarser to the right' if PHI_COARSE_RIGHT else 'finer to the right'}",
             va='center', ha='left')

def save(fig, name):
    fig.savefig(f'{OUT}/{name}.pdf', facecolor='white')     # exact size, no cropping
    fig.savefig(f'{OUT}/{name}_preview.png', dpi=25, facecolor='white')
    plt.close(fig)
    print(f"   saved {name}.pdf  ({fig.get_figwidth() / CM:.1f} x {fig.get_figheight() / CM:.1f} cm)")

print(f"\n4. figures (common depth frame {FRAME[0]}-{FRAME[1]} cm, photo box "
      f"{PHOTO_W_CM:.2f} cm wide)")
H = TOP_CM + FRAME_H + BOTTOM_CM
core_list = list(CORES)

# 4a. both cores + common depth scale between them
W = CORE_W + MID_W_CM + CORE_W
fig = plt.figure(figsize=(W * CM, H * CM))
AX = {}
AX[core_list[0]] = place_core(fig, W, H, 0.0, core_list[0])
AX['common'] = common_scale(fig, W, H, CORE_W)
AX[core_list[1]] = place_core(fig, W, H, CORE_W + MID_W_CM, core_list[1])
legend_note(fig, W, H)
fig.canvas.draw()
LAYOUT = {k: [a.get_position().bounds for a in (v if isinstance(v, list) else [v])]
          for k, v in AX.items()}
save(fig, f'{core_list[0]}_{core_list[1].split("-")[-1]}_realscale')

# 4b. each core alone, same page height / depth frame (pages line up side by side)
for core in core_list:
    W1 = CORE_W
    fig = plt.figure(figsize=(W1 * CM, H * CM))
    place_core(fig, W1, H, 0.0, core)
    legend_note(fig, W1, H)
    save(fig, f'{core}_realscale')

# ----------------------------------------------------------------------------
# 5. self-check of the drawing geometry (cm on the page)
# ----------------------------------------------------------------------------
print("\n5. checks")
for core in core_list:
    b = LAYOUT[core]
    t0, b0 = C[core]['window']
    h_cm = b[0][3] * H
    print(f"   {core}: boxes {h_cm:.2f} cm tall for {b0 - t0:.2f} cm of core "
          f"({h_cm / (b0 - t0) / SCALE * 100:.2f} % of real scale); widths "
          + " | ".join(f"{bb[2] * W:.2f}" for bb in b) + " cm")
    # the same depth must sit at the same page height in the core and on the common scale
    for d in (t0, b0):
        y_core = (b[0][1] + b[0][3] * (b0 - d) / (b0 - t0)) * H
        cb = LAYOUT['common'][0]
        y_com = (cb[1] + cb[3] * (FRAME[1] - d) / (FRAME[1] - FRAME[0])) * H
        print(f"      {d:5.1f} cm: core {y_core:.3f} cm, common scale {y_com:.3f} cm "
              f"from page bottom (difference {abs(y_core - y_com) * 10:.2f} mm)")

# ----------------------------------------------------------------------------
# 6. show + download
# ----------------------------------------------------------------------------
try:
    from IPython.display import Image as Show, display
    for f in sorted(glob.glob(f'{OUT}/*_preview.png')):
        print(os.path.basename(f)); display(Show(f, height=900))
except Exception:
    pass
zip_path = '/content/GUAC_24A_29A_realscale_PDFs.zip'
with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as z:
    for f in sorted(glob.glob(f'{OUT}/*.pdf')):
        z.write(f, os.path.basename(f))
print(f"\n6. {zip_path}: {len(glob.glob(OUT + '/*.pdf'))} PDF(s)")
try:
    from google.colab import files
    files.download(zip_path)
except Exception:
    pass
