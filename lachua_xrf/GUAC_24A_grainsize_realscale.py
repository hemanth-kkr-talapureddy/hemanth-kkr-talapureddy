# =============================================================================
# GUAC-24A ONLY - grain size next to the core photo, at REAL depth scale
# Paste this whole file into ONE Colab cell and run it.
#
# Makes (1 cm of core = 1 cm on the page):
#   GUAC-24A_DATA_realscale.pdf      D10-D50-D90 | Sorting            (no photo)
#   GUAC-24A_realscale.pdf           photo | D10-D50-D90 | Sorting
# The DATA figure is saved FIRST, so it exists even if the photo step fails.
# =============================================================================
import os, re, glob, zipfile
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')                       # draw to files only (robust in Colab)
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.ticker import MultipleLocator, LogLocator, FuncFormatter, NullFormatter
from matplotlib.lines import Line2D

# ----------------------------------------------------------------------------
# SETTINGS
# ----------------------------------------------------------------------------
FOLDER     = '/content'             # where your files are
UPLOAD     = True                   # True = show the upload box first
STATS_FILE = 'stats 24A'            # words in the stats csv name
PHOTO_FILE = '24A'                  # words in the photo name (jpg/png/tif)
CORE       = 'GUAC-24A'
DEPTH_TOP, DEPTH_BOT = 13.0, 80.0   # plotted depth window, cm (photo cropped to it)

# photo calibration: core depth at the image's TOP and BOTTOM edge (from the
# ruler in your photo). rotate=None because the 24A photo is upright.
PHOTO_TOP_CM, PHOTO_BOT_CM, PHOTO_ROTATE = -0.09, 105.65, None
PHOTO_PX_PER_CM = 100               # photo reduced to this before drawing (fast)

# GUAC-24A event depths (cm), copied from
# Lachua_Graine_Siza_Samples_N_Events_with_Ages.xlsx  (E9, E12, E13 = ABSENT)
EVENTS = {1: (11.5, 12.4),  2: (13.0, 17.2),  3: (19.6, 26.0),  4: (27.4, 28.0),
          5: (28.6, 30.6),  6: (30.6, 33.0),  7: (33.0, 47.8),  8: (51.6, 52.8),
         10: (55.6, 56.2), 11: (56.6, 61.6), 14: (63.4, 63.8), 15: (66.2, 66.5),
         16: (66.5, 66.9), 17: (68.6, 69.6), 18: (70.2, 70.6), 19: (71.4, 72.4),
         20: (75.2, 79.6)}

# event colours from your event table (8 colours, repeating every 8 events)
BASE8 = ['#f6a362', '#2a9e8f', '#eac46b', '#6d8fae',
         '#e86f52', '#8ab27d', '#b5838f', '#9c8ab9']
BAND_ALPHA, LINE_DARKEN = 0.45, 0.35  # bands as in the table; curves darker

CLAY_FINE_UM = 1.0                  # clay class runs 3.9 -> 1 µm ("<1 µm removed")
SCALE = 1.0                         # 1.0 = true size
PANEL_W_CM, SORT_W_CM, GAP_CM = 4.0, 3.0, 0.5
MARGIN_CM = dict(left=1.6, right=1.1, top=2.0, bottom=0.6)
CONTACT_TOL = 0.15                  # events closer than this (cm) = touching
OUT = '/content/realscale_24A'

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['font.size'] = 8

# ----------------------------------------------------------------------------
# 0. upload
# ----------------------------------------------------------------------------
if UPLOAD:
    try:
        from google.colab import files
        print("Upload the 24A stats csv and the 24A core photo.")
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
            and not os.path.basename(f).startswith(('GS_STATS_', 'EVENT_TABLE_'))]
    return max(hits, key=os.path.getmtime) if hits else None

# ----------------------------------------------------------------------------
# 1. read the 24A stats csv
# ----------------------------------------------------------------------------
path = find(STATS_FILE, ('.csv', '.xlsx', '.xls'))
if path is None:
    raise SystemExit(f"STOP: no file with all of the words '{STATS_FILE}' in its name "
                     f"in {FOLDER}. Upload GUAC_24A_GSD_Removedless1um_FolkWard_Stats_1.csv "
                     f"(or change STATS_FILE).")
raw = pd.read_excel(path) if path.lower().endswith(('.xlsx', '.xls')) \
      else pd.read_csv(path, encoding='utf-8-sig')
print(f"\n1. stats file: {os.path.basename(path)}  ({len(raw)} rows)")
print("   columns:", list(raw.columns))

# depth of every sample, from the 'Sample' label: '13.0cm', '15.5 - 16cm', ...
def label_depth(s):
    s = str(s)
    m = re.findall(r'(\d+\.?\d*)\s*-\s*(\d+\.?\d*)', s)
    if m: a, b = float(m[-1][0]), float(m[-1][1]); return (a + b) / 2
    m = re.findall(r'(\d+\.?\d*)\s*cm', s, re.I) or re.findall(r'\d+\.\d+', s)
    return float(m[-1]) if m else np.nan

label_col = next((c for c in raw.columns if norm(c) in ('sample', 'samplename', 'sampleid')),
                 raw.columns[0])
depth = raw[label_col].map(label_depth).astype(float)
print(f"   depths from column '{label_col}': {depth.min():.2f}-{depth.max():.2f} cm, "
      f"{depth.notna().sum()} of {len(depth)} rows")
if depth.notna().sum() == 0:
    raise SystemExit(f"STOP: no depth could be read from '{label_col}'. "
                     f"First labels: {raw[label_col].head(3).tolist()}")

def col(*keys):
    return next((c for c in raw.columns if norm(c).startswith(keys)), None)
c_sort = col('sortingphi', 'sorting')
if c_sort is None:
    raise SystemExit(f"STOP: no Sorting column found in {list(raw.columns)}")

# D10 / D50 / D90 from the Udden-Wentworth class % (phi intervals, coarse->fine)
CLASSES = [('granule', -2, -1), ('vcoarsesand', -1, 0), ('coarsesand', 0, 1),
           ('mediumsand', 1, 2), ('finesand', 2, 3), ('vfinesand', 3, 4),
           ('coarsesilt', 4, 5), ('mediumsilt', 5, 6), ('finesilt', 6, 7),
           ('vfinesilt', 7, 8), ('clay', 8, -np.log2(CLAY_FINE_UM / 1000))]
cls = []
for key, lo, hi in CLASSES:
    c = next((c for c in raw.columns if norm(c).startswith(key)), None)
    if c is None:
        raise SystemExit(f"STOP: class column '{key}' not found in {list(raw.columns)}")
    cls.append(c)
edges = [CLASSES[0][1]] + [hi for _, _, hi in CLASSES]
P = raw[cls].apply(pd.to_numeric, errors='coerce').fillna(0).to_numpy(float)

def phi_at(row, pct_coarser):
    if row.sum() <= 0: return np.nan
    cum = np.concatenate([[0], np.cumsum(row) / row.sum() * 100])
    return np.interp(pct_coarser, cum, edges)

pct = {p: np.array([phi_at(r, p) for r in P]) for p in (10, 16, 50, 84, 90)}
data = pd.DataFrame({
    'depth':   depth.values,
    'D10':     1000 * 2.0 ** -pct[90],     # µm (10 % finer = 90 % coarser)
    'D50':     1000 * 2.0 ** -pct[50],
    'D90':     1000 * 2.0 ** -pct[10],
    'Sorting': pd.to_numeric(raw[c_sort], errors='coerce').values,
}).dropna(subset=['depth']).groupby('depth', as_index=False).mean()

c_mean = col('meanphi')
if c_mean:
    mz = (pct[16] + pct[50] + pct[84]) / 3
    d = np.abs(mz - pd.to_numeric(raw[c_mean], errors='coerce').values)
    print(f"   check: Folk & Ward mean from the classes vs '{c_mean}': "
          f"max difference {np.nanmax(d):.3f} phi {'(OK)' if np.nanmax(d) < 0.1 else '(LARGE!)'}")
print("   first rows:\n", data.head(3).round(3).to_string(index=False))

# ----------------------------------------------------------------------------
# 2. events and colours
# ----------------------------------------------------------------------------
COLOR = {e: BASE8[(e - 1) % 8] for e in EVENTS}
DARK = {e: tuple(v * (1 - LINE_DARKEN) for v in mcolors.to_rgb(c)) for e, c in COLOR.items()}
EV = sorted([(e, t, b) for e, (t, b) in EVENTS.items() if b >= DEPTH_TOP and t <= DEPTH_BOT],
            key=lambda x: x[1])
CONTACTS = [(EV[i][2] + EV[i + 1][1]) / 2 for i in range(len(EV) - 1)
            if EV[i + 1][1] - EV[i][2] <= CONTACT_TOL]

def in_event(t, b):
    return data[(data.depth >= t) & (data.depth <= b)]

n = {e: len(in_event(t, b)) for e, t, b in EV}
print(f"\n2. samples per event: " + ", ".join(f"E{e} {n[e]}" for e in n))
print(f"   total plotted: {sum(n.values())} samples   contacts at {CONTACTS} cm")
if sum(n.values()) == 0:
    raise SystemExit("STOP: no sample falls inside any event - check depths above.")

# ----------------------------------------------------------------------------
# 3. drawing
# ----------------------------------------------------------------------------
CM = 1 / 2.54
BOX_H = (DEPTH_BOT - DEPTH_TOP) * SCALE
os.makedirs(OUT, exist_ok=True)

def depth_axis(ax, labels):
    ax.set_ylim(DEPTH_BOT, DEPTH_TOP)
    ax.yaxis.set_major_locator(MultipleLocator(5))
    ax.yaxis.set_minor_locator(MultipleLocator(1))
    if labels: ax.set_ylabel('Depth (cm)')
    else: ax.tick_params(axis='y', which='both', labelleft=False)

STYLE = {'D10': ('o', ':'), 'D50': ('s', '-'), 'D90': ('^', '--')}

def draw_data(ax, vars_):
    vals = []
    for e, t, b in EV:
        ax.axhspan(t, b, color=COLOR[e], alpha=BAND_ALPHA, lw=0, zorder=0)
        seg = in_event(t, b)
        for v in vars_:
            s = seg[seg[v].notna()]
            if s.empty: continue
            vals += list(s[v])
            mk, ls = STYLE.get(v, ('o', '-'))
            ax.plot(s[v], s.depth, color=DARK[e], lw=0.9, ls=ls, marker=mk, ms=2.6,
                    markeredgewidth=0, zorder=2)
    for y in CONTACTS:
        ax.axhline(y, color='0.25', lw=0.8, ls=(0, (4, 2)), zorder=3)
    lo, hi = min(vals), max(vals)
    if vars_[0].startswith('D'):                          # µm, log axis
        ax.set_xscale('log'); ax.set_xlim(lo / 1.15, hi * 1.15)
        dec = np.log10(hi / lo)
        ax.xaxis.set_major_locator(LogLocator(10, subs=(1, 2, 5) if dec < 1.6 else (1, 5) if dec < 2.2 else (1,)))
        ax.xaxis.set_minor_locator(LogLocator(10, subs=np.arange(2, 10)))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f'{x:g}'))
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.set_xlabel('D10, D50, D90 (µm)')
        ax.legend(handles=[Line2D([], [], color='0.3', lw=0.9, ms=3, marker=STYLE[v][0],
                                  ls=STYLE[v][1], label=v) for v in vars_],
                  loc='lower center', ncol=3, fontsize=6.5, frameon=False,
                  bbox_to_anchor=(0.5, 1.0 + 0.9 / BOX_H), borderaxespad=0)
    else:
        pad = 0.06 * (hi - lo or 1); ax.set_xlim(lo - pad, hi + pad)
        ax.set_xlabel('Sorting (phi)')
    ax.xaxis.tick_top(); ax.xaxis.set_label_position('top')
    for sp in ('right', 'bottom'): ax.spines[sp].set_visible(False)

def labels(ax):
    for e, t, b in EV:
        ax.text(1.03, (max(t, DEPTH_TOP) + min(b, DEPTH_BOT)) / 2, f'E{e:02d}',
                color=DARK[e], fontsize=6.5, fontweight='bold', va='center',
                transform=ax.get_yaxis_transform(), clip_on=False)

def page(panels, name):
    W = MARGIN_CM['left'] + sum(w for _, w, _ in panels) + GAP_CM * (len(panels) - 1) + MARGIN_CM['right']
    H = MARGIN_CM['top'] + BOX_H + MARGIN_CM['bottom']
    fig = plt.figure(figsize=(W * CM, H * CM))
    x, axes = MARGIN_CM['left'], []
    for i, (kind, w, arg) in enumerate(panels):
        ax = fig.add_axes([x / W, MARGIN_CM['bottom'] / H, w / W, BOX_H / H])
        if kind == 'photo': draw_photo(ax, arg)
        else: draw_data(ax, arg)
        depth_axis(ax, labels=(i == 0)); axes.append(ax); x += w + GAP_CM
    labels(axes[-1])
    fig.text(0.25 / W, 1 - 0.25 / H, CORE, fontsize=8, fontweight='bold', va='top')
    fig.savefig(f'{OUT}/{name}.pdf', facecolor='white')     # exact size, no cropping
    fig.savefig(f'{OUT}/{name}_preview.png', dpi=40, facecolor='white')
    plt.close(fig)
    print(f"   saved {name}.pdf  ({W:.1f} x {H:.1f} cm, data box {BOX_H:.1f} cm tall)")

# ---- 3a. DATA figure first (no photo) ----------------------------------------
print("\n3. figures")
page([('data', PANEL_W_CM, ['D10', 'D50', 'D90']), ('data', SORT_W_CM, ['Sorting'])],
     f'{CORE}_DATA_realscale')

# ---- 3b. photo figure ---------------------------------------------------------
def load_photo():
    p = find(PHOTO_FILE, ('.jpg', '.jpeg', '.png', '.tif', '.tiff'))
    if p is None or 'realscale' in p:
        raise FileNotFoundError(f"no photo with '{PHOTO_FILE}' in its name in {FOLDER}")
    from PIL import Image
    Image.MAX_IMAGE_PIXELS = None
    img = Image.open(p)
    per_cm = img.size[1] / (PHOTO_BOT_CM - PHOTO_TOP_CM)
    if per_cm > PHOTO_PX_PER_CM:                          # shrink big scans first
        k = PHOTO_PX_PER_CM / per_cm
        img.draft('RGB', (int(img.size[0] * k), int(img.size[1] * k)))
        img = img.convert('RGB')
        k = PHOTO_PX_PER_CM / (img.size[1] / (PHOTO_BOT_CM - PHOTO_TOP_CM))
        if k < 0.95:
            img = img.resize((round(img.size[0] * k), round(img.size[1] * k)), Image.LANCZOS)
    img = img.convert('RGB')
    w, h = img.size
    per_cm = h / (PHOTO_BOT_CM - PHOTO_TOP_CM)
    r0 = int(round((DEPTH_TOP - PHOTO_TOP_CM) * per_cm))
    r1 = int(round((DEPTH_BOT - PHOTO_TOP_CM) * per_cm))
    crop = np.asarray(img.crop((0, max(r0, 0), w, min(r1, h))))
    print(f"   photo {os.path.basename(p)}: {w} x {h} px used, {per_cm:.1f} px/cm, "
          f"true width {w / per_cm * SCALE:.2f} cm")
    return dict(img=crop, top=PHOTO_TOP_CM + max(r0, 0) / per_cm,
                bot=PHOTO_TOP_CM + min(r1, h) / per_cm, width_cm=w / per_cm * SCALE)

def draw_photo(ax, ph):
    ax.imshow(ph['img'], extent=[0, 1, ph['bot'], ph['top']], aspect='auto',
              interpolation='none')
    ax.set_xlim(0, 1); ax.set_xticks([])
    xb = 1 + 0.15 / ph['width_cm']
    for e, t, b in EV:
        ax.plot([xb, xb], [max(t, DEPTH_TOP), min(b, DEPTH_BOT)], color=COLOR[e], lw=3,
                solid_capstyle='butt', clip_on=False, zorder=4)
    for y in CONTACTS:
        ax.plot([xb - 0.08 / ph['width_cm'], xb + 0.08 / ph['width_cm']], [y, y],
                color='k', lw=1.0, clip_on=False, zorder=5)

try:
    ph = load_photo()
    page([('photo', ph['width_cm'], ph), ('data', PANEL_W_CM, ['D10', 'D50', 'D90']),
          ('data', SORT_W_CM, ['Sorting'])], f'{CORE}_realscale')
except Exception as ex:
    print(f"   !!! photo figure NOT made: {type(ex).__name__}: {ex}")
    print(f"       (the DATA figure above was saved)")

# ----------------------------------------------------------------------------
# 4. show + download
# ----------------------------------------------------------------------------
try:
    from IPython.display import Image as Show, display
    for f in sorted(glob.glob(f'{OUT}/*_preview.png')):
        print(os.path.basename(f)); display(Show(f, width=420))
except Exception:
    pass
zip_path = f'/content/{CORE}_realscale_PDFs.zip'
with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as z:
    for f in sorted(glob.glob(f'{OUT}/*.pdf')):
        z.write(f, os.path.basename(f))
print(f"\n4. {zip_path}: {len(glob.glob(OUT + '/*.pdf'))} PDF(s)")
try:
    from google.colab import files
    files.download(zip_path)
except Exception:
    pass
