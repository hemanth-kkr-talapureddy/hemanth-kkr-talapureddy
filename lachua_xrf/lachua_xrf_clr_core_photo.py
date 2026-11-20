# =============================================================================
# Lachua (GUAC-29A, GUAC-24A) - core photo + clr of every XRF element
#
# INPUT FILES (nothing is typed into this code):
#   1. GUAC-29A-1G-1_RepoOG.xlsx            XRF counts, 29A
#   2. GUAC-24A-1G-1_RepoOG.xlsx            XRF counts, 24A
#   3. Lachua_Graine_Siza_Samples_N_Events_with_Ages.xlsx
#                                           event depths, ages, uncertainties
#   4. Lachua_Varve_counts.xlsx             top and bottom depth of each core
#   5. core photo of 29A (.jpg / .png / .tif), with its ruler, as scanned
#   6. core photo of 24A
#   + for each photo, the ruler depth at its first and last pixel: PHOTOS below
#
# OUTPUT, for each core:
#   <core>_photo_clr_events : photo | 10 clr panels, with event bands + ages
#   <core>_photo_clr_lines  : photo | 10 clr panels, lines only
#   <core>_photo_check      : the whole photo with the depth lines the code
#                             uses drawn over it - they must sit on the ruler
#   saved as 300 dpi PNG + vector PDF
#
# THE RULES
#   * clr over the same 10 elements as the event figures, so values match them
#   * every XRF sample in the core is drawn: ONE continuous line per element
#     from the top to the bottom of the core (no break at background or varves)
#     with a dot at each sample (one XRF measurement every 5 mm)
#   * the photo is turned upright if it was scanned lying down, placed by its
#     ruler, and cut to the sediment (ruler and colour card removed)
#   * depth axis of each core: first -> last 'Depth (in core)' in the varve file
#   * ONE common depth scale, same panel widths, same x ranges, 24 pt text
#
# THE STEPS (each is one function below, in this order):
#   STEP 1  event file                  -> read_events()
#   STEP 2  varve file: depth range     -> read_core_range()
#   STEP 3  XRF: clr                    -> read_xrf()
#   STEP 4  core photo                  -> read_photo()
#   STEP 5  draw the figure             -> draw_core()
#   STEP 6  check the drawn figure      -> check_figure()   (stops on any error)
#
# Reference
#   clr: Bertrand et al. (2024) Earth-Sci. Rev. 249, 104639, Eqs. 1-2
# =============================================================================
import os, re, glob
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, MultipleLocator, AutoMinorLocator
from PIL import Image

Image.MAX_IMAGE_PIXELS = None    # core photos are often very tall scans

# =============================================================================
# SETTINGS - the only lines you should need to change
# =============================================================================
DATA_DIR = '/content'            # where the input files are (Colab default)
OUT_DIR  = '/content/figures'    # where the figures are written

# Each file is found by WORDS in its name (case, '-', '_', spaces ignored).
FILE_WORDS = {
    'events':   ['Events', 'Ages'],
    'varves':   ['Varve', 'counts'],
    'GUAC-29A': dict(xrf=['29A', 'RepoOG']),
    'GUAC-24A': dict(xrf=['24A', 'RepoOG']),
}
CORE_IDS = {'GUAC-29A': 'GUAC-29A-1G-1',      # core names as written in the
            'GUAC-24A': 'GUAC-24A-1G-1'}      # varve file

# ---- CORE PHOTOS -------------------------------------------------------------
# words    finds the photo file by name.
# rotate   'auto' : a photo wider than tall is turned upright so that its LEFT
#                   end becomes the top (your 29A photo: core top on the left)
#          'none' : use the photo as it is
#          'left_top' / 'right_top' : turn it, with that end at the top
# first_cm, last_cm : ruler depth (cm) at the photo's first and last pixel
#          along the core (after turning). Measured on your photos from the
#          ruler's white/grey 10 cm blocks (10, 20 ... 100 cm), with a straight
#          line fitted through all ten block edges:
#            29A: 19.90 px per cm, worst misfit 0.04 cm -> 0.00 to 100.52 cm
#            24A: 18.95 px per cm, worst misfit 0.05 cm -> 0.00 to 105.54 cm
#          These depths are fractions of the photo length, so they hold for a
#          larger copy of the SAME photo. Check them on <core>_photo_check.png.
# across   part of the photo's WIDTH to keep, as fractions (0 = the side that
#          is left after turning, 1 = the right side). Keeps the sediment and
#          drops the ruler strip and the core edges. Measured on your photos.
PHOTOS = {
    'GUAC-29A': dict(words=['29A'], rotate='auto', first_cm=0.00, last_cm=100.52,
                     across=(0.22, 0.84)),
    'GUAC-24A': dict(words=['24A'], rotate='auto', first_cm=0.00, last_cm=105.54,
                     across=(0.22, 0.96)),
}
PHOTO_EXT = ('.jpg', '.jpeg', '.png', '.tif', '.tiff')
PHOTO_MAX_ROWS = 6000            # very tall scans are reduced to this many rows

# Column names, exactly as written in your files.
XRF_DEPTH_COL   = 'Section Depth (cm)'
VARVE_DEPTH_COL = 'Depth (in core)'

# The 10 clr elements, in panel order, grouped as in your v2 notebook
# (Bertrand et al. 2024): redox | lithogenic | authigenic.
CLR_GROUPS = [
    ('Redox',                 ['Mn']),
    ('Lithogenic (detrital)', ['Al', 'Si', 'K', 'Ti', 'Fe', 'Zr']),
    ('Authigenic',            ['Ca', 'Sr', 'S']),
]
CLR_ELEMENTS = [el for _, els in CLR_GROUPS for el in els]
# Element colours from your v2 notebook (ELCOL).
ELEMENT_COLORS = {'Si': '#4299e1', 'K': '#dd6b20', 'Ca': '#2b6cb0', 'Ti': '#2c7a7b',
                  'Mn': '#38a169', 'Fe': '#c05621', 'Sr': '#805ad5', 'Zr': '#b83280',
                  'S': '#975a16', 'Al': '#5a67d8'}

# ---- SIZE --------------------------------------------------------------------
FONT_PT       = 24               # ALL text and scale numbers
DEPTH_SCALE   = 1.0              # cm on paper per 1 cm of core (same both cores)
PHOTO_W_CM    = 5.0              # width of the photo column (same both cores)
PANEL_W_CM    = 6.5              # width of every clr panel (same both cores)
PANEL_GAP_CM  = 1.2              # gap between panels
LABEL_W_CM    = 11.0             # column for "E01  2010.5 ± 1.7" (events version)
LEFT_CM       = 4.5              # room for the depth numbers + "Depth (cm)"
TOP_CM        = 8.5              # room for title, group names, panel names
BOTTOM_CM     = 1.5
DEPTH_MAJOR_CM = 5               # labelled depth tick every 5 cm
DEPTH_MINOR_CM = 0.5             # small depth tick every 5 mm (XRF step)
DPI = 300

# Event colours from your event colour table (8 colours that repeat).
_TABLE = ['#f6a362', '#2a9f8f', '#ebc56b', '#6e8fae',
          '#e97053', '#8bb37d', '#b5848f', '#9d8bba']
EVENT_COLORS = {e: _TABLE[(e - 1) % 8] for e in range(1, 21)}
BAND_ALPHA = 0.30

# =============================================================================
# STYLE - every font size is FONT_PT
# =============================================================================
matplotlib.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'Liberation Sans', 'DejaVu Sans'],
    'font.size': FONT_PT, 'axes.labelsize': FONT_PT, 'axes.titlesize': FONT_PT,
    'xtick.labelsize': FONT_PT, 'ytick.labelsize': FONT_PT,
    'axes.linewidth': 2.0,
    'xtick.major.width': 2.0, 'ytick.major.width': 2.0,
    'xtick.minor.width': 1.4, 'ytick.minor.width': 1.4,
    'xtick.major.size': 10, 'ytick.major.size': 10,
    'xtick.minor.size': 5, 'ytick.minor.size': 5,
    'xtick.major.pad': 6, 'ytick.major.pad': 6,
    'xtick.direction': 'out', 'ytick.direction': 'out',
    'pdf.fonttype': 42, 'ps.fonttype': 42,      # PDF text stays editable
    'axes.grid': False,
})
CM = 1 / 2.54                                   # cm -> inch
LINE_W = 2.0                                    # data line width (pt)
DOT_SIZE = 5.0                                  # sample dot diameter (pt)


# =============================================================================
# FINDING AND READING FILES
# =============================================================================
def _norm(s):
    """'GUAC-29A_GSD' -> 'guac29agsd' so names can be compared loosely."""
    return re.sub(r'[^a-z0-9]', '', str(s).lower())


def find_file(words, ext=('.csv', '.xlsx', '.xls')):
    """The newest file in DATA_DIR with one of these extensions whose name
    contains every word."""
    hits = [f for f in glob.glob(os.path.join(DATA_DIR, '*'))
            if f.lower().endswith(ext)
            and all(_norm(w) in _norm(os.path.basename(f)) for w in words)]
    return max(hits, key=os.path.getmtime) if hits else None


def needed():
    """(description, words, extensions) of every input file."""
    out = [('event table', FILE_WORDS['events'], ('.xlsx', '.xls', '.csv')),
           ('varve counts', FILE_WORDS['varves'], ('.xlsx', '.xls', '.csv'))]
    for core in CORE_IDS:
        out.append((f'{core} XRF', FILE_WORDS[core]['xrf'], ('.xlsx', '.xls', '.csv')))
        out.append((f'{core} photo', PHOTOS[core]['words'], PHOTO_EXT))
    return out


def ensure_files():
    """Stop early, with a clear message, if a photo depth or a file is missing."""
    for core, p in PHOTOS.items():
        if p.get('first_cm') is None or p.get('last_cm') is None:
            raise ValueError(f'Fill in first_cm and last_cm for {core} in PHOTOS '
                             f'(ruler depth at the photo\'s first and last pixel).')
        if not p['last_cm'] > p['first_cm']:
            raise ValueError(f'{core}: photo last_cm must be deeper than first_cm')
        a, b = p.get('across', (0, 1))
        if not 0 <= a < b <= 1:
            raise ValueError(f'{core}: photo across must be two fractions 0 <= a < b <= 1')
        if p.get('rotate', 'auto') not in ('auto', 'none', 'left_top', 'right_top'):
            raise ValueError(f'{core}: photo rotate must be auto, none, left_top or right_top')
    if all(find_file(w, e) for _, w, e in needed()):
        return
    try:
        from google.colab import files
        print('Upload: 2 XRF (.xlsx), the event table, the varve counts, 2 core photos')
        for name, data in files.upload().items():
            with open(os.path.join(DATA_DIR, name), 'wb') as fh:
                fh.write(data)
    except ImportError:
        pass
    missing = [d for d, w, e in needed() if not find_file(w, e)]
    if missing:
        raise FileNotFoundError(f'Missing in {DATA_DIR}: {missing}')


def read_table(path, **kw):
    return pd.read_excel(path, **kw) if path.lower().endswith(('.xlsx', '.xls')) \
        else pd.read_csv(path, encoding='utf-8-sig', **kw)


# =============================================================================
# STEP 1 - EVENTS, AGES AND UNCERTAINTIES FROM YOUR EVENT FILE
# =============================================================================
def read_events(path):
    """One row per event, one column block per core:
         age_CE | age_unc | event | in_24A | 24A_top_cm | 24A_base_cm | ...
    Events marked ABSENT (or with an empty depth) are not in that core."""
    df = read_table(path)
    need = ['age_CE', 'age_unc', 'event']
    for tag in ('24A', '29A'):
        need += [f'in_{tag}', f'{tag}_top_cm', f'{tag}_base_cm']
    missing = [c for c in need if c not in df.columns]
    if missing:
        raise KeyError(f'Event file is missing columns {missing}. Found {list(df.columns)}')
    events = {}
    for core in CORE_IDS:
        tag = core.split('-')[1]                        # 'GUAC-24A' -> '24A'
        rows = []
        for _, r in df.iterrows():
            top, base = r[f'{tag}_top_cm'], r[f'{tag}_base_cm']
            if str(r[f'in_{tag}']).strip().lower() != 'yes' or pd.isna(top) or pd.isna(base):
                continue
            rows.append(dict(event=int(re.search(r'\d+', str(r['event'])).group()),
                             top=float(top), base=float(base),
                             age=float(r['age_CE']), unc=float(r['age_unc'])))
        events[core] = pd.DataFrame(rows).sort_values('top').reset_index(drop=True)
        print(f'{core}: {len(rows)} events in the event file')
    return events


# =============================================================================
# STEP 2 - DEPTH RANGE OF EACH CORE FROM THE VARVE FILE
# =============================================================================
def read_core_range(path, core_id):
    """The varve sheet has one block of columns per core:
         row 1:  GUAC-24A-1G-1 |                 |
         row 2:  Year          | Depth (in core) | Depth (adjusted)
    Returns (first, last) 'Depth (in core)' of this core."""
    raw = read_table(path, header=None)
    hit = [(r, c) for r in range(min(5, len(raw))) for c in range(raw.shape[1])
           if str(raw.iat[r, c]).strip() == core_id]
    if not hit:
        raise KeyError(f'{core_id} not found in {os.path.basename(path)}')
    r0, c0 = hit[0]
    header = [str(v).strip() for v in raw.iloc[r0 + 1, c0:c0 + 4]]
    if VARVE_DEPTH_COL not in header:
        raise KeyError(f'No "{VARVE_DEPTH_COL}" column under {core_id}. Found {header}')
    col = c0 + header.index(VARVE_DEPTH_COL)
    depth = pd.to_numeric(raw.iloc[r0 + 2:, col], errors='coerce').dropna()
    if not depth.is_monotonic_increasing:
        raise ValueError(f'{core_id}: varve depths are not in increasing order')
    return float(depth.iloc[0]), float(depth.iloc[-1])


# =============================================================================
# STEP 3 - XRF: clr
# =============================================================================
def clr(counts):
    """clr(x) = ln(x) - mean(ln(x)) for each sample (row)."""
    L = np.log(counts)
    return L - L.mean(axis=1, keepdims=True)


def read_xrf(path):
    df = read_table(path)
    df.columns = [str(c).strip() for c in df.columns]
    missing = [c for c in [XRF_DEPTH_COL] + CLR_ELEMENTS if c not in df.columns]
    if missing:
        raise KeyError(f'{os.path.basename(path)} is missing columns {missing}')
    d = df[[XRF_DEPTH_COL] + CLR_ELEMENTS].apply(pd.to_numeric, errors='coerce')
    d = d.rename(columns={XRF_DEPTH_COL: 'depth'}).dropna()
    d = d.groupby('depth', as_index=False).mean()       # repeat scans -> mean
    # ln(0) does not exist: a zero count is replaced by half the smallest
    # non-zero count of that element (Bertrand et al. 2024). It is printed.
    for e in CLR_ELEMENTS:
        zero = d[e] <= 0
        if zero.any():
            half = d.loc[~zero, e].min() / 2
            print(f'   {e}: {zero.sum()} zero count(s) replaced by {half:g}')
            d.loc[zero, e] = half
    C = clr(d[CLR_ELEMENTS].to_numpy(float))
    assert np.abs(C.sum(axis=1)).max() < 1e-9          # clr rows must sum to 0
    for i, e in enumerate(CLR_ELEMENTS):
        d[f'{e}_clr'] = C[:, i]
    return d.sort_values('depth').reset_index(drop=True)


# =============================================================================
# STEP 4 - CORE PHOTO
# =============================================================================
def read_photo(path, rotate='auto', across=(0.0, 1.0)):
    """The photo as an RGB array with depth running DOWN its rows.
    1. a photo lying on its side is turned upright (see `rotate` in PHOTOS);
    2. only the `across` part of its width is kept (the sediment);
    3. very tall scans are reduced to PHOTO_MAX_ROWS rows.
    Returns (photo, full upright photo for the check image, what was done)."""
    img = Image.open(path).convert('RGB')
    w, h = img.size
    done = []
    if rotate == 'auto':
        rotate = 'left_top' if w > h else 'none'
    if rotate == 'left_top':                 # left end -> top: turn clockwise
        img = img.transpose(Image.ROTATE_270)
        done.append('turned upright (left end at the top)')
    elif rotate == 'right_top':              # right end -> top: turn anticlockwise
        img = img.transpose(Image.ROTATE_90)
        done.append('turned upright (right end at the top)')
    full = img
    w, h = img.size
    x0, x1 = round(across[0] * w), round(across[1] * w)
    img = img.crop((x0, 0, x1, h))
    done.append(f'width cut to {across[0]:.2f}-{across[1]:.2f} (pixels {x0}-{x1} of {w})')
    if h > PHOTO_MAX_ROWS:
        f = PHOTO_MAX_ROWS / h
        img = img.resize((max(1, round(img.size[0] * f)), PHOTO_MAX_ROWS), Image.LANCZOS)
        full = full.resize((max(1, round(full.size[0] * f)), PHOTO_MAX_ROWS), Image.LANCZOS)
    return np.asarray(img), np.asarray(full), '; '.join(done)


def save_photo_check(core, full, photo_depths, across, depth_range):
    """The whole upright photo, with a red line every 10 cm at the depth the code
    believes it is, and the kept width and the plotted depth range outlined.
    If the red lines sit on the ruler's 10 cm marks, the photo is placed right."""
    top, bottom = photo_depths
    h_px, w_px = full.shape[:2]
    fig_h = 40.0
    fig_w = max(6.0, fig_h * w_px / h_px) + 4.0
    fig = plt.figure(figsize=(fig_w * CM, (fig_h + 3) * CM))
    ax = fig.add_axes([3.0 / fig_w, 1.0 / (fig_h + 3), (fig_w - 4.0) / fig_w,
                       fig_h / (fig_h + 3)])
    ax.imshow(full, extent=[0, 1, bottom, top], aspect='auto')
    for d in np.arange(np.ceil(top / 10) * 10, bottom + 1e-9, 10):
        ax.axhline(d, color='red', lw=1.0)
    ax.axvspan(across[0], across[1], ymin=0, ymax=1, fill=False, ec='cyan', lw=2)
    for y in depth_range:
        ax.axhline(y, color='cyan', lw=1.5, ls='--')
    ax.set_ylim(bottom, top)
    ax.set_xticks([])
    ax.yaxis.set_major_locator(MultipleLocator(10))
    ax.tick_params(labelsize=10)
    ax.set_title(f'{core}\nred: every 10 cm\ncyan: part used', fontsize=10)
    name = os.path.join(OUT_DIR, f'{core}_photo_check.png')
    fig.savefig(name, dpi=200)
    plt.close(fig)
    return name


# =============================================================================
# STEP 5 - DRAWING
# =============================================================================
def in_range(depth, top, bottom):
    return (depth >= top - 1e-9) & (depth <= bottom + 1e-9)


def darker(color, f=0.72):
    r, g, b = matplotlib.colors.to_rgb(color)
    return (r * f, g * f, b * f)


def layout(with_events, depth_top, depth_bottom):
    """Figure size and panel positions in cm. Every width and the depth scale
    are fixed numbers, so they are identical in both cores."""
    n = len(CLR_ELEMENTS)
    plot_h = (depth_bottom - depth_top) * DEPTH_SCALE
    photo_x = LEFT_CM
    lefts = [photo_x + PHOTO_W_CM + PANEL_GAP_CM + i * (PANEL_W_CM + PANEL_GAP_CM)
             for i in range(n)]
    fig_w = lefts[-1] + PANEL_W_CM + (LABEL_W_CM if with_events else 1.0)
    fig_h = TOP_CM + plot_h + BOTTOM_CM
    return fig_w, fig_h, plot_h, photo_x, lefts


def draw_core(core, xrf, photo, photo_depths, ev, xlims, depth_range, with_events):
    top, bottom = depth_range
    fig_w, fig_h, plot_h, photo_x, lefts = layout(with_events, top, bottom)
    fig = plt.figure(figsize=(fig_w * CM, fig_h * CM))
    frac = lambda x, y, w, h: [x / fig_w, y / fig_h, w / fig_w, h / fig_h]

    # ---- core photo: placed by the depths of its top and bottom edge -------
    pax = fig.add_axes(frac(photo_x, BOTTOM_CM, PHOTO_W_CM, plot_h))
    pax.imshow(photo, extent=[0, 1, photo_depths[1], photo_depths[0]],
               aspect='auto', interpolation='lanczos')
    pax.set_ylim(bottom, top)                       # depth increases downwards
    pax.set_xlim(0, 1)
    pax.set_xticks([])
    pax.xaxis.set_label_position('top')
    pax.set_xlabel('Core\nphoto', labelpad=12, fontweight='bold')
    pax.spines[['top', 'right', 'bottom']].set_visible(False)
    pax.yaxis.set_major_locator(MultipleLocator(DEPTH_MAJOR_CM))
    pax.yaxis.set_minor_locator(MultipleLocator(DEPTH_MINOR_CM))
    pax.set_ylabel('Depth (cm)', labelpad=10)

    # every XRF sample inside the core's depth range, as ONE continuous line
    d = xrf[in_range(xrf.depth, top, bottom)]
    axes = []
    for i, el in enumerate(CLR_ELEMENTS):
        ax = fig.add_axes(frac(lefts[i], BOTTOM_CM, PANEL_W_CM, plot_h), sharey=pax)
        axes.append(ax)
        if with_events:
            for _, e in ev.iterrows():
                ax.axhspan(e.top, e.base, color=EVENT_COLORS[int(e.event)],
                           alpha=BAND_ALPHA, lw=0, zorder=0)
            for y in sorted(set(np.round(np.r_[ev.top, ev.base], 3))):
                ax.axhline(y, color='0.3', lw=1.2, ls=(0, (5, 3)), zorder=1)
        # one continuous line with a dot at every sample (5 mm XRF step)
        ax.plot(d[f'{el}_clr'], d.depth, '-', color=ELEMENT_COLORS[el], lw=LINE_W,
                marker='o', ms=DOT_SIZE, mfc=ELEMENT_COLORS[el], mec='white',
                mew=0.6, zorder=3, gid='data')

        lo, hi = xlims[el]                          # same range in both cores
        pad = 0.08 * (hi - lo)
        ax.set_xlim(lo - pad, hi + pad)
        ax.xaxis.set_major_locator(MaxNLocator(3, steps=[1, 2, 2.5, 5, 10]))
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position('top')
        ax.set_xlabel(f'{el}\n(clr)', labelpad=12, fontweight='bold')
        ax.spines[['bottom', 'right', 'left']].set_visible(False)
        ax.tick_params(axis='y', which='both', left=False, labelleft=False)

    # group names with a line above the panels
    fig.canvas.draw()
    label_tops = max(a.xaxis.label.get_window_extent().y1 for a in axes + [pax])
    y_line = label_tops / fig.dpi / CM + 0.6
    k = 0
    for name, els in CLR_GROUPS:
        group_axes = axes[k:k + len(els)]
        k += len(els)
        x0 = group_axes[0].get_position().x0 + 0.002
        x1 = group_axes[-1].get_position().x1 - 0.002
        fig.add_artist(plt.Line2D([x0, x1], [y_line / fig_h] * 2, color='k', lw=2,
                                  transform=fig.transFigure))
        fig.text((x0 + x1) / 2, (y_line + 0.3) / fig_h, name, ha='center',
                 va='bottom', fontweight='bold')

    # event labels (events version only)
    if with_events:
        label_x0 = lefts[-1] + PANEL_W_CM + 0.6
        lab = fig.add_axes(frac(label_x0, BOTTOM_CM, LABEL_W_CM - 0.6, plot_h),
                           sharey=pax)
        lab.axis('off')
        min_gap = 1.15 * FONT_PT / 72 * 2.54 / DEPTH_SCALE
        y_prev = -np.inf
        tr = lab.get_yaxis_transform()
        for _, e in ev.iterrows():
            mid = (e.top + e.base) / 2
            y = max(mid, y_prev + min_gap)
            y_prev = y
            if abs(y - mid) > 0.05:
                lab.plot([0.0, 0.05], [mid, y], color='0.4', lw=1.2,
                         transform=tr, clip_on=False)
            lab.text(0.07, y, f'E{int(e.event):02d}  {e.age:.1f} ± {e.unc:.1f}',
                     color=darker(EVENT_COLORS[int(e.event)]), fontweight='bold',
                     va='center', ha='left', transform=tr, clip_on=False,
                     gid=f'E{int(e.event)}')
        lab.text(0.07, top, 'Event  Age (CE)', fontweight='bold', va='bottom',
                 ha='left', transform=tr)

    fig.text(0.5 / fig_w, 1 - 0.4 / fig_h, core.replace('-', ' '),
             fontweight='bold', ha='left', va='top')
    return fig, pax, axes


# =============================================================================
# STEP 6 - CHECK THE DRAWN FIGURE AGAINST THE DATA (stops on any error)
# =============================================================================
def check_figure(core, fig, pax, axes, xrf, photo_depths, ev, depth_range, with_events):
    problems = []
    top, bottom = depth_range
    for t in fig.findobj(matplotlib.text.Text):
        if t.get_text().strip() and t.get_visible() and abs(t.get_fontsize() - FONT_PT) > 1e-6:
            problems.append(f'text "{t.get_text()[:20]}" is {t.get_fontsize()} pt')
    fig_w, fig_h = fig.get_size_inches() / CM
    if abs(pax.get_position().width * fig_w - PHOTO_W_CM) > 1e-6:
        problems.append('photo column width')
    img = pax.get_images()[0]
    if img.get_array().shape[0] <= img.get_array().shape[1]:
        problems.append('photo is not upright (wider than tall)')
    ext = img.get_extent()                           # [left, right, bottom, top]
    if (ext[3], ext[2]) != tuple(photo_depths):
        problems.append(f'photo placed at {ext[3]}-{ext[2]} cm, not {photo_depths}')
    want = xrf[in_range(xrf.depth, top, bottom)]
    for ax, el in zip(axes, CLR_ELEMENTS):
        pos = ax.get_position()
        if abs(pos.width * fig_w - PANEL_W_CM) > 1e-6:
            problems.append(f'{el}: panel width {pos.width * fig_w:.3f} cm')
        scale = pos.height * fig_h / (bottom - top)
        if abs(scale - DEPTH_SCALE) > 1e-6:
            problems.append(f'{el}: depth scale {scale:.4f}')
        lines = [l for l in ax.get_lines() if l.get_gid() == 'data']
        # exactly one line, holding every sample in depth order (continuous)
        if len(lines) != 1:
            problems.append(f'{el}: {len(lines)} data lines, expected 1 continuous line')
            continue
        if lines[0].get_marker() != 'o':
            problems.append(f'{el}: samples are not drawn as dots')
        got = np.c_[lines[0].get_ydata(), lines[0].get_xdata()].astype(float)
        if not np.allclose(got, want[['depth', f'{el}_clr']].to_numpy(float), atol=1e-12):
            problems.append(f'{el}: drawn line differs from the data')
    if tuple(pax.get_ylim()) != (bottom, top):
        problems.append(f'depth axis {pax.get_ylim()} != varve range {depth_range}')
    if with_events:
        labels = {t.get_gid() for t in fig.findobj(matplotlib.text.Text) if t.get_gid()}
        missing = [f'E{e}' for e in ev.event if f'E{e}' not in labels]
        if missing:
            problems.append(f'no label for {missing}')
    if problems:
        raise AssertionError(f'{core}: ' + '; '.join(problems))
    print(f'   CHECK PASSED: {len(want)} XRF samples in one continuous line per '
          f'element, photo at {photo_depths[0]:g}-{photo_depths[1]:g} cm, all text '
          f'{FONT_PT} pt, panels {PANEL_W_CM} cm, depth scale {DEPTH_SCALE} cm/cm'
          + (f', {len(ev)} events labelled' if with_events else ''))


def _in_notebook():
    try:
        from IPython import get_ipython
        return get_ipython() is not None
    except ImportError:
        return False


# =============================================================================
# RUN EVERYTHING
# =============================================================================
def main():
    ensure_files()
    events = read_events(find_file(FILE_WORDS['events']))           # STEP 1
    varve_path = find_file(FILE_WORDS['varves'])

    cores = {}
    for core in CORE_IDS:
        xp = find_file(FILE_WORDS[core]['xrf'])
        pp = find_file(PHOTOS[core]['words'], PHOTO_EXT)
        rng = read_core_range(varve_path, CORE_IDS[core])           # STEP 2
        xrf = read_xrf(xp)                                          # STEP 3
        ph = PHOTOS[core]
        photo, full, done = read_photo(pp, ph.get('rotate', 'auto'),
                                       ph.get('across', (0.0, 1.0)))  # STEP 4
        pdep = (float(ph['first_cm']), float(ph['last_cm']))
        n_in = int(in_range(xrf.depth, *rng).sum())
        print(f'\n{core}\n   XRF   : {os.path.basename(xp)}'
              f'\n   photo : {os.path.basename(pp)} ({photo.shape[1]} x {photo.shape[0]} px),'
              f' placed at {pdep[0]:g}-{pdep[1]:g} cm; {done}'
              f'\n   depth axis : {rng[0]:g}-{rng[1]:g} cm (varve file); '
              f'{n_in} of {len(xrf)} XRF samples fall inside it')
        if pdep[0] > rng[0] or pdep[1] < rng[1]:
            print(f'   ! the photo ({pdep[0]:g}-{pdep[1]:g} cm) does not cover the '
                  f'whole depth axis ({rng[0]:g}-{rng[1]:g} cm): that part stays blank')
        step = np.median(np.diff(xrf.depth))
        print(f'   XRF step   : {step * 10:g} mm (median spacing of the samples)')
        os.makedirs(OUT_DIR, exist_ok=True)
        chk = save_photo_check(core, full, pdep, ph.get('across', (0.0, 1.0)), rng)
        print(f'   photo check: {os.path.basename(chk)} (red lines must sit on the '
              f'ruler\'s 10 cm marks)')
        cores[core] = dict(xrf=xrf, photo=photo, pdep=pdep, ev=events[core], rng=rng)

    # one x range per element, shared by both cores
    xlims = {}
    for el in CLR_ELEMENTS:
        v = pd.concat([c['xrf'].loc[in_range(c['xrf'].depth, *c['rng']), f'{el}_clr']
                       for c in cores.values()])
        xlims[el] = (v.min(), v.max())

    os.makedirs(OUT_DIR, exist_ok=True)
    for with_events, tag in ((True, 'photo_clr_events'), (False, 'photo_clr_lines')):
        for core, c in cores.items():
            fig, pax, axes = draw_core(core, c['xrf'], c['photo'], c['pdep'], c['ev'],
                                       xlims, c['rng'], with_events)        # STEP 5
            check_figure(core, fig, pax, axes, c['xrf'], c['pdep'], c['ev'],
                         c['rng'], with_events)                            # STEP 6
            name = f'{core}_{tag}'
            for ext in ('png', 'pdf'):
                fig.savefig(os.path.join(OUT_DIR, f'{name}.{ext}'), dpi=DPI)
            w, h = fig.get_size_inches() / CM
            print(f'   saved {name}.png / .pdf  ({w:.1f} x {h:.1f} cm)')
            if _in_notebook():
                plt.show()
            plt.close(fig)
    print(f'\nAll figures are in {OUT_DIR}')


if __name__ == '__main__':
    main()
