# =============================================================================
# Lachua (GUAC-29A, GUAC-24A) - core photo + mean grain size + clr of every
# XRF element + carbon and nitrogen (29A)
#
# INPUT FILES (nothing is typed into this code):
#   1. GUAC-29A-1G-1_RepoOG.xlsx            XRF counts, 29A
#   2. GUAC-24A-1G-1_RepoOG.xlsx            XRF counts, 24A
#   3. GUAC_29A_..._FolkWard_Stats.csv      grain-size statistics, 29A
#   4. GUAC_24A_..._FolkWard_Stats_1.csv    grain-size statistics, 24A
#   5. Obrist_Lachua_July_2026_All_results.xlsx
#                                           carbon + nitrogen (29A samples)
#   6. Lachua_Graine_Siza_Samples_N_Events_with_Ages.xlsx
#                                           event depths, ages, uncertainties
#   7. Lachua_Varve_counts.xlsx             top and bottom depth of each core
#   8. core photo of 29A (.jpg / .png / .tif), with its ruler, as scanned
#   9. core photo of 24A
#
# OUTPUT, for each core (photo | grain size | 10 clr | carbon, if the core has
# carbon data), in two line styles x two backgrounds:
#   <core>_V7_continuous_events : one line through all samples, event bands
#   <core>_V7_continuous_lines  : one line through all samples
#   <core>_V7_gaps_events       : lines stop where samples are missing, bands
#   <core>_V7_gaps_lines        : lines stop where samples are missing
#   <core>_photo_check          : the photo with the depth lines the code uses
#   saved as 300 dpi PNG + vector PDF
#
# THE RULES
#   * a dot at every sample (XRF every 5 mm; grain size and carbon as sampled)
#   * 'continuous': ONE line per panel from the first to the last sample
#     'gaps'      : the line stops where a sample is missing (spacing more
#                   than GAP_FACTOR x the usual spacing), so no line crosses
#                   unsampled core
#   * clr over the same 10 elements as the event figures
#   * C/N is the ATOMIC ratio: (TOC / 12.011) / (N / 14.007)   (Meyers 1994)
#   * depth axis of each core: first -> last 'Depth (in core)' in the varve file
#   * ONE common depth scale, same panel widths, same x ranges in both cores,
#     24 pt text
#
# THE STEPS (each is one function below, in this order):
#   STEP 1  event file                  -> read_events()
#   STEP 2  varve file: depth range     -> read_core_range()
#   STEP 3  XRF: clr                    -> read_xrf()
#   STEP 4  grain size                  -> read_grain_size()
#   STEP 5  carbon + nitrogen           -> read_carbon()
#   STEP 6  core photo                  -> read_photo()
#   STEP 7  draw the figure             -> draw_core()
#   STEP 8  check the drawn figure      -> check_figure()   (stops on any error)
#
# References
#   clr       : Bertrand et al. (2024) Earth-Sci. Rev. 249, 104639, Eqs. 1-2
#   grain size: Folk & Ward (1957) J. Sediment. Petrol. 27, 3-26 (mean, phi)
#   C/N       : Meyers (1994) Chem. Geol. 114, 289-302
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
    'carbon':   ['Obrist', 'results'],
    'GUAC-29A': dict(xrf=['29A', 'RepoOG'], gs=['29A', 'FolkWard']),
    'GUAC-24A': dict(xrf=['24A', 'RepoOG'], gs=['24A', 'FolkWard']),
}
CORE_IDS = {'GUAC-29A': 'GUAC-29A-1G-1',      # core names as written in the
            'GUAC-24A': 'GUAC-24A-1G-1'}      # varve and carbon files

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
GS_SAMPLE_COL   = 'Sample'       # e.g. '15.5 - 16cm' (29A) or '13.0cm' (24A)
GS_MEAN_COL     = 'Mean_phi'     # Folk & Ward graphic mean, in phi
VARVE_DEPTH_COL = 'Depth (in core)'
# Carbon file: a value column takes its sample names from the nearest of these
# columns to its LEFT (your file has one sample column per measurement block).
CARBON_SAMPLE_COLS = ('Sample ID', 'Name', 'Sample')

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
GS_COLOR = '#4a5568'

# Carbon and nitrogen panels: (key, column in your file, panel label, colour).
# 'CN' has no column: it is calculated from TOC and N (atomic ratio).
CARBON_PANELS = [
    ('TC',    'wt %T',                     'TC\n(wt %)',             '#1a202c'),
    ('TIC',   '% Inorganic Carbon',        'TIC\n(wt %)',            '#2b6cb0'),
    ('TOC',   '%TOC by difference',        'TOC\n(wt %)',            '#276749'),
    ('CaCO3', '%CaCO3',                    'CaCO$_3$\n(wt %)',       '#3182ce'),
    ('N',     'wt %N',                     'N\n(wt %)',              '#9c4221'),
    ('d15N',  'd15N (permil, vs AIR) raw', '$\\delta^{15}$N\n(‰ AIR)', '#6b46c1'),
    ('CN',    None,                        'C/N\n(atomic)',          '#b7791f'),
]
C_ATOM, N_ATOM = 12.011, 14.007

# ---- LINES -------------------------------------------------------------------
JOIN_MODES = ('continuous', 'gaps')   # both versions are made
GAP_FACTOR = 1.6                 # 'gaps': a gap = spacing > 1.6 x usual spacing

# ---- SIZE --------------------------------------------------------------------
FONT_PT       = 24               # ALL text and scale numbers
DEPTH_SCALE   = 1.0              # cm on paper per 1 cm of core (same both cores)
PHOTO_W_CM    = 5.0              # width of the photo column (same both cores)
PANEL_W_CM    = 6.5              # width of every data panel (same both cores)
PANEL_GAP_CM  = 1.2              # gap between panels
MIN_LABEL_GAP_CM = 0.3           # least space between neighbouring panels' numbers
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
    'mathtext.default': 'regular',               # math text in the same font
    'pdf.fonttype': 42, 'ps.fonttype': 42,       # PDF text stays editable
    'axes.grid': False,
})
CM = 1 / 2.54                                    # cm -> inch
LINE_W = 2.0                                     # data line width (pt)
DOT_SIZE = 5.0                                   # sample dot diameter (pt)


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
    tab = ('.xlsx', '.xls', '.csv')
    out = [('event table', FILE_WORDS['events'], tab),
           ('varve counts', FILE_WORDS['varves'], tab),
           ('carbon results', FILE_WORDS['carbon'], tab)]
    for core in CORE_IDS:
        out.append((f'{core} XRF', FILE_WORDS[core]['xrf'], tab))
        out.append((f'{core} grain size', FILE_WORDS[core]['gs'], tab))
        out.append((f'{core} photo', PHOTOS[core]['words'], PHOTO_EXT))
    return out


def check_settings():
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
    bad = [m for m in JOIN_MODES if m not in ('continuous', 'gaps')]
    if bad:
        raise ValueError(f'JOIN_MODES may hold only "continuous" and "gaps", not {bad}')


def ensure_files():
    """Stop early, with a clear message, if a setting or a file is missing."""
    check_settings()
    if all(find_file(w, e) for _, w, e in needed()):
        return
    try:
        from google.colab import files
        print('Upload: 2 XRF (.xlsx), 2 grain size (.csv), the carbon results, '
              'the event table, the varve counts, 2 core photos')
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


def interval_depth(text):
    """Depth of a sample from its name: the centre of the sampled interval.
       '15.5 - 16cm' -> 15.75 ; '16.5-17.0cm-G2' -> 16.75 ; '13.0cm' -> 13.0
       '15.5-16' (carbon file, after the core name) -> 15.75"""
    numbers = re.findall(r'\d+(?:\.\d+)?', str(text).split('cm')[0])[:2]
    return float(np.mean([float(n) for n in numbers])) if numbers else np.nan


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
# STEP 4 - GRAIN SIZE
# =============================================================================
def read_grain_size(path):
    df = read_table(path)
    for c in (GS_SAMPLE_COL, GS_MEAN_COL):
        if c not in df.columns:
            raise KeyError(f'{os.path.basename(path)} has no column "{c}"')
    d = pd.DataFrame({'depth': df[GS_SAMPLE_COL].map(interval_depth),
                      'phi': pd.to_numeric(df[GS_MEAN_COL], errors='coerce')}).dropna()
    return d.groupby('depth', as_index=False).mean()    # repeat runs -> mean


# =============================================================================
# STEP 5 - CARBON AND NITROGEN
# =============================================================================
def read_carbon(path, core_id):
    """Your results sheet has one block per measurement, each with its own
    sample column ('Sample ID' | 'Name' | 'Sample' | 'Sample'). Each value is
    matched to the sample name in the nearest sample column to its LEFT, so the
    blocks may be in different row orders. Returns a DataFrame with depth and
    one column per CARBON_PANELS key, or None if the core has no samples."""
    raw = read_table(path, header=None)
    hdr = [str(v).strip() for v in raw.iloc[0]]
    body = raw.iloc[1:]
    sample_cols = [i for i, h in enumerate(hdr) if h in CARBON_SAMPLE_COLS]
    if not sample_cols:
        raise KeyError(f'{os.path.basename(path)}: no sample column {CARBON_SAMPLE_COLS}')
    series = {}
    for key, col, _, _ in CARBON_PANELS:
        if col is None:
            continue
        if col not in hdr:
            raise KeyError(f'{os.path.basename(path)} has no column "{col}"')
        i = hdr.index(col)
        left = [s for s in sample_cols if s < i]
        if not left:
            raise KeyError(f'"{col}" has no sample column to its left')
        names = body.iloc[:, left[-1]].astype(str).str.strip()
        vals = pd.to_numeric(body.iloc[:, i], errors='coerce')
        keep = names.str.startswith(core_id + ' ')
        s = pd.Series(vals[keep].values, index=names[keep].values, name=key)
        if s.index.duplicated().any():
            raise ValueError(f'"{col}": the same sample appears twice: '
                             f'{s.index[s.index.duplicated()].tolist()}')
        series[key] = s
    df = pd.DataFrame(series)
    if df.empty:
        return None
    df.insert(0, 'depth', [interval_depth(n[len(core_id):]) for n in df.index])
    if 'TOC' in df and 'N' in df:
        df['CN'] = (df.TOC / C_ATOM) / (df.N / N_ATOM)          # atomic C/N
    if df.depth.isna().any():
        raise ValueError(f'No depth in sample names {df.index[df.depth.isna()].tolist()}')
    return df.sort_values('depth').reset_index(drop=True)


# =============================================================================
# STEP 6 - CORE PHOTO
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
# STEP 7 - DRAWING
# =============================================================================
def in_range(depth, top, bottom):
    return (depth >= top - 1e-9) & (depth <= bottom + 1e-9)


def darker(color, f=0.72):
    r, g, b = matplotlib.colors.to_rgb(color)
    return (r * f, g * f, b * f)


def panel_list(has_carbon):
    """Every data panel, left to right: key, data set, column, label, colour,
    group, reversed axis."""
    panels = [dict(key='phi', data='gs', col='phi', label='Mean grain\nsize (φ)',
                   color=GS_COLOR, group='Grain size', reverse=True)]
    for group, els in CLR_GROUPS:
        for el in els:
            panels.append(dict(key=f'{el}_clr', data='xrf', col=f'{el}_clr',
                               label=f'{el}\n(clr)', color=ELEMENT_COLORS[el],
                               group=group))
    if has_carbon:
        for key, _, label, color in CARBON_PANELS:
            panels.append(dict(key=key, data='carbon', col=key, label=label,
                               color=color, group='Carbon and nitrogen'))
    return panels


def panel_points(df, col, top, bottom):
    """The samples of one panel inside the depth axis, without empty values,
    in depth order (each panel drops only its own empty values)."""
    d = df.loc[in_range(df.depth, top, bottom), ['depth', col]].dropna()
    return d.sort_values('depth')


def runs(depth, value, mode):
    """'continuous': one run. 'gaps': split where the spacing is more than
    GAP_FACTOR x the usual (median) spacing."""
    depth, value = np.asarray(depth, float), np.asarray(value, float)
    if mode == 'continuous' or len(depth) < 3:
        return [(depth, value)]
    step = np.median(np.diff(depth))
    cut = np.where(np.diff(depth) > GAP_FACTOR * step)[0] + 1
    return list(zip(np.split(depth, cut), np.split(value, cut)))


def draw_core(core, data, photo, photo_depths, ev, panels, xlims, depth_range,
              with_events, mode):
    top, bottom = depth_range
    n = len(panels)
    plot_h = (bottom - top) * DEPTH_SCALE
    lefts = [LEFT_CM + PHOTO_W_CM + PANEL_GAP_CM + i * (PANEL_W_CM + PANEL_GAP_CM)
             for i in range(n)]
    fig_w = lefts[-1] + PANEL_W_CM + (LABEL_W_CM if with_events else 1.0)
    fig_h = TOP_CM + plot_h + BOTTOM_CM
    fig = plt.figure(figsize=(fig_w * CM, fig_h * CM))
    frac = lambda x, y, w, h: [x / fig_w, y / fig_h, w / fig_w, h / fig_h]

    # ---- core photo: placed by its ruler depths ----------------------------
    pax = fig.add_axes(frac(LEFT_CM, BOTTOM_CM, PHOTO_W_CM, plot_h))
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

    axes = []
    for i, p in enumerate(panels):
        ax = fig.add_axes(frac(lefts[i], BOTTOM_CM, PANEL_W_CM, plot_h), sharey=pax)
        axes.append(ax)
        if with_events:
            for _, e in ev.iterrows():
                ax.axhspan(e.top, e.base, color=EVENT_COLORS[int(e.event)],
                           alpha=BAND_ALPHA, lw=0, zorder=0)
            for y in sorted(set(np.round(np.r_[ev.top, ev.base], 3))):
                ax.axhline(y, color='0.3', lw=1.2, ls=(0, (5, 3)), zorder=1)
        d = panel_points(data[p['data']], p['col'], top, bottom)
        for dd, vv in runs(d.depth, d[p['col']], mode):
            ax.plot(vv, dd, '-', color=p['color'], lw=LINE_W, marker='o',
                    ms=DOT_SIZE, mfc=p['color'], mec='white', mew=0.6,
                    zorder=3, gid='data')

        lo, hi = xlims[p['key']]                    # same range in both cores
        pad = 0.08 * (hi - lo) if hi > lo else 1.0
        lo, hi = lo - pad, hi + pad
        ax.set_xlim((hi, lo) if p.get('reverse') else (lo, hi))
        ax.xaxis.set_major_locator(MaxNLocator(3, steps=[1, 2, 5, 10]))
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position('top')
        ax.set_xlabel(p['label'], labelpad=12, fontweight='bold')
        ax.spines[['bottom', 'right', 'left']].set_visible(False)
        ax.tick_params(axis='y', which='both', left=False, labelleft=False)

    # group names with a line above the panels
    fig.canvas.draw()
    label_tops = max(a.xaxis.label.get_window_extent().y1 for a in axes + [pax])
    y_line = label_tops / fig.dpi / CM + 0.6
    groups = []
    for ax, p in zip(axes, panels):
        if groups and groups[-1][0] == p['group']:
            groups[-1][1].append(ax)
        else:
            groups.append((p['group'], [ax]))
    for name, axs in groups:
        x0 = axs[0].get_position().x0 + 0.002
        x1 = axs[-1].get_position().x1 - 0.002
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
# STEP 8 - CHECK THE DRAWN FIGURE AGAINST THE DATA (stops on any error)
# =============================================================================
def check_figure(core, fig, pax, axes, panels, data, photo_depths, ev,
                 depth_range, with_events, mode):
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
    if tuple(pax.get_ylim()) != (bottom, top):
        problems.append(f'depth axis {pax.get_ylim()} != varve range {depth_range}')
    counts = []
    for ax, p in zip(axes, panels):
        name = p['key']
        pos = ax.get_position()
        if abs(pos.width * fig_w - PANEL_W_CM) > 1e-6:
            problems.append(f'{name}: panel width {pos.width * fig_w:.3f} cm')
        if abs(pos.height * fig_h / (bottom - top) - DEPTH_SCALE) > 1e-6:
            problems.append(f'{name}: depth scale')
        want = panel_points(data[p['data']], p['col'], top, bottom)
        lines = [l for l in ax.get_lines() if l.get_gid() == 'data']
        if any(l.get_marker() != 'o' for l in lines):
            problems.append(f'{name}: samples are not drawn as dots')
        if mode == 'continuous' and len(want) and len(lines) != 1:
            problems.append(f'{name}: {len(lines)} lines, expected 1 continuous line')
        got = np.concatenate([np.c_[l.get_ydata(), l.get_xdata()] for l in lines]) \
            if lines else np.empty((0, 2))
        same_count = len(got) == len(want)
        same_values = same_count and np.allclose(got.astype(float),
                                                  want.to_numpy(float), atol=1e-12)
        if not same_values:
            problems.append(f'{name}: drawn points differ from the data')
        if mode == 'gaps' and len(want) >= 3:
            step = np.median(np.diff(want.depth))
            for l in lines:
                if len(l.get_ydata()) > 1 and np.diff(l.get_ydata()).max() > GAP_FACTOR * step + 1e-9:
                    problems.append(f'{name}: a line crosses a gap in the samples')
                    break
        counts.append(len(want))
    # tick numbers of neighbouring panels must stay MIN_LABEL_GAP_CM apart
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    def shown(ax):
        lo, hi = sorted(ax.get_xlim())
        return [t.get_window_extent(renderer) for t, v in
                zip(ax.get_xticklabels(), ax.get_xticks())
                if t.get_text() and lo - 1e-9 <= v <= hi + 1e-9]
    for (a1, p1), (a2, p2) in zip(zip(axes, panels), zip(axes[1:], panels[1:])):
        b1, b2 = shown(a1), shown(a2)
        min_gap_px = MIN_LABEL_GAP_CM * CM * fig.dpi
        if b1 and b2 and min(b.x0 for b in b2) - max(b.x1 for b in b1) < min_gap_px:
            problems.append(f'tick numbers of {p1["key"]} and {p2["key"]} are '
                            f'less than {MIN_LABEL_GAP_CM} cm apart')
    if with_events:
        labels = {t.get_gid() for t in fig.findobj(matplotlib.text.Text) if t.get_gid()}
        missing = [f'E{e}' for e in ev.event if f'E{e}' not in labels]
        if missing:
            problems.append(f'no label for {missing}')
    if problems:
        raise AssertionError(f'{core}: ' + '; '.join(problems))
    summary = ', '.join(f'{p["key"]} {c}' for p, c in zip(panels, counts)
                        if p['data'] != 'xrf')
    print(f'   CHECK PASSED ({mode}{", events" if with_events else ""}): every point '
          f'matches the data (XRF {counts[1]} per element; {summary}), all text '
          f'{FONT_PT} pt, panels {PANEL_W_CM} cm, depth scale {DEPTH_SCALE} cm/cm')


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
    carbon_path = find_file(FILE_WORDS['carbon'])
    os.makedirs(OUT_DIR, exist_ok=True)

    cores = {}
    for core in CORE_IDS:
        xp = find_file(FILE_WORDS[core]['xrf'])
        gp = find_file(FILE_WORDS[core]['gs'])
        pp = find_file(PHOTOS[core]['words'], PHOTO_EXT)
        rng = read_core_range(varve_path, CORE_IDS[core])           # STEP 2
        print(f'\n{core}\n   XRF        : {os.path.basename(xp)}'
              f'\n   grain size : {os.path.basename(gp)}')
        xrf = read_xrf(xp)                                          # STEP 3
        gs = read_grain_size(gp)                                    # STEP 4
        carbon = read_carbon(carbon_path, CORE_IDS[core])           # STEP 5
        ph = PHOTOS[core]
        photo, full, done = read_photo(pp, ph.get('rotate', 'auto'),
                                       ph.get('across', (0.0, 1.0)))  # STEP 6
        pdep = (float(ph['first_cm']), float(ph['last_cm']))
        print(f'   photo      : {os.path.basename(pp)} ({photo.shape[1]} x '
              f'{photo.shape[0]} px), placed at {pdep[0]:g}-{pdep[1]:g} cm; {done}'
              f'\n   depth axis : {rng[0]:g}-{rng[1]:g} cm (varve file)'
              f'\n   XRF        : {int(in_range(xrf.depth, *rng).sum())} samples, step '
              f'{np.median(np.diff(xrf.depth)) * 10:g} mm'
              f'\n   grain size : {int(in_range(gs.depth, *rng).sum())} samples')
        if carbon is None:
            print(f'   carbon     : no {CORE_IDS[core]} samples in '
                  f'{os.path.basename(carbon_path)} -> no carbon panels')
        else:
            n = {k: int(carbon[k].notna().sum()) for k, *_ in CARBON_PANELS}
            print(f'   carbon     : {len(carbon)} samples from '
                  f'{os.path.basename(carbon_path)}; values per panel {n}')
        if pdep[0] > rng[0] or pdep[1] < rng[1]:
            print(f'   ! the photo ({pdep[0]:g}-{pdep[1]:g} cm) does not cover the '
                  f'whole depth axis ({rng[0]:g}-{rng[1]:g} cm): that part stays blank')
        chk = save_photo_check(core, full, pdep, ph.get('across', (0.0, 1.0)), rng)
        print(f'   photo check: {os.path.basename(chk)} (red lines must sit on the '
              f'ruler\'s 10 cm marks)')
        cores[core] = dict(data=dict(xrf=xrf, gs=gs, carbon=carbon), photo=photo,
                           pdep=pdep, ev=events[core], rng=rng)

    # one x range per panel, shared by every core that has that panel
    xlims = {}
    for p in panel_list(True):
        v = pd.concat([panel_points(c['data'][p['data']], p['col'], *c['rng'])[p['col']]
                       for c in cores.values() if c['data'][p['data']] is not None])
        xlims[p['key']] = (v.min(), v.max())

    for mode in JOIN_MODES:
        for with_events, bg in ((True, 'events'), (False, 'lines')):
            for core, c in cores.items():
                panels = panel_list(c['data']['carbon'] is not None)
                fig, pax, axes = draw_core(core, c['data'], c['photo'], c['pdep'],
                                           c['ev'], panels, xlims, c['rng'],
                                           with_events, mode)               # STEP 7
                check_figure(core, fig, pax, axes, panels, c['data'], c['pdep'],
                             c['ev'], c['rng'], with_events, mode)          # STEP 8
                name = f'{core}_V7_{mode}_{bg}'
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
