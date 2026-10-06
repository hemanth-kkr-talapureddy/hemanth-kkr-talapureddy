# =============================================================================
# Lachua (GUAC-29A, GUAC-24A) - XRF clr + log-ratio event figures, poster style
#
# INPUT FILES (all six are required; nothing is typed into this code):
#   1. GUAC-29A-1G-1_RepoOG.xlsx            XRF counts, 29A
#   2. GUAC-24A-1G-1_RepoOG.xlsx            XRF counts, 24A
#   3. GUAC_29A_..._FolkWard_Stats.csv      grain-size statistics, 29A
#   4. GUAC_24A_..._FolkWard_Stats_1.csv    grain-size statistics, 24A
#   5. Lachua_Graine_Siza_Samples_N_Events_with_Ages.xlsx
#                                           event depths, ages, uncertainties
#   6. Lachua_Varve_counts.xlsx             top and bottom depth of each core
#
# OUTPUT, for each core:
#   <core>_GS_XRF_ratios : Mean grain size (phi) | Ti (clr) | Zr (clr) |
#                          ln(Mn/Fe) | ln(Ti/Ca) | ln(Zr/K)
#   <core>_XRF_ratios    : the same without the grain-size panel
#   saved as 300 dpi PNG + vector PDF, plus event_check.csv
#
# THE RULES (from your instructions)
#   * depth axis of each core: first -> last 'Depth (in core)' in the varve file
#   * EVERY event present in the core (in_24A / in_29A = Yes) is drawn:
#     shaded band, dashed top and base, label "E01  2010.5 ± 1.7"
#   * XRF points are drawn only in events with >= 3 XRF points
#     (thin turbidites with 1-2 points get the band and label only)
#   * grain-size points are drawn wherever grain-size data exist in an event
#   * ONE common depth scale: 1 cm of core = DEPTH_SCALE cm on paper in both
#     cores, so the longer core gives the taller figure
#   * same panel width, same x-axis ranges and 24 pt text in both cores
#
# THE STEPS (each is one function below, in this order):
#   STEP 1  event file                  -> read_events()
#   STEP 2  varve file: depth range     -> read_core_range()
#   STEP 3  XRF: clr + log-ratios       -> read_xrf()
#   STEP 4  grain size                  -> read_grain_size()
#   STEP 5  count samples per event     -> count_events()
#   STEP 6  draw the figure             -> draw_core()
#   STEP 7  check the drawn figure      -> check_figure()   (stops on any error)
#
# References
#   clr       : Bertrand et al. (2024) Earth-Sci. Rev. 249, 104639, Eqs. 1-2
#   grain size: Folk & Ward (1957) J. Sediment. Petrol. 27, 3-26 (mean, phi)
# =============================================================================
import os, re, glob
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, MultipleLocator, AutoMinorLocator

# =============================================================================
# SETTINGS - the only lines you should need to change
# =============================================================================
DATA_DIR = '/content'            # where the 6 input files are (Colab default)
OUT_DIR  = '/content/figures'    # where the figures are written

# Each file is found by WORDS in its name (case, '-', '_', spaces ignored).
FILE_WORDS = {
    'events':   ['Events', 'Ages'],
    'varves':   ['Varve', 'counts'],
    'GUAC-29A': dict(xrf=['29A', 'RepoOG'], gs=['29A', 'FolkWard']),
    'GUAC-24A': dict(xrf=['24A', 'RepoOG'], gs=['24A', 'FolkWard']),
}
CORE_IDS = {'GUAC-29A': 'GUAC-29A-1G-1',      # core names as written in the
            'GUAC-24A': 'GUAC-24A-1G-1'}      # varve file and XRF files

# Column names, exactly as written in your files.
XRF_DEPTH_COL   = 'Section Depth (cm)'
GS_SAMPLE_COL   = 'Sample'       # e.g. '15.5 - 16cm' (29A) or '13.0cm' (24A)
GS_MEAN_COL     = 'Mean_phi'     # Folk & Ward graphic mean, in phi
VARVE_DEPTH_COL = 'Depth (in core)'

# Elements in the clr (the geometric mean is taken over THESE elements).
CLR_ELEMENTS = ['Mn', 'Al', 'Si', 'K', 'Ti', 'Fe', 'Zr', 'Ca', 'Sr', 'S']
CLR_PANELS   = ['Ti', 'Zr']                              # clr panels drawn
RATIOS       = [('Mn', 'Fe'), ('Ti', 'Ca'), ('Zr', 'K')] # ln(a/b), raw counts

MIN_XRF_POINTS = 3               # XRF is drawn in an event only with >= this

# ---- SIZE: change these to set the poster size ------------------------------
FONT_PT       = 24               # ALL text and scale numbers
DEPTH_SCALE   = 1.0              # cm on paper per 1 cm of core (same both cores)
PANEL_W_CM    = 10.0             # width of every data panel (same both cores)
PANEL_GAP_CM  = 1.6              # gap between panels
LABEL_W_CM    = 11.0             # column for "E01  2010.5 ± 1.7"
LEFT_CM       = 4.5              # room for the depth numbers + "Depth (cm)"
TOP_CM        = 8.5              # room for title, group names, panel names
BOTTOM_CM     = 1.5
DEPTH_MAJOR_CM = 5               # labelled depth tick every 5 cm
DEPTH_MINOR_CM = 1               # small depth tick every 1 cm
DPI = 300

BREAK_AT_GAPS = True             # inside an event, do not join across missing samples
GAP_FACTOR    = 1.6              # "missing" = spacing > 1.6 x the usual spacing

# Event colours from your event colour table (8 colours that repeat:
# E1 = E9 = E17, E2 = E10 = E18, ...).
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
LINE_W, MARKER = 2.6, 9                         # data line width, dot size (pt)


# =============================================================================
# FINDING AND READING FILES
# =============================================================================
def _norm(s):
    """'GUAC-29A_GSD' -> 'guac29agsd' so names can be compared loosely."""
    return re.sub(r'[^a-z0-9]', '', str(s).lower())


def find_file(words):
    """The newest file in DATA_DIR whose name contains every word."""
    hits = [f for f in glob.glob(os.path.join(DATA_DIR, '*'))
            if f.lower().endswith(('.csv', '.xlsx', '.xls'))
            and all(_norm(w) in _norm(os.path.basename(f)) for w in words)]
    return max(hits, key=os.path.getmtime) if hits else None


def all_file_words():
    yield FILE_WORDS['events']
    yield FILE_WORDS['varves']
    for core in CORE_IDS:
        yield FILE_WORDS[core]['xrf']
        yield FILE_WORDS[core]['gs']


def ensure_files():
    """In Colab: ask for an upload if a file is missing. Stop if still missing."""
    if all(find_file(w) for w in all_file_words()):
        return
    try:
        from google.colab import files
        print('Upload the 6 files: 2 XRF (.xlsx), 2 grain size (.csv), '
              'the event table (.xlsx) and the varve counts (.xlsx)')
        for name, data in files.upload().items():
            with open(os.path.join(DATA_DIR, name), 'wb') as fh:
                fh.write(data)
    except ImportError:
        pass
    missing = [w for w in all_file_words() if not find_file(w)]
    if missing:
        raise FileNotFoundError(f'No file in {DATA_DIR} has these words in its name: {missing}')


def read_table(path, **kw):
    return pd.read_excel(path, **kw) if path.lower().endswith(('.xlsx', '.xls')) \
        else pd.read_csv(path, encoding='utf-8-sig', **kw)


# =============================================================================
# STEP 1 - EVENTS, AGES AND UNCERTAINTIES FROM YOUR EVENT FILE
# =============================================================================
def read_events(path):
    """One row per event, one column block per core:
         age_CE | age_unc | event | in_24A | 24A_top_cm | 24A_base_cm | 24A_n_samples
                         | event | in_29A | 29A_top_cm | 29A_base_cm | 29A_n_samples
    Events marked ABSENT (or with an empty depth) are not in that core."""
    df = read_table(path)
    need = ['age_CE', 'age_unc', 'event']
    for tag in ('24A', '29A'):
        need += [f'in_{tag}', f'{tag}_top_cm', f'{tag}_base_cm', f'{tag}_n_samples']
    missing = [c for c in need if c not in df.columns]
    if missing:
        raise KeyError(f'Event file is missing columns {missing}. Found {list(df.columns)}')

    events = {}
    for core in CORE_IDS:
        tag = core.split('-')[1]                        # 'GUAC-24A' -> '24A'
        rows = []
        for _, r in df.iterrows():
            present = str(r[f'in_{tag}']).strip().lower() == 'yes'
            top, base = r[f'{tag}_top_cm'], r[f'{tag}_base_cm']
            if not present or pd.isna(top) or pd.isna(base):
                continue
            rows.append(dict(event=int(re.search(r'\d+', str(r['event'])).group()),
                             top=float(top), base=float(base),
                             age=float(r['age_CE']), unc=float(r['age_unc']),
                             n_file=r[f'{tag}_n_samples']))
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
# STEP 3 - XRF: clr and log-ratios
# =============================================================================
def clr(counts):
    """clr(x) = ln(x) - mean(ln(x)) for each sample (row)."""
    L = np.log(counts)
    return L - L.mean(axis=1, keepdims=True)


def read_xrf(path):
    df = read_table(path)
    df.columns = [str(c).strip() for c in df.columns]
    elements = sorted(set(CLR_ELEMENTS) | {e for pair in RATIOS for e in pair})
    missing = [c for c in [XRF_DEPTH_COL] + elements if c not in df.columns]
    if missing:
        raise KeyError(f'{os.path.basename(path)} is missing columns {missing}')

    d = df[[XRF_DEPTH_COL] + elements].apply(pd.to_numeric, errors='coerce')
    d = d.rename(columns={XRF_DEPTH_COL: 'depth'}).dropna()
    d = d.groupby('depth', as_index=False).mean()       # repeat scans -> mean

    # ln(0) does not exist: a zero count is replaced by half the smallest
    # non-zero count of that element (Bertrand et al. 2024). It is printed.
    for e in elements:
        zero = d[e] <= 0
        if zero.any():
            half = d.loc[~zero, e].min() / 2
            print(f'   {e}: {zero.sum()} zero count(s) replaced by {half:g}')
            d.loc[zero, e] = half

    C = clr(d[CLR_ELEMENTS].to_numpy(float))
    assert np.abs(C.sum(axis=1)).max() < 1e-9          # clr rows must sum to 0
    for i, e in enumerate(CLR_ELEMENTS):
        d[f'{e}_clr'] = C[:, i]
    for a, b in RATIOS:
        d[f'ln_{a}_{b}'] = np.log(d[a] / d[b])
    return d.sort_values('depth').reset_index(drop=True)


# =============================================================================
# STEP 4 - GRAIN SIZE
# =============================================================================
def sample_depth(label):
    """Depth of a grain-size sample from its name.
       '15.5 - 16cm'    -> 15.75  (centre of the sampled interval)
       '13.0cm'         -> 13.0   (24A samples are named by one depth)"""
    numbers = re.findall(r'\d+(?:\.\d+)?', str(label).split('cm')[0])[:2]
    return float(np.mean([float(n) for n in numbers])) if numbers else np.nan


def read_grain_size(path):
    df = read_table(path)
    for c in (GS_SAMPLE_COL, GS_MEAN_COL):
        if c not in df.columns:
            raise KeyError(f'{os.path.basename(path)} has no column "{c}"')
    d = pd.DataFrame({'depth': df[GS_SAMPLE_COL].map(sample_depth),
                      'phi': pd.to_numeric(df[GS_MEAN_COL], errors='coerce')}).dropna()
    return d.groupby('depth', as_index=False).mean()    # repeat runs -> mean


# =============================================================================
# STEP 5 - SAMPLES IN EACH EVENT
# =============================================================================
def in_event(depth, top, base):
    """True for samples with top <= depth <= base (boundaries included)."""
    return (depth >= top - 1e-9) & (depth <= base + 1e-9)


def count_events(core, ev, xrf, gs):
    ev = ev.copy()
    ev['n_xrf'] = [int(in_event(xrf.depth, t, b).sum()) for t, b in zip(ev.top, ev.base)]
    ev['n_gs']  = [int(in_event(gs.depth,  t, b).sum()) for t, b in zip(ev.top, ev.base)]
    ev['xrf_drawn'] = ev.n_xrf >= MIN_XRF_POINTS
    ev['gs_drawn']  = ev.n_gs > 0
    print(f'\n{core}: every event is drawn (band + label); data drawn as below')
    print(ev[['event', 'top', 'base', 'age', 'unc', 'n_xrf', 'xrf_drawn',
              'n_gs', 'gs_drawn']].to_string(index=False))
    return ev


def drawn_samples(df, ev, flag):
    """Every sample inside an event whose data are drawn (flag column)."""
    keep = np.zeros(len(df), bool)
    for t, b in zip(ev.top[ev[flag]], ev.base[ev[flag]]):
        keep |= in_event(df.depth, t, b).to_numpy()
    return df[keep]


# =============================================================================
# STEP 6 - DRAWING
# =============================================================================
def panel_list(with_grain_size):
    panels = []
    if with_grain_size:
        panels.append(dict(col='phi', data='gs', flag='gs_drawn', reverse=True,
                           label='Mean grain\nsize (φ)', group='Grain size'))
    for e in CLR_PANELS:
        panels.append(dict(col=f'{e}_clr', data='xrf', flag='xrf_drawn',
                           label=f'{e} (clr)', group='Lithogenic (detrital)'))
    for a, b in RATIOS:
        panels.append(dict(col=f'ln_{a}_{b}', data='xrf', flag='xrf_drawn',
                           label=f'ln({a}/{b})', group='Log-ratios'))
    return panels


def runs_without_gaps(depth, value):
    """Split a series where samples are missing, so no line bridges a gap."""
    depth, value = np.asarray(depth), np.asarray(value)
    if not BREAK_AT_GAPS or len(depth) < 3:
        return [(depth, value)]
    step = np.median(np.diff(depth))
    cut = np.where(np.diff(depth) > GAP_FACTOR * step)[0] + 1
    return list(zip(np.split(depth, cut), np.split(value, cut)))


def darker(color, f=0.72):
    r, g, b = matplotlib.colors.to_rgb(color)
    return (r * f, g * f, b * f)


def layout(n_panels, depth_top, depth_bottom):
    """Figure size and panel positions in cm. Panel width and the depth scale
    are fixed numbers, so they are identical in every figure."""
    plot_h = (depth_bottom - depth_top) * DEPTH_SCALE
    fig_w = LEFT_CM + n_panels * PANEL_W_CM + (n_panels - 1) * PANEL_GAP_CM + LABEL_W_CM
    fig_h = TOP_CM + plot_h + BOTTOM_CM
    lefts = [LEFT_CM + i * (PANEL_W_CM + PANEL_GAP_CM) for i in range(n_panels)]
    return fig_w, fig_h, plot_h, lefts


def draw_core(core, data, ev, panels, xlims, depth_range):
    top, bottom = depth_range
    n = len(panels)
    fig_w, fig_h, plot_h, lefts = layout(n, top, bottom)
    fig = plt.figure(figsize=(fig_w * CM, fig_h * CM))
    frac = lambda x, y, w, h: [x / fig_w, y / fig_h, w / fig_w, h / fig_h]

    axes = []
    for i in range(n):
        ax = fig.add_axes(frac(lefts[i], BOTTOM_CM, PANEL_W_CM, plot_h),
                          sharey=axes[0] if axes else None)
        axes.append(ax)
    label_x0 = lefts[-1] + PANEL_W_CM + 0.6
    label_ax = fig.add_axes(frac(label_x0, BOTTOM_CM, LABEL_W_CM - 0.6, plot_h),
                            sharey=axes[0])

    boundaries = sorted(set(np.round(np.r_[ev.top, ev.base], 3)))

    for i, (ax, p) in enumerate(zip(axes, panels)):
        df = data[p['data']]
        for _, e in ev.iterrows():
            color = EVENT_COLORS[e.event]
            ax.axhspan(e.top, e.base, color=color, alpha=BAND_ALPHA, lw=0, zorder=0)
            if not e[p['flag']]:
                continue                                  # band only, no data
            inside = df[in_event(df.depth, e.top, e.base)]
            for dd, vv in runs_without_gaps(inside.depth, inside[p['col']]):
                ax.plot(vv, dd, '-', color=color, lw=LINE_W, zorder=3)
            ax.plot(inside[p['col']], inside.depth, 'o', ms=MARKER, mfc=color,
                    mec='k', mew=0.9, zorder=4, clip_on=False, gid='data')
        for y in boundaries:
            ax.axhline(y, color='0.3', lw=1.2, ls=(0, (5, 3)), zorder=2)

        # x axis: same range in both cores, 8 % margin each side
        lo, hi = xlims[p['col']]
        pad = 0.08 * (hi - lo)
        lo, hi = lo - pad, hi + pad
        ax.set_xlim((hi, lo) if p.get('reverse') else (lo, hi))
        ax.xaxis.set_major_locator(MaxNLocator(3, steps=[1, 2, 2.5, 5, 10]))
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position('top')
        ax.set_xlabel(p['label'], labelpad=12, fontweight='bold')
        ax.spines[['bottom', 'right']].set_visible(False)

        # ONE depth scale: only the first panel shows the depth axis
        ax.yaxis.set_major_locator(MultipleLocator(DEPTH_MAJOR_CM))
        ax.yaxis.set_minor_locator(MultipleLocator(DEPTH_MINOR_CM))
        if i == 0:
            ax.set_ylabel('Depth (cm)', labelpad=10)
        else:
            ax.spines['left'].set_visible(False)
            ax.tick_params(axis='y', which='both', left=False, labelleft=False)
    axes[0].set_ylim(bottom, top)                         # depth increases downwards

    # group names with a line above the panels
    fig.canvas.draw()
    label_tops = max(ax.xaxis.label.get_window_extent().y1 for ax in axes)
    y_line = label_tops / fig.dpi / CM + 0.6             # cm, above the panel names
    groups = []
    for ax, p in zip(axes, panels):
        if groups and groups[-1][0] == p['group']:
            groups[-1][1].append(ax)
        else:
            groups.append((p['group'], [ax]))
    for name, axs in groups:
        x0 = axs[0].get_position().x0 + 0.004
        x1 = axs[-1].get_position().x1 - 0.004
        fig.add_artist(plt.Line2D([x0, x1], [y_line / fig_h] * 2, color='k', lw=2,
                                  transform=fig.transFigure))
        fig.text((x0 + x1) / 2, (y_line + 0.3) / fig_h, name, ha='center',
                 va='bottom', fontweight='bold')

    # event labels: name + age ± uncertainty, pushed apart when events are
    # thin and close; a thin line links a moved label to its event
    label_ax.axis('off')
    min_gap = 1.15 * FONT_PT / 72 * 2.54 / DEPTH_SCALE   # one text line, in core cm
    y_prev = -np.inf
    for _, e in ev.iterrows():
        mid = (e.top + e.base) / 2
        y = max(mid, y_prev + min_gap)
        y_prev = y
        tr = label_ax.get_yaxis_transform()
        if abs(y - mid) > 0.05:
            label_ax.plot([0.0, 0.05], [mid, y], color='0.4', lw=1.2,
                          transform=tr, clip_on=False)
        label_ax.text(0.07, y, f'E{e.event:02d}  {e.age:.1f} ± {e.unc:.1f}',
                      color=darker(EVENT_COLORS[e.event]), fontweight='bold',
                      va='center', ha='left', transform=tr, clip_on=False,
                      gid=f'E{e.event}')
    label_ax.text(0.07, top, 'Event  Age (CE)', fontweight='bold', va='bottom',
                  ha='left', transform=label_ax.get_yaxis_transform())

    fig.text(0.5 / fig_w, 1 - 0.4 / fig_h, core.replace('-', ' '),
             fontweight='bold', ha='left', va='top')
    return fig, axes


# =============================================================================
# STEP 7 - CHECK THE DRAWN FIGURE AGAINST THE DATA (stops on any error)
# =============================================================================
def check_figure(core, fig, axes, panels, data, ev, depth_range):
    problems = []
    # every text is FONT_PT
    for t in fig.findobj(matplotlib.text.Text):
        if t.get_text().strip() and t.get_visible() and abs(t.get_fontsize() - FONT_PT) > 1e-6:
            problems.append(f'text "{t.get_text()[:20]}" is {t.get_fontsize()} pt')
    fig_w, fig_h = fig.get_size_inches() / CM
    for ax, p in zip(axes, panels):
        pos = ax.get_position()
        # same panel width and same depth scale in every figure
        if abs(pos.width * fig_w - PANEL_W_CM) > 1e-6:
            problems.append(f'{p["col"]}: panel width {pos.width * fig_w:.3f} cm')
        scale = pos.height * fig_h / (depth_range[1] - depth_range[0])
        if abs(scale - DEPTH_SCALE) > 1e-6:
            problems.append(f'{p["col"]}: depth scale {scale:.4f}')
        # the points drawn are exactly the data inside drawn events
        want = drawn_samples(data[p['data']], ev, p['flag'])
        got = [l for l in ax.get_lines() if l.get_gid() == 'data']
        gx = np.concatenate([l.get_xdata() for l in got]) if got else np.array([])
        gy = np.concatenate([l.get_ydata() for l in got]) if got else np.array([])
        # compare (depth, value) pairs; a sample on a boundary shared by two
        # events is drawn in both, so compare the unique pairs
        a = np.unique(np.round(np.c_[gy, gx], 9), axis=0)
        b = np.unique(np.round(want[['depth', p['col']]].to_numpy(float), 9), axis=0)
        if not np.array_equal(a, b):
            problems.append(f'{p["col"]}: drawn points differ from the data')
    if tuple(axes[0].get_ylim()) != (depth_range[1], depth_range[0]):
        problems.append(f'depth axis {axes[0].get_ylim()} != varve range {depth_range}')
    labels = {t.get_gid() for t in fig.findobj(matplotlib.text.Text) if t.get_gid()}
    missing = [f'E{e}' for e in ev.event if f'E{e}' not in labels]
    if missing:
        problems.append(f'no label for {missing}')
    if problems:
        raise AssertionError(f'{core}: ' + '; '.join(problems))
    print(f'   CHECK PASSED: {len(ev)} events labelled, all text {FONT_PT} pt, '
          f'panels {PANEL_W_CM} cm wide, depth scale {DEPTH_SCALE} cm/cm, '
          f'every drawn point matches the data')


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

    cores, qc = {}, []
    for core in CORE_IDS:
        xp, gp = find_file(FILE_WORDS[core]['xrf']), find_file(FILE_WORDS[core]['gs'])
        rng = read_core_range(varve_path, CORE_IDS[core])           # STEP 2
        print(f'\n{core}\n   XRF        : {os.path.basename(xp)}'
              f'\n   grain size : {os.path.basename(gp)}'
              f'\n   depth axis : {rng[0]:g} - {rng[1]:g} cm (varve file)')
        xrf, gs = read_xrf(xp), read_grain_size(gp)                 # STEPS 3-4
        ev_out = events[core][(events[core].top < rng[0]) |
                                (events[core].base > rng[1])]
        if len(ev_out):
            print(f'   ! events outside the varve depth range: '
                  f'{["E%d" % e for e in ev_out.event]}')
        ev = count_events(core, events[core], xrf, gs)              # STEP 5
        cores[core] = dict(xrf=xrf, gs=gs, ev=ev, rng=rng)
        qc.append(ev.assign(core=core))

    os.makedirs(OUT_DIR, exist_ok=True)
    for with_gs, tag in ((True, 'GS_XRF_ratios'), (False, 'XRF_ratios')):
        panels = panel_list(with_gs)
        # one x range per panel, shared by both cores (drawn data only)
        xlims = {}
        for p in panels:
            v = pd.concat([drawn_samples(c[p['data']], c['ev'], p['flag'])[p['col']]
                           for c in cores.values()])
            xlims[p['col']] = (v.min(), v.max())
        for core, c in cores.items():
            fig, axes = draw_core(core, c, c['ev'], panels, xlims, c['rng'])  # STEP 6
            check_figure(core, fig, axes, panels, c, c['ev'], c['rng'])       # STEP 7
            name = f'{core}_{tag}'
            for ext in ('png', 'pdf'):
                fig.savefig(os.path.join(OUT_DIR, f'{name}.{ext}'), dpi=DPI)
            w, h = fig.get_size_inches() / CM
            print(f'   saved {name}.png / .pdf  ({w:.1f} x {h:.1f} cm)')
            if _in_notebook():
                plt.show()
            plt.close(fig)

    pd.concat(qc)[['core', 'event', 'top', 'base', 'age', 'unc', 'n_xrf', 'xrf_drawn',
                   'n_gs', 'gs_drawn', 'n_file']].to_csv(
        os.path.join(OUT_DIR, 'event_check.csv'), index=False)
    print(f'\nAll figures and event_check.csv are in {OUT_DIR}')


if __name__ == '__main__':
    main()
