# =============================================================================
# Lachua (GUAC-29A, GUAC-24A) - magnetic susceptibility (MSCL), poster style
#
# INPUT FILES (nothing is typed into this code):
#   1. GUAC_MSCL-XYZ_rerun.csv              MSCL data of every GUAC core
#   2. Lachua_Graine_Siza_Samples_N_Events_with_Ages.xlsx
#                                           event depths, ages, uncertainties
#   3. Lachua_Varve_counts.xlsx             top and bottom depth of each core
#
# OUTPUT, for each core:
#   <core>_MS_events : magnetic susceptibility, with event bands + ages
#   <core>_MS_lines  : magnetic susceptibility, line only
#   saved as 300 dpi PNG + vector PDF
#
# THE RULES
#   * only the rows of GUAC-24A-1G and GUAC-29A-1G are used
#   * every MS measurement inside the core's depth range is drawn as ONE
#     continuous line (no break at background, varves or events)
#   * depth axis of each core: first -> last 'Depth (in core)' in the varve file
#   * ONE common depth scale, same panel width, same x range, 24 pt text
#
# THE STEPS (each is one function below, in this order):
#   STEP 1  event file                  -> read_events()
#   STEP 2  varve file: depth range     -> read_core_range()
#   STEP 3  MSCL file: MS of each core  -> read_ms()
#   STEP 4  draw the figure             -> draw_core()
#   STEP 5  check the drawn figure      -> check_figure()   (stops on any error)
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
DATA_DIR = '/content'            # where the input files are (Colab default)
OUT_DIR  = '/content/figures'    # where the figures are written

# Each file is found by WORDS in its name (case, '-', '_', spaces ignored).
FILE_WORDS = {
    'ms':     ['MSCL'],
    'events': ['Events', 'Ages'],
    'varves': ['Varve', 'counts'],
}
CORE_IDS = {'GUAC-29A': 'GUAC-29A-1G-1',      # names in the varve file
            'GUAC-24A': 'GUAC-24A-1G-1'}
MS_CORE_PREFIX = {'GUAC-29A': 'GUAC-29A-1G',  # start of 'Section' in the
                  'GUAC-24A': 'GUAC-24A-1G'}  # MSCL file (e.g. GUAC-29A-1G-1-A)

# Column names, exactly as written in your files.
MS_SECTION_COL  = 'Section'
MS_DEPTH_COL    = 'Section Depth'
MS_COL          = 'Magnetic Susceptibility'
VARVE_DEPTH_COL = 'Depth (in core)'

MS_COLOR = '#2d3748'             # colour of the MS line

# ---- SIZE --------------------------------------------------------------------
FONT_PT       = 24               # ALL text and scale numbers
DEPTH_SCALE   = 1.0              # cm on paper per 1 cm of core (same both cores)
PANEL_W_CM    = 10.0             # width of the MS panel (same both cores)
LABEL_W_CM    = 11.0             # column for "E01  2010.5 ± 1.7" (events version)
LEFT_CM       = 4.5              # room for the depth numbers + "Depth (cm)"
TOP_CM        = 6.5              # room for title and panel name
BOTTOM_CM     = 1.5
DEPTH_MAJOR_CM = 5               # labelled depth tick every 5 cm
DEPTH_MINOR_CM = 1               # small depth tick every 1 cm
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
LINE_W = 2.6                                    # data line width (pt)


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


def ensure_files():
    """In Colab: ask for an upload if a file is missing. Stop if still missing."""
    if all(find_file(w) for w in FILE_WORDS.values()):
        return
    try:
        from google.colab import files
        print('Upload 3 files: the MSCL .csv, the event table (.xlsx), '
              'the varve counts (.xlsx)')
        for name, data in files.upload().items():
            with open(os.path.join(DATA_DIR, name), 'wb') as fh:
                fh.write(data)
    except ImportError:
        pass
    missing = [k for k, w in FILE_WORDS.items() if not find_file(w)]
    if missing:
        raise FileNotFoundError(f'Missing in {DATA_DIR}: {missing} file '
                                f'(names must contain {[FILE_WORDS[k] for k in missing]})')


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
# STEP 3 - MAGNETIC SUSCEPTIBILITY OF ONE CORE FROM THE MSCL FILE
# =============================================================================
def read_ms(path, prefix):
    """The MSCL csv has a header row, then a UNITS row (cm, mm, SI x 10^-5 ...),
    then one row per measurement for every GUAC core. Keeps the rows whose
    Section starts with `prefix`. Returns (DataFrame depth/ms, unit text)."""
    units = pd.read_csv(path, header=None, nrows=2, encoding='utf-8-sig')
    header = [str(v).strip() for v in units.iloc[0]]
    for c in (MS_SECTION_COL, MS_DEPTH_COL, MS_COL):
        if c not in header:
            raise KeyError(f'{os.path.basename(path)} has no column "{c}"')
    unit = str(units.iloc[1, header.index(MS_COL)]).strip()

    df = pd.read_csv(path, skiprows=[1], low_memory=False, encoding='utf-8-sig')
    df.columns = [str(c).strip() for c in df.columns]
    s = df[df[MS_SECTION_COL].astype(str).str.startswith(prefix + '-')]
    if s.empty:
        raise ValueError(f'No rows for {prefix} in {os.path.basename(path)}')
    sections = s[MS_SECTION_COL].unique().tolist()
    if len(sections) != 1:
        # several sections would need their depths joined end to end; stop
        # rather than guess how they join
        raise ValueError(f'{prefix} has {len(sections)} sections {sections}; '
                         f'this script expects one section per core')
    d = pd.DataFrame({'depth': pd.to_numeric(s[MS_DEPTH_COL], errors='coerce'),
                      'ms': pd.to_numeric(s[MS_COL], errors='coerce')})
    empty = d[d.ms.isna()]
    if len(empty):
        print(f'   {prefix}: {len(empty)} depths without an MS value '
              f'({empty.depth.min():g}-{empty.depth.max():g} cm) are left out')
    d = d.dropna().sort_values('depth').reset_index(drop=True)
    return d, unit, sections[0]


# =============================================================================
# STEP 4 - DRAWING
# =============================================================================
def in_range(depth, top, bottom):
    return (depth >= top - 1e-9) & (depth <= bottom + 1e-9)


def darker(color, f=0.72):
    r, g, b = matplotlib.colors.to_rgb(color)
    return (r * f, g * f, b * f)


def draw_core(core, ms, unit, ev, xlim, depth_range, with_events):
    top, bottom = depth_range
    plot_h = (bottom - top) * DEPTH_SCALE
    fig_w = LEFT_CM + PANEL_W_CM + (LABEL_W_CM if with_events else 1.5)
    fig_h = TOP_CM + plot_h + BOTTOM_CM
    fig = plt.figure(figsize=(fig_w * CM, fig_h * CM))
    frac = lambda x, y, w, h: [x / fig_w, y / fig_h, w / fig_w, h / fig_h]
    ax = fig.add_axes(frac(LEFT_CM, BOTTOM_CM, PANEL_W_CM, plot_h))

    if with_events:
        for _, e in ev.iterrows():
            ax.axhspan(e.top, e.base, color=EVENT_COLORS[int(e.event)],
                       alpha=BAND_ALPHA, lw=0, zorder=0)
        for y in sorted(set(np.round(np.r_[ev.top, ev.base], 3))):
            ax.axhline(y, color='0.3', lw=1.2, ls=(0, (5, 3)), zorder=1)

    # every MS value inside the depth range, as ONE continuous line
    d = ms[in_range(ms.depth, top, bottom)]
    ax.plot(d.ms, d.depth, '-', color=MS_COLOR, lw=LINE_W, zorder=3, gid='data')

    lo, hi = xlim                                   # same range in both cores
    pad = 0.08 * (hi - lo)
    ax.set_xlim(lo - pad, hi + pad)
    ax.xaxis.set_major_locator(MaxNLocator(4, steps=[1, 2, 2.5, 5, 10]))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')
    # 'SI x 10^-5' -> SI × 10⁻⁵, drawn with math text (every font has it)
    unit_txt = re.sub(r'x\s*10\^(-?\d+)', r'$\\times\\,10^{\1}$', unit)
    ax.set_xlabel(f'Magnetic susceptibility\n({unit_txt})', labelpad=12,
                  fontweight='bold')
    ax.spines[['bottom', 'right']].set_visible(False)
    ax.set_ylim(bottom, top)                        # depth increases downwards
    ax.yaxis.set_major_locator(MultipleLocator(DEPTH_MAJOR_CM))
    ax.yaxis.set_minor_locator(MultipleLocator(DEPTH_MINOR_CM))
    ax.set_ylabel('Depth (cm)', labelpad=10)

    if with_events:                                 # event name + age ± unc
        lab = fig.add_axes(frac(LEFT_CM + PANEL_W_CM + 0.6, BOTTOM_CM,
                                LABEL_W_CM - 0.6, plot_h), sharey=ax)
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
    return fig, ax


# =============================================================================
# STEP 5 - CHECK THE DRAWN FIGURE AGAINST THE DATA (stops on any error)
# =============================================================================
def check_figure(core, fig, ax, ms, ev, depth_range, with_events):
    problems = []
    top, bottom = depth_range
    for t in fig.findobj(matplotlib.text.Text):
        if t.get_text().strip() and t.get_visible() and abs(t.get_fontsize() - FONT_PT) > 1e-6:
            problems.append(f'text "{t.get_text()[:20]}" is {t.get_fontsize()} pt')
    fig_w, fig_h = fig.get_size_inches() / CM
    pos = ax.get_position()
    if abs(pos.width * fig_w - PANEL_W_CM) > 1e-6:
        problems.append(f'panel width {pos.width * fig_w:.3f} cm')
    if abs(pos.height * fig_h / (bottom - top) - DEPTH_SCALE) > 1e-6:
        problems.append('depth scale')
    if tuple(ax.get_ylim()) != (bottom, top):
        problems.append(f'depth axis {ax.get_ylim()} != varve range {depth_range}')
    lines = [l for l in ax.get_lines() if l.get_gid() == 'data']
    want = ms[in_range(ms.depth, top, bottom)]
    if len(lines) != 1:
        problems.append(f'{len(lines)} MS lines, expected 1 continuous line')
    else:
        got = np.c_[lines[0].get_ydata(), lines[0].get_xdata()].astype(float)
        if np.isnan(got).any():
            problems.append('the MS line contains a gap (NaN)')
        elif not np.allclose(got, want[['depth', 'ms']].to_numpy(float), atol=1e-12):
            problems.append('drawn MS line differs from the data')
    if with_events:
        labels = {t.get_gid() for t in fig.findobj(matplotlib.text.Text) if t.get_gid()}
        missing = [f'E{e}' for e in ev.event if f'E{e}' not in labels]
        if missing:
            problems.append(f'no label for {missing}')
    if problems:
        raise AssertionError(f'{core}: ' + '; '.join(problems))
    print(f'   CHECK PASSED: {len(want)} MS values in one continuous line, '
          f'all text {FONT_PT} pt, panel {PANEL_W_CM} cm, depth scale {DEPTH_SCALE} cm/cm'
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
    varve_path, ms_path = find_file(FILE_WORDS['varves']), find_file(FILE_WORDS['ms'])

    cores = {}
    for core in CORE_IDS:
        rng = read_core_range(varve_path, CORE_IDS[core])           # STEP 2
        ms, unit, section = read_ms(ms_path, MS_CORE_PREFIX[core])  # STEP 3
        d = ms[in_range(ms.depth, *rng)]
        print(f'\n{core}\n   MSCL section : {section} ({os.path.basename(ms_path)})'
              f'\n   depth axis   : {rng[0]:g}-{rng[1]:g} cm (varve file)'
              f'\n   MS values    : {len(d)} inside the depth axis, '
              f'{d.ms.min():.3f} to {d.ms.max():.3f} {unit}')
        cores[core] = dict(ms=ms, unit=unit, ev=events[core], rng=rng)

    # one x range, shared by both cores
    v = pd.concat([c['ms'].loc[in_range(c['ms'].depth, *c['rng']), 'ms']
                   for c in cores.values()])
    xlim = (v.min(), v.max())

    os.makedirs(OUT_DIR, exist_ok=True)
    for with_events, tag in ((True, 'MS_events'), (False, 'MS_lines')):
        for core, c in cores.items():
            fig, ax = draw_core(core, c['ms'], c['unit'], c['ev'], xlim, c['rng'],
                                with_events)                                # STEP 4
            check_figure(core, fig, ax, c['ms'], c['ev'], c['rng'], with_events)  # STEP 5
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
