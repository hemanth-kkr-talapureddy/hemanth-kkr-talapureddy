# =============================================================================
# Lachua (GUAC-29A, GUAC-24A) - XRF clr + log-ratio event figures
#
# INPUT FILES (all five are required; nothing is typed into this code):
#   1. GUAC-29A-1G-1_RepoOG.xlsx            XRF counts, 29A
#   2. GUAC-24A-1G-1_RepoOG.xlsx            XRF counts, 24A
#   3. GUAC_29A_..._FolkWard_Stats.csv      grain-size statistics, 29A
#   4. GUAC_24A_..._FolkWard_Stats_1.csv    grain-size statistics, 24A
#   5. Lachua_Graine_Siza_Samples_N_Events_with_Ages.xlsx
#                                           event depths, ages, uncertainties
#
# OUTPUT, for each core:
#   <core>_GS_XRF_ratios : Mean grain size (phi) | Ti (clr) | Zr (clr) |
#                          ln(Mn/Fe) | ln(Ti/Ca) | ln(Zr/K)
#   <core>_XRF_ratios    : the same without the grain-size panel
#   saved as 600 dpi PNG + vector PDF, plus a QC table of every event.
#
# THE STEPS (each is one function below, in this order):
#   STEP 1  read the event file           -> read_events()
#   STEP 2  read XRF, compute clr+ratios  -> read_xrf()
#   STEP 3  read grain size               -> read_grain_size()
#   STEP 4  count samples in each event,
#           decide which events to plot   -> select_events()
#   STEP 5  draw the figure               -> draw_core()
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
DATA_DIR = '/content'            # where the 5 input files are (Colab default)
OUT_DIR  = '/content/figures'    # where the figures are written

# Each file is found by WORDS in its name (case, '-', '_', spaces ignored).
FILE_WORDS = {
    'events':    ['Events', 'Ages'],
    'GUAC-29A':  dict(xrf=['29A', 'RepoOG'], gs=['29A', 'FolkWard']),
    'GUAC-24A':  dict(xrf=['24A', 'RepoOG'], gs=['24A', 'FolkWard']),
}

# Column names in the XRF files and the event file (exactly as written there).
XRF_DEPTH_COL = 'Section Depth (cm)'
GS_SAMPLE_COL = 'Sample'         # e.g. '15.5 - 16cm' (29A) or '13.0cm' (24A)
GS_MEAN_COL   = 'Mean_phi'       # Folk & Ward graphic mean, in phi

# Elements in the clr. The clr divides each element by the geometric mean of
# THESE elements, so changing the list changes every clr value.
CLR_ELEMENTS = ['Mn', 'Al', 'Si', 'K', 'Ti', 'Fe', 'Zr', 'Ca', 'Sr', 'S']
CLR_PANELS   = ['Ti', 'Zr']                              # clr panels drawn
RATIOS       = [('Mn', 'Fe'), ('Ti', 'Ca'), ('Zr', 'K')] # ln(a/b), raw counts

# Which events are drawn (your rule): the event must have grain-size data
# AND at least this many XRF points inside its top-base interval.
MIN_XRF_POINTS = 3
MIN_GS_POINTS  = 1

# Depth axis (your rule): from the TOP of this event down to the deepest
# data point of the last event that is drawn, separately for each core.
START_EVENT = 2

# Figure size: identical for both cores, so both figures are the same length.
FIG_W_MM, FIG_H_MM = 180, 230    # 180 mm = double-column journal width
DEPTH_MAJOR_CM = 5               # labelled depth tick every 5 cm
DEPTH_MINOR_CM = 1               # small depth tick every 1 cm
SHARE_X_ACROSS_CORES = True      # same x-axis range in 24A and 29A
BREAK_AT_GAPS = True             # inside an event, do not join across missing samples
GAP_FACTOR    = 1.6              # "missing" = spacing > 1.6 x the usual spacing
DPI = 600

# Event colours: read from your event colour table (E1 ... E20).
# The table uses 8 colours that repeat: E1 = E9 = E17, E2 = E10 = E18, ...
_TABLE = ['#f6a362', '#2a9f8f', '#ebc56b', '#6e8fae',
          '#e97053', '#8bb37d', '#b5848f', '#9d8bba']
EVENT_COLORS = {e: _TABLE[(e - 1) % 8] for e in range(1, 21)}
BAND_ALPHA = 0.30                # shading behind each event

# =============================================================================
# STYLE (fonts and line widths for print)
# =============================================================================
matplotlib.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'Liberation Sans', 'DejaVu Sans'],
    'font.size': 7, 'axes.labelsize': 7.5, 'xtick.labelsize': 6.5,
    'ytick.labelsize': 7, 'axes.linewidth': 0.7,
    'xtick.major.width': 0.7, 'ytick.major.width': 0.7,
    'xtick.minor.width': 0.5, 'ytick.minor.width': 0.5,
    'xtick.major.size': 3, 'ytick.major.size': 3.5,
    'xtick.minor.size': 1.6, 'ytick.minor.size': 2,
    'xtick.direction': 'out', 'ytick.direction': 'out',
    'pdf.fonttype': 42, 'ps.fonttype': 42,      # PDF text stays editable
    'axes.grid': False,
})
MM = 1 / 25.4                                   # mm -> inch


# =============================================================================
# FINDING THE FILES
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
    for core in ('GUAC-29A', 'GUAC-24A'):
        yield FILE_WORDS[core]['xrf']
        yield FILE_WORDS[core]['gs']


def ensure_files():
    """In Colab: ask for an upload if a file is missing. Stop if still missing."""
    if all(find_file(w) for w in all_file_words()):
        return
    try:
        from google.colab import files
        print('Upload the 5 files: 2 XRF (.xlsx), 2 grain size (.csv), 1 event table (.xlsx)')
        for name, data in files.upload().items():
            with open(os.path.join(DATA_DIR, name), 'wb') as fh:
                fh.write(data)
    except ImportError:
        pass
    missing = [w for w in all_file_words() if not find_file(w)]
    if missing:
        raise FileNotFoundError(f'No file in {DATA_DIR} has these words in its name: {missing}')


def read_table(path):
    return pd.read_excel(path) if path.lower().endswith(('.xlsx', '.xls')) \
        else pd.read_csv(path, encoding='utf-8-sig')


# =============================================================================
# STEP 1 - EVENTS, AGES AND UNCERTAINTIES FROM YOUR EVENT FILE
# =============================================================================
def read_events(path):
    """Your file has one row per event and one column block per core:
         age_CE | age_unc | event | in_24A | 24A_top_cm | 24A_base_cm | 24A_n_samples
                         | event | in_29A | 29A_top_cm | 29A_base_cm | 29A_n_samples
    Returns {core: DataFrame(event, top, base, age, unc, n_file)}.
    Events marked ABSENT (or with an empty depth) are left out of that core."""
    df = read_table(path)
    need = ['age_CE', 'age_unc', 'event']
    for tag in ('24A', '29A'):
        need += [f'in_{tag}', f'{tag}_top_cm', f'{tag}_base_cm', f'{tag}_n_samples']
    missing = [c for c in need if c not in df.columns]
    if missing:
        raise KeyError(f'Event file is missing columns {missing}. Found {list(df.columns)}')

    events = {}
    for tag, core in (('24A', 'GUAC-24A'), ('29A', 'GUAC-29A')):
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
# STEP 2 - XRF: clr and log-ratios
# =============================================================================
def clr(counts):
    """clr(x) = ln(x) - mean(ln(x)) for each sample (row).
    Same as ln(x / geometric mean), but does not overflow."""
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
# STEP 3 - GRAIN SIZE
# =============================================================================
def sample_depth(label):
    """Depth of a grain-size sample from its name.
       '15.5 - 16cm'    -> 15.75  (centre of the sampled interval)
       '16.5-17.0cm-G2' -> 16.75
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
# STEP 4 - WHICH SAMPLES ARE IN WHICH EVENT, WHICH EVENTS ARE DRAWN
# =============================================================================
def in_event(depth, top, base):
    """True for samples with top <= depth <= base.
    A sample lying exactly on a boundary shared by two events is drawn once,
    at its true depth, on the dashed boundary line."""
    return (depth >= top - 1e-9) & (depth <= base + 1e-9)


def select_events(core, ev, xrf, gs):
    ev = ev.copy()
    ev['n_xrf'] = [int(in_event(xrf.depth, t, b).sum()) for t, b in zip(ev.top, ev.base)]
    ev['n_gs']  = [int(in_event(gs.depth,  t, b).sum()) for t, b in zip(ev.top, ev.base)]
    ev['drawn']  = (ev.n_xrf >= MIN_XRF_POINTS) & (ev.n_gs >= MIN_GS_POINTS)
    ev['reason'] = np.where(ev['drawn'], 'plotted',
                   np.where(ev.n_gs < MIN_GS_POINTS, 'no grain-size data',
                            f'fewer than {MIN_XRF_POINTS} XRF points'))
    print(f'\n{core}: event check (n_file = n_samples written in your event file)')
    print(ev[['event', 'top', 'base', 'age', 'unc', 'n_xrf', 'n_gs', 'n_file',
              'reason']].to_string(index=False))
    return ev


def depth_range(ev, xrf, gs):
    """Top of START_EVENT -> deepest data point of the last plotted event."""
    start = ev.loc[ev.event == START_EVENT, 'top']
    if start.empty:
        raise ValueError(f'E{START_EVENT} is not in this core')
    last = ev[ev['drawn']].iloc[-1]
    deepest = max(xrf.depth[in_event(xrf.depth, last.top, last.base)].max(),
                  gs.depth[in_event(gs.depth, last.top, last.base)].max())
    return float(start.iloc[0]), float(deepest)


# =============================================================================
# STEP 5 - DRAWING
# =============================================================================
def panel_list(with_grain_size):
    panels = []
    if with_grain_size:
        panels.append(dict(col='phi', data='gs', label='Mean grain size (φ)',
                           group='Grain size', reverse=True))   # coarse to the right
    for e in CLR_PANELS:
        panels.append(dict(col=f'{e}_clr', data='xrf', label=f'{e} (clr)',
                           group='Lithogenic (detrital)'))
    for a, b in RATIOS:
        panels.append(dict(col=f'ln_{a}_{b}', data='xrf', label=f'ln({a}/{b})',
                           group='Log-ratios'))
    return panels


def event_samples(df, ev):
    """Every sample that falls inside a plotted event."""
    keep = np.zeros(len(df), bool)
    for t, b in zip(ev.top, ev.base):
        keep |= in_event(df.depth, t, b).to_numpy()
    return df[keep]


def runs_without_gaps(depth, value):
    """Split a series where samples are missing, so no line bridges a gap."""
    depth, value = np.asarray(depth), np.asarray(value)
    if not BREAK_AT_GAPS or len(depth) < 3:
        return [(depth, value)]
    step = np.median(np.diff(depth))
    cut = np.where(np.diff(depth) > GAP_FACTOR * step)[0] + 1
    return list(zip(np.split(depth, cut), np.split(value, cut)))


def darker(color, f=0.75):
    r, g, b = matplotlib.colors.to_rgb(color)
    return (r * f, g * f, b * f)


def draw_core(core, data, ev, panels, xlims, ylims, fname):
    ev = ev[ev['drawn']]
    top, bottom = ylims
    n = len(panels)

    fig = plt.figure(figsize=(FIG_W_MM * MM, FIG_H_MM * MM))
    grid = fig.add_gridspec(1, n + 1, width_ratios=[1] * n + [0.95],
                            left=0.075, right=0.995, top=0.895, bottom=0.03,
                            wspace=0.16)
    axes = [fig.add_subplot(grid[0, i]) for i in range(n)]
    label_ax = fig.add_subplot(grid[0, n], sharey=axes[0])
    for ax in axes[1:]:
        ax.sharey(axes[0])

    # every event top and base gets a dashed line (each depth drawn once)
    boundaries = sorted(set(np.round(np.r_[ev.top, ev.base], 3)))

    for i, (ax, p) in enumerate(zip(axes, panels)):
        df = data[p['data']]
        for _, e in ev.iterrows():
            color = EVENT_COLORS[e.event]
            ax.axhspan(e.top, e.base, color=color, alpha=BAND_ALPHA, lw=0, zorder=0)
            inside = df[in_event(df.depth, e.top, e.base)]
            for dd, vv in runs_without_gaps(inside.depth, inside[p['col']]):
                ax.plot(vv, dd, '-', color=color, lw=0.9, zorder=3)
            ax.plot(inside[p['col']], inside.depth, 'o', ms=2.4, mfc=color,
                    mec='k', mew=0.25, zorder=4, clip_on=False)
        for y in boundaries:
            ax.axhline(y, color='0.3', lw=0.5, ls=(0, (4, 2.5)), zorder=2)

        # x axis: same range for both cores, 8 % margin each side
        lo, hi = xlims[p['col']]
        pad = 0.08 * (hi - lo)
        lo, hi = lo - pad, hi + pad
        ax.set_xlim((hi, lo) if p.get('reverse') else (lo, hi))
        ax.xaxis.set_major_locator(MaxNLocator(4 if p['col'] == 'phi' else 3,
                                               steps=[1, 2, 2.5, 5, 10]))
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position('top')
        ax.set_xlabel(p['label'], labelpad=4, fontweight='bold')
        ax.spines[['bottom', 'right']].set_visible(False)

        # depth axis: ticks on every panel, numbers on the first one only
        ax.yaxis.set_major_locator(MultipleLocator(DEPTH_MAJOR_CM))
        ax.yaxis.set_minor_locator(MultipleLocator(DEPTH_MINOR_CM))
        if i == 0:
            ax.set_ylabel('Depth (cm)')
        else:
            ax.tick_params(axis='y', labelleft=False)
    axes[0].set_ylim(bottom, top)                         # depth increases downwards

    # group names with a line above the panels
    fig.canvas.draw()
    groups = []
    for ax, p in zip(axes, panels):
        if groups and groups[-1][0] == p['group']:
            groups[-1][1].append(ax)
        else:
            groups.append((p['group'], [ax]))
    for name, axs in groups:
        x0 = axs[0].get_position().x0 + 0.004
        x1 = axs[-1].get_position().x1 - 0.004
        fig.add_artist(plt.Line2D([x0, x1], [0.945, 0.945], color='k', lw=0.8,
                                  transform=fig.transFigure))
        fig.text((x0 + x1) / 2, 0.951, name, ha='center', va='bottom',
                 fontsize=8, fontweight='bold')

    # event name + age +/- uncertainty, at the middle of each event
    label_ax.axis('off')
    min_gap = 0.022 * (bottom - top)                      # keep labels apart
    y_prev = -np.inf
    for _, e in ev.iterrows():
        y = max((e.top + e.base) / 2, y_prev + min_gap)
        y_prev = y
        label_ax.text(0.04, y, f'E{e.event:02d}  {e.age:.1f} ± {e.unc:.1f}',
                      color=darker(EVENT_COLORS[e.event]), fontsize=6.5,
                      fontweight='bold', va='center', ha='left',
                      transform=label_ax.get_yaxis_transform(), clip_on=False)
    label_ax.text(0.04, top, 'Event  Age (CE)', fontsize=6.5, fontweight='bold',
                  va='bottom', ha='left', transform=label_ax.get_yaxis_transform())

    fig.text(0.01, 0.995, core.replace('-', ' '), fontsize=11, fontweight='bold',
             ha='left', va='top')
    os.makedirs(OUT_DIR, exist_ok=True)
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(OUT_DIR, f'{fname}.{ext}'), dpi=DPI)
    if _in_notebook():
        plt.show()
    plt.close(fig)
    print(f'   saved {fname}.png / .pdf  (depth {top:g}-{bottom:g} cm)')


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
    events = read_events(find_file(FILE_WORDS['events']))      # STEP 1

    cores, qc = {}, []
    for core in ('GUAC-29A', 'GUAC-24A'):
        xp, gp = find_file(FILE_WORDS[core]['xrf']), find_file(FILE_WORDS[core]['gs'])
        print(f'\n{core}\n   XRF        : {os.path.basename(xp)}'
              f'\n   grain size : {os.path.basename(gp)}')
        xrf = read_xrf(xp)                                     # STEP 2
        gs = read_grain_size(gp)                               # STEP 3
        ev = select_events(core, events[core], xrf, gs)        # STEP 4
        cores[core] = dict(xrf=xrf, gs=gs, ev=ev, ylims=depth_range(ev, xrf, gs))
        qc.append(ev.assign(core=core))

    for with_gs, tag in ((True, 'GS_XRF_ratios'), (False, 'XRF_ratios')):   # STEP 5
        panels = panel_list(with_gs)

        def limits(core_names):
            return {p['col']: (min(event_samples(cores[c][p['data']], cores[c]['ev'][cores[c]['ev']['drawn']])[p['col']].min() for c in core_names),
                               max(event_samples(cores[c][p['data']], cores[c]['ev'][cores[c]['ev']['drawn']])[p['col']].max() for c in core_names))
                    for p in panels}

        for core, c in cores.items():
            xlims = limits(list(cores) if SHARE_X_ACROSS_CORES else [core])
            draw_core(core, c, c['ev'], panels, xlims, c['ylims'], f'{core}_{tag}')

    pd.concat(qc)[['core', 'event', 'top', 'base', 'age', 'unc', 'n_xrf', 'n_gs',
                   'n_file', 'reason']].to_csv(os.path.join(OUT_DIR, 'event_check.csv'),
                                               index=False)
    print(f'\nAll figures and event_check.csv are in {OUT_DIR}')


if __name__ == '__main__':
    main()
