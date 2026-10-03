# =============================================================================
# Lachua (GUAC-29A, GUAC-24A) - XRF clr + log-ratio event figures
#
# Makes, for EACH core, two figures in the layout of the GUAC 29A reference:
#   A) Mean grain size (phi) | Ti (clr) | Zr (clr) | ln(Mn/Fe) | ln(Ti/Ca) | ln(Zr/K)
#   B) the same without grain size (XRF elements + ratios only)
#
# * one colour per event, taken from the event colour table (E1 ... E20)
# * both cores share the same figure size, depth axis and x-axis limits,
#   so the two figures are the same length and can be compared directly
# * saved as 600 dpi PNG + vector PDF (editable text)
#
# Works in Google Colab (asks for the files) or locally (reads DATA_DIR).
#
# clr       : Bertrand et al. (2024) Earth-Sci. Rev. 249, 104639, Eqs. 1-2
# grain size: Folk & Ward (1957) graphic mean, in phi
# =============================================================================
import os, re, glob
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, MultipleLocator, AutoMinorLocator

# =============================================================================
# SETTINGS
# =============================================================================
DATA_DIR = '/content'          # Colab default; change for a local run
OUT_DIR  = '/content/figures'

# File names are matched by WORDS (case, '-', '_' and spaces ignored).
FILES = {
    'GUAC-29A': dict(xrf=['29A', 'RepoOG'], gs=['29A', 'FolkWard']),
    'GUAC-24A': dict(xrf=['24A', 'RepoOG'], gs=['24A', 'FolkWard']),
}

# Elements used for the clr (the denominator = geometric mean of THESE).
CLR_ELEMENTS = ['Mn', 'Al', 'Si', 'K', 'Ti', 'Fe', 'Zr', 'Ca', 'Sr', 'S']
CLR_PANELS   = ['Ti', 'Zr']                      # clr panels that are drawn
RATIOS       = [('Mn', 'Fe'), ('Ti', 'Ca'), ('Zr', 'K')]   # ln(a/b), raw counts
DROP_INVALID = True            # drop XRF rows with validity = 0

# Optional event table (csv/xlsx) with columns:
#   core | event | top_cm | base_cm | age_CE | age_err
# None = use the event depths and ages written below.
EVENTS_FILE = None

MIN_SAMPLES = 4                # an event is drawn if XRF OR grain size has >= this
DEPTH_LIM   = (10.0, 90.0)     # SAME depth axis for both cores (top, bottom cm)
FIG_W_MM    = 180              # double-column journal width
FIG_H_MM    = 230              # same height for both cores
SHADE_EVENTS = True            # pastel band behind each event (colour table)
SHOW_BACKGROUND = False        # True = also draw non-event samples in light grey
BREAK_AT_GAPS = True           # line stops where samples are missing
GAP_FACTOR  = 1.6              # gap = spacing > 1.6 x normal sample spacing
SHARE_X_ACROSS_CORES = True    # same x-axis limits in 24A and 29A
GS_CLIP_PCT = None            # grain-size axis spans the finest -> this percentile;
                               # coarser outliers sit at the edge as a triangle
                               # with their real value written next to them
                               # (None = no clipping, axis fits every point; try 99)
DPI         = 600

# ---- event colours: the event colour table (8 colours, repeating) -----------
# Hex = the solid colour; the table shows it at ~45 % opacity (BAND_ALPHA).
_BASE = ['#f6a362', '#2a9f8f', '#ebc56b', '#6e8fae',
         '#e97053', '#8bb37d', '#b5848f', '#9d8bba']
EVENT_COLORS = {e: _BASE[(e - 1) % 8] for e in range(1, 21)}
BAND_ALPHA = 0.30              # band opacity (table = 0.45; lower keeps lines clear)

# ---- event ages (CE, +/- 1 sigma) and depths (cm) ---------------------------
AGES = {1:(2010.5,1.7), 2:(2000.1,1.4), 3:(1977.9,2.5), 4:(1966.1,3.4),
        5:(1962.0,3.7), 6:(1960.6,3.7), 7:(1960.6,3.7), 8:(1931.5,2.3),
        9:(1921.0,3.9),10:(1902.1,3.7),11:(1895.6,3.9),12:(1893.8,4.3),
       13:(1886.7,4.8),14:(1881.3,5.3),15:(1861.9,6.2),16:(1859.9,6.5),
       17:(1845.3,6.6),18:(1836.7,3.8),19:(1820.6,4.4),20:(1795.5,3.5)}
EVENTS = {
 'GUAC-29A': {1:(13.40,14.10), 2:(15.10,19.10), 3:(22.60,30.15), 4:(31.30,33.20),
              5:(33.55,36.80), 6:(36.80,41.00), 7:(41.00,49.10), 8:(53.25,55.80),
              9:(56.80,58.10),10:(60.20,60.60),11:(61.20,64.10),12:(64.65,66.60),
             13:(67.15,67.50),14:(67.90,68.30),15:(70.80,71.15),16:(71.25,71.50),
             17:(73.15,75.40),18:(76.65,80.00),19:(87.75,87.75)},
 'GUAC-24A': {1:(11.75,11.90), 2:(12.90,17.40), 3:(19.35,26.20), 4:(27.50,27.85),
              5:(28.35,30.80), 6:(30.90,33.00), 7:(33.00,48.00), 8:(51.70,53.00),
             10:(55.70,55.85),11:(56.50,62.00),14:(63.50,63.75),15:(66.25,66.50),
             16:(66.70,67.10),17:(68.80,69.75),18:(70.15,70.55),19:(71.30,72.70),
             20:(75.20,79.80)},
}
CONTIG_TOL = 0.15              # cm: events this close share a dashed boundary

# =============================================================================
# STYLE
# =============================================================================
matplotlib.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'Liberation Sans', 'DejaVu Sans'],
    'font.size': 7, 'axes.labelsize': 7.5, 'axes.titlesize': 8,
    'xtick.labelsize': 6.5, 'ytick.labelsize': 7,
    'axes.linewidth': 0.7, 'xtick.major.width': 0.7, 'ytick.major.width': 0.7,
    'xtick.minor.width': 0.5, 'ytick.minor.width': 0.5,
    'xtick.major.size': 3, 'ytick.major.size': 3,
    'xtick.minor.size': 1.6, 'ytick.minor.size': 1.6,
    'xtick.direction': 'out', 'ytick.direction': 'out',
    'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none',
    'axes.grid': False, 'savefig.dpi': DPI,
})
MM = 1 / 25.4


# =============================================================================
# FILE FINDING / READING
# =============================================================================
def _norm(s):
    return re.sub(r'[^a-z0-9]', '', str(s).lower())


def find_file(words):
    hits = [f for f in glob.glob(os.path.join(DATA_DIR, '*'))
            if f.lower().endswith(('.csv', '.xlsx', '.xls'))
            and all(_norm(w) in _norm(os.path.basename(f)) for w in words)]
    return max(hits, key=os.path.getmtime) if hits else None


def ensure_files():
    """In Colab, ask for uploads if any file is missing."""
    missing = [w for c in FILES.values() for w in c.values() if not find_file(w)]
    if not missing:
        return
    try:
        from google.colab import files
        os.makedirs(DATA_DIR, exist_ok=True)
        print('Upload the 2 XRF (RepoOG .xlsx) and 2 grain-size (FolkWard .csv) files'
              + (' and the event table' if EVENTS_FILE else ''))
        up = files.upload()
        for name, data in up.items():
            with open(os.path.join(DATA_DIR, name), 'wb') as fh:
                fh.write(data)
    except ImportError:
        pass
    still = [w for c in FILES.values() for w in c.values() if not find_file(w)]
    if still:
        raise FileNotFoundError(f'No file in {DATA_DIR} matches {still}')


def read_table(path):
    if path.lower().endswith(('.xlsx', '.xls')):
        return pd.read_excel(path)
    return pd.read_csv(path, encoding='utf-8-sig')


def clr(X):
    """clr(x) = ln(x) - mean(ln(x)), row-wise (Bertrand et al. 2024, Eq. 2)."""
    L = np.log(X)
    return L - L.mean(axis=1, keepdims=True)


def load_xrf(path):
    df = read_table(path)
    df.columns = [str(c).strip() for c in df.columns]
    dcol = next(c for c in df.columns if 'depth' in c.lower())
    if DROP_INVALID and 'validity' in df.columns:
        df = df[pd.to_numeric(df['validity'], errors='coerce') == 1]
    need = sorted(set(CLR_ELEMENTS) | {e for r in RATIOS for e in r})
    missing = [e for e in need if e not in df.columns]
    assert not missing, f'{os.path.basename(path)} has no column(s) {missing}'
    d = df[[dcol] + need].apply(pd.to_numeric, errors='coerce').dropna()
    d = d.rename(columns={dcol: 'depth'}).groupby('depth', as_index=False).mean()
    for e in need:                      # zeros -> half the minimum non-zero
        nz = d.loc[d[e] > 0, e]
        if (d[e] <= 0).any():
            print(f'   {e}: {(d[e] <= 0).sum()} zero counts -> {nz.min() / 2:g}')
            d.loc[d[e] <= 0, e] = nz.min() / 2
    C = clr(d[CLR_ELEMENTS].to_numpy(float))
    assert np.abs(C.sum(axis=1)).max() < 1e-9
    for i, e in enumerate(CLR_ELEMENTS):
        d[f'{e}_clr'] = C[:, i]
    for a, b in RATIOS:
        d[f'ln_{a}_{b}'] = np.log(d[a] / d[b])
    return d.sort_values('depth').reset_index(drop=True)


def _sample_depth(label):
    """'15.5 - 16cm' -> 15.75 ; '16.5-17.0cm-G2' -> 16.75 ; '13.0cm' -> 13.0"""
    nums = re.findall(r'\d+(?:\.\d+)?', str(label).split('cm')[0])
    if not nums:
        return np.nan
    v = [float(n) for n in nums[:2]]
    return float(np.mean(v))


def load_gs(path):
    df = read_table(path)
    scol = df.columns[0]
    pcol = next(c for c in df.columns if _norm(c) in ('meanphi', 'mzphi', 'meanfw'))
    d = pd.DataFrame({'depth': df[scol].map(_sample_depth),
                      'phi': pd.to_numeric(df[pcol], errors='coerce')}).dropna()
    return d.groupby('depth', as_index=False).mean()   # repeats averaged


def load_event_table(path):
    df = read_table(path)
    cols = {_norm(c): c for c in df.columns}
    pick = lambda *ks: next((cols[k] for k in ks if k in cols), None)
    cc, ec = pick('core', 'corename'), pick('event', 'eventid', 'id')
    tc, bc = pick('topcm', 'top', 'depthtop'), pick('basecm', 'base', 'bottomcm', 'bottom')
    ac, uc = pick('agece', 'age'), pick('ageerr', 'err', 'ageunc', 'uncertainty')
    for _, r in df.iterrows():
        core = 'GUAC-24A' if '24a' in _norm(r[cc]) else 'GUAC-29A'
        n = int(re.search(r'\d+', str(r[ec])).group())
        EVENTS.setdefault(core, {})[n] = (float(r[tc]), float(r[bc]))
        if ac and pd.notna(r[ac]):
            AGES[n] = (float(r[ac]), float(r[uc]) if uc and pd.notna(r[uc]) else 0.0)
    print(f'Events read from {os.path.basename(path)}')


# =============================================================================
# EVENTS
# =============================================================================
def in_event(depth, top, base):
    return (depth >= top - 1e-6) & (depth <= base + 1e-6)


def kept_events(core, xrf, gs):
    """Events with >= MIN_SAMPLES in XRF or in grain size."""
    out = {}
    for n, (t, b) in sorted(EVENTS[core].items()):
        nx = int(in_event(xrf.depth, t, b).sum())
        ng = int(in_event(gs.depth, t, b).sum())
        if max(nx, ng) >= MIN_SAMPLES:
            out[n] = (t, b)
        print(f'   E{n:02d} {t:6.2f}-{b:6.2f} cm  XRF n={nx:2d}  GS n={ng:2d}  '
              f'{"plotted" if n in out else "skipped (< MIN_SAMPLES)"}')
    return out


def segments(depth, value):
    """Split into runs without gaps so missing samples are not bridged."""
    depth, value = np.asarray(depth), np.asarray(value)
    if not BREAK_AT_GAPS or len(depth) < 3:
        return [(depth, value)]
    step = np.median(np.diff(depth))
    cut = np.where(np.diff(depth) > GAP_FACTOR * step)[0] + 1
    return list(zip(np.split(depth, cut), np.split(value, cut)))


# =============================================================================
# PLOTTING
# =============================================================================
def panel_specs(with_gs):
    specs = []
    if with_gs:
        specs.append(dict(key='phi', src='gs', label='Mean grain size (φ)',
                          group='Grain size', invert=True))
    for e in CLR_PANELS:
        specs.append(dict(key=f'{e}_clr', src='xrf', label=f'{e} (clr)',
                          group='Lithogenic (detrital)'))
    for a, b in RATIOS:
        specs.append(dict(key=f'ln_{a}_{b}', src='xrf', label=f'ln({a}/{b})',
                          group='Log-ratios'))
    return specs


def event_values(df, key, events):
    """All samples that fall in a plotted event."""
    m = np.zeros(len(df), bool)
    for t, b in events.values():
        m |= in_event(df.depth, t, b).to_numpy()
    return df.loc[m, key]


def _pad(lo, hi, f=0.08):
    r = hi - lo if hi > lo else 1.0
    return lo - f * r, hi + f * r


def plot_core(core, xrf, gs, events, specs, xlims, fname):
    src = {'xrf': xrf, 'gs': gs}
    n = len(specs)
    label_w = 0.95                     # width (in panel units) of the event-label column
    fig = plt.figure(figsize=(FIG_W_MM * MM, FIG_H_MM * MM))
    gsp = fig.add_gridspec(1, n + 1, width_ratios=[1] * n + [label_w],
                           left=0.075, right=0.995, top=0.895, bottom=0.03,
                           wspace=0.16)
    axes = [fig.add_subplot(gsp[0, i]) for i in range(n)]
    lab_ax = fig.add_subplot(gsp[0, n], sharey=axes[0])
    for ax in axes[1:]:
        ax.sharey(axes[0])

    ev_sorted = sorted(events.items(), key=lambda kv: kv[1][0])
    # dashed boundaries where two plotted events touch
    shared_bounds = [b1 for (_, (t1, b1)), (_, (t2, b2))
                     in zip(ev_sorted, ev_sorted[1:]) if abs(t2 - b1) <= CONTIG_TOL]

    for i, (ax, s) in enumerate(zip(axes, specs)):
        df = src[s['src']]
        if SHADE_EVENTS:
            for e, (t, b) in events.items():
                ax.axhspan(t, b, color=EVENT_COLORS[e], alpha=BAND_ALPHA,
                           lw=0, zorder=0)
        if SHOW_BACKGROUND:
            for dd, vv in segments(df.depth, df[s['key']]):
                ax.plot(vv, dd, color='0.75', lw=0.5, zorder=1)
        # grain size: coarse outliers (small phi) are drawn at the axis edge
        clip = None
        if s['key'] == 'phi' and GS_CLIP_PCT:
            clip = np.nanpercentile(xlims['phi'], 100 - GS_CLIP_PCT)
        for e, (t, b) in events.items():
            sub = df[in_event(df.depth, t, b)].copy()
            if clip is not None:
                far = sub[s['key']] < clip
                for _, r in sub[far].iterrows():
                    ax.plot(clip, r.depth, '>', ms=3.6, mfc=EVENT_COLORS[e], mec='k',
                            mew=0.3, zorder=5, clip_on=False)
                    ax.annotate(f'{r[s["key"]]:.1f}', (clip, r.depth),
                                xytext=(-4, 0), textcoords='offset points',
                                ha='right', va='center', fontsize=5.5, color='0.2')
                sub.loc[far, s['key']] = clip
                sub_pts = sub[~far]
            else:
                sub_pts = sub
            for dd, vv in segments(sub.depth, sub[s['key']]):
                ax.plot(vv, dd, '-', color=EVENT_COLORS[e], lw=0.9, zorder=3,
                        solid_joinstyle='round')
            ax.plot(sub_pts[s['key']], sub_pts.depth, 'o', ms=2.4, mfc=EVENT_COLORS[e],
                    mec='k', mew=0.25, zorder=4, clip_on=False)
        for yb in shared_bounds:
            ax.axhline(yb, color='0.3', lw=0.6, ls=(0, (4, 2.5)), zorder=2)

        lo, hi = _pad(np.nanmin(xlims[s['key']]), np.nanmax(xlims[s['key']]))
        ax.set_xlim((hi, lo) if s.get('invert') else (lo, hi))
        ax.xaxis.set_major_locator(MaxNLocator(3 if s['key'] != 'phi' else 4,
                                               steps=[1, 2, 2.5, 5, 10]))
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position('top')
        ax.set_xlabel(s['label'], labelpad=4, fontweight='bold')
        ax.spines['bottom'].set_visible(False)
        ax.spines['right'].set_visible(False)
        if i == 0:
            ax.set_ylabel('Depth (cm)')
            ax.yaxis.set_major_locator(MultipleLocator(10))
            ax.yaxis.set_minor_locator(MultipleLocator(2))
        else:
            ax.tick_params(axis='y', labelleft=False)
            ax.spines['left'].set_visible(True)
    axes[0].set_ylim(DEPTH_LIM[1], DEPTH_LIM[0])

    # ---- group brackets above the panels --------------------------------
    fig.canvas.draw()
    groups = []
    for ax, s in zip(axes, specs):
        if groups and groups[-1][0] == s['group']:
            groups[-1][1].append(ax)
        else:
            groups.append((s['group'], [ax]))
    y_br = 0.945
    for name, axs in groups:
        x0 = axs[0].get_position().x0 + 0.004
        x1 = axs[-1].get_position().x1 - 0.004
        fig.add_artist(plt.Line2D([x0, x1], [y_br, y_br], color='k', lw=0.8,
                                  transform=fig.transFigure))
        fig.text((x0 + x1) / 2, y_br + 0.006, name, ha='center', va='bottom',
                 fontsize=8, fontweight='bold')

    # ---- event labels: "E02  2000.1 ± 1.4" in the event colour ------------
    lab_ax.axis('off')
    mids = [((t + b) / 2, e) for e, (t, b) in ev_sorted]
    min_gap = 0.022 * (DEPTH_LIM[1] - DEPTH_LIM[0])
    ys = []
    for y, _ in mids:                        # push labels apart if they collide
        ys.append(max(y, ys[-1] + min_gap) if ys else y)
    for (y0, e), y in zip(mids, ys):
        age, err = AGES.get(e, (np.nan, np.nan))
        txt = f'E{e:02d}' + (f'  {age:.1f} ± {err:.1f}' if np.isfinite(age) else '')
        lab_ax.text(0.04, y, txt, color=_darker(EVENT_COLORS[e]), fontsize=6.5,
                    fontweight='bold', va='center', ha='left',
                    transform=lab_ax.get_yaxis_transform())
    lab_ax.text(0.04, DEPTH_LIM[0], 'Event  Age (CE)', fontsize=6.5,
                fontweight='bold', va='bottom', ha='left',
                transform=lab_ax.get_yaxis_transform())

    fig.text(0.01, 0.995, core.replace('-', ' '), fontsize=11, fontweight='bold',
             ha='left', va='top')
    os.makedirs(OUT_DIR, exist_ok=True)
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(OUT_DIR, f'{fname}.{ext}'), dpi=DPI)
    plt.show() if _in_notebook() else None
    plt.close(fig)
    print(f'   saved {fname}.png / .pdf')


def _darker(hex_color, f=0.78):
    """Slightly darker shade for text, so light colours stay legible."""
    r, g, b = matplotlib.colors.to_rgb(hex_color)
    return (r * f, g * f, b * f)


def _in_notebook():
    try:
        from IPython import get_ipython
        return get_ipython() is not None
    except ImportError:
        return False


# =============================================================================
# RUN
# =============================================================================
def main():
    ensure_files()
    if EVENTS_FILE:
        load_event_table(os.path.join(DATA_DIR, EVENTS_FILE)
                         if not os.path.isabs(EVENTS_FILE) else EVENTS_FILE)
    data, events = {}, {}
    for core, f in FILES.items():
        xp, gp = find_file(f['xrf']), find_file(f['gs'])
        print(f'\n{core}\n   XRF : {os.path.basename(xp)}\n   GS  : {os.path.basename(gp)}')
        xrf, gs = load_xrf(xp), load_gs(gp)
        print(f'   XRF {len(xrf)} rows {xrf.depth.min()}-{xrf.depth.max()} cm | '
              f'GS {len(gs)} rows {gs.depth.min()}-{gs.depth.max()} cm')
        data[core] = {'xrf': xrf, 'gs': gs}
        events[core] = kept_events(core, xrf, gs)

    for with_gs, tag in ((True, 'GS_XRF_ratios'), (False, 'XRF_ratios')):
        specs = panel_specs(with_gs)
        xl = {}
        for s in specs:
            vals = [data[c][s['src']][s['key']] if SHOW_BACKGROUND
                    else event_values(data[c][s['src']], s['key'], events[c])
                    for c in data]
            xl[s['key']] = pd.concat(vals).to_numpy()
        if with_gs and GS_CLIP_PCT:          # phi axis: finest .. clip percentile
            v = xl['phi']
            xl['phi'] = v[v >= np.nanpercentile(v, 100 - GS_CLIP_PCT)]
        for core in data:
            if not SHARE_X_ACROSS_CORES:
                xl = {s['key']: (data[core][s['src']][s['key']] if SHOW_BACKGROUND
                                 else event_values(data[core][s['src']], s['key'],
                                                   events[core])).to_numpy()
                      for s in specs}
                if with_gs and GS_CLIP_PCT:
                    v = xl['phi']
                    xl['phi'] = v[v >= np.nanpercentile(v, 100 - GS_CLIP_PCT)]
            plot_core(core, data[core]['xrf'], data[core]['gs'], events[core],
                      specs, xl, f'{core}_{tag}')

    # Event-sample table, for checking what was plotted
    rows = []
    for core in data:
        for e, (t, b) in events[core].items():
            sub = data[core]['xrf'][in_event(data[core]['xrf'].depth, t, b)]
            for _, r in sub.iterrows():
                rows.append(dict(core=core, event=f'E{e:02d}', depth_cm=r.depth,
                                 **{f'{x}_clr': r[f'{x}_clr'] for x in CLR_PANELS},
                                 **{f'ln_{a}_{b}': r[f'ln_{a}_{b}'] for a, b in RATIOS}))
    pd.DataFrame(rows).to_csv(os.path.join(OUT_DIR, 'event_samples_xrf.csv'),
                              index=False)
    print(f'\nAll figures in {OUT_DIR}')


if __name__ == '__main__':
    main()
