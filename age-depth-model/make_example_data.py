"""
Write a SYNTHETIC varve-count workbook (example_varve_counts.xlsx) so
plot_age_depth.py can be tried without real core data.

It uses the same layout as the Lake Lachuá counting sheet: all cores side by
side on one sheet, the core ID in row 1, and Year / Depth (in core) /
Depth (adjusted) headers in row 2. An event layer is written as the same year
repeated on consecutive rows, with depth increasing from its top to its base.

A second sheet, "Events", lists the events shared by all cores in the same
layout, two rows per event (its top and base in each core).
"""

import numpy as np
import pandas as pd

rng = np.random.default_rng(42)
TOP_YEAR = 2023
CORES = {  # core ID: (oldest year, core-top offset in cm)
    "CORE-01": (1720, 10.5),
    "CORE-02": (1800, 14.8),
    "CORE-03": (1840, 9.8),
    "CORE-04": (1780, 14.1),
}
SHARED_EVENTS = {1999: 4.0, 1976: 7.0, 1954: 25.0, 1931: 2.0, 1897: 4.0}  # year: cm

blocks, event_blocks = [], []
for name, (base_year, offset) in CORES.items():
    rows, depth, shared = [], 0.0, []
    for year in range(TOP_YEAR, base_year - 1, -1):
        rows.append((year, depth))
        thick = None
        if year in SHARED_EVENTS:
            thick = SHARED_EVENTS[year] * rng.uniform(0.6, 1.4)
        elif rng.random() < 0.05:
            thick = rng.uniform(0.2, 2.5)
        if thick:
            if year in SHARED_EVENTS:
                shared += [(year, depth), (year, depth + thick)]
            depth += thick
            rows += [(year, depth), (year, depth)]  # event base, repeated as in the counts
        depth += max(0.05, round(rng.lognormal(np.log(0.12), 0.4), 2))
    year, adj = np.array(rows).T
    blocks.append(pd.DataFrame({(name, "Year"): year.astype(int),
                                (name, "Depth (in core)"): (adj + offset).round(2),
                                (name, "Depth (adjusted)"): adj.round(2),
                                (name, ""): None}))
    ey, ed = np.array(shared).T
    event_blocks.append(pd.DataFrame({(name, "Year"): ey.astype(int),
                                      (name, "Depth (in core)"): (ed + offset).round(2),
                                      (name, "Depth (adjusted)"): ed.round(2),
                                      (name, ""): None}))


def with_header(blocks):
    """Stack the core-ID row and the Year/Depth header row on top of the data."""
    table = pd.concat(blocks, axis=1)
    header = pd.DataFrame([[c[0] if c[1] == "Year" else None for c in table.columns],
                           [c[1] or None for c in table.columns]])
    return pd.concat([header, pd.DataFrame(table.to_numpy())], ignore_index=True)


with pd.ExcelWriter("example_varve_counts.xlsx") as xw:
    with_header(blocks).to_excel(xw, sheet_name="Sheet1", header=False, index=False)
    with_header(event_blocks).to_excel(xw, sheet_name="Events", header=False, index=False)
print("Wrote example_varve_counts.xlsx (synthetic data)")
