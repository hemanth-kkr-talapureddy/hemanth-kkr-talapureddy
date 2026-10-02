"""
Write a SYNTHETIC varve-count workbook (example_varve_counts.xlsx) in the
one-sheet-per-core layout, so plot_age_depth.py can be tried without real
core data. Shape loosely mimics eight cores with a slow lower section, a
thick event deposit in the 1860s, and fast recent accumulation.
"""

import numpy as np
import pandas as pd

rng = np.random.default_rng(42)
TOP_YEAR = 2023
cores = {  # core: (base year, depth of event top, event thickness)
    "GUAC-22A-1G-1": (1840, 20, 33),
    "GUAC-23A-1G-1": (1785, 17, 31),
    "GUAC-24A-1G-1": (1720, 13, 25),
    "GUAC-25A-1G-1": (1835, 19, 28),
    "GUAC-26A-1G-1": (1810, 25, 23),
    "GUAC-27A-1G-1": (1780, 18, 27),
    "GUAC-28A-1G-1": (1815, 13, 24),
    "GUAC-29A-1G-1": (1800, 15, 23),
}

with pd.ExcelWriter("example_varve_counts.xlsx") as xw:
    for name, (base, ev_top, ev_thick) in cores.items():
        years = np.arange(TOP_YEAR, base - 1, -1)
        thick = np.where(years > 1864, ev_top / (TOP_YEAR - 1864), 0.22)
        thick = thick * rng.lognormal(0, 0.35, len(years))
        # occasional thin event layers
        thick += np.where(rng.random(len(years)) < 0.04, rng.uniform(0.5, 2.5, len(years)), 0)
        thick[years == 1864] += ev_thick  # thick event deposit
        depth = np.concatenate([[0], np.cumsum(thick)[:-1]])
        pd.DataFrame({"Depth (cm)": depth.round(2), "Year (CE)": years}).to_excel(
            xw, sheet_name=name, index=False)

print("Wrote example_varve_counts.xlsx (synthetic data)")
