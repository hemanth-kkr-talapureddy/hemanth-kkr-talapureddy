# Varve age-depth model plot

Plots varve-counted age-depth models for several sediment cores on one figure,
in the same style as the Excel chart: age (yr CE) on a top axis, depth (cm)
increasing downward, and one stepped line per core.

```bash
pip install pandas numpy matplotlib openpyxl
python plot_age_depth.py my_varve_counts.xlsx -o age_depth_model.png --xlim 1700 2050 --ylim 90 0
```

Save as `.pdf` or `.svg` instead of `.png` for a vector figure. Use `--style line`
to connect points directly instead of drawing steps.

## Workbook layouts

The script detects the layout from the column headers (case-insensitive):

| Layout | Columns |
|---|---|
| One sheet per core (sheet name = core ID) | `Depth (cm)`, `Year (CE)` |
| All cores on one sheet | `GUAC-22A-1G-1 Depth`, `GUAC-22A-1G-1 Year`, `GUAC-23A-1G-1 Depth`, ... |
| Varve thickness only | `Thickness (mm)`, one row per varve from the top. Pass `--top-year 2023` |

## Example

`make_example_data.py` writes a **synthetic** workbook so you can try the script
without real core data. `example_age_depth.png` is the result.

![example](example_age_depth.png)
