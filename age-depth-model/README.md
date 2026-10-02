# Varve age-depth model plot

Plots varve-counted age-depth models for several sediment cores as one
three-panel figure with linked age, depth and thickness scales:

- **a** age-depth model: age (yr CE) vs depth (cm), one line per core
- **b** thickness (mm, log) vs depth, on the same depth axis as **a**: background
  varves as points, event layers as bars spanning their depth interval
- **c** thickness (mm, log) through time, one strip per core, on the same age
  axis as **a**, with event layers as markers

Event layers are read from the counts: when the same year is listed on
consecutive rows with increasing depth, that depth interval is an
instantaneous deposit (e.g. a turbidite) within that year.

## Run in Google Colab

Open `age_depth_colab.ipynb` in Colab (**File → Upload notebook**), then
**Runtime → Run all**. It asks you to upload your `.xlsx` file and downloads the
figure as PNG and PDF.

## Run locally

```bash
pip install pandas numpy matplotlib openpyxl
python plot_age_depth.py my_varve_counts.xlsx -o age_depth_model.pdf --xlim 1700 2050 --ylim 90 0
```

Options:

| Option | Meaning |
|---|---|
| `--depth-col "in core"` | Use a different depth column (default: the one whose header contains "adjusted") |
| `--event-mm 10` | Dashed thickness line on panels b and c |
| `--layout simple` | Age-depth curves only, like the original Excel chart |
| `--sheet Sheet1` | Read only one sheet |

## Workbook layouts

Detected from the header cells (case-insensitive):

| Layout | Columns |
|---|---|
| All cores on one sheet | Row 1: core ID above each block. Row 2: `Year`, `Depth (in core)`, `Depth (adjusted)` |
| One sheet per core (sheet name = core ID) | `Year`, `Depth (cm)` |
| Varve thickness only | `Thickness (mm)`, one row per varve from the top. Pass `--top-year 2023` |

If the same core appears on several sheets (e.g. an extra "Events" sheet), only
the first one is used. Depths are in cm.

## Example

`make_example_data.py` writes a **synthetic** workbook in the all-cores-on-one-sheet
layout so you can try the script without real core data. `example_age_depth.png`
is the result.

![example](example_age_depth.png)
