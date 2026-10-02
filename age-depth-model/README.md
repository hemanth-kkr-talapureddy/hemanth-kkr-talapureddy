# Varve age-depth model plot

Plots varve-counted age-depth models for several sediment cores as one figure
with linked age and depth scales:

- **a** age-depth model: age (yr CE) vs depth (cm), one line per core. Every
  event layer is drawn as a thick bar on its core's curve. Each event correlated
  between cores (from an **Events** sheet) is shaded as a **coloured column
  hanging from the age axis**: as wide as the event's age range across the
  cores (each layer spanning its varve year, ±0.5 yr), marked by a coloured bar
  on the age scale, with its lower edge running through each core's event base.
  So the column matches both the age and the depth scale. Columns are labelled E1, E2, …, and a
  table beside the plot gives each event's age range, depth range and number of
  cores (event names coloured as in the plot).
- **b** thickness (mm, log) through time, one strip per core, on the same age
  axis as **a**, with event layers as markers and each event's age range shaded
  in the same colour as in **a**

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
| `--event-mm 10` | Dashed thickness line on panel b |
| `--events-sheet Events` | Sheet with correlated events (default: a sheet whose name contains "event") |
| `--layout simple` | Age-depth curves only, like the original Excel chart |
| `--sheet Sheet1` | Read only one sheet |

## Workbook layouts

Detected from the header cells (case-insensitive):

| Layout | Columns |
|---|---|
| All cores on one sheet | Row 1: core ID above each block. Row 2: `Year`, `Depth (in core)`, `Depth (adjusted)` |
| One sheet per core (sheet name = core ID) | `Year`, `Depth (cm)` |
| Varve thickness only | `Thickness (mm)`, one row per varve from the top. Pass `--top-year 2023` |

If the same core appears on several sheets, only the first one is used for the
counts. Depths are in cm.

**Correlated events (optional):** a sheet named like `Events`, in the same
layout as the counts, with two rows per event: the event's top and base in every
core (rows 1–2 = E1, rows 3–4 = E2, …). Leave a core's cells blank when the
event is missing from that core. Reading stops at the first fully blank row.

## Example

`make_example_data.py` writes a **synthetic** workbook in the all-cores-on-one-sheet
layout, with an Events sheet, so you can try the script without real core data. `example_age_depth.png`
is the result.

![example](example_age_depth.png)
