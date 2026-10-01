# Grain-size down-core plots at true (cm) scale

`grain_size_core_plots.py` plots pre-calculated D10, D50, D90, Mean and Sorting
for two sediment cores, each parameter as a **separate figure**, with the depth
axis at **real scale** (default 1 cm core = 1 cm on the figure), so plots can be
placed directly beside the core photo.

## Colab
1. Paste the script into a cell.
2. Edit the `USER SETTINGS` block: file names, `COLUMN_MAP`, core top/bottom and
   your exact `events` depths `(top_cm, bottom_cm, "label")`.
3. Run — Colab asks you to upload the data files and downloads a zip of results.

Data format: see `core1_grainsize_TEMPLATE.csv` (CSV or Excel; column names are
mapped in `COLUMN_MAP`).

## Output (`grain_size_plots/`)
- `<core>/full_core/` — whole core, event intervals shaded
- `<core>/event_<label>/` — exactly the event top→bottom depth
- SVG (editable), PDF, 600-dpi PNG for each; `plot_scale_summary.csv`

## Aligning with the core photo
The data panel starts `MARGIN_TOP_CM` (1.3 cm) below the top edge of each file
and its height equals (bottom − top) × `PRINT_CM_PER_CORE_CM`. Import the SVG/PDF
at 100 % into Illustrator/Inkscape and scale the photo so that its cm ruler
matches. Do not resize the plots. Set `PRINT_CM_PER_CORE_CM = 0.5` for a 1:2 layout.
