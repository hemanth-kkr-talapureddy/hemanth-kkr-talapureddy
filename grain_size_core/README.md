# GUAC-29A event grain-size panels (E2 flood, E3 earthquake)

Notebook: `GUAC29A_event_grain_size_colab.ipynb` (Colab: open it, then Runtime → Run all, and upload the 3 files when asked).
Script version: `GUAC29A_event_grain_size_colab.py` (same code).

**Colab:** upload the script and run `%run GUAC29A_event_grain_size_colab.py`. When asked, upload the grain-size CSV, the events Excel file and the original core photo. The figures are shown and downloaded as a zip.
**Local:** `python GUAC29A_event_grain_size_colab.py` (it reads the inputs from `data/` and writes to `output/`).

- Panels: core photo on a true depth scale, mean grain size (φ) with a ±1σ Folk–Ward sorting envelope, and class fractions (%).
- Both figures share a 9 cm depth window, the same x-axis limits, the same panel widths and 24 pt Liberation Sans. They are 575 pt tall, matching the cartoon PDFs.
- Built-in checks run before plotting. Depth error is ≤ 0.8 mm (under one photo pixel, 0.5 mm), which is ≥ 99 % of the 9 cm window. Grain-size values are plotted exactly as they appear in the CSV.
