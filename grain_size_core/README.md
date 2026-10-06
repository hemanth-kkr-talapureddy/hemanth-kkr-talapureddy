# GUAC-29A event grain-size panels (E2 flood, E3 earthquake)

`python plot_event_grain_size.py` → `output/GUAC29A_E2_core_grain_size.{png,pdf,svg}` and `..._E3_...`

- The core photo is cropped on a true depth scale using a ruler calibration of 19.887 px/cm (max residual 0.74 mm).
- Panels: mean grain size in µm (log2) and in φ, each with a ±1σ Folk–Ward sorting envelope, plus stacked class fractions (%).
- Both figures share the same 9 cm depth span, the same x-axis limits, the same panel widths (set in inches) and 24 pt Liberation Sans. They are 575 pt tall, matching the cartoon PDFs.
- `run_checks()` runs before plotting: it checks the ruler calibration, re-detects the ruler edges, checks the crop mapping and verifies that fractions sum to 100%. It also checks that µm and φ agree, and it checks the event depths and ages.

Requires: matplotlib, pandas, numpy, pillow, openpyxl.
