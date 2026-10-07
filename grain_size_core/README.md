# GUAC-24A & GUAC-29A event grain-size panels (E2 flood, E3 earthquake)

**Colab:** open `GUAC_event_grain_size_colab.ipynb`, then choose Runtime → Run all and upload the 5 files when asked: the 29A CSV, the 24A CSV, the events Excel file, the 29A core photo (JPG) and `GUAC-24A_realscale.pdf`.
**Local:** `python GUAC_event_grain_size_colab.py` reads `data/` and writes 4 figures to `output/`.

- Panels: core photo on a true depth scale, mean grain size (φ) with a ±1σ Folk–Ward sorting envelope, and class fractions (%).
- All 4 figures share a 9 cm depth window, the same x-axis limits, the same panel widths and 24 pt Liberation Sans. They are 575 pt tall, matching the cartoon PDFs.
- Depth checks run on every execution:
  - 29A uses the 10 cm ruler blocks: error ≤ 0.8 mm (99.1 %).
  - 24A uses 109 cm and half-cm ruler ticks in the image from the real-scale PDF: error ≤ 0.12 mm (99.9 %).
- 24A sample labels are single depths of 2 mm slices and are drawn as depth ± 1 mm. The two replicates at 35.8 cm are averaged. Each figure plots only the samples that belong to its event. The mean grain size is a continuous line through the samples, and fractions are continuous. Missing samples inside an event (24A E2, 14.3–14.9 cm) are filled by linear interpolation; set `SHOW_SORTING = False` to hide the ±1σ envelope.
