"""
Write SYNTHETIC core photos (example_photos/CORE-0x.jpg) and a photo-depth
table (example_photo_depths.csv) for trying plot_core_photos.py without real
photos. Laminations and event layers are drawn at the in-core depths of
example_varve_counts.xlsx (run make_example_data.py first). CORE-03 is saved
sideways (core top on the left) to show that rotated photos are handled.
"""

from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image

from plot_age_depth import read_cores, split_varves_events

PX_PER_CM = 40
WIDTH_CM = 6.5
rng = np.random.default_rng(7)
out = Path("example_photos")
out.mkdir(exist_ok=True)

cores = read_cores("example_varve_counts.xlsx", depth_col="in core")
rows = []
for name, df in cores.items():
    top_cm, bottom_cm = 0.0, float(np.ceil(df["depth"].max() + 3))
    h, w = int((bottom_cm - top_cm) * PX_PER_CM), int(WIDTH_CM * PX_PER_CM)
    depth = np.arange(h) / PX_PER_CM + top_cm
    shade = np.full(h, 0.93)  # empty liner above the sediment
    sed = depth >= df["depth"].min()
    shade[sed] = 0.42 + 0.05 * np.sin(depth[sed] * 2 * np.pi / 0.25) + rng.normal(0, 0.02, sed.sum())
    var, evt = split_varves_events(df)
    for _, e in evt.iterrows():  # event layers: paler, sandier
        shade[(depth >= e["top"]) & (depth < e["base"])] = 0.62 + rng.normal(0, 0.03)
    rgb = np.stack([shade * 0.95, shade * 0.85, shade * 0.70], axis=1)
    img = np.repeat(rgb[:, None, :], w, axis=1) * (1 + rng.normal(0, 0.03, (h, w, 1)))
    img[:, :int(0.04 * w)] *= 0.8  # liner edges
    img[:, -int(0.04 * w):] *= 0.8
    photo = Image.fromarray(np.clip(img * 255, 0, 255).astype(np.uint8))
    if name == "CORE-03":
        photo = photo.rotate(90, expand=True)  # core top now on the left
    photo.save(out / f"{name}.jpg", quality=90)
    rows.append({"Core": name, "Top (cm)": top_cm, "Bottom (cm)": bottom_cm})

pd.DataFrame(rows).to_csv("example_photo_depths.csv", index=False)
print(f"Wrote {len(rows)} synthetic photos to {out}/ and example_photo_depths.csv")
