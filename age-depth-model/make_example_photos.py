"""
Write SYNTHETIC uncropped core photos (example_photos/CORE-0x.jpg) and the
ruler readings for them (example_photo_depths.csv), for trying
plot_core_photos.py without real photos.

Each photo is laid out like a real core photo: background around the core, a
ruler beside the core liner with 1 cm ticks (numbered every 10 cm), and the
ruler's 0 at the top of the liner, some way below the photo's top edge. The
sediment and the event layers are drawn at the in-core depths of
example_varve_counts.xlsx (run make_example_data.py first). CORE-03 is saved
sideways (core top on the left); CORE-04 stores its resolution (DPI).
"""

from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont

from plot_age_depth import read_cores, split_varves_events

PX_PER_CM = 50  # = 127 DPI, a whole number as stored in JPEG files


def draw_core_photo(sediment_top, events, length_cm, px_per_cm=PX_PER_CM, margin_top_cm=6.0,
                    margin_bottom_cm=4.0, seed=0):
    """An upright, uncropped core photo. Depths are in cm from the liner top
    (ruler 0). Returns (PIL image, pixel row of ruler 0)."""
    rng = np.random.default_rng(seed)
    ppc = px_per_cm
    h = int(round((margin_top_cm + length_cm + margin_bottom_cm) * ppc))
    w = int(round(14 * ppc))
    zero = int(round(margin_top_cm * ppc))  # pixel row of ruler 0 / liner top
    bg = np.array([70, 78, 90], float)
    img = np.ones((h, w, 3)) * bg * (1 + rng.normal(0, 0.02, (h, w, 1)))

    # core liner and sediment, 6.5 cm wide, from x = 5 cm
    x0, x1 = int(5 * ppc), int(11.5 * ppc)
    rows = np.arange(zero, zero + int(round(length_cm * ppc)))
    depth = (rows - zero) / ppc
    shade = np.full(len(rows), 0.93)  # empty liner above the sediment
    sed = depth >= sediment_top
    shade[sed] = 0.42 + 0.05 * np.sin(depth[sed] * 2 * np.pi / 0.25) + rng.normal(0, 0.02, sed.sum())
    for t, b in events:  # event layers: paler
        shade[(depth >= t) & (depth < b)] = 0.66
    core = np.stack([shade * 0.95, shade * 0.85, shade * 0.70], axis=1) * 255
    img[rows[0]:rows[-1] + 1, x0:x1] = core[:, None, :] * (1 + rng.normal(0, 0.02, (len(rows), x1 - x0, 1)))
    img[rows[0]:rows[-1] + 1, x0:x0 + 8] *= 0.7  # liner walls
    img[rows[0]:rows[-1] + 1, x1 - 8:x1] *= 0.7
    photo = Image.fromarray(np.clip(img, 0, 255).astype(np.uint8))

    # ruler: white strip from x = 1.5 to 4 cm, ticks every 1 cm from 0
    d = ImageDraw.Draw(photo)
    rx0, rx1 = int(1.5 * ppc), int(4 * ppc)
    d.rectangle([rx0, zero - int(0.5 * ppc), rx1, zero + int((length_cm + 0.5) * ppc)],
                fill=(238, 236, 228))
    font = ImageFont.load_default(size=int(0.55 * ppc))
    for cm in range(0, int(length_cm) + 1):
        y = zero + cm * ppc
        tick = 0.6 if cm % 10 == 0 else 0.4 if cm % 5 == 0 else 0.25
        d.line([rx1 - int(tick * ppc * 2.5), y, rx1, y], fill=(0, 0, 0), width=2)
        if cm % 10 == 0:
            d.text((rx0 + 4, y), str(cm), fill=(0, 0, 0), font=font, anchor="lm")
    return photo, zero


if __name__ == "__main__":
    out = Path("example_photos")
    out.mkdir(exist_ok=True)
    cores = read_cores("example_varve_counts.xlsx", depth_col="in core")
    rows = []
    for i, (name, df) in enumerate(cores.items()):
        _, evt = split_varves_events(df)
        length = float(np.ceil(df["depth"].max() + 3))
        photo, zero = draw_core_photo(df["depth"].min(), zip(evt["top"], evt["base"]), length,
                                      margin_top_cm=4 + 2 * i, seed=i)
        ruler_100 = zero + 100 * PX_PER_CM if length >= 100 else zero + 80 * PX_PER_CM
        cm_2 = 100.0 if length >= 100 else 80.0
        row = {"Core": name, "Ref1 cm": 0.0, "Ref1 px": zero, "Ref2 cm": cm_2,
               "Ref2 px": ruler_100, "Top side": ""}
        save = {"quality": 90}
        if name == "CORE-03":  # sideways: core top on the left, readings are x coordinates
            photo = photo.rotate(90, expand=True)
            row["Top side"] = "left"
        if name == "CORE-04":  # stores its resolution: one ruler reading is enough
            save["dpi"] = (PX_PER_CM * 2.54,) * 2
            row.update({"Ref2 cm": None, "Ref2 px": None})
        photo.save(out / f"{name}.jpg", **save)
        rows.append(row)
    pd.DataFrame(rows).to_csv("example_photo_depths.csv", index=False)
    print(f"Wrote {len(rows)} synthetic uncropped photos to {out}/ and example_photo_depths.csv")
