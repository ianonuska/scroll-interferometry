#!/usr/bin/env python3
"""End-to-end demo: stream one PHerc 1667 cross-section from the Vesuvius
Challenge open-data bucket, solve the winding coordinate, export spiral-fitting
point-collection JSONs, and render the showcase figure.

Runs on a laptop in ~2 minutes; downloads only the CT chunks it touches
(~40 MB), never the full 8.6 TB volume.

    pip install numpy scipy zarr s3fs matplotlib pillow
    python examples/run_demo.py
"""
import sys, time, pathlib

import numpy as np
import zarr
from scipy import ndimage as ndi

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
from winding_phase import winding_coordinate
from winding_to_vc import export

STORE = ("https://vesuvius-challenge-open-data.s3.us-east-1.amazonaws.com/"
         "PHerc1667/volumes/20251217075048-2.399um-0.2m-78keV-masked.zarr")
LEVEL, SCALE = "3", 8            # 19.2 µm working resolution

def main():
    print("streaming slice from S3 …")
    t0 = time.time()
    g = zarr.open_group(STORE, mode="r")
    a = g[LEVEL]
    z = a.shape[0] // 2
    img = np.asarray(a[z]).astype(np.float64)
    print(f"  slice {img.shape} @ level {LEVEL} in {time.time()-t0:.1f}s")

    mask = ndi.binary_erosion(img > 0, iterations=3)
    ys, xs = np.nonzero(mask)
    core = (float(xs.mean()), float(ys.mean()))

    print("solving winding coordinate …")
    t0 = time.time()
    W, q, theta, kmag, info = winding_coordinate(img, mask, core)
    print(f"  solved in {time.time()-t0:.0f}s (cg info={info}); "
          f"winding span ≈ {np.nanmax(W)-np.nanmin(W):.1f}")

    export(W, q, z_index=z, scale=SCALE, out_prefix=f"pherc1667_z{z}")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(8, 8))
        ax.imshow(img, cmap="gray", vmax=np.percentile(img[img > 0], 99))
        ax.contour(np.ma.masked_invalid(W),
                   levels=np.arange(np.floor(np.nanmin(W)), np.ceil(np.nanmax(W)) + 1),
                   cmap="plasma", linewidths=0.8)
        ax.set_title("automated winding coordinate — each contour = one winding")
        ax.axis("off")
        fig.savefig(f"pherc1667_z{z}_windings.png", dpi=130, bbox_inches="tight")
        print(f"wrote pherc1667_z{z}_windings.png")
    except ImportError:
        pass

if __name__ == "__main__":
    main()
