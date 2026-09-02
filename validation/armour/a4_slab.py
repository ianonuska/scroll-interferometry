"""Armour A4 (declared first): adjacent-slice agreement of independent solves vs the
z-coupled slab solve (winding_slab.py), at 9.6 um around z4634 (3 slices, z 4633-4635).
Metric: adjacent_consistency = median |dW| between adjacent slices after gauge alignment.
Report both; the claim 'coupling tightens the residual' holds only if slab median is lower
at both adjacent pairs."""
import sys, time, json, numpy as np, zarr
from scipy import ndimage as ndi
sys.path.insert(0, "/Users/ianonuska/projects/vesuvius/scroll-interferometry"); import winding_slab as ws
t0 = time.time(); log = lambda m: print(f"[{time.time()-t0:5.0f}s] {m}", flush=True)
B = "https://vesuvius-challenge-open-data.s3.us-east-1.amazonaws.com"
g = zarr.open_group(f"{B}/PHerc1667/volumes/20251217075048-2.399um-0.2m-78keV-masked.zarr", mode="r")
imgs, masks, cores = [], [], []
import os
for z in (4633, 4634, 4635):
    f = "img_L2.npy" if z == 4634 else f"img_L2_z{z}.npy"
    img = None
    if os.path.exists(f): img = np.load(f).astype(float)
    else:
        for attempt in range(8):
            try:
                g = zarr.open_group(f"{B}/PHerc1667/volumes/20251217075048-2.399um-0.2m-78keV-masked.zarr", mode="r")
                img = np.asarray(g["2"][z]).astype(float); break
            except Exception as e: log(f"fetch z{z} attempt {attempt} failed: {str(e)[:80]}"); time.sleep(30)
        if img is None: raise SystemExit(f"could not fetch z{z}")
        np.save(f, img.astype(np.float32))
    mask = ndi.binary_erosion(img > 0, iterations=6); ys, xs = np.nonzero(mask)
    imgs.append(img); masks.append(mask); cores.append((float(xs.mean()), float(ys.mean())))
log("slices loaded")
ind = [ws.solve_independent(imgs[s], masks[s], cores[s])/(2*np.pi) for s in range(3)]; log("independent solves done")
ci = ws.adjacent_consistency(ind, masks); log(f"independent adjacent consistency (median, p90 |dW|): {ci}")
slab, sinfo = ws.solve_slab(imgs, masks, cores); slab = [p/(2*np.pi) for p in slab]; log(f"slab info {sinfo}"); log("slab solve done")
cs = ws.adjacent_consistency(slab, masks); log(f"slab adjacent consistency: {cs}")
json.dump({"independent": ci, "slab": cs}, open("a4_slab.json", "w"), indent=1); print("A4 DONE")
