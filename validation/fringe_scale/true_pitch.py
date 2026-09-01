"""Ground-truth pitch per annulus on 1667 z2317 from the ACCEPTED single-winding
segments (2.399um frame; /8 -> L3 px). For each adjacent pair of winding
numbers present in the z-slab, nearest-neighbour distance from wrap w+1 points
to wrap w points; bucket by radius from the L3 mask centroid (in L3 px)."""
import glob, os, re, json, numpy as np, tifffile
from scipy.spatial import cKDTree
ZF = 2317*8; SCALE = 8
mask = np.load("mask_L3.npy"); ys, xs = np.nonzero(mask); core = (xs.mean(), ys.mean())
pts = {}
for d in sorted(glob.glob("meshes/*")):
    m = re.search(r"-w(\d+)", d); 
    if not m: continue
    z = tifffile.imread(d + "/z.tif"); sel = (z > 0) & (np.abs(z - ZF) < 24)
    if sel.sum() < 50: continue
    x = tifffile.imread(d + "/x.tif")[sel]; y = tifffile.imread(d + "/y.tif")[sel]
    pts[int(m.group(1))] = np.stack([x, y], 1) / SCALE   # L3 px
winds = sorted(pts); out = []
for a, b in zip(winds, winds[1:]):
    if b - a != 1: continue
    t = cKDTree(pts[a]); dd, _ = t.query(pts[b])
    r = np.hypot(pts[b][:, 0] - core[0], pts[b][:, 1] - core[1])
    out.append(np.stack([r, dd], 1))
out = np.concatenate(out)
annuli = [(r0, r0+120) for r0 in range(60, 660, 120)]
truth = {}
for r0, r1 in annuli:
    s = (out[:, 0] >= r0) & (out[:, 0] < r1)
    truth[f"{r0}-{r1}"] = float(np.median(out[s, 1])) if s.sum() > 30 else None
    print(f"{r0:4d}-{r1:<4d} true pitch {truth[f'{r0}-{r1}']} px (n={int(s.sum())})")
json.dump(truth, open("true_pitch_L3.json", "w"), indent=1)
print("adjacent-wrap pairs used:", [(a, a+1) for a in winds if a+1 in pts])
