"""Per-wrap fringe scale on 1667 z2317 (L3): true adjacent-wrap spacing from
the accepted segments vs the field's implied pitch 1/|grad W| at the same
points. ratio = W_pitch / true_pitch (1.0 = one fringe per winding)."""
import glob, re, json, numpy as np, tifffile
from scipy.spatial import cKDTree
from scipy import ndimage as ndi
ZF = 2317*8; S = 8
W = np.load("W_L3.npy").astype(float); Q = np.load("Q_L3.npy")
gy, gx = np.gradient(np.nan_to_num(W)); g = np.hypot(gx, gy)
pts = {}
for d in sorted(glob.glob("meshes/*")):
    m = re.search(r"-w(\d+)", d)
    if not m: continue
    z = tifffile.imread(d + "/z.tif"); sel = (z > 0) & (np.abs(z - ZF) < 24)
    if sel.sum() < 50: continue
    pts[int(m.group(1))] = np.stack([tifffile.imread(d + "/x.tif")[sel], tifffile.imread(d + "/y.tif")[sel]], 1) / S
rows = []
print(f"{'pair':>8} {'true px':>8} {'W-pitch':>8} {'ratio':>6} {'n':>6}")
for a in sorted(pts):
    if a + 1 not in pts: continue
    t = cKDTree(pts[a]); dd, _ = t.query(pts[a+1])
    P = pts[a+1]; xi = np.clip(P[:, 0].astype(int), 0, W.shape[1]-1); yi = np.clip(P[:, 1].astype(int), 0, W.shape[0]-1)
    ok = (Q[yi, xi] > 0.1) & (g[yi, xi] > 1e-4) & np.isfinite(W[yi, xi])
    if ok.sum() < 30: continue
    tp = float(np.median(dd[ok])); wp = float(np.median(1.0 / g[yi, xi][ok]))
    rows.append({"pair": f"w{a}->w{a+1}", "true_px": tp, "W_pitch_px": wp, "ratio": wp/tp, "n": int(ok.sum())})
    print(f"w{a:02d}->w{a+1:02d} {tp:8.1f} {wp:8.1f} {wp/tp:6.2f} {ok.sum():6d}")
r = np.array([x["ratio"] for x in rows]); tp = np.array([x["true_px"] for x in rows])
print(f"median ratio {np.median(r):.2f}; ratio where true<12px: {np.median(r[tp<12]):.2f} (n={int((tp<12).sum())}); where true>=12px: {np.median(r[tp>=12]):.2f} (n={int((tp>=12).sum())})")
json.dump(rows, open("fringe_scale_by_wrap_L3.json", "w"), indent=1)
