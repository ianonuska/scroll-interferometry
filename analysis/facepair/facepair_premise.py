"""Item 9.2 PREMISE (declared first), PHerc 1667 z4634, 9.6 um.
Each ridge-adjacency pair (p,q) is either the two FACES of one sheet or two different
SHEETS. Truth from accepted wraps: a ridge within R px of wrap A's mesh points belongs
to A. Pair = 'same' if both ridges belong to the same wrap; 'different' if to two
different wraps. Statistic: intensity along the straight segment p->q in the CT
(min and mean, normalised by local intensity); also the gap length.
Declared: PROMISING if AUC(same vs different) >= 0.80 for some single statistic;
KILL if the best AUC < 0.70 (distributions overlap). Everything reported."""
import sys, glob, re, json, numpy as np, tifffile
from scipy import ndimage as ndi
sys.path.insert(0, "/Users/ianonuska/projects/vesuvius/scroll-interferometry"); import winding_phase as wp
img = np.load("img_L2.npy").astype(float); mask = np.load("mask_L2.npy")
ys, xs = np.nonzero(mask); core = (float(xs.mean()), float(ys.mean()))
W3 = np.load("W_L3.npy").astype(float); Q3 = np.load("Q_L3.npy"); gy, gx = np.gradient(np.nan_to_num(W3)); g3 = np.hypot(gx, gy)
g3 = np.where(np.isfinite(W3)&(Q3>0.05)&(g3>1e-4), g3, np.nan); g3 = np.where(np.isnan(g3), np.nanmedian(g3), g3); g3 = ndi.gaussian_filter(g3, 4.0)
k_prior = 2*np.pi*ndi.zoom(g3, (img.shape[0]/g3.shape[0], img.shape[1]/g3.shape[1]), order=1)/2.0
theta, coh = wp.structure_tensor(img); kmag, amp = wp.local_frequency(img, theta, k_prior=k_prior, lock_tol=0.6); kmag = np.clip(kmag, 2*np.pi/60, 2*np.pi/3)
ny, nx = img.shape; yy, xx = np.mgrid[0:ny, 0:nx]; rx, ry = xx-core[0], yy-core[1]
nxv, nyv = np.cos(theta), np.sin(theta); flip = (nxv*rx + nyv*ry) < 0; nxv = np.where(flip, -nxv, nxv); nyv = np.where(flip, -nyv, nyv)
# use the UNGATED pair set (min_frac=0) so face pairs are present
pairs = wp.ridge_pairs(img, (nxv, nyv), amp, mask, max_gap=45, kmag=kmag, min_frac=0.0)
py, px, qy, qx = [np.asarray(a) for a in pairs[:4]]
print(f"pairs (ungated): {len(py)}")
# wrap id map from accepted meshes
ZF = 2317*8; S = 4; R = 4
wid = np.full((ny, nx), -1, int); dist = np.full((ny, nx), 1e9)
for d in sorted(glob.glob("meshes/*")):
    m = re.search(r"-w(\d+)", d)
    if not m or "merged" in d: continue
    z = tifffile.imread(d+"/z.tif"); s = (z > 0) & (np.abs(z-ZF) < 24)
    if s.sum() < 50: continue
    x = tifffile.imread(d+"/x.tif")[s]/S; y = tifffile.imread(d+"/y.tif")[s]/S
    pt = np.zeros((ny, nx), bool); pt[np.clip(y.astype(int),0,ny-1), np.clip(x.astype(int),0,nx-1)] = True
    dt = ndi.distance_transform_edt(~pt); upd = dt < dist; wid[upd] = int(m.group(1)); dist[upd] = dt[upd]
near = dist <= R
wp_, wq_ = np.where(near[py, px], wid[py, px], -1), np.where(near[qy, qx], wid[qy, qx], -1)
lab = np.where((wp_ >= 0) & (wq_ >= 0), np.where(wp_ == wq_, 1, 0), -1)   # 1 same, 0 different, -1 unknown
# consecutive-wrap 'different' only (adjacent accepted wraps), to avoid gaps across missing wraps
diff_ok = (lab == 0) & (np.abs(wp_ - wq_) == 1)
same = lab == 1
print(f"labelled pairs: same-sheet {same.sum()}, different-adjacent-sheets {diff_ok.sum()}, unknown {(lab==-1).sum()}")
loc = ndi.gaussian_filter(img, 25.0) + 1e-6
def seg_stats(i):
    n = max(int(np.hypot(qy[i]-py[i], qx[i]-px[i])), 2); t = np.linspace(0, 1, n)
    yy_ = py[i] + t*(qy[i]-py[i]); xx_ = px[i] + t*(qx[i]-px[i])
    v = ndi.map_coordinates(img, [yy_, xx_], order=1) / ndi.map_coordinates(loc, [yy_, xx_], order=1)
    return v.min(), v.mean(), v[1:-1].min() if n > 3 else v.min(), n
idx = np.nonzero(same | diff_ok)[0]
S_ = np.array([seg_stats(i) for i in idx]); y = same[idx].astype(int)
def auc(x, y):
    a, b = x[y==1], x[y==0]; a = a[::max(1,len(a)//1500)]; b = b[::max(1,len(b)//1500)]
    return np.mean([(u > v) + 0.5*(u == v) for u in a for v in b])
names = ["min intensity along segment", "mean intensity along segment", "interior min", "gap length (px)"]
res = {}
for j, nm in enumerate(names):
    A_ = auc(S_[:, j], y); res[nm] = float(A_); print(f"  AUC(same-sheet > different | {nm}) = {A_:.3f}")
best = max(res.values()); print(f"best AUC {best:.3f} -> " + ("PROMISING" if best >= 0.80 else ("KILL" if best < 0.70 else "INCONCLUSIVE")))
print(f"gap length medians: same {np.median(S_[y==1,3]):.1f} px, different {np.median(S_[y==0,3]):.1f} px")
json.dump({"n_same": int(same.sum()), "n_diff": int(diff_ok.sum()), "auc": res}, open("facepair_premise.json", "w"), indent=1); print("FACEPAIR PREMISE DONE")
