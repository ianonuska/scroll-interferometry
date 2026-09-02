"""Item 9.2b PREMISE (declared first): striation matching. For a ridge pair (p,q), sample
the CT intensity along each ridge's tangent (perpendicular to the sheet normal) for
+/- L px, high-pass, and correlate the two profiles. Two faces of one sheet should
share fibre texture along the sheet (same fibres seen from both sides); two different
sheets should not. Same labels as facepair_premise.py. Declared: PROMISING if
AUC(same > different | correlation) >= 0.80 for some L in {16, 32, 64}; KILL if < 0.70."""
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
pairs = wp.ridge_pairs(img, (nxv, nyv), amp, mask, max_gap=45, kmag=kmag, min_frac=0.0); py, px, qy, qx = [np.asarray(a) for a in pairs[:4]]
ZF = 2317*8; S = 4; R = 4; wid = np.full((ny, nx), -1, int); dist = np.full((ny, nx), 1e9)
for d in sorted(glob.glob("meshes/*")):
    m = re.search(r"-w(\d+)", d)
    if not m or "merged" in d: continue
    z = tifffile.imread(d+"/z.tif"); s = (z > 0) & (np.abs(z-ZF) < 24)
    if s.sum() < 50: continue
    x = tifffile.imread(d+"/x.tif")[s]/S; y = tifffile.imread(d+"/y.tif")[s]/S
    pt = np.zeros((ny, nx), bool); pt[np.clip(y.astype(int),0,ny-1), np.clip(x.astype(int),0,nx-1)] = True
    dt = ndi.distance_transform_edt(~pt); upd = dt < dist; wid[upd] = int(m.group(1)); dist[upd] = dt[upd]
near = dist <= R; wp_ = np.where(near[py, px], wid[py, px], -1); wq_ = np.where(near[qy, qx], wid[qy, qx], -1)
same = (wp_ >= 0) & (wq_ >= 0) & (wp_ == wq_); diff = (wp_ >= 0) & (wq_ >= 0) & (np.abs(wp_ - wq_) == 1)
idx = np.nonzero(same | diff)[0]; rng = np.random.default_rng(0); idx = rng.choice(idx, size=min(6000, len(idx)), replace=False); y = same[idx].astype(int)
hp = img - ndi.gaussian_filter(img, 3.0)
def profile(y0, x0, L):
    tx, ty = -nyv[y0, x0], nxv[y0, x0]  # tangent (perpendicular to normal)
    t = np.arange(-L, L+1); return ndi.map_coordinates(hp, [y0 + t*ty, x0 + t*tx], order=1)
def auc(x, y):
    a, b = x[y==1], x[y==0]; a = a[::max(1,len(a)//1500)]; b = b[::max(1,len(b)//1500)]
    return float(np.mean([(u > v) + 0.5*(u == v) for u in a for v in b]))
res = {}
for L in (16, 32, 64):
    c = []
    for i in idx:
        a = profile(py[i], px[i], L); b = profile(qy[i], qx[i], L)
        a = a - a.mean(); b = b - b.mean(); c.append((a*b).sum()/(np.sqrt((a*a).sum()*(b*b).sum())+1e-9))
    c = np.array(c); res[L] = auc(c, y)
    print(f"  L={L:2d}: corr medians same {np.median(c[y==1]):+.3f} different {np.median(c[y==0]):+.3f} | AUC(same>diff) = {res[L]:.3f}")
best = max(res.values()); print(f"best AUC {best:.3f} -> " + ("PROMISING" if best >= 0.80 else ("KILL" if best < 0.70 else "INCONCLUSIVE")))
json.dump({str(k): v for k, v in res.items()}, open("facepair_striation.json", "w"), indent=1); print("STRIATION PREMISE DONE")
