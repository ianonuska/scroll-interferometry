"""Item 8 prototype: ambiguity map by constraint resampling on 1667 L2 (locked
config). Ambiguity = per-pixel std of the winding coordinate across B solves
that each drop 30% of the ridge-adjacency constraints (global offset removed).
Declared criteria: (a) KILL if corr(ambiguity, quality map) ~ 1 — same info,
new name; (b) PASS if ambiguity predicts where the field is wrong against the
accepted segments: AUC(ambiguity, |within-wrap deviation| > 0.5 vs <= 0.25) >= 0.7."""
import glob, re, sys, time
import numpy as np, tifffile
from scipy import ndimage as ndi
sys.path.insert(0, "/Users/ianonuska/projects/vesuvius/scroll-interferometry")
import winding_phase as wp
t0 = time.time(); log = lambda m: print(f"[{time.time()-t0:5.0f}s] {m}", flush=True)
img = np.load("img_L2.npy").astype(float); mask = np.load("mask_L2.npy")
ys, xs = np.nonzero(mask); core = (float(xs.mean()), float(ys.mean()))
# prior exactly as in locked2_L2.py
W3 = np.load("W_L3.npy").astype(float); Q3 = np.load("Q_L3.npy")
gy, gx = np.gradient(np.nan_to_num(W3)); g3 = np.hypot(gx, gy)
g3 = np.where(np.isfinite(W3) & (Q3 > 0.05) & (g3 > 1e-4), g3, np.nan); g3 = np.where(np.isnan(g3), np.nanmedian(g3), g3)
g3 = ndi.gaussian_filter(g3, 4.0)
g2 = ndi.zoom(g3, (img.shape[0]/g3.shape[0], img.shape[1]/g3.shape[1]), order=1) / 2.0
k_prior = 2*np.pi*g2
# front end once (mirrors winding_coordinate)
theta, coh = wp.structure_tensor(img)
kmag, amp = wp.local_frequency(img, theta, k_prior=k_prior, lock_tol=0.6)
kmag = np.clip(kmag, 2*np.pi/60, 2*np.pi/3)
ny, nx = img.shape; yy, xx = np.mgrid[0:ny, 0:nx]; rx, ry = xx-core[0], yy-core[1]
nxv, nyv = np.cos(theta), np.sin(theta); flip = (nxv*rx + nyv*ry) < 0
nxv = np.where(flip, -nxv, nxv); nyv = np.where(flip, -nyv, nyv)
kx, ky = kmag*nxv, kmag*nyv
a95 = np.percentile(amp[mask], 95) + 1e-9; w = np.where(mask, coh*np.clip(amp/a95, 0, 1), 0.0)
pairs = wp.ridge_pairs(img, (nxv, nyv), amp, mask, max_gap=45, kmag=kmag, min_frac=0.55)
n = len(pairs[0]); log(f"front end done, {n} pairs")
def solve(keep):
    p = tuple(a[keep] for a in pairs)
    phi, info = wp.integrate_phase_multiscale(kx, ky, w, mask, pairs=p)
    return phi/(2*np.pi)
W_full = solve(np.ones(n, bool)); log("full solve done")
rng = np.random.default_rng(0); B = 8; devs = []
for b in range(B):
    keep = rng.random(n) > 0.30
    Wb = solve(keep); d = Wb - W_full; d -= np.nanmedian(d[mask]); devs.append(d); log(f"bootstrap {b+1}/{B}")
A = np.nanstd(np.stack(devs), axis=0)                     # ambiguity map (windings)
np.save("ambiguity_L2.npy", A.astype(np.float32)); np.save("W_L2_amb_full.npy", W_full.astype(np.float32))
# (a) kill test vs quality
m = mask & np.isfinite(A) & np.isfinite(W_full)
print(f"corr(ambiguity, quality w) = {np.corrcoef(A[m], w[m])[0,1]:+.3f}  (kill if ~1)")
print(f"corr(ambiguity, log quality) = {np.corrcoef(A[m], np.log(w[m]+1e-3))[0,1]:+.3f}")
# (b) spatial prediction vs accepted segments: within-wrap deviation of W_full
ZF = 2317*8; S = 4; amb_err = []
for d in sorted(glob.glob("meshes/*")):
    mm = re.search(r"-w(\d+)", d)
    if not mm or "merged" in d: continue
    z = tifffile.imread(d+"/z.tif"); sel = (z > 0) & (np.abs(z-ZF) < 24)
    if sel.sum() < 50: continue
    x = tifffile.imread(d+"/x.tif")[sel]/S; y = tifffile.imread(d+"/y.tif")[sel]/S
    xi = np.clip(x.astype(int), 0, nx-1); yi = np.clip(y.astype(int), 0, ny-1)
    wv = W_full[yi, xi]; av = A[yi, xi]; ok = np.isfinite(wv) & np.isfinite(av) & mask[yi, xi]
    if ok.sum() < 30: continue
    dev = np.abs(wv[ok] - np.median(wv[ok])); amb_err.append(np.stack([av[ok], dev], 1))
E = np.concatenate(amb_err); a, e = E[:, 0], E[:, 1]
hi, lo = e > 0.5, e <= 0.25
auc = np.mean([(ah > al) + 0.5*(ah == al) for ah in a[hi][::max(1, hi.sum()//400)] for al in a[lo][::max(1, lo.sum()//400)]])
print(f"spatial prediction: n_hi={hi.sum()} n_lo={lo.sum()}  AUC(ambiguity | high-error vs low-error) = {auc:.3f}  (pass >= 0.7)")
print(f"ambiguity percentiles (windings): p50 {np.nanpercentile(A[m],50):.2f} p90 {np.nanpercentile(A[m],90):.2f} p99 {np.nanpercentile(A[m],99):.2f}")
print("AMBIGUITY PROTOTYPE DONE")
