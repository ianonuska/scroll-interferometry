"""10.C follow-on, declared before running: accuracy of band-selection RULES against
true pitch on accepted wraps (1667 z4634, 9.6 um). Correct = k_sel/k_true in [0.7, 1.4].
Rules: (a) dominant amplitude (current default); (b) self-prior: band nearest the
local median (sigma = one fringe) of the dominant k, no external info; (c) L3 lock
(current pitch lock, prior from the 19.2 um solve); (d) ORACLE: best band per point
(upper bound of any per-point band choice). A rule is worth building only if it
closes a real fraction of the gap between (a) and (d) -- declared: >= 1/3 of the gap."""
import sys, json, glob, re, numpy as np, tifffile
from scipy import ndimage as ndi
sys.path.insert(0, "/Users/ianonuska/projects/vesuvius/scroll-interferometry"); import winding_phase as wp
img = np.load("img_L2.npy").astype(float); mask = np.load("mask_L2.npy")
theta, coh = wp.structure_tensor(img)
bands = ((0.8, 3.0), (1.6, 6.0), (3.2, 12.0), (6.0, 24.0)); K, A = [], []
for s1, s2 in bands:
    k, a = wp._band_frequency(img, theta, s1, s2); K.append(k); A.append(ndi.gaussian_filter(a, 2.0))
K = np.stack(K); A = np.stack(A); iy, ix = np.mgrid[0:img.shape[0], 0:img.shape[1]]
sel_dom = np.argmax(A, 0); k_dom = K[sel_dom, iy, ix]
lam = float(np.median(2*np.pi/np.clip(k_dom[mask], 2*np.pi/60, 2*np.pi/3)))
k_med = ndi.median_filter(np.clip(k_dom, 2*np.pi/60, 2*np.pi/3), size=int(2*lam)|1)
with np.errstate(all="ignore"):
    d_self = np.abs(np.log(K) - np.log(k_med)[None]); sel_self = np.argmin(np.where(np.isfinite(d_self), d_self, 9), 0)
# L3 lock prior exactly as locked2
W3 = np.load("W_L3.npy").astype(float); Q3 = np.load("Q_L3.npy"); gy, gx = np.gradient(np.nan_to_num(W3)); g3 = np.hypot(gx, gy)
g3 = np.where(np.isfinite(W3)&(Q3>0.05)&(g3>1e-4), g3, np.nan); g3 = np.where(np.isnan(g3), np.nanmedian(g3), g3); g3 = ndi.gaussian_filter(g3, 4.0)
kp = 2*np.pi*ndi.zoom(g3, (img.shape[0]/g3.shape[0], img.shape[1]/g3.shape[1]), order=1)/2.0
with np.errstate(all="ignore"):
    d_lock = np.abs(np.log(K) - np.log(kp)[None]); ok = np.isfinite(d_lock) & (d_lock < 0.6)
    sc = np.where(ok, A/(1+d_lock), -1); sel_lock = np.where(ok.any(0), np.argmax(sc, 0), sel_dom)
pairs = json.load(open("/Users/ianonuska/projects/vesuvius/scroll-interferometry/validation/fringe_scale/fringe_scale_by_wrap_L2.json"))
kt_by = {int(p["pair"].split("->")[0][1:]): 2*np.pi/p["true_px"] for p in pairs if p["true_px"] < 40}
ZF = 2317*8; S = 4; ny, nx = img.shape; acc = {n: [] for n in ("dominant", "self-prior", "L3 lock", "ORACLE")}; nlist = []
for d in sorted(glob.glob("meshes/*")):
    m = re.search(r"-w(\d+)", d)
    if not m or "merged" in d or int(m.group(1)) not in kt_by: continue
    z = tifffile.imread(d+"/z.tif"); s = (z > 0) & (np.abs(z-ZF) < 24)
    if s.sum() < 50: continue
    x = tifffile.imread(d+"/x.tif")[s]/S; y = tifffile.imread(d+"/y.tif")[s]/S
    xi = np.clip(x.astype(int), 0, nx-1); yi = np.clip(y.astype(int), 0, ny-1); okm = mask[yi, xi]; xi, yi = xi[okm], yi[okm]
    kt = kt_by[int(m.group(1))]; Kp = K[:, yi, xi]
    corr = lambda ks: ((ks/kt >= 0.7) & (ks/kt <= 1.4))
    acc["dominant"].append(corr(Kp[sel_dom[yi, xi], np.arange(len(xi))])); acc["self-prior"].append(corr(Kp[sel_self[yi, xi], np.arange(len(xi))]))
    acc["L3 lock"].append(corr(Kp[sel_lock[yi, xi], np.arange(len(xi))])); acc["ORACLE"].append(corr(Kp).any(0))
res = {n: float(np.concatenate(v).mean()) for n, v in acc.items()}; n = len(np.concatenate(acc["dominant"]))
gap = res["ORACLE"] - res["dominant"]
print(f"n={n} points on {len(kt_by)} wraps; fringe wavelength est {lam:.1f} px")
for k, v in res.items(): print(f"  {k:11s} correct {v:.3f}" + (f"   closes {(v-res['dominant'])/gap:.0%} of the gap to oracle" if k not in ("dominant","ORACLE") and gap>0 else ""))
json.dump(res, open("jointband_rules.json", "w"), indent=1); print("RULES DONE")
