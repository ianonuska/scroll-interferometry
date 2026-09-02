"""Armour B4 (declared first): the 9.6 um fringe-scale claim and the pitch-lock fix rest
on one slice (z4634). Repeat at four more stations spanning the accepted meshes' z range.
Per station: L3 solve (19.2 um) -> prior -> L2 locked2 solve (lock_tol 0.6, ridge_min_frac 0.55)
and the static L2 solve; then per-wrap fringes-per-winding (field pitch / accepted NN spacing)
for static vs locked. Claim holds if the static overcount (>1.2 median) and the locked
correction (closer to 1 than static) reproduce at >= 3 of 4 stations."""
import sys, time, glob, re, json, numpy as np, zarr, tifffile
from scipy import ndimage as ndi
sys.path.insert(0, "/Users/ianonuska/projects/vesuvius/scroll-interferometry"); import winding_phase as wp
t0 = time.time(); log = lambda m: print(f"[{time.time()-t0:5.0f}s] {m}", flush=True)
B = "https://vesuvius-challenge-open-data.s3.us-east-1.amazonaws.com"
g = zarr.open_group(f"{B}/PHerc1667/volumes/20251217075048-2.399um-0.2m-78keV-masked.zarr", mode="r")
meshes = []
for d in sorted(glob.glob("meshes/*")):
    m = re.search(r"-w(\d+)", d)
    if m and "merged" not in d: meshes.append((int(m.group(1)), d))
def per_wrap_scale(W, mask, ZF, S):
    pts = {}
    ny, nx = W.shape
    for w, d in meshes:
        z = tifffile.imread(d+"/z.tif"); sel = (z > 0) & (np.abs(z-ZF) < 24)
        if sel.sum() < 50: continue
        x = tifffile.imread(d+"/x.tif")[sel]/S; y = tifffile.imread(d+"/y.tif")[sel]/S
        xi = np.clip(x.astype(int), 0, nx-1); yi = np.clip(y.astype(int), 0, ny-1); ok = mask[yi, xi] & np.isfinite(W[yi, xi])
        if ok.sum() < 30: continue
        pts[w] = (x[ok], y[ok], np.median(W[yi, xi][ok]))
    ws = sorted(pts); ratios = []
    for a, b in zip(ws, ws[1:]):
        xa, ya, Wa = pts[a]; xb, yb, Wb = pts[b]
        # true pitch: median NN distance from wrap a points to wrap b points
        sub = np.arange(0, len(xa), max(1, len(xa)//400)); d = np.sqrt((xa[sub,None]-xb[None,::max(1,len(xb)//2000)])**2 + (ya[sub,None]-yb[None,::max(1,len(yb)//2000)])**2).min(1)
        true_px = float(np.median(d))
        if true_px >= 40: continue
        ratios.append({"pair": f"w{a}->w{b}", "true_px": true_px, "dW": float(Wb-Wa)})
    return ratios
results = {}
for ZF in (10000, 14000, 24000, 30000):
    z3, z2 = ZF//8, ZF//4
    img3 = np.asarray(g["3"][z3]).astype(float); mask3 = ndi.binary_erosion(img3 > 0, iterations=3); ys, xs = np.nonzero(mask3); core3 = (float(xs.mean()), float(ys.mean()))
    W3, Q3, th3, k3, info3 = wp.winding_coordinate(img3, mask3, core3); log(f"ZF{ZF}: L3 span {np.nanmax(W3)-np.nanmin(W3):.1f}")
    img2 = np.asarray(g["2"][z2]).astype(float); mask2 = ndi.binary_erosion(img2 > 0, iterations=6); ys, xs = np.nonzero(mask2); core2 = (float(xs.mean()), float(ys.mean()))
    W2s, Q2s, _, _, _ = wp.winding_coordinate(img2, mask2, core2); log(f"ZF{ZF}: L2 static span {np.nanmax(W2s)-np.nanmin(W2s):.1f}")
    gy, gx = np.gradient(np.nan_to_num(W3)); g3 = np.hypot(gx, gy); g3 = np.where(np.isfinite(W3)&(Q3>0.05)&(g3>1e-4), g3, np.nan); g3 = np.where(np.isnan(g3), np.nanmedian(g3), g3); g3 = ndi.gaussian_filter(g3, 4.0)
    kp = 2*np.pi*ndi.zoom(g3, (img2.shape[0]/g3.shape[0], img2.shape[1]/g3.shape[1]), order=1)/2.0
    W2l, Q2l, _, _, _ = wp.winding_coordinate(img2, mask2, core2, k_prior=kp, lock_tol=0.6, ridge_min_frac=0.55); log(f"ZF{ZF}: L2 locked span {np.nanmax(W2l)-np.nanmin(W2l):.1f}")
    rs, rl = per_wrap_scale(W2s, mask2, ZF, 4), per_wrap_scale(W2l, mask2, ZF, 4)
    # windings counted per true winding = dW (accepted pair = exactly 1 winding); >1 overcount
    ms = float(np.median([abs(r["dW"]) for r in rs])) if rs else float("nan"); ml = float(np.median([abs(r["dW"]) for r in rl])) if rl else float("nan")
    results[ZF] = {"pairs": len(rs), "static_median_dW": ms, "locked_median_dW": ml, "static": rs, "locked": rl}
    log(f"ZF{ZF}: pairs {len(rs)} | field windings per true winding: static {ms:.2f} locked {ml:.2f}")
    np.save(f"W_L2_static_ZF{ZF}.npy", W2s.astype(np.float32)); np.save(f"W_L2_locked_ZF{ZF}.npy", W2l.astype(np.float32)); np.save(f"img_L2_ZF{ZF}.npy", img2.astype(np.float32)); np.save(f"mask_L2_ZF{ZF}.npy", mask2)
json.dump(results, open("b4_slices.json", "w"), indent=1); print("B4 DONE")
