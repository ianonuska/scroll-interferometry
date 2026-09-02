"""Armour A5 (declared first): repeat the fringe-scale and band-ceiling measurements on a
SECOND scroll with per-wrap ground truth at the PRIZE resolution: PHerc 0172, 7.91 um
(volume 20241024131838), 44 accepted per-winding segments w052-w095. Stations z=6000, 9000.
Per station: level-1 solve (15.8 um) -> prior -> level-0 pitch-locked solve, and level-0
static solve. Metrics: (i) field windings per true winding per adjacent accepted pair
(static vs locked); (ii) band-selection accuracy: dominant / lock / ORACLE against the true
pitch. Claim armoured if the 1667 pattern (static >1.2 overcount, lock closer to 1, oracle
far above the rules) reproduces at both stations."""
import sys, time, glob, re, json, numpy as np, zarr, tifffile
from scipy import ndimage as ndi
sys.path.insert(0, "/workspace/wp_current"); import winding_phase as wp
t0 = time.time(); log = lambda m: print(f"[{time.time()-t0:5.0f}s] {m}", flush=True)
B = "https://vesuvius-challenge-open-data.s3.us-east-1.amazonaws.com"
g = zarr.open_group(f"{B}/PHerc0172/volumes/20241024131838-7.910um-53keV-masked.zarr", mode="r")
meshes = [(int(re.search(r"-w(\d+)", d).group(1)), d) for d in sorted(glob.glob("/workspace/p0172/meshes/*")) if re.search(r"-w(\d+)", d)]
def mesh_points(ZF, S, ny, nx, mask):
    pts = {}
    for w, d in meshes:
        z = tifffile.imread(d+"/z.tif"); sel = (z > 0) & (np.abs(z-ZF) < 6)
        if sel.sum() < 50: continue
        x = tifffile.imread(d+"/x.tif")[sel]/S; y = tifffile.imread(d+"/y.tif")[sel]/S
        xi = np.clip(x.astype(int), 0, nx-1); yi = np.clip(y.astype(int), 0, ny-1); ok = mask[yi, xi]
        if ok.sum() < 30: continue
        pts[w] = (x[ok], y[ok], xi[ok], yi[ok])
    return pts
def true_pitch(pts):
    ws = sorted(pts); out = {}
    for a, b in zip(ws, ws[1:]):
        xa, ya, _, _ = pts[a]; xb, yb, _, _ = pts[b]
        sa = slice(0, None, max(1, len(xa)//400)); sb = slice(0, None, max(1, len(xb)//3000))
        d = np.sqrt((xa[sa,None]-xb[None,sb])**2 + (ya[sa,None]-yb[None,sb])**2).min(1); tp = float(np.median(d))
        if tp < 60: out[a] = (b, tp)
    return out
results = {}
for ZF in (6000, 9000):
    img1 = np.asarray(g["1"][ZF//2]).astype(float); m1 = ndi.binary_erosion(img1 > 0, iterations=3); ys, xs = np.nonzero(m1); c1 = (float(xs.mean()), float(ys.mean()))
    W1, Q1, _, _, _ = wp.winding_coordinate(img1, m1, c1); log(f"z{ZF}: L1 {img1.shape} span {np.nanmax(W1)-np.nanmin(W1):.1f}")
    img0 = np.asarray(g["0"][ZF]).astype(float); m0 = ndi.binary_erosion(img0 > 0, iterations=6); ys, xs = np.nonzero(m0); c0 = (float(xs.mean()), float(ys.mean()))
    gy, gx = np.gradient(np.nan_to_num(W1)); g1 = np.hypot(gx, gy); g1 = np.where(np.isfinite(W1)&(Q1>0.05)&(g1>1e-4), g1, np.nan); g1 = np.where(np.isnan(g1), np.nanmedian(g1), g1); g1 = ndi.gaussian_filter(g1, 4.0)
    kp = 2*np.pi*ndi.zoom(g1, (img0.shape[0]/g1.shape[0], img0.shape[1]/g1.shape[1]), order=1)/2.0
    Ws, Qs, th, ks, _ = wp.winding_coordinate(img0, m0, c0); log(f"z{ZF}: L0 static {img0.shape} span {np.nanmax(Ws)-np.nanmin(Ws):.1f}")
    Wl, Ql, _, kl, _ = wp.winding_coordinate(img0, m0, c0, k_prior=kp, lock_tol=0.6, ridge_min_frac=0.55); log(f"z{ZF}: L0 locked span {np.nanmax(Wl)-np.nanmin(Wl):.1f}")
    ny, nx = img0.shape; pts = mesh_points(ZF, 1, ny, nx, m0); tp = true_pitch(pts); log(f"z{ZF}: {len(pts)} wraps on slice, {len(tp)} adjacent pairs with true pitch")
    dWs = [abs(np.nanmedian(Ws[pts[b][3], pts[b][2]]) - np.nanmedian(Ws[pts[a][3], pts[a][2]])) for a, (b, _) in tp.items()]
    dWl = [abs(np.nanmedian(Wl[pts[b][3], pts[b][2]]) - np.nanmedian(Wl[pts[a][3], pts[a][2]])) for a, (b, _) in tp.items()]
    # band oracle
    theta, coh = wp.structure_tensor(img0); K, A = [], []
    for s1, s2 in ((0.8, 3.0), (1.6, 6.0), (3.2, 12.0), (6.0, 24.0)):
        k, a = wp._band_frequency(img0, theta, s1, s2); K.append(np.clip(k, 2*np.pi/80, 2*np.pi/2)); A.append(ndi.gaussian_filter(a, 2.0))
    K = np.stack(K); A = np.stack(A); sel_dom = np.argmax(A, 0)
    with np.errstate(all="ignore"):
        dl = np.abs(np.log(K) - np.log(np.clip(kp, 2*np.pi/80, 2*np.pi/2))[None]); ok = dl < 0.6; sc = np.where(ok, A/(1+dl), -1); sel_lock = np.where(ok.any(0), np.argmax(sc, 0), sel_dom)
    corr = {"dominant": [], "lock": [], "oracle": []}
    for a, (b, tpx) in tp.items():
        _, _, xi, yi = pts[a]; kt = 2*np.pi/tpx; Kp = K[:, yi, xi]; ar = np.arange(len(xi))
        inr = lambda ks: (ks/kt >= 0.7) & (ks/kt <= 1.4)
        corr["dominant"].append(inr(Kp[sel_dom[yi, xi], ar])); corr["lock"].append(inr(Kp[sel_lock[yi, xi], ar])); corr["oracle"].append(inr(Kp).any(0))
    acc = {k: float(np.concatenate(v).mean()) for k, v in corr.items()}; n = len(np.concatenate(corr["dominant"]))
    results[ZF] = {"pairs": len(tp), "true_pitch_px_median": float(np.median([t for _, t in tp.values()])), "static_median_dW": float(np.median(dWs)), "locked_median_dW": float(np.median(dWl)), "band_points": n, "band_accuracy": acc}
    log(f"z{ZF}: true pitch median {results[ZF]['true_pitch_px_median']:.1f} px | windings per true winding: static {np.median(dWs):.2f} locked {np.median(dWl):.2f} | band correct: dominant {acc['dominant']:.3f} lock {acc['lock']:.3f} ORACLE {acc['oracle']:.3f} (n={n})")
    np.save(f"/workspace/p0172/W_static_z{ZF}.npy", Ws.astype(np.float32)); np.save(f"/workspace/p0172/W_locked_z{ZF}.npy", Wl.astype(np.float32)); np.save(f"/workspace/p0172/mask_z{ZF}.npy", m0)
json.dump(results, open("/workspace/p0172/a5_0172.json", "w"), indent=1); print("A5 DONE")
