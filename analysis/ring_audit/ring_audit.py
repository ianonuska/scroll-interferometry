"""Item 10.A.3 CONFOUND AUDIT (declared first): ring artifacts in the CT are concentric
about the reconstruction centre (image centre), i.e. fake circular 'windings'. Remove them
by the classic polar median destripe (subtract the angular median of the high-passed image
at each radius about the IMAGE centre), re-solve the pitch-locked 9.6 um field on 1667 z4634,
and measure how much the field changes. MATERIAL if median |dW| (after gauge alignment)
> 0.25 winding or the per-wrap windings-per-true-winding median moves by > 0.05; otherwise
the ring confound is immaterial at this station and prior numbers stand."""
import sys, glob, re, json, numpy as np, tifffile
from scipy import ndimage as ndi
sys.path.insert(0, "/Users/ianonuska/projects/vesuvius/scroll-interferometry"); import winding_phase as wp
img = np.load("img_L2.npy").astype(float); mask = np.load("mask_L2.npy"); ny, nx = img.shape
yy, xx = np.mgrid[0:ny, 0:nx]; cy, cx = (ny-1)/2, (nx-1)/2; r = np.hypot(yy-cy, xx-cx); th = np.arctan2(yy-cy, xx-cx)
hp = img - ndi.gaussian_filter(img, 6.0)
rb = np.round(r).astype(int); ring = np.zeros(rb.max()+1)
for i in range(rb.max()+1):
    sel = (rb == i) & (img > 0)
    if sel.sum() > 50: ring[i] = np.median(hp[sel])
ring_img = ring[rb]; ring_img[img <= 0] = 0
print(f"ring profile: rms {np.sqrt(np.mean(ring**2)):.2f}, max |ring| {np.abs(ring).max():.2f}, image high-pass rms {hp[img>0].std():.2f}")
img_d = img - ring_img
ys, xs = np.nonzero(mask); core = (float(xs.mean()), float(ys.mean()))
W3 = np.load("W_L3.npy").astype(float); Q3 = np.load("Q_L3.npy"); gy, gx = np.gradient(np.nan_to_num(W3)); g3 = np.hypot(gx, gy)
g3 = np.where(np.isfinite(W3)&(Q3>0.05)&(g3>1e-4), g3, np.nan); g3 = np.where(np.isnan(g3), np.nanmedian(g3), g3); g3 = ndi.gaussian_filter(g3, 4.0)
kp = 2*np.pi*ndi.zoom(g3, (ny/g3.shape[0], nx/g3.shape[1]), order=1)/2.0
Wd, Qd, _, _, _ = wp.winding_coordinate(img_d, mask, core, k_prior=kp, lock_tol=0.6, ridge_min_frac=0.55)
W0 = np.load("W_L2_locked2.npy").astype(float)
both = mask & np.isfinite(W0) & np.isfinite(Wd); d = (Wd - W0)[both]; d = d - np.median(d)
print(f"field change after destripe: median |dW| {np.median(np.abs(d)):.3f}, p90 {np.percentile(np.abs(d),90):.3f}, frac >0.5 {(np.abs(d)>0.5).mean():.3f}")
pairs = json.load(open("/Users/ianonuska/projects/vesuvius/scroll-interferometry/validation/fringe_scale/fringe_scale_by_wrap_L2.json"))
ZF = 2317*8; S = 4; med = {}
for d_ in sorted(glob.glob("meshes/*")):
    m = re.search(r"-w(\d+)", d_)
    if not m or "merged" in d_: continue
    z = tifffile.imread(d_+"/z.tif"); s = (z > 0) & (np.abs(z-ZF) < 24)
    if s.sum() < 50: continue
    x = tifffile.imread(d_+"/x.tif")[s]/S; y = tifffile.imread(d_+"/y.tif")[s]/S
    xi = np.clip(x.astype(int), 0, nx-1); yi = np.clip(y.astype(int), 0, ny-1); ok = mask[yi, xi]
    med[int(m.group(1))] = (np.nanmedian(W0[yi, xi][ok]), np.nanmedian(Wd[yi, xi][ok]))
ws = sorted(med); s0, sd = [], []
for a, b in zip(ws, ws[1:]):
    tp = [p for p in pairs if p["pair"] == f"w{a}->w{b}"]
    if tp and tp[0]["true_px"] < 40: s0.append(abs(med[b][0]-med[a][0])); sd.append(abs(med[b][1]-med[a][1]))
m0, md = float(np.median(s0)), float(np.median(sd))
print(f"windings per true winding (median over {len(s0)} pairs): original {m0:.3f}, destriped {md:.3f}, change {md-m0:+.3f}")
verdict = "MATERIAL" if (np.median(np.abs(d)) > 0.25 or abs(md-m0) > 0.05) else "IMMATERIAL"
json.dump({"ring_rms": float(np.sqrt(np.mean(ring**2))), "hp_rms": float(hp[img>0].std()), "median_abs_dW": float(np.median(np.abs(d))), "p90_abs_dW": float(np.percentile(np.abs(d),90)), "scale_original": m0, "scale_destriped": md, "verdict": verdict}, open("ring_audit.json","w"), indent=1)
print("VERDICT:", verdict); print("RING AUDIT DONE")
