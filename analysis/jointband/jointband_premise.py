"""Item 10.C PREMISE TEST (declared before running), PHerc 1667 z4634 at 9.6 um.
Question: where the dominant Riesz band is a face-split harmonic (k_dom ~ 2 x true
winding wavenumber), does a lower band with k ~ k_true exist with competitive
amplitude? If yes, joint multi-band integer resolution has the information to fix
the overcount; if the lower band is absent/negligible, the bands do not carry it.
Truth: per adjacent accepted-wrap pair, true pitch = NN spacing (fringe_scale_by_wrap_L2.json),
pairs restricted to plausible adjacency (true_px < 40). Points on wrap A get k_true from A->B.
Declared: PROMISING if >= 50% of face-split points (k_dom/k_true in [1.6, 2.6]) have a band
within 0.35 log-units of k_true with amplitude >= 0.3 x A_dom. KILL if < 25%.
Also reported: same statistic on non-split points (should be high too -> those are fine either way)."""
import sys, json, glob, re, numpy as np, tifffile
from scipy import ndimage as ndi
sys.path.insert(0, "/Users/ianonuska/projects/vesuvius/scroll-interferometry"); import winding_phase as wp
img = np.load("img_L2.npy").astype(float); mask = np.load("mask_L2.npy")
theta, coh = wp.structure_tensor(img)
bands = ((0.8, 3.0), (1.6, 6.0), (3.2, 12.0), (6.0, 24.0))
K, A = [], []
for s1, s2 in bands:
    k, a = wp._band_frequency(img, theta, s1, s2); K.append(k); A.append(ndi.gaussian_filter(a, 2.0))
K = np.stack(K); A = np.stack(A); sel = np.argmax(A, axis=0)
pairs = json.load(open("/Users/ianonuska/projects/vesuvius/scroll-interferometry/validation/fringe_scale/fringe_scale_by_wrap_L2.json"))
ktrue_by_wrap = {int(p["pair"].split("->")[0][1:]): 2*np.pi/p["true_px"] for p in pairs if p["true_px"] < 40}
print("wraps with true pitch:", sorted(ktrue_by_wrap), flush=True)
ZF = 2317*8; S = 4; ny, nx = img.shape; rows = []
for d in sorted(glob.glob("meshes/*")):
    m = re.search(r"-w(\d+)", d)
    if not m or "merged" in d: continue
    w = int(m.group(1))
    if w not in ktrue_by_wrap: continue
    z = tifffile.imread(d+"/z.tif"); s = (z > 0) & (np.abs(z-ZF) < 24)
    if s.sum() < 50: continue
    x = tifffile.imread(d+"/x.tif")[s]/S; y = tifffile.imread(d+"/y.tif")[s]/S
    xi = np.clip(x.astype(int), 0, nx-1); yi = np.clip(y.astype(int), 0, ny-1); ok = mask[yi, xi]
    xi, yi = xi[ok], yi[ok]; kt = ktrue_by_wrap[w]
    kd = K[sel[yi, xi], yi, xi]; ad = A[sel[yi, xi], yi, xi]
    lk = np.log(K[:, yi, xi] + 1e-9); dist = np.abs(lk - np.log(kt)); amp_ok = A[:, yi, xi] >= 0.3*ad[None]
    has_true_band = ((dist < 0.35) & amp_ok).any(axis=0)
    ratio = kd/kt
    rows.append(np.stack([np.full(len(xi), w), ratio, has_true_band.astype(float), ad], 1))
E = np.concatenate(rows); ratio, has = E[:, 1], E[:, 2] > 0
split = (ratio >= 1.6) & (ratio <= 2.6); fine = (ratio >= 0.7) & (ratio <= 1.4); under = ratio < 0.7
print(f"points {len(E)} | face-split {split.mean():.2%} | ~correct {fine.mean():.2%} | under {under.mean():.2%}")
print(f"P(true-pitch band present, amp>=0.3 A_dom | face-split) = {has[split].mean():.3f}   (promising >= 0.50, kill < 0.25)")
print(f"P(true-pitch band present | ~correct) = {has[fine].mean():.3f}")
for w in sorted(set(E[:,0].astype(int))):
    m = E[:,0]==w; print(f"  w{w:02d}: n={m.sum():5d} split {split[m].mean():.2f} correct {fine[m].mean():.2f} | true-band|split {has[m&split].mean() if (m&split).any() else float('nan'):.2f}")
print("JOINTBAND PREMISE DONE")
