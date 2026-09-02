"""Plan item 5 PREMISE (declared first): is the integration step where error lives?
The field is the weighted-LS integral of (k-vector field + ridge-adjacency +1 pairs).
Branch-cut / min-cost-flow integrators only help if the constraint set is INCONSISTENT
in a localised way (dislocations) that LS smears. Measure on 1667 z4634 (locked2 field):
 (a) per ridge pair, residual r = (W[q]-W[p]) - 1 winding; fraction |r| > 0.5;
 (b) plaquette residues of the k-vector field: round(closed-loop integral / 2pi) != 0
     over 8x8-px loops, density on-mask;
 (c) spatial structure of violated pairs: are they clustered (few components holding
     most violations) or diffuse?
Declared: the integrator family is WORTH testing if violated-pair fraction >= 10% AND
the top 10 clusters hold >= 50% of violations (localised inconsistency). KILL the whole
family if violations are diffuse (top-10 clusters < 20%) or rare (< 3%)."""
import sys, numpy as np
from scipy import ndimage as ndi
sys.path.insert(0, "/Users/ianonuska/projects/vesuvius/scroll-interferometry"); import winding_phase as wp
img = np.load("img_L2.npy").astype(float); mask = np.load("mask_L2.npy"); W = np.load("W_L2_locked2.npy").astype(float)
ys, xs = np.nonzero(mask); core = (float(xs.mean()), float(ys.mean()))
W3 = np.load("W_L3.npy").astype(float); Q3 = np.load("Q_L3.npy"); gy, gx = np.gradient(np.nan_to_num(W3)); g3 = np.hypot(gx, gy)
g3 = np.where(np.isfinite(W3)&(Q3>0.05)&(g3>1e-4), g3, np.nan); g3 = np.where(np.isnan(g3), np.nanmedian(g3), g3); g3 = ndi.gaussian_filter(g3, 4.0)
k_prior = 2*np.pi*ndi.zoom(g3, (img.shape[0]/g3.shape[0], img.shape[1]/g3.shape[1]), order=1)/2.0
theta, coh = wp.structure_tensor(img); kmag, amp = wp.local_frequency(img, theta, k_prior=k_prior, lock_tol=0.6); kmag = np.clip(kmag, 2*np.pi/60, 2*np.pi/3)
ny, nx = img.shape; yy, xx = np.mgrid[0:ny, 0:nx]; rx, ry = xx-core[0], yy-core[1]
nxv, nyv = np.cos(theta), np.sin(theta); flip = (nxv*rx + nyv*ry) < 0; nxv = np.where(flip, -nxv, nxv); nyv = np.where(flip, -nyv, nyv)
kx, ky = kmag*nxv, kmag*nyv
pairs = wp.ridge_pairs(img, (nxv, nyv), amp, mask, max_gap=45, kmag=kmag, min_frac=0.55)
py, px, qy, qx = [np.asarray(a) for a in pairs[:4]] if len(pairs) >= 4 else (None,)*4
if py is None:
    print("pairs format:", type(pairs), [getattr(p, "shape", None) for p in pairs]); sys.exit()
r = (W[qy, qx] - W[py, px]) - 1.0; ok = np.isfinite(r); r = r[ok]
viol = np.abs(r) > 0.5
print(f"(a) ridge pairs {len(r)}: residual median {np.median(r):+.3f}, |r|>0.5 on {viol.mean():.1%}, |r|>0.25 on {(np.abs(r)>0.25).mean():.1%}")
print(f"    residual percentiles 5/25/50/75/95: {np.percentile(r,[5,25,50,75,95]).round(2)}")
# (b) plaquette residues of the k field, loop size L
L = 8; res_map = np.zeros((ny//L, nx//L)); 
for i in range(ny//L):
    for j in range(nx//L):
        y0, x0 = i*L, j*L
        if not mask[y0:y0+L, x0:x0+L].all(): continue
        top = kx[y0, x0:x0+L].sum(); right = ky[y0:y0+L, x0+L-1].sum(); bot = -kx[y0+L-1, x0:x0+L].sum(); left = -ky[y0:y0+L, x0].sum()
        res_map[i, j] = np.round((top+right+bot+left)/(2*np.pi))
onm = mask[::L, ::L][:ny//L, :nx//L]
print(f"(b) k-field plaquette residues ({L}px loops): nonzero on {(res_map[onm]!=0).mean():.2%} of on-mask loops")
# (c) clustering of violated pairs
vm = np.zeros((ny, nx), bool); my, mx = ((py+qy)//2)[ok][viol], ((px+qx)//2)[ok][viol]; vm[my, mx] = True
vm = ndi.binary_dilation(vm, iterations=6); lab, n = ndi.label(vm); sizes = ndi.sum(np.ones_like(lab), lab, range(1, n+1))
cnt = np.bincount(lab[my, mx], minlength=n+1)[1:]; top10 = np.sort(cnt)[::-1][:10].sum()/max(1, cnt.sum())
print(f"(c) violated pairs form {n} clusters (6px dilation); top-10 clusters hold {top10:.1%} of violations")
verdict = "WORTH TESTING" if (viol.mean() >= 0.10 and top10 >= 0.5) else ("KILL (diffuse)" if top10 < 0.2 else ("KILL (rare)" if viol.mean() < 0.03 else "INCONCLUSIVE"))
print("VERDICT:", verdict); print("RESIDUE PREMISE DONE")
