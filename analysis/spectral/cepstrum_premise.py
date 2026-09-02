"""Item 10.B ZERO-CROSSING variant, same harness and bar as 10.A.2: period = 2 * window length / number of sign changes of the detrended profile.
(original header follows)
"""
"""Item 10.A.2 PREMISE (declared first): does the cepstrum recover the true winding pitch
where the quadrature/autocorrelation estimators lock to a harmonic? Field-frame test on
PHerc 1667 z4634 (9.6 um): at points on accepted wraps (true pitch known per adjacent pair),
take a 96-px profile along the sheet normal, detrend, Hann, cepstrum = IFFT(log|FFT|^2);
pitch = quefrency of the largest peak in [8, 60] px. Compare per point to (a) the static
dominant-band estimate and (b) the pitch-locked estimate. Correct = within [0.7, 1.4] x true.
Declared: WORTH TESTING FURTHER if cepstrum correct fraction >= lock (0.44) + 0.10 = 0.54;
KILL if below the static dominant band (0.37)."""
import sys, json, glob, re, numpy as np, tifffile
from scipy import ndimage as ndi
sys.path.insert(0, "/Users/ianonuska/projects/vesuvius/scroll-interferometry"); import winding_phase as wp
img = np.load("img_L2.npy").astype(float); mask = np.load("mask_L2.npy"); ny, nx = img.shape
ys, xs = np.nonzero(mask); core = (float(xs.mean()), float(ys.mean()))
theta, coh = wp.structure_tensor(img)
yy, xx = np.mgrid[0:ny, 0:nx]; rx, ry = xx-core[0], yy-core[1]; nxv, nyv = np.cos(theta), np.sin(theta); flip = (nxv*rx + nyv*ry) < 0; nxv = np.where(flip, -nxv, nxv); nyv = np.where(flip, -nyv, nyv)
pairs = json.load(open("/Users/ianonuska/projects/vesuvius/scroll-interferometry/validation/fringe_scale/fringe_scale_by_wrap_L2.json"))
tp = {int(p["pair"].split("->")[0][1:]): p["true_px"] for p in pairs if p["true_px"] < 40}
ZF = 2317*8; S = 4; L = 48; win = np.hanning(2*L+1); t = np.arange(-L, L+1); rows = []
rng = np.random.default_rng(0)
for d in sorted(glob.glob("meshes/*")):
    m = re.search(r"-w(\d+)", d)
    if not m or int(m.group(1)) not in tp or "merged" in d: continue
    z = tifffile.imread(d+"/z.tif"); s = (z > 0) & (np.abs(z-ZF) < 24)
    if s.sum() < 50: continue
    x = tifffile.imread(d+"/x.tif")[s]/S; y = tifffile.imread(d+"/y.tif")[s]/S
    xi = np.clip(x.astype(int), 0, nx-1); yi = np.clip(y.astype(int), 0, ny-1); ok = mask[yi, xi]; xi, yi = xi[ok], yi[ok]
    sel = rng.choice(len(xi), min(300, len(xi)), replace=False); truth = tp[int(m.group(1))]
    for i in sel:
        px, py = xi[i], yi[i]; prof = ndi.map_coordinates(img, [py + t*nyv[py, px], px + t*nxv[py, px]], order=1)
        prof = (prof - prof.mean()) * win; F = np.abs(np.fft.rfft(prof, 4*len(prof)))**2; c = np.fft.irfft(np.log(F + 1e-6))
        q = np.arange(len(c)) / 4.0; band = (q >= 8) & (q <= 60); cep = q[band][np.argmax(c[band])]
        p0 = prof - ndi.gaussian_filter(prof, 2.0); zc = np.count_nonzero(np.diff(np.sign(p0)) != 0); zcp = 2*len(prof)/max(zc, 1)
        rows.append((truth, cep, zcp))
E = np.array(rows); ratio = E[:, 1]/E[:, 0]; correct = ((ratio >= 0.7) & (ratio <= 1.4)).mean()
rz = E[:, 2]/E[:, 0]; zc_correct = ((rz >= 0.7) & (rz <= 1.4)).mean(); print(f"zero-crossing period: correct {zc_correct:.3f} | half pitch {((rz >= 0.4) & (rz < 0.7)).mean():.3f} | above 1.4x {(rz > 1.4).mean():.3f} -> " + ("WORTH" if zc_correct >= 0.54 else ("KILL" if zc_correct < 0.37 else "NOT BETTER")))
harm2 = ((ratio >= 0.4) & (ratio < 0.7)).mean(); sub = (ratio > 1.4).mean()
print(f"n={len(E)} | cepstrum correct {correct:.3f} | at half pitch (harmonic) {harm2:.3f} | above 1.4x {sub:.3f}")
print("reference from jointband_rules: dominant 0.369, lock 0.440, oracle 0.778")
print("VERDICT:", "WORTH TESTING FURTHER" if correct >= 0.54 else ("KILL" if correct < 0.37 else "NOT BETTER THAN LOCK"))
json.dump({"n": int(len(E)), "cepstrum_correct": float(correct), "half_pitch": float(harm2), "above_1.4": float(sub), "zero_crossing_correct": float(zc_correct)}, open("cepstrum_premise.json", "w"), indent=1); print("CEPSTRUM PREMISE DONE")
