"""Item 9.7 PREMISE (declared first): does the field's angular phase advance smoothly with z,
so one station can seed the next? For each 9.6 um pitch-locked station (ZF 10000, 14000,
18536, 24000, 30000), take the fractional part of W around a ring of radius R about the core
and find the angle where W crosses an integer (the 'spoke' direction); unwrap across stations
and fit angle vs z. Also the winding count on the ring (span) per station. VIABLE if the
linear fit of spoke angle vs z has RMS residual < 15 deg over the five stations (so a
station's phase predicts the next within a fraction of a winding); otherwise propagation
needs a per-station solve anyway and the idea is shelved."""
import numpy as np, json
st = {10000: "W_L2_locked_ZF10000.npy", 14000: "W_L2_locked_ZF14000.npy", 18536: "W_L2_locked2.npy", 24000: "W_L2_locked_ZF24000.npy", 30000: "W_L2_locked_ZF30000.npy"}
mk = {10000: "mask_L2_ZF10000.npy", 14000: "mask_L2_ZF14000.npy", 18536: "mask_L2.npy", 24000: "mask_L2_ZF24000.npy", 30000: "mask_L2_ZF30000.npy"}
rows = []
for z in sorted(st):
    W = np.load(st[z]).astype(float); m = np.load(mk[z]); ys, xs = np.nonzero(m); cy, cx = ys.mean(), xs.mean()
    ny, nx = W.shape; yy, xx = np.mgrid[0:ny, 0:nx]; r = np.hypot(yy-cy, xx-cx); th = np.degrees(np.arctan2(yy-cy, xx-cx))
    for R in (150, 250, 350):
        ring = m & np.isfinite(W) & (np.abs(r-R) < 2)
        if ring.sum() < 200: continue
        w = W[ring]; t = th[ring]; o = np.argsort(t); w, t = w[o], t[o]
        # spoke = angle where fractional part crosses 0, taken as circular mean of exp(i*2pi*frac) phase
        frac = w - np.floor(w); ph = np.degrees(np.angle(np.mean(np.exp(1j*2*np.pi*frac)))); coh = abs(np.mean(np.exp(1j*2*np.pi*frac)))
        rows.append((z, R, ph, coh, float(w.max()-w.min())))
        print(f"z{z} R{R}: phase {ph:7.1f} deg, coherence {coh:.2f}, windings around ring {w.max()-w.min():.2f}")
res = {}
for R in (150, 250, 350):
    pts = [(z, ph) for z, r_, ph, c, sp in rows if r_ == R]
    if len(pts) < 4: continue
    zs = np.array([p[0] for p in pts]); ph = np.array([p[1] for p in pts]); ph = np.degrees(np.unwrap(np.radians(ph)))
    A = np.vstack([zs, np.ones_like(zs)]).T; coef, *_ = np.linalg.lstsq(A, ph, rcond=None); resid = ph - A@coef
    res[R] = {"deg_per_1000z": float(coef[0]*1000), "rms_resid_deg": float(np.sqrt(np.mean(resid**2)))}
    print(f"R{R}: twist {coef[0]*1000:+.1f} deg per 1000 z, RMS residual {np.sqrt(np.mean(resid**2)):.1f} deg -> " + ("VIABLE" if np.sqrt(np.mean(resid**2)) < 15 else "not viable"))
json.dump(res, open("twist_premise.json", "w"), indent=1); print("TWIST PREMISE DONE")
