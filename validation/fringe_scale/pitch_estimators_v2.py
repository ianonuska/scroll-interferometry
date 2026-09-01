#!/usr/bin/env python3
"""Item 3 / 10.A.2 gate: pitch estimators on PHerc1667 z2317 (L3, 19.2um).

TRUTH for the gate: the validated W field. In an annulus, pitch = median of
1/|grad W| (px per winding). This is the number the old estimator missed.

Estimators, all textbook:
  A. autocorrelation first-peak (the one that failed: returned 13.9 px)
  B. cepstrum peak (Bogert/Healy/Tukey 1963)
  C. Radon-Wigner / fractional-Fourier concentration: for a polar strip the
     radial profile is a chirp; find the chirp rate + center frequency that
     maximally concentrate energy. Implemented as: dechirp with candidate
     rate alpha, FFT, take peak sharpness; report 1/f_peak as pitch.
Gate (declared): |estimate - truth| / truth < 0.10 in at least 4 of 5 annuli.
"""
import numpy as np
from scipy import ndimage as ndi

W = np.load("W_L3.npy").astype(np.float64); Q = np.load("Q_L3.npy"); img = np.load("img_L3.npy").astype(np.float64)
mask = np.load("mask_L3.npy")
ys, xs = np.nonzero(mask); core = (xs.mean(), ys.mean())
ny, nx = W.shape
yy, xx = np.mgrid[0:ny, 0:nx]
r = np.hypot(xx - core[0], yy - core[1])

# polar image (2048 theta x r)
n_th = 2048; r_max = int(np.percentile(r[mask], 99)); rs = np.arange(8, r_max, 1.0)
ths = np.linspace(0, 2*np.pi, n_th, endpoint=False)
PX = core[0] + rs[None, :]*np.cos(ths[:, None]); PY = core[1] + rs[None, :]*np.sin(ths[:, None])
P = ndi.map_coordinates(img, [PY.ravel(), PX.ravel()], order=1, cval=0.0).reshape(n_th, len(rs))
M = ndi.map_coordinates(mask.astype(float), [PY.ravel(), PX.ravel()], order=1, cval=0.0).reshape(n_th, len(rs)) > 0.5
PW = ndi.map_coordinates(np.nan_to_num(W), [PY.ravel(), PX.ravel()], order=1, cval=np.nan).reshape(n_th, len(rs))
PQ = ndi.map_coordinates(Q, [PY.ravel(), PX.ravel()], order=1, cval=0.0).reshape(n_th, len(rs))

def truth_pitch(r0, r1):
    sl = (rs >= r0) & (rs < r1)
    dW = np.gradient(PW, axis=1)[:, sl]
    ok = M[:, sl] & (PQ[:, sl] > 0.1) & np.isfinite(dW) & (np.abs(dW) > 1e-4)
    return float(np.median(1.0/np.abs(dW[ok]))) if ok.sum() > 100 else np.nan

def est_autocorr(prof):
    prof = prof - prof.mean(); ac = np.correlate(prof, prof, "full")[len(prof)-1:]
    ac /= ac[0] + 1e-12
    # first local max after the first zero crossing (fixes the lag-3 floor bug)
    zc = np.argmax(ac < 0) if np.any(ac < 0) else 3
    seg = ac[zc:zc+200]
    return float(zc + np.argmax(seg))

def est_cepstrum(prof):
    prof = prof - prof.mean()
    spec = np.abs(np.fft.rfft(prof * np.hanning(len(prof))))**2 + 1e-12
    cep = np.fft.irfft(np.log(spec))
    lo, hi = 6, min(200, len(cep)//2)
    return float(lo + np.argmax(cep[lo:hi]))

def est_dechirp(prof, alphas=np.linspace(-0.004, 0.004, 81)):
    """Radon-Wigner-style: multiply by conjugate chirp exp(-i*pi*alpha*n^2),
    FFT, pick alpha with the sharpest peak; pitch = 1/f_peak at that alpha."""
    prof = prof - prof.mean(); n = np.arange(len(prof)); best = (-1, None, None)
    win = np.hanning(len(prof))
    for a in alphas:
        x = prof * win * np.exp(-1j*np.pi*a*n**2)
        S = np.abs(np.fft.fft(x, 4*len(prof)))**2
        f = np.fft.fftfreq(4*len(prof))
        band = (f > 1/200) & (f < 1/5)
        k = np.argmax(S[band]); peak = S[band][k]; sharp = peak / (S[band].mean() + 1e-12)
        if sharp > best[0]:
            best = (sharp, f[band][k], a)
    return float(1.0/best[1]), float(best[2])


def detrend(prof, sigma=40.0):
    return prof - ndi.gaussian_filter1d(prof, sigma)

annuli = [(r0, r0+120) for r0 in range(60, 660, 120)]
N_WEDGE = 64
print(f"{'annulus':>12} {'truth':>7} {'autocorr':>9} {'cepstrum':>9} {'dechirp':>8}  (medians over wedges; n wedges used)")
hits = {"autocorr": 0, "cepstrum": 0, "dechirp": 0}
for r0, r1 in annuli:
    sl = (rs >= r0) & (rs < r1)
    t = truth_pitch(r0, r1)
    ests = {"autocorr": [], "cepstrum": [], "dechirp": []}
    for w in range(N_WEDGE):
        th_sl = slice(w*n_th//N_WEDGE, (w+1)*n_th//N_WEDGE)
        blk = P[th_sl][:, sl]; mb = M[th_sl][:, sl]
        if mb.mean() < 0.95: continue
        prof = detrend(blk.mean(axis=0))
        ests["autocorr"].append(est_autocorr(prof))
        ests["cepstrum"].append(est_cepstrum(prof))
        ests["dechirp"].append(est_dechirp(prof)[0])
    n = len(ests["dechirp"])
    if n < 8: print(f"{r0:4d}-{r1:<4d}   {t:7.1f}   (only {n} wedges inside mask)"); continue
    med = {k: float(np.median(v)) for k, v in ests.items()}
    for k, v in med.items():
        if np.isfinite(t) and abs(v - t)/t < 0.10: hits[k] += 1
    print(f"{r0:4d}-{r1:<4d}   {t:7.1f} {med['autocorr']:9.1f} {med['cepstrum']:9.1f} {med['dechirp']:8.1f}   n={n}")
print("GATE (<10% err in >=4/5 annuli):", {k: f"{v}/5 {'PASS' if v>=4 else 'fail'}" for k, v in hits.items()})
