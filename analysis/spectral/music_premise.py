"""Item 10.A.1 SYNTHETIC PREMISE (declared first). Super-resolution spectral estimation
(MUSIC) claims to separate two sinusoids closer than the Fourier limit given the model
order. Our use: two sheets merged into one fringe at 19.2 um (undercount). Danger: a wrong
model order INVENTS a component. Test on synthetic 1-D profiles along the sheet normal,
window N=64 samples, unit-amplitude sinusoids with random phases, white noise at SNR 10 dB,
plus a mild chirp (spacing drifts 5% across the window, as a real scroll's pitch does):
 (1) resolution: two sinusoids at spacings p and p*r, r in {1.10, 1.25, 1.5, 2.0}; MUSIC with
     order 2 -> fraction of trials where BOTH estimated frequencies are within 8% of truth;
     compared to the FFT peak-pick with the same window (Fourier limit ~ 1/N -> r_min ~ 1.6 at p=20).
 (2) false positives: ONE sinusoid, MUSIC told order 2 -> fraction of trials where a second
     component is reported > 15% away from the true one with pseudo-spectrum peak above half
     the main peak (an 'invented sheet').
Declared: WORTH BUILDING if resolution at r=1.25 >= 0.8 AND invented-sheet rate <= 0.05.
KILL if invented-sheet rate > 0.20 at any SNR tested (10, 20 dB), regardless of resolution."""
import numpy as np, json
rng = np.random.default_rng(0); N = 64; p = 12.0; trials = 400
def profile(spacings, snr_db, chirp=0.05):
    t = np.arange(N); x = np.zeros(N)
    for s in spacings:
        ph = rng.uniform(0, 2*np.pi); k = 2*np.pi/s * (1 + chirp*(t/N - 0.5)); x += np.cos(np.cumsum(k) + ph)
    x += rng.normal(0, np.sqrt(len(spacings)/2)/10**(snr_db/20), N); return x
def music(x, order, M=32, grid=None):
    # covariance from forward-backward averaging over M-length snapshots
    X = np.lib.stride_tricks.sliding_window_view(x, M); R = X.T @ X / len(X); R = 0.5*(R + np.flipud(np.fliplr(R)).conj())
    w, V = np.linalg.eigh(R); En = V[:, :M-order]
    grid = np.linspace(2*np.pi/40, 2*np.pi/5, 1200) if grid is None else grid
    a = np.exp(1j*np.outer(np.arange(M), grid)); P = 1/np.maximum((np.abs(En.conj().T @ a)**2).sum(0), 1e-12)
    return grid, P
def peaks(grid, P, n):
    idx = [i for i in range(1, len(P)-1) if P[i] > P[i-1] and P[i] >= P[i+1]]; idx = sorted(idx, key=lambda i: -P[i])[:n]
    return [(2*np.pi/grid[i], P[i]) for i in idx]
res = {}
for snr in (10, 20):
    res[snr] = {"resolution": {}, "invented_rate": None}
    for r in (1.10, 1.25, 1.5, 2.0):
        ok_m = ok_f = 0
        for _ in range(trials):
            x = profile([p, p*r], snr); g, P = music(x, 2); est = [s for s, _ in peaks(g, P, 2)]
            def matched(est):
                if len(est) != 2: return False
                e = sorted(est); t = sorted((p, p*r)); return abs(e[0]/t[0]-1) < 0.08 and abs(e[1]/t[1]-1) < 0.08
            ok_m += matched(est)
            F = np.abs(np.fft.rfft(x*np.hanning(N), 8*N)); f = np.fft.rfftfreq(8*N); ip = [i for i in range(1, len(F)-1) if F[i] > F[i-1] and F[i] >= F[i+1]]; ip = sorted(ip, key=lambda i: -F[i])[:2]; estf = [1/f[i] for i in ip if f[i] > 0]
            ok_f += matched(estf)
        res[snr]["resolution"][r] = {"music": ok_m/trials, "fft": ok_f/trials}
        print(f"SNR {snr} dB, spacing ratio {r}: MUSIC resolves both {ok_m/trials:.2f} | FFT {ok_f/trials:.2f}", flush=True)
    inv = 0
    for _ in range(trials):
        x = profile([p], snr); g, P = music(x, 2); pk = peaks(g, P, 2)
        if len(pk) == 2 and abs(pk[1][0]/p - 1) > 0.15 and pk[1][1] > 0.5*pk[0][1]: inv += 1
    res[snr]["invented_rate"] = inv/trials; print(f"SNR {snr} dB, ONE sheet, order 2: invented second sheet in {inv/trials:.2f} of trials", flush=True)
worth = all(res[s]["resolution"][1.25]["music"] >= 0.8 for s in res) and all(res[s]["invented_rate"] <= 0.05 for s in res)
kill = any(res[s]["invented_rate"] > 0.20 for s in res)
print("VERDICT:", "KILL (invents sheets)" if kill else ("WORTH BUILDING" if worth else "NOT WORTH IT (fails resolution or FP bar)"))
json.dump(res, open("/workspace/music_premise.json", "w"), indent=1); print("MUSIC PREMISE DONE")
