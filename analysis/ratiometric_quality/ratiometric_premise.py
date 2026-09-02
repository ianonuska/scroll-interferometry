"""Queue item 2 PREMISE TEST (declared before running), 1667 L2 locked2 field.
Current quality: w = coherence * clip(amp / global_p95, 0, 1)  -- amplitude normalised
by ONE global number, so the weight tracks beam/brightness (AoE 'loop gain' symptom).
Candidate: q_ratio = coherence * clip(amp / local_rms, 0, 1) with local_rms = RMS of
the high-passed image over ~one fringe wavelength -> dimensionless local fringe contrast.
KILL if |corr(q_ratio, coherence)| > 0.95 on-mask (it is coherence re-expressed).
PROMISING (worth the full solve + box M3 gauge) if, ranking nodes by quality,
AUC(low quality | high within-wrap error vs low) for q_ratio exceeds that for w by
>= 0.05 and reaches >= 0.60. Otherwise: not worth building, record and stop."""
import glob, re, sys, time
import numpy as np, tifffile
from scipy import ndimage as ndi
sys.path.insert(0, "/Users/ianonuska/projects/vesuvius/scroll-interferometry"); import winding_phase as wp
t0=time.time(); log=lambda m: print(f"[{time.time()-t0:5.0f}s] {m}", flush=True)
img=np.load("img_L2.npy").astype(float); mask=np.load("mask_L2.npy")
W3=np.load("W_L3.npy").astype(float); Q3=np.load("Q_L3.npy")
gy,gx=np.gradient(np.nan_to_num(W3)); g3=np.hypot(gx,gy)
g3=np.where(np.isfinite(W3)&(Q3>0.05)&(g3>1e-4),g3,np.nan); g3=np.where(np.isnan(g3),np.nanmedian(g3),g3)
g3=ndi.gaussian_filter(g3,4.0); g2=ndi.zoom(g3,(img.shape[0]/g3.shape[0],img.shape[1]/g3.shape[1]),order=1)/2.0
k_prior=2*np.pi*g2
theta,coh=wp.structure_tensor(img); kmag,amp=wp.local_frequency(img,theta,k_prior=k_prior,lock_tol=0.6)
kmag=np.clip(kmag,2*np.pi/60,2*np.pi/3); log("front end done")
lam=float(np.median(2*np.pi/kmag[mask])); log(f"median fringe wavelength {lam:.1f} px")
a95=np.percentile(amp[mask],95)+1e-9; w=np.where(mask,coh*np.clip(amp/a95,0,1),0.0)
hp=img-ndi.gaussian_filter(img,lam); local_rms=np.sqrt(ndi.gaussian_filter(hp**2,lam))+1e-9
q_ratio=np.where(mask,coh*np.clip(amp/local_rms/np.percentile((amp/local_rms)[mask],95),0,1),0.0)
m=mask&np.isfinite(w)&np.isfinite(q_ratio)
print(f"corr(q_ratio, coherence) = {np.corrcoef(q_ratio[m],coh[m])[0,1]:+.3f}   KILL if > 0.95")
print(f"corr(q_ratio, w)         = {np.corrcoef(q_ratio[m],w[m])[0,1]:+.3f}")
print(f"corr(w, coherence)       = {np.corrcoef(w[m],coh[m])[0,1]:+.3f}")
print(f"corr(w, local brightness)       = {np.corrcoef(w[m],ndi.gaussian_filter(img,lam)[m])[0,1]:+.3f}")
print(f"corr(q_ratio, local brightness) = {np.corrcoef(q_ratio[m],ndi.gaussian_filter(img,lam)[m])[0,1]:+.3f}")
W=np.load("W_L2_locked2.npy").astype(float); ny,nx=img.shape; ZF=2317*8; S=4; rows=[]
for d in sorted(glob.glob("meshes/*")):
    if not re.search(r"-w(\d+)",d) or "merged" in d: continue
    z=tifffile.imread(d+"/z.tif"); sel=(z>0)&(np.abs(z-ZF)<24)
    if sel.sum()<50: continue
    x=tifffile.imread(d+"/x.tif")[sel]/S; y=tifffile.imread(d+"/y.tif")[sel]/S
    xi=np.clip(x.astype(int),0,nx-1); yi=np.clip(y.astype(int),0,ny-1)
    wv=W[yi,xi]; ok=np.isfinite(wv)&mask[yi,xi]
    if ok.sum()<30: continue
    dev=np.abs(wv[ok]-np.median(wv[ok])); rows.append(np.stack([w[yi,xi][ok],q_ratio[yi,xi][ok],coh[yi,xi][ok],dev],1))
E=np.concatenate(rows); e=E[:,3]; hi,lo=e>0.5,e<=0.25
def auc_low(q):  # prob that a high-error node has LOWER quality than a low-error node
    qh=q[hi][::max(1,hi.sum()//600)]; ql=q[lo][::max(1,lo.sum()//600)]
    return np.mean([(a<b)+0.5*(a==b) for a in qh for b in ql])
print(f"n_hi={hi.sum()} n_lo={lo.sum()}")
for name,col in [("w (current)",0),("q_ratio",1),("coherence alone",2)]:
    q=E[:,col]; print(f"AUC(low {name:16s}| high-error vs low-error) = {auc_low(q):.3f}")
    dec=np.quantile(q,np.linspace(0,1,11)); acc=[np.mean(e[(q>=dec[i])&(q<dec[i+1]+1e-12)]>0.5) for i in range(10)]
    print("   frac high-error by quality decile (low->high):", " ".join(f"{a:.2f}" for a in acc))
np.save("q_ratio_L2.npy",q_ratio.astype(np.float32)); log("PREMISE TEST DONE")
