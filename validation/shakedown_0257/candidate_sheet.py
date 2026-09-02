import sys, csv, json, numpy as np, tifffile
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
P="/workspace/first-light-pherc0826/runbook/out/PHerc0257/shakedown1/predictions"; A="/workspace/first-light-pherc0826/runbook/out/PHerc0257/shakedown1/audit"
FL="/workspace/first-light-pherc0826/runbook/out/PHerc0257/shakedown1/spiral-fit/baseline/meshes/fitted_scoped_w010-065/concat/w010-065_flat"
thr=json.load(open("/workspace/first-light-pherc0826/analysis/control-PHerc0257/calibration.json"))["forward"]["median_ink"]
F=tifffile.imread(P+"/segment.tif"); R=tifffile.imread(P+"/segment_reverse.tif")
if F.ndim==3: F=F[0]
if R.ndim==3: R=R[0]
F=F.astype(np.float32); R=R.astype(np.float32)
if F.max()>1.5: F/=255; R/=255
rows=[r for r in csv.DictReader(open(A+"/w010-065_candidates.csv")) if r["passes_size_floor"]=="True"]
X=np.stack([tifffile.imread(f"{FL}/{a}.tif").astype(float) for a in ("x","y","z")],-1); s_flat=20.0
eu=(X[:,1:]-X[:,:-1])[:-1]; ev=(X[1:,:]-X[:-1,:])[:,:-1]; lu=np.linalg.norm(eu,axis=-1); lv=np.linalg.norm(ev,axis=-1)
area=np.linalg.norm(np.cross(eu,ev),axis=-1)/s_flat**2; ang=np.degrees(np.arccos(np.clip((eu*ev).sum(-1)/(lu*lv+1e-9),-1,1)))
out=[]; n=len(rows); cols=5; rws=int(np.ceil(n/cols))
fig,axs=plt.subplots(rws,cols,figsize=(cols*4.2,rws*2.6),dpi=90); axs=np.array(axs).reshape(-1)
for i,r in enumerate(rows):
    y0,x0,h,w=int(r["y0"]),int(r["x0"]),int(r["h"]),int(r["w"]); cy,cx=y0+h//2,x0+w//2
    both=(F[y0:y0+h,x0:x0+w]>=thr)&(R[y0:y0+h,x0:x0+w]>=thr); ys,xs=np.nonzero(both)
    c=np.cov(np.stack([xs,ys]).astype(float)); ev_,evec=np.linalg.eigh(c); v=evec[:,-1]; theta=abs(np.degrees(np.arctan2(v[1],v[0]))); theta=min(theta,180-theta)
    ci,cj=min(cy//20,area.shape[0]-1),min(cx//20,area.shape[1]-1); ar=float(np.nanmedian(area[max(0,ci-2):ci+3,max(0,cj-2):cj+3])); sh=float(np.nanmedian(np.abs(ang[max(0,ci-2):ci+3,max(0,cj-2):cj+3]-90)))
    seam=bool(sh>20 or ar<0.6 or ar>1.6)
    r.update({"angle_from_horizontal_deg":round(theta,1),"local_area_ratio":round(ar,3),"local_shear_deg":round(sh,1),"seam_or_distorted":seam}); out.append(r)
    pad=120; a=axs[i]; a.imshow(F[max(0,cy-pad):cy+pad, max(0,cx-3*pad):cx+3*pad],cmap="gray",vmin=0,vmax=1,aspect="auto"); a.set_xticks([]); a.set_yticks([])
    a.set_title("#%s %.2fmm asp%.1f ang%.0f p%.2f area%.2f sh%.0f%s" % (r["id"], float(r["long_axis_mm"]), float(r["aspect"]), theta, float(r["mean_p_fwd"]), ar, sh, " SEAM" if seam else ""), fontsize=7, color="#b03a2e" if seam else "black")
for j in range(len(rows),len(axs)): axs[j].axis("off")
plt.tight_layout(); plt.savefig(A+"/w010-065_candidate_sheet.png")
with open(A+"/w010-065_candidates_annotated.csv","w",newline="") as f:
    wtr=csv.DictWriter(f,fieldnames=list(out[0].keys())); wtr.writeheader(); wtr.writerows(out)
ns=sum(1 for r in out if r["seam_or_distorted"]); nh=sum(1 for r in out if float(r["angle_from_horizontal_deg"])<=20)
print("passing %d | flagged seam/distorted %d | within 20 deg of horizontal %d | horizontal AND not flagged %d" % (len(out), ns, nh, sum(1 for r in out if float(r["angle_from_horizontal_deg"])<=20 and not r["seam_or_distorted"])))
print("SHEET DONE")
