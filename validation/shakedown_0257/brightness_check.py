"""Does the target's ink probability track the raw render brightness? (artifact test)
Compare, on the same 1/8 subsample: corr(pred, mean intensity over the 21 central slices),
and the brightness percentile at the 50 candidates vs the whole surface. Do the same for the
control, where real ink exists, for reference."""
import csv, json, numpy as np, tifffile, zarr
P="/workspace/first-light-pherc0826/runbook/out/PHerc0257/shakedown1"; C="/workspace/first-light-pherc0826/analysis/control-PHerc0257"
def load(zpath):
    g=zarr.open_group(zpath,mode="r"); a=g["0"]; return a
def stats(pred_tif, zpath, cands=None, ds=8):
    F=tifffile.imread(pred_tif); F=F[0] if F.ndim==3 else F; F=F.astype(np.float32); F=F/255 if F.max()>1.5 else F
    a=load(zpath); nz=a.shape[0]; mid=slice(nz//2-10, nz//2+11)
    I=np.asarray(a[mid, ::ds, ::ds]).astype(np.float32).mean(0); Fd=F[::ds, ::ds][:I.shape[0], :I.shape[1]]; I=I[:Fd.shape[0], :Fd.shape[1]]
    v=(I>0)&np.isfinite(Fd); print("  corr(pred, brightness) = %.3f on %d px" % (np.corrcoef(Fd[v], I[v])[0,1], v.sum()))
    ranks=None
    if cands:
        pct=[]
        for r in cands:
            cy,cx=(int(r["y0"])+int(r["h"])//2)//ds,(int(r["x0"])+int(r["w"])//2)//ds
            loc=I[max(0,cy-3):cy+4, max(0,cx-8):cx+9]; loc=loc[loc>0]
            if len(loc): pct.append((I[v] < np.median(loc)).mean())
        pct=np.array(pct); print("  brightness percentile at candidates: median %.2f, p25 %.2f, p75 %.2f (0.5 = typical)" % (np.median(pct), np.percentile(pct,25), np.percentile(pct,75)))
    hi=Fd[v]>=np.percentile(Fd[v],99); print("  brightness percentile of the top-1%% prediction pixels: %.2f" % (I[v] < np.median(I[v][hi])).mean())
cands=[r for r in csv.DictReader(open(P+"/audit/w010-065_candidates_annotated.csv")) if r["seam_or_distorted"]=="False"]
print("TARGET (no known ink):"); stats(P+"/predictions/segment.tif", P+"/render/segment.zarr", cands)
print("CONTROL (real ink present):"); stats(C+"/control_prediction.tif", C+"/control.zarr")
print("BRIGHTNESS CHECK DONE")
