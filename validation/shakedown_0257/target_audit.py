"""Mechanical application of the preregistered criteria to the target predictions.
Inputs: forward tif, reverse tif, calibration.json (threshold = forward median on
control ink), pixel size um, output prefix. Outputs: candidates.csv, summary.json,
crops of every candidate that passes the size floor (for the human audit)."""
import sys, json, numpy as np, tifffile
from scipy import ndimage as ndi
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
pf, pr, cal_p, um, out = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4]), sys.argv[5]
thr = json.load(open(cal_p))["forward"]["median_ink"]
F = tifffile.imread(pf).astype(np.float32); R = tifffile.imread(pr).astype(np.float32)
if F.ndim == 3: F = F[0]
if R.ndim == 3: R = R[0]
if F.max() > 1.5: F /= 255; R /= 255
sub = (slice(None, None, 7), slice(None, None, 7)); cc = float(np.corrcoef(F[sub].ravel(), R[sub].ravel())[0, 1]); print("corr(forward, reverse) on 1/49 subsample = %.4f; identical arrays: %s" % (cc, bool(np.array_equal(F[sub], R[sub]))))
both = (F >= thr) & (R >= thr)
lab, n = ndi.label(both)
floor_px = 500.0/um
rows = []; objs = ndi.find_objects(lab)
for i, sl in enumerate(objs, 1):
    if sl is None: continue
    m = lab[sl] == i; ys, xs = np.nonzero(m)
    h, w = sl[0].stop-sl[0].start, sl[1].stop-sl[1].start
    # long axis via PCA of pixel coords
    if len(ys) >= 3:
        c = np.cov(np.stack([xs, ys]).astype(float)); ev = np.linalg.eigvalsh(c); long_ax = 4*np.sqrt(max(ev[-1], 1e-9)); short_ax = 4*np.sqrt(max(ev[0], 1e-9))
    else: long_ax, short_ax = max(h, w), min(h, w)
    passes = long_ax >= floor_px
    rows.append(dict(id=i, y0=sl[0].start, x0=sl[1].start, h=h, w=w, area_px=int(m.sum()), long_axis_px=round(long_ax,1), long_axis_mm=round(long_ax*um/1000,3), aspect=round(long_ax/max(short_ax,1e-6),2), mean_p_fwd=float(F[sl][m].mean()), mean_p_rev=float(R[sl][m].mean()), passes_size_floor=bool(passes), at_array_edge=bool(sl[0].start==0 or sl[1].start==0 or sl[0].stop==F.shape[0] or sl[1].stop==F.shape[1])))
rows.sort(key=lambda r: -r["long_axis_px"])
import csv
with open(out+"_candidates.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()) if rows else ["id"]); w.writeheader(); w.writerows(rows)
summ = dict(threshold=thr, shape=list(F.shape), pixels_total=int(F.size), pixels_above_both=int(both.sum()), components=int(n), pass_size_floor=int(sum(r["passes_size_floor"] for r in rows)), size_floor_px=floor_px)
json.dump(summ, open(out+"_summary.json", "w"), indent=1); print(json.dumps(summ))
for r in [r for r in rows if r["passes_size_floor"]][:40]:
    pad = 150; y0, x0 = max(0, r["y0"]-pad), max(0, r["x0"]-pad); y1, x1 = min(F.shape[0], r["y0"]+r["h"]+pad), min(F.shape[1], r["x0"]+r["w"]+pad)
    fig, ax = plt.subplots(1, 2, figsize=(10, 5), dpi=100)
    ax[0].imshow(F[y0:y1, x0:x1], cmap="gray", vmin=0, vmax=1); ax[0].set_title("cand %d fwd  long %.3f mm  aspect %.2f" % (r["id"], r["long_axis_mm"], r["aspect"]))
    ax[1].imshow(both[y0:y1, x0:x1], cmap="gray"); ax[1].set_title("above threshold in BOTH directions")
    for a in ax: a.axis("off")
    plt.tight_layout(); plt.savefig("%s_cand%04d.png" % (out, r["id"])); plt.close(fig)
print("AUDIT DONE")
