import zarr, numpy as np, json
from scipy.spatial import cKDTree
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
z0,z1,y0,y1,x0,x1 = 18408,18664,2000,3024,6352,7376
sig = lambda x: 1/(1+np.exp(-x.astype(np.float32)))
A = sig(np.asarray(zarr.open_array("/workspace/ink3d/runA_blended.zarr", mode="r")[0]))          # (256,1024,1024)
Bp = sig(np.asarray(zarr.open_array("/workspace/ink3d/runB_blended.zarr", mode="r")[0]))[96:160]  # (64,256,256), padding removed
mA = A >= 0.5; mB = Bp >= 0.5
res = {"runA_frac_p_ge_0.5": float(mA.mean()), "runA_p99": float(np.percentile(A, 99)), "runB_frac_p_ge_0.5": float(mB.mean()), "runB_p99": float(np.percentile(Bp, 99))}
Bup = np.repeat(np.repeat(np.repeat(mB, 4, 0), 4, 1), 4, 2)
Adn = A.reshape(64,4,256,4,256,4).mean((1,3,5)) >= 0.5
dice = lambda a, b: float(2*(a&b).sum()/max(1, a.sum()+b.sum()))
res["dice_B_up_vs_A"] = dice(Bup, mA); res["dice_B_vs_A_down"] = dice(mB, Adn)
print(json.dumps(res, indent=1), flush=True)
# figure: mid slice
blk = np.load("/workspace/ink3d/block_2p4um.npy"); zc = 128
pts = np.load("/workspace/ink3d/block_mesh_points.npz"); px, py, pz, pw = pts["x"]-x0, pts["y"]-y0, pts["z"]-z0, pts["w"]
fig, ax = plt.subplots(1, 3, figsize=(19, 6.6), dpi=100)
ax[0].imshow(blk[zc], cmap="gray"); ax[0].set_title("PHerc 1667, 2.4 µm, block mid-slice (z 18536)")
ax[1].imshow(blk[zc], cmap="gray"); ov = np.ma.masked_where(A[zc] < 0.5, A[zc]); ax[1].imshow(ov, cmap="autumn", vmin=0.5, vmax=1, alpha=0.9); ax[1].set_title("run A: 3D ink p ≥ 0.5 (red/yellow) at 2.4 µm")
sel = np.abs(pz - zc) < 6; ax[1].scatter(px[sel], py[sel], c=pw[sel], cmap="tab20", s=2); 
ax[2].imshow(blk[zc], cmap="gray"); ovb = np.ma.masked_where(Bup[zc] == 0, Bup[zc].astype(float)); ax[2].imshow(ovb, cmap="cool", alpha=0.6); ax[2].set_title("run B: same block at 9.6 µm, p ≥ 0.5 (cyan), upsampled")
for a in ax: a.axis("off")
plt.tight_layout(); plt.savefig("/workspace/ink3d/fig_ink3d_block.png"); print("fig ok", flush=True)
# run C: page assignment
W = np.load("/workspace/gauge1667/W_L2_locked2.npy").astype(float)
iz, iy, ix = np.nonzero(mA); rng = np.random.default_rng(0)
if len(iz) > 400000: k = rng.choice(len(iz), 400000, replace=False); iz, iy, ix = iz[k], iy[k], ix[k]
tree = cKDTree(np.stack([pz, py, px], 1)); d, j = tree.query(np.stack([iz, iy, ix], 1), distance_upper_bound=12)
ok = np.isfinite(d); iz, iy, ix, true_w = iz[ok], iy[ok], ix[ok], pw[j[ok]]
Wv = W[np.clip((iy+y0)//4, 0, W.shape[0]-1), np.clip((ix+x0)//4, 0, W.shape[1]-1)]
ref = pw == 28; Wref = np.nanmedian(W[np.clip(py[ref].astype(int)//4 + y0//4, 0, W.shape[0]-1), np.clip(px[ref].astype(int)//4 + x0//4, 0, W.shape[1]-1)])
assigned = np.round(Wv - Wref); true_off = true_w - 28; good = np.isfinite(assigned)
acc = float((assigned[good] == true_off[good]).mean()); res["runC_ink_voxels_labelled"] = int(good.sum()); res["runC_accuracy"] = acc
per = {}
for w in sorted(set(true_w)):
    m = good & (true_w == w); 
    if m.sum() < 50: continue
    per[int(w)] = {"n": int(m.sum()), "acc": float((assigned[m] == true_off[m]).mean()), "median_assigned_minus_true": float(np.median(assigned[m] - true_off[m]))}
res["runC_per_wrap"] = per
print(f"RUN C: {good.sum()} labelled ink voxels, page-assignment accuracy {acc:.3f} (pass >= 0.80, kill < 0.50)")
for w, r in per.items(): print(f"  w{w:02d} n={r['n']:6d} acc {r['acc']:.2f} median(assigned-true) {r['median_assigned_minus_true']:+.1f}")
json.dump(res, open("/workspace/ink3d/ink3d_results.json", "w"), indent=1); print("INK3D ANALYSIS DONE")
