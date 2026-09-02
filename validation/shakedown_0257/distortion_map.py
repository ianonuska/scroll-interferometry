"""Item 9.5: flattening distortion map for a flattened tifxyz (lasagna/SLIM output).
Each flat cell is nominally a square of side s_flat = 1/scale px in the render; its 3D
footprint has edges e_u, e_v (voxels). Area ratio = |e_u x e_v| / s_flat^2 (1 = no
stretch); shear = angle between e_u and e_v (90 deg = no shear). Reported as maps and
percentiles; every published render from this surface should carry these numbers."""
import sys, json, numpy as np, tifffile
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
d, out = sys.argv[1], sys.argv[2]
X = np.stack([tifffile.imread(f"{d}/{a}.tif").astype(float) for a in ("x","y","z")], -1)
meta = json.load(open(f"{d}/meta.json")); sc = float(meta.get("scale", [0.05,0.05])[0]); s_flat = 1.0/sc
valid = np.isfinite(X[...,0]) & (X[...,0] > 0)
eu = X[:, 1:] - X[:, :-1]; ev = X[1:, :] - X[:-1, :]
eu = eu[:-1]; ev = ev[:, :-1]; v = valid[:-1, :-1] & valid[:-1, 1:] & valid[1:, :-1]
lu, lv = np.linalg.norm(eu, axis=-1), np.linalg.norm(ev, axis=-1)
area = np.linalg.norm(np.cross(eu, ev), axis=-1) / (s_flat**2)
cosang = (eu*ev).sum(-1) / (lu*lv + 1e-9); ang = np.degrees(np.arccos(np.clip(cosang, -1, 1)))
su, sv = lu/s_flat, lv/s_flat
stats = {}
for name, arr in (("area_ratio", area), ("stretch_u", su), ("stretch_v", sv), ("angle_deg", ang)):
    a = arr[v]; stats[name] = {"p5": float(np.percentile(a,5)), "p25": float(np.percentile(a,25)), "median": float(np.median(a)), "p75": float(np.percentile(a,75)), "p95": float(np.percentile(a,95))}
stats["cells_valid"] = int(v.sum()); stats["frac_area_within_10pct"] = float((np.abs(area[v]-1) < 0.1).mean()); stats["frac_area_within_25pct"] = float((np.abs(area[v]-1) < 0.25).mean()); stats["frac_shear_within_10deg"] = float((np.abs(ang[v]-90) < 10).mean())
json.dump(stats, open(out+"_distortion.json","w"), indent=1); print(json.dumps(stats, indent=1))
fig, ax = plt.subplots(2, 1, figsize=(16, 5.5), dpi=110)
im0 = ax[0].imshow(np.where(v, np.log2(np.clip(area, 0.25, 4)), np.nan), cmap="RdBu_r", vmin=-2, vmax=2, aspect="auto", interpolation="nearest"); ax[0].set_title("area distortion, log2(3D area / flat area): blue = compressed in the render, red = stretched, white = faithful"); plt.colorbar(im0, ax=ax[0], fraction=0.02)
im1 = ax[1].imshow(np.where(v, ang-90, np.nan), cmap="PuOr", vmin=-45, vmax=45, aspect="auto", interpolation="nearest"); ax[1].set_title("shear, degrees from a right angle (0 = no shear)"); plt.colorbar(im1, ax=ax[1], fraction=0.02)
for a in ax: a.set_yticks([]); a.set_xlabel("flat cells along the unrolled scroll (w010 at left → w065 at right)")
plt.tight_layout(); plt.savefig(out+"_distortion.png"); print("wrote", out+"_distortion.png")
