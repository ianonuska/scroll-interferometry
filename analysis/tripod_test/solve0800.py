"""Tripod-test stations on PHerc0800: for each inspected window, L2 solve ->
prior -> L1 locked solve (pitch lock + face gate), cored on the community
umbilicus. Preregistered in preregistrations/2026-09-02_verifier_tripod_test.md."""
import json, os, sys, time
import numpy as np
from scipy import ndimage as ndi
sys.path.insert(0, "/workspace/scroll-interferometry")
import zarr
from winding_phase import winding_coordinate
t0 = time.time(); log = lambda m: print(f"[{time.time()-t0:6.0f}s] {m}", flush=True)
OUT = "/workspace/p0800"; os.makedirs(OUT, exist_ok=True)
B = "https://vesuvius-challenge-open-data.s3.us-east-1.amazonaws.com"
g = zarr.open_group(f"{B}/PHerc0800/volumes/20250521135224-8.640um-1.2m-116keV-masked.zarr", mode="r")
L1, L2 = g["1"], g["2"]
log(f"L1 {L1.shape} L2 {L2.shape}")
umb = json.load(open("/workspace/herculaneum-umbilici/PHerc0800_umbilicus.json"))["control_points"]
uz = np.array([p["z"] for p in umb], float); ux = np.array([p["x"] for p in umb], float); uy = np.array([p["y"] for p in umb], float)
o = np.argsort(uz); uz, ux, uy = uz[o], ux[o], uy[o]
WINDOWS = [5664, 6272, 6864, 7472, 8064, 8672, 9264, 9872, 11072, 12272, 13472, 14064, 14672, 15264, 15872, 16464, 17072]
for wz in WINDOWS:
    zs = wz + 400
    if os.path.exists(f"{OUT}/L1lock_z{zs}_W.npy"):
        log(f"z{zs} exists, skip"); continue
    # level 2 solve (prior)
    img2 = np.asarray(L2[zs // 4]).astype(np.float64); m2 = ndi.binary_erosion(img2 > 0, iterations=3)
    core2 = (float(np.interp(zs, uz, ux)) / 4, float(np.interp(zs, uz, uy)) / 4)
    W2, Q2, *_ = winding_coordinate(img2, m2, core2)
    gy, gx = np.gradient(np.nan_to_num(W2)); gg = np.hypot(gx, gy)
    gg = np.where(np.isfinite(W2) & (Q2 > 0.05) & (gg > 1e-4), gg, np.nan); gg = np.where(np.isnan(gg), np.nanmedian(gg), gg)
    gg = ndi.gaussian_filter(gg, 4.0)
    # level 1 locked solve
    img1 = np.asarray(L1[zs // 2]).astype(np.float64); m1 = ndi.binary_erosion(img1 > 0, iterations=4)
    g1 = ndi.zoom(gg, (img1.shape[0] / gg.shape[0], img1.shape[1] / gg.shape[1]), order=1) / 2.0
    core1 = (float(np.interp(zs, uz, ux)) / 2, float(np.interp(zs, uz, uy)) / 2)
    W1, Q1, th, km, info = winding_coordinate(img1, m1, core1, k_prior=2 * np.pi * g1, lock_tol=0.6, ridge_min_frac=0.55)
    np.save(f"{OUT}/L1lock_z{zs}_W.npy", W1.astype(np.float32)); np.save(f"{OUT}/L1lock_z{zs}_Q.npy", Q1.astype(np.float32))
    np.save(f"{OUT}/mask_z{zs}.npy", m1)
    log(f"z{zs}: L2 span {np.nanmax(W2)-np.nanmin(W2):.1f} -> L1 locked cg={info} span {np.nanmax(W1)-np.nanmin(W1):.1f}")
log("P0800 STATIONS DONE")
