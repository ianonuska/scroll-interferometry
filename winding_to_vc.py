#!/usr/bin/env python3
"""Export an automated winding-coordinate field as Vesuvius Challenge
point-collection JSONs (vc_pointcollections_json_version 1), ready for the
spiral-fitting pipeline (ScrollPrize/villa spiral-fitting).

Two outputs, matching the two constraint families the fitter consumes:

  relative_windings.json  — "pearl string" collections: points sampled along
      radial paths at successive integer winding values, each point carrying
      wind_a = 0,1,2,…  (relative-winding constraints)
  same_windings.json      — collections sampled along individual iso-winding
      contours (same-winding constraints; no wind_a)

Coordinates are emitted as [z, y, x] in FULL-RESOLUTION voxels (spiral-input
dataset convention); pass --scale for the pyramid level used (level 3 → 8).
"""
from __future__ import annotations

import argparse
import json

import numpy as np
from scipy import ndimage as ndi


ORDER = "xyz"  # shipped spiral-input files use p=[x,y,z]; set "zyx" for docs-order


def _collection(cid, name, pts_zyx, wind_a=None, color=(1.0, 0.4, 0.0)):
    points = {}
    for i, p in enumerate(pts_zyx):
        z, y, x = float(p[0]), float(p[1]), float(p[2])
        coords = [x, y, z] if ORDER == "xyz" else [z, y, x]
        d = {"p": coords, "creation_time": 0}
        if wind_a is not None:
            d["wind_a"] = float(wind_a[i])
        points[str(i)] = d
    return {
        "name": name,
        "color": list(color),
        "metadata": {"generator": "scroll-interferometry/winding_to_vc",
                     "method": "auto (interferometric winding field)"},
        "points": points,
    }


def export(W, quality, z_index, scale, out_prefix,
           n_rays=24, q_thresh=0.05, contour_step=2.0, pts_per_contour=24):
    mask = np.isfinite(W) & (quality > q_thresh)
    ys, xs = np.nonzero(mask)
    core = (float(xs.mean()), float(ys.mean()))
    ny, nx = W.shape

    # ---- relative-winding pearl strings along GRADIENT STREAMLINES of W.
    # Straight rays assume the scroll is star-convex around the core; squashed
    # scrolls are not, and a ray can re-enter windings, producing collections
    # whose wind_a oscillates — self-contradictory constraints. Streamlines of
    # +grad(W) are monotone in W by construction on any geometry.
    rel = {}
    cid = 0
    # mask-aware (normalized-convolution) smoothing: holes and the outside
    # region must not inject fake values into the gradient field
    maskf = mask.astype(np.float64)
    num = ndi.gaussian_filter(np.where(mask, W, 0.0), 3.0)
    den = ndi.gaussian_filter(maskf, 3.0)
    Ws = np.where(den > 0.2, num / np.maximum(den, 1e-9), np.nan)
    gy, gx = np.gradient(np.nan_to_num(Ws))
    valid_dir = den > 0.35
    gnorm = np.hypot(gx, gy) + 1e-9
    ux_f, uy_f = gx / gnorm, gy / gnorm
    # quality gate: |grad W| is the local winding density (windings/px); where
    # the implied wavelength drops below the resolvable limit the field merges
    # windings and its deltas under-count — assert nothing there.
    min_wavelength_px = 6.0
    resolved = gnorm < (1.0 / min_wavelength_px)
    # seeds: a well-spread grid of HIGH-QUALITY points across the whole field
    # (seeding at the winding minimum lands in the crushed core, where damage
    # kills walks immediately). From each seed, walk both +grad and -grad,
    # then splice into one monotone pearl string.
    step = 2.0
    hole_limit = 24          # coast up to ~3 local wavelengths of damage
    max_steps = int(3 * max(ny, nx) / step)
    n_seeds = max(n_rays * 4, 96)
    qgood = mask & (quality > np.nanpercentile(np.where(mask, quality, np.nan), 60))
    gy2, gx2 = np.nonzero(qgood)
    if gy2.size:
        picks = np.random.default_rng(0).choice(gy2.size, min(n_seeds, gy2.size), replace=False)

        def walk(x, y, sgn):
            out = []
            holes = 0
            last_dx, last_dy = 1.0, 0.0
            w_prev = None
            for _ in range(max_steps):
                xi, yi = int(round(x)), int(round(y))
                if not (0 <= xi < nx and 0 <= yi < ny):
                    break
                if not mask[yi, xi]:
                    holes += 1
                    if holes > hole_limit:
                        break
                    x += sgn * step * last_dx; y += sgn * step * last_dy
                    continue
                holes = 0
                w_here = W[yi, xi]
                if np.isfinite(w_here) and quality[yi, xi] > q_thresh:
                    if not resolved[yi, xi]:
                        # compression-suspect region: break the chain so no
                        # delta is asserted across it
                        w_prev = None
                    elif w_prev is not None and (w_here - w_prev) * sgn > 0:
                        lo, hi = sorted((w_prev, w_here))
                        for m in range(int(np.ceil(lo)), int(np.floor(hi)) + 1):
                            out.append((m, x, y))
                        w_prev = w_here
                    else:
                        w_prev = w_here
                dx_, dy_ = ux_f[yi, xi], uy_f[yi, xi]
                if not valid_dir[yi, xi] or not (np.isfinite(dx_) and np.isfinite(dy_)):
                    break
                last_dx, last_dy = dx_, dy_
                x += sgn * step * dx_
                y += sgn * step * dy_
            return out

        for p_i in picks:
            x0_, y0_ = float(gx2[p_i]), float(gy2[p_i])
            fwd = walk(x0_, y0_, +1)
            bwd = walk(x0_, y0_, -1)
            merged = {}
            for m, x, y in bwd + fwd:      # fwd overwrites duplicates at the seam
                merged[m] = (x, y)
            ms = sorted(merged)
            if len(ms) >= 4:
                base = ms[0]
                pearls = [(z_index * scale, merged[m][1] * scale, merged[m][0] * scale) for m in ms]
                winds = [m - base for m in ms]
                rel[str(cid)] = _collection(cid, f"auto_rel_stream{cid:02d}", pearls, winds)
                cid += 1

    # ---- same-winding contours
    same = {}
    sid = 0
    Wf = np.where(mask, W, np.nan)
    for level in np.arange(np.ceil(np.nanmin(Wf)), np.floor(np.nanmax(Wf)) + 0.1, contour_step):
        near = np.abs(Wf - level) < 0.05
        yy, xx = np.nonzero(near)
        if yy.size < pts_per_contour:
            continue
        sel = np.linspace(0, yy.size - 1, pts_per_contour).astype(int)
        pts = [(z_index * scale, float(yy[s]) * scale, float(xx[s]) * scale) for s in sel]
        same[str(sid)] = _collection(sid, f"auto_same_w{level:+.0f}", pts,
                                     color=(0.2, 0.5, 0.9))
        sid += 1

    for name, coll in (("relative_windings", rel), ("same_windings", same)):
        out = {"vc_pointcollections_json_version": "1", "collections": coll}
        path = f"{out_prefix}_{name}.json"
        with open(path, "w") as f:
            json.dump(out, f)
        npts = sum(len(c["points"]) for c in coll.values())
        print(f"wrote {path}: {len(coll)} collections, {npts} points")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("winding_npy")
    ap.add_argument("quality_npy")
    ap.add_argument("--z", type=int, required=True, help="slice index at the analysis level")
    ap.add_argument("--scale", type=int, default=8, help="voxel scale to full-res (level 3 = 8)")
    ap.add_argument("--out", default="auto")
    ap.add_argument("--order", choices=["xyz", "zyx"], default="xyz",
                    help="coordinate order in output p (shipped datasets: xyz)")
    a = ap.parse_args()
    ORDER = a.order
    W = np.load(a.winding_npy)
    Q = np.load(a.quality_npy)
    export(W, Q, a.z, a.scale, a.out)
