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
    gy, gx = np.gradient(ndi.gaussian_filter(np.nan_to_num(W), 3.0))
    gnorm = np.hypot(gx, gy) + 1e-9
    ux_f, uy_f = gx / gnorm, gy / gnorm
    # seeds: points near the inner end of the field, spread by angle around core
    Wmin = np.nanpercentile(np.where(mask, W, np.nan), 3)
    seed_band = mask & np.isfinite(W) & (W < Wmin + 1.5)
    sy, sx = np.nonzero(seed_band)
    if sy.size:
        angs = np.arctan2(sy - core[1], sx - core[0])
        order_idx = np.argsort(angs)
        picks = order_idx[np.linspace(0, order_idx.size - 1, n_rays).astype(int)]
        step = 2.0
        max_steps = int(3 * max(ny, nx) / step)
        for p_i in picks:
            x, y = float(sx[p_i]), float(sy[p_i])
            pearls, winds = [], []
            base = None
            w_prev = None
            for _ in range(max_steps):
                xi, yi = int(round(x)), int(round(y))
                if not (0 <= xi < nx and 0 <= yi < ny) or not mask[yi, xi]:
                    break
                w_here = W[yi, xi]
                q_here = quality[yi, xi]
                if np.isfinite(w_here) and q_here > q_thresh and w_prev is not None and w_here > w_prev:
                    for m in range(int(np.ceil(w_prev)), int(np.floor(w_here)) + 1):
                        if base is None:
                            base = m
                        pearls.append((z_index * scale, y * scale, x * scale))
                        winds.append(m - base)
                if np.isfinite(w_here):
                    w_prev = w_here
                dx_, dy_ = ux_f[yi, xi], uy_f[yi, xi]
                if not (np.isfinite(dx_) and np.isfinite(dy_)):
                    break
                x += step * dx_
                y += step * dy_
            # enforce strict monotonicity (streamline should guarantee it, but
            # guard against blur-field loops): keep the increasing prefix run
            if len(pearls) >= 4:
                keep_p, keep_w = [pearls[0]], [winds[0]]
                for pt, wv in zip(pearls[1:], winds[1:]):
                    if wv > keep_w[-1]:
                        keep_p.append(pt); keep_w.append(wv)
                if len(keep_p) >= 4:
                    rel[str(cid)] = _collection(cid, f"auto_rel_stream{cid:02d}", keep_p, keep_w)
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
