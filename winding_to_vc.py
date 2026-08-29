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


def _collection(cid, name, pts_zyx, wind_a=None, color=(1.0, 0.4, 0.0)):
    points = {}
    for i, p in enumerate(pts_zyx):
        d = {"p": [float(p[0]), float(p[1]), float(p[2])], "creation_time": 0}
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

    # ---- relative-winding pearl strings along rays
    rel = {}
    cid = 0
    for ang in np.linspace(0, 2 * np.pi, n_rays, endpoint=False):
        ts = np.arange(0, max(ny, nx), 0.5)
        rx = core[0] + ts * np.cos(ang)
        ry = core[1] + ts * np.sin(ang)
        ok = (rx >= 0) & (rx < nx - 1) & (ry >= 0) & (ry < ny - 1)
        rx, ry = rx[ok], ry[ok]
        inb = mask[ry.astype(int), rx.astype(int)]
        if inb.sum() < 20:
            continue
        Wray = ndi.map_coordinates(np.nan_to_num(W), [ry, rx], order=1)
        Qray = ndi.map_coordinates(quality, [ry, rx], order=1)
        Wray = np.where(inb & (Qray > q_thresh), Wray, np.nan)
        # walk outward, drop a pearl at each integer winding crossing
        pearls, winds = [], []
        base = None
        for i in range(1, len(ts[ok])):
            a, b = Wray[i - 1], Wray[i]
            if not (np.isfinite(a) and np.isfinite(b)) or b == a:
                continue
            lo, hi = sorted((a, b))
            for m in range(int(np.ceil(lo)), int(np.floor(hi)) + 1):
                f = (m - a) / (b - a)
                px, py = rx[i - 1] + f * (rx[i] - rx[i - 1]), ry[i - 1] + f * (ry[i] - ry[i - 1])
                if base is None:
                    base = m
                pearls.append((z_index * scale, py * scale, px * scale))
                winds.append(m - base)
        if len(pearls) >= 4:
            rel[str(cid)] = _collection(cid, f"auto_rel_ray{cid:02d}", pearls, winds)
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
    a = ap.parse_args()
    W = np.load(a.winding_npy)
    Q = np.load(a.quality_npy)
    export(W, Q, a.z, a.scale, a.out)
