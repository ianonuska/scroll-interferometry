#!/usr/bin/env python3
"""Encode winding-coordinate fields as a spiral-fitting winding_inference
store ("winding_inference_crossings", format_version 1).

This is the fitter's dense power port: the input that normally comes from a
trained winding model and only exists for scrolls that have one. This encoder
lets the interferometric field serve that port instead, so fit_spiral's
dense_spacing_mode="winding_model" runs on scrolls with no trained model and
no annotations.

Store semantics (mirrors spiral-fitting/winding_supervision.py, the loader):
  - straight rays, point(t) = origin_zyx + t * step_zyx, coordinates in
    reference-zarr scale-0 voxels, zyx order;
  - per ray, ascending crossing_t with integer crossing_level = the winding
    index at each crossing; the fitter samples crossing pairs along a ray and
    supervises with the level difference.

Our pearl streams (collect_pearl_streams) are exactly integer-winding
crossings along gradient streamlines, which curve. Rays must be straight, so
each stream is chord-split greedily: extend a chord while every pearl on it
stays within `chord_tol` full-res voxels of the fitted line and the
projections stay strictly monotone. Deltas are therefore only ever asserted
along locally straight, monotone runs — never across a bend that could alias.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os

import numpy as np

from winding_to_vc import collect_pearl_streams

ARTIFACT_TYPE = "winding_inference_crossings"
FORMAT_VERSION = 1


def _canonical_digest(value) -> str:
    encoded = json.dumps(value, sort_keys=True,
                         separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _sha256_file(path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(8 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def chord_split(points, levels, chord_tol):
    """Split one polyline of crossings into straight monotone rays.

    points: (n, 3) float array (zyx, full-res voxels), levels: (n,) int.
    Returns a list of (origin, step, t, level) with len(t) >= 2, t[0] = 0,
    strictly increasing t, |step| = 1.
    """
    rays = []
    n = len(points)
    i = 0
    while i < n - 1:
        # unit-step guard: consecutive crossings must differ by exactly one
        # winding. Anything else means the stream skipped fringes (a hole,
        # a resolution break) and no delta may be asserted across it.
        if abs(int(levels[i + 1]) - int(levels[i])) != 1:
            i += 1
            continue
        j = i + 1
        best = j
        while j + 1 <= n - 1:
            if abs(int(levels[j + 1]) - int(levels[j])) != 1:
                break
            j += 1
            a, b = points[i], points[j]
            d = b - a
            length = np.linalg.norm(d)
            if length < 1e-6:
                break
            u = d / length
            rel = points[i:j + 1] - a
            t = rel @ u
            if not np.all(np.diff(t) > 1e-6):
                break
            perp = rel - t[:, None] * u
            if np.max(np.linalg.norm(perp, axis=1)) > chord_tol:
                break
            best = j
        a, b = points[i], points[best]
        d = b - a
        length = np.linalg.norm(d)
        if length > 1e-6:
            u = d / length
            rel = points[i:best + 1] - a
            t = rel @ u
            if np.all(np.diff(t) > 1e-6):
                rays.append((a.astype(np.float32), u.astype(np.float32),
                             t.astype(np.float32),
                             np.asarray(levels[i:best + 1], dtype=np.int16)))
        i = best
    return rays


def encode_station(W, quality, z_index, scale, q_thresh=0.05,
                   chord_tol=2.0, n_seeds=96):
    """One z-station -> list of rays in full-res zyx voxels."""
    rays = []
    for g in collect_pearl_streams(W, quality, q_thresh=q_thresh,
                                   n_seeds=n_seeds):
        ms = sorted(g)
        pts = np.array([[float(z_index) * scale,
                         g[m][1] * scale,
                         g[m][0] * scale] for m in ms], dtype=np.float64)
        out = chord_split(pts, ms, chord_tol * scale)
        rays.extend(out)
    return rays


def write_store(stations, out_dir, provenance=None):
    """stations: list of (name, rays) -> one shard per station."""
    os.makedirs(out_dir, exist_ok=True)
    shard_entries = []
    num_rays = 0
    num_crossings = 0
    for name, rays in stations:
        shard_dir = os.path.join(out_dir, name)
        os.makedirs(shard_dir, exist_ok=True)
        origins = np.array([r[0] for r in rays], dtype=np.float32) \
            if rays else np.empty((0, 3), np.float32)
        steps = np.array([r[1] for r in rays], dtype=np.float32) \
            if rays else np.empty((0, 3), np.float32)
        t_flat = np.concatenate([r[2] for r in rays]).astype(np.float32) \
            if rays else np.empty((0,), np.float32)
        lvl_flat = np.concatenate([r[3] for r in rays]).astype(np.int16) \
            if rays else np.empty((0,), np.int16)
        offsets = np.zeros(len(rays) + 1, dtype=np.int64)
        np.cumsum([len(r[2]) for r in rays], out=offsets[1:])
        seed = np.array([r[3][0] for r in rays], dtype=np.int16) \
            if rays else np.empty((0,), np.int16)
        arrays = {}
        for arr_name, arr in [("ray_origin_zyx", origins),
                              ("ray_step_zyx", steps),
                              ("crossing_t", t_flat),
                              ("crossing_level", lvl_flat),
                              ("seed_winding", seed),
                              ("crossing_offsets", offsets)]:
            path = os.path.join(shard_dir, arr_name + ".npy")
            np.save(path, arr)
            arrays[arr_name] = {
                "file": arr_name + ".npy",
                "dtype": np.dtype(arr.dtype).str,
                "shape": list(arr.shape),
                "bytes": os.path.getsize(path),
                "sha256": _sha256_file(path),
            }
        shard_entries.append({
            "name": name,
            "arrays": arrays,
            "num_crossings": int(len(t_flat)),
            "num_retained_rays": int(len(rays)),
        })
        num_rays += len(rays)
        num_crossings += int(len(t_flat))
    manifest = {
        "artifact_type": ARTIFACT_TYPE,
        "format_version": FORMAT_VERSION,
        "coordinate_order": "zyx",
        "coordinate_space": "reference zarr scale-0 voxels",
        "generator": "scroll-interferometry/winding_to_store",
        "num_rays": num_rays,
        "num_crossings": num_crossings,
        "shards": shard_entries,
    }
    if provenance:
        manifest["source_attributes"] = provenance
    manifest["fingerprint"] = _canonical_digest(manifest)
    with open(os.path.join(out_dir, "manifest.json"), "w") as f:
        json.dump(manifest, f, indent=2, sort_keys=True)
    print(f"wrote store {out_dir}: {num_rays} rays, "
          f"{num_crossings} crossings, {len(stations)} shard(s)")
    return manifest


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("pairs", nargs="+",
                    help="stationZ:W.npy:Q.npy triplets, e.g. "
                         "11056:st11056_W.npy:st11056_Q.npy")
    ap.add_argument("--scale", type=int, default=1,
                    help="analysis-level to full-res factor (station solves: 1)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--chord-tol", type=float, default=2.0)
    a = ap.parse_args()
    stations = []
    for spec in a.pairs:
        z_str, wp, qp = spec.split(":")
        W = np.load(wp)
        Q = np.load(qp)
        rays = encode_station(W, Q, int(z_str), a.scale,
                              chord_tol=a.chord_tol)
        n_cross = sum(len(r[2]) for r in rays)
        print(f"  station z={z_str}: {len(rays)} rays, {n_cross} crossings")
        stations.append((f"station_z{z_str}", rays))
    write_store(stations, a.out,
                provenance={"fields": a.pairs, "scale": a.scale,
                            "chord_tol": a.chord_tol})
