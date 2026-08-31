#!/usr/bin/env python3
"""Confound analysis for the listening pilot — the table the prereg
promised alongside the result. POST-HOC and clearly labeled as such: the
preregistered verdict comes from listening_pilot.py alone; this script asks
whether that verdict survives the obvious confounds.

Checks per fragment:
  C1 brightness: ink patches may sit on systematically brighter/darker
     papyrus. Report the ink-vs-control mean-intensity Cliff's delta, then
     REDO each statistic's Mann-Whitney on intensity-matched pairs
     (nearest-neighbor matching on patch mean intensity, one control per
     ink patch). A texture statistic that only proxies brightness dies here.
  C2 edge-distance match quality: median |edge_dist(ink) - edge_dist(ctl)|
     after the main analysis' matching.
  C3 spatial clustering: ink letters cluster; report the median
     nearest-other-ink-patch distance vs the same for controls, so a reader
     can judge how spatially interleaved the comparison really is.

Usage: python analysis/listening_confounds.py FRAGDIR...
"""
from __future__ import annotations

import os
import sys

import numpy as np
from scipy import ndimage as ndi
from scipy import stats

sys.path.insert(0, os.path.dirname(__file__))
from listening_pilot import (LAYERS, PATCH, ALPHA, N_TESTS, otsu,
                             patch_stats, cliffs_delta)


def analyze(fragdir):
    import tifffile
    from PIL import Image
    Image.MAX_IMAGE_PIXELS = None
    name = os.path.basename(fragdir.rstrip("/"))
    stack = np.stack([
        tifffile.imread(os.path.join(fragdir, f)).astype(np.float64)
        for f in LAYERS])
    ir_path = os.path.join(fragdir, "ir.png")
    if os.path.exists(ir_path):
        ir = np.asarray(Image.open(ir_path).convert("L"), dtype=np.float64)
        fmask = np.asarray(Image.open(os.path.join(fragdir, "mask.png"))
                           .convert("L")) > 127
        ink = (ir < otsu(ir[fmask])) & fmask
    else:  # Frag4, amendment 1
        ink = np.asarray(Image.open(os.path.join(fragdir, "inklabels.png"))
                         .convert("L")) > 127
        fmask = np.asarray(Image.open(os.path.join(fragdir, "mask.png"))
                           .convert("L")) > 127
        ink &= fmask
    edge = ndi.distance_transform_edt(fmask)
    ny, nx = fmask.shape
    rows = []          # (kind, cy, cx, edge_med, intensity, S1..S4)
    for y0 in range(0, ny - PATCH + 1, PATCH):
        for x0 in range(0, nx - PATCH + 1, PATCH):
            sl = np.s_[y0:y0 + PATCH, x0:x0 + PATCH]
            if fmask[sl].mean() < 0.95:
                continue
            frac = ink[sl].mean()
            kind = "ink" if frac >= 0.70 else ("ctl" if frac == 0.0 else None)
            if kind is None:
                continue
            sub = stack[(np.s_[:],) + sl]
            rows.append((kind, y0 + 32, x0 + 32, np.median(edge[sl]),
                         float(sub.mean())) + patch_stats(sub))
    kinds = np.array([r[0] for r in rows])
    P = np.array([r[1:] for r in rows], dtype=float)
    ik, ck = kinds == "ink", kinds == "ctl"
    if ik.sum() < 10 or ck.sum() < 10:
        print(f"{name}: underpowered, skipped")
        return
    print(f"\n== {name}  (ink {ik.sum()}, control {ck.sum()})")
    inten_d = cliffs_delta(P[ik, 3], P[ck, 3])
    print(f"  C1 brightness delta (ink vs all controls): {inten_d:+.3f}")

    # intensity-matched controls, one per ink patch, no reuse
    ctl_idx = np.nonzero(ck)[0]
    avail = np.ones(len(ctl_idx), bool)
    pairs = []
    for i in np.nonzero(ik)[0]:
        cand = np.nonzero(avail)[0]
        if not len(cand):
            break
        j = cand[np.argmin(np.abs(P[ctl_idx[cand], 3] - P[i, 3]))]
        avail[j] = False
        pairs.append((i, ctl_idx[j]))
    ii = np.array([p[0] for p in pairs])
    jj = np.array([p[1] for p in pairs])
    resid = np.median(np.abs(P[ii, 3] - P[jj, 3]))
    spread = np.median(np.abs(P[ck, 3] - np.median(P[ck, 3])))
    print(f"  C1 match quality: median |dI| {resid:.1f} "
          f"(control spread {spread:.1f})")
    for k, sname in enumerate(["S1_var", "S2_kurt", "S3_ac1", "S4_hf"]):
        col = 4 + k
        p = stats.mannwhitneyu(P[ii, col], P[jj, col]).pvalue
        d = cliffs_delta(P[ii, col], P[jj, col])
        alive = p < ALPHA / N_TESTS
        print(f"  C1 {sname} after intensity matching: p={p:.2e} "
              f"delta {d:+.3f} -> {'SURVIVES' if alive else 'dies'}")
    print(f"  C2 edge-dist match: ink med {np.median(P[ik, 2]):.0f} vs "
          f"matched-ctl med {np.median(P[jj, 2]):.0f} px")
    from scipy.spatial import cKDTree
    for lab, m in [("ink", ik), ("ctl", ck)]:
        pts = P[m][:, :2]
        if len(pts) > 1:
            t = cKDTree(pts)
            dd, _ = t.query(pts, k=2)
            print(f"  C3 {lab} median NN distance: {np.median(dd[:, 1]):.0f} px")


if __name__ == "__main__":
    for d in sys.argv[1:]:
        analyze(d)
