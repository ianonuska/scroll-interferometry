#!/usr/bin/env python3
"""Listening pilot — the DIRECTION confound test (post-hoc, labeled).

Hole in the published result: `patch_stats` averages the lag-1
autocorrelation along x and y, discarding exactly the information that
separates real microstructure (isotropic, or aligned with papyrus fibers)
from reconstruction artifacts (aligned with the sampling grid / scanner
geometry). This script recomputes the ink-vs-control effect PER DIRECTION.

Criterion, declared before running (2026-09-01):
  The ink-vs-control shift must NOT be confined to a single grid axis. Report
  Cliff's delta for the x-lag and y-lag autocorrelations separately and for
  autocorrelation at 45/135 degrees. If |delta| on one grid axis exceeds three
  times the other AND the diagonal deltas are near zero, the signal is
  grid-locked and the fingerprint claim must be weakened in the results file.
  Caveat to carry: surface volumes are resampled along the surface normal,
  which partially scrambles scanner geometry — a grid-locked signal is
  suspicious, an isotropic one is not proof of physics.

Usage: python analysis/listening_anisotropy.py FRAGDIR... [--mask-source inklabels]
"""
from __future__ import annotations
import argparse, os, sys
import numpy as np
from scipy import ndimage as ndi
from scipy import stats

sys.path.insert(0, os.path.dirname(__file__))
import listening_pilot as lp
from listening_pilot import LAYERS, PATCH, ALPHA, N_TESTS, otsu, cliffs_delta


def dir_stats(stack):
    """Per-layer lag-1 autocorrelation along x, y, and the two diagonals,
    averaged over the 5 layers. Returns (ax, ay, ad1, ad2)."""
    out = np.zeros(4)
    for L in stack:
        L = L - L.mean()
        v = L.var()
        if v <= 0:
            continue
        ax = np.mean(L[:, :-1] * L[:, 1:]) / v
        ay = np.mean(L[:-1, :] * L[1:, :]) / v
        ad1 = np.mean(L[:-1, :-1] * L[1:, 1:]) / v
        ad2 = np.mean(L[:-1, 1:] * L[1:, :-1]) / v
        out += (ax, ay, ad1, ad2)
    return out / stack.shape[0]


def analyze(fragdir, mask_source):
    import tifffile
    from PIL import Image
    Image.MAX_IMAGE_PIXELS = None
    name = os.path.basename(fragdir.rstrip("/"))
    stack = np.stack([tifffile.imread(os.path.join(fragdir, f)).astype(np.float64)
                      for f in LAYERS])
    fmask = np.asarray(Image.open(os.path.join(fragdir, "mask.png")).convert("L")) > 127
    ir_path = os.path.join(fragdir, "ir.png")
    if mask_source == "inklabels" or not os.path.exists(ir_path):
        ink = (np.asarray(Image.open(os.path.join(fragdir, "inklabels.png")).convert("L")) > 127) & fmask
    else:
        ir = np.asarray(Image.open(ir_path).convert("L"), dtype=np.float64)
        ink = (ir < otsu(ir[fmask])) & fmask
    edge = ndi.distance_transform_edt(fmask)
    ny, nx = fmask.shape
    ink_f, ctl_f, ink_e, ctl_e = [], [], [], []
    for y0 in range(0, ny - PATCH + 1, PATCH):
        for x0 in range(0, nx - PATCH + 1, PATCH):
            sl = np.s_[y0:y0 + PATCH, x0:x0 + PATCH]
            if fmask[sl].mean() < 0.95:
                continue
            fr = ink[sl].mean()
            if fr >= 0.70:
                ink_f.append(dir_stats(stack[(np.s_[:],) + sl])); ink_e.append(np.median(edge[sl]))
            elif fr == 0.0:
                ctl_f.append(dir_stats(stack[(np.s_[:],) + sl])); ctl_e.append(np.median(edge[sl]))
    ink_f, ctl_f = np.array(ink_f), np.array(ctl_f)
    ink_e, ctl_e = np.array(ink_e), np.array(ctl_e)
    if len(ink_f) < 10 or len(ctl_f) < 10:
        print(f"{name}: underpowered ({len(ink_f)} ink / {len(ctl_f)} control)")
        return None
    # same edge-distance matching as the main analysis
    avail = np.ones(len(ctl_e), bool); chosen = []
    for i in np.argsort(ink_e):
        cand = np.nonzero(avail)[0]
        if not len(cand):
            break
        j = cand[np.argmin(np.abs(ctl_e[cand] - ink_e[i]))]
        avail[j] = False; chosen.append(j)
    ctl_m = ctl_f[chosen]
    names = ["x-lag", "y-lag", "diag+", "diag-"]
    deltas, ps = [], []
    for k, nm in enumerate(names):
        d = cliffs_delta(ink_f[:, k], ctl_m[:, k])
        p = stats.mannwhitneyu(ink_f[:, k], ctl_m[:, k]).pvalue
        deltas.append(d); ps.append(p)
    dx, dy, d1, d2 = deltas
    axis_ratio = max(abs(dx), abs(dy)) / max(min(abs(dx), abs(dy)), 1e-6)
    diag_mean = 0.5 * (abs(d1) + abs(d2))
    grid_locked = axis_ratio > 3.0 and diag_mean < 0.05
    print(f"{name}: n_ink={len(ink_f)} n_ctl={len(ctl_m)} | "
          + " ".join(f"{nm} d={d:+.3f}(p={p:.1e})" for nm, d, p in zip(names, deltas, ps))
          + f" | axis-ratio {axis_ratio:.2f} diag-mean {diag_mean:.3f} -> "
          + ("GRID-LOCKED (suspicious)" if grid_locked else "not grid-locked"))
    return {"fragment": name, "deltas": dict(zip(names, deltas)),
            "axis_ratio": axis_ratio, "diag_mean": diag_mean, "grid_locked": grid_locked}


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("fragdirs", nargs="+")
    ap.add_argument("--mask-source", choices=["ir", "inklabels"], default="inklabels")
    a = ap.parse_args()
    res = [analyze(d, a.mask_source) for d in a.fragdirs]
    locked = [r["fragment"] for r in res if r and r["grid_locked"]]
    print("VERDICT:", ("grid-locked on " + ", ".join(locked)) if locked else
          "no fragment shows a grid-locked signal; direction confound not supported")
