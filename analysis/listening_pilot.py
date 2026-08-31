#!/usr/bin/env python3
"""The "listening" pilot: does carbon ink leave a microstructure fingerprint
in standard CT texture? Preregistered analysis — see
preregistrations/2026-08-31_listening_pilot.md (committed before any data
was analyzed). This script implements exactly the declared statistics; the
few implementation choices the prereg left open are documented in ANALYSIS
CHOICES below and were fixed before running on any real fragment.

Usage:
    python analysis/listening_pilot.py FRAGDIR...   # one dir per fragment
    python analysis/listening_pilot.py --selftest   # synthetic blind checks

Each FRAGDIR must contain: 30.tif 31.tif 32.tif 33.tif 34.tif (surface
volume middle layer 32 +-2), ir.png (infrared photo), mask.png (fragment
extent). Output: per-fragment CSV + a verdict against the preregistered
criteria.

ANALYSIS CHOICES (fixed 2026-08-31, before any real-fragment run):
  - IR ink mask: Otsu threshold on ir.png inside mask.png; ink = darker side.
  - Patches: non-overlapping 64x64 grid; ink patch = >=70% ink pixels;
    control patch = 0% ink pixels; both require >=95% inside mask.png.
  - Edge matching: controls are subsampled so their distance-to-mask-edge
    distribution matches the ink patches' (nearest-neighbor matching on
    median patch edge-distance, one control per ink patch, no reuse).
  - Statistics per patch: computed per layer on the 5 layers, averaged.
      S1 variance; S2 excess kurtosis (Fisher); S3 mean of lag-1
      autocorrelation along x and along y; S4 fraction of FFT power above
      1/4 Nyquist (radial).
  - Test: two-sided Mann-Whitney U per fragment per statistic;
    significance p < 0.01 / 24 (Bonferroni: 4 stats x 6 fragments).
  - Effect size: Cliff's delta, reported regardless of significance.
"""
from __future__ import annotations

import argparse
import csv
import os
import sys

import numpy as np
from scipy import ndimage as ndi
from scipy import stats

LAYERS = ["30.tif", "31.tif", "32.tif", "33.tif", "34.tif"]
PATCH = 64
MASK_SOURCE = "ir"  # "ir" = preregistered Otsu-on-IR; "inklabels" = post-hoc sensitivity
ALPHA = 0.01
N_TESTS = 24  # 4 statistics x 6 fragments, declared in the prereg


def otsu(values):
    hist, edges = np.histogram(values, bins=256)
    mids = 0.5 * (edges[:-1] + edges[1:])
    w0 = np.cumsum(hist).astype(float)
    w1 = w0[-1] - w0
    m0 = np.cumsum(hist * mids)
    mu0 = np.where(w0 > 0, m0 / np.maximum(w0, 1), 0)
    mu1 = np.where(w1 > 0, (m0[-1] - m0) / np.maximum(w1, 1), 0)
    var_b = w0 * w1 * (mu0 - mu1) ** 2
    return mids[int(np.argmax(var_b))]


def patch_stats(stack):
    """stack: (5, 64, 64) float. Returns S1..S4 averaged over layers."""
    s1 = s2 = s3 = s4 = 0.0
    n = stack.shape[0]
    for L in stack:
        L = L - L.mean()
        v = L.var()
        s1 += v
        s2 += stats.kurtosis(L, axis=None, fisher=True, bias=True)
        if v > 0:
            ax = np.mean(L[:, :-1] * L[:, 1:]) / v
            ay = np.mean(L[:-1, :] * L[1:, :]) / v
        else:
            ax = ay = 0.0
        s3 += 0.5 * (ax + ay)
        F = np.abs(np.fft.rfft2(L)) ** 2
        fy = np.fft.fftfreq(L.shape[0])[:, None]
        fx = np.fft.rfftfreq(L.shape[1])[None, :]
        rad = np.sqrt(fy ** 2 + fx ** 2)
        tot = F.sum()
        s4 += float(F[rad > 0.125].sum() / tot) if tot > 0 else 0.0
    return s1 / n, s2 / n, s3 / n, s4 / n


def cliffs_delta(a, b):
    a, b = np.asarray(a), np.asarray(b)
    gt = sum((x > b).sum() for x in a)
    lt = sum((x < b).sum() for x in a)
    return (gt - lt) / (len(a) * len(b))


def analyze_fragment(fragdir, writer=None):
    import tifffile
    from PIL import Image
    Image.MAX_IMAGE_PIXELS = None   # fragment IR photos exceed PIL's default
    name = os.path.basename(fragdir.rstrip("/"))
    stack = np.stack([
        tifffile.imread(os.path.join(fragdir, f)).astype(np.float64)
        for f in LAYERS])
    fmask = np.asarray(Image.open(os.path.join(fragdir, "mask.png")).convert("L")) > 127
    ir_path = os.path.join(fragdir, "ir.png")
    if MASK_SOURCE == "inklabels":
        # POST-HOC SENSITIVITY MODE (not the preregistered analysis): the
        # declared Otsu-on-IR mask failed on several fragments (ink
        # fractions up to 0.95 — it thresholds illumination, not ink).
        # This mode uses the provider-aligned inklabels for every fragment.
        ink = np.asarray(Image.open(os.path.join(fragdir, "inklabels.png"))
                         .convert("L")) > 127
        ink &= fmask
    elif os.path.exists(ir_path):
        ir = np.asarray(Image.open(ir_path).convert("L"), dtype=np.float64)
        if ir.shape != stack.shape[1:]:
            raise SystemExit(f"{name}: ir.png {ir.shape} != volume {stack.shape[1:]}")
        thr = otsu(ir[fmask])
        ink = (ir < thr) & fmask
    else:
        # Amendment 1: Frag4 ships no aligned ir.png; its provider-derived
        # aligned inklabels.png is the ink mask for that fragment only.
        ink = np.asarray(Image.open(os.path.join(fragdir, "inklabels.png"))
                         .convert("L")) > 127
        ink &= fmask
        print(f"{name}: no ir.png — using inklabels.png per amendment 1")
    edge_dist = ndi.distance_transform_edt(fmask)

    ny, nx = fmask.shape
    ink_feats, ctl_feats, ink_ed, ctl_ed = [], [], [], []
    for y0 in range(0, ny - PATCH + 1, PATCH):
        for x0 in range(0, nx - PATCH + 1, PATCH):
            sl = np.s_[y0:y0 + PATCH, x0:x0 + PATCH]
            if fmask[sl].mean() < 0.95:
                continue
            frac = ink[sl].mean()
            feats = None
            if frac >= 0.70:
                feats = patch_stats(stack[(np.s_[:],) + sl])
                ink_feats.append(feats); ink_ed.append(np.median(edge_dist[sl]))
            elif frac == 0.0:
                feats = patch_stats(stack[(np.s_[:],) + sl])
                ctl_feats.append(feats); ctl_ed.append(np.median(edge_dist[sl]))
    ink_feats = np.array(ink_feats); ctl_feats = np.array(ctl_feats)
    ink_ed = np.array(ink_ed); ctl_ed = np.array(ctl_ed)
    if len(ink_feats) < 10 or len(ctl_feats) < 10:
        print(f"{name}: UNDERPOWERED (ink {len(ink_feats)}, "
              f"control {len(ctl_feats)}) — reported, not tested")
        return None

    # nearest-neighbor edge-distance matching, one control per ink patch
    order = np.argsort(ink_ed)
    avail = np.ones(len(ctl_ed), bool)
    chosen = []
    for i in order:
        cand = np.nonzero(avail)[0]
        if not len(cand):
            break
        j = cand[np.argmin(np.abs(ctl_ed[cand] - ink_ed[i]))]
        avail[j] = False
        chosen.append(j)
    ctl_m = ctl_feats[chosen]

    out = {"fragment": name, "n_ink": len(ink_feats), "n_control": len(ctl_m)}
    for k, sname in enumerate(["S1_variance", "S2_kurtosis",
                               "S3_lag1_autocorr", "S4_hf_fraction"]):
        u = stats.mannwhitneyu(ink_feats[:, k], ctl_m[:, k],
                               alternative="two-sided")
        d = cliffs_delta(ink_feats[:, k], ctl_m[:, k])
        sig = u.pvalue < ALPHA / N_TESTS
        out[sname] = {"p": float(u.pvalue), "cliffs_delta": float(d),
                      "significant": bool(sig),
                      "sign": int(np.sign(d)) if sig else 0}
        if writer:
            writer.writerow([name, sname, len(ink_feats), len(ctl_m),
                             f"{u.pvalue:.3e}", f"{d:+.3f}", sig])
    return out


def verdict(results):
    """Preregistered criteria: FINGERPRINT iff >=1 statistic significant in
    >=4 of 6 fragments with the same sign in all significant cases."""
    stats_names = ["S1_variance", "S2_kurtosis", "S3_lag1_autocorr",
                   "S4_hf_fraction"]
    fingerprint = False
    for sname in stats_names:
        signs = [r[sname]["sign"] for r in results
                 if r is not None and r[sname]["significant"]]
        if len(signs) >= 4 and len(set(signs)) == 1:
            fingerprint = True
            print(f"  {sname}: significant same-sign in {len(signs)} fragments")
    print("VERDICT:", "FINGERPRINT" if fingerprint else "NULL",
          f"({sum(r is not None for r in results)} fragments tested; "
          f"criteria: >=1 stat significant same-sign in >=4 of 6)")
    return fingerprint


def selftest():
    """Blind checks on synthetic data — no real fragments touched.
    (a) Pure-noise null: fabricated ink mask on iid noise -> ~uniform p.
    (b) Planted signal: ink regions get extra high-frequency texture ->
        S4 must detect with positive delta."""
    rng = np.random.default_rng(7)
    n_sig = 0
    ps = []
    for rep in range(20):
        stack = rng.normal(size=(5, 512, 512))
        feats_a = [patch_stats(stack[:, y:y+64, x:x+64])
                   for y in range(0, 512, 64) for x in range(0, 512, 64)
                   if (y // 64 + x // 64) % 2 == 0]
        feats_b = [patch_stats(stack[:, y:y+64, x:x+64])
                   for y in range(0, 512, 64) for x in range(0, 512, 64)
                   if (y // 64 + x // 64) % 2 == 1]
        a, b = np.array(feats_a), np.array(feats_b)
        for k in range(4):
            p = stats.mannwhitneyu(a[:, k], b[:, k]).pvalue
            ps.append(p)
            n_sig += p < ALPHA / N_TESTS
    print(f"null selftest: {n_sig}/{len(ps)} significant "
          f"(expect ~0), p median {np.median(ps):.2f} (expect ~0.5)")
    assert n_sig <= 1

    stack = rng.normal(size=(5, 512, 512))
    hf = rng.normal(size=(5, 512, 512))
    hf -= ndi.gaussian_filter(hf, (0, 2, 2))     # high-pass texture
    planted = stack.copy()
    planted[:, :, :256] += 0.6 * hf[:, :, :256]  # "ink" = left half
    a = np.array([patch_stats(planted[:, y:y+64, x:x+64])
                  for y in range(0, 512, 64) for x in range(0, 256, 64)])
    b = np.array([patch_stats(planted[:, y:y+64, x:x+64])
                  for y in range(0, 512, 64) for x in range(256, 512, 64)])
    p4 = stats.mannwhitneyu(a[:, 3], b[:, 3]).pvalue
    d4 = cliffs_delta(a[:, 3], b[:, 3])
    print(f"planted-signal selftest: S4 p={p4:.2e}, delta={d4:+.2f} "
          f"(must be significant, positive)")
    assert p4 < ALPHA / N_TESTS and d4 > 0
    print("SELFTEST PASS")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("fragdirs", nargs="*")
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--csv", default="listening_pilot_results.csv")
    ap.add_argument("--mask-source", choices=["ir", "inklabels"], default="ir")
    a = ap.parse_args()
    MASK_SOURCE = a.mask_source
    if a.selftest:
        selftest()
        sys.exit(0)
    if not a.fragdirs:
        ap.error("give fragment dirs or --selftest")
    results = []
    with open(a.csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["fragment", "statistic", "n_ink", "n_control",
                    "p_value", "cliffs_delta", "significant"])
        for d in a.fragdirs:
            r = analyze_fragment(d, w)
            results.append(r)
            if r:
                print(r["fragment"], {k: v for k, v in r.items()
                                      if k.startswith("S")})
    verdict(results)
    print("wrote", a.csv)
