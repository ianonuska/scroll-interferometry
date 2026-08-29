#!/usr/bin/env python3
"""Z-coupled slab solve: joint winding coordinate over a stack of nearby slices.

Extends the per-slice interferometric solve to a thin 3D slab: each slice keeps
its in-plane data terms (steered-quadrature k field + ridge-adjacency 2π pairs),
and adjacent slices are coupled by a soft z-continuity term
    wz * (phi[s+1](y,x) - phi[s](y,x))^2
reflecting that windings drift slowly along the scroll axis. One sparse LS
system, warm-started from gauge-aligned per-slice solutions.

Outputs a z-consistent winding volume — the natural precursor to full-3D
constraints for spiral fitting — plus a fair before/after consistency metric
(independent solves vs joint solve on the same slices).
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage as ndi
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import cg

from winding_phase import (structure_tensor, local_frequency, ridge_pairs,
                           integrate_phase_multiscale)


def slice_terms(img, mask, core):
    """Per-slice data terms: kx, ky, weights, ridge pairs."""
    theta, coh = structure_tensor(img)
    kmag, amp = local_frequency(img, theta)
    kmag = np.clip(kmag, 2 * np.pi / 60, 2 * np.pi / 3)
    ny, nx = img.shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    rx, ry = xx - core[0], yy - core[1]
    nxv, nyv = np.cos(theta), np.sin(theta)
    flip = (nxv * rx + nyv * ry) < 0
    nxv = np.where(flip, -nxv, nxv)
    nyv = np.where(flip, -nyv, nyv)
    a95 = np.percentile(amp[mask], 95) + 1e-9
    w = np.where(mask, coh * np.clip(amp / a95, 0, 1), 0.0)
    pairs = ridge_pairs(img, (nxv, nyv), amp, mask, kmag=kmag)
    return kmag * nxv, kmag * nyv, w, pairs


def solve_independent(img, mask, core):
    kx, ky, w, pairs = slice_terms(img, mask, core)
    phi, info = integrate_phase_multiscale(kx, ky, w, mask, pairs=pairs)
    return phi


def solve_slab(imgs, masks, cores, wz_rel=0.4, x0_slices=None,
               tol=1e-7, maxiter=9000):
    """Joint solve over S slices. Returns list of phi arrays (winding * 2π)."""
    S = len(imgs)
    ny, nx = imgs[0].shape
    terms = [slice_terms(imgs[s], masks[s], cores[s]) for s in range(S)]

    # global indexing
    idx = np.full((S, ny, nx), -1, dtype=np.int64)
    n = 0
    for s in range(S):
        ids = np.flatnonzero(masks[s].ravel())
        idx[s].ravel()[ids] = np.arange(n, n + ids.size)
        n += ids.size

    rows, cols, vals = [], [], []
    b = np.zeros(n)

    def add_pair(iA, iB, we, delta):
        rows.extend([iA, iA, iB, iB])
        cols.extend([iA, iB, iB, iA])
        vals.extend([we, -we, we, -we])
        np.add.at(b, iB, we * delta)
        np.add.at(b, iA, -we * delta)

    wz_all = []
    for s in range(S):
        kx, ky, w, pairs = terms[s]
        m = masks[s]
        yy, xx = np.nonzero(m)
        for dy, dx, kcomp in ((0, 1, kx), (1, 0, ky)):
            y2 = yy + dy
            x2 = xx + dx
            ok = (y2 < ny) & (x2 < nx)
            ok &= m[np.minimum(y2, ny - 1), np.minimum(x2, nx - 1)]
            yA, xA = yy[ok], xx[ok]
            yB, xB = yA + dy, xA + dx
            we = 0.5 * (w[yA, xA] + w[yB, xB]) + 1e-4
            ke = 0.5 * (kcomp[yA, xA] + kcomp[yB, xB])
            add_pair(idx[s, yA, xA], idx[s, yB, xB], we, ke)
        ya, xa, yb, xb, wc, delta = pairs
        ok = m[ya, xa] & m[yb, xb]
        if ok.any():
            add_pair(idx[s, ya[ok], xa[ok]], idx[s, yb[ok], xb[ok]],
                     np.asarray(wc)[ok], np.asarray(delta)[ok])
        wz_all.append(np.mean(w[m]))
        # z-coupling to next slice
        if s + 1 < S:
            both = m & masks[s + 1]
            yy2, xx2 = np.nonzero(both)
            wzv = np.full(yy2.size, wz_rel * np.mean(w[m]))
            add_pair(idx[s, yy2, xx2], idx[s + 1, yy2, xx2], wzv,
                     np.zeros(yy2.size))

    A = coo_matrix((np.concatenate([np.asarray(v, float) for v in vals]),
                    (np.concatenate(rows), np.concatenate(cols))),
                   shape=(n, n)).tocsr()
    A[0, 0] += 1.0
    M = diags(1.0 / np.maximum(A.diagonal(), 1e-8))

    x0 = None
    if x0_slices is not None:
        # gauge-align warm starts by chaining median offsets
        aligned = [x0_slices[0]]
        for s in range(1, S):
            both = masks[s - 1] & masks[s] & np.isfinite(x0_slices[s - 1]) & np.isfinite(x0_slices[s])
            off = np.median((aligned[-1] - x0_slices[s])[both])
            aligned.append(x0_slices[s] + off)
        x0 = np.concatenate([np.nan_to_num(aligned[s])[masks[s]] for s in range(S)])

    phi_flat, info = cg(A, b, rtol=tol, maxiter=maxiter, M=M, x0=x0)
    out = []
    for s in range(S):
        phi = np.full((ny, nx), np.nan)
        phi[masks[s]] = phi_flat[idx[s][masks[s]]]
        out.append(phi)
    return out, info


def adjacent_consistency(fields, masks):
    """Median |ΔW| between adjacent slices after per-pair gauge alignment."""
    stats = []
    for s in range(len(fields) - 1):
        both = masks[s] & masks[s + 1] & np.isfinite(fields[s]) & np.isfinite(fields[s + 1])
        d = (fields[s] - fields[s + 1])[both]
        r = np.abs(d - np.median(d))
        stats.append((float(np.median(r)), float(np.percentile(r, 90))))
    return stats
