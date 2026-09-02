#!/usr/bin/env python3
"""Winding coordinate from a scroll CT slice, by interferometric phase methods.

Idea (the ONUSKA wedge): a scroll cross-section is an interferogram. The rolled
papyrus sheet produces a quasi-periodic layered pattern; one winding = one fringe.
So "relative winding number" — the annotation the Vesuvius spiral-fit pipeline
needs and currently gets by hand — is fringe order, and the machinery of RF/SAR
interferometry applies:

  1. structure tensor  -> local sheet orientation + coherence (quality map)
  2. Riesz/monogenic quadrature steered along the sheet normal -> local spatial
     frequency k (rad/px): |k| = 2*pi / local winding spacing
  3. orient k radially (sign fix), then solve the weighted least-squares phase
     integration  min sum w |grad(phi) - k|^2   (the classic Poisson / phase-
     unwrapping step from interferometry)
  4. phi / 2*pi is a continuous WINDING COORDINATE: iso-contours are individual
     windings; differences between any two points are relative winding numbers.

Output: winding-coordinate map + quality map + auto-generated relative-winding
annotations, validated against direct sheet-crossing counts along radial rays.

Stdlib + numpy/scipy only. 2D proof of concept on one slice; the 3D extension
stacks slices with z-continuity constraints.
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage as ndi
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import cg


# ---------------------------------------------------------------- structure tensor
def structure_tensor(img: np.ndarray, sg: float = 1.2, st: float = 5.0):
    """Return sheet-normal angle theta (rad, mod pi) and coherence in [0,1]."""
    gx = ndi.gaussian_filter(img, sg, order=(0, 1))
    gy = ndi.gaussian_filter(img, sg, order=(1, 0))
    Jxx = ndi.gaussian_filter(gx * gx, st)
    Jyy = ndi.gaussian_filter(gy * gy, st)
    Jxy = ndi.gaussian_filter(gx * gy, st)
    theta = 0.5 * np.arctan2(2 * Jxy, Jxx - Jyy)  # dominant gradient direction
    tr = Jxx + Jyy
    det = np.sqrt((Jxx - Jyy) ** 2 + 4 * Jxy ** 2)
    coherence = np.where(tr > 1e-9, det / (tr + 1e-9), 0.0)
    return theta, coherence


# ---------------------------------------------------------------- local frequency
def _band_frequency(img_f: np.ndarray, theta: np.ndarray, s1: float, s2: float):
    f = ndi.gaussian_filter(img_f, s1) - ndi.gaussian_filter(img_f, s2)
    F = np.fft.fft2(f)
    ny, nx = f.shape
    u = np.fft.fftfreq(nx)[None, :]
    v = np.fft.fftfreq(ny)[:, None]
    q = np.sqrt(u * u + v * v)
    q[0, 0] = 1.0
    R1 = np.real(np.fft.ifft2(1j * u / q * F))
    R2 = np.real(np.fft.ifft2(1j * v / q * F))
    nxv, nyv = np.cos(theta), np.sin(theta)
    qd = nxv * R1 + nyv * R2
    amp = np.hypot(f, qd)
    fx = ndi.sobel(f, axis=1) / 8.0
    fy = ndi.sobel(f, axis=0) / 8.0
    qx = ndi.sobel(qd, axis=1) / 8.0
    qy = ndi.sobel(qd, axis=0) / 8.0
    dfn = nxv * fx + nyv * fy
    dqn = nxv * qx + nyv * qy
    k = (f * dqn - qd * dfn) / (f * f + qd * qd + 1e-6)
    return np.abs(k), amp


def local_frequency(img: np.ndarray, theta: np.ndarray,
                    bands=((0.8, 3.0), (1.6, 6.0), (3.2, 12.0), (6.0, 24.0)),
                    k_prior: np.ndarray | None = None, lock_tol: float = 0.6):
    """Multi-band steered-quadrature frequency: per pixel, take the frequency of
    the band with the strongest local fringe amplitude (then median-denoise).
    A single wide band biases k low in densely packed zones; the dominant-band
    pick keeps high-frequency (tight-winding) regions honest."""
    img_f = img.astype(np.float64)
    ks, amps = [], []
    for s1, s2 in bands:
        k_b, a_b = _band_frequency(img_f, theta, s1, s2)
        # smooth amplitude for stable band selection
        ks.append(k_b)
        amps.append(ndi.gaussian_filter(a_b, 2.0))
    K = np.stack(ks)
    A = np.stack(amps)
    if k_prior is None:
        sel = np.argmax(A, axis=0)
    else:
        # Pitch lock (field-frame harmonic suppression). At fine resolution the
        # two faces of one sheet resolve as separate fringes and the strongest
        # band sits at twice the true winding frequency (measured: 1.39
        # fringes/winding at 9.6 um, validation/fringe_scale). Given a prior
        # wavenumber per pixel (e.g. from a coarser solve), prefer the band
        # whose frequency is within `lock_tol` (log units) of the prior,
        # scored by amplitude; fall back to plain argmax where no band is.
        with np.errstate(divide="ignore", invalid="ignore"):
            dist = np.abs(np.log(K) - np.log(k_prior)[None])
        ok = np.isfinite(dist) & (dist < lock_tol)
        score = np.where(ok, A / (1.0 + dist), -1.0)
        sel_lock = np.argmax(score, axis=0)
        sel_amp = np.argmax(A, axis=0)
        sel = np.where(ok.any(axis=0), sel_lock, sel_amp)
    iy, ix = np.mgrid[0:img.shape[0], 0:img.shape[1]]
    kmag = K[sel, iy, ix]
    amp = A[sel, iy, ix]
    kmag = ndi.median_filter(kmag, size=5)
    return kmag, amp


# ------------------------------------------------------------- phase integration
def integrate_phase(kx: np.ndarray, ky: np.ndarray, w: np.ndarray, mask: np.ndarray,
                    tol: float = 1e-7, maxiter: int = 6000, x0=None, pairs=None):
    """Solve the weighted LS phase problem on the masked grid:
         min  sum_edges w |grad(phi) - k|^2
            + sum_pairs wc (phi_b - phi_a - delta)^2      (ridge-adjacency terms)
    The pairwise terms encode 'next sheet outward = +2π' — auto-generated
    relative winding annotations."""
    ny, nx = mask.shape
    idx = -np.ones(mask.shape, dtype=np.int64)
    ids = np.flatnonzero(mask.ravel())
    idx.ravel()[ids] = np.arange(ids.size)
    n = ids.size
    rows, cols, vals = [], [], []
    b = np.zeros(n)
    yy, xx = np.nonzero(mask)
    wk = w
    for dy, dx, kcomp in ((0, 1, kx), (1, 0, ky)):
        y2, x2 = yy + dy, xx + dx
        ok = (y2 < ny) & (x2 < nx)
        ok &= mask[y2 % ny, x2 % nx] if False else mask[np.minimum(y2, ny - 1), np.minimum(x2, nx - 1)]
        yA, xA = yy[ok], xx[ok]
        yB, xB = yA + dy, xA + dx
        iA = idx[yA, xA]
        iB = idx[yB, xB]
        we = 0.5 * (wk[yA, xA] + wk[yB, xB]) + 1e-4
        ke = 0.5 * (kcomp[yA, xA] + kcomp[yB, xB])
        rows += [iA, iA, iB, iB]
        cols += [iA, iB, iB, iA]
        vals += [we, -we, we, -we]
        np.add.at(b, iA, -we * ke)
        np.add.at(b, iB, we * ke)
    if pairs is not None and len(pairs[0]) > 0:
        ya, xa, yb, xb, wc, delta = pairs
        ok = mask[ya, xa] & mask[yb, xb]
        iA = idx[ya[ok], xa[ok]]
        iB = idx[yb[ok], xb[ok]]
        wcv = np.asarray(wc, float)[ok] if np.ndim(wc) else np.full(iA.size, float(wc))
        dv = np.asarray(delta, float)[ok] if np.ndim(delta) else np.full(iA.size, float(delta))
        rows += [iA, iA, iB, iB]
        cols += [iA, iB, iB, iA]
        vals += [wcv, -wcv, wcv, -wcv]
        np.add.at(b, iB, wcv * dv)
        np.add.at(b, iA, -wcv * dv)
    A = coo_matrix((np.concatenate([np.asarray(v, float) for v in vals]),
                    (np.concatenate(rows), np.concatenate(cols))), shape=(n, n)).tocsr()
    # The masked domain can shatter into many disconnected components (debris,
    # erosion specks, coarse-level subsampling). Each component has its own
    # floating gauge, so the full system is singular and CG stalls. Solve only
    # the largest component; everything else gets NaN.
    from scipy.sparse.csgraph import connected_components
    ncomp, labels = connected_components(A, directed=False)
    if ncomp > 1:
        import sys as _sys
        big = np.argmax(np.bincount(labels))
        keep = labels == big
        print(f"  (mask has {ncomp} disconnected components; solving largest "
              f"= {int(keep.sum())}/{n} px)", file=_sys.stderr)
        A = A[keep][:, keep].tocsr()
        b = b[keep]
        if x0 is not None:
            x0 = x0[keep]
    else:
        keep = None
    # anchor gauge
    A[0, 0] += 1.0
    # Preconditioner: algebraic multigrid when available (essential for
    # large domains — Jacobi-PCG stalls beyond ~5M unknowns), Jacobi fallback.
    M = None
    if A.shape[0] > 200_000:
        try:
            import pyamg
            ml = pyamg.smoothed_aggregation_solver(A.tocsr(), max_coarse=500)
            M = ml.aspreconditioner(cycle='V')
        except Exception as e:
            import sys
            print(f"pyamg unavailable ({e}); falling back to Jacobi", file=sys.stderr)
    if M is None:
        from scipy.sparse import diags
        d = A.diagonal()
        M = diags(1.0 / np.maximum(d, 1e-8))
    phi_flat, info = cg(A, b, rtol=tol, maxiter=maxiter, M=M, x0=x0)
    full_flat = np.full(n, np.nan)
    if keep is None:
        full_flat[:] = phi_flat
    else:
        full_flat[keep] = phi_flat
    phi = np.full(mask.shape, np.nan)
    phi.ravel()[ids] = full_flat
    return phi, info


def integrate_phase_multiscale(kx, ky, w, mask, factors=None, pairs=None,
                               finest_maxiter=25000):
    """Coarse-to-fine weighted phase integration: solve small, upsample, refine.
    Ridge-pair constraints are applied at the finest level. For large domains
    (>2M masked pixels) a deeper pyramid is used so the finest solve starts
    from a good warm start; finest_maxiter bounds the final refinement."""
    if factors is None:
        factors = (8, 4, 2, 1) if mask.sum() > 2_000_000 else (4, 2, 1)
    phi0 = None
    for f in factors:
        if f > 1:
            sl = (slice(None, None, f), slice(None, None, f))
            # dilate-then-sample keeps the coarse mask connected where the
            # fine mask is (strided sampling shatters it into islands)
            mdil = ndi.maximum_filter(mask.astype(np.uint8), size=f) > 0
            kxs, kys, ws, ms = kx[sl] * f, ky[sl] * f, w[sl], mdil[sl]
        else:
            kxs, kys, ws, ms = kx, ky, w, mask
        x0 = None
        if phi0 is not None:
            up = ndi.zoom(np.nan_to_num(phi0), (ms.shape[0] / phi0.shape[0],
                                                ms.shape[1] / phi0.shape[1]), order=1)
            x0 = up[ms]
            if not np.all(np.isfinite(x0)):
                x0 = np.nan_to_num(x0)
        phi, info = integrate_phase(kxs, kys, ws, ms, x0=x0,
                                    maxiter=3000 if f > 1 else finest_maxiter,
                                    pairs=pairs if f == 1 else None)
        if f == 1 and info != 0:
            import sys
            print(f"WARNING: finest-level CG not converged (info={info})",
                  file=sys.stderr)
        phi0 = phi
    return phi0, info


def ridge_pairs(img: np.ndarray, theta_out: tuple, amp: np.ndarray, mask: np.ndarray,
                min_gap: int = 3, max_gap: int = 45, amp_q: float = 55.0,
                kmag: np.ndarray = None, min_frac: float = 0.0):
    """Detect sheet ridges and, for each ridge pixel, march outward along the
    sheet normal to the next ridge — emitting a (+2π between these two points)
    constraint. These ARE relative winding annotations, generated automatically."""
    nxv, nyv = theta_out
    f = ndi.gaussian_filter(img.astype(np.float64), 1.0)
    # local maxima along the normal: f greater than both normal-shifted copies
    ny, nx = f.shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    def samp(field, dy, dx):
        return ndi.map_coordinates(field, [np.clip(yy + dy, 0, ny - 1),
                                           np.clip(xx + dx, 0, nx - 1)], order=1)
    fp = samp(f, nyv * 1.5, nxv * 1.5)
    fm = samp(f, -nyv * 1.5, -nxv * 1.5)
    thr = np.percentile(amp[mask], amp_q)
    ridge = (f > fp) & (f > fm) & (amp > thr) & mask
    ry, rx = np.nonzero(ridge)
    if ry.size == 0:
        return (np.array([], int),) * 4 + (np.array([]),) * 2
    # march outward
    dirx = nxv[ry, rx]
    diry = nyv[ry, rx]
    hit_y = np.full(ry.size, -1, int)
    hit_x = np.full(ry.size, -1, int)
    alive = np.ones(ry.size, bool)
    # damage-aware gate: a ridge may only march ~3.5 local wavelengths before we
    # stop trusting that the next ridge found is the adjacent winding
    if kmag is not None:
        lam = 2 * np.pi / np.maximum(kmag[ry, rx], 1e-3)
        tmax = np.clip(3.5 * lam, min_gap + 2, max_gap).astype(int)
        # face gate: with a (locked) local wavelength, a ridge closer than
        # min_frac * lambda is the other face of the SAME sheet, not the next
        # winding — march past it instead of asserting +2pi. Inert at 0.0.
        tmin = np.maximum(min_gap, (min_frac * lam).astype(int)) if min_frac > 0 else np.full(ry.size, min_gap)
    else:
        tmax = np.full(ry.size, max_gap)
        tmin = np.full(ry.size, min_gap)
    for t in range(min_gap, max_gap + 1):
        qy = np.clip((ry + diry * t).round().astype(int), 0, ny - 1)
        qx = np.clip((rx + dirx * t).round().astype(int), 0, nx - 1)
        is_hit = alive & (t >= tmin) & (t <= tmax) & ridge[qy, qx]
        hit_y[is_hit] = qy[is_hit]
        hit_x[is_hit] = qx[is_hit]
        alive &= ~is_hit
        alive &= t < tmax
        if not alive.any():
            break
    got = hit_y >= 0
    wc = np.full(got.sum(), 0.6)
    delta = np.full(got.sum(), 2 * np.pi)
    return ry[got], rx[got], hit_y[got], hit_x[got], wc, delta


# ---------------------------------------------------------------------- pipeline
def winding_coordinate(img: np.ndarray, mask: np.ndarray, core: tuple[float, float],
                       kmin: float = 2 * np.pi / 60, kmax: float = 2 * np.pi / 3,
                       bands=None, ridge_max_gap: int = 45,
                       k_prior: np.ndarray | None = None, lock_tol: float = 0.6,
                       ridge_min_frac: float = 0.0):
    theta, coh = structure_tensor(img)
    kw = {"k_prior": k_prior, "lock_tol": lock_tol}
    kmag, amp = (local_frequency(img, theta, bands=bands, **kw) if bands
                 else local_frequency(img, theta, **kw))
    kmag = np.clip(kmag, kmin, kmax)
    # sign: orient the normal outward from the core
    ny, nx = img.shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    rx, ry = xx - core[0], yy - core[1]
    nxv, nyv = np.cos(theta), np.sin(theta)
    flip = (nxv * rx + nyv * ry) < 0
    nxv = np.where(flip, -nxv, nxv)
    nyv = np.where(flip, -nyv, nyv)
    kx = kmag * nxv
    ky = kmag * nyv
    a95 = np.percentile(amp[mask], 95) + 1e-9
    w = coh * np.clip(amp / a95, 0, 1)
    w = np.where(mask, w, 0.0)
    pairs = ridge_pairs(img, (nxv, nyv), amp, mask, max_gap=ridge_max_gap, kmag=kmag,
                        min_frac=ridge_min_frac)
    n_pairs = len(pairs[0])
    phi, info = integrate_phase_multiscale(kx, ky, w, mask, pairs=pairs)
    print(f"  ridge-adjacency constraints: {n_pairs}")
    return phi / (2 * np.pi), w, theta, kmag, info


# -------------------------------------------------------------------- validation
def ray_crossings(img: np.ndarray, mask: np.ndarray, core, ang, step=0.5):
    """Count sheet crossings (intensity peaks) along a ray from the core."""
    from scipy.signal import find_peaks
    ny, nx = img.shape
    ts = np.arange(0, max(ny, nx), step)
    xs = core[0] + ts * np.cos(ang)
    ys = core[1] + ts * np.sin(ang)
    ok = (xs >= 0) & (xs < nx - 1) & (ys >= 0) & (ys < ny - 1)
    xs, ys = xs[ok], ys[ok]
    inside = mask[ys.astype(int), xs.astype(int)]
    if inside.sum() < 10:
        return None
    last = np.max(np.flatnonzero(inside))
    xs, ys = xs[:last + 1], ys[:last + 1]
    prof = ndi.map_coordinates(img.astype(float), [ys, xs], order=1)
    prof_s = ndi.gaussian_filter1d(prof, 2.0)
    thr = np.percentile(prof_s, 60)
    peaks, _ = find_peaks(prof_s, height=thr, distance=7)
    return xs, ys, peaks, prof_s


if __name__ == "__main__":
    import sys, time
    from PIL import Image
    src = sys.argv[1] if len(sys.argv) > 1 else "v_l3_mid.npy"
    img = np.load(src).astype(np.float64)
    mask = img > 0
    mask = ndi.binary_erosion(mask, iterations=3)
    # core = darkest cavity near centroid of mask… use provided or centroid
    ys, xs = np.nonzero(mask)
    core = (float(np.mean(xs)), float(np.mean(ys)))
    t0 = time.time()
    W, w, theta, kmag, info = winding_coordinate(img, mask, core)
    print(f"winding solve: {time.time()-t0:.1f}s  cg_info={info}  "
          f"winding range {np.nanmin(W):.1f}..{np.nanmax(W):.1f}")
    np.save("winding_coord.npy", W)
    np.save("winding_quality.npy", w)

    # validation along 8 rays
    print("\nray validation (peaks counted vs winding-coordinate span):")
    for ang in np.linspace(0, 2 * np.pi, 8, endpoint=False):
        r = ray_crossings(img, mask, core, ang)
        if r is None:
            continue
        xs_, ys_, peaks, prof = r
        Wray = ndi.map_coordinates(np.nan_to_num(W), [ys_, xs_], order=1)
        span = abs(Wray[-1] - Wray[np.argmin(np.hypot(xs_ - core[0], ys_ - core[1]))])
        print(f"  ang {np.degrees(ang):5.1f}°: peaks={len(peaks):3d}   |ΔW|={span:6.1f}")
