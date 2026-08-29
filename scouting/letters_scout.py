"""Scout candidate First Letters volumes: stream a mid-volume slice at a preview
level, run the winding field zero-shot, and report convergence quality."""
import sys, json, time
import numpy as np
import zarr, s3fs
sys.path.insert(0, 'scroll-interferometry')
from winding_phase import winding_coordinate

BUCKET = 'vesuvius-challenge-open-data'
VOLS = {
    'PHerc1203': '20250820131727-9.362um-1.2m-113keV-masked.zarr',
    'PHerc1218': '20250521120456-8.640um-1.2m-116keV-masked.zarr',
    'PHerc0125': '20250821151825-9.362um-1.2m-113keV-masked.zarr',
    'PHerc0257': '20250821151750-9.362um-1.2m-113keV-masked.zarr',
}
LEVEL = '2'   # ~4x downsample -> ~35-37um preview; fast smoke of scan character

fs = s3fs.S3FileSystem(anon=True)
report = {}
for scroll, vol in VOLS.items():
    t0 = time.time()
    try:
        store = f"https://vesuvius-challenge-open-data.s3.us-east-1.amazonaws.com/{scroll}/volumes/{vol}/{LEVEL}"
        arr = zarr.open(store, mode='r')
        zmid = arr.shape[0] // 2
        sl = arr[zmid].astype(np.float32)
        nz = sl[sl > 0]
        frac = len(nz) / sl.size
        stats = dict(shape=list(arr.shape), nonzero_frac=round(float(frac), 3),
                     p20=float(np.percentile(nz, 20)), p90=float(np.percentile(nz, 90)))
        from scipy import ndimage as ndi
        img = sl.astype(np.float64)
        mask = ndi.binary_erosion(img > 0, iterations=3)
        ys, xs = np.nonzero(mask)
        core = (float(xs.mean()), float(ys.mean()))
        W, q, theta, kmag, info = winding_coordinate(img, mask, core)
        wv = W[np.isfinite(W) & mask]
        span = float(np.nanmax(wv) - np.nanmin(wv)) if wv.size else 0.0
        qmed = float(np.nanmedian(q[mask])) if q is not None else -1
        stats.update(winding_span=round(span, 1), cg_info=int(info),
                     quality_median=round(qmed, 3), secs=round(time.time() - t0, 1))
        report[scroll] = stats
        print(scroll, stats, flush=True)
    except Exception as e:
        report[scroll] = {'error': f'{type(e).__name__}: {e}'}
        print(scroll, 'ERR', type(e).__name__, str(e)[:200], flush=True)
json.dump(report, open('letters_scout_report.json', 'w'), indent=1)
