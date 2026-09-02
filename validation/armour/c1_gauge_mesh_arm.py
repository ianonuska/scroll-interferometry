"""Armour C1 (declared first): score the same field against TWO ground-truth arms of the
constraint gauge: (a) the point file we assembled (gt_accepted_1667_clean.json) and (b) the
gauge's own --gt-mesh loader over the accepted meshes directly (independently derived).
Same adapter (W_L2_locked2, step 9), same pitch-um 224, um-per-vox 2.399. Report M1/M2
from both; the claim is armoured if both arms agree on direction (locked better than static)
and within ~20% on magnitude."""
import json, subprocess, sys, numpy as np
ROOT = "/Users/ianonuska/projects/vesuvius/local1667"; G = "/Users/ianonuska/projects/vesuvius/constraint-gauge"
SCALE = 4; ZL = 4634; step = 9
for tag in ("L2_locked2", "L2"):
    W = np.load(f"{ROOT}/W_{tag}.npy").astype(float); Q = np.load(f"{ROOT}/Q_{tag}.npy").astype(float); mask = np.load(f"{ROOT}/mask_L2.npy")
    good = mask & np.isfinite(W); sel = np.zeros_like(good); sel[::step, ::step] = True; yy, xx = np.nonzero(good & sel)
    pts = np.stack([xx*SCALE, yy*SCALE, np.full(len(xx), ZL*SCALE)], 1).astype(float)
    json.dump({"name": f"windingmeter-1667-z4634-9.6um-{tag}", "points_xyz": pts.tolist(), "winding": W[yy, xx].tolist(), "conf": Q[yy, xx].tolist()}, open(f"{ROOT}/adapter_c1_{tag}.json", "w"))
    for arm, args in (("points", ["--gt", f"{ROOT}/gt_accepted_1667_clean.json"]), ("meshes", ["--gt-mesh", f"{ROOT}/meshes"])):
        r = subprocess.run([sys.executable, f"{G}/run_gauge.py", *args, "--adapter", f"json:{ROOT}/adapter_c1_{tag}.json", "--subject", f"wm-{tag}", "--gt-arm", f"1667-{arm}", "--pitch-um", "224", "--um-per-vox", "2.399", "--max-pairs", "200000", "--out-prefix", f"{ROOT}/c1_{tag}_{arm}"], capture_output=True, text=True, cwd=G)
        print(f"=== {tag} vs {arm} arm (rc {r.returncode})"); print("\n".join(l for l in r.stdout.splitlines() if any(k in l for k in ("M1","M2","M3","M4","tau","density","wrote","pairs","collections"))))
        if r.returncode: print(r.stderr[-600:])
print("C1 DONE")
