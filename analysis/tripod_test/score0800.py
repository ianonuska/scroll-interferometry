"""Tripod test scoring — exactly as preregistered (2026-09-02):
S_switch = fraction of a mesh's sampled vertices whose W deviates from the
mesh median by >0.5 windings; S_mad = MAD of W. Primary: AUC(S_switch,
aprova vs reprova) with bootstrap CI; pass = AUC>=0.75 and 5th pct>0.5.
Secondary: Spearman(S_switch, wrap number) over all 23."""
import csv, json, glob, os
import numpy as np, tifffile
from scipy import stats
IDX = "/workspace/eligible/data/index.csv"; MESH = "/workspace/eligible"; OUT = "/workspace/p0800"
rows = [r for r in csv.DictReader(open(IDX)) if r["scroll"] == "PHerc0800" and r["gate_verdict"]]
res = []
for r in rows:
    wz = int(r["window_z"]); zs = wz + 400
    Wp, Qp = f"{OUT}/L1lock_z{zs}_W.npy", f"{OUT}/L1lock_z{zs}_Q.npy"
    if not os.path.exists(Wp): print("missing station", zs); continue
    W = np.load(Wp); Q = np.load(Qp)
    d = f"{MESH}/{r['path']}"
    z = tifffile.imread(d + "/z.tif"); sel = (z > 0) & (np.abs(z - zs) <= 24)
    if sel.sum() < 30:
        print(f"{r['path']}: only {sel.sum()} vertices in slab"); continue
    x = tifffile.imread(d + "/x.tif")[sel] / 2; y = tifffile.imread(d + "/y.tif")[sel] / 2
    xi = np.clip(x.astype(int), 0, W.shape[1]-1); yi = np.clip(y.astype(int), 0, W.shape[0]-1)
    w = W[yi, xi]; q = Q[yi, xi]; ok = np.isfinite(w) & (q > 0.05)
    if ok.sum() < 30:
        print(f"{r['path']}: only {ok.sum()} usable samples"); continue
    w = w[ok]; med = np.median(w)
    s_switch = float(np.mean(np.abs(w - med) > 0.5)); s_mad = float(np.median(np.abs(w - med)))
    res.append({"path": r["path"], "wrap": int(r["wrap"][1:]), "verdict": r["gate_verdict"],
                "n": int(ok.sum()), "S_switch": s_switch, "S_mad": s_mad})
    print(f"{r['path']:34s} {r['gate_verdict']:8s} n={ok.sum():6d} S_switch={s_switch:.3f} S_mad={s_mad:.3f}")
json.dump(res, open(f"{OUT}/tripod_scores.json", "w"), indent=1)
pos = [x for x in res if x["verdict"] == "aprova"]; neg = [x for x in res if x["verdict"] == "reprova"]
def auc(p, n, key):  # holds-up should have LOWER S_switch -> AUC of (neg score > pos score)
    P = np.array([x[key] for x in p]); N = np.array([x[key] for x in n])
    return float(np.mean([(nn > pp) + 0.5 * (nn == pp) for pp in P for nn in N]))
rng = np.random.default_rng(0)
for key in ["S_switch", "S_mad"]:
    a = auc(pos, neg, key)
    boots = [auc([pos[i] for i in rng.integers(0, len(pos), len(pos))], [neg[i] for i in rng.integers(0, len(neg), len(neg))], key) for _ in range(2000)]
    p5, p95 = np.percentile(boots, [5, 95])
    print(f"PRIMARY {key}: AUC={a:.3f} (n={len(pos)} aprova vs {len(neg)} reprova), bootstrap 5-95%: {p5:.2f}-{p95:.2f} -> "
          + ("PASS" if (key == "S_switch" and a >= 0.75 and p5 > 0.5) else ("pass" if a >= 0.75 and p5 > 0.5 else "fail")))
rho, pv = stats.spearmanr([x["S_switch"] for x in res], [x["wrap"] for x in res])
print(f"SECONDARY Spearman(S_switch, wrap) over {len(res)}: rho={rho:+.3f} p={pv:.3f} (must be positive to reproduce his falloff)")
print("TRIPOD SCORING DONE")
