# scroll-interferometry — the winding meter

A winding meter for Herculaneum scrolls: automated relative winding
annotations, made by treating the CT cross-section as an interferogram.
It measures a winding coordinate at each pixel. It is a measuring
instrument, not a fitter. Its output loads into the community's spiral
fitter as annotations.

Work-in-progress toward Vesuvius Challenge Open Problem #6 (winding-number
automation — "automating relative winding number annotations will boost
scalability by a great extent").

![showcase](figures/winding_showcase.png)

## The idea

A scroll cross-section is a fringe pattern: the rolled papyrus sheet is
quasi-periodic, therefore **one winding can be treated as one fringe**, and
*relative winding number* (the annotation the spiral-fitting pipeline needs and
currently gets from human annotators) **can be treated as fringe order**. This lets us use standard tools from RF/SAR interferometry:

1. **Structure tensor** → local sheet-normal orientation + coherence (a quality
   map for free).
2. **Multi-band steered Riesz quadrature** → local spatial frequency
   |k| = 2π / winding spacing. Per-pixel dominant-band selection avoids the
   low-frequency bias of any single band in densely packed zones.
3. **Ridge-adjacency constraints** — detect sheet ridges, march along the outward
   normal to the next ridge, emit a "+2π between these two points" pairwise term.
   *These are relative winding annotations, generated automatically.* They fix
   the classic isolated-fringe undercount: a lone sheet in a wide gap swings
   quadrature phase by π, not 2π; the discrete term restores the count.
4. **Weighted least-squares phase integration** (coarse-to-fine, Jacobi-PCG):
   `min Σ w|∇φ − k|² + Σ w_c(φ_b − φ_a − 2π)²`. The result φ/2π is a smooth,
   monotonic **winding coordinate**: iso-contours are individual windings, and
   the difference between any two points is their relative winding number.

This can be performed by streaming per-slice CT chunks directly from the
open data without extra GPU time, training data, or manual annotation. To me
it seems pretty lean.

## Day-one validation (PHerc 1667, mid-scroll slice, 19.2 µm)

PHerc 1667 was read end-to-end in June 2026, so ground truth exists. That is
why it was the dev target.

- **83,098** ridge-adjacency constraints generated automatically; solver fully
  converges on ~950k unknowns.
- Ray check (8 radial rays): winding span 29.3–31.9 vs. 33–43 naively counted
  sheet-crossing peaks. The field's ray-to-ray spread is **±4%**, versus ±15%
  for the raw peak counts it's checked against.
- External check: ~30–33 recovered windings at ~36 mm diameter imply
  ~1.3–1.5 m of rolled papyrus. The published unrolled length is 1.4 m.
  This came from a single slice.

## Multi-slice and second-scroll results (day one, continued)

- **Five slices across a 19 mm z-range of PHerc 1667** (z = 1800/2100/2317/2500/
  2800 at level 3): winding spans 33.0 / 35.3 / 33.9 / 37.4 / 35.8, every solve
  fully converged. Adjacent-slice fields (solved *independently*, 5.8 mm apart)
  agree after gauge alignment to a median residual of 1.4–2.7 windings
  (`validation/consistency.json`). These solves share no information, so this
  is an upper bound; the z-coupled solve below tightens it.

  ![multislice](figures/multislice_strip.png)

- Second scroll, no re-tuning: PHerc 172 (the Bodleian scroll; its title
  *On Vices* won the First Title Prize). Different scanner, resolution
  (7.9 µm), and energy (53 keV). The solver converges to a monotonic field of
  ~67 windings, core to rim. Contours are noisier than on PHerc 1667;
  scan-adaptive preprocessing is on the roadmap.

  ![second scroll](figures/second_scroll.png)

## Ground-truth benchmark against the accepted unwrap

PHerc 1667's released per-winding segments (`w011`–`w041`) are the accepted
solution of the fully-read scroll. We sample each segment's mesh near our five
slices and evaluate our field there (90 segment-slice samples across 20
segments, `validation/gt_rows.json`):

![gt](figures/gt_benchmark.png)

- Same-winding coherence: the field is near-constant within a ground-truth
  winding. Intra-segment MAD median 0.90 windings (p90 2.2).
- **Ordering:** winding order matches ground truth with median |residual|
  **0.25–1.2 windings** per slice against a linear map.
- Slope: the naive regression gives 0.65–0.77. About a third of the deficit
  traces to defective ground truth:
  one released segment labeled `w011` contains ~200k points spanning ~5 windings
  internally (an umbrella segment, not a single winding), and three 2024-era
  segments use an older segmentation's indexing. Restricting to the 2025–26
  single-winding segments: **slope +1.02 with median residual 0.20 windings at
  z=1800**, and +0.74–0.82 elsewhere (mean 0.81, residuals 0.13–0.53; one
  anomalous slice at z=2100 under investigation). Remaining caveat at 19.2 µm:
  the clean segments cover w11–12 and w28–41, leaving the dense middle band
  unmeasured — exactly where resolution blur would hide.
- Resolution: rerunning the same slice at
  **9.6 µm** resolves the dense middle windings the coarser level merged — the
  recovered span jumps 33.9 → 54.3, the naive slope moves **0.647 → 0.940**, and
  the clean regression (15 windings, middle band w13–23 now included) gives
  **slope +1.10, median residual 0.94 windings**. The compression at 19.2 µm was
  resolution blur, not an algorithmic defect. Operating envelope: run at
  ≥9.6 µm for production constraints (~6 min/slice on a laptop, no GPU);
  19.2 µm is a fast preview with a known ~0.8 gain.

Round-trip integrity: the exported JSONs load through the official
`spiral-fitting/point_collection.py` loader verbatim (24 collections / 783
points relative + 17 / 408 same-winding). The full spiral fitter requires a CUDA host; benchmark results below.

## Output: spiral-fitting constraint files (format-verified)

`winding_to_vc.py` exports the field as `vc_pointcollections_json_version: 1`
files matching the spiral-input dataset conventions (full-resolution `[z,y,x]`
voxel coordinates):

- `*_relative_windings.json` — radial "pearl strings", each point carrying
  `wind_a` = 0,1,2,… (relative-winding constraints)
- `*_same_windings.json` — iso-winding contour collections (same-winding
  constraints)

Historical note: the round-trip integrity check above (24 relative
collections / 783 points + 17 / 408) was performed with the original
ray-based exporter. The current exporter walks gradient streamlines instead
(the ray version silently broke on non-star-convex scrolls — bug #8) and
emits sparser, strictly monotone strings: e.g. ~40 collections / ~294 points
per PHercParis4 station. Both formats load through the official loader.

## Run it

### Docker (reproducible)

    docker build -t winding-meter .
    mkdir -p out && docker run --rm -v $PWD/out:/out winding-meter

Streams a PHerc1667 slice from the public S3 bucket and writes the winding
field + exported annotation JSONs — no local data needed.


```bash
pip install numpy scipy zarr s3fs matplotlib pillow
python examples/run_demo.py     # streams a PHerc 1667 slice, solves, exports
```

## Status & the benchmark results (updated Aug 31 — negatives included)

The fit_spiral benchmark ran on PHercParis4 (z 8500–9500, production config,
30k steps, scored by the fitter's own satisfaction metrics over the identical
89,237-patch verified set). Failures are listed with the successes. Run A was repeated to measure the
run-to-run noise floor: ±0.6 points on patches (52.0 vs 52.6).

| Run | Annotations | Satisfied patches | Satisfied area |
|---|---|---|---|
| A — manual | 2,173 human points | 52.0% | 75.9% |
| A — repeat (noise floor) | identical inputs | 52.6% | 76.0% |
| **A+Bf** — human + loop-vetted machine | +379 vetted machine points | **51.6%** | **75.3%** |
| C — none | 0 | 48.5% | 71.1% |
| B raw machine (best of 5 attempts) | 1,086 points | 44.2% | 66.7% |

**Confirmatory (held-out) band results — z 11000–12000, untouched by any
debugging decision, single pre-declared run per arm:**

| Held-out band | Satisfied patches | Satisfied area |
|---|---|---|
| confA — manual (run twice) | 54.6% / 55.1% | 77.1% / 77.1% |
| confC — no annotations | 53.6% | 76.8% |
| confA+Bf — manual + vetted machine | 50.0% | 77.1% |

Two results. (1) The redundancy finding replicated: manual annotations are
worth ~1–1.5 points here, against a ±0.5 noise floor. This confirms
pscamillo's result out-of-band with an independent method. (2) The
exploration-band parity claim failed replication: the vetted combination
lands below even the no-annotation control. That is the headline for that
arm. Annotations look like the wrong integration point for this instrument;
see the roadmap for the alternative (encoding the field as a
winding-inference store for scrolls that lack one).

What this shows and does not show: raw machine annotations as a drop-in
replacement hurt the fit (44.2% < C's 48.5%): the field under-counts windings
across compressed regions at this scan resolution (measured slope 0.647 vs.
manual over long spans; local correlation 0.93, n=15 matched pairs at one
station, so indicative rather than precise). A feedback loop that vets
each machine constraint against a fitted spiral and drops the untrustworthy
ones brought the human+machine combination to statistical parity with human
annotations alone, on the exploration band only. That parity did not survive
the pre-declared held-out replication (table above). No claim of parity or
improvement is made. Getting to this table took eight root-caused bugs (thread caps,
non-comparable losses, clobbered checkpoints, silent CG non-convergence, a
resolution-envelope violation, a stitched-volume coordinate frame,
non-star-convex ray sampling, gradient vortices at damage holes). Each fix is
in this repo's history with its diagnosis in the commit message.

Also measured along the way: under the current production config, the manual
annotations contribute only ~3.5 satisfaction points over no annotations at
all (A vs C). This is consistent with pscamillo's earlier, more thorough finding
(July 2026, "Winding evidence, measured": annotations statistically redundant
in patch-dense windows, worth +3–5pp in patch-sparse ones, 2 seeds per arm).
Ours is a single-band, single-seed corroboration, not a first. pscamillo also reported three constraint
generators degrading fits vs. a no-constraint baseline with root causes in
lasagna materialization resolution; our generator reads raw CT rather than
lasagna fields, so our B-run failures have distinct mechanics (documented
above), but the phenomenon has precedent. Next validation step before any
further claim: scoring this generator through pscamillo's constraint-gauge,
the community's shared referee for winding generators.

### Roadmap

1. Score the field through constraint-gauge (community referee).
2. Encode the field as a winding-inference store for scrolls without one.
3. Close the compression gap: resolution escalation in flagged regions
   (2.4 µm partial scans) + fit-feedback constraint refinement.
4. Patch auto-verification via winding-field consistency.
5. 3D-native solve.

## Ablations & extensions (running log)

- **Damage-aware ridge gate** (marching capped at 3.5 local wavelengths so a
  damaged gap can't alias into a skipped winding): validated on the clean GT
  set at 19.2 µm — slope 0.787 → **0.838**, intra-winding MAD 1.19 → **0.97**,
  max residual unchanged. Now the default.
- **Z-coupled slab solve** (`winding_slab.py`): joint least squares over a
  5-slice stack (115 µm spacing) with soft z-continuity, warm-started from
  gauge-aligned per-slice solutions. Result: adjacent-slice disagreement drops
  from median 0.36–0.50 windings (independent solves) to **0.018–0.026** —
  ~20× — with the winding span preserved (35.3 → 34.6, no oversmoothing), at
  **zero ground-truth cost** (slab center slice: slope +0.833, residual median
  0.29, intra-MAD 0.92 — matching the best per-slice result). The joint solve
  of all 5 slices (264 s) is faster than five independents.

  ![slab](figures/slab_result.png)

## Generalization: zero-shot on four untouched Grand-Prize volumes

`scouting/letters_scout.py` streams a mid-volume slice from four 2027
Grand-Prize-eligible volumes with essentially no community segmentation
activity (PHerc1203, PHerc1218, PHerc0125, PHerc0257 — 8.6–9.4 µm masked
scans) and runs the winding solver with zero per-scroll tuning. All four
converge (CG info=0) with winding spans 37–57 at preview resolution
(level 2, ~35 µm — expect undercounting from blur per the resolution
study above; this is a scan-character smoke test, not a count claim).

![four-volume scout](figures/letters_scout_4vol.png)

## How this differs from neighboring tools — stated plainly

Several community tools sound similar. The differences, one line each:

| Tool | What it does | What ours does differently |
|---|---|---|
| **winding-sync** (Balmaceda) | Integer winding labels from structure-tensor orientation, solved discretely | Ours also measures winding **spacing** (a frequency channel), and solves one continuous convex system — a dense coordinate at every voxel, not labels on patches |
| **winding-ruler / constraint-gauge** (pscamillo) | Measures winding pitch along rays; scores any generator against human GT | A ruler and a referee — not a field. We are a **generator**, and we will be scored by their referee |
| **Prior constraint generators** (three, per pscamillo) | Generated constraints from **lasagna fields** (ML predictions); all three degraded fits; root cause = lasagna resolution | Ours reads **raw CT** — no ML predictions anywhere upstream, so lasagna's resolution ceiling does not bind us. Our failures had different causes (documented above) |
| **scroll-truth** (karasukun) | Answers "same wrap or not?" for a PAIR of patches, from raw CT | Same philosophy (raw CT, judge models with model-free references) — but pairwise verdicts, not a global coordinate. Complementary: our field assigns the wrap number both patches share |
| **ScrollAnchor** (olgaiv39) | Per-vertex drift/sheet-switch diagnostics for an existing tifxyz surface | Diagnoses surfaces after they exist; our field exists **before** any surface and can position disconnected pieces globally |
| **ScrollFiesta** (Kyles) / **scrollreading** (Stevens) | Build and assemble surface meshes from ML predictions | Sheet **builders**. Ours is a **map**: the global winding address their assembled pieces need for placement |

Summary: the other tools label, build, check, or score pieces. This one
measures a continuous winding coordinate for every voxel directly from raw
CT. Whether that map is accurate enough to matter is an empirical question,
scored above and next against constraint-gauge.

## Relation to existing community work

Two community projects attack the same problem with different machinery, and this
work should be read alongside them:

- **winding-sync** (Joseph Balmaceda) derives relative winding constraints from
  structure-tensor lamina *orientation* and reconciles contradictions globally as
  L1 integer synchronization.
- **winding-ruler** (pscamillo) measures winding-pitch evidence across scrolls.

What this repo adds beyond orientation-only approaches: a frequency channel
(steered quadrature measures local winding spacing, not just direction); a
continuous dense winding field from a convex least-squares solve, with no
integer programming; quality maps (coherence times fringe amplitude) marking
where the answer is trustworthy; and a cheap external check (winding count
times mean circumference vs. known unrolled length). Whether continuous-field
or integer-sync works better inside fit_spiral is untested. Combining them
(this field as winding-sync's prior) might beat both.

Feedback, failure cases, and prior art welcome, especially from anyone who
has tried a frequency-domain framing before.

## Data & license

Scan data: [Vesuvius Challenge open data](https://scrollprize.org/data)
(`s3://vesuvius-challenge-open-data/`, CC BY-NC 4.0) — streamed, never bulk
downloaded. Code: MIT. Built by Ian Onuska (ONUSKA Industries) with Claude.
