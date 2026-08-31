# Cross-check on pscamillo's constraint-gauge (2026-08-31)

The community constraint-gauge bench scores winding generators per point:
exact agreement on winding-difference-1 pairs (M1), mean absolute error of
the winding difference (M2), coverage (M4), and whether declared confidence
predicts accuracy (M3). We ran our field through it unmodified.

## Setup

- Scorer: `constraint-gauge` at commit `a72c423`, code untouched.
- Ground-truth arm: the 19 accepted single-winding PHerc 1667 segments
  (w011–w041), loaded through the gauge's own tifxyz loader, each vertex
  labeled with its accepted winding number. The merged whole-scroll surface
  (`20260612...merged_v4`) is excluded — it spans the full radius
  (968–6067 voxels at our test slice) and is not a single winding. This is
  the same umbrella segment already excluded from the clean benchmark above.
- Adapter: our z=2317 field (19.2 µm working resolution), 36,502 nodes,
  96 µm node spacing, confidence = the solver's quality channel.
- Same volume frame on both sides (`20251217075048-2.399um`), so no
  cross-frame conversion. Sign checked by correlation (+0.51 raw), not flipped.
- Caveat: this arm is ours (built from public accepted segments), not one of
  the gauge author's arms, so comparisons to his published lines indicate
  scale only.

## Result (`constraint_gauge_summary.json`)

- Density gate: node gap 40 vox vs sheet gap 93 vox — SCORABLE.
- **M1 exact-dw1: 0.190.  M2 MAE: 2.47 windings.  M4 coverage: 0.352.**
- **M3: our confidence channel is uninformative** — accuracy is flat
  (~0.13) across all ten confidence deciles, ECE-rank 0.39. The quality
  channel gates where we emit points; within emitted regions it does not
  rank winding accuracy. Treat exported confidences as coverage flags, not
  probabilities.

## Reading

The two benchmarks measure different things and both are true:

- Medians per segment (the benchmark above): corr 0.93–0.99, slope up to
  +1.10 at 9.6 µm. The field is *systematically* right — the winding
  ordering and scale are correct.
- Per-point pairs (this gauge): MAE 2.47 windings at 19.2 µm. Any single
  point can sit ±2 windings off. The field is a coordinate, not yet a
  per-sheet assignment, at this working resolution.

For scale: on the gauge author's own (harder, different-arm) lines,
winding-sync solvers score M1 0.017–0.050 / M2 21–25, and his tuned
estimators M1 0.23–0.30 / M2 2.6–3.5.

## What this changes

Per-point exactness at higher resolution (9.6 µm halved the median-level
error) is now a measured, prioritized gap, and the confidence channel needs
an accuracy-linked definition before anyone should weight by it.
