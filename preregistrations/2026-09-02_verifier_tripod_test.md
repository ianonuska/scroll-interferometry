# PREREGISTRATION: the verifier's external test ("tripod tablet")
# Committed before any score is computed. 2026-09-02.

## Claim under test
Our winding field can grade a fitted surface — i.e. tell a surface that
"holds up" from one that has wandered off its sheet — on a quiet scroll,
with no reference other than raw CT.

## External labels (not ours)
pscamillo's `vesuvius-eligible-meshes` (published 2026-09-01): 84 spiral-
fitter meshes inspected by eye against a fiber-weave reference, with a
per-mesh verdict `aprova` / `parcial` / `reprova`. Pilot set: **PHerc0800**,
23 inspected meshes over 17 z windows (14 aprova, 4 parcial, 5 reprova).
Replication set, if the pilot passes: PHerc0813 (24), 0211 (25), 0125 (12).

## Our measurement (fixed now)
For each inspected mesh: solve our field at the mesh window's central z
(window_z + 400) at level 1 (17.3 um) with the pitch lock (prior from a
level-2 solve, lock_tol 0.6, ridge_min_frac 0.55), cored on the community
umbilicus. Sample the mesh's vertices within +-24 voxels of the station z,
read W at (x, y)/2. Per mesh, two scores:
  S_switch = fraction of sampled vertices whose W deviates from the mesh's
             median W by more than 0.5 windings;
  S_mad    = median absolute deviation of W across the mesh's samples.
A single-winding surface that holds up should be near-constant in W.

## Statistics and criteria (declared)
- Primary: AUC of S_switch for aprova (positive = holds up) vs reprova,
  parcial excluded. With n = 14 vs 5 the CI is wide; report a bootstrap CI.
  **Pass:** AUC >= 0.75 and the bootstrap 5th percentile > 0.5.
- Secondary: Spearman correlation of S_switch with wrap number across all 23,
  which must be positive if we reproduce his "quality falls with distance
  from the umbilicus" (83 % usable at w020 -> 33 % at w100).
- Failure is reported as failure. A pass on the pilot is a *pilot* result;
  the replication set decides.

## What would embarrass us and is therefore checked
- Frame: mesh coordinates are in the 8.64 um volume's native frame; our
  stations are solved on the same volume and scaled by 2. Verified by
  overlaying mesh vertices on the station slice before scoring.
- Circularity: the labels were made by eye against fiber weave, not with any
  winding field; our field never sees the labels until scoring.
- The labels are one person's judgment (his README says so).
