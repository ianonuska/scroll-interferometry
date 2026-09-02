# Ambiguity map prototype — FAILED its pass criterion (2026-09-02)

**Idea.** Re-solve the winding field with random 30% subsets of the
ridge-adjacency constraints dropped (B = 8 bootstraps) and take the per-node
standard deviation as an "ambiguity" score, on the theory that places where
the constraint set is fragile are the places the field is wrong.

**Declared criteria before running.** Kill if ambiguity is ~perfectly
correlated with the existing quality weight (then it is nothing new). Pass if
ambiguity predicts within-wrap error (accepted-segment residual) with
AUC ≥ 0.7.

**Result (PHerc 1667, z4634, 9.6 µm, 218,275 pairs).**
- corr(ambiguity, quality) = +0.087 → not killed; it is independent of quality.
- AUC(ambiguity | high-error vs low-error nodes) = **0.489** → chance. FAIL.
- ambiguity percentiles: p50 0.05, p90 0.08, p99 0.12 windings.

**Reading.** The field barely moves under random constraint loss, and where it
moves has nothing to do with where it is wrong. The errors in this field are
systematic (the fringe-scale defect documented in ../../validation/fringe_scale),
not variance from a fragile constraint set, so a variance-type ambiguity map
cannot see them. A different design — enumerating alternative *pairings*
rather than dropping constraints at random — is a separate hypothesis and is
not claimed here.

Script: `ambiguity_L2.py` (run as committed). Raw tail of the log:
`ambiguity_L2_results.txt`.
