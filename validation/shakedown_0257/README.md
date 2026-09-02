# Shakedown on PHerc 0257, z 9000–10000 — artefacts

Preregistration: `../../preregistrations/2026-09-02_shakedown_pherc0257.md`.

## Positive control
`control_w035_prediction_vs_labels.png` — PHerc0139 w035 through the
identical pipeline; forward median on labelled ink 0.7294 vs background
0.2431, reverse none. Threshold for the target fixed at 0.7294.

## Flattening distortion map (item 9.5), scoped surface w010–w065
`distortion_map.py` computes, per flat cell, the ratio of its 3D footprint
area to its nominal flat area and the angle between its edges.
`w010-065_distortion.png`, `w010-065_distortion.json` (811,829 valid cells):

| | p5 | median | p95 |
|---|---|---|---|
| area ratio (1 = faithful) | 0.62 | **0.96** | 1.23 |
| stretch along the scroll | 0.71 | 1.01 | 1.32 |
| stretch along z | 0.74 | 0.97 | 1.15 |
| edge angle (90° = no shear) | 73° | 90.6° | 116° |

52 % of cells are within 10 % of true area, 81 % within 25 %; 74 % are
within 10° of unsheared. Distortion concentrates in the innermost wraps
(compressed and sheared) and in periodic diagonal streaks at the seams
between concatenated windings; the body of the surface is faithful. Any
letter called out from this render must quote the local area ratio, and
nothing at a seam streak counts (already excluded by the preregistration).
