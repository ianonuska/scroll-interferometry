# Preregistration — end-to-end shakedown on PHerc 0257, z 9000–10000

Committed 2026-09-02, before any ink inference on this window. The criteria
below are adopted unchanged from the first-light-pherc0826 team's readout
(their prereg/readout.md); they are not ours to soften, and they are not
edited after inference. Sections marked "filled later" are additions, not
edits.

## Framing
This is a **shakedown of the reading chain**, not a letters attempt. Success
means learning where the chain breaks, with letters as upside. Any rendered
image published from this run carries its scale, the ML patch size, and the
label "control" where applicable, and claims nothing beyond this file.

## Input
- Spiral fit: the organizers' `fit_spiral` **baseline** run on PHerc 0257,
  window z 9000–10000, 30,000 iterations, satisfied-track-points fraction
  0.57875 — the arm with no winding field, so nothing from this repository
  is upstream of the render. A chain break here is the chain's.
  Fitted meshes: 120 wraps w010–w129 (`w*_spliced`), coverage proportional
  to circumference, no thinning.
- Volume: `s3://vesuvius-challenge-open-data/PHerc0257/volumes/20250821151750-9.362um-1.2m-113keV-masked.zarr/`
- Ink model: `scrollprize/ink_9um`, `hybrid_3d2d-seed42/step-075000.pth`,
  `--direction both`, hann blend, overlap 0.5 — the first-light recipe
  unchanged.

## Scope — declared before any target render
The full range w010–w129 was flattened first (one 49,794 × 204-cell surface,
3,055,435 valid cells, z 9000–10000 confirmed): 1.22 gigapixels per rendered
slice at the model's resolution, beyond available storage for the 28-slice
recipe. The first pass is therefore **wraps w010–w065**, the inner range,
the same scope the first-light team used. Wraps w066–w129 are **not
evaluated** in this pass — not "no ink." A second pass on the outer range is
a separate run under these same criteria.

## What counts as a positive signal
- Ink probability at or above the control-derived threshold: the median
  probability the identical pipeline assigns to the control segment's
  labeled ink pixels, computed before target inference and recorded in the
  Calibration section at that time; never adjusted afterwards.
- Spatially coherent, contiguous, stroke-like (elongated) region at least
  0.5 mm (~53 px at 9.362 µm) along its long axis.
- Aligned with the fibre direction visible in the same render.
- Present in **both** `--direction` outputs.

## What does not count
Isolated pixels or blobs; anything at tile-blend seams or winding
boundaries; signal in only one direction output; anything failing the
control comparison; anything outside the declared scope.

## Held-out check / positive control
Control: PHerc0139 segment w035 with its published ink labels
(`scrollprize/datasets`, `ink/0139/w035_2026031718`), run first through the
identical pipeline. That segment's coordinates are in the 2.399 µm volume
`20260102150214-2.399um-0.2m-78keV-masked`; the ink_9um model card states
its PHerc0139 training segments were rendered from that volume downsampled
4× (pyramid level 2, 9.596 µm isotropic), so the control is rendered that
way — the model's own training route. Rendered against the 9.362 µm volume
the same segment produces an empty output, which is recorded here as a
finding. Labels (22,640 × 20,000) are nearest-resampled to the level-2
render grid for the calibration medians. Caveat stated up front: w035 is in
the model's training set; this validates the pipeline, not generalisation.
If the control does not show its letters with clearly separated
ink/background medians, the pipeline is declared uncalibrated and no
statement of any kind is made about the target.

## Calibration — filled from the control run, before target inference
(empty until then)

## False-positive statement — written after inference; criteria above never edited
(empty)

## Verdict — one sentence, human-written, after the statement above
(empty)
