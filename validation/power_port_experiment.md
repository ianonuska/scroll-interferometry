# The power-port experiment (2026-08-31)

fit_spiral takes winding information through two very different inputs:

- the **annotation port** (`relative_windings.json` pcls) — sparse pearl
  strings, the port human annotations use. Our earlier benchmark measured
  this port's ceiling at ~3.5 points, and our field's constraints through it
  failed the held-out replication (README feedback-loop section).
- the **power port** (`winding_inference` store + `dense_spacing_mode=
  "winding_model"`) — dense ray crossings, normally produced by a trained
  neural winding model. This is what makes the human-free C configuration
  strong, and it only exists for scrolls that have a trained model.

`winding_to_store.py` encodes our field as that store. Question: can a
raw-CT physics field replace the trained model at the port that matters?

## Setup

- Held-out band (PHercParis4 z 11000–12000), 30k steps, no tracks, no
  manual pcl annotations — identical to the confC configuration except the
  `winding_inference` directory, swapped by symlink. Same dataset otherwise.
- Two stores tested: v1 (pre bug-9 fix; 6,904 rays / 15,962 crossings,
  ~2.6% of rays carried fabricated cross-hole deltas) and v2 (bug-9 fixed,
  audited zero non-unit rays; 3,450 rays / 7,929 crossings).
- The trained model's store: 1,460,171 rays / 26,743,553 crossings.

## Results (satisfaction on the fixed 89,237-patch denominator)

| run | winding source | patches | area |
|---|---|---|---|
| confA (reference) | manual annotations + trained model | 54.6 / 55.1 | 77.1 |
| confC | trained neural model | 53.6 | 76.8 |
| confC seed 2 | trained neural model (replicate) | 53.2 | 76.6 |
| confOurWM (v1) | our field, flawed store | 53.0 | 76.3 |
| confOurWM2 (v2) | our field, clean store | 53.1 | 76.0 |

The fitter's log confirms it consumed our stores (ray/crossing counts
match). Run-repeat noise on this band measured ±0.5–0.6 points.

## Reading

- With the trained model's own replicate in (53.6 / 53.2), the
  between-family gap (~0.35 points, means 53.4 vs 53.05) is smaller than
  the trained model's own seed-to-seed spread (0.4). Our field is
  **statistically indistinguishable from the trained model on this band**
  at the resolution these runs can measure. We still claim comparable, not
  better.
- The supervision budget difference is the point: **~3,400× fewer
  crossings** (7,929 vs 26.7M) for near-identical fit quality. The port
  works; the field's information survives the socket change.
- Contrast with the same field through the annotation port (confABf 50.0):
  the port, not the field, was the bottleneck. This closes the B-run
  question with a mechanism.
- The trained model exists only for this scroll. The 11 quiet Grand-Prize
  scrolls have neither annotations nor a model — there, this store is not
  an alternative but the only dense winding source on offer. That
  benchmark (quiet scroll, C-baseline vs our store) is the next
  experiment.
