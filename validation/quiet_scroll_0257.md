# The market test: PHerc 0257, a scroll with no annotations and no trained
# winding model (2026-08-31)

The question the whole quiet-scroll strategy rests on: on the terrain where
the prize scrolls live — no manual annotations, no trained winding model —
does our winding store help the community's fitter?

## Setup

- Scroll: PHerc 0257 (volume `20250821151750-9.362um`), chosen because it has
  a public tracks dataset, a hash-verified community umbilicus, and our field
  had already converged on it during the letters scout.
- Baseline: the published quiet-scroll recipe from the first-light PHerc 0826
  team — patches disabled, dense-spacing weight zero, tracks only — on
  z 9000–10000, 30,000 steps. In that configuration the fitter's
  winding-model input (the "power port") is empty, because nothing exists to
  fill it.
- Metric: the fitter's own `satisfied_track_points_fraction` over 34.4
  million track points. One variable per arm.
- Winding sense decided empirically by a 1,500-step CW/ACW pair, as the
  recipe prescribes (CW won by 0.25 points — a thin margin, noted).
- Our products: five station fields cored on the community umbilicus, an
  outer-shell tifxyz built from the field's masks, and a winding store
  encoded with the bug-9-fixed encoder (audited: zero non-unit rays).

## Results

| arm | winding source | satisfied track points |
|---|---|---|
| baseline | tracks only (power port empty) | **57.88 %** |
| shell-only | + our outer shell, no store | 57.51 % |
| ours, stations solved at 37.4 µm | + our winding store | **55.75 %** |
| ours, stations re-solved at 18.7 µm | + our winding store | **57.53 %** |

**The first result was negative: our store made the fit worse by 2.1
points.** The fitter's log confirms the store was consumed (41,598 rays). The
baseline sanity-checks against the 0826 team's published numbers (16.2 %
fully-satisfied tracks here versus their 14.4 % on their scroll).

## Why — measured, not asserted

The suspected mechanism was our own fringe-scale finding
(`constraint_gauge_cross_check.md`): stations solved at 37.4 µm have a
measured pitch of 10.3 px, inside the regime where the instrument merges
windings. So the store asserted "one winding" where the truth was about
three, at a trust weight calibrated for a much denser, more accurate
supplier.

To test that, the agreement between our field and the fitter's own 120
fitted wraps was measured directly — five stations, ~235,000 sampled points,
the same protocol as the 1667 ground-truth benchmark
(`quiet_scroll_0257_agreement_L2.json`, `..._L1.json`):

| field resolution | correlation with the fit | slope | median residual |
|---|---|---|---|
| 37.4 µm | +0.990 to +0.999 | **0.30 – 0.39** | 0.26 – 0.98 windings |
| 18.7 µm | +0.994 to +0.999 | **0.65 – 0.79** | 0.54 – 1.76 windings |

Two predictions were written down before the corresponding measurement, and
both held: that re-solving at twice the resolution would roughly double the
slope (it did), and that the refit would recover *part but not all* of the
lost 2.1 points (1.8 of 2.1).

## Reading

- At usable resolution our store is **net-neutral** on a tracks-rich quiet
  scroll: it stopped hurting, and it adds nothing the 500 million track
  points do not already supply. That is the honest state of the "become the
  winding model for quiet scrolls" thesis: not yet a contribution to surface
  quality.
- The ordering agreement (correlation +0.99 at both resolutions) is
  excellent; the *scale* is not yet right, and the ladder says why —
  slope 0.33 harms, 0.70 is neutral, and ~1.0 has not been tested because
  it needs full-resolution stations or the de-chirp module on the roadmap.
- What the field can already do here is **verification**: a public-tree
  search found no published global winding reference of any kind for
  PHerc 0257, and a field that agrees with the fit's ordering at 0.99 is an
  independent check the scroll otherwise lacks. Ordering-verifier now;
  count-verifier only after envelope-compliant solves.

## Two things tried and stopped

**A track filter, killed by its own pre-measurement.** Premise: does our
field disagree with individual tracks often enough to be worth filtering
them? Median disagreement along a track is 0.11–0.14 windings — an
incidental mutual validation of both structures. Only 3–6 % of tracks are
flagged, and flagged tracks are *straighter* than unflagged ones by the
track store's own tortuosity metric (1.08 vs 1.18): a length confound, not a
defect detector (`quiet_scroll_0257_track_flags.json`). Not built.

**A warm-start experiment, resolved by its control (2026-09-01).** Letting
our store shape the first 3,000 steps and then handing over to the pure
baseline for 27,000 more scored 53.6 % against the baseline's 57.9 %. The
control — baseline for 3,000 steps, then the identical resume — scored
**53.6 % as well** (0.5359 vs 0.5362). The entire 4.3-point drop is the
resume mechanics (learning-rate schedule and optimizer state), not the
store-shaped start. So the warm-start question is *unmeasured*, not
negative: a fair test needs a single uninterrupted run whose store weight is
switched off mid-run, which the fitter's run-mutable weights appear to
allow. Both numbers are kept so neither can be mistaken for a result.

## Provenance notes

- Two chirality prototypes (handedness from the polar shift field, then from
  a polar-FFT asymmetry) were built and both failed their own pre-declared
  gates on 1667; no prediction was issued. An earlier claim that the field
  could read handedness directly was wrong — the onion-form coordinate is
  achiral by construction — and was retracted the same hour.
- Every number above regenerates from the public data with the scripts in
  this repository plus the community's `fit_spiral` at the pinned villa
  commit.
