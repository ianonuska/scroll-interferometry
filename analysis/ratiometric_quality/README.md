# Ratiometric (locally-normalised) quality weight — premise test FAILED (2026-09-02)

**Idea (from the normalising-transimpedance-amplifier pattern):** the current
quality weight is `coherence × clip(amp / global_p95)`, i.e. fringe amplitude
normalised by one global number, so the weight follows local brightness
rather than local fringe contrast. Candidate: normalise amplitude by the local
RMS of the high-passed image over one fringe wavelength (dimensionless
contrast). Declared before running: KILL if it is coherence re-expressed
(|corr| > 0.95); worth the full solve + box M3 gauge only if it predicts
within-wrap error better than the current weight by ≥ 0.05 AUC and ≥ 0.60.

**Result (PHerc 1667 z4634, 9.6 µm, locked field, 60,368 labelled nodes on
accepted wraps; script `ratiometric_premise.py`):**

| quality map | corr with current w | corr with local brightness | AUC(low quality ⇒ high error) |
|---|---|---|---|
| current w | 1 | +0.42 | 0.541 |
| ratiometric | +0.97 | +0.29 | 0.534 |
| coherence alone | +0.65 | — | 0.530 |

Fraction of high-error nodes by quality decile (low → high), current weight:
0.80 0.79 0.78 0.77 0.75 0.76 0.77 0.77 0.76 0.73 — flat. Ratiometric: the
same to two decimals.

**Reading.** The brightness dependence drops as predicted, but the candidate
is the current weight in different clothes (0.97) and neither map knows where
the field is wrong. The deeper fact is in the decile row: ~77 % of nodes on
accepted wraps sit more than half a winding from their own wrap's median in
this field. That is the fringe-scale bias (`../../validation/fringe_scale`),
not noise, and no confidence channel can be calibrated against an error that
is a bias shared by every node. Consistent with `../ambiguity` (variance-type
map also blind). Order of operations, then: fix the scale first; re-ask about
calibration (gauge M3) only afterwards. Not built.
