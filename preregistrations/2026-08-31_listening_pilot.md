# PREREGISTRATION (DRAFT for Ian's review & rewording): the "listening" pilot
# Does carbon ink leave a microstructure fingerprint in standard CT texture?
# Status: PREREGISTERED — criteria approved by Ian Onuska 2026-08-31 (all five
# decision points approved verbatim in review). No analysis has been run.
# is touched.

## Hypothesis
Carbon ink differs from papyrus in microstructure (particle granularity),
not density. Reconstruction averages away most scattering information, but
some may survive as local texture statistics. H1: texture statistics at
ink locations differ measurably from no-ink locations in standard scans.
H0: no measurable difference (a publishable null that motivates dark-field
acquisition, where the scattering signal is captured directly).

## Data (all public, none of it ours)
- Fragments with infrared ground truth (dl.ash2txt.org/fragments/): the IR
  photo proves where ink is on the exposed surface; the aligned surface
  volume gives the CT texture at those exact locations.
- Planned set: Frag1-Frag6 (each fragment = one independent replicate).
- Regions: ink mask from IR threshold; no-ink control regions sampled from
  the same fragment at matched depth/distance-from-edge (to control for
  surface position confounds).

## Statistics (declared now, before looking)
Per 64x64-px patch, computed on the middle surface-volume layer ±2:
  S1 local variance; S2 excess kurtosis; S3 lag-1 autocorrelation (x and y);
  S4 high-frequency band energy fraction (above 1/4 Nyquist).
Test: per-fragment Mann-Whitney U between ink and control patch
populations, per statistic; significance threshold p < 0.01 AFTER
Bonferroni correction (4 statistics x 6 fragments = 24 tests).

## Success / failure criteria (committed in advance)
- FINGERPRINT: >= 1 statistic significant in >= 4 of 6 fragments with the
  SAME SIGN of effect in all significant cases. (Same-sign requirement
  prevents cherry-picking mixed noise.)
- NULL: anything less. Reported as the headline either way.
- Effect size reported regardless (Cliff's delta), with per-fragment CSV.

## Known confounds we will check and report
- Ink regions may correlate with surface position/curvature (IR photographs
  the exposed side): control patches matched on depth-layer and local
  surface tilt; a confound analysis table ships with the result.
- Scan-to-scan differences across fragments: no pooling across fragments;
  per-fragment tests only, aggregated by vote as above.

## Outcomes and what each means
- Fingerprint -> a training-free "listening" ink channel exists at reduced
  strength in current data (candidate detector #4), and dark-field
  acquisition is motivated as its full-strength version.
- Null -> the microstructure signal provably does not survive standard
  absorption reconstruction; this is the empirical case FOR a dark-field
  pilot scan (e.g., of the campfire scrolls) and gets written up as such.

## Provenance
Idea: Ian Onuska (2026-08-31), from RF channel-sounding/passive-vs-active
framing. Analysis code will be MIT, one-command, CPU-only. AI-assisted
implementation (Claude), disclosed. Prior art checked: dual-energy
absorption subtraction (PaulG) = different channel, null; no community
proposal of scatter/texture-statistics ink detection found as of this date.

## Amendment 1 (2026-08-31, before any real-fragment analysis)

Declared after inspecting the public directory layouts and BEFORE running
the analysis on any fragment. Data-availability accommodations only; no
statistic, threshold, or criterion changes:

- Frag4 ships no aligned ir.png. Its aligned inklabels.png (derived from
  the IR photograph by the data providers) is the ink mask for Frag4 only.
  All other fragments use the declared Otsu threshold on ir.png.
- Frag5 is scanned at 70 keV / 3.24 um and Frag6 at multiple energies; we
  use Frag6's 53 keV / 3.24 um series (closest to the scrolls' 54 keV).
  Per-fragment scan differences were already handled by the declared
  no-pooling rule.
- "Middle surface-volume layer +-2" is computed from each fragment's actual
  layer count (middle = n//2), not assumed to be layer 32 of 65.
- File layouts differ per fragment; the exact staged paths are recorded in
  the staging script committed alongside this amendment.
