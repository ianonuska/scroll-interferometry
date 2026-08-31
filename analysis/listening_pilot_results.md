# Listening pilot: results (2026-08-31)

Question (preregistered before any data was touched —
`preregistrations/2026-08-31_listening_pilot.md`, commit 3c2199d; amendment
1 committed before analysis): does carbon ink leave a measurable texture
fingerprint in standard absorption CT, where it has ~no density contrast?

## The short answer

Yes, on 4 of 6 fragments, with modest effect sizes — and the result
survives the brightness confound on the two fragments where that check has
statistical power. Two fragments are honest nulls. One declared analysis
choice (the Otsu ink mask) failed on real data and is fully disclosed below.

## What happened, in order

1. **Preregistered run** (Otsu-on-IR masks, exactly as declared):
   verdict FINGERPRINT — S1 variance significant, same sign, in 4/6
   fragments. BUT our own mask diagnostic then showed the declared Otsu
   thresholding is invalid on several fragments: it classified 92%/95% of
   Frag1/Frag3 as "ink" (it splits illumination, not ink). The
   preregistered verdict therefore cannot be interpreted as an ink effect,
   and we do not claim it. (`listening_pilot_results.csv`)
2. **Sensitivity rerun** (post-hoc, labeled: provider-aligned
   `inklabels.png` as the mask for all six fragments; ink fractions a
   plausible 9–18% everywhere): **verdict FINGERPRINT again** — S1
   significant and positive in the same 4/6 fragments (Frag1 +0.22,
   Frag2 +0.14, Frag3 +0.16, Frag6 +0.16, Bonferroni-corrected p ≤ 2e-6;
   Frag4 and Frag5 null). S3 (lag-1 autocorrelation, positive) and S4
   (high-frequency fraction, negative) co-move in 3 of those 4.
   (`listening_sensitivity_inklabels.csv`)
3. **Confound analysis** (promised in the prereg): controls re-matched on
   patch mean intensity, one control per ink patch. On the two fragments
   with plausible masks AND large control pools, the signal survives:
   Frag6 all four statistics (|delta| 0.18–0.24, p ≤ 6e-26 after
   matching), Frag2 three of four. On Frag4/Frag5 the (already-null)
   effects die, consistently. Frag1/Frag3 could not be confound-tested
   under the preregistered masks (79/33 controls); re-testing them under
   inklabels masks is open work.

## Reading

The signature — MORE variance, MORE short-range correlation, LESS
high-frequency power in ink regions — is what you would expect if ink
fills and coats papyrus fibers: fine fiber texture damped, larger-scale
structure added. It is weak (Cliff's delta ~0.15–0.25), far from a
per-pixel ink detector, and absent on 2 of 6 fragments (different scrolls
and scan energies; the prereg's no-pooling rule applies). It is, however,
a training-free signal in a channel the community's ink detection does not
currently use, and it is the empirical motivation for the dark-field
acquisition proposal: if this much microstructure signal survives
absorption reconstruction, a scatter-sensitive modality should see far
more.

## Files

- `listening_pilot.py` — the analysis (preregistered mode and
  `--mask-source inklabels` sensitivity mode)
- `listening_confounds.py` — the confound table generator
- `listening_pilot_results.csv`, `listening_sensitivity_inklabels.csv`
- Raw logs on request; every number regenerates with one command from
  public data.
