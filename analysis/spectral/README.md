# Super-resolution spectral estimation (item 10.A.1) — synthetic premise: not worth building (2026-09-02)

Aim: separate two sheets merged into one fringe (the 19.2 µm undercount)
with a model-order method (MUSIC), and measure the invented-component rate
first. `music_premise.py` (criteria in the header): 64-sample profiles, true
pitch 12 px, unit sinusoids, white noise at 10 and 20 dB, and a 5 % pitch
drift across the window as a real scroll has. Forward–backward spatially
smoothed MUSIC (snapshot 32, order 2) vs a Hann-windowed 8× zero-padded FFT
peak pick; one-to-one matching within 8 %.

| spacing ratio | MUSIC resolves both (10 / 20 dB) | FFT (10 / 20 dB) |
|---|---|---|
| 1.10 | 0.00 / 0.00 | 0.01 / 0.01 |
| 1.25 | **0.00 / 0.00** | 0.12 / 0.17 |
| 1.50 | 0.18 / 0.26 | 0.66 / 0.65 |
| 2.00 | 0.89 / 0.99 | 1.00 / 1.00 |
| one sheet, order 2 → invented second sheet | **0.00 / 0.00** | — |

Bar: worth building if ratio 1.25 ≥ 0.80 and invented rate ≤ 0.05.
**Not worth building.** It does not hallucinate sheets here, but it resolves
nothing the FFT cannot, and less: a 5 % chirp breaks the stationary-sinusoid
assumption MUSIC rests on. A chirp-compensated variant would be a different
item and is not built. First harness attempt had two bugs (snapshot shorter
than the longer period; a scorer that let one peak satisfy both truths);
both fixed before the numbers above, criteria untouched.

## Cepstral pitch (10.A.2) and zero-crossing period (10.B) — both KILLED on real data
`cepstrum_premise.py`, `cepstrum_premise.json`. PHerc 1667 z4634, 9.6 µm,
3000 points on accepted wraps with known true pitch; 96-px profiles along
the sheet normal. Bar (declared): worth testing further if correct ≥ 0.54
(pitch lock + 0.10); kill if below the default dominant band (0.37).

| estimator | correct (0.7–1.4× true) | at half the true pitch |
|---|---|---|
| cepstrum, largest quefrency peak in 8–60 px | **0.359** | 0.570 |
| zero-crossing period | **0.029** | — |
| (reference) dominant band / pitch lock / oracle | 0.369 / 0.440 / 0.778 | |

Both killed. The cepstrum locks to the face spacing (half the winding
pitch) on 57 % of points — the same harmonic the dominant band locks to. At
this resolution the strongest periodicity along the normal *is* the
face-to-face spacing, so any estimator of "the period" inherits the
overcount. This closes the estimator family (autocorrelation, quadrature,
cepstrum, zero-crossing, MUSIC) as a route to the true pitch; the
information is in which band to trust, not in a better period estimate.
