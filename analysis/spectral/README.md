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
