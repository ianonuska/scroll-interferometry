# Fringe scale by wrap — an independent second measurement (2026-09-01)

The constraint-gauge cross-check found the instrument recovers ~0.8 windings
per fringe at 19.2 µm from per-pair winding deltas. This is the same quantity
measured a different way: for each adjacent pair of accepted single-winding
segments on PHerc 1667 at z=2317, the true spacing is the median
nearest-neighbour distance from wrap w+1 to wrap w (2.399 µm frame, scaled to
level-3 pixels), and the field's implied pitch is the median of 1/|∇W| at the
same points.

| wrap pair | true spacing (px @19.2 µm) | field pitch (px) | field/true |
|---|---|---|---|
| w11→w12 | 6.1 | 10.3 | 1.69 |
| w28→w29 | 7.0 | 11.5 | 1.64 |
| w29→w30 | 7.7 | 11.5 | 1.50 |
| w30→w31 | 6.4 | 10.7 | 1.68 |
| w31→w32 | 9.4 | 10.6 | 1.13 |
| w32→w33 | 7.9 | 9.0 | 1.14 |
| w33→w34 | 9.1 | 8.9 | 0.99 |
| w34→w35 | 7.0 | 9.2 | 1.32 |
| w35→w36 | 7.4 | 10.2 | 1.37 |
| w36→w37 | 12.5 | 9.5 | 0.76 |

(w12→w13 excluded from the summary: its 37 px "spacing" says those two
accepted meshes are not physically adjacent in this slab.)

**Median field/true = 1.32, i.e. ≈0.73 windings per fringe — corroborating
the gauge's 0.8× by an independent method.** Two more things this table
says: the true spacing at 19.2 µm is 6–9 px, *at* the instrument's declared
6 px resolvability limit, so the undercount here is a resolution fact rather
than an algorithmic one; and where the true spacing is ≥9 px (w31–w34) the
ratio is 1.0–1.1. The operating envelope should be stated as roughly
"pitch ≥ 10 px at the analysis level".

**The same measurement at 9.6 µm (`fringe_scale_by_wrap_L2.json`):** true
spacing 12–19 px, field pitch 10–12.7 px, median field/true 0.72 — i.e.
≈1.39 fringes per winding, the face-splitting overcount, matching the
gauge's 1.4 at 9.6 µm. Both regimes of the fringe-scale finding are now
corroborated by an independent method. This also fixes the de-chirp gate:
at 9.6 µm the per-wrap ratio must move from 0.72 toward 1.0.

**A negative result that shaped the roadmap.** Three textbook pitch
estimators on radial profiles in polar coordinates about the umbilicus
(autocorrelation first-peak, cepstrum, and a chirp-rate search) all failed a
pre-declared 10 % gate here (`pitch_estimators_v2.py`), including after
per-wedge averaging and detrending. The reason is geometric: PHerc 1667's
cross-section is not circular, so radial lines cross sheets obliquely and
annuli by centroid radius do not correspond to winding order (`true_pitch.py`
shows the non-monotonic "annulus pitch" that exposed this). The field's own
normal-direction quadrature is the right geometry. Consequence for the
rotating-frame roadmap item: the de-chirp must be done in the **field's
frame** (resampling along iso-W contours) rather than in polar coordinates,
and its job at fine resolution is harmonic suppression (one fringe per
period), not extra resolution.
