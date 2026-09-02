# Integrator family (branch-cut, moving-front, min-cost-flow) — KILLED at the premise (2026-09-02)

**Question.** Would a localising integrator (Goldstein branch cuts, Costantini
min-cost flow, or a moving-front solve) beat the weighted least-squares
integration? Only if the constraint set is inconsistent in a localised way
that least squares smears. Criteria declared in `residue_premise.py`:
worth testing if ≥ 10 % of ridge pairs are violated by > 0.5 winding AND the
ten largest clusters hold ≥ 50 % of violations; kill if violations are rare
(< 3 %) or diffuse (top-10 < 20 %).

**Measured, PHerc 1667 z4634, 9.6 µm, pitch-locked field, 218,275 pairs.**

| quantity | value |
|---|---|
| pair residual, median | −0.04 winding |
| pairs violated by > 0.5 winding | **0.7 %** |
| pairs off by > 0.25 winding | 10.5 % |
| k-field plaquette residues (8 px loops) | 0.48 % of loops |
| clustering of violations | 923 clusters; top-10 hold **10.9 %** |

**Verdict: KILL (diffuse and rare).** The solve satisfies its constraints;
the field is wrong because the constraints are wrong (wrong pitch, wrong
pairs — the scale bias), not because of how they are integrated. No
integrator will fix a consistent set of wrong constraints. This also
answers the "do we need exotic optimisation for the integer part" question
in the negative from the other side: the integer part is not where the
error is.
