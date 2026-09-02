# Twist-rate z-propagation (item 9.7) — premise not viable as designed; shelved (2026-09-02)

Idea: estimate one scalar (rotation of the winding phase per unit z) so a
solved station seeds the next. Premise test (`twist_premise.py`, criteria in
the header): on five pitch-locked 9.6 µm stations of PHerc 1667 (z 10000–
30000), take the winding phase around rings of radius 150/250/350 px about
the core and fit its angle against z.

Result: not viable — RMS residual 34–50° against a 15° bar. But the
diagnostic matters more than the verdict: the rings cross **6–17 windings**
at fixed radius, varying by station, and the phase coherence around a ring
is 0.02–0.25. A circle about the centroid is not an iso-winding curve on a
non-circular scroll, so "the phase on a ring" is undefined — the same
failure mode that retired the polar pitch estimators. A field-frame version
(phase tracked along the sheet normal from the core, per sector) would be
the honest design; it is not built. Shelved, not killed.
