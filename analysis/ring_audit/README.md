# Ring-artifact confound audit (item 10.A.3) — IMMATERIAL at this station (2026-09-02)

Ring artifacts are concentric about the reconstruction centre: fake circular
"windings." Audit, declared first in `ring_audit.py`: polar median destripe
about the image centre, re-solve the pitch-locked 9.6 µm field on PHerc 1667
z4634, measure the change. Material if median |ΔW| > 0.25 winding or the
per-wrap scale moves by > 0.05.

| | |
|---|---|
| ring profile RMS / image high-pass RMS | 3.6 / 30.3 (12 %) |
| field change after destripe, median |ΔW| | **0.068** winding (p90 0.21; 0.2 % of pixels > 0.5) |
| windings per true winding, original → destriped | 0.985 → 0.968 |

**Verdict: immaterial.** Ring artifacts do not drive the field or the
fringe-scale numbers at this station. One station, one scroll; the audit
script runs unchanged on any other.
