# Face pairing (item 9.2) — KILLED at the premise, both statistics (2026-09-02)

**Question.** The 9.6 µm overcount comes from a sheet's two faces being read
as two windings. Can a ridge pair be classified as *two faces of one sheet*
versus *two different sheets* from the image alone? If so, the multi-band
ceiling measured in `../jointband` (78 %) becomes reachable.

**Labels.** PHerc 1667 z4634, 9.6 µm; ungated ridge-adjacency pairs
(237,547); a ridge within 4 px of an accepted wrap's mesh points belongs to
that wrap. Same-sheet pairs: 38,058 (median separation 8 px ≈ one papyrus
thickness). Different-adjacent-sheet pairs: 17,120 (median 12 px). The
labels are physically plausible; a 4 px membership radius is a choice and is
stated.

| statistic (criteria declared in each script header) | AUC same vs different |
|---|---|
| minimum CT intensity along the segment between the ridges | 0.649 |
| mean intensity along the segment | 0.592 |
| gap length | 0.372 (i.e. 0.63 reversed: same-sheet gaps shorter) |
| along-ridge texture correlation, ±16 / ±32 / ±64 px (striation matching) | 0.483 / 0.468 / 0.452 |

Bar: promising ≥ 0.80, kill < 0.70. **Both killed.**

**Reading.** Between tightly packed sheets the dark gap is one or two voxels
at 9.6 µm and is not reliably darker than a sheet interior; and the fibre
texture along a face does not correlate more with its own sheet's other
face than with the neighbouring sheet. At this resolution the local image
does not carry the face/sheet distinction. Consistent with the 22 % of
points in `../jointband` that have no correct band at all, and with the
resolution cliff in `../prior_art/topography_2026.md`. Face pairing is
retired as an image-local method at 9.6 µm; a thickness prior (same-sheet
gaps ≈ 8 px) is the only residual signal, weak (0.63), and is noted, not
built.
