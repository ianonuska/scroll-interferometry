# Preregistration — 3D ink segmentation at prize resolution, and the winding field as page assigner

Committed 2026-09-02 before any run. Motivated by the organisers' 2026
open-problems note that the optimised 2.4 µm scan of PHerc. Paris 4 supports
direct volumetric ink segmentation (`scrollprize/ink_3d_dino_guided`), while
all eight 2027-eligible scrolls have only 8.6–9.4 µm scans.

## Block
PHerc 1667, volume `20251217075048-2.399um-0.2m-78keV-masked`, bbox
z 18408:18664, y 2000:3024, x 6352:7376 (2.4 µm voxels; 16 accepted
per-winding segments w018–w041 pass through it). Checkpoint
`ckpt_78k_fullsup.pth` (EMA weights), 256³ patches, the tool's defaults
otherwise.

## Run A — 3D ink at 2.4 µm (the favourable condition)
Expected to show ink. Output: per-voxel probability. Reported: fraction of
voxels above 0.5; a mid-block slice figure with the accepted wraps overlaid.
Human ink labels for 1667 are on Hugging Face (`ink/1667/…`) and are not on
the box at commit time; Dice against them is added when they are, and is
not required for the runs below.

## Run B — the same block downsampled 4× to 9.6 µm, same model
The question that matters for the prize scrolls: how much of the 3D-ink
capability survives at ~9 µm? Metric: Dice between run B's mask (upsampled
back) and run A's mask, treating A as reference. Declared: **useful at prize
resolution if Dice ≥ 0.5; not useful if < 0.25**; between is reported as
partial. A second reading is required regardless: if run A shows no ink in
this block, run B is uninterpretable and both are reported as such.

## Run C — page assignment with the winding field
Take run A's ink voxels (p ≥ 0.5). Truth: each voxel's page = the accepted
wrap whose mesh is nearest (within 12 voxels; others dropped). Our field: the
pitch-locked 9.6 µm solve at this z (`W_L2_locked2`), sampled at (y/4, x/4),
made relative to wrap w028's median. Assigned page = round(ΔW). Metric:
fraction of ink voxels whose assigned page offset equals the true offset.
Declared: **pass ≥ 0.80, kill < 0.50.** Also reported per true wrap, so the
scale bias (which grows with distance from the reference) is visible rather
than averaged away.

## What is not claimed
Nothing about scrolls other than 1667; nothing about letters; run B is a
resolution experiment on one block, not a verdict on the 3D model.

## Readout — written after the runs, 2026-09-02; criteria above unedited
Run A: 2.0 % of voxels p ≥ 0.5, deposits on sheet faces. Run B: Dice 0.087
against run A → **not useful at prize resolution** by the declared bar (one
block, one patch, reflect-padded). Run C: accuracy 0.380 → **kill**; the
failure is the field's scale bias (+2 pages at 5 wraps, +4 at 10). Details
and figure in `../validation/ink3d/`.
