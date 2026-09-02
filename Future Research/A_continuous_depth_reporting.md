# A — Report placement as a depth, not a flag

**Cost:** no docking, no new analysis run. One figure and two paragraphs.
**Type:** reanalysis of committed data, plus rewriting.

## Goal

The placement endpoint is reported as a binary: a pose is inside the slab or outside it. The continuous
depth of every pose is already stored and never shown. Showing it discloses how close the surviving
poses sit to the line.

## Why it matters

The retained usable AutoDock poses hug the boundary on both sides:

| side | poses | axial fraction |
|---|---|---|
| extracellular | 22 | 1.00 to 1.61 |
| cytosolic | 43 | −0.02 to −0.88 |

Fifteen of the 65, or **23%**, lie within 0.3 slab-widths of the cut, which at a 22.6 Å slab is about
7 Å, roughly one ligand radius. The closest cytosolic pose is 0.02 widths below its plane. A reader
currently sees a clean binary and cannot tell that a quarter of the surviving poses would change status
under a small shift of the threshold.

This also pre-empts the obvious examiner question, which is whether the result is an artefact of where
exactly the planes were drawn.

## Prerequisites

None. Everything needed is already committed.

## Data

`posebusters_results/_orai_matched_root/orai_jku/transmembrane_filter/tm_pose_classification.csv`
and the matching `orai_benchmark/...` file for the control panel.

**Use the `_orai_matched_root` tree.** The plain `posebusters_results/orai_jku/` tree is one generation
behind and has already caused one false finding in review.

Columns, all present and populated (1,424 of 1,428 rows):

| column | meaning |
|---|---|
| `axial_t` | axial coordinate in Å along the pore axis, 0 at the Arg91 plane |
| `axial_frac` | `axial_t` divided by slab thickness; **0 = Arg91 plane, 1 = Glu106 plane** |
| `radial_dist` | distance from the pore axis in Å |
| `frac_atoms_in_layer` | fraction of the pose's heavy atoms inside the slab |
| `majority_in_layer` | whether most atoms are inside |
| `in_transmembrane` | the binary verdict actually used, on the centroid |

The slab is `0 <= axial_frac <= 1`. Admissible is outside that interval.

## Steps

1. Write a small script under `Scripts/Analysis/`, or extend
   `orai_transmembrane_exclusion.py`, to emit one panel per tool: the distribution of `axial_frac` for
   retained poses, with the two ring planes marked at 0 and 1, the slab shaded, and the modulator and
   control panels overlaid. A violin or a histogram with a rug both work; the rug is better here because
   the experimental n is 12 units.

2. Add a second panel or an inset showing `frac_atoms_in_layer`, which exposes the poses that straddle
   the plane. The centroid rule calls these in or out; the atom fraction shows they are neither.

3. Compute the threshold sensitivity: usable yield with the slab padded by **±2 Å and ±5 Å** on each
   plane. Report it as a short table. This is the number that answers "is the result an artefact of
   where you drew the line".

4. Register the new output in `Scripts/Analysis/REGENERATE.md` with the exact command.

5. Copy the figure into `thesis_latex/media/media/imageNN.png`, choosing the next free number, and
   verify with `md5sum`. Add the `\includegraphics` block and caption in the Orai1 results section.

6. Add two paragraphs: one describing the distribution, one giving the padding sensitivity.

## Deliverable

One figure, one small sensitivity table, two paragraphs.

## Decision to make before starting

**Keep the binary endpoint as primary.** Every statistic in the chapter is built on it, and switching
the primary endpoint would invalidate the Mann-Whitney, Fisher and Cliff's delta results that were just
corrected under M4. The depth view is a disclosure alongside, not a replacement.

## Risks

- **Low.** The main hazard is scope creep into re-defining the endpoint. Resist it.
- The chapter is page-constrained; the short build sits at 142 pages. A new figure costs roughly a page.
  Decide in advance whether the figure goes in the Results or in Appendix C.8 next to the placement
  geometry section, which is the cheaper home.

## Done when

A reader can see where the surviving poses sit relative to the planes, and knows what happens to the
yield if the planes move by 2 and 5 Å.
