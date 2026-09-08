# E — Write down the recovered MD provenance

> **UPDATE 2026-09-03 — read this first.** The timing caveat below is now resolved by measurement,
> and two of the "keep on the unknown list" items have moved. The frame interval was measured from
> bulk-water displacement at about 0.66 ps, so `bar3` spans roughly 0.33 ns and `fix1` about 38 ps.
> `fix1` is a fixed-protein heating stage, which makes Fr0 the minimised pre-dynamics structure
> rather than an equilibrated one. The homology model is identified as Frischauf et al. 2015
> (ModelArchive `ma-akdjp`), and no publication describes this membrane, so Najjar et al. 2025 must
> not be cited for the setup. Full evidence and drop-in replacement text:
> `Scripts/Analysis/MD_FRAME_PROVENANCE_ANSWERS_2026-09-03.md`.


**Cost:** one paragraph, no computation. **Priority: do this first.**
**Type:** rewriting.

## Goal

Appendix D declares the receptor's provenance unrecorded. Most of it is now recoverable from the files
delivered on 2026-09-02. Replacing the declaration with the facts converts a stated limitation into
reported provenance, and it explains an anomaly the thesis measured but could not account for.

## Why this is the highest value per unit effort

Three separate things close at once.

1. The examiner review flagged the missing provenance under reproducibility and under red-flag item 13,
   "ligand or receptor preparation is not reproducible", where the Orai1 arm was the partial-present case.
2. Appendix D.1 reports that Fr0 lies 5.88 to 6.02 Å from the three later frames and is "not a member of
   the ensemble the later frames sample", and attributes this only to its `START` label. The reason is
   now known and is mechanistic: **Fr0 is not a production frame**. It is frame 0 of a
   protein-restrained equilibration run in which the protein was held while lipid and water relaxed.
3. It removes the awkward framing that the four frames are "four unlabelled conformational samples of one
   membrane simulation rather than timed points along a characterised trajectory". Three of them are now
   identified frame indices of a named production run.

## Prerequisites

None. Every fact below is already established and recorded in
[`README.md`](README.md). No file needs to be opened to write this.

## Steps

1. Open `thesis_latex/body_appendix_short.tex` and find the provenance paragraph at **line 570**,
   beginning "The provenance of the model stops at the coordinates, and the gap is stated here rather
   than left silent."

2. Replace the list of unknowns with what is now known:

   - built with **CHARMM-GUI v3.7**, 6 April 2022, assembling a previously generated lipid bilayer,
     protein and pore water
   - **CHARMM** force field (PSF EXT with CMAP)
   - bilayer of **426 lipids**: 213 cholesterol, 71 sphingomyelin, 71 DLPC, 71 DLPE
   - **40,123 TIP3** waters
   - **109 K+ and 150 Cl-**, and **no calcium**
   - production trajectory `Orai1_WT_bar3.dcd`, 500 frames; restrained equilibration
     `Orai1_WT_fix1.dcd`, 501 frames
   - **Fr300, Fr400, Fr499 are `bar3` frames 300, 400 and 499**; **Fr0 is `fix1` frame 0**, all
     confirmed by CA-RMSD after superposition at 0.0005 to 0.0006 Å

3. Keep two items on the unknown list, and say why:
   - the real time per frame, because the DCD header's `dt = 0.049 ps` is almost certainly unscaled
   - the homology-model template PDB identifiers

4. Add one sentence to **Appendix D.1**, near the Fr0 spread measurement, stating that Fr0 comes from the
   restrained equilibration rather than the production run, which is why it sits outside the ensemble.
   This turns an observation into an explanation and costs one sentence.

5. Add one sentence recording that the simulation contained no calcium. **Frame it correctly.** The
   docking receptor faithfully reflects its source; this is an inherited property of the source
   simulation, not something the preparation removed. Do not present it as a preparation defect.

6. Check whether Chapter 3.1.2 in `body_main_short.tex` needs a matching touch. It currently says the
   snapshots "represented limited receptor flexibility" and cross-references Appendix D.

7. Rebuild and confirm:
   ```
   cd thesis_latex && latexmk -pdf -f Thesis_short.tex
   ```
   Expect no page-count change from a paragraph swap of similar length. Confirm 0 undefined references.

## Deliverable

One rewritten paragraph in Appendix D, one added sentence in D.1, one added sentence on the ion content.

## Risks

- **Overclaiming the timing.** The single real hazard. Do not convert 500 frames into nanoseconds
  without the run input. If a reviewer asks how long the trajectory was, the correct answer today is
  that the frame count is known and the timestep is not.
- **Scope creep.** It is tempting to also re-examine whether Fr0 should have been included at all, given
  it is now known to be a restrained structure. That is a separate argument and belongs in the
  Limitations chapter, not here. The thesis already defends including Fr0 as a deliberate
  minimised-versus-thermal contrast, and that defence is now better supported, not worse.

## Done when

Appendix D states the composition, the force field, the builder, the trajectory files and the frame
indices; Appendix D.1 explains the Fr0 offset by its origin; the timing and the template IDs remain
explicitly unknown; and the build is clean.
