# D — Membrane-aware docking arm

**Cost:** moderate. 12 dockings plus receptor construction. **This is the real research project here.**
**Type:** additional computation.
**Status:** unblocked as of 2026-09-02.

## Goal

Test whether the lipid-facing half of the slab exclusion is reproduced by physics. If a bilayer
sterically excludes poses from the outer band, that half of the rule is validated rather than assumed.

## Why this is now a good project

The frame identification is **exact**, so the environment extracted is the true environment of the
precise receptor conformation that was docked, not a re-solvated approximation:

| snapshot | source frame | CA RMSD |
|---|---|---|
| Fr300 | `Orai1_WT_bar3.dcd` frame 300 | 0.0005 Å |
| Fr400 | `Orai1_WT_bar3.dcd` frame 400 | 0.0005 Å |
| Fr499 | `Orai1_WT_bar3.dcd` frame 499 | 0.0005 Å |
| Fr0 | `Orai1_WT_fix1.dcd` frame 0 | 0.0006 Å |

Most membrane-docking studies have to build and equilibrate a bilayer and then argue it is
representative. This one does not.

## Prerequisites

- `hORAI1_WT_ARN.psf`, `Orai1_WT_bar3.dcd`, `Orai1_WT_fix1.dcd` from the OneDrive `MDFiles` folder.
- MDAnalysis 2.10.0, present in the `vina` env.
- Copy the three files into the project rather than reading across `/mnt/c` for every frame; the DCDs
  are about 1.1 GB each and the WSL mount is slow.

## Steps

1. **Extract.** For each of the four frames, write out protein plus lipid within a shell of the protein
   surface. Start at 6 Å and treat the shell radius as a parameter to vary, not a constant to fix.
   Decide explicitly whether to keep whole lipid molecules or truncate at the shell; whole molecules are
   cleaner and avoid dangling valences.

2. **Sanity-check the extraction** before docking anything. Confirm the protein heavy-atom coordinates
   in the extracted file are identical to the docked receptors already in
   `posebusters_results/_orai_jku_staging/receptors/`. If they are not, the frame mapping or the
   superposition is wrong and everything downstream is void.

3. **Prepare.** Run the extracted receptors through the same ADFRsuite path the reported arm used, and
   confirm the lipid atoms are typed sensibly rather than silently dropped. PDBFixer and
   `prepare_receptor -U nphs_lps_waters_nonstdres` both strip heterogens; this path must be bypassed or
   the lipid will vanish exactly as the ions did.

4. **Dock.** Three modulators into four frames with Vina at the reported settings: exhaustiveness 128,
   energy range 6, 30 modes, seed 42, same box construction. Then the same gnina pass with
   `optimize_rank_by: cnn_affinity`.

5. **Two variants worth separating:** lipid only, and lipid plus the K+ and Cl- that were actually
   present. The second is the faithful one.

6. **Compare** against the protein-only arm on slab occupancy and on the gnina-rank-versus-slab
   correlation. The protein-only numbers are already known: 45.8% of retained experimental AutoDock
   poses inside the slab, Spearman −0.332 experimental and −0.211 control.

7. **Screen** with PoseBusters against the same receptor the search used. Not doing this cost the Fr0
   control panel 733 spurious minimum-distance failures once already.

8. Register the command in `Scripts/Analysis/REGENERATE.md`.

## What it can and cannot conclude

**Can.** Whether the lipid-facing half of the exclusion is physical. If poses stop appearing in the
outer band once lipid is present, that half of the rule is earned.

**Cannot.** The pore half. Lipid does not occlude the lumen. And it cannot validate any binding mode.

## Risks, in order

1. **Out of distribution for the learned tools.** A lipid-bearing receptor is nothing like the soluble
   crystal structures DiffDock and EquiBind were trained on. **Keep this Vina-only** or it destroys the
   information parity the whole study rests on.
2. **Box volume.** The whole-receptor box already spans 91.8 to 102.0 Å per edge. Adding a bilayer shell
   enlarges it and dilutes a fixed search budget. Exhaustiveness may need raising, and if it is raised
   it must be raised in a matched protein-only control, or the comparison confounds membrane with
   search effort.
3. **One frame is one lipid sample.** Lipids are mobile, so the exclusion produced is physical but not
   unique. Reporting all four frames mitigates this; claiming a single definitive answer does not.
4. **Atom typing and charges.** Confirm ADFRsuite assigns lipid atoms types the Vina function
   understands. Gasteiger charges will be written but neither Vina nor Vinardo has a charge term, so
   they will not reach the objective; gnina's CNN does see atom types and will.
5. **The trade being made.** This swaps a stated geometric assumption for an unstated force-field one.
   It is still the right experiment, because lipid occupies the band for physical reasons rather than by
   fiat, but it is not automatically more objective than what the thesis does now. Say so in the write-up.

## Done when

There is a measured answer to whether the outer-band half of the slab rule survives contact with the
bilayer that was actually there.
