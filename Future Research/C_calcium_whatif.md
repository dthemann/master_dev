# C — Calcium what-if

**Cost:** small. 12 dockings plus a gnina pass.
**Type:** additional computation.
**Status: DOWNGRADED on 2026-09-02. Prefer package D. Read the next section before starting.**

## Why this was downgraded

The original motivation was that a calcium channel had been docked without calcium, so the acidic rings
at the E106 filter and at D110/D112/D114 were an artificially attractive electrostatic sink.

The delivered topology removes most of that motivation. The simulation contained **109 K+, 150 Cl- and
zero Ca2+**. So:

1. **The docking receptor faithfully reflects its source.** The absence is inherited, not introduced by
   the preparation. Any write-up must say so, or it misattributes the problem to the thesis's own work.
2. **The side chains never accommodated calcium.** The filter and the accumulating-region aspartates
   reached their conformations in KCl. Inserting Ca2+ into those coordinates post hoc places an ion into
   geometry that never made room for it, which is a worse modelling error than leaving it out.

**Honest framing if it is run at all:** this asks what an affinity ranker does when the acidic rings are
neutralised. It is not "restoring the physiological ion". Label it a what-if.

## What it would still test

The gnina ranker's slab preference is measured and real: Spearman −0.211 with slab occupancy over
12,315 control poses, against −0.004 for Vina. If that preference is driven by an unneutralised acidic
ring, neutralising it should weaken the correlation. If the correlation is unchanged, the preference is
not electrostatic and package B's burial hypothesis gains support by elimination.

That is a falsifiable prediction, which is why this package is kept rather than deleted.

## Steps

1. Place ions at the selectivity filter on the pore axis at the E106 ring. Source the position from
   published Orai structure and function, not by eye, and record the source.
2. Second variant with ions also at the D110/D112/D114 accumulating region, so the sensitivity to the
   placement choice is visible rather than hidden.
3. **Bypass the ion stripping.** PDBFixer removes heterogens and `prepare_receptor -U
   nphs_lps_waters_nonstdres` removes what survives. The ion must be reinserted after preparation,
   directly into the PDBQT, with the correct AutoDock atom type and charge.
4. Dock the three modulators into four frames at the reported settings, then the gnina pass.
5. Screen with PoseBusters against the same receptor the search used.
6. Compare slab occupancy and the gnina-rank-versus-slab correlation, with and without.

## Risks

1. **This tests the ranker, not the search.** Neither Vina nor Vinardo carries an electrostatic term, so
   an ion reaches the search only through sterics and a single small ion displaces very little volume.
   The measurable effect will sit almost entirely in gnina's convolutional rescorer, which does see atom
   types. That is acceptable because gnina is the component under suspicion, but the write-up must not
   present it as a test of the physics-based search.
2. **A known atom-naming landmine.** This project has already been bitten by the protein backbone atom
   name `CA` being parsed as calcium, which broke 25 of 308 template pairs before it was fixed in August
   2026. Introducing genuine calcium re-enters that failure mode from the other side. Every parser
   touching the new receptors must be checked explicitly and the earlier fix re-verified against the new
   inputs. **Budget real time for this; it is the most likely thing to go wrong.**
3. **Placement is a modelling choice.** Running only one variant hides the sensitivity. Run both.
4. **Vina-only by design.** Same reasoning as package D.

## What it can and cannot conclude

**Can.** Whether the learned ranker's slab preference survives neutralising the acidic rings.

**Cannot.** Whether the slab rule is biologically correct. That needs a ligand-bound structure.

## Done when

There is a number for the gnina-rank-versus-slab correlation with the rings neutralised, and a sentence
saying whether the electrostatic explanation survived.
