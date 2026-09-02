# Plan: does the Orai1 slab rule encode physics or assumption?

**Date:** 2026-09-02 · **Status:** proposal, nothing executed

## Context

The M4 correction established that the transmembrane exclusion interacts with the ranking step. On the
control panel, over 12,315 poses from 308 unrelated ligands, the gnina order puts 76.4% of its first
five ranks outside the slab against 89.2% of ranks six to ten, and rank correlates with slab occupancy
at Spearman −0.211 (p = 7e-124). The Vina order moves 0.5 points across the same bands and correlates
at −0.004 (p = 0.63). So the learned affinity ranker prefers membrane-buried poses and the physics-based
search does not.

That leaves the question the chapter cannot currently answer. Do ligands land in the slab because the
receptor model is missing the physics that would keep them out, namely the bilayer and the permeant
calcium, or because the tools genuinely favour that region? The slab rule is an operational hypothesis,
and the thesis says so, but it has never been tested against a more complete receptor.

A second, unrelated improvement is free: placement is reported as a binary in-or-out when the continuous
depth of every pose is already stored.

> **UPDATE, 2026-09-02, same day.** The MD system arrived from the collaborating group while this plan
> was being written. It changes three of the four work packages. Finding 1 below is superseded; see the
> boxed revision after it, and the revised D and C sections.

## Three findings that shape this plan

1. ~~**The membrane and ion coordinates do not exist in this project.**~~ **SUPERSEDED.** They arrived on
   2026-09-02 at `.../Orai1-Docking-MScProjekt-DMann/Orai/MDFiles/`: a CHARMM topology
   `hORAI1_WT_ARN.psf` (182,692 atoms), a production trajectory `Orai1_WT_bar3.dcd` (500 frames) and a
   protein-restrained equilibration `Orai1_WT_fix1.dcd` (501 frames). Work package D is unblocked and
   work package C is weakened. The system holds 426 lipids (213 cholesterol, 71 sphingomyelin, 71 DLPC,
   71 DLPE), 40,123 TIP3 waters, 109 K+ and 150 Cl-, and **no calcium at all**. The four snapshots were
   identified exactly by CA-RMSD after superposition: Fr300, Fr400 and Fr499 are `bar3` frames 300, 400
   and 499 at 0.0005 Å, and Fr0 is `fix1` frame 0 at 0.0006 Å, which is why it sits 5.8 Å outside the
   ensemble the other three sample.
2. **All five depth columns already exist and are populated** in
   `tm_pose_classification.csv`: `axial_t`, `axial_frac`, `radial_dist`, `frac_atoms_in_layer`,
   `majority_in_layer`, with 1,424 of 1,428 rows non-null. Work package A needs no new computation.
3. **The preparation path strips ions twice**, once in PDBFixer and again via
   `prepare_receptor -U nphs_lps_waters_nonstdres`. Any calcium added to a receptor will be removed
   before docking unless that path is deliberately bypassed.

---

## Work package A — report placement as a depth, not a flag

**Cost:** no docking, no new analysis run. A figure and two paragraphs.

**What.** Replace, or rather supplement, the binary usable count with the distribution of `axial_frac`
across retained poses. Zero is the Arg91 plane, one is the Glu106 plane, the slab is the interval
between, and everything outside is admissible under the current rule.

**Why it matters.** The retained usable AutoDock poses sit hard against the boundary:

| side | poses | axial fraction |
|---|---|---|
| extracellular | 22 | 1.00 to 1.61 |
| cytosolic | 43 | −0.02 to −0.88 |

Fifteen of the 65, or 23%, lie within 0.3 slab-widths of the cut, which at a 22.6 Å slab is about 7 Å,
roughly one ligand radius. The closest cytosolic pose is 0.02 widths below its plane. A reader currently
sees a clean binary and has no way to know that a quarter of the surviving poses would change status
under a small shift of the threshold.

**Deliverable.** One panel per tool showing the `axial_frac` distribution of retained poses with the two
ring planes marked, the slab shaded, and the modulator and control panels overlaid. Plus a sensitivity
line in the text giving usable yield as the slab is padded by ±2 and ±5 Å.

**Risk.** Low. The one judgement call is whether to keep the binary endpoint as primary. It should stay
primary, because every statistic in the chapter is built on it; the depth view is a disclosure, not a
replacement.

---

## Work package B — is the ranker's slab preference a burial preference?

**Cost:** no docking. One analysis script over existing poses.

This is the cheapest route to a real answer on "physics or assumption", and it should be done before
any re-docking.

**B1, mechanism.** For every AutoDock pose on both Orai panels, compute a burial count, the number of
receptor heavy atoms within 4.5 Å of any pose heavy atom. Then test whether gnina rank tracks burial,
and whether slab occupancy is explained by burial once it is controlled for. If gnina rank correlates
with burial and slab occupancy is largely a consequence of burial, the mechanism proposed in Appendix
C.8 is confirmed rather than asserted.

**B2, generality.** Run the same burial-versus-gnina-rank test on the 303-complex soluble benchmark,
where the receptors are ordinary crystal structures. Two outcomes, both informative:

- gnina prefers buried poses there too. Then the Orai1 effect is a general property of an affinity
  ranker meeting a receptor whose most buried region happens to be a place ligands should not sit. That
  is the stronger and more transferable claim.
- gnina shows no burial preference on soluble receptors. Then the effect is specific to this receptor or
  to membrane proteins, which is a narrower but sharper claim.

**Deliverable.** A short appendix paragraph and one figure, plus a sentence in the Discussion that
replaces the current mechanistic assertion with a measured one.

**Risk.** Low, and the outcome is informative either way. The one caveat is that burial and slab
occupancy are geometrically entangled on a channel, so the analysis must report the partial association
rather than claim to have separated them cleanly.

---

## Work package C — a calcium-restored Vina arm

> **Downgraded on 2026-09-02.** The delivered topology shows the simulation was run in KCl with **zero
> calcium**: 109 K+, 150 Cl-, no Ca2+. So the docked receptor is not missing an ion its source had, and
> the absence is inherited rather than introduced by the preparation. That removes the strongest
> motivation for this package. Worse, the E106 filter and the D110/D112/D114 ring reached their
> conformations *without* calcium present, so inserting it post hoc places an ion into side-chain
> geometry that never accommodated one. Keep this package only as a deliberate what-if, clearly labelled
> as such, and prefer work package D. If it is run at all, the honest framing is "what would an affinity
> ranker do if the acidic rings were neutralised", not "restoring the physiological ion".

**Cost:** small. Three ligands × four frames = 12 dockings at exhaustiveness 128, plus the gnina pass.
Comparable to a fraction of one benchmark arm.

**What.** Re-dock the three modulators into receptors carrying calcium at the selectivity filter, and in
a second variant also at the D110/D112/D114 accumulating region. Everything else identical to the
reported arm: exhaustiveness 128, energy range 6, 30 modes, seed 42, same box construction, same gnina
rescoring with `optimize_rank_by: cnn_affinity`.

**The comparison.** Slab occupancy and the gnina-rank-versus-slab correlation, with calcium against
without. The prediction is directional and falsifiable. If the ranker's slab preference is driven by an
unneutralised acidic ring, restoring the ion should weaken it. If the correlation is unchanged, the
preference is not about the missing ion and the slab rule survives that particular challenge.

**Implementation notes.**

- Bypass the ion stripping. PDBFixer removes heterogens and `prepare_receptor -U nphs_lps_waters_nonstdres`
  removes what survives. The calcium must be reinserted after preparation, directly into the PDBQT, with
  the correct AutoDock atom type and charge.
- Screen validity against the same receptor the search used. The Fr0 control panel already cost 733
  spurious minimum-distance failures from exactly this mismatch, so the calcium receptor must be carried
  through to PoseBusters unchanged.
- Place the ions from published Orai structure and function rather than by eye, and record the source.

**Risks, in order of seriousness.**

1. **This tests the ranker more than the search.** Neither the Vina nor the Vinardo function carries an
   electrostatic term, so a calcium ion reaches the search only through sterics, and a single small ion
   displaces very little volume. The measurable effect will be concentrated in gnina's convolutional
   rescorer, which does see atom types. That is not a flaw, since gnina is the component under
   suspicion, but the write-up must not present the result as a test of the physics-based search.
2. **A known atom-naming landmine.** This project has already been bitten once by the protein backbone
   atom name `CA` being parsed as calcium, which broke 25 of 308 template pairs before it was fixed in
   August. Introducing genuine calcium into these files re-enters that failure mode from the other side.
   Any parser touching the new receptors must be checked explicitly, and the fix that resolved the
   earlier bug must be re-verified against the new inputs.
3. **Ion placement is itself a modelling choice.** Replacing one assumption with another is only progress
   if the new one is better sourced. Run both the filter-only and the filter-plus-accumulating-region
   variants so the sensitivity to that choice is visible.
4. **Vina-only, by design.** DiffDock and EquiBind were trained on soluble crystal receptors. Feeding
   them an ion-bearing membrane channel is out of distribution and would degrade them unevenly, breaking
   the information parity the study rests on. This arm cannot become a three-tool comparison.

**What it can and cannot conclude.** It can say whether the learned ranker's slab preference survives
neutralising the acidic rings. It cannot say whether the slab rule is biologically correct, because that
needs a ligand-bound structure.

---

## Work package D — the membrane arm, now unblocked

**Status:** feasible. The data arrived on 2026-09-02.

**Extraction is exact, not approximate.** Because the frames are identified to 0.0005 Å, the lipid,
water and ion coordinates belonging to each docked snapshot can be pulled from the trajectory with no
matching ambiguity: `bar3` frames 300, 400 and 499, and `fix1` frame 0. Whatever is built will be the
true environment of that exact receptor conformation, not a re-solvated approximation.

**Plan.** For each of the four frames, rebuild a docking receptor carrying the protein plus the lipid
within a shell of the protein surface, then re-dock the three modulators with Vina at the reported
settings and compare slab occupancy against the protein-only arm. Two variants are worth separating:
lipid-only, and lipid plus the K+ and Cl- that were actually present.

**What it can conclude.** Whether the lipid-facing half of the slab exclusion is reproduced by physics.
If the bilayer sterically excludes poses from the outer band, that half of the rule is validated rather
than assumed. The pore half is a separate question the lipid cannot answer.

**Risks.**

1. **Out of distribution for the learned tools.** A lipid-bearing receptor is nothing like the soluble
   crystal structures DiffDock and EquiBind were trained on. This must stay a Vina-only arm or it will
   destroy the cross-tool parity the study rests on.
2. **Lipid is mobile, so one frame's lipid coordinates are one sample.** The exclusion this produces is
   physical but not unique. Reporting all four frames mitigates it.
3. **Box volume explodes.** The whole-receptor box already spans 91.8 to 102.0 Å per edge; adding a
   bilayer shell enlarges it further and dilutes a fixed search budget. Exhaustiveness may need raising
   for the comparison to be fair, and that must then be matched in the protein-only control.
4. **Atom typing.** ADFRsuite must type the lipid atoms sensibly, and PoseBusters must screen against
   the same receptor the search used, or the Fr0-style spurious-failure problem returns.

**Bonus deliverable, and the cheapest thing in this whole plan.** Appendix D currently states that the
modelling program, force field, lipid composition, solvent and ion treatment, trajectory length and
frame identity are all unrecorded. Most of that can now simply be written down: CHARMM-GUI v3.7 of
6 April 2022, a CHARMM force field, 426 lipids in the composition above, 40,123 TIP3 waters, 109 K+ and
150 Cl-, and the exact frame indices. That converts a stated limitation into reported provenance for the
cost of one paragraph. Two things must stay unstated: the real time per frame, because the DCD header
reports 0.049 ps which is almost certainly unscaled, and the homology-model template identifiers, which
are still not in evidence.

---

## Sequencing and scope

| package | cost | dependency | thesis or future work |
|---|---|---|---|
| A, continuous depth | free | none | thesis, do it |
| **E, write down the recovered provenance** | **free** | **none** | **thesis, do it first** |
| B, burial mechanism and generality | free | none | thesis if time allows |
| D, membrane arm | moderate | none any more | future work, now feasible |
| C, calcium what-if | 12 dockings | prep bypass | optional, downgraded |

**Work package E, added 2026-09-02.** Replace Appendix D's paragraph of unrecorded provenance with the
recovered facts, and state that Fr0 comes from the restrained equilibration rather than the production
run. This is the highest value per unit effort in the plan: it costs one paragraph, it closes a
limitation the examiner review flagged under reproducibility, and it explains the Fr0 outlier the thesis
already measured but could not account for.

A and B use only committed data and would strengthen the chapter that M4 just corrected. C is a genuine
new experiment with a falsifiable prediction and a bounded cost. D is a research proposal.

**Recommendation, revised 2026-09-02.** Do E first: it is one paragraph, it needs no computation, and it
turns a stated limitation into reported fact. Then A, which costs nothing and closes a disclosure gap an
examiner can see. Then B, which converts the mechanistic claim in Appendix C.8 from an assertion into a
measurement, with either outcome usable. D is now the right future-work headline rather than a data
request, and it is a genuinely good next project because the frame identification makes the environment
exact. C is optional and should be relabelled as a what-if if it is run at all.

**One caution that applies to all of it.** None of these packages can validate the Orai1 poses. They
test whether the exclusion rule is defensible, not whether the predicted binding modes are right. That
still needs mutagenesis, a competition assay, or a ligand-bound structure.
