# B — Is the ranker's slab preference a burial preference?

**Cost:** no docking. One analysis script over existing poses.
**Type:** reanalysis.

## Goal

Appendix C.8 now asserts a mechanism: gnina prefers slab poses "because the gnina ranking orders poses
by predicted affinity while the endpoint scores placement, and the two are not aligned on this
receptor". That is an assertion. This package turns it into a measurement, and it is the cheapest
honest route to the physics-or-assumption question. **Do it before any re-docking.**

## What is already established

Over the 12,315 control-panel AutoDock poses:

| ordering | usable, ranks 1-5 | usable, ranks 6-10 | Spearman with slab occupancy |
|---|---|---|---|
| gnina | 76.4% | 89.2% | −0.211 (p = 7e-124) |
| Vina | 82.6% | 83.1% | −0.004 (p = 0.63) |

So the learned ranker prefers slab poses and the physics score is indifferent. **Why** is not measured.

## B1 — mechanism

**Hypothesis.** gnina's CNNaffinity rewards burial and enclosure. The transmembrane band is the most
enclosed region of a channel searched without a bilayer. So the slab preference is a burial preference
expressing itself through the receptor's geometry.

**Steps.**

1. For every AutoDock pose on both Orai panels, compute a burial count: the number of receptor heavy
   atoms within 4.5 Å of any pose heavy atom. Pose files are under
   `Dockings/Orai_JKU_MGLTools_exh128/mgl_tools/` for the experimental panel; receptors are in
   `posebusters_results/_orai_jku_staging/receptors/`.
2. Test whether gnina rank tracks burial (Spearman, per unit and pooled).
3. Test whether slab occupancy still associates with gnina rank **after** controlling for burial
   (partial correlation, or stratify by burial decile and check whether the gnina gradient survives
   within strata).
4. Report the partial association, not a claim of clean separation. Burial and slab occupancy are
   geometrically entangled on a channel and cannot be fully disentangled from observational data.

**Interpretation.**
- gnina rank tracks burial, and the slab association largely disappears once burial is controlled →
  the mechanism is confirmed. The slab rule is partly a proxy for "not buried", and the appendix
  sentence becomes a measured statement.
- The slab association survives controlling for burial → something else is going on, and the current
  appendix sentence is wrong and must be softened.

## B2 — generality

**Question.** Does gnina prefer buried poses on ordinary soluble receptors too, or only here?

**Steps.** Run the same burial-versus-gnina-rank test on the 303-complex calibration benchmark, using
`posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report/per_pose_metrics.csv` for
the rank columns (`autodock_rank`, `optimized_rank`) and the corresponding pose files.

**Interpretation, and both outcomes are useful.**
- **gnina prefers buried poses on soluble receptors too.** Then the Orai1 effect is a general property
  of an affinity ranker meeting a receptor whose most buried region is a place ligands should not sit.
  This is the stronger and more transferable claim, and it generalises well beyond this thesis.
- **No burial preference on soluble receptors.** Then the effect is specific to this receptor or to
  membrane proteins without a bilayer. Narrower, but sharper, and it points directly at package D.

## Deliverable

One appendix paragraph and one figure. A sentence in the Discussion replacing the current mechanistic
assertion with the measured version.

## Risks

- **Entanglement.** Burial and slab occupancy are correlated by the receptor's shape. Report partial
  associations and say plainly that the two cannot be fully separated observationally.
- **Burial metric choice.** A 4.5 Å heavy-atom count is one choice among several. Check the conclusion
  survives at 4.0 and 5.0 Å before reporting.
- **Do not reuse the PandaMap contact counts as a burial proxy.** They were computed only on
  PoseBusters-valid, outside-slab, top-three poses, which excludes exactly the in-slab poses this
  analysis needs. Using them would build the answer into the question.

## Done when

Appendix C.8's mechanism sentence cites a number instead of asserting a mechanism, and the generality
question has a yes or no.
