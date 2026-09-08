# Defence Deck — Every Number in the Speaker Notes, Explained

Source deck: `Mann_MasterThesis_Defence_v3.pptx` (15 slides, notes on 13).
Every figure below was checked line-by-line against `body_main_short.tex` and
`body_appendix_short.tex` on 2026-09-03. Line references are to those two files.

Reading key for each entry:

- **Is** — what the number actually measures.
- **Over** — the denominator. This is where most examiner questions land.
- **Source** — where it lives in the thesis.
- **Trap** — the way the number can be mis-said under pressure.

---

## Slide 2 — Motivation and Research Questions

### `exhaustiveness 128`

- **Is** the Vina search-effort setting carried forward for every reported AutoDock
  number. Vina spreads a fixed number of Monte Carlo runs over whatever box it is
  given, so with a whole-receptor box this had to become a decision rather than a
  default.
- **Over** all 303 analysed complexes, and all 308 prepared ones. Every docking log
  records the same value, so effort is constant within the arm.
- **Source** appendix :163, :169, :1199.
- **Trap** the deck note calls the AutoDock arm "Vina search at exhaustiveness 128
  followed by a gnina CNN re-ranking pass". That is exactly right and it is the
  sentence that defuses the "AI versus physics" framing. Say it in that order.

### `three` (three configured workflows)

- **Is** AutoDock Vina, DiffDock, EquiBind, each as a full pipeline rather than a bare
  model. Two of the three carry a learned component after the deck's own definition,
  because AutoDock is CNN-rescored.

---

## Slide 3 — Tools, Pipeline and Metrics

No hard figures in the note prose. Two definitions carry weight.

### `rank-1` versus `best-of-top-k`

- **Is** precision at 1 versus recall at k, in information-retrieval terms. The deck
  note makes that translation explicitly and it is the single most useful bridge for a
  computer-science examiner.
- **Trap** best-of-top-k is a retrospective oracle. The thesis says so at main :868.
  Never let it be quoted as prospective accuracy.

### `22 physics checks` (slide face, not notes)

- **Is** the PoseBusters battery, in three groups: chemical validity, intramolecular
  geometry, intermolecular geometry.
- **Source** the failure decomposition is appendix table `tab:appendix-pb-decomposition-full`.
- **Trap** of 72,219 poses across eight reported variants, **none** fails a chemical
  check. The discrimination is almost entirely intermolecular, and within that almost
  entirely one check, minimum distance to protein.

---

## Slide 4 — Two Datasets

### `303`

- **Is** the analysed benchmark complexes. 308 were prepared; 303 carry a complete
  record on all three pipelines and are the cohort every accuracy number uses.
- **Trap** 308 appears in the Orai1 control panel (308 ligands × 4 frames = 1,232
  complexes) and in preparation counts. 303 is the accuracy denominator. Do not mix them.

### `24 heavy atoms`, `359 Da`, `5 rotatable bonds`

- **Is** cohort medians of the benchmark ligands.
- **Over** the 303 analysed complexes. Heavy-atom interquartile range 16 to 32, full
  range 6 to 64. Molecular-weight interquartile range 236 to 462 Da, range 100 to 872.
- **Source** main :88, appendix :669 and :695.
- **Trap** the distribution is right-skewed and includes peptides, nucleotides and
  sugars. The median is not a typical drug-like molecule.

### `three ligands` and `four frames`

- **Is** the Orai1 experimental panel: 2abp-NH2, Synta-66, GSK-7975A, docked against
  four molecular-dynamics snapshots of one homology-model trajectory.
- **Over** 3 × 4 = 12 frame-ligand units per tool, 40 poses per ligand after the
  top-ten cap, 120 poses per tool.
- **Trap** the four frames are correlated snapshots of one system, not four independent
  receptors. Every cross-panel contrast is therefore re-tested at ligand level, and the
  note's "almost none survives" is accurate.

---

## Slide 5 — Optimisation Gains

### `over 70 percentage points`

- **Is** the pooled physical-validity gap between the AutoDock arm and the raw learned
  output. AutoDock rescored sits at 98.9 % pooled, raw DiffDock at 24.4 %. The gap is
  74.5 points. Against raw EquiBind at 2.9 % it is 96.0 points.
- **Over** poses, not complexes: 8,976 AutoDock poses, 9,013 raw DiffDock poses, 9,090
  raw EquiBind poses.
- **Source** main :294, :295, :300, :406.
- **Trap** this is the *pooled pose* basis. The slide face quotes the *per-complex
  median* basis instead, which gives different numbers for the same phenomenon. Keep
  the two apart, see the next entry.

### `13 %` and `0 %` raw, `93 %` and `80 %` optimised (slide face)

- **Is** per-complex **median** PoseBusters yield. DiffDock rises from a 13.3 % raw
  median to 93.3 % with smina, a Hodges-Lehmann gain of +62.1 points. Unguided EquiBind
  rises from 0 % to 80.0 % with gnina.
- **Over** the 303 complexes, each contributing one yield fraction.
- **Source** main :319.
- **Also worth having** DiffDock reaches 96.6 % with gnina instead of smina, so the two
  refiners are near-interchangeable for it. EquiBind manages only 56.7 % with smina, so
  gnina is decisive there. That asymmetry is why the two arms carry different refiners.
- **Trap** AutoDock's raw median is already 100 % with 99.6 % pooled validity, and
  gnina rescoring moves the pooled figure slightly the **wrong** way, to 98.9 %
  (Wilcoxon p = 0.081, Hodges-Lehmann difference exactly zero, 293 of 303 complexes
  unchanged). AutoDock is adopted for what rescoring does to *ranking*, not geometry.

### `at most four of 303`

- **Is** the number of complexes that requiring PoseBusters validity removes from the
  headline endpoint, at any reported depth. The exact counts are two at rank-1, four at
  top-5 and three at the deeper pools.
- **Source** main :324, main :777, appendix :1194.
- **Why it matters** this is the mitigation for the circularity charge. If the validity
  filter only ever moves four complexes, then even a partly self-fulfilling filter cannot
  be manufacturing the tool ordering.

### `33 of 8,976`

- Worth having ready though not in the notes: only 33 raw AutoDock poses fail
  PoseBusters at all, split as internal clash 12, minimum distance to protein 20,
  internal energy 1. Source main :318.

---

## Slide 6 — Placement and Form

### `36.6 %`, `34.0 %`, `18.2 %` (rank-1 placement)

- **Is** the share of complexes whose **rank-1** pose is both PoseBusters-valid and
  within 2 Å in-place RMSD of the crystal ligand. In counts, 111, 103 and 55 of 303.
- **Over** 303 complexes each.
- **Source** main :347, :382, :834.
- **Trap** the AutoDock–DiffDock gap here is +2.6 points and is **not** significant
  (Newcombe hybrid-score interval −4.2 to +9.4, exact McNemar p = 0.509). The separation
  appears only from top-5 onwards.

### `52.5 %`, `52.1 %`, `33.3 %` (rank-1 form)

- **Is** the same poses judged on best-fit Kabsch RMSD ≤ 1 Å, which superimposes the two
  molecules first and therefore reports internal conformation alone.
- **Source** main :396.
- **Trap** the AutoDock–DiffDock contrast is a tie at every depth, p_holm = 1.000, with
  47 against 46 discordant complexes at rank-1, 30 against 31 at top-15, 25 against 26 at
  top-30. That symmetry in the discordant counts is the cleanest way to show the tie is
  real rather than underpowered noise, though no equivalence margin was prespecified.

### `five` exhaustiveness values

- **Is** the AutoDock ladder: exhaustiveness 18, 32, 64, 92, 128, all sharing prepared
  ligands, receptors, box and seed, so exhaustiveness is the only variable.
- **Trap** the *selection field* is eight arms, not five. Three of the five rungs were
  additionally gnina-rescored. If the examiner asks how many configurations the winner
  was picked from, the answer is eight.
- **Source** appendix :169, :180.

### `0.4 to 0.8 percentage points` (winner's-curse correction)

- **Is** the standard winner's-curse estimator applied **per family over the three
  refiner arms** each learned pipeline offers.
- **Trap — important** this is *not* the AutoDock correction. The eight-arm AutoDock
  field carries a larger one, near **two points at rank-1 and one and a half at top-15**
  (appendix :248). The thesis flags explicitly that the two are different quantities and
  neither supersedes the other. If you quote 0.4 to 0.8 for AutoDock you will have
  understated your own optimism by a factor of three, and the appendix contradicts you
  in writing.
- **Source** main :866, appendix :248.

### `3.0 percentage points` (largest selection margin)

- **Is** the largest winning margin on the variant-selection endpoint at the top-15
  depth at which the choice was made, expressed as a share of the 303. It bounds how
  much optimism the selection could have introduced in any one family.
- **In counts** eight complexes for AutoDock Vina, one for DiffDock, nine for EquiBind.
  Only the AutoDock margin is separable, and it is the one drawn from the largest field.
- **Source** main :866.

### `9.8 percentage points` (minimum detectable difference)

- **Is** the smallest rank-1 AutoDock–DiffDock difference the study could have detected
  at 80 % power. The observed difference is +2.6 points, so the study is underpowered
  by a factor of nearly four on that one contrast.
- **Over** the **112 discordant complexes** only. A McNemar test draws information from
  discordant pairs alone, which is why the effective sample is 112 rather than 303.
- **Source** main :386.
- **Trap** this bounds the claim, it does not rescue it. No equivalence margin was
  prespecified, so a non-significant result is not evidence of equivalence. Say
  "unresolved", never "the same".

### `1 Å` and `2 Å` thresholds

- **Is** the form threshold and the near-nativeness threshold. 2 Å is the field
  convention for RMSD-to-crystal. 1 Å for form is **not** a published convention. It was
  fixed before analysis as half the in-place threshold, and the appendix resolves the
  endpoint across ten thresholds from 0.25 to 2.5 Å so the choice cannot be doing work.
- **Source** main :380 footnote, appendix table `tab:results-near-native-form`.

### `4 Å` centroid clustering

- **Is** the coarse site criterion used two slides later. A cluster centroid within 4 Å
  of the crystal ligand site.
- **Why the note raises it here** because placement is never measured directly on this
  slide. It is a residual, sqrt(in-place² − form²), and that residual bundles orientation
  with position. A ligand in the right pocket but flipped end-for-end counts as
  placement-limited. The 4 Å centroid criterion is the metric that actually tests location.

---

## Slide 7 — Criteria at Rank-1 and Top-15

### `2.6` widening to `6.6` points

- **Is** the AutoDock-over-DiffDock rank-1 lead at two thresholds. At 2 Å in-place the
  rates are 36.6 % against 34.0 %, a gap of 2.6 points. At 1 Å they are 22.1 % against
  15.5 %, a gap of 6.6 points.
- **Source** appendix table `tab:results-near-native-form`, RMSD block, rank-1 rows.
- **Reading** the leaders are closer at the conventional threshold than at a strict one.
  Tightening the requirement favours AutoDock.

### form `nominally ahead` for DiffDock at 2 Å

- **Is** Kabsch RMSD ≤ 2 Å at rank-1: AutoDock 78.9 %, DiffDock 82.5 %. DiffDock leads
  by 3.6 points. At the primary 1 Å threshold AutoDock leads by 0.4 points, 52.5 against
  52.1. The crossover sits between them: at 1.75 Å the two read 74.6 and 76.2, DiffDock
  already ahead by 1.6.
- **Reading** the form tie is threshold-robust and the sign of the difference flips
  between thresholds, which is the strongest available evidence that it is a genuine tie.

### `1.49` to `5.04 Å` (AutoDock) and `1.47` to `1.93 Å` (DiffDock)

- **Is** median in-place RMSD at rank-1 and at top-15. AutoDock degrades from 1.490 to
  5.043 Å. DiffDock moves only from 1.467 to 1.929 Å.
- **Over** the **near-site cohort** only, see the caveat below. Pose counts are 190 then
  2,200 for AutoDock, 164 then 2,139 for DiffDock.
- **Source** appendix table at :1375–1377.
- **Reading** deeper AutoDock ranks add correctly shaped but badly placed ligands. Its
  ranking is informative only near the top. DiffDock's flatter pool is the mirror image:
  it is stable with depth but its ranker adds little.

### `48 %` to `78 %` placement-limited

- **Is** the share of near-site poses whose squared deviation is mostly positional
  rather than conformational, binned on r = form² / in-place² at cut-points 1/3 and 2/3.
  AutoDock moves 47.9 → 65.1 → 77.8 % across depths 1, 5, 15. Form-limited poses fall
  27 → 8 % over the same span.
- **Source** main :410, appendix :1375–1377.
- **Reading** the same statement as the median-RMSD degradation, in a different currency.
  Quote whichever one the question is phrased in.
- **For contrast** DiffDock's placement-limited share moves only 40.9 → 47.9 %.

### `8 Å` near-site trim

- **Is** the cohort restriction on the last three rows: poses whose centroid lies more
  than 8 Å from the crystal-ligand site are dropped. It means "outside the crystal-site
  neighbourhood", not "off the receptor".
- **Trap — volunteer this** those rows are pose-level and cohort-trimmed, so they are
  **not** comparable with the first five complex-level rows. Coverage runs 62.7 % to
  90.8 % of complexes depending on tool and depth, so the cohort itself changes with depth.

### the symmetry caveat

- The in-place and best-fit values each minimise over symmetry-equivalent atom mappings
  **independently** and may pick different mappings for the same pose, so the placement
  decomposition is approximate. Source main :132, :133.
- Separately, the two RMSD implementations agree to about 0.1 Å: rank-1 medians 1.490
  against 1.543 Å for AutoDock, 1.467 against 1.566 Å for DiffDock. Substituting the
  PoseBusters implementation into the 2 Å rank-1 gate moves **two of 909** tool-by-complex
  decisions. Rank-1 recovery would read 37.0 % and 18.5 % against an unchanged 34.0 %.
  The ordering is unaffected.

---

## Slide 8 — Binding Pocket

### `95 %`

- **Is** the retrospective three-tool oracle at the 4 Å crystal-site threshold, at
  cluster level. Pooling every candidate cluster from every tool, the crystal-compatible
  region is present for 95 % of complexes.
- **Trap — say this before you are asked** it is a *cluster-level oracle* at a *4 Å
  centroid* threshold. That is far coarser than a valid pose within 2 Å. It must never be
  quoted as 95 % pose accuracy. The slide deliberately labels it an oracle.
- **Source** main :474.
- **Related** across the full pose cloud rather than clusters, the pooled oracle reaches
  97.7 % at a median 0.34 Å. AutoDock alone reaches 88.1 %, DiffDock 84.5 %, and the two
  together already reach 97.7 %, so EquiBind adds nothing to pooled reach. Source main :476.

### `57 %` and the `38-point` gap

- **Is** the best inexpensive blind first-cluster rule against that oracle. Consensus
  ranking reaches **56.8 %**, its medoid-centred variant 55.4 %, confidence weighting
  55.8 %. The slide rounds the best to 57 %. 95 − 57 = 38 points.
- **Over** the 303 complexes.
- **Source** main :474.
- **Also have ready** the three rules are not separable from one another. They gain 7.3,
  5.9 and 6.3 points on the 49.5 % legacy pose-count rule, and all three clear Holm
  correction against it (p = 3.6e-4, 0.004, 0.008).
- **Reading** the gap is in *selection*, not in sampling. That reframes the whole
  problem as re-ranking.

### `0.45`, `0.42`, `0.38` and `0.46` (Jaccard)

- **Is** mean typed-interaction Jaccard similarity to the crystal contact fingerprint,
  over the union of each tool's top five valid poses. Precisely 0.450 for DiffDock over
  300 complexes, 0.417 for AutoDock over 301, 0.380 for EquiBind over 266. The
  AutoDock–DiffDock pair sits at 0.463.
- **Trap 1** Jaccard counts shared **and** non-shared contacts. It is a similarity, not
  a recall. The thesis says so in as many words.
- **Trap 2 — the depth artefact** DiffDock's lead exists only for the **union of its top
  five**. At rank-1 the two leaders are not separable on contact F1 either: AutoDock
  0.565 against DiffDock 0.529, with the thesis calling it inconclusive rather than
  equivalent. If you quote 0.45 without the depth qualifier you have overstated DiffDock.
- **Trap 3** the cohorts are overlapping, not identical, 300 / 301 / 266. On the common
  264 complexes the values are 0.474 / 0.431 / 0.381, so the ordering survives, but say
  which basis you are on.
- **Source** main :489, main :514–531.
- **Reading** two tools resembling each other more than either resembles the crystal
  means they hold *different* contact hypotheses from the reference, not better ones.

---

## Slide 9 — Orai1, Validity Overstates Usefulness

### `100 %` and `54 %` (AutoDock), `93 %` and `68 %` (DiffDock)

- **Is** PoseBusters validity and *usable* yield on the experimental panel. Usable means
  valid **and** on the receptor **and** outside the transmembrane exclusion slab.
  Precisely: AutoDock 120/120 valid = 100.0 %, 65/120 usable = 54.2 %. DiffDock 112/120
  = 93.3 % valid, 82/120 = 68.3 % usable. EquiBind 107/120 = 89.2 % valid, 1/120 = 0.8 %
  usable.
- **Over** produced poses, capped at the top ten ranked per frame-ligand unit, 3 ligands
  × 4 frames × 10 = 120 per tool.
- **Source** main :565–577.
- **Reading** all 120 AutoDock poses are geometrically sound and 55 of them are in the
  membrane. Validity is an admissibility filter, not evidence of binding. That is the
  slide's whole point.

### the exclusion slab

- **Is** the set of points between the plane of the six Arg91 Cα atoms and the plane of
  the six Glu106 Cα atoms. A pose is transmembrane when its heavy-atom centroid lies in
  that slab. Axial thickness is 22.4 to 22.8 Å, stable across frames, against receptor
  frames spanning 90.6 to 100.5 Å along the pore axis.
- **Trap — concede this first** the slab is a *modelling hypothesis*, not a measurement.
  It encodes an outer-pore binding hypothesis and may reject valid GSK-7975A or Synta-66
  poses, since both act at or near the pore mouth.

### `28.7` to `17.0` points, and `p` from `5.6e-4` to `above 0.2`

- **Is** the benchmark-to-experimental drop in AutoDock's usable yield, before and after
  matching the two panels' selection rules. Unmatched it is 28.7 points with Holm
  p = 5.6e-4, the only contrast surviving at both inferential units. Prefix-matched it
  falls to 17.0 points and clears no correction at either unit, with the frame-collapsed
  median moving 55.0 → 70.0 %, its adjusted p from 0.012 to above 0.2, and the
  frame-ligand p from 5.6e-4 to 0.054.
- **Why the panels differed** the experimental panel kept thirty Vina modes and then
  took the gnina-best ten. The control panel only ever wrote ten. Because gnina re-orders
  all thirty before the cut, **66 of the 120 retained experimental poses are Vina modes
  eleven to thirty** and would never have been written by a ten-mode run.
- **Source** main :625, :627, appendix :363, :367.
- **Trap — a number collision** `28.7` also appears on slide 12 as AutoDock's rank-1 to
  top-15 gain in percentage points. Two completely different quantities. If both come up
  in one answer, name the endpoint each time.

### `75 %` on the Vina prefix

- **Is** AutoDock's per-unit median usable yield when the top-ten cap is applied to its
  raw Vina order instead of the gnina re-ranking: 60.0 → 75.0 %, and pooled 54.2 → 65.8 %.
  The control figures are untouched.
- **Reading** the between-tool reversal on this slide is a property of the **secondary
  ranker**, not of the search. Concede it as such.
- **Source** appendix :367.

### the Synta-66 concordance

- **Is** the residue-recovery result, and the answer is no, it does not validate the
  published site. A decoy control contacts the same residues at the same rate, so the
  concordance is not specific. Have this ready; it is the strongest single objection to
  the Orai1 chapter and you are better off raising it yourself.

---

## Slide 10 — Orai1 Clustering

### `30 Å` on the controls, `33` to `44 Å` on the modulators

- **Is** median inter-tool consensus-site separation, that is the distance between two
  tools' largest-cluster centroids within one frame. Benchmark pooled median 28.4 Å over
  the 1,198 frame-ligand units carrying more than one tool. Frame-collapsed 29.3 Å.
  Experimental pooled 38.1 Å, frame-collapsed 40.3 Å. By tool pair the largest is
  DiffDock–EquiBind at 43.5 Å experimentally against 30.7 Å on the benchmark
  (Cliff's δ = +0.78).
- **Trap** the experimental estimate is **unstable**, moving 38.1 → 40.3 Å under
  frame-robust collapse. Only the benchmark value is a settled central estimate, and the
  thesis says so. Quote the ranges, not a point value, on the modulator side.
- **Source** main :642, :656, appendix :416, :418.

### `14 %` of control pairs and `none of 11`

- **Is** the share of AutoDock–DiffDock pairs whose consensus sites lie within 5 Å.
  14 % on the benchmark control panel, zero of eleven usable experimental pairs.
- **Over** eleven multi-tool frame-ligand units on the experimental side. Across all
  three pairings, 9 % of benchmark pairs and none of thirteen experimental pairs agree.
  Three-way agreement occurs in two of 450 benchmark units and in none of the single
  experimental unit that has all three tools.
- **Source** main :650, :658, :810.
- **Trap** eleven pairs is descriptive, not a corrected finding. Say "exploratory".
  AutoDock–DiffDock is nonetheless the only material agreement channel; EquiBind pairings
  sit at or below 4 % on the benchmark and zero experimentally.

### `two Angstroms` from symmetry folding

- **Is** the effect of folding every inter-tool distance over Orai1's six-fold symmetry
  about the Arg91→Glu106 pore axis. The pooled benchmark median moves 29.5 → 27.4 Å on a
  hydrogen-free centroid basis, or 29.3 → 27.3 Å on the shipped pipeline basis. Agreement
  within 5 Å moves 9.5 → 10.5 %, and the fold flips the verdict for about twenty of the
  roughly 2,100 constituent tool pairs. Under full continuous C∞ symmetrisation, a
  strictly stronger correction than any C6 treatment can be, the median moves only to
  27.3 Å. The thesis prints the unfolded baseline as a pooled median of 28.4 Å with
  9.1 % of 2,098 pairs agreeing within 5 Å; the recomputation used a slightly different
  pair set, so quote the direction and magnitude, not the pair count.
- **Why so little** each separation splits into an axial offset along the pore, a radial
  offset and an azimuthal chord, and rotation can shorten **only the chord**. Median
  axial separation is 15.4 Å, median chord 3.8 Å, and 76 % of pairs exceed 5 Å on the
  axial component alone. Tool sites sit at a median 2.2 Å from the pore axis, and a
  rotation moves a point at radius r by at most r, which bounds the achievable median
  reduction at 1.3 Å.
- **Source** appendix :416, and the recomputation recorded in the project memory
  `orai-c6-symmetry-refuted`.
- **Trap** the thesis text says all quoted separations are **unfolded** and does not
  print the folded median in the Results. The two-ångström figure is a defence-only
  number from the direct recomputation. Present it as "we recomputed it and folding moves
  the median by about two ångströms", not as a thesis result.
- **One honest caveat** on the *experimental* frame-robust statistic alone, C6 folding
  moves 42.1 → 30.8 Å, because those sites sit at a median radius of 18.3 Å. The thesis
  already discounts that statistic as unstable, so nothing rests on it, but do not be
  caught claiming folding never matters.
- **The closing argument** running cross-tool agreement in the examiner's own preferred
  symmetry-folded residue representation gives 12.6 % pooled agreement against the
  Cartesian metric's 9.5 %. The proposed correction, applied natively, reproduces the
  thesis's verdict.

### `no equivalence test`

- Filtering to valid, outside-slab poses improves bootstrap reproducibility but not
  spread. No equivalence test was run, so this **bounds** the filtering effect rather
  than proving the scatter is intrinsic. Say it that way.

---

## Slide 11 — Cost Per Pose

### the charged basis

- **Is** CPU-core-seconds divided by the 32 hardware threads of the workstation, plus
  GPU-seconds. It charges *device occupancy* rather than wall clock, so work that merely
  ran with more concurrency does not look cheaper.
- **Hardware** Intel Core i9-14900HX, 32 threads, NVIDIA RTX 4070 Laptop, 8 GB VRAM,
  32 GB system memory.
- **Source** main :710, appendix :330, :355.

### `0.74` wall-hours, `sixteen` ways, `10.27` GPU-hours

- **Is** the gnina rescoring pass that completes the AutoDock arm. It ran sixteen poses
  at a time on the GPU, consuming 10.27 GPU-hours of device occupancy in 0.74 wall-hours.
- **Trap** the 0.74 figure is a *scheduling result*, not a property of the method. It is
  the number that would flatter this arm.

### `15.00` charged hours and `68.5 %`

- **Is** AutoDock's total: 4.73 wall-hours of blind Vina search plus 10.27 GPU-hours of
  gnina occupancy. The rescorer is 68.5 % of the bill and is roughly twice as expensive
  as the search it re-ranks.
- **Source** main :710, :756, :850.

### `5.47` hours

- **Is** what AutoDock would total if the rescoring pass were charged at wall clock
  instead of occupancy, with rescoring falling to 13.5 % of the bill. That is the
  flattering accounting, and the thesis declines it.
- **Source** main :710.

### `137.2`, `49.4`, `2.4` seconds and `199`, `169`, `80` complexes

- **Is** median charged seconds per **qualifying pose**, where qualifying means
  PoseBusters-valid and within 2 Å, and the count of complexes contributing at least one
  such pose. AutoDock 137.2 s over 199 complexes, DiffDock 49.4 s over 169, EquiBind
  2.4 s over 80.
- **Source** main :737–739, :848.
- **Worth noticing** these three counts are exactly the best-of-top-30 recovery counts
  from the accuracy table, 199 / 169 / 80 of 303. The cost denominator and the deepest
  accuracy pool are the same quantity, which is why the cost slide and the RQ1 slide can
  be read against each other.
- **Reading** on 53 complexes where every tool succeeds the ordering persists
  (Friedman χ² = 85.2, p = 3.2e-19, Kendall W = 0.80) and every pair clears Holm.
  Median per-complex cost ratios are 0.14 for DiffDock against AutoDock and 0.02 for
  EquiBind.

### `223`

- **Is** the complexes where EquiBind produces no qualifying pose: 303 − 80. Its 2.4 s
  median describes only the 80 successes.
- **Trap — concede immediately** cost conditional on success favours whoever fails most.
  The thesis states it. Alternative denominators are reported and none changes the
  ordering: per attempted docking the medians are 157.7 / 227.6 / 8.1 s; allocating the
  whole campaign to successful complexes gives 271 / 479 / 38 s, on which AutoDock is
  actually *cheaper* than DiffDock in charged and GPU terms.

### `21.7` seconds

- **Is** the median charged cost per qualifying pose for **Vina search alone**, with no
  GPU time, over 202 complexes and 302 qualifying poses. Below DiffDock on the charged
  basis.
- **Trap** Vina alone nevertheless has the highest CPU cost per qualifying pose,
  1,805.0 s against the rescored arm's 1,752.8 s, because rescoring adds nine qualifying
  poses at no CPU cost. Its 32 threads compress 151.42 CPU-core-hours into 4.73
  wall-hours. AutoDock's advantage is coverage and admissibility, not economy.

### `5.3`, `5.3`, `137.2` — the endpoint-tightening series

- Worth having ready: median charged seconds per **generated**, **valid**, and
  **valid-near-native** pose. AutoDock 5.3 / 5.3 / 137.2. DiffDock 7.6 / 9.2 / 49.4.
  EquiBind 0.3 / 0.4 / 2.4. Contributing complexes fall 303 → 301 → 199, 303 → 300 → 169,
  303 → 266 → 80. Tightening from valid to valid-near-native multiplies AutoDock's unit
  cost by 25.9, against 5.4 and 6.0.
- **Trap** these are charged seconds, **not** wall seconds. The wall figures are
  1.4 / 1.4 / 26.6 and quoting those instead breaks the RQ3 argument.

### `a quarter of the time`

- EquiBind succeeds on 80 of 303 complexes, which is 26.4 %. The note's "a quarter" is
  that. Its charged figure is also an idealised floor, because the uniform 32-thread
  basis credits it with concurrency it reaches only by running separate complexes in
  parallel. Its 21.39 CPU-core-hours are under a seventh of AutoDock's 151.42.

---

## Slide 12 — RQ1 and RQ2

### `98.9 %`, `84.5 %`, `62.0 %`

- **Is** pooled physical validity of the three carried-forward arms: AutoDock + gnina,
  DiffDock + smina, EquiBind + gnina.
- **Over** poses, not complexes.
- **Trap** DiffDock reaches 87.1 % with gnina rather than smina, and EquiBind only
  51.3 % with smina. The reported figures are variant-specific.

### `36.6 / 34.0 / 18.2` and `65.3 / 55.1 / 26.4`

- **Is** rank-1 and best-of-top-15 recovery on the headline endpoint, valid and within
  2 Å. In counts, rank-1 gives 111, 103 and 55 of 303; top-15 gives 198, 167 and 80;
  top-30 gives 199, 169 and 80.
- **Depth gains** AutoDock +28.7 points (Wilson 23.9–34.0), DiffDock +21.1 (16.9–26.1),
  EquiBind +8.3 (5.7–11.9). AutoDock adds 87 complexes, DiffDock 64.
- **Separation** AutoDock over DiffDock is +13.2 points at top-5 (p_holm = 7.6e-4),
  +10.2 at top-15 (p_holm = 0.009), +9.9 at top-30 (p_holm = 0.011). Not at rank-1.
- **Source** main :382, :388, :834.
- **Trap** EquiBind's small depth gain is a **sampling** constraint, not a ranking one.
  A single-shot regression model provides few distinct alternatives and rescoring cannot
  recover poses that were never generated.

### `24.4 → 84.5 %` and `2.9 → 62.0 %`

- **Is** the RQ2 repair effect on pooled validity. Rank-biserial effects 0.994 to 1.000,
  adjusted p below 1e-35. These are the largest effects in the thesis.

### `31.0 → 36.6 %`

- **Is** what gnina re-ranking alone buys AutoDock at rank-1, and 61.1 → 65.3 % at
  top-15. It reorders the pool rather than enlarging it, so the gain is confined to the
  top of the list.
- **Trap** over the whole thirty-pose pool the **raw** search is nominally ahead, 202
  complexes against 199, because rescoring minimises every pose and a few cross the
  threshold the wrong way. And that single step does **not** clear Holm correction at
  rank-1 on its own; it reaches significance only from top-5 onwards.

### `0.1` to `4.6` points

- **Is** the change in fixed-rank near-native rates from optimisation, spanning all
  arms: 0.1 to 0.7 points for DiffDock, 4.1 to 4.6 for EquiBind with gnina, 3.5 to 4.4
  with smina.
- **Reading** this is the mechanism claim. Optimisation repairs geometry and re-orders a
  pool. At a *fixed* rank it barely moves accuracy, because it cannot create a binding
  mode that was never sampled. The larger top-k gains accumulate these small
  per-rank changes across ranks.

### `top-5 onwards`

- The one number to protect. Everything separable in this thesis becomes separable at
  top-5, not at rank-1.

---

## Slide 13 — RQ3 and RQ4

No new figures in the notes. The claims rest on slide 11's numbers.

### `three` failure modes

- AutoDock fails by mis-ranking a pool that contains the answer. DiffDock fails by not
  localising. EquiBind fails by not generating anything usable. These are the RQ4 answer
  and they are what makes the conclusion a workflow rather than a winner.
- **Supporting figure** about six in ten rank-1 failures contain no qualifying pose
  anywhere in the pool, so generation and localisation remain the larger constraint.
  Source main :406.

---

## Slide 14 — Conclusions

No figures in the notes. The practical finding to have in one sentence: one cheap local
minimisation pass converts most learned output from physically impossible to admissible,
at a cost negligible next to the generation step. The supporting numbers are slide 5's
13.3 → 93.3 % and 0 → 80.0 % medians, against slide 11's per-pose seconds.

---

## Slide 15 — Limits and Outlook

### `303 complexes`

- The comparison's outer bound. Variant selection and evaluation use the same 303, with
  no held-out split.

### `three ligands and four correlated frames`

- The Orai1 chapter's outer bound. One trajectory, four snapshots, three modulators.

### `unseeded` DiffDock

- The DiffDock benchmark run is a single unseeded draw. Run-to-run variance is not
  characterised.

### `5 Å` and `78 of 303`

- **Is** the metal stratum. A metal ion or metal-containing cofactor lies within 5 Å of
  the crystal ligand in 78 of the 303 analysed complexes, and receptor preparation
  stripped waters, ions, metals and cofactors for the AutoDock and DiffDock arms. A
  quarter of the set is therefore scored without a coordinating partner the native pose
  depends on.
- **Source** main :88, :868, appendix :165.

### `+6.7` against `+20.5` points

- **Is** the AutoDock-over-DiffDock top-15 contrast split by that stratum. On the 225
  metal-free complexes it is +6.7 points (134 against 149 of 225, exact McNemar
  p = 0.142). On the 78 metal-adjacent ones it is +20.5 points (33 against 49 of 78,
  p = 0.020).
- **Reading** the pooled contrast averages a narrow gap with a much wider one. DiffDock
  is disproportionately hurt by the missing metal, which is consistent with a learned
  model trained on holo structures.
- **Trap** the split was made **after the fact**. Present it as a sensitivity analysis,
  not as a prespecified subgroup finding.

### the halogen blind spot

- The interaction profiler cannot see chlorine or bromine in its halogen-bond test, so
  any halogen-bond result is discounted. Not a number, but the same class of concession.

---

## Cross-slide number collisions to guard against

| Value | Meaning A | Meaning B |
|---|---|---|
| `28.7` | AutoDock rank-1 → top-15 gain, percentage points (slide 12) | Orai1 benchmark-to-experimental usable-yield drop (slide 9) |
| `2 Å` | in-place near-nativeness threshold | Kabsch form threshold at the loose end of the ten-threshold table |
| `5 Å` | inter-tool consensus-site agreement radius (slide 10) | metal-to-ligand proximity cut (slide 15) |
| `303` | analysed benchmark complexes | never 308, which is the prepared / Orai-control ligand count |
| `0.4–0.8 pp` | winner's curse over three refiner arms per learned family | **not** the AutoDock correction, which is ~2.0 pp at rank-1 |
| `95 %` | cluster-level 4 Å oracle (slide 8) | **not** a confidence interval, and **not** pose accuracy |
| `5.3 s` | charged seconds per generated pose | **not** wall seconds, which are 1.4 s |

## Three numbers to keep out of your mouth unless asked

1. **97.7 %** pooled full-cloud oracle. It is true and it sounds like an accuracy claim.
2. **0.74 wall-hours** for the gnina pass. It is the accounting you rejected.
3. **75 %** AutoDock on the Vina prefix. Correct, but it is the number that removes your
   own Orai1 headline. Volunteer it only in the caveat, where it reads as rigour.
