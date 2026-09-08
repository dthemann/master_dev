# Examiner Report — MSc Thesis in Computational Chemistry / Molecular Modelling

**Comparative evaluation of AutoDock Vina, DiffDock and EquiBind with post-docking optimisation, on a crystal-referenced calibration benchmark and a reference-free Orai1 application**

Date of report: 2026-09-06
Files examined: `thesis_latex/Thesis_short.tex`, `thesis_latex/body_main_short.tex`, `thesis_latex/body_appendix_short.tex`, `thesis_latex/Literatur.bib`, `thesis_latex/Thesis_short.pdf` (with `.toc`, `.lot`, `.lof`, `.log`), and the supporting repository at `/home/manndo/master_dev`.

---

## Methodology note

This review was produced by a multi-agent examination process rather than a single reading. Each rubric sub-category (research question, literature foundation, methodological design, the four docking-protocol audits, the two validation axes, results and analysis, interpretation, reproducibility, presentation, and the docking red-flag checklist) was assessed by a dedicated primary examiner. Each primary assessment was then attacked by two independent adversarial challengers instructed to refute findings, correct misreadings of the thesis, and surface anything the primary missed. An adjudicator reconciled the three inputs per category, re-verifying every load-bearing quote and, where possible, re-deriving numbers directly from the thesis files and from the supporting repository. Every surviving finding then went through a second, independent validation pass with two verifiers, each of whom could confirm, weaken, downgrade or refute it.

The process was deliberately conservative. **Four findings were refuted outright and removed** at the validation stage, and a fifth was reduced to nil; a substantial further number were downgraded in severity or reclassified after the thesis was shown to state, quantify or disclose the point the finding alleged was missing. Those removals are recorded in Section 9.4 so that the reader can see what was tested and rejected, not only what survived. Prior review documents present in the working directory were excluded from the evidence base to keep this assessment independent.

**Outcome: 79.3 / 100 — strong.** Four major findings, no critical findings.

---

## 1. Executive assessment

The thesis asks four numbered research questions (`body_main_short.tex:11-14`) about physical validity, post-docking optimisation, computational cost and failure modes for three configured docking workflows — AutoDock Vina, DiffDock and EquiBind, each with optional smina/gnina refinement — and answers them on two arms: a 303-complex crystal-referenced calibration benchmark under blind whole-protein search, and a reference-free application to four molecular-dynamics frames of a human Orai1 homology model with three published modulators.

The scientific approach is sound and unusually well controlled for an MSc. Site information is withheld uniformly from all three headline pipelines (`body_appendix_short.tex:148`), recovery is resolved by ranking depth rather than reported from one selected pose, physical validity is screened with a twenty-two-check PoseBusters battery that carries no crystal reference, and the reference-free arm is deliberately denied the calibration arm's warrant. The principal contribution is not a new method but a carefully instrumented, depth-resolved comparison of deployed workflows with an honest transfer test to a target that has no ligand-bound structure.

The strongest aspect is the thesis's willingness to attack its own results and report the damage: it shows that its most eye-catching Orai1 finding is a ranker artefact (`body_main_short.tex:836`), bounds its own selection optimism (`:866`), and concedes that part of its validity advantage is constitutive (`body_appendix_short.tex:1188`). The most important weakness is a cluster of internal inconsistencies concentrated in the Orai1 chapter and its receptor appendix, where the documented receptor is not the receptor that was docked.

Conclusions are broadly supported but overstated in a small number of load-bearing sentences. Overall quality: strong, upper second-class-plus to first-class band for an MSc. Confidence: high — every headline number checked reproduced, and the supporting assertion harness runs clean.

---

## 2. Summary of the research approach

**Design.** Two arms. The *calibration* arm uses 308 post-2021 protein–ligand crystal complexes from the PoseBusters benchmark (`body_main_short.tex:88`), of which 303 are analysed after five DiffDock production failures, handled with an intention-to-treat sensitivity check (`body_appendix_short.tex:1072`). The *application* arm docks three Orai1 modulators (2-APB, Synta-66, GSK-7975A) into four MD frames of a homology model, alongside a 1,232-unit background control panel built from the calibration ligands (`body_main_short.tex:541`).

**Pipelines.** AutoDock Vina (locally patched build `v1.2.7-20-g93cdc3d-mod`), DiffDock-L, and a custom unguided EquiBind harness, each optionally followed by smina or gnina local optimisation. All three search the whole receptor: "AutoDock Vina, DiffDock and unguided EquiBind all search the whole receptor without a crystal-derived restriction, so none is told where the ligand belongs" (`body_appendix_short.tex:148`). Six pocket-guided EquiBind arms exist but are excluded from the three-tool comparison on information-parity grounds and reported as a negative result (`body_appendix_short.tex:265`).

**Endpoints.** Near-nativeness (symmetry-corrected in-place RMSD ≤ 2 Å, `body_main_short.tex:107-113`); Form (Kabsch RMSD ≤ 1 Å, `:122`); PoseBusters validity over twenty-two applicable reference-free checks (`:138`); operational placement for Orai1 via an Arg91–Glu106 transmembrane exclusion slab (`:145-149`); cross-tool clustering with a prespecified 4 Å crystal-site threshold (`:152-159`); PandaMap interaction fingerprints (`:161`); and computational cost on a charged CPU/GPU basis (`:168`).

**Statistics.** Exact McNemar on paired binary outcomes, Cochran's Q, Friedman/Wilcoxon/Mann–Whitney with matched effect sizes, Holm and Benjamini–Hochberg families declared a priori, 2,000-resample bootstraps and 10,000-permutation nulls with fixed seeds (`body_appendix_short.tex:1096-1120`, `:1214`). The Orai1 inferential unit is fixed as the frame–ligand pair and honestly deflated: "The effective number of independent experimental observations is therefore three ligands" (`body_main_short.tex:181`).

**Scope declaration.** "Each question is posed about the three workflows as configured here rather than about physics-based and AI-based docking as method classes" (`body_main_short.tex:16`), repeated at `:851` and in the Abstract (`Thesis_short.tex:65`).

---

## 3. Five principal strengths

1. **Self-refutation of the study's own headline application result.** The Orai1 between-tool reversal is dismantled by the thesis itself: "AutoDock is capped on its gnina re-ranking, and restricting it instead to Vina modes one to ten puts it at 75\% and removes the reversal" and "The five-point deficit is therefore a property of the secondary ranker rather than of the search" (`body_main_short.tex:836`). The mechanism is then measured independently on 12,315 control poses — Spearman −0.211 under the gnina order against −0.004 under the Vina order (`body_appendix_short.tex:369`) — and carried into the Results, moving the adjusted p from 5.6e-4 to 0.054 (`body_main_short.tex:625`).

2. **Protocol reporting at publication rather than MSc standard.** Every binary is version-pinned by invocation, including the patched Vina build identified by upstream commit (`body_appendix_short.tex:152`); both local source patches are disclosed with impact audits (`:155-158`); and five preparation choices that cannot be recovered from output files — protonation, Gasteiger charges, occupancy, alternate locations, symmetry mates — are recorded explicitly (`:150`). Alternate-location and water handling reproduce exactly against the shipped files (zero altloc indicators and zero HOH records across all 428 receptors).

3. **Exhaustiveness treated as a decision, not a default.** An eight-arm, five-rung ladder is tabulated at four ranking depths with exact McNemar tests, Holm correction, measured search hours and a bootstrapped marginal-return interval of "+65.6 per hour with a 95\% interval of +32.1 to +137.2" (`body_appendix_short.tex:232`), then discounted for its own optimism: "Choosing the best of eight configurations on the same 303 complexes that then report its performance carries a winner's curse" (`:248`).

4. **Instrument self-criticism against the author's own interest.** The thesis concedes that "Part of the physics-based advantage reported in the Results is constitutive rather than an independent verdict" (`body_appendix_short.tex:1188`), then bounds it: requiring validity "removes at most four complexes of 303 at any reported depth" (`:1195`). It also documents three PandaMap execution-path defects with their direction of bias (`:545-549`), including halogen blindness that "inflates the halogen-bond contrast reported in the Results in the same direction that contrast runs".

5. **Verifiable end-to-end reproducibility.** The assertion harness `Scripts/Analysis/thesis_assertions.py` was executed during this review and returned "1364/1364 checks reproduce. Every number checked matches the thesis", "Figures 39 included, 39 accounted for", exit 0. Independent regeneration of prepared ligand PDBQT files from the appendix's own printed commands produced byte-identical output for the three test cases, and the RMSD chain is cross-validated against the PoseBusters implementation, flipping only two of 909 rank-1 decisions (`body_appendix_short.tex:1385`).

---

## 4. Five principal weaknesses

1. **The documented Orai1 Fr0 receptor is not the receptor that was docked.** `body_appendix_short.tex:560` states the Fr0 AutoDock search "ran against 1,333 residues and 10,321 heavy atoms"; every staged Fr0 file in both reported panels carries 1,338 residues and 10,368 heavy atoms, and the quoted figures reproduce a retired Meeko-prepared file the same appendix disowns at `:152`. Separately, the Fr0 AutoDock receptor was PDBFixer-repaired and OpenMM-minimised (0.819 Å heavy-atom, 0.636 Å Cα from the delivered frame) while DiffDock and EquiBind read the delivered coordinates at 0.000 Å.

2. **A blind whole-protein box makes an inherited scoring convention into a live, unquantified confound.** 143 of 308 distributed entries carry more than one crystallographic copy of the reference ligand, all inside the searched volume; scoring is against "the single deposited reference instance rather than against the closest copy" (`body_appendix_short.tex:664`), with no count and no sensitivity arm.

3. **A Results recommendation states the wrong sign for one of the thesis's own measurements.** "Refined-score ordering can be considered an unclaimed remedy for DiffDock, worth about three points of rank-1 recovery on the reported coordinates" (`body_main_short.tex:406`), against the thesis's own two measurements at `:197` and `body_appendix_short.tex:962`, where re-ordering "does not help" (32.7% and 32.3% against 33.7%).

4. **The Orai1 chapter carries several internal tensions on exactly the axis its design is meant to isolate.** The exclusion slab discards the lipid-facing interface the four-frame design was built to sample (`body_appendix_short.tex:323` against `:581`); a Conclusions-level neutrality claim about DiffDock's pose cap (`body_main_short.tex:836`) is declared unmeasurable at `body_appendix_short.tex:365`; a declared "controlled contrast" between Fr0 and the thermal frames (`:601`, `:611`) is never reported; and two enumerations of residual panel differences (`:361`, `:363`) omit the start-conformer non-parity their own table prints (`:389`).

5. **The Literature Review chapter is a tool description rather than a foundation.** It cites no empirical evaluation study at all — including `ref035`, which supplies the dataset and the primary endpoint — and the thesis nowhere introduces molecular recognition, binding thermodynamics or the scoring-function and search taxonomies its own contrast depends on. The student's genuine synthesis exists, but it is filed in the Motivation, Appendix B and the Discussion.

---

## 5. Detailed scoring table

Scores below are the panel's adjudicated values and are reproduced exactly; contributions are the panel's stated weighted contributions.

| # | Category | Score /5 | Sub-scores | Weight | Contribution | Confidence |
|---|---|---|---|---|---|---|
| 1 | Research question and significance | 4.00 | B = 4 | 10% | 8.00 | high |
| 2 | Literature and theoretical foundation | 3.50 | C = 3.5 | 10% | 7.00 | high |
| 3 | Overall methodological design | 4.00 | D = 4 | 10% | 8.00 | high |
| 4 | Docking protocol and technical execution | 3.88 | E1 = 3.5, E2 = 3.5, E3 = 4, E4 = 4.5 | 20% | 15.50 | high |
| 5 | Validation and benchmarking | 4.25 | F1 = 4, F2 = 4.5 | 15% | 12.75 | high |
| 6 | Results and data analysis | 4.00 | G = 4 | 10% | 8.00 | high |
| 7 | Interpretation, discussion and limitations | 4.00 | H = 4 | 10% | 8.00 | high |
| 8 | Reproducibility and transparency | 4.00 | I = 4 | 10% | 8.00 | high |
| 9 | Writing and presentation | 4.00 | J = 4 | 5% | 4.00 | high |

**Weight check:** 10 + 10 + 10 + 20 + 15 + 10 + 10 + 10 + 5 = **100%** ✓

**Arithmetic:** 8.00 + 7.00 + 8.00 + 15.50 + 12.75 + 8.00 + 8.00 + 8.00 + 4.00 = **79.25**, reported as **79.3 / 100**.

A supplementary docking red-flag category (K) was scored at 4.0 with high confidence. It carries no rubric weight and is reported in full at Section 8.

### 5.1 Per-category justification, evidence and key deficiency

**1. Research question and significance — 4.0 (contribution 8.0), high confidence.**
Four numbered research questions each receive a dedicated Conclusions section (`:832`, `:841`, `:846`, `:853`), and every headline number traces to a Results float (36.6/34.0/18.2 rank-1; 65.3/55.1/26.4 top-15; 137.2/49.4/2.4 s charged). *Strongest evidence:* `body_main_short.tex:16` — "Each question is posed about the three workflows as configured here rather than about physics-based and AI-based docking as method classes" — a restriction honoured at `:851` and in the Abstract. *Key deficiency:* the reference-free Orai1 arm rests on two claims the thesis's own pages undercut — a radially unbounded usable-yield endpoint (`body_appendix_short.tex:323`) that excludes a region App. D:581 names as "where small-molecule Orai modulators are proposed to bind", and a Synta-66 concordance (`:695`) never discounted against the decoy control ten lines earlier (`:685`, "No residue differs after Holm correction").

**2. Literature and theoretical foundation — 3.5 (contribution 7.0), high confidence.**
Split foundation. Tool scholarship is at the top of the MSc range: `body_appendix_short.tex:66` and `:77` give the complete Vina functional form, the six fitted weights and the full Vinardo diff, with a footnote recording that the constants "were read from the source of the patched Vina build ... rather than transcribed from the publication". *Key deficiency:* the chapter that carries this category's name is a tool description. `grep -ci posebusters` over lines 20–73 returns 0; `enthalp` returns 0 across all three files; "van der Waals", "knowledge-based" and "Monte Carlo" return 0 in the main body. The chapter is 1,151 words against an 18,539-word body.

**3. Overall methodological design — 4.0 (contribution 8.0), high confidence.**
The inferential unit is prespecified and deflated (`body_main_short.tex:181`), the application arm is denied the calibration arm's warrant (`body_appendix_short.tex:864`), a 1,232-unit background panel is correctly disclaimed (`:541`), and winner's curse is measured rather than confessed (`:866`). *Key deficiency:* the Orai1 endpoint and its receptor-ensemble rationale are calibrated against different biological hypotheses and never reconciled, compounded by a Conclusions neutrality claim the appendix declares unmeasurable, a panel-parity enumeration that omits a non-parity its own table prints, and a "controlled contrast" (`:601`) no result reports.

**4. Docking protocol and technical execution — 3.88 (contribution 15.5), high confidence.**
E1 protein preparation 3.5; E2 ligand preparation 3.5; E3 binding-site definition 4.0; E4 engines, settings and pose selection 4.5. *Strongest evidence:* `body_appendix_short.tex:152` version pinning verified by invoking each binary, plus the eight-variant exhaustiveness ladder with Holm-adjusted McNemar tests. *Key deficiency:* Orai1 receptor preparation is not internally consistent and not parity-matched on Fr0; the searched, rescored and analysed AutoDock ligand are not the same molecule; and the whole-protein box places genuine alternate copies of the correct site inside the searched volume for 143 of 308 entries.

**5. Validation and benchmarking — 4.25 (contribution 12.75), high confidence.**
F1 pose-prediction validation 4.0; F2 ranking, affinity and statistics 4.5. *Strongest evidence:* the RMSD chain is validated against an independent implementation and the disagreement counted — `body_appendix_short.tex:1385`, "Substituting the PoseBusters implementation into the 2\angstrom{} rank-1 gate moves two of the 909 tool-by-complex decisions". *Key deficiency:* the thesis never establishes that a pose it counts as near-native reproduces the native interactions, in either the visual or quantitative register; and the primary endpoint carries Holm-corrected p-values whose family Appendix G declares excluded (`:869`).

**6. Results and data analysis — 4.0 (contribution 8.0), high confidence.**
The cost chain closes end to end and independently of the prose (151.42 CPU-core-hours / 32 = 4.73 wall-hours; +10.27 GPU-hours = the printed 15.00 charged hours, 68.5% for rescoring), and the depth ladder reproduces. *Key deficiency:* the Orai1 interaction analysis is over-read against an under-specified instrument; two of fourteen BH-clearing cells are structurally impossible for the three modulators under the executed atom-typing rule yet are given a chemotype explanation (`body_main_short.tex:683`).

**7. Interpretation, discussion and limitations — 4.0 (contribution 8.0), high confidence.**
Scope is restricted (`:770`, `:774`), nulls are never converted into equivalence (`:794`), no docking score is presented as an affinity (`:870`), and limitations are quantified rather than gestured at (`:866`, `:868`). *Key deficiency:* two load-bearing sentences fail — the DiffDock cap neutrality claim (`:836` against `body_appendix_short.tex:365`) and the filter-residue recovery claim at `:810` against `:694`. The panel's third item under this heading, that the physics-versus-learning qualification never states the AutoDock rank-1 lead is produced by a convolutional rescorer, was refuted at validation and is recorded in Section 9.4 as H-3; `body_main_short.tex:384` states it.

**8. Reproducibility and transparency — 4.0 (contribution 8.0), high confidence.**
1364/1364 asserted numbers reproduce; 39/39 figures accounted for; every version pin at `body_appendix_short.tex:152` matches `Conda_Env/vina.yml`; `Scripts/Analysis/REGENERATE.md` (1,177 lines) maps each float to its command and determinism class. *Key deficiency:* Appendix D, the section whose purpose is to make the Orai1 receptor reproducible, describes the wrong file.

**9. Writing and presentation — 4.0 (contribution 4.0), high confidence.**
Abstract and Kurzfassung are numerically faithful across twelve checked figures; the build carries zero undefined references and one oversized-float warning. *Key deficiency:* Table 1's Near-Nativeness column reads backwards against the thesis's own conclusion with an undefined "Comp." header and an empty footnote row, and "benchmark" names more than one set across the Results chapter with the bridging sentence only in the appendix.

---

## 6. Docking-protocol audit

### 6.1 Protein preparation (E1 — 3.5 / 5)

**Adequately reported.** Alternate-location, occupancy and biological-assembly handling is disclosed to a level almost never seen at MSc and reproduces exactly: `body_appendix_short.tex:150` records that no atom record in any of the 308 benchmark receptors or four Orai1 frames carries an alternate-location indicator, that fractional occupancy survives as a label on 2.1% of atoms across 258 files and is never weighted on, and that entries were docked as supplied so no biological assembly was generated. Water is stated correctly rather than assumed: "Crystallographic water is absent from the distribution itself". The protonation model actually produced is reported residue class by residue class with counts (`:163`: histidine doubly protonated in all 4,485 intact rings; 316 of 2,282 cysteines carrying an inherited thiol hydrogen; `prepare_receptor` "accepts no pH argument and performs no pK_a estimate"). Metal and cofactor stripping is quantified, stratified on 78 of 303 complexes, tested at 4 and 6 Å, and carried into Limitations (`body_main_short.tex:868`). Software is named with versions verified by invoking each binary (`body_appendix_short.tex:152`). Residue renumbering by the whole-protein preparation is anticipated and neutralised (`:551`, `body_main_short.tex:163`).

**Inadequately justified.** Metals, ions and cofactors are stripped for two of three arms with no stated scientific rationale — `body_appendix_short.tex:163` reports the operation as a pipeline fact, and `:533` concedes "This is nonetheless a difference from the source benchmark, which ships receptors retaining their cofactors". The blanket cationic histidine model is likewise stated but not defended, and on the Orai1 frames it discards force-field-assigned per-residue states (33 HSP, 23 HSE, 4 HSD), altering 27 of 60 residues including His113 in two of six subunits. Terminal treatment of the truncated Orai1 66–288 construct is not reported.

**Entirely missing.** Missing heavy atoms, incomplete side chains and chain breaks are neither repaired nor reported for any benchmark receptor (measured 1.20% of residues over 217 files). Receptor preparation for the DiffDock arm is never named on either dataset. No stereochemical or model-quality assessment of the Orai1 homology model is reported. The absence of calcium from the Orai1 receptor is never stated, although one of the two candidate sites is defined by its calcium-binding function. And, most seriously, the Fr0 AutoDock search receptor's PDBFixer repair plus amber14-all/GBn2 minimisation is absent from Appendix D and contradicted by two statements there and at `body_main_short.tex:97`.

### 6.2 Ligand preparation (E2 — 3.5 / 5)

**Adequately reported.** Disclosure is verification-grade and reproduces independently: `body_appendix_short.tex:150` states "It adds 353 such atoms across 161 of the 308 benchmark ligands, 352 of them on nitrogen", and comparing each `*_ligand_start_conf.sdf` with its staged PDBQT gives exactly 353 added polar hydrogens on exactly 161 ligands. The net-formal-charge census (299/7/1/1), the InChIKey verification of all three Orai modulators against published identifiers (`:793`), the 30-versus-1 conformer asymmetry (`:148`), the boron atom-type patch down to its radius and well depth (`:154`), and the PoseBusters `energy_ratio` patch (`:158`) are all recorded. Conformer-generation difficulty is measured rather than asserted (median 1.56 Å best-fit start-conformer-to-crystal RMSD, `:782`). Stereochemistry is handled: the reference-based identity checks were measured post hoc and cost DiffDock one complex at rank-1, with "AutoDock and EquiBind remain unchanged" (`body_main_short.tex:140`).

**Inadequately justified.** One protomer per ligand, states inherited from the source benchmark, and no pK_a model at any stage is a proportionate MSc choice applied identically to all three engines — but the consequence is never sized. No ionisable-group census accompanies the net-charge counts, ligand and receptor sit on inconsistent protonation models, and ligand preparation is absent from a Limitations chapter that enumerates five other parity axes.

**Entirely missing / misstated.** The AutoDock arm searches an N-protonated ligand and analyses a neutral one, and the write-up states the opposite in three places ("the docked poses carry them", `:150`; "The three compounds were docked as neutral species", `body_main_short.tex:97`; "as Docked and Analysed", `body_appendix_short.tex:807`). The boron reassurance at `:156` — "Every docked 2-APB file in the study records the atom types A, C, HD, NA and OA and not one records B" — is false for the reported Orai1 arm, where all four docked 2-APB files record thirty B atoms. Macrocyclic ligands are never mentioned although eight are present with frozen ring conformations. An anionic GSK-7975A arm was docked on all four frames and is silently absent from every result.

### 6.3 Binding-site and search-space definition (E3 — 4.0 / 5)

**Adequately reported.** The box is a deterministic, reproducible rule: "Its centre is the midpoint of the axis-aligned bounding box over every atom record of the prepared receptor, and each edge is that axis's span plus 1\angstrom{} on either side. No further padding is applied" (`body_appendix_short.tex:167`), with measured edges 29.5–167.4 Å, volumes 44,430–2,479,703 Å³, median 263,874 Å³ against the source benchmark's 15,625 Å³ crystal-centred cube, and an explicit warning at `:1199` that "a single edge figure must not be cubed". Site information is withheld uniformly (`:148`); no crystal information enters the refiners (`:169`, `:253`) or the validity screen (`body_main_short.tex:140`). The Orai1 placement rule is given as explicit geometry with its own limits stated (`body_appendix_short.tex:311-325`).

**Inadequately justified.** The exclusion slab's "membrane band" is asserted rather than determined — no bilayer position was computed, and the 22.4–22.8 Å slab is narrower than a bilayer while being radially unbounded to the 32–36 Å outer shell. The usable criterion is purely exclusionary, and the stricter atom-level rule the appendix says was recorded (`:321`) is never exercised. The threshold of the protein–ligand maximum-distance check is never quoted.

**Entirely missing.** No radially bounded or lumen-restricted sensitivity of the usable-yield endpoint. No quantification of the multiple-copy confound (143 of 308 entries). No stratification of the benchmark comparison by box volume despite a 56-fold range at fixed exhaustiveness. No internal site-boxed or predicted-pocket Vina arm — the thesis says so itself: "A run against a predicted pocket was not performed, and the comparison with a site-specific protocol therefore remains inferred from the literature rather than established here" (`:1207`). No sampling-adequacy check on the Orai1 boxes, from which the calibration ladder result is explicitly withdrawn (`:1205`).

### 6.4 Docking settings and pose selection (E4 — 4.5 / 5)

**Adequately reported.** This is the strongest sub-category. All binaries version-pinned by invocation including a Vina build by upstream commit; both local patches disclosed with impact audits showing neither moves a reported number; exhaustiveness earned through a five-rung, eight-variant ladder with exact McNemar/Holm tests, measured search hours, a bootstrap marginal-return interval and an explicit winner's-curse correction (`body_appendix_short.tex:206-250`); the gnina ranking head subjected to a sensitivity analysis that runs against the author's interest (CNNaffinity 111, CNNscore 112, minimised Vina 90 at rank-1, `:173`); pose selection fully automated and depth-resolved with an oracle ceiling and three benchmarked consensus rules (`body_main_short.tex:474`); a genuine refiner non-parity found and repaired by re-running under `--minimize` (`body_appendix_short.tex:276`).

**Inadequately justified.** DiffDock's twenty denoising steps are stated without citation, default declaration or sensitivity check. The 6 kcal/mol energy window is unjustified. Exhaustiveness 128 was transferred unchanged to the Orai1 boxes with no convergence check. The primary DiffDock variant was chosen on a 133-versus-132 margin the thesis itself calls unstable and then carried into the reference-free arm.

**Entirely missing.** No replication at any seed for any stochastic arm on either dataset; the 303-complex DiffDock calibration run is "one unseeded draw" (`body_main_short.tex:864`) and the Orai1 control panel records "Ctrl. none recorded" (`body_appendix_short.tex:394`). The Vina pose-redundancy threshold governing which of thirty modes are distinct is never reported. The exhaustiveness-92 rung was never rescored, leaving the selection grid incomplete at its cheapest cell. Version identifiers are absent for the two learned engines and the DiffDock checkpoint.

---

## 7. Validation assessment

### 7.1 Pose-prediction validation (F1 — 4.0 / 5)

**Grade: unusually strong for an MSc, with one conspicuous hole.**

The reference is genuinely held out — the benchmark supplies receptor, all crystallographic instances, the chosen ground-truth instance and a separate RDKit ETKDGv3 + UFF start conformer (`body_appendix_short.tex:662`). Symmetry handling is explicit, its one limitation named, and its bias shown to be conservative (`body_main_short.tex:113`). The arbitrary 1 Å Kabsch bound is declared arbitrary and recovery resolved over ten thresholds from 0.25 to 2.5 Å at three depths (`body_appendix_short.tex:1261-1340`). Sampling is separated from ranking with a formal rank-1-to-top-15 decomposition (`:1224`). The exclusion of five complexes carries an intention-to-treat check that moves the top-15 AutoDock margin only from 10.2 to 10.1 points (`:1072`). The validity half of the conjoint endpoint costs at most four complexes of 303 at any reported depth (`:1195`).

The hole is that the thesis never establishes that a pose it counts as near-native reproduces the native contacts. The contact analysis is computed on a pool with median 7.57 Å deviation from the crystal ligand and 21.3% of poses within 2 Å (`:543`), and the one distance-conditioned sensitivity deliberately stops short of the gate. No figure anywhere superposes a predicted pose on its crystal reference. There is also no positive control screening a deposited crystal pose through the twenty-two reference-free checks, which is the only non-circular calibration available given the thesis's own concession at `:1188`.

### 7.2 Ranking and affinity claims (F2 — 4.5 / 5)

**Grade: unusually strong for an MSc.**

Rank ordering and absolute affinity are kept strictly apart: "No endpoint reported here resolves binding free energy. Recovery, physical validity and RMSD are geometric criteria, and no score used for ranking was compared against measured affinities" (`body_main_short.tex:870`). No docking score is presented anywhere as a measured binding energy. Virtual-screening over-claiming is refused by declaration: "The control panel is only a background comparator, and no enrichment claim is made" (`:541`); ROC-AUC, EF and BEDROC appear only in the literature comparison. Score-to-quality concordance is measured per complex with Kendall tau and a one-sample Wilcoxon, BH-corrected across nine cells (`body_appendix_short.tex:1403`), reporting AutoDock at tau −0.20 to −0.32 (q < 0.001) and DiffDock's native confidence at tau ≈ 0. Prospective ranking rules are benchmarked against an oracle ceiling and a naive prior rule with corrected paired tests (`body_main_short.tex:474`).

### 7.3 Statistics

Test-to-design mapping is tabulated and correct for every data type (`body_appendix_short.tex:1096-1108`). Multiplicity families are declared a priori with a stated non-pooling rationale (`:869`). The headline null is explained rather than hidden: 112 discordant complexes, minimum detectable difference 9.8 points, observed power 0.12, ~4,200 pairs required (`:1258`). Equivalence is claimed once and its post-hoc construction confessed: "No equivalence margin was prespecified for this work. The margins quoted here were located after the fact by scanning the four candidate values 5, 10, 12 and 15 percentage points" (`:1258`). Non-significance is never converted into equivalence (`body_main_short.tex:183`, `:637`). Resampling is fully specified with fixed seeds and a stated p-floor (`:1214`).

Residual statistical defects are three and all correctable: the primary endpoint carries Holm labels whose family Appendix G excludes (`:869`); two pose-level p-values are pseudoreplicated against the declared inferential unit (`:369`, `body_main_short.tex:627`); and `:1214` promises bootstrap intervals for effect sizes that no Cliff's delta or rank-biserial in the thesis carries.

### 7.4 Practical consequence — which claims the results support

The results **support**: (i) relative ranking of these three configured workflows on this benchmark under blind whole-protein search, at declared ranking depths; (ii) the claim that raw learned output is largely physically inadmissible and that local optimisation repairs geometry without recovering unsampled binding modes; (iii) a cost ordering on a stated charged basis, conditional on the optional gnina pass; (iv) pose *hypotheses* for Orai1, explicitly framed as such.

The results **do not support**: (i) any statement about physics-based versus AI-based docking as method classes — a restriction the thesis itself imposes; (ii) virtual-screening prioritisation or enrichment; (iii) any binding-affinity or thermodynamic claim; (iv) a validated Orai1 binding mode or binding site; (v) a resolved rank-1 difference between AutoDock and DiffDock, which the thesis correctly reports as underpowered.

---

## 8. Docking-specific red-flag checklist

| # | Flag | Status | Evidence |
|---|---|---|---|
| 1 | Docking scores treated as experimental binding energies | **absent** | `body_main_short.tex:870` "no score used for ranking was compared against measured affinities"; `body_appendix_short.tex:66` Vina is "a fast and well-behaved surrogate objective rather than a physical free-energy decomposition"; `body_main_short.tex:17` "None of the endpoints is a claim about thermodynamic or time-dependent stability, which docking does not test." The only loose phrasing is the generic field description at `:4` and `:32`, corrected at `:39`. |
| 2 | Small score differences presented as decisive | **partially present** | Two instances. `body_main_short.tex:197` "DiffDock native ranking outperforms the smina ranking" on 33.7% vs 32.7% (three complexes of 303) with no interval or test; DiffDock + smina carried on a 133-to-132 margin. More seriously `:406` asserts refined-score ordering is "worth about three points of rank-1 recovery" where the thesis measured minus one (`body_appendix_short.tex:962`). Mitigated by "the ordering is unstable" at `:197` and by the refusal to claim the exh92→128 step, "too small to separate from the noise" (`body_appendix_short.tex:249`). |
| 3 | Docking protocol not validated | **absent** | 308 crystal complexes as the calibration arm (`body_main_short.tex:88`); eight-arm exhaustiveness ladder with exact McNemar, Holm correction, measured hours and bootstrapped marginal return (`body_appendix_short.tex:213`); transfer limit stated at `body_main_short.tex:868`, "cognate self-docking is the easier of the two regimes, so the calibration figures bound the Orai1 expectation from above rather than predicting it." |
| 4 | Receptor structure selected without scientific rationale | **absent** | `body_appendix_short.tex:581` "A single static conformer would over-weight one arbitrary geometry of exactly the two regions the Orai1 chapter asks about"; `body_main_short.tex:95` gives the rigid-receptor reason for four states; measured frame geometry and the absence of any conformation-selection heuristic are both disclosed. |
| 5 | Protonation, tautomerism or stereochemistry ignored | **partially present** | `body_appendix_short.tex:150` "Every ligand entered preparation in the tautomer and the ionisation state its source file carried, and no protomer or tautomer enumeration and no pK_a model ran at any stage"; `:163` histidine cationic in all 4,485 intact rings against bare imidazole nitrogens on 3,898 of 4,844 rings in the EquiBind path, so the arms differ on receptor chemistry; 2-APB docked neutral despite an expected protonated fraction (`:653`). Documented in unusual detail but not modelled and not matched. Stereochemistry is handled (`body_main_short.tex:140`). |
| 6 | Conserved waters, metals, cofactors or ions removed without justification | **partially present** | `body_appendix_short.tex:163` "Receptors are cleaned with PDBFixer, which strips waters, ions, metals and cofactors" — stated with no scientific reason; `:533` concedes the difference from the source benchmark. Only the rationale is missing: consequences are quantified at `body_main_short.tex:88` ("apo-like rates"), in a full stratified appendix paragraph (`:165`) and in Limitations (`:868`, +6.7 against +20.5 points). |
| 7 | Binding site or grid insufficiently documented | **absent** | `body_appendix_short.tex:167` gives the box construction rule verbatim with edges 29.5–167.4 Å over 303 receptors; `:1197` adds the edge-ratio distribution and the anti-cubing warning; the Orai1 placement envelope and pore-axis construction are specified at `:311`. |
| 8 | Only the best-looking poses or compounds reported | **absent** | `body_main_short.tex:192` "Table~\ref{tab:results-pose-production} reports up to thirty poses per complex for every variant and thus covers the full pool rather than one selected pose"; failures counted at `:406`; pocket-guided EquiBind reported as worse than no guidance. One undisclosed exception is documentary only — the Form column silently drops 42 exploded DiffDock poses, which works against DiffDock. |
| 9 | Visual inspection used without transparent criteria | **absent** | Both PyMOL renders are captioned "is not the output of any analysis script, so it carries no measurement and none of the quantities in the surrounding text is read from it" (`body_appendix_short.tex:411`, `:567`). Every selection rule is numeric and prespecified. |
| 10 | A favourable pose treated as proof of binding | **absent** | `body_main_short.tex:810` "Neither partial contacts nor their combination validates a complete pose"; `:586` "The slab remains a modelling choice rather than proof that every excluded pose is biologically wrong"; `:695` "evidence from different tools cannot be combined to validate one binding mode." |
| 11 | Docking results treated as proof of biological activity or mechanism | **partially present** | `body_main_short.tex:695` offers His113, Tyr115, Pro201 and Leu202 as "characteristic of both tools' modulator poses", while `:685` reports the inactive control panel concentrating on the same residues with "No residue differs after Holm correction of the fifteen within-tool tests", and `:702` recommends prioritising poses that "agree with independent pharmacology". Heavily hedged at `:541`, `:695`, `:700` and `:810`, but the discriminative consequence of the control result is never drawn. |
| 12 | Experimental validation implied but not performed | **absent** | `body_main_short.tex:876` "The Orai1 poses remain computational hypotheses. Validation requires site-directed mutagenesis, competition or displacement assays and ideally a ligand-bound structure. This work continues at JKU. The reference-free arm contributes a documented, physically screened workflow for testable hypotheses, not a validated complex." |
| 13 | Ligand or receptor preparation not reproducible | **partially present** | The benchmark arm is fully reproducible with exact binaries and flags (`body_appendix_short.tex:163`). The Orai1 arm is not: "There is no record of the modelling program or the template Protein Data Bank identifiers used to build the homology model, of the force field name or version, of the lipid composition" (`:572`), and "The externally supplied optimised geometries lacked records of both the optimisation program and its level of theory" (`body_main_short.tex:97`). Both gaps are disclosed and both originate with the supplying group. |
| 14 | The same data used to tune and to validate the protocol | **present** | `body_main_short.tex:866` "Variant selection and evaluation use the same 303 complexes, with no held-out split", declared in the Abstract itself (`Thesis_short.tex:61`). Optimism bounded per family at 0.4–0.8 percentage points and at about two points at rank-1 for the eight-arm AutoDock field (`body_appendix_short.tex:250`), with the tuned-versus-untuned asymmetry flagged in the Results at `:384`. |
| 15 | Conclusions exceed what docking alone can establish | **absent** | `body_main_short.tex:17` scopes the endpoints before any result; the Conclusions restrict every claim to configured pipelines; Limitations at `:870` states "Recovery, physical validity and RMSD are geometric criteria"; the Abstract's Orai1 reversal is immediately attributed to a ranking artefact (`Thesis_short.tex:65`). |

**Summary: 10 absent, 4 partially present, 1 present, 0 impossible to determine.** The single fully present flag (#14) is declared by the thesis in its own Abstract and bounded numerically, which is the correct handling.

---

## 9. Issues

### 9.1 Critical issues

**None.** No finding survived validation at a severity that would invalidate a reported result or require the thesis to be withdrawn and re-run.

### 9.2 Major issues (4)

---

#### M1 — The Fr0 Orai1 AutoDock search receptor was repaired and energy-minimised, breaking receptor parity with DiffDock and EquiBind on that frame, while two other passages deny that any pH-based protonation occurred
*Finding E1-1 · scientific_flaw · CONFIRMED by both verifiers*

**Where.** `body_appendix_short.tex:558` (Appendix D, Orai1 Receptor Model); contradicted by `:489` (Appendix C.10, PoseBusters Validity Checks); the same denial repeated at `body_main_short.tex:97` (Methods, Experimental Dataset).

**Verbatim evidence.** "Polar hydrogens are then re-added by the PDBQT converter under its own geometric and atom-type rules rather than by a titration model, and no pH-based protonation was applied at any stage."

**Why it matters.** The Orai1 chapter compares four receptor frames and three tools and treats receptor conformation as an explicit controlled factor — "Four states are what turn receptor conformation into an explicit factor rather than an uncontrolled nuisance variable" (`:611`). On Fr0 the AutoDock arm searched coordinates that DiffDock and EquiBind never saw, so every Fr0 between-tool statement (validity, transmembrane placement, usable yield, cross-tool clustering) carries an undisclosed preparation difference on top of the conformational one. Appendix D also builds an argument at `:609` on Fr0 being "a controlled contrast between a minimised and a thermally sampled receptor geometry", attributing the minimised character to the MD starting structure when a second, preparation-time minimisation is present in the docked file. `body_main_short.tex:156` further asserts that "All three tools dock into the same receptor file for a given frame", which is untrue on Fr0.

**Evidence missing or problematic.** `Scripts/Utilities/prep_docking.py:993-1006` runs a fallback whenever `prepare_receptor`'s hydrogen addition crashes; the embedded script at `:239-286` calls PDBFixer `findMissingResidues`, `findMissingAtoms`, `addMissingAtoms`, `addMissingHydrogens(7.0)` and then an amber14-all + implicit/gbn2 OpenMM minimisation of up to 500 steps. It fired on exactly one receptor. The staged Fr0 PDBQT in both reported panels carries 1,338 residues, 10,368 heavy atoms and 6 OXT, and is 0.819 Å heavy-atom, 0.636 Å Cα and 3.63 Å maximum from `Data/Receptors/Orai1WT-START-Fr0.pdb`, whereas staged Fr300/Fr400/Fr499 are 0.000 Å from their delivered files and DiffDock's prepared Fr0 is 0.000 Å. Verification confirmed by SHA-256 provenance: the run's own JSON for a Fr0 pose records the digest of the staged file. What is missing is why the minimisation happened, how far it moved the receptor, that six OXT atoms were added on that frame alone, and that DiffDock and EquiBind read the unminimised coordinates. The thesis is not silent — `:489` states the search read "the energy-minimised one the search actually read" — but that disclosure is in the validity-screening appendix and is contradicted by Appendix D.

**Proportionate correction.** Add a paragraph to Appendix D and a clause to the preparation summary at `:150` stating that the Fr0 AutoDock search receptor alone was passed through a PDBFixer repair and an OpenMM amber14-all/GBn2 minimisation because the converter's hydrogen addition failed on that frame; report the resulting displacement (0.82 Å heavy-atom, 0.64 Å Cα, 3.63 Å maximum) and the six added OXT atoms; state that DiffDock and EquiBind read the unminimised Fr0 so the Fr0 cross-tool comparison is not receptor-matched; add a one-clause Fr0 exception to `body_main_short.tex:156`; and reword `:558` and `body_main_short.tex:97` so they are true of the Fr0 path.

**Work required: rewriting.** No re-docking; all measurements exist on disk. Magnitude bounds the severity below critical: 0.636 Å Cα against the 5.88–6.02 Å separating Fr0 from the later frames (`:583`), on one frame of four, with the largest downstream consequence already detected and corrected (`:489`, 733 spurious minimum-distance failures rescreened away).

---

#### M2 — The blind whole-protein box places additional genuine copies of the correct binding site inside the search volume for 143 of 308 benchmark entries, and the confound is stated once but never quantified
*Finding E3-2 · omitted_from_write_up · CONFIRMED by both verifiers*

**Where.** `body_appendix_short.tex:664` (Appendix F, Composition and Provenance); protocol at `body_appendix_short.tex:167` (Appendix C, AutoDock Vina).

**Verbatim evidence.** "Where a structure contains multiple crystallographic copies of the same ligand, the comparison is made against the single deposited reference instance rather than against the closest copy."

**Why it matters.** The source benchmark gave Vina a 25 Å cube centred on the crystal ligand, so alternative copies of the site lay outside the searched space. This study's whole-receptor box puts them inside it, converting an inert convention into a live confound. A pose that correctly identifies a bona-fide copy of the binding site is scored as a total miss on both the RMSD ≤ 2 Å gate and the prespecified 4 Å crystal-site reach. Of the 308 ids in `posebusters_pdb_ccd_ids.txt`, 143 have more than one record in their `*_ligands.sdf` — nearly half the cohort, and the affected quantity is the denominator of the primary accuracy endpoint.

**Evidence missing or problematic.** Independent verification went further than the finding claimed. For all 143 entries at least one alternate copy lies more than 2 Å from the reference centroid (median nearest-other-copy separation 36.0–40.4 Å), inside the axis-aligned protein bounding box + 1 Å that `:167` defines, and within 5 Å of a protein atom; for 139 of 143 the alternate copy's 4.5 Å residue environment is ≥ 0.6 Jaccard-identical to the reference pocket on a different chain. Using the reported arm's own poses, 60 of 143 rank-1 poses land within 4 Å of an alternate copy and nearer it than the reference; 58% of the rank-1 poses in this stratum that miss the reference under the thesis's own 4 Å rule are sitting on a bona-fide copy of the correct pocket; and on a symmetry-corrected RMSD, 35–38 of 303 complexes would flip from miss to near-native at rank-1 under a nearest-copy convention (≈ +11.6 to +12.5 percentage points). DiffDock places nearly twice AutoDock's share of poses on alternate copies (20.1% vs 11.5%), so the confound is not certainly non-differential. Nothing in the main body, appendix or front matter reports the count, the change of regime, or a sensitivity arm; `grep` on "copies", "copy", "asymmetric unit", "alternate", "duplicat" locates the convention only at `:664`.

**Proportionate correction.** Report the affected count (143 of 308) in Appendix F beside the existing convention; add one sentence noting that the whole-protein box places those alternate sites inside the searched volume whereas the source benchmark's crystal-centred cube did not; and recompute near-native recovery and 4 Å crystal-site reach in a sensitivity arm scoring each pose against the nearest crystallographic copy.

**Work required: additional computation** (re-scoring of existing poses; all inputs already on disk, and the rank-1 half was reproduced during verification in minutes). Severity is major rather than critical because all three headline arms search the whole receptor and are diluted alike, so no reported conclusion is shown to invert, and the direction deflates the study's own absolutes.

---

#### M3 — Appendix D's Fr0 receptor description does not match the receptor the reported AutoDock arms docked, and the caveat built on it is stale
*Finding I-1 · scientific_flaw · CONFIRMED by both verifiers*

**Where.** `body_appendix_short.tex:558` and `:560` (Appendix D, Orai1 Receptor Model).

**Verbatim evidence.** "That search ran against 1,333 residues and 10,321 heavy atoms rather than the full assembly."

**Why it matters.** This paragraph is the thesis's only stated qualification on one of four Orai1 receptor states, and it concludes that "The reported Fr0 AutoDock validity is therefore marginally conservative for those poses". If no residues were dropped from the receptor actually docked, the caveat attaches to nothing: the claim that 42 of 240 Fr0 AutoDock poses lie within 8 Å of a dropped residue, and that two Synta-66 poses fail only against an atom of a deleted Phe253, are unsupported for the reported arm. A reader auditing Fr0 validity is sent to the wrong structure, and the qualification the arm does carry is missing from the section that exists to supply it. The error also runs opposite to its stated direction — Fr0 gained six repaired atoms rather than losing 41 — so no conservative correction is owed.

**Evidence missing or problematic.** Both reported Orai panels dock the same staged file (`Dockings/Orai_Benchmark_MGLTools_exh128` and `Dockings/Orai_JKU_MGLTools_exh128`, `_staging/receptors/pdbqt/Orai1WT-START-Fr0.pdbqt`, md5 `9fd4928d`), which carries 1,338 residues and 10,368 heavy atoms with no residues dropped. The 1,333/10,321 figures and the five named residues (Ser93 B/D, Leu248 C, Phe253 A, His256 C) reproduce `Data/Receptors/pdbqt/Orai1WT-START-Fr0_meeko_retry1.pdbqt` exactly — a per-chain alignment recovers the dropped set name-for-name and atom-for-atom — and `:152` states of the Meeko trees "neither of which is carried into the reported calibration or Orai1 results". `:363` adds "The receptor files no longer differ by converter." The paragraph is stale text from the retired Meeko arm. A further internal contradiction: `:560` asserts "PoseBusters screened all frames including Fr0 against the intact heavy-atom structure" while `:489` says the Fr0 screening file "was the deposited structure rather than the energy-minimised one the search actually read".

**Proportionate correction.** Delete the paragraph at `:560` in its entirety, including the 42-of-240 pose count and the "marginally conservative" claim; verify — as this review did, cleanly — that no Results or Discussion sentence cites the retracted caveat; and replace it with the correct disclosure under M1, naming the force field, the pH-7 hydrogen placement, the iteration cap, the resulting heavy-atom count and the displacement from the supplied frame. Reconcile `:560` with `:489` on which file was screened. One sentence at `:558` should be retained rather than deleted: "purely heavy-atom structure of 10,362 atoms" is correct for the cleaned frame.

**Work required: rewriting.** Major rather than critical because no reported figure changes and the caveat errs against the thesis's own interest.

---

#### M4 — A Results summary sentence credits DiffDock refined-score ordering with about three points of rank-1 recovery that the thesis's own two measurements report as a loss
*Finding K-1 · scientific_flaw · CONFIRMED by both verifiers*

**Where.** `body_main_short.tex:406` (Results, Calibration Dataset, Near-Nativeness and Form, closing summary); contradicted by the footnote at `body_main_short.tex:197` and by `body_appendix_short.tex:962` (Appendix G, Ranking-Recovery Statistics, table note).

**Verbatim evidence.** "Refined-score ordering can be considered an unclaimed remedy for DiffDock, worth about three points of rank-1 recovery on the reported coordinates."

**Why it matters.** This sentence closes the paragraph delivering the sampling-versus-scoring verdict for the whole calibration chapter, and it recommends an operational intervention to the reader. The thesis measured that intervention twice and both times it lost. `body_appendix_short.tex:962` states: "Re-ordering was measured separately on the same 303 complexes and does not help, since ranking the unrefined poses by their refined affinity rather than by confidence nominates a near-native rank-1 pose for 32.3% of complexes with gnina and 32.7% with smina, against the 33.7% of the confidence order." The sign is wrong, not merely the emphasis. A reader who trusts `:406` would draw the opposite operational conclusion from the one the data support. This is the cleanest instance of red flag 2 in the document.

**Evidence missing or problematic.** No number in the thesis supports "about three points" for any DiffDock ranking operation — "three points", "unclaimed" and "refined-score" each occur exactly once in the whole thesis, only in this sentence, with no cross-reference to the appendix that measured the quantity. The nearest measured quantity is the geometry-refinement gain at fixed rank, 33.7 → 35.3 and 33.7 → 35.6, summarised in the same appendix as "DiffDock starts at 33.7% and gains two". The two operations are conflated. A secondary discrepancy should be resolved while fixing this: `:197` describes "Re-ranking the refined poses by the smina score" while `:962` describes "ranking the unrefined poses by their refined affinity", reporting identical numbers; only one can be what was computed.

**Proportionate correction.** Replace the sentence with the measured result, for example: "Re-ordering the refined DiffDock poses by their refined score does not improve rank-1 recovery, at 32.7% with smina and 32.3% with gnina against 33.7% under the confidence order (Appendix G); the two-point gain at rank-1 comes from geometry refinement at fixed rank."

**Work required: rewriting.** Major rather than critical because the error does not propagate: the Discussion (`:790`), the Conclusions (`:843`) and the variant-selection paragraph (`:197`) all state the correct, opposite direction, and no headline number depends on the claim.

### 9.3 Minor issues

#### 9.3.1 Findings filed as major by a primary examiner and downgraded at validation

Each item below was verified as factually correct in its core but reduced to minor on grounds of scope, disclosure elsewhere in the thesis, or proportionality to an MSc standard. The **Verdict** column reports the two independent verifiers' own verdicts, which in a few cases sit below the panel's summary label. The **Work** column reports the remedy actually required after validation: where both verifiers reduced a filed *additional computation* to a text repair — E1-3, E1-4, E2-1, E2-2, F1-3, G-2 and I-2 — the reduced remedy is shown, and the original filing is recoverable from the panel record.

| ID | Issue | Location | Verbatim anchor | Correction | Work | Verdict |
|---|---|---|---|---|---|---|
| B-1 | The novelty clause is refuted by the thesis's own appendix for both works it cites | `body_main_short.tex:8` vs `body_appendix_short.tex:274`, `body_main_short.tex:819` | ":8 the absence of post-docking optimisation together with explicit physical-validity screening \cite{ref012}, \cite{ref035}" vs App. ":274 mirrors the post-prediction energy minimisation reported alongside the PoseBusters benchmark \cite{ref035}" | Add the missing qualifier (not applied *uniformly to every compared arm*) or relocate the two citations | rewriting | CONFIRMED both |
| B-2 | The usable-yield endpoint excludes by construction one of two named modulator regions; excluded-pose split recorded but not reported | `body_appendix_short.tex:323` vs `:581`, `:607`, `body_main_short.tex:95` | ":323 It therefore excludes a ligand resting against the lipid-facing outer surface of the hexamer just as it excludes one in the conduction pathway. This is deliberate" | One reconciling sentence in App. D; optionally a lumen-restricted sensitivity yield | reanalysis | PLAUSIBLE (1 dissent) |
| B-3 | The only biological inference is not discounted against the study's own decoy control | `body_main_short.tex:695` vs `:685`, `:701` | ":685 Both panels concentrate contacts on Asp110, Gln108, Tyr115, His113, Pro201 and Leu202. No residue differs after Holm correction of the fifteen within-tool tests." | One clause at `:695` citing `:685`, phrased as low power rather than demonstrated equivalence | rewriting | PLAUSIBLE (1 dissent) |
| C-1 | No molecular recognition, binding thermodynamics or scoring-function/search taxonomy anywhere in the thesis | `body_main_short.tex:23-41` | `grep -ci enthalp` returns 0 in all three files; "van der Waals" returns 0 in the main body | 400–600 words in Molecular Docking Process using `ref015`, `ref016`, already cited | rewriting | PLAUSIBLE (1 dissent) |
| C-2 | The Literature Review chapter cites no empirical evaluation study, including `ref035` | `body_main_short.tex:20-73` | citation set: ref009 ref010 ref014 ref015 ref016 ref022 ref023 ref024 ref027 ref055–ref059 ref067 ref068 ref073 ref075 ref076 | Move or duplicate 3–4 sentences from Motivation and §7.4 into the chapter; add a cross-reference to Table 1 | rewriting | PLAUSIBLE (1 dissent) |
| D-1 | The exclusion slab scores as unusable one of the two regions the four-state design exists to sample | `body_appendix_short.tex:323`, `:581`; `body_main_short.tex:149` | ":323 It therefore excludes a ligand resting against the lipid-facing outer surface of the hexamer just as it excludes one in the conduction pathway. This is deliberate" | One reconciling clause; optional radially bounded sensitivity | additional computation | PLAUSIBLE (1 dissent) |
| D-2 | Conclusions assert the DiffDock cap is "close to neutral", which the appendix says cannot be measured | `body_main_short.tex:836` vs `body_appendix_short.tex:365` | ":836 whereas the DiffDock cap is close to neutral on this endpoint. The five-point deficit is therefore a property of the secondary ranker rather than of the search." | Quote the measured depth sensitivity, or replace with the appendix's actual position | rewriting | PLAUSIBLE (1 dissent) |
| D-3 | The declared Fr0-versus-thermal "controlled contrast" is never reported; no outcome is resolved by frame | `body_appendix_short.tex:601`, `:611`; `body_main_short.tex:95` | ":601 Fr0 against the three later frames is accordingly a controlled contrast between a minimised and a thermally sampled receptor geometry." | Reword to "differs from", or add a per-frame control-panel validity table | additional computation | PLAUSIBLE ×2 |
| D-4 | Both enumerations of residual panel differences omit the start-conformer non-parity their own table records | `body_appendix_short.tex:361`, `:363` vs `:389` | ":361 They are therefore matched on search settings rather than on every input, and one difference remains, namely the DiffDock sampler seed." | Name the start-conformer axis in both enumerations; restate the shared item as ligand PDBQT conversion | rewriting | CONFIRMED both |
| E1-2 | The Fr0 five-dropped-residue caveat describes a Meeko file no reported arm used | `body_appendix_short.tex:560` vs `:363`, `:152`, `:489` | ":560 Its AutoDock search receptor alone required the permissive residue filter of the PDBQT converter, which dropped five residues carrying 41 heavy atoms" | Delete the paragraph and the two derived numbers (see M3) | rewriting | CONFIRMED both |
| E1-3 | The converter types 6–11 Asn/Gln/Arg nitrogens as acceptors in the thermal frames and none in Fr0, including Gln108 | `body_appendix_short.tex:558`, `:163`; `body_main_short.tex:685` | ":558 Polar hydrogens are then re-added by the PDBQT converter under its own geometric and atom-type rules rather than by a titration model" | One reporting sentence in App. D naming the per-frame counts, Gln108 and Arg77 NH1 | rewriting | PLAUSIBLE (1 dissent) |
| E1-4 | The simulation's own per-residue histidine states were discarded for a blanket cationic model, flipping 27 of 60 residues incl. His113 | `body_appendix_short.tex:558`, `:163`, `:646`; absent from Limitations | ":558 The preparation step removes every hydrogen and renames HSD, HSE and HSP to plain HIS" | Give the 33/23/4 split, state that 27 residues change, name His113; one Limitations sentence | rewriting | CONFIRMED both |
| E2-1 | The boron-patch reassurance is false for the reported Orai1 AutoDock arm | `body_appendix_short.tex:156` | ":156 Every docked 2-APB file in the study records the atom types A, C, HD, NA and OA and not one records B." | State that the reported arm docked 2-APB under the patched type (X-Score radius 1.87 Å, non-hydrophobic, neither donor nor acceptor) | rewriting | CONFIRMED both |
| E2-2 | The AutoDock arm searches an N-protonated ligand and analyses a neutral one; three statements say otherwise | `body_appendix_short.tex:150`; `body_main_short.tex:97`; `body_appendix_short.tex:807` | ":150 It adds 353 such atoms across 161 of the 308 benchmark ligands, 352 of them on nitrogen, and the docked poses carry them." | Say the hydrogens protonate aromatic *acceptor* nitrogens, are present in the searched PDBQT and absent from the template rebuild; fix the `:807` caption | rewriting (+ optional stratification) | PLAUSIBLE ×2 |
| E3-3 | The Orai1 slab is radially unbounded and no radially bounded sensitivity is reported behind the yield ordering | `body_appendix_short.tex:323`; `:581`; `body_main_short.tex:95` | ":323 it is bounded only along the pore axis and unbounded radially." | One reconciling sentence; optionally a positive extracellular-face sensitivity | reanalysis | CONFIRMED both |
| E4-1 | The DiffDock calibration run is a single unseeded draw; no stochastic arm was replicated at any seed | `body_main_short.tex:864`; `body_appendix_short.tex:301`, `:394` | ":864 The 303-complex calibration run is one unseeded draw made before the local seed patch. Rerunning it may change the pose sample. ... No arm was assessed across multiple seeds." | Optional two-seed spot check on ~60 complexes, or a one-sentence bound | new experimental work (optional) | PLAUSIBLE (1 dissent) |
| E4-2 | The panel-parity statement lists ligand preparation as shared while the adjacent table records a different start conformer | `body_appendix_short.tex:361`, `:363` vs `:389` | ":361 the receptor and ligand preparation, the AutoDock seed and the post-docking variant carried forward for each tool ... one difference remains" | Add the start-conformer axis; restate the shared item precisely | rewriting | PLAUSIBLE (1 dissent) |
| F1-1 | No test of whether a pose counted as near-native reproduces the native interactions, visually or quantitatively | `body_appendix_short.tex:543`; `body_main_short.tex:491-493`; `Thesis_short.lof` | ":543 The profiled pool is the retained valid pose set rather than its near-native subset, being 4,253 poses over 303 complexes with a median deviation of 7.57\angstrom{} from the crystal ligand" | Report contact P/R/F1 on the ≤ 2 Å subset; add one pose-versus-crystal overlay figure | additional computation | CONFIRMED both |
| F1-2 | No positive control establishing what the battery does to a deposited crystal pose under this configuration | `body_main_short.tex:140`; `body_appendix_short.tex:485-489`, `:1188` | ":140 Every pose was therefore evaluated with the PoseBusters ``dock'' battery and counted as PB-valid only if it passed all twenty-two applicable checks" | Screen each `*_ligand.sdf` through the same twenty-two checks against the same prepared receptor; add a ceiling row with a Wilson interval | additional computation | CONFIRMED both |
| F1-3 | The Results report the interaction F1 ordering without the two appendix corrections that reframe it | `body_main_short.tex:493` vs `body_appendix_short.tex:543`, `:549` | ":543 Standardising rank-1 contact recovery to a common placement distribution moves the F1 triple ... from 0.565, 0.529 and 0.471 to 0.532, 0.527 and 0.513." | Two sentences at `:493` giving the pose basis and standardised triple; add metal coordination to the `:870` profiler-defect list | rewriting | CONFIRMED both |
| F2-1 | The primary endpoint carries Holm-corrected p-values although Appendix G declares that endpoint excluded from the family | `body_appendix_short.tex:869` vs `body_main_short.tex:382`, `:834`, `Thesis_short.tex:63` | ":869 The validity-aware endpoint is not part of the family ... Its gap is reported as an effect size ... rather than tested." | One clarifying sentence in Appendix G naming the family behind the Results contrasts, reconciled with H.3:1116 | rewriting | PLAUSIBLE (1 dissent) |
| G-1 | Two reported Orai1 contrasts are structurally impossible for the three modulators under the executed atom-typing rule | `body_main_short.tex:681`; `body_appendix_short.tex:470-471`, `:545` | ":681 and all ionic, salt-bridge and attractive-charge contacts. These charged classes are absent from both tools' modulator poses." | State that no modulator carries an atom the executed criterion classifies as charged; discount those cells as the halogen cells are discounted; correct `:545` | rewriting | CONFIRMED both |
| G-2 | Interaction-class definitions do not match the executed profiler; the stated family size of thirteen is contradicted by the code | `body_appendix_short.tex:545`, `:547` | ":545 The aromatic family, namely pi-stacking, cation-pi, pi-cation, carbon-pi, donor-pi, amide-pi and alkyl-pi, uses 5.5\angstrom{}." | Correct the ligand aromaticity description (a carbon-neighbour count, not a ring test), the pi-cation/cation-pi nesting, and the repulsion threshold (effectively 5.5 Å) | rewriting | PLAUSIBLE ×2 (both weakened the family-size claim) |
| G-3 | No figure shows a docked pose with a labelled residue or a drawn contact | `body_appendix_short.tex:411`; `Thesis_short.lof` (39 figures) | ":411 The render carries no key, so none of those classes is identified by the figure itself" | One labelled render of the retained top-ranked Synta-66 pose; a colour key for Figure 22 | rewriting (one new render) | CONFIRMED (1 weakened) |
| G-5 | The Orai prefix-matching sensitivity is quantified for usable yield alone; clustering, agreement and interaction sit on the same capped set unquantified | `body_appendix_short.tex:363`, `:367`; `body_main_short.tex:640-665`, `:671-701` | ":367 raises its usable yield from 54.2 to 65.8\% pooled and its per-unit median from 60.0 to 75.0\%, while leaving the control figures untouched" | One sentence stating the bound is endpoint-limited; delete or soften "and is quantified below" | rewriting | PLAUSIBLE ×2 |
| H-1 | Conclusions assert a comparative DiffDock-cap fact with no supporting measurement | `body_main_short.tex:836` vs `body_appendix_short.tex:371` | ":836 whereas the DiffDock cap is close to neutral on this endpoint" | Quote the depth sensitivity, or state that no comparable counterfactual is available so the corrected AutoDock figure bounds rather than removes the reversal | rewriting / additional computation | PLAUSIBLE (1 dissent) |
| H-2 | The one biology-facing sentence claims recovery of a filter residue the Results say was not recovered | `body_main_short.tex:810` vs `:694`, `:685` | ":694 Each tool produces one modulator pose contacting Glu106, so neither consistently recovers the filter residue." | Name the loop residues actually recovered, drop or qualify "filter", add a half-sentence on the control panel | rewriting | PLAUSIBLE (1 dissent) |
| I-2 | The Fr0 minimisation leaves AutoDock searching a different receptor, and the frame-geometry justification does not hold for the arm it justifies | `body_appendix_short.tex:588`, `:609` | ":609 differing by only 0.3 to 0.8\angstrom{} at the gate and the filter. The compaction is confined to the loop rather than distributed over the fold." | Qualify `:609` for the supplied coordinates; cross-reference `:489`; the recomputed AutoDock Fr0 offsets are ≈ 1.1 and 1.6 Å at gate and filter | rewriting | PLAUSIBLE ×2 |
| J-1 | Table 1 is not self-contained; its Near-Nativeness column inverts the headline ordering with the denominator unstated | `body_main_short.tex:260`, `:289-295`, `:302` | ":293 AutoDock (raw) & ... 8,943 & 99.6\% & 202 & 307 & 3.4\%" against ":295 DiffDock (raw) & ... 2,198 & 24.4\% & 152 & 1,683 & 18.7\%" | Populate the empty footnote row: define "Comp.", state the pose-share basis and its divergence from complex counts, explain "Setting" | rewriting | CONFIRMED both (denominator claim weakened) |
| J-2 | "Benchmark" names more than one set; the reconciling sentence is appendix-only | `body_main_short.tex:606` and ×10 more, vs `:88`, §4.1; bridge at `body_appendix_short.tex:402` | ":606 EquiBind changes most, from a benchmark median of 60\% and mean of 52.5\% to an experimental median of 100\% and mean of 89.2\%." | One declarative sentence after `:541`; disambiguate the sign conventions at `:625` and `:671` | rewriting | CONFIRMED both |
| K-3 | The Synta-66 concordance is asserted from residues the inactive control contacts indistinguishably, and agreement is then given decisional weight | `body_main_short.tex:695`; `:685`; `:702` | ":695 Its defining extracellular loop1 and loop3 residues (His113, Tyr115, Pro201 and Leu202) lie above the filter and are characteristic of both tools' modulator poses." | One or two sentences at `:695`; the `:702` rewording is optional since "prioritise" already reads as a soft filter | rewriting | CONFIRMED both |

**Withdrawn at validation.** Finding G-4 ("The Results report the DiffDock interaction lead without the appendix's own finding that roughly four fifths of it is site-reaching") was reduced to nil. Verification established that the four-fifths statement decomposes the *rank-1 F1* separation, not the top-five Jaccard lead, and that the appendix's 10 Å restriction — the direct control for site-reaching — *strengthens* the DiffDock lead from p = 0.014 to p = 5.7e-5, holding "at every threshold from 8 to 20\angstrom{}". Adopting the proposed correction would have introduced an error.

#### 9.3.2 Remaining minor findings

The following were filed and validated as minor. "confirmed" means the verification pass reproduced the finding as stated; "weakened" means it survived in narrowed form.

**Research question and significance (B).**
- B-4 (confirmed) — Two endpoint families in the Results are covered by no research question, and the cross-tool clustering conclusion never reaches the Conclusions. `body_main_short.tex:701`: "Second, cross-tool consensus provides no reliable confidence signal." *Rewriting.*
- B-5 (weakened) — The Motivation criticises unmatched search conditions while this study equalises information, not effort or input-conformer pool (`:8` vs `body_appendix_short.tex:257`). *Rewriting.*
- B-6 (weakened) — `body_main_short.tex:539` announces a transfer test that `:593` says two of three endpoint families cannot deliver. *Rewriting.*
- B-7 (weakened) — The opening framing presents docking scores as binding-affinity estimates (`:4`, `:32`). *Rewriting.*
- B-8 (weakened) — The homology model carrying the reference-free arm has no external structural validation (`body_appendix_short.tex:572`). *Additional computation (one superposition).*

**Literature and theoretical foundation (C).**
- C-3 (confirmed) — No literature supports docking a polytopic membrane protein without a bilayer; `grep -ciE "membrane|lipid|bilayer" Literatur.bib` → 0; no Orai1 prior work in the review chapter. `body_main_short.tex:20-73 vs body_appendix_short.tex:323, :369, :581`. *Rewriting.*
- C-4 (confirmed) — The rigid-receptor approximation is never previewed in the review chapter, and `ref036` was deleted as uncited (`Thesis_short.tex:126`). *Rewriting.*
- C-5 (confirmed) — `body_main_short.tex:88` states the strong no-leakage claim without qualification; training corpora documented for gnina only. *Rewriting.*
- C-6 (confirmed) — DiffDock-L is what was run, but the review and Appendix B describe the 2023 release. `body_main_short.tex:47-50; body_appendix_short.tex:113-128, :106-107`. *Rewriting.*
- C-7 (confirmed) — `:37` asserts AutoDock Vina needs no post-docking optimisation, 34 lines before describing the pass its own headline arm runs. *Rewriting.*
- C-8 (confirmed) — Five AI docking methods named with a classification claim and no primary citation (`body_appendix_short.tex:82`). *Rewriting — the cheapest fix in the report.*
- C-9 (confirmed) — Vinardo's published advantage is enumerated incompletely, dropping the scoring-power result (`:77`). *Rewriting.*

**Methodological design (D).**
- D-6 (confirmed) — The custom thirty-conformer EquiBind harness is appendix-only while the body calls EquiBind single-forward-pass (`body_main_short.tex:54`, `:414` vs `body_appendix_short.tex:258`). *Rewriting.*
- D-7 (confirmed) — Training corpora and cutoffs of the two learned checkpoints are never stated, although `body_appendix_short.tex:143` does this properly for the rescorer. *Rewriting.*
- D-8 (confirmed) — No crystallographic quality, completeness or construct metadata for the 308 receptors, which are called "high-quality". `body_main_short.tex:88; body_appendix_short.tex:657, :771`. *Additional computation.*
- D-9 (confirmed) — The winner's-curse estimator prices the arm field but not the selection-depth degree of freedom. `body_main_short.tex:866, :195; body_appendix_short.tex:246`. *Reanalysis.*
- D-10 (confirmed) — The Orai1 homology model is never compared with any experimental structure. `body_appendix_short.tex:558, :583; body_main_short.tex:591`. *Additional computation.*
- D-11 (confirmed) — The Orai1 control panel is not property-matched to the three modulators (`body_appendix_short.tex:793`). *Reanalysis.*
- D-12 (confirmed) — The Methods never state that the Orai1 receptor contains no bilayer, ions or water, although the empty band is why the slab exists. `body_main_short.tex:145-150; body_appendix_short.tex:558`. *Rewriting.*

**Protein preparation (E1).**
- E1-5 (confirmed) — Missing heavy atoms, incomplete side chains and chain breaks neither repaired nor reported. `body_appendix_short.tex:163, :150`. *Rewriting.*
- E1-6 (weakened) — DiffDock receptor preparation never named on either dataset; Orai1 EquiBind path undescribed. `body_appendix_short.tex:251, :375-399`. *Rewriting.*
- E1-7 (confirmed) — Terminal treatment of the truncated Orai1 66–288 construct unreported; CHARMM caps stripped; frames not mutually consistent. `body_appendix_short.tex:558`. *Rewriting.*
- E1-8 (confirmed) — The Orai1 receptor contains no calcium and this is never stated, although a candidate site is defined by calcium binding. `body_appendix_short.tex:581, :558; body_main_short.tex:868`. *Rewriting.*
- E1-9 (confirmed) — The later frames are called an "equilibrium ensemble" although the appendix disclaims any record of trajectory length or frame times (`body_appendix_short.tex:613` vs `:566`). *Rewriting.*
- E1-10 (confirmed) — No stereochemical or model-quality assessment of the homology model. `body_appendix_short.tex:558, :566`. *Additional computation.*

**Ligand preparation (E2).**
- E2-3 (confirmed) — `body_appendix_short.tex:720` tells the reader the validity battery probes protonation, contradicted by `body_main_short.tex:783`; the charge census is presented as chemistry rather than a preparation artefact. *Rewriting.*
- E2-4 (confirmed) — No ionisation model and no sizing of the consequence; ligand preparation absent from Limitations. `body_appendix_short.tex:150; absent from body_main_short.tex:860-876`. *Rewriting.*
- E2-5 (confirmed) — The 30-versus-1 conformer non-parity never leaves Appendix C. `body_appendix_short.tex:148`. *Rewriting.*
- E2-6 (confirmed) — An anionic GSK-7975A arm was docked on all four frames and is silently absent from every result. `body_appendix_short.tex:793`. *Rewriting.*
- E2-7 (confirmed) — The Orai1 protocol table states all EquiBind conformers were MMFF-minimised, untrue for 2-APB (no MMFF94 boron parameters). `body_appendix_short.tex:389, :258`. *Rewriting.*
- E2-8 (confirmed) — Macrocyclic ligands never mentioned although eight are present with frozen rings; whole document, `grep -cin "macrocycl|ring-open"` returns 0 in `body_main_short.tex`, `body_appendix_short.tex` and `Thesis_short.tex`. Fix belongs in Appendix C (Docking Protocols). *Rewriting.*
- E2-9 (confirmed) — The dominant physiological species of 2-APB was not docked; the cheap sensitivity run is deferred (`body_appendix_short.tex:652`). *Additional computation.*
- E2-10 (confirmed) — The strongest defence of the preparation protocol sits in a dataset footnote 500 lines from the admission it justifies. `body_appendix_short.tex:662 footnote vs :150`. *Rewriting.*

**Binding-site definition (E3).**
- E3-4 (confirmed) — The "membrane band" was never determined, and the 22.4 Å slab is narrower than a bilayer while radially unbounded. `body_appendix_short.tex:323, :572`. *Rewriting.*
- E3-5 (confirmed) — The success criterion is purely exclusionary; the recorded atom-level rule is never exercised. `body_main_short.tex:147; body_appendix_short.tex:321`. *Additional computation.*
- E3-6 (weakened) — The Orai1 concordance evidence is computed only on poses that already passed the placement hypothesis it supports. `body_main_short.tex:671; body_appendix_short.tex:543`. *Additional computation.*
- E3-7 (confirmed) — The Synta-66 concordance is presented as support although the matched decoy panel contacts the same residues indistinguishably. `body_main_short.tex:695 vs :685`. *Rewriting.*
- E3-8 (confirmed) — The numeric threshold of the protein–ligand maximum-distance check is never stated. `body_appendix_short.tex:527, :311`. *Rewriting.*
- E3-9 (confirmed) — The receptor-parity argument covers only AutoDock and DiffDock; the metal stratification omits the one arm that retains metals. `body_appendix_short.tex:165; body_main_short.tex:870`. *Additional computation.*
- E3-10 (confirmed) — Sampling adequacy on the Orai1 box is unmeasured and the calibration result is explicitly withdrawn from that panel (`body_appendix_short.tex:1209`). *Additional computation.*
- E3-11 (confirmed) — No stratification by search-box volume despite a 56-fold range at fixed exhaustiveness. `body_appendix_short.tex:1199`. *Additional computation.*
- E3-12 (confirmed) — No internal site-boxed or predicted-pocket Vina arm; the blind-docking penalty is imported from a different benchmark. `body_appendix_short.tex:1207; body_main_short.tex:874`. *Additional computation.*
- E3-13 (confirmed) — Withholding the box is described as equalising site information, although the learned pipelines carry an implicit pocket prior. `body_appendix_short.tex:148; body_main_short.tex:872`. *Rewriting.*
- E3-14 (confirmed) — The 1 Å padding is never shown to enclose the reference poses. `body_appendix_short.tex:167`. *Rewriting.*

**Settings and pose selection (E4).**
- E4-3 (weakened) — `body_appendix_short.tex:175` says the AutoDock gnina pass re-orders "without driving their geometry" while it ran with `--minimize`. *Rewriting.*
- E4-4 (confirmed) — DiffDock's twenty denoising steps stated without citation, default declaration or sensitivity check. `body_appendix_short.tex:253, :298, :391`. *Rewriting.*
- E4-5 (confirmed) — "Exhaustiveness is the only variable" is contradicted by the same table's thread-allocation footnote. `body_appendix_short.tex:169 vs :227`. *Rewriting.*
- E4-6 (confirmed) — Version identifiers exact for every classical tool, absent for the two learned engines and their weights. `body_appendix_short.tex:152`. *Rewriting.*
- E4-7 (confirmed) — EquiBind's ordering rests on an unnamed gnina score field, in a thesis proving the head choice moves rank-1 by seven points. `body_appendix_short.tex:300, :398`. *Rewriting.*
- E4-8 (confirmed) — The Vina pose-redundancy threshold is never reported. `body_appendix_short.tex:163, :28`. *Rewriting.*
- E4-9 (confirmed) — The 6 kcal/mol window is unjustified and the pool falls 114 modes short of thirty per complex without explanation. `body_appendix_short.tex:163, :305`. *Rewriting.*
- E4-10 (confirmed) — The AutoDock ligand hydrogen-bond donor non-parity never reaches the body or Limitations. `body_appendix_short.tex:150; body_main_short.tex:868`. *Rewriting.*
- E4-11 (confirmed) — Depth-resolved recovery is reported for raw EquiBind although that arm has no ranking; the order used is named in one footnote. `body_appendix_short.tex:997, :911`. *Rewriting.*
- E4-12 (confirmed) — The Results attribute EquiBind's flat depth curve to the model while the pool is an author-built harness. `body_main_short.tex:388; body_appendix_short.tex:258`. *Rewriting.*
- E4-13 (confirmed) — The exhaustiveness-92 rung was never rescored, leaving the grid incomplete at its cheapest cell. `body_appendix_short.tex:169, :206-227`. *Additional computation.*
- E4-14 (confirmed) — Exhaustiveness 128 transferred unchanged to the Orai1 boxes with no convergence check. `body_appendix_short.tex:361, :388`. *Additional computation.*
- E4-15 (confirmed) — The primary DiffDock variant was chosen on a one-complex margin the thesis calls unstable, then carried into the reference-free arm. `body_main_short.tex:197; body_appendix_short.tex:361`. *Additional computation.*
- E4-16 (weakened) — The prespecified between-tool family names the raw and gnina-refined DiffDock arms, not the DiffDock + smina arm the Results headline. `body_appendix_short.tex:869 vs body_main_short.tex:382-384`. *Rewriting.*

**Pose-prediction validation (F1).**
- F1-4 (confirmed) — No predicted-pocket or ligand-centred Vina arm, so blind search and under-sampling separate only along the sampling axis. `body_appendix_short.tex:1207; body_main_short.tex:874`. *Additional computation.*
- F1-5 (confirmed) — One third of the benchmark comparison rests on a single unseeded stochastic draw. `body_main_short.tex:864`. *Additional computation.*
- F1-6 (confirmed) — `body_main_short.tex:113` reports three recovered complexes for the symmetrised-RMSD check where `body_appendix_short.tex:1385` reports two. *Rewriting.*
- F1-7 (confirmed) — The Form column applies an undisclosed exploded-pose filter contradicting the stated definition; 42 DiffDock poses removed. `body_main_short.tex:236-300, :205, :133`. *Rewriting.*
- F1-8 (confirmed) — No recovery rate is placed numerically beside a published value for the same tool on the same benchmark. `body_main_short.tex:797-802, :868`. *Rewriting.*
- F1-9 (confirmed) — The count of entries with multiple crystallographic copies is never reported (see M2). `body_appendix_short.tex:664, :662`. *Additional computation.*
- F1-10 (confirmed) — The body attributes EquiBind's deficit to the model class without recording the input non-parity. `body_main_short.tex:404, :398; body_appendix_short.tex:148, :258`. *Rewriting.*

**Statistics (F2).**
- F2-2 (confirmed) — `body_appendix_short.tex:1214` states bootstrap intervals are computed for effect sizes; no Cliff's delta or rank-biserial carries one. *Rewriting.*
- F2-3 (weakened) — The design-to-test table names a marginal Wilson interval as the paired-binary effect size. `body_appendix_short.tex:1102`. *Additional computation.*
- F2-4 (confirmed) — A Spearman pooled over 12,315 nested poses is reported with p = 7e-124, against the declared inferential unit. `body_appendix_short.tex:369`. *Reanalysis.*
- F2-5 (confirmed) — A pooled pose-level Fisher p at `body_main_short.tex:627` is reported as corrected while neighbouring pose-level contrasts are labelled descriptive. *Rewriting.*
- F2-6 (confirmed) — Run-to-run variation is disclosed but never quantified, and the statistics chapter never states its intervals are a lower bound. `body_main_short.tex:864; body_appendix_short.tex:1212-1216`. *Rewriting.*
- F2-7 (confirmed) — The pharmacology-concordance sentence is never connected to the same-receptor null eight lines earlier. `body_main_short.tex:695 vs :687`. *Rewriting.*
- F2-9 (confirmed) — Two omnibus statistics reported without df or p; a Kendall W of 0.055 called significant without comment. `body_appendix_short.tex:1392`. *Rewriting.*
- F2-10 (confirmed) — Cross-panel comparisons are not property-matched although the appendix documents the modulators as chemical outliers. `body_appendix_short.tex:793; body_main_short.tex:625, :703`. *Rewriting.*
- F2-11 (confirmed) — Ranking rules are benchmarked against an oracle ceiling and a legacy rule but not a random-selection floor. `body_main_short.tex:474`. *Additional computation.*
- F2-12 (confirmed) — The study-wide number of hypotheses is never stated, so aggregate false-positive exposure cannot be judged. `body_appendix_short.tex:1116-1120`. *Rewriting.*

**Results and data analysis (G).**
- G-6 (confirmed) — The control-panel null is never linked to the Synta-66 concordance, and `:701` converts it into a positive claim. *Rewriting.*
- G-7 (confirmed) — `body_main_short.tex:409` states AutoDock returns a valid pose for all 303 complexes, contradicting Table 4 (301). *Rewriting.*
- G-8 (confirmed) — The Form endpoint applies an undocumented exploded-pose guard. `body_main_short.tex:135, :293-300, :206-213`. *Rewriting.*
- G-9 (confirmed) — `:591` sends the reader to a figure that does not show the exclusion slab. *Rewriting.*
- G-10 (confirmed) — Two different non-reach counts four paragraphs apart, with no depth qualifier on the first. `body_main_short.tex:427 vs :469, :472`. *Rewriting.*
- G-11 (confirmed) — Interaction P/R/F1 printed to three decimals across 45 cells with no dispersion, interval or within-table test. `body_main_short.tex:498-527`. *Additional computation.*
- G-12 (confirmed) — The AutoDock top-15 panel of Figure 4 loses the centroid and drift arrow the text emphasises. `body_main_short.tex:409, Figure 4`. *Rewriting.*
- G-13 (confirmed) — Single-pose EquiBind series plotted at the same visual weight as thousands-of-pose distributions. `body_main_short.tex:671`. *Rewriting.*
- G-14 (weakened) — The slab removes the region the Methods cite as implicated; the two passages are never reconciled. `body_main_short.tex:95, :148; body_appendix_short.tex:323`. *Rewriting.*
- G-15 (confirmed) — No protonation sensitivity analysis, although both ligand and receptor protonation are declared modelling choices. `body_main_short.tex:97; body_appendix_short.tex:150, :163`. *Additional computation.*
- G-16 (confirmed) — No run-to-run replication of any stochastic arm while single-draw margins are reported and in places tested. `body_main_short.tex:864, :489`. *Additional computation.*
- G-17 (weakened) — Cost currencies carry no uncertainty band. `body_main_short.tex:756`. *Rewriting.*

**Interpretation and limitations (H).**
- H-4 (weakened) — The tool-selection advice at `:824` contradicts the RQ3 finding that removing the gnina pass inverts the cost ordering. *Rewriting.*
- H-5 (confirmed) — The between-tool ranking conclusion rests on a difference in significance; the standardisation that nearly erases the contact gap is not carried into the Discussion. `body_main_short.tex:794, :493; body_appendix_short.tex:543`. *Rewriting.*
- H-6 (confirmed) — The mechanism for the Orai1 ranker artefact (bilayer-free receptor, affinity objective rewarding burial) never reaches the Discussion. `body_main_short.tex:806; body_appendix_short.tex:369`. *Rewriting.*
- H-7 (weakened) — The Orai1 interpretation never carries the neutral-ligand, no-titration and ion-free choices into the reading of contacts with charged residues. `body_main_short.tex:806-812, :97, :39`. *Rewriting.*
- H-8 (weakened) — The Discussion omits the appendix figure showing most of the validity gap disappears once poses are conditioned on near-nativeness. `body_main_short.tex:777, :836; body_appendix_short.tex:1186`. *Rewriting.*
- H-9 (confirmed) — EquiBind's result is read as a property of one-shot architecture while the custom harness appears only in the appendix. `body_main_short.tex:801; body_appendix_short.tex:258, :148, :297`. *Rewriting.*
- H-10 (confirmed) — Limitations attribute the rigid-receptor limitation to AutoDock Vina alone although all three pipelines hold the receptor rigid. `body_main_short.tex:862, :95; body_appendix_short.tex:579`. *Rewriting.*
- H-11 (confirmed) — The Discussion generalises the validity-cost bound beyond the depths at which it was computed ("at any depth" for "at any reported depth"). `body_main_short.tex:777; body_appendix_short.tex:1191`. *Rewriting.*
- H-12 (confirmed) — The Abstract's "almost categorically" overstates a repair leaving pooled validity at 62.0%. `Thesis_short.tex:63`. *Rewriting.*
- H-13 (confirmed) — No research question or Conclusions section names the reference-free arm; its deliverable is stated only in Limitations. `body_main_short.tex:11-16, :829-859, :876`. *Rewriting.*
- H-14 (weakened) — Solvation is acknowledged only qualitatively, with no counterpart to the metal-adjacency sensitivity analysis. `body_main_short.tex:868`. *Rewriting.*

**Reproducibility (I).**
- I-3 (confirmed) — Code and data obtainable only by author request; no licence, DOI or commit reference; `body_appendix_short.tex:260` describes run records as "committed" that are not in version control. *Rewriting (plus an archival deposit).*
- I-4 (confirmed) — The DiffDock local patch is named but never specified; the two learned tools and the DiffDock checkpoint are unpinned. `body_main_short.tex:864; body_appendix_short.tex:152`. *Rewriting.*
- I-5 (confirmed) — Appendix A's summary of what could not be reconstructed omits the items that matter most, including the Fr0 receptor. `body_appendix_short.tex:13`. *Rewriting.*
- I-6 (confirmed) — The thesis never names its own verification harness or reproduction notebook. `body_appendix_short.tex:13`. *Rewriting.*
- I-7 (confirmed) — No operating system, kernel or GPU driver version is reported anywhere. `body_appendix_short.tex:330`. *Rewriting.*
- I-8 (confirmed) — No retrieval date or distribution version for the benchmark set. `body_main_short.tex:88`. *Rewriting.*

**Writing and presentation (J).**
- J-3 (confirmed) — `:591` names the wrong figure and mis-describes it. *Rewriting.*
- J-4 (confirmed) — "Fifteen interaction classes" contradicts the twenty-six-cell family four sentences later; the reconciliation is appendix-only. `body_main_short.tex:671, :683; body_appendix_short.tex:549`. *Rewriting.*
- J-5 (confirmed) — The Abstract and Conclusions headline a 75% figure that appears nowhere in the Results, which quotes 70.0% on a different inferential unit. `Thesis_short.tex:65, body_main_short.tex:844 vs :625; body_appendix_short.tex:367`. *Rewriting.*
- J-6 (confirmed) — Pro201 is named as a concentrated contact with a pointer to a figure whose fifteen-residue axis excludes it. `body_main_short.tex:685; body_appendix_short.tex:479`. *Rewriting.*
- J-7 (confirmed) — The arm titled "Experimental Dataset" is the one defined by its lack of an experimental reference, and the Discussion uses the adjective oppositely. `body_main_short.tex:91, :537 vs :768, :806`. *Rewriting.*
- J-8 (weakened) — Main-body figure captions are bare noun phrases, below the standard the thesis sets in its own appendix. `body_main_short.tex:419, :432, :603, :620 vs body_appendix_short.tex:423, :479`. *Rewriting.*
- J-9 (confirmed) — No Results figure or table is cited anywhere in the Discussion, Conclusions or Limitations. `body_main_short.tex:763-876`. *Rewriting.*
- J-10 (confirmed) — Structural and ligand nomenclature is inconsistent across chapters, including in the primary Orai1 results table. `body_main_short.tex:149, :695, :808, :565-571; body_appendix_short.tex:558`. *Rewriting.*
- J-11 (confirmed) — Copy-editing residue: four grammatical slips, mixed spelling variant, inconsistent tool-name and panel-letter case, three detached footnote marks. `body_main_short.tex:81, :384, :408, :591, :197`. *Rewriting.*
- J-12 (confirmed) — Figure 15's annotations are overprinted and its title carries raw statistical output; Figure 26 overruns the page (`Thesis_short.log:1799`). *Additional computation.*
- J-13 (confirmed) — The Orai1 receptor description sits in Results after its own yield table and largely repeats Methods §3.1.2. `body_main_short.tex:589-593 vs :93-95`. *Rewriting.*
- J-14 (confirmed) — The bibliography is 94 `@misc` entries with hand-typed note fields, so the declared IEEE style formats nothing. `Literatur.bib; Thesis_short.tex:7, :37, :149`. *Rewriting (optional for submission).*
- J-15 (confirmed) — The two molecular renders carry no key, labels or per-subunit colouring, and the Results chapter has no structural figure. `body_appendix_short.tex:411, :567`. *Additional computation.*
- J-16 (confirmed) — `:32` presents docking scoring as an assessment of binding affinity. *Rewriting.*

**Red-flag checklist (K).**
- K-4 (weakened) — No protonation, protomer or tautomer modelling at any stage, and the arms do not dock into the same receptor chemistry. `body_appendix_short.tex:150, :163, :653; body_main_short.tex:97`. *Additional computation.*
- K-5 (confirmed) — A three-complex margin is described as one ranking outperforming another, with no interval or test. `body_main_short.tex:197, :192`. *Rewriting.*
- K-6 (confirmed) — Metals, ions and cofactors are stripped for two of three arms with no stated scientific rationale. `body_appendix_short.tex:163, :533`. *Rewriting.*
- K-7 (confirmed) — Variants selected and evaluated on the same 303 complexes with no held-out split. `body_main_short.tex:866; body_appendix_short.tex:213`. *Additional computation.*
- K-8 (confirmed) — The Orai1 ligand geometries and receptor provenance cannot be reproduced from the record supplied. `body_main_short.tex:97; body_appendix_short.tex:572`. *Rewriting.*
- K-9 (confirmed) — The Form column applies an in-place RMSD guard the printed Methods definition denies. `body_main_short.tex:205, :133`. *Rewriting.*

### 9.4 Findings refuted and removed (4)

These were filed by primary examiners, survived adversarial challenge in some form, and were then **refuted at the independent validation stage**. They are recorded so the reader can see what was tested and rejected.

1. **D-5 — "The pharmacological concordance is never adjudicated against the study's own two controls."** Refuted. Two modulator-versus-control contact-frequency statistics sit inside the concordance paragraph itself (`body_main_short.tex:695`): "Leu109 is characteristic of both tools on the benchmark panel but of neither experimentally, at 9.4\% for AutoDock Vina and 5.6\% for DiffDock", and Phe199 at 19.4% versus 6.3%. The control result appears at `:685` and the thesis adjudicates explicitly at `:701`: "Neither cytosolic cluster supports the Synta-66 concordance." The slab-filtered pose basis is stated three times, and the proposed mechanism is contradicted by the retained cytosolic clusters. The finding also converted an underpowered Holm null into positive evidence of non-specificity, which the thesis itself forbids (`body_appendix_short.tex:482`).

2. **E3-1 — "Per-tool crystal-site reach is scored inside a pooled partition, costing AutoDock ten complexes, and no uncoupled measure is reported anywhere."** Refuted on its stated defect. The mechanism and arithmetic are correct (pooled p = 0.4704 reproducing the printed 0.470; uncoupled p = 0.0567), but the uncoupled measure *is* reported, at `body_main_short.tex:476`: "AutoDock alone reaches 88.1\% at 0.53A and DiffDock 84.5\% at 0.56A". The pooling is disclosed at `:154`, the table is captioned "Crystal-*Cluster* Reach", `:159` states "The cluster nearest the crystal is not automatically counted as recovered", and no conclusion inverts on either basis.

3. **H-3 — "The paragraph qualifying the physics-versus-learning framing never states that the AutoDock rank-1 lead is produced by a convolutional rescorer."** Refuted. `body_main_short.tex:384` states it: "The carried-forward AutoDock arm is itself CNN-rescored by gnina. The contrast is therefore between two configured workflows on this benchmark and not between physics-based and AI-based docking." `:770` adds the class-level disclaimer; `body_appendix_short.tex:184` isolates the ranker on fixed poses. The proposed correction would also have introduced an asymmetric raw-versus-refined comparison, and the thesis does not claim the rank-1 lead the finding attacked (`:382`: "not statistically resolved").

4. **K-2 — "The Orai1 usable-yield endpoint excludes by construction the implicated modulator site, and this is never stated."** Refuted. The exclusion is stated in the line the finding itself quotes (`body_appendix_short.tex:323`), with a rationale, and is disclaimed at `:325`, `body_main_short.tex:613`, `:586`, `:808` and `:870`. The endpoint is anchored to compound-specific pharmacology (`body_appendix_short.tex:654`, `body_main_short.tex:99`), and for the three compounds actually docked the proposed sites lie above the Glu106 plane and are retained. Related concerns survive in narrowed form as B-2, D-1 and E3-3.

---

## 10. Are the conclusions supported?

**Verdict: broadly supported but overstated in places.**

The four research-question answers reproduce from the Results, the scope of inference is restricted in advance and honoured, negative results are reported as negative, nulls are never converted into equivalence, and no docking score is presented as a measured affinity. The thesis is also unusually willing to weaken its own claims: it labels the head-to-head "a tuned pipeline against two untuned ones" (`body_main_short.tex:384`), bounds its selection optimism at 0.4–0.8 percentage points (`:866`), concedes that part of its validity advantage is constitutive (`body_appendix_short.tex:1188`), and dismantles its own most quotable Orai1 result. That is why the verdict is not "only partially supported".

Overstatement is confined to a small number of load-bearing sentences, three of which are identified above as M4, H-1/D-2 and H-2, together with the novelty clause at `:8` and the Abstract's "almost categorically".

**Most important overclaim, quoted verbatim** (`body_main_short.tex:836`, Conclusions, RQ1):

> "whereas the DiffDock cap is close to neutral on this endpoint. The five-point deficit is therefore a property of the secondary ranker rather than of the search."

This is the clause that licenses correcting AutoDock's Orai1 usable yield from 60% to 75% while leaving DiffDock's 65% untouched, in the chapter with the least hedging. The appendix says of DiffDock and EquiBind: "Neither tool admits a prefix counterfactual, because a ten-sample run draws its own poses rather than truncating a thirty-sample list, so the size of the bias cannot be measured from these runs" (`body_appendix_short.tex:365`), and adds only the weaker, directional statement that the bias "runs against the direction of the reported result". A dissenting verifier established that a depth counterfactual *is* computable for DiffDock and lands near-neutral, which vindicates the substance — but that number appears nowhere in the thesis, so as written the Conclusions assert an unmeasured comparative fact.

**Defensible replacement sentence:**

> "AutoDock is capped on its gnina re-ranking, and restricting it instead to Vina modes one to ten puts it at 75% and removes the reversal; DiffDock's cap is applied on its own confidence order, which carries no comparable re-ranking step. The five-point deficit is therefore a property of the secondary ranker rather than of the search, on twelve correlated units from three ligands and descriptively rather than as a tested result."

If the neutrality clause is to be retained rather than deleted, the depth counterfactual the dissenting verifier showed to be computable — DiffDock's usable yield over its full thirty-pose pool against its top-ten confidence cap — must first be run and its value printed, since no number in the thesis currently supports it.

---

## 11. Ten prioritised revision recommendations

Ordered by return on effort — the first six cost hours and remove the report's most serious findings; the last four cost days and strengthen rather than repair.

1. **Rewrite the Fr0 receptor paragraph in Appendix D and reconcile it with Appendix C.10.** Delete `body_appendix_short.tex:560` including the 42-of-240 and Phe253 claims; state the PDBFixer repair and OpenMM amber14-all/GBn2 minimisation, the 0.82 Å / 0.64 Å displacement, the six added OXT, and that DiffDock and EquiBind read the unminimised Fr0; add a one-clause exception at `body_main_short.tex:156`; reword the two "no pH-based protonation at any stage" sentences. *Removes M1 and M3.* **Work: rewriting.**

2. **Correct the wrong-signed refined-score sentence at `body_main_short.tex:406`** with the measured 32.7 / 32.3 versus 33.7 figures and a cross-reference to Appendix G, and resolve the `:197` versus `:962` refined/unrefined discrepancy. *Removes M4; closes red flag 2's clearest instance.* **Work: rewriting.**

3. **Fix the three false or unsupported claims about the author's own data**: the boron reassurance at `body_appendix_short.tex:156`, the "docked poses carry them" / "docked as neutral species" / "as Docked and Analysed" triple (`:150`, `body_main_short.tex:97`, `body_appendix_short.tex:807`), and the "close to neutral" clause at `body_main_short.tex:836`. **Work: rewriting.**

4. **Repair the Orai1 interaction analysis.** State that no modulator carries an atom the executed criterion classifies as charged, so the charged and repulsion cells are structurally zero; discount them as the halogen cells are already discounted; correct the class definitions at `:545` (ligand aromaticity is a carbon-neighbour count, pi-cation nests inside cation-pi, repulsion is truncated at 5.5 Å); and link `:695` to the control null at `:685`. **Work: rewriting.**

5. **Quantify the multiple-copy confound and add a nearest-copy sensitivity arm.** Report 143 of 308 in Appendix F, state that the whole-protein box places those sites inside the searched volume where the source benchmark's cube did not, and recompute rank-1 and top-15 recovery plus 4 Å reach against the nearest crystallographic copy. *Removes M2.* **Work: additional computation on existing poses.**

6. **Close the enumeration and cross-reference defects in one editing pass**: the start-conformer axis at `:361`/`:363`; the endpoint-limited prefix bound at `:363`; the `:591` figure reference; the "Comp." header, denominator and empty footnote row of Table 1; the "benchmark" naming collision after `:541`; the three-versus-two symmetrised-RMSD count at `body_main_short.tex:113`; the "all 303 complexes" at `:409`; and the exploded-pose guard on the Form column. **Work: rewriting.**

7. **Add the near-native-conditioned contact statistic and one structural figure.** Report contact precision, recall and F1 on the ≤ 2 Å subset per tool, carry the placement-standardised triple (0.532 / 0.527 / 0.513) into `:493`, and add one figure superposing a recovered pose on its crystal reference plus a labelled Orai1 pocket render. *Closes the largest gap in the validation chapter.* **Work: additional computation.**

8. **Screen the deposited crystal poses through the twenty-two applicable checks** against the same prepared receptors, and add the resulting pass rate with a Wilson interval as a ceiling row to the validity table. **Work: additional computation.**

9. **Add roughly 800–1,000 words to the Literature Review**: molecular recognition and the enthalpic/entropic components of binding, the non-covalent interaction classes with their physical basis, the scoring-function families placing Vina, Vinardo, the gnina CNNs and EquiBind's objective, an Orai1 target subsection summarising Appendices D and E, and citations for membrane-protein docking without a bilayer and for `ref035` and the four tabulated evaluations. All sources are already in the bibliography. **Work: rewriting.**

10. **Run the two cheap missing controls**: a two- or three-seed DiffDock replicate on a fixed subset of about sixty complexes to bound draw-to-draw variance, and a lumen-restricted or positive extracellular-face variant of the Orai1 usable-yield endpoint reported beside the current one. Optionally add the predicted-pocket Vina arm the thesis itself names as the missing intermediate control (`body_appendix_short.tex:1207`). **Work: new experimental work (bounded).**

---

## 12. Five technically demanding viva questions

1. Your Appendix D states the Fr0 AutoDock search ran "against 1,333 residues and 10,321 heavy atoms", while Appendix C.10 describes that same search receptor as "the energy-minimised one the search actually read". The staged file in both reported panels carries 1,338 residues and 10,368 heavy atoms and sits 0.636 Å Cα from the delivered frame, while DiffDock's prepared Fr0 sits at 0.000 Å. Which file did the reported Fr0 AutoDock arm search, how did it acquire six OXT atoms, and what does that do to the within-frame cross-tool agreement metric that `body_main_short.tex:156` justifies by asserting all three tools read the same file?

2. `body_main_short.tex:406` recommends refined-score ordering as "worth about three points of rank-1 recovery" for DiffDock. Appendix G reports that the same operation "does not help", at 32.7% and 32.3% against 33.7%. Which measurement is the sentence describing, and — since `:197` says "re-ranking the refined poses by the smina score" and `:962` says "ranking the unrefined poses by their refined affinity" while reporting identical numbers — which of those two operations was actually computed?

3. Your whole-protein box is the receptor bounding box plus 1 Å, with a median volume of 263,874 Å³ against the source benchmark's 15,625 Å³ crystal-centred cube. 143 of the 308 distributed entries carry more than one crystallographic copy of the reference ligand, and Appendix F scores against "the single deposited reference instance rather than against the closest copy". How many of your rank-1 misses are poses sitting on a genuine alternate copy of the correct pocket, and does the AutoDock-versus-DiffDock ordering survive a nearest-copy convention given that DiffDock places nearly twice AutoDock's share of poses on those copies?

4. `body_appendix_short.tex:545` describes the aromatic family as pi-stacking, cation-pi, pi-cation, carbon-pi, donor-pi, amide-pi and alkyl-pi, and states the correction family "holds thirteen distinct hypotheses rather than fifteen". In the executed PandaMap 4.1.0 path, the ligand aromaticity flag is a carbon-neighbour count that any internal chain carbon satisfies, pi-cation is a strict subset of cation-pi, and the 6.0 Å repulsion criterion is truncated at 5.5 Å by the neighbour search. Which of your reported interaction classes are measuring what their names say, and what happens to your Benjamini–Hochberg verdicts if the family is redefined on the executed criteria?

5. Appendix G states "The validity-aware endpoint is not part of the family, because the validity-aware set is a strict subset of the near-native set. Its gap is reported as an effect size ... rather than tested." The Results report that same validity-aware endpoint with Holm-corrected p-values (`:382`), and the Abstract states it as the primary result. Which family were 7.6e-4, 0.009 and 0.011 corrected within, and how would you demonstrate to an examiner that the multiplicity correction behind your headline ordering is auditable?

---

## 13. Five broader conceptual viva questions

1. Your design withholds the search box from all three pipelines so that "none is told where the ligand belongs". But DiffDock and EquiBind carry a pocket prior learned from training structures that no protocol choice removes. In what sense is site information actually equalised, and what would a genuinely information-matched comparison of a physics-based and a learned docking pipeline look like?

2. The thesis restricts every claim to "the three workflows as configured here rather than ... method classes". Given that your AutoDock arm's rank-1 performance depends on a convolutional rescorer and your DiffDock arm on a classical minimiser, what would have to change in the design for a method-class claim to become legitimate — and is that a question docking benchmarks can answer at all?

3. Your reference-free endpoint is a geometric admissibility rule that scores as unusable everything in the Arg91–Glu106 band at any radius, including the lipid-facing interface Appendix D gives as the reason for carrying four receptor states. What is the correct epistemic status of "usable yield" when the success criterion encodes the hypothesis under test, and how would you construct a placement endpoint that could in principle be wrong?

4. You dock a polytopic membrane channel with no bilayer, no permeant ion and no water, and then measure that your secondary ranker promotes membrane-buried poses across 12,315 control poses. Is that a failure of the scoring function, of the receptor preparation, or of the framing that treats a channel lumen as a docking site — and what does the answer imply for the whole practice of membrane-protein docking in vacuo?

5. PoseBusters validity discriminates the tools on checks that read coordinates a Vina-driven objective has already minimised, which you concede is "partly constitutive". If physical plausibility is partly a property of the objective that produced the pose rather than an independent verdict on it, what would an objective-independent physical-validity criterion have to be, and would it still be cheap enough to use as a screening gate?

---

## 14. Final weighted score

| Category | Score /5 | Weight | Contribution |
|---|---|---|---|
| Research question and significance | 4.00 | 10% | 8.00 |
| Literature and theoretical foundation | 3.50 | 10% | 7.00 |
| Overall methodological design | 4.00 | 10% | 8.00 |
| Docking protocol and technical execution | 3.88 | 20% | 15.50 |
| Validation and benchmarking | 4.25 | 15% | 12.75 |
| Results and data analysis | 4.00 | 10% | 8.00 |
| Interpretation, discussion and limitations | 4.00 | 10% | 8.00 |
| Reproducibility and transparency | 4.00 | 10% | 8.00 |
| Writing and presentation | 4.00 | 5% | 4.00 |
| **Total** | | **100%** | **79.25** |

Weights: 10 + 10 + 10 + 20 + 15 + 10 + 10 + 10 + 5 = 100 ✓
Contributions: 8.00 + 7.00 + 8.00 + 15.50 + 12.75 + 8.00 + 8.00 + 8.00 + 4.00 = 79.25

### **Final score: 79.3 / 100**

---

## 15. Overall classification

### **Strong.**

A well-designed, honestly reported and largely reproducible comparative docking study that answers all four of its research questions within a scope it declares in advance and honours. Four major findings, none critical; every one is correctable by rewriting or by re-analysis of data already on disk, and none invalidates a reported result.

---

## 16. Confidence in this assessment

**High.**

Every load-bearing quotation in this report was located at the file and line cited. Headline numbers were re-derived independently rather than accepted: the cost chain (151.42 CPU-core-hours / 32 threads + 10.27 GPU-hours = 15.00 charged hours, 68.5% for rescoring), the depth-recovery ladder, the pose-production table, the pooled and uncoupled crystal-site reach statistics (p = 0.4704 reproducing the printed 0.470), the 353-hydrogen / 161-ligand preparation audit, the 143-of-308 multiple-copy count, the Fr0 receptor atom and residue counts with SHA-256 provenance, and the four-frame pore-axis geometry table all reproduced. The supporting assertion harness was executed and returned 1364/1364 checks reproducing with exit 0.

Confidence is high rather than very high for three reasons. Two of the four major findings (M1, M3) rest on repository state that could in principle have drifted since the reported runs, although checksum provenance in the run JSON argues against it. Several minor findings carry one dissenting verifier and are marked PLAUSIBLE rather than CONFIRMED. And the Orai1 arm's biological interpretation cannot be independently adjudicated from within the thesis, because the target has no ligand-bound structure — which is the thesis's own stated position.

---

## 17. Examiner recommendation

This is a strong master's thesis that should be accepted subject to minor corrections. The candidate has designed and executed a genuinely controlled comparison of three docking workflows under blind whole-protein search, instrumented it with a depth-resolved recovery endpoint, a reference-free physical-validity battery and a properly specified statistical apparatus, and then applied it to a target where no ground truth exists while consistently refusing the inferences that absence forbids. The quality of self-criticism is the standout feature: the thesis measures its own winner's curse, concedes that part of its validity advantage is built into the objective, discovers that its most quotable application result is an artefact of its own secondary ranker, and reports all three. Reproducibility is verifiable rather than asserted. Against that, the Orai1 chapter and its receptor appendix contain a cluster of internal inconsistencies — the documented Fr0 receptor is not the receptor that was docked, a stale caveat describes a retired arm, a "controlled contrast" is declared but never reported, and a Conclusions-level comparative claim rests on a quantity the appendix calls unmeasurable — and a small number of sentences elsewhere state the wrong sign or the wrong species for the author's own data. The Literature Review chapter also does not discharge its title, filing the student's real synthesis in three other places. None of this requires new docking. Corrections 1–6 of Section 11 are hours of careful editing plus one re-scoring pass, and they would remove every major finding in this report.

---

## Closing question

**Does this thesis demonstrate that the student can independently plan, execute, critically evaluate, and communicate a scientifically sound protein-ligand docking project at master's level?**

### **Yes, with minor reservations.**

Planning is demonstrated by a design that withholds site information uniformly, prespecifies its inferential unit and denies the application arm the calibration arm's warrant. Execution is demonstrated by version-pinned, patch-disclosed pipelines whose numbers reproduce 1364/1364 under the candidate's own harness. Critical evaluation is the strongest evidence of all: the thesis quantifies its own selection optimism, concedes constitutive coupling in its validity metric, and dismantles its headline Orai1 reversal as a ranker artefact. Communication is competent, with a numerically faithful abstract and a clean build. The reservations are that the Orai1 receptor appendix documents a file that was not docked, that three sentences state the wrong sign or species for the author's own data, and that the literature chapter does not establish the theoretical foundation its title promises. All are correctable by rewriting; none reflects a failure of scientific judgement.
