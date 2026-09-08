# Examiner Review — MSc Thesis

**Title:** *Whole-Protein Docking with AutoDock Vina, DiffDock and EquiBind: Physical Validity and Near-Nativeness Across Ranking Depth*
**Candidate:** Dominik Mann · **Programme:** MSc Artificial Intelligence Engineering, FH Technikum Wien
**Artefact reviewed:** short build as of 2026-09-03 01:27 — `Thesis_short.tex`, `body_main_short.tex` (876 ll.), `body_appendix_short.tex` (1,411 ll.); 133 numbered TOC pages, main body 1–49, appendices A–H 64–133.
**Standard applied:** strong MSc, *not* PhD · **Date:** 2026-09-03

> **Relationship to the 2026-09-01 review.** An earlier pass against the 09-01 build is archived at `obsolete/EXAMINER_REVIEW_2026-09-01_full_msc_rubric.md` and scored 85/100. This is a fresh assessment of the current build, which has moved materially. Where a finding has been fixed, closed or refined since, that is stated explicitly. One of my earlier findings is **corrected** below (M4c).

---

## 1. Executive assessment

The thesis asks how three configured blind whole-protein docking workflows — AutoDock Vina + gnina rescoring, DiffDock-L + smina, EquiBind + gnina — differ in physical validity, near-nativeness, ranking depth, response to post-docking optimisation and computational cost (RQ1–RQ4, p.1). It answers with a 303-complex crystal-referenced calibration benchmark and a reference-free case study docking three Orai1 modulators into four supplied MD snapshots. Docking is used for exactly what it can deliver — geometric pose hypotheses plus an admissibility screen — and is never converted into an affinity, potency or mechanism claim.

The principal contribution is a decomposition rather than a tool ranking. Pairing a 22-check PoseBusters battery with symmetry-corrected RMSD across ranking depth separates four failure modes — generation, ranking, local geometry, localisation — and shows generation dominates: about six in ten rank-1 failures of the two leading tools contain no qualifying pose at any depth (p.16).

The strongest aspect remains validation and statistical discipline, and it has strengthened since 09-01. Appendix H.9 (p.128) converts the non-significant primary contrast into quantitative information — Newcombe interval, 112 informative discordant pairs, minimum detectable difference 9.8 pp, realised power 0.12, ~4,200 complexes required. The multiplicity architecture has since been repaired: a single pooled 24-test family was split into two families of sixteen after the candidate established that adding a within-tool optimisation row could flip an unrelated between-tool verdict (p.121, `:1118`). I re-derived Table 4 and the headline recovery figures from primary data; all reproduce exactly.

The most important weakness is that a superseded analysis generation still reaches print in the Orai1 chapter. The candidate's own 2026-09-03 coverage audit establishes that Table 12's AutoDock column does not reproduce — 18 of 35 cells differ, all in that column — and that the Figure 26 caption and one body sentence trace to superseded runs. These remain in the document.

Conclusions are supported. Quality is very strong and has improved. Confidence: **high**.

*(≈300 words)*

---

## 2. Summary of the research approach

Three end-to-end workflows run blind over the entire receptor; none is given the binding site.

| | Search / generation | Refinement | Ranking used |
|---|---|---|---|
| **AutoDock\*** | Vina 1.2.7-mod, exhaustiveness 128, seed 42, whole-receptor box, ≤30 modes | gnina 1.3.2 rescoring (`--minimize`) | gnina CNNaffinity (replaces Vina rank) |
| **DiffDock\*** | DiffDock-L, 20 denoising steps, ≤30 samples | smina local minimisation, 4 Å box (`--minimize`) | native confidence on refined coordinates |
| **EquiBind\*** | 30 × ETKDGv3+MMFF conformers, unguided, uncropped receptor | gnina (`--minimize`) | gnina affinity (EquiBind emits no score) |

**Calibration arm.** 308 PoseBusters complexes, 303 analysed after DiffDock produced no output for five (Appendix H.1, p.121, with an intention-to-treat sensitivity check). Primary endpoint: PoseBusters-valid **and** ≤2 Å symmetry-corrected heavy-atom RMSD, at rank-1, top-5, top-15, top-30. Secondary: Kabsch form ≤1 Å, cross-tool centroid clustering against the crystal pocket, PandaMap interaction fingerprints, cost per qualifying pose on a uniform charged basis.

**Experimental arm.** 2-APB (neutral), GSK-7975A, Synta-66 docked into four Orai1 MD snapshots supplied by a collaborating group. With no crystallographic ligand, near-nativeness is replaced by an operational endpoint — PoseBusters-valid, on the receptor, outside an R91–E106 transmembrane slab. A control panel docks the 308 benchmark ligands into the same frames as a background comparator; no enrichment claim is made.

**Design class.** Correctly self-described as **exploratory and comparative**. Variant selection and evaluation share the same 303 complexes with no held-out split; disclosed in the abstract, Methods and Chapter 7, with the winner's curse quantified.

---

## 3. Five principal strengths

**S1 — Validation is unusually strong for an MSc, and its limits are computed rather than asserted.**
The calibration arm is a 303-complex redocking study against experimental reference poses. Appendix H.9 (p.128): "*At rank-1 both pipelines recover 51 of the 303 complexes and neither recovers 140. Only the remaining 112 discordant complexes, 37.0 % of the set, resolve the contrast … minimum detectable difference at 9.8 percentage points for 80 % power … a power of 0.12 … roughly 4,200 paired complexes.*" Appendix H.11 (p.131) re-runs the primary gate under an independent RMSD implementation, moving two of 909 tool-by-complex decisions.

**S2 — The cost analysis was rebuilt on a single defensible currency, and it is now correct.**
On 09-01, §4.3 quoted a Friedman statistic from a *mixed* vector — AutoDock charged, DiffDock and EquiBind elapsed — while pointing at a wall-basis figure that marked the AutoDock–DiffDock pair "ns". Every cost number is now on one uniform basis, CPU-core-seconds / 32 + GPU-seconds, defined in a new Methods paragraph (§3.2.7, p.11) and printed in the figure's own subtitle. I recomputed the whole block from `docking_effort_gnina_v2_charged/per_complex_effort.csv`: AutoDock 137.2 s (IQR 75.6–167.4, n=199), DiffDock 49.4 (14.1–127.3, n=169), EquiBind 2.4 (1.0–7.3, n=80); Friedman χ² = 85.2, p = 3.2e-19, Kendall W = 0.80 on n = 53; ratios 0.14 and 0.02. **Every value matches the thesis exactly**, all three pairs now separate under Holm, and the fix propagated to §6.3 ("twentyfold") and Appendix Table 10 (6.02 / 9.01 / 0.34). This is a model correction: it did not reword the claim, it made the claim true.

**S3 — Protocol documentation is at publication standard, including the parts most authors hide.**
Every binary is version-pinned by invocation (Appendix C, p.71): `AutoDock Vina v1.2.7-20-g93cdc3d-mod`, gnina 1.3.2, ADFRsuite 1.0, PDBFixer 1.12.0/OpenMM 8.3.1, PoseBusters 0.6.3, RDKit 2025.09.1. Two local source patches are declared and bounded — a boron atom type shown to be latent because "*not one [2-APB file] records B*", and a `posebusters/modules/energy_ratio.py` handler fix. The search box is given as a construction rule with its realised distribution (edges 29.5–167.4 Å, median 63.5 Å; volumes to 2,479,703 Å³). Receptor preparation is reported to residue counts: histidine cationic in all 4,485 intact rings; 316 of 2,282 cysteines carrying an inherited thiol hydrogen.

**S4 — Confounds are hunted, quantified and reported against the author's own interest.**
Metal stripping: 78 of 303 complexes carry a metal within 5 Å of the crystal ligand, and the top-15 contrast is re-reported split by stratum (−6.7 pp metal-free vs −20.5 pp metal-adjacent, p.75), with the interaction declared unresolved. The PoseBusters cofactor asymmetry is audited to the pose: "*Removing all four checks moves 12 of its 80,100 poses across the validity threshold … The asymmetry therefore penalises EquiBind rather than inflating the reported ordering*" (p.93). Appendix H.4 (p.124) volunteers that part of the physics-based validity advantage is constitutive. Since 09-01 the multiplicity family was split after the candidate found that extending a pooled family from 24 to 32 tests cost an unrelated cross-tool contrast its significance — a self-inflicted finding acted on rather than buried.

**S5 — Numerical integrity is demonstrable, and the verification infrastructure is now exceptional.**
I independently recomputed Table 4 from `per_pose_metrics.csv`: all ten variant rows reproduce exactly (AutoDock raw 8,943/8,976 = 99.6 %, 303 complexes; EquiBind + gnina 5,633 = 62.0 %, 266; combined 164/160/83/138/138/2/48/55). Headline recovery reproduces: AutoDock\* 36.6 / 60.7 / 65.3 / 65.7 %, DiffDock\* 34.0 / 47.5 / 55.1 / 55.8 %, including the +13.2 pp top-5 gap. Cross-reference machinery is clean: 39 figures, 27 tables, zero unreferenced labels, zero broken `\ref`, zero duplicates, `\nocite` diffing empty against 82 cited keys both ways. Beyond the thesis, the accompanying harness now carries **1,282 assertions covering 23 of 27 tables and 37 of 39 figure assets**, with a completeness assertion that fails if a printed figure is neither paired nor declared, and — unusually — **adversarial mutation testing of the checks themselves**, on the principle that "*a check that still passes after its expected value is changed is vacuous and worse than no check*".

---

## 4. Five principal weaknesses

**W1 — A superseded analysis generation still reaches print in the Orai1 chapter.**
Established by the candidate's own 2026-09-03 coverage audit and still live in the document: Table 12's AutoDock column does not reproduce (18 of 35 cells differ, all in that column; e.g. AutoDock hydrogen bond benchmark 3.85 printed against 4.15 canonical); the Figure 26 caption cites 7,231 poses and Asp110 at 30.2 % against a canonical 7,458 and 28.8 %, and ends its axis on the wrong residue; and Appendix C.8 still carries the pre-Fr0 cluster (29.5 Å, 42.2 Å, "seventeen pairs") against the Results' 28.4 Å, 40.3 Å and thirteen.

**W2 — One panel-parity claim underpinning the Orai1 chapter's strongest result is factually wrong.**
Appendix C.8 (p.86, `:363`) still asserts that the analysed experimental AutoDock poses are the Vina-ranked prefix of a thirty-mode run, and exempts AutoDock from the selection-bias concession granted to the other two tools. Table 19 on the facing page records the cap as "gnina CNNaffinity order". Recomputing from `posebusters_filtered_results.csv`, **66 of the 120 analysed poses carry Vina ranks 11–30** (ρ between the two orderings = 0.35).

**W3 — EquiBind is evaluated through a custom harness on non-matched ligand inputs, and this never reaches the body.**
EquiBind natively emits one pose; the study builds a 30-pose pool from 30 ETKDGv3+**MMFF** conformers, where AutoDock and DiffDock share one ETKDGv3+**UFF** start conformer. Stated once in Appendix C.4 — "*that expansion is a custom pipeline rather than native EquiBind behaviour*" — and in a table row. I confirmed that `MMFF`, `ETKDG`, `custom pipeline` and `native EquiBind` return **zero hits in `body_main_short.tex`**, including Chapter 7, where the failure is attributed to the model.

**W4 — The receptor-state factor is declared load-bearing and then never analysed.**
Appendix D.1 (p.104, `:609`) still states that four states are "*what permits the frame-level tests the Results chapter reports*". No per-frame Orai1 result, figure or test exists in Chapter 4; the only frame mentions are naming and figure captions. Fr0 is measured at 5.88–6.02 Å from the other three (which sit 0.93–1.19 Å apart) and is pooled with them throughout.

**W5 — The thesis understates its own reproducibility, and omits basic availability hygiene.**
Appendix A (p.64) still describes only a regeneration guide and concedes "*the tables whose provenance it did not audit*". That sentence is now false in the candidate's favour: 23 of 27 tables are asserted. Meanwhile the repository is request-gated with no licence, no DOI and no archived snapshot; no release tag or retrieval date is given for the PoseBusters distribution; and DiffDock-L and EquiBind carry no commit hash or checkpoint identifier.

---

## 5. Detailed scoring table

| # | Category | Wt | Score /5 | Weighted | Confidence |
|---|---|---:|---:|---:|---|
| 1 | Research question and significance | 10 % | 4 | 8.0 | High |
| 2 | Literature and theoretical foundation | 10 % | 4 | 8.0 | High |
| 3 | Overall methodological design | 10 % | 4 | 8.0 | High |
| 4 | Docking protocol and technical execution | 20 % | 4.5 | 18.0 | High |
| 5 | Validation and benchmarking | 15 % | 5 | 15.0 | High |
| 6 | Results and data analysis | 10 % | 4 | 8.0 | High |
| 7 | Interpretation, discussion, limitations | 10 % | 4 | 8.0 | High |
| 8 | Reproducibility and transparency | 10 % | 4.5 | 9.0 | High |
| 9 | Writing and presentation | 5 % | 4 | 4.0 | High |
| | **Total** | **100 %** | | **86.0** | **High** |

Arithmetic: 8.0 + 8.0 + 8.0 + 18.0 + 15.0 + 8.0 + 8.0 + 9.0 + 4.0 = **86.0**. Weights sum to 100 %.
*Change from 09-01 (85.0): reproducibility 4 → 4.5, on the 1,282-assertion harness, the completeness assertion and the mutation testing of the checks. Every other category is unchanged; the cost fix repaired a defect rather than adding capability.*

### Per-category justification

**1 · Research question and significance — 4 (High).** Four explicit, numbered, answerable questions mapped one-to-one onto four Conclusions sections, scoped in the same breath to configured workflows rather than method classes (p.2). *Deficiency:* RQ4 is answered with a prescriptive staged workflow ending in an "endpoint-specific selector" defined nowhere; the Orai1 arm occupies a third of the Results without being a numbered question and carries two mutually inverted purpose statements (p.5 vs p.26).

**2 · Literature and theoretical foundation — 4 (High).** Appendix B.1.2 (pp.65–67) reproduces the Vina function completely — both Gaussians, repulsion parabola, hydrophobic and H-bond terms, rotatable-bond divisor, 8 Å truncation, all six weights — then the Vinardo re-parametrisation with equal precision. §2.1.1 correctly notes the absence of explicit electrostatic and desolvation terms. 45 of 94 bibliography entries are 2023 or later. *Deficiency:* Chapter 2 contains no Orai1, CRAC-channel or membrane-protein-docking literature at all, though Orai1 carries half the thesis; a body-only reader meets the target first in the Methods. The physical reason lipophilic ligands land in a bilayer-stripped membrane band is never stated, though the resulting loss is reported as a tool property.

**3 · Overall methodological design — 4 (High).** Appendix C.2 (p.80) quantifies the winner's curse with the standard estimator ("*near two percentage points at rank-1*"); Appendix F.10 (p.116) states the calibration-to-Orai1 transfer as "*an upper bound … a stated assumption and not a measured result*". *Deficiency:* W4. Secondarily the selection field is unequal (8 AutoDock arms vs 3 refiner arms per learned pipeline), so the differential optimism is of the same order as the +2.6 pp headline rank-1 difference, and §6.1 reports that difference without the caveat.

**4 · Docking protocol and technical execution — 4.5 (High).** Reporting standard exceeds most published work; deductions for two gaps and one non-parity. *Deficiency:* handling of missing residues, missing side-chain atoms and chain breaks is absent from the entire thesis — conspicuous because every other preparation decision is reported to the residue count. DiffDock's receptor-preparation recipe is never stated. Realised per-complex pose budgets are never reported, so the top-15/top-30 pools are not demonstrably matched across tools.

**5 · Validation and benchmarking — 5 (High).** See §7. Strengthened since 09-01: the candidate found that Table 21 had been tested against the raw Vina arm rather than the declared gnina arm "in no commit" of the script's history, established that **the thesis was right and the script was wrong**, and pinned the corrected contrast in the harness. *Deficiency:* the refinement-failure fallback is still described two incompatible ways (Appendix C.3 "*failed outputs excluded*" vs C.5 "*retained the unoptimised pose under the optimised label*") and its magnitude is never reported for the three carried-forward arms.

**6 · Results and data analysis — 4 (High).** Distributions rather than winners; negative results printed; denominators normally explicit; effect sizes beside every p-value; a declared and now-repaired multiplicity architecture. *Deficiency:* W1 — and it is worse than an inconsistency, because the candidate has established that one printed table does not reproduce from any artefact on disk and has left it printed. Plus a scatter of local defects (§9).

**7 · Interpretation, discussion and limitations — 4 (High).** The overclaim hunt again came back essentially empty; every candidate over-reading is hedged in place. §5.1 (p.42): "*Physics-based is itself an operational label here … not a molecular-mechanics or free-energy simulation.*" Chapter 7 runs to ~25 limitations, several damaging to the author's own headline. Improved since 09-01: §6.2 now discloses that the rank-1 re-ranking step "*does not clear Holm correction on its own*". *Deficiency:* Appendix H.4's finding that the validity spread "*narrows to roughly four percentage points once the pose is in the right place*" still never reaches §5.2, §6.1 or the abstract, all of which present the raw spread.

**8 · Reproducibility and transparency — 4.5 (High).** Procedural reproducibility is now excellent. The harness covers 23 of 27 tables and 37 of 39 figure assets, records the two it cannot generate (both one-off PyMOL renders whose captions already disclaim measurement), and validates its own checks by adversarial mutation. It also caught two silent-input defects — a stage passing one of nine paths and defaulting the other eight to superseded trees, and a diagnostics script defaulting to a non-canonical tree. *Deficiency:* W5. **Conclusion reproducibility remains weaker than procedural reproducibility**: no arm was replicated across seeds and the 303-complex DiffDock run is "*one unseeded draw made before the local seed patch*" (p.48), so no run-to-run component enters any interval.

**9 · Writing and presentation — 4 (High).** Prose is precise and economical; abstract and Kurzfassung agree and track the body; cross-reference machinery verified clean. *Deficiency:* §4.2.1 (p.28) still sends the reader to Figure 27 "*with the transmembrane exclusion slab*", but that figure shows no slab, no ligands and no labels — at the one point where the reference-free endpoint is defined. Main-body captions remain one-line and not self-contained. The printed Methods still carries no protocol section.

---

## 6. Docking-protocol audit

Legend: **R** reported adequately · **I** reported but inadequately justified · **A** absent.

### 6.1 Protein preparation

| Item | | Evidence / assessment |
|---|:-:|---|
| Crystallographic waters | **R** | "*Crystallographic water is absent from the distribution itself, so every arm docks a water-free receptor*" (p.71). §3.1 and Ch.7 scope removal to "the AutoDock and DiffDock arms", under-describing the appendix. |
| Conserved waters | **A** | No count of natives with a bridging water and no stratification of recovery on it, though the analogous metal stratification is done in full. Partly inherent to docking; the analysis was within reach. |
| Cofactors, ions, metals | **R** | Stripped for AutoDock/DiffDock, retained for EquiBind; asymmetry quantified twice (pp.75, 93), affected stratum counted (78/303) and the endpoint re-reported. Exemplary. |
| Hydrogen addition | **R** | `prepare_receptor -A hydrogens`, non-polar merged; EquiBind uses `reduce 4.15.250408` without amide/His flips. |
| Protonation-state assignment | **I** | Uniform, pH-free, fully reported — not justified, not sensitivity-tested. Histidine cationic in **all 4,485** intact rings; EquiBind's `reduce` leaves both imidazole nitrogens bare on 3,898 of 4,844. Three arms, three receptor protonation models, changing the donor/acceptor set Vina's H-bond term reads. |
| His/Asp/Glu/Lys/Cys/termini | **R** | Enumerated with counts, including "*Cysteine is the exception, since the preparation builds no thiol hydrogen*" and unrepaired termini. Disulfide-bridged vs free cysteines not split. |
| Missing atoms / residues / loops / chain breaks | **A** | **Nowhere in the thesis.** Clearest single omission in the audit. |
| Structural repair / minimisation | **I** | PDBFixer for stripping and modified-residue mapping only; the Orai1 Fr0 receptor was additionally energy-minimised, stated in C.8 but absent from Appendix D where the frame is characterised. |
| Alternate conformations | **R** | "*No atom record … carries an alternate-location indicator*" (p.71); occupancy handling stated. |
| Biological assembly vs ASU | **R** | Docked as supplied; the benchmark authors had already withdrawn 120 entries whose ligand contacts a symmetry mate. |

### 6.2 Ligand preparation

| Item | | Evidence / assessment |
|---|:-:|---|
| Structure verification | **R** | Heavy-atom conformer identical to source in all 308 cases; formal-charge tally given. |
| Stereochemistry | **R** | The four reference-based identity checks excluded by design but measured post hoc: applying them "*reduces DiffDock recovery … by one complex at rank-1*" (p.9 footnote). |
| Protonation / ionisation | **I** | "*no protomer or tautomer enumeration and no pKa model ran at any stage*" (p.71). Disclosed, uniform, non-differential — but Appendix F.4 then reads the source files' neutral formal charges as a chemical property of the ligand set, which it is not. |
| Tautomer generation | **A** | None run; disclosed. Acceptable at MSc level for a fixed benchmark; consequence not weighed. |
| Salt removal | **R** | Not applicable — benchmark ships single desalted instances. |
| 3D geometry / conformers | **I** | **W3.** AutoDock and DiffDock share one ETKDGv3+UFF conformer; EquiBind uses 30 × ETKDGv3+MMFF. Stated once, never discussed as a confound. |
| Partial charges | **R** | Gasteiger, written to both PDBQTs — and correctly noted never to reach the objective, since neither Vina nor Vinardo carries a charge-dependent term. |
| Energy minimisation | **R** | UFF for the shared start conformer; the EquiBind in-protein UFF stage switched off for every reported run, with the reason given. |
| Macrocycles / metals / covalent / unusual chemistry | **R** | Boron patch with demonstrated latency; UFF-unparameterisable ligands handled by the `energy_ratio` patch. |
| Multiple states considered | **A** | One state per ligand. For 2-APB the appendix concedes the docked neutral form may be the minor species. |

### 6.3 Binding-site definition

| Item | | Evidence / assessment |
|---|:-:|---|
| Site rationale | **R** | Deliberately blind; no site supplied to any arm, and the reason argued. |
| Co-crystal ligand / known residues used? | **R** | No — the absence is the design; leakage from a ligand-centred box explicitly rejected and quantified against the benchmark's 25 Å cube. |
| Grid coordinates and dimensions | **R** | Construction rule plus realised distribution (p.76). Orai1: per-frame boxes 91.8–102.0 Å per edge, 0.375 Å grid. |
| Box too small / large / biased? | **I** | Deliberately very large, cost acknowledged and cited, exhaustiveness raised in response. Not verified that every crystal ligand lies inside its own box (1 Å margin). |
| Multiple / allosteric sites | **R** | Central to the Orai1 design. **But** the radially unbounded slab discards the lipid-facing M3 interface that Appendix D.1 itself names as a candidate site. |

### 6.4 Docking settings

| Item | | Assessment |
|---|:-:|---|
| Software + exact versions | **R** | Verified by invocation; two local patches declared. **Except** DiffDock-L and EquiBind, which carry no commit hash or checkpoint id. |
| Engine and scoring function | **R** | Vina function throughout; Vinardo arm excluded with its reason. |
| Parameter choices + justification | **R** | Five-rung exhaustiveness ladder with paired McNemar tests and a dedicated appendix section. |
| Number of poses generated | **I** | Budget stated (≤30); **realised** per-complex counts never reported. |
| Random seeds / repeats | **I** | Seed 42 for Vina and the gnina pass; several other seeds described as fixed but unstated; **no arm replicated across seeds**; benchmark DiffDock run unseeded. |
| Receptor rigidity | **R** | All three rigid; stated. |
| Flexible residues / constraints | **A** | None used; correctly not claimed. |
| Symmetry / equivalent poses | **R** | Graph-isomorphism mapping, residual limitation stated and bounded (adds three complexes, removes none). |
| Final-pose selection criteria | **R** | Explicit per arm; the ranking-head choice quantified against alternatives, including the head that would have scored better (CNNscore 37.0 % vs adopted CNNaffinity 36.6 %). |
| Rescoring / consensus / clustering | **R** | All three used, each with thresholds and rationale. |
| Visual inspection | **R** (absent by design) | No pose selected by eye; both renders disclaim measurement. |

---

## 7. Validation assessment

**Verdict: unusually strong for an MSc thesis**, with one gap limiting generalisation rather than internal validity.

**Pose-prediction validation.** A 303-complex redocking study against experimental reference poses. Ligands are absent from the docked receptor coordinates; RMSD is symmetry-corrected by graph isomorphism, computed in place without superposition, with the residual limitation stated and bounded; the 2 Å threshold is conventional and cited; form is separated from placement by Kabsch superposition and the decomposition honestly labelled approximate. Recovery is resolved over ten thresholds (H.10), four depths, and split by metal adjacency, and the primary gate is re-run under an independent RMSD implementation.

Two caveats bound it. It is **self-docking** throughout — every complex redocked into its own cognate crystal receptor — and no cross-docking arm exists, so the cognate advantage is asserted rather than estimated. And **no visual confirmation of key interactions is presented**: there is not one figure showing a docked calibration pose beside its crystal reference. Interaction agreement is measured rather than shown, which is defensible, but one illustrative panel would let an examiner sanity-check the pipeline at a glance.

**Virtual-screening validation.** Correctly **not applicable**. No enrichment, ROC-AUC, BEDROC or decoy claim anywhere; the control panel is explicitly "*only a background comparator, and no enrichment claim is made*" (p.26). Their absence is **not** a deficiency and I do not penalise it. I checked for implied screening utility and found none.

**Ranking / affinity claims.** Correctly **refused**. No score is correlated with any measured affinity, and Chapter 7 says so. The Orai1 IC₅₀ data in Appendix E frames the placement rule and is never correlated with scores; its footnote states plainly that "*Proposed sites are model-derived and filter-proximal rather than experimentally resolved ligand-bound poses.*"

**Practical consequence.** The evidence supports: (i) **pose hypotheses** at the stated rates; (ii) **relative ranking of the three configured workflows** on this benchmark at the depths where intervals separate (top-5 onward), but **not** at rank-1, where the study is explicitly underpowered; (iii) **descriptive comparison** on Orai1; (iv) **exploratory, testable site hypotheses**. It does **not** support a method-class ranking, a prospective accuracy estimate on unseen targets, a virtual-screening prioritisation, or any statement about Orai1 binding modes. The thesis claims (i)–(iv) and disclaims the rest.

---

## 8. Docking-specific red-flag checklist

| # | Red flag | Verdict | Deciding evidence |
|---|---|---|---|
| 1 | Docking scores treated as experimental binding energies | **ABSENT** | "*no score used for ranking was compared against measured affinities*" (p.48). "kcal" occurs zero times in main body and abstract. |
| 2 | Small score differences presented as decisive | **ABSENT** | The 2.6 pp gap is labelled unresolved everywhere with interval, exact McNemar p and power; the refiner ordering is called "*unstable*". |
| 3 | Docking protocol not validated | **ABSENT** | 303-complex reference-based validation plus a five-rung ladder with paired tests. |
| 4 | Receptor selected without scientific rationale | **ABSENT** | Four-state design argued from receptor rigidity and two literature-implicated regions; "*No conformation-selection heuristic was applied*". |
| 5 | Protonation, tautomerism or stereochemistry ignored | **PARTIALLY PRESENT, fully disclosed** | Nothing ignored — every choice recorded — but no enumeration ran, and three arms carry three receptor protonation models without a comparability statement. |
| 6 | Waters / metals / cofactors removed without justification | **ABSENT** | Justified on arm-parity grounds, cost quantified (78/303), endpoint re-reported by stratum. |
| 7 | Binding site or grid insufficiently documented | **ABSENT** | Construction rule, realised distribution, grid spacing, and a dedicated appendix section (H.5). |
| 8 | Only best-looking poses or compounds reported | **ABSENT** | Full 30-pose pools tabulated; all eight AutoDock rungs printed; failing pocket-guided EquiBind arms reported; the better-scoring CNN head reported and not adopted. |
| 9 | Visual inspection without transparent criteria | **ABSENT** (exemplary) | Both molecular renders explicitly disclaim measurement; every placement statement defined by an equation (C.6). |
| 10 | Favourable pose treated as proof of binding | **ABSENT** | "*evidence from different tools cannot be combined to validate one binding mode*" (p.34). |
| 11 | Docking treated as proof of biological activity or mechanism | **ABSENT** | Mechanism imported from cited pharmacology to frame the placement rule, never inferred from poses. |
| 12 | Experimental validation implied but not performed | **ABSENT** | Named as required and absent: "*Validation requires site-directed mutagenesis, competition or displacement assays and ideally a ligand-bound structure*" (p.49). |
| 13 | Ligand or receptor preparation not reproducible | **ABSENT** | Command lines, flags, versions, seeds and patches given; the one irreproducible element (supplied Orai1 model provenance) stated as unavailable. |
| 14 | Same data used to tune and validate the protocol | **PRESENT — disclosed, quantified, bounded** | "*Variant selection and evaluation use the same 303 complexes, with no held-out split*" (p.48), flagged in the abstract, corrected with the standard estimator, shown robust to substituting the losing refiner. Residual: the *differential* optimism is never subtracted from the headline margin. |
| 15 | Conclusions exceed what docking alone can establish | **ABSENT** | "*The result is a workflow, not a universal winner*" (p.47). |

**Net:** 12 decisively absent, items 9 and 13 handled above MSc level, item 5 partially present but recorded rather than ignored, item 14 present but disclosed, quantified and shown not to change the ordering.

---

## 9. Critical, major and minor issues

### Critical issues

**None.** No defect I could verify invalidates the main results or conclusions. The calibration findings (RQ1, RQ2) rest on data I reproduced cell-for-cell; RQ3 now reproduces from first principles on a single currency; and the Orai1 arm is framed as descriptive and exploratory throughout, so W1 and W2 damage claims the thesis has already declined to lean on.

*(The one item labelled "critical" in the candidate's own `FINDINGS_2026-09-03.md` — Table 21 tested against the wrong AutoDock arm — was a **script** defect, not a thesis defect. The printed table was correct and the pipeline was wrong. It is fixed and pinned.)*

### Major issues

---

**M1 — A superseded analysis generation still reaches print in the Orai1 chapter.**

*Where.* Table 12 and its footnote, Appendix C.8 p.92 (`:471`); the Figure 26 caption, p.90; Appendix C.8 pp.87 (`:414`, `:418`); body `:683`; §4.1.4 p.24 vs Appendix C.10 p.99.

*Evidence.* Confirmed live in the current build, and independently established by the candidate's own `REPRODUCTION_COVERAGE_2026-09-03.md`:

| Quantity | Printed | Canonical |
|---|---|---|
| Table 12 AutoDock column | 18 of 35 cells stale (e.g. H-bond benchmark 3.85, Cliff's δ +0.42) | 4.15, +0.45 — all 14 DiffDock cells reproduce exactly |
| Figure 26 caption | 7,231 poses; Asp110 30.2 %; axis ends Asp114 | 7,458 poses; 28.8 %; axis ends Pro146 |
| Table 12 footnote | eleven of twenty-six cells; 0.142 for AutoDock | fourteen; 0.089 |
| Body `:683` | smallest ligand-unit adjusted values 0.086 and 0.106 | **0.089** and 0.106 |
| Appendix C.8 | 29.5 → 27.3 Å; 42.2 → 30.9 Å; "seventeen pairs" | 28.4 Å; 40.3 Å; thirteen |
| Appendix C.10 | baseline Jaccard p = 0.024 on "the shared cohort the Results chapter reports" | §4.1.4's cohort is 264 complexes at p = 0.014; C.10's is a 247-complex EquiBind-limited pool |

*Correction to my 2026-09-01 review.* I previously wrote that for the 0.086-vs-0.142 pair "the body is right and the appendix is stale". That is **wrong in one direction**: the candidate's audit isolates three generations, and the canonical value is **0.089**. The body's 0.086 is itself from a superseded run — a different one from the footnote's 0.142. DiffDock's 0.106 agrees across all three, which is why the discrepancy is invisible unless the AutoDock arm is isolated.

*Why it matters.* These are the descriptive statistics and multiplicity-corrected inferences behind the Orai1 clustering and interaction sections. An examiner cross-checking body against appendix finds two or three answers and no way to tell which supersedes. The Table 12 case is the sharper one: a printed table that the candidate has established cannot be produced from anything in the repository.

*Proportionate correction.* Regenerate Table 12's AutoDock column, the Figure 26 caption, the C.8 symmetry-folding decomposition and the C.10 sensitivity audits from the canonical `_matched` trees, then propagate one set of numbers and correct body `:683` to 0.089. Where a genuinely different cohort is intended (C.10), say so and give its size. The harness now pins these values, so the fix is verifiable.

*Effort.* **Reanalysis over existing outputs, plus rewriting.** No re-docking.

---

**M2 — The Orai1 panel-parity claim is factually wrong, and it protects the chapter's strongest result.**

*Where.* Appendix C.8, p.86 (`:363`), against Table 19 p.87 (`:392`) and Methods p.13 (`:199`). Result affected: §4.2.3 p.30 and §4.2.5 p.37.

*What it says.* "*For AutoDock Vina the mode count and the energy window govern only how many of the searched modes are written out, so the first ten modes of a thirty-mode run are the same ranked prefix that a ten-mode run reports … For DiffDock and EquiBind the prefix argument is weaker.*"

*Why it matters.* This sentence is the entire defence of comparing a 30-mode experimental panel against a 10-mode control panel, and it explicitly exempts AutoDock from the selection-bias concession granted to the other two. AutoDock's 28.7-pp benchmark-to-experimental usable-yield loss is "*the only result significant at both units and supported by all three ligands*" (p.37).

*Evidence.* The carried-forward variant is Vina + gnina, and Table 19 records the cap as "*gnina CNNaffinity order*". Recomputing from `orai_jku_mgltools_exh128_fr0corrected/dock/posebusters_filtered_results.csv`, the top ten by `optimized_rank` within each of the twelve frame-ligand units gives exactly 120 poses, of which **66 carry `autodock_rank` 11–30** (25 / 16 / 25 for 2-APB / Synta-66 / GSK-7975A). Spearman(`autodock_rank`, `optimized_rank`) = 0.35.

*Proportionate correction.* Rewrite `:363` to state the cap is applied in gnina CNNaffinity order over all thirty written modes; extend the direction-of-bias paragraph at `:365` to AutoDock; add a one-paragraph sensitivity recomputing experimental AutoDock usable yield restricted to Vina ranks 1–10. Note that the correction **preserves the conclusion** — the direction-of-bias argument already made for the other two tools extends to AutoDock — which is why this is major rather than critical.

*Effort.* **Rewriting plus one recomputation over data already held.**

---

**M3 — EquiBind is evaluated through a custom harness on non-matched ligand inputs, and this never reaches the body.**

*Where.* Appendix C.4, p.81 (`:258`) and the "Start conformer" row of Table 9, p.83 (`:297`).

*Evidence.* "*To create a comparable retrospective pool, this study docks thirty fixed RDKit ETKDGv3 conformers, minimised with MMFF … That expansion is a custom pipeline rather than native EquiBind behaviour.*" Table 9: AutoDock and DiffDock "RDKit ETKDGv3 + UFF"; EquiBind "30 × ETKDGv3 + MMFF". Verified 2026-09-03: `MMFF`, `ETKDG`, `custom pipeline` and `native EquiBind` return **zero hits** in `body_main_short.tex`, including Chapter 7.

*Why it matters.* Two consequences, neither carried forward. (i) The comparison is carefully matched on the receptor and search-space axes — argued at length — but **not** on the ligand axis for EquiBind, so an unknown part of its deficit is the conformer generator and force field. (ii) The reported "EquiBind" is a 30-conformer ensemble ranked by an external CNN, so 18.2 % rank-1 corresponds to no deployment anyone would run. Meanwhile the body attributes the failure to the model: "*A single-shot regression model provides few distinct alternatives*" (p.17).

*Proportionate correction.* Two sentences in §3.1 or §4.1 and one clause in Chapter 7; then report, as a one-line baseline, what EquiBind's own single forward pass on the shared start conformer achieves. The poses exist.

*Effort.* **Rewriting, plus one tally over existing output.**

---

**M4 — The receptor-state factor is declared load-bearing and then never analysed.**

*Where.* Appendix D.1, p.104 (`:609`), against all of Chapter 4.2.

*Evidence.* "*Four states are what turn receptor conformation into an explicit factor rather than an uncontrolled nuisance variable, which is what permits the frame-level tests the Results chapter reports.*" No per-frame Orai1 result, figure or test exists in Chapter 4; grepping the frame names returns only the Methods introduction, the receptor description and two figure captions. Fr0 sits 5.88–6.02 Å from the other three (which sit 0.93–1.19 Å apart) and its AutoDock search receptor additionally lost five residues (41 heavy atoms) to a permissive filter.

*Proportionate correction.* Either add a per-frame breakdown of the usable-yield table with a Friedman test across frames within tool plus a leave-Fr0-out sensitivity, or delete the claim at `:609` and state that the four states are pooled to represent receptor variability with no frame-level inference attempted.

*Effort.* **Reanalysis (cheap) or rewriting.** Either removes the contradiction.

---

**M5 — Appendix H.4's placement-versus-validity decomposition never reaches the Discussion, Conclusions or abstract.**

*Where.* Appendix H.4, p.124, against §5.2 p.42, §6.1 p.46 and the abstract. Verified: "four percentage points" occurs once in the appendix and **zero times** in the body.

*Evidence.* "*Intermolecular failure for DiffDock falls from 14.6 % over all produced poses to 5.4 % over its 1,715 near-native poses, and for EquiBind from 37.0 % to 6.2 % … The 98.9 % against 62.0 % spread of Table 4 narrows to roughly four percentage points once the pose is in the right place.*"

*Why it matters.* This is the most interesting negative result in the thesis — most of the headline validity separation is a *placement* effect, not a geometry-quality effect. The Discussion offers only the unquantified "*partly constitutive*"; the Conclusions print the raw 98.9 / 84.5 / 62.0; the abstract prints 99.6 / 24.4 / 2.9. A reader of the abstract alone would substantially over-read the validity gap.

*Proportionate correction.* One sentence in §5.2 and one clause in §6.1 carrying the four-point figure and its conditioning.

*Effort.* **Rewriting only.** The analysis is done.

---

**M6 — Handling of missing residues, side chains and chain breaks is absent from the entire thesis.**

*Where.* Appendix C.1, p.74 (`:163`); nowhere else. Verified: zero hits for "missing residue", "missing side", "chain break" or "findMissing" across both body files.

*Why it matters.* Under a whole-protein blind search an unmodelled surface loop is a cavity the search can occupy, so incompleteness interacts directly with the design. The effect is likely small and non-differential — the same receptor file feeds AutoDock and DiffDock — but the reader cannot tell. It is conspicuous precisely because every other preparation decision is reported to the residue count.

*Proportionate correction.* One paragraph in Appendix C.1 stating the PDBFixer configuration actually used (`findMissingResidues` / `findMissingAtoms` on or off) and cohort-level counts of missing residues and incomplete side chains, over files already held.

*Effort.* **Rewriting plus a tally.**

---

**M7 — The printed Methods chapter contains no docking protocol, no software version and no seed.**

*Where.* Chapter 3, pp.5–11; all protocol material in Appendix C, pp.71–100. Verified: Chapter 3 still contains only Datasets and Evaluation Metrics.

*Why it matters.* A defect of the artefact rather than the science. An examiner reading the 49-page main body cannot state what software, version, box, exhaustiveness or seed produced any number, and Chapter 3 contains no forward pointer saying Appendix C holds it. It also drives the mass imbalance: 49 pages of body against 70 of appendix.

*Partially improved since 09-01:* §3.2.7 (p.11) gained a paragraph defining the charged cost basis, which is exactly the right pattern — applied to one setting only.

*Proportionate correction.* A one-page §3.3 "Docking protocol" carrying a compact table (tool · preparation · search settings · pose budget · seed · refiner · ranking) with one forward reference to Appendix C, and a home in the printed body for the software-version manifest.

*Effort.* **Rewriting only.** The table can be condensed from Table 9.

---

**M8 — Appendix A understates the project's own reproducibility, and omits basic availability hygiene.**

*Where.* Appendix A, p.64 (`:13`).

*Evidence.* Appendix A still describes only a regeneration guide and concedes it marks "*the tables whose provenance it did not audit*". As of 2026-09-03 the accompanying harness carries **1,282 assertions covering 23 of 27 tables and 37 of 39 figure assets**, with a completeness assertion over printed figures and adversarial mutation testing of the checks. None of this is mentioned in the thesis. Separately: the repository is request-gated and not publicly resolvable, with no licence file, no DOI and no archived snapshot; no release tag or retrieval date is given for the PoseBusters distribution; and DiffDock-L and EquiBind carry no commit hash or model-checkpoint identifier, unlike every other binary.

*Why it matters.* This is the rare defect that costs the candidate marks for work actually done. An examiner assessing reproducibility from the thesis alone would score it materially lower than the evidence warrants.

*Proportionate correction.* Rewrite Appendix A to state the coverage (assertions, tables, figure assets), name the notebook, and record the two figures that have no generator and why. Add a licence, deposit a tagged snapshot with a DOI, pin the repository commit, and add the two missing version identifiers.

*Effort.* **Rewriting, plus repository hygiene.**

---

### Minor issues

Grouped; each real and specific, none changing a conclusion. All verified live in the current build unless marked.

**Numerical and cross-reference.** Crystal-site non-reach count is 16 in the text (p.21) and 18 in the figure caption and panel-B sentence (p.23). Orai1 pooled yields printed to two decimals on a 120-pose denominator where one pose is 0.83 pp (p.30), inconsistently with Table 10's one decimal. §5.7 (p.45) quotes 71.1 % for DiffDock beside two per-unit means (the per-unit value is 70.9 %). §4.2.5 (p.38) reports a "56"-contrast family defined nowhere, in a sentence missing its noun ("*15 of 56 contact have adjusted p below 0.05*").

**Statistics.** "*It is significantly best on both axes*" (p.19) is asserted with no test, statistic or cross-reference, across cohorts of unequal and informative composition. "*gnina is decisive for EquiBind*" (p.15) has no within-family test and §5.2 states the opposite. Several Results p-values are unlabelled as to correction status, against the convention Appendix H.3 sets. Every Cliff's δ in the Orai1 chapter is a bare point estimate although Appendix H.6 declares bootstrap intervals. No Orai1 proportion carries a confidence interval although the calibration chapter uses Wilson intervals throughout. No minimum-detectable-effect statement anywhere in the Orai1 chapter, although H.9 supplies the template.

**Protocol reporting.** DiffDock's receptor-preparation recipe is never stated. Appendix C.1 (p.75) says the AutoDock gnina pass does not move geometry, contradicting C.5, §4.1.1 and the Literature Review. The exhaustiveness ladder claims exhaustiveness is "the only variable" while its own table footnote records a different thread allocation on one rung. The 1 Å box margin is never verified to contain every crystal ligand. The Orai1 pocket-detection result (p.104) names neither tool nor version and reports no numbers. The Orai1 histidine state is never given after the supplied CHARMM HSD/HSE/HSP assignments are discarded — relevant because His113 is load-bearing in the contact analysis.

**Chemistry and instrument caveats.** The PandaMap profiler cannot assign a charge to a nitrogen or oxygen bonded to carbon; since all three modulators are neutral with carbon-bonded heteroatoms, the reported absence of ionic, salt-bridge and attractive-charge contacts in the experimental panel is structurally forced rather than observed, yet is read as a panel chemistry difference (p.35). Water-mediated contacts cannot be recorded at all, which is not in C.10's otherwise candid defect list. The covalent contact class matches ligand hydrogens at H-bond distance. Appendix E describes cationic 2-APB as "a protonated fraction" where a primary aliphatic amine would be the majority species at pH 7.4. Appendix F.4 reads the source SDFs' neutral formal charges as a property of chemical space; 141 of 308 ligands carry an acid group written neutral.

**Presentation.** §4.2.1 (p.28) still sends the reader to Figure 27 "*with the transmembrane exclusion slab*"; that figure shows no slab, no ligands, no labels. Main-body captions are one-line and not self-contained. Declared naming conventions are not consistently honoured. Several abbreviations (RMSD, SE(3), STIM1, PDBQT) are used before or without expansion. §4.1 runs four pages before its first subsection, the landscape Table 4 forcing near-empty pages. The Orai1 control panel is defined first inside Results §4.2, after prior use in §3.2.6. Chapter 2 carries no Orai1 literature. Chapter 7 attributes the rigid-receptor limitation to AutoDock alone though the Methods state all three pipelines are rigid. No Acknowledgements section, and the supplier of the three Orai1 modulator geometries is never named.

---

## 10. Are the conclusions supported?

**Verdict: broadly supported. The one place where presentation outran the evidence has been repaired since 09-01.**

**RQ1 — supported, with the right hedges in place.** AutoDock ≥ DiffDock ≫ EquiBind on validity-aware recovery is established at top-5 and deeper (+13.2 pp, p_holm = 7.6e-4) and correctly declared unresolved at rank-1 (+2.6 pp, 95 % CI −4.2 to +9.4). The generation-versus-ranking decomposition is the strongest claim in the thesis and is well supported. §4.1.2 already states that AutoDock is the only tuned arm; it should add that the differential winner's curse is of the same order as the rank-1 gap.

**RQ2 — fully supported, and now more carefully stated.** Local optimisation lifts median per-complex validity from 13.3 % to 93.3 % (DiffDock) and 0 % to 80.0 % (EquiBind) with rank-biserial effects of 0.994–1.000, while fixed-rank near-nativeness moves 0.1–4.6 pp. §6.2 now adds that the AutoDock rank-1 re-ranking step "*does not clear Holm correction on its own, and the gain reaches significance only from top-5 onwards*" — a disclosure the 09-01 build lacked, prompted by the candidate's own table extension.

**RQ3 — now supported and internally consistent.** This was my strongest criticism of the 09-01 build and it has been properly resolved, not papered over. Every cost figure is on one charged basis defined in the Methods and printed in the figure subtitle; I reproduced the entire block from primary data. All three pairs separate under Holm, so "every pair passes Holm correction" is literally true for the first time. The basis-dependence is still disclosed (§4.3 p.40 notes the ordering depends on the denominator).

**RQ4 — partly a recommendation rather than a finding.** The staged workflow is a sensible synthesis, but was never run as a composite and its final component ("an endpoint-specific selector") is defined nowhere. It should be framed as a proposal.

**Most important overclaim.** Not in the Discussion, where I looked hardest — every candidate over-reading is hedged in place. It is Appendix C.8, p.86:

> "*For AutoDock Vina the mode count and the energy window govern only how many of the searched modes are written out, so the first ten modes of a thirty-mode run are the same ranked prefix that a ten-mode run reports and the twenty deeper experimental modes are discarded before analysis.*"

**Defensible replacement:**

> "For all three tools the top-ten cap is applied in the ranking each analysed variant carries, which for AutoDock Vina is the gnina CNNaffinity order rather than the Vina mode order. Sixty-six of the 120 analysed experimental AutoDock poses are Vina modes 11–30, which a ten-mode control run would not have written. The experimental panel therefore enters every comparison with the more heavily selected pose set for all three tools. That bias runs against the direction of the reported result, since the control panel returns the higher usable yield in every case and a deeper control pool would if anything have narrowed the gap."

---

## 11. Ten prioritised revision recommendations

1. **Regenerate the stale Orai1 material (M1)** — Table 12's AutoDock column, the Figure 26 caption, the C.8 folding decomposition, the C.10 audits — and correct body `:683` from 0.086 to the canonical 0.089. *Reanalysis over existing outputs.* Highest priority: one printed table is known not to reproduce from anything on disk.
2. **Correct the Orai1 panel-parity sentence (M2)** and add the Vina-rank-1-to-10 sensitivity. *Rewriting + one recomputation.*
3. **Carry Appendix H.4's four-percentage-point result into §5.2, §6.1 and the abstract (M5).** *Rewriting.* Cheapest high-value change available: it makes the thesis's most interesting negative result visible.
4. **Rewrite Appendix A to state the actual reproduction coverage (M8)**, name the notebook, and record the two figures with no generator. Add a licence, a DOI-bearing snapshot, and the DiffDock-L / EquiBind commit hashes. *Rewriting + repository hygiene.*
5. **Disclose the EquiBind harness and ligand-input non-parity in the body (M3)**, and report the native single-forward-pass baseline. *Rewriting + one tally.*
6. **Resolve the receptor-state contradiction (M4)** — add a per-frame breakdown with a leave-Fr0-out sensitivity, or delete the claim that frame-level tests are reported. *Either.*
7. **Add a one-page §3.3 "Docking protocol" to the printed Methods (M7)**, following the pattern §3.2.7 already sets for the cost basis. *Rewriting.*
8. **Report receptor structural completeness (M6)** — the PDBFixer configuration used, plus cohort counts of missing residues and incomplete side chains. *Rewriting + tally.*
9. **Add one illustrative structural figure** showing a recovered calibration complex (crystal ligand plus each tool's rank-1 pose, contact residues labelled) and one form-correct/placement-wrong near-miss; and **replace or relabel Figure 27** so the figure cited for the exclusion slab actually shows it. *Rewriting + two renders.*
10. **Clear the numerical minors** — the 16-vs-18 non-reach count, the undefined "56" family and its broken sentence, the two-decimal Orai percentages, the 71.1/70.9 basis mix, and the unlabelled p-values. *Rewriting.*

*(Deliberately **not** recommended: re-running with multiple seeds. It would strengthen the work and Chapter 7 already names the gap, but at MSc scale on a ~15-hour-per-arm campaign it is disproportionate. A 100-complex three-seed DiffDock spread, reported as a bound, is the proportionate version if time allows.)*

---

## 12. Five technically demanding viva questions

1. **On the parity defect.** Appendix C.8 argues the top-ten cap restores matching between the 30-mode experimental and 10-mode control panels for AutoDock. But the cap is applied in gnina CNNaffinity order, and 66 of the 120 analysed poses are Vina modes 11–30. Walk me through what that does to the 28.7-point usable-yield loss you call your strongest cross-panel result. Which way does the bias run, and can you bound it without re-docking?

2. **On the validity endpoint.** Appendix H.4 states the 98.9 %-versus-62.0 % spread "*narrows to roughly four percentage points once the pose is in the right place*", and that part of the physics-based advantage is constitutive because PoseBusters' intermolecular checks read the coordinate Vina minimises. Given both, what does PoseBusters validity measure that RMSD does not — and why should a reader of your abstract not conclude AutoDock produces physically better geometry than DiffDock?

3. **On the receptor model.** Your AutoDock receptors carry histidine protonated on both ring nitrogens in all 4,485 intact imidazole rings; the EquiBind receptors carry both nitrogens bare on 3,898 of 4,844. Vina's hydrogen-bond term is typed from the donor/acceptor set. Explain how this affects what each arm's objective sees at a histidine-containing site, whether it is differential across arms, and what you would run to bound the effect.

4. **On the cost fix.** You rebuilt the cost analysis on CPU-core-seconds / 32 + GPU-seconds after finding that Table 6's "charged" column was not one quantity. Defend dividing CPU-core-seconds by the thread count as a *time* rather than a resource-second, and tell me what that basis would do to a comparison run on a machine with a different core count or a second GPU.

5. **On the EquiBind harness.** EquiBind natively returns one pose from one forward pass. You report rank-1 and top-15 recovery for a 30-conformer pool built with a different conformer pipeline from the shared start conformer the other two arms receive, ranked by an external CNN. What exactly is the object your 18.2 % rank-1 figure describes, and how much of EquiBind's deficit can you attribute to the model rather than the harness?

---

## 13. Five broader conceptual viva questions

1. Your primary endpoint is a conjunction of PoseBusters validity and RMSD ≤ 2 Å. Both are geometric. What binding-relevant failure modes can a pose pass while still being wrong — and where do desolvation, conformational entropy and induced fit sit in that list?

2. You made the search blind and whole-protein to equalise information across three tools that localise the ligand differently. Argue the opposite case: that equalising *information* while leaving *difficulty* unequal produces a comparison that answers no question a practitioner has. What would a fairer design look like?

3. The post-2021 release criterion excludes these complexes as direct training entries but not their homologues, pockets or scaffolds. If you had to defend one claim about generalisation to genuinely novel pockets from this dataset alone, which would it be — and what would you have had to measure to defend more?

4. Your Orai1 exclusion slab is radially unbounded, so it removes not only the conduction pathway but the lipid-facing M3 interface that Appendix D.1 names as a proposed modulator site. Justify that. Then tell me what a positive site-plausibility criterion would look like for a membrane channel docked without its bilayer.

5. Rank-1 recovery is 36.6 %, top-15 is 65.3 %, and about six in ten rank-1 failures contain no qualifying pose at any depth. With a fixed compute budget to raise prospective rank-1 accuracy on this benchmark, where would you spend it — sampling, scoring, refinement or ensemble receptors — and what does your own data say about the expected return of each?

---

## 14. Final weighted score

**86.0 / 100**

| Category | Weight | Score | Contribution |
|---|---:|---:|---:|
| Research question and significance | 10 % | 4/5 | 8.0 |
| Literature and theoretical foundation | 10 % | 4/5 | 8.0 |
| Overall methodological design | 10 % | 4/5 | 8.0 |
| Docking protocol and technical execution | 20 % | 4.5/5 | 18.0 |
| Validation and benchmarking | 15 % | 5/5 | 15.0 |
| Results and data analysis | 10 % | 4/5 | 8.0 |
| Interpretation, discussion, limitations | 10 % | 4/5 | 8.0 |
| Reproducibility and transparency | 10 % | 4.5/5 | 9.0 |
| Writing and presentation | 5 % | 4/5 | 4.0 |
| **Sum** | **100 %** | | **86.0** |

Weights verified to sum to 100 %; contributions verified to sum to 86.0.

## 15. Overall classification

**Very strong (80–89).** Deliberately not converted to an institutional grade, as no formal scale was supplied.

## 16. Confidence in this assessment

**High.** Basis: I read the complete main body and the load-bearing appendix sections in the current source; I re-derived Table 4 and the four headline recovery figures from `per_pose_metrics.csv` (all reproduce exactly); I recomputed the entire revised cost block from `docking_effort_gnina_v2_charged/per_complex_effort.csv` and confirmed every printed value including the Friedman statistic and Kendall W; I re-verified the 66-of-120 Vina-rank finding from `posebusters_filtered_results.csv`; I read the regenerated Figure 15 image directly; and I audited float and citation integrity programmatically. Each of my eight majors was re-tested against the current text rather than carried forward. Confidence is lowest on the precise scope of M1, where I rely partly on the candidate's own coverage audit for the Table 12 cell-by-cell comparison, though I confirmed the disputed values are still printed.

## 17. Examiner recommendation

This is a very strong MSc thesis that should pass without reservation, subject to minor corrections. Its distinguishing quality is intellectual honesty operating at a level rarely seen at this stage, and the two days since my previous reading demonstrate it in action rather than merely in prose: presented with a figure that contradicted its own caption, the candidate did not reword the sentence but rebuilt the entire cost analysis on a single defensible currency, then verified that the repaired claim was true; presented with a request to extend two tables, he discovered that the extension would flip an unrelated between-tool verdict under the pooled family, and split the family rather than accept the convenient result; and an independent coverage pass found that one appendix table had been tested against the wrong arm since its first commit, concluding — correctly — that the thesis was right and the script was wrong. The defects that remain are of a specific and benign kind: a superseded analysis generation still reaching print, a parity argument stated more strongly than the data support, an analysis whose most interesting negative result never travels from the appendix to the abstract, and a reproducibility statement that undersells work already done. None invalidates a conclusion; all are correctable by rewriting and re-running scripts over data already held, and the new assertion harness makes the corrections verifiable. I would require M1, M2 and M5 before final deposit, recommend M3, M4, M6, M7 and M8, and defer the rest to the candidate's judgement.

---

## Closing question

> *Does this thesis demonstrate that the student can independently plan, execute, critically evaluate, and communicate a scientifically sound protein–ligand docking project at master's level?*

### **Yes, with minor reservations.**

Planning shows in a design that anticipates its own confounds — blind whole-protein search for information parity, an exhaustiveness ladder, refiner arms, a control panel, four receptor states. Execution spans 303 complexes across three pipelines with version-pinned tooling, declared source patches and reported seeds. Critical evaluation is the standout and is demonstrably ongoing: the candidate computes the power his primary contrast lacks, quantifies his own selection optimism, split a multiplicity family after finding it could flip a verdict in his favour, and rebuilt his cost analysis rather than reword a claim it did not support. Communication is precise and the abstract does not outrun the body. The reservations are that a superseded analysis generation still reaches print, one parity argument is stated more strongly than the data support, and the printed Methods carries no protocol — defects of upkeep and placement, not of scientific judgement.
