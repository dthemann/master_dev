# Examiner Review — MSc Thesis

**Title.** *Whole-Protein Docking with AutoDock Vina, DiffDock and EquiBind: Physical Validity and Near-Nativeness Across Ranking Depth*
**Degree course.** Artificial Intelligence Engineering (FH Master)
**Documents assessed.** `Thesis_short.tex` (abstract, Kurzfassung, front/back matter), `body_main_short.tex` (Ch. 1–7), `body_appendix_short.tex` (App. A–H), `Literatur.bib`; printed page numbers resolved from `.toc/.lof/.lot`. 140 pages; 6 main tables + 22 appendix tables; 16 main figures + 23 appendix figures.
**Review date.** 2026-09-09
**Standard applied.** Strict but proportionate MSc standard. Not assessed as a PhD dissertation.

---

## 1. Executive assessment

The thesis asks how three configured blind whole-protein docking workflows — AutoDock Vina with gnina rescoring, DiffDock-L with smina minimisation, and EquiBind with gnina optimisation — differ in physical validity, near-nativeness, ranking depth, cost and blind-docking behaviour (RQ1–RQ4, `body_main_short.tex:11-16`). The approach is a 303-complex crystal-referenced calibration benchmark under a compound endpoint that ANDs PoseBusters validity with heavy-atom RMSD ≤ 2 Å, evaluated at rank-1, top-5, top-15 and top-30, plus a reference-free Orai1 case study (3 modulators × 4 MD snapshots, with a 308-ligand control panel).

The principal contribution is not a tool ranking but a **decomposition of failure**: the thesis separates generation from ranking, placement from internal form, and geometric accuracy from physical admissibility, and it shows that generation dominates — about six to seven in ten rank-1 failures of the two leading tools contain no qualifying pose at any depth (`:433`, `:811`). That is a genuinely useful and non-obvious result.

The strongest aspect is **evidential discipline**. The document consistently reports against its own interest: the rank-1 lead is declared unresolved with its interval and a power calculation (`:405`, App. H.9); the alternative CNN ranking head that would have given a better number is printed (152 vs the reported 149, App. C.1); the reference-convention sensitivity arm is tabulated in full (Table 23); the pocket-guided EquiBind arms are reported "because the answer they give is negative" (App. C.4.1); the selection-on-the-evaluation-set optimism is quantified rather than merely confessed (`:889`).

The most important weakness is structural and has one root cause: **Chapter 3 "Methods" contains no docking-protocol section at all.** It has exactly two sections — Datasets and Evaluation Metrics — so receptor preparation, ligand preparation, search settings, exhaustiveness, refiner configuration and pose-selection rules are absent from the main body and reachable only through Appendix C. Everything downstream inherits that displacement: a **body/appendix asymmetry of disclosure** in which qualifications that materially change how a result should be read live only in the appendix or only in Ch. 7, while Ch. 4, Ch. 5 and the abstract state the unqualified version. The ligand-input non-parity between arms, the magnitude of the nearest-copy reference rule, and the metal-stripping stratification are the three clearest cases. Layered on this are a small number of **hard textual defects** — one claim inverted in sign against the thesis's own measurement (`:429`), one factually false premise in the RQ3 conclusion (`:871`), and one Discussion overclaim contradicted by the Results two chapters earlier (`:833` vs `:718`).

Conclusions are **broadly supported but overstated in three identifiable places**, all correctable by rewriting. The headline calibration result itself is sound and has survived more sensitivity testing than most published docking papers receive.

Quality level: **very strong for an MSc**, with appendix-level rigour that in places exceeds the published-paper norm and a main body that does not always carry that rigour forward.

Confidence in this assessment: **high** for the calibration chapter, the protocol audit, the statistics and the reproducibility assessment; **medium** for the Orai1 chapter, where the receptor model's provenance gap limits what any reader can verify.

---

## 2. Summary of the research approach

| Element | As reported |
|---|---|
| Design | Exploratory/comparative, retrospective, crystal-referenced calibration + reference-free application |
| Calibration set | PoseBusters benchmark, 308 post-2021 PDB complexes; 303 analysed after 5 DiffDock production failures (App. H.1) |
| Application set | Orai1 homology model (Hou 2012 *Drosophila* architecture, human numbering), 4 MD snapshots Fr0/Fr300/Fr400/Fr499; 3 modulators (2-APB as neutral 2abp-NH2, GSK-7975A, Synta-66) |
| Control | 308 calibration ligands × 4 frames = 1,232 units (1,215 for DiffDock), decoy-style background comparator; no enrichment claim made (`:564`) |
| Search | Blind whole-protein for all three arms; Vina box built per receptor from the atom bounding box + 1 Å, median edge 63.5 Å, median volume 263,874 Å³ (App. C.1) |
| Primary endpoint | PB-valid (22 reference-free checks) **AND** symmetry-corrected heavy-atom RMSD ≤ 2 Å to the **nearest deposited copy**, in place, no superposition |
| Secondary endpoints | Form (Kabsch RMSD ≤ 1 Å), 8 Å centroid clustering + 4 Å crystal-site reach, PandaMap interaction fingerprints, charged device-occupancy cost, Orai1 operational placement (R91→E106 exclusion slab) |
| Statistics | Distribution-free throughout; exact McNemar, Wilcoxon, Cochran's Q, Friedman, Mann–Whitney/Cliff's δ, Fisher; Newcombe paired and Wilson intervals; declared Holm and Benjamini–Hochberg families (App. H.3); inferential unit = complex (benchmark) / frame–ligand pair (Orai1) |
| Headline result | Rank-1 valid-near-native 49.2 % / 43.6 % / 18.8 %; top-15 70.3 % / 59.1 % / 27.4 %; AutoDock–DiffDock rank-1 gap +5.6 pp (−1.7 to +12.8, unresolved), resolved from top-5 onwards |

---

## 3. Five principal strengths

1. **The primary endpoint is defensible and multiply sensitivity-tested.** Requiring validity *and* accuracy jointly is the right response to the failure mode the source benchmark documents, and the thesis then tests the endpoint against itself: ten distance thresholds (Table 27), both reference conventions at four depths (Table 23), an intention-to-treat arm restoring the five excluded complexes (App. H.1), and a decomposition showing the validity term removes at most five of 303 complexes (App. H.4). Very few MSc theses interrogate their own endpoint this hard.

2. **Protocol reporting is at or above publication standard.** Exact ADFRsuite invocation (`prepare_receptor -A hydrogens -U nphs_lps_waters_nonstdres`), Vina seed 42, energy range 6 kcal/mol, 30 modes, per-receptor box construction rule with the full edge and volume distribution, gnina 1.3.2 pinned to `crossdock_default2018`, CNNaffinity rather than default CNNscore ordering with the alternative reported, `--minimize` for all three refined arms, 20 denoising steps and 30 samples for DiffDock-L, 30 ETKDGv3+MMFF conformers for EquiBind, and versions "verified by invoking each binary" (App. C).

3. **Negative, unflattering and against-interest results are reported in full.** The pocket-guided EquiBind arms (worse than unguided, App. C.4.1); the raw exhaustiveness-128 search beating its own rescored sibling over the whole pool (217 vs 214); the CNNscore head that would have scored better than the one reported; the Orai1 prefix counterfactual that removes the study's one surviving cross-panel significance (App. C.8); and the winner's-curse correction on its own selection (`:889`).

4. **The scope of inference is policed, and policed consistently.** The scope-limiting sentence at `:18` is restated in the Discussion (`:793`), applied reflexively to the AutoDock arm ("Physics-based is itself an operational label here", `:796`), and closed in Ch. 7: "no score used for ranking was compared against measured affinities" (`:893`). The docking-score/affinity boundary is never crossed.

5. **Statistical craftsmanship.** Declared multiplicity families with the design defended rather than asserted (App. H.3); the inferential unit fixed by design and the Orai1 units explicitly declared correlated with an effective n of three ligands (`:187`); correct refusal to read a non-significant rank-1 contrast as equivalence, backed by a power calculation showing the minimum detectable difference is 10.5 pp against an observed 5.6 pp (App. H.9).

---

## 4. Five principal weaknesses

1. **The Methods chapter has no docking-protocol section** (M0). Ch. 3 runs Datasets → Evaluation Metrics and stops; no receptor prep, ligand prep, search setting, exhaustiveness, refiner or pose-selection rule appears in the main body. This is the root cause of weakness 2 and of much of M4, M6 and M7.

2. **Material qualifications are stranded in the appendix while the body states the unqualified result** — ligand-input non-parity (M4), nearest-copy magnitude (M7), metal-stratum heterogeneity (M6).

3. **One claim runs opposite in sign to the thesis's own measurement** (`:429`, M1). This is the single hardest textual defect in the document.

4. **The Literature Review does not support the study's own design.** PoseBusters (`ref035`), the source of the benchmark and of the entire validity axis, is never cited in Ch. 2; blind whole-protein docking — the defining design and the explicit subject of RQ4 — is never reviewed there either; and no physical account of molecular recognition (enthalpy/entropy, desolvation, induced fit, a non-covalent interaction taxonomy) is given before Ch. 3 and Ch. 4 begin reasoning with those concepts.

5. **The Orai1 chapter rests on a receptor whose provenance and quality cannot be established.** No modelling program, template PDB IDs, force field, lipid composition, ion treatment or trajectory length is recorded (App. D); the one hard piece of evidence about geometric quality — Cα–Cβ bonds stretched to 1.85–2.06 Å at Ser93 in every subunit of Fr0 — was treated as a tooling obstacle and never swept across Fr300/Fr400/Fr499, which DiffDock and EquiBind read as delivered.

6. **Reproducibility of the numbers, as opposed to the procedure, is not achievable by a third party** — and for the DiffDock arm not by the author either. Code and data are "available from the author on request" with no licence, no persistent identifier and no commit; the benchmark DiffDock run is "one unseeded draw made before the local seed patch" (`:887`); no arm was replicated across seeds.

---

## 5. Detailed scoring table

Scale: 0 absent/invalid · 1 very weak · 2 below MSc standard · 3 competent MSc · 4 strong · 5 excellent/publication-quality for MSc.

| # | Category | Weight | Score | Weighted | Strongest supporting evidence | Most important deficiency | Confidence |
|---|---|---:|---:|---:|---|---|---|
| 1 | Research question and significance | 10 % | 4.0 | 0.40 | RQ1–RQ4 explicit and numbered (`:11-16`); each answered in a dedicated Conclusions section whose numbers reproduce Table 2 exactly; scope-limiting sentence at `:18` honoured throughout | The Motivation's third claimed gap is contradicted by the thesis's own Table 7 and §5.6 (M11); no RQ covers the Orai1 arm, which is half the Results | high |
| 2 | Literature and theoretical foundation | 10 % | 3.5 | 0.35 | Vina and Vinardo scoring functions described exactly, with weights (App. B.1.2); DiffDock's two learned components correctly separated; pose-prediction / ranking / affinity distinction held from Ch. 2 to Ch. 6 | `ref035` never cited in Ch. 2; blind whole-protein docking never reviewed; no physical foundation for molecular recognition; Orai1 pharmacology absent from Ch. 2; Table 7 lists four prior studies without stating what any concluded | high |
| 3 | Overall methodological design | 10 % | 4.0 | 0.40 | Cohort reduction 308→303 documented and tested by an intention-to-treat arm; control panel adopts decoy practice and explicitly declines an enrichment claim; selection optimism quantified per family | Selection and evaluation on the same 303 complexes with no held-out split, and the tuning budget is asymmetric — an eight-arm AutoDock ladder against three refiner arms per learned tool (disclosed, `:889`, but not neutralised) | high |
| 4 | Docking protocol and technical execution | 20 % | 4.0 | 0.80 | Item-level protocol disclosure (§3 above); protonation model reported residue class by residue class with population counts; Fr0 non-parity quantified atom by atom; refiner-mode parity actively repaired by re-running EquiBind under `--minimize` | Ligand input not matched across arms and the body never says so (M4); no titration model anywhere and all 4,485 intact histidines are doubly protonated and cationic, with no sensitivity check; metals/cofactors stripped for two arms but retained for the third | high |
| 5 | Validation and benchmarking | 15 % | 4.0 | 0.60 | 303-complex crystal-referenced validation; Table 23 reference-convention sensitivity; ten-threshold resolution (Table 27); metal-stratified recovery; exhaustiveness ladder with measured cost | Nearest-copy convention supplies 38 of AutoDock's 149 and 29 of DiffDock's 132 rank-1 recoveries and doubles the rank-1 lead, and site-equivalence of the crediting copies is never established (M7); no cross-docking arm bounds the regime the Orai1 application actually faces; deposited crystal poses never passed through the study's own battery as a positive control | high |
| 6 | Results and data analysis | 10 % | 4.5 | 0.45 | Distributions rather than headline points throughout; per-tool per-rank denominators printed (298/282/266 at k=1) with the overlapping-cohort caveat; declared multiplicity families; correct refusal to resolve the underpowered rank-1 contrast | The same three p-values are labelled `p_holm` in §4.1.2 and "exact McNemar" in App. H.9; Ch. 5 introduces three quantitative counts with no printed source anywhere (M5) | high |
| 7 | Interpretation, discussion and limitations | 10 % | 4.0 | 0.40 | Ch. 7 is unusually complete — stochasticity, selection optimism, denominators, leakage, receptor-axis mismatch, profiler defects, slab arbitrariness, cost conditionality all stated | The §5.5 filter/TM1–TM3 "recovery" sentence is contradicted by §4.2.5 (M3); the Orai1 reversal is attributed causally on a comparator the appendix says cannot be measured (M9); the metal stratification never reaches Ch. 5 or Ch. 6 (M6) | high |
| 8 | Reproducibility and transparency | 10 % | 4.0 | 0.40 | Versions verified by binary invocation; seeds, box rule, flags, CNN ensemble and hardware all printed; a regeneration guide exists and names what it could not reconstruct | Repository private and "on request"; no licence, DOI or commit; the guide itself is inside the private repo and admits unnamed tables with unaudited provenance; the primary DiffDock arm is unseeded and unrepeatable; Orai1 receptor frames have no named source or availability | high |
| 9 | Writing and presentation | 5 % | 4.0 | 0.20 | Disciplined, concise register; abstract's quantitative claims trace to the body; German Kurzfassung numerically faithful; tables carry substantive footnotes | Eleven of sixteen main-body figure captions are bare titles with no panel key, axis definition, units or cohort size; one wrong figure cross-reference (§4.2.1 cites Fig. 27 for the slab; the slab is Fig. 22); `2-APB` in prose vs `2abp-NH2` in every table | high |

**Weighted total.**
0.40 + 0.35 + 0.40 + 0.80 + 0.60 + 0.45 + 0.40 + 0.40 + 0.20 = **4.00 / 5**
Weights check: 10 + 10 + 10 + 20 + 15 + 10 + 10 + 10 + 5 = 100 %.
**4.00 × 20 = 80.0 / 100.**

---

## 6. Docking-protocol audit

Legend: **A** adequately reported · **J** reported but inadequately justified · **M** missing.

### 6.1 Protein preparation

| Item | Status | Evidence |
|---|---|---|
| Crystallographic waters | **A** | Absent from the benchmark distribution itself, so every arm docks a water-free receptor (App. C) |
| Cofactors, ions, metals | **J** | "PDBFixer, which strips waters, ions, metals and cofactors" (App. C.1). Consequence *measured* (78/303 metal-adjacent; rank-1 51.1/49.3 % metal-free vs 43.6/26.9 % metal-adjacent) but the *choice* is never justified, and the source benchmark ships a cofactor-retaining file that was not used |
| Hydrogen addition | **A** | ADFRsuite `prepare_receptor -A hydrogens`, non-polar merged |
| Protonation assignment | **J** | Explicitly no titration model, no pKa estimate, no pH argument. All 4,485 intact histidines doubly protonated and cationic; Asp/Glu bare carboxylates; Lys/Arg charged; 316 of 2,282 Cys carry an inherited thiol H. Uniform across the cohort, but non-physical and untested for sensitivity |
| Missing atoms / residues / loops | **M** | No statement anywhere on unmodelled loops or incomplete side chains in the 308 benchmark receptors |
| Chain selection | **A** | Not applicable — the whole cleaned receptor is used |
| Repair / minimisation | **A** (with a parity break) | Only the AutoDock copy of Fr0: PDBFixer + OpenMM amber14 / GBn2, quantified at 0.82 Å heavy-atom and 0.64 Å Cα, 149 atoms > 2 Å (App. D). Note: implicit *aqueous* solvent applied to a bilayer-stripped membrane channel, which is not acknowledged as inappropriate |
| Alternate conformations / occupancy | **A** | Resolved with counts rather than assertion (App. C) |
| Biological assembly vs ASU | **A**, but see M7 | The whole deposited assembly is retained; multi-copy handling is the nearest-copy rule |
| Orai1 homology model | **M** on provenance | No modelling program, template PDB IDs, force field, lipid composition, ion treatment or trajectory length recorded; the gap is stated honestly and nothing is inferred (App. D) — but it is never carried into Ch. 7 |

### 6.2 Ligand preparation

| Item | Status | Evidence |
|---|---|---|
| Structure verification | **A** | InChIKeys verified against published identifiers for all three Orai1 compounds, including the neutral-phenol vs phenolate distinction |
| Stereochemistry | **A** | Demonstrated rather than assumed: the four reference-dependent identity checks (formula, bond graph, tetrahedral chirality, double-bond stereo) run post hoc and remove nothing on any arm |
| Protonation / tautomers / protomers | **J** | No pKa, tautomer or protomer treatment anywhere; the three Orai1 compounds "docked as neutral species" as a stated modelling choice. Disclosed only in App. C and never listed as a limitation |
| Salt removal | **A** | Benchmark ligands supplied as single entities |
| Conformer generation | **A** in appendix, **M** in body | RDKit ETKDGv3, UFF (AutoDock/DiffDock, one shared conformer) vs MMFF (EquiBind, thirty). See M4 |
| Partial charges | **A** | Gasteiger, and correctly reasoned inert — neither Vina nor Vinardo carries a charge-dependent term |
| Macrocycles / covalent / metals in ligand | **M** | Never mentioned, although two of three arms cannot alter a ring conformation they are handed |
| Boron | **A** | Stock Vina defines no boron type; the added type, radius (1.92 Å) and well depth are documented at parameter level |
| Multiple relevant states | **M** | Only one state per ligand considered |
| Orai1 ligand geometry provenance | **M** | "lacked records of both the optimisation program and its level of theory" (`:97`) |

### 6.3 Binding-site / search-space definition

| Item | Status | Evidence |
|---|---|---|
| Site rationale | **A** | Blind by design, so that Vina "faces the same problem as the two methods that localise the ligand themselves" |
| Co-crystallised ligand used to place the box | **A** — deliberately *not* used | Explicitly contrasted with the source benchmark's 25 Å cube (15,625 Å³) against this study's median 263,874 Å³ |
| Grid coordinates and dimensions | **A** | Reproducible construction rule + full empirical distribution (edges 29.5–167.4 Å, volumes 44,430–2,479,703 Å³, edge-ratio median 1.29); Orai1 per-frame boxes 91.8–102.0 Å per edge; 0.375 Å grid spacing |
| Box adequacy | **A** | The dilution penalty is named, cited (CASF-2016 whole-protein loss) and answered empirically with the exhaustiveness ladder |
| Multiple / allosteric sites | **A** | Handled by clustering and by the nearest-copy rule; the Orai1 slab encodes an explicit outer-pore hypothesis |
| Orai1 slab geometry | **A** on geometry, **J** on scope | Pore axis, axial coordinate and region fully specified and rebuilt per frame; but radially unbounded — see M9 |
| Orai1 box volume / search adequacy | **M** | Orai1 boxes are ~3.2–3.6× the benchmark median volume at the identical exhaustiveness 128, and no convergence evidence is offered for that regime |

### 6.4 Docking settings and pose selection

| Item | Status | Evidence |
|---|---|---|
| Software and versions | **A** (mostly) | Verified by binary invocation; gnina 1.3.2, fpocket 4.0, P2Rank 2.5, PandaMap 4.1.0 named. **Gap:** the DiffDock-L and EquiBind trained checkpoints are identified by model family and date rather than by release or checksum |
| Engine and scoring function per arm | **A** | Stated per arm, including which score orders which pool |
| Exhaustiveness / samples / steps | **A** | Vina exh 128 (ladder 18–128); DiffDock 20 steps, 30 samples; EquiBind 30 conformers |
| Seeds | **Partial** | Vina seed 42 and gnina seed 42 stated; the benchmark DiffDock run is unseeded (`:887`); Orai1 experimental panel uses seed 42, control panel unseeded; resampling "a fixed seed" without a value |
| Repeat runs / stochastic variance | **M** | "No arm was assessed across multiple seeds" (`:888`) |
| Receptor rigidity / flexible residues | **A** | All three rigid; flexibility represented only by the four Orai1 states, with the limits of that stated |
| Symmetry / equivalent poses | **A** | Graph-isomorphism mapping; the residual conjugated-terminal-group limitation is disclosed *and bounded* ("adds two recovered complexes across the three tools and removes none") |
| Final-pose criteria | **A**, but asymmetric | AutoDock re-ordered by gnina CNNaffinity; DiffDock keeps confidence order on refined coordinates; EquiBind ordered by gnina affinity. Asymmetry stated at `:203` but its magnitude (the re-ranker is worth +6.3 pp at rank-1, larger than the +5.6 pp lead it wins) is only in App. C.1 |
| Rescoring / consensus / clustering | **A** | Fully specified, including the pose-redundancy question left open |
| Visual inspection | **A** — none enters the evidence chain | Both molecular renders carry captions declaring they carry no measurement |
| Tuning on the evaluation set | **A** (disclosed) | `:889`, App. C.2; not neutralised by a split |

**Verdict on the protocol.** Reporting is exceptional; execution carries four avoidable and unexamined weak points — no titration model, uniformly cationic histidines, asymmetric receptor composition across arms, and unmatched ligand inputs. None invalidates the calibration result; all are correctable by re-analysis of data already held.

---

## 7. Validation assessment

**Pose-prediction validation: adequate, in places unusually strong for an MSc.**
The study re-docks 303 co-crystallised complexes with the ligand removed, scores symmetry-corrected heavy-atom RMSD against the deposited pose, and pairs it with a 22-check physical-validity battery. Atom mapping is by graph isomorphism with the residual limitation disclosed and quantified. The 2 Å threshold is conventional and cited; the 1 Å Kabsch "Form" threshold is explicitly declared as the author's own choice with no published convention claimed — and, importantly, Table 27 resolves it over ten thresholds so the reader can test the choice. Cross-docking into different receptor conformations is **not** performed on the calibration set; the appendix states plainly that the cognate case is the easier regime and that the calibration figures therefore *bound* the Orai1 expectation from above rather than predicting it. That is the right framing, but a cross-docking arm would have measured the penalty rather than assumed it.

**Two qualifications on the primary endpoint.**
(i) The nearest-copy rule is a minimum over deposited copies and can only lower RMSD. It is justified by appeal to the source benchmark's own convention, and it is tested (Table 23) — but it supplies 38 of AutoDock's 149 and 29 of DiffDock's 132 rank-1 recoveries and doubles the rank-1 lead from +2.6 to +5.6 pp, and the crediting copies are never classified for binding-site equivalence. (ii) The deposited crystal poses were never passed through the study's own battery, so the achievable validity ceiling on these prepared receptors is unmeasured — a cheap positive control that would have strengthened every validity comparison in the thesis.

**Virtual-screening validation: absent, and appropriately so.**
No enrichment factor, ROC-AUC, BEDROC or early-recognition metric is reported. This is **not** a deficiency: none of RQ1–RQ4 is a screening question, the Orai1 control panel is explicitly described as "only a background comparator" with "no enrichment claim", and demanding enrichment statistics here would be demanding an experiment the thesis did not set out to run.

**Ranking and affinity claims: correctly scoped.**
No correlation with measured affinity is attempted or implied, and Ch. 7 closes the question explicitly. Rank ordering and absolute affinity are kept distinct throughout — including in the treatment of gnina's CNNaffinity head, which is used and described as a *ranking* device rather than an affinity estimate.

**Classification: adequate, approaching unusually strong for an MSc**, held back from the higher grade by the absent cross-docking arm, the absent crystal-pose positive control, and the un-replicated stochastic arms.

**Practical consequence.** The evidence supports: (a) comparative statements about *these configured pipelines* on *this benchmark* regarding validity, near-nativeness and ranking depth; (b) the generation-versus-ranking decomposition, which is robust across depths, thresholds and both reference conventions; (c) cost statements conditional on the declared occupancy basis. It does **not** support: a physics-vs-AI method-class verdict (the thesis says so itself), transfer of the ordering to unseen targets (selection optimism, no held-out split), or any binding-mode, site or activity claim for the three Orai1 modulators. The Orai1 arm supports **hypothesis generation and workflow comparison only** — which is precisely what the thesis claims for it.

---

## 8. Docking-specific red-flag checklist

| # | Red flag | Status | Evidence |
|---|---|---|---|
| 1 | Docking scores treated as experimental binding energies | **Absent** | `:893` "no score used for ranking was compared against measured affinities"; `:806` PoseBusters validity "establishes neither favourable binding free energy…" |
| 2 | Small score differences presented as decisive | **Absent** | The 5.6-pp rank-1 gap is declared unresolved with its interval and a power calculation; the one-complex DiffDock refiner margin is called "unstable" at `:201` |
| 3 | Protocol not validated | **Absent** | 303 crystal-referenced complexes, four depths, ten thresholds, two reference conventions |
| 4 | Receptor selected without rationale | **Absent** | App. D §"Choice of the Four Receptor States" argues the case structurally and measures the spread the four states span |
| 5 | Protonation/tautomerism/stereochemistry ignored | **Absent** (but simplified) | Not ignored — reported residue class by residue class, with stereochemistry survival demonstrated. Simplified: no titration model, no tautomer/protomer enumeration |
| 6 | Waters/metals/cofactors removed without justification | **PRESENT** | Removal documented and its cost measured, but never *justified*; the source benchmark's cofactor-retaining file was not used, and App. C.9 concedes the difference |
| 7 | Binding site / grid insufficiently documented | **Absent** | Construction rule + full empirical distribution + per-frame Orai1 boxes + slab equations |
| 8 | Only best-looking poses/compounds reported | **Absent** | Table 1 covers up to 30 poses per complex for every variant; negative arms reported in full |
| 9 | Visual inspection without transparent criteria | **Absent** | Both renders carry captions stating they carry no measurement |
| 10 | Favourable pose treated as proof of binding | **Absent** | `:899` "The Orai1 poses remain computational hypotheses" |
| 11 | Docking treated as proof of biological activity/mechanism | **Absent** | `:834` "Neither partial contacts nor their combination validates a complete pose" |
| 12 | Experimental validation implied but not performed | **Absent** | Required assays named explicitly as future work, located at JKU |
| 13 | Ligand/receptor preparation not reproducible | **PRESENT** | App. A names the Fr0 receptor copy and the prepared ligand files among items the regeneration guide "could not reconstruct"; Orai1 ligand geometry and receptor model provenance both unrecorded |
| 14 | Same data used to tune and validate | **PRESENT** | `:889` "Variant selection and evaluation use the same 303 complexes, with no held-out split"; disclosed, quantified (winner's curse 0.4–1.0 pp per family), not neutralised |
| 15 | Conclusions exceed what docking can establish | **PRESENT, localised** | §5.5 `:833` filter/TM1–TM3 "recovery" (M3) and §6.1 `:859` causal attribution of the Orai1 reversal (M9). Three sentences out of a consistently careful document |

Eleven of fifteen absent; all four present items are disclosed by the candidate rather than concealed. This is a very clean checklist by MSc standards.

---

## 9. Critical, major and minor issues

### Critical issues
**None.** No defect identified in this review invalidates the main calibration results or the conclusions drawn from them. The headline ordering survives both reference conventions, the intention-to-treat cohort, ten distance thresholds and four ranking depths.

### Major issues

**M0 — Chapter 3 "Methods" contains no docking-protocol section.**
*Where.* Ch. 3, printed pp. 5–11, `body_main_short.tex:74-189`. The chapter has exactly two sections: `\section{Datasets}` (`:79`) and `\section{Evaluation Metrics}` (`:102`).
*Why it matters.* Receptor preparation, ligand preparation, search-space construction, exhaustiveness, refiner mode and pose-selection rules — the entire docking protocol — appear nowhere in the main body. A reader of Chapters 1–7 cannot state what was actually run. Everything the protocol audit in §6 rates as "adequately reported" is reported in Appendix C, and Ch. 3 does not signpost it for preparation at all. This single structural gap is the root cause of M4 (ligand non-parity invisible in the body), and it explains why M6 and M7 read as concealment when the underlying analyses were in fact done and published.
*Evidence.* Grepping `body_main_short.tex:73-190` for `exhaustiv|gnina|smina|refin` returns three incidental hits and no protocol statement; §4.1 (`:195-203`) then has to open the Results chapter with nine paragraphs of variant-selection methodology that the Methods should have carried.
*Correction.* Add a §3.2 "Docking protocol" of roughly two pages summarising, per arm, receptor preparation, ligand preparation, search space, engine settings, refiner mode and ranking key, with a table pointer to Appendix C Tables 9 and 11; renumber Evaluation Metrics to §3.3; move the §4.1 selection methodology into it. *Effort: rewriting (1 day); no new analysis.*

**M1 — A Results claim is inverted in sign against the thesis's own measurement.**
*Where.* §4.1.2, printed p. 17, `body_main_short.tex:429`: "Refined-score ordering can be considered an unclaimed remedy for DiffDock, worth between three and four points of rank-1 near-native rate on the reported coordinates before the validity requirement is applied."
*Why it matters.* It is placed at exactly the point where the AutoDock–DiffDock rank-1 contrast is unresolved, so it implies DiffDock has recoverable headroom that would close the gap.
*Evidence.* The chapter's own footnote at `:201` measures 41.6 % under smina ordering and 41.6 % under gnina against 43.6 % under confidence ranking "on those same coordinates" — a **2.0-point loss**. Appendix G (`:964`) repeats the measurement and states outright that re-ordering "does not help". No number anywhere in the thesis supports "three to four points".
*Correction.* Delete the sentence, or replace it with: "Re-ordering DiffDock's poses by their refined score does not help: it nominates a near-native rank-1 pose for 41.6 % of complexes under either refiner against 43.6 % under the native confidence order, so the confidence ranking is retained." *Effort: rewriting.*

**M2 — The RQ3 cost defence rests on a factually false premise.**
*Where.* §6.3, printed p. 48, `body_main_short.tex:871`: "AutoDock Vina and DiffDock each produced a single configuration, and for them the charged cost and the whole run are the same number."
*Why it matters.* The sentence exists to justify charging every tool only for the configuration carried forward, in the section that concludes AutoDock is the most expensive route. It is false for AutoDock and misleading for DiffDock, and it makes the accounting rule look fairer than the campaign was.
*Evidence.* App. C.1 (`:169`): "Eight AutoDock variants were run on this set, five independent docking passes"; Table 8 records those five passes at 2.14 + 2.56 + 3.14 + 4.84 + 4.73 = **17.41 Vina search hours**, of which only 4.73 h is charged. DiffDock's single raw run underlies three arms (Table 1). The thesis's own Ch. 7 (`:889`) calls them "the eight configurations of Table 8 against three refiner arms for each learned pipeline". The EquiBind half of the same sentence is also wrong: six of its nine configurations came from two further cropped-receptor inference runs, not from one pass (App. C.4.1).
*Correction.* Replace with an accurate statement — all three pipelines discarded uncharged selection work; report the campaign totals alongside the charged ones. *Effort: rewriting (the numbers already exist).*

**M3 — The single most important overclaim, in §5.5.**
*Where.* §5.5, printed p. 45, `body_main_short.tex:833`: "Recovery of literature-associated filter and TM1/TM3 residues therefore generates hypotheses."
*Why it matters.* It is the one place where the reference-free chapter reads as partial confirmation of published pharmacology.
*Evidence.* §4.2.5 (`:718`) says the opposite about the filter: "Each tool produces one modulator pose contacting Glu106, so neither consistently recovers the filter residue." And the residues actually recovered — His113, Tyr115, Pro201, Leu202 — are among the most contacted residues on the **control panel of 308 ligands with no known Orai1 activity**, with "no residue differ[ing] after Holm correction of the fifteen within-tool tests" (`:708`).
*Proposed replacement (paste-ready).* "Single poses from each tool contact the selectivity-filter residue Glu106, and the loop residues His113, Tyr115, Pro201 and Leu202 recur in the modulator poses. Neither observation is discriminative: the same residues are contacted at a statistically indistinguishable frequency by 308 control ligands with no known Orai1 activity, so these contacts reflect the accessible surface of the outer pore rather than modulator-specific recognition. They are therefore not evidence for the published sites."
*Effort: rewriting.*

**M4 — The three arms did not receive the same ligand, and the body never says so.**
*Where.* Disclosed at App. C (`:148`), App. C.4 (`:258`) and Table 9 row "Start conformer"; **absent** from Ch. 3, Ch. 4, Ch. 5 and Ch. 7.
*Why it matters.* Chapter 3 contains no ligand-preparation section at all. §4.1.2 (`:427`) then reads EquiBind's form deficit as evidence that "direct regression trails diffusion and hybrid methods", while §3.2.2 (`:135`) has already told the reader that "a form-limited pose implicates conformer generation or refinement". The reader cannot separate the model from its harness.
*Evidence.* AutoDock and DiffDock receive one shared ETKDGv3 + **UFF** conformer; EquiBind receives thirty ETKDGv3 + **MMFF** conformers through "a custom pipeline rather than native EquiBind behaviour". Additionally, the AutoDock ligand path alone rebuilds hydrogens on 161 of 308 ligands (353 atoms, 352 on nitrogen), changing the donor set the Vina H-bond term reads.
*Correction.* State conformer source, generator, force field and pose budget per arm in the new §3.2 (M0); add a cross-reference from §4.1.2 and one clause to Ch. 7's comparability paragraph. *Effort: rewriting.*

**M5 — The Discussion introduces quantities that appear nowhere else.**
*Where.* §5.2, printed pp. 42–43, `body_main_short.tex:804`.
*Evidence.* "179 against 177" (form term dropped), "146 against 145" (validity term dropped) and EquiBind + smina's top-15 near-native count of 64 have no printed source in Ch. 4 or in the 1,454-line appendix; the attached claim that the EquiBind refiner difference "survives Holm correction" is stated without a test, a p-value or a family. (Two of the reviewer-flagged counts — 186/186 and 236/235 — *are* traceable to `:201` and Table 1.)
*Why it matters.* The refiner-selection argument the Discussion builds cannot be audited, and it breaches the Methods/Results/Discussion separation.
*Correction.* Move the three counts and the test into §4.1 or Appendix G with their families. *Effort: rewriting (the analysis appears to exist).*

**M6 — The metal-stripping stratification never reaches the interpretation chapters.**
*Where.* Disclosed at `:88` (Methods), `:891` (Ch. 7) and App. C.1 (`:165`). The word "metal" occurs **nowhere** in Ch. 5 or Ch. 6, and §4.1.2 (`:405`) states the pooled +13.2 / +11.2 / +10.2-point contrasts unqualified.
*Why it matters.* On the 225 metal-free complexes — 74 % of the analysed set — the top-15 AutoDock-over-DiffDock contrast falls to **+8.4 pp at an uncorrected exact p = 0.056**, against +19.2 pp (p = 0.028) on the 78 complexes from which the study's own preparation deleted a coordinating metal. The headline RQ1 ordering is therefore partly a preparation artefact.
*Correction.* One cross-reference sentence in §4.1.2 and one qualifying sentence each in §5.3 and §6.1, retaining the appendix's own caveat that the stratum difference is itself unresolved (Fisher p = 0.440). *Effort: rewriting.*

**M7 — The magnitude of the nearest-copy convention is invisible outside the appendix.**
*Where.* §4.1.2, §6.1 and both abstracts quote 49.2 / 43.6 / 18.8 %; Table 23 (App. H.1, p. 125) supplies the single-instance twins **36.6 / 34.0 / 18.2 %** (111 / 103 / 55 of 303) and shows the rank-1 lead halving from +5.6 to +2.6 pp.
*Why it matters.* About a quarter of the leading arm's rank-1 recoveries depend on a rule that takes a minimum over deposited copies. The convention is declared, justified by appeal to the source benchmark, and honestly tested — but a reader of Ch. 4 or the abstract cannot see its size. Compounding this, the crediting copies are never classified for binding-site equivalence; the only pocket-level statement in the document is negative ("Four alternates do not share the reference pocket", App. F).
*Correction.* Print the single-instance twins in parentheses at `:405` and in the abstract; classify the ~38 AutoDock and ~29 DiffDock rank-1 complexes that the rule supplies by whether the crediting copy shares the reference pocket. *Effort: rewriting + cheap re-analysis.*

**M8 — The primary DiffDock arm is a single unseeded draw and no arm was replicated.**
*Where.* Ch. 7 (`:887-888`), disclosed.
*Why it matters.* Variant selection turned on a **one-complex** margin for DiffDock (smina 144 vs gnina 143), and the headline rank-1 contrast has an interval (−1.7 to +12.8) whose width is comparable to plausible run-to-run variance that is nowhere estimated. Neither the author nor a third party can regenerate the benchmark DiffDock numbers.
*Correction.* Run 3–5 seeds on the three selected arms and report the between-seed spread of the rank-1 and top-15 endpoints; or state explicitly that the DiffDock refiner choice is arbitrary within noise. *Effort: additional computation (modest — DiffDock inference is the cheap arm).*

**M9 — The Orai1 placement endpoint and its causal interpretation both overreach.**
*(a) Slab scope.* Methods (`:95`) and App. D (`:583`) name the "mobile non-pore-lining transmembrane interfaces" — the lipid-facing M3 interface — as a proposed small-molecule modulator site and use that as a design rationale for four MD frames. The exclusion slab is **radially unbounded** and therefore deletes that band along with the conduction pathway. The thesis does disclose this in the body (`:636`) and defends it in App. C.6 as an artefact filter for a bilayer-free receptor, but the two statements are never reconciled, and the excluded poses are never split into pore-lumen and lipid-facing even though the radial coordinate is already computed for the placement envelope.
*(b) Caveat names the wrong compounds.* Ch. 7 (`:891`) says the slab "may exclude valid GSK-7975A or Synta-66 poses", whereas App. E (`:652`) states that the evidence for those two "converges on the extracellular pore mouth and the selectivity filter, which is what the placement criterion of the Results chapter encodes", and that it is **2-APB** for which the criterion "is a uniform study assumption rather than an established fact".
*(c) Causal attribution.* §6.1 (`:859`) and the abstract attribute the yield reversal to gnina re-ranking "rather than of the search", with the comparator clause "whereas the DiffDock cap is close to neutral on this endpoint". App. C.8 states that for DiffDock and EquiBind "Neither tool admits a prefix counterfactual… so the size of the bias cannot be measured from these runs."
*Correction.* (a) split the excluded poses by radial band — the data exist; (b) swap the compound names in the Ch. 7 caveat; (c) restate as: "…so the deficit is at least in part an effect of the gnina re-ranking applied before the cap. No equivalent prefix counterfactual can be computed for DiffDock, so the contribution of its own confidence cap is unknown." *Effort: rewriting + cheap re-analysis.*

**M10 — The Literature Review does not underwrite the study's own design.**
*Where.* Ch. 2, printed pp. 2–4.
*Evidence.* The 19 citation keys in Ch. 2 are all tool or scoring-function papers; **`ref035` (PoseBusters) is not among them**, and first appears in Methods §3.1.1. Blind whole-protein docking — the defining design and the explicit subject of RQ4 — is never reviewed; its entire literature basis (the CASF-2016 whole-protein penalty, the box-dilution citations) sits in App. C.1. No physical foundation for molecular recognition is given: no enthalpy/entropy decomposition, no desolvation, no non-covalent interaction taxonomy, no induced fit or conformational selection — although §3.2.2 and §4.1.4 reason with hydrogen bonds, hydrophobic burial, π-stacking and salt bridges. Orai1 receives no coverage in Ch. 2 at all. Table 7 catalogues four prior comparisons by tools, datasets and metrics but never states what any of them concluded.
*Correction.* Promote the App. C.1 blind-docking material into a §2.2; add ~1.5 pages on the physical basis of binding and a short Orai1 subsection; add a conclusions column to Table 7. *Effort: rewriting.* (Note: the Kurzfassung page has no slack, so German additions would need matching cuts.)

**M11 — The Motivation's third claimed gap is falsifiable from inside the document.**
*Where.* `body_main_short.tex:8`: "the absence of post-docking optimisation together with explicit physical-validity screening \cite{ref012}, \cite{ref035}".
*Evidence.* Table 7 lists "Physical Validity" among `ref012`'s metrics, and §5.6 (`:843`) states that `ref012` "reports that relaxation removes many AI steric failures". Two of the three claimed limitations also cite works (`ref011`, `ref035`) that Table 7 does not contain.
*Correction.* Restate the gap as what the evidence supports — that no prior comparison applies matched blind whole-protein search conditions *and* post-docking optimisation *and* validity screening *and* ranking-depth resolution to these three tools. *Effort: rewriting.*

**M12 — Neither the procedure's inputs nor the conclusions are independently reproducible.**
*Where.* App. A, printed p. 65.
*Evidence.* Repository "available from the author on request"; several hundred thousand output files also on request; no licence, no persistent identifier, no commit hash pinning the analysis behind Tables 1–6. The regeneration guide that would map floats to commands sits **inside** that private repository and admits an unnamed set of "tables whose provenance it did not audit". Among the items it "could not reconstruct" are the Fr0 receptor copy and the prepared ligand files — the former being the direct measurement substrate for the quantitative block of App. D. The four Orai1 receptor frames have no named source and no stated availability.
*Correction.* Archive code and the derived result tables (not the raw pose files) under a DOI with a licence and a commit hash; move the regeneration guide into the thesis appendix or the public archive; name the tables with unaudited provenance. *Effort: rewriting + archiving; no new computation.*

### Minor issues

1. **Statistical labelling.** §4.1.2 (`:405`) reports the top-5/15/30 contrasts as `p_holm = 7.6e-4 / 0.004 / 0.008`; App. H.9 (`:1299`) reports the identical values as "exact McNemar p = 7.6e-4, 0.0036 and 0.0075". One label is wrong. The conclusion survives either way (Holm over a family of three leaves all three below 0.05), but the correction status of the primary inferential claim should not be ambiguous.
2. **Figure captions.** Eleven of sixteen main-body figure captions are one-line titles with no panel key, axis definition, units or cohort size (e.g. Fig. 12 "Cross-tool agreement on Orai", Fig. 13 "Interaction-type profile on Orai"). Figures should be readable without the surrounding prose. Figs. 1, 6, 11, 15 and 16 show the standard the rest should meet.
3. **Wrong figure cross-reference.** §4.2.1 sends the reader to the receptor render (Fig. 27) for the transmembrane exclusion slab; the slab figure is Fig. 22.
4. **Compound naming.** `2-APB` is declared and used in prose, but every ligand table prints `2abp-NH2`; the two appendix tables also disagree with each other.
5. **Residue nomenclature.** One-letter (R91, V102, E106) and three-letter (Arg91, Glu106) codes alternate for the same residues across Methods, Results, Discussion and appendix.
6. **Spelling.** `body_main_short.tex:201` breaks the document's British convention twice in one sentence pair ("minimization", "minimizing"); the body otherwise uses `-is-` forms throughout.
7. **Abstract scope wording.** Both abstracts describe "a calibration benchmark of 303 co-crystallised complexes"; the benchmark contains 308 and 303 is the analysed cohort.
8. **Abstract traceability.** The 75 % Vina-prefix figure carried by the abstract and §6.1 never appears in Ch. 4, which quotes 70.0 % on a different inferential unit.
9. **Orai control-panel basis mixing.** §5.7 quotes DiffDock's control usable yield as 71.1 % (pooled, pose-level) beside two per-unit means; §4.2.3 gives 70.9 % for the same quantity.
10. **`217 against 214`** in §4.1.2 is not readable from Table 1, whose near-nativeness column prints 217 for both AutoDock arms.
11. **Fr0 description.** Called "the equilibrated starting structure" in §4.2.1 and "one pre-dynamics starting geometry" in §3.1.2, while App. D.1 shows it lies 5.88–6.02 Å outside the ensemble the other three frames sample.
12. **Appendix C.1 internal contradiction.** "the CNN scores and re-orders the poses without driving their geometry" contradicts App. C.5, Table 9's "Refiner mode: Energy minimisation" and three main-body statements that the AutoDock gnina pass ran `--minimize`.
13. **Undocumented Form guard.** Table 1's Form block applies an exploded-pose filter (`rmsd < 1000`) that removes 42 DiffDock coordinate-explosion poses; the Methods define Form as placement-independent, and the guard is not disclosed in the printed text.
14. **Model checkpoints.** DiffDock-L and EquiBind weights are identified by model family and date, while every other component carries an exact version.
15. **Discarded CHARMM protonation states.** The supplied Orai1 MD coordinates carried HSD/HSE/HSP labels — physically assigned per-residue tautomer states — which the preparation discards in favour of geometric re-protonation. The loss is documented but not remarked on.
16. **Ionisation caveat.** The deliberate absence of any pKa, tautomer or protomer treatment appears only in App. C and never as a limitation in Ch. 7.
17. **Macrocycles / covalent ligands / ring flexibility** are never mentioned, although two of three arms cannot alter a ring conformation they are handed.
18. **Vina pose redundancy.** The filter deciding which Vina modes are written is not reported, so the AutoDock top-15 pool may be deduplicated where the DiffDock and EquiBind pools are not.
19. **Resampling seed** described as "a fixed seed" without its value or generator.
20. **Benchmark file source.** Methods attribute the 308 complexes to the PDB; the appendix shows the docked files come from the PoseBusters distribution.
21. **Unqualified leakage claim.** §3.1.1 (`:88`) states without qualification that "The post-2021 release criterion excludes published training data of AI docking models" — a claim the thesis itself retracts three times (§5.4 `:822`, Ch. 7 `:891`, App. F.1). It should be stated as a *reduction* of temporal overlap from the outset. Relatedly, the training corpus and cut-off are named for only one learned component (gnina, CrossDocked2020); DiffDock-L's and EquiBind's are not.
22. **Boron handling documented for one arm only.** Boron typing is documented at parameter level for AutoDock (radius 1.92 Å, well depth), but neither DiffDock's featuriser coverage for boron nor the force field that produced 2-APB's conformers is reported for the learned arms — although 2-APB is one of the three experimental ligands.

---

## 10. Are the conclusions supported?

> **Status note added 2026-09-09.** This section records the assessment *as first filed*. All three overstatements below have since been corrected in the thesis — see Addendum 5. The section is left in its original form so the finding and the fix can be read against each other.

**Broadly supported, with three identifiable overstatements.**

*Fully supported.* RQ1's central finding — that generation rather than ranking dominates failure, that AutoDock and DiffDock are indistinguishable at rank-1 but separate from top-5 onwards, and that EquiBind is decisively behind both — is supported by paired tests, intervals, two reference conventions, ten thresholds and four depths. RQ2 — that local optimisation repairs geometry but does not supply a binding mode — is supported by the fixed-rank decomposition (0.3–1.0 pp for DiffDock against a 62.1-pp validity gain) and is one of the cleanest results in the thesis. RQ4's synthesis is proportionate.

*Overstated in three places.* (i) `:429`, the DiffDock refined-score claim, which is contradicted in sign by the thesis's own two measurements (M1). (ii) `:833`, the Orai1 filter/TM1–TM3 "recovery" sentence, contradicted by §4.2.5 and undercut by the study's own control panel (M3). (iii) `:859` and the abstract, the causal attribution of the Orai1 reversal to gnina re-ranking, whose comparator clause rests on a measurement the appendix says cannot be made (M9c).

*Understated in one place.* §5.7 says only that "Their rank-1 near-nativeness gap is unresolved", declining to state the practical conclusion that its own Holm-resolved top-5 and top-15 contrasts support — at the very depth its selection footnote (`:199`) argues is the more informative one.

*One conclusion carries an unstated qualification.* RQ1's ordering claim is materially weaker on the 74 % of the benchmark from which no metal was stripped (M6).

**The most important overclaim** is M3 (`:833`), with a paste-ready replacement given above.

---

## 11. Prioritised revision recommendations

| # | Recommendation | Addresses | Effort |
|---:|---|---|---|
| 1 | **Add a §3.2 "Docking protocol"** to Chapter 3, summarising per arm the receptor prep, ligand prep, search space, engine settings, refiner mode and ranking key, pointing to Appendix C Tables 9 and 11; move the §4.1 selection methodology into it | M0, and the root cause of M4/M6/M7 | Rewriting (1 day) |
| 2 | Delete or invert the refined-score sentence at `:429`; state the measured −2.0-point result | M1 | Rewriting (30 min) |
| 3 | Replace the §5.5 filter/TM1–TM3 sentence with the control-anchored version supplied in M3 | M3, RF15 | Rewriting (30 min) |
| 4 | Correct the false premise at `:871` and report campaign totals (17.41 Vina search hours) beside the 4.73 charged hours | M2 | Rewriting (1 h) |
| 5 | In the new §3.2, state the conformer source, generator, force field and pose budget **per arm**, and add a matching clause to Ch. 7's comparability paragraph | M4 | Rewriting (1 h) |
| 6 | Print the single-instance twins (36.6 / 34.0 / 18.2 %) in parentheses at `:405` and in both abstracts; classify the copies that supply the nearest-copy gain by pocket identity | M7 | Rewriting + ½-day re-analysis |
| 7 | Carry the metal stratification into §4.1.2, §5.3 and §6.1, with the Fisher p = 0.440 caveat | M6 | Rewriting (1 h) |
| 8 | Split the Orai1 slab-excluded poses into pore-lumen and lipid-facing bands (the radial coordinate is already computed); swap the compound names in the Ch. 7 caveat; soften the causal clause at `:859` | M9 | ½-day re-analysis + rewriting |
| 9 | Locate or remove the three unsourced counts in §5.2, with their test and family | M5 | Rewriting (1 h) |
| 10 | Run 3–5 seeds on the three selected arms and report the between-seed spread of rank-1 and top-15 recovery; failing that, state that the DiffDock refiner choice is arbitrary within noise | M8 | Additional computation (1–2 days) |
| 11 | Promote blind whole-protein docking and the PoseBusters framework into Ch. 2; add the physical basis of binding and a short Orai1 subsection; add a conclusions column to Table 7; archive code + derived tables under a DOI and licence | M10, M11, M12 | Rewriting (2–3 days) + archiving |

Eleven recommendations are listed because M0 was added after the initial pass; recommendations 2, 3, 4, 7 and 9 are same-day text fixes that close the hardest defects, and recommendation 1 is the structural change that prevents them recurring. Nothing in this list requires new experimental work, and only recommendation 10 requires new computation.

---

## 12. Five technically demanding viva questions

1. **On the primary endpoint.** Your headline rank-1 recovery scores each pose against the nearest of all deposited copies of the ligand. Table 23 shows the rule supplies 38 of AutoDock's 149 rank-1 recoveries and doubles the AutoDock–DiffDock lead from +2.6 to +5.6 points. For those 38 complexes, is the crediting copy in a pocket equivalent to the reference — and if some are crystal-packing or surface sites, in what sense has the tool "recovered" the binding mode? Given App. F reports four alternates that do not share the reference pocket, why is the single-instance convention the sensitivity arm rather than the primary?

2. **On the histidine model.** ADFRsuite protonated all 4,485 intact histidines on both ring nitrogens, making every one cationic. In Vina's PDBQT typing that removes the imidazole acceptor and leaves two donors, so the hydrogen-bond term sees a different receptor from the physical one at pH 7. You argue the model is uniform across the cohort. Uniform bias is still bias in a *between-tool* comparison — EquiBind never sees receptor hydrogens at all. What would you expect a HIS-tautomer-resolved re-run to do to the AutoDock and DiffDock validity and recovery figures, and could you bound it from data you already hold?

3. **On the re-ranker.** AutoDock's rank-1 recovery before gnina rescoring is 42.9 %, below DiffDock's 43.6 %; after it, 49.2 %. The +6.3 points the learned re-ranker supplies is larger than the +5.6-point lead it wins. DiffDock's refined scores were computed and written to the optimiser log but never used to order poses (App. C.3). Why is the asymmetry defensible, and what would the fair counterfactual be — a CNN-re-ranked DiffDock arm, or an AutoDock arm at raw Vina order?

4. **On the Orai1 slab.** Your slab is bounded axially between the Arg91 and Glu106 Cα planes and unbounded radially, so it deletes the lipid-facing M3 interface that §3.1.2 and App. D name as a proposed modulator site. Appendix E further says the outer-pore premise is best supported for GSK-7975A and Synta-66 and weakest for 2-APB, while Ch. 7's caveat names the first two as the compounds at risk. Take the poses you excluded: what fraction lie inside the ~7 Å pore lumen and what fraction against the outer shell, and does the DiffDock-over-AutoDock usable-yield reversal survive a lumen-only rule?

5. **On cost and selection.** RQ3 charges each pipeline for the configuration carried forward, defended at `:871` on the ground that AutoDock produced a single configuration — but Table 8 records five search passes totalling 17.41 Vina hours. Under a campaign-cost accounting that charges each tool for everything it ran to reach its selected arm, does the ordering EquiBind < DiffDock < AutoDock survive, and what does that say about the practical recommendation in §5.7?

## 13. Five broader conceptual viva questions

1. Your primary endpoint ANDs geometric accuracy with physical validity. Appendix H.4 shows the validity term removes at most five of 303 complexes, and App. H.4 also shows the discriminating intermolecular checks inspect coordinates Vina has already minimised. If the validity term is near-inert at the endpoint and partly constitutive of the tool that wins it, what work is it actually doing — and would a validity-*and*-strain-*and*-desolvation criterion be a better admissibility filter?

2. A pose that is PoseBusters-valid and within 2 Å of a crystal ligand is neither a binding free energy nor a residence time. Where, in a real structure-based design campaign, does the endpoint you optimised stop being the right thing to optimise — and what would you measure instead?

3. Your calibration set is cognate self-docking; your application is non-cognate cross-docking into a homology model. You treat the calibration figures as an upper bound rather than a prediction, but never measure the gap. What experiment inside the resources of this thesis would have quantified it, and would you now consider that a better use of compute than the exhaustiveness ladder?

4. The post-2021 release criterion is offered as a leakage control. `ref065`'s TEMPL result — 22.1 % overall rising to 67.3 % where related templates exist — suggests time is a weak proxy for novelty. If you re-ran this study with protein-, pocket- and ligand-novel splits, which of your four conclusions would you expect to survive, and which is most exposed?

5. All three pipelines hold the receptor rigid, and you represent flexibility with four MD snapshots that App. D shows keep the pore closed throughout and vary only in the extracellular loop. Is docking into a static ensemble a scientifically different claim from docking into one structure, or only a more expensive one — and what would it take to make the four-frame design an inference rather than a sensitivity statement?

---

## 14. Final weighted score

**80.0 / 100** (weighted mean 4.00 / 5; see §5 for the per-category calculation and the weight check).

## 15. Overall classification

**Very strong** (80–89 band). At the lower edge of that band: the appendix work is 85–90 quality and the main body's disclosure asymmetry, the Literature Review's gaps and a small number of hard textual defects hold the whole down.

## 16. Confidence in the final assessment

**High.** Fourteen independent reviews were run across the rubric's dimensions and every finding was adversarially re-tested against the verbatim source; 5 of 91 findings were refuted outright and a further 59 were downgraded or narrowed on verification, which is the expected profile for a document that discloses heavily in its appendix. The three highest-severity items (M1, M2, M3) and the two structural ones (M5, M9b) were additionally verified by hand against the LaTeX source and the cited tables. Confidence is lower for the Orai1 chapter only because the receptor's provenance gap places a ceiling on what *any* reader can verify — a limit the thesis itself states rather than conceals.

## 17. Examiner recommendation

This is a strong MSc thesis that does the difficult thing well. The candidate designed a comparison that most published work in this area does not attempt — matched blind whole-protein conditions, a compound accuracy-and-validity endpoint, ranking depth treated as a variable rather than a footnote, and a cost basis defined once and applied identically — and then spent an unusual amount of effort trying to break their own result rather than defend it. The appendix is the intellectual centre of the work and is, in places, better than the published literature it calibrates against: Table 23, the exhaustiveness ladder, the metal stratification and the Orai1 prefix counterfactual are all analyses that a careful referee would have had to ask for, and each is already there, reported with the result rather than merely with its existence. The defects that remain are almost entirely defects of *transmission*, not of *execution*: qualifications stranded in the appendix while the body states the unqualified result, a Literature Review that does not underwrite the design the rest of the thesis defends, and three sentences that say more than the data behind them allow — one of which (`:429`) says the opposite. All twelve major issues are correctable by rewriting or by re-analysis of data already held; only one (M8, seed replication) needs new computation, and that is a day or two on the cheapest arm. I recommend acceptance subject to minor-to-moderate revision, with recommendations 1–4, 6 and 8 treated as required and the remainder as strongly advised.

---

### Does this thesis demonstrate that the student can independently plan, execute, critically evaluate, and communicate a scientifically sound protein–ligand docking project at master's level?

**Yes, with minor reservations.**

Planning and execution are demonstrated beyond doubt: three pipelines, eight AutoDock variants, nine EquiBind configurations, two Orai1 panels and 72,219 screened poses, all under a protocol documented to a standard that permits re-implementation. Critical evaluation is the thesis's outstanding quality — the candidate quantifies their own selection optimism, prints the ranking head that would have scored better, and repeatedly declines to resolve contrasts their data cannot resolve. The reservations are confined to communication and are specific rather than diffuse: material qualifications sit in the appendix while the body states the unqualified result, the Literature Review does not support the study's own design, and three sentences overstate — one of them against the thesis's own measurement. These are revisions to a well-executed study, not evidence of a weakness in the science or in the candidate's judgement.

---

## Addendum — disposition of M0–M6 (2026-09-09, after author response)

The findings above are left as filed. This section records what was ruled on, what was applied and what was re-validated.

| ID | Disposition |
|---|---|
| **M0** — no docking-protocol section in Ch. 3 | **Not a defect.** The author confirms the protocol was moved to Appendix C deliberately. Withdrawn. The consequential parts of M4, M6 and M7 stand on their own and do not depend on it. |
| **M1** — `:429` inverted in sign | **Applied.** Sentence deleted. No other passage in the body, appendix or abstracts referred to it. |
| **M2** — `:871` false premise | **Applied.** The two false sentences are replaced by one that states the accounting rule symmetrically: "Every tool is charged for the configuration carried forward and for nothing else, so the selection work each pipeline discarded on the way to that configuration is excluded on all three arms alike." Only dominant-variant timings are now referenced. |
| **M3** — `:833` overclaim | **Applied**, in the thesis's register rather than the draft wording. |
| **M4** — ligand non-parity | **Closed, no action needed.** Already the first paragraph of Appendix C (`body_appendix_short.tex:148`), repeated at `:258` and tabulated in the "Start conformer" row of Tables 9 and 11. With M0 withdrawn, the information is where the author intends protocol detail to live. |
| **M5** — unsourced Discussion counts | **Re-validated, holds (narrowed). Rewrite applied.** See below. |
| **M6** — metal stratification absent from interpretation | **Re-validated, holds. Both insertions applied.** See below. |
| **M11** — Motivation's third claimed gap | **Applied.** See below. |

Build after all edits (M1, M2, M3, M5, M6 and the §4.2.5 seam): `latexmk -pdf` exit 0, no errors, no undefined references, converged, **146 pages** — identical to the committed PDF at HEAD. (Note for the record: plain `latexmk` without `-pdf` takes the DVI route on this project, mis-handles the PNG figures and inflates the document to 155 pages. Always build with `latexmk -pdf`.)

### M5 re-validation — holds, for four quantities

§5.2 (`body_main_short.tex:804`) contains seven paired counts. Three are traceable and four are not.

| Quantity | Status |
|---|---|
| 144 to 143 (top-15 selection depth) | Traceable — `:201` |
| 186 against 186 (pooled near-native) | Traceable — Table 1 |
| 236 against 235 (pooled form) | Traceable — Table 1 |
| nine-complex EquiBind margin | Traceable — `:199` (57 against 48) |
| **179 against 177** (form term dropped) | **Not printed anywhere.** The only other `179` in the appendix are DiffDock's 179 of 308 in the intention-to-treat arm and the "179 to 199" Form range of the guided EquiBind arms. Neither is this quantity. |
| **146 against 145** (validity term dropped) | **Not printed anywhere.** |
| **64 for smina** | **Not printed anywhere.** Table 1 and App. C.4.1 give EquiBind + smina 67 near-native complexes over the full pool. |
| **86 "at top-15"** | Printed, but as the **full-pool** count (Table 1; App. C.4.1 "19, 67 and 86 complexes unguided"), not at top-15. |
| "that difference survives Holm correction" | **No test, p-value or family given** anywhere for this contrast. |

**Rewrite of §5.2, second half — APPLIED** (uses only printed quantities and drops the two significance claims that have no printed test). Table 19 was checked first and tests raw-to-optimised recovery within each tool, not smina against gnina, so no printed test backs the deleted clause:

> Refiner selection remains an empirical pipeline choice, and the two families behave differently. For DiffDock the refiners are interchangeable. smina is carried forward on the selection criterion by 144 complexes to 143 at the top-15 selection depth, while Table~\ref{tab:results-pose-production} leaves the two level on pooled near-native recovery at 186 against 186, puts smina ahead on pooled form at 236 against 235 and gnina ahead on pooled validity at 87.1\% against 84.5\%. Exact McNemar separates them at no depth. The one-complex margin therefore rests on which components the endpoint combines rather than on a measurable difference between the refiners. For unguided EquiBind the choice is not immaterial. Over the whole pool gnina recovers 86 near-native complexes against 67 for smina and reaches 62.0\% pooled validity against 51.3\%, and at the top-15 selection depth its margin is nine complexes, at 57 against 48. The study measured neither optimisation trajectories nor refiner performance on an independent set, so proposed mechanisms such as gentle preservation or aggressive repair remain interpretations.

What this changes: it drops "179 against 177" and "146 against 145", restates the EquiBind contrast on the full-pool counts Table 1 actually prints (86 against **67**, not 64), and removes the unsupported "survives Holm correction" clause. If the 179/177 and 146/145 figures do exist in the analysis outputs, the alternative is to keep the original sentence and add the two endpoint-variant columns to Appendix G, where the reader can audit them.

### M6 re-validation — holds

`metal` occurs **zero times** in Chapter 5 (`:786-851`) and **zero times** in Chapter 6 (`:852-882`). Its single occurrence in Chapter 4 (`:694`) is the interaction class "metal coordination" on Orai1, not the stratification. The stratification is therefore confined to Methods `:88`, Ch. 7 `:891` and Appendix C.1 `:165`, while §4.1.2 `:405` and §6.1 `:857` state the pooled contrasts unqualified. Both now carry the qualification.

**Insertion 1 — §5.3 — APPLIED**, after "…while AutoDock led on near-nativeness from top-5 onwards.":

> That lead is not uniform across the benchmark. On the 225 complexes from which no metal was removed the top-15 contrast narrows to 8.4 percentage points at an uncorrected exact p of 0.056, against 19.2 points on the 78 metal-adjacent complexes, although the difference between the two strata is itself unresolved (Appendix~\ref{autodock-vina-2}).

**Insertion 2 — §6.1 — APPLIED**, after "That ordering is conditional on search effort.":

> It is also conditional on receptor preparation, since the metal-adjacent quarter of the benchmark carries a wider AutoDock margin than the metal-free remainder that supplies three quarters of the set (Appendix~\ref{autodock-vina-2}).

Both use the `autodock-vina-2` label the body already cites at `:88` and `:891`, and both carry the appendix's own caveat that the stratum difference is unresolved (Fisher p = 0.440), so neither overstates the heterogeneity.

### Consequence of the M3 edit — APPLIED

§4.2.5 `:718` opened "The surviving modulator poses partially agree with published pharmacology" and closes "These contacts support biological interpretation and hypotheses for mutagenesis or ligand-bound simulation". The Discussion now says the same contacts are not discriminative against the control panel. The Results sentence is not wrong — it is about the poses matching the literature description — but a reader meeting `:718` before `:833` will feel the change of footing. The closing sentence of `:718` now carries the control comparison:

> These contacts are not specific to the modulators, because the control ligands engage the same residues at rates that do not differ after correction. They therefore support hypotheses for mutagenesis or ligand-bound simulation rather than a binding-mode assignment, and evidence from different tools cannot be combined to validate one binding mode.

Results and Discussion are now on the same footing on this point.


---

## Addendum 2 — M11 applied (2026-09-09)

**Verified against Table 7 before rewriting.** Three of the four tabulated studies report a validity metric (`ref012` "Physical Validity", `ref013` "PB-valid", `ref032` "PB-valid" plus a **combined RMSD-and-PB-valid criterion** — the very endpoint this thesis uses), and §5.6 records that `ref012` relaxes poses. The claim that this literature lacks "post-docking optimisation together with explicit physical-validity screening" is therefore not defensible from the table it points at. Limitations one and two do hold and are unchanged.

**Restated at `body_main_short.tex:8`**, from an absence claim to a conjunction claim:

> …carry recurring limitations, namely evaluation from a selected output pose [12], [13], [33] and search conditions that are not matched across the competing methods [11], [33]. **Validity screening is established in that literature [12], [13], [31], but no tabulated study combines it with matched blind whole-protein search and recovery resolved across ranking depth.** None of the tabulated evaluations addresses Orai1.

This is now consistent with §5.6 rather than contradicting it, and it concedes the point the table actually makes while stating the gap the design genuinely fills.

**Pagination note.** The Motivation page has *zero slack*, exactly like the Kurzfassung — at HEAD the RQ1–RQ4 list ended on the final line of printed page 1. A first, fuller restatement pushed the document to 147 pages and orphaned RQ4 onto page 2. To hold 146 pages and keep the research questions with the paragraph that motivates them, one adjacent sentence in the same paragraph was trimmed:

> …while the fourth benchmarks classical and ligand-guided cross-docking ~~strategies~~ without an AI arm and ~~so~~ enters the table as a non-AI reference ~~rather than as a physics-versus-AI comparison~~ [32].

No meaning is lost — the excised clause restated what "non-AI reference" already says. Revert it if you would rather carry the extra page. **Treat Motivation page 1 as page-budget-locked in future edits.**

Build: `latexmk -pdf` exit 0, no errors, 146 pages, page 1 ends exactly as at HEAD.

### Open after this round

Majors: **M7** (single-instance twins absent from body and abstracts, copy site-equivalence unestablished), **M8** (unseeded DiffDock, no seed replication), **M9** (slab lumen-vs-lipid split; Ch. 7 caveat omits 2-APB; `:859` unmeasurable comparator), **M10** (Literature Review), **M12** (archiving). All 22 minors.


---

## Addendum 3 — minor issues round (2026-09-09)

### Minor 1 (statistical labelling) — **REFUTED, my error. No edit made.**

Both labels are correct and the values agree because Holm leaves the largest p in its family unadjusted. App. H.3 declares the family as "the three tool-versus-tool contrasts within one metric", so at each depth the family is {AutoDock–DiffDock, AutoDock–EquiBind, DiffDock–EquiBind}. §4.1.2 records the two EquiBind contrasts at `p_holm < 3e-7`, so the AutoDock–DiffDock contrast is the **largest** of the three. Holm's step-down multiplies the largest by 1, and the two tiny contrasts stay far below the threshold after their ×3 and ×2 factors, so monotonicity never binds. Therefore `p_holm` = raw exact McNemar p exactly, and §4.1.2's `7.6e-4 / 0.004 / 0.008` and App. H.9's `7.6e-4 / 0.0036 / 0.0075` are the same numbers under two correct labels, differing only in rounding. **Withdraw this finding. Do not re-file it.**

### Applied

| Item | Change |
|---|---|
| Minor 3 — wrong figure cross-reference | §4.2.1 now reads "Figure 27 illustrates the Fr300 snapshot, and Figure 22 in Appendix C.8 adds the transmembrane exclusion slab." |
| Minor 4 — compound naming | All **7 printed** occurrences of `2abp-NH2` replaced by `2-APB` (`:97` declaration, Table 5 ×3, App. E, App. F, Table 17). The `:97` declaration now reads "These are 2-APB, docked in the neutral protonation state its appendix entry records". Verified first that **no figure image carries a ligand label** — Figs. 9 and 13 facet by tool and interaction type — so no raster/text mismatch is created. The three remaining `2abp` strings are LaTeX provenance comments recording the source filename `2abp-nh2-OPT.sdf` and were deliberately left. |
| Minor 6 — spelling | `:201` "minimization"→"minimisation", "minimizing"→"minimising". The two surviving `minimiz` strings are the literal `\texttt{-{}-minimize}` flag and must keep US spelling. |
| Minor 9 — Orai basis mixing | §5.7 `71.1\%`→`70.9\%`, matching the per-unit mean basis of the 82.8 % and 3.8 % it stands beside (§4.2.3). |
| Minor 10 — 217 against 214 | §4.1.2 now names the gate and both tables: "…nominally ahead **on the validity-aware gate**, 217 complexes against 214 in Table 8, because rescoring minimises every pose and a few cross the validity threshold the wrong way. **The two arms are level at 217 on near-nativeness alone in Table 1, so the difference sits in the validity term rather than in placement.**" |
| Minor 12 — App. C.1 contradiction | Reworded to match App. C.5, Table 9 and the three body statements: "…so the empirical objective relaxes each pose into its nearest minimum and the CNN then scores and re-orders the pool **without itself driving the geometry**." |

Build: `latexmk -pdf` exit 0, no errors, no undefined references, **146 pages**.

### Open after this round

Majors: **M7**, **M8**, **M9**, **M10**, **M12**. Minors: 2, 5, 7, 8, 11, 13–22 (minor 1 withdrawn; 3, 4, 6, 9, 10, 12 closed).

---

## Addendum 4 — abstract and Kurzfassung audit (2026-09-09)

Both versions were re-verified claim by claim against the updated body, and the German was checked sentence by sentence against the English. **No number in either version was wrong.** Every figure reproduces exactly somewhere in the body, and the numeral sequences of the two versions are identical apart from English "top-5 / top-15" against German "fünf / fünfzehn". The defects were of *basis*, *qualifier* and *strength*, and six were fixed.

### Fixed — English

| Was | Now | Why |
|---|---|---|
| "the gaps at top-5 and top-15 are statistically resolved." | "…statistically resolved, though conditional on search effort and receptor preparation." | §5.3 and §6.1 now qualify the lead by the metal stratum (M6). The abstract carried the pooled verdict unconditionally over the 74 % of the benchmark where the top-15 contrast is +8.4 pp at p = 0.056. |
| "the validity ordering held" | "the validity ordering held **on the panel means**" | On the per-unit *median* basis that governs the 65/60/0 in the same sentence, Orai validity is a three-way tie at 100 % (§4.2.2). The ordering holds on means (100.0 / 93.3 / 89.2). |
| "That reversal is a property of…" | "That **untested** reversal is a property of…" | §6.1 calls it "This untested reversal … descriptive"; the abstract stated it flat. |

### Fixed — German

| Was | Now | Why |
|---|---|---|
| "war der **AutoDock-Ablauf** am teuersten" | "war der **konfigurierte** AutoDock-Ablauf am teuersten" | Most serious of the German issues. The body explicitly denies the unrestricted reading twice — "This does not apply to Vina alone" (`:873`), and Vina alone is 17.8 s, *below* DiffDock. The German as written asserted something the thesis contradicts. |
| "lag die **Wiederfindung gültiger Posen** … bei 49,2 %" | "lag der **Anteil der Komplexe** mit einer gültigen Pose … bei 49,2 %" | Basis slip present only in German: it read as a share of poses, where the endpoint is complex-level (149 of 303). |
| "möglicherweise zu optimistisch" | "**für unbekannte Zielproteine** möglicherweise zu optimistisch" | German made the optimism claim unbounded; the body bounds it to transfer beyond the benchmark and quantifies the within-benchmark winner's curse at 0.4–1.0 pp. |
| "Rangfolge der Gültigkeit erhalten" | "Rangfolge der Gültigkeit **im Mittel** erhalten" | Same basis fix as the English. |
| "Diese Umkehr stammt aus…" | "Diese **ungeprüfte** Umkehr stammt aus…" | Same as the English. |
| "gemeinsame **Beurteilung** von Genauigkeit und Gültigkeit" | "gemeinsames **Filtern nach** Genauigkeit und Gültigkeit" | "Beurteilung" is assessment; §5.7 recommends a filter that "rejects invalid or biologically inadmissible structures". |
| "statistisch gesichert." | "statistisch gesichert, **wenn auch abhängig von Suchaufwand und Rezeptorvorbereitung**." | Same as the English. |

### Considered and rejected

**"against 99.6 % for AutoDock" is the raw arm, not the configured 98.9 %.** Correct as written and *not* to be changed. The sentence is explicitly about raw output — "Raw poses from the learned workflows … at 24.4 % … and 2.9 % … against 99.6 % for AutoDock" — so it is a raw-vs-raw comparison, and that contrast is the point being made. Substituting 98.9 % would introduce an error.

### Still carried by both abstracts (open, matching the body)

The causal clause "a property of AutoDock's learned re-ranking rather than of its search" is **M9c** and remains open in `:859` as well. The abstracts were left aligned with the body rather than fixed unilaterally; fixing M9c should change all three together.

### Pagination — a live trap, now recorded

Both abstract pages are page-budget-locked, not just the German one. The additions above initially overflowed **both**, orphaning `Keywords` and `Schlagworte` onto pages of their own and taking the build from 146 to 148 pages. The binding constraint is **not** the visible text block: the frame runs to ~749 pt, but `\paragraph*` needs ~25 pt of headroom, so prose must end by about **700 pt**. The English abstract still overflowed at 710 pt with ~130 pt of apparent white space below it. Ten compensating trims (six English, four German) restored 146 pages with both keyword lines back in place.

Final state: `latexmk -pdf` exit 0, no errors, **146 pages**, English abstract p. 3 and Kurzfassung p. 4, each complete on one page.


---

## Addendum 5 — the three overstatements are closed (2026-09-09)

Section 10 filed three overstatements. All three are now corrected in the thesis.

**(i) `:429`, the DiffDock refined-score claim — DELETED.** The sentence claiming refined-score ordering was "worth between three and four points" is gone. Zero occurrences remain, and nothing elsewhere in the body, appendix or either abstract referred to it. (Applied in Addendum 1.)

**(ii) `:833`, the Orai1 filter/TM1–TM3 "recovery" sentence — REPLACED.** Now reads "The contacts recovered at the literature-associated filter and TM1/TM3 residues are **not discriminative** either…", with the control-panel equality stated in the same breath. §4.2.5 `:718` received a matching clause so Results and Discussion no longer sit on different footings. (Applied in Addendum 1 and Addendum 2.)

**(iii) M9c, the causal attribution of the Orai1 reversal — REPLACED in all three places.**

The defect was two clauses. "*whereas the DiffDock cap is close to neutral on this endpoint*" asserted something App. C.8 states cannot be measured — "a ten-sample run draws its own poses rather than truncating a thirty-sample list, so the size of the bias cannot be measured from these runs" — and the following "*therefore a property of*" leaned on it, since attributing a between-tool difference to one tool's ranker requires knowing the other's cap contributed nothing.

The fix substitutes a measured mechanism for the unmeasurable comparator. App. C.8 already establishes the gnina ranker's bias directly on 12,315 control poses (77.3 % of first-five ranks outside the slab against 90.1 % of ranks six to ten, Spearman −0.211 at p = 7e-124, where the Vina order moves 0.5 points at p = 0.63) — a chemotype-independent result the Conclusions were not using.

*§6.1 now reads:*
> AutoDock is capped on its gnina re-ranking, and restricting it instead to Vina modes one to ten puts it at 75 % and removes the reversal. **That ranker measurably promotes membrane-buried poses on the control panel as well. The deficit is therefore at least in part an effect of the secondary ranker rather than of the search, though DiffDock admits no equivalent counterfactual, so the contribution of its own cap is unmeasured.**

*English abstract:* "is a property of … rather than of its search" → "**tracks** … rather than its search" (identity claim to association, and two words shorter).

*Kurzfassung:* "**stammt aus**" (causal) → "**folgt**" (associative), matching the English and one word shorter.

Discussion §5.5 `:829` needed no change — it already says "That ordering depends on the ranking each pipeline applies before the cap" without asserting DiffDock neutrality or exclusive causation.

Build: `latexmk -pdf` exit 0, no errors, **146 pages**, English abstract p. 3 and Kurzfassung p. 4 each complete on one page.

### Section 10 restated as it now stands

Conclusions are **broadly supported**. No overstatement identified in this review remains in the thesis. What remains open is *reporting* rather than *overclaiming*: the single-instance twins of the headline endpoint are absent from the body and abstracts (M7), the primary DiffDock arm is one unseeded draw (M8), the Orai1 slab excludes without decomposition (M9a) and its Ch. 7 caveat omits 2-APB (M9b), the Literature Review does not underwrite the design (M10), and the artefacts are not archived (M12).

---

## Addendum 6 — every percentage now names its tool (2026-09-09)

Requested change: in both abstracts, no percentage should be left for the reader to attach to a tool by inference.

**Two places were unattributed**, in both languages, and both relied on carrying a tool order across a sentence boundary:

| Was | Now |
|---|---|
| "lifting the median per-complex validity yield from 13.3 % to 93.3 % and from 0 % to 80.0 %" | "lifting **DiffDock's** median per-complex validity yield from 13.3 % to 93.3 % and **EquiBind's** from 0 % to 80.0 %" |
| "raised recovery to 70.3 % and 59.1 %" | "raised recovery to 70.3 % **for AutoDock** and 59.1 % **for DiffDock**" |
| "restricting **it** to Vina modes one to ten puts it at 75 %" | "restricting **AutoDock** to Vina modes one to ten puts it at 75 %" |

German equivalents applied identically — "hob **DiffDocks** medianen Gültigkeitsanteil … und **EquiBinds** von 0 % auf 80,0 %", "Wiederfindung auf 70,3 % **für AutoDock** und 59,1 % **für DiffDock**".

Source check for the first pair, §4.1.1 (`body_main_short.tex:331`): "DiffDock rises from a 13.3 % raw median to 93.3 % **with smina** … Unguided EquiBind rises from 0 % to 80.0 % **with gnina**." Both are the carried-forward refiners, so the figures belong to the two workflows the abstract names in its second sentence.

**One residual, deliberately left:** "leaving the leading pair unresolved at 5.6 points (95 % interval −1.7 to +12.8)". The 5.6 is a *difference between* AutoDock and DiffDock, named in the immediately preceding clause, and the 95 % is a confidence level rather than a result. Naming tools again there would be redundant.

**A basis juxtaposition worth knowing about.** The sentence before gives *pooled* validity for the same raw arms (24.4 % DiffDock, 2.9 % EquiBind); this one gives *median per-complex* yield (13.3 %, 0 %). Same tools, same raw state, different numbers, because pooled ≠ median on a skewed per-complex distribution. Both are labelled, and both are correct, but the two readings sit two clauses apart.

**Selection-optimism clause.** "and possibly optimistic for unseen targets" is deliberately absent from the English. It was briefly restored in error during this round and has been removed again, and the German was mirrored ("ist der Vergleich deskriptiv.") so the two versions no longer diverge on it — the earlier German-only breadth was itself a finding in Addendum 4.

Fourteen compensating trims were needed to hold the page budget (nine English, five German), all wording-level with no loss of content. Verified after rebuild: result numerals identical between the two versions (the only difference being English "top-5 / top-15" against German "fünf / fünfzehn"), `latexmk -pdf` exit 0, no errors, **146 pages**, abstract p. 3 and Kurzfassung p. 4 each complete on one page.

---

## Addendum 7 — DiffDock confidence-model claim validated and extended (2026-09-09)

The §2.1.2 sentence "A confidence model ranks poses by whether they fall below a 2\angstrom{} RMSD threshold rather than by an affinity-like energy. Performance therefore depends heavily on the training data and on confidence calibration \cite{ref024}" was checked on three axes — internal consistency, primary literature, adversarial — each adversarially verified, consulting both papers and the installed DiffDock source tree.

**Verdict: substantially correct.** Decomposed:

| | Sub-claim | Verdict |
|---|---|---|
| C1 | Two learned components | Accurate |
| C2 | Score model steers denoising and does not rank | Accurate — verified in the paper *and* in the source tree the run loads (the score model's output enters `sampling()` and never the sort) |
| C3 | Confidence model ranks by a 2 Å RMSD threshold | Accurate but imprecise |
| C4 | Not an affinity-like energy | Accurate, and load-bearing across the thesis |
| C5 | "depends heavily on training data and confidence calibration \cite{ref024}" | Imprecise, wrong citation |

**The 2 Å number survives scrutiny**, which was the main risk. It is the literal positive-class label of the confidence model's binary cross-entropy objective, and it holds for DiffDock-L rather than only for the 2023 release. The one route by which it could have been wrong — the 2 Å / 4 Å two-sided labelling in ref091 — belongs to DiffDock-S under Confidence Bootstrapping, and under either scheme the positive class is RMSD < 2 Å.

**C3** described an oracle: knowing whether a pose falls below 2 Å requires the crystal pose, which is unavailable at inference. Appendix B.2 already had the right register ("scores each generated complex **for whether** it falls below"), as did the thesis's own gnina sentence at `:66`. Only `:48` slipped. **C5** attributed a training-data claim to ref024 while ref091 is both the model actually run and the paper whose headline result is generalisation; "confidence calibration" appeared exactly once in the thesis, undefined and unmeasured, and is not something ref024 measures.

Both corrected at `:48`, which now reads "…ranks poses by **its predicted probability that each falls below** a 2\angstrom{} RMSD threshold…" and "…on **how well that confidence ranking separates good poses from bad** \cite{ref024}, \cite{ref091}", with "Full details are given in the Appendix" softened to "provides further details" (the run configuration sits in a different chapter than the pointer names).

### A new finding, and it runs in the thesis's favour — now stated at §5.3

DiffDock's confidence model is trained on **exactly the criterion the primary endpoint measures**, RMSD < 2 Å. So is gnina's CNNscore, by the thesis's own description at `:66`. But the AutoDock arm is ranked by **CNNaffinity**, a pK regression that is *not* endpoint-aligned — and App. C.1 records that CNNscore would have recovered 152 complexes at rank-1 against CNNaffinity's 149. The comparison therefore hands DiffDock a ranker optimised for the scoring criterion and hands AutoDock one that is not, and AutoDock still leads from top-5 onwards. A grep confirmed no sentence anywhere noted this.

Added to §5.3 ¶1, bounded by the one measurement available so it states a quantified small effect rather than a speculative advantage:

> The two leading rankers are also not equally aligned with the endpoint. DiffDock orders its poses with a confidence model trained to separate poses within 2\angstrom{} of the reference from poses outside it, which is the criterion scored here. The AutoDock pool is ordered instead by gnina CNNaffinity, a regression onto a binding constant. The alignment favours DiffDock by construction. Where it can be measured it is worth little, since substituting gnina's own 2\angstrom{} discrimination head for CNNaffinity moves AutoDock rank-1 recovery from 149 complexes to 152 on an unresolved paired contrast (Appendix~\ref{autodock-vina-2}). AutoDock leads from top-5 onwards under the less favourable head.

Placed in §5.3 rather than Ch. 7 because it is an interpretation that reads as a strength, and §6.1 already carries two conditionals.

**Deliberately excluded, and worth closing separately:** EquiBind. The thesis says only "gnina affinity" for its ranking key — at `:200`, `:513`, App. C.5 `:274`, `:553`, `:964` and in the "Analysis ranking" row of both protocol tables — which never resolves to CNNaffinity or to the empirical minimised affinity, two different quantities. AutoDock's key is unambiguous ("gnina CNNaffinity", Table 9 and App. C.1), so the claim was restricted to the two tools whose ranker can be named exactly. One word in Table 9's EquiBind column would close it, and an examiner reading the new §5.3 sentence will ask.

Build: `latexmk -pdf` exit 0, no errors, no undefined references, **146 pages**, both abstracts still complete on one page each. Reference [80] resolves to the DiffDock-L paper.

---

## Addendum 8 — EquiBind's ranking key disambiguated (2026-09-09)

Adding the §5.3 ranker-alignment sentence exposed a latent ambiguity: the thesis described AutoDock's ranker precisely ("gnina CNNaffinity") but EquiBind's only as "gnina affinity", seven times, never resolving which quantity that is.

**It is not the same head, and it does not even point the same way.** Verified against the pipeline rather than inferred:

- **AutoDock** sorts `cnn_affinity` **descending** — CNNaffinity, the learned pK head that §2.1.3 describes.
- **EquiBind** sorts `gnina_affinity` **ascending, most-negative first**, which the code calls "the gnina **minimisation energy**" (`Scripts/Analysis/orai_transmembrane_exclusion.py:740`, `Scripts/Analysis/run_pandamap.py:768`).

Most-negative-first is an empirical energy in kcal/mol, not a pK. A reader taking "gnina affinity" to mean CNNaffinity — the natural reading after §2.1.3 — would conclude that EquiBind and AutoDock share a ranking criterion. They do not, and the difference is a sign convention as well as a quantity.

All seven occurrences corrected. Both protocol tables (Tables 9 and 11) now read "gnina minimised affinity, **most negative first**", carrying the direction; the five prose mentions read "the gnina minimised affinity" (`body_main_short.tex:198`, `:511`; `body_appendix_short.tex:274`, `:553`, `:964`). §4.1.4 `:511` was the most exposed, since it lists all three rankers side by side, and now reads "AutoDock uses gnina **CNNaffinity** re-ranking, EquiBind **the gnina minimised affinity** and DiffDock native confidence retained through smina refinement." No bare "gnina affinity" remains in the rendered document.

**This retroactively justifies the §5.3 restriction.** EquiBind was excluded from the ranker-alignment sentence because its key could not be named exactly. Now that it can, the exclusion is the right call on stronger grounds: EquiBind is not a third learned head at all but an empirical energy, which is the ordering criterion furthest of the three from the 2 Å endpoint the study scores.

Build: `latexmk -pdf` exit 0, no errors, no undefined references, **146 pages**, both abstracts complete on one page each.

**Process note for future edits:** the shell working directory reset to the repository root mid-session and one `latexmk` invocation consequently ran outside `thesis_latex/`, reporting three errors and no PDF. The committed PDF was untouched. Always confirm the working directory before reading a build result as a failure.

---

## Addendum 9 — where the two quantities are explained, and a correction to Addendum 8 (2026-09-09)

**The difference is explained, and it is even measured.** Three places:

1. **§2.1.3, `body_main_short.tex:66` (p. 4)** defines both. The empirical quantity: "gnina is a fork of smina, thus share the Vina docking algorithm and **Vina energy, i.e. quantity optimised during search**". The learned heads: CNNscore, "trained to distinguish poses within 2\angstrom{} RMSD of the native mode from poses outside it", and CNNaffinity, "**a regression head that predicts a binding constant on a pK scale**".
2. **§2.1.3, `:68`** gives the mechanism that separates them: "the inherited quasi-Newton minimiser first drives the pose to a local minimum of the empirical function. Afterwards the CNNs are evaluated and the poses reordered by the learned score."
3. **Appendix C.1, `body_appendix_short.tex:173` (p. 84)** quantifies all three orderings on one pose pool: CNNaffinity recovers 149 complexes at rank-1, CNNscore 152, and **the minimised Vina energy 128**. The two convolutional heads are unresolved against each other (exact McNemar p = 0.720) while both beat the empirical energy by seven to eight points.

**What was missing** was any link from EquiBind's ranker to that contrast — and Addendum 8 made it worse rather than better. "Gnina minimised affinity" was a **fourth name** for a quantity the thesis already calls "the minimised Vina energy" at `:173`, so the disambiguation named the quantity precisely while disconnecting it from the two passages that define it.

**Corrected.** All seven occurrences now read **"the minimised Vina energy"**, the thesis's own term, which is defined by `:66` and `:68` and already carries a measured comparison against CNNaffinity at `:173`. The protocol tables read "minimised Vina energy, most negative first", where the sign convention now follows from the word *energy* rather than needing to be taken on trust.

`body_main_short.tex:198`, where EquiBind's ranker is first named, now closes the loop explicitly:

> Since EquiBind provides no native confidence score, poses are ranked by the minimised Vina energy, **the empirical objective of the search rather than the learned head that reorders the AutoDock pool**.

A reader meeting that sentence is pointed straight at the §2.1.3 distinction, and Appendix C.1 already tells them what the choice is worth.

Build: `latexmk -pdf` exit 0, no errors, **146 pages**, both abstracts complete on one page each.

---

## Addendum 10 — multi-copy convention made explicit in Methods §3.2.1 (2026-09-09)

The near-nativeness definition named the rule for choosing among deposited ligand copies before telling the reader that multiple copies exist. "The 165 single-copy complexes are unaffected" arrived with no antecedent, so the multi-copy population had to be inferred by subtraction from 303.

Replaced at `body_main_short.tex:110`:

> **Many entries deposit the ligand more than once. This is true of 138 of the 303 analysed complexes, and over all multi-copy entries of the source set the nearest other copy lies a median 36.0\angstrom{} from the reference instance. Each pose is therefore scored against every deposited copy and keeps the lowest value, without superposition onto that copy. That copy also supplies the centroid, form, contact and interaction references for the pose, so one reference geometry serves every endpoint. The choice is made per pose rather than per complex.** This is the convention of the source benchmark, whose validity battery loads every deposited instance and reports the lowest RMSD \cite{ref035}. **The remaining 165 complexes hold one copy and are unaffected.**

What the change buys, beyond stating the fact before the rule:

- **The 36.0 Å median** tells the reader the copies are separate chains or sites rather than near-duplicates in one pocket, which is what makes the choice consequential. Both figures already appear at `:888` and in Appendix F, so nothing new is asserted.
- **"per pose rather than per complex"** states the point most readers get wrong. Two poses of the same complex can be scored against different copies, which was previously carried only by the words "for that pose".
- **"so one reference geometry serves every endpoint"** supplies the reason the same copy governs centroid, form, contact and interaction. As written it read as an arbitrary stipulation rather than a coherence requirement.
- **"The remaining 165"** lands as the complement of 138 instead of arriving unannounced.
- A duplicated "computed in place without superposition", already stated before the equation, was absorbed into the new scoring sentence.

This addresses the readability half of **M7**. The evidential half is unchanged and still open: the crediting copies are classified for site equivalence only for the four alternates named in Appendix F, and the single-instance twins (36.6 / 34.0 / 18.2 per cent) still appear nowhere in the body or either abstract.

Build: `latexmk -pdf` exit 0, no errors, no undefined references, **146 pages** — the roughly four added lines absorbed without a page shift. Both abstracts complete on one page each.

---

## Addendum 11 — the two-tier gate now carries its reason (2026-09-09)

Investigated why the three-part gate (PB-valid AND in-place ≤ 2 Å AND Kabsch ≤ 1 Å) selects the variants while the two-part gate (PB-valid AND ≤ 2 Å) reports every downstream number. Two independent sweeps, each adversarially verified, agreed: **the reason was nowhere stated**. `:199` gave only the consequence, "which admits more complexes", and §3.2.2 stated the practice but not the rationale.

**Two sentences that answered it had been deleted from §3.2.2 since HEAD** — "That bound was fixed before analysis as a conformer-level counterpart to the cited 2\angstrom{} in-place threshold, half its value and **applied to only one of the two components of the in-place deviation**" and "**No published convention is claimed for it.**" The first is the non-orthogonality point, which is the single clearest answer to the question; the second is the candour that makes the two-tier design read as principled rather than opportunistic; and a third clause, "so the sensitivity to this choice can be read directly", had gone with them. All three restored.

**Reason added at `:199`:**

> Form is left out of that endpoint because it is a component of the in-place deviation rather than an independent axis, so a conjunction would gate the same error twice and would attach an uncited threshold to a cited one. The three-part twin is printed at four depths in Table~\ref{tab:reference_convention_sensitivity}, where the ordering between tools is unchanged.

The second sentence is the load-bearing defence: Table 23's Combined column gives AutoDock 122/164/175/175, DiffDock 98/130/144/148, EquiBind 38/52/57/57, so the gate choice moves the levels but not the ranking.

### The strongest surviving objection, verified

Sharper than "the switch flatters every number" (it does: 175→213, 144→179, 57→83) is that **the two gates disagree about which variant was carried forward**. DiffDock + smina was selected on the three-part gate at top-15, 144 to 143. But Appendix G `:871` reports the same two arms on the *reporting* gate at rank-1 as **43.6 % smina against 44.2 % gnina**, and Table 1's three-part column over the full pool gives gnina 149 to smina 148. The arm supplying every reported DiffDock number is therefore not the arm the reporting endpoint would have chosen. That defeats the usual defence that a stricter screen is a conservative version of the same ordering, and it needs no imputation of motive.

Secondary, also verified: the winner's-curse correction is computed "on the three-part selection endpoint at top-15" (App. C.2) and then carried into Limitations as the caveat on two-part numbers. Where it was recomputed on the validity-aware gate it roughly doubled, to 1.4 and 1.6 points against 0.4–1.0. Neither point is currently answered in the document.

**The strongest available defence is the thesis's own rule, one tier down**, and it is still not invoked for form. Appendix G `:1062`: "The validity-aware set is a strict subset of the near-native set, so the gap is one-directional and is reported as an effect size rather than tested." *Combined* stands to *Recovery* in exactly that relation, so the rule prescribes the two-tier design and Table 23's Combined column is the effect size it demands. One cross-reference would convert the design from unexplained to principled by the thesis's own stated standard.

### Build note

The `.aux` file had been left corrupted by an earlier interrupted run ("File ended while scanning use of `\@newl@bel`"), producing hundreds of spurious undefined references. Cleaning the aux family and rebuilding resolved it. Final state verified on the last pass rather than on the multi-pass log: **0 undefined references, 0 undefined citations, no LaTeX errors, 146 pages, and no "??" on any page.** The `!`-prefixed lines that remain in a full log are `tocbasic` number-width warnings, not errors — do not count them with `grep -c '^!'`.

---

## Addendum 12 — the gate disagreement reframed as a defence (2026-09-09)

Addendum 11 filed the two gates disagreeing about the carried-forward variant as the strongest objection to the two-tier design. That framing was wrong, or at least half wrong, and the correction came from the author.

**Selecting on one metric and reporting on another decouples the reported number from the selection criterion.** The winner's curse transfers only to the extent that the two metrics rank the candidates identically. The disagreement is therefore the *evidence* that they do not, which is exactly the condition under which the design protects the headline rather than inflating it.

For DiffDock it demonstrably did. On the reporting gate at rank-1 the arms stand at **43.6 % smina against 44.2 % gnina** (App. G), and smina was carried forward. The reported DiffDock figure is thus **not** the maximum of its family on the metric it is reported on — it sits 0.6 points below it, and carries no winner's curse on that endpoint at all.

**The residual is narrower than what Addendum 11 filed.** The protection is asymmetric. Checked against Table 8, the eight AutoDock arms give 92, 109, 129, 123, 141, 123, 130 and **149** at rank-1 on the reporting gate, and 149 is the selected arm. So AutoDock's carried arm *is* its reporting-gate maximum while DiffDock's is not, and the asymmetry widens the headline gap: 49.2 against 43.6 is +5.6, where a reporting-gate-optimal DiffDock arm would give +5.0.

Not "the gates disagree, so the design is suspect" but "the decoupling protects two arms of three, and the one it does not protect is the one that wins". Added to Ch. 7's selection-optimism paragraph:

> That correction transfers only in part, because the two gates do not rank the variants identically. DiffDock + smina was carried forward on the selection gate, whereas DiffDock + gnina stands higher on the reporting gate at rank-1, at 44.2\% against 43.6\% (Appendix~\ref{ranking-recovery-statistics}), so the reported DiffDock figure is not its family maximum on the endpoint it is reported on. The selected AutoDock arm is the maximum of its eight on both gates. The decoupling therefore protects the DiffDock and EquiBind figures and not the AutoDock one, and the rank-1 gap would be 5.0 rather than 5.6 points against a reporting-gate DiffDock arm.

This states the defence and concedes the residual in the same breath, and it costs 0.6 of a percentage point on a gap the thesis already reports as unresolved.

**Page count: 147, up one.** The addition runs about six lines and Chapter 7 now ends one page later. The extra page is not sparse — pages 55 to 58 are all fully set at 708 to 711 pt, the Bibliography still starts immediately after, and no heading is orphaned. Both abstracts remain complete on one page each. Holding 146 would mean moving the detail into Appendix C.2 and leaving a pointer in Ch. 7, which buries a defence where an examiner is least likely to read it. Recommendation: keep the page.

Final pass: 0 undefined references, 0 undefined citations, no LaTeX errors.
