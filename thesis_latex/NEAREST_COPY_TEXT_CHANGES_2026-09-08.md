# Nearest-copy endpoint: proposed thesis text changes (2026-09-08)

**Purpose.** The thesis author has forbidden every in-place edit of the thesis. This file therefore collects each prose change that `PLAN_nearest_copy_primary_endpoint_2026-09-07_v2.md` calls for as a PROPOSED edit, with the exact current text quoted verbatim from the live file and its line number as of 2026-09-08, and the exact proposed replacement, for the author to apply or reject. Nothing in `body_main_short.tex`, `body_appendix_short.tex` or `Thesis_short.tex` has been touched.

**Provenance of the quotations.** Every "Current" block below was extracted mechanically from the named file and line at generation time by a throw-away generator script (`gen_text_changes.py`, kept in the session scratchpad and deliberately not added to the repository). The generator aborts if a quotation is not a substring of the named line, so each block is verbatim by construction. Line numbers are those of the working tree on 2026-09-08 and agree with the plan's anchors (the plan's line numbers are "as of 22:30 on 2026-09-07"); one anchor drifted and is flagged where it occurs (the interaction values the plan places at App. :545-553 sit at :543 in the live file). Locate every edit by the quoted text, never by the line.

**Plan sections read.** Section 0 (D1-D17, in particular D11 and D14), Section 1 (1.1-1.11), Phase 0.3 and 0.4, Phase 5 (5.1-5.6 and both drafts), Section 3 (3.1-3.3), Section 7, Section 8.

**Conventions in this file.**
- Part A holds corrections that are true today, independent of whether the endpoint is switched.
- Part B holds edits that belong to the endpoint switch and are to be applied only if the switch is adopted.
- A number in square brackets, e.g. `[7,565]`, is a preview derived from the probe CSVs (plan D16). Where the plan says recompute, the text says *value from the Phase 3 rebuild*. No bracketed number is to be transcribed into the thesis before the rebuild has reproduced it.
- House style of every proposed sentence: no mid-sentence semicolons, colons or clause-dashes, bare tool names in the body after the declaration at `body_main_short.tex:199`, full variant names where the surrounding appendix paragraph uses them, "rank-1", "top-15", "ranking depth".
- File abbreviations in headings: M = `body_main_short.tex`, A = `body_appendix_short.tex`, T = `Thesis_short.tex`.

---

## Part A. Corrections valid today, independent of the endpoint switch

> **Status 2026-09-08 14:05.** Part A (A1, A1b, A2, A3, A4, A6) was applied verbatim to the shipped `thesis_latex/` on the user's instruction (decision 10, "apply exactly Part A to thesis_latex/"); A5 stays a flag. The shipped PDF was rebuilt (142 pages, 0 errors) and probed for the new strings. No yaml value changes.

Source: plan Phase 0.3, Phase 0.4 and Section 8 ("Pre-existing defects in the current thesis found by the audits, independent of the switch"). Each item was re-verified on 2026-09-08 against the live files and the read-only canonical table before being written down.

#### A1. M:113 — "adds three recovered complexes" → "two"

*Status:* apply today

**Current (`body_main_short.tex:113`):**
```latex
Substituting the fully symmetrised value into the rank-1 gate adds three recovered complexes across the three tools and removes none (Appendix~\ref{supplementary-statistics}).
```

**Proposed:**
```latex
Substituting the fully symmetrised value into the rank-1 gate adds two recovered complexes across the three tools and removes none (Appendix~\ref{supplementary-statistics}).
```

**Notes:** Re-verified 2026-09-08 on the read-only canonical table `posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report/per_pose_metrics.csv` (241,913 rows). Rank-1 PB-valid poses with `rmsd` > 2 and `pb_rmsd` ≤ 2: AutoDock Vina + gnina (`autodock_mgltools_exh128_gnina`, `rank == 1`) gains 7PGX_FMN only, DiffDock + smina (`diffdock_smina`, `rank == 1`) gains none, EquiBind + gnina (`equibind_unguided_gnina`, rank-1 taken as the lowest `gnina_affinity` per complex because the `rank` column is the 999 sentinel for that arm) gains 6XBO_5MC only. No pose runs the other way in any arm. Total two, none removed, matching plan 0.3. The same table shows `pb_rmsd` > `rmsd` on 0 of 241,913 rows (largest excess 4.7e-10 Å), which is the "never smaller" claim of the preceding sentence. The count is identical under both conventions (plan 0.3), so this edit belongs in Part A. Two further observations. (i) The appendix already prints the correct count: the Table 27 footnote at `body_appendix_short.tex:1387` says "moves two of the 909 tool-by-complex decisions, one for AutoDock Vina + gnina and one for EquiBind + gnina, with none for DiffDock + smina", so today the body contradicts its own appendix. (ii) Deviation from the plan: the plan says the `\ref` "points at nothing". In the live file `\label{supplementary-statistics}` exists at `body_appendix_short.tex:1067` and Table 27 (`:1377-1387`) lies inside that chapter, so the reference resolves to the chapter that carries the supporting footnote. What is missing is a sentence in running prose, which A1b supplies. Registering the `pb_rmsd` gate as a yaml prose check (plan 0.3) is a harness task outside this file.

#### A1b. A:1196 — new supporting sentence in the paragraph "What that requirement costs the headline endpoint"

*Status:* apply today

**Current (`body_appendix_short.tex:1196`):**
```latex
That is a span of 0.66 to 1.32 percentage points. The recovery ordering, the paired differences and the intervals reported in the Results are therefore set almost entirely by the 2\angstrom{} placement criterion.
```

**Proposed:**
```latex
That is a span of 0.66 to 1.32 percentage points. The recovery ordering, the paired differences and the intervals reported in the Results are therefore set almost entirely by the 2\angstrom{} placement criterion. The choice of RMSD implementation moves the same endpoint even less. Substituting the fully symmetrised PoseBusters RMSD for the study\textquotesingle s own value in the rank-1 gate recovers two further complexes, 7PGX\_FMN for AutoDock Vina + gnina and 6XBO\_5MC for unguided EquiBind + gnina, and removes none, because the PoseBusters value never exceeds the study\textquotesingle s own on any of the 241,913 scored poses of the canonical table.
```

**Notes:** Insertion of two sentences after the quoted second sentence; the rest of the paragraph is unchanged. The two names and the zero-excess fact were verified as described under A1. "241,913" is the row count of the canonical table over all 27 method keys (it includes Uni-Dock2 and the ladder arms); if the author prefers a figure restricted to the three headline arms it must be re-counted. Under Part B (nearest-copy convention) the plan states the count stays two with the same complexes, but the names must be re-read from the Phase 3 rebuild before this sentence is kept in a switched build. The sentence uses the paragraph's own variant names ("AutoDock Vina + gnina", "unguided EquiBind + gnina"). If the author would rather keep the Table 27 footnote as the sole carrier, the alternative is to change the `\ref` at M:113 to point at Table~\ref{tab:results-placement-form} instead of the chapter.

#### A2. M:396 — "p_holm ≤ 4e-10" → "≤ 9e-9"

*Status:* apply today

**Current (`body_main_short.tex:396`):**
```latex
AutoDock and DiffDock exceed EquiBind at every depth (p\_holm \ensuremath{\leq} 4e-10).
```

**Proposed:**
```latex
AutoDock and DiffDock exceed EquiBind at every depth (p\_holm \ensuremath{\leq} 9e-9).
```

**Notes:** Evidence: the companion report of the Kabsch depth figure, `posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report/18_topn_within_thresholds_kabsch_pbvalid_depths_report.txt`, lines 49-50, prints `AutoDock* > EquiBind: discordant +79/−21, p_holm=8.997e-09` and `DiffDock* > EquiBind: discordant +76/−19, p_holm=8.997e-09` at top-1 (1 Å, PB-valid). The deeper depths are smaller (top-15 1.634e-11 for both, top-30 2.427e-11 and 3.556e-11), so the largest corrected p over the six tool-by-depth cells is 8.997e-09 and the bound "≤ 4e-10" is false today. Under Part B the plan lists this cell as recompute after the rebuild (Section 3.2 row 396), so the 9e-9 is the Part A value only.

#### A3. M:414 — "median in-place RMSD of 2.9" → "3.0"

*Status:* apply today

**Current (`body_main_short.tex:414`):**
```latex
Its rank-1 median in-place RMSD of 2.9\angstrom{} already exceeds the near-native threshold.
```

**Proposed:**
```latex
Its rank-1 median in-place RMSD of 3.0\angstrom{} already exceeds the near-native threshold.
```

**Notes:** Evidence inside the thesis itself: Table 27 at `body_appendix_short.tex:1383` prints the EquiBind + gnina rank-1 near-site median as 2.977 Å, which rounds to 3.0, not 2.9 (plan 0.3 and Section 1.6 call it a pre-existing truncation). Under Part B the same cell becomes [3.2] (preview 3.239 Å), value from the Phase 3 rebuild.

#### A4. A:871 — lead-in "validity-aware recovery (34.0 % and 34.3 %)" → name the gate each table prints

*Status:* apply today

**Current (`body_appendix_short.tex:871`):**
```latex
These tables carry both DiffDock refiners, so the smina arm the Results select is printed beside the gnina one, at 34.0\% and 34.3\% rank-1 validity-aware recovery respectively. The following tables give depth-resolved statistics for the variants displayed below.
```

**Proposed:**
```latex
These tables carry both DiffDock refiners, so the smina arm the Results select is printed beside the gnina one. Table~\ref{tab:appendix-top-k-recovery} prints both gates, near-native placement alone and the validity-aware conjunction the Results use, and the two arms stand at 35.3\% against 35.6\% rank-1 on the first and at 34.0\% against 34.3\% on the second. Every paired test in this chapter is taken on the near-native gate alone. The following tables give depth-resolved statistics for the variants displayed below.
```

**Notes:** Deviation from the plan's description, recorded here rather than silently absorbed. The plan (Section 8, Section 3.3 row 871) says the lead-in "calls its tables validity-aware while they print the near-native gate". Checked against the live file: the first table of the chapter, `tab:appendix-top-k-recovery` at `:881-921`, prints TWO columns, "Near-native% [95% CI]" and "Validity-aware% [95% CI]", and the 34.0 % and 34.3 % are exactly its validity-aware rank-1 cells for DiffDock + smina (`:905`) and DiffDock + gnina (`:909`), whose near-native cells are 35.3 % and 35.6 %. The two McNemar tables (`tab:appendix-refinement-mcnemar` at `:935`, `tab:appendix-cross-tool-mcnemar` at `:1014`) are tested on the near-native gate, which the same paragraph already states further down ("each taken on the near-native endpoint", "The validity-aware endpoint is not part of the family"). So the defect is narrower than the plan states: the opening sentence names one gate and one table's cells for a chapter whose tests run on the other gate. The proposal names both and says which one the tests use. Under Part B the four percentages become [45.2] % and *value from the Phase 3 rebuild* on the near-native gate and [43.6] % and [44.2] % on the validity-aware gate (plan Section 3.3 row 871 and Section 3.2 row 191-199).

#### A5. A:165 — cofactor / free-ion split 15 / 63, flag only

*Status:* flag, do not edit yet

**Current (`body_appendix_short.tex:165`):**
```latex
Fifteen of the 78 involve a metal-bearing cofactor such as haem, an iron-sulfur cluster or cobalamin, and the remaining 63 involve a free ion.
```

**Proposed:** no replacement text yet (see notes).

**Notes:** To reconcile with the registered rule (16 / 62). Plan D4 and Phase 0.4: the residue-level rule that reproduces the printed 73 / 78 / 81 stratum sizes yields 16 metal-bearing-cofactor complexes and 62 free-ion complexes, whereas the sentence prints 15 / 63. No number is proposed here. The registered generator `Scripts/Analysis/metal_stratum.py` (Phase 0.4, not yet written) is to decide the split, and the sentence is then edited to whatever that generator asserts. The two candidate outcomes are (a) the sentence reads "Sixteen of the 78 ... and the remaining 62" or (b) the rule is amended to reproduce 15 / 63 and the sentence stands. The 78 itself is not in question.

#### A6. A:664 — footnote roles of the four benchmark files

*Status:* apply today

**Current (`body_appendix_short.tex:664`):**
```latex
\footnote{For every entry the benchmark provides four files, namely the receptor stripped of solvent but retaining its cofactors (*\_protein.pdb), all crystallographic instances of the ligand of interest (*\_ligands.sdf), the single chosen crystal instance that serves as the ground-truth answer pose (*\_ligand.sdf), and one computer-generated starting conformer produced with RDKit's ETKDGv3 algorithm followed by UFF minimisation (*\_ligand\_start\_conf.sdf).}
```

**Proposed:**
```latex
\footnote{For every entry the benchmark provides four files, namely the protein structure without the ligand of interest, without solvent and with all cofactors (*\_protein.pdb), all crystallographic instances of the ligand of interest (*\_ligands.sdf), one of those instances, which the benchmark supplies to mark the binding site for docking methods that require one and which this study scores against (*\_ligand.sdf), and one computer-generated starting conformer produced with RDKit's ETKDGv3 algorithm followed by UFF minimisation (*\_ligand\_start\_conf.sdf).}
```

**Notes:** Documented roles taken from the dataset README, `Data/PoseBuster Benchmark Set/README.txt` lines 51-63 (read 2026-09-08): (1) `_protein.pdb` "The protein structure without the ligand of interest without solvents and with all cofactors", (2) `_ligands.sdf` "All instances of the ligand of interest", (3) `_ligand.sdf` "One of the instances of the ligand of interest. This crystal pose marks the binding site for those docking methods that require a binding site", (4) `_ligand_start_conf.sdf` the ETKDGv3 + UFF conformer. "Without the ligand of interest" means every copy is removed from the receptor (plan, External convention). The current footnote's "single chosen crystal instance that serves as the ground-truth answer pose" is the study's own use of the file, not the benchmark's documented role, so the proposal keeps the study's use as a separate clause ("and which this study scores against"), which is true today. Under Part B that clause becomes "and which is the first record of the multi-instance file" (the reference is record 0 of `_ligands.sdf` in 308 / 308 entries, plan Data facts). The following sentence of the same line, "The crystallographic pose is the experimentally validated reference for scoring every docked pose", is true under both conventions and needs no change.


---

## Part B. Edits that belong to the endpoint switch (apply only if the switch is adopted)

Everything below presumes plan decisions D1-D17 are adopted and the Phase 3 rebuild has run. Previews in square brackets come from the probe CSVs (plan Sections 1.1-1.11) and are NOT to be transcribed before the rebuild reproduces them (D16).

### B1. Methods redefinition (plan D11, Phase 5.1) at M:109-115

#### B1a. M:109 — definition sentence

**Current (`body_main_short.tex:109`):**
```latex
Near-nativeness is the symmetry-corrected heavy-atom root-mean-square deviation (RMSD) between a docked pose and the crystallographic ligand:
```

**Proposed:**
```latex
Near-nativeness is the symmetry-corrected heavy-atom root-mean-square deviation (RMSD) between a docked pose and the nearest of all crystallographic copies of the ligand deposited with the structure, computed in place without superposition:
```

**Notes:** First sentence of the plan's D11 draft, adapted so that the displayed equation at `:111` still follows a colon that introduces it (the colon before a displayed formula is not prose punctuation). The remaining sentences of the draft go into `:113` (B1b) because the equation sits between the two lines.

#### B1b. M:113 — copy rule, external convention, sensitivity arm, and the "two" of A1

**Current (`body_main_short.tex:113`):**
```latex
RMSD is computed in place without superposition onto the crystal ligand. A pose is accurate at the conventional threshold of RMSD \ensuremath{\leq} 2\angstrom{} \cite{ref022}, \cite{ref050}.
```

**Proposed:**
```latex
RMSD is computed in place without superposition onto that copy. The copy with the lowest RMSD is chosen for each pose and supplies the centroid, form, contact and interaction references for that pose. This is the convention of the source benchmark, whose validity battery loads every deposited instance and reports the lowest RMSD \cite{ref035}. The 165 single-copy complexes are unaffected. The single-instance value is retained and reported as a sensitivity arm in Appendix~\ref{cohort-and-exclusions}. A pose is accurate at the conventional threshold of RMSD \ensuremath{\leq} 2\angstrom{} \cite{ref022}, \cite{ref050}.
```

**Notes:** Sentences two to five are the D11 draft verbatim apart from the citation and the appendix label. The label `cohort-and-exclusions` (`body_appendix_short.tex:1072`) is where the plan places the new D9 sensitivity table ("beside :1074"); if the table is placed elsewhere the `\ref` follows it. External-convention sources per D11: README :55-59, `redock.yml:307 load_all: True`, Chem. Sci. 2024 Section 2.3 "closest crystallographic ligand"; `\cite{ref035}` is the PoseBusters paper already used at `:113`. The same line also carries the A1 change ("three" → "two"), which applies in both parts. The sentence "The reported value is therefore never smaller than a fully symmetrised RMSD on any of the poses analysed here" stays true under the nearest-copy convention because both implementations minimise over the same copies (PoseBusters `rmsd.py:52-70` argmin over conformers), but it must be re-verified on the rebuilt table exactly as A1 verified it today.

#### B1c. M:115 — headline criterion names the copy

**Current (`body_main_short.tex:115`):**
```latex
The headline crystal-referenced criterion, valid-near-native, requires both RMSD \ensuremath{\leq} 2\angstrom{} and PB-validity.
```

**Proposed:**
```latex
The headline crystal-referenced criterion, valid-near-native, requires both RMSD \ensuremath{\leq} 2\angstrom{} to that nearest copy and PB-validity.
```

**Notes:** Optional if B1b is adopted, since the definition already carries the copy rule. Listed because plan Phase 5.1 names `:115`.


### B2. Convention paragraph at A:666 (plan Phase 5 draft, D3, Section 1.3)

#### B2. A:666 — replace the two single-instance sentences

**Current (`body_appendix_short.tex:666`):**
```latex
Where a structure contains multiple crystallographic copies of the same ligand, the comparison is made against the single deposited reference instance rather than against the closest copy. The symmetry correction therefore resolves symmetry within that one instance and not equivalence between copies.
```

**Proposed:**
```latex
In 143 of the 308 entries the deposited file holds more than one copy of the ligand, 138 of them among the 303 analysed complexes, and the nearest other copy lies a median 36.0\angstrom{} from the reference over those 143. Every pose is scored against the nearest deposited copy, and the form, centroid and interaction references follow that copy. The source benchmark scored every method this way. Its 25\angstrom{} cube centred on the reference kept almost every other copy outside Vina\textquotesingle s search, whereas the whole-protein box used here contains every one of them, so the convention is live for all three tools in this study. The choice of copy can lower a form count where the copies differ in conformation, and it does so for one complex in each of [five] printed variants of Table~\ref{tab:results-pose-production}. Four alternates do not share the reference pocket (7A9E\_R4W, 7TUO\_KL9, 7VKZ\_NOJ, 7Z1Q\_NIO) and no complex gains recovery through any of them.
```

**Notes:** The first two sentences of the current paragraph ("The experimentally resolved ligand conformations provide the reference ..." and the leakage sentence) stay as they are; only the quoted third and fourth sentences are replaced. The replacement is the plan's Phase 5 draft for `:666` with three adjustments. (i) LaTeX: `\angstrom{}`, `\textquotesingle`, escaped underscores. (ii) Deviation from the draft: the draft says "one complex in each of three reported variants", but plan Section 1.3 and the Section 8 refutation ledger ("three arms lose a Form complex → five printed rows") establish that the Form block loses a complex in FIVE printed Table 1 rows (DiffDock raw, DiffDock + smina, EquiBind + gnina, fpocket raw, P2Rank raw), so the draft was not updated after its own audit. Written here as "[five]" and tied to the table, *value from the Phase 3 rebuild*. (iii) D3 wording is as the plan requires ("no complex gains recovery through any of them", not the refuted "no reported pose recovers any of them", since three PB-valid ladder poses sit 1.70-1.75 Å from 7Z1Q_NIO's alternate). The 36.0 Å median is over the 143 multi-copy ids; 36.3 over the 138 analysed and 40.4 over all 211 alternates are the other denominators and must not be mixed (plan Data facts). "Almost every other copy" per D11 (8 of 211 alternate centroids lie inside a 25 Å cube on the reference).


### B3. Limitations sentence at M:868 (plan Phase 5.3, D4)

#### B3. M:868 — metal stratum wording, stratum contrasts, and the convention disclosure

**Current (`body_main_short.tex:868`):**
```latex
A metal ion or metal-containing cofactor lies within 5\angstrom{} of the crystal ligand in 78 of the 303 calibration complexes, so recovery across that quarter of the set is measured without a coordinating partner the crystal pose depends on. The AutoDock-over-DiffDock top-15 contrast is +6.7 percentage points on the 225 metal-free complexes against +20.5 on the 78 metal-adjacent ones, so the pooled contrast averages a narrow gap with a much wider one (Appendix~\ref{autodock-vina-2}).
```

**Proposed:**
```latex
A residue carrying a metal ion or a metal-containing cofactor lies within 5\angstrom{} of the deposited reference instance of the ligand in 78 of the 303 calibration complexes, so recovery across that quarter of the set is measured without a coordinating partner the crystal pose depends on. The AutoDock-over-DiffDock top-15 contrast is [+8.4] percentage points on the 225 metal-free complexes (p = [0.056]) against [+19.2] on the 78 metal-adjacent ones (p = [0.028]), so the pooled contrast averages a narrow gap with a much wider one (Appendix~\ref{autodock-vina-2}). Near-nativeness is scored against the nearest of the deposited copies of the ligand, which matters in 143 of the 308 entries and 138 of the 303 analysed complexes, where the nearest other copy lies a median 36.0\angstrom{} from the reference. The source benchmark scored every method the same way, and the single-instance value is kept as a sensitivity arm in Appendix~\ref{cohort-and-exclusions}.
```

**Notes:** Three changes in one sentence pair plus two new sentences. (i) D4 residue-level rule wording ("a residue carrying ... of the deposited reference instance"), since membership is a complex-level property of the reference instance and does not follow the scored copy. (ii) Stratum contrasts [+8.4] (p [0.056]) and [+19.2] (p [0.028]) are previews from plan Section 3.2 row 868 and Phase 5.4, *value from the Phase 3 rebuild*. (iii) The convention disclosure sentence the plan asks for at `:868` (Phase 5.3: 143 / 308, 138 / 303, median 36.0 Å over the 143, source benchmark scored the same way, single-instance kept as sensitivity). The rest of the line is unchanged. The same D4 wording applies at M:88 (see B5).


### B4. Abstract and Kurzfassung (plan Section 1.11, Phase 5.5, Section 3.1)

The Kurzfassung page is full (memory: zero slack), so the German swaps need a compensating cut of at least the added length. Character deltas below are counted on the LaTeX source; the plan's check is `wc -m` on the page and a pymupdf raster showing no overflow.

#### B4a. T:61 — Abstract endpoint definition (W optional)

**Current (`Thesis_short.tex:61`):**
```latex
A calibration benchmark of 303 co-crystallised complexes permits assessment against known poses under a primary endpoint combining PoseBusters validity with a heavy-atom root-mean-square deviation of at most 2\angstrom{}.
```

**Proposed:**
```latex
A calibration benchmark of 303 co-crystallised complexes permits assessment against known poses under a primary endpoint combining PoseBusters validity with a heavy-atom root-mean-square deviation of at most 2\angstrom{} to the nearest deposited copy of the ligand.
```

**Notes:** Plan Section 3.1 row 61 marks this "W optional" (the one optional convention clause of Section 1.11). The English abstract page has more slack than the German one, but the author should raster-check it too.

#### B4b. T:63 — Abstract numbers

**Current (`Thesis_short.tex:63`):**
```latex
Raw poses from the learned workflows were largely inadmissible, at 24.4\% pooled validity for DiffDock and 2.9\% for EquiBind against 99.6\% for AutoDock. Local optimisation repaired this almost categorically, lifting the median per-complex validity yield from 13.3\% to 93.3\% and from 0\% to 80.0\%. It corrected local geometry but recovered no binding mode that was unsampled or grossly misplaced. Rank-1 recovery of a pose both valid and within 2\angstrom{} reached 36.6\% of complexes for AutoDock, 34.0\% for DiffDock and 18.2\% for EquiBind, leaving the leading pair unresolved at 2.6 points with a 95\% interval of \ensuremath{-}4.2 to +9.4. Generation rather than ranking dominated failure, with about six in ten rank-1 failures for the leading pair containing no qualifying pose, though inspecting fifteen poses instead of one still raised recovery to 65.3\% and 55.1\%. Unlike the rank-1 contrast, the gaps at top-5 and top-15 are statistically resolved. On a single device-occupancy basis the configured AutoDock workflow was the most expensive per qualifying pose, and EquiBind's low cost was conditional on recovering only 80 of 303 complexes.
```

**Proposed:**
```latex
Raw poses from the learned workflows were largely inadmissible, at 24.4\% pooled validity for DiffDock and 2.9\% for EquiBind against 99.6\% for AutoDock. Local optimisation repaired this almost categorically, lifting the median per-complex validity yield from 13.3\% to 93.3\% and from 0\% to 80.0\%. It corrected local geometry but recovered no binding mode that was unsampled or grossly misplaced. Rank-1 recovery of a pose both valid and within 2\angstrom{} reached [49.2]\% of complexes for AutoDock, [43.6]\% for DiffDock and [18.8]\% for EquiBind, leaving the leading pair unresolved at [5.6] points with a 95\% interval of \ensuremath{-}[1.7] to +[12.8]. Generation rather than ranking dominated failure, with about six to seven in ten rank-1 failures for the leading pair containing no qualifying pose, though inspecting fifteen poses instead of one still raised recovery to [70.3]\% and [59.1]\%. Unlike the rank-1 contrast, the gaps at top-5 and top-15 are statistically resolved. On a single device-occupancy basis the configured AutoDock workflow was the most expensive per qualifying pose, and EquiBind's low cost was conditional on recovering only [83] of 303 complexes.
```

**Notes:** Five number swaps and one wording swap per plan Section 1.11 (36.6 / 34.0 / 18.2 → 49.2 / 43.6 / 18.8; +2.6 → +5.6; −4.2..+9.4 → −1.7..+12.8; 65.3 / 55.1 → 70.3 / 59.1; 80 → 83; "six in ten" → "six to seven in ten", since rank-1 failures with no qualifying pose anywhere become [57.8] % and [70.2] % for the leading pair). All bracketed values are previews, *value from the Phase 3 rebuild*. Two clauses are kept pending checks the plan names. (i) "It corrected local geometry but recovered no binding mode that was unsampled or grossly misplaced" is to be re-verified (DiffDock raw → smina rank-1 near-native 132 → 137 under nearest, plan 1.11), and stays as written unless the rebuild contradicts it. (ii) "the configured AutoDock workflow was the most expensive per qualifying pose" is pending plan Section 1.9 (the median succeeding AutoDock complex yields two qualifying poses instead of one, so every cost median is re-derived); the clause changes only if the ordering flips. "Unlike the rank-1 contrast, the gaps at top-5 and top-15 are statistically resolved" holds (top-5 identical, top-15 p 0.004).

#### B4c. T:91-100 — Kurzfassung numbers (quoted as the full lines)

**Current (`Thesis_short.tex:91-100`):**
```latex
Posen innerhalb von 2\angstrom{} bei 36,6\% für AutoDock, 34,0\% für DiffDock und 18,2\%
für EquiBind, womit das führende Paar bei 2,6 Prozentpunkten unaufgelöst blieb,
95\%-Intervall
\ensuremath{-}4,2 bis +9,4. Nicht die Reihung, sondern die Posenerzeugung dominierte das
Versagen, denn etwa sechs von zehn Fehlschlägen am Rang 1 enthielten beim führenden Paar
keine qualifizierende Pose. Fünfzehn statt einer Pose erhöhten die Wiederfindungsrate
dennoch auf 65,3\% und 55,1\%. Anders als am Rang 1 sind die Abstände bei fünf und
fünfzehn Posen statistisch gesichert. In einheitlicher Belegungsrechnung war der
AutoDock-Ablauf am teuersten pro qualifizierender Pose, EquiBinds niedrige Kosten galten
nur für 80 von 303 Komplexen.
```

**Proposed:**
```latex
Posen innerhalb von 2\angstrom{} bei [49,2]\% für AutoDock, [43,6]\% für DiffDock und [18,8]\%
für EquiBind, womit das führende Paar bei [5,6] Prozentpunkten unaufgelöst blieb,
95\%-Intervall
\ensuremath{-}[1,7] bis +[12,8]. Nicht die Reihung, sondern die Posenerzeugung dominierte das
Versagen, denn etwa sechs bis sieben von zehn Fehlschlägen am Rang 1 enthielten beim führenden Paar
keine qualifizierende Pose. Fünfzehn statt einer Pose erhöhten die Wiederfindungsrate
dennoch auf [70,3]\% und [59,1]\%. Anders als am Rang 1 sind die Abstände bei fünf und
fünfzehn Posen statistisch gesichert. In einheitlicher Belegungsrechnung war der
AutoDock-Ablauf am teuersten pro qualifizierender Pose, EquiBinds niedrige Kosten galten
nur für [83] von 303 Komplexen.
```

**Notes:** Swaps per plan Section 1.11 and Section 3.1 rows 91 / 92 / 94 / 95 / 97 / 100. Length accounting on the source text (brackets removed): 36,6 / 34,0 / 18,2 → 49,2 / 43,6 / 18,8 (0 chars), 2,6 → 5,6 (0), "−4,2 bis +9,4" → "−1,7 bis +12,8" (+1), "sechs von zehn" → "sechs bis sieben von zehn" (+11), 65,3 / 55,1 → 70,3 / 59,1 (0), 80 → 83 (0). Net +12 characters, all on the two lines `:94-95`, so the page needs a cut of at least 12 characters that does not change meaning (B4d). The cost clause at `:98-99` ("war der AutoDock-Ablauf am teuersten pro qualifizierender Pose") is W only if plan Section 1.9 flips the ordering. Line breaks inside the paragraph are a property of the source and the author may re-flow them.

#### B4d. T:89-90 — compensating cut candidate for the Kurzfassung page

**Current (`Thesis_short.tex:89-90`):**
```latex
auf 93,3\% sowie von 0\% auf 80,0\%. Sie korrigierte lokale Geometrie, nicht jedoch
fehlende oder grob fehlplatzierte Bindungsmodi. Am Rang 1 lag die Wiederfindung gültiger
```

**Proposed:**
```latex
auf 93,3\% sowie von 0\% auf 80,0\%. Sie korrigierte lokale Geometrie, nicht
fehlende oder grob fehlplatzierte Bindungsmodi. Am Rang 1 lag die Wiederfindung gültiger
```

**Notes:** The plan names the optimisation clause at `:89-90` as the candidate ("length-neutral in meaning"). The clause "Sie korrigierte lokale Geometrie, nicht jedoch fehlende oder grob fehlplatzierte Bindungsmodi." can be shortened without loss by "nicht jedoch" → "nicht aber" (−2 characters) or by dropping "jedoch" (−7). Neither alone reaches the required 12, so a second meaning-neutral cut is needed. Candidate on `:96-97`: "Fünfzehn statt einer Pose erhöhten die Wiederfindungsrate dennoch auf" → "Fünfzehn statt einer Pose hoben die Wiederfindung dennoch auf" (−8). Dropping "jedoch" (−7) plus the `:96-97` cut (−8) gives −15 ≥ 12. Alternative single cut with a small change of emphasis: drop "grob" as well (−5 more), which the author may not want because "grossly misplaced" mirrors the English abstract. Whatever is chosen, the page must be re-rastered (plan 5.5) because line breaks, not character counts, decide overflow. The Kurzfassung definition at `:80` is "no change" per plan Section 3.1.


### B5. Every Section 3 row marked W (wording), with current text and proposed replacement

Index first, details below. \"W-if\" means the plan makes the wording change conditional. Rows whose Section 3 entry carries only numbers are listed at the end of Part B as value-only rows and are not quoted.

| # | Location | Plan gate | What changes | Replacement given? |
|---|---|---|---|---|
| W1 | M:88 | stratum | D4 residue-rule wording | yes |
| W2 | M:109-115 | def. | D11 redefinition, \"three\" → \"two\" | yes (B1, A1) |
| W3 | M:119 | def. | clustering pocket follows any copy (D5) | yes |
| W4 | M:130 (row 124-133) | def. | form value taken on the copy the in-place minimum selects (D2) | yes |
| W5 | M:158 | def. | 4 Å reach to the nearest copy, any copy reaches (D5) | yes |
| W6 | M:165 | def. | fingerprint reference per copy, union rule (1.8) | yes |
| W7 | M:170 | def. | cost definition names the copy | yes |
| W8 | M:193, 195, 197 | triple / double / pooled | \"level\" not \"gnina higher\", counts | yes, previews |
| W9 | M:264, 277-283 | — | Table 1 optional footnote | yes (optional) |
| W10 | M:324 | valid ∧ ≤ 2 | \"at most four\" → \"five\" | yes |
| W11 | M:326 | def. | Table 2 lead-in names the copy | yes |
| W12 | M:380 | def. | \"as-placed crystal RMSD\" | yes |
| W13 | M:393, 401 | — | Fig 2 / 3 captions | yes for 393, 401 holds |
| W14 | M:396 | valid ∧ Kabsch | tie wording, discordant counts | yes, previews |
| W15 | M:408 | trim | footnote \"crystal-ligand site\", counts | yes, previews |
| W16 | M:422 | trim | near-site interpretive paragraph | no, re-read after rebuild |
| W17 | M:481 (row 481-532) | interaction | copy and union rule stated | yes |
| W18 | M:697 | valid ∧ ≤ 2 | \"at rank-1\" | yes |
| W19 | M:708 | def. | cost definition | yes |
| W20 | M:744, 758 (row 744-760) | cost | \"yields one\" → \"two\", \"nine\" → \"eleven\", W-if ordering | yes / W-if |
| W21 | M:768 | valid ∧ ≤ 2 / cluster | interval, reach clause | yes / rebuild |
| W22 | M:777 | valid ∧ ≤ 2 | \"at most four\" → \"five\" | yes |
| W23 | M:781 | triple / double / near | \"smina ahead by two\", \"at the top-15 selection depth\" | yes |
| W24 | M:788, 857 | valid ∧ ≤ 2 | \"six to seven in ten\", counts | yes |
| W25 | M:834 | valid ∧ ≤ 2 | \"eight points\" → \"twelve points\", numbers | yes |
| W26 | M:843 | valid ∧ ≤ 2 (D14) | Holm family named | yes |
| W27 | M:848-850 | cost | W-if ordering flips | W-if |
| W28 | M:866 | triple | \"eight\" → \"nine\" | yes |
| W29 | M:868 | stratum | D4 wording, disclosure | yes (B3) |
| W30 | A:165 | stratum | crossover sentence, residue rule, numbers | yes, previews |
| W31 | A:167 | def. | every copy inside the box | yes |
| W32 | A:182 | valid ∧ ≤ 2 | depth lead, six of seven, 149 vs 141 | yes, previews |
| W33 | A:195 | — | Fig 19 caption top step | yes, preview |
| W34 | A:227 (row 218-227) | valid ∧ ≤ 2 | Table 8 footnote gate | yes |
| W35 | A:230 | triple (pool) | D10 sentence with nearest values | yes, previews |
| W36 | A:246 | near | plateau argument rewritten | yes, previews |
| W37 | A:267 (row 267-269) | near | P2Rank raw recovers one | yes |
| W38 | A:276 | def. | \"to the crystal ligand\" | yes |
| W39 | A:487 | — | optional `load_all` clause | yes (optional) |
| W40 | A:543, 551 (plan row 545-553) | interaction | chain-keyed fingerprint, per-copy reference | yes / rebuild |
| W41 | A:664 | — | footnote roles | yes (A6, Part B variant) |
| W42 | A:666 | — | convention paragraph | yes (B2) |
| W43 | A:722 | — | 46 % multi-copy share | yes |
| W44 | A:784 | — | D13 clause | yes |
| W45 | A:871 | — | gate lead-in | yes (A4, Part B variant) |
| W46 | A:1181 (row 1172-1181) | near | footnote near-native definition | yes |
| W47 | A:1188, 1196 (row 1188-1196) | near / valid | strata, \"at most five\", span | yes, previews |
| W48 | A:1203 | — | regime clause \"almost every other copy\" | yes |
| W49 | A:1260 (row 1258-1260) | valid ∧ ≤ 2 | \"within 12 points\", numbers | yes, previews |
| W50 | A:1265 | — | Kabsch criterion names the copy | yes |
| W51 | A:1387 (row 1377-1387) | trim | footnote site wording, seven numbers | yes, previews |
| W52 | A:1399 (row 1394-1399) | cluster | Fig 38 caption | yes |
| W53 | A:1412 (row 1407-1415) | interaction | Fig 39 caption | yes |
| W54 | T:61, 98-99 | — | abstract clause optional, Kurzfassung cost W-if | yes (B4a) / W-if |

#### W1. M:88 — stratum sentence in the Calibration Dataset paragraph

**Current (`body_main_short.tex:88`):**
```latex
78 of the 303 analysed complexes place a metal within 5\angstrom{} of the crystal ligand
```

**Proposed:**
```latex
78 of the 303 analysed complexes place a residue carrying a metal within 5\angstrom{} of the deposited reference instance of the ligand
```

**Notes:** D4: the registered rule is residue-level (any atom of a residue containing a metal element) and is measured to the reference instance, not to the scored copy. The 78 is unchanged.

#### W3. M:119 — binding-site recovery sentence

**Current (`body_main_short.tex:119`):**
```latex
It is assessed by clustering poses across tools and against the crystallographic pocket, a resolution also applicable to reference-free Orai1.
```

**Proposed:**
```latex
It is assessed by clustering poses across tools and against the crystallographic pocket, where a cluster reaches the pocket if it lies within the threshold of any deposited copy of the ligand, a resolution also applicable to reference-free Orai1.
```

**Notes:** D5 (any copy counts for cluster reach). The Orai1 clause is unaffected.

#### W4. M:130 — form value names the copy (plan row 124-133)

**Current (`body_main_short.tex:130`):**
```latex
The value is obtained from the PoseBusters implementation \cite{ref035}.
```

**Proposed:**
```latex
The value is obtained from the PoseBusters implementation \cite{ref035} and is taken against the deposited copy that the in-place minimum selects for the pose.
```

**Notes:** D2: one argmin-RMSD copy j* per pose drives centroid, Kabsch, torsions, contacts and PLIF. Nothing else in `:124-133` needs a wording change; the decomposition caveat at `:133` stays valid.

#### W5. M:158 — crystal-site recovery threshold

**Current (`body_main_short.tex:158`):**
```latex
Crystal-site recovery is assessed separately at the prespecified 4\angstrom{} centroid threshold. The cluster nearest the crystal is not automatically counted as recovered.
```

**Proposed:**
```latex
Crystal-site recovery is assessed separately at the prespecified 4\angstrom{} centroid threshold, measured to the nearest deposited copy of the ligand, and a cluster reaches the site when it lies within that distance of any copy. The cluster nearest the crystal is not automatically counted as recovered.
```

**Notes:** D5. Note the plan's own consequence: the cluster report's pooled 4 Å oracle follows D5 (any copy, [302 / 303 = 99.7] %) while the hub's per-pose `oracle_centroid_le_4A_%` follows D2 ([301 / 303 = 99.3] %); `:476` prints the cluster-report value.

#### W6. M:165 — fingerprint reference and the union rule

**Current (`body_main_short.tex:165`):**
```latex
On the calibration benchmark, each tool is compared with the crystal-ligand fingerprint using mean Jaccard overlap and rank-resolved precision, recall and F1 for native-contact recovery.
```

**Proposed:**
```latex
On the calibration benchmark, each pose is compared with the fingerprint of the deposited ligand copy it lies nearest to, and a tool\textquotesingle s top-5 union is compared with the union of the fingerprints of the copies its retained poses are nearest to, using mean Jaccard overlap and rank-resolved precision, recall and F1 for native-contact recovery.
```

**Notes:** Plan Section 1.8: adopt and state the pose-set rule for the top-5 union. Values (Jaccard, F1, Table 4, Fig 7) are recompute after per-copy crystal fingerprints exist.

#### W7. M:170 — cost definition

**Current (`body_main_short.tex:170`):**
```latex
Computational cost is reported per pose that is PoseBusters-valid and within 2\angstrom{} of the crystal ligand, charging each pipeline against successful benchmark poses.
```

**Proposed:**
```latex
Computational cost is reported per pose that is PoseBusters-valid and within 2\angstrom{} of the nearest deposited copy of the crystal ligand, charging each pipeline against successful benchmark poses.
```

**Notes:** Definition only. Qualifying poses become [442 / 2,068 / 483] (plan 1.9).

#### W8a. M:193 — pooled counts and the DiffDock refiner agreement

**Current (`body_main_short.tex:193`):**
```latex
Rescoring cannot extend a pool, so at that variant raw search leads its own rescored sibling, 164 complexes against 160. The full-pool endpoint cannot express the ranking gain that motivates rescoring. The two DiffDock refiners each recover 138 complexes. They agree on 133 and differ on five in each direction (exact McNemar p = 1.000), although gnina has higher pooled physical validity than smina (87.1\% versus 84.5\%). For unguided EquiBind, gnina also exceeds smina in pooled validity (62.0\% versus 51.3\%) and complex count (55 versus 48).
```

**Proposed:**
```latex
Rescoring cannot extend a pool, so at that variant raw search leads its own rescored sibling, [178] complexes against [175]. The full-pool endpoint cannot express the ranking gain that motivates rescoring. The two DiffDock refiners recover [148] and [149] complexes. They agree on [145] and differ on three and four (exact McNemar p = [rebuild]), although gnina has higher pooled physical validity than smina (87.1\% versus 84.5\%). For unguided EquiBind, gnina also exceeds smina in pooled validity (62.0\% versus 51.3\%) and complex count ([57] versus [50]).
```

**Notes:** Plan Section 3.2 row 191-199: 164 vs 160 → 178 vs 175; 138 / 138 agree 133, 5 / 5 → 148 / 149 agree 145, 3 / 4; 55 vs 48 → 57 vs 50. The "each recover" wording must go because the two counts no longer coincide. Validity percentages are unchanged.

#### W8b. M:195 — top-15 selection counts

**Current (`body_main_short.tex:195`):**
```latex
AutoDock Vina at exhaustiveness 128 with gnina rescoring reaches 160 complexes and is considered the leading configuration. Appendix~\ref{exhaustiveness-selection} elaborates on the examined AutoDock configurations. DiffDock with smina reaches 133, compared with 132 for gnina. Unguided EquiBind with gnina reaches 55, compared with 46 for smina. Pocket guidance is worse than no guidance at all for EquiBind, and no fpocket- or P2Rank-guided arm exceeds 18 at any depth.
```

**Proposed:**
```latex
AutoDock Vina at exhaustiveness 128 with gnina rescoring reaches [175] complexes and is considered the leading configuration. Appendix~\ref{exhaustiveness-selection} elaborates on the examined AutoDock configurations. DiffDock with smina reaches [144], compared with [143] for gnina. Unguided EquiBind with gnina reaches [57], compared with [48] for smina. Pocket guidance is worse than no guidance at all for EquiBind, and no fpocket- or P2Rank-guided arm exceeds [rebuild] at any depth.
```

**Notes:** Also on the same line: "Under that endpoint the same AutoDock variant stands at 198 complexes at top-15 rather than 160." → "[213] ... rather than [175]". The guided-arm ceiling "18" is recompute (plan 1.5: P2Rank-gnina sits at 18 today and the probe covers five of nine configurations).

#### W8c. M:197 — smina against gnina, "level" not "gnina higher"

**Current (`body_main_short.tex:197`):**
```latex
DiffDock + smina is carried forward on the selection criterion, leading gnina by 133 to 132 complexes and by 834 to 828 qualifying poses at the top-15 ranking depth. Since the arms disagree on nine complexes, including four recovered only by gnina the ordering is unstable. Moreover, Table~\ref{tab:results-pose-production} gives gnina higher pooled validity and near-native recovery (176 versus 173 complexes), but smina higher form recovery (237 versus 235), i.e. components also conflict.
```

**Proposed:**
```latex
DiffDock + smina is carried forward on the selection criterion, leading gnina by [144] to [143] complexes and by [1,031] to [1,022] qualifying poses at the top-15 ranking depth. Since the arms disagree on seven complexes, including three recovered only by gnina the ordering is unstable. Moreover, Table~\ref{tab:results-pose-production} gives gnina higher pooled validity, leaves the two level on near-native recovery ([186] against [186] complexes) and gives smina higher form recovery ([236] against [235]), i.e. components also conflict.
```

**Notes:** Plan D12 and Section 3.2 row 191-199 ("W: level, not gnina higher"). The footnote on the same line (32.7 / 33.7 / 32.3 %) is recompute; its no-validity basis moves 35.3 → [45.2] % for the confidence ranking.

#### W9. M:277-283 — Table 1 header (optional footnote)

**Current (`body_main_short.tex:278`):**
```latex
RMSD \ensuremath{\leq} 2\angstrom{}\strut}} & \multicolumn{3}{c}{%
```

**Proposed:**
```latex
% No header change. Optional footnote text for the table follows.
% Near-nativeness and form are taken to the deposited copy of the ligand nearest to each pose (Methods, Evaluation Metrics).
```

**Notes:** Plan Section 3.2 row 264, 277-283 says "W optional footnote". The header cells themselves need no change. If a footnote is added it goes in the table's `\multicolumn` note the way Tables 8 and 27 carry theirs.

#### W10. M:324 — validity cost "at most four"

**Current (`body_main_short.tex:324`):**
```latex
The validity requirement removes at most four of 303 complexes from the headline endpoint at any depth.
```

**Proposed:**
```latex
The validity requirement removes at most five of 303 complexes from the headline endpoint at any depth.
```

**Notes:** Plan 1.1: DiffDock + smina validity-gate cost by depth becomes [5 / 4 / 3 / 3], AutoDock and EquiBind unchanged. Same change at M:777 (W22) and A:1196 (W47).

#### W11. M:326 — Table 2 lead-in

**Current (`body_main_short.tex:326`):**
```latex
A complex is recovered if any inspected pose is PoseBusters-valid and within 2\angstrom{}.
```

**Proposed:**
```latex
A complex is recovered if any inspected pose is PoseBusters-valid and within 2\angstrom{} of the nearest deposited copy of the ligand.
```

**Notes:** Definition sentence for Table 2. The table rows (`:347-349`) are value-only.

#### W12. M:380 — "as-placed crystal RMSD"

**Current (`body_main_short.tex:380`):**
```latex
against as-placed crystal RMSD (i.e. near-nativeness) and best-fit Kabsch RMSD (i.e. form)
```

**Proposed:**
```latex
against as-placed RMSD to the nearest deposited ligand copy (i.e. near-nativeness) and best-fit Kabsch RMSD to the same copy (i.e. form)
```

**Notes:** Both axes are taken on the copy the in-place minimum selects (D2).

#### W13. M:393 — Fig 2 caption (Fig 3 caption at M:401 holds)

**Current (`body_main_short.tex:393`):**
```latex
\caption{Best-of-top-N accuracy against the as-placed crystal RMSD}\label{fig:results-r2a}
```

**Proposed:**
```latex
\caption{Best-of-top-N accuracy against the as-placed RMSD to the nearest deposited ligand copy}\label{fig:results-r2a}
```

**Notes:** Caption at `:401` ("Best-of-top-N accuracy against the best-fit (Kabsch) RMSD") needs no wording change; both figures (image3, image4) regenerate in Phase 6.

#### W14. M:396 — form tie wording and discordant counts

**Current (`body_main_short.tex:396`):**
```latex
In contrast, the AutoDock--DiffDock comparison is a tie throughout, at p\_holm = 1.000 with 47 against 46 discordant complexes at rank-1, 30 against 31 at top-15 and 25 against 26 at top-30.
```

**Proposed:**
```latex
In contrast, the AutoDock--DiffDock comparison is a tie throughout, not significant at any depth, with [48] against [47] discordant complexes at rank-1, [30] against [32] at top-15 and [25] against [26] at top-30.
```

**Notes:** Plan 1.4: the top-15 Holm p becomes [0.899] within its three-pair family, so the literal "p_holm = 1.000" no longer holds at every depth. Alternative wording the plan allows: "at p_holm = 1.000, 0.899 and 1.000". The EquiBind clause on the same line carries the A2 correction (≤ 9e-9) today and is recompute after the switch. The 74.3 / 74.6, 95.4 / 97.0 and 52.5 / 52.1 / 33.3 percentages are unchanged at printed precision.

#### W15. M:408 — Fig 4 footnote

**Current (`body_main_short.tex:408`):**
```latex
\footnote{Poses whose docked centroid lies more than 8\angstrom{} from the crystal-ligand site were removed leaving 6,169 PoseBusters-valid poses for scoring. Axes are clipped at 5\angstrom{} keeping dense near-native core legible, yet removing 2,329 poses from the graph frame.}
```

**Proposed:**
```latex
\footnote{Poses whose docked centroid lies more than 8\angstrom{} from the nearest deposited copy of the ligand were removed leaving [7,565] PoseBusters-valid poses for scoring. Axes are clipped at 5\angstrom{} keeping dense near-native core legible, yet removing [rebuild] poses from the graph frame.}
```

**Notes:** Plan 1.6: near-site pool 6,169 → [7,565]; the clipped count is taken from the regenerated `20d_…caption.txt` (the printed 2,329 is that script's own rule on the `pb_rmsd` axis and is NOT a defect today).

#### W16. M:422 — near-site interpretive paragraph

*Status:* wording pending rebuild

**Current (`body_main_short.tex:422`):**
```latex
In the near-site cohort, deeper ranks add valid but misplaced AutoDock and EquiBind poses and few for DiffDock. In-place RMSD deteriorates with depth. DiffDock appears more form-limited only because its low in-place RMSD inflates the form-to-in-place ratio. Its absolute form RMSD is lowest.
```

**Proposed:** no replacement text yet (see notes).

**Notes:** Plan Section 3.2 row 422: W, to be re-read against the regenerated 20d sidecar. No replacement is proposed before the rebuild because every claim in the paragraph ("few for DiffDock", "form-limited only because") is a reading of Table 27, all 117 cells of which change.

#### W17. M:481 — interaction section opening (plan row 481-532)

**Current (`body_main_short.tex:481`):**
```latex
Crystal-pose interactions are compared with those of PoseBusters-valid pipeline poses.
```

**Proposed:**
```latex
Crystal-pose interactions are compared with those of PoseBusters-valid pipeline poses, each pose against the deposited ligand copy it lies nearest to and each top-5 union against the union of the fingerprints of those copies.
```

**Notes:** Plan 1.8. All values at `:489`, `:493`, `:514-525` and `:530` (0.450 / 0.417 / 0.380, F1 0.57 / 0.53 / 0.47, Table 4, 21.0 / 19.4 / 15.3 / 16.4) are recompute; preview F1 [0.532 / 0.527 / 0.513].

#### W18. M:697 — Orai1 usable-yield paragraph, calibration clause

**Current (`body_main_short.tex:697`):**
```latex
Calibration point estimates show the same ordering, although not a statistically separable difference.
```

**Proposed:**
```latex
Calibration point estimates at rank-1 show the same ordering, although not a statistically separable difference.
```

**Notes:** Plan Phase 5.3: true at rank-1 only under both conventions (from top-5 onward AutoDock leads and the contrast is resolved). This is the one Section 4.2 sentence that changes (plan Section 7).

#### W19. M:708 — cost definition in the Results

**Current (`body_main_short.tex:708`):**
```latex
Cost is normalised per qualifying pose, defined as PoseBusters-valid and within 2\angstrom{} of the crystal ligand.
```

**Proposed:**
```latex
Cost is normalised per qualifying pose, defined as PoseBusters-valid and within 2\angstrom{} of the nearest deposited copy of the crystal ligand.
```

**Notes:** Definition only.

#### W20a. M:744 — median yield per succeeding complex and the ordering sentence

**Current (`body_main_short.tex:744`):**
```latex
The median successful complex yields one, although 199 of 303 complexes now yield at least one, far more than the 169 of DiffDock.
```

**Proposed:**
```latex
The median successful complex yields two, although [214] of 303 complexes now yield at least one, far more than the [183] of DiffDock.
```

**Notes:** Plan 1.9: AutoDock's median succeeding complex yields two qualifying poses instead of one. The opening sentence of the same line, "The configured AutoDock pipeline is the most expensive route to a qualifying pose under blind whole-protein search.", is W-if: it changes only if the re-derived medians flip the ordering. Friedman 85.2 / 3.2e-19 / 0.80, the 53 shared complexes (→ [54]) and the 0.14 / 0.02 ratios are recompute.

#### W20b. M:758 — "adds nine qualifying poses"

**Current (`body_main_short.tex:758`):**
```latex
because rescoring adds nine qualifying poses at no CPU cost.
```

**Proposed:**
```latex
because rescoring adds eleven qualifying poses at no CPU cost.
```

**Notes:** Plan 1.1: 442 − 431 = 11 (preview). The 1,805.0 vs 1,752.8 CPU-seconds on the same line are recompute. `:760` (271 / 479 / 38 s) recompute.

#### W21. M:768 — principal findings interval and reach clause

**Current (`body_main_short.tex:768`):**
```latex
Under the blind whole-protein search, AutoDock had the highest physical validity and led DiffDock on near-native recovery and native-pocket reach at every depth. The near-native lead is inconclusive at rank-1 and separable from top-5 onwards, and reach never separated them. The rank-1 interval spans a 4.2-point DiffDock lead to a 9.4-point AutoDock lead.
```

**Proposed:**
```latex
Under the blind whole-protein search, AutoDock had the highest physical validity and led DiffDock on near-native recovery at every depth [and on native-pocket reach, value from the Phase 3 rebuild under D5]. The near-native lead is inconclusive at rank-1 and separable from top-5 onwards[, and reach never separated them, value from the Phase 3 rebuild]. The rank-1 interval spans a [1.7]-point DiffDock lead to a [12.8]-point AutoDock lead.
```

**Notes:** Plan Section 3.2 row 768. The interval swap is a preview (Newcombe [−1.7, +12.8]); the two reach clauses are recompute under D5 and are marked in brackets rather than rewritten.

#### W22. M:777 — validity cost in the Discussion

**Current (`body_main_short.tex:777`):**
```latex
Requiring validity removes at most four of 303 complexes at any depth.
```

**Proposed:**
```latex
Requiring validity removes at most five of 303 complexes at any depth.
```

**Notes:** As W10.

#### W23. M:781 — refiner selection paragraph

**Current (`body_main_short.tex:781`):**
```latex
smina is carried forward on the selection criterion by 133 complexes to 132, but that one-complex margin survives only in the three-way conjunction of validity, near-nativeness and form. Dropping either the validity or the form term leaves the two arms exactly level, and gnina is level or ahead on every component taken alone.
```

**Proposed:**
```latex
smina is carried forward on the selection criterion by [144] complexes to [143] at the top-15 selection depth, and that one-complex margin does not depend on the three-way conjunction of validity, near-nativeness and form. Dropping either the validity or the form term leaves smina ahead by two, the two arms are level on pooled near-native recovery ([186] against [186]) and smina stays ahead on pooled form ([236] against [235]).
```

**Notes:** Plan 1.5 gives the replacement wording ("smina ahead by two on both drops; level on pooled near-native recovery 186 against 186; smina ahead on pooled form 236 against 235") and requires "at the top-15 selection depth" because `:781` is at the selection depth whereas `:197` is Table 1's pooled basis. Same line, later sentences: "gnina recovers 83 near-native complexes at top-15 against 61 for smina" → "[86] ... against [64]"; "where the margin is nine complexes" → "seven complexes" (plan Section 3.2 row 781). "Exact McNemar separates them at no depth" is to be re-checked from the rebuild.

#### W24a. M:788 — ranking-depth paragraph

**Current (`body_main_short.tex:788`):**
```latex
Moving from rank-1 to the best of fifteen substantially increased recovery for DiffDock and AutoDock but added only about eight percentage points for EquiBind (Results, Near-Nativeness and Form). Useful DiffDock and AutoDock poses were therefore often under-ranked. Selection was not the main loss. About six in ten rank-1 failures of the two tools contained no qualifying pose at any depth.
```

**Proposed:**
```latex
Moving from rank-1 to the best of fifteen substantially increased recovery for DiffDock and AutoDock but added only about nine percentage points for EquiBind (Results, Near-Nativeness and Form). Useful DiffDock and AutoDock poses were therefore often under-ranked. Selection was not the main loss. About six to seven in ten rank-1 failures of the two tools contained no qualifying pose at any depth.
```

**Notes:** Plan Section 3.2 row 788, 857: "eight pp" → +8.6 (which rounds to nine, so "about eight" would misstate a preview of +8.6; the author may prefer "between eight and nine"), "six in ten" → "six to seven in ten" ([57.8] % and [70.2] %). "nine in ten" for EquiBind holds ([89.4] %).

#### W24b. M:857 — RQ4 counts

**Current (`body_main_short.tex:857`):**
```latex
About six in ten AutoDock and DiffDock rank-1 failures, and nine in ten EquiBind failures, contained no qualifying pose. Selection still recovered 87 and 64 additional complexes at top-15 for the leading tools, compared with 25 for EquiBind.
```

**Proposed:**
```latex
About six to seven in ten AutoDock and DiffDock rank-1 failures, and nine in ten EquiBind failures, contained no qualifying pose. Selection still recovered [64] and [47] additional complexes at top-15 for the leading tools, compared with [26] for EquiBind.
```

**Notes:** Plan 1.1: rank-1 → top-15 gains +87 / +64 / +25 → [+64 / +47 / +26].

#### W25. M:834 — RQ1 numbers and the compatibility interval

**Current (`body_main_short.tex:834`):**
```latex
Rank-1 recovery was 36.6\% for AutoDock, 34.0\% for DiffDock and 18.2\% for EquiBind. Retrospective top-15 recovery was 65.3\%, 55.1\% and 26.4\%. The AutoDock--DiffDock rank-1 difference was +2.6 percentage points (95\% interval \ensuremath{-}4.2 to +9.4) and was not statistically resolved, although AutoDock led from top-5 onwards, by +13.2 points there (p\_holm = 7.6e-4) and +10.2 at top-15.
```

**Proposed:**
```latex
Rank-1 recovery was [49.2]\% for AutoDock, [43.6]\% for DiffDock and [18.8]\% for EquiBind. Retrospective top-15 recovery was [70.3]\%, [59.1]\% and [27.4]\%. The AutoDock--DiffDock rank-1 difference was [+5.6] percentage points (95\% interval \ensuremath{-}[1.7] to +[12.8]) and was not statistically resolved, although AutoDock led from top-5 onwards, by +13.2 points there (p\_holm = 7.6e-4) and [+11.2] at top-15.
```

**Notes:** Same line, later: "A post-hoc 90\% compatibility interval places the rank-1 difference within about eight points." → "within about twelve points." (plan 1.2: 90 % interval [−0.5, +11.7], equivalence within 12 not 10). The top-5 contrast is identical under both conventions.

#### W26. M:843 — RQ2, AutoDock re-ranking step and its Holm family (D14)

**Current (`body_main_short.tex:843`):**
```latex
AutoDock gained by re-ranking, from 31.0\% to 36.6\% at rank-1. That single step does not clear Holm correction on its own, and the gain reaches significance only from top-5 onwards.
```

**Proposed:**
```latex
AutoDock gained by re-ranking, from [42.9]\% to [49.2]\% at rank-1. That single step does not clear Holm correction within the sixteen-test within-tool family of Table~\ref{tab:appendix-refinement-mcnemar}, although it does within the seven-contrast ladder family of Appendix~\ref{exhaustiveness-selection}, and the gain is resolved at every depth from top-5 onwards.
```

**Notes:** D14: the two families must be named, otherwise `:843` ("does not clear Holm") and App. C `:182` ("six of seven survive") contradict each other under nearest. Previews: PB-valid ∧ ≤ 2 Å basis 42.9 → 49.2 %, 20 / 39 discordant, p 0.0183, Holm over App. C's seven contrasts 0.037 (survives); Table 19's sixteen-test family gives Holm 0.060 (does not). The rate bands earlier on the line (0.1-0.7, 4.1-4.6, 3.5-4.4) are recompute (DiffDock + smina [0.4-0.8] derived).

#### W27. M:848-850 — RQ3 (W only if the cost ordering flips)

*Status:* W-if

**Current (`body_main_short.tex:850`):**
```latex
Configured AutoDock was most expensive per qualifying pose in charged, CPU-core and GPU seconds.
```

**Proposed:** no replacement text yet (see notes).

**Notes:** Plan Section 3.2 row 848-850: medians 2.4 / 49.4 / 137.2 s recompute; "DiffDock succeeded on 169 complexes, AutoDock on 199 and EquiBind on 80" (`:848`) → [183 / 214 / 83]; "twentyfold" recompute; "quarter" holds. The quoted ordering sentence, and its sibling at `:744`, `:758`, `:817` and in the Abstract and Kurzfassung, change wording only if the re-derived per-qualifying-pose medians flip the AutoDock-most-expensive ordering. No text is proposed until plan Section 1.9 is re-derived from the timing files.

#### W28. M:866 — winning margins

**Current (`body_main_short.tex:866`):**
```latex
the winning margins were eight complexes for AutoDock Vina, one for DiffDock and nine for EquiBind
```

**Proposed:**
```latex
the winning margins were nine complexes for AutoDock Vina, one for DiffDock and nine for EquiBind
```

**Notes:** Plan 1.5: AutoDock 175 vs raw 166 = nine. "3.0 percentage points" holds; "between 0.4 and 0.8 percentage points" is recompute.

#### W30. A:165 — metal stratum paragraph (D4)

**Current (`body_appendix_short.tex:165`):**
```latex
In 78 of the 303 analysed entries a metal ion or a metal-containing cofactor lies within 5\angstrom{} of a crystal-ligand heavy atom, most often magnesium, zinc or iron.
```

**Proposed:**
```latex
In 78 of the 303 analysed entries a residue carrying a metal ion or a metal-containing cofactor has an atom within 5\angstrom{} of a heavy atom of the deposited reference instance of the ligand, most often magnesium, zinc or iron. Membership is a property of the complex and does not follow the copy a pose is scored against, so one AutoDock rank-1 success (7JHQ\_VAJ, scored against a copy with no metal within 20\angstrom{}) counts as metal-adjacent and one DiffDock rank-1 and top-15 success (7TB0\_UD1, scored against the one copy that has a metal within 5\angstrom{}) counts as metal-free.
```

**Notes:** D4 crossover sentence and residue-rule wording. Value swaps on the same line, all previews from plan Phase 5.4 and Section 3.3 row 165, *value from the Phase 3 rebuild*: "Rank-1 recovery is 37.8\% for AutoDock and 39.6\% for DiffDock on the 225 metal-free complexes, against 33.3\% and 17.9\% on the 78 metal-adjacent ones." → [51.1] / [49.3] and [43.6] / [26.9]; "\ensuremath{-}6.7 percentage points on the metal-free stratum (134 against 149 of 225, exact McNemar p = 0.142) and \ensuremath{-}20.5 points on the metal-adjacent stratum (33 against 49 of 78, p = 0.020)" → [−8.4] ([141] against [160], p [0.056]) and [−19.2] ([38] against [53], p [0.028]); "Fisher\textquotesingle s exact test on the discordant pairs p = 0.256" → [0.440]; failure shares "53.6\% of AutoDock and 65.4\% of DiffDock" → [58.2] / [70.2] and "55.8\% and 70.3\%" → [56.8] / [70.2]; "\ensuremath{-}7.4 and \ensuremath{-}19.2 points and of \ensuremath{-}6.3 and \ensuremath{-}21.0 points" → [−9.1 / −17.8] and [−8.1 / −19.8]. The 73 / 78 / 81 stratum sizes hold (registered rule). The cofactor split is A5. The 7JHQ_VAJ alternate distances (21.8-39.5 Å) and the 7TB0_UD1 copy-1 distance (K at 3.17 Å) are from D4.

#### W31. A:167 — box construction

**Current (`body_appendix_short.tex:167`):**
```latex
No further padding is applied.
```

**Proposed:**
```latex
No further padding is applied, and every deposited copy of the ligand lies inside the box.
```

**Notes:** Plan Phase 5.4 ("add that every deposited copy lies inside this box"). Holds because the box is the receptor's bounding box plus 1 Å and every copy is bound to the receptor.

#### W32. A:182 — exhaustiveness selection paragraph

**Current (`body_appendix_short.tex:182`):**
```latex
At top-15 the selected arm recovers 198 of 303 complexes and beats each of the seven alternatives on an exact McNemar test, with Holm-adjusted p-values running from 2e-21 against exhaustiveness 18 to 0.007 against its own raw search. It leads at every depth from top-3 to top-22. At rank-1 only two of the seven contrasts survive correction and the arm ties with exhaustiveness 64 rescored, at 111 against 110. Over the full thirty-pose pool the raw exhaustiveness-128 arm is nominally ahead, 202 against 199, because rescoring minimises poses and a few cross the threshold the wrong way. The selected arm is therefore the best configuration at the depth a user inspects rather than at every depth. On the three-part endpoint the Results use for selection the same ordering holds, at 160 complexes against 152 for its own raw search and 95 at exhaustiveness 18, with Holm-adjusted p-values from 6e-18 to 0.039 over the same family. That family is the seven contrasts against the selected arm, and the step over raw exhaustiveness 128 is the one contrast that would not survive a correction taken over all twenty-eight pairs of the ladder.
```

**Proposed:**
```latex
At top-15 the selected arm recovers [213] of 303 complexes and beats each of the seven alternatives on an exact McNemar test, with Holm-adjusted p-values running from [7e-20] against exhaustiveness 18 to [0.005] against its own raw search and against exhaustiveness 64 rescored. It leads strictly at every depth from rank-1 to top-25 and ties at top-26. At rank-1 six of the seven contrasts survive correction within this seven-contrast family. The one that does not is the contrast with exhaustiveness 64 rescored, where the arm leads nominally at [149] against [141] on [16] against [8] discordant complexes (exact p = [0.15]) and the ordering is unresolved. Over the full thirty-pose pool the raw exhaustiveness-128 arm is nominally ahead, [217] against [214], because rescoring minimises poses and a few cross the threshold the wrong way. The selected arm is therefore the best configuration at every depth a user is likely to inspect, and the raw search overtakes it only over the whole pool. On the three-part endpoint the Results use for selection the same ordering holds, at [175] complexes against [166] for its own raw search and [rebuild] at exhaustiveness 18, with Holm-adjusted p-values from [rebuild] over the same family. That family is the seven contrasts against the selected arm, and [the steps over raw exhaustiveness 128 and over exhaustiveness 64 rescored are the contrasts that would not survive a correction taken over all twenty-eight pairs of the ladder, value from the Phase 3 rebuild].
```

**Notes:** Plan 1.5 and Section 3.3 row 182. Three wording changes: the depth lead ("top-3 to top-22" → strictly 1-25, ties at 26, so the claim is strengthened), "only two of the seven contrasts survive" → "six of seven" (only 64 + gnina fails, p [0.152]) and the rank-1 tie → nominal lead 149 against 141 that stays unresolved (16 / 8 discordant, exact p 0.15). The twenty-eight-pair sentence is recompute (uncorrected 0.0026 × 28 = 0.073 for the raw-128 contrast and 0.0044 × 28 = 0.12 for the 64 + gnina contrast, so both likely fail, but the rebuild decides). The D14 family statement ("within this seven-contrast family") is what lets `:843` (W26) quote the other family without contradiction. Holm range preview 7.3e-20 to 5.2e-3 (plan 1.5).

#### W33. A:195 — Fig 19 caption

**Current (`body_appendix_short.tex:195`):**
```latex
and the top step is worth three complexes at rank-1 but does not become individually resolvable until about depth 22}
```

**Proposed:**
```latex
and the top step is worth [seven] complexes at rank-1, where it is already individually resolvable at an unadjusted p of [0.039]}
```

**Notes:** Plan 1.5: exh92 → exh128 becomes +7 at rank-1 with uncorrected p 0.039 already at rank-1. The `:235` caption ("Exhaustiveness 92 ... buys nothing at rank-1") is re-checked from the rebuild (91 vs 94 today, [123] vs [123] under nearest, so the clause likely holds).

#### W34. A:227 — Table 8 footnote gate

**Current (`body_appendix_short.tex:227`):**
```latex
Complexes of 303 with at least one pose that is PoseBusters-valid and within 2\angstrom{} among the first d ranks.
```

**Proposed:**
```latex
Complexes of 303 with at least one pose that is PoseBusters-valid and within 2\angstrom{} of the nearest deposited copy of the ligand among the first d ranks.
```

**Notes:** Table 8 rows (`:218-225`) are value-only (plan 1.5 ladder table). The rest of the footnote (search hours, thread allocation) is unchanged.

#### W35. A:230 — D10 confirmatory cell kept on the single instance, nearest values added

**Current (`body_appendix_short.tex:230`):**
```latex
In absolute terms the five rungs recover 98, 128, 149, 159 and 164 complexes, so the first step adds thirty and the last adds five.
```

**Proposed:**
```latex
In absolute terms the five rungs recover 98, 128, 149, 159 and 164 complexes, so the first step adds thirty and the last adds five. Those counts and the confirmatory interval above are kept on the single deposited instance on which the cell was fixed in advance. Scored against the nearest deposited copy the same rungs recover [115, 143, 166, 172 and 178] complexes, so the first step adds [28] and the last [6], and the per-hour returns are [rebuild].
```

**Notes:** D10 as adopted: the whole paragraph stays on the instance convention and one sentence gives the nearest values, because Figures 20-21 and `:182` ("160 against 152" → [175 against 166]) are on the same pool. Previews from plan Section 3.3 row 230.

#### W36. A:246 — plateau argument (rewritten, not renumbered)

**Current (`body_appendix_short.tex:246`):**
```latex
The ceiling that remains is also a ranking ceiling rather than a sampling one. Measured on placement alone, without the validity requirement the ladder table imposes, raw rank-1 recovery tops out at 95 complexes anywhere above exhaustiveness 64, while rescoring lifts it to 113 and doubling the pool beneath the rescorer leaves it at exactly 113. Extra sampling keeps finding poses that the empirical function then fails to promote.
```

**Proposed:**
```latex
The ceiling that remains is only partly a ranking ceiling. Measured on placement alone, without the validity requirement the ladder table imposes, raw rank-1 recovery is [124] complexes at exhaustiveness 64 and 92 and [131] at 128, while rescoring lifts the two rescored rungs to [144] and [151], so doubling the pool beneath the rescorer still adds [seven] complexes at rank-1. Extra sampling keeps finding poses that the empirical function fails to promote, and the rescorer promotes a part of what the harder search finds.
```

**Notes:** Plan 1.5: the plateau argument COLLAPSES under nearest (95 / 95 / 95 becomes [124 / 124 / 131], "exactly 113" becomes [144] vs [151]), so the sentence is rewritten. The opening "Two things follow, and a third does not." must then read "One thing follows, one follows in part, and a third does not." or the author's equivalent, because the second claimed consequence is no longer clean. The surrounding numbers at `:244` (five, three, 18.7, −0.110 h) are recompute except the wall-time delta.

#### W37. A:267 — pocket-guided EquiBind

**Current (`body_appendix_short.tex:267`):**
```latex
both raw guided arms recovering no near-native complex at all
```

**Proposed:**
```latex
the fpocket raw arm recovering no near-native complex and the P2Rank raw arm one
```

**Notes:** Plan 1.5 and Section 3.3 row 267-269: P2Rank raw recovers one complex over the pool (2 poses) under nearest. All counts on the line (2 / 46 / 55, 0 / 7 / 14, 0 / 5 / 18, 37, 19 / 64 / 83, 0 / 14 / 24, 0 / 18 / 35, "between 180 and 200") are recompute; previews 2 / [48] / [57], 19 / [67] / [86], form [179-199] on five of nine configurations. Validity percentages hold. `:269` (261 / 282 / 266) unchanged.

#### W38. A:276 — accuracy scoring sentence

**Current (`body_appendix_short.tex:276`):**
```latex
and geometric accuracy was scored as the symmetry-corrected heavy-atom RMSD to the crystal ligand.
```

**Proposed:**
```latex
and geometric accuracy was scored as the symmetry-corrected heavy-atom RMSD to the nearest deposited copy of the crystal ligand.
```

**Notes:** Definition only.

#### W39. A:487 — optional clause on the reference-based RMSD check

**Current (`body_appendix_short.tex:487`):**
```latex
It is kept off the validity axis here because near-nativeness is gated separately.
```

**Proposed:**
```latex
It is kept off the validity axis here because near-nativeness is gated separately. That test loads every deposited instance of the ligand and reports the lowest RMSD, which is the copy rule the Methods adopt.
```

**Notes:** Plan Section 3.3 row 487: optional clause on `load_all` (installed `redock.yml:301-307`, `modules/rmsd.py:52-70`). The M:140 footnote on the post-hoc reference run is recompute and must state which file was passed as `mol_true`.

#### W40a. A:543 — profiled pool and standardised F1 (plan anchor :545-553 drifted here)

**Current (`body_appendix_short.tex:543`):**
```latex
being 4,253 poses over 303 complexes with a median deviation of 7.57\angstrom{} from the crystal ligand and 21.3\% of poses within 2\angstrom{}.
```

**Proposed:**
```latex
being 4,253 poses over 303 complexes with a median deviation of [rebuild]\angstrom{} from the nearest deposited copy of the crystal ligand and [rebuild]\% of poses within 2\angstrom{}.
```

**Notes:** Anchor note: plan Section 3.3 row 545-553 lists 4,253, 7.57 Å, 21.3 %, F1 0.565 / 0.529 / 0.471, 750 and 0.0003; in the live file the first four sit on `:543` and the 750 / 0.0003 on `:549`. The 4,253 profiled poses do not change (PB-valid top-5 is convention-free); the median deviation, the 2 Å share, the F1 triple "0.565, 0.529 and 0.471 to 0.532, 0.527 and 0.513" and the 10 Å-restriction figures are recompute against the nearest copy.

#### W40b. A:551 — canonical receptor and per-copy crystal fingerprint

**Current (`body_appendix_short.tex:551`):**
```latex
All poses of the calibration benchmark are therefore profiled against a canonical receptor set that carries author numbering and was verified heavy-atom coordinate-compatible with each tool\textquotesingle s own docking receptor, which relabels residues without moving any atom.
```

**Proposed:**
```latex
All poses of the calibration benchmark are therefore profiled against a canonical receptor set that carries author numbering and was verified heavy-atom coordinate-compatible with each tool\textquotesingle s own docking receptor, which relabels residues without moving any atom. The crystal fingerprint is computed for every deposited copy of the ligand on that receptor, and each pose is scored against the copy it lies nearest to, because a pose placed at an alternate copy would otherwise be scored against the residues of a site it does not occupy.
```

**Notes:** Plan 1.8: 105 / 119 / 11 top-5 PB-valid poses of AutoDock* / DiffDock* / EquiBind* sit within 2 Å of an alternate copy and today score Jaccard ≈ 0 because the fingerprint keys carry the reference chain.

#### W43. A:722 — name the multi-copy share beside Figure fig:appendix-dataset-3

**Current (`body_appendix_short.tex:722`):**
```latex
Roughly a quarter of the ligands contain phosphorus (25\%) or sulfur (24\%), 22\% carry at least one halogen, and metal atoms are essentially absent from the ligand records.
```

**Proposed:**
```latex
Roughly a quarter of the ligands contain phosphorus (25\%) or sulfur (24\%), 22\% carry at least one halogen, and metal atoms are essentially absent from the ligand records. In 143 of the 308 entries, 46\%, the deposited structure holds more than one copy of the ligand of interest, and 138 of those are among the 303 analysed complexes.
```

**Notes:** Plan Section 3.3 row 722 ("name the 46 % (143 of 308) share"). The figure the plan calls Fig 30 is `fig:appendix-dataset-3` (`:727`). 143 / 308 = 46.4 %.

#### W44. A:784 — start-conformer difficulty (D13)

**Current (`body_appendix_short.tex:784`):**
```latex
The best-fit (Kabsch) symmetry-corrected RMSD between this start conformer and the crystallographic pose is computed after optimal superposition, so it isolates conformational difference from placement.
```

**Proposed:**
```latex
The best-fit (Kabsch) symmetry-corrected RMSD between this start conformer and the deposited reference instance of the crystallographic pose is computed after optimal superposition, so it isolates conformational difference from placement. It is a property of the start conformer and is not re-taken against the nearest copy.
```

**Notes:** D13: kept on the single instance; 1.56 Å and 65 % unchanged.

#### W46. A:1181 — Table 24 footnote near-native definition

**Current (`body_appendix_short.tex:1181`):**
```latex
Near-native poses are those within 2\angstrom{} of the crystal ligand by symmetry-corrected heavy-atom RMSD without superposition
```

**Proposed:**
```latex
Near-native poses are those within 2\angstrom{} of the nearest deposited copy of the crystal ligand by symmetry-corrected heavy-atom RMSD without superposition
```

**Notes:** Rows `:1172-1179` near-native half: [437 / 98.6; 451 / 98.0; 2,118; 2,169 / 95.3; 2,171; 515 / 93.8], per-check columns recompute. All-poses half unchanged.

#### W47a. A:1188 — validity-objective coupling, near-native strata

**Current (`body_appendix_short.tex:1188`):**
```latex
Intermolecular failure for DiffDock falls from 14.6\% over all produced poses to 5.4\% over its 1,715 near-native poses, and for EquiBind from 37.0\% to 6.2\% over 498 near-native poses.
```

**Proposed:**
```latex
Intermolecular failure for DiffDock falls from 14.6\% over all produced poses to [rebuild]\% over its [2,169] near-native poses, and for EquiBind from 37.0\% to [rebuild]\% over [515] near-native poses.
```

**Notes:** Same line: "EquiBind still fails the intermolecular group on 31 of its 498 near-native poses against 5 of 319 for AutoDock." → "[rebuild] of its [515] ... against [rebuild] of [451]".

#### W47b. A:1196 — validity cost of the headline endpoint

**Current (`body_appendix_short.tex:1196`):**
```latex
On the three variants carried through the Results it removes at most four complexes of 303 at any reported depth. The counts are two at rank-1, four at top-5 and three at the deeper pools for AutoDock Vina + gnina, three to four for DiffDock + smina and two to three for unguided EquiBind + gnina. That is a span of 0.66 to 1.32 percentage points.
```

**Proposed:**
```latex
On the three variants carried through the Results it removes at most five complexes of 303 at any reported depth. The counts are two at rank-1, four at top-5 and three at the deeper pools for AutoDock Vina + gnina, [five] at rank-1, [four] at top-5 and [three] at the deeper pools for DiffDock + smina and two to three for unguided EquiBind + gnina. That is a span of 0.66 to [1.65] percentage points.
```

**Notes:** Plan 1.1: DiffDock* 4 / 4 / 3 / 4 → [5 / 4 / 3 / 3]; span 0.66-1.32 → [0.66-1.65]. The A1b insertion follows this passage in the same paragraph.

#### W48. A:1203 — 25 Å cube regime clause

**Current (`body_appendix_short.tex:1203`):**
```latex
The benchmark this study calibrates against gave Vina a 25\angstrom{} cube centred on the crystal-ligand heavy atoms \cite{ref035}, a volume of 15,625\angstrom{}\textsuperscript{3}.
```

**Proposed:**
```latex
The benchmark this study calibrates against gave Vina a 25\angstrom{} cube centred on the crystal-ligand heavy atoms \cite{ref035}, a volume of 15,625\angstrom{}\textsuperscript{3}. That cube kept almost every other deposited copy of the ligand outside Vina\textquotesingle s search, since only 8 of the 211 alternate copies have a centroid inside it, whereas every box used here contains every copy, so the nearest-copy convention is live for all three tools in this study.
```

**Notes:** D11: "kept almost every other copy outside" (8 of 211 alternate centroids, 16 ids by any heavy atom, lie inside the cube), softened per plan Phase 5.4 ("regime clause, softened per D11"). The inertness of the PUBLISHED Vina anchor is inferred from the thesis's own boxed arm (0 flips), not measured on the paper's unreleased poses.

#### W49. A:1260 — two-pipeline contrast, equivalence margin

**Current (`body_appendix_short.tex:1260`):**
```latex
On that grid only rank-1 is equivalent, within 10 points, and its smallest passing value of 8.4 points is read off the observed interval rather than taken from the grid.
```

**Proposed:**
```latex
On that grid only rank-1 is equivalent, within 12 points and no longer within 10, and its smallest passing value of [11.7] points is read off the observed interval rather than taken from the grid.
```

**Notes:** Plan 1.2 and Section 3.3 row 1258-1260. Value swaps on `:1258` and `:1260`, all previews: +2.6 [−4.2, +9.4] → [+5.6 [−1.7, +12.8]]; +10.2 [+2.8, +17.5] → [+11.2 [+3.9, +18.4]]; +9.9 [+2.5, +17.1] → [+10.2 [+3.0, +17.3]]; "both pipelines recover 51 ... neither recovers 140. Only the remaining 112 discordant complexes, 37.0\%" → [76 / 98 / 129 (42.6 %)]; MDD 9.8 → [10.5]; power 0.12 → [0.32]; "roughly 4,200" → ["roughly 1,100"] (harness rounds to the nearest hundred); 8.2 → [8.8]; "19.4, 16.4 and 16.0" → [19.3, 17.2, 16.2]; 90 % interval "−3.1 to +8.3" → [−0.5 to +11.7]; p 0.009 / 0.011 → [0.0036 / 0.0075] (exact McNemar, uncorrected, as the line prints them). Top-5 identical.

#### W50. A:1265 — threshold-resolved recovery lead-in

**Current (`body_appendix_short.tex:1265`):**
```latex
for both the as-placed and the best-fit (Kabsch) criterion.
```

**Proposed:**
```latex
for both the as-placed and the best-fit (Kabsch) criterion, each taken to the deposited copy of the ligand nearest to the pose.
```

**Notes:** Table 26 RMSD rows are derived and 47 of 90 Kabsch cells move by 0.3-0.7 pp (plan 1.4).

#### W51. A:1387 — Table 27 footnote

**Current (`body_appendix_short.tex:1387`):**
```latex
trimmed of poses whose centroid lies more than 8\angstrom{} from the crystal-ligand site, that is outside the crystal-site neighbourhood rather than off the receptor
```

**Proposed:**
```latex
trimmed of poses whose centroid lies more than 8\angstrom{} from the nearest deposited copy of the ligand, that is outside the crystal-site neighbourhood rather than off the receptor
```

**Notes:** Same footnote, seven numbers (plan 1.6), all previews: "1.490 against 1.543\angstrom{} for AutoDock Vina + gnina and 1.467 against 1.566\angstrom{} for DiffDock + smina" → [1.376 against 1.491] and [1.532 against 1.585]; "two of the 909 tool-by-complex decisions" holds at two; "Rank-1 recovery would then read 37.0\% for AutoDock Vina + gnina and 18.5\% for EquiBind + gnina against an unchanged 34.0\% for DiffDock + smina." → [49.5] / [19.1] / [43.6]. All 117 table cells (`:1377-1385`) change.

#### W52. A:1399 — Fig 38 caption

**Current (`body_appendix_short.tex:1399`):**
```latex
(n = 303, 8\angstrom{} centroid-clustering cut and 4\angstrom{} crystal-site recovery threshold)
```

**Proposed:**
```latex
(n = 303, 8\angstrom{} centroid-clustering cut and 4\angstrom{} crystal-site recovery threshold to any deposited copy of the ligand)
```

**Notes:** D5. The statistics at `:1394` are recompute AFTER the Phase 0.3 loader fix in `pose_cluster_crystal_pocket_report.py:133` (hydrogens counted in pose centroids), which moves every cluster-geometry number even today.

#### W53. A:1412 — Fig 39 caption

**Current (`body_appendix_short.tex:1412`):**
```latex
\caption[Typed contacts versus the crystal by pose rank]{Matched, missed and spurious typed contacts versus the crystal, by pose rank, whole-protein search}
```

**Proposed:**
```latex
\caption[Typed contacts versus the crystal by pose rank]{Matched, missed and spurious typed contacts versus the nearest deposited copy of the crystal ligand, by pose rank, whole-protein search}
```

**Notes:** Kendall τ values at `:1407` and the 15.3 vs 14.1 counts at `:1415` are recompute; image11 regenerates.

#### W54. T:98-99 — Kurzfassung cost clause (W only if the ordering flips)

*Status:* W-if

**Current (`Thesis_short.tex:98-99`):**
```latex
fünfzehn Posen statistisch gesichert. In einheitlicher Belegungsrechnung war der
AutoDock-Ablauf am teuersten pro qualifizierender Pose, EquiBinds niedrige Kosten galten
```

**Proposed:** no replacement text yet (see notes).

**Notes:** See W27. The German clause "war der AutoDock-Ablauf am teuersten pro qualifizierender Pose" changes only if plan Section 1.9 flips the per-qualifying-pose ordering, and any change must be length-neutral on the full page (B4d). T:61 is B4a.


### B6. Value-only rows of Section 3 (no wording change, numbers from the Phase 3 rebuild and the yaml)

Listed for completeness so the author can see that nothing in Section 3 was skipped. None is quoted because the plan prescribes no wording for them; every number comes from the rebuild (D16).

- `body_main_short.tex`: 140 (footnote −1 / −2 / −2, recompute and name the `mol_true` file), 205-214 (provenance comments), 240 (unchanged), 286-304 (Table 1 rows, plan 1.3), 347-349 (Table 2 rows), 353-372 (footnote and comments), 382, 386, 388, 404, 406, 410-414 (values; `:414` also carries A3), 419 (Fig 4 regenerate), 427-435, 455-461, 469-474, 476 (cluster, after the Phase 0.3 loader fix), 489, 493, 514-525, 530 (interaction), 715, 737-741, 751 (cost), 760, 779, 790, 792, 806, 817, 824, 836, 870.
- `body_appendix_short.tex`: 77 (plan row dropped, no such text), 173, 184, 189 / 235 / 241 (captions re-check, images regenerate), 218-225 (Table 8 rows), 244, 248, 269, 549, 893-920 (Wilson table), 947-964 (Table 19 with the star flips at DiffDock + smina and + gnina k = 5), 992-999, 1026-1031, 1056-1062 (DiffDock raw k = 1 / 5 only), 1074 (ITT), 1172-1179, 1207 (148 → 185 becomes [161 → 198]), 1247-1250 (Table 25), 1258 (values, see W49), 1303-1340 (Table 26), 1377-1385 (Table 27 cells), 1394, 1407, 1415.
- New material the plan requires that has no current text to quote: the D9 sensitivity table (single-instance convention with Form and Combined columns plus the single-copy stratum row [88 / 80 / 53 of 165 = 53.3 / 48.5 / 32.1] % rank-1 on the PB-valid ∧ ≤ 2 Å gate) beside A:1074, and the ITT sentence rewrite at A:1074 (rank-1 none within 2 Å, closest [6.3] Å; top-15 none but closest [2.1] Å, 7FRX_O88; [149 / 308] and [213 / 308] against [132] and [179]; margins [5.5] and [11.2 → 11.0]).

---

## Deviations from the plan, and what was verified

1. **A1 `\ref` claim.** The plan says the `\ref{supplementary-statistics}` at M:113 \"points at nothing\". The label exists (A:1067) and the chapter contains the Table 27 footnote (A:1387) that already states the correct count of two. Recorded under A1; the proposal keeps the `\ref` and adds a running-prose sentence (A1b).
2. **A4 gate description.** The plan says A:871 calls tables \"validity-aware\" that print the near-native gate. The Wilson table prints both gates and the quoted 34.0 / 34.3 are its validity-aware cells; the McNemar tables test the near-native gate. The proposal names both gates and says which the tests use.
3. **B2 Form-loss count.** The plan's `:666` draft says \"three reported variants\"; its own Section 1.3 and Section 8 say five printed Table 1 rows. Written as [five], value from the Phase 3 rebuild.
4. **W24a rounding.** The plan's +8.6 pp rounds to \"nine\", not the current \"eight\"; both readings are offered.
5. **Anchor drift.** Plan row A:545-553 (interaction values) sits at A:543 and A:549 in the live file. Every other plan anchor matched the live line.
6. **W36.** The plan calls for a rewrite of the plateau sentence; a draft is supplied, but it rests on previews (124 / 124 / 131, 144 / 151) and must be re-read after the rebuild.

Verification performed on 2026-09-08 (all read-only):
- Every \"Current\" block was matched as a substring of the named line of the live file by the generator (a miss aborts generation). Ranged quotations (T:89-90, T:91-100, T:98-99) are the full lines.
- Canonical table check for A1: `per_pose_metrics.csv`, 241,913 rows; rank-1 PB-valid own-RMSD > 2 with `pb_rmsd` ≤ 2: 7PGX_FMN (AutoDock exh128 + gnina, `rank == 1`), 6XBO_5MC (EquiBind unguided + gnina, lowest `gnina_affinity`), none for DiffDock + smina; reverse direction empty in all three; `pb_rmsd` > `rmsd` + 1e-3 on 0 rows.
- A2 source: `18_topn_within_thresholds_kabsch_pbvalid_depths_report.txt` lines 47-58 (rank-1 8.997e-09 for both contrasts, deeper depths ≤ 3.556e-11).
- A3 source: A:1383 prints 2.977.
- A6 source: `Data/PoseBuster Benchmark Set/README.txt` lines 51-63.
- A4 source: A:881-921 column headers and rows `:905`, `:909`.

No file other than this one was created or modified.
