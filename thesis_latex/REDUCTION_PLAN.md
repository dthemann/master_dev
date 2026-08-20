# Main-Body Reduction Plan
**Target file:** `/home/manndo/master_dev/thesis_latex/body_main.tex` (1,224 lines, 37,570 words)
**Destinations:** `/home/manndo/master_dev/thesis_latex/body_appendix.tex` (Chapters B, C, D, E, H)
**Status:** plan only. No file was edited.

---

## 1. VERDICT

About **8,600 words, or 23 per cent of the main body, is genuinely reducible**, and roughly 7 to 8 printed pages of floats can leave the Results chapter on top of that. The claimed 25 to 35 per cent is not reachable without either damaging the caveat spine or moving evidence away from the prose that reads it.

The bloat has exactly three shapes, in descending mass.

**Float walkthrough is the largest single source.** The Results chapters narrate tables and figure panels cell by cell. Lines 802 to 804 spend 627 words on four panels of one figure, line 835 spends 425 on a nine-cell Kendall tau matrix, line 997 spends 373 reciting counts that panel (B) already prints, and lines 704 to 710 re-read the quartiles of Table 9. The figures already draw all of it.

**Fourfold restatement is second.** Nearly every headline quantity is reported at a Results primary site, again in a "taken together" paragraph closing the same subsection, again in the Discussion, again in the Conclusions, and for six of them again in the abstract. The Orai chapter carries three consecutive summary paragraphs (L1060, L1062, L1064) that each re-report the same funnel, and Discussion L1144 and L1164 are full re-derivations of Results L825 to L847 and L1064 rather than references to them.

**Literature-review duplication against Appendix B is third and is the cleanest free money.** Main L53 and L55 are byte-identical to appendix L135 and L137, and main L73 is an 83 per cent verbatim copy of appendix L185 with the appendix strictly the superset. The same text is printed twice in the submitted document.

**Where the body is already tight and must be left alone:** the Motivation (L2-15, 224 words, carries RQ1-RQ4 verbatim), the definitional core of Evaluation Metrics (L234-262), the caveat spine (L286, L288, L316, L436-442, L964, L976, L1058, L1108, L1160), RQ2 to RQ4 (L1195-1209), and the Limitations chapter apart from its duplicated Orai provenance. Nine separate proposals against these ranges were tested and refuted, most of them because the real saving was 10 to 30 words for a whole-paragraph rewrite of audit-hardened text.

---

## 2. THREE TIERS

Line ranges are deduplicated across all tiers. Collisions are resolved and stated.

### TIER 1 — SAFE (risk none/low, no meaning loss) — **5,360 words**

| Lines | Action | What it is | Words |
|---|---|---|---|
| 712-754 | move → App H | Landscape Table 9 `tab:results-placement-form`, nine rows by fourteen columns of in-place-versus-form distributions plus footnote | 652 |
| 1007-1013 | move → App D | Figure `fig:results-r18` plus its whole filtering-effect analysis (benchmark-only robustness check, changes no conclusion) | 313 |
| 42-46 | compress | AutoDock Vina / Vina score / Vinardo exposition, duplicated at appendix L22-130 with equations and weight table | 205 |
| 430-434 | compress + move BH sweep → App C | PoseBusters decomposition read cell by cell, plus the per-check Benjamini-Hochberg family | 200 |
| 835 | compress, tau matrix → App H | Figure r8 rank-resolved F1 recital plus nine BH-corrected Kendall tau cells | 200 |
| 216-218 | compress | Receptor-protonation block (twin of App D L378) plus the GSK-7975A neutral-form story told twice | 200 |
| 144 (footnote) | move → App H | 238-word intention-to-treat sensitivity footnote plus superseded-run identifiers | 166 |
| 547-574 | move → App H | Table 8 `tab:results-depth-gain` plus its framing and trailing reader paragraphs | 164 |
| 98-120 | move → App B | Longtable `tab:literature-comparisons`, four rows of other people's comparative studies | 163 |
| 310 | move roster → App H | The family roster only (eighteen cluster-quality cells, four consensus, fifteen residue, fifteen type, nine tau). **L312 and L314 untouched.** | 160 |
| 912-918 | move → App D | Figure `fig:results-r10`, the Fr300 render whose own caption says it carries no measurement | 150 |
| 897-903 | compress | Four paragraphs walking the Orai yield table row by row | 145 |
| 130-132 | compress | Chapter-opening precis of every metric defined again 100 lines later | 145 |
| 1144 | compress | Discussion interaction paragraph, full re-derivation of L825/L835 plus a doubled hedge sentence pair | 139 |
| 57-65 | compress | DiffDock and EquiBind subsections against appendix L142-177 | 130 |
| 73-78 | compress | smina and gnina subsections against appendix L185-196 | 127 |
| 587 | move → App H | Power, minimum-detectable-difference and equivalence-margin derivation | 121 |
| 146-174 | delete | Table 2 `tab:methods-calibration-characteristics`, row-for-row duplicate of `tab:appendix-ligand-distribution` | 111 |
| 83-91 | compress | gnina rescoring and refinement, duplicated at App C L216-220 | 110 |
| 265-267 | compress | Three-group PoseBusters check walkthrough of `tab:appendix-posebusters-checks` | 110 |
| 837-847 | move fig r9 → App H, compress | Matched / missed / spurious contact decomposition plus section-closing restatement | 110 |
| 962 | compress | Figure r13 usable-yield boxplot recital, eighteen distributional statistics | 110 |
| 1082-1084 | move → App C | Workstation specification plus ratio-of-totals reconciliation | 107 |
| 1071-1077 | compress | Dataset Summary, fourth site for every number it holds | 100 |
| 50-55 | compress | AI-docking landscape survey, byte-identical to appendix L135/L137 | 102 |
| 1118-1120 | compress | Discussion opener stating the endpoint split and the hedge twice in four lines | 95 |
| 1138-1142 | compress | Three Discussion paragraphs each opening by re-reading a triplet out of a Results table | 84 |
| 696 | compress | Figure 5 lead-in footnote, plotting provenance | 80 |
| 1171 | delete | "Current work points toward hybrid rather than exclusive workflows", third statement of one recommendation | 76 |
| 220-222 | compress | SOCE framing and site preview, restated at L288 where the criterion is defined | 75 |
| 908-910 | compress | Orai receptor description, twin of App D L378-382 | 75 |
| 585 | compress | Four-depth Newcombe interval recital plus the estimator function name | 70 |
| 694 | compress | Headroom summary, third statement of the 62-against-64 gain | 70 |
| 502 | move → App H | Depth-family McNemar / Wilson / Holm procedure recital | 67 |
| 1024-1030 | move → App C | Figure r19 total-contacts null result plus its four effect sizes | 65 |
| 759 | compress | Cross-Tool Clustering opener restating the Methods 4 Å rule | 62 |
| 668-672 | compress | Stranded "band between the curves" plus AutoDock and DiffDock walk | 60 |
| 817 (footnote) | move → App C | Receptor-numbering footnote, twin of App C L371 | 60 |
| 808-812 | move → App H | Figure `fig:results-r6` six-panel cluster-quality assessment | 50 |
| 972-974 | compress | Transmembrane-loss rates reported twice in consecutive paragraphs | 47 |
| 1188 | compress | RQ1 inline interval endpoints, cohort size and Wilcoxon p | 46 |
| 806 | compress | Cluster section opening restatement (pooled oracle figures retained) | 40 |

### TIER 2 — MODERATE (low/medium risk, texture loss but no claim) — **2,760 words**

| Lines | Action | What it is | Words |
|---|---|---|---|
| 1060-1064 | compress | Three closing Orai synthesis paragraphs, 861 words re-reporting L897/L964/L974/L1015 | 320 |
| 204-213 | move → App E | Pharmacology rows of Table 3 plus the 130-word potency footnote. **Descriptor rows and the label stay in the body**, so the `\ref` at L218 and "The descriptors above" at L216 both survive | 250 |
| 946-954 | move → App D/C | Placement opener third restatement plus Figure `fig:results-r12` with its 134-word disclaimer caption | 218 |
| 274-284 | move → App C | Placement-envelope thresholds, pore-axis vector construction, two display equations | 195 |
| 802-804 | compress, omnibus → App H | Figure r5 panel-by-panel walkthrough. **Rescued form**: keep the 232/204/138 denominators and the McNemar p = 0.040 in the body, move only the Friedman / Cochran-Q / rank-biserial / Kruskal-Wallis block | 190 |
| 854-858 | compress | Orai chapter roadmap, per-frame drop breakdown, pose totals | 182 |
| 308 | move → App H | Design-to-test catalogue, as a three-column table. Equivalence provision stays in the body | 173 |
| 326-332 | compress, p-values → App H | Variant-selection protocol sitting in Results, three table walkthroughs and four p-values | 170 |
| 1216 | move → App C | Orai run provenance inside the Limitations paragraph, twice-stated pose accounting | 166 |
| 1164 | compress | Discussion cross-panel ledger, near-verbatim restatement of L1064. Future-work clause relocates to L1224 | 150 |
| 1038 | move → new App table | Interaction-type per-cell enumeration, ten mean pairs and eight effect sizes | 150 |
| 933-937 | compress | Figure r11 recital plus ceiling-and-floor paragraph plus six-effect-size statistics recital | 200 |
| 1110 | move triplets → App C table | Alternative-denominator cost arithmetic. **Prerequisite: create the appendix cost table first** | 120 |
| 406-411 | move rows → App C | Six pocket-guided EquiBind rows of `tab:results-pose-production`, never read anywhere | 110 |
| 426-428 | compress | Figure 1 Friedman / Kendall W / Hodges-Lehmann / rank-biserial recital | 90 |
| 767 | compress | Reach rates reported twice over plus seven phi cells | 75 |

### TIER 3 — AGGRESSIVE (medium risk, only under a hard limit) — **510 words**

| Lines | Action | What it is | Words |
|---|---|---|---|
| 997 | compress | Figure r16 per-pair recital. Keeps the delta +0.91, the two-ligand caveat, the 15 per cent rate and the three-way result | 120 |
| 1005 | compress | Figure r17 cluster-quality recital. Calinski-Harabasz artefact flag must stay named, not anonymised | 100 |
| 1046 | compress | Figure r21 residue hot-spot recital. **Retain one DiffDock reversal pair (Asp110 or Leu109)** or the retained claim loses its only evidence | 100 |
| 596-665 | trim columns → App H | Reduce Table 8 to the 0.5, 1, 1.5, 2 and 2.5 Å columns. The other six thresholds are cited nowhere | 75 |
| 1169 | compress | Discussion efficiency medians, third of four sites | 45 |
| 1158 | compress | Discussion Orai yield restatement, fourth site | 40 |
| 1127-1129 | compress | Validity triple plus the four per-rank refinement bands, third statement | 32 |

**Contingency, not recommended.** Moving the whole of Table 8 (L593-666) to Appendix H yields 663 words but was refuted on three grounds, namely that Appendix H's preamble declares k = 1, 5, 10, 15 against this table's k = 1, 15, 30, that its variant naming differs, and that five body sites read it. It is mutually exclusive with the column trim above. Moving Table 6 (L444-500) yields only about 131 words of actual prose against 549 raw tokens and would strand four discussion sites, so take it only if the metric really is `wc -w` on the `.tex`.

---

### Collisions resolved

- **L146-174** had three separate proposals (77, 111, 157 words) for the same delete. Counted once at 111.
- **L1110** had three proposals. The naive compressions save 36 and 60. The 145-word version destroys nine numbers that exist nowhere else, namely the per-generated-pose currency triplets 0.89/5.73/9.53, 0.07/4.85/9.00 and 0.53/28.20/28.50. Resolved by requiring the appendix cost table to exist first, then counting 120.
- **L587** had two identical move proposals plus a partly-broken variant. Counted once at 121, using the version that keeps the 37.0 per cent discordant share and the `p_holm` 0.109/0.161 values.
- **L1060-1064** collided with three per-line proposals totalling 119 words. Resolved in favour of the single coherent rewrite at 320.
- **L802-804** collided with **L804-806**. Resolved by splitting at the line boundary, 802-804 in Tier 2 and 806 in Tier 1, so the pooled oracle figures at L806 (95 per cent at 0.35 Å, 94.7 per cent at 0.40 Å, the 1.30 Å cluster centre) are never inside a moved range.
- **L827-835** collided. The 827-833 figure move was refuted because it deletes the EquiBind "always starts at k = 1" clause. Resolved by taking the L835-only compression at 200 and leaving Figure r8 and its lead-in in the body.
- **L962-964** collided. Resolved by compressing L962 only, because L964 is the sole site of the pooled 71.07, 2.86 and 79.17 that RQ1 quotes as 71.1, 2.9 and 79.2.
- **L181-214** collided (whole-table move at 293 against half-table move at 250). Resolved in favour of the half move, which keeps the label and the `\ref` in the body.
- **L98-120** was counted at 163 rather than 260, and **L308/L310** were separated so the roster move and the catalogue move do not double-count.

---

## 3. APPENDIX MOVES BY DESTINATION

### Appendix B — Tools and Functions (L17)
| From | What | Label / ref work |
|---|---|---|
| L98-120 | `tab:literature-comparisons` → after the AI-Based section at appendix L137 | Sole `\ref` is body L96. Label travels with the table, cross-chapter `\ref` resolves. No other edit. |

Note that Appendix B is a compromise destination, since it documents the tools this thesis ran rather than the literature. Take it only if that is acceptable.

### Appendix C — Docking Protocols (L203)
| From | What | Label / ref work |
|---|---|---|
| L274-284 | Placement-envelope thresholds, pore-axis construction, two display equations | No `\label` in range. Nothing breaks. |
| L406-411 | Six pocket-guided EquiBind rows out of `tab:results-pose-production` | Row move inside a landscape longtable. Table keeps its label and all nine `\ref` sites. |
| L432 | Per-check Benjamini-Hochberg sweep → beside `tab:appendix-posebusters-checks` | No `\label` in range. |
| L817 footnote | Receptor-numbering explanation → Interaction Fingerprints | Footnote only. |
| L854-858 | Per-frame drop breakdown, DiffDock 89-nine-pose / ten-one-pose split, pose totals → Orai1 Panels | **Prerequisite.** The replacement adds a `\ref{orai1-panels}` for the totals 12,010 / 12,150 / 11,673, which Appendix C does not currently carry. Add them or the pointer is false. |
| L948-952 | Figure `fig:results-r12` → Orai1 Panels | Sole `\ref` at L954, inside the range. |
| L1024-1030 | Figure `fig:results-r19` → Interaction Fingerprints | Sole `\ref` at L1030, inside the range. **Rename the label** to `fig:appendix-orai-contacts` and repoint the surviving body sentence. |
| L1082-1084 | Workstation specification, ratio-of-totals reconciliation, 7WCF_ACP exclusion → new "Hardware and Timing Basis" section | No `\label` in range. |
| L1110 | Three per-currency cost triplets → new per-denominator table in the same section | **Create the table first.** These nine values have no other site, and `fig:results-r24` carries a different axis. |
| L1216 | Two DiffDock run defects and the panel pose accounting → Orai1 Panels, beside `tab:appendix-orai-protocols` | The defect disclosure exists nowhere else, so it must move rather than be deleted. |

### Appendix D — Orai1 Receptor Model (L376)
| From | What | Label / ref work |
|---|---|---|
| L912-918 | Figure `fig:results-r10` and the pre-superposition 7.3 to 7.5 Å figures | Sole `\ref` at L918, which the replacement rewrites. Keep the Kabsch 5.90 / 5.88 / 6.02 Å values in the body, since Discussion L1164 leans on them. |
| L1007-1013 | Figure `fig:results-r18` and its whole analysis paragraph | Sole `\ref` at L1013, inside the range. Self-contained. |
| L946-954 (part) | The placement-render half | See Appendix C row for `fig:results-r12` if the render is preferred beside the protocol. |

### Appendix E — Experimental Ligand Pharmacology (L391)
| From | What | Label / ref work |
|---|---|---|
| L204-213 | Eight pharmacology rows of `tab:methods-experimental-ligands` plus the potency footnote | Label and descriptor rows stay in the body. `\ref` at L218 unaffected, "The descriptors above" at L216 stays true. |

### Appendix H — Ranking-Recovery Statistics (L765)
| From | What | Label / ref work |
|---|---|---|
| L144 footnote | Intention-to-treat sensitivity check and excluded identifiers | Footnote only. |
| L308 | Design-to-test-to-effect-size catalogue, as a table | Equivalence provision stays in the body for RQ1. |
| L310 | Family roster only | **L312 and L314 must not move.** L312 is the sole statement of the add-one permutation convention and its 1e-4 floor, quoted seven times in Results. L314 holds the listwise-deletion and "non-constant values" rules, which occur once in the thesis. |
| L502 | Depth-family McNemar / Wilson / Holm procedure and the nested-pool one-sidedness argument | No `\label`. |
| L547-574 | Table `tab:results-depth-gain` | Keep the label unchanged. Both `\ref` sites resolve cross-chapter. |
| L587 | Power, MDD, 1,400-complex projection, per-depth TOST margins | Twelve-point bound stays in the body for RQ1. |
| L712-754 | Table `tab:results-placement-form`, inside its `\begin{landscape}` wrappers | Contains `\ref` to `tab:results-pose-production`, `tab:results-near-native-form`, `fig:results-r3`, `fig:results-r2a`, `fig:results-r2b`. All resolve from the appendix. Its own label is `\ref`'d from surviving body prose at L696 and L704-710, which becomes a forward reference. |
| L802-804 (part) | Friedman, Cochran-Q, rank-biserial and Kruskal-Wallis omnibus block | Keep the McNemar p = 0.040 and the 232/204/138 denominators in the body. |
| L808-812 | Figure `fig:results-r6` | **Three `\ref` sites**, at L295 and L299 in Methods and L804 in Results. The two Methods sites become forward references into the appendix and must be reworded. |
| L835 | Nine-cell Kendall tau matrix | No `\label`. |
| L839-843 | Figure `fig:results-r9` | Sole `\ref` at L837, inside the range. |
| L1038 | New table `tab:appendix-orai-interaction`, thirty tool-by-type cells at both inferential units | **Label does not exist.** Create it or the new body `\ref` dangles. |

**Blocking prerequisite for Appendix H.** Its preamble at `body_appendix.tex:766-767` declares that recovery is reported for k = 1, 5, 10 and 15, names variants as "DiffDock (raw)" and "DiffDock (gnina-opt)", and asserts that every corrected p value in the chapter refers to a variant printed there. Material arriving at k = 1/15/30, under "+ smina" naming, and carrying clustering and interaction statistics falsifies all three statements. Rewrite that preamble before any move lands.

---

## 4. THE FIVE HIGHEST-VALUE SINGLE ACTIONS

Together these five edits remove **1,536 words** and one landscape page from the body.

### Action 1 — move Table 9 to Appendix H (652 words)
**Cut** `body_main.tex` L712-754 entire, including the `\begin{landscape}` and `\end{landscape}` wrappers, and paste unchanged into Appendix H. **Replace in the body with:**

> Table~\ref{tab:results-placement-form} in Appendix~\ref{ranking-recovery-statistics} gives the underlying distributions by tool and ranking depth.

The cells the argument actually uses (58.4 per cent coverage, the 1.9 to 5.6 Å drift, the 51 to 81 per cent placement-limited swing, the 2.34 against 5.74 Å diversity spreads) are all quoted in surviving prose at L696 and L704-710. This is a word saving only, since the landscape page reappears in the appendix.

### Action 2 — move the filtering-effect analysis to Appendix D (313 words)
**Cut** L1007-1013 entire, figure and paragraph together. **Replace with:**

> Filtering the benchmark pose clouds to poses that are PoseBusters-valid and outside the transmembrane exclusion slab does not tighten them. Within-cluster spread is essentially unmoved for all three tools, while bootstrap stability rises significantly for all three and DiffDock also gains on silhouette, so what filtering improves is reproducibility rather than tightness. The paired per-unit tests, the reduced EquiBind subset on which its contrast rests and the degenerate single-cluster clouds that inflate its stability gain are reported in Appendix~\ref{orai1-receptor-model}. Filtering therefore does not detectably tighten the wide scatter along the pore, so on this evidence it is not an artefact of invalid or mislocated poses that filtering could clean away. No equivalence test was run, so this bounds the effect of filtering rather than establishing that the scatter is intrinsic to how the tools sample.

Nothing in the Conclusions, the abstract, the Discussion or the subsection summary at L1015 depends on any number in the moved block. `fig:results-r18` is referenced only at L1013, inside the range. This is a word saving **and** a page saving of close to one page.

### Action 3 — compress the AutoDock Vina exposition (205 words)
**Replace** L42-46 with:

> AutoDock Vina couples an efficient global search engine to an interchangeable empirical scoring function \cite{ref022}, \cite{ref023}. Engine and objective are worth keeping apart, because the sampler is largely independent of the function that ranks its output, and two functions are relevant here, namely the original Vina function and the later Vinardo reparametrisation, both selectable within the same program \cite{ref023}. The engine generates random poses inside a user-defined search box, refines them by iterated local search over precomputed interaction grids, clusters the survivors by RMSD and returns the lowest-energy representatives. Because the chosen function is evaluated together with its gradient, each pose is minimised against it during the search, so no separate post-docking optimisation is needed.
>
> The Vina function itself is a weighted sum of two attractive Gaussian steric terms, a repulsive parabola, a hydrophobic term and a hydrogen-bond term, divided by a factor growing with the number of rotatable bonds, with its six weights fitted once against the PDBbind set rather than derived from first principles \cite{ref022}. Vinardo keeps those families of terms but drops the second, longer-range Gaussian that tended to seat ligands too deeply in the pocket, re-optimises the atomic radii and re-fits the remaining weights, and its parameters were selected on redocking behaviour rather than on correlation with measured affinities \cite{ref073}.

Everything removed is present at greater length in Appendix B L22-130, in Equations eq:vina-1 to eq:vina-9 and in `tab:appendix-vina-weights`. The Vinardo introduction that Results L330 needs is preserved.

### Action 4 — merge the Orai protonation paragraphs (200 words)
**Replace** L216-218 with:

> Limited receptor flexibility was represented with four supplied Orai1 molecular-dynamics snapshots, namely Fr0, Fr300, Fr400 and Fr499 \cite{ref038}, \cite{ref072}. No pH-based protonation was applied to them, and the polar hydrogens the docked receptors carry were added by the PDBQT converter under its own geometric rules rather than by a titration model, as detailed in Appendix~\ref{orai1-receptor-model}. The ligands were prepared for a nominal pH of 7.0 and all three were docked as neutral species, namely native 2-APB (2-aminoethoxydiphenyl borate, labelled 2abp-NH2 in the per-compound results tables), GSK-7975A \cite{ref041}, \cite{ref070} and Synta-66 \cite{ref071}, so the reported 2-APB pharmacology \cite{ref040} characterises the docked molecule directly. A deprotonated phenolate form of GSK-7975A was also prepared during ligand setup and was corrected during this work, its earlier preparation having produced a neutral radical with a dearomatised ring rather than the phenolate, and it enters no analysis reported here, so Table~\ref{tab:methods-experimental-ligands} and Table~\ref{tab:appendix-orai-descriptors} both describe the neutral species (C\textsubscript{18}H\textsubscript{12}F\textsubscript{5}N\textsubscript{3}O\textsubscript{2}, 397.3~Da) that underlies every Orai1 result. No cognate Orai1 crystal pose is available for any of the docked compounds, so this set supports a reference-free case study rather than a conventional accuracy benchmark.

The radical-file correction record and the reference-free scope sentence both survive. The CHARMM-hydrogen and HSD/HSE/HSP detail lives in Appendix D L378 already.

### Action 5 — move the intention-to-treat footnote to Appendix H (166 words)
**Cut** the 238-word footnote hanging off "the 302 cases with comparable timing records" at L144 and paste it verbatim into Appendix H. **Replace in the body with:**

> Five of the 308 source entries produced no DiffDock output and were dropped before the three-tool comparison, leaving the 303 complexes analysed throughout, and one further complex carries no usable wall-clock record and is excluded from the timing comparison alone. An intention-to-treat sensitivity check that retains the five and scores them as DiffDock failures, together with the identifiers of every excluded entry, is reported in Appendix~\ref{ranking-recovery-statistics}. On that basis the DiffDock margin over AutoDock Vina + gnina narrows at both depths without changing direction.

Move the footnote verbatim rather than trimming it, since it was added to answer an audit. The abstract and the Conclusions quote the 303-basis rates, which stay.

---

## 5. WHAT MUST NOT BE CUT

### Protected because a research question or the abstract depends on it
- **Table 4 (`tab:results-pose-production`, L334-415).** Primary site for 98.2, 99.5, 24.4, 2.9, 84.5, 87.1 and 51.8 per cent, six of which the abstract quotes.
- **Table 7 (`tab:results-depth-recovery`, L506-545).** Sole body site for 34.0 / 55.1 / 55.8, 29.4 / 49.8 / 49.8, 15.2 / 20.1 / 20.5 and the +62 / +64 / +15 gains that RQ1 quotes.
- **L440-442.** The per-variant, per-depth enumeration is the derivation of the 0.00 to 1.32 point span that the abstract states as "at most 1.3 percentage points". Note the trap, the 15.5 at L442 is a **complex-level** span for raw DiffDock and is a numerical coincidence with the abstract's pose-level 15.5 per cent. Never merge them.
- **L964.** Sole site of the pooled 71.07, 2.86 and 79.17, which RQ1 and Discussion L1176 quote as 71.1, 2.9 and 79.2, and the only cross-panel contrast that survives correction at both units.
- **L1094.** Cost medians 8.5 / 50.9 / 119.7 s and the 168 / 151 / 62 counts that RQ3 quotes.
- **L585 rank-1 and top-15 intervals** and the twelve-point TOST bound. Only the top-5 and top-30 intervals may leave.

### Protected because they are deliberately added hedges
- **L286 and L288**, the interpretation caveats and the operational-hypothesis statement on the placement rule.
- **L316**, the two overarching caveats that the abstract's optimism disclosure rests on.
- **L436-438**, the circularity qualification. Rewriting it nets 72 words and hollows out a rebuttal into a content-free promissory sentence.
- **L976, L1058, L1108, L1160**, the pseudoreplication, power and "not ground truth" statements.
- Every occurrence of "inconclusive rather than positively equivalent", "no equivalence test was run", "informative rather than random" and "absence of detected difference".

### Tempting but refuted, with the reason
- **L6, the Motivation history paragraph.** Real net saving is 10 words, not 42, because the replacement keeps all five citations and the closing sentence verbatim.
- **L91, the gnina design rationale.** Real saving 26 words. It is also the only body statement of the allocation rule that rescoring goes to the learned pipelines and refinement to the AutoDock variants, which is the premise for the six AutoDock variants in Table 4.
- **L299, the cluster-quality axes.** Real saving 25 words, and "adjusted Rand" and "nearest-cluster separation" occur exactly once in the thesis, here. No caption carries them.
- **L312 and L314.** L312 is the sole statement of the add-one permutation convention and the 1e-4 floor, which is quoted as a result seven times in the Results. L314 is the sole statement of the listwise-deletion rule and of "non-constant values". Move the L310 roster only.
- **L444-500, Table 6.** Four surviving body discussions read it, at L430, L434, L436 and L901. About 131 of its 549 tokens are prose.
- **L704-710.** Contains "by top-5 and top-15 it is significantly the best tool on both axes", which appears nowhere else and is the only assertion that DiffDock is statistically best on the continuous medians.
- **L589-591.** The Wilson intervals on the depth gains and the per-tool Cochran's Q values 124 / 128 / 30 exist nowhere else. Table 7 carries no intervals.
- **L806.** The 1.30 Å pooled-cluster-centre sentence with its "consequence of averaging rather than a cost of pooling" clarification, and the 94.7 per cent AutoDock-plus-DiffDock anchor, are both unique. Without the latter, "adding EquiBind moves reach by 0.3 percentage points" has no antecedent.
- **L827-833.** "EquiBind's gnina-affinity rank is assigned over its retained valid poses and so always starts at k = 1" is the only record of that cohort asymmetry, and it is what qualifies EquiBind's 0.46 F1.
- **L1176.** "retain uncertainty when methods disagree" appears in neither of its two twins at L1205 and L1212, and it is the stage the near-zero agreement result demands.
- **L1212.** "crystal-blind selector" occurs once in the thesis. RQ4's twin says only "validated for the intended endpoint", which is weaker.
- **L1222.** The six named future directions belong in the Future Outlook, not as a back-reference to a citation list one chapter earlier.
- **L1110 currency triplets.** Nine values with no other site. `fig:results-r24` carries the per-qualifying-pose axis, not the per-generated-pose axis.

---

## 6. BOTTOM LINE

| Stage | Words removed | Body word count | Reduction |
|---|---|---|---|
| Start | — | **37,570** | — |
| After Tier 1 | 5,360 | **32,210** | 14.3 % |
| After Tier 1 + 2 | 8,120 | **29,450** | 21.6 % |
| After all three tiers | 8,630 | **28,940** | 23.0 % |

**Page savings, counted separately.** Deletions that remove a page from the whole document, not merely from the body, are Table 2 (L146-174, about one third of a page) and the eight pharmacology rows plus footnote of Table 3 (about half a page), so the document shortens by roughly one page. Floats that leave the body but reappear in the appendix are Table 1, Table 8, Table 9 landscape, Figures r6, r9, r10, r12, r18 and r19, which together release about **7 to 8 pages from the Results and Methods chapters** while leaving the total page count broadly unchanged. If total pages rather than body words is the real constraint, the only levers are the two deletions and the Tier 3 column trim of Table 8.

**Honest assessment of the gap to 35 per cent.** Reaching it would require the two refuted contingencies (the full Table 8 move at 663 words and the Table 6 move at 131 prose words) plus roughly 1,500 further words that can only come out of the caveat spine, the Conclusions or the Limitations chapter. Every one of those was tested and refuted. Twenty-three per cent is what this body gives up without losing an argument.