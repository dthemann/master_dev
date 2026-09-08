# Adversarial revalidation of the five hand-verified defects

Date: 2026-09-03. Source under test: `thesis_latex/obsolete/ABSTRACT_REVIEW_2026-09-03.md`, section G,
"Verified by hand, new".

Method: each claim was independently re-derived from primary artefacts, then attacked by three
skeptics with distinct lenses (arithmetic/provenance, semantics/scope, consequence/fix-correctness),
then adjudicated on the facts. 25 agents, 0 errors. Every number below was also re-derived by hand
against the generating data files, not against the prose.

## Verdicts

| Claim | Report said | Verdict | Severity |
|---|---|---|---|
| main :683 vs appendix :471 | direct contradiction, 14 vs 11 and 0.086 vs 0.142 | **RECLASSIFIED** | major |
| main :591 slab figure | cites `fig:results-r10`, slab is `fig:results-r12` | **CONFIRMED** | minor |
| main :97 appendix letter | cites Appendix E, statement is in Appendix F | **CONFIRMED** | minor |
| appendix :1192 | "55 to 58" and "19.1" wrong | **CONFIRMED** | major |
| main :474 | "between six and seven points" excludes 7.3 | **REFUTED** | trivial |

Four of five survive. One is refuted. One is real but was mis-diagnosed, and the mis-diagnosis
pointed at the wrong fix.

---

## C1 RECLASSIFIED — not a contradiction, a stale appendix table

The report framed this as two passages disagreeing, which invites an edit that aligns the sentences.
That fix would have been wrong. The body was right and the appendix was a whole generation stale.

All seven AutoDock cells of Table C.12 matched
`obsolete/pandamap_results/orai_interaction_compare_PRE_FR0/` exactly, digit for digit, including
Cliff's delta. That generation independently reproduces every disputed appendix number:

- BH-clearing cells in PRE_FR0: AutoDock 5 + DiffDock 6 = **11**, the appendix's figure.
- AutoDock ligand-unit minimum in PRE_FR0: **0.142**, the appendix's figure.
- AutoDock repulsion in PRE_FR0: BH 0.052, which does not clear, so the footnote's "AutoDock clears
  neither" was true then.

Current matched data (`pandamap_results/orai_interaction_compare_matched/`, and
`posebusters_results/_orai_matched_root/.../orai_ligand_level_contrasts.txt`):

- BH-clearing cells: AutoDock 8 + DiffDock 6 = **14**. The body's "fourteen of twenty-six" was correct,
  and so was "twelve interpretable cells" after discounting the two halogen cells.
- AutoDock ligand-unit minimum: **0.089** (halogen bonds). DiffDock: **0.106**, tied across four types.
- AutoDock repulsion: BH 0.047, which does clear.

The denominator 26 is correct in both passages. The generating script collapses the triplicated
charged classes to one hypothesis, giving 13 types per tool across the two mature tools. Every
steelman for 11 was tested and failed. No assignment of the charged triplicate produces 14 and 11
against a common denominator of 26.

The body's "0.086" matched neither generation. It is the value from the non-matched run at n=308.

Two further appendix sentences were false under current data and the report did not flag them.
"Every one of them appears above" failed, because AutoDock carbon-pi at 0.018 and pi-pi stacking at
0.037 both clear but had no printed row. "The DiffDock effect sizes run in the same direction as the
AutoDock ones on every row" failed on pi-pi stacking, where AutoDock is -0.41 against DiffDock +0.07.

**Status: already repaired on disk by a concurrent session, and the repair is correct.** The table now
carries nine rows and all eighteen AutoDock and DiffDock cells match the matched run. The body reads
0.089 and the footnote reads fourteen. The repair is uncommitted, so a rebuild from HEAD would
reintroduce the defect.

## C2 CONFIRMED — wrong figure for the exclusion slab

Line 591 promises "the Fr300 snapshot with the transmembrane exclusion slab" and cites
`fig:results-r10`, which resolves to Figure 27 on page 103. I opened both renders.

- Figure 27 (`media/media/image12.png`) is a bare pale-yellow cartoon. No ligands, no slab.
- Figure 22 (`media/media/image14.png`, `fig:results-r12`) shows the tan cartoon, the dense yellow
  sphere layer, the magenta dot surface and green and blue ligand sticks.

The half of the sentence naming Fr300 is right, since Figure 27 is Fr300. Only the slab clause is wrong.

Preferred fix is to delete the false clause rather than re-target the reference, because line 613
introduces `fig:results-r12` as a first mention twenty-two lines later and swapping it in at 591 would
force that sentence to be rewritten too.

## C3 CONFIRMED — Appendix E cited for a statement in Appendix F

The label `experimental-ligand-pharmacology` resolves to Appendix E, page 105. I read that chapter in
full. It carries potencies, assays, proposed sites and selectivity liabilities, and says nothing about
the optimisation record. The statement lives at `body_appendix_short.tex:791`, inside
"Contrast with the Functionally Tested Ligands", which resolves to Appendix F.9 on page 115.

Fix is to cite `contrast-with-the-functionally-tested-ligands`, which renders as F.9. Section-level
appendix references are already the house form in the main body.

## C4 CONFIRMED — and the same sentence carries two more errors

Table C.x, "Accurate-but-Invalid Share by Ranking Depth", gives the raw DiffDock row as
15.5 / 19.1 / 19.5 / 18.8 percent of N=303. Each percentage maps to a unique integer count, so the
gaps are forced:

| depth | near-native | validity-aware | complexes removed |
|---|---|---|---|
| k=1 | 102 | 55 | 47 |
| k=5 | 135 | 77 | 58 |
| k=10 | 146 | 87 | 59 |
| k=15 | 151 | 94 | 57 |

The prose says "47 complexes at rank-1 and 55 to 58 at the deeper pools, a span of 15.5 to 19.1".
"Deeper pools" is fixed by the same paragraph's AutoDock sentence, which separates top-5 from the
deeper pools and is exactly right at 2 / 4 / 3. Under either reading the DiffDock range is **57 to 59**,
and the maximum share is 19.5 at k=10, not 19.1. The stray 55 is the k=1 validity-aware count, taken
from the wrong column.

Two further errors in the same sentence, not in the report:

- "one to three for unguided EquiBind + gnina" — the gaps are 2, 3, 3, 3. It should read two to three.
- "a span of 0.33 to 1.32 percentage points" — the smallest gap across the three variants is 2 of 303,
  which is 0.66 points, not 0.33. The upper bound of 1.32 is right.

## C5 REFUTED — the band is defensible under the paragraph's own convention

All three skeptics refuted. I agree.

Exact rates from `posebusters_results/cluster_crystal_pocket_matched_equibind/.../ranking_ablation.csv`:
legacy 150/303, consensus 172/303, confidence 169/303, medoid 168/303, oracle 287/303. Gains are
7.26, 6.27 and 5.94 points.

The same paragraph writes "The 38-point oracle gap" for an exact 37.95, which establishes a
whole-point convention. Under it the gains are 7, 6 and 6, and "between six and seven points" holds.
The four underlying rates are printed verbatim one clause earlier, so a reader can compute the exact
gains. Nothing downstream depends on the magnitude.

One incidental note. The report's own arithmetic used the thesis's printed 55.5 percent for the
medoid rule to get a gain of 6.0. The data give 55.4455 percent, so the true gain is 5.94. The printed
55.5 is itself a double-rounding of 55.45. This does not change the verdict.

---

## Open items

The seven candidates listed as "not yet adversarially verified" in the review report were not
examined in this pass.

---

# Edits applied, 2026-09-03

Rebuilt with `latexmk -f -pdf -interaction=nonstopmode -file-line-error Thesis_short.tex`.
Still 142 pages, zero undefined references. Every new value verified present in the PDF and every
superseded value verified absent.

## C1 — no edit needed, verified complete

The table and its footnote had already been regenerated from the matched run by a concurrent session.
I verified rather than redid. Table 12 now carries nine rows and all eighteen AutoDock and DiffDock
cells match `orai_interaction_compare_matched`. The row set is exactly the union of types clearing in
at least one tool, five clear in both, amide-pi is DiffDock only, and carbon-pi, pi-pi stacking and
charged classes are AutoDock only, which sums to the stated fourteen cells. A sweep for the retired
markers 0.086, 0.142 and "eleven of the twenty-six" returns nothing. The Figure 26 caption already
carries the current 7,458 poses and Asp110 at 28.8 percent. No further prose is stale.

## C4 — four corrections plus two scope repairs

Validated twice more before editing, once by decoding each printed percentage to its unique integer
count out of 303, and once by an independent agent that counted the per-complex boolean columns
directly. Both agree with the table.

| was | now |
|---|---|
| 55 to 58 at the deeper pools | **57 to 59** |
| a span of 15.5 to 19.1 percentage points | **15.5 to 19.5** |
| one to three for unguided EquiBind + gnina | **two to three** |
| a span of 0.33 to 1.32 percentage points | **0.66 to 1.32** |
| at most four complexes of 303 at any depth | at any **reported** depth |
| three to four at every depth for DiffDock + smina | three to four for DiffDock + smina |

The stray 0.33 is one complex of 303, and one complex is the raw AutoDock Vina figure, an arm the
sentence explicitly excludes. The stray 55 is the rank-1 validity-aware count read from the wrong
column.

The last two rows are a scope repair, not a transcription fix. I re-ran the project's own recovery
script across all thirty depths. Over depths 1 to 30 the DiffDock plus smina gap reaches 6 at depth 2
and 5 at depths 3 and 4, so an unbounded "at any depth" was false. Over the four depths the tables
actually report it is true, so both sentences are now scoped to the reported depths.

## C5 — replaced the band with explicit one-decimal gains

The band itself was defensible under the paragraph's whole-point convention, but the user asked for
one decimal so the error cannot recur. Two changes at `body_main_short.tex:474`:

- The medoid rate 55.5 percent was a double-rounding of the CSV's already-rounded 0.5545. The exact
  rate is 168 of 303, which is 55.4455 percent, so it now reads **55.4**.
- "each gains between six and seven points" now reads **"they gain 7.3, 5.9 and 6.3 points
  respectively"**, in the same order as the rates clause.

The printed rates and the printed gains are now mutually consistent: 56.8 minus 49.5 is 7.3,
55.4 minus 49.5 is 5.9, and 55.8 minus 49.5 is 6.3. Every other number in the passage was checked and
is correct, including the 38-point oracle gap, which resolves uniquely as the oracle minus consensus
ranking at 37.95 points.
