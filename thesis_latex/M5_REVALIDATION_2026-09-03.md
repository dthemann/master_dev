# Re-validation of examiner finding M5 (2026-09-03)

**Finding under test** — "Appendix H.4's placement-versus-validity decomposition never reaches the
Discussion, Conclusions or abstract."
**Source** — `thesis_latex/EXAMINER_REVIEW_2026-09-03_full_msc_rubric.md:305-315`
**Method** — 58-agent adversarial workflow: 8 independent audit dimensions, 2 refuters per material
claim, a completeness critic, a final adjudication. Every number below was then reproduced by hand
from the primary pose data.
**Adjudicated against** — the live short build, `Thesis_short.pdf` (142 pp), and
`posebusters_results/benchmark_matched_equibind/`.

---

## Verdict: PARTIAL — act, but not on the remedy M5 proposes

The literal claim is exactly true and the effort estimate is exactly right. The rationale is false,
and the abstract limb is a category error. M5 has also been filed before, verbatim, and never
actioned.

| # | Sub-claim | Verdict |
|---|---|---|
| 1 | "four percentage points" occurs once in the appendix, zero times in the body | **CONFIRMED** |
| 2 | The decomposition reaches neither §5.2, §6.1, the abstract nor the Kurzfassung | **CONFIRMED** |
| 3 | §5.2 "offers only the unquantified *partly constitutive*" | **PARTIAL** — "only" is false |
| 4 | §6.1 "prints the raw 98.9 / 84.5 / 62.0" | **PARTIAL** — not raw, and not bare |
| 5 | The abstract makes a reader over-read the validity gap | **REFUTED** — category error |
| 6 | "Most of the separation is a placement effect, not a geometry-quality effect" | **REFUTED** by the thesis's own data |
| 7 | "Effort. Rewriting only. The analysis is done." | **CONFIRMED** |

---

## 1. What M5 gets right

The string count is exact.

| File | occurrences |
|---|---|
| `Thesis_short.tex` | 0 |
| `body_main_short.tex` | 0 |
| `body_appendix_short.tex` | 1 (line 1184) |

One occurrence in the 142-page PDF, printed page 126. No German equivalent in the Kurzfassung.
The conditioned trio 97.5 / 94.6 / 93.6 appears nowhere outside the appendix table.

The operands are machine-asserted at `Scripts/Analysis/thesis_expected_values.yaml:656,662` and
recomputed by `check_table_24`. I re-ran the join independently: 241,913 rows, 0 unjoined, all
sixteen table cells reproduce to the digit. Chapters 5 to 7 contain zero floats, so no pagination
work is possible. "Rewriting only" holds.

---

## 2. The rationale is refuted by the thesis's own data

M5's load-bearing sentence is "most of the headline validity separation is a placement effect, not a
geometry-quality effect." Appendix H.4's own topic sentence commits the same error ("most of the
apparent gap turns out to be placement").

Validity by RMSD band, reproduced from the primary data:

| Arm | ≤2 Å | 2–5 Å | 5–10 Å | 10–20 Å | >20 Å |
|---|---|---|---|---|---|
| AutoDock + gnina | 97.49 | 99.00 | 99.31 | 99.13 | 98.69 |
| DiffDock + smina | 94.58 | 81.63 | 76.71 | 77.84 | 85.90 |
| EquiBind + gnina | 93.57 | 78.22 | 60.83 | 50.44 | 56.62 |

AutoDock is flat. Its near-native band is in fact its *worst* band. EquiBind collapses.

Stratified decomposition of the 36.95 pp AutoDock-versus-EquiBind gap:

| Component | Value | Share of gap |
|---|---|---|
| Within-stratum (geometry quality) | +37.59 pp | 101.7 % |
| Composition (placement mix) | −0.64 pp | −1.7 % |

Across 2-band, 5-band and symmetric specifications the composition term runs from −4.7 % to +0.1 %.
The mechanism is plain: near-native shares are 3.55 % for AutoDock against 5.48 % for EquiBind.
EquiBind produces the *larger* fraction of near-native poses, so reweighting cannot transfer
validity to AutoDock. Over the 95.5 % of poses that are misplaced the gap is **38.83 pp**, wider
than the pooled figure.

The narrowing is an **interaction**, not a **mediation**. Near-native poses are validity-clean for
every arm, so the physics-based advantage does not apply inside that stratum. That is a real and
publishable observation. It is not a statement that placement explains the gap.

---

## 3. The abstract limb is a category error

The four-point figure conditions the *optimised* trio. The abstract prints the *raw* trio.

| Trio | Unconditioned spread | Near-native spread |
|---|---|---|
| Optimised 98.9 / 84.5 / 62.0 | 36.9 pp | **3.9 pp** |
| Raw 99.6 / 24.4 / 2.9 (the abstract's) | 96.7 pp | **78.0 pp** |

Under the same conditioning the raw arms go 99.6→98.4, 24.4→50.3, 2.9→20.4. The gap does not
narrow to four points, it stays enormous. Appendix H.4 says so itself: "Physical validity
discriminates decisively between raw and minimised learned output and only marginally between the
minimised pipelines." The abstract also already qualifies its numbers in the next sentence.

Carrying the four-point figure into the abstract, as M5's priority item 3 asks, would attach a
minimised-arm statistic to raw-arm numbers.

---

## 4. Four factual errors in M5 as written

1. **"Zero times in the body" is true of the string, not the finding.** `body_main_short.tex:324`
   (§4.1.1, printed p.15) already states: "Conditioning on poses already within 2 Å of the crystal
   indicates that placement, not chemistry, explains most of the remaining gap." The body cites
   `\ref{validity-objective-coupling}` three times, not the two M5's narrative implies. What is
   missing is the number, not the idea. That sentence carries the same attribution error, and its
   "not chemistry" comparator is vacuous because the chemical-validity column is 0.0 in all sixteen
   cells.
2. **§5.2 is not unquantified.** Line 777 carries "Requiring validity removes at most four of 303
   complexes at any depth" alongside the concession and the cross-reference. That is a different
   H.4 result, so the gap on M5's axis is real, but "only" is wrong.
3. **§6.1 is neither raw nor bare.** Line 836 reads "at 98.9 %, 84.5 % and 62.0 %, although part of
   this separation is constitutive (Appendix H.4)". In house vocabulary "raw" means pre-minimisation
   and names 99.6 / 24.4 / 2.9. The hedge is nonetheless the wrong one, since "constitutive" names
   H.4's Vina-repulsion coupling, which the appendix says "is not a quantity this study estimates".
4. **The locators are wrong.** The evidence quote renders "spread of **Table 4**"; the source is
   `\ref{tab:results-pose-production}`, which resolves to **Table 1**, p.13. The decomposition table
   is Table 24 on p.125, inside **H.3 Multiplicity Families**, not H.4. H.4 starts on p.126 and the
   four-point sentence prints on p.126, not p.124.

---

## 5. M5 is a re-file, not a new finding

It is near-verbatim identical to M6 of `obsolete/EXAMINER_REVIEW_2026-09-01_full_msc_rubric.md:349`.
`diff` of the two heading lines returns one difference, "M6" versus "M5". It has never been refuted,
survived that round's adversarial pass, and was simply not actioned. Its August ancestor ("the
constitutive-advantage concession never leaves Appendix H", `obsolete/EXAMINER_REVIEW_2026-08-23.md:282`)
*was* filed, upheld, fixed and closed — and the fix is the very "partly constitutive" clause M5 now
calls unquantified. The examiner is not reading stale text: the §6.1 clause predates git tracking of
the short build, and the target lines 324, 777 and 836 are untouched by the working diff.

---

## 6. Recommendation

**Decline** M5's remedy as written. Do not carry "roughly four percentage points" into §5.2, §6.1 or
the abstract with the placement-effect framing.

**Act instead** on the defect M5 accidentally surfaced: Appendix H.4 ¶1 makes an attribution claim
its own data refute.

### Edit 1 (required) — `body_appendix_short.tex:1184`

Replace the opening clause:

```latex
Restricting the decomposition of Table~\ref{tab:appendix-pb-decomposition-full} to poses already
within 2\angstrom{} of the crystal ligand separates validity from placement. The spread is confined
to misplaced output rather than explained by it.
```

Insert after "…once the pose is in the right place." and before "A residual gap survives that
conditioning.":

```latex
That narrowing is not a compositional effect. Near-native poses are 3.6\% of AutoDock Vina + gnina
output against 5.5\% of unguided EquiBind + gnina output, so reweighting either arm to the other's
placement distribution moves its pooled validity by under one percentage point. Over the 95\% of
poses that are not near-native the spread is 38.8 percentage points, wider than the pooled figure.
Placement and validity therefore interact, and the conditioned figure marks where the physics-based
advantage does not apply rather than measuring how much of it is placement.
```

### Edit 2 (optional, defensible) — `body_main_short.tex:836`

Append after the existing clause. Do not print "four percentage points":

```latex
That separation is confined to misplaced output. Among poses already within 2\angstrom{} of the
crystal ligand the three rates are 97.5\%, 94.6\% and 93.6\%, on small and unequal subsets of 319,
1,715 and 498 poses.
```

Same pose-pooled unit as the 98.9 / 84.5 / 62.0 it qualifies, so the inferential-unit objection does
not bite. States the interaction, not a mediation. No "four" collision with §5.2.

### Edit 3 — decline

No insertion in §5.2 (its margin is the raw seventy-point one, and "four of 303 complexes" already
sits two sentences away), and none in the abstract or Kurzfassung.

### Also fix

`body_main_short.tex:324` carries the same attribution error and a vacuous "not chemistry"
comparator.

### Housekeeping

Pin 97.5 / 93.6 / 3.9 / 319 / 498 / 31 / 5 in a `prose: appendix_h4:` registry entry beside
`check_prose_appendix_h9`. `table_24` asserts the operands; the derived prose is unasserted. If M5
is logged, correct its locators.

---

## 7. Unverified

- The decomposition is pose-level and point-estimated. A paired analysis on the 69 complexes where
  both arms produce near-native poses gives +4.53 pp, and a cluster-adjusted interval on the 3.92 pp
  residual includes zero. The *direction* — composition share near zero, misplaced-stratum gap at or
  above the pooled gap — is stable across every specification run.
- The full `thesis_assertions.py` harness was not executed; `check_table_24`'s join was
  reimplemented and all sixteen cells reproduced.
- One of 58 agents died on an API safeguard error. Its dimension was covered by three other agents.
- `obsolete/` was not searched for whether H.4's attribution error has been filed before. It appears
  in neither the 09-01 nor the 09-03 round.
