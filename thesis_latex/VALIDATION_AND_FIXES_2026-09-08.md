# Double-validated number sweep of the promoted thesis, and the eight fixes it produced

**Scope.** Every number a reader sees in the short build (`Thesis_short.tex`, `body_main_short.tex`,
`body_appendix_short.tex`) and every shipped figure, validated after the nearest-copy endpoint was promoted
into the canonical locations on 2026-09-08 (tag `promoted-nearest-copy-20260908`). The question asked of every
value was not only "does it match a sidecar" but "does it match the *promoted* data, and is it the right
quantity for the sentence that prints it".

**Standing rule.** Nothing in this sweep was fixed on one confirmation. Every fix below was established on at
least two independent routes, and the last two were additionally put to three adversarial agents whose brief
was to refute them.

---

## 1. How the validation was run

### 1.1 The two routes

Every text chunk was validated twice, by two agents that were forbidden to share work:

| route | what it was allowed to use | what it establishes |
|---|---|---|
| **A, sidecar** | `thesis_expected_values.yaml`, the stage sidecars, `REGENERATE.md` | that the printed value equals what the pipeline wrote |
| **B, recomputation** | `per_pose_metrics.csv` and the other primary tables, `stats_utils`, own scripts | that the pipeline's value is itself correct |

Route B was explicitly barred from reading the yaml or the txt sidecars except to learn what a number means,
so a value that both routes confirm has been checked against the data twice by two different derivations.
Each route returned its confirmed values as `L<line>:<value>` strings, which were intersected to separate
values confirmed twice from values confirmed once.

### 1.2 The gates and orders every recomputation had to use

Mis-stating any of these silently changes a count, so they were fixed in the brief given to every agent:

- near-native = `rmsd <= 2`, where `rmsd` is now the distance to the **nearest deposited copy**
- double gate = `pb_valid AND rmsd <= 2` (the Table 2 endpoint)
- triple gate = double gate `AND bestfit_rmsd <= 1` (the variant-selection endpoint)
- Form = `bestfit_rmsd <= 1 AND rmsd < 1000` (the exploded-pose guard)
- depth-`d` recovery = **existence** of a qualifying pose within the tool's own top-`d` order:
  AutoDock by `optimized_rank` on `autodock_mgltools_exh128_gnina`, DiffDock by `rank` on `diffdock_smina`,
  EquiBind by **ascending `gnina_affinity`** on `equibind_unguided_gnina`, because EquiBind's `rank`
  column is a 999 sentinel
- cohort = the 303 ids in `analysed_cohort_ids.txt`

### 1.3 The stale-value trap the sweep was built to catch

The promotion moved the primary endpoint from the single deposited ligand instance to the nearest deposited
copy. A number that still carried a retired instance value would look perfectly self-consistent, so
`Scripts/Analysis/nearest_copy_program/harness_signatures/harness_signature_A_oldvalues_newtree.csv` was given
to every agent: it lists, per harness key, the OLD instance value, the NEW nearest value and the thesis
location. Any printed value matching the OLD column where the two differ was to be reported as
`STALE-INSTANCE`.

**Result: not one such value was found anywhere in the validated text.** The Abstract and the Kurzfassung were
additionally checked against each other and agree number for number.

### 1.4 Coverage achieved

| part | independent routes | outcome |
|---|---|---|
| Abstract and Kurzfassung | 3 | clean |
| Methods, Results 4.1 / 4.2 / 4.3, Discussion, Conclusions, Limitations | 3 to 4 per section | 8 fixes below |
| 39 shipped figures | 2 (provenance and content) | 37 confirmed, 2 declared hand renders unverifiable, 0 stale, 0 inconsistent |
| Appendix, 1,454 lines in 7 chunks | 2 plus refuters | running at the time of writing |

The 14 figures regenerated for the switch are exactly `image3, 4, 5, 6, 7, 8, 9, 11, 24, 25, 44, 45, 46, 47`
(Figures 2, 3, 4, 5, 6, 38, 7, 39, 15, 16, 18, 19, 20, 21). Each differs from its pre-promotion twin, matches
`figure_copy_manifest.csv` md5 for md5, and matches the file the harness pairs it with. The other 25 are
byte-identical to the backup and each is convention-independent.

---

## 2. What was fixed

### 2.1 Round one, six values (commit `2446c658`)

| # | location | was | now | route to the correction |
|---|---|---|---|---|
| 1 | main :333 | rank-biserial `0.994--1.000` | `0.998--1.000` | `09f_pbvalid_yield_report.txt` prints 0.998, 0.998, 1.000, 1.000; both routes recomputed the same four Wilcoxon contrasts from per-complex yields. No estimator variant produces 0.994 |
| 2 | main Table 4, six cells | 0.641 / 0.671 / 0.632 / 0.508 / 0.459 / 0.499 | 0.640 / 0.670 / 0.631 / 0.507 / 0.458 / 0.498 | the table had been rounded from the 4-dp CSV, which double-rounds. Recomputing the means from `recovery_detail_per_pose.csv` gives 0.640476, 0.670448, 0.631477, 0.507464, 0.458474, 0.498483. The other 39 cells are exact. The appendix twin at :543 and the yaml pins were corrected with them |
| 3 | main :733 footnote | 0.74 h, 5.47 h, 13.5 %, "fourteen times" | 0.72 h, 5.45 h, 13.2 %, "fifteen times" | all four derived from a 2,677.7 s span that `REGENERATE_NOTES.md` itself calls retired. `exhaustiveness_arm_status.json` records `gnina_optimizer_wall_clock_s` 2,587.3 s and a serial sum of 37,693.7 s, so the ratio is 14.6 |
| 4 | main :779 | "pipeline success differs ... Q = 61.1" | "the share of complexes with at least one PoseBusters-valid pose differs ..." | `effort_stats.json` tests the PB-valid indicator (301, 300, 266 complexes). The qualifying-success indicator the sentence implied would give Q = 132.7, p = 1.6e-29 |
| 5 | main :648 | "p = 0.667", "to above 0.2" | "Holm p = 0.667", "to about 0.2" | the sidecar prints raw 0.333 and Holm 0.667; the neighbouring clauses already said "Holm p". The ligand-unit value recomputes to 0.1997 |
| 6 | main :429 | "three and four points of rank-1 recovery" | "... rank-1 near-native rate ... before the validity requirement is applied" | the 3.0 and 3.6 pp come from the rerank summaries on the **ungated** 2 Å rate, whereas every other rank-1 figure in the section is validity-gated |

### 2.2 Round two, two values (this commit)

Both were put to three adversarial agents, each on a different route and each instructed to **refute** the
defect and to state the strongest case against its own verdict. All three angles and the judge confirmed both.

#### Fix 7, main :433, "all 303" to "301 of the 303"

```
- AutoDock returns a valid pose for all 303 complexes, but only 82.2\% remain after the rank-1
  crystal-site trim.
+ AutoDock returns a valid pose for 301 of the 303 complexes, but only 82.2\% of the cohort remain
  after the rank-1 validity and crystal-site trim.
```

**Route 1, primary data.** Restricted to the 303-id cohort, `autodock_mgltools_exh128_gnina` has at least one
PoseBusters-valid pose for **301** complexes and a valid rank-1 pose for 298. The raw arm
`autodock_mgltools_exh128` reaches 303. Reproducing the trim (`pb_valid AND centroid_dist <= 8`) gives 249 of
303 at rank-1 for the gnina arm, that is 82.18 %, against 237 (78.2 %) for the raw arm.

**Route 2, generator and sidecar.** `filmstrip_stats__per_tool_depth.csv` carries the row
`AutoDock, 1, 249, 249, 303, 82.178` under the header `n_poses, n_valid_complexes, universe_complexes,
coverage_pct`, and the panel sidecar states the tool collapse `AutoDock = autodock_mgltools_exh128_gnina`.
So the paragraph is unambiguously the gnina arm, and 82.2 % is 249/303 over the benchmark universe.

**Route 3, document-internal.** Table 1 prints `AutoDock (raw)` 303 and `AutoDock + gnina` **301** on the same
page range. Line :332 attributes "every one of the 303 complexes has at least one valid pose" explicitly to
the raw search. The parallel sentence at :437, "EquiBind produces at least one valid pose for 266 of 303
complexes", uses the same at-least-one-valid convention with the gnina arm's own value. Line :433 was the only
place attaching 303 valid complexes to the carried-forward AutoDock arm, and it contradicted Table 1.

**Why the percentage was not touched.** The strongest counterargument raised was that 82.2 % is literally
249/303, so replacing the denominator naively would invite a reader to compute 249/301 = 82.7 %. The fix
therefore keeps 82.2 % and names "the cohort" as its base. The same edit repairs a second, smaller error in the
clause: the rank-1 drop from 303 to 249 is a two-step (303 poses, 298 valid, 249 also within 8 Å), not a pure
site trim, which is why "validity and" was added.

**Nature of the defect.** Pre-existing and convention-independent. PoseBusters validity does not move with the
nearest-copy switch, so this was not a promotion miss. The two complexes the gnina arm loses are `7NUT_GLP` and
`7OZC_G6S`, which fail internal-energy checks on all thirty poses after rescoring.

#### Fix 8, main :694, "fifteen" to "sixteen ... fifteen of which"

```
- PandaMap assigns residue-level contacts to fifteen interaction classes.
+ PandaMap assigns residue-level contacts to sixteen interaction classes, fifteen of which occur on
  Orai1 because the receptor carries no metal ion.
```

**Route 1, the third-party library.** The installed PandaMap 4.1.0 initialises its interaction dictionary with
exactly **sixteen** keys, twice and identically, at `pandamap/core.py:717` and `:797`, and every one is a live
detector; metal coordination is emitted at `core.py:1072-1089`. The project wrapper says so in its own words at
`Scripts/Analysis/run_pandamap.py:83`, "The 16 PandaMap interaction types".

**Route 2, the data.** The calibration-benchmark panel exhibits all sixteen classes over 232,024 interaction
rows. The Orai benchmark panel exhibits fifteen over 384,552 rows, with `metal_coordination` the sole absentee;
the Orai JKU panel exhibits ten. The mechanism was confirmed rather than assumed: the four Orai receptor files
contain no HETATM record at all, so no metal ion exists for the class to fire on.

**Route 3, the document.** The sentence is the reader's only statement of the profiler's vocabulary, is in
generic present tense with no "here" or "on this receptor" qualifier, and cross-references an appendix section
that enumerates all sixteen classes and then devotes a paragraph to metal coordination as a class PandaMap
emitted on the calibration panel. As written it undercounted the tool by one and contradicted its own appendix.

**Why both numbers are printed.** Fifteen is not simply wrong, it is the number realised on Orai1, and the
appendix's false-discovery-rate family over the Orai types depends on it. Swapping fifteen for sixteen would
have broken that; naming both keeps every downstream statement true.

---

## 3. Values that were challenged and survived

| location | challenge | why the printed value stands |
|---|---|---|
| main :516, `Wilcoxon p = 0.064` | claimed to be 0.063 from unrounded F1 | recomputing F1 exactly from the integer contingency counts, `2tp/(2tp+fp+fn)`, over the 280 shared rank-1 complexes gives p = 0.06379, which prints 0.064. Worth knowing for a defence: the 26 exact ties drive it, and the Pratt and z-split variants give 0.043 |
| main :871, "sixteenfold" | claimed to be tenfold | 39.0/2.4 = 16.25 follows from the two medians the same sentence prints. The cost sidecar's **paired** estimator on the same tier gives 10.45 with a 95 % interval of 7.14 to 14.20. Both are defensible and they are different estimators; naming which one is used would be an improvement, not a correction |

## 4. Values that cannot be validated from this repository

These are transcribed literature figures with no repository source. They are unverifiable by construction and
must be checked against the papers, not against the pipeline.

| location | value | source |
|---|---|---|
| main :614 | 73 % Orai1 transmembrane sequence identity | Hou et al. 2012, `ref037` |
| main :822 | 22.1 % and 67.3 % TEMPL success rates | Fülöp, Šícho and Dehaen 2025, `ref065` |

---

## 5. Verification after the fixes

| check | result |
|---|---|
| `Scripts/Analysis/thesis_assertions.py` | **1398 / 1398 checks reproduce** |
| `latexmk -pdf Thesis_short.tex` | 0 errors, 145 pages |
| PDF probes | all new strings present, all superseded strings absent |
| Kurzfassung page | keywords still on page 4 |

Neither round-two fix changes a computed quantity, so no figure, table or statistics file needed regenerating
and no harness pin moved.

## 6. Known gap this sweep exposed

No harness assertion covers prose clauses. That is precisely how a raw-arm count survived inside a
dominant-arm paragraph at :433. An equivalent raw-versus-gnina swap could recur silently elsewhere. Closing it
would mean adding prose-level checks for the coverage rows, which is recorded here as an open item rather than
done.
