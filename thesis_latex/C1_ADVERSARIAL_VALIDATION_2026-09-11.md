# C1 — adversarial validation of the charged-cost basis finding

**Date** 2026-09-11 · **Branch** `nearest-copy-endpoint` · **Verdict: C1 is CONFIRMED**, with three
corrections to how the examiner states it and one correction to its proposed remedy.

Method: seven independent lenses (steelman/defence, physics, accounting symmetry, consequence,
fix-critic, independent reproduction, history) plus adversarial verification of each, over the
primary artefacts — not the thesis's summary CSVs.

---

## 1. What reproduces

Everything the thesis prints reproduces exactly. This is not a computation error.

| quantity | thesis | reproduced |
|---|---|---|
| Vina search wall | 4.73 h | 4.7319 h |
| CPU-core hours | 151.42 | 151.4192 (= search wall × 32 exactly, all 303 complexes) |
| gnina GPU-hours | 10.27 | 10.2700 (8,976 cohort sidecars) |
| charged total | 15.00 h | 15.0017 h |
| rescoring share | 68.5 % | 68.46 % |
| Table 6 charged medians | 91.8 / 39.0 / 2.4 | 91.828 / 39.041 / 2.413 |
| Table 6 CPU-core / GPU | 1,233.3 / 2.30 / 159.4 · 83.6 / 39.1 / 1.4 | exact |
| Friedman (n=54) | χ²=87.1, W=0.81 | 87.11, 0.807, all pairs pass Holm |
| App Table 10 | 6.02 / 4.12 / 60.73 | 6.017 / 4.1189 / 60.7296 |
| App C.7 alternative | 5.45 h, 13.2 % | 5.4506 h, 13.19 % |

## 2. Why the finding holds

### 2.1 The physical ceiling (needs no concurrency model)

One 32-thread CPU, one GPU. The GPU was idle throughout the 4.7319 h serial Vina phase
(`max_workers: 1`). The gnina phase spans 0.7450 h end to end. The arithmetic maximum of
`CPU-core-s/32 + GPU-s` for this campaign is therefore

    4.7319 (CPU) + 0.7450 (CPU) + 0.7450 (GPU) = 6.22 h

Published: **15.0017 h — 2.41× above the ceiling.** A quantity named "device occupancy" reports
more occupancy than the devices could physically supply. Even the loosest bound, 2 × campaign
wall = 10.95 h, is exceeded.

### 2.2 The measurement

From 9,126 `optimized_gnina/*.provenance.json` sidecars (`created_at` + `optimizer_elapsed_time_s`):

    SUM of per-invocation elapsed   10.47 h
    UNION of the intervals           0.7254 h  (start-anchored) / 0.7187 h (end-anchored)
                                     0.7116 h  restricted to the 303-complex cohort
    campaign span                    0.7450 h
    inflation                       ~14.2 – 14.6 ×
    max concurrency                  16, exactly `optimize_workers`

`gpu_max_concurrent: 16` equals `optimize_workers: 16`, so the flock semaphore
(`_interprocess_gpu_lock`, `run_autodock.py:2262`) never blocks: the sixteen processes genuinely
shared one device. The thesis's 0.72 h was the 308-complex, end-anchored union (0.7187 h). The
cohort-consistent, start-anchored value is 0.7116 h, and the prose remedy on branch
`c1-cost-basis-remedy` prints it as 0.71 h (with 5.44 h and 13.1 %) at `body_appendix_short.tex:322`
and `:326` and in the Results and Conclusions.

### 2.3 The asymmetry is same-tool and it is 32×

`_equibind_time` (`docking_effort_comparison.py:803-855`) and `_autodock_effort` (`:420-448`)
treat the *same binary doing the same job* differently, purely on the `gnina_use_gpu` flag:

| | AutoDock arm | EquiBind arm |
|---|---|---|
| gnina time | 10.27 h, 8,976 poses | 20.21 h, 9,088 poses (**94.5 % of its whole CPU bill**) |
| routed to | `gpu_s` | `cpu_core_s` |
| divisor in the charged formula | **1** | **32** |
| measured per pose | **4.12 s** | 8.01 s |
| charged per pose | **4.12 s** | 0.250 s |

The same tool, measured **1.94× faster** on the AutoDock arm, is billed **16.5× more expensive**.
Amortised over their own concurrency the two passes cost the machine almost identically —
**EquiBind 0.632 h, AutoDock 0.642 h**.

### 2.4 The comparator really is serial, so "one basis" fails in substance

DiffDock, 303-complex cohort, reading `timestamp` as the completion time `run_diffdock.py:472`
actually writes: SUM 22.4609 h, UNION 22.3930 h, **inflation 1.003×**, mean concurrency 1.003.
The identical formula charges DiffDock at 1.00× its device occupancy and AutoDock at 14.4×.

### 2.5 The divisor 32 is three different operations wearing one name

Proved by algebraic identity over all 303 complexes:

    AutoDock   wall_s ≡ cpu_core_s/32 + gpu_s/16     (residual 1.1e-13)
    EquiBind   wall_s ≡ (cpu_core_s + gpu_s)/32      (residual 1.8e-15)
    DiffDock   wall_s ≡ cpu_core_s + gpu_s           (residual 1.1e-13)

So `/32` is an algebraic inverse of an assumed thread count for Vina's search, a real amortisation
across 32 concurrent processes for EquiBind, and a 32× discount on a documented-serial
**single-core** smina stage for DiffDock (`optimize_cpu: 1`). In every arm the GPU stream is
exempted from the concurrency division its own CPU stream receives.

---

## 3. Three corrections to how C1 is stated

1. **The A1 identity is definitional, not evidence.** `cpu/32 + gpu/16 == wall_s` to 1e-13 holds
   because the code *builds* `wall_s = dock + opt/16` and `cpu_core_s = dock × 32`
   (ratio exactly 32.000000 on all 303 rows). It cannot independently corroborate anything. The
   physics argument in §2.1–2.2 stands on its own and should carry the finding.

2. **"Reverses the conclusion" is true of one statistic, not all three.**

   | statistic | published | examiner's fix | symmetric fix |
   |---|---|---|---|
   | **median per qualifying pose** (Table 6, §6.3, abstracts) | 91.8 / 39.0 / 2.4 | 22.8 / 39.0 / 2.4 | 28.8 / 38.9 / 1.7 |
   | **paired, n=54** (Figure 15, "every pair passes Holm") | 115.6 / 17.7 / 1.9, Holm p=3.9e-05 | 21.2 / 17.7 / 1.9, **p=0.43** | 27.4 / 17.6 / 1.3, **p=0.87** |
   | **pooled ratio-of-totals** | 122 / 39 / 6 | 44 / 39 / 6 | **50 / 39 / 5** |

   The reversal is robust on the median — the statistic the headline claims actually use — under
   every correction. On the paired test the pair becomes **not separable**, with AutoDock still
   nominally *higher*; the examiner's 22.8-vs-39.0 compares two different cohorts (214 vs 183
   complexes). On a pooled statistic AutoDock stays highest. **The defensible replacement claim is
   "AutoDock and DiffDock are not separable per qualifying pose", not "DiffDock is dearer."**

3. **The "Mitigating" paragraph is unfair as written.** `body_appendix_short.tex:326` gives *two*
   grounds, not one: the directional one ("would flatter this arm") leads, but is followed by a
   concurrency-invariance argument. The correct rebuttal is that the invariance is bought by making
   the GPU term not a time — contradicting §3.2.7's own definition of the unit — and that it is
   asymmetric, because `/32` *does* credit EquiBind's inter-complex concurrency, which caveat four
   at `:774` concedes without ever connecting the two.

---

## 4. Defects the review did not find

1. **The gnina stage's CPU cost is charged to nobody.** `_autodock_effort:447` returns
   `dock_cpu_s` unchanged on the GPU branch. Every invocation ran at `optimize_cpu: 4`, and the
   config's own measurement (`…exh128_gnina.yaml:88`) is "~2.1 cores" each — 16 × 2.1 = 33.6 cores
   against 32 hardware threads, i.e. the pass **saturated the processor**. The omitted term is
   ~21.6 core-hours = 0.674 charged hours, almost exactly the stage's GPU occupancy of 0.72 h.
   `body_main_short.tex:772` states "The subsequent GPU rescoring pass adds device occupancy
   without adding CPU cost", and `:776` "at no CPU cost" — both false by the project's own
   measurement. Correcting this makes the surviving CPU-core finding **stronger**.

2. **DiffDock's charged total is below its own wall clock** — 22.5022 h charged against 23.7829 h
   elapsed. A device-occupancy figure smaller than the time the device was occupied is impossible
   on any reading of the term. Cause: its single-core smina stage is divided by 32.

3. **EquiBind carries the same defect in miniature.** Its `gpu_s` (0.1817 h) is an undivided sum
   over the same 32 workers by which its `cpu_core_s` is divided; charged 0.8502 h exceeds its own
   modelled wall 0.6742 h by 26 %. Its published 2.41 s should be 1.72 s.

4. **Applying one consistent divisor per resource stream reproduces the elapsed basis exactly**
   (5.3737 / 23.7829 / 0.6742 h). "Charged device occupancy" is therefore not a distinct currency:
   it is wall clock plus three arm-specific errors.

5. **The headline swings 9 h on a boolean.** Flipping `--autodock-gnina-gpu` changes no measured
   quantity (`wall_s` is identical on both branches) but re-routes the same 10.27 h log to
   `cpu_core_s × optimize_cpu / 32`, giving 6.02 charged hours at `cpu=4` — against the published
   15.00 h.

6. **Table 6 mixes estimators inside one row, and labels only one column.** The "Charged median"
   column is a per-complex median; the CPU-core and GPU columns are ratios of totals. Per-complex
   medians would be 556.1 / 1.2 / 54.2 CPU-core-s and 65.7 / 39.0 / 0.7 GPU-s. A reader comparing
   91.8 with 1,233.3 in the adjacent cell is comparing two different statistics of the same cohort.
   **This also affects the replacement headline the review proposes.**

7. **The thesis already publishes an allocation on which its headline reverses.**
   `body_main_short.tex:778`: "On this denominator AutoDock is cheaper than DiffDock in charged and
   GPU terms, with 173 versus 442 GPU-seconds per success" — then concludes "Neither denominator
   changes the per-qualifying-pose result." The opening sentence of the same paragraph, "The
   ordering depends on the denominator rather than on the currency", is false under C1: within the
   single charged currency the ordering flips on how the GPU term is charged.

8. **Self-contradiction across 450 lines.** Caveat three (`:774`) says each invocation ran "on a GPU
   that a single-pose evaluation leaves substantially idle". App C.7 (`:326`) says the sixteen-way
   sum "changes elapsed time without changing the work done". An idle device cannot be doing
   sixteen calls' worth of work.

9. **Disclosure regression in the uncommitted working tree.** At HEAD the Table 10 footnote read
   "…10.27 hours of gnina device occupancy over 8,976 poses **rather than the 0.72 hours the
   sixteen-way concurrent run actually took. It is a device-occupancy time and not a measured
   elapsed time on any row.**" That is deleted in the working tree (`body_appendix_short.tex:349`).
   The disclosure survives at `:322` and `:326`, but is one site weaker than when this area was
   last adjudicated.

10. **Root cause, named by no prior review.** The un-divided convention was *correct* for the run it
    was written against — the pre-exh128 gnina pass was proven serial (`max_workers: 1`, exclusive
    flock, zero interval overlap). The exh128 migration of 2026-08-29/30 moved to sixteen workers
    and **silently invalidated the accounting**; no round re-derived it. Same failure pattern as the
    known `exh128-migration-left-discussion-behind` case.

---

## 5. Retractions from this pass

Two of my own intermediate measurements did **not** survive verification and must not be quoted:

- **"Per-invocation elapsed does not inflate with concurrency (4.26 s at concurrency 1–2 vs 3.93 s
  at 12–16, ρ = −0.259)" — RETRACTED.** Anchoring concurrency at interval midpoints, **9,101 of
  9,126 invocations sit at concurrency 12–16 and only 3 below 6**: there is no low-concurrency
  stratum in this campaign. The sign of ρ flips with the anchor (+0.53 midpoint, +0.48 start,
  −0.45 end, −0.40 within-complex). *Whether the GPU was the contended resource is not answerable
  from this run.* The primary finding does not depend on it.
- **"DiffDock inflation 1.28×" — RETRACTED.** I read `timestamp` as a start time; `_build_log_row`
  is called with a completed result, so it is a completion time. Correct figure is **1.003×**,
  which strengthens the finding.

The AutoDock union is convention-independent (0.7187–0.7254 h; span 0.7450 h either way), and only
the start-anchored reading yields max concurrency = 16, matching `optimize_workers` exactly.

---

## 6. What to change

### 6.1 Code — fix `gpu_s` at source, not in the cost numerator

**Do not patch `_cost_numerator`.** The inflated `gpu_s` also writes `total_gpu_h` (10.27),
`gpu_s_per_generated` (4.12 → Table 10), `gpu_s_per_near2` (83.6 → Table 6) and all five
`resource_*` panels including Figure 16. Correcting only the charged numerator would leave Table 6's
own GPU column, Table 10 and Figure 16 asserting 83.6 / 4.12 / 10.27 beside a charged median
computed from the corrected values — one page contradicting itself. The bug is that `gpu_s` *means*
device-seconds for DiffDock and EquiBind and process-seconds for AutoDock; the fix is to make it
mean one thing.

`Scripts/Analysis/docking_effort_comparison.py::_autodock_effort:447`

```python
# now
return dock + opt_wall, opt, dock_cpu_s
# should be: charge the stage what it occupied, and stop dropping its CPU
return dock + opt_wall, opt_occupancy, dock_cpu_s + opt * max(int(optimizer_cpu), 1)
```

**The identical defect sits at `_collect_rows_orai:1262`** (`gpu_s = opt`, and unlike the benchmark
path the wall is not divided either: `wall = w + opt`). No thesis artefact depends on it, but it
must move in the same commit or the two dataset paths diverge.

**Gate it behind `--gpu-accounting {process,device}` defaulting to `process`**, so today's behaviour
is bit-preserved and the byte-parity regression gate still holds — that gate is currently **23 of 23
files identical**, tighter than the 21/23 recorded in `REGENERATE_NOTES.md` on 2026-09-01.

**Decide the divisor first; it is a scientific choice, and all three are defensible:**
`/16` configured workers (0.642 h), `/14.6` measured mean concurrency (0.705 h), or the per-complex
union of each complex's own sidecar intervals (0.7116 h over the 303 cohort). The union is the
measured quantity and needs no flag. Complexes ran strictly sequentially (0 of 303 spans overlap
another), so the sum of the per-complex unions equals the global union to 1e-6: the per-complex
union is an exact additive decomposition of device occupancy and remains a true occupancy figure
under any cohort restriction. Say which one is used. (An earlier draft claimed the per-complex
unions "sum to more than" the global union and gave 0.755 h; both were wrong.)

**Assertions that actually hold.** The tempting one, charged total equals `total_wall_h`, is
FALSE once the gnina stage's CPU is billed; it holds only for `/16` with the CPU still dropped. Use
instead: (A) per complex, the sidecar elapsed sum agrees with the CSV sum to within 0.01 s × n
(measured max 0.0435 s); (B) per complex, longest single interval ≤ union ≤ summed process time,
and union ≥ summed / `optimizer_workers`; (C) Σ per-complex unions == global union (measured
difference 0.0 s); (D) the corrected charged total (≈6.1 h with CPU billed) stays below twice the
campaign span (2 × 5.48 h), the ceiling the published 15.00 h exceeded by 2.41×.

For symmetry the same pass should stop dividing DiffDock's single-core smina stage by 32 and should
divide EquiBind's `gpu_s` by its own `n_parallel_workers`. Both corrections move the comparators
*against* AutoDock (DiffDock 39.04 → 39.84 s, EquiBind 2.41 → 1.72 s), so the package is not
one-sided.

### 6.2 Data, figures and generated provenance

Re-run `bench_effort_charged` and `bench_effort_elapsed` as `Scripts/Analysis/REGENERATE.md` §3.1b
gives them, then copy the regenerated panels to `media/media/image24.png` (Figure 15) and
`image25.png` (Figure 16). Expect 18 of 23 files to move in the charged tree and 9 of 23 in the
elapsed tree — all of the latter GPU-bearing; the five `wall_*` panels and both stats JSONs must be
unchanged. Anything outside those sets moving means the change leaked.

Figure 15 prints the basis on the panel itself and its AutoDock–DiffDock bracket goes `***` → `ns`.
`thesis_assertions.check_figures` asserts byte parity against `thesis_expected_values.yaml`, so
skipping the copy fails loudly rather than silently.

Also update: `thesis_expected_values.yaml` (`table_6`, `table_10`), `validate_regeneration.py`
Case '15', `_build_values_registry_main.py` (the charged recipes), and rebuild the two generated
files — `_build_reproduction_notebook.py` then `regenerate_guide.py`. Hand-edit only
`REGENERATE_NOTES.md` §3.1b, whose claim that "the charged basis is not sensitive to
`--autodock-optimizer-workers`" becomes the defect's own epitaph. **While there, fix a pre-existing
error in the same paragraph:** `REGENERATE_NOTES.md:686` records "Friedman chi2 = 85.2, p = 3.2e-19,
Kendall W = 0.80, n = 53 paired" against the committed 87.1 / 1.2e-19 / 0.81 / n = 54 — stale since
the nearest-copy promotion.

**Blast radius:** `--basis charged` appears in exactly one regeneration stage, and no Orai1 artefact
depends on the twin defect. No new docking is required. Build with `pdflatex` twice, **not
`latexmk`** — it deletes the tracked `.bbl`.

### 6.3 Prose — fourteen sites

| file:line | what breaks |
|---|---|
| `body_main_short.tex:176` | Methods definition — the GPU term is not "a time rather than a resource-second" |
| `:728` | 10.27 "GPU-hours of device occupancy"; "roughly twice as expensive as the search" (true ratio ≈ 0.15×, **inverted**); 15.00 h; 68.5 % |
| `:733` | Figure 15 caption |
| `:755` | Table 6 AutoDock row: 91.8, IQR, 1,233.3, 83.6 |
| `:762` | "most expensive route to a qualifying pose"; χ²=87.1 / W=0.81 / ratio 0.15; **"every pair passes Holm" becomes false** |
| `:764` | "5.3, 5.3 and 91.8"; multiplier 17.3 |
| `:772` | "most expensive on every axis"; "adds device occupancy without adding CPU cost"; "highest GPU cost per qualifying pose" |
| `:774` | caveat one (4.73 of 15.00; 10.27 on a shared GPU) |
| `:776` | "most expensive in charged, CPU and GPU terms"; "at no CPU cost" |
| `:778` | "ordering depends on the denominator rather than on the currency"; "Neither denominator changes the per-qualifying-pose result" |
| `:835` | §5.6 "charged on one basis … comparison of like with like" |
| `:842` | §5.7 "DiffDock costs less per qualifying pose" — becomes unsupported (the pair is not separable) |
| `:864` | §6.3 medians; "comparable in basis as well as in scope" |
| `:866` | §6.3 "most expensive … in charged, CPU-core and GPU seconds"; "Removing gnina puts its charged and GPU cost below DiffDock" |
| `body_appendix_short.tex:322` | "consuming 10.27 GPU-hours of device occupancy in 0.72 wall-hours" |
| `:326` | the whole App C.7 justification paragraph — **rewrite, do not patch** (see below) |
| `:345` | Table 10 AutoDock row 6.02 / 4.12 → ≈2.15 / 0.26 (60.73 unchanged) |
| `:349` | restore the footnote disclosure deleted from HEAD |
| `Thesis_short.tex:63` (EN), `:99–:100` (DE) | "On a single device-occupancy basis … most expensive per qualifying pose". **Two traps:** the DE sentence *wraps*, with "In einheitlicher Belegungsrechnung" on `:99` and "am teuersten" on `:100` — fixing one line alone leaves the basis claim standing; and both abstract pages are page-locked at ~700 pt, so the replacement must be length-neutral or the Keywords/Schlagworte block orphans |

**Cleared as SAFE** (checked, not assumed): Appendix Table 8's "Search h" column and Figures 30/31
(`image46.png`, `image47.png`) are measured Vina search time only, with no GPU term. RQ3 itself
(`:13`) carries no number. Table 6's DiffDock, EquiBind and "Vina search only" rows are unchanged.

**App C.7 `:326` needs rewriting rather than patching.** It is the thesis's own defence of the
defect, and its central premise — that charging the pass at its union "would flatter this arm
against the pipelines it is compared with" — is refutable on the same cohort: DiffDock's GPU
intervals give SUM 22.4609 h against UNION 22.3930 h, so the comparator is *already* charged at its
occupancy. Occupancy is the symmetric treatment, not a concession.

### 6.4 What the headline should become

Two of the three currencies in `:866` flip; one survives.

- **charged** — AutoDock and DiffDock are **not separable** (Holm p = 0.43–0.87 depending on the
  correction). Say that, rather than swapping the winner.
- **GPU-seconds** — flips hard: 83.6 s → ~5.2 s per qualifying pose, from highest of three to
  second-lowest, against DiffDock's 39.1. `:772`'s attribution of this to low qualifying yield is an
  artefact of the undivided worker count.
- **CPU-core-seconds** — **survives, and gets stronger** once the omitted ~21.6 gnina core-hours are
  billed. AutoDock remains the most processor-intensive engine by an order of magnitude on either
  estimator (1,233 vs 2.3 vs 159 pooled; 556 vs 1.2 vs 54 as per-complex medians). Promote it — but
  label the estimator, and state that AutoDock's CPU-core-seconds are search wall × 32 *requested*
  threads rather than a sampled occupancy.

### 6.5 Cheapest defensible alternative

If the schedule will not take a regeneration: **adopt the thesis's own 5.45 h figure**, already
computed and printed at `body_appendix_short.tex:326`, and state plainly that the ordering depends
on whether the concurrent GPU stage is charged at its summed process-elapsed time or at its measured
device occupancy. That converts a false claim into a disclosed limitation without touching the data
— but it does not repair Figure 15, which prints the contested basis on the panel.

---

## 7. Is this double jeopardy?

No. The reversal has been in the project's records since 2026-08-30, but every prior round (08-30
twelve-refuter, 08-31 external, both 09-01 examiners, 09-03, 09-06, 09-09) litigated a *different*
proposition — that the accounting was non-commensurate **across tools**. That was closed on 09-01 by
unifying onto the charged basis. C1 attacks the unified basis **itself**, a target that did not exist
before 09-01, using the Methods sentence that the 09-01 repair itself wrote. No round has
adjudicated it.
