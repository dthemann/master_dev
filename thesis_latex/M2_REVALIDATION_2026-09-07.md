# Re-validation of examiner finding M2 / E3-2 (2026-09-07)

**Finding under test** — "The blind whole-protein box places additional genuine copies of the correct binding
site inside the search volume for 143 of 308 benchmark entries, and the confound is stated once but never
quantified."
**Source** — `thesis_latex/EXAMINER_REVIEW_2026-09-06_multiagent.md:253-262`
**Method** — every pose of six whole-protein arms (53,651 poses, 303 complexes) and the rank-1 poses of the
old boxed AutoDock run were re-scored against every record of `<ID>_ligands.sdf` with the thesis's own
loaders, template reassignment and `symmetry_rmsd`, imported from `posebusters_pose_comparison.py`. The
reference-copy value reproduces the `rmsd` column of `per_pose_metrics.csv` to 4.7e-10 Å over all rows, and
the recomputed centroid distance reproduces `centroid_dist` to 4.5e-13 Å, so the engine is the thesis's
engine and the nearest-copy value is the only new quantity.
**Adjudicated against** — `posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report/per_pose_metrics.csv`
(canonical, Table 5 source) and `Data/PoseBuster Benchmark Set/`.
**Probe** — `Scripts/Analysis/nearest_copy_probe/` (script, rank-1 CSVs, per-depth summary, README).

---

## Verdict: CONFIRMED, with the flip counts understated

Every structural claim reproduces exactly. The two "on an alternate copy" counts reproduce exactly. The
rank-1 flip counts are four higher per arm than the report prints, and the direction claim (the correction
widens AutoDock's margin) is confirmed and holds at every depth. The proportionate remedy the report proposes
is the right one.

| # | Sub-claim | Verdict |
|---|---|---|
| 1 | 143 of 308 ids have >1 record in `*_ligands.sdf` | **CONFIRMED** — 111 with two, 7 three, 22 four, one each with five, six and twelve |
| 2 | The reference `*_ligand.sdf` is one of those records; the convention is "single instance" | **CONFIRMED** — reference centroid matches a record in all 308 |
| 3 | For all 143 an alternate copy lies >2 Å from the reference, inside the protein bounding box + 1 Å, within 5 Å of protein | **CONFIRMED** — 143/143 on all three |
| 4 | Median nearest-other-copy centroid separation 36.0 Å | **CONFIRMED** — 36.0 Å (min 4.9, max 119.5) |
| 5 | 139 of 143 alternate copies share ≥0.6 Jaccard of the 4.5 Å residue environment | **CONFIRMED** — 139/143, median Jaccard 0.96; the four exceptions are 7A9E_R4W, 7TUO_KL9, 7VKZ_NOJ, 7Z1Q_NIO |
| 6 | 138 of the 143 survive into the analysed 303 | **CONFIRMED** |
| 7 | Rank-1 poses within 4 Å of an alternate copy and nearer it than the reference: AutoDock* 58/138, DiffDock* 53/138 | **CONFIRMED** — exactly |
| 8 | Rank-1 near-native flips under nearest-copy scoring: AutoDock* 34 (+11.2 pp), DiffDock* 26 (+8.6 pp) | **UNDERSTATED** — 38 (+12.5 pp) and 30 (+9.9 pp). Four flips per arm have nearest-copy RMSD in (1.5, 2.0] Å, consistent with the verifier's symmetry handling being incomplete |
| 9 | "The correction moves AutoDock further ahead" | **CONFIRMED** — margin +8 → +17 complexes at rank-1, +31 → +34 at top-15 |
| 10 | The convention is stated only at `body_appendix_short.tex:664`; count, regime change and sensitivity arm are absent | **CONFIRMED** — grep on copies/copy/asymmetric unit/alternate/duplicat hits :664 and nothing relevant elsewhere |
| 11 | "The source benchmark's crystal-centred cube kept the copies outside the searched space" | **CONFIRMED and now measured** — see the regime contrast below |

AutoDock* is `autodock_mgltools_exh128_gnina` ranked by CNNaffinity, DiffDock* is `diffdock_smina`,
EquiBind* is `equibind_unguided_gnina` ranked by gnina affinity. Rank-1 PB-valid and near-native counts
reproduce Table 5 (111, 103, 55) before any change.

---

## The regime change, measured

The old boxed AutoDock run (`posebusters_results/benchmark/dock/`, 25 Å cube on the crystal ligand) and the
whole-protein run share the same DiffDock arm, so the two AutoDock arms can be compared on the same 303
complexes.

| Arm | Box | Rank-1 near-native, reference | Rank-1 near-native, nearest copy | Flips | On an alternate copy (of 138) |
|---|---|---|---|---|---|
| AutoDock raw Vina, boxed run | 25 Å crystal-centred cube | 175 | 175 | **0** | 2 |
| DiffDock-smina | none (blind) | 107 | 137 | 30 | 53 |
| AutoDock* exh128 + gnina | whole protein | 113 | 151 | **38** | 58 |
| AutoDock raw exh128 | whole protein | 95 | 131 | 36 | 59 |

Under the source benchmark's box the convention is inert for Vina, exactly as the finding states. Under the
whole-protein box AutoDock reaches alternate copies slightly more often than DiffDock (42.0 % against 38.4 %
of multi-copy rank-1 poses). A note in the project memory dated 2026-08-22 recorded the tool asymmetry the
other way round (DiffDock 11.7 %, AutoDock 0.9 %). That measurement was taken on the boxed arm and has been
annotated, not deleted, because it is the same regime effect seen from the other side.

---

## Sensitivity arm: recovery under both conventions

Existence gate over the top-d poses, n = 303, PB-valid and within 2 Å (the Table 5 endpoint). Near-native
without the validity gate is in the summary CSV and differs by at most two complexes.

| Arm | d | Reference copy | Nearest copy | Gain | Gain, pp |
|---|---|---|---|---|---|
| AutoDock* | 1 | 111 (36.6 %) | 149 (49.2 %) | +38 | +12.5 |
| DiffDock* | 1 | 103 (34.0 %) | 132 (43.6 %) | +29 | +9.6 |
| EquiBind* | 1 | 55 (18.2 %) | 57 (18.8 %) | +2 | +0.7 |
| AutoDock* | 15 | 198 (65.3 %) | 213 (70.3 %) | +15 | +5.0 |
| DiffDock* | 15 | 167 (55.1 %) | 179 (59.1 %) | +12 | +4.0 |
| EquiBind* | 15 | 80 (26.4 %) | 83 (27.4 %) | +3 | +1.0 |
| AutoDock* | 30 | 199 (65.7 %) | 214 (70.6 %) | +15 | +5.0 |
| DiffDock* | 30 | 169 (55.8 %) | 183 (60.4 %) | +14 | +4.6 |
| EquiBind* | 30 | 80 (26.4 %) | 83 (27.4 %) | +3 | +1.0 |

No complex loses recovery under the nearest-copy convention at any depth. Uni-Dock2 gains 14 at every depth
and raw DiffDock gains 12 / 9 / 9 on the valid-and-near endpoint.

AutoDock* against DiffDock*, exact McNemar on the discordant complexes:

| d | Convention | AutoDock* | DiffDock* | Margin | Discordant | p |
|---|---|---|---|---|---|---|
| 1 | reference | 111 | 103 | +8 (+2.6 pp) | 60 / 52 | 0.51 |
| 1 | nearest | 149 | 132 | +17 (+5.6 pp) | 73 / 56 | 0.16 |
| 15 | reference | 198 | 167 | +31 (+10.2 pp) | 82 / 51 | 0.009 |
| 15 | nearest | 213 | 179 | +34 (+11.2 pp) | 82 / 48 | 0.004 |
| 30 | reference | 199 | 169 | +30 (+9.9 pp) | 81 / 51 | 0.011 |
| 30 | nearest | 214 | 183 | +31 (+10.2 pp) | 79 / 48 | 0.008 |

The crystal-site reach at 4 Å moves the same way. Rank-1 reach goes 183 → 241 for AutoDock*, 172 → 225 for
DiffDock*, 137 → 144 for EquiBind*. The pooled three-tool oracle at `body_main_short.tex:476` goes from
98.0 % to 99.7 % of complexes (the sentence prints 97.7 %, one complex off on the present canonical file, a
pose-basis detail unrelated to M2 and not pursued here).

Single-copy stratum (165 complexes, where the confound cannot arise): AutoDock* 54.5 %, DiffDock* 50.9 %,
EquiBind* 32.7 % rank-1 near-native. Multi-copy stratum (138): 16.7 %, 16.7 %, 2.2 % against the reference
and 44.2 %, 38.4 %, 3.6 % against the nearest copy. The ordering holds in the unconfounded stratum.

Flipped poses sit 12.7 to 72.0 Å from the reference (median 38.4 Å) and 0.98 Å (median) from the copy they
recover. These are the same pocket on another chain, not near-misses of the reference.

---

## Options to close the finding

**A. Disclose the count and the regime change (text only).** Two sentences after `body_appendix_short.tex:664`
and one in the Limitations paragraph at `body_main_short.tex:868` that already lists denominator limits. Cost
under an hour. On its own it leaves the effect unsized, which is the part of M2 the examiners upheld.

**B. Add a nearest-copy sensitivity table (recommended).** One appendix table in "Supplementary Statistics
and Extended Results", beside the intention-to-treat check at `body_appendix_short.tex:1072`, printing the
first table above for the three headline arms at d = 1, 15 and 30, with one paragraph stating that no
complex loses, that AutoDock gains most, and that the margin widens. The primary endpoint stays on the
single deposited instance, which is the conservative reading. Inputs already exist. Cost: the probe, one
registered reproduction stage, three yaml assertions, one table, about half a day including a rebuild.

**C. Switch the primary endpoint to nearest-copy.** Rejected on two grounds. It credits a pose for a site
the study did not set out to predict, and it moves every ≤2 Å number in the thesis (Tables 1, 5, 7, the
appendix depth tables, the discordance tests, the recovery figures, the Abstract and Kurzfassung rates and
the 194 harness assertions) for a change the sensitivity arm shows alters no conclusion. The upstream tool does score this way (installed PoseBusters `redock.yml:307` sets `load_all:
True` on `mol_true` with the comment "important to set if there are multiple ligands that could be
correct"), which is worth one clause in the appendix as justification for offering the sensitivity arm,
not for adopting it.

**D. Report the single-copy stratum.** Cheap complement to B, needs no re-scoring, one row in the same
table. Shows the ordering survives where the confound is absent.

**E. Re-dock with the alternate copies masked or the other chains removed.** Rejected as disproportionate.
It is a new campaign on 138 receptors, it changes the receptor the three tools share, and B already answers
the question the examiners asked.

### Draft text for B (register and punctuation follow the thesis house style)

Appendix F, after the sentence at :664:

> That convention is not idle in this study. In 143 of the 308 entries the deposited file holds more than
> one copy of the ligand, 138 of them among the 303 analysed complexes, and the nearest other copy lies a
> median 36.0 Å from the reference. The source benchmark's 25 Å cube centred on the reference excluded those
> copies from the search, whereas the whole-protein box used here contains every one of them. A pose that
> recovers an alternate copy is therefore counted as a miss on every crystal-referenced endpoint.
> Appendix~\ref{supplementary-statistics} sizes that effect.

Appendix I, new paragraph and table beside the intention-to-treat check:

> Scoring each pose against the nearest of all deposited copies rather than the single reference instance
> raises rank-1 recovery by 38 complexes for AutoDock, 29 for DiffDock and 2 for EquiBind, and by 15, 12 and
> 3 at top-15. No complex loses recovery under that convention. The AutoDock margin over DiffDock widens
> from 8 to 17 complexes at rank-1 and from 31 to 34 at top-15, so the single-instance convention understates
> every tool and does not produce the reported ordering.

Limitations, one sentence in the denominators paragraph:

> Scoring against the single deposited ligand instance while searching a box that contains every other copy
> understates rank-1 recovery by 9.6 to 12.5 points and top-15 recovery by at most 5 points, and does so
> more for AutoDock than for DiffDock (Appendix~\ref{supplementary-statistics}).

### Reproduction chain for B

- Move `Scripts/Analysis/nearest_copy_probe/nearest_copy_rescoring.py` to `Scripts/Analysis/nearest_copy_sensitivity.py`,
  write the per-depth summary into `pose_comparison_report/nearest_copy_sensitivity.csv`.
- Register stage `bench_nearest_copy_sensitivity` in `_build_reproduction_notebook.py` (section "3. Benchmark
  analysis", determinism bitexact, cost cheap, needs `bench_pose_comparison`) and mirror it in REGENERATE.md.
- Add `nearest_copy` assertions to `thesis_expected_values.yaml`: multi_copy_ids 143, multi_copy_analysed 138,
  rank1 gains 38 / 29 / 2, top15 gains 15 / 12 / 3.

---

## Not part of M2, noticed in passing

- `body_main_short.tex:476` prints a pooled three-tool 4 Å reach of 97.7 % at median 0.34 Å and AutoDock
  88.1 % at 0.53 Å. The canonical file gives 98.0 % at 0.27 Å and 88.4 % at 0.44 Å over all thirty poses.
  One complex and the pose basis separate them. `pocket_localization.csv` is not its source (6 Å cutoff, raw arms), so
  trace the generating script before touching the sentence.
