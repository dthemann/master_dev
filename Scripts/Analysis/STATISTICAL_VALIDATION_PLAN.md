# Statistical validation plan — Benchmark & Orai × JKU analyses

Audit of every figure-producing analysis script for the **Benchmark Set** (has crystal
reference → RMSD-to-native) and the **Orai × JKU set** (no crystal → geometry / validity /
region only). For each analysis: the claim the figure makes, the data structure, the test
that would establish it, and whether a test already exists.

Scope note: this catalogs what *would* validate the presented results and flags where nothing
does yet. It does not (yet) implement anything.

---

## 0. Five cross-cutting principles (apply to almost every figure)

These matter more than any single per-figure test. Getting them right is what actually makes
the results defensible; getting them wrong makes even a "significant" p-value meaningless.

1. **Unit of analysis / pseudoreplication — the #1 issue.**
   Each tool emits ~30 (DiffDock/AutoDock) to ~135 (EquiBind) *correlated* poses per complex.
   Bar charts that pool raw poses imply n = tens of thousands when the true independent n is
   **~303 complexes (Benchmark)** or **12 (frame,ligand) pairs (Orai × JKU, 3 ligands × 4 frames)**.
   Any test must aggregate to the complex / (frame,ligand) level first, **or** use a
   cluster (receptor / complex) bootstrap. Everything below assumes the correct unit.

2. **The design is paired, so use paired tests.**
   All tools dock the *same* complexes → cross-tool comparisons are paired at the complex level
   (with informative missingness where a tool yields no valid/near-native pose). Paired tests
   (McNemar / Cochran's Q for proportions; Wilcoxon / Friedman for continuous) are both correct
   and more powerful than independent-sample tests. Use independent-sample tests only for the
   genuinely disjoint splits (Ro5 pass/fail, best-N vs worst-N, flexibility tertiles).

3. **Multiplicity.** Comparisons run across 3 tools × 20 PoseBusters checks × many
   thresholds/ranks/descriptors. Control the error rate: **Holm** (FWER) for a small pairwise
   family, **Benjamini–Hochberg** (`scipy.stats.false_discovery_control`) across large families
   (the 20 checks, the 14 descriptors, per-residue hotspots).

4. **Report effect size + CI, not just p.** At n ≈ 303 a 1-pp difference can be "significant"
   yet irrelevant. Pair every test with an effect size (risk difference, rank-biserial,
   Cliff's δ, Kendall's W) and a CI (Wilson for proportions, bootstrap/Hodges–Lehmann for
   continuous). "Validity of the presented results" is really a question about effect size + CI.

5. **Orai × JKU small-n honesty.** With **3 ligands**, formal inference is severely limited.
   Use exact / permutation tests on the 12 (frame,ligand) units, label conclusions exploratory,
   and never manufacture significance by treating hundreds of poses as independent.

**Environment:** `vina` env has scipy 1.16 / numpy / pandas / sklearn only —
**no statsmodels / scikit-posthocs / pingouin**. Every test below is buildable from scipy
primitives + the hand-rolled helpers (`_holm`, McNemar via `binomtest`, Cochran's Q, Dunn,
Jonckheere–Terpstra, permutation, Wilson CI, Cliff's δ) that **already exist** in
`pose_cluster_crystal_pocket_report.py` and `pocket_comparison_report.py`. Recommendation:
factor those into a shared `stats_utils.py` and import everywhere.

---

## 1. What already has a statistical test (the good news)

| Analysis | Figure | Test in place |
|---|---|---|
| Pose clustering (Benchmark) | `crystal_cluster_homogeneity.png` (`--stats`) | Friedman+Kendall W+Wilcoxon/Holm (A); permutation (B); Kruskal–Wallis+Dunn/Holm (C); Jonckheere–Terpstra+Spearman (D) |
| Pose clustering (Benchmark) | `cluster_quality_metrics.png` panel F | Cochran's Q + pairwise exact McNemar/Holm + Spearman enrichment |
| Pose comparison (Benchmark) | `20d…depth_filmstrip_pbvalid__stats.png` | Friedman + Wilcoxon/Holm + Cliff's δ + cluster-bootstrap CI + Spearman (only panel C draws a star) |
| gnina re-rank (Benchmark) | `ranking_stability_*.png` panel B | per-complex Spearman ρ distribution (Kendall τ computed, CSV-only) |
| gnina re-rank (Benchmark) | `chemotype_benefit_*.png` panels B, C | Spearman r + p (Mann–Whitney U computed, CSV-only) |
| Pocket comparison (Benchmark) | `crystal_success_vs_{topn,threshold}.png` | Wilson 95% CI bands |
| Pocket comparison (Benchmark) | `pocket_attr_ranking_*.png` | Spearman ρ (p not shown) |
| PandaMap interactions (both) | `14_rmsd_vs_recovery.png` | per-method Spearman ρ (in legend) |
| Ligand/receptor difficulty (Benchmark) | `*_difficulty_delta_heatmap.png` | Cliff's δ effect size (no p / CI) |

Everything else in the suite has **no inferential test**. The three biggest gaps are the
**validity report (both datasets, all figures)**, the **~50 pose-comparison accuracy figures**,
and **every Orai × JKU figure**.

---

## 2. PoseBusters validity report — `posebusters_validity_report.py`  (BOTH datasets, currently 0 tests)

Runs on Benchmark and Orai × JKU. Boolean `pb_valid` per pose, grouped by tool/variant. This is
arguably the headline result of the whole study and carries no statistics.

| Figure | Claim it makes | Structure | Recommended test |
|---|---|---|---|
| `01_per_tool_validity`, `03_grouped_bars`, `04_per_pair_distribution` | "Tool X yields more valid poses" | per-complex valid-pose **count**, 3 paired tools | **Friedman + Kendall's W**, **Wilcoxon** post-hoc (Holm), rank-biserial — identical to homogeneity panel A; reuse it. For "≥1 valid pose per complex" (binary): **Cochran's Q + McNemar/Holm** + Wilson CIs |
| `05_per_check_passrate` | "EquiBind fails bond-angle/length checks; others don't" | 20 checks × tools, paired proportions per complex | **per-check McNemar / Cochran's Q**, **BH across the 20 checks**; report risk difference. Directly substantiates the geometry-failure claim |
| `06_equibind_variant_validity`, `07_diffdock_variant_validity`, `07b_variant_matrix` | "gnina refinement beats smina / raw" | same poses minimized → **paired** binary | **Cochran's Q + McNemar/Holm**. (A paired smina-vs-gnina test exists elsewhere in the project — p≈6e-17 — but is *not* on these figures; port it here) |
| `08a/08b_validity_vs_rmsd` (Benchmark only) | "valid poses are also near-native" | cumulative valid-≤τ curves per tool | pointwise **Wilson CI** band; compare tools at fixed τ via **McNemar**, or paired **permutation** on curve area |
| `09_failure_waterfall` | attrition per check | per-method Pareto (descriptive) | low priority; optional Wilson CI on each drop |

**Orai × JKU caveat:** correct unit = 12 (frame,ligand) pairs, not the pooled poses. Cochran's Q /
exact McNemar are computable but underpowered — report exact CIs and state n = 12 explicitly.

---

## 3. Pose comparison vs crystal — `posebusters_pose_comparison.py`  (Benchmark only; ~50 figs, only `20d…__stats` tested)

Every figure is RMSD-to-crystal driven (no-op on Orai × JKU). Same 303 complexes → paired.

- **Oracle RMSD distributions** — `01_cdf`, `02_boxplot`, `14_twist_turn` (6 geometry measures):
  3 paired continuous → **Friedman + Kendall's W**, **Wilcoxon** post-hoc (Holm) on per-complex
  values + Hodges–Lehmann median-difference CI. For CDFs, test the *paired per-complex difference*
  with Wilcoxon — **not** a 2-sample KS (the samples are paired, not independent). For `14`, BH
  across the 6 measures; top-1 vs oracle is a within-tool paired Wilcoxon.
- **Success-rate bars** — `03`, `09a/b/c`, `10`, `12`, `13`, `19`, `20c(A)`: proportion of complexes
  ≤2 Å (± PB-valid), same complexes → **McNemar / Cochran's Q** (Holm/BH) + Wilson CIs.
  Oracle-vs-top1 and raw-vs-optimized are within-tool paired McNemar.
- **`11_vs_posebusters_paper`**: our rate vs a fixed external % → one-sample **binomial test**
  (`binomtest`) + Wilson CI. Caveat inline that complex sets differ (not strictly comparable).
- **Rank / threshold curves** — `04_rank_success`, `05_cumulative`, `15b`, `16a`, `18_*`:
  success vs rank/threshold → pointwise **Wilson CI**; compare tools at chosen depths via
  **McNemar**; whole-curve via paired **permutation on AUC**.
- **`15c/15d` optimization benefit**: Δpp raw→refined per rank, paired → **McNemar** per rank bin
  (BH) or bootstrap CI on Δpp.
- **`06`, `08` scatter (top1-vs-oracle, cross-tool)**: paired continuous → **Wilcoxon** on the
  paired difference + concordance (Spearman/Lin's CCC).
- **`07`, `15` oracle-rank histograms**: descriptive; optionally test rank distribution vs a
  "confidence is uninformative" **uniform null** via a permutation / multinomial test.
- **20-series form fidelity** (many panels): the `__stats` variant is the template — extend the
  same Friedman + Wilcoxon/Holm + Cliff's δ + bootstrap-CI machinery to the untested siblings
  (`20_form_fidelity_AB/CD`, `20b`, the `rank_depth_*` families). `20e` clustering is a
  no-cluster diagnostic → a **permutation test on ARI / silhouette vs label-shuffled null**
  makes the "no real clusters" claim rigorous.
- **`17*` PoseBusters waterfalls**: descriptive cascades; the meaningful test is raw-vs-refined
  recovery (`17d/e/h`) → **McNemar** on "complex recovered by refinement" (paired).

---

## 4. Rank / selection success — `rank_success_analysis.py`  (Benchmark only, 0 tests)

- `rank_success_curves`, `rank_validity_curves` (paired by complex): pointwise **Wilson CI**;
  compare tools at fixed k via **McNemar**; whole-curve via paired **permutation on success@k AUC**.
- `rank_success_by_flexibility` (disjoint rigid/mid/flexible bins, **unpaired**, ordered):
  **Cochran–Armitage trend test** for a monotone success-vs-flexibility trend (hand-roll as a
  linear-by-linear permutation), or Jonckheere–Terpstra on the continuous rot-bond axis.

---

## 5. Pose clustering — crystal pocket — `pose_cluster_crystal_pocket_report.py`  (Benchmark; 2 figs tested, ~30 not)

Homogeneity (`--stats`) and quality panel F are done. Untested figures worth a test:

- `summary_hit_rate`, `summary_triangulation`, `topN_crystal_reach_curves`, `rank_in_crystal_cluster_by_rank`,
  `topN_crystal_cluster_matrix`, `rank1_cluster_agreement_matrix`: "reach the crystal cluster"
  is a **paired binary** indicator across tools → **Cochran's Q + McNemar/Holm** at each depth
  (BH across depths) + Wilson CI. Agreement matrices → **exact McNemar** per cell, or an
  independence-null **permutation** ("do tools agree above chance?").
- Oracle-distance ECDFs (`summary_oracle_distance_ecdf`, `oracle_ecdf_split`, `top5_*_by_rank`,
  `cross_tool_agreement`): paired continuous distances → **Friedman / Wilcoxon** per-complex.
- `top5_pose_concentration`, `top5_distinct_sites`: per-tool paired distributions → Friedman/Wilcoxon.
- Descriptor panels (`descriptor_accuracy_vs_ro5` = unpaired 2-group → **Mann–Whitney + Cliff's δ**;
  `descriptor_success_vs_flexibility` = ordered bins → **Jonckheere–Terpstra**;
  `descriptor_property_correlations` already computes Spearman ρ but **discards p** → surface p + BH).
- `placement_vs_centroid` ARI: **permutation** of observed ARI vs a label-shuffle null.

---

## 6. Pose clustering — Orai × JKU — `orai_pose_cluster_report.py`  (JKU; 0 tests, small n)

Unit = 12 (frame,ligand) pairs. Reference axis = consensus site (no crystal).

- `orai_pb_validity_per_frame`: validity trend across the 4 ordered MD frames →
  **Cochran–Armitage / Jonckheere–Terpstra** (report as exploratory — n small).
- `orai_tool_agreement_matrix`, `orai_cross_tool_agreement`, `orai_tool_divergence_matrix`:
  cross-tool site distances paired by (frame,ligand) → **Wilcoxon** per tool-pair, or an
  independence-null **permutation** for above-chance agreement.
- `orai_valid_vs_consensus` (PB survivors vs failures, **unpaired** poses): distance-to-consensus
  → **Mann–Whitney + Cliff's δ**, but resample by (frame,ligand) — a **cluster permutation** — to
  avoid pose pseudoreplication.
- `orai_tool_dispersion`, `orai_axis_depth`: per-tool spread → Kruskal–Wallis/Dunn (or Friedman
  if paired by pair).
- `orai_descriptor_quality` panel C: Spearman ρ with p discarded → surface p + BH.

---

## 7. Consensus binding-region — `orai_ligand_region_consensus.py`  (Benchmark + JKU projection; 0 tests)

- `01_region_overview` panels C/D (region × method, region × frame): contingency tables →
  **χ² / G-test of independence** + **standardized residuals** (which cells drive it) + **Cramér's V**.
  Guard pseudoreplication: aggregate to per-ligand occupancy, or permute tool/frame labels *within
  ligand*.
- `03_ligandcluster_region_consensus`: chemo-cluster × region association → χ² / G-test.
- `04_ligand_consistency` purity histogram: test observed dominant-region purity vs a **random
  region-assignment permutation null** to show consistency exceeds chance.
- `06c_jku_vs_similar`: JKU vs chemically-similar-benchmark region occupancy → **permutation /
  χ²** on the two occupancy distributions (with n = 3 JKU ligands, strictly exploratory).
- `02_ligand_chemo_clusters` (KMeans + silhouette for k): add **gap statistic** or **bootstrap
  cluster-stability** to justify the chosen k.

---

## 8. Transmembrane-pore exclusion — `orai_transmembrane_exclusion.py`  (both Orai; 0 tests)

- `fig_tool_tm_loss` (toolchain × fate counts): **χ² / G-test of independence** ("does the
  fraction lost to the TM pore differ by tool?") + standardized residuals + Cramér's V; per-tool
  TM fraction with Wilson CI. Aggregate to (frame,ligand) or cluster-permute.
- `fig_pbvalid_lost_matrix` (frame × ligand): test frame effect on loss via
  Cochran–Armitage; otherwise descriptive.
- `fig_rank_survival`: survival vs rank per tool → descriptive; compare tools via permutation on AUC.

---

## 9. PandaMap interactions — `pandamap_interaction_report.py`  (both; fig 14 Spearman, 06b Pearson w/o p)

- `01_interaction_profile`, `02_type_heatmap`, `02b_total_boxplot`: interaction counts per pose,
  paired by complex → **Friedman + Wilcoxon** per interaction type (BH across the 16 types);
  totals → Friedman/Wilcoxon.
- `03_residue_hotspots`, `07_residue_class`, `08_ligand_element`: categorical contingency
  (residue/class/element × tool or × interaction type) → **χ² / G-test** + standardized residuals;
  per-residue contact fraction across tools → **McNemar** (BH across residues).
- Native recovery `04/04b/04c/09/10/12/13` (Benchmark): per-pose precision/recall/F1 vs crystal,
  paired across tools → **Friedman / Wilcoxon** on per-complex F1 + bootstrap CI on the mean bars.
  `12_typed_vs_loose` and `13_native_f1_oracle` are within-pose / within-complex paired → **Wilcoxon**.
- `05_fingerprint_similarity` (Jaccard): compare tool↔crystal Jaccard across tools →
  Friedman/Wilcoxon; test "above chance" vs a contact-shuffle **permutation** null.
- `06b_descriptor_corr` (Pearson, no p): surface p + BH; prefer **Spearman** (interaction counts
  are non-normal).

---

## 10. Docking difficulty — `ligand_docking_difficulty.py` / `receptor_docking_difficulty.py`  (Benchmark; Cliff's δ only)

The δ heatmaps (best-N vs worst-N, unpaired) show effect size but no p / CI:
- Add **Mann–Whitney U** p per descriptor with **BH** across the 14 (ligand) / 13 (receptor)
  rows × 3 tools, and a **bootstrap CI on Cliff's δ**.
- Stronger than the arbitrary N=60 extreme-group split: **Spearman ρ (+ BH)** of each descriptor
  vs continuous oracle RMSD across *all* complexes, and/or a **logistic regression** of success ~
  descriptor (per-tool). Report both so the difficulty drivers rest on the full sample, not the tails.

---

## 11. Docking effort — `docking_effort_comparison.py`  (Benchmark; 0 tests)

- `wall_per_complex_distribution` (boxplot, same complexes timed → paired): **Friedman + Wilcoxon**
  (Holm) on per-complex wall-clock + median-ratio with bootstrap CI.
- `wall_validity_rate`, `wall_poses_generated_vs_valid`: paired proportions → **McNemar / Cochran's Q**
  + Wilson CI.
- Aggregate total bars: descriptive; add bootstrap CI on per-pose cost. (Effort is a
  resource claim, not a hypothesis — the main win is pairing the per-complex time comparison.)

---

## 12. gnina re-ranking — `diffdock_gnina_rerank_analysis.py` / `_chemotype.py`  (Benchmark; partly tested)

- `selection_benefit_*.png` panels A & C: the "+X pp from full pipeline" claim currently rests on
  descriptive proportions. Add **McNemar** (Holm) across the selection strategies (native / rerank /
  minimize / full), paired by complex — this is the test that substantiates the headline benefit.
- `selection_benefit_*.png` panel B (paired RMSD scatter): **Wilcoxon signed-rank** native vs full
  + Hodges–Lehmann CI.
- Chemotype Mann–Whitney U is computed but CSV-only → **surface on panel A** and BH across predictors.

---

## 13. Pocket & receptor comparison — `pocket_comparison_report.py` / `receptor_comparison_report.py`  (Benchmark; partly tested)

- Pocket success curves (18, 19) already carry Wilson CI. Add **McNemar** (paired by receptor) to
  compare fpocket vs p2rank at matched top-N / threshold.
- Pocket top-1 distance ECDF / paired scatter (10, 12): paired continuous → **Wilcoxon signed-rank**
  (fpocket vs p2rank DCC per receptor) + Hodges–Lehmann CI.
- Attribute-ranking Spearman (8, 9): surface p + BH.
- `receptor_comparison_report.py` (Orai n=4 vs benchmark, different proteins, unpaired): the
  descriptive "Orai sits at the Xth percentile" is honestly the right call at n=4 — a formal test
  would be underpowered and is not recommended. Keep descriptive; state n=4.

---

## 14. Priority order (highest validity payoff first)

1. **Validity report (both datasets)** — Cochran's Q / McNemar per tool & per check (§2). Headline
   result, zero tests, directly substantiates the geometry-failure and gnina>smina claims.
2. **Pose-comparison success bars & oracle RMSD (Benchmark)** — Friedman/Wilcoxon + McNemar (§3).
   The core accuracy claims.
3. **gnina-rerank selection benefit (Benchmark)** — McNemar on the strategies (§12). Makes the
   advertised "+pp" honest.
4. **Fix the "computed-but-hidden" p-values** — surface the discarded Spearman/Mann–Whitney p's
   already sitting in CSVs (difficulty, region/cluster descriptor panels, pocket attr, PandaMap 06b),
   add BH. Near-zero effort.
5. **Transmembrane & region contingency tables (Orai)** — χ²/G-test + Cramér's V + residuals (§7,§8).
6. **Rank/threshold curves everywhere** — Wilson CI bands + matched-depth McNemar (§4, §5, §9).
7. **Orai × JKU cluster/validity** — exact/permutation on 12 units, flagged exploratory (§6).

**Reusable scaffolding:** the Holm, McNemar, Cochran's Q, Dunn, Jonckheere–Terpstra, permutation,
Wilson-CI and Cliff's-δ helpers already exist in `pose_cluster_crystal_pocket_report.py` and
`pocket_comparison_report.py`. Extract to a shared `stats_utils.py` before rolling these out.
