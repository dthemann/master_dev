# Implemented statistical tests — what was added to each graph, and what it means

Every test below is computed by the shared `Scripts/Analysis/stats_utils.py` (scipy-only,
23/23 self-tests pass), annotated on the figure, and written in full to a `*_stats.json`
sidecar next to it. Numbers quoted are from the verification runs (Benchmark n≈303 complexes;
Orai×JKU n≈12–16 (frame,ligand) pairs). Three rules were enforced everywhere: poses are
**aggregated to one value per complex / per (frame,ligand) pair** before testing (never pooled
as independent); **paired** tests are used across tools on the same complexes; families of tests
are **BH-FDR corrected**. Orai×JKU results are labelled **exploratory** (tiny n).

Legend: `***` p<.001 · `**` p<.01 · `*` p<.05 · `ns` not significant.

---

## A. PoseBusters validity report — `posebusters_validity_report.py` (Benchmark **and** Orai×JKU)

| Graph | Test added | What it implies |
|---|---|---|
| **01_per_tool_validity** / **04_valid_per_pair_distribution** | Friedman + Kendall's W on per-complex valid-pose **count**; Wilcoxon signed-rank post-hoc (Holm); companion Cochran's Q on "≥1 valid pose" | Tools differ enormously in how many valid poses they yield: **Friedman p≈0, Kendall's W=0.66, n=303** — a strong, consistent ranking. AutoDock ≫ raw EquiBind (rank-biserial=0.99); AutoDock vs raw-DiffDock rb=0.87 `***`. Establishes the headline "tools are not interchangeable on physical validity." |
| **03_valid_grouped_bars** | Paired-proportions (Cochran's Q + pairwise McNemar/Holm + Wilson CIs) on per-complex "≥1 valid pose" | Cochran's **Q=2131, p≈0, df=12**. AutoDock produces a valid pose on far more complexes than raw DiffDock (56 vs 6 discordant complexes, p_holm=9×10⁻¹⁰) or raw EquiBind (262 vs 2, p=2×10⁻⁷³). |
| **05_per_check_passrate** | Per PB check: Cochran's Q on per-complex "any pose passes"; **BH-FDR across the 20 checks**; significant checks starred | Pinpoints *where* tools diverge: **bond_lengths, bond_angles, internal_energy, minimum_distance_to_protein, volume_overlap_with_protein** (+cofactor/water clashes) are the BH-significant checks. This is the statistical backing for "EquiBind's failures are covalent-geometry (bond angles/lengths), plus protein clashes." |
| **06_equibind_variant_validity** | Paired-proportions + Friedman on valid-count across the 9 EquiBind variants | Variants differ (Friedman **p≈4×10⁻²⁹⁶, W=0.58**): smina/gnina refinement lifts validity massively over raw (raw→smina rb=−1.0), and gnina beats smina within every pocket mode. |
| **07_diffdock_variant_validity** | Paired-proportions + Friedman across raw / smina / gnina | Optimisation is decisive: **Friedman p=3×10⁻⁹⁷, W=0.73**; raw→gnina rb=−0.99; and critically **gnina > smina, rb=0.55, p_holm=4×10⁻¹⁰** `***` — the minimiser choice is real, not noise. |
| **07b_variant_validity_matrix** | Cochran's Q per tool's optimiser axis | Confirms the raw→smina→gnina lift is significant for both DiffDock and EquiBind on the "≥1 valid pose" level. |
| **Orai×JKU** (same figures) | Same tests, guarded, **exploratory** | Even at n≈12 the geometry signal survives: valid-count Friedman **p=9×10⁻⁴, W=0.92**, and the per-check BH still flags **bond_lengths, bond_angles** — the same failure mode as the Benchmark, on independent ligands. |

---

## B. Pose accuracy vs crystal — `posebusters_pose_comparison.py` (Benchmark)

| Graph | Test added | What it implies |
|---|---|---|
| **01_oracle_rmsd_cdf** / **02_oracle_rmsd_boxplot** | Friedman + Kendall's W + Wilcoxon/Holm on per-complex **oracle RMSD** | Accuracy ordering is unambiguous: median oracle RMSD **AutoDock 1.01 Å < DiffDock 1.56 Å < EquiBind 6.74 Å**; AutoDock vs DiffDock rb=−0.48 `***` (p=7×10⁻¹³), both ≫ EquiBind (rb=−0.95 `***`). EquiBind rarely reaches a near-native pose even at oracle depth. |
| **03_oracle_vs_top1_success** / **10_pb_valid_rmsd2_success_bars** | Cochran's Q + McNemar (cross-tool) and McNemar within-tool (oracle vs top-1); Wilson CIs on bars | Tests both "who is more accurate" (cross-tool) and "how much ranking costs you" (oracle vs top-1). CIs make the ~2 Å success bars comparable; the within-tool McNemar quantifies the top-1-vs-oracle gap that ranking imperfection imposes. |
| **11_vs_posebusters_paper** | One-sample exact **binomial** vs the paper's published rate; Wilson CI on our rate | Places our AutoDock rate (171/303 ≤2 Å) against the published PoseBusters number with an honest CI — flagged **descriptive** because the complex sets differ (not a like-for-like test). |
| **19_within2_validity_dumbbell** | Per method: within-method **Wilcoxon signed-rank** of the oracle (min-RMSD) pose vs the near-native-pool PB-validity, per complex; **cluster (complex) bootstrap** 95% CI on the circle, **Wilson** CI on the diamond | Tests whether the closest-to-native pose is systematically more/less physically valid than a typical near-native pose. AutoDock and optimised DiffDock show **no gap** (ns); raw DiffDock `*`, and EquiBind-unguided smina (rb=0.63, **p=0.009**) / gnina (`*`) show the oracle pose is *more* valid than the pool. CI whiskers + a right-hand significance column are drawn on the figure. |
| **17f / 17g / 17h_pb_top30_recovery_by_test** | Per panel, on the per-complex "passes all PoseBusters tests" outcome: exact **McNemar** for rank-1 vs best-of-top-30 (the coloured recovery bars) and, on refined panels, raw vs refined at **both** rank-1 and top-30 depths; **Wilson 95% CI** on every pass-all rate; all p-values in a figure form one **Holm-corrected** family. Annotated in each panel's note box and written to the same `posebusters_pose_comparison_stats.json` sidecar (`top30_recovery_{refined,raw,raw_vs_refined}`). | Turns the "+k recovered" bars and the "raw → refined +pp" callouts into tested claims: whether offering the ranker top-30 (or refining the geometry) recovers *significantly* more complexes that pass everything, or sits within sampling noise. McNemar is two-sided — the top-N pick minimises RMSD, not validity, so a complex can pass at rank-1 yet fail the top-N pick — and the unit is one representative pose per complex (no pose-level pseudoreplication). |
| **15b_topk_recovery_validity** | On the per-complex "≥1 top-k pose is near-native (RMSD ≤ 2 Å)" outcome: exact **McNemar** between pre-specified tool pairs at **pre-specified depths k ∈ {1, 5, 10, 15}** — one **Holm-corrected** family of 24 tests. The 30 cumulative curve points are nested (recovered by k ⟹ by k+1) and auto-correlated, so a fixed small depth-set is tested rather than every point. **Wilson 95% CI** bands are shaded on every curve (a proportion CI — stays valid at the EquiBind ≈1% / AutoDock ≈89% extremes where a Gaussian/Wald interval would spill past 0/100%). The near-native-vs-validity-aware gap is reported as an **effect size** (accurate-but-invalid share + Wilson CI), *not* a p-value, because the validity set is a strict subset of near-native so McNemar there is degenerate. The full numeric results (each rate + Wilson CI, McNemar p / Holm p / star) are written to `topk_recovery_stats.csv` + the `topk_recovery` block of the sidecar; the figure itself shows the CI as shaded bands only. | Turns the recovery curves into tested claims at a non-pseudoreplicated depth-set. Headline: **gnina refinement of DiffDock is not significant at top-1 (p_holm=0.11) but is from k=5 on (p_holm≤0.02)** — the reranking benefit needs depth. EquiBind raw→gnina and every cross-tool gap are `***`. The nested gap's honest statistic is its size: raw DiffDock loses **15.8 pp [12.2, 20.4]** of its near-native recoveries to invalidity at k=1, vs only ~3 pp once refined. |

| **15c / 15d_optimization_benefit_by_rank** (standalone `optimization_benefit_stats.py`; reads `per_pose_metrics.csv`, the same rows fig 15c is built from) | Two pre-specified questions on the per-pose 0/1 outcomes (the plotted CSV holds only marginal %+n, which cannot recover the paired discordant cells). **(a) "does optimization help?"** — per (tool, endpoint ∈ {near-native, PB-valid}) at **rank 1 (primary, clustering-proof) + pre-specified ranks {5,15,30}**: **Cochran's Q** across {raw, smina, gnina} + pairwise **exact McNemar** (Holm within the triple) + **Wilson CIs** and the paired risk-difference (pp); omnibus family **BH-FDR** corrected. `both_%` is deterministic (near AND pbv) → reported descriptively, never tested. **(b) "does the benefit depend on rank?"** — per (tool, endpoint, contrast smina/gnina-vs-raw): logistic `logit P(pass) = b0 + b1·opt + b2·rank + b3·(opt×rank)`; **b3** is the single rank-dependence coefficient, inferred by **cluster (complex) bootstrap** (no statsmodels/GEE in the env). Raw-only rank slope **b2** is reported separately — DiffDock = confidence gradient, EquiBind = pre-specified **negative control** (generation order). Results → `optimization_benefit_stats_pairwise.csv` / `_interaction.csv` / `.json`, and **figures 15c/15d are annotated** by the script (it reuses the pipeline's own plot functions via an optional `stats=` payload, so the plotting stays single-sourced): 15c gets per-rank McNemar stars at the tested depths + a rank-1 summary box; 15d gets the interaction result, plain-language trend, and the raw-slope negative-control line. A plain report run (`posebusters_pose_comparison.py`) still emits un-annotated 15c/15d; running `optimization_benefit_stats.py` afterwards overwrites them with the annotated versions (use `--no-annotate` to skip). | Separates physics from correctness. **(a) PB-valid: massive and universal** — DiffDock +45→67 pp, EquiBind +40→52 pp, **McNemar p_holm<1e-4 `***`** at every rank, near-one-directional flips (e.g. rank-1 DiffDock 140 gain / 3 lose; many depths 0 losses); gnina ≥ smina. **Near-native: essentially null** — DiffDock +0–2 pp, all `ns` after Holm (rank-1 Cochran p=0.045 fails BH); EquiBind a marginal +2–3 pp on a ~1% base (gnina reaches Holm `*` only at rank 15). Minimisation fixes geometry, it does not move a wrong pose toward native. **(b) The PB-valid benefit does depend on rank — but grows with depth, not shrinks**: DiffDock b3=+0.021/+0.029 per rank (**BH p=0.008 `**`**), because raw PB-validity **decays** with confidence rank (raw slope −0.039 log-odds/rank, all-boots one-sided) while the optimised rate stays high/flat → the gap widens from ~+51 pp (ranks 1–10) to ~+59–63 pp (21–30). It rescues the low-confidence tail hardest. Near-native interaction is null (nothing to depend on rank). **Negative control validated**: DiffDock raw slopes are significantly negative (confidence rank genuinely orders quality), EquiBind raw slopes are `ns` (generation order carries none), and EquiBind's PB-valid interaction is `ns` after BH. |

| **20_form_fidelity__rank_depth_gate_vs_rank** (`_stats_form_fidelity_gate_vs_rank`; annotated on all three panels + the `form_fidelity_gate_vs_rank` block of the sidecar) | Four per-complex paired tests, each tied to a panel claim (poses aggregated to one value per complex first — no pose-level pseudo-replication). **Gate effect (Panel C):** the near-native cohort is a strict SUBSET of PB-valid, so rather than the degenerate subset-vs-superset contrast, each tool's PB-valid poses are PARTITIONED into near-native (≤2 Å) vs far (>2 Å); per complex holding both, their median best-fit (Kabsch) form RMSD is compared with a paired **Wilcoxon signed-rank + Hodges–Lehmann median-difference CI** (far − near). **Cross-tool form (Panels A/B):** **Friedman + Kendall's W + Wilcoxon/Holm** on the SHALLOWEST-ranked PB-valid (A) / near-native (B) pose per complex, listwise-complete across tools. **Rank trend:** per-complex **Kendall τ**(form, effective rank), then a one-sample **Wilcoxon** of the τ vs 0 (repeated ranks within a complex are correlated → one slope per complex, never pooled; EquiBind gnina-rank = pre-specified negative control). **Form⊥placement coupling:** per-complex **Spearman ρ**(form, in-place), one-sample Wilcoxon vs 0. Each per-tool family is Holm-corrected. | Turns the figure's headline "the ≤2 Å gate improves form fidelity" into a tested claim. **The gate genuinely removes worse-form poses:** far−near median-form Δ = **AutoDock +0.64 Å [0.55, 0.74], DiffDock* +0.44 [0.29, 0.61], EquiBind* +0.36 [0.19, 0.54], all `***`**. It can do so because form and placement are only **weakly-to-moderately coupled** (per-complex Spearman **ρ = 0.26 / 0.20 / 0.10, all `***`**) — best-fit RMSD removes rigid-body placement by construction, yet a residual coupling survives, so the in-place gate still pulls median form down. Cross-tool: **AutoDock ≈ DiffDock* (ns), both beat EquiBind* `***`** (Panel A **W=0.06** — the ordering is almost entirely the EquiBind gap). Form error **degrades with rank for every tool** (**τ = 0.16 / 0.10 / 0.05, all `***`**), and the EquiBind negative control is **not** null → gnina-affinity rank carries some form-ordering signal. Panel B is EquiBind-limited (**n=41** complete blocks). |

*(Scope note: the remaining pose-comparison figures — the other waterfalls, 09-series, 16-series and the rest of the 20-series — were left for a later pass; the `20d…__stats` figure already had its own tests, the 17f/17g/17h top-30 recovery figures plus 15b carry the McNemar/Wilson tests above, fig 18 carries the within-threshold paired tests, and **20_…rank_depth_gate_vs_rank** now carries the per-complex paired tests above. The 15c/15d optimization-benefit tests live in the standalone `optimization_benefit_stats.py`, which writes the sidecars **and** overwrites 15c/15d with the annotated panels.)*

---

## C. Rank / selection success — `rank_success_analysis.py` (Benchmark)

| Graph | Test added | What it implies |
|---|---|---|
| **rank_success_curves** | Pointwise **Wilson 95% CI** bands on success@k; **McNemar** across tools at top-1 and top-5 (Holm) | Success@1 now carries CIs, e.g. **AutoDock 56.4% [50.8, 61.2], n=303**, rising to ~89% by top-30. Matched-depth McNemar says whether the tool gaps at top-1/top-5 are real rather than eyeballed from overlapping curves. |
| **rank_validity_curves** | Wilson CI bands + McNemar at top-1/top-5 on "a valid pose within top-k" | Separates "can rank a pose" from "can rank a *valid* pose"; the CIs expose where the validity curve is genuinely below the accuracy curve. |
| **rank_success_by_flexibility** | **Cochran–Armitage** trend across the ordered rigid/mid/flexible rot-bond tertiles | Tests the claim "flexible ligands dock worse" as a monotone trend rather than three disconnected bars. |

---

## D. Computational effort — `docking_effort_comparison.py` (Benchmark)

| Graph | Test added | What it implies |
|---|---|---|
| **wall_per_complex_distribution** | Friedman + Kendall's W + Wilcoxon/Holm on per-complex wall-clock; median-ratio bootstrap CI | Runtime ordering is essentially deterministic: **Friedman p=5×10⁻⁹⁸, W=0.94, n=239 paired complexes**; AutoDock vs DiffDock rb=−1.0 (AutoDock faster on *every* common complex). The cost differences are not anecdotal. |
| **wall_validity_rate** / **wall_poses_generated_vs_valid** | Cochran's Q + McNemar on per-complex "produced ≥1 valid pose"; Wilson CIs | Validity-yield differs by method (**Cochran Q=40.3, p=2×10⁻⁹**), so cost-per-valid-pose comparisons rest on a real yield gap, not sampling noise. |
| aggregate total / per-pose-cost bars | left descriptive (noted in sidecar) | These are resource totals, not paired samples — no test is appropriate. |

---

## E. gnina re-ranking — `diffdock_gnina_rerank_analysis.py` + `_chemotype.py` (Benchmark)

| Graph | Test added | What it implies |
|---|---|---|
| **selection_benefit panel A** (strategy bars) | Cochran's Q + McNemar/Holm across native / rerank / minimize / full; Wilson CIs | **Honest negative-ish result:** omnibus **Q=7.11, p=0.068** (n.s.), and full-vs-native McNemar 36/24, p_holm=0.62 `ns`. The advertised "+pp from the full pipeline" is **not** a significant top-1 selection gain — consistent with it being a modest recovery effect, not a headline win. Full *does* beat rerank-only (20/4, p=0.009 `**`). |
| **selection_benefit panel B** (native vs full RMSD) | Wilcoxon signed-rank + Hodges–Lehmann median-diff CI (paired) | Quantifies the per-complex RMSD shift from the full pipeline with a CI instead of a raw win/lose count. |
| **selection_benefit panel C** (decomposition) | McNemar of rerank-only / minimize-only / full each vs native | Attributes any change to the minimise vs rerank component; both are individually n.s. vs native. |
| **chemotype panel A** (PCA) | Mann–Whitney U + Cliff's δ on high- vs low-communality benefit | High- and low-communality ligands don't differ in benefit (**U, p=0.16, δ=−0.09 `ns`**): the benefit is **not** chemotype-clustered. |
| **chemotype panels B/C** | Surface the previously-discarded Spearman p; **BH-FDR** across 22 predictors | The one real driver: **benefit vs DiffDock native RMSD, ρ=0.48, p=8×10⁻¹⁹** — re-ranking helps precisely when the native top-1 was already poor. It's a **recovery** effect, statistically, not a chemistry-predictable one. |

---

## F. Docking difficulty — `ligand_docking_difficulty.py` + `receptor_docking_difficulty.py` (Benchmark)

| Graph | Test added | What it implies |
|---|---|---|
| **ligand_difficulty_delta_heatmap** | Mann–Whitney U per descriptor×tool (worst-60 vs best-60) with **BH-FDR** + Cliff's-δ bootstrap CI; **plus a full-sample Spearman** of each descriptor vs oracle RMSD over all 303 complexes | Moves difficulty from an arbitrary extreme-group δ to a tested, whole-sample correlation. Strongest AutoDock drivers: **TPSA ρ=0.41 (p=5×10⁻¹⁴), HBA ρ=0.39, rotatable bonds ρ=0.39, HBD ρ=0.37, MW ρ=0.36** — bigger, more polar, more flexible ligands are harder, and now with p-values + FDR control. |
| **receptor_difficulty_delta_heatmap** | Same design over 13 receptor/pocket properties × tool | Identifies which receptor/pocket properties (size, druggability, pocket count) actually track difficulty, FDR-controlled, rather than reading colour off an untested heatmap. |

---

## G. Pocket predictor comparison — `pocket_comparison_report.py` (Benchmark; fpocket vs p2rank)

| Graph | Test added | What it implies |
|---|---|---|
| **crystal_top1_distance_ecdf** / **crystal_paired_top1_scatter** | Wilcoxon signed-rank + Hodges–Lehmann CI (paired by receptor) | The two predictors' top-1→crystal distances differ significantly (**rb=0.37, p=3×10⁻¹¹**) — p2rank places its top pocket closer to the ligand. |
| **crystal_success_vs_topn** / **crystal_success_vs_threshold** | McNemar (paired by receptor) at matched N/threshold (keeps existing Wilson bands) | At rank-1 p2rank finds the true site far more often: **59.6% vs 30.4%, McNemar p=2×10⁻²³** `***` (n=421). Quantifies the p2rank rank-1 advantage the Wilson bands only hinted at. |
| **pocket_attr_ranking_fpocket / _p2rank** | Surface the Spearman p + **BH-FDR** across attributes | fpocket's own score drives its ranking strongly (**ρ=−0.54, p=1×10⁻⁹⁵**), alpha-spheres and volume next — the ranking is explained by these attributes, now significance-tested. |

---

## H. Consensus binding-region — `orai_ligand_region_consensus.py` (Orai×Benchmark + JKU projection)

| Graph | Test added | What it implies |
|---|---|---|
| **01_region_overview** (region×method) | G-test of independence + Cramér's V + residuals, on a distinct-ligand table | Tools occupy **different regions**: **G=1413, p≈10⁻²⁹¹, V=0.46** — a large, real association (the docking method strongly conditions where the ligand lands). |
| **01_region_overview** (region×frame) | G-test | MD frame matters far less: **G=345, p=2×10⁻⁵⁸ but V=0.11** — statistically detectable, practically small. |
| **03_ligandcluster_region_consensus** | G-test (chemo-cluster × region) | **Not** significant: **G=21, p=0.17, V=0.03** — ligand chemotype does **not** determine binding region here. |
| **04_ligand_consistency** | Label-shuffle **permutation** null (2000 iters) for dominant-region purity | Ligands are more region-consistent than chance: **observed purity 0.29 vs null 0.21, p=5×10⁻⁴** `***` — the consensus regions are real structure, not an artefact of clustering. |
| **06c_jku_vs_similar** | G-test on JKU vs similar-benchmark region occupancy | JKU blockers occupy **different regions** than their nearest benchmark analogues (**G=153, p=5×10⁻²⁹, V=0.45**) — exploratory (3 JKU ligands) but a clear separation. |

---

## I. Transmembrane-pore exclusion — `orai_transmembrane_exclusion.py` (both Orai datasets)

| Graph | Test added | What it implies |
|---|---|---|
| **fig_tool_tm_loss** | G-test on toolchain × fate + Cramér's V + residuals; per-tool TM-loss fraction with Wilson CI; **per-(frame,ligand) aggregate reported alongside** | Which tool dumps poses in the transmembrane pore differs strongly (**G=431, p=6×10⁻⁹⁰, V=0.49**). The sidecar carries an explicit pseudoreplication caveat (pose counts overstate significance) and the honest **n=14 (frame,ligand)** unit. |
| **fig_pbvalid_lost_matrix** | Cochran–Armitage trend in TM-loss across the ordered MD frames | Tests whether later MD frames lose more/fewer poses to the pore as a monotone trend. |
| **fig_rank_survival** | left descriptive | A tool-vs-tool AUC test would need far more independent units than 3 ligands provide. |

---

## J. Pose clustering, Orai×JKU — `orai_pose_cluster_report.py` (no crystal, **n≈12, all exploratory**)

| Graph | Test added | What it implies |
|---|---|---|
| **orai_pb_validity_per_frame** | Cochran–Armitage across the 4 ordered MD frames, per tool | Trend test for validity vs MD frame (exploratory at this n). |
| **orai_tool_agreement / cross_tool / divergence** | Wilcoxon signed-rank of per-(frame,ligand) inter-tool site distance vs the 5 Å agreement threshold (Holm across pairs) | Tools essentially **never** converge on the same site: AutoDock vs DiffDock median inter-tool distance **30 Å, agreement fraction 0, p=1×10⁻³** — the cross-tool disagreement on Orai is real, not small-sample scatter. |
| **orai_valid_vs_consensus** | Mann–Whitney + Cliff's δ, with **cluster-bootstrap CI** by (frame,ligand) | PB-valid and PB-invalid poses are **not** at different distances from consensus (**U, p=0.17, δ=−0.15 `ns`**): validity doesn't explain consensus proximity here. |
| **orai_descriptor_quality panel C** | Surface Spearman p + BH-FDR across properties | Recovers which ligand properties track Orai docking quality with actual p-values. |

---

## K. PandaMap interactions — `pandamap_interaction_report.py` (both datasets)

| Graph | Test added | What it implies |
|---|---|---|
| **01_interaction_profile** / **02b_total_boxplot** | Friedman + Wilcoxon on per-complex mean count per tool, per interaction type; **BH-FDR across 16 types** | Tests, type by type, whether tools form different interaction spectra (FDR-controlled). On Orai×JKU only 2 complete units exist → **correctly skipped and labelled exploratory** (n=2 too small). |
| **03_residue_hotspots** | G-test (residue × tool) + per-residue McNemar with BH | Which contact residues are tool-dependent hotspots vs shared. |
| **04_native_recovery** | Friedman on per-complex native-F1 + cluster-bootstrap CIs on the precision/recall/F1 bars (Benchmark) | Puts CIs on the recovery bars and tests whether tools recover native contacts differently. |
| **06b_descriptor_corr** | Switched to Spearman, **surfaced the p**, BH across the grid | Marks which descriptor↔interaction correlations are real (was an untested Pearson heatmap). |
| **07_residue_class** / **08_ligand_element** | G-test of independence + Cramér's V | Whether interaction type is associated with residue class / ligand element. |
| **12_typed_vs_loose** / **13_native_f1_oracle** | Wilcoxon signed-rank + BH across tools | Tests the typed-vs-loose recall gap and the oracle-vs-top1 F1 gap within each tool. |
| **14_rmsd_vs_recovery** | Kept per-method Spearman, added BH across methods | The RMSD↔recovery correlation, now multiplicity-controlled. |
| **04b_native_recovery_cdf** | Per-tool **total-miss** rate (fraction of poses with F1≈0 — the CDF's y-intercept) with a **Wilson 95% CI**, plus a **per-complex miss-fraction Friedman** across tools (the CDF is pose-level/pseudoreplicated, so the paired Friedman is the inference-grade comparison). | Turns the CDF's zero-recovery jump into a tested claim: AutoDock essentially never totally misses (**0% [0.0, 0.3]**), while DiffDock* and EquiBind* leave **~22%** of poses recovering *zero* native contacts (**22.3% [20.2, 24.5]** / **22.0% [19.8, 24.5]**); the rate differs across tools (**Friedman W=0.20, n=236, p<1e-4 `***`**). |
| **04c_native_recovery_by_rank** / **10b_contact_decomposition_by_rank** | Per-complex **Kendall τ**(pose rank, score) → **one-sample Wilcoxon vs 0**, per tool, **BH-FDR** across the tool×metric family (one slope per complex — correlated poses never pooled; EquiBind gnina-affinity rank = pre-specified negative-control-ish baseline). A negative τ on a quality metric ⇒ the tool's own ranking puts better poses first. | Tests whether each tool's ranking actually orders pose quality. AutoDock does strongly (F1 **τ=−0.40 `***`**; matched **−0.36 `***`**, spurious **+0.32 `***`**); EquiBind moderately (F1 **τ=−0.20 `***`**); **DiffDock*'s confidence rank does NOT** (F1 **τ=0.00 `ns`**, matched only −0.11 `*`) — the statistical backing for DiffDock's large oracle-vs-top-1 gap in fig 13. |
| **09_type_resolved_recovery** (+ `_top1`) | Per-cell **Wilson 95% CI** on every recall/precision (recall k=tp n=tp+fn; precision n=tp+fp), **printed in each heatmap cell** (value with its CI) and written to the CSV + sidecar; low-n (native n<200) and saturated (recall=precision=1.00 for all tools) types recorded in the sidecar as `flagged_low_n` / `flagged_saturated`. | Makes the per-type heatmap read honestly and self-contained: the AutoDock > DiffDock* > EquiBind* ordering holds per type, but the in-cell CIs show which cells are trustworthy — **halogen_bonds** (n=64) is visibly wide (e.g. recall 0.54 [0.42, 0.66]) and the **metal_coordination** row is saturated (1.00 [0.99, 1.00] for all — trivially recovered), so neither should be over-read. |

---

## L. Pose clustering, crystal pocket — `pose_cluster_crystal_pocket_report.py` (Benchmark)

*(This script already had the richest tests — `crystal_cluster_homogeneity` under `--stats` and
`cluster_quality_metrics` panel F — which were left untouched. Four additions:)*

| Graph | Test added | What it implies |
|---|---|---|
| **summary_hit_rate** | Cochran's Q + pairwise McNemar/Holm across the 6 oracle sources (per-complex "oracle ≤ thr") | Tests whether docking tools + pocket predictors differ in how often they hit the crystal site. |
| **topN_crystal_reach_curves** | Paired-proportions across tools at each depth N∈{1,5,10,15} | Tests whether tools reach the crystal cluster at equal rates at matched depth. |
| **descriptor_accuracy_vs_ro5** | Mann–Whitney + Cliff's δ on oracle RMSD, Ro5 pass vs fail | Whether Rule-of-5 compliance separates docking accuracy. |
| **descriptor_property_correlations** | Recompute Spearman **keeping the p** (was discarded), BH across descriptors | Marks which descriptor↔accuracy correlations are significant. |

---

## M. Orai × JKU crystal-free completeness pass (2026-07-19)

A dedicated sweep to close every **produced** Orai×JKU figure that was still descriptive
but where a crystal-free test is feasible. All new tests aggregate to the ~9–16
`(frame,ligand)` unit (no pose pseudoreplication), are guarded, and are labelled
**exploratory**; full results land in each script's existing `*_stats.json` sidecar.

| Figure (script) | Test added | Result on Orai×JKU |
|---|---|---|
| **divergence-matrix RIGHT panel** (`orai_pose_cluster_report`) | Wilcoxon signed-rank of the per-pair CLOSEST cross-tool pose distance vs the agree threshold, Holm across tool-pairs (`_stats_cross_pose_touch`); RIGHT column headers now starred like the LEFT | pose *clouds* graze but don't converge: all pairs `ns` after Holm (p_holm≈0.09) |
| **cluster_structure / Panel C** (`orai_pose_cluster_report`) | Jonckheere–Terpstra trend across ordered MD frames on per-pair #clusters + dominant fraction (`_stats_cluster_structure_trend`); trend p in the legend | #clusters trend z=−1.78 **ns** (p=0.076); dominant-frac ns |
| **tool_cluster_contribution** (`orai_pose_cluster_report`) | Paired Friedman + Kendall's W + Wilcoxon/Holm across tools on per-pair distinct-cluster count (`_stats_tool_cluster_contribution`) | Friedman `*` p=0.047, W=0.51, n=6; pairwise ns after Holm |
| **tool_dispersion / Panel E** (`orai_pose_cluster_report`) | Paired Friedman/Wilcoxon across tools on per-pair pose spread (`_stats_tool_dispersion`) | Friedman ns p=0.31 |
| **fig_tool_tm_loss** (`orai_transmembrane_exclusion`) | Pseudoreplication-**free** companion to the pooled-pose G-test: paired Friedman across tools on the per-(frame,ligand) TM-loss fraction + bootstrap CI on each tool's mean per-unit loss | The honest contrast: pooled-pose G=562.7 `***` vs per-unit Friedman **ns** p=0.14 (n=7); means AutoDock 0.32 [0.18,0.47], DiffDock\* 0.15 [0.07,0.23], EquiBind\* 0.70 [0.40,0.99] |
| **06a_jku_where / 06b_similar_where** (`orai_ligand_region_consensus`) | Adaptive: distinct-ligand region×tool **G-test** when ≥5 distinct ligands, else an explicit *"descriptive — n too small, exploratory"* label (mirrors 06c) | 06a (n=4 JKU ligands) → descriptive-labelled; 06b (similar set) → G=25.5, p=0.062 `ns`, V on distinct-ligand counts |

**Deliberately left test-free** (crystal-free but no meaningful test / no host graph):
PandaMap has no crystal-free *by-rank* figure to attach a Kendall-τ(rank, interaction-count)
test to (the τ helper is wired only to the crystal-gated 04c/10b); adding a whole new figure
was out of scope for "tests in the **produced** graphs". Every produced Orai×JKU PandaMap
graph (01/02b Friedman, 03 McNemar+G, 06b Spearman, 07/08 G-test) already carries its
feasible test. **Known quality caveat (not a missing test):** the pooled-pose G-tests
(PandaMap 03-title/07/08; TM Fig 1/2) still stamp stars from pseudoreplicated poses — honest
only in the JSON `unit` field; the TM per-unit companion above is the inference-grade view.

---

### Still descriptive (no test — by design)
3-D scene figures, example-cluster panels, per-pair count heatmaps (validity 02/02b), failure
waterfalls (validity 09; pose-comparison 17-series), the ~44 deferred pose-comparison figures,
and aggregate resource totals — these are illustrative or not paired samples. See
`STATISTICAL_VALIDATION_PLAN.md` §3 and §14 for the deferred pose-comparison list.
