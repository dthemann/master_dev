# Plan: switch the primary near-native endpoint to the nearest crystallographic copy (2026-09-07)

**Scope.** Replace "symmetry-corrected in-place heavy-atom RMSD to the single deposited instance
`<ID>_ligand.sdf`" by "the same RMSD to the nearest of all deposited copies in `<ID>_ligands.sdf`",
per pose, with every crystal-referenced metric of that pose derived from the copy it is nearest to.
The single-instance value is kept as a secondary column and becomes the appendix sensitivity arm.
**Status.** PLAN ONLY. Nothing in the thesis, the scripts or the harness has been changed.
**Evidence base.** M2 re-validation (`thesis_latex/M2_REVALIDATION_2026-09-07.md`), the probe
(`Scripts/Analysis/nearest_copy_probe/`), and a nine-dimension read-only inventory run on 2026-09-07
(470 file:line sites; dossier and two-convention CSVs archived under
`Scripts/Analysis/nearest_copy_probe/wf_plan_2026-09-07/`). Every number below was recomputed from
the canonical `per_pose_metrics.csv` with the thesis's own loaders, `symmetry_rmsd`, `posebusters_rmsd`
and `stats_utils`; each reference-convention value was reproduced to the printed digit before the
nearest-copy counterpart was taken.
**Verification status, stated plainly.** The inventory sweeps completed for all dimensions except
"numbers and selection", which was recomputed by hand afterwards (Section 1). Of the eighteen planned
adversarial audits only the main-body evidence audit ran before the session quota was exhausted; it
refuted nothing, corrected five sites and added five (all folded in below). The other seven dimensions
are inventoried but not independently audited. Phase 0 therefore includes the audit pass.

---

## 0. Decision record

Each decision names the options considered and the one this plan adopts. All are upstream-consistent
with PoseBusters' own redock path, which is the convention the source benchmark published under.

| # | Decision | Adopted | Rejected alternatives and why |
|---|---|---|---|
| D1 | Definition of "nearest" | Per pose, minimum over all records of `_ligands.sdf` of the thesis's own in-place `symmetry_rmsd`, NaN as +inf, record 0 (the reference) as tie-breaker | PoseBusters `check_rmsd` argmin on a combined multi-conformer `mol_true` gives the same ≤2 Å verdict on 414/414 rank-1 poses tested but symmetrises terminal groups, which the thesis discloses it does not (`body_main_short.tex:113`). Centroid-nearest can disagree with RMSD-nearest only for six adjacent-copy ids and RMSD is the gate. |
| D2 | Which copy drives centroid, Kabsch form, torsions, contact and PLIF recovery | The single argmin-RMSD copy j* for every metric of that pose (PoseBusters `rmsd.py:61-70` rule: `choose_by="rmsd"`, kabsch and centroid follow) | Per-metric minima break `placement² = inplace² − form²` (`_form_components`, hub :7043-7051) and `autodock_exhaustiveness_returns.py:570` asserts `bestfit_rmsd == pb_kabsch_rmsd`. Leaving contacts on the reference makes a pose on chain B's copy score ~0 recovery against chain A while being counted near-native. Measured cost of "follow": Kabsch median shift 0.001 Å, form verdict changes on 0/2/0 of the 38/29/2 flipped rank-1 poses; 4 Å reach differs from the per-metric minimum by ≤2 complexes (239 vs 241 AutoDock*). |
| D3 | Four alternates that do not share the reference pocket (7A9E_R4W, 7TUO_KL9, 7VKZ_NOJ, 7Z1Q_NIO) | Include every deposited instance, as PoseBusters and the dataset README do, and name the four in the appendix | A Jaccard ≥0.6 restriction adds a parameter for zero effect: no headline arm gains a complex from these four at d = 1, 15 or 30. 7TUO's alternate is a genuine second site 26 Å away that no pose approaches within 17 Å; the other three are side-by-side fragment copies 4.9-6.2 Å from the reference with minimum alternate RMSD 2.1-2.6 Å. |
| D4 | Metal-adjacent stratum (78 of 303) | Keep as a complex-level property of the reference instance; verify once that membership agrees when tested against every copy; register its currently unregistered generator | Per-pose or any-copy definitions turn a receptor-preparation property into a pose property. Only the per-stratum recovery rates are recomputed. |
| D5 | Crystal-cluster reach when several crystal sites exist | A cluster reaches if its centroid lies within 4 Å of ANY copy centroid; the crystal-closest cluster is defined by the nearest copy; the six ids whose adjacent copy sits 4.9-7.8 Å from the reference (inside the 8 Å clustering radius) count once | Reference-only is the status quo being replaced. Per-site multi-recovery is a new analysis the thesis does not need. |
| D6 | Column naming in `per_pose_metrics.csv` | Keep `rmsd`, `pb_rmsd`, `pb_kabsch_rmsd`, `bestfit_rmsd`, `centroid_dist`, `contact_recovery`, `plif_recovery` as the primary names carrying the nearest-copy value; add `rmsd_ref_instance`, `centroid_dist_ref_instance`, `pb_rmsd_ref_instance`, `bestfit_rmsd_ref_instance`, `n_copies`, `ref_copy_index`, `nearest_copy_index`, `nearest_copy_is_ref`, `reference_convention` | Renaming the primary forces edits in 14 consumer scripts, two test files and seven gate lines in `thesis_assertions.py`; keeping the names leaves all of them untouched and keeps the yaml column roles true. |
| D7 | Where the new cache lives | Back up the current report directory as `pose_comparison_report.bak-instance-<date>` (repo precedent `.bak-preexh128rerun-20260830`) and write the rebuild to the canonical path | A new directory name would force path edits in the builder, REGENERATE.md and every notebook cell that reads `BENCH_REPORT`. |
| D8 | CLI and cache safety | Add `--reference-convention {instance,nearest}` to the hub; land it with default `instance`, prove the rebuilt cache reproduces today's table bit-for-bit, then flip the default to `nearest` in the same commit that edits the stage command. Bump the cache schema 5→6, add `rmsd_ref_instance` to `_CACHE_REQUIRED_COLS`, record the convention in the manifest and make `--reuse-cache` refuse a cache built under the other convention | A warn-only path leaves the documented `--reuse-cache --collapse-plots-only` re-render convention-blind, which is the Gen A/B trap in a new guise. |
| D9 | Sensitivity arm | The single-instance convention becomes the appendix sensitivity table (the M2 table with its columns inverted) plus a single-copy stratum row (165 complexes: 54.5 / 50.9 / 32.7 % rank-1) | No sensitivity arm would discard the conservative bound M2 valued. |
| D10 | Pre-specified confirmatory cell (App. C :230, +65.6 complexes per hour, 95 % +32.1 to +137.2) | Keep the confirmatory reading on the single-instance column and print the nearest-copy value beside it as descriptive | Re-declaring a pre-specified cell after seeing the data is what the paragraph says it avoids. |
| D11 | Methods wording | Redefine the endpoint at `body_main_short.tex:109-115` and cite the external convention (dataset README "All instances of the ligand of interest"; `redock.yml:307 load_all: True`; paper Section 2.3 "closest crystallographic ligand") | Describing it as a repair for the blind box misstates it: the paper scored every method this way and its Vina anchor was inert only because the 25 Å cube excluded the copies. |
| D12 | Variant selection | Re-run the top-15 triple gate under the new convention and report the outcome. Measured: AutoDock* 160→175, DiffDock smina 133→144 against gnina 132→143 (1,031 vs 1,022 qualifying poses, seven discordant, three gnina-only), EquiBind gnina 55→57 against smina 46→48. All three selections hold | Freezing the selection under the old convention invites the same examiner objection twice. |
| D13 | Start-conformer difficulty descriptor (App. F :782, Fig 37) | Keep on the single instance and say it is a start-conformer property, not a docked-pose metric | Copies differ only conformationally; regenerating image40/image42 buys nothing. |

**Data facts the decisions rest on** (308 ids, `Data/PoseBuster Benchmark Set`): 143 multi-copy ids
(111×2, 7×3, 22×4, one each ×5, ×6, ×12), 138 in the 303, 211 alternate copies. The reference is record 0
of `_ligands.sdf` in 308/308. All 211 alternates have the same heavy-atom count, formula, element sequence,
formal charge and all-atom count as the reference and substructure-match it in both directions, so no
copy can silently drop out of the minimum. Eight nucleotide alternates differ only in perceived stereo
tags (irrelevant with `useChirality=False`). Alternate centroid separation from the reference: none within
4 Å, six within 8 Å (7A9E_R4W 4.9, 7CD9_FVR 5.0, 7JY3_VUD 5.8, 7VKZ_NOJ 6.0, 7Z1Q_NIO 6.2, 6YYO_Q1K 7.8),
203 beyond 12 Å, median 36.0 Å. Seven ids carry the second copy on the same chain id. All five dropped
complexes (7B2C_TP7, 7D6O_MTE, 7FRX_O88, 7M31_TDR, 8F4J_PHO) are multi-copy (2, 2, 4, 4, 4 records).

**External convention, confirmed.** Chem. Sci. 2024 full text (PMC10901501), Section 2.3: "PoseBusters
calculates the minimum heavy-atom symmetry-aware root-mean-square deviation (RMSD) between the predicted
ligand binding mode and the closest crystallographic ligand". Installed `posebusters/config/redock.yml:307`
sets `load_all: True` on `mol_true` ("important to set if there are multiple ligands that could be correct");
`modules/rmsd.py:52-70` takes argmin over conformers and reports kabsch and centroid from that copy; the
v0.2.2 paper-era code already did so. Dataset `README.txt:55-59` documents `_ligands.sdf` as "All instances
of the ligand of interest" and `_ligand.sdf` as "One of the instances ... marks the binding site for those
docking methods that require a binding site". The thesis footnote at `body_appendix_short.tex:662` calling
`_ligand.sdf` "the ground-truth answer pose" is therefore not what the source documents. Recomputed from the
Zenodo results CSV on the 308 ids: Vina 59.7 % RMSD ≤ 2 Å, DiffDock 38.0 %; under that same nearest-copy
scoring the published DiffDock rate is 31.5 % on the 143 multi-copy ids against 43.6 % on the 165 single-copy
ids. One wording conflict to note if quoted: the paper says GetBestRMS (superposing) while every code tag uses
CalcRMS for the in-place value.

---

## 1. What changes numerically

Convention "nearest" below means D1 + D2 (form and centroid follow the in-place copy). n = 303 throughout.
"Derived" values come from the probe CSVs; "recompute" means the quantity needs the rebuilt pipeline.

### 1.1 Headline recovery (PB-valid AND ≤ 2 Å, existence gate; Table 5, Abstract, Kurzfassung, Conclusions)

| Arm | d | Reference | Nearest | Gain |
|---|---|---|---|---|
| AutoDock* (exh128 + gnina) | 1 | 111 (36.6 %) | 149 (49.2 %) | +38 |
| DiffDock* (smina) | 1 | 103 (34.0 %) | 132 (43.6 %) | +29 |
| EquiBind* (unguided gnina) | 1 | 55 (18.2 %) | 57 (18.8 %) | +2 |
| AutoDock* | 15 | 198 (65.3 %) | 213 (70.3 %) | +15 |
| DiffDock* | 15 | 167 (55.1 %) | 179 (59.1 %) | +12 |
| EquiBind* | 15 | 80 (26.4 %) | 83 (27.4 %) | +3 |
| AutoDock* | 30 | 199 (65.7 %) | 214 (70.6 %) | +15 |
| DiffDock* | 30 | 169 (55.8 %) | 183 (60.4 %) | +14 |
| EquiBind* | 30 | 80 (26.4 %) | 83 (27.4 %) | +3 |

Rank-1 to top-15 gains: +87 → +64 (+28.7 → +21.1 pp, Wilson 16.9-26.1), +64 → +47 (+21.1 → +15.5, 11.9-20.0),
+25 → +26 (+8.3 → +8.6, 5.9-12.3). Top-15 to top-30: +1 / +4 / +0 (DiffDock b = 4, c = 0, p = 0.125, no star),
so "at most 0.7 points" at `:404` becomes "at most 1.3". Near-native without the validity gate at rank-1:
113 / 107 / 57 → 151 / 137 / 59, so the validity gate removes at most FIVE complexes (DiffDock rank-1), not four
(`:324`, `:777`, App. H :1194). Rank-1 failures with no qualifying pose anywhere: 54.2 / 67.0 / 89.9 % →
57.8 / 70.2 / 89.4 % ("about six in ten" at `:406`, `:788`, `:857`, Abstract `:63`, Kurzfassung `:95` becomes
"six to seven in ten"). Raw exh128: 31.0 / 61.1 / 66.7 % → 42.9 / 65.3 / 71.6 %; the raw rung still leads the
rescored one over the full pool, 217 against 214 (was 202 against 199).

### 1.2 AutoDock* against DiffDock* (paired; `:382-388`, `:834`, App. H.9 :1256-1258, glossary)

| d | Convention | Margin | Newcombe 95 % | Discordant | Exact McNemar p |
|---|---|---|---|---|---|
| 1 | reference | +8 (+2.6 pp) | −4.2 to +9.4 | 60 / 52 (112) | 0.509 |
| 1 | nearest | +17 (+5.6 pp) | −1.7 to +12.8 | 73 / 56 (129) | 0.159 |
| 5 | both | +13.2 pp | +5.7 to +20.4 | 88 / 48 | 7.6e-4 (identical under both conventions) |
| 15 | reference | +10.2 pp | +2.8 to +17.5 | 82 / 51 | 0.009 |
| 15 | nearest | +11.2 pp | +3.9 to +18.4 | 82 / 48 | 0.004 |
| 30 | reference | +9.9 pp | +2.5 to +17.1 | 81 / 51 | 0.011 |
| 30 | nearest | +10.2 pp | +3.0 to +17.3 | 79 / 48 | 0.008 |

Power statement (`:386`, App. H.9, glossary :191): discordant 112 → 129, minimum detectable difference
9.8 → 10.5 pp at 80 % power, observed power 0.12 → 0.32, complexes needed for the observed effect 4,162 → 1,062.
Equivalence (App. H.9 :1258, yaml :1042): the 90 % interval moves from −3.1..+8.3 ("about eight points") to
−0.5..+11.7, so the rank-1 contrast is no longer equivalent within 10 points, only within 12. Cochran's Q at
rank-1 p ≈ 2e-17 (still < 1e-7). Metal-free versus metal-adjacent stratum margins (+6.7 on 225, +20.5 on 78 at
`:868`, `App. C :165`) recompute; the stratum itself stays 78 (D4).

### 1.3 Table 1 pose-level columns (`body_main_short.tex:286-304`, yaml table_1 :71-80)

Validity columns are unchanged. Near-native, Form and Combined (D2 "follow"):

| Row | Near-native cplx / poses / % | Form cplx / poses / % | Combined cplx / poses / % |
|---|---|---|---|
| AutoDock (raw exh128) | 202 / 307 / 3.4 → 217 / 437 / 4.9 | 225 / 2,848 / 31.7 → 225 / 2,829 / 31.5 | 164 / 221 / 2.5 → 178 / 324 / 3.6 |
| AutoDock + gnina | 202 / 319 / 3.6 → 217 / 451 / 5.0 | 227 / 2,836 / 31.6 → 227 / 2,823 / 31.5 | 160 / 221 / 2.5 → 175 / 325 / 3.6 |
| DiffDock (raw) | 152 / 1,683 / 18.7 → 161 / 2,118 / 23.5 | 225 / 4,012 / 44.5 → 224 / 4,016 / 44.6 | 83 / 722 / 8.0 → 92 / 869 / 9.6 |
| DiffDock + smina | 173 / 1,715 / 19.1 → 186 / 2,169 / 24.1 | 237 / 3,903 / 43.4 → 236 / 3,887 / 43.2 | 138 / 1,227 / 13.6 → 148 / 1,516 / 16.9 |
| DiffDock + gnina | 176 / 1,716 / 19.1 → 186 / 2,171 / 24.1 | 235 / 3,873 / 43.1 → 235 / 3,854 / 42.9 | 138 / 1,219 / 13.6 → 149 / 1,513 / 16.8 |
| EquiBind (raw, unguided) | 19 / 103 / 1.1 → unchanged | 185 / 2,479 / 27.3 → 186 / 2,487 / 27.4 | 2 / 21 / 0.2 → unchanged |
| EquiBind + smina | 64 / 453 / 5.0 → 67 / 465 / 5.1 | 191 / 2,610 / 28.7 → 192 / 2,617 / 28.8 | 48 / 249 / 2.7 → 50 / 253 / 2.8 |
| EquiBind + gnina | 83 / 498 / 5.5 → 86 / 515 / 5.7 | 200 / 2,661 / 29.3 → 199 / 2,661 / 29.3 | 55 / 282 / 3.1 → 57 / 289 / 3.2 |
| EquiBind fpocket raw / P2Rank raw | 0 → 0 / 0 → 1 (2 poses) | 180 / 2,449 → 179 / 2,441; 181 / 2,336 → 180 / 2,339 | 0 → 0 |
| Uni-Dock2 | 81 / 117 / 1.4 → 95 / 162 / 1.9 | 196 / 2,729 / 31.7 → 196 / 2,715 / 31.6 | 65 / 85 / 1.0 → 82 / 125 / 1.5 |

The Form block can LOSE a complex under "follow" (AutoDock raw none, DiffDock raw 225 → 224, DiffDock + smina
237 → 236, EquiBind + gnina 200 → 199), because the copy a pose sits on may be a slightly different conformer.
The M2 statement "no complex loses" is true for near-nativeness and must not be carried over to Form.

### 1.4 Form axis (Kabsch; `:396`, `:404`, Fig 3, Table 26 Kabsch rows, glossary :151-159)

Unchanged at the printed precision: 74.3 / 74.6 % at 1 Å best-of-top-30, 95.4 / 97.0 % at 2 Å, rank-1
52.5 / 52.1 / 33.3 %. EquiBind best-of-top-30 54.1 → 53.8 % (gain 20.8 → 20.5 points). 26 of the 90 Table 26
Kabsch cells move by 0.3-0.7 pp. Form McNemar AutoDock vs DiffDock: 47/46, 30/31, 25/26 → 48/47, 30/32, 25/26,
Holm p stays 1.000, so "indistinguishable on form" (`:790`) holds.

### 1.5 Variant selection and the exhaustiveness ladder (`:191-199`, `:781`, `:866`, App. C :178-248, Table 8)

Top-15 triple gate (PB-valid & ≤ 2 Å & Kabsch ≤ 1 Å): AutoDock* 160 → 175; DiffDock smina 133 → 144 against
gnina 132 → 143 (poses 834 / 828 → 1,031 / 1,022; discordant nine → seven, gnina-only four → three); EquiBind
gnina 55 → 57 against smina 46 → 48. Double gate top-15: smina 167 → 179 against gnina 167 → 177, so "exactly
level without the form term" (`:781`) becomes "smina ahead". Raw exh128 against its rescored sibling on the
triple gate over the pool: 164 vs 160 → 178 vs 175.

Ladder (PB-valid & ≤ 2 Å; d = 1 / 15 / 30; reference → nearest):

| Rung | Reference | Nearest |
|---|---|---|
| Exhaustiveness 18 | 72 / 121 / 125 | 92 / 138 / 143 |
| Exhaustiveness 32 | 85 / 148 / 154 | 109 / 161 / 170 |
| Exhaustiveness 64 | 94 / 175 / 182 | 123 / 190 / 201 |
| Exhaustiveness 92 | 91 / 182 / 194 | 123 / 194 / 209 |
| Exhaustiveness 128 | 94 / 185 / 202 | 130 / 198 / 217 |
| 32 + gnina | 100 / 154 / 154 | 129 / 170 / 170 |
| 64 + gnina | 110 / 182 / 182 | 141 / 200 / 200 |
| 128 + gnina (selected) | 111 / 198 / 199 | 149 / 213 / 214 |

The selected rung still beats each of the seven alternatives at top-15 on exact McNemar (uncorrected p from
1.1e-20 against exh18 to 2.6e-3 against its own raw search and 4.4e-3 against 64 + gnina; Holm to recompute)
and still leads at every depth from 3 to 22. The rank-1 "tie with exhaustiveness 64 rescored, 111 against
110" (`App. C :182`) becomes a lead, 149 against 141. "Adding rescoring moves top-15 from 185 to 198, gaining
17 and losing 4" (`:184`) becomes 198 to 213, gaining 19 and losing 4. "Raising exhaustiveness 32 → 128 moves
top-15 from 154 to 198" becomes 170 to 213. The AutoDock gnina-versus-raw rank-1 step (`:843`, App. H :945):
31.4 → 37.3 % becomes 43.2 → 49.8 % on the top-k table basis, discordant 39/20, uncorrected p 0.018 (was 0.040);
whether it clears Holm is recompute. Winner's-curse estimate (`App. C :248`, `:866` "3.0 points", "0.4-0.8")
recompute on the new ladder.

### 1.6 Placement, cohort trims and the H.4 decomposition (`:408-414`, Fig 4, Table 27, App. H :1186-1192, :1375-1385)

Near-site (8 Å) PB-valid top-15 pool 6,169 → 7,565 poses; clipped-beyond-5 Å 2,329 → 2,854 (the figure
script's clipping rule reproduces 2,337 on the reference basis, so take the new value from the script, not from
this table). Rank-1 coverage after the trim 62.7 / 54.1 / 47.2 % → 82.2 / 73.9 / 51.2 %. AutoDock in-place
median 1.490 → 1.376 Å at rank-1 and 5.043 → 4.922 Å at top-15, drift +3.55 Å unchanged; DiffDock drift
+0.46 → +0.53; EquiBind +1.09 → +0.96. EquiBind rank-1 near-site median 2.9 → 3.3 Å (`:414`). All 117 Table 27
cells change. Near-native strata 319 / 1,715 / 498 → 451 / 2,169 / 515 poses; all-check validity within them
97.5 / 94.6 / 93.6 → 98.0 / 95.3 / 93.8 %; the per-check intermolecular rates (14.6 → 5.4 %, 37.0 → 6.2 %) need
the per-check PoseBusters columns and are recompute. Depth-gain decomposition (Table 25): [88, 2, 1, 87] /
[63, 1, 2, 64] / [26, 2, 1, 25] → [65, 2, 1, 64] / [45, 0, 2, 47] / [27, 2, 1, 26].

### 1.7 Crystal-site reach and clusters (`:158`, `:427-476`, `:792`, Table 6, Figs 5, 6, 38, App. H :1392-1397)

Pose-level rank-1 4 Å reach 183 / 172 / 137 → 239 / 225 / 143 under D2 (241 / 225 / 144 under a per-metric
minimum). Any-pose reach per tool 88.4 → 92.4 % (AutoDock*), 85.1 → 90.1 % (DiffDock*); pooled three-tool
oracle 98.0 → 99.7 % (302/303). The printed sentence at `:476` (97.7 / 88.1 / 84.5 %) is already one to two
complexes off the canonical file; trace its generating script before writing 99.7 in. Cluster-level reach,
co-reach φ, "sixteen complexes with no qualifying cluster", "18 contribute no pose", precision-at-one 56.8 /
55.4 / 55.8 % and the 95 % oracle ceiling are recompute under D5.

### 1.8 Native-interaction recovery (`:165`, `:481-532`, Fig 7, Table 7, App. H :543-551, :1405-1413, Fig 39)

Everything is recompute after per-copy crystal fingerprints exist. Direction: 105 / 119 / 11 top-5 PB-valid
poses of AutoDock* / DiffDock* / EquiBind* sit within 2 Å of an alternate copy and today score Jaccard ≈ 0
because the fingerprint keys carry the reference chain. Reconstruction on the probe CSV gives crude F1
0.565 / 0.529 / 0.471 → about 0.532 / 0.527 / 0.513 and the 4,253-pose audit cohort 7.57 Å median, 21.3 % within
2 Å → recompute. The poses beyond 20 Å with no metal contact (750, F1 0.0003) shrink because alternate-copy
poses leave that band.

### 1.9 Cost per qualifying pose (`:708-760`, `:817`, `:848`, `:850`, Table 8, Figs 15, 16, Abstract `:63`, Kurzfassung `:98`)

Qualifying poses 311 / 1,622 / 466 → 442 / 2,068 / 483 (Vina-only 302 → 431); qualifying complexes
199 / 169 / 80 → 214 / 183 / 83; the all-three shared set 53 → 54. The AutoDock median succeeding complex now
yields two qualifying poses instead of one, so every charged / CPU / GPU median, IQR, Friedman, ratio and the
"most expensive per qualifying pose" ordering (`:744`, `:754`, `:850`, Abstract, Kurzfassung) must be re-derived
from the timing files. Campaign-per-success values (271 / 479 / 38 s) recompute over 214 / 183 / 83.

### 1.10 Intention-to-treat check (`App. I :1072`)

Re-scored on 2026-09-07 with the probe engine over the 150 AutoDock exh128 + gnina poses of the five dropped
complexes, ranked by CNNaffinity: rank-1 none within 2 Å (closest 6.3 Å, 7FRX_O88 to an alternate copy, was
80.7 Å to the reference); top-15 none within 2 Å but the closest pose is now 2.1 Å (7FRX_O88) instead of 6.4 Å;
over all thirty poses 8F4J_PHO reaches 1.5 Å of an alternate copy at a rank beyond 15. Base ITT counts become
149 / 308 and 213 / 308 against 132 and 179. The 6.4 Å and 150-pose figures have no registered generator today;
Phase 2 registers one.

### 1.11 What the Abstract and Kurzfassung print

Abstract `:61-63`: five number swaps (36.6 / 34.0 / 18.2 → 49.2 / 43.6 / 18.8; +2.6 → +5.6; −4.2..+9.4 →
−1.7..+12.8; 65.3 / 55.1 → 70.3 / 59.1; 80 of 303 → 83), the "six in ten" clause, the cost-ordering clause
pending 1.9, and one optional convention clause. Kurzfassung `:91-100`: 36,6 / 34,0 / 18,2 → 49,2 / 43,6 / 18,8
(same length), 2,6 → 5,6, "−4,2 bis +9,4" → "−1,7 bis +12,8" (+1 character on a full page), 65,3 / 55,1 →
70,3 / 59,1, 80 → 83, "sechs von zehn" → "sechs bis sieben von zehn" (+11 characters, needs a matching cut on
the same page), the cost clause at `:98-99` pending 1.9. No German convention clause (page is full).

---

## 2. Step-by-step plan

Effort classes: T trivial (minutes), S small (< 1 h), M medium (half a day), L large (1-2 days).
Every step names its verification check. Order is a dependency order.

### Phase 0. Decide, audit, freeze (M, 0.5 day)

0.1 Adopt the decision record (Section 0) in a dated note under `thesis_latex/`. Check: every later step cites a D-number.
0.2 Run the seven missing adversarial audits over the archived inventories (`wf_plan_2026-09-07/dossier/`), or re-run the workflow with `resumeFromRunId wf_9b6bad77-5da` once quota allows. Check: refuted sites removed from Sections 3-4 before any edit.
0.3 Freeze baselines: copy `posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report/` to `pose_comparison_report.bak-instance-<date>`; copy `thesis_latex/media/media/` and `Scripts/Analysis/thesis_expected_values.yaml` beside it; note `git rev-parse HEAD`. Check: `md5sum` manifest of the backup written next to it.
0.4 Trace two pre-existing discrepancies before they are overwritten: the `:476` reach sentence (97.7 / 88.1 / 84.5 vs canonical 98.0 / 88.4 / 85.1) and the `:408` clip count (2,329 vs 2,337 by the stated rule). Check: generating script and pose basis named in the note.

### Phase 1. Hub code, `Scripts/Analysis/posebusters_pose_comparison.py` (L, 1 day)

1.1 `:1149` add `load_all_mols(sdf_path)` beside `load_first_mol`, same sanitize fallback per record (port of `nearest_copy_rescoring.py:19-28`).
1.2 `:1520-1528` load `copies_sdf = cdir / f"{pdb_id}_ligands.sdf"`, fall back to the single file if absent; `copies_h = [RemoveHs(reassign_template(m, crystal_h) or m) for m in load_all_mols(...)]`; `ref_copy_index` by centroid identity with `crystal_h` (0 on this dataset).
1.3 `:1531-1532` compute `native_contacts[j]`, `n_native[j]` for every copy.
1.4 `:1560-1562` and `compute_plif_recovery` (`:1383-1415`): pass all copies as reference ligands in the one ProLIF run, return an (n_poses × n_copies) matrix, select column j* per pose.
1.5 `:1566` `r = [symmetry_rmsd(pose, c) for c in copies_h]`, `j_star = nanargmin(r)` with NaN → +inf and record 0 as tie-breaker; `rmsd = r[j_star]`, `rmsd_ref_instance = r[ref_copy_index]`.
1.6 `:1567`, `:1572`, `:1581`, `:1584`, `:1585`, `:1588` evaluate `centroid_distance`, `posebusters_rmsd`, `rigid_body_fit`, `bestfit`, `torsion_metrics` and contact recovery at `copies_h[j_star]` (D2). Keep `*_ref_instance` twins for rmsd, centroid, pb_rmsd, bestfit.
1.7 `:1459-1511` `PoseRecord`: append `reference_convention`, `n_copies`, `ref_copy_index`, `nearest_copy_index`, `nearest_copy_is_ref`, the four `*_ref_instance` floats, at the END so existing column positions hold. `:1592` wire them.
1.8 `:1513-1514` and the work-tuple builder at `:13494`: thread `reference_convention` through `process_pair`; under `instance` set j* = ref_copy_index so the old table is reproduced exactly.
1.9 `:13346` argparse `--reference-convention {instance,nearest}` (default `instance` until 1.12); `:13224` schema 5 → 6 and add the convention to `_per_pose_signature`; `:13187` add `rmsd_ref_instance`, `nearest_copy_index`, `n_copies`, `reference_convention` to `_CACHE_REQUIRED_COLS`; `:13249-13267` make `--reuse-cache` read the manifest and refuse a cache whose convention differs from the flag (D8).
1.10 `:9313-9330` mechanism exemplars: load record `nearest_copy_index` of `_ligands.sdf` instead of `_ligand.sdf`. `:46-55`, `:661`, `:723-732` docstrings and comments name the convention. Sidecar text writers (`:4827`, `:5625`, `:5690`, `:6199`, `:6402`, `:6508`, `:10667`, `:12043`) print the convention string from the manifest.
1.11 Regression: rebuild with `--reference-convention instance --force` into a scratch out-dir; assert `rmsd`, `pb_rmsd`, `bestfit_rmsd`, `centroid_dist` equal the frozen table to 1e-9 Å on all 241,914 rows and `oracle_summary.csv` is byte-identical apart from the new columns. Add a unit test (`Scripts/Docking/tests/test_downstream_variants.py:601-624` style) that a single-copy id yields `rmsd == rmsd_ref_instance` and `nearest_copy_is_ref` True. Check: both pass.
1.12 Flip the default to `nearest`; record the multi-conformer PoseBusters result once as `pb_rmsd_loadall` during the first rebuild to show parity with `redock.yml:307`, then drop it if identical (expected: ≤ 2 Å verdicts agree 100 %, values differ only by terminal-group symmetrisation).
1.13 Pocket localisation (`:10998-11044`, `POCKET_CENTROID_CUTOFF` `:733`): keep `centroid_dist` at j* (D2); document that "validated pocket" now means the pose's nearest copy.

### Phase 2. Downstream scripts with their own crystal logic (L, 1.5 days)

2.1 `Scripts/Analysis/pose_cluster_crystal_pocket_report.py`: `:921` `crystal_centroid` returns the list of copy centroids; `:1109-1111` a cluster is correct if within `--match-thr` of ANY copy (D5) with the six adjacent-copy ids counted once; `:562` `_precision_at_1`, `:963` `cdist`, `:1020`, `:4064` min over copies; `:3697` and `:980` follow the hub column automatically. Check: on the 165 single-copy complexes every reach / co-reach / purity value equals the frozen run.
2.2 `Scripts/Analysis/run_pandamap.py:963`: emit one crystal task per record of `_ligands.sdf` (`pose_name .../crystal_k`), profiled against the canonical `hetatm_cleaned/<ID>_protein.pdb`. OPEN CHECK FIRST: confirm that receptor retains every chain so alternate copies can be profiled; if it does not, the fallback is to restrict Fig 7 / Table 7 / Fig 39 to poses whose nearest copy is the reference and disclose the cohort.
2.3 `Scripts/Analysis/pandamap_interaction_report.py`: `:1181` `native_recovery`, `:1716` `build_recovery_detail` and the per-protein natives join `per_pose_metrics` on the pose key to obtain `nearest_copy_index` and select `crystal_k`; `:2251-2255` RMSD join is automatic; `:599` descriptors are copy-invariant. `Scripts/Analysis/interaction_pose_basis_audit.py:206` automatic. Check: F1 on single-copy complexes equals the frozen run.
2.4 `Scripts/Analysis/autodock_exhaustiveness_returns.py:570` (asserts `bestfit_rmsd == pb_kabsch_rmsd`): unchanged under D2; run it as the first consumer to prove the coupling holds.
2.5 Pure re-readers need no code change and re-run automatically: `posebusters_validity_report.py:748`, `docking_effort_comparison.py:256`, `topk_recovery_gnina_arm.py:130`, `thesis_endpoint_diagnostics.py:77`, `optimization_benefit_stats.py:323`, `rank_success_analysis.py:139/:286`, `filmstrip_rank1_top5_top15.py:133`, `ligand_/receptor_docking_difficulty`, `diffdock_gnina_rerank_analysis`, `rebust_hydrogen_normalised.py:254-267`. `boxed308_raw_arm_metrics.py:37` optional parity branch (0 flips on the boxed run).
2.6 Register two currently unregistered generators as harness stages: the 78-complex metal-adjacent stratum (no script under `Scripts/` produces it) and the intention-to-treat re-score of the five dropped complexes (Section 1.10; the probe's `work()` over `Dockings/vina_results_full_protein_vina_scoring_mgltools_exh128/<ID>/mgl_tools/docking/optimized_gnina/`).
2.7 Register `bench_reference_convention` (section 3, BITEXACT, CHEAP, needs `bench_pose_comparison`) writing `reference_convention_sensitivity.csv` from the `*_ref_instance` columns, the inverted M2 table plus the single-copy stratum row (D9). Move `nearest_copy_probe/nearest_copy_rescoring.py` under it as the stage script or retire it.

### Phase 3. Rebuild (compute ≈ 1-2 h wall, attention S)

3.1 `bench_pose_comparison --force` with the stage command at `_build_reproduction_notebook.py:386-400` plus `--reference-convention nearest` (cost class moderate, "minutes to about an hour" on 24 workers; the multi-copy work adds RMSD, Kabsch, torsion and contact evaluations for up to 12 copies on 138 complexes).
3.2 Then, in order: `bench_validity_report`, `bench_filmstrip`, `bench_topk`, `bench_endpoint_diagnostics`, `bench_optimization_benefit`, `bench_exhaustiveness_returns`, `bench_clusters`, `bench_effort_charged`, `bench_effort_elapsed`, PandaMap crystal tasks then `bench_pandamap_report`, `bench_interaction_audit`, `bench_diffdock_rerank`, the two new stages of 2.6 and `bench_reference_convention` (14 of the 19 section-3 stages plus 3 new). Not re-run: `bench_figure1`, all PoseBusters, docking, Orai and dataset stages.
3.3 Re-render collapsed figures with `--reuse-cache --collapse-plots-only` only after 1.9 is in place (the flag pair is REQUIRED, and the manifest check now protects it).

### Phase 4. Selection re-check (S)

4.1 Confirm from the rebuilt sidecars: exh128 + gnina still selected at top-15 (expected 213, beats all seven); DiffDock smina still carried (144 vs 143 triple, 179 vs 177 double); EquiBind unguided gnina still the winner under `BEST_VARIANT_METRIC` (`:2136`). If any flips, stop and re-plan the affected arm's text before Phase 5.
4.2 Confirm `_TOPK_TEST_PAIRS` (`:9764-9797`, AutoDock gnina-opt vs DiffDock raw / gnina-opt / EquiBind gnina-opt) significance pattern; derived expectation: vs DiffDock raw 77 / 58 at k = 1 (uncorrected 0.12, "genuine tie" wording holds), significant from k = 5.

### Phase 5. Thesis text (L, 2-3 days; short build only)

5.1 Methods: `body_main_short.tex:109` redefinition (D11 draft below), `:113` recompute the "adds three recovered complexes" sentence, `:115`, `:119`, `:124-133` form definition names the copy rule, `:158`, `:165` interaction reference, `:170` cost definition.
5.2 Results: every site in Section 3 (Tables 1, 5, 6, 7, 8; `:191-199`, `:324-326`, `:380-414`, `:427-435`, `:455-476`, `:481-532`, `:708-760`); captions of Figs 2-7, 10, 11, 15, 16.
5.3 Discussion, Conclusions, Limitations: `:768`, `:777-781`, `:788-794`, `:806`, `:817`, `:824`, `:834-836`, `:843-850`, `:857`, `:866-870`; add the one-sentence convention disclosure at `:868` (143 / 308, 138 / 303, median 36.0 Å, source benchmark scored the same way, single-instance kept as sensitivity).
5.4 Appendix: `body_appendix_short.tex:77` (Vinardo, needs probe run on the vinardo arms), `:165`, `:173-184`, `:218-248` ladder and figures 18-21 captions, `:267`, `:276`, `:487` optional clause, `:543-551` interaction, `:662-664` rewrite (draft below), `:869`, `:891-918`, `:945-962`, `:990-997`, `:1024-1029`, `:1054-1060`, `:1072` ITT, `:1170-1194`, `:1201` regime clause, `:1205`, `:1245-1258`, `:1263`, `:1301-1338`, `:1375-1385`, `:1392-1413`; new sensitivity table (D9) beside `:1072`.
5.5 Abstract and Kurzfassung per Section 1.11; pay the +12 German characters with a cut on the same page.
5.6 House style: no mid-sentence semicolons, colons or clause-dashes; bare tool names after the dominant-variant declaration; "ranking depth", "rank-1", "top-15".

Draft for `:109` (D11): "Near-nativeness is the symmetry-corrected heavy-atom RMSD between a docked pose and the
nearest of all crystallographic copies of the ligand deposited with the structure, computed in place without
superposition. The copy with the lowest RMSD is chosen for each pose and supplies the centroid, form, contact
and interaction references for that pose. This is the convention of the source benchmark, whose validity
battery loads every deposited instance and reports the lowest RMSD. The 165 single-copy complexes are
unaffected. The single-instance value is retained and reported as a sensitivity arm in Appendix~\ref{...}."

Draft for `:664`: "In 143 of the 308 entries the deposited file holds more than one copy of the ligand, 138 of
them among the 303 analysed complexes, and the nearest other copy lies a median 36.0 Å from the reference.
Every pose is scored against the nearest deposited copy, and the form, centroid and interaction references
follow that copy. The source benchmark scored every method this way. Its 25 Å cube centred on the reference
kept the other copies outside Vina's search, whereas the whole-protein box used here contains every one of
them, so the convention is live for all three tools in this study. Four alternates do not share the
reference pocket (7A9E_R4W, 7TUO_KL9, 7VKZ_NOJ, 7Z1Q_NIO) and no reported pose recovers any of them."

### Phase 6. Figures and media (M, 0.5 day)

6.1 Regenerate the 14 convention-dependent figures: image3, image4 (hub 18_topn curves), image5 (filmstrip),
image6, image7, image8 (cluster report), image9, image11 (PandaMap), image24, image25 (effort), image44-47
(exhaustiveness). Hand-copy each into `thesis_latex/media/media/` and confirm by md5 against the producing
report file, never by timestamp. The other 25 figures are unchanged (Fig 1 validity-only, 12 Orai figures,
flow chart, two PyMOL renders, 10 dataset-descriptor figures; Fig 30 panel A already plots the 46 % multi-instance share).
6.2 Update every caption number listed in Section 3; `validate_regeneration.py` re-renders 9 generators and only
DECLARES Figs 2, 3, 5, 6, 38, so those five have no byte gate and rely on the md5 hand-check.

### Phase 7. Harness, notebook, documentation (L, 1 day)

7.1 `thesis_expected_values.yaml`: bump `meta.recorded`, add `meta.reference_convention: nearest-copy`, state the
column convention in the gate comment (`:56-58`), then re-transcribe FROM THE PRINTED THESIS (not from the
outputs) every changing block: table_1 (90 of 130 cells), table_2, table_3, table_4, table_6, table_8, table_18,
table_19, table_20 (24 of 36), table_21, table_22 (DiffDock raw row), table_24 near half, table_25, table_26,
table_27, prose.refiner_contrasts, prose.appendix_h9, prose.figure_39, prose.audit, figures.pairs md5s; add a
`reference_convention` block (143 / 138, gains 38 / 29 / 2 and 15 / 12 / 3) and blocks for the two new stages.
About 850 of the 1,364 current result rows move; cohort, determinism, Tables 5, 9, 10, 11, 12, 14, 16, 17 and all
Orai rows stay.
7.2 `thesis_assertions.py`: docstring `:12-30` (202 / 199 → 217 / 214, convention sentence); no gate-line change
under D6 (`:122-124`, `:393`, `:426`, `:611-615`, `:1432-1457`); add `check_reference_convention` and the two new
checks to `run_all` (`:808-841`).
7.3 `_build_reproduction_notebook.py`: `:386-400` add the flag to the stage command and a dated CHANGED note;
`:345` section text (convention sentence, "fourteen downstream stages"); `:484` register the three new stages;
`:1138-1168` section 10 paragraph; run `validate_regeneration.py` once after 6.1 so `:1231-1245` records the
re-based result; rebuild and execute the notebook so `REGENERATE.md` regenerates (it is generated, never hand-edited).
7.4 Documentation: new dated `REPRODUCTION_COVERAGE_<date>.md` and a new glossary version beside
`DEFENCE_NOTES_NUMBER_GLOSSARY.md` (deliverables are never overwritten); edit the living files `README.md:15`
(stage and assertion counts are already stale: 60 stages, 1,364 checks), `THESIS_REPRODUCTION.md:57/:76/:380`,
`thesis_latex/README.md:12` page count; add the convention to `nearest_copy_probe/README.md` or retire the probe.

### Phase 8. Build and QA (M, 1 day)

8.1 `latexmk -f Thesis_short.tex`; rasterise with pymupdf and inspect every regenerated figure page and the
Kurzfassung page for overflow.
8.2 Registry sweep: grep the short build for every OLD value in Sections 1 and 3 (111, 103, 55, 198, 167, 80, 199,
169, 36.6, 34.0, 18.2, 65.3, 55.1, 26.4, 2.6, 9.4, 4.2, 9.8, 112, 311, 1,622, 466, 6,169, 2,329, 62.7, 133, 132,
160, 164, 97.7, 88.1, 84.5, 154, 185, 110 ...) and account for every hit as changed or legitimately retained.
8.3 Run `thesis_assertions.py`; zero mismatches or each mismatch adjudicated in writing.
8.4 Re-read Discussion and Conclusions for claims whose direction could have moved (cost ordering, six-in-ten,
equivalence within ten points, "adds three complexes", "ties with exhaustiveness 64").
8.5 Update project memory: the convention switch, the new stage names, the traps met.

---

## 3. Every affected thesis location

W = wording must change (not only numbers). Values are derived unless marked recompute.

### 3.1 `thesis_latex/Thesis_short.tex`

| Line | Current | New |
|---|---|---|
| 61 | Abstract endpoint definition | W optional: "to the nearest deposited copy of the ligand" |
| 63 | 36.6 / 34.0 / 18.2 %; +2.6; −4.2..+9.4; "six in ten"; 65.3 / 55.1; 80 of 303; "most expensive per qualifying pose" | 49.2 / 43.6 / 18.8; +5.6; −1.7..+12.8 (p 0.159); "six to seven in ten"; 70.3 / 59.1; 83; cost clause pending 1.9 |
| 80 | Kurzfassung endpoint definition | no change (page full) |
| 91 | 36,6 / 34,0 / 18,2 | 49,2 / 43,6 / 18,8 |
| 92 | 2,6 Prozentpunkten | 5,6 |
| 94 | −4,2 bis +9,4 | −1,7 bis +12,8 (+1 char) |
| 95 | etwa sechs von zehn | sechs bis sieben von zehn (+11 chars, cut elsewhere on the page) |
| 97 | 65,3 % und 55,1 % | 70,3 % und 59,1 % |
| 98-99 | AutoDock am teuersten pro qualifizierender Pose | W only if the cost ordering flips (1.9) |
| 100 | 80 von 303 | 83 von 303 |

### 3.2 `thesis_latex/body_main_short.tex`

| Line | Current | New |
|---|---|---|
| 88 | 78 of 303 place a metal within 5 Å of the crystal ligand | 78 (D4); W "of the deposited reference instance" |
| 109 | Methods definition against "the crystallographic ligand" | W: D11 draft |
| 113 | "adds three recovered complexes across the three tools and removes none" | recompute |
| 115, 119 | valid-near-native criterion; "against the crystallographic pocket" | W: nearest deposited copy |
| 124-133 | form definition and mapping caveat | W: form taken against the copy fixed by the in-place minimum |
| 140 | footnote: reference-supplied battery −1 / −2 / −2 DiffDock | recompute (state which file was passed as mol_true) |
| 158 | 4 Å crystal-site reach definition | W |
| 165 | interaction fingerprint against the crystal | W plus D2 method |
| 170 | cost per pose within 2 Å of the crystal ligand | W |
| 191-199 | selection text: 164 vs 160; 138 / 138, 133, 5 / 5; 55 vs 48; 160; 133 vs 132; 834 vs 828; nine, four; 176 vs 173; 237 vs 235; footnote 32.7 / 33.7 / 32.3 | 178 vs 175; 148 / 149 agree on 145, 3 / 4; 57 vs 50; 175; 144 vs 143; 1,031 vs 1,022; seven, three; 186 vs 186; 236 vs 235; footnote recompute (basis: without validity 35.3 → 45.2 %) |
| 205-214 | provenance comments (gates, 4,012 / 3,903 / 3,873; 225 / 237 / 235) | 4,016 / 3,887 / 3,854; 224 / 236 / 235 |
| 264, 277-283 | Table 1 headers | W optional footnote |
| 286-304 | Table 1 rows | Section 1.3 |
| 324 | "at most four" | five |
| 326 | Table 5 lead-in | W |
| 347-349 | Table 5 rows | 49.2 % (149) / 70.3 % (213) / 70.6 % (214) / +64 (+21.1 %) *** / +1 (+0.3 %); 43.6 % (132) / 59.1 % (179) / 60.4 % (183) / +47 (+15.5 %) *** / +4 (+1.3 %); 18.8 % (57) / 27.4 % (83) / 27.4 % (83) / +26 (+8.6 %) *** / +0 |
| 353-372 | footnote and provenance comments | update counts |
| 380 | "as-placed crystal RMSD" | W |
| 382 | Q p < 1e-7; 36.6 / 34.0 / 18.2; +2.6 (−4.2, +9.4) p 0.509; +13.2 | holds; 49.2 / 43.6 / 18.8; +5.6 (−1.7, +12.8) p 0.159; +13.2 unchanged |
| 386 | 112 discordant; 9.8 pp; +2.6 | 129; 10.5; +5.6 (power 0.32) |
| 388 | +28.7 (23.9-34.0); +21.1 (16.9-26.1); +8.3 (5.7-11.9); 36.6 → 65.3; 202 vs 199 | +21.1 (16.9-26.1); +15.5 (11.9-20.0); +8.6 (5.9-12.3); 49.2 → 70.3; 217 vs 214 |
| 393, 401 | Fig 2 / Fig 3 captions | W; regenerate image3 / image4 |
| 396 | 74.3 / 74.6; 95.4 / 97.0; 52.5 / 52.1 / 33.3; 47 / 46, 30 / 31, 25 / 26 | unchanged; unchanged; unchanged; 48 / 47, 30 / 32, 25 / 26 |
| 404 | 33.3 → 54.1; 20.8; "at most 0.7" | 33.3 → 53.8; 20.5; "at most 1.3" |
| 406 | 87 and 64; six in ten | 64 and 47; six to seven in ten |
| 408 | 6,169; 2,329; "crystal-ligand site" | 7,565; script value (≈ 2,854); W |
| 410-414 | 62.7 %; 1.5 → 5.0; +3.55; 48 → 78; 27 → 8; 5.72; +0.46 / +3.55 / +1.09; 137; 2.34; 2.9 Å; 26.4 / 26.4 | 82.2 %; 1.4 → 4.9; +3.55; recompute; recompute; recompute; +0.53 / +3.55 / +0.96; recompute; recompute; 3.3 Å; 27.4 / 27.4 |
| 419 | Fig 4 caption | regenerate image5 |
| 427-435 | sixteen non-reach; Fig 5; 26 / 24 / 4 pp; φ +0.23..+0.28; 2.8e-4 | recompute (D5); regenerate image6 |
| 455-461 | Table 6 reach and co-reach | recompute |
| 469-474 | Fig 6 caption "18"; median six / five poses, 81 / 83 / 47 % reach; 95 %; 56 %; 56.8 / 55.4 / 55.8 | recompute; regenerate image7 |
| 476 | 97.7 % / 0.34; 88.1 / 0.53; 84.5 / 0.56; 97.7 / 0.37; 1.23 | 99.7 %; 92.4; 90.1; medians recompute (trace generator first) |
| 481-532 | interaction methods, Fig 7, 0.450 / 0.417 / 0.380, F1 0.57 / 0.53 / 0.47, Table 7, 21.0 / 19.4 / 15.3 / 16.4 | W plus recompute; regenerate image9 |
| 708 | cost definition | W |
| 715, 751 | Fig 10 / Fig 11 captions: 199, 169, 80; 53 shared | 214, 183, 83; 54; regenerate image24 / image25 |
| 737-741 | Table 8 rows | medians / IQR recompute; complexes 214 / 183 / 83 / 217; poses 442 / 2,068 / 483 / 431 |
| 744-760 | Friedman 85.2 / 3.2e-19 / 0.80; 53; ratios; "yields one"; 311 / 1,622 / 466; 21.7; 1,805.0 vs 1,752.8; 271 / 479 / 38; "most expensive" | recompute; 54; recompute; "two"; 442 / 2,068 / 483; recompute; recompute; recompute; W if ordering flips |
| 768 | 4.2 / 9.4 | 1.7-point DiffDock lead to 12.8-point AutoDock lead |
| 777 | at most four | five |
| 781 | 133 to 132; "exactly level"; 83 vs 61; nine | 144 to 143; "smina ahead" (179 to 177); 86 vs 64; seven |
| 788, 857 | eight pp; six in ten; nine in ten; 87 / 64 / 25 | +8.6; six to seven in ten; nine in ten holds; 64 / 47 / 26 |
| 790 | indistinguishable on form; lead from top-5 | holds |
| 792 | 95 %; 57 % | recompute (D5) |
| 794, 806, 824, 870 | qualitative | hold (cost ordering to confirm) |
| 817 | 80 of 303 vs 169 | 83 vs 183 |
| 834 | 36.6 / 34.0 / 18.2; 65.3 / 55.1 / 26.4; +2.6 (−4.2, +9.4); "eight points" | 49.2 / 43.6 / 18.8; 70.3 / 59.1 / 27.4; +5.6 (−1.7, +12.8); "twelve points" |
| 836 | rank-1 interaction recovery did not separate | re-test |
| 843 | 0.1-0.7; 4.1-4.6; 3.5-4.4; 31.0 → 36.6; Holm clause | recompute bands; 42.9 → 49.2; p 0.018 uncorrected, Holm recompute |
| 848-850 | 2.4 / 49.4 / 137.2 s; 169 / 199 / 80; twentyfold; quarter; RQ3 conclusion | recompute; 183 / 214 / 83; recompute; holds; W if ordering flips |
| 866 | margins eight / one / nine; 3.0; 0.4-0.8 | recompute on the new triple gate |
| 868 | 78; +6.7 on 225; +20.5 on 78; add convention sentence | 78; recompute; recompute; W |

### 3.3 `thesis_latex/body_appendix_short.tex`

| Line | Current | New |
|---|---|---|
| 77 | Vinardo 37 vs 80; 42 vs 90 | recompute (probe on vinardo arms) |
| 165 | 37.8 / 39.6 / 33.3 / 17.9 %; 134 vs 149 of 225; 33 vs 49 of 78; p 0.142 / 0.020 | recompute (stratum stays 78) |
| 173 | 111 (36.6 %) / 112 (37.0 %) / 90 (29.7 %); 31 = 15 + 16; 198 / 197 / 190; 199 | 149 (49.2 %) / recompute / recompute; recompute; 213 / recompute; 214 |
| 182 | 198; Holm 2e-21..0.007; top-3 to top-22; 111 vs 110; 202 vs 199; 160 / 152 / 95 | 213; recompute Holm (uncorrected 1.1e-20..2.6e-3); holds; 149 vs 141 (no longer a tie, W); 217 vs 214; 175 / recompute |
| 184 | 154 → 198; 185 → 198 (+17 / −4); 93 / 183 vs 94 / 185 | 170 → 213; 198 → 213 (+19 / −4); recompute |
| 189, 195, 235, 241 | Figs 18-21 captions ("three complexes at rank-1", "depth 22", "buys nothing at rank-1") | recompute; regenerate image44-47 |
| 218-227 | ladder table and footnote | Section 1.5; W footnote |
| 230 | +71.6 / +35.9 / +5.9 per hour; +65.6 (+32.1..+137.2) | keep confirmatory on instance (D10), print nearest beside |
| 244-248 | five; three; 18.7 pp; ~2 pp; ~1.5 pp | recompute |
| 246 | 95; 113; "exactly 113" | 131 raw exh128, 151 rescored (near-native without validity) |
| 267-269 | 2 / 46 / 55; 0 / 7 / 14; 0 / 5 / 18; 37; 19 / 64 / 83; 180-200; 261 / 282 / 266 | 2 / 48 / 57; recompute; recompute; recompute; 19 / 67 / 86; 179-199; unchanged |
| 276 | "to the crystal ligand" | W |
| 487 | fifth reference-based check | optional clause on load_all |
| 543-551 | 4,253; 7.57 Å; 21.3 %; F1 0.565 / 0.529 / 0.471; 750; 0.0003; fingerprint keys with chain | recompute after per-copy fingerprints; W |
| 662 | footnote roles of the four files | W: documented roles (`_ligands.sdf` all instances, `_ligand.sdf` site marker) |
| 664 | convention sentence | W: draft in Phase 5 |
| 782 | 1.56 Å; 65 % | unchanged (D13), W one clause |
| 869 | 34.0 % and 34.3 % | 43.6 % and recompute (DiffDock + gnina rank-1 134 valid-near = 44.2 %) |
| 891-918 | top-k Wilson table | derived for 5 of 7 variants (e.g. AutoDock raw k = 1 43.2 [37.8, 48.9]); two recompute |
| 945-962 | optimisation table and footnote | AutoDock k = 1 43.2 → 49.8, 19 / 39, p 0.012; others derived / recompute |
| 990-997 | band table and footnote | DiffDock smina near-native +0.6 / +0.4 / +0.8 derived; rest recompute |
| 1024-1029 | cross-tool table and footnote | row 1: 49.8 vs 43.6, 67.3 vs 50.2, 70.0 vs 52.1, 71.3 vs 52.8; k = 1 77 / 58 uncorrected 0.12; rest recompute |
| 1054-1060 | accurate-but-invalid table | AutoDock + gnina unchanged (2 / 4 / 3 / 3); DiffDock raw 21.5 [17.2, 26.4] at k = 1; footnote unchanged |
| 1072 | ITT: none within 2 Å; closest 6.4 Å; 111 / 308, 198 / 308 vs 103, 167 | Section 1.10: rank-1 none (6.3 Å); top-15 none but closest 2.1 Å; 149 / 308, 213 / 308 vs 132, 179 |
| 1170-1179 | PB decomposition near-native half; footnote | 437 / 98.6; 451 / 98.0; 2,118; 2,169 / 95.3; 2,171; 515 / 93.8 ...; per-check split recompute; W footnote |
| 1186-1194 | 5.4 %; 1,715; 6.2 %; 498; 31 of 498; 5 of 319; "at most four"; 2 / 4 / 3 | strata 2,169 / 515 / 451; rates recompute; "at most five"; AutoDock 2 / 4 / 3 / 3 / 3 unchanged, DiffDock 5 / 4 / 3 / 3 / 3 |
| 1201 | 25 Å cube volume argument | W: add the regime clause (cube excluded the copies) |
| 1205 | 148 → 185 | recompute → 198 |
| 1245-1248 | +88 / −2 / +1 / +87; +63 / −1 / +2 / +64; +26 / −2 / +1 / +25 | +65 / −2 / +1 / +64; +45 / 0 / +2 / +47; +27 / −2 / +1 / +26 |
| 1256-1258 | +2.6 [−4.2, +9.4]; +13.2; +10.2 [+2.8, +17.5]; +9.9; 51 / 140 / 112 (37.0 %); MDD 9.8; power 0.12; 4,200; 8.2; "within 10 points" | +5.6 [−1.7, +12.8]; +13.2 unchanged; +11.2 [+3.9, +18.4]; +10.2 [+3.0, +17.3]; 76 / 98 / 129 (42.6 %); 10.5; 0.32; ~1,060; 11.7; "within 12 points" (W) |
| 1263 | Kabsch criterion sentence | W: best-fit to the same copy |
| 1301-1318 | Table 26 RMSD rows | derived (AutoDock + gnina d = 1: 0.0 8.6 18.2 30.0 38.0 41.6 44.9 49.2 49.8 53.5 ...) |
| 1319-1338 | Table 26 Kabsch rows | 26 of 90 cells move 0.3-0.7 pp |
| 1375-1385 | Table 27 and footnote | all 117 cells (e.g. AutoDock (1) 249 / 249 / 82.2 % / 1.376 / 0.682 / 3.734 / 0.824 / 0.402 / 1.513 / 0.368 / 46.6 / 25.7 / 27.7); W footnote |
| 1392-1397 | cluster statistics; Fig 38 caption | recompute (D5); regenerate image8; W |
| 1405-1413 | Kendall τ; matched / missed / spurious 15.3 vs 14.1; Fig 39 | recompute; regenerate image11; W |

---

## 4. Affected figures, tables, sidecars and harness assertions

**Figures to regenerate (14 of 39):** Fig 2 (image3), 3 (image4), 4 (image5), 5 (image6), 6 (image7), 7 (image9),
10 (image24), 11 (image25), 18 (image44), 19 (image45), 20 (image46), 21 (image47), 38 (image8), 39 (image11).
Producing stages: hub (`18_topn_within_thresholds_*`), `bench_filmstrip`, `bench_clusters`, `bench_pandamap_report`,
`bench_effort_charged` / `bench_effort_elapsed`, `bench_exhaustiveness_returns`. Two shipped assets (image12, image14)
are PyMOL renders with no source match by design.

**Tables:** 1, 5, 6, 7, 8 (main); ladder (App. C), top-k Wilson, optimisation effect, band gains, cross-tool,
accurate-but-invalid, PB decomposition near-native half, depth-gain, threshold-resolved RMSD and Kabsch, placement /
form (App. H, I); new sensitivity table (D9). Unchanged: 2, 3, 4 (dataset), protocol tables, Orai tables, cost totals
(Table 10), PoseBusters checks, design-to-test mapping.

**Hub sidecars that change (all regenerated by one `--force` run):** `per_pose_metrics.csv` and manifest,
`oracle_summary*.csv`, `top1_summary*.csv`, `per_rank_metrics.csv`, `per_rank_ifp_recovery.csv`,
`topn_within_thresholds*.csv` and the `18_*_report.txt` files, `topk_recovery_validity*.csv`, `rank1_vs_topn.csv`
(Gen A, never diff against Tables 5-7), `optimization_raw_vs_best.csv`, `optimization_benefit_by_rank.csv`,
`pose_validity_cascade*.csv`, `within2_validity_comparison.csv`, `form_fidelity_summary*.csv`,
`oracle_selection_comparison.csv`, `filmstrip_stats__*.csv`, `20d_*` tables, `twist_turn_summary.csv`,
`rank_quality_*.csv`, `pb_test_waterfall*.csv`, `pbvalid_filter_influence*`, `posebusters_pose_comparison_stats.json`,
the `pb_valid/` re-derivations. Unchanged: `pbvalid_yield_*` (validity only).

**Harness:** ~850 of 1,364 assertion rows; yaml blocks listed in 7.1; 14 md5 pairs; three new stages
(`bench_reference_convention`, metal stratum, ITT) taking the registry from 60 to 63.

---

## 5. Risks and traps

1. **Selection flips.** The DiffDock refiner margin is one complex under both conventions (144 vs 143). Phase 4 gates
   Phase 5. If smina falls behind, every "DiffDock*" number in the thesis moves to the gnina arm and `:197`, `:781`,
   `:962`, `:869` change direction.
2. **"No complex loses" is false for Form.** Three arms lose one form complex under D2. Do not copy the M2 sentence
   into the Form or Combined columns.
3. **Gen A / Gen B.** `rank1_vs_topn.csv` (min-RMSD pose then validity) will again differ from Tables 5-7 (existence
   gate); the divergence pattern changes because the argmin pose itself moves. Never diff them.
4. **`rmsd` versus `pb_rmsd`.** `pb_rmsd` symmetrises terminal groups and is ≤ `rmsd`; the ≤ 2 Å verdicts agreed on
   414/414 rank-1 poses tested but the disclosure at `:113` must be re-derived at all depths.
5. **Cache convention blindness.** `--reuse-cache --collapse-plots-only` today reads any cache; without D8 a
   post-switch re-render on a pre-switch cache prints single-instance numbers with no error.
6. **Media by md5, not timestamp.** Figures 2, 3, 5, 6, 38 have no byte gate in `validate_regeneration.py`.
7. **Kurzfassung.** +1 character (interval) and +11 (six to seven) on a page with zero slack; pay with a cut.
8. **PandaMap receptor chains.** Whether `hetatm_cleaned/<ID>_protein.pdb` keeps every chain is unverified; if not,
   the interaction family needs the disclosed-cohort fallback (2.2).
9. **Adjacent copies inside the 8 Å clustering radius.** Six ids; D5's count-once rule must be implemented and stated.
10. **Unregistered generators.** The 78-complex metal stratum and the ITT 6.4 Å figure have no script today; they must
    be registered before the numbers can be re-derived reproducibly.
11. **Pre-existing drift.** `:476` (97.7 vs 98.0 %, DiffDock 84.5 vs 85.1 %) and `:408` (2,329 vs 2,337) are already
    off; fix the provenance before overwriting, or the new numbers inherit an unexplained offset.
12. **Pre-specified confirmatory cell** (App. C :230): keep on instance (D10) or the paragraph contradicts itself.
13. **EquiBind ranking ties.** 333 poses in 133 same-complex groups tie on `gnina_affinity`; Table 27 depth-5 cells
    depend on the tie-break and must come from the regenerated sidecar, not from a hand recomputation.
14. **Never overwrite deliverables.** Glossary and coverage documents get new dated versions; `REGENERATE.md` is
    generated; the full build (`obsolete/`) is not touched.
15. **Yaml transcription direction.** Values are transcribed from the printed thesis, then checked against outputs.
    Copying outputs into the yaml makes the harness a tautology (how Table 21 drifted).
16. **Stale counts already in the repo.** README says 257 numbers and 60 stages; the harness runs 1,364 checks;
    THESIS_REPRODUCTION.md says 58 stages. Correct them in the same pass.
17. **Comparability wording.** Once the scoring axis matches the paper, only the search-space axis separates the
    studies; `:868` and `App. C :1201` must say exactly that and no more.
18. **Unaudited inventories.** Seven of eight dimensions were not adversarially audited (Phase 0.2).

---

## 6. Effort estimate and critical path

| Phase | Effort | Wall |
|---|---|---|
| 0 Decide, audit, freeze | 0.5 day | 0.5 day |
| 1 Hub code + regression | 1 day | 1 day |
| 2 Downstream code + new stages | 1.5 days | 1.5 days |
| 3 Rebuild | 1-2 h compute, 0.5 day attention | 0.5 day |
| 4 Selection re-check | 0.5 day | 0.5 day |
| 5 Thesis text | 2-3 days | 2.5 days |
| 6 Figures and media | 0.5 day | 0.5 day |
| 7 Harness and docs | 1 day | 1 day |
| 8 Build and QA | 1 day | 1 day |
| **Total** | **8.5-9.5 working days** | ≈ 2 weeks with review |

Critical path: hub code → rebuild → PandaMap per-copy fingerprints → selection re-check → text → harness → QA.
The PandaMap per-copy work (2.2-2.3) is the single largest and least certain item; if its receptor check fails,
the fallback costs a day less but adds a disclosed exemption to `:165`, `:481` and Table 7.

---

## 7. What does NOT change

- **Orai1 panels** (Figs 8-14, 22-27, Tables 9-15, Sections 4.2 and 5.2): no crystal ligand exists, nothing is
  crystal-referenced.
- **X-ray crystal control** (`XRay_PoseBusters_Control.ipynb`, `_build_xray_control_notebook.py:283/:503`): busts the
  deposited ligand with no reference and no RMSD.
- **The PoseBusters validity battery and Fig 1**: `pb_valid` is reference-free (`run_posebusters.py` passes `None` as
  `mol_true`); Table 1 poses and valid columns, the 09f yield outputs and Table 24's all-poses half are unchanged.
- **Cost totals and hardware currencies** (Table 10): only the per-qualifying-pose normalisation moves.
- **Dataset descriptors** (Figs 28-37, Tables 2-3): Fig 30 panel A already shows the 46 % multi-instance share;
  the start-conformer RMSD stays on the single instance (D13).
- **Top-5 AutoDock-DiffDock contrast**: identical under both conventions (88 / 48, +13.2 pp, p_holm 7.6e-4).
- **Kabsch headline**: 74.3 / 74.6 % at 1 Å and 95.4 / 97.0 % at 2 Å; form tie p_holm 1.000.
- **Validity-gate cost for AutoDock**: 2 / 4 / 3 / 3 complexes at d = 1 / 5 / 15 / 30.
- **Table 22 rows** for AutoDock + gnina, EquiBind + gnina and the raw-AutoDock footnote coincide numerically
  (confirm by harness, not by inspection).
- **Metal-adjacent stratum size** 78 (D4), the five-complex exclusion, the 303 denominator, all determinism evidence.
- **Docking outputs, PoseBusters runs, receptor and ligand preparation**: nothing is re-docked or re-busted.
