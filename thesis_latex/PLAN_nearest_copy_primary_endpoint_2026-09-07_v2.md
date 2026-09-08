# Plan v2: switch the primary near-native endpoint to the nearest crystallographic copy (2026-09-07, audited)

**Scope.** Replace "symmetry-corrected in-place heavy-atom RMSD to the single deposited instance `<ID>_ligand.sdf`"
by "the same RMSD to the nearest of all deposited copies in `<ID>_ligands.sdf`", per pose, with every crystal-referenced
metric of that pose derived from the copy it is nearest to. The single-instance value is kept as a secondary column and
becomes the appendix sensitivity arm.
**Status.** PLAN ONLY. Nothing in the thesis, the scripts or the harness has been changed. v1 is kept unchanged beside
this file; Section 8 lists what v2 changed and why.
**Verification status.** All 17 outstanding adversarial audits ran on 2026-09-07 (22:26-23:10): an evidence auditor and a
completeness auditor for eight dimensions, two auditors over the hand-computed numerical section, then a plan critic and
six gap-fill agents. Totals: 7 refutations, 80 line or value corrections, 91 additional sites, 117 plan errors
(1 critical, 36 major, 80 minor). The critic accepted 28 of 29 distinct critical-or-major findings and rejected one.
Every headline number in v1 Sections 1.1, 1.2 and D12 survived every auditor's recomputation; the defects were
sequencing (pins, PandaMap gate, cache schemas), undecided or self-contradicting decisions, mixed gate bases, stale line
anchors and table numbers. The audit JSONs are archived under
`Scripts/Analysis/nearest_copy_probe/wf_plan_2026-09-07/dossier/` (`audit__*.json`, `critic.json`, `gapfill.json`).
**Anchors.** `body_appendix_short.tex` and `body_main_short.tex` were edited at 22:03 on 2026-09-07 (two lines inserted at
appendix :558-560), so every appendix anchor from :560 onward is two lines later than in v1. Line numbers below are as of
22:30 on 2026-09-07 and every edit is to be located by the quoted text, never by the line.
**Table and figure numbering.** Thesis numbering from `Thesis_short.lot` / `.lof`: Table 1 pose production, Table 2 depth
recovery, Table 3 crystal-cluster reach, Table 4 native-interaction recovery, Table 5 Orai1 yield, Table 6 cost per
qualifying pose, Table 8 exhaustiveness ladder, Tables 16-17 ligand descriptors, Tables 18-27 appendix statistics.
Figures 15 and 16 are the cost figures (image24, image25). v1 used a different numbering in places; v2 uses only this one,
and it matches the yaml block keys `table_N`.

---

## 0. Decision record

| # | Decision | Adopted | Notes from the audit |
|---|---|---|---|
| D1 | Definition of "nearest" | Per pose, argmin over all records of `_ligands.sdf` of the thesis's own in-place `symmetry_rmsd`, NaN as +inf, record 0 (the reference) as tie-breaker | Matches PoseBusters' `check_rmsd` argmin on 414/414 rank-1 verdicts tested. Centroid-nearest and RMSD-nearest disagree on 317 of 24,486 multi-copy rows across 77 ids (not "six ids" as v1 said), which is why one explicit index must drive every metric. |
| D2 | Which copy drives centroid, Kabsch, torsions, contacts, PLIF | The single argmin-RMSD copy j* for every metric of that pose (PoseBusters `rmsd.py:61-70`, `choose_by="rmsd"`) | Cost measured on the flipped rank-1 poses: median Kabsch shift 0.013 Å (v1 said 0.001), max 0.25 / 0.68 / 0.03 Å; form verdict changes on 0 / 2 / 0 of the 38 / 30 / 2 flipped poses (30 not 29: the 29 is the PB-valid subset). 4 Å reach at j* differs from the per-metric minimum by ≤ 2 complexes (Section 1.7). |
| D3 | Alternates that do not share the reference pocket (7A9E_R4W, 7TUO_KL9, 7VKZ_NOJ, 7Z1Q_NIO) | Include every deposited instance, as PoseBusters and the dataset README do; name the four in the appendix | No complex GAINS recovery through them at any depth in any arm. The v1 draft sentence "no reported pose recovers any of them" is false: three PB-valid ladder poses sit 1.70-1.75 Å from 7Z1Q_NIO's alternate (exh32 rank 30, exh18 rank 28, exh32+gnina rank 21), all in complexes already recovered on the reference. Write "no complex gains recovery through any of them". |
| D4 | Metal-adjacent stratum (78 of 303) | Define by the rule that reproduces the printed 73 / 78 / 81: minimum heavy-atom distance from the reference instance to ANY atom of a residue containing a metal element, on `Data/PoseBuster Benchmark Set/<ID>/<ID>_protein.pdb`, ≤ 5 Å (4 Å → 73, 6 Å → 81). Membership is a complex-level property of the reference instance and does not follow the scored copy | The atom-level rule gives 65 / 77 / 80 and contradicts the print. Membership is NOT copy-invariant: 7JHQ_VAJ is in by the reference only (Na at 2.42 Å; alternates 21.8-39.5 Å), 7TB0_UD1 is in by copy 1 only (K at 3.17 Å; reference 9.04 Å). Under the switch one AutoDock* rank-1 success (7JHQ_VAJ, 1.49 Å to a metal-free copy) counts as metal-adjacent and one DiffDock* rank-1 / top-15 success (7TB0_UD1, 1.51 / 1.46 Å to the metal-adjacent copy) counts as metal-free. Disclose both. The registered rule yields 16 cofactor / 62 free-ion, where App. C :165 prints 15 / 63; assert whichever the rule yields. |
| D5 | Crystal-cluster reach with several crystal sites | A cluster reaches if its centroid lies within 4 Å of ANY copy centroid; the crystal-closest cluster is defined by the nearest copy; the six ids whose adjacent copy sits 4.9-7.8 Å from the reference (inside the 8 Å clustering radius) count once | Consequence for :476 (Section 1.7): the cluster report's pooled 4 Å oracle follows D5 (any copy, 302/303 = 99.7 %) while the hub's `oracle_centroid_le_4A_%` follows D2 (301/303 = 99.3 %). State that :476 prints the cluster-report value under D5. |
| D6 | Column naming | Keep `rmsd`, `pb_rmsd`, `pb_kabsch_rmsd`, `pb_rmsd_within_2A`, `bestfit_rmsd`, `centroid_dist`, `contact_recovery`, `plif_recovery` as primary names carrying the nearest-copy value; add `rmsd_ref_instance`, `centroid_dist_ref_instance`, `pb_rmsd_ref_instance`, `bestfit_rmsd_ref_instance`, `n_copies`, `ref_copy_index`, `nearest_copy_index`, `nearest_copy_is_ref`, `reference_convention` | 14 consumers, two test files and seven gate lines in `thesis_assertions.py` stay untouched. `pb_rmsd_within_2A` is also crystal-referenced (v1 omitted it). |
| D7 | Where the new cache lives | Back up the current report directory as `pose_comparison_report.bak-instance-<date>` and write the rebuild to the canonical path | Pattern precedent only: `.bak-preexh128rerun-20260830` exists under `benchmark_full_protein_vina_scoring/`, not under `benchmark_matched_equibind/dock/`. |
| D8 | CLI and cache safety | Add `--reference-convention {instance,nearest}`; keep the default `instance` until Phase 7; pass the flag explicitly at every invocation in Phases 1-6; bump the cache schema 5 → 6; add the new columns to `_CACHE_REQUIRED_COLS`; record the convention in the manifest; `--reuse-cache` refuses a cache built under the other convention; the non-reuse path REFUSES (not recomputes) when only the convention differs unless `--force` is given | v1 flipped the default in Phase 1 and edited the stage command in Phase 7, leaving a five-phase window in which any bare call rebuilt silently. Today a signature mismatch at hub `:13274-13276` returns None and triggers a full recompute; the only refusal is the `--reuse-cache` branch at `:13478-13482`, which loads without a signature check (`:13261-13266`). The schema bump alone invalidates every live manifest, so the window opens at the schema bump, not at the default flip. |
| D9 | Sensitivity arm | The single-instance convention becomes the appendix sensitivity table (the M2 table with columns inverted) with a Form and a Combined column, plus a single-copy stratum row on the SAME PB-valid ∧ ≤ 2 Å gate: 88 / 80 / 53 of 165 = 53.3 / 48.5 / 32.1 % rank-1 | v1's 54.5 / 50.9 / 32.7 % were RMSD-only (90 / 84 / 54 of 165) and would have sat beside validity-gated columns. The Form column is required because the Form block LOSES complexes under D2 (Section 1.3). |
| D10 | Pre-specified confirmatory cell (App. C :230, +65.6 complexes per hour, 95 % +32.1 to +137.2) | Keep the confirmatory reading on the single-instance column and print the nearest-copy value beside it | Scope stated explicitly: the same paragraph prints descriptive pool counts "98, 128, 149, 159 and 164" and three per-hour returns. Either the whole paragraph stays on instance with one sentence giving the nearest values, or the descriptive numbers move and only the confirmatory cell stays. Adopt the first, because Figures 20-21 and :182 ("160 against 152") are on the same pool and would otherwise print 164 (instance) at :230 beside 178 (nearest) at :182. |
| D11 | Methods wording | Redefine the endpoint at `body_main_short.tex:109-115` and cite the external convention (README :47-49 and :55-59; `redock.yml:307 load_all: True`; paper Section 2.3 "closest crystallographic ligand") | Draft in Phase 5. The inertness of the PUBLISHED Vina anchor is inferred from the thesis's own boxed arm (0 flips), not measured on the paper's poses, which were never released; and 8 of 211 alternate centroids (16 ids by any heavy atom) lie inside a 25 Å cube on the reference, so write "kept almost every other copy outside". |
| D12 | Variant selection | Re-run the top-15 triple gate under the new convention. Measured: AutoDock* 160 → 175; DiffDock smina 133 → 144 against gnina 132 → 143 (1,031 vs 1,022 poses, seven discordant, three gnina-only); EquiBind gnina 55 → 57 against smina 46 → 48. All three selections hold | Confirmed by three auditors. |
| D13 | Start-conformer difficulty descriptor (App. F :784, Fig 37) | Keep on the single instance; say it is a start-conformer property | Unchanged. |
| D14 (new) | Holm family for the AutoDock raw → gnina rank-1 step | Name the family the body quotes. In App. C's seven-contrast family the rank-1 raw-exh128 vs selected contrast has Holm p 0.037 (survives) under nearest; in Table 19's sixteen-test within-tool family the same contrast has Holm p 0.060 (does not). Today `:843` says "does not clear Holm on its own" and App. C :182 says two of seven survive; under nearest six of seven survive, so the two sentences contradict each other unless the family is named | Adopt: `:843` quotes Table 19's family and says so; App. C :182 quotes its own seven-contrast family. |
| D15 (new) | Rescore scope | Rescore the canonical `benchmark_matched_equibind` tree AND the `benchmark_full_protein_vina_scoring` tree that `Master_Docking_AD_Full_Protein.ipynb` cells 33 / 35 / 37 / 46 / 121-124 read and that four scripts default to (`thesis_endpoint_diagnostics.py:41`, `optimization_benefit_stats.py:348`, `autodock_exhaustiveness_returns.py:222`, `rebust_hydrogen_normalised.py:67`); leave the boxed `benchmark/dock` tree with a note (0 flips) | v1 never mentioned the second tree. A mixed-convention notebook is the alternative. |
| D16 (new) | Single source of truth for printed numbers | Probe CSVs are previews used to write expectations in advance; every printed number comes from the Phase 3 rebuild. The Vinardo and guided-gnina arms need no probe because the rebuild scores all 27 method keys | The probe already disagrees with a printed EquiBind-raw ordering at k = 5 / 10 (3.3 / 5.3 vs 3.0 / 4.0), which shows why previews must not be transcribed. |
| D17 (new) | `pocket_comparison_report.py` (fpocket / P2Rank pocket accuracy to the crystal, `:1141`) | Out of scope; none of its crystal outputs appears in the short build or the yaml. List in Section 7 | Verified by grep. |

**Data facts (308 ids, all verified twice).** 143 multi-copy ids (111×2, 7×3, 22×4, one each ×5, ×6, ×12), 138 in the 303,
211 alternate copies. The reference is record 0 of `_ligands.sdf` in 308/308 and its coordinates match `_ligand.sdf` to
1e-4 Å in all 428 folders. All 211 alternates have the same heavy-atom count, formula, element sequence, formal charge and
all-atom count as the reference, so no copy can silently drop out of the minimum; the files carry no hydrogens (519/519
records on the 308). Alternate centroid separation from the reference: none within 4 Å, six within 8 Å (7A9E_R4W 4.9,
7CD9_FVR 5.0, 7JY3_VUD 5.8, 7VKZ_NOJ 6.0, 7Z1Q_NIO 6.2, 6YYO_Q1K 7.8), eight within 12.5 Å, 203 beyond 12 Å; median 36.0 Å
over the 143 ids (36.3 Å over the 138 analysed ids, 40.4 Å over all 211 alternates; pair the denominator with the
statistic). Six ids carry the second copy on the same chain id (7A9E_R4W, 7P2I_MFU, 7TUO_KL9, 7VKZ_NOJ, 7WL4_JFU,
7Z1Q_NIO). All five dropped complexes are multi-copy (2, 2, 4, 4, 4 records).

**External convention.** Chem. Sci. 2024 (PMC10901501) Section 2.3: "the minimum heavy-atom symmetry-aware RMSD between
the predicted ligand binding mode and the closest crystallographic ligand". Installed `redock.yml:301-307` sets
`load_all: True` on both `mol_pred` and `mol_true`; `modules/rmsd.py:52-70` takes argmin over conformers and reports
kabsch and centroid from that copy; paper-era v0.2.2 did the same. Two loader facts matter for implementation:
`tools/loading.py:158-166` combines records "without checking identity or atom order", and because `next(supplier)`
restarts the iteration it adds record 0 twice (n + 1 conformers). Build the combined `mol_true` from the thesis's own
per-record loader instead. Dataset `README.txt:47-49` and `:55-59`: `_protein.pdb` is "without the ligand of interest"
(every copy removed), `_ligands.sdf` is "All instances of the ligand of interest", `_ligand.sdf` "marks the binding
site for those docking methods that require a binding site". The thesis footnote at `body_appendix_short.tex:664`
calling `_ligand.sdf` the ground-truth answer pose misdescribes it. Zenodo results on the 308, `post-processing == none`
rows only (the CSV also carries `energy minimization` rows: Vina 54.9 %, DiffDock 39.6 %): Vina 59.7 % RMSD ≤ 2 Å,
DiffDock 38.0 %. Under that nearest-copy scoring the paper's DiffDock rate is 31.5 % on multi-copy against 43.6 % on
single-copy ids, and its Vina rate 57.3 % against 61.8 %, so part of the multi / single gap is a difficulty confound
rather than the convention.

---

## 1. What changes numerically

Convention "nearest" means D1 + D2. n = 303. Every value marked derived was recomputed by at least two auditors from the
probe CSVs after reproducing the printed reference value to the digit. Every table row below now carries its GATE.
Per D16 these are previews for writing yaml expectations in advance; the printed numbers come from the Phase 3 rebuild.

### 1.1 Headline recovery (gate: PB-valid ∧ ≤ 2 Å, existence over top-d; Table 2, Abstract, Kurzfassung, Conclusions)

| Arm | d | Reference | Nearest | Gain |
|---|---|---|---|---|
| AutoDock* | 1 | 111 (36.6 %) | 149 (49.2 %) | +38 |
| DiffDock* | 1 | 103 (34.0 %) | 132 (43.6 %) | +29 |
| EquiBind* | 1 | 55 (18.2 %) | 57 (18.8 %) | +2 |
| AutoDock* | 15 | 198 (65.3 %) | 213 (70.3 %) | +15 |
| DiffDock* | 15 | 167 (55.1 %) | 179 (59.1 %) | +12 |
| EquiBind* | 15 | 80 (26.4 %) | 83 (27.4 %) | +3 |
| AutoDock* | 30 | 199 (65.7 %) | 214 (70.6 %) | +15 |
| DiffDock* | 30 | 169 (55.8 %) | 183 (60.4 %) | +14 |
| EquiBind* | 30 | 80 (26.4 %) | 83 (27.4 %) | +3 |

Rank-1 → top-15 gains: +87 → +64 (+28.7 → +21.1 pp, Wilson 16.9-26.1), +64 → +47 (+15.5, 11.9-20.0), +25 → +26 (+8.6,
5.9-12.3). Top-15 → top-30: +1 / +4 / +0 (DiffDock b = 4, c = 0, p = 0.125). Near-native without the validity gate at
rank-1: 113 / 107 / 57 → 151 / 137 / 59. Validity-gate cost by depth 1 / 5 / 15 / 30: AutoDock* 2 / 4 / 3 / 3
(unchanged), DiffDock* 4 / 4 / 3 / 4 → 5 / 4 / 3 / 3, EquiBind* 2 / 3 / 3 / 3, so "at most four" becomes "at most five"
and the span 0.66-1.32 pp becomes 0.66-1.65 pp (`:324`, `:777`, App. H :1196). Rank-1 failures with no qualifying pose
anywhere: 54.2 / 67.0 / 89.9 % → 57.8 / 70.2 / 89.4 %. Raw exh128: 31.0 / 61.1 / 66.7 % → 42.9 / 65.3 / 71.6 %, and it
still leads the rescored arm over the pool, 217 against 214. Rescoring adds 442 − 431 = 11 qualifying poses at no CPU cost
(`:758` says nine).

### 1.2 AutoDock* against DiffDock* (gate: PB-valid ∧ ≤ 2 Å; `:382-388`, `:834`, App. H.9 :1258-1260, glossary)

| d | Convention | Margin | Newcombe 95 % | Discordant | Exact McNemar p |
|---|---|---|---|---|---|
| 1 | reference | +8 (+2.6 pp) | −4.2 to +9.4 | 60 / 52 (112) | 0.509 |
| 1 | nearest | +17 (+5.6 pp) | −1.7 to +12.8 | 73 / 56 (129) | 0.159 |
| 5 | both | +13.2 pp | +5.7 to +20.4 | 88 / 48 | 7.6e-4 (identical) |
| 15 | reference | +10.2 pp | +2.8 to +17.5 | 82 / 51 | 0.009 |
| 15 | nearest | +11.2 pp | +3.9 to +18.4 | 82 / 48 | 0.004 |
| 30 | reference | +9.9 pp | +2.5 to +17.1 | 81 / 51 | 0.011 |
| 30 | nearest | +10.2 pp | +3.0 to +17.3 | 79 / 48 | 0.008 |

Power (`:386`, App. H.9): discordant 112 → 129 (37.0 → 42.6 % of the set), MDD 9.8 → 10.5 pp at 80 % power, observed
power 0.12 → 0.32, complexes needed 4,162 → 1,062 (the harness compares to the nearest hundred, so the yaml reads 1100 and
the prose "roughly 1,100"), the 428-entry MDD 8.2 → 8.8 (v1 wrongly wrote 11.7 there). Equivalence: 90 % interval
−3.1..+8.3 → −0.5..+11.7, so rank-1 is no longer equivalent within 10 points, only within 12; smallest passing margins
[8.4, 19.4, 16.4, 16.0] → [11.7, 19.3, 17.2, 16.2]. On the intention-to-treat 308 basis the margins are 5.5 pp at rank-1
and 11.0 at top-15 (App. I :1074 prints 2.6 and "10.2 to 10.1").

### 1.3 Table 1 pose-level columns (`body_main_short.tex:286-304`; yaml table_1 :71-80; D2 "follow")

Validity columns unchanged. Table 1 has ten rows; Uni-Dock2 is deliberately omitted (`:240`), so its values below are
for the sidecars only.

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
| EquiBind fpocket raw | 0 → 0 | 180 / 2,449 → 179 / 2,441 | 0 → 0 |
| EquiBind P2Rank raw | 0 → 1 (2 poses) | 181 / 2,336 → 180 / 2,339 | 0 → 0 |
| (Uni-Dock2, sidecars only) | 81 / 117 / 1.4 → 95 / 162 / 1.9 | 196 / 2,729 → 196 / 2,715 | 65 / 85 / 1.0 → 82 / 125 / 1.5 |

The Form block LOSES a complex in FIVE printed rows (DiffDock raw, DiffDock + smina, EquiBind + gnina, fpocket raw, P2Rank
raw), and 7 of the 11 scored arms lose at least one form complex (two regain one elsewhere). Pose level: 149 lost / 91
gained over the six headline-family arms, 234 / 167 over all eleven. The M2 sentence "no complex loses" holds for
near-nativeness only and must appear nowhere near the Form or Combined columns; the D9 table carries a Form column and
the :666 draft discloses the loss. Under the yaml's exact tolerances 68 of the 130 `table_1` rows change under D2
(70 under an independent Kabsch minimum), not "90 of 130" as v1 said; `poses` and `valid_*` rows never move.

### 1.4 Form axis (Kabsch; gate: PB-valid ∧ Kabsch ≤ 1 Å; `:396`, `:404`, Fig 3, Table 26 Kabsch rows, glossary :151-159)

Unchanged at the printed precision: 74.3 / 74.6 % at 1 Å best-of-top-30, 95.4 / 97.0 % at 2 Å, rank-1
52.5 / 52.1 / 33.3 %. EquiBind best-of-top-30 54.1 → 53.8 %. 47 of the 90 Table 26 Kabsch cells move by 0.3-0.7 pp (v1
said 26). Form McNemar AutoDock vs DiffDock: 47 / 46, 30 / 31, 25 / 26 → 48 / 47, 30 / 32, 25 / 26; the tie verdict holds
but the top-15 Holm p becomes 0.899 within its three-pair family, so the literal "p_holm = 1.000 at every depth" at
`:396` changes to "1.000, 0.899 and 1.000" or "ns at every depth". The same sentence's "exceed EquiBind at every depth
(p_holm ≤ 4e-10)" is already wrong today (companion report: 8.997e-09 at rank-1; Section 0.4).

### 1.5 Variant selection and the exhaustiveness ladder (App. C :178-248, Table 8, `:191-199`, `:781`, `:866`)

Top-15 triple gate (PB-valid ∧ ≤ 2 Å ∧ Kabsch ≤ 1 Å): AutoDock* 160 → 175; DiffDock smina 133 → 144 against gnina
132 → 143 (poses 834 / 828 → 1,031 / 1,022; discordant nine → seven, gnina-only four → three); EquiBind gnina 55 → 57
against smina 46 → 48. Double gate top-15: smina 167 → 179 against gnina 167 → 177. Near-native only (no validity, no
form) top-15: smina 182 against gnina 180. So `:781` "dropping either term leaves the two arms exactly level, and gnina
is level or ahead on every component" becomes "smina ahead by two on both drops; the arms are level on pooled
near-native recovery (186 against 186) and smina stays ahead on pooled form (236 against 235)". Note `:781` is at the
top-15 selection depth and `:197` is Table 1's pooled basis; both are true today and the sentence at `:781` should say
"at the top-15 selection depth". Winning margins at `:866`: nine (AutoDock, 175 vs raw 166) / one / nine, largest 3.0 pp,
unchanged in wording.

Ladder (gate: PB-valid ∧ ≤ 2 Å; d = 1 / 15 / 30; reference → nearest; d = 5 and the Holm column from the rebuild):

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

The selected rung beats each of the seven alternatives at top-15 (uncorrected 1.1e-20 to 2.6e-3, Holm 7.3e-20 to 5.2e-3
against its own raw search and 64 + gnina; App. C :182 prints "2e-21 to 0.007"). Depth lead: reference leads at 1 and
3-22 (ties at 2 and 23); nearest leads strictly at every depth 1-25 and ties at 26, so "from top-3 to top-22" is
strengthened, not merely retained. Rank-1: the "tie with 64 rescored, 111 against 110" becomes 149 against 141, but that
is 16 / 8 discordant, exact p 0.15, a nominal lead that stays unresolved; and "only two of the seven contrasts survive
correction" becomes SIX of seven (only 64 + gnina fails, p 0.152) under App. C's seven-contrast family (D14). "Adding
rescoring at 128 moves top-15 from 185 to 198, gaining 17 and losing 4" → 198 to 213, gaining 19 and losing 4. "Raising
exhaustiveness 32 → 128 moves top-15 from 154 to 198" → 170 to 213. Parent-Vina-order carry (`:184`): 93 / 183 vs raw
94 / 185 → 129 / 196 vs raw 130 / 198, direction holds. Fig 19 caption (`:195`, "top step worth three at rank-1,
resolvable only at about depth 22"): the exh92 → exh128 step becomes +7 at rank-1 with uncorrected p 0.039 already at
rank-1. The plateau argument at `:246` COLLAPSES rather than renumbers: raw rank-1 "tops out at 95 above exhaustiveness
64" becomes 124 / 124 / 131 for 64 / 92 / 128 (a +7 rise), and "doubling the pool beneath the rescorer leaves it at
exactly 113" becomes 144 vs 151 (a 7-complex gap), so the "ranking ceiling, not sampling ceiling" sentence must be
rewritten, not re-numbered. The `:182` clause "the step over raw exhaustiveness 128 is the one contrast that would not
survive a correction over all twenty-eight pairs" is recompute (uncorrected 0.0026 × 28 = 0.073, so it likely still
fails, but the 64 + gnina contrast at 0.0044 × 28 = 0.12 now fails too).

AutoDock raw → gnina rank-1 step, two bases that must not be mixed (v1 mixed them):
- Table 19 basis (App. H :947, near-native WITHOUT validity): 31.4 → 37.3 %, 21 / 39, p 0.0273 becomes 43.2 → 49.8 %,
  19 / 39, p 0.0119; Holm over Table 19's sixteen tests 0.060 (ns).
- `:843` basis (PB-valid ∧ ≤ 2 Å): 31.0 → 36.6 %, 22 / 39, p 0.0396 becomes 42.9 → 49.2 %, 20 / 39, p 0.0183; Holm over
  App. C's seven contrasts 0.037 (survives). D14 decides which sentence quotes which.

Table 19 star flips (Holm over sixteen within-tool tests, near-native gate): DiffDock + smina k = 5 44.6 → 48.8 (5 / 18,
0.0468 *) becomes 50.2 → 54.1 (5 / 17, p 0.0169, Holm 0.0676, ns); DiffDock + gnina k = 5 likewise → ns; DiffDock gnina
k = 10 keeps its star (0.0245); AutoDock k = 1 stays ns (0.0596); all other marks hold.

Guided EquiBind (App. C :267, `:195`): the probe covers five of nine configurations (fpocket raw 179, P2Rank raw 180,
unguided raw 186, smina 192, gnina 199 form complexes), so "179-199" is a preview and the row is recompute; the probed
guided arms gain 3-5 complexes (fpocket-smina triple 7 → 10, P2Rank-smina 5 → 9, P2Rank-smina double 14 → 19), and
P2Rank-gnina sits at the "18" ceiling of `:195` today, so that clause is re-checked from the rebuild. "Both raw guided
arms recovering no near-native complex at all" (`:267`) becomes false: P2Rank raw recovers one over the pool.

### 1.6 Placement, cohort trims and the H.4 decomposition (`:408-422`, Fig 4, Table 27 :1377-1387, App. H :1188-1196)

Near-site (8 Å to the copy) PB-valid top-15 pool 6,169 → 7,565 poses; clipped-beyond-5 Å count is taken from the
regenerated `20d_…caption.txt` (the printed 2,329 is that script's own rule on the `pb_rmsd` axis; v1's "2,337" was the
`rmsd` column and is withdrawn). Rank-1 coverage after the trim 62.7 / 54.1 / 47.2 % → 82.2 / 73.9 / 51.2 %. AutoDock
in-place median 1.490 → 1.376 Å at rank-1, 5.043 → 4.922 at top-15, drift +3.55 unchanged; DiffDock +0.46 → +0.53;
EquiBind +1.09 → +0.96. EquiBind rank-1 near-site median: printed 2.9, generator 2.977 (rounds to 3.0, a pre-existing
truncation), nearest 3.239 → 3.2 Å (v1 wrote 3.3). Table 27 footnote (`:1387`): 1.490 / 1.543 and 1.467 / 1.566 →
1.376 / 1.491 and 1.532 / 1.585, EquiBind 3.239 / 3.276; "two of the 909 decisions" still two; 37.0 / 18.5 / 34.0 % →
49.5 / 19.1 / 43.6 %. All 117 Table 27 cells change. Near-native strata 319 / 1,715 / 498 → 451 / 2,169 / 515; all-check
validity within them 97.5 / 94.6 / 93.6 → 98.0 / 95.3 / 93.8 %; per-check intermolecular rates recompute. Depth-gain
decomposition (Table 25): [88, 2, 1, 87] / [63, 1, 2, 64] / [26, 2, 1, 25] → [65, 2, 1, 64] / [45, 0, 2, 47] /
[27, 2, 1, 26]. The interpretive paragraph at `:422` ("deeper ranks add valid but misplaced poses ... few for DiffDock")
is a wording site v1 missed.

### 1.7 Crystal-site reach and clusters (`:158`, `:427-476`, `:768`, `:792`, Table 3, Figs 5, 6, 38, App. H :1394-1399)

Pose-level rank-1 4 Å reach 183 / 172 / 137 → 239 / 225 / 143 under D2 (241 / 225 / 144 under a per-metric minimum).
Any-pose reach per headline arm on the canonical table: 88.4 → 92.4 % (AutoDock*), 85.1 → 90.1 % (DiffDock*). Pooled
three-tool oracle: 98.0 → 99.7 % under D5 (any copy), 99.3 % under D2 (copy j*); :476 is produced by the cluster report
and follows D5. The printed :476 values (97.7 / 88.1 / 84.5, medians 0.34 / 0.53 / 0.56) equal their generator exactly;
the generator differs from the canonical table because `pose_cluster_crystal_pocket_report.py:133` loads poses with
`Chem.MolFromMolFile(..., removeHs=True, sanitize=False)`, which does NOT strip hydrogens when `sanitize=False`, so pose
centroids include the explicit hydrogens the refined pose files carry (28 / 27 / 30 of 30 sampled AutoDock* / DiffDock* /
EquiBind* files) while `_ligand.sdf` is heavy-only; 215 of 303 AutoDock* centroids differ by > 0.05 Å (max 0.58) and
four complexes cross the 4 Å boundary. This generator defect drives every cluster-geometry number of the current thesis
(Table 3, Figs 5 / 6 / 38, `:474` precision-at-one, the 8 Å radius) and is fixed in Phase 0 before any regeneration.
Cluster-level reach, co-reach φ, "sixteen complexes with no qualifying cluster", "18 contribute no pose",
precision-at-one and the 95 % oracle ceiling are recompute under D5. The `:768` reach clauses ("led DiffDock on
native-pocket reach at every depth") join the recompute list.

### 1.8 Native-interaction recovery (`:165`, `:481-532`, Table 4, Fig 7, App. H :545-553, :1407-1415, Fig 39)

Recompute after per-copy crystal fingerprints exist. Direction: 105 / 119 / 11 top-5 PB-valid poses of AutoDock* /
DiffDock* / EquiBind* sit within 2 Å of an alternate copy and today score Jaccard ≈ 0 because fingerprint keys carry the
reference chain. Preview F1 0.565 / 0.529 / 0.471 → about 0.532 / 0.527 / 0.513. The union fingerprints of the top-5
(Fig 7 top-5 union, audit `:104-118`) raise a pose-SET question D2 does not settle: when a method's five poses are nearest
to different copies, the union is scored against the union of those copies' fingerprints. Adopt that rule and state it.

### 1.9 Cost per qualifying pose (`:708-760`, `:817`, `:848-850`, Table 6, Figs 15-16, Abstract `:63`, Kurzfassung `:98`)

Qualifying poses 311 / 1,622 / 466 → 442 / 2,068 / 483 (Vina-only 302 → 431); qualifying complexes 199 / 169 / 80 →
214 / 183 / 83; all-three shared set 53 → 54 (rank-1 all-three 27 → 28). AutoDock's median succeeding complex yields two
qualifying poses instead of one, so every charged / CPU / GPU median, IQR, Friedman, ratio and the "most expensive per
qualifying pose" ordering (`:744`, `:754`, `:850`, Abstract, Kurzfassung `:98-99`) is re-derived from the timing files.
Table 6 yaml `poses_near2` targets can be entered in advance (442 / 2,068 / 483).

### 1.10 Intention-to-treat (App. I :1074)

Re-scored with the probe engine over the 150 AutoDock exh128 + gnina poses of the five dropped complexes, ranked by
CNNaffinity: rank-1 none within 2 Å (closest 6.3 Å, 7FRX_O88 to an alternate copy); top-15 none within 2 Å but the
closest pose is 2.1 Å (7FRX_O88) instead of 6.4 Å; over all thirty poses 8F4J_PHO reaches 1.5 Å of an alternate copy at
a rank beyond 15. Base counts 149 / 308 and 213 / 308 against 132 and 179; margins 5.5 pp and 11.0 pp on 308. The 6.4 Å
and 150-pose figures have no registered generator; Phase 2 registers one.

### 1.11 Abstract and Kurzfassung

Abstract `:61-63`: five number swaps (36.6 / 34.0 / 18.2 → 49.2 / 43.6 / 18.8; +2.6 → +5.6; −4.2..+9.4 → −1.7..+12.8;
65.3 / 55.1 → 70.3 / 59.1; 80 of 303 → 83), "six in ten" → "six to seven in ten", the cost-ordering clause pending 1.9,
the "corrected local geometry but recovered no binding mode" clause re-verified (DiffDock raw → smina rank-1 near-native
132 → 137 under nearest), one optional convention clause. Kurzfassung `:91-100`: 36,6 / 34,0 / 18,2 → 49,2 / 43,6 / 18,8
(same length), 2,6 → 5,6, "−4,2 bis +9,4" → "−1,7 bis +12,8" (+1 char), 65,3 / 55,1 → 70,3 / 59,1, 80 → 83, "sechs von
zehn" → "sechs bis sieben von zehn" (+11 chars), cost clause `:98-99` pending 1.9. The compensating cut must be named
before the step is executable (Phase 5.5).

---

## 2. Step-by-step plan

Effort: T trivial, S small (< 1 h), M medium (half a day), L large (1-2 days). Every step, including Phase 5, names its
check. Order is a dependency order. All edits are located by quoted text.

### Phase 0. Decide, freeze, fix pre-existing defects (L, 1 day)

0.1 Adopt the decision record (D1-D17) in a dated note. Check: every later step cites a D-number.
0.2 Freeze baselines: copy `posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report/` to
`pose_comparison_report.bak-instance-<date>`; freeze beside it the six other harness inputs that change:
`docking_effort_gnina_v2/`, `docking_effort_gnina_v2_charged/`,
`cluster_crystal_pocket_matched_equibind/autodock_mgltools_exh128_gnina__diffdock_smina_allposes/`,
`pandamap_results/benchmark_matched_equibind/report/`, `posebusters_results/autodock_exhaustiveness_returns/`,
`thesis_latex/media/media/`, plus the WORKING-TREE `thesis_expected_values.yaml` (2,213 uncommitted lines against HEAD;
never `git show HEAD:`). Run `thesis_assertions.py` once through a 12-line dump wrapper (import `run_all`, write
`Report.results` as `key, expected, actual, ok, source` CSV) and store `harness_signature_0_instance.csv`. Check: 1,364
rows, 0 failures, exit 0; md5 manifest of the backup written beside it.
0.3 Fix the pre-existing defects that would otherwise be overwritten with their offset carried forward:
- `body_main_short.tex:113` "adds three recovered complexes": no generator exists and the `\ref` points at nothing;
  recompute gives two (AutoDock* +7PGX_FMN, EquiBind* +6XBO_5MC, DiffDock* 0), none removed, identical under both
  conventions. Write "two", add the supporting appendix sentence, register the `pb_rmsd` gate as a prose check.
- `:396` "p_holm ≤ 4e-10": companion report gives 8.997e-09 at rank-1. Write "≤ 9e-9"; add to yaml.
- `:414` "2.9 Å": generator 2.977. Write "3.0 Å".
- `:476` and every cluster-geometry number: fix the loader at `pose_cluster_crystal_pocket_report.py:133` (sanitize then
  `RemoveHs`, or `RemoveHs(sanitize=False)`), bump `_CACHE_SCHEMA` at `:3927`, regenerate the cluster report on the
  CURRENT convention, and re-verify :474 / :476 / Table 3 / Figs 5, 6, 38 before Phase 5. Fix the double-rounding
  formatter at `:2852` / `:2914` (sidecar prints 55.5 % from a 4-dp value where 168/303 = 55.4 %).
- NOT defects (do not touch): `:408` 2,329 (the script's own `pb_rmsd` rule) and `:781` vs `:197` (different depths).
Check: prose corrected at :113 / :396 / :414; regenerated cluster sidecars reproduce :476 from the canonical table
(268 / 258 / 297); yaml prose entries added for :113, :396, :474, :476.
0.4 Register the metal-stratum generator: `Scripts/Analysis/metal_stratum.py` (D4 rule; RDKit SDMolSupplier
`sanitize=False`, cKDTree, metal element set) and stage `bench_metal_stratum` (BITEXACT, CHEAP, needs
`bench_pose_comparison` for the 303-id cohort) writing `metal_stratum.csv` (protein, copy, is_ref, distances, nearest
element / residue, membership at 4 / 5 / 6 Å) and `metal_stratum_contrasts.csv` (rank-1 and top-15 by stratum under both
conventions, McNemar and Fisher). Yaml asserts 78 / 73 / 81, the two crossover complexes, and the cofactor / ion split
the rule yields (16 / 62; App. C :165 prints 15 / 63, to be reconciled). Check: reference-convention output reproduces
App. C :165 verbatim (37.8 / 39.6 / 33.3 / 17.9 %; 134 vs 149 of 225, p 0.142; 33 vs 49 of 78, p 0.020).
0.5 Inspect the eight rank-1 flips with nearest-copy RMSD in (1.5, 2.0] Å (four per headline arm; M2 note :35) against
`pb_rmsd_nearest` and a manual overlay, and record whether any is a symmetry artefact. Check: a written verdict per pose;
if any is an artefact the derived 149 / 132 drop accordingly before Phase 5.
0.6 Adjudicate the remaining minor audit items (Section 8, `audit__*.json`) and fold accepted ones into Phases 1-7.

### Phase 1. Hub code, `Scripts/Analysis/posebusters_pose_comparison.py` (L, 1 day)

1.1 `:1149` add `load_all_mols(sdf_path)` beside `load_first_mol`, same sanitize fallback per record.
1.2 `:1520-1528` load `copies_sdf = cdir / f"{pdb_id}_ligands.sdf"` (fall back to the single file if absent);
`copies_h = [RemoveHs(reassign_template(m, crystal_h) or m) for m in load_all_mols(...)]`; `ref_copy_index` by
coordinate identity with `crystal_h` (0 on this dataset).
1.3 `:1531-1532` native contacts per copy: `native_contacts[j]`, `n_native[j]`. Note the hub's contact and PLIF
receptor is `cdir / f"{pdb_id}_protein.pdb"` (`:1521`, `:1530`, `:1534`), the PoseBusters-shipped file with every chain,
so alternate copies are always covered (a different file from the PandaMap receptor).
1.4 `:1560-1562` and `compute_plif_recovery` (`:1383-1415`): pass all copies as the leading reference ligands of one
ProLIF run, return an (n_poses × n_copies) matrix, select column j*; a copy failing `_ligand_for_plif` (`:1391`) gets a
NaN column instead of zeroing the pair.
1.5 `:1566` `r = [symmetry_rmsd(pose, c) for c in copies_h]`; `j_star = int(np.argmin(np.nan_to_num(r, nan=np.inf)))`;
`rmsd = r[j_star]`; `rmsd_ref_instance = r[ref_copy_index]`; in the `except` branch at `:1568` set
`j_star = ref_copy_index`.
1.6 `:1567`, `:1572`, `:1581`, `:1584`, `:1585`, `:1588` evaluate centroid, `posebusters_rmsd` (incl.
`pb_rmsd_within_2A`), `rigid_body_fit`, `bestfit`, `torsion_metrics` and contact recovery at `copies_h[j_star]` (D2).
Keep `*_ref_instance` twins for rmsd, centroid, pb_rmsd, bestfit.
1.7 `:1459-1511` `PoseRecord`: append the D6 fields at the END. `:1592` wire them.
1.8 `:1513-1514` and the work-tuple builder at `:13507-13508` (not :13494): thread `reference_convention`; under
`instance` set j* = ref_copy_index so the old table is reproduced exactly.
1.9 `:13346` argparse `--reference-convention {instance,nearest}`, default `instance` (kept until 7.3); `:13225` schema
5 → 6 and add the convention to `_per_pose_signature` (`sig` is computed at `:13475`); `:13187-13191` add
`rmsd_ref_instance`, `nearest_copy_index`, `n_copies`, `reference_convention` to `_CACHE_REQUIRED_COLS`; `:13243` print
WHICH required columns are missing; `:13261-13276` make `--reuse-cache` read the manifest and refuse a cache whose
convention differs, and make the non-reuse signature mismatch REFUSE when only the convention differs unless `--force`
(schema-only or pb_csv-only mismatches keep today's recompute). Add a unit test for both refusals.
1.10 Second crystal-load site `:9323` and `:9351-9352` (fig 20g exemplars and receptor-context PDBs `:9470-9540`): load
record `nearest_copy_index` of `_ligands.sdf`; legend "nearest deposited copy" instead of "crystal (original)".
Docstrings `:46-55`, `:503-528`, comments `:661`, `:723-732`, argparse help `:13387-13390`; the eight sidecar text
writers print the convention string from the manifest.
1.11 Regression under `--reference-convention instance --force` into a scratch out-dir: EVERY PoseRecord float / bool
column equals the frozen table (rmsd, pb_rmsd, pb_kabsch_rmsd, bestfit_rmsd, centroid_dist to 1e-9 Å; contact_recovery,
plif_recovery to 1e-12; torsion columns and `pb_rmsd_within_2A` exactly) on all 241,913 data rows (`wc -l` 241,914
includes the header); `oracle_summary.csv` byte-identical (it gains no columns). Unit test: a single-copy id yields
`rmsd == rmsd_ref_instance` and `nearest_copy_is_ref` True. Check: both pass.
1.12 Add `--reference-convention nearest --force` to the four hand-maintained invocations now:
`_build_reproduction_notebook.py:393-400` (stage cmd), `Master_Docking_AD_Full_Protein.ipynb` cell 33
(`:10658-10701`, also switch its out-dir per D15), `REGENERATE_NOTES.md:935 / :954` (also correct their stale
`autodock_gnina` pin), and the hub docstring examples `:514-524`. Record `pb_rmsd_loadall` once during the first rebuild
as a parity column (expected: ≤ 2 Å verdicts agree 100 %, values differ only by terminal-group symmetrisation; CalcRMS
argmin and symmetry_rmsd argmin pick different copies on 3 of 24,486 multi-copy poses), drop it if so.
1.13 Pocket localisation (`:10998-11044`, cutoff `:733`): `centroid_dist` at j* (D2); the switch point for an
any-copy variant would be `:11025-11026`.

### Phase 2. Downstream scripts with their own crystal logic (L, 1.5 days)

2.1 `pose_cluster_crystal_pocket_report.py`: `:921` `crystal_centroid` returns the list of copy centroids;
`:1107-1109` a cluster is correct if within `--match-thr` of ANY copy (D5) with the six adjacent-copy ids counted once;
`:553` `_cluster_purity`, `:562` `_precision_at_1`, `:614 / :624` `_ensemble_stats` (the pooled oracle behind :476),
`:944` `_pocket_crystal_stats`, `:963`, `:1020`, `:1314-1315` example scatter (draw every copy), `:4064`: min over
copies. Bump `_CACHE_SCHEMA` (`:3927`) again or pass `--force`, because the analysis signature fingerprints the CSV but
not the code and keys the crystal directory by path only. Check: on the 165 single-copy ids every per-complex column
of `cluster_summary` (reach, purity, rank_in_correct_cluster, has_crystal) equals the Phase 0.3 run; co-reach by a
scratch recomputation.
2.2 `run_pandamap.py:963`: one crystal task per record of `_ligands.sdf` (`pose_name .../crystal_k`). The receptor
check is DONE: for all 395 alternate copies of the 138 analysed multi-copy ids the canonical
`hetatm_cleaned/<ID>_protein.pdb` keeps every hosting chain (no copy has < 10 receptor atoms within 4.5 Å). The gate at
`:1174-1176` recomputes crystal references only when `crystal_interactions.csv` is absent or `overwrite` is true, and
config 07 sets `overwrite: false` with the file dated 2026-08-14: add an explicit pre-step that moves
`crystal_interactions.csv` and `crystal_pose_summary.csv` aside (or a crystal-only refresh switch; `--overwrite` would
re-profile every docked pose at `:1161`). Add both crystal CSVs to the stage outputs (builder `:528`).
2.3 `pandamap_interaction_report.py`: every `crystal.groupby(["protein","ligand"])` (`:1185`, `:1484` similarity
matrix, `:1729`, `:2044`, `:2435` per-complex figure) selects `crystal_k` by the pose's `nearest_copy_index` joined from
`per_pose_metrics`; top-5 union scored against the union of the selected copies' fingerprints (1.8).
`interaction_pose_basis_audit.py:227` and `:406` likewise (they read `crystal_interactions.csv` directly at `:184`).
Check: F1 on single-copy complexes equals the frozen run.
2.4 Pure re-readers needing no code change: `posebusters_validity_report.py:748`, `docking_effort_comparison.py:256`,
`topk_recovery_gnina_arm.py:130` (but its report embeds raw-arm literals at `:305-306` and `:331-333` and a docstring
`:36-40` that must be regenerated or removed), `thesis_endpoint_diagnostics.py:77`, `optimization_benefit_stats.py:115`
(the gate; `:323` is an import), `rank_success_analysis.py:139 / :286`, `ligand_docking_difficulty.py:99 / :112`,
`receptor_docking_difficulty.py:93`, `diffdock_gnina_rerank_analysis`. `filmstrip_rank1_top5_top15.py:313` reads the
cache through the hub validator and therefore inherits the tightened required columns; move its `DEFAULT_REPORT`
(`:38-39`) to the matched tree. `boxed308_raw_arm_metrics.py:123` optional parity branch. `rebust_hydrogen_normalised.py:67`
is pinned to the full-protein tree (D15). `Scripts/Analysis/tests/test_diffdock_sweep_gate.py:29` points at a notebook
that now lives under `obsolete/` and cannot load; fix or retire.
2.5 Re-pin `autodock_exhaustiveness_returns.py`: `CASCADE_PINS` (`:214-217`) hard-codes the published cascade rows for
exh18 and exh92 and `main()` raises `SystemExit("CASCADE REGRESSION ...")` at `:1373-1375` on any 1-pose difference, so
the stage ABORTS after the rebuild. Key the pins by the manifest's `reference_convention` (keep the instance pins as
`CASCADE_PINS["instance"]`) and fill the nearest pins from the rebuilt `pose_validity_cascade.csv` in step 3.1a. Note the
`:570` coupling check `bestfit_rmsd == pb_kabsch_rmsd` is warn-only (written to `exh_gate_input_audit.csv`), not an
assert; read the CSV. Never invoke the script bare: `DEFAULT_PER_POSE` (`:222-223`) points at the full-protein tree.
2.6 Register the intention-to-treat re-score of the five dropped complexes as a stage (probe `work()` over
`Dockings/vina_results_full_protein_vina_scoring_mgltools_exh128/<ID>/mgl_tools/docking/optimized_gnina/`, ranked by
CNNaffinity) with yaml asserts (rank-1 closest 6.3 Å, top-15 closest 2.1 Å, one pose ≤ 2 Å beyond rank 15).
2.7 Register `bench_reference_convention` (section 3, BITEXACT, CHEAP, needs `bench_pose_comparison`) writing the D9
table from the `*_ref_instance` columns (no re-scoring from SDFs); archive the probe and the five wf_plan helper
scripts read-only under it; rewrite `nearest_copy_probe/README.md:5` (after the switch `rmsd_ref_re` reproduces
`rmsd_ref_instance`).

### Phase 3. Rebuild (compute ≈ 1-2 h estimated, no wall-time record exists; attention M)

3.1 Regenerate `Thesis_Reproduction.ipynb` / `REGENERATE.md` from the edited builder, then run `bench_pose_comparison`
in harness force mode (`RUN_MODE="force"`, `FORCE_STAGES`, builder `:94-95`); the `--force` and the convention flag live
in the stage command (1.12), because `repro_harness.py:391-406` runs `st.cmd` verbatim. Then rebuild the
`benchmark_full_protein_vina_scoring` report tree the same way (D15). Check: both manifests carry
`reference_convention: nearest`, `schema: 6`; `find posebusters_results -name per_pose_metrics.manifest.json | xargs grep
-L reference_convention` lists only the boxed tree.
3.1a Re-pin the ladder check: read rows `autodock_mgltools_exh18` and `autodock_mgltools_exh92` of the rebuilt
`pose_validity_cascade.csv` into `CASCADE_PINS["nearest"]`. Check: pins 1-5 unchanged (303, 8986, 303, 8974, 99.9 and
303, 8984, 303, 8954, 99.7); pins 6-8 equal 144 / 219 / 2.4 (exh18) and 209 / 407 / 4.5 (exh92), pre-computed; pins 9-15
differ from the instance values by no more than the D2 Kabsch effect measured on exh128 (kabsch1 poses −19 of 2,848,
≤ 47 verdict flips). If pins 6-8 do not match, the rebuild did not apply the rule and 3.2 must not start.
3.1b Run `bench_exhaustiveness_returns` with its stage command. Check: `exh_gate_input_audit.csv` row
`bestfit_rmsd == pb_kabsch_rmsd` is True with max difference 0.0; stdout prints the cascade self-check success line.
3.2 Then, in order: `bench_validity_report`, `bench_filmstrip`, `bench_topk`, `bench_endpoint_diagnostics`,
`bench_optimization_benefit`, `bench_ligand_difficulty`, `bench_receptor_difficulty`, `bench_clusters` (after 2.1's
schema bump), `bench_effort_charged`, `bench_effort_elapsed`, the crystal-CSV move-aside then FORCE `bench_pandamap`
(the crystal fingerprints are produced INSIDE this stage; docked poses resume-skip at `:1157-1161`), `bench_pandamap_report`,
`bench_interaction_audit`, `bench_diffdock_rerank`, then the new `bench_metal_stratum`, `bench_itt`,
`bench_reference_convention`. That is 16 of the 19 section-3 stages plus 3 new (only `bench_posebusters`,
`bench_posebusters_matched`, `bench_figure1` stay); `FORCE_STAGES` carries all of them.
3.3 Re-render collapsed figures with `--reuse-cache --collapse-plots-only --reference-convention nearest` (the flag pair
is REQUIRED; the manifest check now protects it).
3.4 Change signature (new). Run the UNCHANGED harness (old yaml) against the rebuilt tree through the 0.2 wrapper →
`harness_signature_A_oldyaml_newtree.csv`. Acceptance: (a) every failing key lies in a yaml block Section 4 names as
changing; (b) zero failures in `cohort`, `determinism`, `table_5`, `table_9`, `table_10`, `table_11`, `table_12`,
`table_14`, `table_16`, `table_17`, `prose.fig23`, `prose.fig24`, `figures[coverage]` and the 25 figure rows outside the
14 regenerated ones (about 290 rows with a hard zero requirement; the changing rows are bounded above by ≈ 1,074, v1's
"~850" was an estimate); (c) `table_1` shows exactly 68 failures with `poses` / `valid_*` never among them; `table_22`
exactly 6 (DiffDock raw, k = 1 and k = 5: 47 → 65, 58 → 60); `table_8` exh128 rows show 130 / 198 / 217 and
149 / 213 / 214 at d = 1 / 15 / 30; (d) the `actual` column of every failure equals the value Section 3 prints as new.
Any failure outside (a) stops the plan. Run this BEFORE 6.1 so the 14 figure rows fail for the right reason.

### Phase 4. Selection re-check (S)

4.1 From the rebuilt sidecars confirm exh128 + gnina still selected (213 at top-15, beats all seven), smina still
carried (144 vs 143 triple, 179 vs 177 double), unguided gnina still the EquiBind winner under `BEST_VARIANT_METRIC`
(`:2136`). If any flips, stop and re-plan the arm's text.
4.2 Confirm `_TOPK_TEST_PAIRS` (`:9764-9797`, which also contains DiffDock gnina-opt vs EquiBind gnina-opt at `:9788`):
vs DiffDock raw 77 / 58 at k = 1 (uncorrected 0.12, "genuine tie" holds), significant from k = 5. Confirm the App. C
seven-contrast rank-1 family (six of seven survive) and the Table 19 sixteen-test family (star flips of 1.5) and apply D14.

### Phase 5. Thesis text (L, 2-3 days; short build only; locate by quoted text)

5.1 Methods: `:109` redefinition (draft below), `:113` "two" plus the generator (0.3), `:115`, `:119`, `:124-133` form
definition names the copy rule, `:158`, `:165` interaction reference and the union rule, `:170` cost definition.
Check: each edited sentence re-read against D1 / D2 / D11; `grep -c 'crystallographic ligand'` on the Methods lines
returns only the retained hits recorded in the 8.2 table.
5.2 Results: every row of Section 3.2; captions of Figs 2-7, 15, 16. Check: after the phase, grep the short build for
each old literal of the phase's rows; the hits must equal the RETAINED rows of the 8.2 table and nothing else.
5.3 Discussion, Conclusions, Limitations: `:768` (interval and reach clauses), `:777-781` (five; "at the top-15
selection depth"), `:788-794`, `:806`, `:817`, `:824`, `:834-836`, `:843-850` (D14 wording), `:857`, `:866-870`; add the
convention disclosure sentence at `:868` (143 / 308, 138 / 303, median 36.0 Å over the 143, source benchmark scored the
same way, single-instance kept as sensitivity). `:697` (Section 4.2) "Calibration point estimates show the same
ordering, although not a statistically separable difference" is true at rank-1 only under both conventions and must
say "at rank-1". Check: as 5.2.
5.4 Appendix: `:165` (D4 numbers: 51.1 / 43.6 and 49.3 / 26.9 % rank-1; top-15 −8.4 pp p 0.056 and −19.2 pp p 0.028;
Fisher 0.440; 4 Å −9.1 / −17.8, 6 Å −8.1 / −19.8; failure shares 58.2 / 56.8 and 70.2 / 70.2 %; the two crossover
complexes; the residue-level rule wording; 16 / 62), `:167` (add "every deposited copy of the ligand lies inside this
box"), `:173-184` (Section 1.5 values incl. six of seven, 1-25, 129 / 196), `:189-195 / :235-241` captions, `:218-227`
ladder and footnote, `:230` (D10 scope), `:244-248` recompute, `:246` rewrite (plateau collapses), `:267-269`, `:276`,
`:487` optional clause, `:545-553` interaction, `:664` footnote roles, `:666` convention paragraph (draft below), `:722`
(name the 46 % share beside Fig 30), `:784` (D13 clause), `:871` lead-in must name the gate its tables use (near-native
WITHOUT validity; today it says "validity-aware"), `:893-920` Wilson table, `:947-964` Table 19 incl. star flips,
`:992-999` bands, `:1026-1031` cross-tool, `:1056-1062` Table 22, `:1074` ITT, `:1172-1181` PB decomposition,
`:1188-1196`, `:1203` regime clause (softened per D11), `:1207`, `:1247-1250` Table 25, `:1258-1260` H.9 (11.7 margin,
8.8, 42.6 %, 1,100, "within 12 points"), `:1265`, `:1303-1340` Table 26, `:1377-1387` Table 27 and footnote,
`:1394-1399` cluster statistics, `:1407-1415` interaction. New D9 sensitivity table with Form and Combined columns and
the single-copy row (53.3 / 48.5 / 32.1 %) beside `:1074`. Check: as 5.2, plus the regenerated sidecar for every
table value.
5.5 Abstract and Kurzfassung per 1.11. Name the compensating German cut before editing (candidate: shorten the
optimisation clause at `:89-90`, which is length-neutral in meaning). Check: `wc -m` of the Kurzfassung page unchanged
within the cut; pymupdf raster shows no overflow.
5.6 House style. Check: `grep -n '; \|: \|—'` over the edited lines returns no prose hits.

Draft for `:109` (D11): "Near-nativeness is the symmetry-corrected heavy-atom RMSD between a docked pose and the
nearest of all crystallographic copies of the ligand deposited with the structure, computed in place without
superposition. The copy with the lowest RMSD is chosen for each pose and supplies the centroid, form, contact and
interaction references for that pose. This is the convention of the source benchmark, whose validity battery loads every
deposited instance and reports the lowest RMSD. The 165 single-copy complexes are unaffected. The single-instance value
is retained and reported as a sensitivity arm in Appendix~\ref{...}."

Draft for `:666`: "In 143 of the 308 entries the deposited file holds more than one copy of the ligand, 138 of them
among the 303 analysed complexes, and the nearest other copy lies a median 36.0 Å from the reference over those 143.
Every pose is scored against the nearest deposited copy, and the form, centroid and interaction references follow that
copy. The source benchmark scored every method this way. Its 25 Å cube centred on the reference kept almost every other
copy outside Vina's search, whereas the whole-protein box used here contains every one of them, so the convention is
live for all three tools in this study. The choice of copy can lower a form count where the copies differ in
conformation, and it does so for one complex in each of three reported variants. Four alternates do not share the
reference pocket (7A9E_R4W, 7TUO_KL9, 7VKZ_NOJ, 7Z1Q_NIO) and no complex gains recovery through any of them."

### Phase 6. Figures and media (M, 0.5 day)

6.1 Regenerate the 14 convention-dependent figures: image3 (Fig 2), image4 (Fig 3), image5 (Fig 4), image6 (Fig 5),
image7 (Fig 6), image8 (Fig 38), image9 (Fig 7), image11 (Fig 39), image24 (Fig 15), image25 (Fig 16), image44-47
(Figs 18-21). Producing stages: hub (`18_topn_*`), `bench_filmstrip`, `bench_clusters`, `bench_pandamap` +
`bench_pandamap_report`, `bench_effort_charged` / `bench_effort_elapsed` (image25 has two byte-identical sources; the
yaml pairs it to the elapsed dir), `bench_exhaustiveness_returns`. Hand-copy into `thesis_latex/media/media/` and
confirm by md5 against the producing file. Check: a second harness signature run after the copy flips exactly the 14
`figure[imageN]` rows from fail to pass and nothing else.
6.2 Captions per Section 3. `validate_regeneration.py` byte-gates nine generators against the CURRENT canonical table
and only DECLARES Figs 2, 3, 5, 6, 38 (its stated reason, an out-dir-keyed cache, is false); add a Case for
`pose_cluster_crystal_pocket_report` with `--force` so those three regain a byte gate. The convention itself is asserted
only through the manifest check (1.9 / 3.1) and `check_reference_convention` (7.2); say so in Risk 6.

### Phase 7. Harness, notebook, documentation (L, 1.5 days)

7.1 `thesis_expected_values.yaml`: bump `meta.recorded`, add `meta.reference_convention`, state the column convention in
the gate comment (`:56-58`); re-transcribe FROM THE PRINTED THESIS every changing block (`table_1` 68 rows, `table_2`,
`table_3`, `table_4`, `table_6`, `table_8`, `table_18`, `table_19`, `table_20` (24 of 36), `table_21`, `table_22`
(6 rows), `table_24` near half, `table_25`, `table_26`, `table_27`, `prose.refiner_contrasts`, `prose.appendix_h9`
(`n_needed: 1100`), `prose.figure_39`, `prose.audit`); add blocks for the three new stages and the four Phase 0.3 prose
values; `figures.pairs` holds paths, not md5s, so nothing there changes. Then run the wrapper (new yaml, new tree) →
`harness_signature_B.csv`; the set of keys whose `expected` changed must equal the failing set of signature A and no
other key may change. Add `ROOT = Path(os.environ.get("THESIS_ASSERT_ROOT", ...))` at `thesis_assertions.py:42` and run
the new yaml against a symlink-farm root pointing at the 0.2 backups → `harness_signature_C_newyaml_oldtree.csv`, whose
failing set must equal A's. Transcription DIRECTION is proven separately: for every changed key, grep the new `expected`
literal in the tex file its `source:` names; a value present in the outputs but absent from the tex is rejected.
Check: A, B, C and the tex-presence list archived under `nearest_copy_probe/wf_.../harness_signatures/`.
7.2 `thesis_assertions.py`: docstring `:12-30` (202 / 199 → 217 / 214, convention sentence); comments `:1583` (0.12),
`:1596-1609` (4,200; 8.2; "only rank-1 is equivalent"); no gate-line change under D6; add `check_reference_convention`
asserting `manifest.signature.reference_convention == yaml meta.reference_convention` and `schema == 6` BEFORE the
numeric block (today nothing reads `meta.*` and verify mode checks presence only), plus the checks for the new stages;
register them in `run_all` (`:808-841`).
7.3 Flip the argparse default at `:13346` to `nearest` and delete the explicit flag from the same four invocation sites
in the same commit. `_build_reproduction_notebook.py`: `:386-400` dated CHANGED note; `:345` and `:392` "nine downstream
stages" → "fourteen" (thirteen declare the need, the audit depends transitively); `:484` register the three new stages
(registry 60 → 63; REGENERATE.md §3 19 → 22); `:528` add the crystal CSVs to `bench_pandamap` outputs; `:1138-1168`
section 10 paragraph; run `validate_regeneration.py` after 6.1 so `:1231-1245` records the re-based result; rebuild and
execute the notebook so `REGENERATE.md` regenerates (never hand-edited). `regenerate_guide.py:323` indexes only
`table_*` keys: key the D9 block `table_28` or extend the filter; re-check the hard-coded manual-table sentence at
`:341`. Check: the invocation grep of Section 8 lists exactly the known sites and none carries the flag.
7.4 Documentation: new dated `REPRODUCTION_COVERAGE_<date>.md` and glossary version (entries :121-129, :141-159,
:170-191, :205-214, :224-240, :276-299, :314-331, :496-507, :561-598, :647-662 re-derived, not copied); living files
`README.md:15 / :67` (the assertion count 257 is stale today, the harness runs 1,364; the stage count 60 is current
until 7.3), `THESIS_REPRODUCTION.md:57 / :76 / :380-399` ("225 of 225", "58 stages", the gate table's 202 / 199),
`thesis_latex/README.md:12`; dated banners (house pattern) in `REGENERATE_NOTES.md:69-71`; no edit to `CLEANUP_PLAN.md`
(dated record) and say so for the 8.2 sweep.

### Phase 8. Build and QA (M, 1 day)

8.1 `latexmk -f Thesis_short.tex`; pymupdf raster of every regenerated figure page and the Kurzfassung page.
8.2 Registry sweep, generated by script from the Current column of Sections 3.1-3.3 and Section 1, plus the appendix
descriptive values 98 / 128 / 149 / 159 / 164 (`:230`), "only two of the seven" and "top-3 to top-22" (`:182`),
93 / 183 (`:184`), 1.490 / 1.543 / 1.467 / 1.566 / 37.0 / 18.5 (`:1387`), 31.4 / 37.3 / 21 / 39 (`:893`, `:897`,
`:947`, `:1026-1031`). Every hit in the three short-build files is entered as CHANGED (with the new value) or RETAINED
(with the reason; Orai chapter, Vina weights, PoseBusters check table, dataset descriptors and LaTeX float fractions
account for nearly all retained hits). Two-digit literals cannot be swept and are covered by the quoted-anchor checks
of Phase 5.
8.3 `thesis_assertions.py` zero mismatches, accepted only together with signatures A, B, C and the tex-presence list.
8.4 Re-read Abstract, Kurzfassung, Results near-site paragraph (`:422`), Discussion and Conclusions for claims whose
direction could have moved (cost ordering, six-in-ten, equivalence within ten points, "adds two complexes", "ties with
64", the `:246` ceiling argument, `:697`).
8.5 Invocation check: `git grep -nE '(VINA_PY,|"run",|python)[^#]*posebusters_pose_comparison\.py' -- '*.py' '*.ipynb' '*.md'`
lists exactly `_build_reproduction_notebook.py:393`, `Thesis_Reproduction.ipynb:430`, `REGENERATE.md:359`,
`Master_Docking_AD_Full_Protein.ipynb:10658`, `REGENERATE_NOTES.md:935 / :954` and the hub docstring; none carries the
flag after 7.3.
8.6 Update project memory: the convention switch, the new stages, the traps met, the Phase 0.3 defects.

---

## 3. Every affected thesis location

W = wording must change. Values are derived previews (D16) unless marked recompute. Lines as of 22:30 on 2026-09-07;
locate by quoted text.

### 3.1 `thesis_latex/Thesis_short.tex`

| Line | Gate | Current | New |
|---|---|---|---|
| 61 | def. | Abstract endpoint definition | W optional |
| 63 | valid ∧ ≤ 2 | 36.6 / 34.0 / 18.2 %; +2.6; −4.2..+9.4; six in ten; 65.3 / 55.1; 80 of 303; "corrected local geometry ..."; "most expensive per qualifying pose" | 49.2 / 43.6 / 18.8; +5.6; −1.7..+12.8; six to seven in ten; 70.3 / 59.1; 83; re-verify; pending 1.9 |
| 80 | def. | Kurzfassung definition | no change |
| 89-90 | — | optimisation clause | verify; candidate for the compensating cut |
| 91 / 92 / 94 / 95 / 97 / 100 | valid ∧ ≤ 2 | 36,6 / 34,0 / 18,2; 2,6; −4,2 bis +9,4; sechs von zehn; 65,3 / 55,1; 80 | 49,2 / 43,6 / 18,8; 5,6; −1,7 bis +12,8; sechs bis sieben von zehn; 70,3 / 59,1; 83 |
| 98-99 | cost | AutoDock am teuersten | W only if 1.9 flips the ordering |

### 3.2 `thesis_latex/body_main_short.tex`

| Line | Gate | Current | New |
|---|---|---|---|
| 88 | stratum | 78 of 303 within 5 Å of the crystal ligand | 78 (D4); W "of the deposited reference instance, measured to the metal-bearing residue" |
| 109-115 | def. | Methods definition | W: D11 draft; `:113` "three" → "two" (0.3) |
| 119, 124-133 | def. | clustering pocket; form definition | W: nearest copy; copy fixed by the in-place minimum |
| 140 | — | footnote −1 / −2 / −2 DiffDock | recompute; state which file was passed as `mol_true` (no run_posebusters path sets it) |
| 158, 165, 170 | def. | 4 Å reach; fingerprint; cost definition | W; union rule at :165 |
| 191-199 | triple / double / pooled | 164 vs 160; 138 / 138, 133, 5 / 5; 55 vs 48; 160; 133 vs 132; 834 vs 828; nine, four; 176 vs 173; 237 vs 235; footnote 32.7 / 33.7 / 32.3 | 178 vs 175; 148 / 149 agree 145, 3 / 4; 57 vs 50; 175; 144 vs 143; 1,031 vs 1,022; seven, three; 186 vs 186 (W: "level", not "gnina higher"); 236 vs 235; footnote recompute (no-validity basis 35.3 → 45.2 %) |
| 205-214 | — | provenance comments | 4,016 / 3,887 / 3,854; 224 / 236 / 235; name the nearest-copy columns |
| 240 | — | "unidock2 row deliberately omitted" | unchanged |
| 264, 277-283 | — | Table 1 headers | W optional footnote |
| 286-304 | pose-level | Table 1 rows | Section 1.3 |
| 324 | valid ∧ ≤ 2 | at most four | five |
| 326 | def. | Table 2 lead-in | W |
| 347-349 | valid ∧ ≤ 2 | Table 2 rows | 49.2 % (149) / 70.3 % (213) / 70.6 % (214) / +64 (+21.1 %) *** / +1 (+0.3 %); 43.6 % (132) / 59.1 % (179) / 60.4 % (183) / +47 (+15.5 %) *** / +4 (+1.3 %); 18.8 % (57) / 27.4 % (83) / 27.4 % (83) / +26 (+8.6 %) *** / +0 |
| 353-372 | — | footnote and comments | update counts |
| 380 | def. | "as-placed crystal RMSD" | W |
| 382 | valid ∧ ≤ 2 | Q p < 1e-7; 36.6 / 34.0 / 18.2; +2.6 (−4.2, +9.4) p 0.509; +13.2 | holds; 49.2 / 43.6 / 18.8; +5.6 (−1.7, +12.8) p 0.159; +13.2 unchanged |
| 386 | valid ∧ ≤ 2 | 112; 9.8; +2.6 | 129; 10.5; +5.6 (power 0.32) |
| 388 | valid ∧ ≤ 2 | +28.7 (23.9-34.0); +21.1 (16.9-26.1); +8.3 (5.7-11.9); 36.6 → 65.3; 202 vs 199 | +21.1 (16.9-26.1); +15.5 (11.9-20.0); +8.6 (5.9-12.3); 49.2 → 70.3; 217 vs 214 |
| 393, 401 | — | Fig 2 / 3 captions | W; regenerate image3 / image4 |
| 396 | valid ∧ Kabsch | 74.3 / 74.6; 95.4 / 97.0; 52.5 / 52.1 / 33.3; 47 / 46, 30 / 31, 25 / 26; "p_holm = 1.000 at every depth"; "p_holm ≤ 4e-10" | unchanged; unchanged; unchanged; 48 / 47, 30 / 32, 25 / 26; "1.000, 0.899, 1.000" or "ns at every depth"; ≤ 9e-9 now (0.3), recompute after |
| 404 | valid ∧ Kabsch / ≤ 2 | 33.3 → 54.1; 20.8; "at most 0.7"; 2.3 / 2.3 / 3.0 | 33.3 → 53.8; 20.5; "at most 1.3"; 2.6 / 2.3 / 2.3 |
| 406 | valid ∧ ≤ 2 | 87 and 64; six in ten | 64 and 47; six to seven in ten |
| 408 | trim | 6,169; 2,329; "crystal-ligand site" | 7,565; regenerated caption value; W |
| 410-414 | trim | 62.7 %; 1.5 → 5.0; +3.55; 48 → 78; 27 → 8; 5.72; +0.46 / +3.55 / +1.09; 1.5-1.9; 137; 2.34; 2.9 Å; 26.4 / 26.4 | 82.2 %; 1.4 → 4.9; +3.55; recompute; recompute; recompute; +0.53 / +3.55 / +0.96; 1.5-2.1; recompute; recompute; 3.0 now (0.3) → 3.2; 27.4 / 27.4 |
| 419 | — | Fig 4 caption | regenerate image5 |
| 422 | trim | near-site interpretive paragraph ("few for DiffDock", "form-limited only because") | W: re-read against the regenerated 20d sidecar |
| 427-435 | cluster | sixteen; Fig 5; 26 / 24 / 4 pp; φ +0.23..+0.28; 2.8e-4 | recompute (D5); regenerate image6 |
| 455-461 | cluster | Table 3 reach and co-reach | recompute |
| 469-474 | cluster | "18"; six / five poses, 81 / 83 / 47 %; 95 %; 56 %; 56.8 / 55.4 / 55.8 | recompute; regenerate image7 (55.4 is right today, sidecar rounding fixed in 0.3) |
| 476 | cluster (D5) | 97.7 / 0.34; 88.1 / 0.53; 84.5 / 0.56; 97.7 / 0.37; 1.23 | loader fixed first (0.3); then 99.7 %; 92.4; 90.1; medians recompute |
| 481-532 | interaction | Fig 7; 0.450 / 0.417 / 0.380; F1 0.57 / 0.53 / 0.47; Table 4; 21.0 / 19.4 / 15.3 / 16.4 | W plus recompute; regenerate image9 |
| 697 | valid ∧ ≤ 2 | "Calibration point estimates show the same ordering, although not a statistically separable difference" (Section 4.2) | W: "at rank-1" |
| 708 | def. | cost definition | W |
| 715, 751 | cost | Fig 15 / Fig 16 captions: 199, 169, 80; 53 | 214, 183, 83; 54; regenerate image24 / image25 |
| 737-741 | cost | Table 6 rows | medians / IQR recompute; complexes 214 / 183 / 83 / 217; poses 442 / 2,068 / 483 / 431 |
| 744-760 | cost | Friedman 85.2 / 3.2e-19 / 0.80; 53; ratios; "yields one"; 311 / 1,622 / 466; 21.7; "nine"; 1,805.0 vs 1,752.8; 271 / 479 / 38; "most expensive" | recompute; 54; recompute; "two"; 442 / 2,068 / 483; recompute; "eleven"; recompute; recompute; W if ordering flips |
| 768 | valid ∧ ≤ 2 / cluster | 4.2 / 9.4; "led on native-pocket reach at every depth" | 1.7-point DiffDock lead to 12.8-point AutoDock lead; reach recompute (D5) |
| 777 | valid ∧ ≤ 2 | at most four (Section 5.2) | five |
| 779 | — | "fixed-rank near-native rate barely changed" | holds (+0.4 to +0.8 pp for DiffDock smina) |
| 781 | triple / double / near | 133 to 132; "exactly level"; "gnina level or ahead"; 83 vs 61; nine | 144 to 143; "smina ahead by two"; W as in 1.5, add "at the top-15 selection depth"; 86 vs 64; seven |
| 788, 857 | valid ∧ ≤ 2 | eight pp; six in ten; nine in ten; 87 / 64 / 25 | +8.6; six to seven in ten; nine in ten holds; 64 / 47 / 26 |
| 790 | Kabsch | indistinguishable on form; lead from top-5 | holds |
| 792 | cluster | 95 %; 57 % | recompute (D5) |
| 806, 824, 870 | — | qualitative | hold (cost ordering to confirm) |
| 817 | valid ∧ ≤ 2 | 80 vs 169 | 83 vs 183 |
| 834 | valid ∧ ≤ 2 | 36.6 / 34.0 / 18.2; 65.3 / 55.1 / 26.4; +2.6 (−4.2, +9.4); "eight points" | 49.2 / 43.6 / 18.8; 70.3 / 59.1 / 27.4; +5.6 (−1.7, +12.8); "twelve points" |
| 836 | interaction | did not separate | re-test |
| 843 | valid ∧ ≤ 2 (D14) | 0.1-0.7; 4.1-4.6; 3.5-4.4; 31.0 → 36.6; "does not clear Holm" | bands recompute (DiffDock smina 0.4-0.8 derived); 42.9 → 49.2; Holm per D14 family |
| 848-850 | cost | 2.4 / 49.4 / 137.2 s; 169 / 199 / 80; twentyfold; quarter; RQ3 ordering | recompute; 183 / 214 / 83; recompute; holds; W if ordering flips |
| 866 | triple | eight / one / nine; 3.0; 0.4-0.8 | nine / one / nine; 3.0 unchanged; recompute |
| 868 | stratum | 78; +6.7 on 225; +20.5 on 78; add convention sentence | 78; +8.4 (p 0.056); +19.2 (p 0.028); W |

### 3.3 `thesis_latex/body_appendix_short.tex`

| Line | Gate | Current | New |
|---|---|---|---|
| 77 | — | (v1 listed "Vinardo 37 vs 80; 42 vs 90" here; no such text exists in the short build) | drop the row |
| 165 | valid ∧ ≤ 2, stratum | 37.8 / 39.6 / 33.3 / 17.9 %; 134 vs 149 of 225, p 0.142; 33 vs 49 of 78, p 0.020; 73 / 81; −7.4 / −19.2; −6.3 / −21.0; 15 / 63; failure shares 53.6 / 65.4 / 55.8 / 70.3 | 51.1 / 49.3 / 43.6 / 26.9 %; 141 vs 160, p 0.056; 38 vs 53, p 0.028; 73 / 81 (rule kept); −9.1 / −17.8; −8.1 / −19.8; 16 / 62 or reconcile; 58.2 / 70.2 / 56.8 / 70.2; W crossover sentence and residue-rule wording |
| 167 | def. | box construction | W: add that every deposited copy lies inside the box |
| 173 | valid ∧ ≤ 2 | 111 (36.6 %) / 112 / 90; 31 = 15 + 16; 198 / 197 / 190; 199 | 149 (49.2 %) / recompute / recompute; recompute; 213 / recompute; 214 |
| 182 | valid ∧ ≤ 2 | 198; 2e-21..0.007; top-3 to top-22; "two of the seven survive"; 111 vs 110; 202 vs 199; 160 / 152 / 95; "one contrast over all twenty-eight pairs" | 213; 7.3e-20..5.2e-3; strictly 1-25 (W); six of seven (W); 149 vs 141, p 0.15 unresolved (W); 217 vs 214; 175 / 166 / recompute; recompute |
| 184 | valid ∧ ≤ 2 | 154 → 198; 185 → 198 (+17 / −4); 93 / 183 vs 94 / 185 | 170 → 213; 198 → 213 (+19 / −4); 129 / 196 vs 130 / 198 |
| 189, 195, 235, 241 | — | Figs 18-21 captions | `:195` W (+7 at rank-1, p 0.039); `:235` re-check; regenerate image44-47 |
| 218-227 | valid ∧ ≤ 2 | Table 8 and footnote | Section 1.5; W footnote |
| 230 | triple (pool) | 98 / 128 / 149 / 159 / 164; +71.6 / +35.9 / +5.9; +65.6 (+32.1..+137.2) | D10: keep on instance, add one sentence with nearest values (115 / 143 / 166 / 172 / 178 preview, first step +28, last +6; per-hour recompute) |
| 244-248 | — | five; three; 18.7; ~2; ~1.5 | recompute |
| 246 | near (no validity) | 95; 113; "exactly 113"; ranking-ceiling sentence | 124 / 124 / 131 and 144 vs 151: W, argument rewritten |
| 267-269 | triple / near / form | 2 / 46 / 55; 0 / 7 / 14; 0 / 5 / 18; 37; 19 / 64 / 83; 180-200; "no near-native complex at all"; 261 / 282 / 266 | 2 / 48 / 57; recompute (probed smina arms 10 / 9); recompute; recompute; 19 / 67 / 86; recompute (preview 179-199 on five of nine); W (P2Rank raw recovers one); unchanged |
| 276 | def. | "to the crystal ligand" | W |
| 487 | — | fifth reference-based check | optional clause on `load_all` |
| 545-553 | interaction | 4,253; 7.57 Å; 21.3 %; F1 0.565 / 0.529 / 0.471; 750; 0.0003; chain-keyed fingerprint | recompute; W |
| 664 | — | footnote file roles | W: documented roles |
| 666 | — | convention sentence | W: draft (Phase 5) |
| 722 | — | Fig 30 prose | W: name the 46 % (143 of 308) share |
| 784 | — | 1.56 Å; 65 % | unchanged (D13), W one clause |
| 871 | — | lead-in "validity-aware recovery (34.0 % and 34.3 %)" | W: name the near-native gate the tables use; 43.6 % and 44.2 % |
| 893-920 | near / valid | Wilson table | derived for five of seven variants; two recompute |
| 947-964 | near | Table 19 and footnote | AutoDock k = 1 43.2 → 49.8, 19 / 39, p 0.012, Holm 0.060 ns; DiffDock smina and gnina k = 5 lose their star; rest derived / recompute |
| 992-999 | near | bands | DiffDock smina +0.6 / +0.4 / +0.8 derived; rest recompute |
| 1026-1031 | near | cross-tool table and footnote | row 1: 49.8 vs 43.6, 67.3 vs 50.2, 70.0 vs 52.1, 71.3 vs 52.8; k = 1 77 / 58 uncorrected 0.12; rest recompute |
| 1056-1062 | accurate-invalid | Table 22 | AutoDock + gnina unchanged; DiffDock raw k = 1 / 5 change (47 → 65, 58 → 60); footnote unchanged |
| 1074 | ITT | none within 2 Å; 6.4 Å; 111 / 308, 198 / 308 vs 103, 167; 2.6; 10.2 → 10.1 | rank-1 none (6.3 Å); top-15 none but closest 2.1 Å; 149 / 308, 213 / 308 vs 132, 179; 5.5; 11.2 → 11.0 |
| 1172-1181 | near | PB decomposition near half; footnote | 437 / 98.6; 451 / 98.0; 2,118; 2,169 / 95.3; 2,171; 515 / 93.8; per-check recompute; W |
| 1188-1196 | near / valid | 5.4 %; 1,715; 6.2 %; 498; 31 of 498; 5 of 319; "at most four"; 2 / 4 / 3; 0.66-1.32 | strata 2,169 / 515 / 451; rates recompute; "at most five"; AutoDock 2 / 4 / 3 / 3 unchanged, DiffDock 5 / 4 / 3 / 3; 0.66-1.65 |
| 1203 | — | 25 Å cube argument | W: regime clause, "almost every other copy" |
| 1207 | valid ∧ ≤ 2 | 148 → 185 | 161 → 198 |
| 1247-1250 | Table 25 | +88 / −2 / +1 / +87; +63 / −1 / +2 / +64; +26 / −2 / +1 / +25 | +65 / −2 / +1 / +64; +45 / 0 / +2 / +47; +27 / −2 / +1 / +26 |
| 1258-1260 | valid ∧ ≤ 2 | +2.6 [−4.2, +9.4]; +13.2; +10.2 [+2.8, +17.5]; +9.9; 51 / 140 / 112 (37.0 %); MDD 9.8; power 0.12; 4,200; 8.2; 8.4; "within 10 points"; p 0.009 / 0.011 | +5.6 [−1.7, +12.8]; unchanged; +11.2 [+3.9, +18.4]; +10.2 [+3.0, +17.3]; 76 / 98 / 129 (42.6 %); 10.5; 0.32; 1,100; 8.8; 11.7; "within 12 points" (W); 0.0036 / 0.0075 |
| 1265 | — | Kabsch criterion sentence | W |
| 1303-1320 | valid ∧ ≤ thr | Table 26 RMSD rows | derived |
| 1321-1340 | valid ∧ Kabsch | Table 26 Kabsch rows | 47 of 90 cells move 0.3-0.7 pp |
| 1377-1387 | trim | Table 27 and footnote | all 117 cells; footnote seven numbers (Section 1.6); W |
| 1394-1399 | cluster | statistics; Fig 38 caption | recompute (D5); regenerate image8; W |
| 1407-1415 | interaction | Kendall τ; 15.3 vs 14.1; Fig 39 | recompute; regenerate image11; W |

---

## 4. Affected figures, tables, sidecars and harness assertions

**Figures (14 of 39):** Fig 2 (image3), 3 (image4), 4 (image5), 5 (image6), 6 (image7), 7 (image9), 15 (image24),
16 (image25), 18-21 (image44-47), 38 (image8), 39 (image11). The other 25 are unchanged: Fig 1 (validity-only script),
twelve Orai figures (Fig 24, image19, is Orai despite its "benchmark ligands" caption), the flow chart, two PyMOL renders
(image12, image14; orphan by design), Figs 28-37 dataset descriptors (Fig 30 panel A prints "Multi-instance ligands: 143
(46 % of set)" from the dataset notebook).

**Tables:** 1, 2, 3, 4, 6 (main); 8, 18, 19, 20 (near rows), 21, 22 (DiffDock raw rows), 24 (near half), 25, 26, 27
(appendix); new sensitivity table (D9, keyed `table_28`). Unchanged: 5, 7, 9-17, 23.

**Hub sidecars (one `--force` run):** `per_pose_metrics.csv` + manifest, `oracle_summary*.csv`, `top1_summary*.csv` (no
EquiBind row by design), `per_rank_metrics.csv`, `per_rank_ifp_recovery.csv`, `topn_within_thresholds*.csv` + `18_*_report.txt`,
`rank1_vs_topn.csv` (Gen A), `optimization_raw_vs_best.csv`, `optimization_benefit_by_rank.csv`, `pose_validity_cascade*.csv`,
`within2_validity_comparison.csv`, `form_fidelity_summary*.csv`, `oracle_selection_comparison.csv`, `filmstrip_stats__*.csv`,
`20d_*` tables and `caption.txt`, `twist_turn_summary.csv`, `rank_quality_*.csv`, `pb_test_waterfall*.csv`,
`comparison_vs_posebusters_paper.csv`, `posebusters_pose_comparison_stats.json`, the `pb_valid/` re-derivations.
**Non-hub sidecars:** `topk_recovery_validity_gnina_arm.csv` (bench_topk), `topN_crystal_cluster_matrix_stats.txt`,
`crystal_cluster_homogeneity_stats.*`, `cluster_quality_metrics_stats.txt`, `summary.json`, `ensembles_summary.csv`
(cluster report), `effort_summary.csv`, `effort_stats.json`, `effort_by_quality_stats.json` (effort), every PandaMap
`<stem>.txt` companion, `native_recovery_by_rank.csv`, `10b_contact_decomposition_by_rank.csv`,
`interaction_pose_basis_audit.csv`, `exh_*_stats.txt` (their `source:` line changes to the matched path), `exh_returns_report.txt`,
`exh_*.csv`, `validity_gate_cost.csv`, `metal_stratum*.csv`, `reference_convention_sensitivity.csv`. Unchanged:
`pbvalid_yield_*`.

**Harness:** changing blocks `table_1` (68 rows), `table_2`, `table_3`, `table_4`, `table_6`, `table_8`, `table_18`,
`table_19`, `table_20` (24 of 36), `table_21`, `table_22` (6 of 24), `table_24` (near half), `table_25`, `table_26`,
`table_27`, `prose.refiner_contrasts`, `prose.appendix_h9`, `prose.figure_39`, `prose.audit`, 14 `figure[imageN]` rows;
hard-zero blocks listed in 3.4; three new stages (registry 60 → 63) and four new prose checks (0.3).

---

## 5. Risks and traps

1. **Selection flips** (D12 margin one complex for DiffDock). Phase 4 gates Phase 5.
2. **Form loses complexes** in five printed Table 1 rows; "no complex loses" is near-nativeness only.
3. **Gen A / Gen B.** `rank1_vs_topn.csv` will again differ from Table 2; never diff them.
4. **`rmsd` versus `pb_rmsd`.** Verdicts agreed on 414 / 414 rank-1 poses; the `:113` count is two under both.
5. **Cache convention blindness** in BOTH paths (reuse: no signature check; non-reuse: silent recompute) until 1.9;
   the schema bump alone opens the window at every live manifest, and harness force mode does not pass `--force`.
6. **Figure gates are convention-blind.** md5 hand-check and `validate_regeneration.py` prove media == re-render of the
   current table; only the manifest check and `check_reference_convention` assert the convention. Figs 5, 6, 38 gain a
   byte gate only with the new Case (6.2).
7. **Kurzfassung** +1 and +11 characters; the cut must be named first.
8. **PandaMap crystal gate**: `overwrite: false` keeps the single-instance natives silently unless the two crystal CSVs
   are moved aside; `--overwrite` re-profiles every pose. Five `groupby` sites union all copies if any is missed.
9. **CASCADE_PINS abort** in `autodock_exhaustiveness_returns.py`; re-pin in 3.1a from the rebuilt cascade.
10. **Two cache schemas**: hub (`:13225`) and cluster report (`:3927`); bump both.
11. **Pre-existing defects** (0.3): `:113`, `:396`, `:414`, and the hydrogen-counting centroid loader behind `:476`,
    Table 3 and Figs 5 / 6 / 38. Fix before writing new numbers on top.
12. **Pre-specified confirmatory cell** (D10) and its paragraph's mixed conventions.
13. **EquiBind gnina_affinity ties** (333 poses in 133 groups) affect Table 27 depth-5 cells; take them from the sidecar.
14. **Never overwrite deliverables**; `REGENERATE.md` is generated; `CLEANUP_PLAN.md` and dated coverage files are
    records.
15. **Yaml transcription direction**: proven by the tex-presence check (7.1), not by the frozen-tree run, which cannot
    distinguish copy-from-outputs from transcribe-from-thesis. (Table 21 drifted because it had NO assertion, not
    because outputs were copied; the risk stands, the precedent does not.)
16. **Stale counts in the repo**: README 257 numbers (harness runs 1,364), THESIS_REPRODUCTION "225 of 225" and "58
    stages" (registry has 60).
17. **Comparability wording**: the cube kept "almost every other copy" outside; inertness of the published Vina anchor
    is inferred from the thesis's boxed arm, not measured.
18. **Two conventions in one paragraph**: App. C :230 / :182 / Figs 20-21 (D10), `:781` vs `:197` (depth basis),
    Table 19 vs `:843` (D14 family), the ladder's near-native (`:246`) vs validity-aware (`:182`) sentences.
19. **Metal stratum**: rule now registered; two complexes cross strata; 16 / 62 vs printed 15 / 63.
20. **Line anchors move** (the appendix shifted two lines during this very session); locate by quoted text.
21. **Probe previews vs rebuilt truth** (D16): EquiBind-raw k = 5 / 10 already disagrees with the print.
22. **Symmetry artefacts**: eight rank-1 flips in (1.5, 2.0] Å inspected in 0.5 before any number is written.

---

## 6. Effort estimate and critical path

| Phase | Effort |
|---|---|
| 0 Decide, freeze, fix pre-existing defects, metal stratum stage, flip inspection | 1 day |
| 1 Hub code, refusal semantics, regression | 1 day |
| 2 Downstream code, PandaMap gate, re-pin, new stages | 1.5 days |
| 3 Rebuild (two trees), re-pin, 19 stages, change signature | 0.5 day attention, 1-2 h compute (estimate) |
| 4 Selection and family re-check | 0.5 day |
| 5 Thesis text with per-step checks | 2.5-3 days |
| 6 Figures, media, byte gate for Figs 5 / 6 / 38 | 0.5 day |
| 7 Harness signatures A / B / C, yaml, notebook, docs | 1.5 days |
| 8 Build, sweep, QA | 1 day |
| **Total** | **10-10.5 working days** (v1: 8.5-9.5) |

Critical path: Phase 0.3 loader fix and cluster regeneration → hub code → rebuild → re-pin → PandaMap per-copy
fingerprints → change signature → text → harness signatures → QA. The PandaMap per-copy work remains the largest item;
its receptor question is now answered (all hosting chains present), so its fallback branch is gone.

---

## 7. What does NOT change

- **Orai1 panels** (Results Section 4.2 except the `:697` sentence, Discussion Section 5.5, Figs 8-14 and 22-27, Tables
  5, 7, 9-15, 23): no crystal ligand exists. (v1 wrote "Sections 4.2 and 5.2"; 5.2 is the validity section and DOES
  change at `:777-781`.)
- **X-ray crystal control**: busts `<ID>_ligand.sdf` with no reference and no RMSD. Its 307/308 ceiling covers one
  instance per entry; optionally bust the 211 alternates too, since they now count as correct.
- **PoseBusters validity battery and Fig 1**; Table 1 poses and valid columns; the 09f yield outputs; Table 24's
  all-poses half; Table 20's validity rows.
- **Cost totals and hardware currencies** (Table 10); only the per-qualifying-pose normalisation moves.
- **Dataset descriptors** (Figs 28-37, Tables 16-17); start-conformer RMSD (D13).
- **Top-5 AutoDock-DiffDock contrast**: identical (88 / 48, +13.2 pp, p_holm 7.6e-4).
- **Kabsch headline**: 74.3 / 74.6 % at 1 Å and 95.4 / 97.0 % at 2 Å; the form tie stays ns at every depth (the top-15
  Holm p becomes 0.899, so the literal 1.000 does move).
- **AutoDock* validity-gate cost** 2 / 4 / 3 / 3 and EquiBind* 2 / 3 / 3 / 3; DiffDock* changes at rank-1 and top-30.
- **Table 22** rows other than DiffDock raw k = 1 / 5 (confirm by harness).
- **Metal-adjacent stratum size** 78 (D4, registered rule); the five-complex exclusion; the 303 denominator;
  determinism evidence.
- **`pocket_comparison_report.py`** crystal outputs (D17, not in the short build); docking, PoseBusters runs, receptor and
  ligand preparation.

---

## 8. Audit ledger and changes from v1

**Audit run.** Workflow `wf_e3513a2f-47b`, 24 agents, 41 minutes, 0 agent errors. Per dimension (evidence / completeness):
hub-script 2 refuted, 17 corrected, 4 + 6 missing, 8 + 4 plan errors; downstream-scripts 0 / 4 / 5 + 11 / 6 + 11
(1 critical); main-body (completeness only) 0 / 0 / 5 / 5; appendix 0 / 7 + 1 / 3 + 3 / 9 + 8; harness-docs 1 / 14 + 1 /
7 + 9 / 6 + 8; figures-media 0 / 15 / 5 + 7 / 5 + 7; form-axis 0 / 9 / 2 + 4 / 5 + 6; methodology 1 / 8 / 4 + 2 / 7 + 7;
numbers-selection 3 / 4 / 6 + 8 / 7 + 8. Critic: 28 of 29 distinct critical-or-major findings accepted, one rejected
(the demand to probe the two gnina-guided EquiBind arms before Phase 5; the rebuild scores all 27 keys, so the residue
is to label the "179-199" preview as recompute, which v2 does). Six of the critic's eight gaps received gap-fill agents;
gaps 7 (gate column and Holm family) and 8 (single source of truth) were closed from the auditors' own findings as D14 and
D16.

**Refutations of v1 content, all accepted:** "26 of 90 Kabsch cells" → 47; the D9 stratum row was RMSD-only → validity
gate applied; Section 1.5 mixed the Table 19 basis with the `:843` basis; "2,337 by the script's rule" → the script's rule
gives 2,329 on the `pb_rmsd` axis; EquiBind near-site median 3.3 → 3.2; the 428-entry MDD 11.7 → 8.8; Uni-Dock2 is not a
Table 1 row; "Sections 4.2 and 5.2 unchanged" → 5.2 changes; table numbers (v1's Tables 5 / 6 / 7 / 8 are Tables 2 / 3 /
4 / 6); Figures 10 / 11 → 15 / 16; "three arms lose a Form complex" → five printed rows; "six adjacent ids" for the
centroid-vs-RMSD disagreement → 317 rows across 77 ids; Kabsch median shift 0.001 → 0.013 Å; "seven same-chain ids" → six
in the 308; "the cube excluded the copies" → almost all; the pooled reach under D2 is 99.3 %, under D5 99.7 %; "90 of 130
table_1 rows" → 68; "~850 rows" → bounded above by ≈ 1,074 with ≈ 290 hard-zero rows; "14 of 19 stages" → 16 of 19;
step 2.4 had no valid slot and its target aborts on `CASCADE_PINS`; the `:570` check is warn-only, not an assert; D8's
window opened at the schema bump and the plan flipped the default five phases before editing the stage command;
`figures.pairs` holds paths, not md5s; the Vinardo row at App. C :77 had no matching text; the `:664` draft's "no
reported pose" was false for 7Z1Q_NIO; the `:476` offset is a generator defect (hydrogens counted in pose centroids), not
a pose-basis detail; `:474` is not drift (sidecar double-rounding); `:781` vs `:197` is not a contradiction (different
depths).

**Pre-existing defects in the current thesis found by the audits, independent of the switch:** `:113` "three" should be
"two" and its `\ref` points at nothing; `:396` "p_holm ≤ 4e-10" should be "≤ 9e-9"; `:414` "2.9 Å" should be "3.0 Å";
`pose_cluster_crystal_pocket_report.py:133` counts hydrogens in pose centroids, moving every cluster-geometry number
(Table 3, Figs 5 / 6 / 38, `:474`, `:476`) by up to 0.58 Å per pose with four complexes crossing 4 Å; the cluster
sidecar prints 55.5 % for 168 / 303 by double rounding; App. C :165 prints a 15 / 63 cofactor split where the reproducing
rule gives 16 / 62; App. H :871 calls its tables "validity-aware" while they print the near-native gate;
`test_diffdock_sweep_gate.py` points at a notebook now under `obsolete/`; README and THESIS_REPRODUCTION quote stale
assertion counts. These are worth fixing whether or not the endpoint is switched.
