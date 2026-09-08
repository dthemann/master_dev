# REGENERATE.md — how every reported figure and table is produced

> **This file is generated. Do not edit it by hand.**
> It is written by `Scripts/Analysis/regenerate_guide.py`, which the last cell
> of `Thesis_Reproduction.ipynb` runs. The ordering, commands and determinism
> classes come from the same stage registry the notebook executes, so they
> cannot disagree with what actually runs. The figure index is rebuilt by
> checksum on every run, so a figure whose source moved shows as moved rather
> than as a stale path.
>
> The previous hand-written guide is kept as [`REGENERATE_NOTES.md`](REGENERATE_NOTES.md) for its
> per-command commentary and its trap notes, which are not derivable from the
> registry. It drifted: its figure numbers came from a build that has since
> been retired and its paths predated the matched-EquiBind migration. That
> drift is what this generator exists to prevent. See [`FINDINGS_2026-09-02.md`](FINDINGS_2026-09-02.md).

Generated 2026-09-08 15:25 from 66 registered stages.

Float numbers are the SHORT build's, read from the figure environments of
`body_main_short.tex` and `body_appendix_short.tex` in document order.

Scope note. The analysis outputs referenced below live in the working tree,
not in the git index. `posebusters_results/`, `pandamap_results/` and
`PoseBusters_Benchmark_Analysis/figures/` are untracked, so byte-identical
here means identical to the file on disk, not to one recoverable from history.

---

## 1. Environments and binaries

Every command below runs from `/home/manndo/master_dev`.

| Environment | Interpreter | Used for |
| --- | --- | --- |
| `vina` | `/home/manndo/anaconda3/envs/vina/bin/python` | screening, statistics and every figure; Python 3.12.12, RDKit 2025.09.1 |
| `diffdock` | `/home/manndo/anaconda3/envs/diffdock/bin/python` | DiffDock inference, torch 2.4.0 / CUDA 12.4; also carries smina |
| `equibind` | `/home/manndo/anaconda3/envs/equibind/bin/python` | EquiBind inference, torch 2.4.1 / CUDA 12.4; also carries smina |

| Binary | Path | Version |
| --- | --- | --- |
| AutoDock Vina | `/home/manndo/AutoDock-Vina/build/linux/release/vina` | v1.2.7-20-g93cdc3d-mod, the boron-patched build. /usr/bin/vina is stock 1.2.5 and must not be used. |
| gnina | `/home/manndo/docking_tools/gnina` | 1.3.2 |
| smina | `<env>/bin/smina` | 9 November 2017 build, on the Vina 1.1.2 scoring base |
| ADFRsuite | `/home/manndo/ADFRsuite-1.0/bin/` | prepare_ligand and prepare_receptor |

Two patches live outside the repository and are load-bearing. The PoseBusters
`energy_ratio.py` exception-handler fix, without which `bust()` aborts on any
ligand UFF cannot parameterise, and the DiffDock seed hook in the local
checkout, without which its sampling is irreproducible. The notebook asserts
both before running anything.

---

## 2. Canonical trees

The trap this table exists to prevent is that the directories whose names read
like scratch are the current ones. Generation order, oldest first, is
`*_PRE_FR0` and `*_PRE_EXH128` and `*_UFFON_backup_*`, then the plain name,
then `_matched`.

| Key | Directory |
| --- | --- |
| `benchmark_pb` | `posebusters_results/benchmark_matched_equibind/dock` |
| `benchmark_report` | `posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report` |
| `benchmark_clusters` | `posebusters_results/cluster_crystal_pocket_matched_equibind` |
| `benchmark_pandamap` | `pandamap_results/benchmark_matched_equibind` |
| `benchmark_effort_charged` | `posebusters_results/benchmark_matched_equibind/docking_effort_gnina_v2_charged` |
| `benchmark_effort_elapsed` | `posebusters_results/benchmark_matched_equibind/docking_effort_gnina_v2` |
| `exhaustiveness` | `posebusters_results/autodock_exhaustiveness_returns` |
| `orai_root` | `posebusters_results/_orai_matched_root` |
| `orai_control` | `posebusters_results/_orai_matched_root/orai_benchmark` |
| `orai_experimental` | `posebusters_results/_orai_matched_root/orai_jku` |
| `orai_yield` | `posebusters_results/_orai_matched_root/orai_pbvalid_yield_compare` |
| `orai_tm_share` | `posebusters_results/_orai_matched_root/orai_pbvalid_tm_share_compare` |
| `orai_pandamap_control` | `pandamap_results/orai_benchmark_matched` |
| `orai_pandamap_experimental` | `pandamap_results/orai_jku_matched` |
| `orai_pandamap_compare` | `pandamap_results/orai_interaction_compare_matched` |
| `dataset_figures` | `PoseBusters_Benchmark_Analysis/figures` |

Present on disk and deliberately NOT read:

| Superseded | Replaced by |
| --- | --- |
| `posebusters_results/benchmark_full_protein_vina_scoring` | `benchmark_pb` |
| `posebusters_results/benchmark` | `benchmark_pb` |
| `posebusters_results/cluster_crystal_pocket_full_protein` | `benchmark_clusters` |
| `pandamap_results/benchmark_full_protein_mgltools` | `benchmark_pandamap` |
| `posebusters_results/orai_benchmark` | `orai_control` |
| `posebusters_results/orai_jku` | `orai_experimental` |
| `posebusters_results/orai_pbvalid_yield_compare` | `orai_yield` |
| `posebusters_results/orai_pbvalid_tm_share_compare` | `orai_tm_share` |
| `pandamap_results/orai_interaction_compare` | `orai_pandamap_compare` |

---

## 3. Order

Declaration order below is a verified topological order of the dependency
graph, so running the sections in sequence can never read a file that a later
step writes. That failure mode is why this section is generated rather than
described.

**2. Benchmark docking** — 14 stages: `bench_prepared_inputs`, `bench_arm_status`, `bench_autodock_exh18_raw`, `bench_autodock_exh32_raw`, `bench_autodock_exh64_raw`, `bench_autodock_exh92_raw`, `bench_autodock_exh128_raw`, `bench_autodock_exh32_gnina`, `bench_autodock_exh64_gnina`, `bench_autodock_exh128_gnina`, `bench_diffdock`, `bench_pockets`, `bench_equibind`, `bench_equibind_minimize`

**3. Benchmark analysis** — 25 stages: `bench_posebusters`, `bench_posebusters_matched`, `bench_pose_comparison`, `bench_validity_report`, `bench_figure1`, `bench_filmstrip`, `bench_topk`, `bench_optimization_benefit`, `bench_endpoint_diagnostics`, `bench_metal_stratum`, `bench_selection_bias`, `bench_reference_identity`, `bench_alternate_copy_containment`, `bench_itt`, `bench_reference_convention`, `bench_clusters`, `bench_pandamap`, `bench_pandamap_report`, `bench_interaction_audit`, `bench_effort_charged`, `bench_effort_elapsed`, `bench_exhaustiveness_returns`, `bench_ligand_difficulty`, `bench_receptor_difficulty`, `bench_diffdock_rerank`

**4. Orai staging** — 2 stages: `orai_receptors`, `orai_ligands`

**5. Orai experimental** — 7 stages: `orai_exp_autodock`, `orai_exp_diffdock`, `orai_exp_equibind`, `orai_exp_posebusters`, `orai_exp_posebusters_fr0`, `orai_exp_pandamap`, `orai_exp_posebusters_equibind_minimize`

**6. Orai control** — 8 stages: `orai_ctl_autodock`, `orai_ctl_diffdock`, `orai_ctl_equibind`, `orai_ctl_posebusters`, `orai_ctl_posebusters_fr0`, `orai_ctl_posebusters_fr0_gnina`, `orai_ctl_posebusters_equibind_minimize`, `orai_ctl_pandamap`

**7. Orai cross-panel** — 9 stages: `orai_tm_exclusion`, `orai_exp_clusters`, `orai_ctl_clusters`, `orai_yield_compare`, `orai_tm_share_compare`, `orai_xtool_agreement`, `orai_pandamap_compare`, `orai_ligand_contrasts`, `orai_region_consensus`

**8. Dataset chapter** — 1 stages: `dataset_figures`

---

## 4. Stages

Determinism is a property of the recorded run, not a policy. `bitexact` means
a re-run must reproduce the stored bytes. `verify` means the reported run was
an unseeded single draw and can be checked but not redrawn. `never` means
re-running destroys the provenance of a reported number.

### Benchmark docking

#### `bench_prepared_inputs` — Prepared ligand and receptor PDBQT (ADFRsuite)

- **determinism** never · **cost** expensive
- **supports** App. B: prepare_ligand -A hydrogens, prepare_receptor -A hydrogens -U nphs_lps_waters_nonstdres
- **writes** `Dockings/vina_results_full_protein_vina_scoring_mgltools/_staging/ligands/pdbqt/*.pdbqt`

No command. This stage refuses to run under any mode.

Shared by every ladder rung. A newer Meeko changes atom ordering, the torsion-tree root and TORSDOF, so these are passed forward with --prepared-inputs-from rather than re-derived.

#### `bench_arm_status` — Exhaustiveness arm status (timing provenance)

- **determinism** never · **cost** cheap
- **supports** App. B hardware and timing basis: the 2,677.7 s gnina wall clock
- **writes** `Dockings/vina_results_full_protein_vina_scoring_mgltools_exh128/exhaustiveness_arm_status.json`

No command. This stage refuses to run under any mode.

A cached re-run overwrites the measured elapsed time with the cost of the cache lookup, which is near zero. No backup exists.

#### `bench_autodock_exh18_raw` — AutoDock Vina search, exhaustiveness 18

- **determinism** bitexact · **cost** expensive
- **supports** Table 8 and Figures 18-21, the appendix ladder
- **needs** `bench_prepared_inputs`
- **writes** `Dockings/vina_results_full_protein_vina_scoring_mgltools_exh18/*/mgl_tools/docking`

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Docking/run_autodock_exhaustiveness_arm.py \
    -c Scripts/Docking/configs_thesis/10_ladder_exh18.yaml \
    --ids-file Dockings/vina_results_full_protein_vina_scoring_mgltools_exh18/all308_ids.txt \
    --expect-exhaustiveness 18 \
    --prepared-inputs-from Dockings/vina_results_full_protein_vina_scoring_mgltools
```

Vina seed 42, energy range 6 kcal/mol, up to 30 modes, whole-receptor box.

#### `bench_autodock_exh32_raw` — AutoDock Vina search, exhaustiveness 32

- **determinism** bitexact · **cost** expensive
- **supports** Table 8 and Figures 18-21, the appendix ladder
- **needs** `bench_prepared_inputs`
- **writes** `Dockings/vina_results_full_protein_vina_scoring_mgltools/*/mgl_tools/docking`

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Docking/run_autodock_exhaustiveness_arm.py \
    -c Scripts/Docking/configs_thesis/11_ladder_exh32_raw.yaml \
    --ids-file Dockings/vina_results_full_protein_vina_scoring_mgltools/all308_ids.txt \
    --expect-exhaustiveness 32 \
    --prepared-inputs-from Dockings/vina_results_full_protein_vina_scoring_mgltools
```

Vina seed 42, energy range 6 kcal/mol, up to 30 modes, whole-receptor box.

#### `bench_autodock_exh64_raw` — AutoDock Vina search, exhaustiveness 64

- **determinism** bitexact · **cost** expensive
- **supports** Table 8 and Figures 18-21, the appendix ladder
- **needs** `bench_prepared_inputs`
- **writes** `Dockings/vina_results_full_protein_vina_scoring_mgltools_exh64/*/mgl_tools/docking`

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Docking/run_autodock_exhaustiveness_arm.py \
    -c Scripts/Docking/configs_thesis/13_ladder_exh64_raw.yaml \
    --ids-file Dockings/vina_results_full_protein_vina_scoring_mgltools_exh64/all308_ids.txt \
    --expect-exhaustiveness 64 \
    --prepared-inputs-from Dockings/vina_results_full_protein_vina_scoring_mgltools
```

Vina seed 42, energy range 6 kcal/mol, up to 30 modes, whole-receptor box.

#### `bench_autodock_exh92_raw` — AutoDock Vina search, exhaustiveness 92

- **determinism** bitexact · **cost** expensive
- **supports** Table 8 and Figures 18-21, the appendix ladder
- **needs** `bench_prepared_inputs`
- **writes** `Dockings/vina_results_full_protein_vina_scoring_mgltools_exh92/*/mgl_tools/docking`

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Docking/run_autodock_exhaustiveness_arm.py \
    -c Scripts/Docking/configs_thesis/15_ladder_exh92.yaml \
    --ids-file Dockings/vina_results_full_protein_vina_scoring_mgltools_exh92/all308_ids.txt \
    --expect-exhaustiveness 92 \
    --prepared-inputs-from Dockings/vina_results_full_protein_vina_scoring_mgltools
```

Vina seed 42, energy range 6 kcal/mol, up to 30 modes, whole-receptor box.

#### `bench_autodock_exh128_raw` — AutoDock Vina search, exhaustiveness 128

- **determinism** bitexact · **cost** expensive
- **supports** Table 1 and the Results chapter
- **needs** `bench_prepared_inputs`
- **writes** `Dockings/vina_results_full_protein_vina_scoring_mgltools_exh128/*/mgl_tools/docking`

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Docking/run_autodock_exhaustiveness_arm.py \
    -c Scripts/Docking/configs_thesis/01_benchmark_autodock_exh128_raw.yaml \
    --ids-file Dockings/vina_results_full_protein_vina_scoring_mgltools_exh128/all308_ids.txt \
    --expect-exhaustiveness 128 \
    --prepared-inputs-from Dockings/vina_results_full_protein_vina_scoring_mgltools
```

Vina seed 42, energy range 6 kcal/mol, up to 30 modes, whole-receptor box.

#### `bench_autodock_exh32_gnina` — gnina rescoring of the exhaustiveness 32 poses

- **determinism** bitexact · **cost** expensive
- **supports** Table 8, the appendix ladder
- **needs** `bench_autodock_exh32_raw`
- **writes** `Dockings/vina_results_full_protein_vina_scoring_mgltools/*/mgl_tools/docking/optimized_gnina`

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Docking/run_autodock_exhaustiveness_arm.py \
    -c Scripts/Docking/configs_thesis/12_ladder_exh32_gnina.yaml \
    --ids-file Dockings/vina_results_full_protein_vina_scoring_mgltools/all308_ids.txt
```

Rescoring reorders a rung's own poses; it does not search again. Ranked by CNNaffinity on the crossdock_default2018 ensemble, seed 42.

#### `bench_autodock_exh64_gnina` — gnina rescoring of the exhaustiveness 64 poses

- **determinism** bitexact · **cost** expensive
- **supports** Table 8, the appendix ladder
- **needs** `bench_autodock_exh64_raw`
- **writes** `Dockings/vina_results_full_protein_vina_scoring_mgltools_exh64/*/mgl_tools/docking/optimized_gnina`

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Docking/run_autodock_exhaustiveness_arm.py \
    -c Scripts/Docking/configs_thesis/14_ladder_exh64_gnina.yaml \
    --ids-file Dockings/vina_results_full_protein_vina_scoring_mgltools_exh64/all308_ids.txt
```

Rescoring reorders a rung's own poses; it does not search again. Ranked by CNNaffinity on the crossdock_default2018 ensemble, seed 42.

#### `bench_autodock_exh128_gnina` — gnina rescoring of the exhaustiveness 128 poses

- **determinism** bitexact · **cost** expensive
- **supports** the DOMINANT arm, Tables 1/2/6 and Figures 1-3, 5-7
- **needs** `bench_autodock_exh128_raw`
- **writes** `Dockings/vina_results_full_protein_vina_scoring_mgltools_exh128/*/mgl_tools/docking/optimized_gnina`

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Docking/run_autodock_exhaustiveness_arm.py \
    -c Scripts/Docking/configs_thesis/02_benchmark_autodock_exh128_gnina.yaml \
    --ids-file Dockings/vina_results_full_protein_vina_scoring_mgltools_exh128/all308_ids.txt
```

Rescoring reorders a rung's own poses; it does not search again. Ranked by CNNaffinity on the crossdock_default2018 ensemble, seed 42.

#### `bench_diffdock` — DiffDock-L, 30 samples, 20 denoising steps

- **determinism** verify · **cost** expensive
- **supports** Tables 1, 2, 6 and Figures 1-3, 5-7, 15, 16
- **writes** `Dockings/Benchmark_DiffDock/*/docking_log.csv`

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Docking/run_diffdock.py \
    -c Scripts/Docking/configs_thesis/03_benchmark_diffdock.yaml
```

UNSEEDED. None of the stored logs carries the seed marker, so re-running draws a new sample and cannot reproduce the reported poses. The driver refuses to run this stage; restore the tree from backup instead.

#### `bench_pockets` — fpocket and P2Rank site detection

- **determinism** bitexact · **cost** moderate
- **supports** App. B pocket-guided EquiBind variants; also the cluster-to-pocket analysis
- **writes** `pocket_results/fpocket_results/*`, `pocket_results/p2rank_results/*`

No committed command; the outputs are staged by the arm that produced them.

Prerequisite for the guided EquiBind arms and for the crystal-pocket cluster report. In the old notebook this sat under an EquiBind heading, which hid the dependency.

#### `bench_equibind` — EquiBind, 30 conformers, unguided and pocket-guided

- **determinism** bitexact · **cost** expensive
- **supports** Tables 1, 2, 6; the guided-arm paragraph of App. B
- **needs** `bench_pockets`
- **writes** `Dockings/Benchmark_Equibind/*/pipeline_summary.json`

No committed command; the outputs are staged by the arm that produced them.

Config 04 replaces the deprecated equibind_docking_config.yaml. uff_minimize is false, which App. B records as the setting for every reported run; the old notebook rewrote it to true.

#### `bench_equibind_minimize` — EquiBind re-refinement with --minimize (matched arm)

- **determinism** bitexact · **cost** expensive
- **supports** examiner item M5: every EquiBind refinement had run --local_only while AutoDock and DiffDock ran --minimize
- **needs** `bench_equibind`
- **writes** `Dockings/Benchmark_Equibind_minimize/*`

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Docking/rerun_equibind_refine.py
```

Re-refines from the committed __refRAW.sdf, which is the exact input the original refinement received, so the flag is the only variable. This is what makes benchmark_matched_equibind the canonical tree.

### Benchmark analysis

#### `bench_posebusters` — PoseBusters validity screen, all benchmark arms

- **determinism** bitexact · **cost** expensive
- **supports** App. B validity screening; PoseBusters 0.6.3
- **needs** `bench_autodock_exh128_gnina`, `bench_diffdock`, `bench_equibind`
- **writes** `posebusters_results/benchmark_full_protein_vina_scoring/dock/posebusters_filtered_results.csv`

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Docking/Posebusters/run_posebusters.py \
    -c Scripts/Docking/configs_thesis/05_benchmark_posebusters.yaml
```

#### `bench_posebusters_matched` — PoseBusters, matched-EquiBind arm (CANONICAL)

- **determinism** bitexact · **cost** expensive
- **supports** every benchmark number in the Results chapter
- **needs** `bench_equibind_minimize`
- **writes** `posebusters_results/benchmark_matched_equibind/dock/posebusters_filtered_results.csv`

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Docking/Posebusters/run_posebusters.py \
    -c Scripts/Docking/configs_thesis/06_benchmark_posebusters_matched_equibind.yaml
```

#### `bench_pose_comparison` — Pose-comparison hub: per-pose metrics and recovery oracle

- **determinism** bitexact · **cost** moderate
- **supports** Table 1; the source table for nine downstream stages
- **needs** `bench_posebusters_matched`
- **writes** `posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report/per_pose_metrics.csv`, `posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report/oracle_summary.csv`

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Analysis/posebusters_pose_comparison.py \
    --pb-csv posebusters_results/benchmark_matched_equibind/dock/posebusters_filtered_results.csv \
    --out-dir posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report \
    --top-n 15 \
    --diffdock-variant all \
    --reference-convention nearest \
    --collapse-plots-only \
    --collapse-diffdock-variant diffdock_smina \
    --collapse-autodock-variant autodock_mgltools_exh128_gnina \
    --exclude-preset meeko \
    --workers 24
```

PROMOTED 2026-09-08: --reference-convention nearest is the primary endpoint (every near-native, form, centroid, contact and PLIF metric is taken against the deposited copy of the ligand nearest to the pose; record 0 of <ID>_ligands.sdf is the single-instance sensitivity arm, kept in the *_ref_instance columns). The cache manifest records the convention and the hub refuses to rebuild a cache of another convention or schema without --force. Program record: Scripts/Analysis/nearest_copy_program/PROGRAM_LOG.md.
--top-n 15 and --diffdock-variant all shaped the cached table. A re-render with --reuse-cache omits them, so dropping the cache after changing the input silently rebuilds at --top-n 5 on the raw diffdock key.
CORRECTED 2026-09-02: the pin was --collapse-autodock-variant autodock_gnina, which the meeko preset cancels, because autodock_gnina is IN that preset. The two flags together dropped the whole AutoDock family and Figures 2 and 3 rendered with two tools instead of three. No error was raised. The dominant arm must be named in full.
CORRECTED 2026-09-03: the cross-tool baseline in _TOPK_TEST_PAIRS was raw Vina, while Appendix C declares the family is computed on the gnina-rescored arm and Table 21 prints that arm's rates. The script carried raw Vina in every commit of its history, so it tested 31.4 against DiffDock's 33.7 instead of the printed 37.3, and reported ns at three of the four depths where the thesis reports significance. Table 21 had no assertion behind it, which is why the drift survived. check_table_21 now pins it.

#### `bench_validity_report` — Validity report per tool and complex

- **determinism** bitexact · **cost** cheap
- **supports** Table 1 validity block
- **needs** `bench_pose_comparison`
- **writes** `posebusters_results/benchmark_matched_equibind/dock/validity_report_mgltools`

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Analysis/posebusters_validity_report.py \
    --csv posebusters_results/benchmark_matched_equibind/dock/posebusters_filtered_results.csv \
    --out-dir posebusters_results/benchmark_matched_equibind/dock/validity_report_mgltools \
    --per-pose-metrics posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report/per_pose_metrics.csv \
    --collapse-plots-only \
    --exclude-preset meeko \
    --exclude-methods autodock_mgltools,autodock_mgltools_gnina,autodock_mgltools_exh18,autodock_mgltools_exh64,autodock_mgltools_exh64_gnina,autodock_mgltools_exh92
```

#### `bench_figure1` — Figure 1: post-hoc optimisation and PoseBusters validity

- **determinism** bitexact · **cost** cheap
- **supports** Figure 1 (media/media/image2.png)
- **needs** `bench_posebusters_matched`
- **writes** `posebusters_results/benchmark_matched_equibind/dock/validity_report_mgltools/00_figure2_validity_yield.png`

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Analysis/figure2_validity_yield.py \
    --csv posebusters_results/benchmark_matched_equibind/dock/posebusters_filtered_results.csv \
    --out posebusters_results/benchmark_matched_equibind/dock/validity_report_mgltools/00_figure2_validity_yield.png \
    --exclude-preset meeko \
    --exclude-methods autodock_mgltools,autodock_mgltools_gnina,autodock_mgltools_exh18,autodock_mgltools_exh64,autodock_mgltools_exh64_gnina,autodock_mgltools_exh92
```

The double exclusion is load-bearing. --exclude-preset meeko alone leaves five ladder rungs to be folded onto a bare 'autodock' key.

#### `bench_filmstrip` — Figure 4: form against placement across ranking depth

- **determinism** bitexact · **cost** cheap
- **supports** Figure 4 (media/media/image5.png); Table 27
- **needs** `bench_pose_comparison`
- **writes** `posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report/20d_form_vs_placement_by_family__depth_filmstrip_pbvalid__rank1_top5_top15__thesis.png`

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Analysis/filmstrip_rank1_top5_top15.py \
    --report-dir posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report \
    --exclude-preset meeko
```

--report-dir was added 2026-09-02. The module constant still names the pre-matched tree, and the shipped figure comes from the matched one.

#### `bench_topk` — Top-k recovery sidecar for the gnina arm

- **determinism** bitexact · **cost** cheap
- **supports** Tables 2, 18, 25
- **needs** `bench_pose_comparison`
- **writes** `posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report/topk_recovery_validity_gnina_arm.csv`

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Analysis/topk_recovery_gnina_arm.py \
    --report-dir posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report \
    --autodock-variant autodock_mgltools_exh128_gnina \
    --diffdock-variant diffdock_smina \
    --equibind-variant equibind_unguided_gnina \
    --exclude-preset meeko
```

#### `bench_optimization_benefit` — Optimisation benefit, paired categorical inference

- **determinism** bitexact · **cost** cheap
- **supports** Tables 19, 20
- **needs** `bench_pose_comparison`
- **writes** `posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report/optimization_benefit_by_rank.csv`

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Analysis/optimization_benefit_stats.py \
    --report-dir posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report \
    --exclude-preset meeko
```

#### `bench_endpoint_diagnostics` — Endpoint diagnostics: validity gate cost and bounded claim

- **determinism** bitexact · **cost** cheap
- **supports** Table 22 and the Results prose numbers
- **needs** `bench_pose_comparison`
- **writes** `posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report/validity_gate_cost.csv`

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Analysis/thesis_endpoint_diagnostics.py \
    --metrics posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report/per_pose_metrics.csv \
    --csv posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report \
    --exclude-preset meeko
```

#### `bench_metal_stratum` — Metal-adjacent stratum: membership, cofactor/ion split and stratum contrasts

- **determinism** bitexact · **cost** cheap
- **supports** App. C metal-adjacent stratum (78 of 303) and body_main_short.tex:88/:868
- **needs** `bench_pose_comparison`
- **writes** `posebusters_results/metal_stratum/metal_stratum_summary.json`, `posebusters_results/metal_stratum/metal_stratum.csv`, `posebusters_results/metal_stratum/metal_stratum_contrasts.csv`

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Analysis/metal_stratum.py \
    --benchmark-dir "Data/PoseBuster Benchmark Set" \
    --ids-file posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report/analysed_cohort_ids.txt \
    --per-pose-csv posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report/per_pose_metrics.csv \
    --out-dir posebusters_results/metal_stratum \
    --arms autodock_mgltools_exh128_gnina diffdock_smina equibind_unguided_gnina \
    --cutoffs 4 5 6 \
    --overwrite
```

ADDED 2026-09-08 (plan D4 / Phase 0.4). No generator existed before this date: the 78-complex list behind App. C :165 and the 73 / 81 sensitivity counts were never registered and could not be asserted. The rule that reproduces the printed 73 / 78 / 81 at 4 / 5 / 6 A is the RESIDUE-level one: minimum heavy-atom distance from the reference ligand instance (record 0 of <ID>_ligands.sdf == <ID>_ligand.sdf) to ANY atom of a residue containing at least one metal element (alkali, alkaline-earth, transition, post-transition), measured on the PoseBusters-shipped <ID>_protein.pdb, not on the docked receptor. The atom-level rule (distance to the metal atom itself) gives 65 / 77 / 80 and does NOT reproduce the print; the single complex separating 77 from 78 is 7TSF_H4B (haem ring at 2.78 A, Fe at 9.25 A). Membership is a complex-level property of the reference instance and does not follow the scored copy; it differs across deposited copies for exactly two complexes (7JHQ_VAJ in by the reference only, 7TB0_UD1 in by copy 1 only), which the summary JSON lists. The rule yields 16 cofactor / 62 free-ion, where :165 prints 15 / 63. The reference-convention contrasts reproduce :165 verbatim (37.8 / 39.6 / 33.3 / 17.9 % rank-1; 134 vs 149 of 225 p 0.142; 33 vs 49 of 78 p 0.020; Fisher 0.256). Pass --nearest-csv <per-pose CSV with rmsd_nearest_copy> to add the nearest-copy convention block. --overwrite only permits replacing this script's own three output files; any other content in --out-dir is refused.

#### `bench_selection_bias` — Variant-selection winner's curse: per-family paired estimate

- **determinism** bitexact · **cost** cheap
- **supports** App. C 'One caveat applies to the selected arm itself' and the Limitations winner's-curse range
- **needs** `bench_pose_comparison`
- **writes** `posebusters_results/selection_bias/selection_bias.csv`, `posebusters_results/selection_bias/selection_bias_summary.json`

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Analysis/nearest_copy_program/selection_bias.py \
    --per-pose-csv posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report/per_pose_metrics.csv \
    --out-dir posebusters_results/selection_bias \
    --overwrite
```

ADDED 2026-09-08 (nearest-copy program, decision 4). Before this date the winner's-curse figures in App. C and the Limitations had no generator. The estimator is the standard c_k * sigma_d / sqrt(2) with c_k the expected maximum of k standard normals and sigma_d the PAIRED (McNemar) standard error of the recovery difference between the family winner and its runner-up on the same complexes, evaluated on the selection endpoint (top-15, PB-valid AND rmsd <= 2 AND bestfit_rmsd <= 1). Families: the two exh128 AutoDock arms (k = 2), the three DiffDock arms and the three unguided EquiBind arms (k = 3), plus the eight-rung AutoDock exhaustiveness ladder (k = 8; the base rung is keyed autodock_mgltools[_gnina]). The double gate and rank-1 are printed for comparison. Output is labelled by the table's reference_convention.

#### `bench_reference_identity` — Reference-dependent PoseBusters identity checks re-run against every deposited ligand copy

- **determinism** bitexact · **cost** cheap
- **supports** Results footnote on the cost of a co-crystallised reference (-1 / -2 / -2 DiffDock complexes)
- **needs** `bench_pose_comparison`
- **writes** `posebusters_results/reference_identity/identity_cost_by_depth.csv`, `posebusters_results/reference_identity/identity_cost_summary.json`

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Analysis/nearest_copy_program/reference_identity_rerun.py \
    --per-pose-csv posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report/per_pose_metrics.csv \
    --out-dir posebusters_results/reference_identity \
    --arms autodock_mgltools_exh128_gnina autodock_mgltools_exh128 diffdock_smina diffdock_gnina equibind_unguided_gnina \
    --overwrite
```

ADDED 2026-09-08 (nearest-copy program, decision 5). The footnote's -1 / -2 / -2 was measured post hoc on 2026-08-22 with no registered generator. This stage recomputes the four InChI-based identity checks (formula, connections, tetrahedral, double-bond; inchi_options 'w' as in the redock config) per pose against the single instance, the multi-copy file loaded with load_all and the nearest copy of each pose, then counts the complexes the identity requirement removes from the PB-valid AND rmsd <= 2 recovery set per arm and depth. On 2026-09-08 all three references give identical per-pose outcomes (0.95 % of poses fail; no complex has copies with different InChI) and the cost is -1 / -2 / -2 (7ZXV_45D; +7XQZ_FPF) for both DiffDock arms under both conventions, 0 for AutoDock and EquiBind. About 1 minute on 8 workers.

#### `bench_alternate_copy_containment` — Alternate deposited ligand copies inside a 25 A PoseBusters-style cube on the reference instance

- **determinism** bitexact · **cost** cheap
- **supports** App. C box-volume paragraph ('only N of the 211 alternate copies have a centroid inside it')
- **writes** `posebusters_results/reference_convention/alternate_copies_in_25A_cube.json`

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Analysis/nearest_copy_program/alternate_copy_box_containment.py \
    --benchmark-dir "Data/PoseBuster Benchmark Set" \
    --ids-file "Data/PoseBuster Benchmark Set/posebusters_pdb_ccd_ids.txt" \
    --out posebusters_results/reference_convention/alternate_copies_in_25A_cube.json
```

ADDED 2026-09-08 after the appendix numbers review: the printed '8 of the 211' was a plan preview with no generator and was the 12.5 A SPHERE count; the stated rule (alternate copy's heavy-atom centroid inside the 25 A cube centred on the reference heavy-atom centroid) gives 11 of 211 over the 308 entries (16 by any heavy atom). 143 multi-copy entries, 211 alternates.

#### `bench_itt` — Intention-to-treat check of the five dropped complexes

- **determinism** bitexact · **cost** cheap
- **supports** App. C cohort paragraph (closest AutoDock pose of the dropped complexes)
- **needs** `bench_pose_comparison`
- **writes** `posebusters_results/itt/itt_dropped_complexes.csv`, `posebusters_results/itt/itt_summary.json`

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Analysis/nearest_copy_program/itt_dropped_complexes.py \
    --out-dir posebusters_results/itt \
    --per-pose-csv posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report/per_pose_metrics.csv
```

ADDED 2026-09-08 (nearest-copy program, plan 2.6). Scores the raw AutoDock poses of the five complexes dropped from the analysed cohort against every deposited copy, so the cohort paragraph can say how close the excluded complexes came instead of asserting it.

#### `bench_reference_convention` — Reference-convention sensitivity table (nearest copy against single instance)

- **determinism** bitexact · **cost** cheap
- **supports** Table 23 (App. C) and the Methods sensitivity sentence; 143 of 308 and 138 of 303 multi-copy counts
- **needs** `bench_pose_comparison`
- **writes** `posebusters_results/reference_convention/reference_convention_sensitivity.csv`, `posebusters_results/reference_convention/reference_convention_sensitivity.tex`, `posebusters_results/reference_convention/reference_convention_margins.csv`, `posebusters_results/reference_convention/reference_convention_summary.json`

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Analysis/nearest_copy_program/reference_convention_table.py \
    --per-pose-csv posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report/per_pose_metrics.csv \
    --out-dir posebusters_results/reference_convention
```

ADDED 2026-09-08 (nearest-copy program, plan 2.7). Reads the rmsd / rmsd_ref_instance pair of the per-pose table and prints both conventions at every depth for the three selected arms, with the AutoDock-DiffDock margin and its Newcombe interval under each.

#### `bench_clusters` — Pose clusters against the crystal pocket

- **determinism** bitexact · **cost** moderate
- **supports** Table 3; Figures 5, 6, 38
- **needs** `bench_pose_comparison`, `bench_pockets`
- **writes** `posebusters_results/cluster_crystal_pocket_matched_equibind/autodock_mgltools_exh128_gnina__diffdock_smina_allposes/topN_crystal_cluster_matrix.png`, `posebusters_results/cluster_crystal_pocket_matched_equibind/autodock_mgltools_exh128_gnina__diffdock_smina_allposes/cluster_quality_metrics.png`

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Analysis/pose_cluster_crystal_pocket_report.py \
    --ids-file "Data/PoseBuster Benchmark Set/posebusters_pdb_ccd_ids.txt" \
    --per-pose-csv posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report/per_pose_metrics.csv \
    --out-dir posebusters_results/cluster_crystal_pocket_matched_equibind/autodock_mgltools_exh128_gnina__diffdock_smina_allposes \
    --autodock-variant autodock_mgltools_exh128_gnina \
    --diffdock-variant diffdock_smina \
    --equibind-variant equibind_unguided_gnina \
    --site-cluster threshold \
    --pocket-radius 8.0 \
    --crystal-copies any \
    --rank-by consensus \
    --stability-boot 25 \
    --workers 30 \
    --stats
```

The out-dir tag reads 'allposes' because --pb-valid-only is OFF for the shipped figures: the cluster geometry is measured over every pose.

#### `bench_pandamap` — PandaMap interaction fingerprints

- **determinism** bitexact · **cost** moderate
- **supports** Table 4; Figures 7, 39
- **needs** `bench_posebusters_matched`
- **writes** `pandamap_results/benchmark_matched_equibind/pandamap_interactions.csv`

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Analysis/run_pandamap.py \
    -c Scripts/Docking/configs_thesis/07_benchmark_pandamap.yaml
```

#### `bench_pandamap_report` — Interaction-difference report

- **determinism** bitexact · **cost** cheap
- **supports** Table 4; Figure 7
- **needs** `bench_pandamap`, `bench_pose_comparison`
- **writes** `pandamap_results/benchmark_matched_equibind/report/05_fingerprint_similarity_top5.png`, `pandamap_results/benchmark_matched_equibind/report/04c_native_recovery_by_rank.png`

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Analysis/pandamap_interaction_report.py \
    -c Scripts/Docking/configs_thesis/07_benchmark_pandamap.yaml \
    --in-dir pandamap_results/benchmark_matched_equibind \
    --per-pose-metrics posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report/per_pose_metrics.csv \
    --exclude-methods unidock2,autodock_gnina
```

#### `bench_interaction_audit` — Interaction pose-basis audit

- **determinism** bitexact · **cost** cheap
- **supports** App. B "Interaction Fingerprints" quotes this output verbatim
- **needs** `bench_pandamap_report`
- **writes** `pandamap_results/benchmark_matched_equibind/report/interaction_pose_basis_audit.txt`, `pandamap_results/benchmark_matched_equibind/report/interaction_pose_basis_audit.csv`

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Analysis/interaction_pose_basis_audit.py \
    --config Scripts/Docking/configs_thesis/07_benchmark_pandamap.yaml \
    --in-dir pandamap_results/benchmark_matched_equibind \
    --metrics posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report/per_pose_metrics.csv \
    --ids-file posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report/analysed_cohort_ids.txt \
    --exclude-methods unidock2,autodock_gnina
```

--ids-file is REQUIRED and is the whole difficulty of this stage. The PandaMap arm profiles 306 complexes while the analysed cohort is 303, so an unrestricted run leaves nine EquiBind poses in three excluded complexes without an RMSD and the join guard aborts. The guard is right: an unguarded gap would bin those poses into a nan band and contaminate the standardisation reference. The cohort file is derived from per_pose_metrics rather than written by hand.

#### `bench_effort_charged` — Cost per qualifying pose, charged basis

- **determinism** bitexact · **cost** cheap
- **supports** Table 6; Figure 15
- **needs** `bench_pose_comparison`
- **writes** `posebusters_results/benchmark_matched_equibind/docking_effort_gnina_v2_charged/effort_summary.csv`, `posebusters_results/benchmark_matched_equibind/docking_effort_gnina_v2_charged/panels/effort_by_quality_near2.png`

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Analysis/docking_effort_comparison.py \
    --dataset benchmark \
    --autodock-dir Dockings/vina_results_full_protein_vina_scoring_mgltools_exh128 \
    --autodock-prep mgl_tools \
    --autodock-refine gnina \
    --autodock-gnina-gpu \
    --autodock-method autodock_mgltools_exh128_gnina \
    --autodock-optimizer-workers 16 \
    --equibind-dir Dockings/Benchmark_Equibind_cputimed \
    --unidock2-dir  \
    --unidock-dir  \
    --per-pose-csv posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report/per_pose_metrics.csv \
    --out-dir posebusters_results/benchmark_matched_equibind/docking_effort_gnina_v2_charged \
    --basis charged \
    --cpu-threads 32
```

The charged basis divides cpu_core_s/32 + gpu_s uniformly. The default elapsed basis divides wall_s, which is not one quantity across arms and is sensitive to --autodock-optimizer-workers; charged is not.

#### `bench_effort_elapsed` — Hardware resource per near-native valid pose

- **determinism** bitexact · **cost** cheap
- **supports** Table 10; Figure 16
- **needs** `bench_pose_comparison`
- **writes** `posebusters_results/benchmark_matched_equibind/docking_effort_gnina_v2/panels/resource_per_near_native_valid_pose.png`

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Analysis/docking_effort_comparison.py \
    --dataset benchmark \
    --autodock-dir Dockings/vina_results_full_protein_vina_scoring_mgltools_exh128 \
    --autodock-prep mgl_tools \
    --autodock-refine gnina \
    --autodock-gnina-gpu \
    --autodock-method autodock_mgltools_exh128_gnina \
    --autodock-optimizer-workers 16 \
    --equibind-dir Dockings/Benchmark_Equibind_cputimed \
    --unidock2-dir  \
    --unidock-dir  \
    --per-pose-csv posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report/per_pose_metrics.csv \
    --out-dir posebusters_results/benchmark_matched_equibind/docking_effort_gnina_v2
```

#### `bench_exhaustiveness_returns` — Exhaustiveness ladder: marginal return per hour

- **determinism** bitexact · **cost** moderate
- **supports** Table 8; Figures 18-21
- **needs** `bench_pose_comparison`, `bench_autodock_exh18_raw`, `bench_autodock_exh92_raw`
- **writes** `posebusters_results/autodock_exhaustiveness_returns/figures/exh_raw_vs_rescored.png`, `posebusters_results/autodock_exhaustiveness_returns/figures/exh_depth_sweep.png`, `posebusters_results/autodock_exhaustiveness_returns/figures/exh_yield_vs_cost.png`, `posebusters_results/autodock_exhaustiveness_returns/figures/exh_marginal_return.png`

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Analysis/autodock_exhaustiveness_returns.py \
    --per-pose-csv posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report/per_pose_metrics.csv \
    --out-dir posebusters_results/autodock_exhaustiveness_returns \
    --n-boot 5000 \
    --verify-input-parity \
    --cascade-pins-from posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report/pose_validity_cascade.csv
```

5,000 paired complex-level bootstrap resamples at a fixed seed; numerator and denominator recomputed on the same resample.

#### `bench_ligand_difficulty` — Docking difficulty by ligand attributes

- **determinism** bitexact · **cost** cheap
- **supports** Discussion; the appendix difficulty analysis
- **needs** `bench_pose_comparison`
- **writes** `PoseBusters_Benchmark_Analysis/ligand_difficulty`

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Analysis/ligand_docking_difficulty.py \
    --per-pose-metrics posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report/per_pose_metrics.csv \
    --features PoseBusters_Benchmark_Analysis/ligand_protein_features.csv \
    --out-dir PoseBusters_Benchmark_Analysis/ligand_difficulty \
    --diffdock-variant smina \
    --exclude-preset meeko \
    --exclude-methods autodock_mgltools,autodock_mgltools_gnina,autodock_mgltools_exh18,autodock_mgltools_exh64,autodock_mgltools_exh64_gnina,autodock_mgltools_exh92
```

#### `bench_receptor_difficulty` — Docking difficulty by receptor attributes

- **determinism** bitexact · **cost** cheap
- **supports** Discussion; the appendix difficulty analysis
- **needs** `bench_pose_comparison`
- **writes** `PoseBusters_Benchmark_Analysis/receptor_difficulty`

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Analysis/receptor_docking_difficulty.py \
    --per-pose-metrics posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report/per_pose_metrics.csv \
    --features PoseBusters_Benchmark_Analysis/ligand_protein_features.csv \
    --out-dir PoseBusters_Benchmark_Analysis/receptor_difficulty \
    --diffdock-variant smina \
    --exclude-preset meeko \
    --exclude-methods autodock_mgltools,autodock_mgltools_gnina,autodock_mgltools_exh18,autodock_mgltools_exh64,autodock_mgltools_exh64_gnina,autodock_mgltools_exh92
```

#### `bench_diffdock_rerank` — DiffDock re-ranking counterfactual

- **determinism** bitexact · **cost** cheap
- **supports** the appendix note that smina does not re-rank DiffDock
- **needs** `bench_pose_comparison`
- **writes** `PoseBusters_Benchmark_Analysis/smina_rerank`, `PoseBusters_Benchmark_Analysis/gnina_rerank`

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Analysis/diffdock_gnina_rerank_analysis.py \
    --tool smina \
    --per-pose-metrics posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report/per_pose_metrics.csv \
    --out-dir PoseBusters_Benchmark_Analysis/smina_rerank
```

The --per-pose-metrics default still names the pre-whole-protein table, so it must be passed explicitly.

### Orai staging

#### `orai_receptors` — Four Orai1 receptor frames (COPY ONLY)

- **determinism** never · **cost** moderate
- **supports** App. D Orai1 receptor model; App. B Orai1 panels
- **writes** `Data/Receptors/Orai1WT-START-Fr0.pdb`, `Data/Receptors/Orai1WT-MDSnap-Fr300.pdb`, `Data/Receptors/Orai1WT-MDSnap-Fr400.pdb`, `Data/Receptors/Orai1WT-MDSnap-Fr499.pdb`

No command. This stage refuses to run under any mode.

Fr0 passed through an unguarded OpenMM minimisation and is not bit-reproducible (0.394 A between two preps). run_orai_mgltools_arm.py sha256-pins it to 6b3ab996... Never re-prepare; copy.

#### `orai_ligands` — Prepared Orai ligand PDBQT (COPY ONLY)

- **determinism** never · **cost** cheap
- **supports** App. B: the three modulators were docked neutral in every reported arm
- **writes** `Data/Ligands/JKU/pdbqt/*.pdbqt`, `Dockings/Orai_Benchmark_MGLTools_exh128/_staging/ligands/pdbqt/*.pdbqt`

No command. This stage refuses to run under any mode.

Same Meeko-version hazard as the benchmark ligands.

### Orai experimental

#### `orai_exp_autodock` — AutoDock Vina exh128 + gnina, 30 modes

- **determinism** bitexact · **cost** expensive
- **supports** Table 5; Figures 8-14
- **needs** `orai_receptors`, `orai_ligands`
- **writes** `Dockings/Orai_JKU_MGLTools_exh128/mgl_tools`

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Docking/run_orai_mgltools_arm.py \
    -c Scripts/Docking/configs_thesis/20_orai_jku_autodock_exh128_gnina.yaml
```

run_orai_mgltools_arm.py is the only legal driver. run_autodock.py main() re-prepares Fr0 and forces Meeko ligands.

#### `orai_exp_diffdock` — DiffDock-L, 30 samples (SEEDED)

- **determinism** bitexact · **cost** expensive
- **supports** Table 5; Figures 8-14
- **needs** `orai_receptors`, `orai_ligands`
- **writes** `Dockings/diffdock_results`

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Docking/run_diffdock.py \
    -c Scripts/Docking/configs_thesis/21_orai_jku_diffdock.yaml
```

24 stored run logs carry the seed marker, so this arm alone can be redrawn.

#### `orai_exp_equibind` — EquiBind unguided, 30 conformers

- **determinism** bitexact · **cost** expensive
- **supports** Table 5; Figures 8-14
- **needs** `orai_receptors`, `orai_ligands`
- **writes** `Dockings/equibind_results_uffoff`

No committed command; the outputs are staged by the arm that produced them.

Config 22. uff_minimize false and output_dir repointed to the uffoff tree.

#### `orai_exp_posebusters` — PoseBusters validity screen

- **determinism** bitexact · **cost** moderate
- **supports** Table 5; Figures 8-10
- **needs** `orai_exp_autodock`, `orai_exp_diffdock`, `orai_exp_equibind`
- **writes** `posebusters_results/_orai_matched_root/orai_jku/dock/posebusters_filtered_results.csv`

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Docking/Posebusters/run_posebusters.py \
    -c Scripts/Docking/configs_thesis/23_orai_jku_posebusters.yaml
```

variant_filter is load-bearing: raw ADFRsuite PDBQT carries no REMARK SMILES, so an unfiltered run falls back to Open Babel and fails silently on about 63 per cent of poses.

#### `orai_exp_posebusters_fr0` — Fr0 receptor correction, AUTODOCK ROWS ONLY

- **determinism** bitexact · **cost** moderate
- **supports** the corrected Fr0 validity figures
- **needs** `orai_exp_posebusters`
- **writes** `posebusters_results/orai_jku_mgltools_exh128_fr0corrected`

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Docking/Posebusters/run_posebusters.py \
    -c Scripts/Docking/configs_thesis/24_orai_jku_posebusters_fr0corrected.yaml
```

AutoDock only. Extending it to DiffDock or EquiBind drops EquiBind validity from 58.2 to 15.8 per cent, which is the mismatch, not the fix.

#### `orai_exp_pandamap` — PandaMap interaction fingerprints

- **determinism** bitexact · **cost** moderate
- **supports** Table 12; Figures 13, 14, 25, 26
- **needs** `orai_exp_posebusters`
- **writes** `pandamap_results/orai_jku_matched`

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Analysis/run_pandamap.py \
    -c Scripts/Docking/configs_thesis/25_orai_jku_pandamap.yaml
```

PandaMap 4.1.0 cannot see chlorine or bromine, and a renumbered receptor yields zero fingerprints.

### Orai control

#### `orai_ctl_autodock` — AutoDock Vina exh128 + gnina, 10 modes

- **determinism** bitexact · **cost** expensive
- **supports** Table 5; Figures 8-14
- **needs** `orai_receptors`, `orai_ligands`
- **writes** `Dockings/Orai_Benchmark_MGLTools_exh128/mgl_tools`

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Docking/run_orai_mgltools_arm.py \
    -c Scripts/Docking/configs_thesis/30_orai_benchmark_autodock_exh128_gnina.yaml
```

#### `orai_ctl_diffdock` — DiffDock-L, 10 samples

- **determinism** verify · **cost** expensive
- **supports** Table 5; Figures 8-14
- **needs** `orai_receptors`, `orai_ligands`
- **writes** `Dockings/Orai_Benchmark_DiffDock`

No committed command; the outputs are staged by the arm that produced them.

UNSEEDED, like the benchmark run. App. B names this as the one remaining input difference between the two Orai panels.

#### `orai_ctl_equibind` — EquiBind unguided, 10 conformers

- **determinism** bitexact · **cost** expensive
- **supports** Table 5; Figures 8-14
- **needs** `orai_receptors`, `orai_ligands`
- **writes** `Dockings/Orai_Benchmark_Equibind`

No committed command; the outputs are staged by the arm that produced them.

Config 32. This is the arm the old notebook already read strictly, with a cfg_req() helper that raises on a missing key.

#### `orai_ctl_posebusters` — PoseBusters validity screen

- **determinism** bitexact · **cost** expensive
- **supports** Table 5; Figures 8-10
- **needs** `orai_ctl_autodock`, `orai_ctl_diffdock`, `orai_ctl_equibind`
- **writes** `posebusters_results/_orai_matched_root/orai_benchmark/dock/posebusters_filtered_results.csv`

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Docking/Posebusters/run_posebusters.py \
    -c Scripts/Docking/configs_thesis/33_orai_benchmark_posebusters.yaml
```

#### `orai_ctl_posebusters_fr0` — Fr0 correction, AUTODOCK ROWS ONLY

- **determinism** bitexact · **cost** moderate
- **supports** Fr0 validity 73.9 -> 97.7 per cent; the AutoDock control yield
- **needs** `orai_ctl_posebusters`
- **writes** `posebusters_results/orai_benchmark_fr0corrected`

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Docking/Posebusters/run_posebusters.py \
    -c Scripts/Docking/configs_thesis/34_orai_benchmark_posebusters_fr0corrected.yaml
```

#### `orai_ctl_posebusters_fr0_gnina` — Fr0 correction (gnina variant), AUTODOCK ROWS ONLY

- **determinism** bitexact · **cost** moderate
- **supports** Fr0 validity 73.9 -> 97.7 per cent; the AutoDock control yield
- **needs** `orai_ctl_posebusters`
- **writes** `posebusters_results/orai_benchmark_gnina_fr0corrected`

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Docking/Posebusters/run_posebusters.py \
    -c Scripts/Docking/configs_thesis/35_orai_benchmark_posebusters_gnina_fr0corrected.yaml
```

#### `orai_ctl_posebusters_equibind_minimize` — PoseBusters, matched-EquiBind arm (control panel)

- **determinism** bitexact · **cost** moderate
- **supports** the EquiBind rows of Table 5 and of every Orai figure
- **needs** `orai_ctl_posebusters`
- **writes** `posebusters_results/orai_benchmark_equibind_minimize`

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Docking/Posebusters/run_posebusters.py \
    -c Scripts/Docking/configs_thesis/37_orai_benchmark_posebusters_equibind_minimize.yaml
```

ADDED 2026-09-02. This arm had no stage, so the canonical _orai_matched_root trees could not be rebuilt from the declared stages. Configs 23 and 33 screen the PRE-matched EquiBind trees; this screens the --minimize re-refinement whose poses the canonical tables actually cite, 24,640 references on the control panel and 1,440 on the experimental one. Both screens are needed and are not alternatives.

### Orai experimental

#### `orai_exp_posebusters_equibind_minimize` — PoseBusters, matched-EquiBind arm (experimental panel)

- **determinism** bitexact · **cost** moderate
- **supports** the EquiBind rows of Table 5 and of every Orai figure
- **needs** `orai_exp_posebusters`
- **writes** `posebusters_results/orai_jku_equibind_minimize`

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Docking/Posebusters/run_posebusters.py \
    -c Scripts/Docking/configs_thesis/38_orai_jku_posebusters_equibind_minimize.yaml
```

ADDED 2026-09-02. This arm had no stage, so the canonical _orai_matched_root trees could not be rebuilt from the declared stages. Configs 23 and 33 screen the PRE-matched EquiBind trees; this screens the --minimize re-refinement whose poses the canonical tables actually cite, 24,640 references on the control panel and 1,440 on the experimental one. Both screens are needed and are not alternatives.

### Orai control

#### `orai_ctl_pandamap` — PandaMap interaction fingerprints

- **determinism** bitexact · **cost** expensive
- **supports** Table 12; Figures 13, 14, 25, 26
- **needs** `orai_ctl_posebusters`
- **writes** `pandamap_results/orai_benchmark_matched`

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Analysis/run_pandamap.py \
    -c Scripts/Docking/configs_thesis/36_orai_benchmark_pandamap.yaml
```

### Orai cross-panel

#### `orai_tm_exclusion` — Transmembrane pore exclusion, BOTH panels in one run

- **determinism** bitexact · **cost** moderate
- **supports** Table 5; Figures 9, 10
- **needs** `orai_exp_posebusters_fr0`, `orai_ctl_posebusters_fr0`
- **writes** `posebusters_results/_orai_matched_root/orai_jku/transmembrane_filter/tm_pose_classification.csv`, `posebusters_results/_orai_matched_root/orai_benchmark/transmembrane_filter/tm_pose_classification.csv`

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Analysis/orai_transmembrane_exclusion.py \
    --autodock-variant gnina \
    --diffdock-variant smina \
    --equibind-variant gnina \
    --out-root posebusters_results/_orai_matched_root
```

The pore slab runs R91 to E106. In the old notebook a JKU-only copy of this step sat in the experimental section and wrote the same output, so which one won depended on execution order.

#### `orai_exp_clusters` — Reference-free pose clustering (experimental)

- **determinism** bitexact · **cost** moderate
- **supports** Table 3 analogue for Orai; Figures 11, 12, 23, 24
- **needs** `orai_tm_exclusion`
- **writes** `posebusters_results/_orai_matched_root/orai_jku/pose_clusters/per_pair.csv`, `posebusters_results/_orai_matched_root/orai_jku/pose_clusters/cluster_quality_per_tool.csv`

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Analysis/orai_pose_cluster_report.py \
    --per-pose-csv posebusters_results/_orai_matched_root/orai_jku/dock/posebusters_filtered_results.no_tm.csv \
    --out-dir posebusters_results/_orai_matched_root/orai_jku/pose_clusters \
    --autodock-variant gnina \
    --diffdock-variant smina \
    --site-cluster threshold \
    --pocket-radius 8.0 \
    --cluster-quality \
    --top-n-poses 10
```

Clusters are built on the pooled three-tool cloud, so changing one tool's pose set moves the other two tools' cluster values as well.

#### `orai_ctl_clusters` — Reference-free pose clustering (control)

- **determinism** bitexact · **cost** moderate
- **supports** Table 3 analogue for Orai; Figures 11, 12, 23, 24
- **needs** `orai_tm_exclusion`
- **writes** `posebusters_results/_orai_matched_root/orai_benchmark/pose_clusters/per_pair.csv`, `posebusters_results/_orai_matched_root/orai_benchmark/pose_clusters/cluster_quality_per_tool.csv`

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Analysis/orai_pose_cluster_report.py \
    --per-pose-csv posebusters_results/_orai_matched_root/orai_benchmark/dock/posebusters_filtered_results.no_tm.csv \
    --out-dir posebusters_results/_orai_matched_root/orai_benchmark/pose_clusters \
    --autodock-variant gnina \
    --diffdock-variant smina \
    --site-cluster threshold \
    --pocket-radius 8.0 \
    --cluster-quality \
    --top-n-poses 10
```

Clusters are built on the pooled three-tool cloud, so changing one tool's pose set moves the other two tools' cluster values as well.

#### `orai_yield_compare` — Figure 8: PoseBusters-valid yield, control against experimental

- **determinism** bitexact · **cost** cheap
- **supports** Figure 8 (media/media/image13.png)
- **needs** `orai_tm_exclusion`
- **writes** `posebusters_results/_orai_matched_root/orai_pbvalid_yield_compare/09f_pbvalid_yield_dominant_orai_benchmark_vs_experimental.png`

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Analysis/orai_pbvalid_yield_compare.py \
    --results-root posebusters_results/_orai_matched_root
```

Extracted 2026-09-02 from an inline notebook cell. Until then this was the only float in the Results chapter with no committed generator. Verified byte-identical to the shipped asset.

#### `orai_tm_share_compare` — Validity and transmembrane share, control against experimental

- **determinism** bitexact · **cost** cheap
- **supports** Table 5; Figures 9, 10
- **needs** `orai_tm_exclusion`
- **writes** `posebusters_results/_orai_matched_root/orai_pbvalid_tm_share_compare/pbvalid_tm_share_pooled.csv`, `posebusters_results/_orai_matched_root/orai_pbvalid_tm_share_compare/pbvalid_tm_share_per_unit.csv`, `posebusters_results/_orai_matched_root/orai_pbvalid_tm_share_compare/fig_pbvalid_outside_tm_whisker.png`, `posebusters_results/_orai_matched_root/orai_pbvalid_tm_share_compare/fig_tm_loss_relative_compare.png`

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Analysis/orai_pbvalid_tm_share_compare.py \
    --autodock-variant gnina \
    --diffdock-variant smina \
    --equibind-variant gnina \
    --results-root posebusters_results/_orai_matched_root \
    --out-dir posebusters_results/_orai_matched_root/orai_pbvalid_tm_share_compare \
    --top-n-poses 10
```

#### `orai_xtool_agreement` — Cross-tool agreement and reference-free cluster quality

- **determinism** bitexact · **cost** cheap
- **supports** Figures 11, 12, 23
- **needs** `orai_exp_clusters`, `orai_ctl_clusters`
- **writes** `posebusters_results/_orai_matched_root/orai_jku/pose_clusters/panels/orai_cross_tool_agreement_compare.png`, `posebusters_results/_orai_matched_root/orai_jku/pose_clusters/panels/orai_tool_agreement_compare.png`

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Analysis/orai_cross_tool_agreement_compare.py \
    --exp-csv posebusters_results/_orai_matched_root/orai_jku/pose_clusters/per_pair.csv \
    --exp-label "Exp. Ligands" \
    --bench-csv posebusters_results/_orai_matched_root/orai_benchmark/pose_clusters/per_pair.csv \
    --bench-label "Benchmark ligands x Orai" \
    --out-dir posebusters_results/_orai_matched_root/orai_jku/pose_clusters/panels
```

This script takes no --autodock-variant. Its arm is fixed by whichever run wrote the two per_pair.csv inputs, which is why they are named explicitly.

#### `orai_pandamap_compare` — Interaction compare: totals, type profile, hot-spots, overlap

- **determinism** bitexact · **cost** cheap
- **supports** Table 12; Figures 13, 14, 25, 26
- **needs** `orai_exp_pandamap`, `orai_ctl_pandamap`
- **writes** `pandamap_results/orai_interaction_compare_matched/fig_type_profile_compare.png`, `pandamap_results/orai_interaction_compare_matched/fig_fingerprint_overlap_compare.png`

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Analysis/orai_pandamap_interaction_compare.py \
    --exp-dir pandamap_results/orai_jku_matched \
    --bench-dir pandamap_results/orai_benchmark_matched \
    --out-dir pandamap_results/orai_interaction_compare_matched \
    --autodock-variant gnina \
    --diffdock-variant smina \
    --equibind-variant gnina \
    --top-n-poses 10 \
    --exp-posebusters-csv posebusters_results/_orai_matched_root/orai_jku/dock/posebusters_filtered_results.csv \
    --bench-posebusters-csv posebusters_results/_orai_matched_root/orai_benchmark/dock/posebusters_filtered_results.csv
```

The two --*-posebusters-csv paths are REQUIRED and their defaults are wrong for this tree. EquiBind's PandaMap pose_rank is a 999 placeholder, so the top-ten cap has to read the ranks from the PoseBusters CSV. With the default paths, which name the plain trees, the lookup fails and all 324 EquiBind pose rows are dropped without an error. The four figures then render without an EquiBind series and differ from the published ones.

#### `orai_ligand_contrasts` — Ligand-level sensitivity re-test (LAST in the chain)

- **determinism** bitexact · **cost** cheap
- **supports** the Methods commitment that every Orai cross-panel contrast is re-tested at ligand level in a sensitivity analysis reported with its result
- **needs** `orai_tm_share_compare`, `orai_exp_clusters`, `orai_ctl_clusters`, `orai_pandamap_compare`
- **writes** `posebusters_results/_orai_matched_root/orai_pbvalid_tm_share_compare/orai_ligand_level_contrasts.csv`, `posebusters_results/_orai_matched_root/orai_pbvalid_tm_share_compare/orai_ligand_level_contrasts.txt`

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Analysis/orai_ligand_level_contrasts.py \
    --autodock-variant gnina \
    --diffdock-variant smina \
    --equibind-variant gnina \
    --per-unit-csv posebusters_results/_orai_matched_root/orai_pbvalid_tm_share_compare/pbvalid_tm_share_per_unit.csv \
    --exp-quality-csv posebusters_results/_orai_matched_root/orai_jku/pose_clusters/cluster_quality_per_tool.csv \
    --bench-quality-csv posebusters_results/_orai_matched_root/orai_benchmark/pose_clusters/cluster_quality_per_tool.csv \
    --exp-pair-csv posebusters_results/_orai_matched_root/orai_jku/pose_clusters/per_pair.csv \
    --bench-pair-csv posebusters_results/_orai_matched_root/orai_benchmark/pose_clusters/per_pair.csv \
    --pandamap-totals-csv pandamap_results/orai_interaction_compare_matched/pandamap_pose_totals.csv \
    --pandamap-exp-dir pandamap_results/orai_jku_matched \
    --pandamap-bench-dir pandamap_results/orai_benchmark_matched \
    --exp-posebusters-csv posebusters_results/_orai_matched_root/orai_jku/dock/posebusters_filtered_results.csv \
    --bench-posebusters-csv posebusters_results/_orai_matched_root/orai_benchmark/dock/posebusters_filtered_results.csv \
    --out-dir posebusters_results/_orai_matched_root/orai_pbvalid_tm_share_compare
```

The effective number of independent experimental observations is three ligands, so this is what keeps the frame-ligand tests from being read as more powered than they are.
CORRECTED 2026-09-03: the stage passed only --per-unit-csv and the three variant flags. The script takes nine input paths and the other eight default to the superseded plain trees, so the registered command silently re-derived the sensitivity analysis from one generation back. It did not error and the output landed under the canonical out-dir, which is what made it invisible. The divergence is material rather than cosmetic: EquiBind's cluster-quality contrasts move from 6 frame-ligand units to 4, Compactness/EquiBind flips from -0.09 ns to +0.39, and adjusted p-values move across the whole file. This is the same defect already documented for orai_pandamap_compare, which was fixed there and left standing here. Verified 2026-09-03: with the paths above the script reproduces the shipped sidecar exactly, byte for byte apart from its own generated-on timestamp line.

#### `orai_region_consensus` — Consensus binding-region analysis

- **determinism** bitexact · **cost** moderate
- **supports** the Orai1 interpretation section
- **needs** `orai_tm_exclusion`
- **writes** `posebusters_results/_orai_matched_root/orai_benchmark/region_consensus`

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Analysis/orai_ligand_region_consensus.py \
    --pb-csv posebusters_results/_orai_matched_root/orai_benchmark/dock/posebusters_filtered_results.no_tm.csv \
    --out-dir posebusters_results/_orai_matched_root/orai_benchmark/region_consensus
```

Label-shuffle permutation null, 2,000 iterations at a fixed seed.

### Dataset chapter

#### `dataset_figures` — Chemical-space figures and descriptor tables

- **determinism** bitexact · **cost** moderate
- **supports** Tables 16, 17; Figures 28-37
- **writes** `PoseBusters_Benchmark_Analysis/figures/01_ligand_univariate_distributions.png`, `PoseBusters_Benchmark_Analysis/figures/09_pca_scree_loadings.png`, `PoseBusters_Benchmark_Analysis/summary_statistics.csv`

```bash
/home/manndo/anaconda3/envs/vina/bin/python \
    -m jupyter nbconvert \
    --to notebook \
    --execute \
    --inplace \
    --ExecutePreprocessor.timeout=1800 PoseBusters_DataSet_Analysis.ipynb
```

Descriptors are computed on a hydrogen-suppressed molecule; two of the sixteen are basis-dependent. Figure-to-file mapping is not sequential: image39 is the scree plot, image41 the receptor profile.

---

## 5. Figure index

Rebuilt by checksum on every generation. The source column is the file whose
md5 equals the shipped asset's, preferring the canonical tree where several
byte-identical copies exist.

| Fig | Asset | Subject | Source |
| --- | --- | --- | --- |
| 1 | `image2` | Post-hoc optimisation and PoseBusters validity. Significance… | `posebusters_results/benchmark_matched_equibind/dock/validity_report_mgltools/00_figure2_validity_yield.png` |
| 2 | `image3` | Best-of-top-N accuracy against the as-placed RMSD to the… | `posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report/18_topn_within_thresholds_pbvalid_depths.png` |
| 3 | `image4` | Best-of-top-N accuracy against the best-fit (Kabsch) RMSD | `posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report/18_topn_within_thresholds_kabsch_pbvalid_depths.png` |
| 4 | `image5` | Form versus in-place RMSD across ranking depth | `posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report/20d_form_vs_placement_by_family__depth_filmstrip_pbvalid__rank1_top5_top15__thesis.png` |
| 5 | `image6` | Top-N crystal-cluster co-recovery | `posebusters_results/cluster_crystal_pocket_matched_equibind/autodock_mgltools_exh128_gnina__diffdock_smina_allposes/topN_crystal_cluster_matrix.png` (+1 identical) |
| 6 | `image7` | Composition and geometry of the crystal-closest cluster. Rates… | `posebusters_results/cluster_crystal_pocket_matched_equibind/autodock_mgltools_exh128_gnina__diffdock_smina_allposes/crystal_cluster_homogeneity.png` (+1 identical) |
| 7 | `image9` | Mean Jaccard similarity of interaction fingerprints… | `pandamap_results/benchmark_matched_equibind/report/05_fingerprint_similarity_top5.png` |
| 8 | `image13` | PoseBusters-valid yield of produced poses on Orai | `posebusters_results/_orai_matched_root/orai_pbvalid_yield_compare/09f_pbvalid_yield_dominant_orai_benchmark_vs_experimental.png` |
| 9 | `image43` | Usable pose yield per Orai frame-ligand unit | `posebusters_results/_orai_matched_root/orai_pbvalid_tm_share_compare/fig_pbvalid_outside_tm_whisker.png` |
| 10 | `image15` | Transmembrane loss of PoseBusters-valid poses | `posebusters_results/_orai_matched_root/orai_pbvalid_tm_share_compare/fig_tm_loss_relative_compare.png` |
| 11 | `image16` | Distribution of inter-tool consensus-site distances on Orai… | `posebusters_results/_orai_matched_root/orai_jku/pose_clusters/panels/orai_cross_tool_agreement_compare.png` |
| 12 | `image17` | Cross-tool agreement on Orai | `posebusters_results/_orai_matched_root/orai_jku/pose_clusters/panels/orai_tool_agreement_compare.png` |
| 13 | `image21` | Interaction-type profile on Orai | `pandamap_results/orai_interaction_compare_matched/fig_type_profile_compare.png` |
| 14 | `image23` | Contact-fingerprint overlap between the two ligand sets | `pandamap_results/orai_interaction_compare_matched/fig_fingerprint_overlap_compare.png` |
| 15 | `image24` | Charged seconds per qualifying pose | `posebusters_results/benchmark_matched_equibind/docking_effort_gnina_v2_charged/panels/effort_by_quality_near2.png` |
| 16 | `image25` | Hardware resource per qualifying pose. The two currencies are… | `posebusters_results/benchmark_matched_equibind/docking_effort_gnina_v2_charged/panels/resource_per_near_native_valid_pose.png` (+1 identical) |
| 17 | `image1` | Physics-based docking flow chart | `docking_workflow.png` |
| 18 | `image44` | Raw Vina ranking against gnina re-ranking at each rung | `posebusters_results/autodock_exhaustiveness_returns/figures/exh_raw_vs_rescored.png` |
| 19 | `image45` | What each rung of the ladder buys against ranking depth | `posebusters_results/autodock_exhaustiveness_returns/figures/exh_depth_sweep.png` |
| 20 | `image46` | Docking success against the compute it cost | `posebusters_results/autodock_exhaustiveness_returns/figures/exh_yield_vs_cost.png` |
| 21 | `image47` | Marginal return per additional hour of search | `posebusters_results/autodock_exhaustiveness_returns/figures/exh_marginal_return.png` |
| 22 | `image14` | Orai Fr300 with docked ligands and transmembrane exclusion | **no source on disk** |
| 23 | `image18` | Reference-free cluster quality per tool on Orai | `posebusters_results/_orai_matched_root/orai_jku/pose_clusters/panels/orai_cluster_quality_compare.png` |
| 24 | `image19` | Effect of validity and placement filtering on cluster quality | `posebusters_results/_orai_matched_root/orai_benchmark/pose_clusters/orai_cluster_quality_filtering.png` |
| 25 | `image20` | Total contacts per pose on Orai | `pandamap_results/orai_interaction_compare_matched/fig_total_interactions_compare.png` |
| 26 | `image22` | Orai residue hot-spots | `pandamap_results/orai_interaction_compare_matched/fig_residue_hotspots_compare.png` |
| 27 | `image12` | The docking receptor | **no source on disk** |
| 28 | `image33` | Univariate descriptor distributions | `PoseBusters_Benchmark_Analysis/figures/01_ligand_univariate_distributions.png` |
| 29 | `image34` | Drug-likeness of the benchmark ligands | `PoseBusters_Benchmark_Analysis/figures/02_druglikeness_ro5_qed.png` |
| 30 | `image35` | Elemental and charge composition | `PoseBusters_Benchmark_Analysis/figures/03_elemental_charge_composition.png` |
| 31 | `image36` | Descriptor correlation matrix | `PoseBusters_Benchmark_Analysis/figures/04_ligand_correlation_heatmap.png` |
| 32 | `image37` | Two-Dimensional Attribute Mappings Coloured by QED | `PoseBusters_Benchmark_Analysis/figures/05_2d_attribute_mappings.png` |
| 33 | `image38` | Joint distributions of six core descriptors | `PoseBusters_Benchmark_Analysis/figures/06_pairplot_core_descriptors.png` |
| 34 | `image39` | Principal Component Analysis of the Ligand Descriptors | `PoseBusters_Benchmark_Analysis/figures/09_pca_scree_loadings.png` |
| 35 | `image40` | PCA biplot and score plot | `PoseBusters_Benchmark_Analysis/figures/10_pca_biplot.png` |
| 36 | `image41` | Receptor profile | `PoseBusters_Benchmark_Analysis/figures/07_receptor_profile.png` |
| 37 | `image42` | Conformer-generation difficulty | `PoseBusters_Benchmark_Analysis/figures/08_startconf_rmsd.png` |
| 38 | `image8` | Cluster quality on the calibration benchmark | `posebusters_results/cluster_crystal_pocket_matched_equibind/autodock_mgltools_exh128_gnina__diffdock_smina_allposes/cluster_quality_metrics.png` (+1 identical) |
| 39 | `image11` | Typed contacts versus the crystal by pose rank | `pandamap_results/benchmark_matched_equibind/report/10b_contact_decomposition_by_rank.png` |

37 of 39 figures resolve to a file on disk.

---

## 6. Table index

Tables whose printed values are recomputed and asserted on every run by
`thesis_assertions.py`. A failure there is reported, never silently fixed.

| Table | Source | Location in the thesis |
| --- | --- | --- |
| 1 | `posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report/per_pose_metrics.csv` | body_main_short.tex:261, tab:results-pose-production |
| 10 | `posebusters_results/benchmark_matched_equibind/docking_effort_gnina_v2_charged/effort_summary.csv` | body_appendix_short.tex, tab:appendix-cost-per-pose |
| 11 | `` | body_appendix_short.tex, tab:appendix-orai-protocols and the two paragraphs after it |
| 12 | `pandamap_results/orai_interaction_compare_matched/orai_pandamap_interaction_compare_stats.txt` | body_appendix_short.tex, tab:appendix-orai-interaction |
| 14 | `` | body_appendix_short.tex, tab:appendix-frame-geometry |
| 16 | `PoseBusters_Benchmark_Analysis/summary_statistics.csv` | body_appendix_short.tex, tab:appendix-ligand-distribution |
| 17 | `` | body_appendix_short.tex, tab:appendix-orai-descriptors |
| 18 | `posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report/topk_recovery_validity.csv` | body_appendix_short.tex, tab:appendix-top-k-recovery |
| 19 | `posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report/topk_recovery_stats.csv` | body_appendix_short.tex, tab:appendix-refinement-mcnemar |
| 2 | `posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report/topk_recovery_validity_gnina_arm.csv` | body_main_short.tex:336, tab:results-depth-recovery |
| 20 | `posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report/optimization_benefit_by_rank.csv` | body_appendix_short.tex, tab:appendix-refinement-bands |
| 21 | `posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report/topk_recovery_stats.csv` | body_appendix_short.tex, tab:appendix-cross-tool-mcnemar |
| 22 | `posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report/topk_recovery_validity.csv` | body_appendix_short.tex, tab:appendix-accurate-invalid |
| 24 | `` | body_appendix_short.tex, tab:appendix-pb-decomposition-full |
| 25 | `posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report/per_pose_metrics.csv` | body_appendix_short.tex, tab:results-depth-gain |
| 26 | `posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report/per_pose_metrics.csv` | body_appendix_short.tex, tab:results-near-native-form |
| 27 | `posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report/20d_form_vs_placement_by_family__depth_filmstrip_pbvalid__rank1_top5_top15__stats_per_tool_depth.csv` | body_appendix_short.tex, tab:results-placement-form |
| 3 | `posebusters_results/cluster_crystal_pocket_matched_equibind/autodock_mgltools_exh128_gnina__diffdock_smina_allposes/topN_crystal_cluster_matrix_stats.txt` | body_main_short.tex, tab:results-cluster-recovery |
| 4 | `pandamap_results/benchmark_matched_equibind/report/native_recovery_by_rank.csv` | body_main_short.tex, tab:results-native-recovery-rank |
| 5 | `posebusters_results/_orai_matched_root/orai_pbvalid_tm_share_compare/pbvalid_tm_share_pooled.csv` | body_main_short.tex, tab:results-orai-yield |
| 6 | `posebusters_results/benchmark_matched_equibind/docking_effort_gnina_v2_charged/effort_summary.csv` | body_main_short.tex, tab:results-cost-per-qualifying-pose |
| 8 | `posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report/per_pose_metrics.csv` | body_appendix_short.tex, tab:appendix-exhaustiveness-ladder |
| 9 | `` | body_appendix_short.tex, tab:appendix-protocols |

Tables 7, 9, 11, 13, 15 and 23 are manual, compiled from the literature or
from `DOCKING_PROTOCOL.md`. The remainder are not individually asserted; their
sources are given per stage in section 4.

---

## 7. What cannot be reconstructed

Stated explicitly, as Appendix A requires.

- **Figure 22** (`image14`), Orai Fr300 with docked ligands and transmembrane exclusion. Hand-composed PyMOL render, present
  nowhere in the repository. `pymol_pose_cluster_scene.py` builds the scene
  files, but the camera, colouring and render were interactive.
- **Figure 27** (`image12`), The docking receptor. Hand-composed PyMOL render, present
  nowhere in the repository. `pymol_pose_cluster_scene.py` builds the scene
  files, but the camera, colouring and render were interactive.
- **Figure 17**, the docking flow chart, is author-drawn in `Flow Charts.ipynb`.
- The **two unseeded DiffDock runs**, the benchmark and the Orai control. Their
  poses can be verified against the stored trees but cannot be redrawn.
- The **Fr0 receptor**, the **prepared ligand PDBQTs** and
  **`exhaustiveness_arm_status.json`**, for the reasons in their stage notes.

---

Findings and evidence: [`FINDINGS_2026-09-02.md`](FINDINGS_2026-09-02.md).  
Legacy commentary: [`REGENERATE_NOTES.md`](REGENERATE_NOTES.md).  
Pipeline design: `../../THESIS_REPRODUCTION.md`.
