# Nearest-copy endpoint: implementation log (branch `nearest-copy-endpoint`, 2026-09-08)

**Governing rule (from the author, 2026-09-08):** no data from the thesis is overwritten. Nothing under
`thesis_latex/` (tex, media, pdf), no existing directory under `posebusters_results/` or `pandamap_results/`,
and neither `thesis_expected_values.yaml` nor `thesis_assertions.py` nor any notebook has been modified.
Every new output lives in a NEW directory (`_nearest_copy_dev/`, `*_hfix`, `metal_stratum_dev/`) or under
this archive. Scripts changed keep their existing command lines bit-identical, and a cache they would have
rebuilt implicitly is now refused instead. Proposed thesis edits are collected, not applied, in
`thesis_latex/NEAREST_COPY_TEXT_CHANGES_2026-09-08.md`.

Verification of the rule at each checkpoint: `sed 's#pose_comparison_report.bak-instance-20260908#pose_comparison_report#' harness_signatures/backup_md5_20260908.txt | md5sum -c --quiet -`
(83 canonical CSVs identical), `git diff --stat -- thesis_latex Scripts/Analysis/thesis_expected_values.yaml
Scripts/Analysis/thesis_assertions.py` (empty), newest mtime under the canonical report dir 2026-09-03 and under
the canonical cluster dir 2026-08-31.

## Commits

| Commit | Content |
|---|---|
| `f4b17888` | Baseline snapshot of the working tree (the author's uncommitted edits, review and plan documents) |
| `c1ee68d8` | Phase 0.3, 0.4 and Phase 1 code, pre-review-fix checkpoint |
| (this) | Adversarial review fixes; implementation log |

## What was implemented (plan v2 references)

**Phase 0.2 freeze.** Seven harness inputs copied to `*.bak-instance-20260908` siblings (gitignored), the
working-tree yaml and `media/` copied under `harness_signatures/*.backup-instance-20260908` (the `.backup`
spelling keeps them outside the orphan-figure search of `thesis_assertions.py`), md5 list of the 83 canonical
CSVs, and `harness_signature_0_instance.csv` (1,364 rows, 0 failures).

**Phase 0.3, cluster-report defect (`pose_cluster_crystal_pocket_report.py`).** `load_heavy_atom_mol` now strips
hydrogens explicitly (RDKit ignores `removeHs=True` when `sanitize=False`), mirroring the hub's `Chem.RemoveHs`
so both generators share one heavy-atom basis; precision-at-one is printed from the exact fraction (the 4-dp
JSON value stays); `_CACHE_SCHEMA` 5 → 6. Run into the NEW directory
`cluster_crystal_pocket_matched_equibind/autodock_mgltools_exh128_gnina__diffdock_smina_allposes_hfix`
(70 s): any-pose 4 Å reach AutoDock 268, DiffDock 258, EquiBind 151, pooled 297, equal to the canonical
per-pose table (the canonical run gave 267 / 256 / 152 / 296); per-complex centroid distances agree with the
table to the 3-dp rounding of the CSV (max 5e-4 Å). Consequences for the printed thesis (Table 3, Figures 5, 6,
38, `:474`, `:476`) are listed in the plan, Section 1.7, and NOT applied.

**Phase 0.4, metal stratum (`metal_stratum.py`, new; `bench_metal_stratum` stage registered in
`_build_reproduction_notebook.py`, builder NOT executed).** Residue-level 5 Å rule on the reference instance:
73 / 78 / 81 at 4 / 5 / 6 Å; cofactor / ion 16 / 62 (App. C :165 prints 15 / 63); membership differs across
copies for 7JHQ_VAJ and 7TB0_UD1 only; the reference-convention contrasts reproduce App. C :165 verbatim.
Output in `posebusters_results/metal_stratum_dev/`.

**Phase 1, hub (`posebusters_pose_comparison.py`).** `--reference-convention {instance,nearest}`, default
`instance`; `load_all_mols`; per-pose argmin copy j* driving RMSD, centroid, PoseBusters RMSD and Kabsch,
rotation, torsions, contact recovery and PLIF (one ProLIF run, n_poses × n_copies matrix); nine provenance
columns appended (`reference_convention`, `n_copies`, `ref_copy_index`, `nearest_copy_index`,
`nearest_copy_is_ref`, `rmsd_ref_instance`, `centroid_dist_ref_instance`, `pb_rmsd_ref_instance`,
`bestfit_rmsd_ref_instance`); cache schema 6 with the convention in the signature and the manifest; the
mechanism exemplars and receptor-context PDBs load the pose's copy; sidecar wording routed through
`_ref_noun()` so instance output is byte-identical.

**Refusal semantics (D8, tightened after review).** A cache whose convention OR schema differs from the run is
refused with exit 2 in both the `--reuse-cache` path and the implicit path unless `--force`; legacy schema-5
manifests count as `instance`; other input changes (PB CSV content, id filter, variant flags) keep the
long-standing recompute. Consequence: the existing stage command against the canonical directory now REFUSES
instead of rebuilding it (the schema changed), which is the behaviour the author's rule requires; rebuilding
requires `--force` after a backup, or a new `--out-dir`.

## Tests and results

- T1 (instance smoke, 6 pairs, 4,840 poses, ids 5SAK_ZRY 6TW5_9M2 7A9E_R4W 7FHA_ADX 7JHQ_VAJ 7TUO_KL9; 5S8I_2LY
  and 7FRX_O88 have no rows in the PB CSV): every column equals the canonical table except `plif_recovery`
  (46 rows, max 0.125). The UNMODIFIED HEAD script differs from the canonical table on 55 plif rows, so this is
  pre-existing run-to-run nondeterminism of `obabel -h` in the ProLIF receptor load, not the change.
- T2 (nearest smoke, same ids): `rmsd <= rmsd_ref_instance` on every row, the `*_ref_instance` twins equal the
  canonical columns to 0.0, single-copy complex identical to the instance run; 2,189 of 4,840 poses have an
  alternate copy nearest; re-run after the review fixes identical except plif (78 rows).
- T3 refusals: (i) `--reuse-cache` under the other convention → exit 2; (ii) other convention without `--force`
  → exit 2; (ii-b) with a second differing key → exit 2; (iii) unchanged instance command → reuse (or the
  long-standing recompute when an input such as the ids path really differs); (iv) legacy schema-5 cache,
  default flags → exit 2, cache byte-identical afterwards; (v) legacy cache with `--reuse-cache` → exit 2.
- Full-cohort instance regression (`_nearest_copy_dev/full_instance/`, 303 pairs, 24 workers, stage flags,
  `--force --reference-convention instance`): RUNNING at the time of writing; result appended below.
- Reviews: two adversarial reviewers on the hub, one each on the cluster and metal changes; all major findings
  fixed (refusal breadth, legacy manifests, full regression scheduled); recorded minors that are NOT yet done:
  five sidecar strings still say "crystal" under `nearest` (instance output unaffected); no unit test yet for the
  refusals (Phase 7); `plif_recovery` nondeterminism root cause (deterministic protonation or a cached
  receptor) is a Phase 0 defect to schedule; `test_diffdock_sweep_gate.py` points at a notebook now under
  `obsolete/`.

## Not started (plan v2 Phases 2-8)

Phase 2 downstream scripts (cluster D5 any-copy rule, PandaMap per-copy crystal tasks and the five `groupby`
sites, exhaustiveness `CASCADE_PINS` re-pin, ITT stage, `bench_reference_convention`), Phase 3 rebuild under
`nearest` into a NEW directory, Phase 4 selection re-check, Phase 5 thesis text (proposals only), Phase 6
figures, Phase 7 harness and notebook, Phase 8 QA. The `PandaMap` pipeline already profiles only the top-5
PB-valid poses per method and complex (`07_benchmark_pandamap.yaml`: `poses_per_combo: 5`, `pb_valid_only:
true`), so the per-copy change adds only crystal tasks.
