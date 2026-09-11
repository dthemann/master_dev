# REGENERATE_NOTES.md — the hand-written regeneration guide

> **Renamed 2026-09-02.** `REGENERATE.md` is now generated from the pipeline's
> own stage registry by `regenerate_guide.py`, so its ordering, commands, paths
> and figure index cannot drift from what actually runs. This file keeps the
> per-command commentary and the trap notes, which are not derivable from the
> registry and are still worth reading. Where the two disagree on a PATH, the
> generated one is right; where this one adds a WARNING, it is still worth heeding.

## Original heading: how every reported figure and table is produced

This file records, for each figure and table reported in the thesis, the committed
output it was taken from and the exact command that rebuilds that output. It exists
because several thesis outputs are produced by a script invoked directly rather than
by a notebook cell, so the notebooks alone are not a complete build recipe.

Everything below was checked against the repository on 2026-08-16. Provenance marked
**verified** means the committed PNG is byte-identical (MD5) to the file embedded in
`thesis_latex/media/media/`. Anything that could not be reconstructed is marked
**UNRESOLVED** rather than given a guessed command.

Scope note: the analysis outputs referenced below live in the working tree, not in the git
index. `posebusters_results/`, `pandamap_results/` and `PoseBusters_Benchmark_Analysis/figures/`
are untracked, so "stored" and "byte-identical" throughout this guide mean identical to the
file on disk in the author's tree, not to a file recoverable from the repository history.

Thesis figure and table numbers are the automatic LaTeX numbers, read from
`thesis_latex/Thesis.lof` and `thesis_latex/Thesis.lot`.

---

> ## CORRECTION, 2026-09-02 — read this before using Section 2
>
> **Two things below are out of date, and both will send you to the wrong file.**
>
> **The numbering is from a retired document.** Section 2's figure numbers were read
> from `Thesis.lof`, the FULL build, which was retired on 2026-09-02 to
> `thesis_latex/obsolete/`. The thesis is now the SHORT build, which numbers its
> floats differently and has 39 figures against this guide's 36.
>
> **The paths predate the matched-EquiBind migration.** Every published figure has
> since moved to a `_matched` tree. Section 2 still points at the plain ones. This
> was established by md5-matching all 39 shipped assets against every PNG under the
> output roots, so it is measured rather than inferred.
>
> **A corrected, md5-verified figure index for the short build is in
> [`FINDINGS_2026-09-02.md`](FINDINGS_2026-09-02.md), section F.** Use that. The
> commands in Sections 3 and 5 below remain correct; only their paths need the
> substitution in the table.
>
> | Wherever this guide says | Read instead |
> | --- | --- |
> | `posebusters_results/benchmark_full_protein_vina_scoring` | `posebusters_results/benchmark_matched_equibind` |
> | `posebusters_results/cluster_crystal_pocket_full_protein` | `posebusters_results/cluster_crystal_pocket_matched_equibind` |
> | `pandamap_results/benchmark_full_protein_mgltools` | `pandamap_results/benchmark_matched_equibind` |
> | `posebusters_results/orai_benchmark` and `orai_jku` | `posebusters_results/_orai_matched_root/orai_benchmark` and `orai_jku` |
> | `posebusters_results/orai_pbvalid_tm_share_compare` | `posebusters_results/_orai_matched_root/orai_pbvalid_tm_share_compare` |
> | `posebusters_results/orai_pbvalid_yield_compare` | `posebusters_results/_orai_matched_root/orai_pbvalid_yield_compare` |
> | `pandamap_results/orai_interaction_compare` | `pandamap_results/orai_interaction_compare_matched` |
>
> The trap is that the directories reading like scratch are the current ones.
> Generation order, oldest first, is `*_PRE_FR0` and `*_PRE_EXH128` and
> `*_UFFON_backup_*`, then the plain name, then `_matched`.
>
> **Section 5's gap list is also out of date.** Four of its items were closed on
> 2026-09-02 and one figure that had no generator at all now has one. See
> `FINDINGS_2026-09-02.md`, sections C and D.
>
> **The pipeline that supersedes this guide** is `Thesis_Reproduction.ipynb`, which
> runs all 66 stages in dependency order and asserts 1,364 checks against the
> canonical trees. This guide remains useful for the per-command detail it carries.

## Analysis changes of 2026-08-28

**Parallel gnina optimiser (`run_autodock.py`).** gnina is invoked once per pose, and each
call pays CUDA context creation plus loading the 5-model `crossdock_default2018_ensemble`, so
the GPU sat ~10 % busy while a single process held an exclusive `flock`. The cost is startup,
not compute — measured on the exh128 arm, p75/p25 of per-pose time is **1.24×** and the
per-complex median varies only **1.20×** across wildly different ligands.

`_interprocess_gpu_lock(path, slots=N)` is now a counting semaphore over N `.slotN` lock files
rather than one exclusive lock, and the per-pose loop body moved into `_optimize_one_pose()`,
which returns counter *deltas* that are merged **in rank order** after the join — so neither
the optimisation log nor the re-ranking can depend on completion order. Two new config keys,
both defaulting to **1 = the previous strictly-serial behaviour**, so every existing arm is
unaffected:

```yaml
optimize_workers: 4      # poses of one complex optimised concurrently
gpu_max_concurrent: 4    # semaphore slots; must be >= optimize_workers
optimize_cpu: 6          # lower when raising workers (4 x 6 = 24 of 32 logical CPUs)
```

Validated by A/B on 4 complexes / 120 poses: **259 s serial → 76 s at 4 workers (3.41×)**, with
the two optimisation logs **bit-identical** on `optimized_rank` and every score column. Applies
to `gnina` and `gnina_refinement`; smina takes no GPU lock and is unchanged.

> **Cost-accounting caveat.** The published arms were produced at `optimize_workers: 1`, where
> optimiser intervals provably never overlap and `sum(elapsed_time_s)` therefore equals wall
> clock. At `optimize_workers > 1` that sum is *compute* time, roughly N× the wall clock. Do
> not compare optimiser cost across arms run at different worker counts without saying so.

**New analysis: marginal return on compute for the exhaustiveness ladder.**
`autodock_exhaustiveness_returns.py` answers what each step up the raw AutoDock Vina
exhaustiveness ladder buys per hour of search, separately on each of the three gates and on
their conjunction. It is not (yet) a thesis figure, so it has no Figure-index row; it is
recorded here because it is a script invoked directly rather than from a notebook cell.

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Analysis/autodock_exhaustiveness_returns.py \
    --per-pose-csv posebusters_results/benchmark_full_protein_vina_scoring/dock/pose_comparison_report/per_pose_metrics.csv \
    --out-dir posebusters_results/autodock_exhaustiveness_returns \
    --n-boot 5000 --verify-input-parity
```

The ladder is **not hardcoded**. Arms are read from the `docking_directories` block of
`Scripts/Docking/Posebusters/posebusters_benchmark_full_protein_config.yaml`, intersected with
the method keys present in `per_pose_metrics.csv` and with the trees that carry
`docking_log_*.csv`. exh128 was registered in that config but unscored when the script was
written, and joined the analysis with no code change the moment its poses were busted on
2026-08-28. Each arm's `cpus_per_worker`, `exhaustiveness` and `optimization` are read at run
time from the docking config whose `output_dir` **is** that arm's tree, so a renamed config is
still found and a config pointing elsewhere is never misattributed.

`--verify-input-parity` sha256s the staged ligand PDBQT and box per complex per arm and reports
whether exhaustiveness really is the only variable (today: identical, 303/303, five arms). Add
`--full-parity` to hash the ~1 MB receptors as well.

**The script hard-fails if it stops reproducing the published cascade.** `CASCADE_PINS` holds the
exh18 and exh92 rows of `pose_validity_cascade_report.txt` and `assert_cascade()` compares all 15
columns, rounding once at comparison time exactly as the cascade rounds once at display time. A
mismatch raises rather than warns, because every number in the report below it would then be
unmoored from the thesis.

**Palette departure, deliberate.** These figures do **not** use the olive exhaustiveness ramp
from `posebusters_pose_comparison.py:850-854`. That ramp exists to avoid colliding with the
other engine families drawn in the same axes in Figure 09f; these figures carry AutoDock alone,
so the constraint does not apply, and the olive ramp's lightest rung sits at **1.19:1 contrast
on white** — invisible as a line. The replacement is viridis sampled over `[0.02, 0.58]`
(`ARM_COLOR`), measured at minimum **3.16:1** on white, with better minimum adjacent-rung
separation under deuteranopia (0.124 → 0.151) and protanopia (0.124 → 0.176). Colour encodes
the **rung** in all three figures and nothing else; gates are carried by marker shape, line
style and hatch. Orange `#d95f02` is reserved for the "not cost-comparable" flag. Rung colours
are pinned per exhaustiveness value, so a subset run (`--min-exhaustiveness 64`) does not
recolour the rungs it keeps.

Outputs land in `posebusters_results/autodock_exhaustiveness_returns/`: ten CSVs, a
self-describing `exh_returns_report.txt` with a `SOURCES` block, `exh_run_manifest.json`, and
three figures under `figures/` each with a `*_stats.txt` sidecar carrying the inferential detail
that the panels deliberately omit. See `IMPLEMENTED_STATS_TESTS.md` §N for the tests, and for the
list of tests the script **refuses** to run and says so in its own output.

## Analysis changes of 2026-08-29 (thesis migration, later the same day)

**The thesis Benchmark chapter now reports the MGLTools/ADFRsuite exhaustiveness ladder.**
`thesis_latex/body_main_short.tex` and `body_appendix_short.tex` were migrated to
`autodock_mgltools_exh128_gnina`; backups are `*.bak-premgltools-20260829`. Seven figures were
re-copied into `thesis_latex/media/media/`: image2 from `validity_report_mgltools/`, image3 and
image4 from `pose_comparison_report/18_topn_within_thresholds_{,kabsch_}pbvalid_depths.png`, image5
from the regenerated filmstrip, image6 and image7 from
`cluster_crystal_pocket_full_protein/autodock_mgltools_exh128_gnina__diffdock_smina_allposes/`, and
image9 from the new PandaMap report. That media folder is **shared with the full thesis**, whose
captions are now stale by design.

**BUG FIXED — `run_pandamap.py` could not pin any MGLTools arm.** `autodock_label()` folded every
`autodock_mgltools*` tree onto plain `autodock`, so `autodock_variant: autodock_mgltools_exh128_gnina`
printed "no AutoDock variant matches" and silently fell back to keeping **all six Meeko arms**, with
poses from several trees pooled under one method label. The fold now keeps the tree base verbatim.
Trees that were already `autodock` / `autodock_vinardo` keep exactly the keys they had, so published
runs are unaffected. Whole-protein fingerprints for the new arm:

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Analysis/run_pandamap.py \
    --config Scripts/Analysis/pandamap_benchmark_full_protein_mgltools_config.yaml
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Analysis/pandamap_interaction_report.py \
    --config Scripts/Analysis/pandamap_benchmark_full_protein_mgltools_config.yaml \
    --in-dir pandamap_results/benchmark_full_protein_mgltools \
    --exclude-methods unidock2,autodock_gnina
# THIRD STEP, required — the appendix "Interaction Fingerprints" quotes its output verbatim:
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Analysis/interaction_pose_basis_audit.py \
    --config Scripts/Analysis/pandamap_benchmark_full_protein_mgltools_config.yaml \
    --exclude-methods unidock2,autodock_gnina
```

`interaction_pose_basis_audit.py` writes `report/interaction_pose_basis_audit.{txt,csv}` and is the
sole provenance for four things the report script cannot compute, because it never joins the per-pose
RMSD table: the composition of the profiled pool (median 7.62 Å, 20.7 % within 2 Å), the placement
share of the between-tool ordering (80.9 %), the site-filtered top-5 union Jaccard (rmsd < 10 Å), and
the `metal_coordination` position-independence defect. **Every number in appendix section
"Interaction Fingerprints" comes from that sidecar**, so re-run it whenever this arm is regenerated
or those numbers go stale silently. It aborts rather than guessing if the pandamap arm and the
pose-comparison report are out of step. Notebook cell id `pandamap-pose-basis-audit`.

The output dir is seeded from a copy of `benchmark_full_protein/`, so `overwrite: false` recomputes
only the ~1,500 new AutoDock poses and reuses the DiffDock/EquiBind/Uni-Dock2 fingerprints. The
`autodock_gnina` rows carried over from that copy are the superseded Meeko arm and must be excluded
at report time, which is what the flag above does.

**`filmstrip_rank1_top5_top15.py` also hard-pinned the Meeko arm.** `FORCED_AUTODOCK` now reads
`posebusters_pose_comparison._DOMINANT_AUTODOCK_OPT`. The assignment had to move BELOW the
`import posebusters_pose_comparison as P` line — at module top it raises `NameError`. Re-run with
`--exclude-preset meeko`; it writes the `__stats_*.csv` sidecars the thesis placement/form table reads.

**Figure 2 needs a double exclusion, not just the preset.** `posebusters_validity_report.load_and_score`
folds every surviving `autodock_mgltools*` key onto plain `autodock`, so `--exclude-preset meeko`
alone would draw a panel pooling five ladder rungs. Exclude the other rungs explicitly so the fold
sees exactly one raw arm and one rescored arm:

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Analysis/figure2_validity_yield.py \
    --csv posebusters_results/benchmark_full_protein_vina_scoring/dock/posebusters_filtered_results.csv \
    --out posebusters_results/benchmark_full_protein_vina_scoring/dock/validity_report_mgltools/00_figure2_validity_yield.png \
    --exclude-preset meeko \
    --exclude-methods autodock_mgltools,autodock_mgltools_gnina,autodock_mgltools_exh18,autodock_mgltools_exh64,autodock_mgltools_exh64_gnina,autodock_mgltools_exh92
```

**`docking_effort_comparison.py` can express the new arm as of 2026-08-30.** It previously could
not: `_per_pose_methods()` derived the AutoDock per-pose label as `f"autodock_{ad_refine}"`, so only
the four Meeko labels were reachable and `--autodock-prep mgl_tools` changed the directory segment
only, which would have joined exh128 timings to Meeko pose counts. Two flags fix that, and BOTH are
required against the exh128 tree:

* `--autodock-method autodock_mgltools_exh128_gnina` — names the per-pose arm outright instead of
  deriving it. Without it the join silently returns the wrong pose count (8,976 vs 9,043).
* `--autodock-optimizer-workers 16` — a concurrency divisor. `_autodock_effort` adds Σ per-pose
  optimiser elapsed to the search wall, an identity that holds only for a serial optimiser, and this
  arm ran `optimize_workers: 16`. Measured over the 9,126 `optimized_gnina/*.provenance.json`
  sidecars: serial sum **37,693.7 s** against a true wall of **2,587.3 s** (union of the per-pose
  intervals; 2,611.5 s under the opposite created_at convention, first-to-last span **2,681.9 s** —
  that span is the ~2,677.7 s figure this document used to quote). So the sum overstates the wall by
  **~14.5×**. Dividing by 16 lands ~10 % below the measured wall because the speed-up is sublinear;
  it is an approximation, not a measurement. Resource-seconds are deliberately NOT divided — sixteen
  concurrent gnina calls still bill sixteen calls' worth of GPU-s and CPU-core-s.

Both flags are wired into notebook cell `4821046f`, so `docking_effort_gnina/` is now regenerated
from the notebook rather than recomputed by hand. Do NOT read the `timestamp` column of
`optimization_log.csv` for wall clock — it is a per-complex batch write.

`Dockings/vina_results_full_protein_vina_scoring_mgltools_exh128/exhaustiveness_arm_status.json`
carries these timings again. It had been overwritten by an abandoned Variant-B `gnina_refinement`
pass that left it reading `status: running, processed: 5`; that pass is quarantined under
`_ABANDONED_gnina_refinement_20260830/` and the status file was rewritten from the sidecars.

## Analysis changes of 2026-08-29

**Single-point method exclusion across the Benchmark suite (`method_filter.py`).**
Every Benchmark analysis script now takes the same three flags, wired at exactly one
place per script:

    --exclude-methods PAT[,PAT...]   drop these method keys
    --exclude-preset  NAME           drop a named, provenance-documented group
    --exclude-strict                 make a pattern that matches nothing fatal

Unset, all three are inert: a run with none of them produces byte-identical output to
the pre-change code (verified 2026-08-29 — 91 of 91 comparable outputs plus the
`pb_valid/` subtree identical, 92 PNGs each, on the whole-protein pose-comparison
report). So every committed command line still reproduces what it used to.

Why a shared filter instead of editing the keep-sets: `posebusters_pose_comparison.py`
alone enumerates variants at two independent points (`_select_presentation_tools` and
`_all_variant_yield_specs`), and the per-family spec tables add more. Excluding an arm
by editing those one at a time reliably leaves it alive in some figure nobody
re-checked. Filtering the frame itself cannot.

Three rules the module enforces, each a bug this suite has already shipped:

* **Exact keys, never bare prefixes.** `autodock` matches `autodock` and nothing else —
  not `autodock_gnina`, not `autodock_mgltools`. A family needs an explicit glob
  (`autodock_vinardo*`) so the breadth is visible at the call site.
* **A pattern that matches nothing is reported**, never silently ignored.
* **The filter runs before any key-folding step.** This matters most in
  `posebusters_validity_report.py`, whose `_autodock_scoring_base` rewrites every
  `autodock_mgltools*` key onto plain `autodock`. Filtering after that fold would take
  the ADFRsuite arms out along with the Meeko ones. `load_and_score` therefore applies
  the exclusion itself, before its splits, and `Scripts/Analysis/tests/test_method_filter.py`
  asserts that the wrong order actively fails.

The `meeko` preset drops the arms docked from Meeko-prepared ligands: `autodock`,
`autodock_smina`, `autodock_gnina`, `autodock_gnina_refinement`, `autodock_vinardo*`.
Verify it against the trees themselves rather than trusting the list:

```bash
python Scripts/Analysis/method_filter.py --verify-preset meeko \
    --config Scripts/Docking/Posebusters/posebusters_benchmark_full_protein_config.yaml
```

That fingerprints the docked poses in every tree the config reads and exits non-zero if
a Meeko tree is not covered or a non-Meeko tree is. **Note the preset removes the
Vinardo scoring contrast entirely**, because the only Vinardo tree on this dataset is
Meeko-ligand. Receptors are not part of it: all four whole-protein AutoDock trees were
prepared with ADFRsuite, and DiffDock / EquiBind / Uni-Dock2 hold no PDBQT at all.

Each run writes `method_filter.json` beside its figures recording the preset, its
provenance string, and the dropped and kept keys — so a missing arm can be told apart
from an arm that was never run.

Two wiring details worth knowing. `pose_cluster_crystal_pocket_report.py` keys its
pickle cache on `_analysis_signature`, which now includes the exclusion; without that a
cached run would silently return the unfiltered analysis. And `docking_effort_comparison.py`
filters inside `build_table`, not in `main()`, because the cross-dataset path re-enters
`build_table` with a cloned namespace that a `main()`-level filter would miss.

`pandamap_interaction_report.py` had its own `--exclude-methods`; it now uses the shared
one. The flag name and its exact-match semantics are unchanged, so committed PandaMap
commands still work.

**Dominant AutoDock arm is now `autodock_mgltools_exh128_gnina`** (ADFRsuite ligands,
exhaustiveness 128, gnina-rescored). It replaced the Meeko `autodock_gnina`, which the
exclusion above removes. On rank-1 poses that are near-native and PoseBusters-valid it
reaches 36.6% (111/303) against 29.4% (89/303), and 65.3% against 49.8% at top-15.

Two constants in `posebusters_pose_comparison.py` are the single source of truth —
`_DOMINANT_AUTODOCK_RAW` and `_DOMINANT_AUTODOCK_OPT`. The pin allowlist
(`_PINNABLE_AUTODOCK_ARMS`), the 09b/09c/09d/09f and 21/22 spec tables and
`_oracle_variant_specs` all read them, so moving the arm again is a two-line change.
`pose_cluster_crystal_pocket_report.py` carries its own copy in `_STARRED_VARIANT`.

**`--collapse-plots-only` is REQUIRED for the pin to do anything.** `_select_autodock_arm`
and `_select_presentation_tools` run only under that flag (or `--best-variants-only`).
Cell 33 used to pass `--collapse-autodock-variant` without it, so the cell's own comment
about collapsing to primary variants had never actually been true. It now passes both.

Three things had to be repaired for the AutoDock family to survive the switch, because
each named the Meeko keys and nothing else:
* the 09d, 09f, 09c/09e and 21/22 spec tables — figure 09d drew no AutoDock at all and
  figures 21/22 lost their cross-tool statistics entirely (Friedman needs three tools);
* `_oracle_variant_specs`, which hard-coded the literal `"autodock"` and so dropped the
  family from 09b/09e;
* `stats_utils.paired_continuous`, which called Friedman unconditionally and RAISED on a
  two-condition family. The dominant arm has only raw + gnina at its ladder point (no
  smina, no CNN-refine), and that exception was taking the whole 09f payload down with
  it, DiffDock and EquiBind included. It now returns a NaN omnibus for k=2 and still runs
  the pairwise Wilcoxon, which is the actual raw-vs-optimised question.

**Stale outputs are NOT cleared by `--force`.** After the first Meeko-free run the report
directory held 28 files from other commands (the `*_all_variants` tables, the
`*_gnina_arm.*` outputs, the `20d_*filmstrip*` set, `bounded_claim.csv`,
`validity_gate_cost.csv`) that still listed the excluded arms, with nothing marking them
as older. They were deleted on 2026-08-29; the pre-change state is in
`pose_comparison_report.bak-premeeko-20260829`.

**Cell 37 follows cell 33.** `pose_cluster_crystal_pocket_report.py` pinned the Meeko
`autodock_gnina` and would have drawn Figures 6 and 7 on an excluded arm WITHOUT error,
because `per_pose_metrics.csv` deliberately retains every variant. It now pins the
dominant arm and writes to
`cluster_crystal_pocket_full_protein/autodock_mgltools_exh128_gnina__diffdock_smina_allposes`
(the out-dir name is built from the variant, so the old run survives beside it). New run:
mean 11.09 sites/complex, precision@1 0.571, purity 0.899 — against 12.17 / 0.554 / 0.889.

**TRAP the validation exposed — the validity report POOLS whatever survives.** Running
`posebusters_validity_report.py --exclude-preset meeko` removes the two Meeko trees
correctly, but the fold then collapses the five surviving ADFRsuite arms into a row
labelled `autodock` (44,933 poses = exh18 + exh32 + exh64 + exh92 + exh128) and the
three gnina-bearing ones into `autodock_gnina` (26,963 poses). The label names one arm;
the row is five. Nothing was changed to alter that pooling, because doing so would move
published numbers — but `_apply_autodock_split` now prints a WARNING naming every arm it
pooled and the pose count, so the row can no longer be read as a single arm by accident.
The warning fires on ordinary runs too, since the pooling predates this work. Read
`summary_per_tool.csv` from that script as "AutoDock family, pooled" until the fold is
fixed; `posebusters_pose_comparison.py` keeps the arms separate and is the place to
compare them.

## Analysis changes of 2026-08-21

**Cluster-report figure labels are now symmetric (Figures 6 and 7).**
`pose_cluster_crystal_pocket_report.py` relabelled only the AutoDock slot
(`_TOOL_DISPLAY["autodock"] = _autodock_display(ad_variant)` -> "AutoDock Vina + gnina") while
the equally optimised DiffDock and EquiBind slots stayed bare on every axis. AutoDock therefore
read as the sole optimised arm on figures where all three were optimised. `_diffdock_display()`
and `_equibind_display()` now sit beside `_autodock_display()`, and all three slots are set
together, so the run pinned to `autodock_gnina` / `diffdock_smina` / `equibind_unguided_gnina`
draws **AutoDock\* / DiffDock\* / EquiBind\***, the same star convention `posebusters_pose_comparison.py`
adopted on 2026-08-20 for Figures 3 and 4 (see [§5.1](#51-figures-2-3-and-4-pin-the-autodock-arm-or-you-rebuild-raw-vina)).
A starred label means exactly that arm; any other variant is still spelled out, so
`autodock_gnina_refinement` remains "AutoDock Vina + gnina (refine)" and can never be read as
the starred default.

Re-render only — the run reuses `analysis_cache.pkl`, so nothing is recomputed:

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Analysis/pose_cluster_crystal_pocket_report.py \
    --ids-file "Data/PoseBuster Benchmark Set/posebusters_pdb_ccd_ids.txt" \
    --per-pose-csv posebusters_results/benchmark_full_protein_vina_scoring/dock/pose_comparison_report/per_pose_metrics.csv \
    --diffdock-variant diffdock_smina \
    --out-dir posebusters_results/cluster_crystal_pocket_full_protein/autodock_gnina__diffdock_smina_allposes \
    --site-cluster threshold --pocket-radius 8.0 --rank-by consensus --stability-boot 25 \
    --workers 30 --stats \
    --autodock-variant autodock_gnina --equibind-variant equibind_unguided_gnina
```

Regression check: every CSV and JSON in the out-dir comes back byte-identical **except**
`rank1_quality_crosstab.csv` and `crystal_cluster_homogeneity_stats.json`, which carry the tool
label as a display string and change only there. `cluster_quality_metrics.png` (**Figure 8**) is
unchanged, because its panels carry no tool axis. Nine PNGs move; only
`topN_crystal_cluster_matrix.png` and `crystal_cluster_homogeneity.png` are thesis figures and
were copied to `thesis_latex/media/media/image6.png` and `image7.png`.

Both thesis documents share that media folder, so the captions of Figures 6 and 7 were updated in
`body_main_short.tex` **and** `body_main.tex`; the latter previously claimed the axis spells out
"AutoDock Vina + gnina", which the reprinted figure no longer does. `body_main.tex`'s Figure 3
caption also gained the star legend it had been missing since the 2026-08-20 relabel.

## Analysis changes of 2026-08-19

Three corrections from the pre-submission audit. Any output produced before this date is superseded.

1. `orai_pandamap_interaction_compare.py` — the **residue hot-spot axis is now pose-weighted**.
   `_top_shared_residues` previously summed the six per-(dataset, tool) contact fractions without
   weighting by group size, so a residue prominent in a one-pose group scored as highly as one
   contacted in every pose of a 3,524-pose group. Asp287 entered the axis at rank 7 on a single
   experimental EquiBind pose. The default is now `mode="pooled"`, ranking by the pooled contact
   fraction over all 7,231 poses, exposed as `--residue-axis {pooled,macro}`; `macro` reproduces the
   superseded axis exactly. The new axis drops Asp287 and Glu106 and adds Glu173 and Lys85.
   Affects **Figure 23** and `pandamap_residue_hotspots.csv` only — re-running leaves
   `fig_total_interactions_compare.png`, `fig_type_profile_compare.png` and
   `fig_fingerprint_overlap_compare.png` byte-identical, which is the isolation control.

2. `PoseBusters_DataSet_Analysis.ipynb` — **`_load_jku_mol` hydrogen basis fixed**. It passed
   `sanitize=False`, which suppresses RDKit's implicit `RemoveHs`, so the three Orai ligands were
   described on a hydrogen-explicit basis while the 308 benchmark ligands were hydrogen-suppressed.
   Exactly two of the sixteen descriptors are basis-dependent: `rot_bonds` (2abp-NH2 6→5,
   Synta-66 7→5, GSK-7975A 5→4) and `n_element_types` (5→4, printed nowhere). The loader now returns
   `Chem.RemoveHs(m)`, with stereochemistry still assigned on the explicit-H molecule first.
   Affects **Table 3**, **Table 17**, the appendix "four to five rotatable bonds" range and the
   **Figure 31 caption** (`thesis_latex/body_appendix_short.tex:740`; the full build's `body_appendix.tex:491` was retired to `thesis_latex/obsolete/` on 2026-09-02), see item 3 below.
   The benchmark statistics and the PCA fit are unaffected, since the PCA is fit on the benchmark
   alone — PC1/PC2/PC3/PC4 remain 45.0 / 23.3 / 6.7 / 6.4 %.

3. **The appendix chemical-space figures were re-rendered** from the corrected notebook and re-embedded as
   `image33` to `image42`. Note the mapping is **not** sequential: `image39`→`09_pca_scree_loadings`,
   `image40`→`10_pca_biplot`, `image41`→`07_receptor_profile`, `image42`→`08_startconf_rmsd`.
   Derive it by MD5 against `PoseBusters_Benchmark_Analysis/figures/`, never by filename order.
   All ten changed byte-wise, but only the JKU overlay in `01`, the JKU points in `05` and the
   JKU points in `10` changed for data reasons. The `05` change is the largest and was missed on
   the first pass: dropping Synta-66 from 7 to 5 and GSK-7975A from 5 to 4 rotatable bonds
   collapsed their vertical gap in panel (A) to `|dy|/yspan = 0.0351`, under the `tol = 0.04`
   coincidence tolerance in `overlay_jku_points`, so the pair now clusters and is fanned apart by
   `0.035 * xspan = 2.93` heavy atoms each way. Because that exceeds half their true 2-heavy-atom
   gap, the two markers cross over and panel (A) draws them in reversed left-to-right order. The
   Figure 31 caption had named panel (A) as exact and was corrected on 2026-08-19. Note both
   panels (A) and (B) sit within 12 % and 2 % of the `tol` cut, so any future re-render that
   nudges the axis limits can flip a panel in or out of the fan and silently re-stale that caption. The rest is matplotlib rendering non-determinism, roughly 1 to 2.5 % of pixels
   and a 3-pixel canvas change on `04`. The three Orai points moved 0.36 to 0.43 PC units against
   PC1/PC2 standard deviations of 2.68 and 1.93, so all three remain at negative PC1 and positive
   PC2 and no interpretation changes.

---

## Analysis changes of 2026-08-18

Two scripts were corrected after an audit found the prose describing something the
code did not do. Both were re-run and their figures re-embedded, so any output
produced before this date is superseded.

1. `pose_cluster_crystal_pocket_report.py` — crystal-cluster **reach is now gated at
   the 4 Å centroid threshold**. Previously `correct_label` was the crystal-*nearest*
   cluster computed unconditionally, and `correct_cluster_is_hit` was stored but never
   used, so 21 of 303 complexes counted as reached with their nearest cluster 4.01 to
   21.82 Å away. A complex failing the gate now carries no correct cluster at all, so it
   contributes neither reach nor composition. Affects Table 10 and Figures 6 and 7
   (reach, co-reach φ, all-three counts, panels A to D). It does **not** affect Figure 8,
   whose ranking-ablation numbers (oracle 93.07 %, consensus 55.45 %, size 48.84 %) were
   already computed at the same threshold and are byte-identical after the change.
   Re-run with `--force`, otherwise the cached analysis is reused.

2. `orai_pandamap_interaction_compare.py` — the **typed Jaccard now uses the same
   ≥10 % characteristic-contact threshold as the residue-level one** (`_char_typed_set`).
   Previously the residue set was thresholded and the typed set was not, so the
   residue-to-typed drop conflated a change of contact type with a change of basis. The
   superseded unthresholded value is retained in `pandamap_fingerprint_overlap.csv` as
   `jaccard_typed_unthresholded` for provenance. Affects Figure 24 only.

---

## 1. Environments

| Env | Interpreter | Used for |
|---|---|---|
| `vina` | `/home/manndo/anaconda3/envs/vina/bin/python` | every script and notebook in this file. It is the only env carrying pandas, matplotlib, RDKit, PoseBusters and PandaMap together. |
| `diffdock`, `equibind`, `unidock`, `unidock2` | see `Scripts/Docking/` | re-docking only. No analysis in this file needs them. |

All commands are run from the repository root `/home/manndo/master_dev`. Activate with

```bash
source ~/anaconda3/etc/profile.d/conda.sh && conda activate vina
```

The notebooks use the same env through the literal path `VINA_PY`, so a notebook cell
and its command-line equivalent are interchangeable.

---

## 2. Figure index

> **SUPERSEDED 2026-09-02.** These numbers are the retired full build's and these
> paths predate the matched-EquiBind migration. The corrected, md5-verified index
> for the short build is in [`FINDINGS_2026-09-02.md`](FINDINGS_2026-09-02.md),
> section F. The rows below are kept for their per-figure commentary only.

| Fig | Subject | Output file | Produced by | Provenance |
|---|---|---|---|---|
| 1 | Physics-based docking flow chart | `docking_workflow.png` | `Flow Charts.ipynb` (author-drawn schematic) | verified |
| 2 | Post-hoc optimisation and PoseBusters validity | `posebusters_results/benchmark_full_protein_vina_scoring/dock/validity_report_mgltools/00_figure2_validity_yield.png` | `figure2_validity_yield.py`, see [§5.8](#58-figure-2--resolved-2026-08-18-generator-added) | verified, byte-identical to the shipped `image2.png` on 2026-08-30. Needs the double exclusion, and mind the **`_mgltools`** suffix — the copy in plain `validity_report/` was the superseded Meeko arm and has been retired |
| 3 | Best-of-top-N accuracy, as-placed RMSD | same report dir | see [§5.1](#51-figures-2-3-and-4-pin-the-autodock-arm-or-you-rebuild-raw-vina) | OK, needs `--collapse-autodock-variant autodock_gnina` |
| 4 | Best-of-top-N accuracy, Kabsch RMSD | same report dir | see [§5.1](#51-figures-2-3-and-4-pin-the-autodock-arm-or-you-rebuild-raw-vina) | OK, needs `--collapse-autodock-variant autodock_gnina` |
| 5 | Form versus placement filmstrip | `…/pose_comparison_report/20d_form_vs_placement_by_family__depth_filmstrip_pbvalid__rank1_top5_top15__thesis.png` | [§3.3](#33-form-versus-placement-filmstrip-figure-5) | verified |
| 6 | Top-N crystal-cluster co-recovery | `posebusters_results/cluster_crystal_pocket_full_protein/autodock_mgltools_exh128_gnina__diffdock_smina_allposes/topN_crystal_cluster_matrix.png` | `Master_Docking_AD_Full_Protein.ipynb` cell `6b9fe31d` (`pose_cluster_crystal_pocket_report.py`, `PB_VALID_ONLY` off so the out-dir tag is `allposes`) | verified, relabelled `AutoDock*`/`DiffDock*`/`EquiBind*` 2026-08-21. **PATH CORRECTED 2026-08-30 — the old `autodock_gnina__…` path is the SUPERSEDED arm; verified byte-identical to the shipped image** |
| 7 | Crystal-closest cluster composition | same dir, `crystal_cluster_homogeneity.png` | same cell | verified, relabelled 2026-08-21 |
| 8 | Cluster-quality assessment | same dir, `cluster_quality_metrics.png` | same cell | verified, no tool axis so the 2026-08-21 relabel leaves it byte-identical |
| 9 | Mean Jaccard fingerprint similarity | `pandamap_results/benchmark_full_protein_mgltools/report/05_fingerprint_similarity_top5.png` | [§3.2](#32-whole-protein-pandamap-interaction-fingerprints-figures-9-10-and-11) | verified. **PATH CORRECTED 2026-08-30 — mind the `_mgltools` suffix, as for Figure 2; the plain `benchmark_full_protein/report/` is the superseded arm** |
| 10 | Native-interaction recovery by rank | same dir, `04c_native_recovery_by_rank.png` | [§3.2](#32-whole-protein-pandamap-interaction-fingerprints-figures-9-10-and-11) | verified. **Figure 10 was STALE until 2026-08-30 — it still carried the plain (superseded) arm while Figures 9 and 11 came from `_mgltools`, a mixed-provenance bug within one figure family** |
| 11 | Matched, missed and spurious contacts by rank | same dir, `10b_contact_decomposition_by_rank.png` | [§3.2](#32-whole-protein-pandamap-interaction-fingerprints-figures-9-10-and-11) | verified |
| 12 | Orai1 docking receptor (frame 300) | not in repository | hand-composed PyMOL render, see [§5.4](#54-hand-composed-pymol-renders-figures-12-and-14) | **UNRESOLVED** |
| 13 | PoseBusters-valid yield, benchmark against experimental | `posebusters_results/orai_pbvalid_yield_compare/09f_pbvalid_yield_dominant_orai_benchmark_vs_experimental.png` | [§3.5](#35-orai-posebusters-valid-yield-compare-figure-13) | verified |
| 14 | Orai Fr300 with docked ligands and TM exclusion | not in repository | hand-composed PyMOL render, see [§5.4](#54-hand-composed-pymol-renders-figures-12-and-14) | **UNRESOLVED** |
| 15 | Usable pose yield per frame-ligand unit | `posebusters_results/orai_pbvalid_tm_share_compare/fig_pbvalid_outside_tm_whisker.png` | [§3.4](#34-orai-validity-and-transmembrane-share-compare-figures-15-and-16) | verified |
| 16 | Transmembrane loss of PoseBusters-valid poses | same dir, `fig_tm_loss_relative_compare.png` | [§3.4](#34-orai-validity-and-transmembrane-share-compare-figures-15-and-16) | verified |
| 17 | Inter-tool consensus-site distance histograms | `posebusters_results/orai_jku/pose_clusters/panels/orai_cross_tool_agreement_compare.png` | `Master_Docking.ipynb` cell `orai-xtool-compare-code` (`orai_cross_tool_agreement_compare.py`). **Not** the same-id cell in `Master_Docking_AD_Full_Protein.ipynb`, see [§5.7](#57-the-two-orai-xtool-compare-code-cells-read-different-inputs) | verified |
| 18 | Cross-tool agreement on Orai | same dir, `orai_tool_agreement_compare.png` | same cell, and the same warning in [§5.7](#57-the-two-orai-xtool-compare-code-cells-read-different-inputs) | verified |
| 19 | Reference-free cluster quality per tool | same dir, `orai_cluster_quality_compare.png` | same cell, and the same warning in [§5.7](#57-the-two-orai-xtool-compare-code-cells-read-different-inputs) | verified |
| 20 | Validity and placement filtering versus cluster quality | `posebusters_results/orai_benchmark/pose_clusters/orai_cluster_quality_filtering.png` | `Master_Docking_AD_Full_Protein.ipynb` cell `577d6b41` (`orai_pose_cluster_report.py --cluster-quality --top-n-poses 10 --autodock-variant gnina`) | dimension match only, not byte-identical. The arm is gnina, see [§5.5](#55-the-orai-autodock-arm-all-published-figures-are-gnina) |
| 21 | Total contacts per pose on Orai | `pandamap_results/orai_interaction_compare/fig_total_interactions_compare.png` | `Master_Docking_AD_Full_Protein.ipynb` cell `60097bb6` (`orai_pandamap_interaction_compare.py`) | verified |
| 22 | Interaction-type profile on Orai | same dir, `fig_type_profile_compare.png` | same cell | verified |
| 23 | Orai residue hot-spots | same dir, `fig_residue_hotspots_compare.png` | same cell | verified |
| 24 | Contact-fingerprint overlap between ligand sets | same dir, `fig_fingerprint_overlap_compare.png` | same cell | verified |
| 25 | Charged seconds per qualifying pose | `posebusters_results/benchmark_matched_equibind/docking_effort_gnina_v2_charged/panels/effort_by_quality_near2.png` | [§3.1b](#31b-charged-basis-cost-chapter-matched-equibind-arm-figures-25-and-26) | verified 2026-09-01, md5 `d5ab5f948537dca2bc19bc4548fd2d10` |
| 26 | Hardware resource per near-native valid pose | `posebusters_results/benchmark_matched_equibind/docking_effort_gnina_v2/panels/resource_per_near_native_valid_pose.png` | [§3.1b](#31b-charged-basis-cost-chapter-matched-equibind-arm-figures-25-and-26) | verified 2026-09-01, md5 `8ca0604e5711b9443d430778fbe2737f` |
| 27 | Univariate descriptor distributions | `PoseBusters_Benchmark_Analysis/figures/01_ligand_univariate_distributions.png` | [§3.6](#36-appendix-chemical-space-figures-figures-27-to-36) | byte-identical to FIG_DIR, see §5.3 |
| 28 | Drug-likeness | same dir, `02_druglikeness_ro5_qed.png` | [§3.6](#36-appendix-chemical-space-figures-figures-27-to-36) | byte-identical to FIG_DIR, see §5.3 |
| 29 | Elemental and charge composition | same dir, `03_elemental_charge_composition.png` | [§3.6](#36-appendix-chemical-space-figures-figures-27-to-36) | byte-identical to FIG_DIR, see §5.3 |
| 30 | Descriptor correlation heat map | same dir, `04_ligand_correlation_heatmap.png` | [§3.6](#36-appendix-chemical-space-figures-figures-27-to-36) | byte-identical to FIG_DIR, see §5.3 |
| 31 | Two-dimensional attribute mappings | same dir, `05_2d_attribute_mappings.png` | [§3.6](#36-appendix-chemical-space-figures-figures-27-to-36) | byte-identical to FIG_DIR, see §5.3 |
| 32 | Pairplot of six core descriptors | same dir, `06_pairplot_core_descriptors.png` | [§3.6](#36-appendix-chemical-space-figures-figures-27-to-36) | byte-identical to FIG_DIR, see §5.3 |
| 33 | PCA scree and loadings | same dir, `09_pca_scree_loadings.png` | [§3.6](#36-appendix-chemical-space-figures-figures-27-to-36) | byte-identical to FIG_DIR, see §5.3 |
| 34 | PCA biplot | same dir, `10_pca_biplot.png` | [§3.6](#36-appendix-chemical-space-figures-figures-27-to-36) | byte-identical to FIG_DIR, see §5.3 |
| 35 | Receptor profile | same dir, `07_receptor_profile.png` | [§3.6](#36-appendix-chemical-space-figures-figures-27-to-36) | byte-identical to FIG_DIR, see §5.3 |
| 36 | Conformer-generation difficulty | same dir, `08_startconf_rmsd.png` | [§3.6](#36-appendix-chemical-space-figures-figures-27-to-36) | byte-identical to FIG_DIR, see §5.3 |

---

## 2b. Table index

> **Extended 2026-09-02.** This index audited four rows. Tables 1, 2, 5 and 6 are
> now recomputed and asserted on every run by `thesis_assertions.py`, driven by
> `thesis_expected_values.yaml`, which records each value's location in the thesis.

Only the tables whose source could be confirmed by re-deriving their printed values are
listed. The remainder were not audited in this pass and are marked as such rather than
assigned a guessed source.

Table numbers below are the THESIS numbers (corrected 2026-08-20; this index previously ran one to two behind).

| Table | Subject | Source file | Produced by | Provenance |
|---|---|---|---|---|
| 6 | PoseBusters pass-all recovery by ranking depth | `…/pose_comparison_report/topk_recovery_validity_gnina_arm.{csv,txt}` | [§5.1](#51-figures-2-3-and-4-no-committed-generator-for-the-published-arm) | value-checked |
| 7 | Decomposition of the rank-1 to best-of-top-15 gain | same file | [§5.1](#51-figures-2-3-and-4-no-committed-generator-for-the-published-arm) | value-checked |
| 8 | Near-nativeness and form recovery by ranking depth | same file | [§5.1](#51-figures-2-3-and-4-no-committed-generator-for-the-published-arm) | value-checked, and see the warning in [§5.2](#52-topn_within_thresholds_pbvalid_depthscsv-is-the-wrong-arm-do-not-regenerate-table-8-from-it) |
| 11 | Orai1 validity and placement yield | `posebusters_results/orai_pbvalid_tm_share_compare/pbvalid_tm_share_pooled.csv` and `pbvalid_tm_share_per_unit.csv` | [§3.4](#34-orai-validity-and-transmembrane-share-compare-figures-15-and-16) | value-checked. The pooled rows 120 / 119 / 99.2% / 66 / 55.0%, 120 / 112 / 93.3% / 82 / 68.3% and 120 / 29 / 24.2% / 1 / 0.8% reproduce the printed totals row exactly. **Re-verified 2026-08-20.** The AutoDock cell was 84 / 70.0% until the `optimized_rank` fix to `_CARRY` in `orai_transmembrane_exclusion.py`; the CSV was regenerated 2026-08-19 11:46 and now gives 66 / 55.0%, which is what the thesis prints. Do not restore the old pair |
| 1 to 5, 9, 10, 12 to 22 | — | — | — | not audited in this pass |

---

## 3. Commands with no notebook cell

> **Mostly closed 2026-09-02.** `Thesis_Reproduction.ipynb` now carries every
> command below as a declared stage, which is why this section existed. One entry
> is new rather than moved: the Orai PoseBusters-valid yield figure had no
> generator in any form and lived as an inline notebook cell. It is now
> `orai_pbvalid_yield_compare.py` and reproduces the shipped asset byte-identically.
>
> The interaction pose-basis audit in the 2026-08-29 block needs one more argument
> than it shows: `--ids-file` with the 303-id analysed cohort. The PandaMap arm
> profiles 306 complexes, so an unrestricted run aborts on a partial join. The
> cohort file is `analysed_cohort_ids.txt`, beside `per_pose_metrics.csv`.

These are the outputs the thesis depends on that no notebook cell produces. Each was
reconstructed from the script defaults plus the stored output's own provenance
header, and each reconstruction was checked by re-deriving the numbers the output
contains.

### 3.1 Whole-protein cost chapter, gnina arm (Figures 25 and 26)

Output directory `posebusters_results/benchmark_full_protein_vina_scoring/docking_effort_gnina/`.

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Analysis/docking_effort_comparison.py \
    --dataset benchmark \
    --autodock-dir Dockings/vina_results_full_protein_vina_scoring \
    --autodock-prep mgl_tools \
    --autodock-refine gnina --autodock-gnina-gpu \
    --unidock2-dir "" --unidock-dir "" \
    --per-pose-csv posebusters_results/benchmark_full_protein_vina_scoring/dock/pose_comparison_report/per_pose_metrics.csv \
    --out-dir posebusters_results/benchmark_full_protein_vina_scoring/docking_effort_gnina
```

`--diffdock-refine smina` and `--eq-refine gnina` are the `benchmark` dataset defaults
and were left unset, which is why the committed `effort_summary.csv` labels the arms
`DiffDock (smina)` and `EquiBind (unguided+gnina)`.

**TRAP 2026-08-20: the two empty-string flags are new and are load-bearing.** The cost
chapter compares THREE tools, and the original run produced three because Uni-Dock2 had
no whole-protein results yet — that campaign landed 5-10 August, after this run. The
`benchmark` dataset default for `--unidock2-dir` is
`Dockings/unidock2_results_full_protein_vina_scoring`, so re-running the command as it
previously stood now silently admits Uni-Dock2 as a fourth arm. That is not a cosmetic
extra bar. It changes the paired statistics underneath the figures, taking the wall
omnibus from df 2 to df 3, Kendall's W from 0.835 to 0.648 and the validity Cochran Q
from 102.5 to 154.8, none of which the chapter's text describes. Passing an empty string
is falsy in the reader, which drops the tool and reproduces the published three-arm
basis. `--unidock-dir ""` does the same for tiled Uni-Dock 1, which is currently inert
(no validated timing) but would behave the same way once it has any.

How the four non-default flags were recovered.

* `--autodock-refine gnina --autodock-gnina-gpu` follow from the committed summary row
  `AutoDock Vina + gnina` with device `CPU+GPU`. The arithmetic closes exactly. Raw
  Vina wall is 2.2065 h and the gnina pass is 12.1334 h, summing to the recorded
  14.3398 h, while CPU-core hours stay at 70.6067 h = 2.2065 × 32 threads because the
  whole gnina pass was charged to the GPU.
* `--autodock-dir` and `--autodock-prep mgl_tools` match the sibling `docking_effort`
  cell (`Master_Docking_AD_Full_Protein.ipynb` cell `4821046f`) and the on-disk layout
  `Dockings/vina_results_full_protein_vina_scoring/<id>/mgl_tools/`.
* `--per-pose-csv` is the flag that distinguishes this run from its sibling. Under the
  308-id allowlist the whole-protein per-pose table gives `autodock_gnina` 9,043 poses
  of which 8,995 are PoseBusters-valid, `diffdock_smina` 8,993 of which 7,597, and
  `equibind_unguided_gnina` 9,074 of which 4,702. The committed summary carries 9,013 /
  8,965, 8,963 / 7,567 and 9,060 / 4,702, that is exactly 30 poses fewer per method,
  which is the three complexes dropped to reach the recorded common timed set of 302.
  The dataset-default per-pose table (`posebusters_results/benchmark/dock/…`) gives
  6,956, 6,979 and 4,734 valid poses instead and reproduces only the older
  `docking_effort/` directory.

Beware that the sibling directory `docking_effort/` was produced by the notebook cell
without `--per-pose-csv`, so it pairs whole-protein AutoDock timing with pose counts
from a different run. It is a hybrid and is not the arm the cost chapter reports.

### 3.1b Charged-basis cost chapter, matched-EquiBind arm (Figures 25 and 26)

**This supersedes §3.1 as the source of the reported cost figures.** §3.1 documents the
`benchmark_full_protein_vina_scoring` tree, which pairs the exh128 AutoDock timings with the
**retired** `--local_only` EquiBind arm (62 qualifying complexes, 8.51 s median). The thesis
reports the `--minimize` EquiBind arm (80 complexes), so both cost figures and Table 6 come
from `benchmark_matched_equibind` instead. Two directories exist there:

| directory | AutoDock timing | what it is |
|---|---|---|
| `docking_effort_gnina` | 9,043 poses, 217 near2 | **mispaired** — an ADFRsuite tree joined to the derived Meeko method label. Do not use. |
| `docking_effort_gnina_v2` | 8,976 poses, 311 near2 | correct pairing, elapsed basis |
| `docking_effort_gnina_v2_charged` | 8,976 poses, 311 near2 | correct pairing, **charged basis — the reported arm** |

Both v2 runs use the same command; only `--basis` differs.

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Analysis/docking_effort_comparison.py \
    --dataset benchmark \
    --autodock-dir Dockings/vina_results_full_protein_vina_scoring_mgltools_exh128 \
    --autodock-prep mgl_tools --autodock-refine gnina --autodock-gnina-gpu \
    --autodock-method autodock_mgltools_exh128_gnina \
    --autodock-optimizer-workers 16 \
    --equibind-dir Dockings/Benchmark_Equibind_cputimed \
    --unidock2-dir "" --unidock-dir "" \
    --per-pose-csv posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report/per_pose_metrics.csv \
    --out-dir posebusters_results/benchmark_matched_equibind/docking_effort_gnina_v2_charged \
    --basis charged --cpu-threads 32
```

Drop `--basis charged --cpu-threads 32` and change the `--out-dir` to `..._v2` to rebuild the
elapsed tree.

**`--basis` (added 2026-09-01).** `elapsed` divides `wall_s` as each reader built it, which is
**not one quantity across arms**: AutoDock's is vina elapsed plus its gnina sum, DiffDock's is
measured, EquiBind's is best-config compute over its own `n_parallel_workers`. `charged` divides
`cpu_core_s / --cpu-threads + gpu_s` for every arm alike. Default is `elapsed`, and under it every
previously committed effort directory reproduces byte-for-byte apart from two new provenance keys
(`cost_basis`, `cpu_threads`) in the two stats JSONs — verified 2026-09-01, 21 of 23 files
md5-identical and the numeric payload of both JSONs unchanged.

Because neither `cpu_core_s` nor `gpu_s` is ever divided by `--autodock-optimizer-workers`, the
**charged basis is invariant to that flag**; the elapsed basis is not. That is why the same command
serves both.

**Why the basis was unified.** Table 6 previously spliced an AutoDock row charged at device
occupancy onto DiffDock and EquiBind rows charged at elapsed time, and Figure 25 was a third thing
again — a uniform elapsed panel whose AutoDock-versus-DiffDock bracket read `ns` directly beneath
prose claiming every pair cleared Holm. On the uniform charged basis the panel and the text agree:
Friedman chi2 = 87.1, p = 1.2e-19, Kendall W = 0.81, n = 54 paired, all three pairs `***`.
The AutoDock cells are unchanged by the unification (151.42 CPU-core-h / 32 = 4.73 h, its measured
search wall); only DiffDock (50.5 -> 49.4 s) and EquiBind (1.7 -> 2.4 s) move.

**Retired.** `posebusters_results/benchmark_full_protein_vina_scoring/docking_effort_gnina_serial/`
was the undocumented source of the old Table 6 AutoDock row. No command, notebook or guide entry
ever referenced it. It is superseded by `docking_effort_gnina_v2_charged` and should not be used.

### 3.2 Whole-protein PandaMap interaction fingerprints (Figures 9, 10 and 11)

Output directory `pandamap_results/benchmark_full_protein/`. The two steps are recorded
in the header of `Scripts/Analysis/pandamap_benchmark_full_protein_config.yaml`, which
also pins every variant used.

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Analysis/run_pandamap.py \
    --config Scripts/Analysis/pandamap_benchmark_full_protein_config.yaml

/home/manndo/anaconda3/envs/vina/bin/python Scripts/Analysis/pandamap_interaction_report.py \
    --config Scripts/Analysis/pandamap_benchmark_full_protein_config.yaml \
    --in-dir pandamap_results/benchmark_full_protein \
    --exclude-methods unidock2
```

The config header shows the report step invoked with `--in-dir` alone. Passing
`--config` as well is preferred and does not change the result, because the report
reads its variant pins, `benchmark_dir`, `poses_per_combo`, `oracle_summary` and
`per_pose_metrics` from that same file, and any of those left to the CLI defaults would
silently fall back to a different, superseded run.

Two settings in the config are load-bearing and must not be dropped. `poses_per_combo:
5` is what the top-5 fingerprint figure reports, and `canonical_receptor_dir:
posebusters_results/benchmark/dock/hetatm_cleaned` is what stops the renumbered
whole-protein receptors from driving measured native recovery to near zero.

The notebook PandaMap cells (`Master_Docking_AD_Full_Protein.ipynb` cells `3c3a757d`
and `68f54daf`) use `Scripts/Analysis/pandamap_config.yaml`, which points at a
superseded run. They do not produce the figures in the Results chapter.

### 3.3 Form versus placement filmstrip (Figure 5)

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Analysis/filmstrip_rank1_top5_top15.py
```

The script takes no arguments. Every choice is pinned in the file itself, namely the
report directory `posebusters_results/benchmark_full_protein_vina_scoring/dock/pose_comparison_report`,
the depths 1 / 5 / 15, the 5 Å axis clip, and the collapse to `autodock_gnina`,
`diffdock_smina` and `equibind_unguided_gnina`. The thesis figure is the `__thesis`
stem. The script also writes a `__compact_titles` variant and the per-pose CSV behind
every plotted point.

### 3.4 Orai validity and transmembrane share compare (Figures 15 and 16)

Output directory `posebusters_results/orai_pbvalid_tm_share_compare/`.

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Analysis/orai_pbvalid_tm_share_compare.py \
    --autodock-variant gnina --diffdock-variant smina --equibind-variant gnina \
    --top-n-poses 10
```

All four values are the script defaults except `--top-n-poses`, which defaults to 0.
They are stated explicitly because the Orai classification tables now carry both the
raw and the gnina AutoDock arm and a silent default is not a record.

Recovered from the committed sidecar as follows. `pbvalid_tm_share_compare_stats.txt`
records `POSE CAP: kept only each tool's top-10 ranked poses per (frame × ligand)`,
which fixes `--top-n-poses 10`. It reports 12,010 produced AutoDock poses on the control
panel, and the classification table splits its 24,160 AutoDock rows into 12,010 with
`variant == gnina` and 12,150 raw, so the gnina arm is the one selected. On the
experimental panel the same variant gives 360 rows over 12 frame-ligand units, capped
to the recorded 120. DiffDock 11,673 and EquiBind 12,150 match the `smina` and `gnina`
selections.

This directory also holds the ligand-level sensitivity analysis, which the Results and
Discussion chapters cite repeatedly but which has no figure or table number. It must be
run after the command above, because it reads that run's per-unit CSV.

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Analysis/orai_ligand_level_contrasts.py \
    --autodock-variant gnina --diffdock-variant smina --equibind-variant gnina
```

Everything else is a default, and the defaults are the shipped run. The generated
sidecar `orai_ligand_level_contrasts.txt` names its own generator and lists every input
CSV in a `SOURCES` block, so the reconstruction is self-checking.

### 3.5 Orai PoseBusters-valid yield compare (Figure 13)

Output directory `posebusters_results/orai_pbvalid_yield_compare/`.

This one does have a cell. It is `Master_Docking_AD_Full_Protein.ipynb` cell id
`ee79ca6e` (index 101), duplicated as cell index 91 in `Master_Docking.ipynb`. The cell
is self-contained inline code rather than a `subprocess.run` of a script, which is why
a plain grep for a script name does not find it.

Its pins are `AUTODOCK_KEY = "autodock_gnina"`, `DIFFDOCK_KEY = "diffdock_smina"`, a
top-10 pose cap through `pose_topn.top_n_allowlist`, and `EXCLUDE_LIGANDS =
{"gsk7975a-deprot-OPT"}`.

**Corrected 2026-08-16.** An earlier draft of this guide recorded the AutoDock arm here
as raw Vina and reported an arm asymmetry against the transmembrane compare in §3.4.
That was wrong. The stored output is the gnina arm and agrees with the thesis: the
sidecar `09f_pbvalid_yield_dominant_stats.txt` reads `AutoDock Vina Benchmark(gnina)
... Experimental(gnina)`, and every row of `pbvalid_yield_dominant_per_complex.csv`
carries `method_key = autodock_gnina`, 1,201 units on the benchmark panel and 12 on the
experimental one. Both Orai panels are therefore on the same gnina arm, exactly as
`body_main.tex` states. There is a real arm asymmetry among the Orai figures, but it is
in Figures 17, 18 and 19, not here. §5.5 records it.

The confusion came from a genuine defect, now fixed. The two notebooks carried the
SAME cell id `ee79ca6e` with DIFFERENT pins. `Master_Docking.ipynb` had the correct
`autodock_gnina` with an explanatory comment, while `Master_Docking_AD_Full_Protein.ipynb`
still had the older bare `autodock`. The committed figure came from the former.

**The repair is complete as of 2026-08-16 and was re-checked for this entry.** The two
copies of `ee79ca6e` are now byte-identical, 14,832 characters with MD5
`353b36cbfff763323c4781ebc46efad4` in both notebooks (cell index 101 in
`Master_Docking_AD_Full_Protein.ipynb`, index 91 in `Master_Docking.ipynb`). The
successor copy carries all three pieces the pin needs, namely
`AUTODOCK_KEY, DIFFDOCK_KEY = "autodock_gnina", "diffdock_smina"` plus the matching
lookup entries `TOOL_OF["autodock_gnina"] = "autodock"` and
`_VLABEL["autodock_gnina"] = "gnina"`. Both dictionaries also keep their bare
`"autodock"` entries, so either arm can be selected. Re-running either notebook now
reproduces the published figure.

Note for anyone bisecting the history: an intermediate state existed in which the pin
had been changed to `autodock_gnina` while `TOOL_OF` and `_VLABEL` still held only the
bare `"autodock"` key. That state raises `KeyError: 'autodock_gnina'` at plot time
rather than drawing the wrong arm. A cell that fails this way is mid-repair, not
broken by the fix.

**Watch this class of defect.** 21 of the 104 cell ids shared between the two notebooks
have diverged. Most divergences are deliberate, because `Master_Docking_AD_Full_Protein.ipynb`
is the whole-protein successor and points at different result directories. Before
trusting either notebook as the producer of a given output, diff the specific cell
against its twin rather than assuming the two agree.

### 3.6 Appendix chemical-space figures (Figures 27 to 36)

Open `PoseBusters_DataSet_Analysis.ipynb` on the `vina` kernel and run all cells. It
writes ten PNGs to `PoseBusters_Benchmark_Analysis/figures/` through its own `savefig`
helper. See §5.3 for why those ten files are not the ten figures in the appendix.

### 3.7 DiffDock re-ranking counterfactual (appendix top-k table note, and the Methods claim that smina does not re-rank)

Supplies the 32.3 / 32.7 / 33.7 / 35.6 percentages quoted in the note to the top-k recovery
table and in the DiffDock variant paragraph of the Methods. Reads the committed optimiser
logs plus a per-pose metrics table. Re-docks nothing.

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Analysis/diffdock_gnina_rerank_analysis.py \
    --tool smina \
    --per-pose-metrics posebusters_results/benchmark_full_protein_vina_scoring/dock/pose_comparison_report/per_pose_metrics.csv \
    --out-dir PoseBusters_Benchmark_Analysis/smina_rerank
# gnina arm (already committed under gnina_rerank/):
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Analysis/diffdock_gnina_rerank_analysis.py \
    --tool gnina \
    --per-pose-metrics posebusters_results/benchmark_full_protein_vina_scoring/dock/pose_comparison_report/per_pose_metrics.csv \
    --out-dir PoseBusters_Benchmark_Analysis/gnina_rerank
```

Read `selection_strategy_summary_<tool>.csv`. The four strategies are the 2x2 of
*which pose you pick* against *whose coordinates you keep*:

| strategy | meaning | smina | gnina |
|---|---|---|---|
| `A_native` | DiffDock confidence rank-1, raw coordinates | 33.7% | 33.7% |
| `C_rerank` | **re-ranking alone** — affinity rank-1, raw coordinates | **32.7%** | **32.3%** |
| `D_minimize` | confidence rank-1, minimised coordinates | 35.3% | 35.6% |
| `B_full` | affinity rank-1, minimised coordinates | 38.3% | 37.6% |

`A_native` and `D_minimize` reproduce the appendix table's own `DiffDock (raw)` and
`DiffDock (gnina-opt)` k = 1 rows exactly (33.7 and 35.6), which is what anchors the new
`C_rerank` cell to that table. **Re-ranking alone is slightly worse than DiffDock's own
confidence order**, so the benefit of refinement is the coordinate move, not the score —
this is the evidence behind §4.1 of `DOCKING_PROTOCOL.md`.

Two traps. The `--per-pose-metrics` default still points at the pre-whole-protein
`posebusters_results/benchmark/...` table; pass the `benchmark_full_protein_vina_scoring`
path explicitly as above. It happens not to matter (DiffDock is blind, so its poses and
RMSDs are identical in both tables and the gnina numbers come out bit-identical either
way), but do not rely on that silently. Second, part A of the script reports 423 complexes
(every complex with a successful minimisation) while part B reports 303 (those with a
crystal reference); the thesis quotes only the 303-complex selection numbers.

---

## 4. Ordering

Nothing in §3 re-docks anything. Every step reads committed CSVs, so the order matters
only where one step consumes another's output.

1. **Docking**, then **PoseBusters** (`Scripts/Docking/`, `run_posebusters.py`),
   producing `posebusters_results/<dataset>/dock/posebusters_filtered_results.csv`.
2. **Pose comparison** for the whole-protein benchmark
   (`Master_Docking_AD_Full_Protein.ipynb` cell `24f6dd5d`), producing
   `per_pose_metrics.csv` and `oracle_summary.csv` in
   `benchmark_full_protein_vina_scoring/dock/pose_comparison_report/`. Required by
   §3.1, §3.2, §3.3 and by the Table 6 to 8 sidecar in §5.1.
3. **Transmembrane exclusion** for both Orai panels
   (`Master_Docking_AD_Full_Protein.ipynb` cell `tm-both-groups-compare-code`), producing
   `posebusters_results/orai_{jku,benchmark}/transmembrane_filter/tm_pose_classification.csv`.
   Required by §3.4.
4. **Orai pose clustering** (cells `d2a207d2`, `3a10a8ca` and `577d6b41`, the last of
   which carries `--cluster-quality`) and **Orai PandaMap
   compare** (cell `60097bb6`), producing `pose_clusters/cluster_quality_per_tool.csv`,
   `pose_clusters/per_pair.csv` and `pandamap_results/orai_interaction_compare/pandamap_pose_totals.csv`.
   Required by the ligand-level contrasts in §3.4 and by Figures 17 to 19, which read
   the two `per_pair.csv` files rather than any pose table. See the input warning in
   [§5.7](#57-the-two-orai-xtool-compare-code-cells-read-different-inputs).
5. §3.1, §3.2, §3.3 and §3.5 are independent of one another and can run in any order
   once step 2 is done.
6. §3.4 in two parts. `orai_pbvalid_tm_share_compare.py` first, because
   `orai_ligand_level_contrasts.py` reads its `pbvalid_tm_share_per_unit.csv`. The
   contrasts script is the last thing to run in the whole chain, since it also consumes
   steps 3 and 4.

---

## 5. Known gaps and traps

> **Partly closed 2026-09-02.** Four items here were resolved by regenerating the
> missing outputs on the canonical tree, and three thesis figures that shipped from
> superseded trees were refreshed. What was closed, what moved and what remains is
> in [`FINDINGS_2026-09-02.md`](FINDINGS_2026-09-02.md), sections C, D and E.
> The hand-composed PyMOL renders in 5.4 remain unresolved and always will be.

### 5.1 Figures 2, 3 and 4, pin the AutoDock arm or you rebuild raw Vina

RESOLVED 2026-08-17 for Figures 3 and 4, whose problem was a missing flag rather than a
missing script. **Correction 2026-08-18:** this section used to say Figure 2 "always
rebuilt from the notebook cell of §3.5". That is wrong. §3.5 is the Orai yield figure,
and no committed cell rebuilds Figure 2 at all. See [§5.8](#58-figure-2-is-an-orphaned-render-values-reproduce-panel-set-does-not).

`posebusters_pose_comparison.py` collapses each tool family to one presentation variant
before it draws the headline figures. DiffDock had `--collapse-diffdock-variant` to pin
which one. AutoDock had no counterpart, so its slot always fell back to **raw Vina**
while the Results chapter reports the **gnina** arm. That is the whole discrepancy.

`--collapse-autodock-variant` now supplies the missing pin. It mirrors the DiffDock
flag, relabels the collapsed slot (`AutoDock Vina + gnina`) and leaves Vinardo alone.
Unset, it keeps the old raw-Vina behaviour, so every earlier command still reproduces
what it used to.

```bash
# NOTE 2026-08-20: this block previously read "--report-dir", which the script does not
# accept (argparse exits 2). The real flags are --pb-csv and --out-dir.
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Analysis/posebusters_pose_comparison.py \
    --pb-csv  posebusters_results/benchmark_full_protein_vina_scoring/dock/posebusters_filtered_results.csv \
    --out-dir posebusters_results/benchmark_full_protein_vina_scoring/dock/pose_comparison_report \
    --reuse-cache --collapse-plots-only \
    --collapse-diffdock-variant diffdock_smina \
    --collapse-autodock-variant autodock_gnina
```

**TRAP 2026-08-20: that command re-renders, it does NOT rebuild.** `--reuse-cache` loads
`per_pose_metrics.csv` and skips pair scoring entirely, so the two flags that shaped the
cached table are absent from it. Drop `--reuse-cache` after changing the input CSV and the
defaults take over, silently rebuilding the cache at `--top-n 5` with the raw
`--diffdock-variant diffdock` instead of the published `15` / `all`. Nothing errors; the
figures just quietly rest on a different basis. `per_pose_metrics.manifest.json` is what
catches it, since it records both values in its signature.

To REBUILD the cache (after the pb CSV changes), the two flags must be restored:

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Analysis/posebusters_pose_comparison.py \
    --pb-csv  posebusters_results/benchmark_full_protein_vina_scoring/dock/posebusters_filtered_results.csv \
    --out-dir posebusters_results/benchmark_full_protein_vina_scoring/dock/pose_comparison_report \
    --top-n 15 --diffdock-variant all \
    --collapse-plots-only \
    --collapse-diffdock-variant diffdock_smina \
    --collapse-autodock-variant autodock_gnina \
    --workers 24
```

`split_equibind: true`, `best_diffdock_only: true` and `limit_pairs: 0` need no flags:
the first is the default and the second is forced by `--collapse-plots-only`. Check the
rebuilt manifest against the committed signature — every field except the `pb_csv`
size/sha256 and `n_rows` must be unchanged. A full rebuild takes roughly 15 minutes at
`--workers 24`.

**Label change 2026-08-20.** The pinned AutoDock slot is now drawn as `AutoDock*`, matching
the `DiffDock*` / `EquiBind*` convention, so a collapsed legend names all three engines the
same way. Previously it read `AutoDock Vina + gnina`, leaving AutoDock as the only unstarred
entry. The change is in the `_ad_label` map of `posebusters_pose_comparison.py`. Re-running the
command above rewrites `18_topn_within_thresholds{,_kabsch}_pbvalid_depths.{png,csv}`; the CSVs
come back byte-identical to the `_gnina_arm` sidecars, which is the regression check that only
the legend moved. The gnina-refinement arm keeps its explicit name and is never starred.

Verified: that run reproduces Table 8 exactly, AutoDock 29.37 / 49.83 / 49.83, DiffDock
33.99 / 55.12 / 55.78, EquiBind 15.18 / 20.13 / 20.46 at 2 Å across depths 1 / 15 / 30,
and the legend reads `AutoDock Vina + gnina` as the published figures do. Dropping the
new flag reproduces the committed raw-arm CSV byte for byte (md5
`8a4ec9c43d47a38c7c36f44dd9229644`), which is the regression check for the change.

The gnina-arm outputs are committed alongside the raw ones as `…_gnina_arm.{png,csv}`
sidecars so both arms stay inspectable from one report directory. `thesis_latex/media/media/image3.png`
and `image4.png` are the gnina pair.

`Scripts/Analysis/topk_recovery_gnina_arm.py` remains the CSV/text source for the
numbers behind Tables 5, 6 and 7. It draws no figure, and it is no longer the only way
to reach the gnina arm.

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Analysis/topk_recovery_gnina_arm.py
```

### 5.2 `topn_within_thresholds_pbvalid_depths.csv` is the WRONG ARM. Do not regenerate Table 8 from it.

File:
`posebusters_results/benchmark_full_protein_vina_scoring/dock/pose_comparison_report/topn_within_thresholds_pbvalid_depths.csv`

It contains only `[autodock, autodock_vinardo, diffdock, equibind, unidock2]`, and at
the 2 Å threshold its AutoDock series across depths 1 / 15 / 30 reads

```
25.742574  47.854785  49.504950     (= 78 / 145 / 150 of 303 complexes)
```

That is **raw Vina**, the consequence of the hard-coded keep-set in §5.1. Table 8 and
Figures 3 and 4 report the gnina arm

```
29.4%  49.8%  49.8%                 (= 89 / 151 / 151 of 303 complexes)
```

The printed thesis values are the correct ones. They were verified against
`per_pose_metrics.csv` and are reproduced exactly by
`topk_recovery_validity_gnina_arm.txt`, which prints
`AutoDock gnina*   89 (29.4%)   151 (49.8%)   151 (49.8%)`. The thesis itself states
both arms side by side in the Results chapter, giving 25.7% and 47.9% for the
unoptimised search against 29.4% and 49.8% for the rescored one.

It is the committed sidecar CSV that is the raw arm, not the table. Nobody should
regenerate Table 8 from that CSV.

### 5.3 Appendix chemical-space figures — RESOLVED 2026-08-18

`PoseBusters_DataSet_Analysis.ipynb` writes ten PNGs to
`PoseBusters_Benchmark_Analysis/figures/`. Until 2026-08-18 those files were dated
2026-06-20 and none matched any thesis image, because the appendix carried the
notebook's **inline cell output** from a later re-run that was never written back to
`FIG_DIR`. Five of the ten reproduced from the stored notebook outputs (Figures 29, 30,
32, 35 and 36) and **five were unresolved** (Figures 27, 28, 31, 33 and 34). An earlier
revision of this section said four, and credited `pca-fit`/Figure 33 as matching; it did
not, in the working tree or at HEAD.

The gap is now closed. The notebook was re-executed end to end
(`jupyter nbconvert --to notebook --execute --inplace`, exit 0), which refreshed both
`FIG_DIR` and the stored inline outputs from the same run, and the ten published images
were then replaced with the **`FIG_DIR` files themselves** rather than the inline
outputs. `FIG_DIR` renders at `savefig.dpi = 150` against the inline `figure.dpi = 110`,
so the two were never going to be byte-identical, and publishing the saved artifact is
what makes the provenance exact. The published images are correspondingly 1.364x larger
in each dimension at identical aspect ratio and content.

Every appendix chemical-space figure is therefore now byte-identical to a stored `FIG_DIR` file:

| Thesis fig | `media/media/` | `PoseBusters_Benchmark_Analysis/figures/` | notebook cell |
|---|---|---|---|
| 27 | `image33.png` | `01_ligand_univariate_distributions.png` | `0f1d42fb` |
| 28 | `image34.png` | `02_druglikeness_ro5_qed.png` | `d31ff00e` |
| 29 | `image35.png` | `03_elemental_charge_composition.png` | `ea2dcd61` |
| 30 | `image36.png` | `04_ligand_correlation_heatmap.png` | `a22996c6` |
| 31 | `image37.png` | `05_2d_attribute_mappings.png` | `87264e87` |
| 32 | `image38.png` | `06_pairplot_core_descriptors.png` | `1b4b17d0` |
| 33 | `image39.png` | `09_pca_scree_loadings.png` | `pca-fit` |
| 34 | `image40.png` | `10_pca_biplot.png` | `pca-biplot` |
| 35 | `image41.png` | `07_receptor_profile.png` | `37f35e86` |
| 36 | `image42.png` | `08_startconf_rmsd.png` | `1ee29749` |

Note that the `FIG_DIR` numbering follows save order, not cell order, so 07/08 sit after
09/10 in the thesis sequence. Rebuild with

```bash
/home/manndo/anaconda3/envs/vina/bin/python -m jupyter nbconvert --to notebook \
    --execute --inplace --ExecutePreprocessor.timeout=1800 \
    PoseBusters_DataSet_Analysis.ipynb
cp PoseBusters_Benchmark_Analysis/figures/01_*.png thesis_latex/media/media/image33.png
# ...and so on per the table above.
```

### 5.8 Figure 2 — RESOLVED 2026-08-18, generator added

`thesis_latex/media/media/image2.png` used to match **no file anywhere in the repository
and no stored notebook inline output**. It came from a run that was never committed, so
the published copy could not be rebuilt byte for byte. §5.1 previously claimed Figure 2
"always rebuilt from the notebook cell of §3.5"; that was wrong twice, since §3.5 is the
Orai yield figure and no committed cell rebuilt Figure 2 at all.

The obstacle was the **panel set, not the data**. Figure 2 shows a curated nine series in
three families, while `posebusters_validity_report.py` draws every variant in the CSV,
which since then also gained AutoDock Vinardo and Uni-Dock2. The uncollapsed
`01_per_tool_validity.png` is correspondingly 3807 px wide with fourteen variants, and
the report script has no engine keep-set, only `--split-equibind`, `--diffdock-variant`
and `--equibind-variant`.

`Scripts/Analysis/figure2_validity_yield.py` supplies that keep-set. It scores the CSV
through `posebusters_validity_report.load_and_score`, so validity and the variant splits
mean exactly what they mean in the report, and it takes the panel as an explicit
`--series` spec so the nine are a stated choice rather than whatever the CSV holds.
Statistics follow the thesis Methods, namely a paired Wilcoxon signed-rank per optimised
arm against its family's raw arm, a Hodges-Lehmann estimate of the median paired
difference, and Holm correction across whichever comparisons the panel actually draws —
six on the old Meeko panel, five on the current one, because Holm's family is the set of
brackets present and shrinks with the family that `_prune_families` trims.

**UPDATED 2026-08-30, and the command below is the only correct one.** The block that
used to stand here took no exclusion flags and wrote into `validity_report/`. Both halves
are now wrong. It drew the superseded Meeko arm, and it put the result among the twelve
outputs of `posebusters_validity_report.py`, where a stale render is indistinguishable
from that report's own products. Figure 2's home is `validity_report_mgltools/`. See the
double-exclusion note under *Analysis changes of 2026-08-29 (thesis migration)* for why
the preset alone is not enough.

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Analysis/figure2_validity_yield.py \
    --csv posebusters_results/benchmark_full_protein_vina_scoring/dock/posebusters_filtered_results.csv \
    --out posebusters_results/benchmark_full_protein_vina_scoring/dock/validity_report_mgltools/00_figure2_validity_yield.png \
    --exclude-preset meeko \
    --exclude-methods autodock_mgltools,autodock_mgltools_gnina,autodock_mgltools_exh18,autodock_mgltools_exh64,autodock_mgltools_exh64_gnina,autodock_mgltools_exh92
cp posebusters_results/benchmark_full_protein_vina_scoring/dock/validity_report_mgltools/00_figure2_validity_yield.png \
   thesis_latex/media/media/image2.png
```

Current arm values, `autodock_mgltools_exh128` (8,976 produced poses, 303 complexes) and
its gnina rescore: medians 100.0 / 100.0 %, labels 303/8943 and 301/8880, means 99.6 and
98.9 %, and a single AutoDock bracket of +0.0 pp at p_holm = 0.0809, i.e. not significant.
Only five brackets are drawn, not six: `autodock_gnina_refinement` does not exist on this
arm, so `_prune_families` drops that member and the family keeps two boxes instead of
three. Re-running the command on 2026-08-30 reproduced the committed PNG **byte for byte**
(md5 `3ffacb8d5031bceaa8ad652ee082c8a4`), which is also the md5 of the shipped
`thesis_latex/media/media/image2.png`.

**The 2026-08-18 validation below is history and describes the Meeko arm.** Against the
superseded published render all 24 annotated values reproduced exactly: nine medians
(100.0, 100.0, 93.3, 13.3, 93.3, 96.6, 0.0, 26.7, 56.7 %), nine `complexes/valid-poses`
labels (303/8883, 303/8995, 303/8372, 242/2198, 300/7597, 300/7831, 35/263, 218/3565,
247/4702) and six Hodges-Lehmann deltas (+0.0, -5.0, +62.1, +65.0, +36.7, +48.3 pp). It
established that the generator is faithful. It is NOT a description of any figure the
thesis now ships, and those numbers must not be quoted as current.

A `.txt` sidecar carrying every plotted number is written beside the PNG, per the house
convention that inferential detail lives off the panel.

### 5.4 Hand-composed PyMOL renders (Figures 12 and 14)

Neither image exists anywhere in the repository. `Scripts/Analysis/pymol_pose_cluster_scene.py`
builds the per-frame scene files (`poses.pdb` plus `.pml`) that these views were
composed from, but the camera, colouring and final render were set interactively in
PyMOL and were not scripted. Both are therefore **UNRESOLVED** as automatic rebuilds.

### 5.5 The Orai AutoDock arm: all published figures are gnina

Every Orai per-pose table now carries BOTH the raw Vina and the gnina-rescored AutoDock
arm, on both panels. In `posebusters_results/orai_benchmark/dock/posebusters_filtered_results.no_tm.csv`
the AutoDock rows split as `optimizer == original` 10,611 and `optimizer == gnina`
10,477, and `pandamap_results/orai_benchmark/pandamap_pose_summary.csv` carries
`autodock` 6,699 poses next to `autodock_gnina` 3,406. Selecting an arm is therefore a
choice that every Orai run has to make, and the shipped runs did not all make the same
one.

**Every published Orai figure sits on the gnina arm.** Re-derived 2026-08-17 from the
committed artifacts themselves, by matching each thesis image to its source file by md5
and then reading the AutoDock pose provenance out of that source.

| Fig | Output | Exp. (JKU) panel | Orai × Benchmark panel | Evidence |
|---|---|---|---|---|
| 13 | `orai_pbvalid_yield_compare/09f_…` | gnina | gnina | inline cell `ee79ca6e`, `AUTODOCK_KEY = "autodock_gnina"` |
| 15, 16 | `orai_pbvalid_tm_share_compare/` | gnina | gnina | `--autodock-variant gnina` (§3.4) |
| 17, 18, 19 | `orai_jku/pose_clusters/panels/` | gnina | gnina | `per_pose.csv` AutoDock rows are 67/67 and 10,477/10,477 from `optimized_gnina/` |
| 20 | `orai_benchmark/pose_clusters/orai_cluster_quality_filtering.png` | not shown | gnina | same benchmark `per_pose.csv`, 10,477/10,477 gnina |
| 21 to 24 | `pandamap_results/orai_interaction_compare/` | gnina | gnina | stats sidecar line 7, `AutoDock Vina=autodock_gnina` |

There is no cross-panel arm mismatch in any published figure. An earlier revision of
this section claimed one for Figures 17 to 20 and for 21 to 24. That claim was derived
from the `--autodock-variant original` pins in the notebook cells rather than from the
outputs, and the outputs disagree with the pins.

**The stale pins were the error, not the figures.** Cells `577d6b41`, `3a10a8ca` and
`60097bb6` pinned `original` and justified it with `this arm carries only the raw Vina
poses`. That was true when written and stopped being true on 2026-08-15, when the
gnina arm was added to the Orai × Benchmark panel. The pins therefore did **not**
reproduce the committed figures. All three are now pinned to `gnina`, which does, and
their comments record the counts instead of the obsolete claim. The legacy
`Master_Docking.ipynb` twins of those cells carried no `--autodock-variant` at all,
which pools both arms and doubles every AutoDock count. They are pinned to `gnina` too.

**One cell could still emit a mixed comparison.** `tm-both-groups-compare-code` pinned
its Benchmark group to `original` and its JKU group to `gnina`. Its four combined
figures are not in the thesis and its numbers are not cited, so nothing published was
affected, and re-running both passes on `gnina` reproduces all six of its artifacts
byte for byte. The reason is that the combined figures are built in the last pass and
read both datasets under that pass's flag, which was always `gnina`. The Benchmark
group's own flag reached only its single-dataset outputs. Both groups are now pinned to
`gnina` so the cell cannot express a mixed comparison at all.

Note that `group_toolchain_tm_comparison.csv` in that folder is written by the notebook
cell, not by the script, so a script-only re-run leaves it untouched. Its AutoDock
`generated` count of 12,150 is the raw-arm figure. The gnina arm holds 12,010 poses.

**Which scripts can even express the choice.** Checked by
`grep -n "autodock-variant" Scripts/Analysis/orai_*.py`.

| Script | `--autodock-variant` | Default |
|---|---|---|
| `orai_transmembrane_exclusion.py` | yes | `all` |
| `orai_pose_cluster_report.py` | yes | `all` |
| `orai_pbvalid_tm_share_compare.py` | yes | `gnina` |
| `orai_ligand_level_contrasts.py` | yes | `gnina` |
| `orai_pandamap_interaction_compare.py` | yes | unset, takes whichever single variant each dataset carries |
| `orai_cross_tool_agreement_compare.py` | **no such flag** | arm is whatever the `--exp-csv` and `--bench-csv` runs used |
| `pymol_pose_cluster_scene.py` | no | — |

The yield compare has no script at all. It is the inline cell `ee79ca6e` of §3.5 and
sets its arm through the `AUTODOCK_KEY` constant.

State the flag explicitly wherever it exists, including where the default is already
right, because a default is not a record. Where it does not exist, namely
`orai_cross_tool_agreement_compare.py`, the only way to state the arm is to name the
`per_pair.csv` inputs and the runs that produced them.

### 5.6 Cost figures come from `docking_effort_gnina/`, not `docking_effort/`

See the closing note of §3.1. The `docking_effort/` directory pairs whole-protein
AutoDock timing with pose counts from a different run and is not what the cost chapter
reports.

### 5.7 The two `orai-xtool-compare-code` cells read DIFFERENT inputs

Figures 17, 18 and 19 are reproduced by the cell in **`Master_Docking.ipynb`** (index
68), which hard-codes

```python
EXP_CSV   = "posebusters_results/orai_jku/pose_clusters/per_pair.csv"
```

The same-id cell in `Master_Docking_AD_Full_Protein.ipynb` (index 77) does not. It reads

```python
_EXP_MATCHED  = "posebusters_results/orai_jku/pose_clusters_advina_matched/per_pair.csv"
_EXP_HEADLINE = "posebusters_results/orai_jku/pose_clusters/per_pair.csv"
EXP_CSV = _EXP_MATCHED if Path(_EXP_MATCHED).exists() else _EXP_HEADLINE
```

and the matched directory now exists, written by the second `subprocess.run` that the
successor copy of cell `d2a207d2` added. That cell re-runs the JKU clustering with
`--autodock-variant original` into `pose_clusters_advina_matched/`. Running the
successor notebook therefore silently swaps the Exp. input and overwrites all three
published PNGs in place, with no warning printed, because the warning is on the
fallback branch and the fallback is not taken.

The committed figures are the headline gnina run. Three independent checks agree.

* The sidecar `posebusters_results/orai_jku/pose_clusters/panels/orai_tool_agreement_compare_stats.txt`
  names `posebusters_results/orai_jku/pose_clusters/per_pair.csv` as its Exp. source.
* Its printed tool-pair medians are AutoDock ↔ DiffDock 30.1 Å over n = 11,
  AutoDock ↔ EquiBind 27.2 Å over n = 3 and DiffDock ↔ EquiBind 47.6 Å over n = 3. The
  headline `per_pair.csv` reproduces all three exactly. The matched one gives 33.0 Å
  over n = 12, 47.1 Å over n = 3 and 47.6 Å over n = 3.
* For Figure 19, `orai_cluster_quality_compare_stats.txt` prints an Exp. AutoDock
  silhouette of 0.67 and a Calinski-Harabasz of 169.33 over n = 11. The headline
  `cluster_quality_per_tool.csv` gives exactly that once the default
  `--exclude-ligand gsk7975a-deprot` is applied. The matched one gives 0.83 and 650.56
  over n = 12.

All three committed PNGs are byte-identical to the thesis images
`thesis_latex/media/media/image16.png`, `image17.png` and `image18.png`.

The command-line equivalent, which is the safer way to rebuild them, is

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Analysis/orai_cross_tool_agreement_compare.py \
    --exp-csv posebusters_results/orai_jku/pose_clusters/per_pair.csv \
    --exp-label "Exp. Ligands" \
    --bench-csv posebusters_results/orai_benchmark/pose_clusters/per_pair.csv \
    --bench-label "Benchmark ligands × Orai" \
    --out-dir posebusters_results/orai_jku/pose_clusters/panels
```

**Corrected 2026-08-18.** The strict all-three-tools-agree figure annotation and the
stats sidecar both used to assert a hard-coded `0 %` for both panels. That was never
measured. Computed over the units that actually carry all three tools it is
**0/3 for the experimental ligands and 1/467 = 0.2% for the benchmark control panel**,
the single agreeing unit being `7OZC_G6S` on frame `Orai1WT-MDSnap-Fr300`, which is what
the Results prose has always reported. `strict_all3()` in
`orai_cross_tool_agreement_compare.py` now derives it from `au_di_dist`, `au_eq_dist` and
`di_eq_dist`, counting only units where all three distances exist, so two-tool units are
excluded from the denominator rather than silently counted as disagreements. Figure 16
(`image17.png`) was regenerated from this fix.

`--bench-label` must be passed, since the script default is the shorter
`Benchmark ligands`. Both clustering runs must exist first, that is cell `d2a207d2` for
the JKU side and cell `577d6b41` for the benchmark side, which is step 4 of §4. The
`_poreblockers` suffix run that both cells fire afterwards writes a separate ligand
subset and is not a thesis figure.

The successor cell was written to remove a raw-versus-gnina asymmetry that an earlier
revision of §5.5 reported. That asymmetry does not exist. Both published panels are
gnina, so the matched run does not fix anything.

It is worse than redundant. `pose_clusters_advina_matched/` is the JKU side re-clustered
on `--autodock-variant original`, while the Benchmark side it is compared against is
gnina. Adopting the matched input would therefore put raw Vina on one panel and gnina
on the other, which is precisely the mismatch it was meant to remove. The headline
input is the correct one on the merits as well as being the one the thesis shows.

Do not let the successor cell take the `_MATCHED` branch. Either delete
`pose_clusters_advina_matched/` or rebuild the three PNGs with the command above, which
names its inputs explicitly and cannot silently swap them. The gnina cluster quality is
also the conservative reading: the Exp. AutoDock silhouette is 0.67 on the published
arm against 0.83 on the matched one.

### 5.13 7WCF_ACP was re-docked 2026-08-20 (whole-protein EquiBind only)

The 2026-07-10 EquiBind campaign was interrupted (SIGINT) during this complex's phase-3
post-processing, which truncated its pose set to 248 of 270 and lost its
`pipeline_summary.json`. The launcher's resume test only asks whether `__refGNINA` SDFs
exist, so the restart skipped it and the gap became permanent. A 2026-08-04 stub manifest
restored `uff_minimize` but not the clocks.

It was re-docked under `equibind_benchmark_config.yaml` (driver
`Scripts/Docking/rerun_equibind_7wcf_timing.py`, `--gnina-gpu false`), re-busted in
isolation, and its 270 rows spliced into the whole-protein
`posebusters_filtered_results.csv`. The truncated poses are preserved at
`Dockings/_7wcf_truncated_backup_20260820/`, the pre-splice CSV at
`posebusters_filtered_results.csv.bak-7wcf-repropagate-20260820`, and the calibration
evidence at `Dockings/Benchmark_Equibind_timing_rerun/PROVENANCE.md`.

Consequence: the whole-protein cost chapter now runs on 303 complexes rather than 302,
and EquiBind's `unguided+gnina` arm gains 16 poses (9,074 to 9,090), one valid pose
(4,702 to 4,703) and one complex with a valid pose (247 to 248).

The superseded report under `posebusters_results/benchmark/dock/` was deliberately NOT
re-propagated. It is a frozen artifact and holds 274 EquiBind rows for this complex,
matching neither the truncated 248 nor the complete 270 on disk, so it was already on
its own snapshot before this change. Nothing reported in the thesis is drawn from it.
