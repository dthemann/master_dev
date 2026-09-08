# Nearest-copy program log

Rule: no data from the thesis is overwritten. All outputs in `_nearest` locations or under this folder;
`thesis_latex_nearest/` is the working thesis copy. See README.md for the layout.

| When (2026-09-08) | Step | Result |
|---|---|---|
| 10:45 | Author decision: run the full endpoint switch in an independent layout | program started |
| 10:50 | `thesis_latex_nearest/` created (sources, bib, class, media, PICs; build artefacts and review docs excluded) | `latexmk -pdf -f Thesis_short.tex` in the copy: exit 0, PDF 8.7 MB (the class declares dvips drivers but the canonical build is pdflatex too, per its `.fdb_latexmk`) |
| 10:52 | Hub rebuild under `--reference-convention nearest --force` into `pose_comparison_report_nearest/` (24 workers) | running |
| 10:53 | Phase 2 workflow launched: cluster any-copy rule, PandaMap per-copy fingerprints and report joins, exhaustiveness pin keying, ITT + reference-convention stages + Phase 3 driver; two reviewers per change | running |
| 10:57 | `make_harness_copies.py` → `thesis_expected_values_nearest.yaml` (39 path remaps, values still canonical) and `thesis_assertions_nearest.py` (SPEC and two media references remapped) | written; both parse |
| 10:58 | `selection_check.py` (Phase 4) written and exercised on the eight-id smoke table | mechanics verified; ties reported, not broken |
| 11:08 | Hub nearest table written: `pose_comparison_report_nearest/per_pose_metrics.csv` (241,913 rows, manifest `nearest`, schema 6); aggregation and figures still writing | validated against the canonical table: the four `*_ref_instance` twins equal the canonical `rmsd`, `pb_rmsd`, `bestfit_rmsd`, `centroid_dist` with max difference 0.0 on all rows; `rmsd <= rmsd_ref_instance` everywhere; 138 multi-copy complexes, 60,539 poses nearest to an alternate copy. PB-valid & <= 2 A existence gate, nearest vs single instance: AutoDock* 149 / 213 / 214 vs 111 / 198 / 199, DiffDock* 132 / 179 / 183 vs 103 / 167 / 169, EquiBind* 57 / 83 / 83 vs 55 / 80 / 80 at d = 1 / 15 / 30, exactly the probe previews of plan Section 1.1 |
