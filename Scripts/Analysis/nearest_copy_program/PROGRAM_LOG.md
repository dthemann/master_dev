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
