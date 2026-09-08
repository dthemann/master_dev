# Nearest-copy endpoint program (independent layout)

Executes plan v2 (`thesis_latex/PLAN_nearest_copy_primary_endpoint_2026-09-07_v2.md`) Phases 2-8 without touching
any existing finding. Rule from the author (2026-09-08): no data from the thesis is overwritten. Every canonical
directory, the shipped thesis under `thesis_latex/`, the harness yaml and `thesis_assertions.py` are READ-ONLY.
Every output of this program lives in a NEW location carrying the suffix `_nearest`:

| Canonical (read-only) | Program output |
|---|---|
| `posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report/` | `.../pose_comparison_report_nearest/` (hub, `--reference-convention nearest --force`) |
| `.../dock/validity_report_mgltools/` | `.../dock/validity_report_mgltools_nearest/` |
| `posebusters_results/cluster_crystal_pocket_matched_equibind/..._allposes/` | `..._allposes_nearest/` |
| `pandamap_results/benchmark_matched_equibind/` and `report/` | `pandamap_results/benchmark_matched_equibind_nearest/` and `report/` (pose fingerprints copied in, crystal fingerprints per copy recomputed) |
| `posebusters_results/benchmark_matched_equibind/docking_effort_gnina_v2{,_charged}/` | `..._nearest` |
| `posebusters_results/autodock_exhaustiveness_returns/` | `posebusters_results/autodock_exhaustiveness_returns_nearest/` |
| `PoseBusters_Benchmark_Analysis/{ligand,receptor}_difficulty/`, `smina_rerank/` | `..._nearest` |
| (new) | `posebusters_results/metal_stratum_nearest/`, `itt_nearest/`, `reference_convention_nearest/` |
| `thesis_latex/` | `thesis_latex_nearest/` (independent copy of the short build; all text edits, figures and builds happen there) |
| `Scripts/Analysis/thesis_expected_values.yaml`, `thesis_assertions.py` | copies in this folder with the `_nearest` paths and `thesis_latex_nearest/` |

Contents of this folder (filled in as the program advances): `run_program.py` (every Phase 3 command with its
`_nearest` output path), `PROGRAM_LOG.md` (what ran, when, result), `harness_signatures/` (change signature
A / B / C), `thesis_expected_values_nearest.yaml`, `thesis_assertions_nearest.py`.

Phases: 2 downstream code (cluster any-copy rule, PandaMap per-copy crystal fingerprints, interaction-report
joins, exhaustiveness pin keying, ITT and reference-convention stages) → 3 runs into the `_nearest` directories
→ 4 selection re-check → 5 text edits in `thesis_latex_nearest/` → 6 figures into its `media/` → 7 harness copies →
8 build and QA of the copy.
