# Validation — Benchmark accuracy vs. validity table (rank-1 / first-pose / oracle)

**Date:** 2026-07-19
**Verdict:** NOT fully correct. 8 of 10 rows validate against the raw data; the two
`DiffDock* (gnina-opt)` rows carry the *smina*-optimized variant's accuracy values.

## Data source & method

All numbers recomputed independently from the raw per-pose table:

- `posebusters_results/benchmark/dock/pose_comparison_report/per_pose_metrics.csv` (114,814 poses, 14 variants)

Selection rules (replicating `Scripts/Analysis/posebusters_pose_comparison.py`):

- **oracle** — minimum-RMSD pose per (protein, ligand) complex.
- **rank-1** — pose with `rank == 1` (tie-broken by `pose_name`). Used for AutoDock, DiffDock, DiffDock*.
- **first pose** — lowest EquiBind generation index (`_<NNN>__ref`). Used for both EquiBind rows.
- **Accurate (≤2 Å)** — selected pose has `rmsd ≤ 2.0`.
- **Valid (≤2 Å & PB)** — the *same* selected pose has `rmsd ≤ 2.0` AND `pb_valid == True`.
- **Validity gap** — accuracy − validity (percentage points).
- **Total poses** — count of that variant's rows (`n_poses`), cross-checked against
  `oracle_summary_all_variants.csv`.

Cross-checked against report CSVs: `oracle_summary_all_variants.csv`,
`top1_summary_all_variants.csv`, `optimization_raw_vs_best.csv`.

## Per-cell result (submitted → recomputed)

| Variant | Poses | Submitted acc / val / gap | Recomputed acc / val / gap | OK? |
|---|---|---|---|---|
| AutoDock Vina — rank-1 | 7,407 | 56.4 / 53.1 / 3.3 | 56.44 / 53.14 / 3.30 | ✓ |
| AutoDock Vina — oracle | 7,407 | 89.4 / 85.5 / 4.0 | 89.44 / 85.48 / 3.96 | ✓ |
| DiffDock (raw) — rank-1 | 9,316 | 33.7 / 17.8 / 15.8 | 33.66 / 17.82 / 15.84 | ✓ |
| DiffDock (raw) — oracle | 9,316 | 50.2 / 23.4 / 26.7 | 50.17 / 23.43 / 26.74 | ✓ |
| DiffDock* (gnina-opt) — rank-1 | 8,993 | 35.3 / 32.3 / 3.0 | **35.64** / 32.34 / **3.30** | ✗ acc & gap |
| DiffDock* (gnina-opt) — oracle | 8,993 | 57.1 / 50.5 / 5.6 | **58.09** / 50.50 / **7.59** | ✗ acc & gap |
| EquiBind (raw) — first pose | 9,089 | 1.3 / 0.7 / 0.7 | 1.32 / 0.66 / 0.66 | ✓ |
| EquiBind (raw) — oracle | 9,089 | 6.3 / 0.7 / 5.6 | 6.27 / 0.66 / 5.61 | ✓ |
| EquiBind* (gnina-opt) — first pose | 9,074 | 3.3 / 2.6 / 0.7 | 3.30 / 2.64 / 0.66 | ✓ |
| EquiBind* (gnina-opt) — oracle | 9,074 | 21.5 / 19.5 / 2.0 | 21.45 / 19.47 / 1.98 | ✓ |

## Root cause of the 2 errors

The submitted `DiffDock*` accuracy figures are the **smina** variant's, not gnina's:

| | gnina-opt (correct for label) | smina-opt (submitted values) |
|---|---|---|
| rank-1 accuracy | 35.64% | **35.31% → 35.3** |
| oracle accuracy | 58.09% | **57.10% → 57.1** |
| smina oracle validity | — | 51.49% (note: submitted validity 50.5 is gnina's, not smina's) |

The rank-1 row (35.3 / 32.3 / 3.0) is a clean match to smina (35.31 / 32.34 / 2.97).
The oracle row is mixed: accuracy 57.1 and gap 5.6 are smina's, but validity 50.5 is gnina's
(smina oracle validity is 51.49). Since the label reads "(gnina-opt)" and the rest of the
study treats DiffDock* = gnina, the accuracy cells are what to fix.

## Validated / corrected table

| Variant | Total poses | Accurate (≤2 Å) | Valid (≤2 Å & PB) | Validity gap |
|---|---|---|---|---|
| AutoDock Vina — rank-1 | 7,407 | 56.4% | 53.1% | 3.3% |
| AutoDock Vina — oracle | 7,407 | 89.4% | 85.5% | 4.0% |
| DiffDock (raw) — rank-1 | 9,316 | 33.7% | 17.8% | 15.8% |
| DiffDock (raw) — oracle | 9,316 | 50.2% | 23.4% | 26.7% |
| DiffDock* (gnina-opt) — rank-1 | 8,993 | **35.6%** | 32.3% | **3.3%** |
| DiffDock* (gnina-opt) — oracle | 8,993 | **58.1%** | 50.5% | **7.6%** |
| EquiBind (raw) — first pose | 9,089 | 1.3% | 0.7% | 0.7% |
| EquiBind (raw) — oracle | 9,089 | 6.3% | 0.7% | 5.6% |
| EquiBind* (gnina-opt) — first pose | 9,074 | 3.3% | 2.6% | 0.7% |
| EquiBind* (gnina-opt) — oracle | 9,074 | 21.5% | 19.5% | 2.0% |

Only the four **bolded** cells changed from the submitted table. If `DiffDock*` was
meant to be smina rather than gnina, the alternative single-cell fix is the reverse:
keep 35.3 / 57.1 accuracy but change the oracle validity 50.5 → 51.5 (smina's value).
