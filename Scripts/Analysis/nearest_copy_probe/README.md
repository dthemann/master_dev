# Nearest-copy sensitivity probe (examiner finding M2 / E3-2, 2026-09-07)

Re-scores every benchmark pose against every crystallographic copy in
`<ID>_ligands.sdf`, using the thesis's own loaders and `symmetry_rmsd`
(imported from `posebusters_pose_comparison.py`). The reference-copy value
reproduces the `rmsd` column of `per_pose_metrics.csv` to 1e-9 Å, so the
nearest-copy column is the only new quantity.

Files
- `nearest_copy_rescoring.py`     the probe; `--rank1-only` for the fast pass, `--ppm` to point at another per_pose_metrics.csv
- `nearest_copy_rank1.csv`        rank-1 poses, seven whole-protein arms (2,121 rows)
- `nearest_copy_boxed_rank1.csv`  rank-1 poses of the OLD boxed AutoDock run + DiffDock-smina (regime contrast)
- `nearest_copy_sensitivity_summary.csv`  per-arm, per-depth existence-gate recovery under both conventions

Regenerate (vina env, ~1 min on 16 cores for all depths):

    /home/manndo/anaconda3/envs/vina/bin/python Scripts/Analysis/nearest_copy_probe/nearest_copy_rescoring.py \
        --procs 16 --arms autodock_mgltools_exh128_gnina autodock_mgltools_exh128 diffdock_smina diffdock \
        equibind_unguided_gnina unidock2 --out /tmp/nearest_copy_alldepths.csv

Columns added per pose: n_copies, rmsd_ref_re (must equal rmsd), rmsd_nearest_copy,
nearest_copy_is_ref, cd_ref, cd_min_alt (centroid distance to the nearest alternate copy),
rmsd_min_alt. EquiBind's rank==999 sentinel is replaced by eff_rank from gnina_affinity.
See thesis_latex/M2_REVALIDATION_2026-09-07.md for the adjudication.
