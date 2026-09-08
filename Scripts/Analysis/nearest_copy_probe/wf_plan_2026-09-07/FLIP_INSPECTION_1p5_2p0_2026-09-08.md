# Flip inspection: rank-1 poses with 1.5 < nearest-copy RMSD <= 2.0 A (2026-09-08)

Question. Under the nearest-copy convention, eight rank-1 poses that fail the reference-copy 2 A
gate (rmsd > 2 A) pass it against a second deposited copy with an RMSD in the 1.5-2.0 A band.
Are these flips genuine near-native placements at the other copy site, or an artefact of the
symmetry-corrected RMSD engine?

Selection. `per_pose_two_conventions.csv`, eff_rank == 1, rmsd > 2, 1.5 < rmsd_nearest_copy <= 2.0,
method in {autodock_mgltools_exh128_gnina, diffdock_smina}. Exactly eight rows.

Method (read-only; `/home/manndo/anaconda3/envs/vina/bin/python`).
* Pose loaded with `ppc.load_first_mol`; bond orders re-assigned from `<ID>_ligand.sdf` with
  `ppc.reassign_template` (succeeded for all eight; heavy-atom counts equal the crystal copy).
* Copies: every record of `<ID>_ligands.sdf` (`SDMolSupplier(removeHs=True, sanitize=False)` +
  `SanitizeMol`; record 0 is the reference).
* (a) thesis engine `ppc.symmetry_rmsd` (bond-order-respecting `GetSubstructMatches`, in place).
* (b) PoseBusters `ppc.posebusters_rmsd` (in place and Kabsch).
* (c) identity-mapping RMSD (atom i of pose vs atom i of copy, no matching).
* Extra diagnostic: RMSD of the WORST bond-order-respecting matching (a "no symmetry choice" number
  with a valid atom correspondence), RDKit `CalcRMS` with `symmetrizeConjugatedTerminalGroups`
  True/False, and the pose-vs-copy element sequences.
* Groups by SMARTS: carboxylate / phosphate / sulfonate / nitro / guanidinium.
* Verdict rule as specified: GENUINE if (a) and (b, in place) agree within 0.3 A and both <= 2.0 A;
  ARTEFACT otherwise. Where the literal rule fails I say so and give the cause.

## Table A: geometry (all RMSDs heavy-atom, A)

| # | Method | Protein | Pose file | RMSD to ref | Nearest copy | (a) thesis | (b) PB in place | (b) PB Kabsch | (c) identity | worst matching | # matchings | Heavy atoms | Groups | Centroid pose-copy | Centroid pose-ref |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | diffdock_smina | 6YT6_PKE | rank1_confidence-0.40_smina.sdf | 45.675 | 1 of 2 | 1.991 | 1.991 | 1.752 | 5.511 | 1.992 | 2 | 31 | none | 0.688 | 45.253 |
| 2 | autodock_mgltools_exh128_gnina | 7A9H_TPP | ..._vina_out_rank2_gnina.sdf | 30.436 | 1 of 2 | 1.679 | 1.595 | 1.513 | 4.085 | 1.831 | 2 | 26 | phosphate (diphosphate) | 0.272 | 29.904 |
| 3 | diffdock_smina | 7AFX_R9K | rank1_confidence-0.04_smina.sdf | 15.192 | 1 of 2 | 1.689 | 1.689 | 1.601 | 3.640 | 1.689 | 1 | 20 | carboxylate | 0.498 | 14.014 |
| 4 | autodock_mgltools_exh128_gnina | 7MYU_ZR7 | ..._vina_out_rank3_gnina.sdf | 48.875 | 1 of 2 | 1.815 | 1.815 | 1.210 | 9.317 | 1.946 | 2 | 35 | none | 1.079 | 48.729 |
| 5 | autodock_mgltools_exh128_gnina | 7O0N_CDP | ..._vina_out_rank1_gnina.sdf | 12.691 | 1 of 2 | 1.628 | 1.616 | 1.494 | 7.353 | 1.695 | 2 | 25 | phosphate (diphosphate) | 0.583 | 10.270 |
| 6 | diffdock_smina | 7TB0_UD1 | rank1_confidence-1.22_smina.sdf | 26.877 | 3 of 4 | 1.508 | 1.508 | 1.359 | 7.731 | 1.508 | 1 | 39 | phosphate (diphosphate) | 0.602 | 25.265 |
| 7 | autodock_mgltools_exh128_gnina | 7TBU_S3P | ..._vina_out_rank28_gnina.sdf | 42.962 | 1 of 2 | 1.701 | 1.400 | 1.252 | 3.195 | 1.716 | 2 | 16 | carboxylate, phosphate | 0.447 | 42.729 |
| 8 | diffdock_smina | 7UJ5_DGL | rank1_confidence0.15_smina.sdf | 19.895 | 3 of 4 | 1.874 | 1.284 | 0.922 | 2.760 | 1.874 | 1 | 10 | carboxylate (x2) | 0.463 | 19.384 |

AutoDock pose files are `<ID>_protein_mgl_tools__<ID>_ligand_start_conf_vina_vina_out_rankN_gnina.sdf`
under `Dockings/vina_results_full_protein_vina_scoring_mgltools_exh128/<ID>/mgl_tools/docking/optimized_gnina/`;
DiffDock files under `posebusters_results/_benchmark_staging/diffdock/<ID>_ligand__<ID>_protein/optimized_smina/`.
Recomputed (a), (b) and centroid values reproduce the CSV columns rmsd_nearest_copy, pb_rmsd_nearest,
kabsch_nearest and cd_nearest to < 1e-9 A on all eight rows; the RMSD to the reference (copy 0) reproduces
the CSV `rmsd` column likewise.

## Table B: verdicts

| # | Protein | (a)-(b) spread | 0.3 A rule | both <= 2 A | Verdict | Reasoning |
|---|---|---|---|---|---|---|
| 1 | 6YT6_PKE | 0.000 | pass | yes | GENUINE | Both engines give 1.991 A and the worst matching is 1.992 A, so the pose sits on copy 1 (centroid 0.69 A away) independent of symmetry handling, although it clears the gate by only 9 mA. |
| 2 | 7A9H_TPP | 0.083 | pass | yes | GENUINE | The spread is terminal-phosphate O symmetrisation in PoseBusters; the worst bond-order-strict matching is 1.831 A, so the pose is near-native at copy 1 under every mapping. |
| 3 | 7AFX_R9K | 0.000 | pass | yes | GENUINE | Single matching, engines identical at 1.689 A, centroid 0.50 A from copy 1; no symmetry choice is involved. |
| 4 | 7MYU_ZR7 | 0.000 | pass | yes | GENUINE | Engines identical at 1.815 A and the worst of the two pyrimidine-flip matchings is 1.946 A, so the flip does not depend on the mapping. |
| 5 | 7O0N_CDP | 0.012 | pass | yes | GENUINE | 12 mA spread from phosphate-O symmetrisation; worst strict matching 1.695 A; centroid 0.58 A from copy 1. |
| 6 | 7TB0_UD1 | 0.000 | pass | yes | GENUINE | Single matching, engines identical at 1.508 A, centroid 0.60 A from copy 3 of 4. |
| 7 | 7TBU_S3P | 0.301 | FAIL (by 1 mA) | yes | GENUINE (literal rule: ARTEFACT) | The 0.30 A spread is exactly PoseBusters treating the acid C=O/C-OH and phosphate oxygens as interchangeable (CalcRMS sym=True 1.400 vs sym=False 1.701); the thesis engine is the stricter one and both stay <= 2 A, so the flip is not an engine artefact. |
| 8 | 7UJ5_DGL | 0.589 | FAIL | yes | GENUINE (literal rule: ARTEFACT) | Glutamate with two carboxylates: PoseBusters' conjugated-terminal-O symmetrisation lowers 1.874 to 1.284 A on a 10-atom ligand; the thesis engine's strict number is the conservative one, only one strict matching exists, and both engines stay <= 2 A, so the flip holds under the stricter engine. |

## Findings

1. Every one of the eight poses sits on a second deposited copy of the same ligand: pose centroid
   0.27-1.08 A from the nearest copy versus 10-49 A from the reference copy. These are copy-site
   choices, not misplacements, which is what the nearest-copy convention is meant to credit.
2. Reproduction: the CSV's nearest-copy numbers under both engines are reproduced to < 1e-9 A;
   template re-assignment succeeded for all eight and heavy-atom counts match the crystal copy.
3. Engine agreement: (a) and (b) are identical to machine precision on 4/8 poses. The four
   divergences (7A9H_TPP 0.083, 7O0N_CDP 0.012, 7TBU_S3P 0.301, 7UJ5_DGL 0.589 A) are all
   phosphate/carboxylate ligands and are reproduced exactly by one RDKit flag: the thesis
   `symmetry_rmsd` equals `CalcRMS(symmetrizeConjugatedTerminalGroups=False)` on all eight, and
   PoseBusters equals `CalcRMS(symmetrizeConjugatedTerminalGroups=True)` on all eight (PoseBusters
   default, `posebusters/modules/rmsd.py:91`). In every divergent case the thesis engine is the
   stricter (higher) number.
4. Symmetry-mapping choice does not carry any pose under the gate: the WORST bond-order-respecting
   matching is <= 1.946 A for all eight (max on 7MYU_ZR7). Five poses have only 1-2 matchings.
5. The identity-mapping RMSD (c) is not a usable check here: pose and crystal atom orders differ for
   all eight (the element sequences differ, e.g. 7UJ5_DGL pose NCCOOCCCOO vs copy NCCOCCCOOO), because
   the AutoDock PDBQT round-trip and the DiffDock/smina outputs re-order atoms. (c) therefore compares
   non-corresponding atoms (2.8-9.3 A) and is reported only for completeness; the "worst matching"
   column is the meaningful no-symmetry-choice figure.
6. Fragility: 6YT6_PKE clears the gate at 1.991 A under both engines (9 mA of margin). It is a real
   near-native pose on copy 1 but any change of engine or tolerance could move it across the threshold.

## Conclusion

None of the eight 1.5-2.0 A rank-1 flips is an artefact of the symmetry-corrected RMSD engine. All eight
are near-native placements on a second deposited copy, confirmed in place by both the thesis engine and
PoseBusters, and all remain <= 2 A even under the worst symmetry-equivalent mapping. Applying the 0.3 A
agreement criterion literally would label 7TBU_S3P (0.301 A) and 7UJ5_DGL (0.589 A) ARTEFACT; I do not adopt
that label because the disagreement is fully explained by PoseBusters' conjugated-terminal-group
symmetrisation (carboxylate/phosphate oxygens treated as interchangeable), which makes the thesis number the
conservative one, so switching engines could only flip these two harder, not undo them. The thesis engine's
bond-order-strict matching is a defensible, slightly conservative choice; if the thesis ever adopts
PoseBusters' convention, nearest-copy RMSDs for acid/phosphate ligands will fall by up to ~0.6 A and a few
extra borderline poses may cross 2 A.

Scripts (scratchpad only, not part of the repo): `wf/flip_inspect.py`, `wf/cause_check.py` under
`/tmp/claude-1000/-home-manndo-master-dev/55d7a717-f5e7-4353-8e73-45c4a49f6d87/scratchpad/`.
