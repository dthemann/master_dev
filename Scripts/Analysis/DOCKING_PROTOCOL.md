# Docking Protocol — Full Methods Reference

*Auto-generated 2026-07-19 from the actual configs, runner scripts, and installed binaries
(`Scripts/Docking/*`, `Scripts/Utilities/prep_docking.py`, `equibind_pipeline/*`, `Master_Docking.ipynb`).
Every value below was read from source or obtained by running the binary; items that could not be
verified are flagged NOT FOUND / inferred.*

Two docking campaigns are the focus here:

1. **PoseBusters Benchmark set** — 428 protein–ligand complexes, **self-docking / re-docking** (each ligand
   docked back into its own crystal protein; a native pose exists → RMSD-to-native available).
2. **Orai × JKU** — 4 Orai1 wild-type conformations × 4 JKU ligands = **16 complexes**, ensemble cross-docking
   (no crystal reference → validity + interaction analysis only, no RMSD).

A third campaign, **Orai × Benchmark** (4 Orai frames × benchmark ligands, cross-docking), shares the same
code and is noted where it differs but is not the focus.

Each complex is docked by **three independent engines** — AutoDock Vina, DiffDock (DiffDock-L), and a custom
EquiBind pipeline — over a shared file-prep front-end and a shared PoseBusters/PandaMap QC back-end.

---

## 1. Software & versions (ground truth, verified on host)

Host: WSL2, Linux 5.15.167.4-microsoft-standard-WSL2. GPU: **NVIDIA RTX 4070 Laptop, 8 GB**, driver 581.80,
CUDA 12.4 (torch runtime).

| Tool | Version (verbatim) | Path / env |
|---|---|---|
| AutoDock Vina | **v1.2.7-20-g93cdc3d-mod** (self-built, *patched* — not stock 1.2.7) | `/home/manndo/AutoDock-Vina/build/linux/release/vina` |
| AutoGrid | AutoGrid 4.2.6 | `/usr/bin/autogrid4` (configs point at a **dead** `/usr/local/bin/autogrid4`; unused — vina scoring) |
| DiffDock | **DiffDock-L**, commit `85c49b6` (`git describe` = v1.1.3-2-g85c49b6), gcorso/DiffDock | `/home/manndo/docking_tools/DiffDock`, env `diffdock` |
| EquiBind | commit `4cb1b4c` (2025-02-19), HannesStark/EquiBind (no tag) | `/home/manndo/tools/EquiBind`, env `equibind` |
| gnina | **v1.3.2** (master:f23dd2b, built Jul 29 2025) — a **bash wrapper** injecting CUDA libs | `/home/manndo/docking_tools/gnina` |
| smina | Nov 9 2017 (based on AutoDock Vina 1.1.2) | `envs/diffdock/bin/smina` and `envs/equibind/bin/smina` |
| Meeko | 0.7.1 | env `vina` (`mk_prepare_receptor.py`, `mk_prepare_ligand.py`) |
| MGLTools | ADFRsuite-1.0 `prepare_receptor` (primary); MGLTools 1.5.7 (fallback) | `~/ADFRsuite-1.0/bin`; env `mgltools-py27` (Py 2.7.3) |
| Open Babel | Open Babel 3.1.1 (Mar 31 2024) | system PATH |
| reduce | reduce.4.16.250520 | `tools/reduce/reduce_src/reduce` (EquiBind protein prep only) |
| fpocket | fpocket 4.0 | `/usr/local/bin/fpocket` (EquiBind pockets) |
| P2Rank | 2.5 (release) / 2.5.2-dev.6 (source build) | `tools/p2rank_2.5/prank` / `docking_tools/p2rank/distro/prank` |
| PoseBusters | 0.6.3 | env `vina` |
| PandaMap | 4.1.0 | env `vina` |

Per-env interpreters: `diffdock` Py 3.9.23 / torch 2.4.0+cu124 / rdkit 2025.3.5; `equibind` Py 3.8.20 /
torch 2.4.1+cu124 / rdkit 2024.3.5; `vina` Py 3.12 / rdkit 2025.9.1 (CPU, no torch).

---

## 2. File-preparation tool chain (shared)

Orchestrated by `run_workflow()` in `Scripts/Utilities/prep_docking.py`, driven from the notebook.

### 2.1 Receptor preparation
1. **PDBFixer validation** (unless `skip_pdb_validation`, which is `false` everywhere): find/add missing
   residues + atoms, replace nonstandard residues, `removeHeterogens(keepWater=False)` (strips **all** HETATM +
   waters). Hydrogen addition at pH 7.4 is gated on `add_hydrogens`, which **no caller sets → no pH protonation
   at this stage**.
2. **Orai receptor cleaning** (`prepare_receptor_pdb.clean_pdb_like_pymol`, notebook): strip H; `HSD/HSE/HSP →
   HIS`; strip CHARMM caps (`CAY/CY/OY/HY1-3`, `NT/CAT/HNT/HT1-3`); **map CHARMM SEGID `MONA…MONF` → chains
   A–F**; force occupancy 1.00; add TER + `CRYST1`. Result: `Data/Receptors/*.pdb` (6 chains A–F, res 66–288 =
   1338 residues, 10 362 heavy atoms). Raw MD frames kept in `Data/Receptors/Original/`.
   Benchmark receptors (`*_protein.pdb`) are PoseBusters crystal structures used as-is.
3. **PDBQT conversion** (`prep_tool`: `mgltools` | `meeko` | `both`):
   - **MGLTools**: `prepare_receptor -r in.pdb -o out.pdbqt -A hydrogens -U nphs_lps_waters_nonstdres`
     (no pH flag). On H-addition crash → OpenMM/PDBFixer minimise (amber14 + gbn2, ≤500 iters) → retry,
     producing `*_minimised.*`.
   - **Meeko**: `mk_prepare_receptor.py -i in.pdb -o base -p -v -g --box_enveloping in.pdb --padding 2.0`,
     retry cascade adding `--allow_bad_res` then dropping bad residues (`*_cleaned.pdb`, `*_retry1`).
     Atom-type sanitize for stock Vina: `CG0/CG1/CG2/G0/G1/G2 → C`, `Si → S`, **`B → C`**.
4. **Cofactor/metal re-injection** (`inject_hetatms.py`, gated by `retain_hetatm_residues: true`) on the
   `run_autodock.py` CLI path — re-adds curated cofactors (HEM/NAD/FAD/ATP/…) + metals (MG/ZN/CA/FE/…),
   partial charge 0.0. Orai is apo → no-op. The benchmark **self-dock notebook cell does not wire this in**.

### 2.2 Search box
- **Orai (JKU + Orai×Benchmark)** — **whole-protein "blind" box**: bounding box of all atoms + **2.0 Å**
  padding, written to `<rec>.box.txt` (`center_*`, `size_*`) and passed to Vina via `--config`.
- **Benchmark self-dock** — **crystal-ligand-centred box**: `box_from_sdf(<id>_ligand.sdf, padding=10.0)`
  overwrites every `*.box.txt`. The box is centred on the crystal pose, but the molecule docked is the
  *generated start conformer*, not the crystal pose. No `box.pdb` generator; AutoGrid never runs (vina scoring).

### 2.3 Ligand preparation
- **JKU**: QM-optimised `*-OPT.xyz` → `obabel … -h` (adds H, perceives bonds) → `.pdb/.sdf/.mol2`
  → Meeko `mk_prepare_ligand.py -i lig.sdf -o lig.pdbqt` → atom-type sanitize (`B → C`). Boron handling:
  `boron-silicon-atom_par.dat` supplies AD4 params for B/Si (AutoGrid path); B→C for the vina path actually
  docked. (`xyz2mol` is installed but **not** used — Open Babel does the conversion.)
- **Benchmark**: PoseBusters `*_ligand_start_conf.sdf` (RDKit 3D start conformer) → Meeko → pdbqt.

---

## 3. Docking engines & parameters

### 3.1 AutoDock Vina
```
vina --receptor rec.pdbqt --config rec.box.txt --ligand lig.pdbqt --out out.pdbqt \
     --exhaustiveness 32 --num_modes 10 --energy_range 3 --cpu 32 --seed 42
```
| Param | Benchmark self-dock | Orai × JKU | Orai × Benchmark |
|---|---|---|---|
| scoring | vina | vina | vina |
| exhaustiveness | 32 | 32 | 32 |
| num_modes | 10 | 10 | 10 |
| energy_range | 3 kcal/mol | 3 | 3 |
| seed | 42 | 42 | 42 |
| --cpu | 32 | 32 | 32 |
| batch | per-ligand (1 lig/complex) | batch, 10 lig/call | batch, 50 lig/call |
| prep_tool | both | both | mgltools |
| timeout/complex | 600 s | 600 s | 600 s |

Output: multi-model PDBQT, one `MODEL` per pose with `REMARK VINA RESULT: <affinity> <rmsd_lb> <rmsd_ub>`,
≤10 poses/complex. No post-docking optimisation (Vina's own local optimiser is intrinsic).

### 3.2 DiffDock (DiffDock-L)
```
envs/diffdock/bin/python -m inference --config custom_inference_args.yaml \
    --protein_ligand_csv batch.csv --out_dir <out> \
    --samples_per_complex 10 --no_final_step_noise
```
- Blind (whole-protein) docking; **no user box**. Score model coarse-grained (receptor_radius 15 Å),
  ESM2 `esm2_t33_650M_UR50D` embeddings generated at inference from the PDB sequence.
- Input: protein PDB (rewritten `*_prepared.pdb`, HIS names normalised); ligand SDF/MOL2 (the `.pdbqt` in the
  JKU folder is ignored — SDF/MOL2 used).
- **samples_per_complex 10, inference_steps 20** (actual_steps 19), batch_size 15, group_by_receptor true,
  device cuda:0 (OOM → CPU fallback). **No sampling seed passed** (stochastic per run).
- Checkpoints: score `workdir/v1.1/score_model/best_ema_inference_epoch_model.pt`; confidence
  `.../confidence_model/best_model_epoch75.pt`. Poses re-ranked by the **confidence model** →
  `rank<k>_confidence<value>.sdf`, all 10 kept.
- *(A working-tree edit set `default_inference_args.yaml` samples 10→30, but CLI + config override force the
  effective value back to 10.)*

### 3.3 EquiBind (custom pipeline, `equibind_pipeline/`)
Three phases: CPU prep → GPU EquiBind inference (in-process, batch 8) → CPU post-processing. Runs in the
`equibind` env; model checkpoint `runs/flexible_self_docking/best_checkpoint.pt` (n_lays 8, 30 att heads,
noise disabled at inference).

**Pose generation per complex:**
- **Guided** poses from pockets: fpocket + p2rank pockets concatenated, top-`n_top_pockets` each,
  `poses_per_pocket` poses per pocket; ligand centroid translated to pocket centre, protein optionally
  cropped (radius 8 + 3 Å buffer = 11 Å, pocket at origin).
- **Unguided** poses: `n_unguided_poses`, ligand placed freely against the full protein.
- Ligand conformers: **RDKit ETKDGv3** (useRandomCoords, maxIters 500) + MMFF optimise; one conformer per
  seed; seed list of 30 (JKU) or 10 (Orai×Benchmark).
- Post-EquiBind geometry fix: torsion/von-Mises dihedral fit + rigid **Kabsch** alignment.

| Param | Orai × JKU | Orai × Benchmark |
|---|---|---|
| n_top_pockets | 3 (fpocket+p2rank) | 0 (guided OFF) |
| poses_per_pocket | 3 | — |
| n_unguided_poses | 30 | 10 |
| rdkit_seeds | 30 | 10 |
| poses/complex | (3·3·2)+30 = **48** | **10** |
| clamp_mode | off | off |
| **uff_minimize** | **true** | **false** (UFF mangles bond angles; gnina relieves clash → UFF off ~2× PB-valid) |
| **refine_mode** | **both** (`__refRAW` + refined) | on (single refined) |
| **refine_tool** | **gnina** | `smina` *(literal value — all its own comments say gnina; committed value is smina)* |

EquiBind emits only coordinates (no score) → poses are **ranked downstream by the smina/gnina refinement
energy** attached in refinement (§4).

---

## 4. Post-docking optimisation / refinement chain

Both ML methods lack a physics/excluded-volume term, so poses are locally relaxed against the receptor
(box around the pose, so it stays in the same pocket). Vina needs none.

### 4.1 DiffDock — post-hoc gnina minimisation (non-destructive; writes `optimized_gnina/`)
```
gnina --receptor <prot>_prepared.pdb --ligand rankK_confidence-X.sdf \
      --autobox_ligand rankK_confidence-X.sdf --autobox_add 4.0 \
      --out optimized_gnina/rankK_..._gnina.sdf --cpu 1 --seed 0 --num_modes 1 --minimize
```
`optimization: gnina`, `optimize_search: minimize` → `--minimize`; `gnina_use_gpu: true` → `--no_gpu`
**omitted** (CNN on GPU). Each of the 10 poses minimised individually; scores logged
(`minimizedAffinity`, `CNNscore`, `CNNaffinity`) to `optimization_log.csv`. Original DiffDock SDFs untouched.

### 4.2 EquiBind — inline refinement (Phase 3)
```
<smina|gnina> --receptor prot.pdb --ligand pose.sdf --autobox_ligand pose.sdf \
      --autobox_add 4.0 --out out.sdf --cpu 1 --seed 0 --num_modes 1 --local_only [--no_gpu]
```
`smina_search: local_only`, autobox_add 4.0, seed 0, cpu 1, timeout 300 s. For **Orai × JKU**: `refine_mode:
both` + `refine_tool: gnina` → each of 48 poses emits a `__refRAW` (unrefined) and a `__refGNINA` variant;
`gnina_use_gpu: false` → `--no_gpu` (CNN on CPU, GPU busy with EquiBind). A UFF minimisation (protein held
fixed, `vdwThresh 10`, 200 iters) precedes refinement when `uff_minimize` is on (JKU: on). On any refine
failure the pre-refine pose is copied through — **a pose is never lost**.

---

## 5. Dataset details

### 5.1 PoseBusters Benchmark set
- **428 complexes** (`Data/Ligands/PoseBuster_Benchmark_Set/*_ligand_start_conf.sdf` = 428;
  `Data/Receptors/PoseBuster_Benchmark_Set/*_protein.pdb` = 428; crystal refs in
  `Data/PoseBuster Benchmark Set/<PDBID>_<LIG>/` = 428). Naming `<PDBID>_<LIGCODE>` (e.g. `5S8I_2LY`).
- **Self-docking**: each ligand → its own crystal protein. Docking input is the **RDKit start conformer**
  (topology/stereo from the deposit, 3D pose decoyed), so the crystal answer cannot be read off.
- Native crystal pose exists → **RMSD-to-native** and near-native (≤2 Å) metrics available.
- **Coverage actually docked**: Vina **428**, DiffDock **428**, EquiBind **309** (the 308 unique
  `posebusters_pdb_ccd_ids.txt` subset). *(There is no dedicated committed DiffDock config for the benchmark;
  the same parameter-driven engine was pointed at benchmark inputs.)*
- Background: standard PoseBusters Benchmark set (428 structures, PDB depositions 2021+), Buttenschoen,
  Morris & Deane, *Chem. Sci.* 2024. *(Citation is background, not repo-verified.)*

### 5.2 Orai × JKU
- **Receptors (4)**: `Orai1WT-START-Fr0` (starting frame) + `MDSnap-Fr300/Fr400/Fr499` (**MD snapshots**) —
  ensemble docking across conformations of the same protein. Orai1 = pore-forming subunit of the Ca²⁺
  release-activated Ca²⁺ (CRAC) channel; **hexamer** (6 subunits, chains A–F, res 66–288). Origin CHARMM MD
  (inferred from `HSD/HSE/HSP` naming).
- **Ligands (4)**: `2abp-nh2-OPT` (boron-containing 2-APB / 2-aminoethoxydiphenylborate analogue, exactly one
  B), `Synta-66-OPT-Singlet`, `gsk7975a-prot-OPT` (**protonated**), `gsk7975a-deprot-OPT` (**deprotonated**) —
  the GSK-7975A compound treated as two protonation states = two ligands. `-OPT` = QM geometry-optimised.
  (`gsk7975a, Synta-66].pdb` is a mis-saved Open Babel duplicate, not a docking ligand.) All are known
  Orai/CRAC modulators *(background)*.
- **Design**: 4 receptors × 4 ligands = **16 complexes** per engine (confirmed: 16 DiffDock folders).
  **No crystal reference** → **no RMSD**; analysis is PoseBusters validity + PandaMap interactions + pose
  clustering + transmembrane-pore exclusion only.

---

## 6. Post-docking QC & interaction analysis (back-end)

- **PoseBusters 0.6.3** (`config_mode: dock`, receptor-aware): ~20 canonical checks — loading/chemistry,
  bond lengths/angles, internal clash/energy (50-conformer UFF ensemble), ring flatness, and (dock mode)
  min-distance + volume-overlap vs protein/cofactors/waters. AutoDock PDBQT→SDF via Meeko `mk_export` (SMILES
  header) or obabel + crystal-ligand bond-order template. Receptors HETATM-cleaned to match what the docker
  saw. **PoseBusters itself computes no RMSD** (`mol_true=None`); the ≤2 Å near-native metric (symmetry-
  corrected heavy-atom RMSD, no superposition) is computed downstream in `posebusters_pose_comparison.py`
  — **Benchmark only**.
- **PandaMap 4.1.0** (must run in `vina` env): residue-level protein–ligand interaction fingerprints, 16
  interaction types, top-N poses/complex. `benchmark_dir` set for Benchmark (native crystal fingerprint),
  `null` for Orai. EquiBind re-ranked by gnina affinity before selection.

---

## 7. Known config inconsistencies (report as-is)
- `autogrid_bin: /usr/local/bin/autogrid4` is a **dead path** (real binary `/usr/bin/autogrid4`); harmless
  because vina scoring never calls AutoGrid.
- EquiBind Orai×Benchmark config: `refine_tool: smina` while every comment says gnina; `skip_existing: true`
  while its comment says "force full re-dock". Literal values win (smina / skip). *(Orai × JKU is unaffected —
  it uses gnina, refine_mode both.)*
