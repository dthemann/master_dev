# Future Research

Proposed next steps arising from the Orai1 chapter, written 2026-09-02.

Everything here descends from one open question. The transmembrane exclusion slab is an **operational
hypothesis**, not ground truth, and the thesis says so. The M4 correction then showed the rule interacts
with the ranking step: over 12,315 control poses the gnina order correlates with slab occupancy at
Spearman −0.211 while the Vina order sits at −0.004. So poses land in the slab partly because a learned
affinity ranker puts them there. What the chapter cannot currently say is whether they would still land
there if the receptor carried the physics that was stripped from it.

The rationale, the geometry and the risk analysis live in
[`thesis_latex/PLAN_orai_slab_physics_test.md`](../thesis_latex/PLAN_orai_slab_physics_test.md).
This folder holds the **executable instructions**, one file per package.

## Order of work

| # | package | cost | status | file |
|---|---|---|---|---|
| 1 | Write down the recovered MD provenance | one paragraph | **do first** | [E](E_recovered_provenance.md) |
| 2 | Report placement as a depth, not a flag | no docking | do next | [A](A_continuous_depth_reporting.md) |
| 3 | Is the slab preference a burial preference? | no docking | if time allows | [B](B_burial_mechanism.md) |
| 4 | Membrane-aware docking arm | moderate | future work, feasible | [D](D_membrane_aware_docking.md) |
| 5 | Calcium what-if | 12 dockings | optional, downgraded | [C](C_calcium_whatif.md) |

E and A are thesis work and cost almost nothing. B is thesis work if time allows. D is the genuine
research project. C is kept only because it is cheap and someone will ask about it.

## The 2026-09-02 delivery, which changed all of this

The MD system arrived from the collaborating group and is at

```
/mnt/c/Users/domin/OneDrive/00_Endless_Learning/FH Technikum/AI Engineering/Master Project/Orai1-Docking-MScProjekt-DMann/Orai/MDFiles/
```

| file | what it is |
|---|---|
| `hORAI1_WT_ARN.psf` | CHARMM topology, 182,692 atoms, CHARMM-GUI v3.7, 6 April 2022 |
| `Orai1_WT_bar3.dcd` | production trajectory, 500 frames |
| `Orai1_WT_fix1.dcd` | protein-restrained equilibration, 501 frames |
| `FH Wien.7z` | 1.77 GB archive, **not yet opened** |

**Frame identity, established by CA-RMSD after Kabsch superposition.** This is exact, not approximate:

| supplied snapshot | source | CA RMSD |
|---|---|---|
| Fr300 | `bar3` frame 300 | 0.0005 Å |
| Fr400 | `bar3` frame 400 | 0.0005 Å |
| Fr499 | `bar3` frame 499 | 0.0005 Å |
| Fr0 | **`fix1` frame 0** | 0.0006 Å |

`fix1` moves 0.029 Å across its whole length, so it is a restrained-protein equilibration in which
lipid and water relax around a held protein. **That explains the Fr0 anomaly the thesis already
measured**: Fr0 sits 5.9 Å outside the ensemble the other three sample because it is not a production
frame at all.

**System composition.** 426 lipids (213 cholesterol, 71 sphingomyelin, 71 DLPC, 71 DLPE), 40,123 TIP3
waters, 109 K+, 150 Cl-, and **zero calcium**. Segments `IONS MEMB MONA MONB MONC MOND MONE MONF TIP3`.

## Two things that must NOT be claimed

1. **The per-frame simulation time.** The DCD header reports `dt = 0.049 ps`, which is almost certainly
   unscaled. Do not quote a trajectory length or a frame spacing until it is confirmed against the run
   input files, which may be inside `FH Wien.7z`.
2. **The homology-model template PDB IDs.** Still not in evidence anywhere.

## The standing caveat

None of these packages can validate the Orai1 poses. They test whether the exclusion rule is defensible,
not whether the predicted binding modes are right. That still needs site-directed mutagenesis, a
competition or displacement assay, or a ligand-bound structure.

## Conventions

- Python: `/home/manndo/anaconda3/envs/vina/bin/python`, conda env `vina`, run from `/home/manndo/master_dev`.
- Docking outputs go under `posebusters_results/<tree>/`; nothing there is tracked by git.
- Any new analysis must get a row in `Scripts/Analysis/REGENERATE.md` giving the exact command, or it
  will be unreproducible within a week. This has already bitten the project once.
- Any new figure is hand-copied into `thesis_latex/media/media/imageNN.png` and verified by `md5sum`.
  `media/` is **shared with the full build**, so replacing an image changes both documents.
