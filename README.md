# master_dev

Code and data for the master's thesis on blind protein-ligand docking: a
calibration benchmark of AutoDock Vina, DiffDock and EquiBind over 308
PoseBusters complexes, and an application of the same three pipelines to the
Orai1 calcium channel.

## Start here

```bash
conda activate vina
jupyter lab Thesis_Reproduction.ipynb
```

`Thesis_Reproduction.ipynb` is the pipeline. It declares 60 stages, walks them in
dependency order, and finishes by recomputing 257 numbers the thesis prints and
comparing them to the document. The default run mode executes nothing expensive:
it checks that every stage's outputs are present and asserts the numbers, which
takes a few minutes. Re-running the docking itself is days of GPU and CPU time
and is opt-in per stage.

## Layout

| Path | What it holds |
| --- | --- |
| `Thesis_Reproduction.ipynb` | the pipeline |
| `XRay_PoseBusters_Control.ipynb` | standalone control: the PoseBusters battery run on the deposited crystal ligands |
| `Scripts/Analysis/` | analysis, statistics, the stage registry and the assertion harness |
| `Scripts/Docking/` | the docking drivers |
| `Scripts/Docking/configs_thesis/` | one documented config per reported arm |
| `Data/` | benchmark set, receptors, ligands |
| `Dockings/` | docked poses |
| `posebusters_results/`, `pandamap_results/`, `pocket_results/` | screening and analysis outputs |
| `results_current/` | symlink index onto the canonical outputs |
| `thesis_latex/` | the thesis sources and the built PDF |
| `obsolete/` | superseded material, gitignored, nothing deleted |

## Documentation

| File | Covers |
| --- | --- |
| `THESIS_REPRODUCTION.md` | how the pipeline was built and why each decision was made |
| `Scripts/Analysis/REGENERATE.md` | generated from the registry: every float, its source, the command that rebuilds it, the environment it needs and what cannot be reconstructed |
| `Scripts/Analysis/FINDINGS_2026-09-02.md` | defects found and their evidence, plus an md5-verified index of all 39 figures |
| `Scripts/Analysis/IMPLEMENTED_STATS_TESTS.md` | the statistical design |
| `CLEANUP_PLAN.md` | the repository cleanup and its outcome |

## Two things to know before changing anything

**The canonical trees are not the obviously-named ones.** Generation order,
oldest first, is `*_PRE_FR0` and `*_PRE_EXH128` and `*_UFFON_backup_*`, then the
plain name, then `_matched`. The current results live in
`benchmark_matched_equibind`, `_orai_matched_root` and the `_matched` PandaMap
trees. `repro_harness.py` holds the authoritative table and every stage reads its
paths from it.

**Some dependencies live in data, not in code.** The result tables store absolute
paths to the pose and receptor files behind each row, so a directory can be
load-bearing without any script naming it. Before moving or deleting a results
tree, check it with `scratchpad/data_refs.py`. The pose loader skips a missing
file silently, so the failure mode is a table that is quietly short of rows
rather than an error.

## Verifying

```bash
python Scripts/Analysis/thesis_assertions.py       # 257 checks against the thesis
python Scripts/Analysis/validate_regeneration.py   # re-runs 9 figure generators, ~10 min
python Scripts/Analysis/build_results_current.py --check
```

Each exits non-zero on failure.

## Environment

`Conda_Env/vina.yml` builds the screening and analysis environment: Python
3.12.12, RDKit 2025.09.1, PoseBusters 0.6.3. The RDKit pin matters, because most
PoseBusters verdicts depend on that exact release. Two patches live outside this
repository and are load-bearing, one to PoseBusters and one to the local DiffDock
checkout; the notebook asserts both before running anything.
