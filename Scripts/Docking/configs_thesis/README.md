# Thesis configuration set

One documented config per arm reported in the thesis. Read by
`Thesis_Reproduction.ipynb`, which never writes them back.

Each file is its current counterpart in `Scripts/Docking/` copied verbatim with a
provenance header prepended. The generator asserts that every output parses to the
same dictionary as its source, so behaviour cannot drift. The originals stay where
they are, so nothing that currently runs breaks.

The header of each file names the thesis section it implements, the floats it
supports, the tree it writes, the analysis tree that reads it, any deviation from
the source, and its determinism class.

## Index

| File | Arm | Determinism |
| --- | --- | --- |
| `01_benchmark_autodock_exh128_raw.yaml` | AutoDock Vina search, exhaustiveness 128 | bit-exact |
| `02_benchmark_autodock_exh128_gnina.yaml` | **The dominant arm.** gnina rescoring of those poses | bit-exact |
| `03_benchmark_diffdock.yaml` | DiffDock-L, 30 samples | verify only |
| `04_benchmark_equibind.yaml` | EquiBind, unguided and pocket-guided | bit-exact |
| `05_benchmark_posebusters.yaml` | PoseBusters, all benchmark arms | bit-exact |
| `06_benchmark_posebusters_matched_equibind.yaml` | **The canonical benchmark screen** | bit-exact |
| `07_benchmark_pandamap.yaml` | Interaction fingerprints | bit-exact |
| `10`–`15_ladder_*.yaml` | Exhaustiveness ladder at 18, 32, 64, 92, plus two rescored rungs | bit-exact |
| `20_orai_jku_autodock_exh128_gnina.yaml` | Orai experimental AutoDock, 30 modes | bit-exact |
| `21_orai_jku_diffdock.yaml` | Orai experimental DiffDock, **the one seeded DiffDock arm** | bit-exact |
| `22_orai_jku_equibind.yaml` | Orai experimental EquiBind | bit-exact |
| `23_orai_jku_posebusters.yaml` | Orai experimental validity screen | bit-exact |
| `24_orai_jku_posebusters_fr0corrected.yaml` | Fr0 correction, **AutoDock rows only** | bit-exact |
| `25_orai_jku_pandamap.yaml` | Orai experimental fingerprints | bit-exact |
| `30_orai_benchmark_autodock_exh128_gnina.yaml` | Orai control AutoDock, 10 modes | bit-exact |
| `31_orai_benchmark_diffdock.yaml` | Orai control DiffDock, 10 samples | verify only |
| `32_orai_benchmark_equibind.yaml` | Orai control EquiBind | bit-exact |
| `33_orai_benchmark_posebusters.yaml` | Orai control validity screen | bit-exact |
| `34`, `35_orai_benchmark_posebusters_*fr0corrected.yaml` | Fr0 correction, **AutoDock rows only** | bit-exact |
| `36_orai_benchmark_pandamap.yaml` | Orai control fingerprints | bit-exact |

## Three things worth knowing before editing any of these

**The EquiBind configs replace a deprecated file.** Until now the benchmark and
Orai experimental EquiBind arms both read
`Scripts/Docking/equibind_docking_config.yaml`, which is marked deprecated and
which a notebook cell rewrote in place, setting `uff_minimize` to true. Appendix B
records that the in-protein UFF relaxation was switched off for every reported
run. Files 04, 22 and 32 carry `uff_minimize: false` and are never written to.

**The Fr0 correction is AutoDock-only.** Files 24, 34 and 35 exist because
AutoDock docked the Fr0 frame into the OpenMM-minimised receptor while the main
screen validates against the raw PDB. DiffDock and EquiBind docked the raw
receptor, so extending the correction to them creates the mismatch rather than
repairing it. The `.WRONG-RECEPTOR-quarantined-1854` directories under
`posebusters_results/` are that mistake, made once.

**Two DiffDock arms record `inference_seed: null` on purpose.** Files 03 and 31
describe runs that were unseeded single draws, established by counting the seed
patch's marker in the stored run logs: zero for the benchmark, zero for the Orai
control, twenty-four for the Orai experimental. Writing `42` there would claim a
reproducibility the stored trees do not have.

## Regenerating this directory

```bash
/home/manndo/anaconda3/envs/vina/bin/python Scripts/Docking/configs_thesis/_generate.py
```

It copies each source verbatim, prepends the header, applies any declared bake,
then asserts `yaml.safe_load(output) == yaml.safe_load(source)` with the bake
applied. If that assertion fails the copy is wrong and must not be used, which is
why it runs on every generation rather than once at authoring time.

The notebook itself is built the same way, by
`Scripts/Analysis/_build_reproduction_notebook.py`. Edit the builder rather than
the `.ipynb`, so the cell text stays reviewable in a diff.

Rationale for the whole set is in `THESIS_REPRODUCTION.md`, section 3.4.
