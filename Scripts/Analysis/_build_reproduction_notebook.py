#!/usr/bin/env python3
"""Build Thesis_Reproduction.ipynb from a declarative cell list."""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path("/home/manndo/master_dev")
OUT = ROOT / "Thesis_Reproduction.ipynb"

cells: list[dict] = []


def _cid() -> str:
    return f"cell-{len(cells):03d}"


def md(text: str) -> None:
    cells.append({"cell_type": "markdown", "id": _cid(), "metadata": {},
                  "source": text.strip("\n").splitlines(keepends=True)})


def code(text: str) -> None:
    cells.append({"cell_type": "code", "id": _cid(), "execution_count": None,
                  "metadata": {}, "outputs": [],
                  "source": text.strip("\n").splitlines(keepends=True)})


# =============================================================================
md(r"""
# Thesis reproduction pipeline

Regenerates the results reported in the master thesis, from the docked poses on
disk through to the figures and tables the document prints.

**Run it top to bottom.** The default mode runs nothing expensive: it walks the
whole pipeline, checks that every stage's outputs are present, and finishes by
recomputing the headline numbers and comparing them to the values in the thesis.
That takes a few minutes. Executing the pipeline for real is days of GPU and CPU
time and is opt-in per stage.

### What this notebook replaces

`Master_Docking_AD_Full_Protein.ipynb` produced these results and is kept as the
historical record, but it is not runnable as a reproduction. It carries 125 cells
of which about a third are abandoned arms, one-off repairs or dead code; its page
order is not its execution order, so three Orai analysis cells read a file that a
cell forty positions later writes; and three of its cells rewrite YAML configs on
disk. One of those flips `uff_minimize` back to `true`, which is the setting the
reported EquiBind arm depends on being `false`.

Here, configs are read-only, every stage declares what it needs, and the driver
refuses to run a stage before its inputs exist.

### Three determinism classes

Not everything in this pipeline can be reproduced, and the difference is a fact
about the recorded runs rather than a choice made here.

| Class | Meaning |
| --- | --- |
| `bitexact` | Seeded end to end. A re-run must reproduce the stored bytes. |
| `verify` | The reported run was an unseeded single draw. It can be checked, not redrawn. |
| `never` | Not reproducible even in principle, and re-running destroys the provenance of a reported number. |

The `verify` class is exactly two stages, both DiffDock. The evidence is printed
below rather than asserted: the local seed patch prints a marker when it fires,
and that marker appears in none of the benchmark run logs, none of the Orai
control logs, and all of the Orai experimental ones. Appendix B of the thesis
records the same asymmetry.
""")

code(r'''
import os, sys, subprocess, json, hashlib
from pathlib import Path

ROOT = Path("/home/manndo/master_dev")
os.chdir(ROOT)
sys.path.insert(0, str(ROOT / "Scripts/Analysis"))

import repro_harness as H
from repro_harness import (Stage, Registry, Driver, CANONICAL, SUPERSEDED,
                           BITEXACT, VERIFY, NEVER, CHEAP, MODERATE, EXPENSIVE,
                           VINA_PY, ROOT as HROOT)

# ── run mode ────────────────────────────────────────────────────────────────
# "verify"      check every stage's outputs exist; execute nothing.  DEFAULT.
# "run-missing" additionally execute stages whose outputs are absent.
# "force"       additionally re-execute the stages named in FORCE_STAGES.
#
# Forcing everything is deliberately not possible. It would re-run days of
# docking and would overwrite the two unseeded DiffDock trees with a different
# sample, which is the one thing that cannot be undone.
RUN_MODE     = "verify"
FORCE_STAGES = []

CFG   = ROOT / "Scripts/Docking/configs_thesis"
ANA   = ROOT / "Scripts/Analysis"
DOCK  = ROOT / "Scripts/Docking"
PB    = DOCK / "Posebusters"
IDS   = ROOT / "Data/PoseBuster Benchmark Set/posebusters_pdb_ccd_ids.txt"

print(f"root      : {ROOT}")
print(f"run mode  : {RUN_MODE}")
print(f"configs   : {CFG.relative_to(ROOT)}  ({len(list(CFG.glob('*.yaml')))} files)")
''')

# =============================================================================
md(r"""
## 0. Environment and provenance guard

Appendix B names every version this pipeline depends on. Two of them are patches
that live outside the repository and would otherwise be invisible.

The **PoseBusters patch** fixes the handler guarding the UFF parameter check. As
shipped, it indexes the second argument of an assertion that carries only one,
raising `IndexError` and aborting the entire `bust()` call. Every pose of any
ligand UFF cannot parameterise, such as one carrying hypervalent sulfur, a
phosphate or a coordinated metal, would be discarded with no error reported.

The **DiffDock patch** adds a seed hook. Upstream DiffDock exposes no seed option
at all, so its sampling differs on every run.

A failure here means the numbers below will not match, so it is worth reading the
table rather than skipping it.
""")

code(r'''
checks = H.check_environment()
width = max(len(c.name) for c in checks)
for c in checks:
    print(f"  {'ok ' if c.ok else 'DRIFT':5s} {c.name:{width}s}  {c.found[:46]:48s} want: {c.expected}")

failed = [c for c in checks if not c.ok]
print()
print(f"{len(checks) - len(failed)}/{len(checks)} checks pass")
if failed:
    print("\nDRIFT DETECTED - the reported numbers are tied to these versions:")
    for c in failed:
        print(f"  {c.name}: found {c.found!r}, expected {c.expected!r} {c.note}")
''')

code(r'''
# The evidence behind the `verify` determinism class, recomputed rather than
# asserted. The DiffDock seed patch prints "[repro] seeded" when it fires, so
# counting that marker in the stored run logs says which arms were seeded.
seeded = H.diffdock_seed_evidence()
print("DiffDock runs carrying the seed marker:")
for panel, n in seeded.items():
    verdict = "seeded, reproducible" if n > 0 else "UNSEEDED single draw - verify only"
    print(f"  {panel:20s} {n:4d} logs   {verdict}")
print()
print("Appendix B states the same: the reproducibility hook appears in the")
print("experimental run logs and is absent from the control ones.")
''')

# =============================================================================
md(r"""
## 1. Canonical trees

There are 57 directories under `posebusters_results/` and 22 under
`pandamap_results/`. Most are superseded. The trap is that the ones whose names
read like scratch are the current ones: generation order is
`*_PRE_FR0` and `*_PRE_EXH128` and `*_UFFON_backup_*`, then the plain name, then
`_matched`.

This was settled by checksum rather than by naming convention. The shipped
Figure 9 asset, `thesis_latex/media/media/image43.png`, matches the copy under
`_orai_matched_root` and not the plain one.

Every stage below reads its input path from this one table, so an analysis
cannot silently pick up the previous generation.
""")

code(r'''
print("CANONICAL")
for k, v in CANONICAL.items():
    mark = "ok " if v.exists() else "MISS"
    print(f"  {mark} {k:28s} {v.relative_to(ROOT)}")

print("\nSUPERSEDED - present on disk, deliberately NOT read")
for old, replacement in SUPERSEDED.items():
    if (ROOT / old).exists():
        print(f"      {old:58s} -> {replacement}")
''')

code(r'''
# The cohort. 308 ligands were prepared and docked; 303 carry a complete record
# on all three pipelines and are the analysed set. The distinction matters
# because every paired test uses the 303.
ids = [ln.strip() for ln in IDS.read_text().splitlines() if ln.strip()]
assert len(ids) == 308, f"expected 308 benchmark ids, found {len(ids)}"
print(f"benchmark ids           : {len(ids)}  ({ids[0]} ... {ids[-1]})")
print(f"analysed cohort         : 303   (the five without a complete three-tool record are dropped)")
print(f"Orai receptor frames    : 4     (START-Fr0, MDSnap-Fr300, Fr400, Fr499)")
print(f"Orai experimental panel : 3 modulators x 4 frames = 12 units")
print(f"Orai control panel      : 308 ligands x 4 frames = 1,232 units")
''')

# =============================================================================
md(r"""
## 2. Benchmark docking

The calibration benchmark: 308 PoseBusters ligands docked blind against the whole
receptor by three pipelines.

The AutoDock arm is an exhaustiveness ladder. Five independent searches at 18, 32,
64, 92 and 128 share prepared ligands, receptors, box and seed, so exhaustiveness
is the only variable and the rungs are directly comparable. Three of them were
additionally rescored with gnina, which reorders a rung's existing poses rather
than searching again. The reported arm is exhaustiveness 128 with gnina
rescoring; the rest support the appendix ladder.

Two stages here refuse to run under any mode. The prepared inputs are the shared
root of every rung, and re-deriving them with a newer Meeko changes atom order and
`TORSDOF`, which is why the driver script takes `--prepared-inputs-from` instead.
The arm status file is the only surviving record of the gnina wall clock the
thesis cites, and a cached re-run overwrites it with a near-zero value.
""")

code(r'''
reg = Registry()

# ---- inputs that must never be regenerated ---------------------------------
reg.add(Stage(
    name="bench_prepared_inputs", section="2. Benchmark docking",
    title="Prepared ligand and receptor PDBQT (ADFRsuite)",
    determinism=NEVER, cost=EXPENSIVE,
    outputs=["Dockings/vina_results_full_protein_vina_scoring_mgltools/_staging/ligands/pdbqt/*.pdbqt"],
    min_matches=300,
    thesis="App. B: prepare_ligand -A hydrogens, prepare_receptor -A hydrogens "
           "-U nphs_lps_waters_nonstdres",
    notes="Shared by every ladder rung. A newer Meeko changes atom ordering, the "
          "torsion-tree root and TORSDOF, so these are passed forward with "
          "--prepared-inputs-from rather than re-derived."))

reg.add(Stage(
    name="bench_arm_status", section="2. Benchmark docking",
    title="Exhaustiveness arm status (timing provenance)",
    determinism=NEVER, cost=CHEAP,
    outputs=["Dockings/vina_results_full_protein_vina_scoring_mgltools_exh128/exhaustiveness_arm_status.json"],
    thesis="App. B hardware and timing basis: the 2,677.7 s gnina wall clock",
    notes="A cached re-run overwrites the measured elapsed time with the cost of "
          "the cache lookup, which is near zero. No backup exists."))

# ---- the exhaustiveness ladder ---------------------------------------------
LADDER = [
    ("exh18",  "10_ladder_exh18.yaml",       "vina_results_full_protein_vina_scoring_mgltools_exh18",  18),
    ("exh32",  "11_ladder_exh32_raw.yaml",   "vina_results_full_protein_vina_scoring_mgltools",        32),
    ("exh64",  "13_ladder_exh64_raw.yaml",   "vina_results_full_protein_vina_scoring_mgltools_exh64",  64),
    ("exh92",  "15_ladder_exh92.yaml",       "vina_results_full_protein_vina_scoring_mgltools_exh92",  92),
    ("exh128", "01_benchmark_autodock_exh128_raw.yaml",
     "vina_results_full_protein_vina_scoring_mgltools_exh128", 128),
]
for tag, cfg, tree, exh in LADDER:
    reg.add(Stage(
        name=f"bench_autodock_{tag}_raw", section="2. Benchmark docking",
        title=f"AutoDock Vina search, exhaustiveness {exh}",
        determinism=BITEXACT, cost=EXPENSIVE,
        needs=["bench_prepared_inputs"],
        outputs=[f"Dockings/{tree}/*/mgl_tools/docking"], min_matches=300,
        thesis=("Table 1 and the Results chapter" if exh == 128
                else "Table 8 and Figures 18-21, the appendix ladder"),
        cmd=[VINA_PY, DOCK / "run_autodock_exhaustiveness_arm.py", "-c", CFG / cfg,
             "--ids-file", f"Dockings/{tree}/all308_ids.txt",
             "--expect-exhaustiveness", str(exh),
             "--prepared-inputs-from",
             "Dockings/vina_results_full_protein_vina_scoring_mgltools"],
        notes="Vina seed 42, energy range 6 kcal/mol, up to 30 modes, whole-receptor box."))

RESCORED = [("exh32", "12_ladder_exh32_gnina.yaml", "vina_results_full_protein_vina_scoring_mgltools"),
            ("exh64", "14_ladder_exh64_gnina.yaml", "vina_results_full_protein_vina_scoring_mgltools_exh64"),
            ("exh128", "02_benchmark_autodock_exh128_gnina.yaml",
             "vina_results_full_protein_vina_scoring_mgltools_exh128")]
for tag, cfg, tree in RESCORED:
    reg.add(Stage(
        name=f"bench_autodock_{tag}_gnina", section="2. Benchmark docking",
        title=f"gnina rescoring of the exhaustiveness {tag[3:]} poses",
        determinism=BITEXACT, cost=EXPENSIVE,
        needs=[f"bench_autodock_{tag}_raw"],
        outputs=[f"Dockings/{tree}/*/mgl_tools/docking/optimized_gnina"], min_matches=300,
        thesis=("the DOMINANT arm, Tables 1/2/6 and Figures 1-3, 5-7"
                if tag == "exh128" else "Table 8, the appendix ladder"),
        cmd=[VINA_PY, DOCK / "run_autodock_exhaustiveness_arm.py", "-c", CFG / cfg,
             "--ids-file", f"Dockings/{tree}/all308_ids.txt"],
        notes="Rescoring reorders a rung's own poses; it does not search again. Ranked "
              "by CNNaffinity on the crossdock_default2018 ensemble, seed 42."))
print(f"{len(reg)} stages registered")
''')

code(r'''
# ---- DiffDock and EquiBind --------------------------------------------------
reg.add(Stage(
    name="bench_diffdock", section="2. Benchmark docking",
    title="DiffDock-L, 30 samples, 20 denoising steps",
    determinism=VERIFY, cost=EXPENSIVE,
    outputs=["Dockings/Benchmark_DiffDock/*/docking_log.csv"], min_matches=300,
    thesis="Tables 1, 2, 6 and Figures 1-3, 5-7, 15, 16",
    cmd=[VINA_PY, DOCK / "run_diffdock.py", "-c", CFG / "03_benchmark_diffdock.yaml"],
    notes="UNSEEDED. None of the stored logs carries the seed marker, so re-running "
          "draws a new sample and cannot reproduce the reported poses. The driver "
          "refuses to run this stage; restore the tree from backup instead."))

reg.add(Stage(
    name="bench_pockets", section="2. Benchmark docking",
    title="fpocket and P2Rank site detection",
    determinism=BITEXACT, cost=MODERATE,
    outputs=["pocket_results/fpocket_results/*", "pocket_results/p2rank_results/*"],
    min_matches=100,
    thesis="App. B pocket-guided EquiBind variants; also the cluster-to-pocket analysis",
    notes="Prerequisite for the guided EquiBind arms and for the crystal-pocket "
          "cluster report. In the old notebook this sat under an EquiBind heading, "
          "which hid the dependency."))

reg.add(Stage(
    name="bench_equibind", section="2. Benchmark docking",
    title="EquiBind, 30 conformers, unguided and pocket-guided",
    determinism=BITEXACT, cost=EXPENSIVE,
    needs=["bench_pockets"],
    outputs=["Dockings/Benchmark_Equibind/*/pipeline_summary.json"], min_matches=300,
    thesis="Tables 1, 2, 6; the guided-arm paragraph of App. B",
    notes="Config 04 replaces the deprecated equibind_docking_config.yaml. "
          "uff_minimize is false, which App. B records as the setting for every "
          "reported run; the old notebook rewrote it to true."))

reg.add(Stage(
    name="bench_equibind_minimize", section="2. Benchmark docking",
    title="EquiBind re-refinement with --minimize (matched arm)",
    determinism=BITEXACT, cost=EXPENSIVE,
    needs=["bench_equibind"],
    outputs=["Dockings/Benchmark_Equibind_minimize/*"], min_matches=100,
    thesis="examiner item M5: every EquiBind refinement had run --local_only "
           "while AutoDock and DiffDock ran --minimize",
    cmd=[VINA_PY, DOCK / "rerun_equibind_refine.py"],
    notes="Re-refines from the committed __refRAW.sdf, which is the exact input the "
          "original refinement received, so the flag is the only variable. This is "
          "what makes benchmark_matched_equibind the canonical tree."))
print(f"{len(reg)} stages registered")
''')

# =============================================================================
md(r"""
## 3. Benchmark validation and analysis

PoseBusters first, then the pose-comparison hub, then everything that reads its
per-pose table. The ordering here is load-bearing: the hub writes
`per_pose_metrics.csv`, and nine downstream stages join against it.

Two flags recur and both change results.

`--collapse-plots-only` is what makes the variant pins take effect. Without it,
`--collapse-autodock-variant` and `--collapse-diffdock-variant` are accepted and
silently ignored, because the selection step only runs under the collapse flag.

`--exclude-preset meeko` drops the retired Meeko-prepared arms. On its own it is
not enough for the validity figure, because the report folds every surviving
`autodock_mgltools*` key onto a bare `autodock`, so the five unreported ladder
rungs have to be excluded by name as well or one panel silently pools five arms.
""")

code(r'''
BENCH_REPORT = CANONICAL["benchmark_report"].relative_to(ROOT)
BENCH_PB     = CANONICAL["benchmark_pb"].relative_to(ROOT)
LADDER_EXCL  = ("autodock_mgltools,autodock_mgltools_gnina,autodock_mgltools_exh18,"
                "autodock_mgltools_exh64,autodock_mgltools_exh64_gnina,autodock_mgltools_exh92")

reg.add(Stage(
    name="bench_posebusters", section="3. Benchmark analysis",
    title="PoseBusters validity screen, all benchmark arms",
    determinism=BITEXACT, cost=EXPENSIVE,
    needs=["bench_autodock_exh128_gnina", "bench_diffdock", "bench_equibind"],
    outputs=["posebusters_results/benchmark_full_protein_vina_scoring/dock/posebusters_filtered_results.csv"],
    thesis="App. B validity screening; PoseBusters 0.6.3",
    cmd=[VINA_PY, PB / "run_posebusters.py", "-c", CFG / "05_benchmark_posebusters.yaml"]))

reg.add(Stage(
    name="bench_posebusters_matched", section="3. Benchmark analysis",
    title="PoseBusters, matched-EquiBind arm (CANONICAL)",
    determinism=BITEXACT, cost=EXPENSIVE,
    needs=["bench_equibind_minimize"],
    outputs=[f"{BENCH_PB}/posebusters_filtered_results.csv"],
    thesis="every benchmark number in the Results chapter",
    cmd=[VINA_PY, PB / "run_posebusters.py",
         "-c", CFG / "06_benchmark_posebusters_matched_equibind.yaml"]))

reg.add(Stage(
    name="bench_pose_comparison", section="3. Benchmark analysis",
    title="Pose-comparison hub: per-pose metrics and recovery oracle",
    determinism=BITEXACT, cost=MODERATE,
    needs=["bench_posebusters_matched"],
    outputs=[f"{BENCH_REPORT}/per_pose_metrics.csv", f"{BENCH_REPORT}/oracle_summary.csv"],
    thesis="Table 1; the source table for nine downstream stages",
    cmd=[VINA_PY, ANA / "posebusters_pose_comparison.py",
         "--pb-csv", f"{BENCH_PB}/posebusters_filtered_results.csv",
         "--out-dir", str(BENCH_REPORT),
         "--top-n", "15", "--diffdock-variant", "all",
         "--collapse-plots-only",
         "--collapse-diffdock-variant", "diffdock_smina",
         "--collapse-autodock-variant", "autodock_gnina",
         "--exclude-preset", "meeko", "--workers", "24"],
    notes="--top-n 15 and --diffdock-variant all shaped the cached table. A re-render "
          "with --reuse-cache omits them, so dropping the cache after changing the "
          "input silently rebuilds at --top-n 5 on the raw diffdock key."))
print(f"{len(reg)} stages registered")
''')

code(r'''
reg.add(Stage(
    name="bench_validity_report", section="3. Benchmark analysis",
    title="Validity report per tool and complex",
    determinism=BITEXACT, cost=CHEAP,
    needs=["bench_pose_comparison"],
    outputs=[f"{BENCH_PB}/validity_report_mgltools"],
    thesis="Table 1 validity block",
    cmd=[VINA_PY, ANA / "posebusters_validity_report.py",
         "--csv", f"{BENCH_PB}/posebusters_filtered_results.csv",
         "--out-dir", f"{BENCH_PB}/validity_report_mgltools",
         "--per-pose-metrics", f"{BENCH_REPORT}/per_pose_metrics.csv",
         "--collapse-plots-only",
         "--exclude-preset", "meeko", "--exclude-methods", LADDER_EXCL]))

reg.add(Stage(
    name="bench_figure1", section="3. Benchmark analysis",
    title="Figure 1: post-hoc optimisation and PoseBusters validity",
    determinism=BITEXACT, cost=CHEAP,
    needs=["bench_posebusters_matched"],
    outputs=[f"{BENCH_PB}/validity_report_mgltools/00_figure2_validity_yield.png"],
    thesis="Figure 1 (media/media/image2.png)",
    cmd=[VINA_PY, ANA / "figure2_validity_yield.py",
         "--csv", f"{BENCH_PB}/posebusters_filtered_results.csv",
         "--out", f"{BENCH_PB}/validity_report_mgltools/00_figure2_validity_yield.png",
         "--exclude-preset", "meeko", "--exclude-methods", LADDER_EXCL],
    notes="The double exclusion is load-bearing. --exclude-preset meeko alone leaves "
          "five ladder rungs to be folded onto a bare 'autodock' key."))

reg.add(Stage(
    name="bench_filmstrip", section="3. Benchmark analysis",
    title="Figure 4: form against placement across ranking depth",
    determinism=BITEXACT, cost=CHEAP,
    needs=["bench_pose_comparison"],
    outputs=[f"{BENCH_REPORT}/20d_form_vs_placement_by_family__depth_filmstrip_pbvalid__rank1_top5_top15__thesis.png"],
    thesis="Figure 4 (media/media/image5.png); Table 27",
    cmd=[VINA_PY, ANA / "filmstrip_rank1_top5_top15.py",
         "--report-dir", str(BENCH_REPORT), "--exclude-preset", "meeko"],
    notes="--report-dir was added 2026-09-02. The module constant still names the "
          "pre-matched tree, and the shipped figure comes from the matched one."))

reg.add(Stage(
    name="bench_topk", section="3. Benchmark analysis",
    title="Top-k recovery sidecar for the gnina arm",
    determinism=BITEXACT, cost=CHEAP,
    needs=["bench_pose_comparison"],
    outputs=[f"{BENCH_REPORT}/topk_recovery_validity_gnina_arm.csv"],
    thesis="Tables 2, 18, 25",
    cmd=[VINA_PY, ANA / "topk_recovery_gnina_arm.py",
         "--report-dir", str(BENCH_REPORT),
         "--autodock-variant", "autodock_mgltools_exh128_gnina",
         "--diffdock-variant", "diffdock_smina",
         "--equibind-variant", "equibind_unguided_gnina",
         "--exclude-preset", "meeko"]))

reg.add(Stage(
    name="bench_optimization_benefit", section="3. Benchmark analysis",
    title="Optimisation benefit, paired categorical inference",
    determinism=BITEXACT, cost=CHEAP,
    needs=["bench_pose_comparison"],
    outputs=[f"{BENCH_REPORT}/optimization_benefit_by_rank.csv"],
    thesis="Tables 19, 20",
    cmd=[VINA_PY, ANA / "optimization_benefit_stats.py",
         "--report-dir", str(BENCH_REPORT), "--exclude-preset", "meeko"]))

reg.add(Stage(
    name="bench_endpoint_diagnostics", section="3. Benchmark analysis",
    title="Endpoint diagnostics: validity gate cost and bounded claim",
    determinism=BITEXACT, cost=CHEAP,
    needs=["bench_pose_comparison"],
    outputs=[f"{BENCH_REPORT}/validity_gate_cost.csv"],
    thesis="Table 22 and the Results prose numbers",
    cmd=[VINA_PY, ANA / "thesis_endpoint_diagnostics.py",
         "--metrics", f"{BENCH_REPORT}/per_pose_metrics.csv",
         "--csv", str(BENCH_REPORT), "--exclude-preset", "meeko"]))
print(f"{len(reg)} stages registered")
''')

code(r'''
CLUST = CANONICAL["benchmark_clusters"].relative_to(ROOT)
CLUST_ARM = f"{CLUST}/autodock_mgltools_exh128_gnina__diffdock_smina_allposes"
PMAP = CANONICAL["benchmark_pandamap"].relative_to(ROOT)

reg.add(Stage(
    name="bench_clusters", section="3. Benchmark analysis",
    title="Pose clusters against the crystal pocket",
    determinism=BITEXACT, cost=MODERATE,
    needs=["bench_pose_comparison", "bench_pockets"],
    outputs=[f"{CLUST_ARM}/topN_crystal_cluster_matrix.png",
             f"{CLUST_ARM}/cluster_quality_metrics.png"],
    thesis="Table 3; Figures 5, 6, 38",
    cmd=[VINA_PY, ANA / "pose_cluster_crystal_pocket_report.py",
         "--ids-file", str(IDS),
         "--per-pose-csv", f"{BENCH_REPORT}/per_pose_metrics.csv",
         "--out-dir", CLUST_ARM,
         "--autodock-variant", "autodock_mgltools_exh128_gnina",
         "--diffdock-variant", "diffdock_smina",
         "--equibind-variant", "equibind_unguided_gnina",
         "--site-cluster", "threshold", "--pocket-radius", "8.0",
         "--rank-by", "consensus", "--stability-boot", "25",
         "--workers", "30", "--stats"],
    notes="The out-dir tag reads 'allposes' because --pb-valid-only is OFF for the "
          "shipped figures: the cluster geometry is measured over every pose."))

reg.add(Stage(
    name="bench_pandamap", section="3. Benchmark analysis",
    title="PandaMap interaction fingerprints",
    determinism=BITEXACT, cost=MODERATE,
    needs=["bench_posebusters_matched"],
    outputs=[f"{PMAP}/pandamap_interactions.csv"],
    thesis="Table 4; Figures 7, 39",
    cmd=[VINA_PY, ANA / "run_pandamap.py", "-c", CFG / "07_benchmark_pandamap.yaml"]))

reg.add(Stage(
    name="bench_pandamap_report", section="3. Benchmark analysis",
    title="Interaction-difference report",
    determinism=BITEXACT, cost=CHEAP,
    needs=["bench_pandamap", "bench_pose_comparison"],
    outputs=[f"{PMAP}/report/05_fingerprint_similarity_top5.png",
             f"{PMAP}/report/04c_native_recovery_by_rank.png"],
    thesis="Table 4; Figure 7",
    cmd=[VINA_PY, ANA / "pandamap_interaction_report.py",
         "-c", CFG / "07_benchmark_pandamap.yaml", "--in-dir", str(PMAP),
         "--per-pose-metrics", f"{BENCH_REPORT}/per_pose_metrics.csv",
         "--exclude-methods", "unidock2,autodock_gnina"]))

reg.add(Stage(
    name="bench_interaction_audit", section="3. Benchmark analysis",
    title="Interaction pose-basis audit",
    determinism=BITEXACT, cost=CHEAP,
    needs=["bench_pandamap_report"],
    outputs=[f"{PMAP}/report/interaction_pose_basis_audit.txt",
             f"{PMAP}/report/interaction_pose_basis_audit.csv"],
    thesis='App. B "Interaction Fingerprints" quotes this output verbatim',
    cmd=[VINA_PY, ANA / "interaction_pose_basis_audit.py",
         "--config", CFG / "07_benchmark_pandamap.yaml", "--in-dir", str(PMAP),
         "--metrics", f"{BENCH_REPORT}/per_pose_metrics.csv",
         "--ids-file", f"{BENCH_REPORT}/analysed_cohort_ids.txt",
         "--exclude-methods", "unidock2,autodock_gnina"],
    notes="--ids-file is REQUIRED and is the whole difficulty of this stage. The "
          "PandaMap arm profiles 306 complexes while the analysed cohort is 303, so "
          "an unrestricted run leaves nine EquiBind poses in three excluded "
          "complexes without an RMSD and the join guard aborts. The guard is right: "
          "an unguarded gap would bin those poses into a nan band and contaminate "
          "the standardisation reference. The cohort file is derived from "
          "per_pose_metrics rather than written by hand."))
print(f"{len(reg)} stages registered")
''')

code(r'''
EFF_C = CANONICAL["benchmark_effort_charged"].relative_to(ROOT)
EFF_E = CANONICAL["benchmark_effort_elapsed"].relative_to(ROOT)
EFFORT_COMMON = [
    "--dataset", "benchmark",
    "--autodock-dir", "Dockings/vina_results_full_protein_vina_scoring_mgltools_exh128",
    "--autodock-prep", "mgl_tools", "--autodock-refine", "gnina", "--autodock-gnina-gpu",
    "--autodock-method", "autodock_mgltools_exh128_gnina",
    "--autodock-optimizer-workers", "16",
    "--equibind-dir", "Dockings/Benchmark_Equibind_cputimed",
    "--unidock2-dir", "", "--unidock-dir", "",
    "--per-pose-csv", f"{BENCH_REPORT}/per_pose_metrics.csv",
]

reg.add(Stage(
    name="bench_effort_charged", section="3. Benchmark analysis",
    title="Cost per qualifying pose, charged basis",
    determinism=BITEXACT, cost=CHEAP,
    needs=["bench_pose_comparison"],
    outputs=[f"{EFF_C}/effort_summary.csv",
             f"{EFF_C}/panels/effort_by_quality_near2.png"],
    thesis="Table 6; Figure 15",
    cmd=[VINA_PY, ANA / "docking_effort_comparison.py", *EFFORT_COMMON,
         "--out-dir", str(EFF_C), "--basis", "charged", "--cpu-threads", "32"],
    notes="The charged basis divides cpu_core_s/32 + gpu_s uniformly. The default "
          "elapsed basis divides wall_s, which is not one quantity across arms and "
          "is sensitive to --autodock-optimizer-workers; charged is not."))

reg.add(Stage(
    name="bench_effort_elapsed", section="3. Benchmark analysis",
    title="Hardware resource per near-native valid pose",
    determinism=BITEXACT, cost=CHEAP,
    needs=["bench_pose_comparison"],
    outputs=[f"{EFF_E}/panels/resource_per_near_native_valid_pose.png"],
    thesis="Table 10; Figure 16",
    cmd=[VINA_PY, ANA / "docking_effort_comparison.py", *EFFORT_COMMON,
         "--out-dir", str(EFF_E)]))

reg.add(Stage(
    name="bench_exhaustiveness_returns", section="3. Benchmark analysis",
    title="Exhaustiveness ladder: marginal return per hour",
    determinism=BITEXACT, cost=MODERATE,
    needs=["bench_pose_comparison", "bench_autodock_exh18_raw", "bench_autodock_exh92_raw"],
    outputs=["posebusters_results/autodock_exhaustiveness_returns/figures/exh_raw_vs_rescored.png",
             "posebusters_results/autodock_exhaustiveness_returns/figures/exh_depth_sweep.png",
             "posebusters_results/autodock_exhaustiveness_returns/figures/exh_yield_vs_cost.png",
             "posebusters_results/autodock_exhaustiveness_returns/figures/exh_marginal_return.png"],
    thesis="Table 8; Figures 18-21",
    cmd=[VINA_PY, ANA / "autodock_exhaustiveness_returns.py",
         "--per-pose-csv", f"{BENCH_REPORT}/per_pose_metrics.csv",
         "--out-dir", "posebusters_results/autodock_exhaustiveness_returns",
         "--n-boot", "5000", "--verify-input-parity"],
    notes="5,000 paired complex-level bootstrap resamples at a fixed seed; numerator "
          "and denominator recomputed on the same resample."))

for what, script, out in [
        ("ligand", "ligand_docking_difficulty.py", "ligand_difficulty"),
        ("receptor", "receptor_docking_difficulty.py", "receptor_difficulty")]:
    reg.add(Stage(
        name=f"bench_{what}_difficulty", section="3. Benchmark analysis",
        title=f"Docking difficulty by {what} attributes",
        determinism=BITEXACT, cost=CHEAP,
        needs=["bench_pose_comparison"],
        outputs=[f"PoseBusters_Benchmark_Analysis/{out}"],
        thesis="Discussion; the appendix difficulty analysis",
        cmd=[VINA_PY, ANA / script,
             "--per-pose-metrics", f"{BENCH_REPORT}/per_pose_metrics.csv",
             "--features", "PoseBusters_Benchmark_Analysis/ligand_protein_features.csv",
             "--out-dir", f"PoseBusters_Benchmark_Analysis/{out}",
             "--diffdock-variant", "diffdock_smina",
             "--exclude-preset", "meeko", "--exclude-methods", LADDER_EXCL]))

reg.add(Stage(
    name="bench_diffdock_rerank", section="3. Benchmark analysis",
    title="DiffDock re-ranking counterfactual",
    determinism=BITEXACT, cost=CHEAP,
    needs=["bench_pose_comparison"],
    outputs=["PoseBusters_Benchmark_Analysis/smina_rerank",
             "PoseBusters_Benchmark_Analysis/gnina_rerank"],
    thesis="the appendix note that smina does not re-rank DiffDock",
    cmd=[VINA_PY, ANA / "diffdock_gnina_rerank_analysis.py", "--tool", "smina",
         "--per-pose-metrics", f"{BENCH_REPORT}/per_pose_metrics.csv",
         "--out-dir", "PoseBusters_Benchmark_Analysis/smina_rerank"],
    notes="The --per-pose-metrics default still names the pre-whole-protein table, "
          "so it must be passed explicitly."))
print(f"{len(reg)} stages registered")
''')

# =============================================================================
md(r"""
## 4. Orai1 staging

Both Orai panels dock the same four receptor frames. None of this can be
regenerated.

`Orai1WT-START-Fr0` went through an unguarded OpenMM minimisation with a Langevin
integrator. Two preparations of it differ by 0.394 Å all-atom, four conformers
of it already exist on disk, and a fifth would destroy the provenance chain
linking the docked poses to a specific receptor. The AutoDock driver pins its
sha256 for exactly this reason.

The prepared ligand PDBQTs carry the same Meeko-version hazard as the benchmark
ones.
""")

code(r'''
reg.add(Stage(
    name="orai_receptors", section="4. Orai staging",
    title="Four Orai1 receptor frames (COPY ONLY)",
    determinism=NEVER, cost=MODERATE,
    outputs=["Data/Receptors/Orai1WT-START-Fr0.pdb",
             "Data/Receptors/Orai1WT-MDSnap-Fr300.pdb",
             "Data/Receptors/Orai1WT-MDSnap-Fr400.pdb",
             "Data/Receptors/Orai1WT-MDSnap-Fr499.pdb"],
    thesis="App. D Orai1 receptor model; App. B Orai1 panels",
    notes="Fr0 passed through an unguarded OpenMM minimisation and is not "
          "bit-reproducible (0.394 A between two preps). run_orai_mgltools_arm.py "
          "sha256-pins it to 6b3ab996... Never re-prepare; copy."))

reg.add(Stage(
    name="orai_ligands", section="4. Orai staging",
    title="Prepared Orai ligand PDBQT (COPY ONLY)",
    determinism=NEVER, cost=CHEAP,
    outputs=["Data/Ligands/JKU/pdbqt/*.pdbqt",
             "Dockings/Orai_Benchmark_MGLTools_exh128/_staging/ligands/pdbqt/*.pdbqt"],
    min_matches=4,
    thesis="App. B: the three modulators were docked neutral in every reported arm",
    notes="Same Meeko-version hazard as the benchmark ligands."))
print(f"{len(reg)} stages registered")
''')

# =============================================================================
md(r"""
## 5. Orai1 experimental panel

Three modulators against four frames, at the benchmark's sampling budget: thirty
Vina modes in a 6 kcal/mol window, thirty DiffDock samples, thirty unguided
EquiBind conformers.

This is the one DiffDock arm that is reproducible. Its run logs carry the seed
marker; the benchmark and control ones do not.

The Fr0 correction below applies to **AutoDock rows only**. AutoDock docked Fr0
into the minimised receptor while the main screen validates against the raw PDB,
which costs Fr0 its validity. DiffDock and EquiBind docked the raw receptor, so
extending the correction to them creates the mismatch instead of repairing it.
The quarantined `.WRONG-RECEPTOR` directories on disk are that mistake, made once.
""")

code(r'''
OEXP = CANONICAL["orai_experimental"].relative_to(ROOT)

reg.add(Stage(
    name="orai_exp_autodock", section="5. Orai experimental",
    title="AutoDock Vina exh128 + gnina, 30 modes",
    determinism=BITEXACT, cost=EXPENSIVE,
    needs=["orai_receptors", "orai_ligands"],
    outputs=["Dockings/Orai_JKU_MGLTools_exh128/mgl_tools"],
    thesis="Table 5; Figures 8-14",
    cmd=[VINA_PY, DOCK / "run_orai_mgltools_arm.py",
         "-c", CFG / "20_orai_jku_autodock_exh128_gnina.yaml"],
    notes="run_orai_mgltools_arm.py is the only legal driver. run_autodock.py main() "
          "re-prepares Fr0 and forces Meeko ligands."))

reg.add(Stage(
    name="orai_exp_diffdock", section="5. Orai experimental",
    title="DiffDock-L, 30 samples (SEEDED)",
    determinism=BITEXACT, cost=EXPENSIVE,
    needs=["orai_receptors", "orai_ligands"],
    outputs=["Dockings/diffdock_results"],
    thesis="Table 5; Figures 8-14",
    cmd=[VINA_PY, DOCK / "run_diffdock.py", "-c", CFG / "21_orai_jku_diffdock.yaml"],
    notes="24 stored run logs carry the seed marker, so this arm alone can be redrawn."))

reg.add(Stage(
    name="orai_exp_equibind", section="5. Orai experimental",
    title="EquiBind unguided, 30 conformers",
    determinism=BITEXACT, cost=EXPENSIVE,
    needs=["orai_receptors", "orai_ligands"],
    outputs=["Dockings/equibind_results_uffoff"],
    thesis="Table 5; Figures 8-14",
    notes="Config 22. uff_minimize false and output_dir repointed to the uffoff tree."))

reg.add(Stage(
    name="orai_exp_posebusters", section="5. Orai experimental",
    title="PoseBusters validity screen",
    determinism=BITEXACT, cost=MODERATE,
    needs=["orai_exp_autodock", "orai_exp_diffdock", "orai_exp_equibind"],
    outputs=[f"{OEXP}/dock/posebusters_filtered_results.csv"],
    thesis="Table 5; Figures 8-10",
    cmd=[VINA_PY, PB / "run_posebusters.py", "-c", CFG / "23_orai_jku_posebusters.yaml"],
    notes="variant_filter is load-bearing: raw ADFRsuite PDBQT carries no REMARK "
          "SMILES, so an unfiltered run falls back to Open Babel and fails silently "
          "on about 63 per cent of poses."))

reg.add(Stage(
    name="orai_exp_posebusters_fr0", section="5. Orai experimental",
    title="Fr0 receptor correction, AUTODOCK ROWS ONLY",
    determinism=BITEXACT, cost=MODERATE,
    needs=["orai_exp_posebusters"],
    outputs=["posebusters_results/orai_jku_mgltools_exh128_fr0corrected"],
    thesis="the corrected Fr0 validity figures",
    cmd=[VINA_PY, PB / "run_posebusters.py",
         "-c", CFG / "24_orai_jku_posebusters_fr0corrected.yaml"],
    notes="AutoDock only. Extending it to DiffDock or EquiBind drops EquiBind "
          "validity from 58.2 to 15.8 per cent, which is the mismatch, not the fix."))

reg.add(Stage(
    name="orai_exp_pandamap", section="5. Orai experimental",
    title="PandaMap interaction fingerprints",
    determinism=BITEXACT, cost=MODERATE,
    needs=["orai_exp_posebusters"],
    outputs=[str(CANONICAL["orai_pandamap_experimental"].relative_to(ROOT))],
    thesis="Table 12; Figures 13, 14, 25, 26",
    cmd=[VINA_PY, ANA / "run_pandamap.py", "-c", CFG / "25_orai_jku_pandamap.yaml"],
    notes="PandaMap 4.1.0 cannot see chlorine or bromine, and a renumbered receptor "
          "yields zero fingerprints."))
print(f"{len(reg)} stages registered")
''')

# =============================================================================
md(r"""
## 6. Orai1 control panel

All 308 benchmark ligands against the same four frames, 1,232 units, cut back to
ten Vina modes, ten DiffDock samples and ten unguided EquiBind conformers.

The panel exists to answer whether the experimental result is a property of the
three modulators or of the receptor. It shares the four frames, the whole-receptor
boxes, the exhaustiveness, the energy window, the preparation and the AutoDock
seed with the experimental panel. Two differences remain and both are disclosed:
the pose-selection rule, and this panel's DiffDock sampler being unseeded.
""")

code(r'''
OCTL = CANONICAL["orai_control"].relative_to(ROOT)

reg.add(Stage(
    name="orai_ctl_autodock", section="6. Orai control",
    title="AutoDock Vina exh128 + gnina, 10 modes",
    determinism=BITEXACT, cost=EXPENSIVE,
    needs=["orai_receptors", "orai_ligands"],
    outputs=["Dockings/Orai_Benchmark_MGLTools_exh128/mgl_tools"],
    thesis="Table 5; Figures 8-14",
    cmd=[VINA_PY, DOCK / "run_orai_mgltools_arm.py",
         "-c", CFG / "30_orai_benchmark_autodock_exh128_gnina.yaml"]))

reg.add(Stage(
    name="orai_ctl_diffdock", section="6. Orai control",
    title="DiffDock-L, 10 samples",
    determinism=VERIFY, cost=EXPENSIVE,
    needs=["orai_receptors", "orai_ligands"],
    outputs=["Dockings/Orai_Benchmark_DiffDock"],
    thesis="Table 5; Figures 8-14",
    notes="UNSEEDED, like the benchmark run. App. B names this as the one remaining "
          "input difference between the two Orai panels."))

reg.add(Stage(
    name="orai_ctl_equibind", section="6. Orai control",
    title="EquiBind unguided, 10 conformers",
    determinism=BITEXACT, cost=EXPENSIVE,
    needs=["orai_receptors", "orai_ligands"],
    outputs=["Dockings/Orai_Benchmark_Equibind"],
    thesis="Table 5; Figures 8-14",
    notes="Config 32. This is the arm the old notebook already read strictly, with a "
          "cfg_req() helper that raises on a missing key."))

reg.add(Stage(
    name="orai_ctl_posebusters", section="6. Orai control",
    title="PoseBusters validity screen",
    determinism=BITEXACT, cost=EXPENSIVE,
    needs=["orai_ctl_autodock", "orai_ctl_diffdock", "orai_ctl_equibind"],
    outputs=[f"{OCTL}/dock/posebusters_filtered_results.csv"],
    thesis="Table 5; Figures 8-10",
    cmd=[VINA_PY, PB / "run_posebusters.py", "-c", CFG / "33_orai_benchmark_posebusters.yaml"]))

for tag, cfg_name, out in [
        ("", "34_orai_benchmark_posebusters_fr0corrected.yaml",
         "posebusters_results/orai_benchmark_fr0corrected"),
        ("_gnina", "35_orai_benchmark_posebusters_gnina_fr0corrected.yaml",
         "posebusters_results/orai_benchmark_gnina_fr0corrected")]:
    reg.add(Stage(
        name=f"orai_ctl_posebusters_fr0{tag}", section="6. Orai control",
        title=f"Fr0 correction{' (gnina variant)' if tag else ''}, AUTODOCK ROWS ONLY",
        determinism=BITEXACT, cost=MODERATE,
        needs=["orai_ctl_posebusters"],
        outputs=[out],
        thesis="Fr0 validity 73.9 -> 97.7 per cent; the AutoDock control yield",
        cmd=[VINA_PY, PB / "run_posebusters.py", "-c", CFG / cfg_name]))

for panel, tag, cfg_name, out in [
        ("ctl", "control", "37_orai_benchmark_posebusters_equibind_minimize.yaml",
         "posebusters_results/orai_benchmark_equibind_minimize"),
        ("exp", "experimental", "38_orai_jku_posebusters_equibind_minimize.yaml",
         "posebusters_results/orai_jku_equibind_minimize")]:
    reg.add(Stage(
        name=f"orai_{panel}_posebusters_equibind_minimize",
        section="6. Orai control" if panel == "ctl" else "5. Orai experimental",
        title=f"PoseBusters, matched-EquiBind arm ({tag} panel)",
        determinism=BITEXACT, cost=MODERATE,
        needs=[f"orai_{panel}_posebusters"],
        outputs=[out],
        thesis="the EquiBind rows of Table 5 and of every Orai figure",
        cmd=[VINA_PY, PB / "run_posebusters.py", "-c", CFG / cfg_name],
        notes="ADDED 2026-09-02. This arm had no stage, so the canonical "
              "_orai_matched_root trees could not be rebuilt from the declared stages. "
              "Configs 23 and 33 screen the PRE-matched EquiBind trees; this screens the "
              "--minimize re-refinement whose poses the canonical tables actually cite, "
              "24,640 references on the control panel and 1,440 on the experimental one. "
              "Both screens are needed and are not alternatives."))

reg.add(Stage(
    name="orai_ctl_pandamap", section="6. Orai control",
    title="PandaMap interaction fingerprints",
    determinism=BITEXACT, cost=EXPENSIVE,
    needs=["orai_ctl_posebusters"],
    outputs=[str(CANONICAL["orai_pandamap_control"].relative_to(ROOT))],
    thesis="Table 12; Figures 13, 14, 25, 26",
    cmd=[VINA_PY, ANA / "run_pandamap.py", "-c", CFG / "36_orai_benchmark_pandamap.yaml"]))
print(f"{len(reg)} stages registered")
''')

# =============================================================================
md(r"""
## 7. Orai1 cross-panel analysis

The transmembrane exclusion runs once over both panels, which is the fix for a
real bug in the old notebook: the JKU-only version of this cell sat in the
experimental section and wrote the same file, so whichever ran last won.

Ordering here is fixed by data, not preference. The transmembrane classification
feeds everything. The share comparison must precede the ligand-level contrasts,
which read its per-unit table and are the last thing in the whole chain.

Every published Orai figure uses the gnina AutoDock arm. The flag is stated
explicitly on every stage below, including where the default is already right,
because a default is not a record.
""")

code(r'''
OROOT = CANONICAL["orai_root"].relative_to(ROOT)
TMSHARE = CANONICAL["orai_tm_share"].relative_to(ROOT)
VARIANTS = ["--autodock-variant", "gnina", "--diffdock-variant", "smina",
            "--equibind-variant", "gnina"]

reg.add(Stage(
    name="orai_tm_exclusion", section="7. Orai cross-panel",
    title="Transmembrane pore exclusion, BOTH panels in one run",
    determinism=BITEXACT, cost=MODERATE,
    needs=["orai_exp_posebusters_fr0", "orai_ctl_posebusters_fr0"],
    outputs=[f"{OEXP}/transmembrane_filter/tm_pose_classification.csv",
             f"{OCTL}/transmembrane_filter/tm_pose_classification.csv"],
    thesis="Table 5; Figures 9, 10",
    cmd=[VINA_PY, ANA / "orai_transmembrane_exclusion.py", *VARIANTS,
         "--out-root", str(OROOT)],
    notes="The pore slab runs R91 to E106. In the old notebook a JKU-only copy of "
          "this step sat in the experimental section and wrote the same output, so "
          "which one won depended on execution order."))

for panel, key, per_pose in [("exp", "orai_experimental", OEXP), ("ctl", "orai_control", OCTL)]:
    reg.add(Stage(
        name=f"orai_{panel}_clusters", section="7. Orai cross-panel",
        title=f"Reference-free pose clustering ({'experimental' if panel == 'exp' else 'control'})",
        determinism=BITEXACT, cost=MODERATE,
        needs=["orai_tm_exclusion"],
        outputs=[f"{per_pose}/pose_clusters/per_pair.csv",
                 f"{per_pose}/pose_clusters/cluster_quality_per_tool.csv"],
        thesis="Table 3 analogue for Orai; Figures 11, 12, 23, 24",
        cmd=[VINA_PY, ANA / "orai_pose_cluster_report.py",
             "--per-pose-csv", f"{per_pose}/dock/posebusters_filtered_results.no_tm.csv",
             "--out-dir", f"{per_pose}/pose_clusters",
             "--autodock-variant", "gnina", "--diffdock-variant", "smina",
             "--site-cluster", "threshold", "--pocket-radius", "8.0",
             "--cluster-quality", "--top-n-poses", "10"],
        notes="Clusters are built on the pooled three-tool cloud, so changing one "
              "tool's pose set moves the other two tools' cluster values as well."))

reg.add(Stage(
    name="orai_yield_compare", section="7. Orai cross-panel",
    title="Figure 8: PoseBusters-valid yield, control against experimental",
    determinism=BITEXACT, cost=CHEAP,
    needs=["orai_tm_exclusion"],
    outputs=[f"{CANONICAL['orai_yield'].relative_to(ROOT)}/"
             "09f_pbvalid_yield_dominant_orai_benchmark_vs_experimental.png"],
    thesis="Figure 8 (media/media/image13.png)",
    cmd=[VINA_PY, ANA / "orai_pbvalid_yield_compare.py", "--results-root", str(OROOT)],
    notes="Extracted 2026-09-02 from an inline notebook cell. Until then this was the "
          "only float in the Results chapter with no committed generator. Verified "
          "byte-identical to the shipped asset."))

reg.add(Stage(
    name="orai_tm_share_compare", section="7. Orai cross-panel",
    title="Validity and transmembrane share, control against experimental",
    determinism=BITEXACT, cost=CHEAP,
    needs=["orai_tm_exclusion"],
    outputs=[f"{TMSHARE}/pbvalid_tm_share_pooled.csv",
             f"{TMSHARE}/pbvalid_tm_share_per_unit.csv",
             f"{TMSHARE}/fig_pbvalid_outside_tm_whisker.png",
             f"{TMSHARE}/fig_tm_loss_relative_compare.png"],
    thesis="Table 5; Figures 9, 10",
    cmd=[VINA_PY, ANA / "orai_pbvalid_tm_share_compare.py", *VARIANTS,
         "--results-root", str(OROOT), "--out-dir", str(TMSHARE), "--top-n-poses", "10"]))

reg.add(Stage(
    name="orai_xtool_agreement", section="7. Orai cross-panel",
    title="Cross-tool agreement and reference-free cluster quality",
    determinism=BITEXACT, cost=CHEAP,
    needs=["orai_exp_clusters", "orai_ctl_clusters"],
    outputs=[f"{OEXP}/pose_clusters/panels/orai_cross_tool_agreement_compare.png",
             f"{OEXP}/pose_clusters/panels/orai_tool_agreement_compare.png"],
    thesis="Figures 11, 12, 23",
    cmd=[VINA_PY, ANA / "orai_cross_tool_agreement_compare.py",
         "--exp-csv", f"{OEXP}/pose_clusters/per_pair.csv", "--exp-label", "Exp. Ligands",
         "--bench-csv", f"{OCTL}/pose_clusters/per_pair.csv",
         "--bench-label", "Benchmark ligands x Orai",
         "--out-dir", f"{OEXP}/pose_clusters/panels"],
    notes="This script takes no --autodock-variant. Its arm is fixed by whichever run "
          "wrote the two per_pair.csv inputs, which is why they are named explicitly."))

reg.add(Stage(
    name="orai_pandamap_compare", section="7. Orai cross-panel",
    title="Interaction compare: totals, type profile, hot-spots, overlap",
    determinism=BITEXACT, cost=CHEAP,
    needs=["orai_exp_pandamap", "orai_ctl_pandamap"],
    outputs=[f"{CANONICAL['orai_pandamap_compare'].relative_to(ROOT)}/fig_type_profile_compare.png",
             f"{CANONICAL['orai_pandamap_compare'].relative_to(ROOT)}/fig_fingerprint_overlap_compare.png"],
    thesis="Table 12; Figures 13, 14, 25, 26",
    cmd=[VINA_PY, ANA / "orai_pandamap_interaction_compare.py",
         "--exp-dir", str(CANONICAL["orai_pandamap_experimental"].relative_to(ROOT)),
         "--bench-dir", str(CANONICAL["orai_pandamap_control"].relative_to(ROOT)),
         "--out-dir", str(CANONICAL["orai_pandamap_compare"].relative_to(ROOT)),
         *VARIANTS, "--top-n-poses", "10",
         "--exp-posebusters-csv", f"{OEXP}/dock/posebusters_filtered_results.csv",
         "--bench-posebusters-csv", f"{OCTL}/dock/posebusters_filtered_results.csv"],
    notes="The two --*-posebusters-csv paths are REQUIRED and their defaults are wrong "
          "for this tree. EquiBind's PandaMap pose_rank is a 999 placeholder, so the "
          "top-ten cap has to read the ranks from the PoseBusters CSV. With the default "
          "paths, which name the plain trees, the lookup fails and all 324 EquiBind pose "
          "rows are dropped without an error. The four figures then render without an "
          "EquiBind series and differ from the published ones."))

reg.add(Stage(
    name="orai_ligand_contrasts", section="7. Orai cross-panel",
    title="Ligand-level sensitivity re-test (LAST in the chain)",
    determinism=BITEXACT, cost=CHEAP,
    needs=["orai_tm_share_compare", "orai_exp_clusters", "orai_ctl_clusters",
           "orai_pandamap_compare"],
    outputs=[f"{TMSHARE}/orai_ligand_level_contrasts.csv",
             f"{TMSHARE}/orai_ligand_level_contrasts.txt"],
    thesis="the Methods commitment that every Orai cross-panel contrast is re-tested "
           "at ligand level in a sensitivity analysis reported with its result",
    cmd=[VINA_PY, ANA / "orai_ligand_level_contrasts.py", *VARIANTS,
         "--per-unit-csv", f"{TMSHARE}/pbvalid_tm_share_per_unit.csv",
         "--out-dir", str(TMSHARE)],
    notes="The effective number of independent experimental observations is three "
          "ligands, so this is what keeps the frame-ligand tests from being read as "
          "more powered than they are."))

reg.add(Stage(
    name="orai_region_consensus", section="7. Orai cross-panel",
    title="Consensus binding-region analysis",
    determinism=BITEXACT, cost=MODERATE,
    needs=["orai_tm_exclusion"],
    outputs=[f"{OCTL}/region_consensus"],
    thesis="the Orai1 interpretation section",
    cmd=[VINA_PY, ANA / "orai_ligand_region_consensus.py",
         "--pb-csv", f"{OCTL}/dock/posebusters_filtered_results.no_tm.csv",
         "--out-dir", f"{OCTL}/region_consensus"],
    notes="Label-shuffle permutation null, 2,000 iterations at a fixed seed."))
print(f"{len(reg)} stages registered")
''')

# =============================================================================
md(r"""
## 8. Dataset chapter

The chemical-space figures come from a separate notebook, executed in place. Its
rendering is not quite deterministic: re-running moves roughly 1 to 2.5 per cent
of pixels and changes one canvas by three pixels. The underlying statistics are
stable; only the raster differs.
""")

code(r'''
reg.add(Stage(
    name="dataset_figures", section="8. Dataset chapter",
    title="Chemical-space figures and descriptor tables",
    determinism=BITEXACT, cost=MODERATE,
    outputs=["PoseBusters_Benchmark_Analysis/figures/01_ligand_univariate_distributions.png",
             "PoseBusters_Benchmark_Analysis/figures/09_pca_scree_loadings.png",
             "PoseBusters_Benchmark_Analysis/summary_statistics.csv"],
    thesis="Tables 16, 17; Figures 28-37",
    cmd=[VINA_PY, "-m", "jupyter", "nbconvert", "--to", "notebook", "--execute",
         "--inplace", "--ExecutePreprocessor.timeout=1800",
         "PoseBusters_DataSet_Analysis.ipynb"],
    notes="Descriptors are computed on a hydrogen-suppressed molecule; two of the "
          "sixteen are basis-dependent. Figure-to-file mapping is not sequential: "
          "image39 is the scree plot, image41 the receptor profile."))
print(f"{len(reg)} stages registered\n")

for section, stages in reg.by_section().items():
    print(f"{section}: {len(stages)} stages")
''')

# =============================================================================
md(r"""
## 9. Run the pipeline

The driver walks the registry in declaration order, which it first verifies is a
valid topological order. For each stage it reports whether the outputs are
present and what it did.

In `verify` mode nothing is executed. A stage marked `MISSING` is a real gap, not
a run that has yet to happen.
""")

code(r'''
driver = Driver(reg, mode=RUN_MODE, only=FORCE_STAGES)
results = driver.walk()
print(driver.table())
''')

code(r'''
missing = [r for r in results if r["status"] == "MISSING"]
if not missing:
    print("Every stage's outputs are present.")
else:
    print(f"{len(missing)} stage(s) with missing outputs:\n")
    for r in missing:
        print(f"  {r['stage']}")
        print(f"      {r['title']}")
        print(f"      missing: {reg.stages[r['stage']].satisfied()[1]}")
        if reg.stages[r['stage']].notes:
            print(f"      note   : {reg.stages[r['stage']].notes}")
        print()
''')

code(r'''
manifest = H.write_manifest(ROOT / "posebusters_results/_reproduction/run_manifest.json",
                            reg, checks, results)
print(f"run manifest -> {manifest.relative_to(ROOT)}")
print(f"  {len(reg)} stages, {sum(1 for r in results if r['status'] == 'present')} satisfied")
''')

# =============================================================================
md(r"""
## 10. Assert the thesis numbers

The point of the whole exercise. Each entry in
`Scripts/Analysis/thesis_expected_values.yaml` names a number the thesis prints,
where in the document it appears, and how to recompute it from the canonical
CSVs. This cell recomputes them all and prints a pass or fail per row.

Ranking rules are the verified ones and differ per tool, which is easy to get
wrong: AutoDock by `optimized_rank` on the gnina arm, DiffDock by `rank` on the
smina variant, EquiBind by ascending `gnina_affinity`, because EquiBind emits no
confidence of its own.

A failure here is reported, not corrected. If a number does not reproduce, that
is a finding about the thesis and belongs in a conversation, not in a silent edit.
""")

code(r'''
import thesis_assertions as TA

report = TA.run_all(verbose=True)
print()
print(report.summary())
''')

code(r'''
if report.failures:
    print("NUMBERS THAT DID NOT REPRODUCE\n")
    for f in report.failures:
        print(f"  {f.key}")
        print(f"      thesis   : {f.expected}   ({f.source})")
        print(f"      recomputed: {f.actual}")
        print(f"      {f.note}")
        print()
else:
    print("Every checked number reproduces from the canonical trees.")
''')

# =============================================================================
md(r"""
## 11. Byte-level regeneration check

Section 10 proves PROVENANCE: the image printed in the thesis is byte-identical to
a file on disk. That is not the same as proving the pipeline still produces it.

This cell proves REGENERATION. It re-runs each figure's generator from the
canonical inputs, into a scratch directory that is thrown away, and compares the
bytes to the published asset. Nothing it does can touch a canonical output.

It is off by default because it takes about ten minutes. When skipped, the cell
reads back the last recorded run and re-fingerprints its inputs, so it reports
either that the result still holds or exactly what changed since. A cached verdict
cannot outlive the code it describes.

Set `RUN_REGEN_CHECK` to `True` to run it, or run the module directly:

```
python Scripts/Analysis/validate_regeneration.py           # about 10 minutes
python Scripts/Analysis/validate_regeneration.py --quick   # about 4 minutes
```

This check earns its keep. It caught a real defect in the Orai interaction-compare
stage: two arguments whose defaults name the plain trees, without which the
top-ten cap cannot rank EquiBind, because EquiBind's PandaMap pose rank is a
placeholder. All 324 EquiBind pose rows were dropped in silence and four figures
rendered without an EquiBind series. Nothing errored. Only the byte comparison
found it.
""")

code(r'''
# About ten minutes. Set True to run, or invoke the module directly.
RUN_REGEN_CHECK = False

if RUN_REGEN_CHECK:
    import validate_regeneration
    rc = validate_regeneration.main([])
    print("\nall regenerated assets match" if rc == 0 else "\nSOME ASSETS DIFFER, see above")
else:
    import validate_regeneration as VR
    # Reads the recorded result and re-fingerprints the inputs, so a cached
    # verdict cannot outlive the code it describes.
    print(VR.describe_result())
    print()
    print("Figures declared rather than tested, each with its reason:")
    for figs, what, why in VR.DECLARED:
        print(f"  Figure {figs:9s} {what}")
''')

# =============================================================================
md(r"""
## 12. Regenerate the guide

Appendix A of the thesis promises a guide at `Scripts/Analysis/REGENERATE.md` that
maps every reported float to the output it came from and to the command that
rebuilds it, names the environment each command needs, gives the order dependent
steps run in, and states explicitly what could not be reconstructed.

That guide used to be maintained by hand, and it drifted. Its figure numbers came
from a build that has since been retired and its paths predated a migration of the
output trees, so following it landed you on the wrong file for nearly every
figure. Nothing announced the drift, because a hand-written index cannot notice
that the tree beneath it moved.

So the guide is written from this registry instead. The ordering, the commands and
the determinism classes come from the same stage objects the driver above just
executed, and the figure index is rebuilt by checksumming each shipped asset
against every candidate output on disk. A figure whose source moves shows up as
moved rather than as a stale path.

The old hand-written guide is kept as `REGENERATE_NOTES.md` for its per-command
commentary and its traps, which are not derivable from the registry.
""")

code(r'''
import regenerate_guide

guide = regenerate_guide.build(reg)
print()
print("Appendix A promises this file. It now contains, in order:")
for section in ("Environments and binaries", "Canonical trees", "Order", "Stages",
                "Figure index", "Table index", "What cannot be reconstructed"):
    print(f"  - {section}")
''')

# =============================================================================
nb = {
    "cells": cells,
    "metadata": {
        "kernelspec": {"display_name": "vina", "language": "python", "name": "python3"},
        "language_info": {"name": "python", "version": "3.12.12"},
    },
    "nbformat": 4,
    "nbformat_minor": 5,
}
OUT.write_text(json.dumps(nb, indent=1))
print(f"wrote {OUT} with {len(cells)} cells "
      f"({sum(1 for c in cells if c['cell_type'] == 'code')} code, "
      f"{sum(1 for c in cells if c['cell_type'] == 'markdown')} markdown)")
