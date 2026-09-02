#!/usr/bin/env python3
"""Generate Scripts/Docking/configs_thesis/ from the current thesis configs.

Each output is the source file copied VERBATIM with a provenance header prepended,
except where a `bake` mapping records an override the notebook used to apply
in memory. Parsed equality is asserted for every file.
"""
from __future__ import annotations

import shutil
import sys
from pathlib import Path

import yaml

ROOT = Path("/home/manndo/master_dev")
DOCK = ROOT / "Scripts/Docking"
PB = DOCK / "Posebusters"
ANA = ROOT / "Scripts/Analysis"
OUT = DOCK / "configs_thesis"

HDR = "# " + "=" * 76


def header(**kw) -> str:
    lines = [HDR, f"# THESIS REPRODUCTION CONFIG - {kw['arm']}", "#"]
    for label in ("thesis", "supports", "dataset", "driver", "writes", "canonical"):
        val = kw.get(label)
        if not val:
            continue
        pretty = {
            "thesis": "Thesis section",
            "supports": "Supports",
            "dataset": "Dataset",
            "driver": "Driver",
            "writes": "Writes",
            "canonical": "Analysis tree",
        }[label]
        first, *rest = val.split("\n")
        lines.append(f"#   {pretty:<16}: {first}")
        lines.extend(f"#   {'':<16}  {r}" for r in rest)
    lines += ["#", f"#   {'Copied from':<16}: {kw['source']}"]
    dev = kw.get("deviations") or ["none - byte-identical body"]
    lines.append(f"#   {'Deviations':<16}: {dev[0]}")
    lines.extend(f"#   {'':<16}  {d}" for d in dev[1:])
    det_first, *det_rest = kw["determinism"].split("\n")
    lines += ["#", f"#   {'Determinism':<16}: {det_first}"]
    lines.extend(f"#   {'':<16}  {r}" for r in det_rest)
    lines += [
        "#",
        "#   This file is READ-ONLY to the pipeline. Thesis_Reproduction.ipynb loads it",
        "#   and never writes it back. Rationale in THESIS_REPRODUCTION.md.",
        HDR,
        "",
    ]
    # Safety net: every header line must be a YAML comment or empty.
    for ln in lines:
        assert ln == "" or ln.startswith("#"), f"uncommented header line: {ln!r}"
    return "\n".join(lines)


BITEXACT = "bit-exact - seeded, a re-run must reproduce the stored tree"
VERIFY = "verify-only - the reported run was an UNSEEDED single draw; a re-run is a new sample"
DERIVED = "bit-exact - deterministic post-processing of stored poses"

BENCH = "Calibration benchmark, 308 PoseBusters ligands (303 analysed)"
OJKU = "Orai1 experimental panel, 3 modulators x 4 receptor frames"
OBEN = "Orai1 control panel, 308 benchmark ligands x 4 receptor frames = 1,232 units"

SPECS = [
    # ---------------------------------------------------------------- benchmark
    dict(
        name="01_benchmark_autodock_exh128_raw.yaml",
        source=DOCK / "autodock_vina_docking_config_full_protein_search_vina_mgltools_exh128.yaml",
        arm="AutoDock Vina, exhaustiveness 128, raw search (benchmark)",
        thesis='App. B "AutoDock Vina"; the reported arm is exhaustiveness 128',
        supports="Table 1, Table 8 (ladder top rung), Figures 18-21",
        dataset=BENCH,
        driver="run_autodock_exhaustiveness_arm.py -c <this> --ids-file <all308>\n"
        "--expect-exhaustiveness 128 --prepared-inputs-from <exh32 tree>",
        writes="Dockings/vina_results_full_protein_vina_scoring_mgltools_exh128",
        canonical="posebusters_results/benchmark_matched_equibind",
        determinism=BITEXACT + " (Vina seed 42)",
    ),
    dict(
        name="02_benchmark_autodock_exh128_gnina.yaml",
        source=DOCK / "autodock_vina_docking_config_full_protein_search_vina_mgltools_exh128_gnina.yaml",
        arm="AutoDock Vina exh128 + gnina rescoring (benchmark, DOMINANT ARM)",
        thesis='App. B "AutoDock Vina"; gnina 1.3.2 in rescoring mode, ranked by CNNaffinity',
        supports="Tables 1, 2, 6; Figures 1-3, 5-7, 15, 16",
        dataset=BENCH,
        driver="run_autodock_exhaustiveness_arm.py -c <this> (second pass over the same output_dir)",
        writes="Dockings/vina_results_full_protein_vina_scoring_mgltools_exh128/*/optimized_gnina",
        canonical="posebusters_results/benchmark_matched_equibind",
        determinism=BITEXACT + " (Vina seed 42, optimize_seed 42)",
    ),
    dict(
        name="03_benchmark_diffdock.yaml",
        source=DOCK / "diffdock_docking_config.yaml",
        arm="DiffDock-L (benchmark)",
        thesis='App. B "DiffDock"; 20 denoising steps, 30 samples, blind whole receptor',
        supports="Tables 1, 2, 6; Figures 1-3, 5-7, 15, 16",
        dataset=BENCH,
        driver="run_diffdock.py -c <this>",
        writes="Dockings/Benchmark_DiffDock",
        canonical="posebusters_results/benchmark_matched_equibind",
        determinism=VERIFY + ".\n"
        "0 of the stored run logs carry the '[repro] seeded' marker, so the\n"
        "reported benchmark poses cannot be redrawn. Verified by checksum only.",
        bake={
            "output_dir": "Dockings/Benchmark_DiffDock",
            "log_dir": "Dockings/Logs/benchmark_diffdock_logs",
            "optimization": "all",
            "inference_seed": None,
        },
        deviations=[
            "output_dir, log_dir and optimization baked in. The old notebook",
            "applied these in memory (cell 22), so the YAML on disk never",
            "described the run that produced the reported poses.",
            "inference_seed set to null to RECORD the truth: the reported run",
            "was unseeded. Setting 42 here would silently claim reproducibility",
            "the stored tree does not have.",
        ],
    ),
    dict(
        name="04_benchmark_equibind.yaml",
        source=DOCK / "equibind_benchmark_config.yaml",
        arm="EquiBind, unguided and pocket-guided (benchmark)",
        thesis='App. B "EquiBind" and "Pocket-Guided Variants"; 30 ETKDGv3+MMFF conformers',
        supports="Tables 1, 2, 6; Figures 1-3, 5-7, 15, 16; the guided-arm paragraph",
        dataset=BENCH,
        driver="run_equibind.py, EQ_* environment mapped from this file",
        writes="Dockings/Benchmark_Equibind",
        canonical="posebusters_results/benchmark_matched_equibind",
        determinism=BITEXACT + " (30 fixed RDKit seeds, per-job model seeding)",
        deviations=[
            "none in value. This file REPLACES equibind_docking_config.yaml,",
            "which is marked DEPRECATED and which notebook cell 68 rewrote in",
            "place, flipping uff_minimize to true. uff_minimize is false here,",
            "which is the reported setting (App. B: the UFF stage was switched",
            "off for every run in the thesis).",
        ],
    ),
    dict(
        name="05_benchmark_posebusters.yaml",
        source=PB / "posebusters_benchmark_full_protein_config.yaml",
        arm="PoseBusters validity screen (benchmark, all arms)",
        thesis='App. B "Post-Docking Optimisation and Validity Screening"; PoseBusters 0.6.3',
        supports="Table 1 validity block, Table 24; Figure 1",
        dataset=BENCH,
        driver="Posebusters/run_posebusters.py -c <this>",
        writes="posebusters_results/benchmark_full_protein_vina_scoring/dock",
        canonical="superseded by the matched-EquiBind pass, config 06",
        determinism=DERIVED,
    ),
    dict(
        name="06_benchmark_posebusters_matched_equibind.yaml",
        source=PB / "posebusters_benchmark_equibind_minimize_config.yaml",
        arm="PoseBusters, matched-EquiBind arm (benchmark, CANONICAL)",
        thesis="App. B; the --minimize refinement that matches EquiBind to the\n"
        "AutoDock and DiffDock refiner setting (examiner item M5)",
        supports="every benchmark number in the Results chapter",
        dataset=BENCH,
        driver="Posebusters/run_posebusters.py -c <this>",
        writes="posebusters_results/benchmark_equibind_minimize",
        canonical="posebusters_results/benchmark_matched_equibind - THE canonical benchmark tree",
        determinism=DERIVED,
    ),
    dict(
        name="07_benchmark_pandamap.yaml",
        source=ANA / "pandamap_benchmark_matched_equibind_config.yaml",
        arm="PandaMap interaction fingerprints (benchmark)",
        thesis='App. B "Interaction Fingerprints"; PandaMap 4.1.0',
        supports="Table 4; Figure 7",
        dataset=BENCH,
        driver="Analysis/run_pandamap.py -c <this>, then pandamap_interaction_report.py",
        writes="pandamap_results/benchmark_matched_equibind",
        canonical="pandamap_results/benchmark_matched_equibind",
        determinism=DERIVED,
    ),
    # ------------------------------------------------------- exhaustiveness ladder
    dict(
        name="10_ladder_exh18.yaml",
        source=DOCK / "autodock_vina_docking_config_full_protein_search_vina_mgltools_exh18.yaml",
        arm="AutoDock ladder rung, exhaustiveness 18 (raw)",
        thesis='App. B "Choice of Exhaustiveness and Rescoring"',
        supports="Table 8; Figures 18-21",
        dataset=BENCH,
        driver="run_autodock_exhaustiveness_arm.py -c <this> --expect-exhaustiveness 18",
        writes="Dockings/vina_results_full_protein_vina_scoring_mgltools_exh18",
        canonical="posebusters_results/autodock_exhaustiveness_returns",
        determinism=BITEXACT,
    ),
    dict(
        name="11_ladder_exh32_raw.yaml",
        source=DOCK / "autodock_vina_docking_config_full_protein_search_vina_mgltools.yaml",
        arm="AutoDock ladder rung, exhaustiveness 32 (raw)",
        thesis='App. B "Choice of Exhaustiveness and Rescoring"',
        supports="Table 8; Figures 18-21. Also the --prepared-inputs-from source\n"
        "for every other rung, which is why it must never be re-prepared.",
        dataset=BENCH,
        driver="run_autodock_exhaustiveness_arm.py -c <this> --expect-exhaustiveness 32",
        writes="Dockings/vina_results_full_protein_vina_scoring_mgltools",
        canonical="posebusters_results/autodock_exhaustiveness_returns",
        determinism=BITEXACT,
    ),
    dict(
        name="12_ladder_exh32_gnina.yaml",
        source=DOCK / "autodock_vina_docking_config_full_protein_search_vina_mgltools_gnina.yaml",
        arm="AutoDock ladder rung, exhaustiveness 32 + gnina rescoring",
        thesis='App. B "Choice of Exhaustiveness and Rescoring"',
        supports="Table 8; Figure 18",
        dataset=BENCH,
        driver="run_autodock_exhaustiveness_arm.py -c <this>",
        writes="Dockings/vina_results_full_protein_vina_scoring_mgltools/*/optimized_gnina",
        canonical="posebusters_results/autodock_exhaustiveness_returns",
        determinism=BITEXACT,
    ),
    dict(
        name="13_ladder_exh64_raw.yaml",
        source=DOCK / "autodock_vina_docking_config_full_protein_search_vina_mgltools_exh64.yaml",
        arm="AutoDock ladder rung, exhaustiveness 64 (raw)",
        thesis='App. B "Choice of Exhaustiveness and Rescoring"',
        supports="Table 8; Figures 18-21",
        dataset=BENCH,
        driver="run_autodock_exhaustiveness_arm.py -c <this> --expect-exhaustiveness 64",
        writes="Dockings/vina_results_full_protein_vina_scoring_mgltools_exh64",
        canonical="posebusters_results/autodock_exhaustiveness_returns",
        determinism=BITEXACT,
    ),
    dict(
        name="14_ladder_exh64_gnina.yaml",
        source=DOCK / "autodock_vina_docking_config_full_protein_search_vina_mgltools_exh64_gnina.yaml",
        arm="AutoDock ladder rung, exhaustiveness 64 + gnina rescoring",
        thesis='App. B "Choice of Exhaustiveness and Rescoring"',
        supports="Table 8; Figure 18. Ties the selected arm at rank-1 (110 vs 111).",
        dataset=BENCH,
        driver="run_autodock_exhaustiveness_arm.py -c <this>",
        writes="Dockings/vina_results_full_protein_vina_scoring_mgltools_exh64/*/optimized_gnina",
        canonical="posebusters_results/autodock_exhaustiveness_returns",
        determinism=BITEXACT,
    ),
    dict(
        name="15_ladder_exh92.yaml",
        source=DOCK / "autodock_vina_docking_config_full_protein_search_vina_mgltools_exh92.yaml",
        arm="AutoDock ladder rung, exhaustiveness 92 (raw, not rescored)",
        thesis='App. B "Choice of Exhaustiveness and Rescoring"',
        supports="Table 8; Figures 18-21",
        dataset=BENCH,
        driver="run_autodock_exhaustiveness_arm.py -c <this> --expect-exhaustiveness 92",
        writes="Dockings/vina_results_full_protein_vina_scoring_mgltools_exh92",
        canonical="posebusters_results/autodock_exhaustiveness_returns",
        determinism=BITEXACT,
    ),
    # ------------------------------------------------------------- Orai x JKU
    dict(
        name="20_orai_jku_autodock_exh128_gnina.yaml",
        source=DOCK / "orai_jku_autodock_mgltools_exh128_gnina_config.yaml",
        arm="AutoDock Vina exh128 + gnina (Orai1 experimental panel)",
        thesis='App. B "Orai1 Panels"; 30 modes, 6 kcal/mol window',
        supports="Table 5; Figures 8-14",
        dataset=OJKU,
        driver="run_orai_mgltools_arm.py -c <this>   <-- the ONLY legal driver",
        writes="Dockings/Orai_JKU_MGLTools_exh128/mgl_tools",
        canonical="posebusters_results/_orai_matched_root/orai_jku",
        determinism=BITEXACT,
        deviations=[
            "none. run_autodock.py main() must NOT be used here: it re-prepares",
            "the Fr0 receptor, which is not bit-reproducible, and forces Meeko",
            "ligands. run_orai_mgltools_arm.py sha256-pins Fr0 against that.",
        ],
    ),
    dict(
        name="21_orai_jku_diffdock.yaml",
        source=DOCK / "diffdock_docking_config.yaml",
        arm="DiffDock-L (Orai1 experimental panel)",
        thesis='App. B "Orai1 Panels"; 30 samples',
        supports="Table 5; Figures 8-14",
        dataset=OJKU,
        driver="run_diffdock.py -c <this>",
        writes="Dockings/diffdock_results",
        canonical="posebusters_results/_orai_matched_root/orai_jku",
        determinism=BITEXACT + ".\n"
        "This is the ONE DiffDock arm that is reproducible: 24 stored run logs\n"
        "carry the '[repro] seeded' marker. App. B records the same asymmetry.",
        bake={"inference_seed": 42},
        deviations=[
            "inference_seed pinned to 42, which is what the stored run logs show.",
            "Requires the local seed patch in docking_tools/DiffDock/inference.py.",
        ],
    ),
    dict(
        name="22_orai_jku_equibind.yaml",
        source=DOCK / "equibind_orai_jku_config.yaml",
        arm="EquiBind unguided (Orai1 experimental panel)",
        thesis='App. B "Orai1 Panels"; 30 unguided conformers, UFF stage off',
        supports="Table 5; Figures 8-14",
        dataset=OJKU,
        driver="run_equibind.py, EQ_* environment mapped from this file",
        writes="Dockings/equibind_results_uffoff",
        canonical="posebusters_results/_orai_matched_root/orai_jku",
        determinism=BITEXACT,
    ),
    dict(
        name="23_orai_jku_posebusters.yaml",
        source=PB / "posebusters_orai_jku_config.yaml",
        arm="PoseBusters validity screen (Orai1 experimental panel)",
        thesis='App. B "PoseBusters Validity Checks"',
        supports="Table 5; Figures 8-10",
        dataset=OJKU,
        driver="Posebusters/run_posebusters.py -c <this>",
        writes="posebusters_results/orai_jku",
        canonical="posebusters_results/_orai_matched_root/orai_jku",
        determinism=DERIVED,
        deviations=[
            "none. variant_filter {autodock: gnina, diffdock: smina,",
            "equibind_guided: gnina} is load-bearing: raw ADFRsuite PDBQT has no",
            "REMARK SMILES, so an unfiltered run falls back to Open Babel and",
            "fails silently on roughly 63 per cent of poses.",
        ],
    ),
    dict(
        name="24_orai_jku_posebusters_fr0corrected.yaml",
        source=PB / "posebusters_orai_jku_mgltools_exh128_fr0corrected_config.yaml",
        arm="PoseBusters Fr0 correction, AUTODOCK ROWS ONLY (experimental panel)",
        thesis="App. B; AutoDock docked Fr0 into the OpenMM-minimised receptor\n"
        "while the main screen validates against the raw PDB",
        supports="the corrected Fr0 validity figures in Table 5",
        dataset=OJKU,
        driver="Posebusters/run_posebusters.py -c <this>",
        writes="posebusters_results/orai_jku_mgltools_exh128_fr0corrected",
        canonical="merged into _orai_matched_root/orai_jku for AutoDock rows only",
        determinism=DERIVED,
        deviations=[
            "none. DO NOT extend this pass to DiffDock or EquiBind. Those tools",
            "docked the RAW receptor, so the correction creates the mismatch",
            "instead of fixing it (EquiBind 58.2 -> 15.8 per cent). The",
            "*.WRONG-RECEPTOR-quarantined-1854 trees are that mistake on disk.",
        ],
    ),
    dict(
        name="25_orai_jku_pandamap.yaml",
        source=ANA / "pandamap_orai_jku_matched_config.yaml",
        arm="PandaMap fingerprints (Orai1 experimental panel)",
        thesis='App. B "Interaction Fingerprints"',
        supports="Table 12; Figures 13, 14, 25, 26",
        dataset=OJKU,
        driver="Analysis/run_pandamap.py -c <this>",
        writes="pandamap_results/orai_jku_matched",
        canonical="pandamap_results/orai_jku_matched",
        determinism=DERIVED,
    ),
    # -------------------------------------------------------- Orai x benchmark
    dict(
        name="30_orai_benchmark_autodock_exh128_gnina.yaml",
        source=DOCK / "orai_benchmark_autodock_mgltools_exh128_gnina_config.yaml",
        arm="AutoDock Vina exh128 + gnina (Orai1 control panel)",
        thesis='App. B "Orai1 Panels"; cut back to 10 modes at control parity',
        supports="Table 5; Figures 8-14",
        dataset=OBEN,
        driver="run_orai_mgltools_arm.py -c <this>   <-- the ONLY legal driver",
        writes="Dockings/Orai_Benchmark_MGLTools_exh128/mgl_tools",
        canonical="posebusters_results/_orai_matched_root/orai_benchmark",
        determinism=BITEXACT,
    ),
    dict(
        name="31_orai_benchmark_diffdock.yaml",
        source=DOCK / "diffdock_docking_config.yaml",
        arm="DiffDock-L (Orai1 control panel)",
        thesis='App. B "Orai1 Panels"; 10 samples',
        supports="Table 5; Figures 8-14",
        dataset=OBEN,
        driver="run_diffdock.py -c <this>",
        writes="Dockings/Orai_Benchmark_DiffDock",
        canonical="posebusters_results/_orai_matched_root/orai_benchmark",
        determinism=VERIFY + ".\n"
        "0 stored run logs carry the seed marker. App. B names this as the one\n"
        "remaining difference between the two Orai panels.",
        bake={
            "output_dir": "Dockings/Orai_Benchmark_DiffDock",
            "log_dir": "Dockings/Logs/orai_benchmark_diffdock_logs",
            "batch_mode": True,
            "overwrite_existing": False,
            "optimization": "gnina",
            "num_samples": 10,
            "inference_seed": None,
        },
        deviations=[
            "output_dir, log_dir, batch_mode, overwrite_existing and optimization",
            "baked in; the old notebook applied them in memory (cell 95).",
            "num_samples 30 -> 10, which is the control-panel budget App. B",
            "records and which the stored tree actually has.",
            "inference_seed null, recording that this run was unseeded.",
        ],
    ),
    dict(
        name="32_orai_benchmark_equibind.yaml",
        source=DOCK / "equibind_orai_benchmark_config.yaml",
        arm="EquiBind unguided (Orai1 control panel)",
        thesis='App. B "Orai1 Panels"; 10 unguided conformers',
        supports="Table 5; Figures 8-14",
        dataset=OBEN,
        driver="run_equibind.py, EQ_* environment mapped from this file",
        writes="Dockings/Orai_Benchmark_Equibind",
        canonical="posebusters_results/_orai_matched_root/orai_benchmark",
        determinism=BITEXACT,
        deviations=[
            "none. This is the one arm the old notebook already read strictly,",
            "via a cfg_req() helper that raises on a missing key. The new",
            "notebook applies that pattern to all three EquiBind arms.",
        ],
    ),
    dict(
        name="33_orai_benchmark_posebusters.yaml",
        source=PB / "posebusters_orai_benchmark_config.yaml",
        arm="PoseBusters validity screen (Orai1 control panel)",
        thesis='App. B "PoseBusters Validity Checks"',
        supports="Table 5; Figures 8-10",
        dataset=OBEN,
        driver="Posebusters/run_posebusters.py -c <this>",
        writes="posebusters_results/orai_benchmark",
        canonical="posebusters_results/_orai_matched_root/orai_benchmark",
        determinism=DERIVED,
    ),
    dict(
        name="34_orai_benchmark_posebusters_fr0corrected.yaml",
        source=PB / "posebusters_orai_benchmark_config_fr0corrected.yaml",
        arm="PoseBusters Fr0 correction, AUTODOCK ROWS ONLY (control panel)",
        thesis="App. B; the same receptor mismatch as config 24",
        supports="Fr0 validity 73.9 -> 97.7 per cent; the AutoDock control yield",
        dataset=OBEN,
        driver="Posebusters/run_posebusters.py -c <this>",
        writes="posebusters_results/orai_benchmark_fr0corrected",
        canonical="merged into _orai_matched_root/orai_benchmark for AutoDock rows only",
        determinism=DERIVED,
        deviations=[
            "none. AutoDock rows only, for the reason given in config 24.",
        ],
    ),
    dict(
        name="35_orai_benchmark_posebusters_gnina_fr0corrected.yaml",
        source=PB / "posebusters_orai_benchmark_gnina_config_fr0corrected.yaml",
        arm="PoseBusters Fr0 correction, gnina variant (control panel)",
        thesis="App. B; companion pass to config 34 for the gnina-rescored poses",
        supports="the AutoDock* rows of Table 5",
        dataset=OBEN,
        driver="Posebusters/run_posebusters.py -c <this>",
        writes="posebusters_results/orai_benchmark_gnina_fr0corrected",
        canonical="merged into _orai_matched_root/orai_benchmark",
        determinism=DERIVED,
    ),
    dict(
        name="36_orai_benchmark_pandamap.yaml",
        source=ANA / "pandamap_orai_benchmark_matched_config.yaml",
        arm="PandaMap fingerprints (Orai1 control panel)",
        thesis='App. B "Interaction Fingerprints"',
        supports="Table 12; Figures 13, 14, 25, 26",
        dataset=OBEN,
        driver="Analysis/run_pandamap.py -c <this>",
        writes="pandamap_results/orai_benchmark_matched",
        canonical="pandamap_results/orai_benchmark_matched",
        determinism=DERIVED,
    ),
]


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    problems = []
    written = 0

    for spec in SPECS:
        src: Path = spec["source"]
        if not src.exists():
            problems.append(f"MISSING SOURCE {src}")
            continue
        body = src.read_text()
        target = OUT / spec["name"]

        bake = spec.get("bake")
        if bake:
            cfg = yaml.safe_load(body)
            cfg.update(bake)
            body = yaml.safe_dump(cfg, default_flow_style=False, sort_keys=False)

        rel = src.relative_to(ROOT)
        target.write_text(header(source=str(rel), **{k: v for k, v in spec.items()
                                                     if k not in {"name", "source", "bake"}}) + body)
        written += 1

        # verify parsed equality
        got = yaml.safe_load(target.read_text())
        want = yaml.safe_load(src.read_text())
        if bake:
            want.update(bake)
        if got != want:
            diff = {k for k in set(got) | set(want) if got.get(k) != want.get(k)}
            problems.append(f"VALUE DRIFT {spec['name']}: {sorted(diff)}")

    print(f"wrote {written} configs to {OUT}")
    if problems:
        print("\nPROBLEMS:")
        for p in problems:
            print("  -", p)
        return 1
    print("all configs parse identically to their sources (with declared bakes)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
