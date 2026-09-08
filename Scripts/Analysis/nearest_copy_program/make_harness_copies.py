#!/usr/bin/env python
"""Create the program's harness copies from the canonical ones (plan v2 Phase 7, independent layout).

Writes, next to this script:
  thesis_expected_values_nearest.yaml  -- canonical yaml with every canonical result path and the thesis
                                          directory remapped to the program's _nearest counterparts
  thesis_assertions_nearest.py         -- canonical harness with SPEC pointing at that yaml and the two
                                          hard-wired thesis_latex/media references remapped

The canonical files are read only. Existing copies are NOT overwritten unless --force is given, because
Phase 7 re-transcribes values into the yaml copy by hand and a regeneration would discard that work.
Values are NOT changed here; only paths. Every replacement is counted and printed.
"""
from __future__ import annotations
import argparse, re, sys
from pathlib import Path

REPO = Path("/home/manndo/master_dev")
HERE = Path(__file__).resolve().parent
CANON_YAML = REPO / "Scripts/Analysis/thesis_expected_values.yaml"
CANON_ASSERT = REPO / "Scripts/Analysis/thesis_assertions.py"

# Longest keys first so a prefix never shadows a longer path.
PATH_MAP = [
    ("posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report/",
     "posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report_nearest/"),
    ("posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report\"",
     "posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report_nearest\""),
    ("posebusters_results/benchmark_matched_equibind/dock/validity_report_mgltools/",
     "posebusters_results/benchmark_matched_equibind/dock/validity_report_mgltools_nearest/"),
    ("posebusters_results/cluster_crystal_pocket_matched_equibind/autodock_mgltools_exh128_gnina__diffdock_smina_allposes/",
     "posebusters_results/cluster_crystal_pocket_matched_equibind/autodock_mgltools_exh128_gnina__diffdock_smina_allposes_nearest/"),
    ("posebusters_results/benchmark_matched_equibind/docking_effort_gnina_v2_charged/",
     "posebusters_results/benchmark_matched_equibind/docking_effort_gnina_v2_charged_nearest/"),
    ("posebusters_results/benchmark_matched_equibind/docking_effort_gnina_v2/",
     "posebusters_results/benchmark_matched_equibind/docking_effort_gnina_v2_nearest/"),
    ("posebusters_results/autodock_exhaustiveness_returns/",
     "posebusters_results/autodock_exhaustiveness_returns_nearest/"),
    ("pandamap_results/benchmark_matched_equibind/report/",
     "pandamap_results/benchmark_matched_equibind_nearest/report/"),
    ("pandamap_results/benchmark_matched_equibind/",
     "pandamap_results/benchmark_matched_equibind_nearest/"),
    ("PoseBusters_Benchmark_Analysis/ligand_difficulty/", "PoseBusters_Benchmark_Analysis/ligand_difficulty_nearest/"),
    ("PoseBusters_Benchmark_Analysis/receptor_difficulty/", "PoseBusters_Benchmark_Analysis/receptor_difficulty_nearest/"),
    ("PoseBusters_Benchmark_Analysis/smina_rerank/", "PoseBusters_Benchmark_Analysis/smina_rerank_nearest/"),
    ("thesis_latex/media/media", "thesis_latex_nearest/media/media"),
    ("thesis_latex/body_main_short.tex", "thesis_latex_nearest/body_main_short.tex"),
    ("thesis_latex/body_appendix_short.tex", "thesis_latex_nearest/body_appendix_short.tex"),
    ("thesis_latex/Thesis_short.tex", "thesis_latex_nearest/Thesis_short.tex"),
]

def remap(text: str, label: str) -> str:
    counts = []
    for old, new in PATH_MAP:
        n = text.count(old)
        if n:
            text = text.replace(old, new); counts.append((old, n))
    print(f"{label}: {sum(n for _, n in counts)} replacements")
    for old, n in counts:
        print(f"   {n:3d}  {old}")
    return text

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--force", action="store_true", help="overwrite existing copies (discards hand edits)")
    a = ap.parse_args()
    out_yaml = HERE / "thesis_expected_values_nearest.yaml"
    out_py = HERE / "thesis_assertions_nearest.py"
    for p in (out_yaml, out_py):
        if p.exists() and not a.force:
            sys.exit(f"REFUSING: {p.name} exists (hand-edited in Phase 7?). Pass --force to regenerate.")

    y = CANON_YAML.read_text()
    y = remap(y, "yaml")
    header = ("# ---------------------------------------------------------------------------\n"
              "# PROGRAM COPY (nearest-copy endpoint, independent layout). Generated from\n"
              "# Scripts/Analysis/thesis_expected_values.yaml by make_harness_copies.py with the\n"
              "# canonical result paths and the thesis directory remapped to their _nearest\n"
              "# counterparts. VALUES ARE STILL THE CANONICAL ONES until Phase 7 re-transcribes\n"
              "# them from thesis_latex_nearest/. The canonical yaml is untouched.\n"
              "# ---------------------------------------------------------------------------\n")
    y = y.replace("meta:\n", "meta:\n  reference_convention: nearest-copy   # program copy; see header\n", 1)
    out_yaml.write_text(header + y)

    s = CANON_ASSERT.read_text()
    old_spec = 'SPEC = Path(__file__).resolve().parent / "thesis_expected_values.yaml"'
    assert s.count(old_spec) == 1, "SPEC line not found in thesis_assertions.py"
    s = s.replace(old_spec, 'SPEC = Path(__file__).resolve().parent / "thesis_expected_values_nearest.yaml"  # program copy')
    old_skip = 'skip=("thesis_latex/media", "obsolete", ".backup")'
    assert s.count(old_skip) == 1, "skip tuple not found"
    s = s.replace(old_skip, 'skip=("thesis_latex/media", "thesis_latex_nearest/media", "obsolete", ".backup")')
    old_media = '(ROOT / "thesis_latex/media/media")'
    assert s.count(old_media) == 1, "media glob not found"
    s = s.replace(old_media, '(ROOT / "thesis_latex_nearest/media/media")')
    s = ('# PROGRAM COPY of Scripts/Analysis/thesis_assertions.py (nearest-copy endpoint). Only SPEC and the\n'
         '# two hard-wired thesis_latex/media references differ; the canonical harness is untouched.\n' + s)
    out_py.write_text(s)
    print(f"wrote {out_yaml.name} and {out_py.name} under {HERE}")

if __name__ == "__main__":
    main()
