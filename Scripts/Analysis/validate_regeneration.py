#!/usr/bin/env python
"""Re-run each figure's generator and compare the bytes to the published asset.

WHY THIS IS SEPARATE from `thesis_assertions.py`. That module proves PROVENANCE:
the image in the thesis is byte-identical to a file sitting on disk. This module
proves REGENERATION: running the generator again, from the canonical inputs,
produces that same file. They are different claims, and only the second one shows
the pipeline still works.

Each generator writes into a scratch directory, never into the canonical tree, so
a failed or partial run cannot damage a published output.

Figures whose generator is not deterministic, or whose source is not a script at
all, are declared rather than silently skipped. The dataset chapter's figures are
the honest exception: their statistics are stable but matplotlib moves roughly one
to two per cent of pixels between renders, so a byte comparison there would report
a failure that is not one.

    python Scripts/Analysis/validate_regeneration.py            # cheap + moderate
    python Scripts/Analysis/validate_regeneration.py --quick    # cheap only
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass, field
from pathlib import Path

ROOT = Path("/home/manndo/master_dev")
PY = "/home/manndo/anaconda3/envs/vina/bin/python"
ANA = ROOT / "Scripts/Analysis"
CFG = ROOT / "Scripts/Docking/configs_thesis"
MEDIA = ROOT / "thesis_latex/media/media"
RESULT = ROOT / "posebusters_results/_reproduction/regeneration_check.json"

BENCH = ROOT / "posebusters_results/benchmark_matched_equibind"
REPORT = BENCH / "dock/pose_comparison_report"
ORAI = ROOT / "posebusters_results/_orai_matched_root"
PMAP_B = ROOT / "pandamap_results/benchmark_matched_equibind"
EXH = ROOT / "posebusters_results/autodock_exhaustiveness_returns"

LADDER_EXCL = ("autodock_mgltools,autodock_mgltools_gnina,autodock_mgltools_exh18,"
               "autodock_mgltools_exh64,autodock_mgltools_exh64_gnina,autodock_mgltools_exh92")
VARIANTS = ["--autodock-variant", "gnina", "--diffdock-variant", "smina",
            "--equibind-variant", "gnina"]


@dataclass
class Case:
    """One generator invocation and the assets it should reproduce."""
    figures: str                      # thesis figure numbers, for the report
    name: str
    cost: str                         # cheap | moderate
    argv: list                        # {out} is replaced with the scratch dir
    # published asset -> filename the generator writes inside the scratch dir
    assets: dict[str, str]
    stage: list[tuple[Path, str]] = field(default_factory=list)  # files copied in first
    note: str = ""


CASES = [
    Case("1", "figure2_validity_yield", "cheap",
         [PY, ANA / "figure2_validity_yield.py",
          "--csv", BENCH / "dock/posebusters_filtered_results.csv",
          "--out", "{out}/00_figure2_validity_yield.png",
          "--exclude-preset", "meeko", "--exclude-methods", LADDER_EXCL],
         {"image2": "00_figure2_validity_yield.png"},
         note="the double exclusion is load-bearing; the preset alone pools five ladder rungs"),

    Case("8", "orai_pbvalid_yield_compare", "cheap",
         [PY, ANA / "orai_pbvalid_yield_compare.py",
          "--results-root", ORAI, "--out-dir", "{out}"],
         {"image13": "09f_pbvalid_yield_dominant_orai_benchmark_vs_experimental.png"},
         note="extracted from an inline notebook cell on 2026-09-02"),

    Case("9, 10", "orai_pbvalid_tm_share_compare", "cheap",
         [PY, ANA / "orai_pbvalid_tm_share_compare.py", *VARIANTS,
          "--results-root", ORAI, "--out-dir", "{out}", "--top-n-poses", "10"],
         {"image43": "fig_pbvalid_outside_tm_whisker.png",
          "image15": "fig_tm_loss_relative_compare.png"}),

    Case("11, 12, 23", "orai_cross_tool_agreement_compare", "cheap",
         [PY, ANA / "orai_cross_tool_agreement_compare.py",
          "--exp-csv", ORAI / "orai_jku/pose_clusters/per_pair.csv",
          "--exp-label", "Exp. Ligands",
          "--bench-csv", ORAI / "orai_benchmark/pose_clusters/per_pair.csv",
          "--bench-label", "Benchmark ligands × Orai",
          "--exp-quality-csv", ORAI / "orai_jku/pose_clusters/cluster_quality_per_tool.csv",
          "--bench-quality-csv", ORAI / "orai_benchmark/pose_clusters/cluster_quality_per_tool.csv",
          "--out-dir", "{out}"],
         {"image16": "orai_cross_tool_agreement_compare.png",
          "image17": "orai_tool_agreement_compare.png",
          "image18": "orai_cluster_quality_compare.png"},
         note="takes no --autodock-variant; its arm is fixed by which run wrote per_pair.csv"),

    Case("13, 14, 25, 26", "orai_pandamap_interaction_compare", "cheap",
         [PY, ANA / "orai_pandamap_interaction_compare.py",
          "--exp-dir", ROOT / "pandamap_results/orai_jku_matched",
          "--bench-dir", ROOT / "pandamap_results/orai_benchmark_matched",
          "--out-dir", "{out}", *VARIANTS, "--top-n-poses", "10",
          "--exp-posebusters-csv", ORAI / "orai_jku/dock/posebusters_filtered_results.csv",
          "--bench-posebusters-csv", ORAI / "orai_benchmark/dock/posebusters_filtered_results.csv"],
         {"image21": "fig_type_profile_compare.png",
          "image23": "fig_fingerprint_overlap_compare.png",
          "image20": "fig_total_interactions_compare.png",
          "image22": "fig_residue_hotspots_compare.png"},
         note="the two --*-posebusters-csv paths are REQUIRED. Their defaults name the "
              "plain trees, and EquiBind's PandaMap pose_rank is a 999 placeholder, so "
              "without them --top-n-poses cannot rank EquiBind and drops all 324 of its "
              "pose rows in silence. Caught by this harness on 2026-09-02."),

    Case("4", "filmstrip_rank1_top5_top15", "cheap",
         [PY, ANA / "filmstrip_rank1_top5_top15.py",
          "--report-dir", "{out}", "--exclude-preset", "meeko"],
         {"image5": "20d_form_vs_placement_by_family__depth_filmstrip_pbvalid"
                    "__rank1_top5_top15__thesis.png"},
         stage=[(REPORT / "per_pose_metrics.csv", "per_pose_metrics.csv"),
                (REPORT / "per_pose_metrics.manifest.json", "per_pose_metrics.manifest.json")],
         note="--report-dir was added 2026-09-02; the module constant still names the "
              "pre-matched tree"),

    Case("15", "docking_effort_comparison (charged)", "moderate",
         [PY, ANA / "docking_effort_comparison.py", "--dataset", "benchmark",
          "--autodock-dir", "Dockings/vina_results_full_protein_vina_scoring_mgltools_exh128",
          "--autodock-prep", "mgl_tools", "--autodock-refine", "gnina", "--autodock-gnina-gpu",
          "--autodock-method", "autodock_mgltools_exh128_gnina",
          "--autodock-optimizer-workers", "16",
          "--equibind-dir", "Dockings/Benchmark_Equibind_cputimed",
          "--unidock2-dir", "", "--unidock-dir", "",
          "--per-pose-csv", REPORT / "per_pose_metrics.csv",
          "--out-dir", "{out}", "--basis", "charged", "--cpu-threads", "32"],
         {"image24": "panels/effort_by_quality_near2.png",
          "image25": "panels/resource_per_near_native_valid_pose.png"},
         note="the charged basis is invariant to the optimiser worker count; elapsed is not"),

    Case("18, 19, 20, 21", "autodock_exhaustiveness_returns", "moderate",
         [PY, ANA / "autodock_exhaustiveness_returns.py",
          "--per-pose-csv", REPORT / "per_pose_metrics.csv",
          "--out-dir", "{out}", "--n-boot", "5000", "--verify-input-parity"],
         {"image44": "figures/exh_raw_vs_rescored.png",
          "image45": "figures/exh_depth_sweep.png",
          "image46": "figures/exh_yield_vs_cost.png",
          "image47": "figures/exh_marginal_return.png"},
         note="5,000 paired complex-level bootstrap resamples at a fixed seed"),

    Case("7, 39", "pandamap_interaction_report", "moderate",
         [PY, ANA / "pandamap_interaction_report.py",
          "-c", CFG / "07_benchmark_pandamap.yaml", "--in-dir", PMAP_B,
          "--out-dir", "{out}",
          "--per-pose-metrics", REPORT / "per_pose_metrics.csv",
          "--exclude-methods", "unidock2,autodock_gnina"],
         {"image9": "05_fingerprint_similarity_top5.png",
          "image11": "10b_contact_decomposition_by_rank.png"}),
]

# Declared rather than tested, with the reason.
DECLARED = [
    ("2, 3", "posebusters_pose_comparison re-render",
     "writes ~92 PNGs into its own report directory and reads a cache keyed to that "
     "directory; a scratch run is not a faithful reproduction of the published call"),
    ("5, 6, 38", "pose_cluster_crystal_pocket_report",
     "same, and its analysis cache is keyed to the out-dir signature"),
    ("24", "orai_pose_cluster_report",
     "writes into the pose_clusters directory the Orai analyses read; a scratch run "
     "would not exercise the published path"),
    ("17", "Flow Charts.ipynb", "author-drawn schematic, not a data product"),
    ("22, 27", "hand-composed PyMOL renders",
     "present nowhere in the repository; the scene builder exists, the render was interactive"),
    ("28-37", "PoseBusters_DataSet_Analysis.ipynb",
     "statistics are stable but matplotlib moves 1 to 2.5 per cent of pixels between "
     "renders, so a byte comparison would report a failure that is not one"),
]


def md5(p: Path) -> str | None:
    try:
        return hashlib.md5(p.read_bytes()).hexdigest()
    except Exception:
        return None


# =============================================================================
# Result file
# =============================================================================
# A run records what it saw, so a later reader is told the result AND whether it
# still applies. Without that, a cached "20 of 20" outlives the code it describes
# and becomes a claim nobody has checked. Two fingerprints decide staleness: the
# mtime of each generator script, which catches an edit, and the md5 of each
# published asset, which catches a figure being replaced.


def _fingerprint(cases) -> dict:
    scripts, assets = {}, {}
    for c in cases:
        for a in c.argv:
            s = str(a)
            if s.endswith(".py") and Path(s).exists():
                scripts[str(Path(s).relative_to(ROOT))] = round(Path(s).stat().st_mtime, 3)
        for asset in c.assets:
            f = MEDIA / f"{asset}.png"
            if f.exists():
                assets[asset] = md5(f)
    return {"scripts": scripts, "assets": assets}


def write_result(cases, rows, scope: str) -> Path:
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "written": time.strftime("%Y-%m-%d %H:%M:%S"),
        "scope": scope,
        "n_generators": len(cases),
        "n_assets": len(rows),
        "n_identical": sum(1 for r in rows if r[1] == "IDENTICAL"),
        "cases": [c.name for c in cases],
        "rows": [{"asset": a, "verdict": v, "detail": d} for a, v, d in rows],
        "fingerprint": _fingerprint(cases),
    }
    RESULT.write_text(json.dumps(payload, indent=2))
    return RESULT


def load_result() -> dict | None:
    try:
        return json.loads(RESULT.read_text())
    except Exception:
        return None


def describe_result(cases=None) -> str:
    """One block of text for the notebook: the last result and whether it holds.

    Staleness is judged only against the cases the recorded run actually covered.
    A --quick run legitimately leaves the moderate generators unchecked, and
    reporting that as staleness would cry wolf. Anything it did not cover is
    reported separately, as unchecked rather than as stale.
    """
    cases = cases if cases is not None else CASES
    n_assets = sum(len(c.assets) for c in cases)
    head = (f"{len(cases)} generators covering {n_assets} published assets.\n"
            f"Run with:  python Scripts/Analysis/validate_regeneration.py")
    res = load_result()
    if not res:
        return head + "\n\nNo recorded run. The check has never been run on this tree."

    ran_names = set(res.get("cases") or [])
    ran = [c for c in cases if c.name in ran_names] or cases
    not_ran = [c for c in cases if c.name not in ran_names] if ran_names else []

    now = _fingerprint(ran)
    old = res.get("fingerprint", {})
    changed = []
    for kind in ("scripts", "assets"):
        for k, v in old.get(kind, {}).items():
            if k in now[kind] and now[kind][k] != v:
                changed.append(f"{kind[:-1]} {k}")

    verdict = (f"{res['n_identical']} of {res['n_assets']} assets byte-identical"
               f" ({res['scope']} run, {res['written']})")
    lines = [head, "", f"Last run: {verdict}."]
    if changed:
        lines.append(f"STALE: {len(changed)} input(s) changed since, so that result "
                     f"no longer applies:")
        lines += [f"  {c}" for c in changed[:8]]
        if len(changed) > 8:
            lines.append(f"  and {len(changed) - 8} more")
        lines.append("Re-run the check.")
    else:
        lines.append("Nothing it covered has changed since, so it still holds.")
    if not_ran:
        n = sum(len(c.assets) for c in not_ran)
        lines.append(f"NOT COVERED by that run: {len(not_ran)} generator(s), {n} asset(s) "
                     f"— {', '.join(c.name for c in not_ran)}")
    return "\n".join(lines)


def run_case(c: Case, keep: Path | None = None) -> list[tuple[str, str, str]]:
    """Returns (asset, verdict, detail) per published asset."""
    out = Path(tempfile.mkdtemp(prefix=f"regen_{c.name.split()[0]}_"))
    try:
        for src, dst in c.stage:
            if src.exists():
                shutil.copy2(src, out / dst)
        argv = [str(a).replace("{out}", str(out)) for a in c.argv]
        t0 = time.time()
        proc = subprocess.run(argv, cwd=ROOT, capture_output=True, text=True)
        elapsed = time.time() - t0
        if proc.returncode != 0:
            tail = (proc.stderr or proc.stdout or "").strip().splitlines()[-1:] or [""]
            return [(a, "RUN FAILED", f"exit {proc.returncode}: {tail[0][:90]}")
                    for a in c.assets]
        rows = []
        for asset, fname in c.assets.items():
            got, want = md5(out / fname), md5(MEDIA / f"{asset}.png")
            if got is None:
                rows.append((asset, "NOT WRITTEN", fname))
            elif got == want:
                rows.append((asset, "IDENTICAL", f"{elapsed:.0f}s"))
            else:
                rows.append((asset, "DIFFERS", f"regen {got[:8]} vs published {want[:8]}"))
        if keep:
            shutil.copytree(out, keep / c.name.replace(" ", "_"), dirs_exist_ok=True)
        return rows
    finally:
        shutil.rmtree(out, ignore_errors=True)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--quick", action="store_true", help="cheap cases only")
    ap.add_argument("--only", default=None, help="substring match on the case name")
    ap.add_argument("--keep", default=None, help="copy scratch output here for inspection")
    args = ap.parse_args(argv)

    cases = [c for c in CASES
             if (not args.quick or c.cost == "cheap")
             and (not args.only or args.only in c.name)]
    keep = Path(args.keep) if args.keep else None

    print(f"Re-running {len(cases)} generator(s) into scratch directories.\n")
    ok = bad = 0
    rows_all = []
    for c in cases:
        print(f"  Figure {c.figures} — {c.name}")
        for asset, verdict, detail in run_case(c, keep):
            rows_all.append((asset, verdict, detail))
            flag = "ok " if verdict == "IDENTICAL" else "** "
            ok += verdict == "IDENTICAL"
            bad += verdict != "IDENTICAL"
            print(f"    {flag}{asset:9s} {verdict:12s} {detail}")
        if c.note:
            print(f"       note: {c.note}")
        print()

    scope = "quick" if args.quick else ("partial" if args.only else "full")
    if not args.only:
        path = write_result(cases, rows_all, scope)
        print(f"recorded -> {path.relative_to(ROOT)}")
    print(f"{ok} assets byte-identical to the published figure, {bad} not.\n")
    if DECLARED:
        print("Declared rather than tested:")
        for figs, what, why in DECLARED:
            print(f"  Figure {figs:9s} {what}")
            print(f"                  {why}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
