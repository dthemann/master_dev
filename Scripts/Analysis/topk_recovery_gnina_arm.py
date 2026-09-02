#!/usr/bin/env python3
"""Emit the gnina/smina-arm rank-depth recovery sidecar (thesis Tables 5 / 6 / 7).

WHY THIS SCRIPT EXISTS
----------------------
The report generator collapses every engine to ONE arm before it writes the
presentation sidecars, and WHICH arm that is depends on the run flags
(``--collapse-plots-only``, ``--collapse-diffdock-variant``, ``--best-equibind-only``)
and on the pin inside ``_select_autodock_arm``. The sidecars behind the thesis
Tables 5-7 — ``18_topn_within_thresholds_pbvalid_depths_report.txt`` and
``topn_within_thresholds_pbvalid_depths.csv`` — therefore report the surviving arm
under the anonymous ``autodock`` / ``diffdock`` keys, with nothing in the file
recording which variant actually produced the row. Read back later, those files cannot
be attributed to an arm at all.

This script pins the three arms BY NAME instead, so the Tables 5-7 family is
reproducible from a flag-free invocation and the emitted ``method_key`` column states
the arm on every row. It touches neither the report generator nor the existing
sidecars. It imports ``posebusters_pose_comparison`` and calls the very same helpers
the thesis numbers came from:

    * ``_topn_within_frames``   — per-complex best-of-top-d frames (identical de-dup of
      DiffDock's duplicated rank-1 pose, identical gnina-affinity ranking for EquiBind,
      identical PoseBusters gate);
    * ``_stats_within_by_depth`` — Cochran's Q + pairwise exact McNemar (Holm) between
      tools at each depth, and the within-tool Cochran Q + Wilson-CI depth headroom;
    * ``aggregate_within_thresholds_by_depth`` / ``_write_within_by_depth_report`` — the
      full threshold sweep and the report text, byte-format-identical to the raw-arm
      report so the two arms can be diffed line for line.

The only thing this script does differently is WHICH per-method rows are handed to
those helpers: the AutoDock slot is filled with the dominant arm
``pc._DOMINANT_AUTODOCK_OPT`` (whose ``rank`` column already IS ``optimized_rank``,
i.e. the gnina re-ranking) and the DiffDock slot with ``diffdock_smina``, each
relabelled to the canonical ``autodock`` / ``diffdock`` key the helpers key off.
EquiBind is passed as ``eq_df`` exactly as ``main()`` does.

The committed sidecar ``topk_recovery_validity_gnina_arm.csv`` (n = 303 complexes)
records what those defaults produce: PB-valid near-native complexes at ranking depths
1 / 15 / 30 are 111 / 198 / 199 for the AutoDock arm, 103 / 167 / 169 for
``diffdock_smina`` and 46 / 61 / 62 for ``equibind_unguided_gnina``.

THE GATE
--------
A complex counts at ranking depth d when ANY pose among the tool's first d ranked poses
is BOTH PoseBusters-valid AND within 2 A heavy-atom RMSD of the crystal ligand (an
EXISTENCE gate, not "the rank-1 pose"). Depth sets are nested (top-1 c top-15 c top-30),
so depth steps carry an effect size (recovered share + Wilson 95 % CI) rather than a
degenerate McNemar p-value; Cochran's Q across the three depths is the omnibus.

USAGE
-----
    conda activate vina
    python Scripts/Analysis/topk_recovery_gnina_arm.py          # defaults = the thesis run

Writes (into the report dir, new files only — nothing existing is overwritten):
    topk_recovery_validity_gnina_arm.csv
    topk_recovery_validity_gnina_arm.txt
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

import posebusters_pose_comparison as pc          # noqa: E402
import method_filter as mf  # noqa: E402  (shared single-point method exclusion)

DEFAULT_REPORT_DIR = Path(
    "/home/manndo/master_dev/posebusters_results/benchmark_full_protein_vina_scoring"
    "/dock/pose_comparison_report")

DEFAULT_DEPTHS = (1, 15, 30)
DEFAULT_THRESHOLD = 2.0
DEFAULT_STEM = "topk_recovery_validity_gnina_arm"

# The three printed arms. (slot key the helpers use, per_pose_metrics method key, label)
# The AutoDock entry is read from posebusters_pose_comparison — which must already be
# imported above, or the dict body raises NameError at import time — so this file cannot
# drift from the arm the rest of the suite reports. The literal is only a fallback for an
# older sibling module that predates the constant. It was previously the Meeko-ligand key
# 'autodock_gnina', which silently rebuilt the Tables 5/6/7 family on the wrong arm for
# anyone who ran this script without --autodock-variant.
DEFAULT_ARM = {
    "autodock": getattr(pc, "_DOMINANT_AUTODOCK_OPT", "autodock_mgltools_exh128_gnina"),
    "diffdock": "diffdock_smina",
    "equibind": "equibind_unguided_gnina",
}
# Kept to 15 characters — the reused report writer slices column headers to 15.
ARM_LABEL = {
    "autodock": "AutoDock gnina*",
    "diffdock": "DiffDock smina*",
    "equibind": "EquiBind gnina*",
}


def _build_arm_frames(df_full: pd.DataFrame, arm: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return ``(df, eq_df)`` in the shape ``main()`` hands to the fig-18 helpers.

    ``df`` holds the two ranking tools relabelled to their canonical keys; ``eq_df``
    holds the EquiBind variant (no native ranking — the helpers rank it by gnina
    affinity). Raises if a requested variant is absent from per_pose_metrics.csv.
    """
    present = set(df_full["method"].astype(str))
    for key, mkey in arm.items():
        if mkey not in present:
            raise SystemExit(f"variant {mkey!r} (slot {key}) not in per_pose_metrics.csv; "
                             f"present: {sorted(present)}")
    parts = []
    for slot in ("autodock", "diffdock"):
        sub = df_full[df_full["method"].astype(str) == arm[slot]].copy()
        sub["method"] = slot            # canonical key the fig-18 helpers group on
        parts.append(sub)
    df = pd.concat(parts, ignore_index=True)
    eq_df = df_full[df_full["method"].astype(str) == arm["equibind"]].copy()
    return df, eq_df


def _counts(df: pd.DataFrame, eq_df: pd.DataFrame, depths, thr: float) -> dict:
    """{(slot, depth): (n, near_native_k, pbvalid_near_native_k)} via the module helper."""
    out: dict = {}
    for d in depths:
        gated = pc._topn_within_frames(df, int(d), eq_df=eq_df, pb_valid_only=True,
                                       rmsd_col="rmsd")
        plain = pc._topn_within_frames(df, int(d), eq_df=eq_df, pb_valid_only=False,
                                       rmsd_col="rmsd")
        for slot, fr in gated.items():
            vk = int((pd.to_numeric(fr["best_rmsd"], errors="coerce") <= thr).sum())
            pf = plain.get(slot)
            nk = (int((pd.to_numeric(pf["best_rmsd"], errors="coerce") <= thr).sum())
                  if pf is not None else float("nan"))
            out[(slot, int(d))] = (len(fr), nk, vk)
    return out


def _fmt_p(p) -> str:
    if p is None or p != p:
        return "n/a"
    return "<1e-300" if p <= 0 else (f"{p:.3e}" if p < 1e-3 else f"{p:.4f}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--report-dir", type=Path, default=DEFAULT_REPORT_DIR,
                    help="pose_comparison_report dir holding per_pose_metrics.csv")
    ap.add_argument("--depths", type=int, nargs="+", default=list(DEFAULT_DEPTHS))
    ap.add_argument("--threshold", type=float, default=DEFAULT_THRESHOLD,
                    help="near-native RMSD threshold in Angstrom (default 2.0)")
    # The defaults are printed: which arm filled each slot is the one thing a reader of
    # Tables 5/6/7 has to be able to check without opening the source.
    ap.add_argument("--autodock-variant", default=DEFAULT_ARM["autodock"],
                    help="per_pose_metrics method key for the AutoDock slot "
                         "(default: %(default)s)")
    ap.add_argument("--diffdock-variant", default=DEFAULT_ARM["diffdock"],
                    help="per_pose_metrics method key for the DiffDock slot "
                         "(default: %(default)s)")
    ap.add_argument("--equibind-variant", default=DEFAULT_ARM["equibind"],
                    help="per_pose_metrics method key for the EquiBind slot "
                         "(default: %(default)s)")
    ap.add_argument("--out-stem", default=DEFAULT_STEM)
    mf.add_method_filter_args(ap)
    args = ap.parse_args()

    depths = sorted({int(d) for d in args.depths})
    thr = float(args.threshold)
    arm = {"autodock": args.autodock_variant,
           "diffdock": args.diffdock_variant,
           "equibind": args.equibind_variant}

    per_pose = args.report_dir / "per_pose_metrics.csv"
    if not per_pose.is_file():
        raise SystemExit(f"missing {per_pose}")
    df_full = pd.read_csv(per_pose, low_memory=False)
    # Before _build_arm_frames, which overwrites each arm's method key with its
    # slot name ('autodock_gnina' -> 'autodock') and validates the requested
    # variants against what is present.
    df_full = mf.apply_method_filter(df_full, "method", args, label="topk-recovery",
                                     out_dir=args.report_dir)
    df, eq_df = _build_arm_frames(df_full, arm)

    # Labels so the reused report writer names the ACTUAL variants, not the slots.
    for slot, label in ARM_LABEL.items():
        pc._LABEL_OVERRIDES[slot] = f"{label}"

    counts = _counts(df, eq_df, depths, thr)
    stats = pc._stats_within_by_depth(
        df, depths, thr, eq_df=eq_df, pb_valid_only=True, rmsd_col="rmsd",
        figure_name=f"{args.out_stem} (no figure — sidecar only)", metric_label="RMSD")
    if not stats:
        raise SystemExit("stats failed — fewer than two tools carry per-complex frames")

    headroom = {h["method"]: h for h in stats.get("depth_headroom", [])}
    between = {int(b["depth"]): b for b in stats.get("between_tools_by_depth", [])}
    slots = [s for s in ("autodock", "diffdock", "equibind") if s in headroom]

    # ── CSV (tidy long form; row_type keeps the three blocks machine-separable) ──
    rows = []
    for slot in slots:
        h = headroom[slot]
        q = h["omnibus_across_depths"]
        steps = {(s["from_depth"], s["to_depth"]): s for s in h["steps"]}
        prev = None
        for d in depths:
            n, near_k, valid_k = counts[(slot, d)]
            step = steps.get((prev, d)) if prev is not None else None
            lo, hi = (step["recovered_ci"] if step else (float("nan"), float("nan")))
            rows.append({
                "row_type": "recovery",
                "arm": "gnina/smina (post-hoc optimised)",
                "slot": slot,
                "method_key": arm[slot],
                "label": ARM_LABEL[slot],
                "depth": d,
                "n_complexes": n,
                "near_native_k": near_k,
                "near_native_pct": round(100.0 * near_k / n, 4),
                "pbvalid_near_native_k": valid_k,
                "pbvalid_near_native_pct": round(100.0 * valid_k / n, 4),
                "gain_from_prev_depth_k": (step["recovered_k"] if step else ""),
                "gain_from_prev_depth_pp": ("" if prev is None else
                                            round(100.0 * (valid_k - counts[(slot, prev)][2]) / n, 4)),
                "gain_recovered_share_pct": ("" if step is None else
                                             round(100.0 * step["recovered_share"], 4)),
                "gain_wilson_lo_pct": ("" if step is None else round(100.0 * lo, 4)),
                "gain_wilson_hi_pct": ("" if step is None else round(100.0 * hi, 4)),
                "within_tool_cochran_q": round(float(q["Q"]), 4),
                "within_tool_cochran_df": int(q["df"]),
                "within_tool_cochran_p": q["p"],
            })
            prev = d
        # end-to-end step (top-1 -> deepest) as its own row
        s = steps.get((depths[0], depths[-1]))
        if s:
            lo, hi = s["recovered_ci"]
            rows.append({
                "row_type": "depth_step",
                "arm": "gnina/smina (post-hoc optimised)",
                "slot": slot, "method_key": arm[slot], "label": ARM_LABEL[slot],
                "depth": f"{depths[0]}->{depths[-1]}",
                "n_complexes": s and h["n"],
                "gain_from_prev_depth_k": s["recovered_k"],
                "gain_recovered_share_pct": round(100.0 * s["recovered_share"], 4),
                "gain_wilson_lo_pct": round(100.0 * lo, 4),
                "gain_wilson_hi_pct": round(100.0 * hi, 4),
                "within_tool_cochran_q": round(float(h["omnibus_across_depths"]["Q"]), 4),
                "within_tool_cochran_df": int(h["omnibus_across_depths"]["df"]),
                "within_tool_cochran_p": h["omnibus_across_depths"]["p"],
            })
    for d in depths:
        b = between.get(d)
        if not b or "omnibus" not in b:
            continue
        om = b["omnibus"]
        for pr in b.get("pairwise", []):
            rows.append({
                "row_type": "between_tools_mcnemar",
                "arm": "gnina/smina (post-hoc optimised)",
                "depth": d,
                "n_complexes": b["n"],
                "pair_a": ARM_LABEL.get(pr["a"], pr["a"]),
                "pair_b": ARM_LABEL.get(pr["b"], pr["b"]),
                "pair_a_key": arm.get(pr["a"], pr["a"]),
                "pair_b_key": arm.get(pr["b"], pr["b"]),
                "discordant_a_wins": pr["a_wins"],
                "discordant_b_wins": pr["b_wins"],
                "mcnemar_p_exact": pr.get("p_raw"),
                "mcnemar_p_holm": pr.get("p_holm"),
                "star": pr.get("star"),
                "between_tools_cochran_q": round(float(om["Q"]), 4),
                "between_tools_cochran_df": int(om["df"]),
                "between_tools_cochran_p": om["p"],
            })
    csv_path = args.report_dir / f"{args.out_stem}.csv"
    pd.DataFrame(rows).to_csv(csv_path, index=False)

    # ── TXT: the module's own report writer (same format as the raw-arm report) ──
    txt_path = args.report_dir / f"{args.out_stem}.txt"
    within_df = pc.aggregate_within_thresholds_by_depth(
        df, depths, pc.FINE_RMSD_THRESHOLDS, eq_df=eq_df, pb_valid_only=True,
        rmsd_col="rmsd")
    pc._write_within_by_depth_report(
        txt_path, within_df, f"{args.out_stem} (sidecar — no figure)", depths, slots,
        thr, "RMSD = symmetry-corrected heavy-atom, no superposition.", True, stats=stats)

    # Prepend the arm banner + append the printed-table block and the two numbers the
    # thesis currently quotes from the RAW arm.
    W = 90
    head = ["=" * W,
            "ARM: gnina/smina (post-hoc optimised) — the arm PRINTED in thesis Tables 5/6/7",
            "=" * W,
            f"AutoDock slot : {arm['autodock']}   (rank == optimized_rank, i.e. gnina re-ranking)",
            f"DiffDock slot : {arm['diffdock']}   (confidence rank, smina-minimised poses)",
            f"EquiBind slot : {arm['equibind']}   (ranked by gnina affinity)",
            f"Source        : {per_pose}",
            f"Generated by  : Scripts/Analysis/topk_recovery_gnina_arm.py",
            "Gate          : a complex counts at depth d when ANY of the tool's first d ranked",
            f"                poses is BOTH PoseBusters-valid AND within {thr:g} A (existence gate).",
            "NOTE          : the committed 18_topn_within_thresholds_pbvalid_depths* files are",
            "                the RAW AutoDock arm (78/145/150) — they are left untouched.",
            "", ""]

    tail = ["", "=" * W,
            "PRINTED TABLE VALUES (thesis Tables 5 / 6 / 7)",
            "=" * W,
            f"  {'variant':<24s}" + "".join(f"{'top-' + str(d):>14s}" for d in depths),
            "  " + "-" * (24 + 14 * len(depths))]
    for slot in slots:
        cells = "".join(
            f"{counts[(slot, d)][2]:>7d} ({100.0 * counts[(slot, d)][2] / counts[(slot, d)][0]:>4.1f}%)"
            for d in depths)
        tail.append(f"  {ARM_LABEL[slot]:<24s}" + cells)
    tail += ["", "  (cells = complexes with a PB-valid pose within "
             f"{thr:g} A among the first d ranked poses; n = "
             f"{counts[(slots[0], depths[0])][0]})", ""]

    # (a) within-tool Cochran Q for the AutoDock slot; (b) largest top-15 -> top-30 step.
    ad_q = headroom.get("autodock", {}).get("omnibus_across_depths")
    tail += ["=" * W,
             "THESIS CORRECTIONS — numbers currently quoted from the RAW AutoDock arm",
             "=" * W]
    if ad_q:
        tail.append(f"(a) within-tool Cochran Q, AutoDock slot ({arm['autodock']}):")
        tail.append(f"      Q = {float(ad_q['Q']):.2f}, df = {int(ad_q['df'])}, "
                    f"p = {_fmt_p(ad_q['p'])}")
        tail.append("      body_main.tex currently prints 'Vina Q = 135' — that is the RAW arm")
        tail.append("      (Q = 134.69, see 18_topn_within_thresholds_pbvalid_depths_report.txt).")
    if len(depths) >= 2:
        d0, d1 = depths[-2], depths[-1]
        steps_pp = []
        for slot in slots:
            n = counts[(slot, d1)][0]
            pp = 100.0 * (counts[(slot, d1)][2] - counts[(slot, d0)][2]) / n
            steps_pp.append((pp, slot))
        pp_max, slot_max = max(steps_pp)
        tail.append("")
        tail.append(f"(b) largest top-{d0} -> top-{d1} step across the three printed variants:")
        for pp, slot in sorted(steps_pp, reverse=True):
            k = counts[(slot, d1)][2] - counts[(slot, d0)][2]
            tail.append(f"      {ARM_LABEL[slot]:<24s} +{k:>2d} complexes = +{pp:.2f} pp")
        tail.append(f"      MAX = +{pp_max:.1f}% ({ARM_LABEL[slot_max]}).")
        tail.append("      body_main.tex currently prints 'at most 1.7%' — that is the RAW arm")
        tail.append("      (AutoDock raw +5 complexes = +1.65 pp).")
    tail.append("")
    tail.append("=" * W)

    txt_path.write_text("\n".join(head) + txt_path.read_text(encoding="utf-8")
                        + "\n".join(tail) + "\n", encoding="utf-8")

    print(f"wrote {csv_path}")
    print(f"wrote {txt_path}")
    for slot in slots:
        print(f"  {ARM_LABEL[slot]:<24s} " + "  ".join(
            f"top-{d}: {counts[(slot, d)][2]}/{counts[(slot, d)][0]}" for d in depths)
            + f"   Q={float(headroom[slot]['omnibus_across_depths']['Q']):.2f}")


if __name__ == "__main__":
    main()
