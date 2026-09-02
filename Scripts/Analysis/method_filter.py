"""Single-point docking-method exclusion for the Benchmark analysis suite.

Every Benchmark script that reads a per-pose / per-variant frame gained the same
two flags, wired at exactly ONE place per script — immediately after the frame is
assembled and before any consumer touches it:

    --exclude-methods  PAT[,PAT...]   drop these method keys
    --exclude-preset   NAME           drop a named, provenance-documented group

The point of routing every script through this module rather than editing each
keep-set is that the suite already carries several INDEPENDENT variant
enumeration points (``_select_presentation_tools`` and ``_all_variant_yield_specs``
in posebusters_pose_comparison.py are two; the per-family spec tables are more).
Excluding an arm by editing those one at a time reliably leaves it alive in some
figure nobody re-checked. A filter applied to the frame itself cannot.

Design rules this module enforces, each of which is a bug the suite has already
shipped at least once:

  * **Exact keys, never bare prefixes.** ``autodock`` matches the method
    ``autodock`` and nothing else. It does NOT match ``autodock_gnina`` or
    ``autodock_mgltools``. Prefix matching is what let ``_is_separate_autodock_arm``
    swallow every MGLTools arm and what makes ``_engine_family`` order-dependent.
    Where a family really is meant, ask for it with an explicit glob
    (``autodock_vinardo*``) so the intent is visible at the call site.
  * **A pattern that matches nothing is reported, never silently ignored.**
    The suite's recurring failure mode is a selector that prints
    "leaving unchanged" and carries on with the wrong arm.
  * **The on-disk per-pose cache is never touched.** Callers must apply this
    AFTER writing their cache, so exclusion stays a presentation-time decision and
    a later ``--reuse-cache`` run still has every variant available.
  * **Filtering happens before any key-rewriting step.** posebusters_validity_report
    folds every ``autodock_mgltools*`` key onto plain ``autodock``; a filter applied
    after that fold cannot distinguish the arms it is supposed to separate.

Run this file directly to verify a preset against the docking trees themselves:

    python Scripts/Analysis/method_filter.py --verify-preset meeko \\
        --config Scripts/Docking/Posebusters/posebusters_benchmark_full_protein_config.yaml
"""
from __future__ import annotations

import argparse
import fnmatch
import json
import sys
from pathlib import Path

__all__ = [
    "PRESETS",
    "add_method_filter_args",
    "resolve_patterns",
    "match_methods",
    "apply_method_filter",
    "apply_patterns",
    "write_filter_provenance",
    "fingerprint_tree",
]


# --------------------------------------------------------------------------- #
# presets
# --------------------------------------------------------------------------- #
# Each preset is (patterns, provenance). The provenance string is written into
# the sidecar JSON so a figure can be traced back to the reason its arms are
# missing. Keep patterns EXACT unless a family is genuinely meant, in which case
# use a trailing '*' so the breadth is legible.
PRESETS: dict[str, tuple[tuple[str, ...], str]] = {
    "meeko": (
        (
            "autodock",
            "autodock_smina",
            "autodock_gnina",
            "autodock_gnina_refinement",
            "autodock_vinardo*",
        ),
        "Arms whose ligands were converted to PDBQT by Meeko. Verified 2026-08-29 by "
        "fingerprinting the docked poses (not the staging directories, which are known "
        "to disagree with the poses beside them): Dockings/"
        "vina_results_full_protein_vina_scoring and .../vinardo_scoring carry "
        "'REMARK SMILES' Meeko headers, while every .../mgltools[_exh*] arm carries the "
        "'N active torsions' ADFRsuite header. Receptors are NOT part of this preset on "
        "the Benchmark whole-protein set: all four AutoDock trees there were prepared "
        "with ADFRsuite prepare_receptor (files named *_protein_mgl_tools.pdbqt). "
        "DiffDock, EquiBind and Uni-Dock2 hold no PDBQT at all and are unaffected. "
        "NOTE this preset removes the Vinardo scoring contrast entirely, because the "
        "only Vinardo tree on this dataset is Meeko-ligand.",
    ),
}


# --------------------------------------------------------------------------- #
# argparse wiring
# --------------------------------------------------------------------------- #
def add_method_filter_args(ap: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """Add the shared exclusion flags to a script's parser.

    Unset, both flags are inert and the script behaves exactly as before, so
    every committed command line still reproduces what it used to.
    """
    g = ap.add_argument_group("method exclusion (shared: Scripts/Analysis/method_filter.py)")
    g.add_argument(
        "--exclude-methods", default=None, metavar="PAT[,PAT...]",
        help="Comma-separated docking-method keys to drop from the analysis frame "
             "before anything consumes it. Matching is EXACT unless the pattern ends "
             "in '*' (e.g. 'autodock_vinardo*'), so 'autodock' never swallows "
             "'autodock_gnina' or 'autodock_mgltools'. Applied after the per-pose "
             "cache is written, so the cache keeps every variant.")
    g.add_argument(
        "--exclude-preset", default=None, choices=sorted(PRESETS),
        help="Drop a named, provenance-documented group of methods. "
             + "; ".join(f"'{k}': {v[0][0]}…" for k, v in sorted(PRESETS.items()))
             + ". Combines with --exclude-methods (union).")
    g.add_argument(
        "--exclude-strict", action="store_true",
        help="Treat a pattern that matches no method in this frame as a fatal error "
             "instead of a warning. Use it in scripts whose frame is expected to "
             "carry every arm.")
    return ap


def resolve_patterns(args) -> list[str]:
    """Union of --exclude-preset and --exclude-methods, order-preserving."""
    pats: list[str] = []
    preset = getattr(args, "exclude_preset", None)
    if preset:
        pats.extend(PRESETS[preset][0])
    raw = getattr(args, "exclude_methods", None)
    if raw:
        pats.extend(p.strip() for p in str(raw).split(",") if p.strip())
    seen: set[str] = set()
    return [p for p in pats if not (p in seen or seen.add(p))]


# --------------------------------------------------------------------------- #
# matching & application
# --------------------------------------------------------------------------- #
def match_methods(methods, patterns) -> tuple[list[str], list[str]]:
    """Split ``methods`` by ``patterns``. Returns (matched, unmatched_patterns).

    A pattern without '*' must equal a method key exactly. A pattern with '*' is
    an fnmatch glob. Both comparisons are case-sensitive, matching how the method
    keys are written throughout the suite.
    """
    methods = [str(m) for m in methods]
    matched: list[str] = []
    unmatched: list[str] = []
    for pat in patterns:
        if "*" in pat or "?" in pat or "[" in pat:
            hits = [m for m in methods if fnmatch.fnmatchcase(m, pat)]
        else:
            hits = [m for m in methods if m == pat]
        if hits:
            matched.extend(hits)
        else:
            unmatched.append(pat)
    seen: set[str] = set()
    return [m for m in matched if not (m in seen or seen.add(m))], unmatched


def apply_method_filter(df, column: str, args, *, label: str = "",
                        out_dir: Path | None = None, strict: bool | None = None):
    """Drop excluded methods from ``df`` and report exactly what went.

    Call this at ONE point per script: right after the frame is assembled and
    after any per-pose cache has been written, but BEFORE any step that renames
    or folds method keys.

    Returns ``df`` unchanged (and prints nothing) when no exclusion was
    requested, so existing command lines are unaffected.
    """
    return apply_patterns(
        df, column, resolve_patterns(args), label=label, out_dir=out_dir,
        strict=getattr(args, "exclude_strict", False) if strict is None else strict,
        preset=getattr(args, "exclude_preset", None))


def apply_patterns(df, column: str, patterns, *, label: str = "",
                   out_dir: Path | None = None, strict: bool = False,
                   preset: str | None = None):
    """Pattern-list form of :func:`apply_method_filter`.

    Library loaders (``posebusters_validity_report.load_and_score`` and the
    scripts that import it) take a plain list rather than an argparse namespace,
    so they stay usable from a notebook or another module.
    """
    if not patterns:
        return df
    if column not in df.columns:
        raise SystemExit(
            f"method-filter: column '{column}' is absent from the frame"
            f"{f' ({label})' if label else ''}; available: {list(df.columns)[:12]}")

    present = sorted(df[column].astype(str).unique())
    matched, unmatched = match_methods(present, patterns)

    tag = f"method-filter{f' [{label}]' if label else ''}"
    print(f"\n{tag}: patterns {list(patterns)}"
          f"{f' (preset {preset!r})' if preset else ''}")

    if unmatched:
        # Two benign causes, and one that is not: (1) a frame read BEFORE a
        # script's optimizer split carries base tree keys only, so variant-level
        # patterns legitimately miss; (2) a preset names a variant this campaign
        # never produced. The one that matters is a typo'd key silently leaving
        # the arm it meant to drop in the analysis, so always say which patterns
        # missed and let the reader decide.
        msg = (f"{tag}: {len(unmatched)} pattern(s) matched nothing in this frame: "
               f"{unmatched}\n{' ' * len(tag)}  present: {present}\n"
               f"{' ' * len(tag)}  (benign if this frame predates the optimizer split, "
               f"or if the campaign never produced that variant — check for a typo "
               f"otherwise)")
        if strict:
            raise SystemExit(msg + "\n  (--exclude-strict is set)")
        print(f"  note: {msg}", file=sys.stderr)

    if not matched:
        print(f"{tag}: nothing dropped — frame unchanged ({len(present)} methods).")
        return df

    n0, c0 = len(df), len(present)
    out = df[~df[column].astype(str).isin(matched)].reset_index(drop=True)
    kept = sorted(out[column].astype(str).unique())
    if out.empty:
        raise SystemExit(
            f"{tag}: the exclusion removed EVERY row ({n0:,} → 0). "
            f"Dropped {matched}; nothing left to analyse.")

    print(f"{tag}: dropped {len(matched)} of {c0} methods, "
          f"{n0 - len(out):,} of {n0:,} rows.")
    print(f"  dropped: {matched}")
    print(f"  kept   : {kept}")
    if out_dir is not None:
        write_filter_provenance(Path(out_dir), patterns, matched, kept, preset,
                                label=label)
    return out


def write_filter_provenance(out_dir: Path, patterns, dropped, kept,
                            preset: str | None, *, label: str = "") -> Path | None:
    """Record the exclusion beside the figures it shaped.

    Without this a reader cannot tell an arm that was excluded from an arm that
    was never run, which is the same ambiguity that made the stale-provenance
    audits expensive.
    """
    try:
        out_dir.mkdir(parents=True, exist_ok=True)
        path = out_dir / "method_filter.json"
        payload = {
            "label": label,
            "preset": preset,
            "preset_provenance": PRESETS[preset][1] if preset in PRESETS else None,
            "patterns": list(patterns),
            "dropped_methods": list(dropped),
            "kept_methods": list(kept),
        }
        path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        print(f"  provenance → {path}")
        return path
    except OSError as exc:                      # never fail an analysis over a sidecar
        print(f"  ! could not write method_filter.json: {exc}", file=sys.stderr)
        return None


# --------------------------------------------------------------------------- #
# preset verification against the docking trees
# --------------------------------------------------------------------------- #
_MEEKO_LIGAND_MARK = "REMARK SMILES"
_MGL_LIGAND_MARK = "active torsions"


def fingerprint_tree(tree: Path, sample: int = 8) -> dict:
    """Fingerprint a docking tree's LIGAND converter from its docked poses.

    Reads the poses themselves rather than the ``_staging`` directory, because
    staging is overwritten by any later run that shares the input path and is
    known to disagree with the poses sitting beside it.

    Receptor provenance is read from the prepared-receptor filename token, which
    is reliable on the whole-protein trees. AD4 atom-type profiling is reported
    alongside but not used as a verdict: the HD/NA direction separates the two
    converters only in relative terms, so a single file cannot decide it.
    """
    out = {"tree": str(tree), "exists": tree.is_dir(),
           "ligand": "unknown", "n_meeko": 0, "n_mgl": 0,
           "receptor": "unknown", "n_receptors": 0}
    if not out["exists"]:
        return out

    poses = []
    for p in tree.rglob("*out*.pdbqt"):
        poses.append(p)
        if len(poses) >= sample:
            break
    for p in poses:
        try:
            head = p.read_text(errors="ignore")[:4096]
        except OSError:
            continue
        if _MEEKO_LIGAND_MARK in head:
            out["n_meeko"] += 1
        elif _MGL_LIGAND_MARK in head:
            out["n_mgl"] += 1
    if out["n_meeko"] and not out["n_mgl"]:
        out["ligand"] = "meeko"
    elif out["n_mgl"] and not out["n_meeko"]:
        out["ligand"] = "mgltools"
    elif out["n_meeko"] or out["n_mgl"]:
        out["ligand"] = "mixed"
    elif poses:
        out["ligand"] = "indeterminate"
    else:
        out["ligand"] = "no-pdbqt"

    recs = []
    for r in tree.rglob("*protein*.pdbqt"):
        if "out" in r.name:
            continue
        recs.append(r)
        if len(recs) >= sample:
            break
    out["n_receptors"] = len(recs)
    if recs:
        tokens = {("meeko" if "_meeko" in r.name else
                   "mgltools" if "_mgl_tools" in r.name else "untokenised")
                  for r in recs}
        out["receptor"] = tokens.pop() if len(tokens) == 1 else "mixed"
    return out


def verify_preset(preset: str, config_path: Path, root: Path) -> int:
    """Check a preset against the trees the PoseBusters config actually reads.

    Exits non-zero when the preset and the fingerprints disagree, so this can be
    used as a regression gate rather than a report to be eyeballed.
    """
    import yaml                                   # only needed for verification

    patterns = PRESETS[preset][0]
    cfg = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    dirs = cfg.get("docking_directories") or {}
    if not dirs:
        print(f"no docking_directories in {config_path}", file=sys.stderr)
        return 2

    print(f"verify-preset '{preset}' against {config_path}")
    print(f"patterns: {list(patterns)}\n")
    hdr = f"{'config key':28s}{'ligand prep':>14}{'receptor prep':>15}   covered?"
    print(hdr)
    print("-" * len(hdr))

    bad = []
    for key in sorted(dirs):
        fp = fingerprint_tree(root / str(dirs[key]))
        # A base key is "covered" when the preset drops the base key itself; the
        # optimizer variants share that base and are listed alongside it.
        covered = bool(match_methods([key], patterns)[0])
        flag = ""
        if fp["ligand"] == "meeko" and not covered:
            flag = "  <-- MEEKO BUT NOT EXCLUDED"
            bad.append((key, fp["ligand"], covered))
        elif fp["ligand"] in ("mgltools", "no-pdbqt") and covered:
            flag = "  <-- EXCLUDED BUT NOT MEEKO"
            bad.append((key, fp["ligand"], covered))
        print(f"{key:28s}{fp['ligand']:>14}{fp['receptor']:>15}   "
              f"{'yes' if covered else 'no':<4}{flag}")

    if bad:
        print(f"\nFAIL: {len(bad)} tree(s) disagree with the preset: "
              f"{[b[0] for b in bad]}", file=sys.stderr)
        return 1
    print("\nOK: every Meeko-prepared tree is covered by the preset, and no "
          "non-Meeko tree is.")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Verify a method-exclusion preset against the docking trees.")
    ap.add_argument("--verify-preset", required=True, choices=sorted(PRESETS))
    ap.add_argument("--config", type=Path, required=True,
                    help="PoseBusters run config whose docking_directories map "
                         "method keys to trees.")
    ap.add_argument("--root", type=Path, default=Path.cwd(),
                    help="Repo root the config's relative paths resolve against.")
    a = ap.parse_args(argv)
    return verify_preset(a.verify_preset, a.config, a.root)


if __name__ == "__main__":
    raise SystemExit(main())
