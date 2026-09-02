#!/usr/bin/env python
"""Move stale material into obsolete/, in verifiable stages.

Every move is a rename within one filesystem, so it is instant and reversible.
Paths are MIRRORED: `X` becomes `obsolete/X`, so restoring is `mv obsolete/X X`.
Nothing is ever deleted.

    python Scripts/Analysis/cleanup_move.py --dry-run          # print, change nothing
    python Scripts/Analysis/cleanup_move.py --stage pandamap   # move one stage
    python Scripts/Analysis/cleanup_move.py --all              # all stages in order

Stages exist so verification can run between them. If the pipeline breaks after a
stage, that stage alone gets reverted:

    python Scripts/Analysis/cleanup_move.py --revert pandamap
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path("/home/manndo/master_dev")
OBS = ROOT / "obsolete"
LIST = Path("/tmp/claude-1000/-home-manndo-master-dev/8b0dc595-fb2f-487b-b009-c6f2dbe18a4f"
            "/scratchpad/move_list.tsv")

# Ordered least- to most-entangled, so a failure surfaces on cheap ground first.
STAGES = {
    "pandamap": lambda p: p.startswith("pandamap_results/"),
    "standalone": lambda p: p.split("/")[0] in {"posebuster_proved", "logs", ".git.corrupt-backup"},
    "dockings": lambda p: p.startswith("Dockings/"),
    "posebusters": lambda p: p.startswith("posebusters_results/"),
    "pocket": lambda p: p.startswith("pocket_results/"),
    "notebooks": lambda p: p.endswith(".ipynb") or ".ipynb.bak" in p or ".ipynb.snapshot" in p,
}
ORDER = ["pandamap", "standalone", "dockings", "posebusters", "pocket", "notebooks"]


def load() -> list[tuple[str, int, str]]:
    rows = []
    for line in LIST.read_text().splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        rows.append((parts[0], int(parts[1]) if len(parts) > 1 else 0,
                     parts[2] if len(parts) > 2 else ""))
    return rows


def stage_of(path: str) -> str | None:
    for name in ORDER:
        if STAGES[name](path):
            return name
    return None


def move_one(rel: str, dry: bool) -> tuple[bool, str]:
    src = ROOT / rel
    dst = OBS / rel
    if not src.exists():
        return False, "already gone"
    if dst.exists():
        return False, f"destination exists: obsolete/{rel}"
    if dry:
        return True, f"would move -> obsolete/{rel}"
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.move(str(src), str(dst))
    return True, f"moved -> obsolete/{rel}"


def revert_one(rel: str, dry: bool) -> tuple[bool, str]:
    src = OBS / rel
    dst = ROOT / rel
    if not src.exists():
        return False, "not in obsolete/"
    if dst.exists():
        return False, "destination already exists"
    if dry:
        return True, f"would restore <- obsolete/{rel}"
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.move(str(src), str(dst))
    return True, f"restored <- obsolete/{rel}"


def ensure_gitignored() -> str:
    gi = ROOT / ".gitignore"
    txt = gi.read_text()
    if any(l.strip() in {"obsolete/", "/obsolete/"} for l in txt.splitlines()):
        return "obsolete/ already ignored"
    gi.write_text(txt.rstrip("\n") +
                  "\n\n# Stale material parked by Scripts/Analysis/cleanup_move.py\nobsolete/\n")
    return "added obsolete/ to .gitignore"


def write_manifest(done: list[tuple[str, int, str]]) -> Path:
    OBS.mkdir(parents=True, exist_ok=True)
    p = OBS / "MANIFEST.md"
    prior = p.read_text() if p.exists() else ""
    if not prior:
        prior = (
            "# Parked material\n\n"
            "Moved out of `master_dev` by `Scripts/Analysis/cleanup_move.py`. This "
            "directory is gitignored.\n\n"
            "**Paths are mirrored.** A file that was at `X` is now at `obsolete/X`. "
            "Restore one with `mv obsolete/X X`, or a whole stage with "
            "`python Scripts/Analysis/cleanup_move.py --revert <stage>`.\n\n"
            "Nothing here was deleted. The keep rule was derived from the 58-stage "
            "registry, every stage command, every thesis config, and the "
            "unoverridden path defaults of the live analysis scripts, then attacked "
            "by six independent adversarial lenses before any move ran.\n\n"
            "| Path | MB | Why it moved |\n| --- | --- | --- |\n")
    rows = "".join(f"| `{p_}` | {n} | {w} |\n" for p_, n, w in done)
    p.write_text(prior + rows)
    return p


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--stage", choices=ORDER)
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--revert", choices=ORDER)
    args = ap.parse_args(argv)

    rows = load()
    if args.revert:
        sel = [r for r in rows if stage_of(r[0]) == args.revert]
        print(f"REVERT stage '{args.revert}': {len(sel)} entries")
        for rel, n, _ in sel:
            ok, msg = revert_one(rel, args.dry_run)
            print(f"  {'ok ' if ok else '** '}{rel}  {msg}")
        return 0

    if not (args.stage or args.all or args.dry_run):
        ap.error("choose --dry-run, --stage <name> or --all")

    sel = [r for r in rows if (args.stage is None or stage_of(r[0]) == args.stage)]
    total = sum(n for _, n, _ in sel)
    print(f"{'DRY RUN: ' if args.dry_run else ''}{len(sel)} entries, {total} MB "
          f"({total/1024:.1f} GB)\n")

    if not args.dry_run:
        print(f"  {ensure_gitignored()}\n")

    done, failed = [], []
    for name in ORDER:
        group = [r for r in sel if stage_of(r[0]) == name]
        if not group:
            continue
        print(f"  [{name}] {len(group)} entries, {sum(n for _, n, _ in group)} MB")
        for rel, n, why in group:
            ok, msg = move_one(rel, args.dry_run)
            if ok:
                done.append((rel, n, why))
            else:
                failed.append((rel, msg))
            print(f"    {'ok ' if ok else '** '}{n:6d} MB  {rel}")
            if not ok:
                print(f"           {msg}")
        print()

    if done and not args.dry_run:
        print(f"  manifest -> {write_manifest(done).relative_to(ROOT)}")
    moved = sum(n for _, n, _ in done)
    print(f"\n{len(done)} {'would move' if args.dry_run else 'moved'}, "
          f"{moved} MB ({moved/1024:.1f} GB); {len(failed)} skipped")
    for rel, msg in failed:
        print(f"  skipped {rel}: {msg}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
