#!/usr/bin/env python
"""Emit REGENERATE.md from the live stage registry.

WHY THIS IS GENERATED. The previous guide was maintained by hand and drifted. Its
figure index was numbered against a build that has since been retired, and its
output paths predated the matched-EquiBind migration, so following it landed you
on the wrong file for nearly every figure. Nothing announced that, because a
hand-written index has no way to notice the tree beneath it moved.

This version is derived. The ordering, the commands, the environments and the
determinism classes come from the same registry the notebook executes, so they
cannot disagree with what actually runs. The figure index is rebuilt by
checksumming each shipped asset against every candidate output on disk, so a
figure whose source moved shows up as moved rather than as a stale path.

Appendix A of the thesis describes what this file must contain: the output each
reported float was taken from, the exact command that rebuilds it, the conda
environment that command needs, the order dependent steps must run in, and an
explicit list of what could not be reconstructed. Each of those is a section below.

The narrative commentary of the old guide is worth keeping and is not
reproducible from the registry, so it stays alongside as REGENERATE_NOTES.md.
"""
from __future__ import annotations

import hashlib
import re
import time
from pathlib import Path

import yaml

ROOT = Path("/home/manndo/master_dev")
TEX = ROOT / "thesis_latex"
MEDIA = TEX / "media/media"
SPEC = ROOT / "Scripts/Analysis/thesis_expected_values.yaml"
NOTES = "REGENERATE_NOTES.md"
FINDINGS = "FINDINGS_2026-09-02.md"

# Roots searched when resolving a shipped figure back to the file that wrote it.
SEARCH_ROOTS = ("posebusters_results", "pandamap_results", "PoseBusters_Benchmark_Analysis")
# Preferred when one image has several byte-identical copies. These are the
# current generation; see the canonical-tree table the guide prints.
CANON_HINTS = ("_matched", "matched_equibind", "_orai_matched_root",
               "autodock_exhaustiveness_returns", "PoseBusters_Benchmark_Analysis")

ENVIRONMENTS = [
    ("vina", "/home/manndo/anaconda3/envs/vina/bin/python",
     "screening, statistics and every figure; Python 3.12.12, RDKit 2025.09.1"),
    ("diffdock", "/home/manndo/anaconda3/envs/diffdock/bin/python",
     "DiffDock inference, torch 2.4.0 / CUDA 12.4; also carries smina"),
    ("equibind", "/home/manndo/anaconda3/envs/equibind/bin/python",
     "EquiBind inference, torch 2.4.1 / CUDA 12.4; also carries smina"),
]

BINARIES = [
    ("AutoDock Vina", "/home/manndo/AutoDock-Vina/build/linux/release/vina",
     "v1.2.7-20-g93cdc3d-mod, the boron-patched build. /usr/bin/vina is stock 1.2.5 "
     "and must not be used."),
    ("gnina", "/home/manndo/docking_tools/gnina", "1.3.2"),
    ("smina", "<env>/bin/smina", "9 November 2017 build, on the Vina 1.1.2 scoring base"),
    ("ADFRsuite", "/home/manndo/ADFRsuite-1.0/bin/", "prepare_ligand and prepare_receptor"),
]


def _md5(p: Path) -> str | None:
    try:
        return hashlib.md5(p.read_bytes()).hexdigest()
    except Exception:
        return None


def _figure_order() -> list[tuple[str | None, str]]:
    """Every figure environment of the short build, in document order."""
    out: list[tuple[str | None, str]] = []
    for f in ("body_main_short.tex", "body_appendix_short.tex"):
        src = (TEX / f).read_text()
        for m in re.finditer(r"\\begin\{figure\}.*?\\end\{figure\}", src, re.S):
            blk = m.group(0)
            img = re.search(r"media/media/(image\d+)\.png", blk)
            cap = re.search(r"\\caption(?:\[([^\]]*)\])?\{", blk)
            short = cap.group(1) if cap and cap.group(1) else ""
            if not short:
                c2 = re.search(r"\\caption\{(.{0,70})", blk, re.S)
                short = (c2.group(1) if c2 else "").replace("\n", " ")
            short = re.sub(r"\\label.*", "", short).strip().rstrip("}").strip()
            if len(short) > 64:
                short = short[:64].rsplit(" ", 1)[0].rstrip(" ,.;") + "…"
            out.append((img.group(1) if img else None, short))
    return out


def _asset_index() -> dict[str, list[str]]:
    idx: dict[str, list[str]] = {}
    for base in SEARCH_ROOTS:
        d = ROOT / base
        if d.exists():
            for p in d.rglob("*.png"):
                h = _md5(p)
                if h:
                    idx.setdefault(h, []).append(str(p.relative_to(ROOT)))
    for p in ROOT.glob("*.png"):
        h = _md5(p)
        if h:
            idx.setdefault(h, []).append(str(p.relative_to(ROOT)))
    return idx


def _fmt_cmd(cmd) -> str:
    """One flag per continuation line, so the command reads as a recipe.

    Paths under the repo root are shortened, because every command runs from
    there and the absolute prefix is noise.
    """
    parts = []
    for c in cmd:
        s = str(c)
        if s.startswith(str(ROOT) + "/"):
            s = s[len(str(ROOT)) + 1:]
        parts.append('"' + s + '"' if (" " in s and not s.startswith('"')) else s)
    if not parts:
        return ""
    # The interpreter and whatever it runs stay together on the first line.
    head, rest = parts[0], parts[1:]
    if rest and not rest[0].startswith("-"):
        head, rest = head + " " + rest[0], rest[1:]
    lines, cur = [head], ""
    for part in rest:
        if part.startswith("-"):
            if cur:
                lines.append(cur)
            cur = "    " + part
        else:
            cur = (cur + " " + part) if cur else "    " + part
    if cur:
        lines.append(cur)
    return " \\\n".join(lines)


def build(reg, out: Path | None = None, verbose: bool = True) -> Path:
    """Write the guide. `reg` is the Registry the notebook built."""
    out = Path(out) if out else ROOT / "Scripts/Analysis/REGENERATE.md"
    L: list[str] = []
    e = L.append

    e("# REGENERATE.md — how every reported figure and table is produced")
    e("")
    e("> **This file is generated. Do not edit it by hand.**")
    e("> It is written by `Scripts/Analysis/regenerate_guide.py`, which the last cell")
    e("> of `Thesis_Reproduction.ipynb` runs. The ordering, commands and determinism")
    e("> classes come from the same stage registry the notebook executes, so they")
    e("> cannot disagree with what actually runs. The figure index is rebuilt by")
    e("> checksum on every run, so a figure whose source moved shows as moved rather")
    e("> than as a stale path.")
    e(">")
    e(f"> The previous hand-written guide is kept as [`{NOTES}`]({NOTES}) for its")
    e("> per-command commentary and its trap notes, which are not derivable from the")
    e("> registry. It drifted: its figure numbers came from a build that has since")
    e("> been retired and its paths predated the matched-EquiBind migration. That")
    e(f"> drift is what this generator exists to prevent. See [`{FINDINGS}`]({FINDINGS}).")
    e("")
    e(f"Generated {time.strftime('%Y-%m-%d %H:%M')} from {len(reg)} registered stages.")
    e("")
    e("Float numbers are the SHORT build's, read from the figure environments of")
    e("`body_main_short.tex` and `body_appendix_short.tex` in document order.")
    e("")
    e("Scope note. The analysis outputs referenced below live in the working tree,")
    e("not in the git index. `posebusters_results/`, `pandamap_results/` and")
    e("`PoseBusters_Benchmark_Analysis/figures/` are untracked, so byte-identical")
    e("here means identical to the file on disk, not to one recoverable from history.")
    e("")
    e("---")
    e("")

    # ---- 1. environments ----------------------------------------------------
    e("## 1. Environments and binaries")
    e("")
    e("Every command below runs from `/home/manndo/master_dev`.")
    e("")
    e("| Environment | Interpreter | Used for |")
    e("| --- | --- | --- |")
    for name, path, why in ENVIRONMENTS:
        e(f"| `{name}` | `{path}` | {why} |")
    e("")
    e("| Binary | Path | Version |")
    e("| --- | --- | --- |")
    for name, path, ver in BINARIES:
        e(f"| {name} | `{path}` | {ver} |")
    e("")
    e("Two patches live outside the repository and are load-bearing. The PoseBusters")
    e("`energy_ratio.py` exception-handler fix, without which `bust()` aborts on any")
    e("ligand UFF cannot parameterise, and the DiffDock seed hook in the local")
    e("checkout, without which its sampling is irreproducible. The notebook asserts")
    e("both before running anything.")
    e("")
    e("---")
    e("")

    # ---- 2. canonical trees -------------------------------------------------
    from repro_harness import CANONICAL, SUPERSEDED  # local import; needs Scripts/Analysis on path

    e("## 2. Canonical trees")
    e("")
    e("The trap this table exists to prevent is that the directories whose names read")
    e("like scratch are the current ones. Generation order, oldest first, is")
    e("`*_PRE_FR0` and `*_PRE_EXH128` and `*_UFFON_backup_*`, then the plain name,")
    e("then `_matched`.")
    e("")
    e("| Key | Directory |")
    e("| --- | --- |")
    for k, v in CANONICAL.items():
        e(f"| `{k}` | `{v.relative_to(ROOT)}` |")
    e("")
    e("Present on disk and deliberately NOT read:")
    e("")
    e("| Superseded | Replaced by |")
    e("| --- | --- |")
    for old, new in SUPERSEDED.items():
        if (ROOT / old).exists():
            e(f"| `{old}` | `{new}` |")
    e("")
    e("---")
    e("")

    # ---- 3. ordering --------------------------------------------------------
    stages = reg.topo()
    e("## 3. Order")
    e("")
    e("Declaration order below is a verified topological order of the dependency")
    e("graph, so running the sections in sequence can never read a file that a later")
    e("step writes. That failure mode is why this section is generated rather than")
    e("described.")
    e("")
    for section, group in reg.by_section().items():
        e(f"**{section}** — {len(group)} stages: " + ", ".join(f"`{s.name}`" for s in group))
        e("")
    e("---")
    e("")

    # ---- 4. the stages ------------------------------------------------------
    e("## 4. Stages")
    e("")
    e("Determinism is a property of the recorded run, not a policy. `bitexact` means")
    e("a re-run must reproduce the stored bytes. `verify` means the reported run was")
    e("an unseeded single draw and can be checked but not redrawn. `never` means")
    e("re-running destroys the provenance of a reported number.")
    e("")
    current = None
    for s in stages:
        if s.section != current:
            current = s.section
            e("### " + re.sub(r"^[0-9]+\. ", "", current))
            e("")
        e(f"#### `{s.name}` — {s.title}")
        e("")
        e(f"- **determinism** {s.determinism} · **cost** {s.cost}")
        if s.thesis:
            e(f"- **supports** {s.thesis}")
        if s.needs:
            e("- **needs** " + ", ".join(f"`{d}`" for d in s.needs))
        e("- **writes** " + ", ".join(f"`{o}`" for o in s.outputs))
        if s.cmd:
            e("")
            e("```bash")
            e(_fmt_cmd(s.cmd))
            e("```")
        elif s.determinism == "never":
            e("")
            e("No command. This stage refuses to run under any mode.")
        else:
            e("")
            e("No committed command; the outputs are staged by the arm that produced them.")
        if s.notes:
            e("")
            e(f"{s.notes}")
        e("")
    e("---")
    e("")

    # ---- 5. figure index ----------------------------------------------------
    e("## 5. Figure index")
    e("")
    e("Rebuilt by checksum on every generation. The source column is the file whose")
    e("md5 equals the shipped asset's, preferring the canonical tree where several")
    e("byte-identical copies exist.")
    e("")
    e("| Fig | Asset | Subject | Source |")
    e("| --- | --- | --- | --- |")
    idx = _asset_index()
    canon_roots = [str(v.relative_to(ROOT)) for v in CANONICAL.values()]
    unresolved = []
    for i, (img, cap) in enumerate(_figure_order(), 1):
        if not img:
            e(f"| {i} | — | {cap} | no image |")
            continue
        h = _md5(MEDIA / f"{img}.png")
        cands = idx.get(h, [])
        # Rank by how canonical the location is: a directory the registry names
        # beats a byte-identical sibling, which beats anything else.
        def _rank(c: str) -> tuple[int, int]:
            for i, root in enumerate(canon_roots):
                if c.startswith(root):
                    return (0, i)
            return (1 if any(k in c for k in CANON_HINTS) else 2, 0)
        pref = sorted(cands, key=_rank)
        if pref:
            extra = f" (+{len(pref)-1} identical)" if len(pref) > 1 else ""
            e(f"| {i} | `{img}` | {cap} | `{pref[0]}`{extra} |")
        else:
            e(f"| {i} | `{img}` | {cap} | **no source on disk** |")
            unresolved.append((i, img, cap))
    e("")
    e(f"{len(_figure_order()) - len(unresolved)} of {len(_figure_order())} figures "
      f"resolve to a file on disk.")
    e("")
    e("---")
    e("")

    # ---- 6. table index -----------------------------------------------------
    e("## 6. Table index")
    e("")
    spec = yaml.safe_load(SPEC.read_text()) if SPEC.exists() else {}
    rows = [(k, v) for k, v in spec.items() if k.startswith("table_")]
    if rows:
        e("Tables whose printed values are recomputed and asserted on every run by")
        e("`thesis_assertions.py`. A failure there is reported, never silently fixed.")
        e("")
        e("| Table | Source | Location in the thesis |")
        e("| --- | --- | --- |")
        for k, v in sorted(rows):
            num = k.split("_", 1)[1]
            e(f"| {num} | `{v.get('input', '')}` | {v.get('source', '')} |")
        e("")
    e("Tables 7, 9, 11, 13, 15 and 23 are manual, compiled from the literature or")
    e("from `DOCKING_PROTOCOL.md`. The remainder are not individually asserted; their")
    e("sources are given per stage in section 4.")
    e("")
    e("---")
    e("")

    # ---- 7. what cannot be reconstructed ------------------------------------
    e("## 7. What cannot be reconstructed")
    e("")
    e("Stated explicitly, as Appendix A requires.")
    e("")
    for i, img, cap in unresolved:
        e(f"- **Figure {i}** (`{img}`), {cap}. Hand-composed PyMOL render, present")
        e("  nowhere in the repository. `pymol_pose_cluster_scene.py` builds the scene")
        e("  files, but the camera, colouring and render were interactive.")
    e("- **Figure 17**, the docking flow chart, is author-drawn in `Flow Charts.ipynb`.")
    e("- The **two unseeded DiffDock runs**, the benchmark and the Orai control. Their")
    e("  poses can be verified against the stored trees but cannot be redrawn.")
    e("- The **Fr0 receptor**, the **prepared ligand PDBQTs** and")
    e("  **`exhaustiveness_arm_status.json`**, for the reasons in their stage notes.")
    e("")
    e("---")
    e("")
    e(f"Findings and evidence: [`{FINDINGS}`]({FINDINGS}).  ")
    e(f"Legacy commentary: [`{NOTES}`]({NOTES}).  ")
    e("Pipeline design: `../../THESIS_REPRODUCTION.md`.")
    e("")

    out.write_text("\n".join(L))
    if verbose:
        print(f"wrote {out.relative_to(ROOT)}")
        print(f"  {len(reg)} stages, {len(_figure_order())} figures, "
              f"{len(unresolved)} without a source on disk")
    return out
