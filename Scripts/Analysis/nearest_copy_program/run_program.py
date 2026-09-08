#!/usr/bin/env python
"""Phase 3 driver of the nearest-copy endpoint program (plan v2, Phase 3.1b and 3.2).

Reads the canonical stage commands from Scripts/Analysis/REGENERATE.md section 4
(and bench_metal_stratum from _build_reproduction_notebook.py, which REGENERATE.md
does not carry yet), rewrites every canonical output path to its ``_nearest``
counterpart (README.md of this folder), points every per-pose or report-dir INPUT
at the nearest hub table, adds the flags the Phase 2 tasks introduce
(--crystal-copies any, the program PandaMap config with --refresh-crystal,
--cascade-pins-from) and appends the two new stages of this folder.

Safety, enforced in code before anything runs (also under --dry-run):
  * every output argument must contain ``_nearest`` (or the driver refuses)
  * no emitted output path may equal, contain or lie inside any path that
    REGENERATE.md lists as a canonical output (**writes** lines) or as a
    canonical tree (section 2)
  * no canonical directory is ever passed as an output
  * the nearest hub table must exist with a manifest saying reference_convention nearest
  * under a real run the target script must already accept every added flag

Stages run in the order of plan Phase 3.1b then 3.2. Every command is logged with
its wall time and exit code to PROGRAM_LOG.md in this folder (append). The driver
stops at the first non-zero exit unless --keep-going is given.

Usage:
  run_program.py --dry-run                 print every command, run nothing
  run_program.py --dry-run --only bench_clusters bench_itt
  run_program.py --list                    stage names in run order
  run_program.py                           run everything, in order
  run_program.py --only bench_topk         run one stage

The stage bench_pandamap expects the program config
Scripts/Analysis/nearest_copy_program/07_benchmark_pandamap_nearest.yaml to exist
(Phase 2.2 writes it); its output_dir must carry ``_nearest``.
"""
from __future__ import annotations

import argparse
import ast
import json
import os
import re
import shlex
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

ROOT = Path("/home/manndo/master_dev")
HERE = Path(__file__).resolve().parent
ANA = ROOT / "Scripts" / "Analysis"
REGEN = ANA / "REGENERATE.md"
BUILDER = ANA / "_build_reproduction_notebook.py"
LOG = HERE / "PROGRAM_LOG.md"
VINA_PY = "/home/manndo/anaconda3/envs/vina/bin/python"

CANON_REPORT = "posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report"
NEAR_REPORT = CANON_REPORT + "_nearest"
NEAR_PPM = f"{NEAR_REPORT}/per_pose_metrics.csv"
NEAR_MANIFEST = f"{NEAR_REPORT}/per_pose_metrics.manifest.json"
NEAR_IDS = f"{NEAR_REPORT}/analysed_cohort_ids.txt"
NEAR_CASCADE = f"{NEAR_REPORT}/pose_validity_cascade.csv"
CANON_PANDAMAP_YAML = "Scripts/Docking/configs_thesis/07_benchmark_pandamap.yaml"
NEAR_PANDAMAP_YAML = "Scripts/Analysis/nearest_copy_program/07_benchmark_pandamap_nearest.yaml"
NEAR_PANDAMAP_DIR = "pandamap_results/benchmark_matched_equibind_nearest"

# canonical prefix -> nearest prefix (README table); applied per whole token,
# longest canonical first so a parent never shadows a child
PATH_MAP = {
    CANON_REPORT: NEAR_REPORT,
    "posebusters_results/benchmark_matched_equibind/dock/validity_report_mgltools":
        "posebusters_results/benchmark_matched_equibind/dock/validity_report_mgltools_nearest",
    "posebusters_results/cluster_crystal_pocket_matched_equibind/autodock_mgltools_exh128_gnina__diffdock_smina_allposes":
        "posebusters_results/cluster_crystal_pocket_matched_equibind/autodock_mgltools_exh128_gnina__diffdock_smina_allposes_nearest",
    "pandamap_results/benchmark_matched_equibind": NEAR_PANDAMAP_DIR,
    CANON_PANDAMAP_YAML: NEAR_PANDAMAP_YAML,
    "posebusters_results/benchmark_matched_equibind/docking_effort_gnina_v2_charged":
        "posebusters_results/benchmark_matched_equibind/docking_effort_gnina_v2_charged_nearest",
    "posebusters_results/benchmark_matched_equibind/docking_effort_gnina_v2":
        "posebusters_results/benchmark_matched_equibind/docking_effort_gnina_v2_nearest",
    "posebusters_results/autodock_exhaustiveness_returns": "posebusters_results/autodock_exhaustiveness_returns_nearest",
    "PoseBusters_Benchmark_Analysis/ligand_difficulty": "PoseBusters_Benchmark_Analysis/ligand_difficulty_nearest",
    "PoseBusters_Benchmark_Analysis/receptor_difficulty": "PoseBusters_Benchmark_Analysis/receptor_difficulty_nearest",
    "PoseBusters_Benchmark_Analysis/smina_rerank": "PoseBusters_Benchmark_Analysis/smina_rerank_nearest",
    "posebusters_results/metal_stratum": "posebusters_results/metal_stratum_nearest",
}
_MAP_ORDER = sorted(PATH_MAP, key=len, reverse=True)

# Stage table in run order (plan 3.1b, then 3.2). Fields:
#   flags   : argument names whose VALUE is an output location (checked for _nearest)
#   extra   : arguments appended to the canonical command
#   outputs : output paths the stage writes (after rewriting), for the canonical-collision assert
#   inputs  : files that must exist before a real run
#   needs_flags : {script: [flag]} the target script must accept (real run only)
STAGES: list[dict] = [
    dict(name="bench_exhaustiveness_returns", plan="3.1b", flags=["--out-dir"],
         extra=["--cascade-pins-from", NEAR_CASCADE],
         outputs=["posebusters_results/autodock_exhaustiveness_returns_nearest"],
         inputs=[NEAR_PPM, NEAR_CASCADE],
         needs_flags={"autodock_exhaustiveness_returns.py": ["--cascade-pins-from"]}),
    dict(name="bench_validity_report", plan="3.2", flags=["--out-dir"],
         outputs=["posebusters_results/benchmark_matched_equibind/dock/validity_report_mgltools_nearest"],
         inputs=[NEAR_PPM]),
    dict(name="bench_filmstrip", plan="3.2", flags=["--report-dir"], outputs=[NEAR_REPORT], inputs=[NEAR_PPM]),
    dict(name="bench_topk", plan="3.2", flags=["--report-dir"], outputs=[NEAR_REPORT], inputs=[NEAR_PPM]),
    dict(name="bench_endpoint_diagnostics", plan="3.2", flags=["--csv"], outputs=[NEAR_REPORT], inputs=[NEAR_PPM]),
    dict(name="bench_optimization_benefit", plan="3.2", flags=["--report-dir"], outputs=[NEAR_REPORT],
         inputs=[NEAR_PPM]),
    dict(name="bench_ligand_difficulty", plan="3.2", flags=["--out-dir"],
         outputs=["PoseBusters_Benchmark_Analysis/ligand_difficulty_nearest"], inputs=[NEAR_PPM]),
    dict(name="bench_receptor_difficulty", plan="3.2", flags=["--out-dir"],
         outputs=["PoseBusters_Benchmark_Analysis/receptor_difficulty_nearest"], inputs=[NEAR_PPM]),
    dict(name="bench_clusters", plan="3.2", flags=["--out-dir"], extra=["--crystal-copies", "any", "--force"],
         outputs=["posebusters_results/cluster_crystal_pocket_matched_equibind/"
                  "autodock_mgltools_exh128_gnina__diffdock_smina_allposes_nearest"],
         inputs=[NEAR_PPM],
         needs_flags={"pose_cluster_crystal_pocket_report.py": ["--crystal-copies", "--force"]}),
    dict(name="bench_effort_charged", plan="3.2", flags=["--out-dir"],
         outputs=["posebusters_results/benchmark_matched_equibind/docking_effort_gnina_v2_charged_nearest"],
         inputs=[NEAR_PPM]),
    dict(name="bench_effort_elapsed", plan="3.2", flags=["--out-dir"],
         outputs=["posebusters_results/benchmark_matched_equibind/docking_effort_gnina_v2_nearest"],
         inputs=[NEAR_PPM]),
    dict(name="bench_pandamap", plan="3.2", flags=["-c"], extra=["--refresh-crystal"],
         outputs=[NEAR_PANDAMAP_DIR], inputs=[NEAR_PPM, NEAR_PANDAMAP_YAML], yaml_check=True,
         needs_flags={"run_pandamap.py": ["--refresh-crystal"]}),
    dict(name="bench_pandamap_report", plan="3.2", flags=["-c", "--in-dir"],
         outputs=[f"{NEAR_PANDAMAP_DIR}/report"],
         inputs=[NEAR_PPM, NEAR_PANDAMAP_YAML, f"{NEAR_PANDAMAP_DIR}/pandamap_interactions.csv"], yaml_check=True),
    dict(name="bench_interaction_audit", plan="3.2", flags=["--config", "--in-dir"],
         outputs=[f"{NEAR_PANDAMAP_DIR}/report"],
         inputs=[NEAR_PPM, NEAR_PANDAMAP_YAML, NEAR_IDS, f"{NEAR_PANDAMAP_DIR}/pandamap_interactions.csv"],
         yaml_check=True),
    dict(name="bench_diffdock_rerank", plan="3.2", flags=["--out-dir"],
         outputs=["PoseBusters_Benchmark_Analysis/smina_rerank_nearest"], inputs=[NEAR_PPM]),
    dict(name="bench_metal_stratum", plan="3.2 (new)", flags=["--out-dir"], source="builder",
         outputs=["posebusters_results/metal_stratum_nearest"], inputs=[NEAR_PPM, NEAR_IDS]),
    dict(name="bench_itt", plan="3.2 (new, 2.6)", flags=["--out-dir"], source="program",
         cmd=[VINA_PY, "Scripts/Analysis/nearest_copy_program/itt_dropped_complexes.py",
              "--out-dir", "posebusters_results/itt_nearest", "--per-pose-csv", NEAR_PPM],
         outputs=["posebusters_results/itt_nearest"], inputs=[NEAR_PPM]),
    dict(name="bench_reference_convention", plan="3.2 (new, 2.7)", flags=["--out-dir"], source="program",
         cmd=[VINA_PY, "Scripts/Analysis/nearest_copy_program/reference_convention_table.py",
              "--per-pose-csv", NEAR_PPM, "--out-dir", "posebusters_results/reference_convention_nearest"],
         outputs=["posebusters_results/reference_convention_nearest"], inputs=[NEAR_PPM]),
]
INPUT_PATH_FLAGS = {"--per-pose-csv", "--per-pose-metrics", "--metrics", "--ids-file", "--report-dir", "--csv",
                    "--in-dir", "-c", "--config", "--out-dir", "--out", "--cascade-pins-from"}


# ---------------------------------------------------------------------------
# REGENERATE.md parsing
# ---------------------------------------------------------------------------
EMPTY_VALUE_OPTIONS = {"--unidock2-dir", "--unidock-dir"}


def parse_regenerate(text: str) -> tuple[dict, list[str], list[str]]:
    """Return ({stage: {"cmd": tokens, "writes": [paths]}}, all canonical writes, canonical trees)."""
    stages: dict[str, dict] = {}
    all_writes: list[str] = []
    blocks = re.split(r"^#### `", text, flags=re.M)[1:]
    for blk in blocks:
        name = blk.split("`", 1)[0]
        writes = []
        m = re.search(r"^- \*\*writes\*\* (.+)$", blk, flags=re.M)
        if m:
            writes = re.findall(r"`([^`]+)`", m.group(1))
        cm = re.search(r"```bash\n(.*?)```", blk, flags=re.S)
        cmd = shlex.split(cm.group(1).replace("\\\n", " ")) if cm else []
        stages[name] = {"cmd": cmd, "writes": writes}
        all_writes.extend(writes)
    trees = []
    sec2 = text.split("## 2. Canonical trees", 1)[1].split("## 3.", 1)[0]
    for line in sec2.splitlines():
        m = re.match(r"^\| `[^`]+` \| `([^`]+)` \|", line)
        if m:
            trees.append(m.group(1))
    return stages, sorted(set(all_writes)), trees


def parse_builder_stage(name: str) -> dict:
    """Extract cmd and outputs of a stage from _build_reproduction_notebook.py by text,
    without importing or executing the builder (which writes the notebook)."""
    text = BUILDER.read_text()
    start = text.find(f'name="{name}"')
    if start < 0:
        sys.exit(f"ERROR: stage {name} not found in {BUILDER}")
    end = text.find("reg.add(Stage(", start)
    body = text[start:end if end > 0 else None]
    ns = {"VINA_PY": VINA_PY, "ANA": Path("Scripts/Analysis"), "BENCH_REPORT": CANON_REPORT}

    def grab(field: str):
        m = re.search(rf"\b{field}=(\[.*?\])\s*,\s*\n", body, flags=re.S)
        if not m:
            sys.exit(f"ERROR: field {field} of {name} not found in the builder")
        return eval(m.group(1), {"__builtins__": {}}, ns)  # noqa: S307 (our own file, fixed names)

    cmd = [str(t) for t in grab("cmd")]
    outputs = [str(t) for t in grab("outputs")]
    return {"cmd": cmd, "writes": outputs}


# ---------------------------------------------------------------------------
# rewriting and checks
# ---------------------------------------------------------------------------
def rewrite_token(tok: str) -> str:
    for canon in _MAP_ORDER:
        if tok == canon or tok.startswith(canon + "/"):
            return PATH_MAP[canon] + tok[len(canon):]
    return tok


def _norm(p: str) -> str:
    return os.path.normpath(p).rstrip("/")


def is_canonical_output(path: str, canon_writes: list[str], canon_trees: list[str]) -> str | None:
    """Return the offending canonical path when ``path`` would touch one.

    Canonical OUTPUTS (REGENERATE **writes**): refused when equal, when ``path``
    lies inside one that is a directory, or when one lies inside ``path``.
    Canonical TREES (REGENERATE section 2): refused when equal or when the tree
    lies inside ``path``. A path INSIDE a tree is allowed only when the first
    component below the tree carries ``_nearest`` (the README places
    ``dock/pose_comparison_report_nearest`` inside the ``benchmark_pb`` tree);
    anything else below a canonical tree is refused.
    """
    p = _norm(path)
    for c in canon_writes:
        c = _norm(c)
        if p == c or p.startswith(c + "/") or c.startswith(p + "/"):
            return c
    for c in canon_trees:
        c = _norm(c)
        if p == c or c.startswith(p + "/"):
            return c
        if p.startswith(c + "/"):
            leaf = p[len(c) + 1:].split("/")[0]
            if "_nearest" not in leaf:
                return c
    return None


def build_command(stage: dict, regen: dict, canon_writes: list[str], canon_trees: list[str]) -> list[str]:
    src = stage.get("source", "regenerate")
    if src == "regenerate":
        if stage["name"] not in regen:
            sys.exit(f"ERROR: {stage['name']} not in {REGEN}")
        base = regen[stage["name"]]["cmd"]
    elif src == "builder":
        base = parse_builder_stage(stage["name"])["cmd"]
    else:
        base = list(stage["cmd"])
    cmd = [rewrite_token(t) for t in base] + list(stage.get("extra", []))
    # REGENERATE.md renders empty-string option values (builder EFFORT_COMMON passes
    # "--unidock2-dir", "" and "--unidock-dir", "") as bare flags; re-insert the empty
    # value whenever one of these value-taking options is followed by another option
    # or ends the command, otherwise argparse rejects the stage.
    repaired = []
    for i, tok in enumerate(cmd):
        repaired.append(tok)
        if tok in EMPTY_VALUE_OPTIONS and (i + 1 == len(cmd) or cmd[i + 1].startswith("-")):
            repaired.append("")
    cmd = repaired

    # (1) every value of an output flag must carry _nearest and must not be canonical
    for i, tok in enumerate(cmd[:-1]):
        if tok in stage["flags"]:
            val = _norm(cmd[i + 1])
            if "_nearest" not in val:
                sys.exit(f"REFUSED {stage['name']}: output {tok} {val!r} does not contain _nearest")
            hit = is_canonical_output(val, canon_writes, canon_trees)
            if hit:
                sys.exit(f"REFUSED {stage['name']}: output {tok} {val!r} collides with canonical {hit!r}")
    # (2) declared outputs: _nearest and no collision with any REGENERATE canonical output
    for out in map(_norm, stage["outputs"]):
        if "_nearest" not in out:
            sys.exit(f"REFUSED {stage['name']}: declared output {out!r} does not contain _nearest")
        hit = is_canonical_output(out, canon_writes, canon_trees)
        assert hit is None, f"{stage['name']}: emitted output {out!r} equals or overlaps canonical {hit!r}"
    # (3) no canonical path survives anywhere among path-like arguments that name a
    #     report dir or per-pose table (inputs must be the nearest tree too)
    for i, tok in enumerate(cmd[:-1]):
        if tok in INPUT_PATH_FLAGS:
            val = cmd[i + 1]
            if val == CANON_REPORT or val.startswith(CANON_REPORT + "/"):
                sys.exit(f"REFUSED {stage['name']}: {tok} still points at the canonical report {val!r}")
            if val == CANON_PANDAMAP_YAML:
                sys.exit(f"REFUSED {stage['name']}: canonical PandaMap config passed")
    return cmd


def check_yaml(stage: dict, strict: bool) -> list[str]:
    """The program PandaMap config must write into a _nearest directory and read the nearest table."""
    notes = []
    p = ROOT / NEAR_PANDAMAP_YAML
    if not p.exists():
        msg = f"{NEAR_PANDAMAP_YAML} does not exist yet (Phase 2.2 writes it)"
        if strict:
            sys.exit(f"REFUSED {stage['name']}: {msg}")
        return [f"WARNING: {msg}"]
    try:
        import yaml  # type: ignore
        cfg = yaml.safe_load(p.read_text()) or {}
    except Exception as exc:  # noqa: BLE001
        cfg = {}
        notes.append(f"WARNING: could not parse {NEAR_PANDAMAP_YAML}: {exc}")
        for line in p.read_text().splitlines():
            if line.startswith("output_dir:"):
                cfg["output_dir"] = line.split(":", 1)[1].split("#")[0].strip()
            if line.startswith("per_pose_metrics:"):
                cfg["per_pose_metrics"] = line.split(":", 1)[1].split("#")[0].strip()
    out_dir = str(cfg.get("output_dir", ""))
    if "_nearest" not in out_dir:
        sys.exit(f"REFUSED {stage['name']}: {NEAR_PANDAMAP_YAML} output_dir {out_dir!r} lacks _nearest")
    if out_dir.rstrip("/") != NEAR_PANDAMAP_DIR:
        notes.append(f"WARNING: yaml output_dir {out_dir!r} differs from the README layout {NEAR_PANDAMAP_DIR!r}")
    ppm = str(cfg.get("per_pose_metrics", ""))
    if ppm and ppm != NEAR_PPM:
        notes.append(f"WARNING: yaml per_pose_metrics {ppm!r} is not the nearest table {NEAR_PPM!r}")
    return notes


def check_flags(stage: dict) -> list[str]:
    missing = []
    for script, flags in (stage.get("needs_flags") or {}).items():
        text = (ANA / script).read_text() if (ANA / script).exists() else ""
        for f in flags:
            if f'"{f}"' not in text and f"'{f}'" not in text:
                missing.append(f"{script} lacks {f}")
    return missing


def check_hub_table(strict: bool) -> list[str]:
    notes = []
    ppm, man = ROOT / NEAR_PPM, ROOT / NEAR_MANIFEST
    if not ppm.exists() or not man.exists():
        msg = f"nearest hub table or manifest missing ({NEAR_PPM}, {NEAR_MANIFEST})"
        if strict:
            sys.exit(f"REFUSED: {msg}")
        return [f"WARNING: {msg}; the rebuild has not finished"]
    try:
        conv = json.loads(man.read_text()).get("reference_convention")
    except Exception as exc:  # noqa: BLE001
        sys.exit(f"REFUSED: cannot read {NEAR_MANIFEST}: {exc}")
    if conv != "nearest":
        sys.exit(f"REFUSED: {NEAR_MANIFEST} reference_convention is {conv!r}, expected 'nearest'")
    return notes


def ensure_cohort_ids() -> None:
    """analysed_cohort_ids.txt is derived from the hub table (REGENERATE: 'derived from
    per_pose_metrics rather than written by hand'); write it into the NEAREST report dir only."""
    target = ROOT / NEAR_IDS
    assert "_nearest" in str(target)
    if target.exists():
        return
    import pandas as pd
    ids = sorted(pd.read_csv(ROOT / NEAR_PPM, usecols=["protein"])["protein"].unique())
    target.write_text("\n".join(ids) + "\n")
    print(f"  derived {NEAR_IDS} ({len(ids)} ids)")


def fmt_cmd(cmd: list[str]) -> str:
    """REGENERATE.md style: interpreter and script first, then one flag (with all
    of its values) per continuation line."""
    q = [shlex.quote(t) for t in cmd]
    lines, cur = [" ".join(q[:2])], []
    for tok in q[2:]:
        if tok.startswith("-") and cur:
            lines.append(" ".join(cur))
            cur = []
        cur.append(tok)
    if cur:
        lines.append(" ".join(cur))
    return " \\\n    ".join(lines)


def log_line(stage: str, cmd: list[str], wall: float, code: int) -> None:
    """Append one row to PROGRAM_LOG.md in its existing three-column table
    (When | Step | Result); the header is written only when the file is new."""
    new = not LOG.exists()
    now = datetime.now()
    result = ("REFUSED, inputs missing" if code == -1 else f"exit {code} after {wall:.1f} s")
    with LOG.open("a") as fh:
        if new:
            fh.write("# Nearest-copy program log\n\n| When | Step | Result |\n|---|---|---|\n")
        fh.write(f"| {now.strftime('%Y-%m-%d %H:%M:%S')} | run_program.py `{stage}`: "
                 f"`{' '.join(shlex.quote(t) for t in cmd)}` | {result} |\n")


# ---------------------------------------------------------------------------
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dry-run", action="store_true", help="print the commands, run nothing, write no log")
    ap.add_argument("--only", nargs="+", default=None, metavar="STAGE", help="restrict to these stages")
    ap.add_argument("--list", action="store_true", help="print the stage names in run order")
    ap.add_argument("--keep-going", action="store_true", help="continue after a failing stage")
    args = ap.parse_args(argv)
    os.chdir(ROOT)

    names = [s["name"] for s in STAGES]
    if args.list:
        print("\n".join(names))
        return 0
    selected = STAGES
    if args.only:
        unknown = [n for n in args.only if n not in names]
        if unknown:
            sys.exit(f"ERROR: unknown stage(s) {unknown}; known: {names}")
        selected = [s for s in STAGES if s["name"] in set(args.only)]

    regen, canon_writes, canon_trees = parse_regenerate(REGEN.read_text())
    # the builder-only stage's canonical outputs count as canonical too
    canon_writes = sorted(set(canon_writes) | set(parse_builder_stage("bench_metal_stratum")["writes"]))
    strict = not args.dry_run
    notes = check_hub_table(strict)
    for n in notes:
        print(n)

    plan = []
    for st in selected:
        cmd = build_command(st, regen, canon_writes, canon_trees)
        st_notes = []
        if st.get("yaml_check"):
            st_notes += check_yaml(st, strict)
        miss = check_flags(st)
        if miss:
            if strict:
                sys.exit(f"REFUSED {st['name']}: target script does not accept the required flag(s): {miss}")
            st_notes += [f"WARNING: {m} (Phase 2 task pending)" for m in miss]
        plan.append((st, cmd, st_notes))

    # every emitted output, across all stages, against every canonical output (belt and braces)
    emitted = [_norm(o) for st, _, _ in plan for o in st["outputs"]]
    for o in emitted:
        assert is_canonical_output(o, canon_writes, canon_trees) is None, o
        assert "_nearest" in o, o
    print(f"# {len(plan)} stage(s); {len(canon_writes)} canonical outputs and {len(canon_trees)} canonical trees "
          f"checked; no emitted output touches any of them")

    if args.dry_run:
        for st, cmd, st_notes in plan:
            print(f"\n## {st['name']}  (plan {st['plan']}; outputs: {', '.join(st['outputs'])})")
            for n in st_notes:
                print(f"#   {n}")
            missing_in = [p for p in st.get("inputs", []) if not (ROOT / p).exists()]
            if missing_in:
                print(f"#   inputs not yet present: {', '.join(missing_in)}")
            print(fmt_cmd(cmd))
        return 0

    failures = 0
    for st, cmd, st_notes in plan:
        for n in st_notes:
            print(n)
        missing_in = [p for p in st.get("inputs", []) if p != NEAR_IDS and not (ROOT / p).exists()]
        if missing_in:
            print(f"REFUSED {st['name']}: inputs missing: {missing_in}")
            log_line(st["name"], cmd, 0.0, -1)
            failures += 1
            if not args.keep_going:
                return 2
            continue
        if NEAR_IDS in st.get("inputs", []):
            ensure_cohort_ids()
        print(f"\n== {st['name']} ==\n{fmt_cmd(cmd)}", flush=True)
        t0 = time.time()
        proc = subprocess.run(cmd, cwd=ROOT)
        wall = time.time() - t0
        log_line(st["name"], cmd, wall, proc.returncode)
        print(f"== {st['name']} exit {proc.returncode} after {wall:.1f} s")
        if proc.returncode != 0:
            failures += 1
            if not args.keep_going:
                return proc.returncode
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
