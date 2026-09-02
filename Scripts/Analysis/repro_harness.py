#!/usr/bin/env python
"""Stage registry and environment guard for Thesis_Reproduction.ipynb.

WHY THIS EXISTS. The notebook that produced the thesis results,
``Master_Docking_AD_Full_Protein.ipynb``, is 125 free-standing cells whose order
on the page is not their order of execution: three Orai analysis cells read a CSV
that a cell forty positions later writes. Three of its cells rewrite YAML configs
on disk, and one of them flips ``uff_minimize`` back to true, which is the setting
the reported EquiBind arm depends on being false. Running it top to bottom
corrupts its own inputs.

This module replaces that with a declared dependency graph. Each stage names what
it produces, what it needs, what it costs, and whether it can be reproduced at
all. The driver topologically sorts them, so a stage can never run before its
inputs exist, and the ordering bug cannot recur.

THREE DETERMINISM CLASSES, which are a property of the recorded runs rather than
a policy choice:

``BITEXACT``
    Seeded end to end, so a re-run must reproduce the stored bytes. AutoDock Vina
    at seed 42, the gnina and smina refinements at their own seeds, EquiBind's
    thirty fixed RDKit seeds and per-job model seeding, the Orai experimental
    DiffDock run, PoseBusters, and every analysis and statistic. The statistics
    are seeded too: the bootstraps take 2,000 resamples and the permutation tests
    10,000, both at fixed seeds.

``VERIFY``
    The reported run was an unseeded single draw and cannot be redrawn. This is
    the benchmark DiffDock run and the Orai control DiffDock run. Established by
    evidence rather than assumption: the seed patch in the local DiffDock
    checkout prints ``[repro] seeded`` when it fires, and that marker appears in
    0 of the stored benchmark logs, 0 of the Orai control logs and 24 of the Orai
    experimental logs. Appendix B records the same asymmetry. These stages verify
    the stored tree and refuse to overwrite it.

``NEVER``
    Not reproducible even in principle, and destructive to re-run. The Fr0
    receptor passed through an unguarded OpenMM minimisation that differs by
    0.394 A between two preparations; the prepared ligand PDBQTs come from a
    Meeko version whose successor changes atom order and TORSDOF; and
    ``exhaustiveness_arm_status.json`` is the sole surviving provenance for the
    gnina wall clock the thesis cites. These refuse to run under any mode.
"""
from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Iterable, Sequence

ROOT = Path("/home/manndo/master_dev")
VINA_PY = Path("/home/manndo/anaconda3/envs/vina/bin/python")

BITEXACT = "bitexact"
VERIFY = "verify"
NEVER = "never"

CHEAP, MODERATE, EXPENSIVE = "cheap", "moderate", "expensive"

_COST_HINT = {
    CHEAP: "seconds to a few minutes",
    MODERATE: "minutes to about an hour",
    EXPENSIVE: "hours to days",
}


# =============================================================================
# Environment guard
# =============================================================================

# Appendix B states these exactly. The RDKit pin matters most: the appendix notes
# that its distance-geometry bounds set the bond-length and bond-angle tolerances
# and its UFF supplies the internal-energy reference ensemble, so "most
# PoseBusters verdicts are tied to this exact version".
EXPECTED_PY_PACKAGES = {
    "rdkit": "2025.09.1",
    "numpy": "2.3.4",
    "scipy": "1.16.2",
    "pandas": "2.3.3",
    "sklearn": "1.8.0",
    "posebusters": "0.6.3",
    "pandamap": "4.1.0",
}
EXPECTED_PYTHON = "3.12.12"

# Version STRINGS to look for in each binary's own output. The Vina entry is the
# locally patched build that adds a boron atom type; stock /usr/bin/vina reports
# 1.2.5 and must never be used, which is why the path is pinned rather than
# resolved through PATH.
EXPECTED_BINARIES = {
    "vina": (Path("/home/manndo/AutoDock-Vina/build/linux/release/vina"),
             ["--version"], "v1.2.7-20-g93cdc3d-mod"),
    "gnina": (Path("/home/manndo/docking_tools/gnina"), ["--version"], "1.3.2"),
    "smina": (Path("/home/manndo/anaconda3/envs/diffdock/bin/smina"),
              ["--version"], "Nov  9 2017"),
    "obabel": (Path(shutil.which("obabel") or "/usr/bin/obabel"), ["-V"], "3.1.1"),
    "prepare_receptor": (Path("/home/manndo/ADFRsuite-1.0/bin/prepare_receptor"), [], None),
}

# Both live outside the repository and both are described in Appendix B.
POSEBUSTERS_PATCH = Path("/home/manndo/anaconda3/envs/vina/lib/python3.12"
                         "/site-packages/posebusters/modules/energy_ratio.py")
DIFFDOCK_PATCH = Path("/home/manndo/docking_tools/DiffDock/inference.py")


@dataclass
class Check:
    name: str
    ok: bool
    found: str
    expected: str
    note: str = ""


def _run_text(cmd: Sequence[str]) -> str:
    try:
        p = subprocess.run([str(c) for c in cmd], capture_output=True, text=True, timeout=60)
        return (p.stdout or "") + (p.stderr or "")
    except Exception as exc:  # pragma: no cover - environment probe
        return f"<error: {exc}>"


def check_environment() -> list[Check]:
    """Assert the environment matches what Appendix B reports.

    Returns a list of checks rather than raising, so the notebook can print the
    whole table and the reader can see which single item drifted.
    """
    out: list[Check] = []

    found_py = ".".join(map(str, sys.version_info[:3]))
    out.append(Check("python", found_py == EXPECTED_PYTHON, found_py, EXPECTED_PYTHON))

    for mod, want in EXPECTED_PY_PACKAGES.items():
        try:
            m = __import__(mod)
            got = getattr(m, "__version__", "?")
        except Exception as exc:
            got = f"<import failed: {exc}>"
        # RDKit ships 2025.09.1 under the string 2025.09.1 but pip records
        # 2025.9.1; accept either spelling of the same release.
        ok = got == want or got.replace(".0", ".") == want.replace(".0", ".")
        out.append(Check(mod, ok, str(got), want))

    for name, (path, args, want) in EXPECTED_BINARIES.items():
        if not path.exists():
            out.append(Check(name, False, "<not found>", want or "present",
                             note=f"expected at {path}"))
            continue
        if want is None:
            out.append(Check(name, True, str(path), "present"))
            continue
        text = _run_text([path, *args])
        out.append(Check(name, want in text, text.strip().splitlines()[0] if text.strip() else "",
                         want, note=str(path)))

    # The PoseBusters patch. As shipped, the handler guarding the UFF parameter
    # check indexes e.args[1] on an assertion carrying a single argument, raising
    # IndexError and aborting the whole bust() call. Every pose of any ligand UFF
    # cannot parameterise would be dropped with no error reported.
    pb_ok = False
    pb_note = "energy_ratio.py not found"
    if POSEBUSTERS_PATCH.exists():
        src = POSEBUSTERS_PATCH.read_text()
        block = re.search(r"UFFHasAllMoleculeParams.*?except Exception as e:\s*\n(.*?)return _empty_results",
                          src, re.S)
        pb_ok = bool(block) and "e.args[1]" not in (block.group(1) if block else "")
        pb_note = "UFF handler logs the exception" if pb_ok else "UFF handler still indexes e.args[1]"
    out.append(Check("posebusters energy_ratio patch", pb_ok, pb_note, "patched"))

    # The DiffDock seed hook. Upstream DiffDock exposes no seed option at all.
    dd_ok = DIFFDOCK_PATCH.exists() and "DIFFDOCK_SEED" in DIFFDOCK_PATCH.read_text()
    out.append(Check("diffdock seed patch", dd_ok,
                     "DIFFDOCK_SEED hook present" if dd_ok else "absent", "patched",
                     note=str(DIFFDOCK_PATCH)))
    return out


def diffdock_seed_evidence() -> dict[str, int]:
    """Count stored DiffDock runs that actually fired the seed hook.

    This is the evidence behind the VERIFY classification, and it is worth
    recomputing rather than trusting: it is the difference between a run that can
    be reproduced and one that can only be checked.
    """
    trees = {
        "benchmark": ROOT / "Dockings/Benchmark_DiffDock",
        "orai_control": ROOT / "Dockings/Orai_Benchmark_DiffDock",
        "orai_experimental": ROOT / "Dockings/diffdock_results",
    }
    counts = {}
    for label, tree in trees.items():
        if not tree.exists():
            counts[label] = -1
            continue
        p = subprocess.run(["grep", "-rl", r"\[repro\] seeded", str(tree)],
                           capture_output=True, text=True)
        counts[label] = len([ln for ln in p.stdout.splitlines() if ln.strip()])
    return counts


# =============================================================================
# Canonical trees
# =============================================================================

# The trap this table exists to prevent. There are 57 directories under
# posebusters_results/ and 22 under pandamap_results/, and the ones whose names
# read like scratch are the CURRENT ones. Generation order, oldest first:
# *_PRE_FR0 / *_PRE_EXH128 / *_UFFON_backup_* -> plain -> _matched. Verified by
# checksum, not by naming convention: the shipped thesis Figure 9 asset
# (media/media/image43.png, md5 80b9a1b9...) matches the _orai_matched_root copy
# and not the plain one.
CANONICAL = {
    "benchmark_pb": ROOT / "posebusters_results/benchmark_matched_equibind/dock",
    "benchmark_report": ROOT / "posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report",
    "benchmark_clusters": ROOT / "posebusters_results/cluster_crystal_pocket_matched_equibind",
    "benchmark_pandamap": ROOT / "pandamap_results/benchmark_matched_equibind",
    "benchmark_effort_charged": ROOT / "posebusters_results/benchmark_matched_equibind/docking_effort_gnina_v2_charged",
    "benchmark_effort_elapsed": ROOT / "posebusters_results/benchmark_matched_equibind/docking_effort_gnina_v2",
    "exhaustiveness": ROOT / "posebusters_results/autodock_exhaustiveness_returns",
    "orai_root": ROOT / "posebusters_results/_orai_matched_root",
    "orai_control": ROOT / "posebusters_results/_orai_matched_root/orai_benchmark",
    "orai_experimental": ROOT / "posebusters_results/_orai_matched_root/orai_jku",
    "orai_yield": ROOT / "posebusters_results/_orai_matched_root/orai_pbvalid_yield_compare",
    "orai_tm_share": ROOT / "posebusters_results/_orai_matched_root/orai_pbvalid_tm_share_compare",
    "orai_pandamap_control": ROOT / "pandamap_results/orai_benchmark_matched",
    "orai_pandamap_experimental": ROOT / "pandamap_results/orai_jku_matched",
    "orai_pandamap_compare": ROOT / "pandamap_results/orai_interaction_compare_matched",
    "dataset_figures": ROOT / "PoseBusters_Benchmark_Analysis/figures",
}

# Superseded trees, kept so the notebook can say plainly why it is NOT reading
# them. Each maps to the canonical key that replaced it.
SUPERSEDED = {
    "posebusters_results/benchmark_full_protein_vina_scoring": "benchmark_pb",
    "posebusters_results/benchmark": "benchmark_pb",
    "posebusters_results/cluster_crystal_pocket_full_protein": "benchmark_clusters",
    "pandamap_results/benchmark_full_protein_mgltools": "benchmark_pandamap",
    "posebusters_results/orai_benchmark": "orai_control",
    "posebusters_results/orai_jku": "orai_experimental",
    "posebusters_results/orai_pbvalid_yield_compare": "orai_yield",
    "posebusters_results/orai_pbvalid_tm_share_compare": "orai_tm_share",
    "pandamap_results/orai_interaction_compare": "orai_pandamap_compare",
}


# =============================================================================
# Stage registry
# =============================================================================

@dataclass
class Stage:
    name: str
    title: str
    section: str
    determinism: str
    cost: str
    outputs: list[str]                       # paths or globs, relative to ROOT
    thesis: str = ""                         # what it supports
    needs: list[str] = field(default_factory=list)
    cmd: list[str] | None = None             # argv, or None when run is a callable
    run: Callable[[], None] | None = None
    notes: str = ""
    min_matches: int = 1                     # for glob outputs
    env: dict[str, str] = field(default_factory=dict)

    def resolve(self) -> list[Path]:
        found: list[Path] = []
        for pat in self.outputs:
            if any(ch in pat for ch in "*?["):
                found.extend(sorted(ROOT.glob(pat)))
            else:
                p = ROOT / pat
                if p.exists():
                    found.append(p)
        return found

    def satisfied(self) -> tuple[bool, str]:
        missing = []
        for pat in self.outputs:
            if any(ch in pat for ch in "*?["):
                n = len(list(ROOT.glob(pat)))
                if n < self.min_matches:
                    missing.append(f"{pat} ({n} < {self.min_matches})")
            elif not (ROOT / pat).exists():
                missing.append(pat)
        return (not missing), ("; ".join(missing) if missing else "all outputs present")


class Registry:
    """Ordered stage collection with dependency checking."""

    def __init__(self) -> None:
        self.stages: dict[str, Stage] = {}
        self._order: list[str] = []

    def add(self, stage: Stage) -> Stage:
        if stage.name in self.stages:
            raise ValueError(f"duplicate stage {stage.name!r}")
        for dep in stage.needs:
            if dep not in self.stages:
                raise ValueError(
                    f"stage {stage.name!r} needs {dep!r}, which is not registered yet. "
                    "Stages must be declared in dependency order; this is the check "
                    "that prevents the old notebook's read-before-write bug.")
        self.stages[stage.name] = stage
        self._order.append(stage.name)
        return stage

    def __len__(self) -> int:
        return len(self.stages)

    def topo(self) -> list[Stage]:
        """Declaration order, verified to be a valid topological order."""
        seen: set[str] = set()
        for name in self._order:
            st = self.stages[name]
            unmet = [d for d in st.needs if d not in seen]
            if unmet:
                raise ValueError(f"{name} declared before its dependencies {unmet}")
            seen.add(name)
        return [self.stages[n] for n in self._order]

    def by_section(self) -> dict[str, list[Stage]]:
        out: dict[str, list[Stage]] = {}
        for st in self.topo():
            out.setdefault(st.section, []).append(st)
        return out


def md5(path: Path, limit: int = 64 * 1024 * 1024) -> str:
    h = hashlib.md5()
    with path.open("rb") as fh:
        read = 0
        while chunk := fh.read(1 << 20):
            h.update(chunk)
            read += len(chunk)
            if read >= limit:
                break
    return h.hexdigest()


class Driver:
    """Walk the registry under one run mode.

    Modes
    -----
    ``verify``       check outputs exist; run nothing. The default, because the
                     full pipeline is days of GPU and CPU time and every stage's
                     output is already on disk.
    ``run-missing``  additionally execute stages whose outputs are absent.
    ``force``        additionally re-execute the stages named in ``only``.
                     NEVER-class stages still refuse.
    """

    def __init__(self, registry: Registry, mode: str = "verify",
                 only: Iterable[str] | None = None, dry_run: bool = False) -> None:
        if mode not in {"verify", "run-missing", "force"}:
            raise ValueError(f"unknown mode {mode!r}")
        if mode == "force" and not only:
            raise ValueError("force mode requires an explicit stage list in `only`; "
                             "forcing everything would re-run days of docking and "
                             "would destroy the unseeded DiffDock trees")
        self.reg = registry
        self.mode = mode
        self.only = set(only or [])
        self.dry_run = dry_run
        self.results: list[dict] = []

    def _should_run(self, st: Stage, satisfied: bool) -> tuple[bool, str]:
        # The two refusals are checked BEFORE the force branch, deliberately.
        # Putting force first would let `force` bypass exactly the protections
        # that exist for the irreversible cases, which is the one mistake this
        # class is here to make impossible.
        if st.determinism == NEVER:
            return False, ("REFUSED: not reproducible, and a re-run destroys the "
                           "provenance of the reported result")
        if st.determinism == VERIFY:
            return False, ("REFUSED: the reported run was an unseeded single draw. "
                           "Re-running draws a different sample and would not "
                           "reproduce the thesis. To redraw deliberately, invoke the "
                           "driver script by hand; the registry will not do it for you.")
        if self.mode == "force" and st.name in self.only:
            return True, "forced"
        if self.mode in {"run-missing", "force"} and not satisfied:
            return True, "outputs missing"
        return False, "skipped"

    def execute(self, st: Stage) -> tuple[bool, str]:
        if self.dry_run:
            return True, "dry-run"
        started = time.time()
        if st.run is not None:
            st.run()
        elif st.cmd:
            env = dict(os.environ)
            env.update(st.env)
            proc = subprocess.run([str(c) for c in st.cmd], cwd=ROOT, env=env)
            if proc.returncode != 0:
                return False, f"exit {proc.returncode}"
        else:
            return False, "no command or callable declared"
        return True, f"ok in {time.time() - started:.1f}s"

    def walk(self) -> list[dict]:
        self.results = []
        for st in self.reg.topo():
            satisfied, detail = st.satisfied()
            unmet = [d for d in st.needs
                     if not self.reg.stages[d].satisfied()[0]]
            run_it, why = self._should_run(st, satisfied)
            if run_it and unmet:
                run_it, why = False, f"BLOCKED: inputs missing from {unmet}"
            status = "present" if satisfied else "MISSING"
            action = why
            if run_it:
                ok, msg = self.execute(st)
                action = f"ran ({msg})" if ok else f"FAILED ({msg})"
                satisfied, detail = st.satisfied()
                status = "present" if satisfied else "MISSING"
            self.results.append(dict(stage=st.name, section=st.section, title=st.title,
                                     determinism=st.determinism, cost=st.cost,
                                     status=status, action=action, detail=detail))
        return self.results

    def table(self) -> str:
        rows = [("STAGE", "STATUS", "DET", "ACTION")]
        rows += [(r["stage"], r["status"], r["determinism"][:4], r["action"][:52])
                 for r in self.results]
        w = [max(len(r[i]) for r in rows) for i in range(4)]
        out = []
        for i, r in enumerate(rows):
            out.append("  ".join(v.ljust(w[j]) for j, v in enumerate(r)))
            if i == 0:
                out.append("  ".join("-" * w[j] for j in range(4)))
        n_missing = sum(1 for r in self.results if r["status"] == "MISSING")
        out.append("")
        out.append(f"{len(self.results)} stages, {n_missing} with missing outputs")
        return "\n".join(out)


def write_manifest(path: Path, registry: Registry, checks: Sequence[Check],
                   results: Sequence[dict]) -> Path:
    """Record what this run saw, so a later reader can tell drift from difference."""
    payload = {
        "written": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "root": str(ROOT),
        "python": ".".join(map(str, sys.version_info[:3])),
        "environment": [vars(c) for c in checks],
        "diffdock_seed_evidence": diffdock_seed_evidence(),
        "canonical_trees": {k: str(v) for k, v in CANONICAL.items()},
        "stages": [
            dict(name=s.name, section=s.section, determinism=s.determinism,
                 cost=s.cost, thesis=s.thesis, outputs=s.outputs, needs=s.needs)
            for s in registry.topo()
        ],
        "results": list(results),
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2))
    return path
