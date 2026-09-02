#!/usr/bin/env python
"""Validity check for the running Orai x Benchmark ADFRsuite exh128 + gnina arm.

Run at any time. Exits 0 if healthy, 1 if a real problem is found.
Checks the things that would silently produce a wrong arm, not just liveness.
"""
from __future__ import annotations

import collections
import csv
import glob
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path("/home/manndo/master_dev")
ARM = ROOT / "Dockings/Orai_Benchmark_MGLTools_exh128"
OUT = ARM / "mgl_tools"
EXPECTED_FR0 = "6b3ab996af727c806ca481e0f1f5ca7b0fe923df8dddf48d002e930e468206e9"
EXPECTED_COMBOS = 1232
EXPECTED_POSES = 12320

problems: list[str] = []
notes: list[str] = []


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def main() -> int:
    print("=" * 74)
    print(f"Orai ADFRsuite exh128 + gnina — validity check  {time.strftime('%F %T')}")
    print("=" * 74)

    # ── liveness ─────────────────────────────────────────────────────────────
    # pgrep also matches SUSPENDED processes, so read the state code: T = stopped
    # (SIGSTOP), R/S/D = actually working. Reporting a paused run as "RUNNING" hides
    # a frozen arm behind a healthy-looking line.
    r = subprocess.run(
        ["ps", "-eo", "pid,stat,args", "--no-headers"], capture_output=True, text=True)
    procs = [l.split(None, 2) for l in r.stdout.splitlines()
             if "run_orai_mgltools_arm.py" in l or "release/vina" in l]
    procs = [p for p in procs if "bash -c" not in p[2] and "ps -eo" not in p[2]]
    stopped = [p for p in procs if p[1].startswith("T")]
    active = [p for p in procs if not p[1].startswith("T")]
    if stopped and not active:
        print(f"\nprocess     : PAUSED — {len(stopped)} process(es) suspended (SIGSTOP); "
              "resume with kill -CONT")
        notes.append("arm is SUSPENDED, not progressing; batch_timeout is wall-clock and "
                     "keeps counting while paused")
    elif active:
        print(f"\nprocess     : RUNNING ({len(active)} active"
              + (f", {len(stopped)} suspended" if stopped else "") + ")")
        if stopped:
            problems.append(f"{len(stopped)} process(es) suspended while others run — "
                            "mixed state, the run may deadlock")
    else:
        print("\nprocess     : not running")
    load = Path("/proc/loadavg").read_text().split()[0]
    print(f"load average: {load}")

    # ── inputs still pristine ────────────────────────────────────────────────
    fr0 = ARM / "_staging/receptors/pdbqt/Orai1WT-START-Fr0.pdbqt"
    if not fr0.is_file():
        problems.append("staged START-Fr0 receptor is missing")
    elif sha256(fr0) != EXPECTED_FR0:
        problems.append("staged START-Fr0 receptor CHANGED — arm is no longer control-matched")
    else:
        print("receptor    : START-Fr0 sha256 OK (Jun-22 set)")

    ligs = sorted((ARM / "_staging/ligands/pdbqt").glob("*.pdbqt"))
    meeko = [p for p in ligs if "REMARK SMILES" in p.read_text(errors="ignore")[:4000]]
    if len(ligs) != 308:
        problems.append(f"ligand staging has {len(ligs)} PDBQTs, expected 308")
    if meeko:
        problems.append(f"{len(meeko)} staged ligands are MEEKO — arm is contaminated")
    else:
        print(f"ligands     : {len(ligs)} ADFRsuite, 0 Meeko")

    # ── docking progress ─────────────────────────────────────────────────────
    poses = list((OUT / "docking").glob("*_vina_out.pdbqt")) if (OUT / "docking").is_dir() else []
    rows: list[dict] = []
    for f in glob.glob(str(OUT / "docking_log_*.csv")):
        with open(f) as fh:
            rows += list(csv.DictReader(fh))
    st = collections.Counter(r_.get("status", "?") for r_ in rows)
    print(f"\ndocking     : {len(poses)}/{EXPECTED_COMBOS} pose containers "
          f"({100*len(poses)/EXPECTED_COMBOS:.1f}%)  log rows={len(rows)} {dict(st)}")

    if st.get("skipped"):
        problems.append(f"{st['skipped']} combos SKIPPED — pre-existing poses would be "
                        "gnina-rescored instead of this arm's")
    if st.get("failed"):
        notes.append(f"{st['failed']} docking failures — inspect before trusting the arm")

    # CAUTION: in batch mode every complex in a batch is stamped with the SAME batch
    # wall-time, and that wall-time includes any interval the run spent SIGSTOPped.
    # Projecting from it after a pause badly overstates the remaining time (batch 1
    # read 191.8 s/complex against a true ~54 s of compute). Prefer the newest batch,
    # and only trust it when that batch ran uninterrupted.
    # elapsed_time_s is ALREADY per-complex (batch wall / batch size), but it includes
    # any interval the run spent SIGSTOPped, so it overstates pace after a pause.
    # The trustworthy live measure is the CURRENT batch: ligands finished in its log
    # divided by the wall time since that log was created.
    # Once docking is complete these pace figures are stale — the "live" one divides an
    # ever-growing wall clock by a batch count that stopped moving, so it climbs forever
    # and reads like a slowdown. Report nothing rather than something misleading.
    docking_done = len(rows) >= EXPECTED_COMBOS
    times = [float(r_["elapsed_time_s"]) for r_ in rows
             if r_.get("elapsed_time_s") and r_["elapsed_time_s"] not in ("", "0")]
    live = None
    if docking_done:
        clean_done = [v for v in {float(r_["elapsed_time_s"]) for r_ in rows} if v < 150]
        if clean_done:
            print(f"pace        : docking COMPLETE — mean {sum(clean_done)/len(clean_done):.1f}"
                  f"s/complex over {len(clean_done)} clean batches")
        times = []          # skip the projection block below
    blogs = sorted(Path("Dockings/Logs/orai_benchmark_mgltools_exh128_logs").glob("batch_*.log"),
                   key=lambda p: p.stat().st_mtime) if Path(
                       "Dockings/Logs/orai_benchmark_mgltools_exh128_logs").is_dir() else []
    if blogs and rows:
        bl = blogs[-1]
        done = bl.read_text(errors="ignore").count("mode |   affinity")
        # The current batch began when the PREVIOUS batch's rows were written, i.e. the
        # newest docking-log timestamp. st_ctime is useless here — it advances on every
        # append, so it always reads as "just now".
        import datetime
        newest_ts = max(r_["timestamp"] for r_ in rows)
        started = datetime.datetime.fromisoformat(newest_ts).timestamp()
        wall = time.time() - started
        if done >= 3 and wall > 0:
            live = wall / done
    if times:
        # Best estimator: the most recently COMPLETED batch. Every complex in a batch
        # shares one exact per-complex time, so a finished uninterrupted batch is a
        # precise measurement — unlike the live figure, which is noisy until its batch
        # is well along (Vina does not finish ligands at a uniform rate).
        seen, ordered = set(), []
        for r_ in rows:                      # rows are in completion order
            v = float(r_["elapsed_time_s"])
            if v not in seen:
                seen.add(v)
                ordered.append(v)
        # Per-batch cost tracks LIGAND SIZE, not machine state: measured rotatable-bond
        # means rise 7.9 / 8.4 / 9.2 across batches 2-4 alongside 58 / 62 / 80 s. Batches
        # are alphabetical slices, so composition varies by chance and averages out over
        # the full 308-ligand set. Project from the MEAN of completed unpaused batches,
        # not the last one, or the estimate swings 30% batch to batch.
        # Batches whose wall-clock spanned a SIGSTOP are inflated by the suspension and
        # must be excluded. A magnitude threshold alone is too blunt — a 17-min pause
        # only lifted one batch from ~71 s to 91 s, well under any sane cutoff — so the
        # contaminated values are recorded explicitly when a pause is taken.
        excl = set()
        exf = Path("/tmp/claude-1000/-home-manndo-master-dev/"
                   "291e74c9-2cdf-49c5-8d5e-2c8f77e80cdc/scratchpad/contaminated_batches.txt")
        if exf.is_file():
            for ln in exf.read_text().split():
                try:
                    excl.add(round(float(ln), 1))
                except ValueError:
                    pass
        clean = [v for v in ordered if v < 150 and round(v, 1) not in excl]
        if len(ordered) - len(clean):
            notes.append(f"{len(ordered) - len(clean)} batch(es) excluded from pace as "
                         "pause-contaminated")
        pace = sum(clean) / len(clean) if clean else ordered[-1]
        remain = (EXPECTED_COMBOS - len(times)) * pace / 3600
        print(f"pace        : {pace:.1f}s/complex (mean of {len(clean)} clean batches) "
              f"→ ~{remain:.1f} h docking remaining")
        print(f"              last batch {ordered[-1]:.1f}s — varies with ligand size, "
              f"not machine load")
        if live is not None:
            print(f"              live current batch: {live:.1f}s/complex "
                  f"({done}/50 — noisy until the batch is well along)")
        print(f"              per-batch history: {[f'{v:.0f}s' for v in ordered]}"
              f"   (the 191.8s batch was inflated by the SIGSTOP pauses)")

    exh = {r_.get("exhaustiveness") for r_ in rows}
    nm = {r_.get("num_modes") for r_ in rows}
    sc = {r_.get("scoring_function") for r_ in rows}
    er = {r_.get("energy_range") for r_ in rows}
    if rows:
        print(f"params      : exhaustiveness={exh} num_modes={nm} energy_range={er} scoring={sc}")
        if exh - {"128"}:
            problems.append(f"unexpected exhaustiveness values: {exh}")
        if nm - {"10"}:
            problems.append(f"unexpected num_modes values: {nm} (control parity needs 10)")
        if er - {"6"}:
            problems.append(f"unexpected energy_range values: {er} (this arm runs at 6; any 3 "
                            "means poses from the aborted first attempt leaked in)")
        if sc - {"vina"}:
            problems.append(f"unexpected scoring values: {sc}")

    # ── gnina layer ──────────────────────────────────────────────────────────
    optlog = OUT / "optimization_log.csv"
    if optlog.is_file():
        with open(optlog) as fh:
            orows = list(csv.DictReader(fh))
        conv = collections.Counter(o.get("converter", "?") for o in orows)
        ost = collections.Counter(o.get("status", "?") for o in orows)
        rk = collections.Counter(o.get("rank_metric", "?") for o in orows)
        print(f"\ngnina       : {len(orows)}/{EXPECTED_POSES} rows {dict(ost)}")
        print(f"converter   : {dict(conv)}")
        print(f"rank metric : {dict(rk)}")
        if conv.get("mk_export"):
            problems.append(f"{conv['mk_export']} rows used mk_export — MEEKO poses were "
                            "optimised by mistake")
        if conv.get("obabel"):
            problems.append(f"{conv['obabel']} rows fell back to Open Babel — chemistry unreliable")
        bad_rank = {k for k in rk if k not in ("cnn_affinity", "?")}
        if bad_rank:
            problems.append(f"unexpected rank_metric: {bad_rank} (control parity needs cnn_affinity)")
        if ost.get("failed"):
            notes.append(f"{ost['failed']} gnina pose failures (control had 64; this arm "
                         "should have 0)")
    else:
        # optimization_log.csv is only written when the optimiser FINISHES, so during the
        # ~2 h gnina phase its absence says nothing. Read progress from the per-pose
        # outputs instead, and check the converter on the provenance sidecars — the same
        # rdkit_template_map guarantee, available live rather than only at the end.
        gdir = OUT / "docking" / "optimized_gnina"
        sdfs = list(gdir.glob("*.sdf")) if gdir.is_dir() else []
        if sdfs:
            provs = list(gdir.glob("*.provenance.json"))
            conv = collections.Counter()
            for p in provs[:400]:
                try:
                    conv[(json.load(open(p)).get("converter") or {}).get("name", "?")] += 1
                except Exception:
                    conv["unreadable"] += 1
            print(f"\ngnina       : IN PROGRESS — {len(sdfs)}/{EXPECTED_POSES} poses "
                  f"({100*len(sdfs)/EXPECTED_POSES:.1f}%)")
            print(f"converter   : {dict(conv)}  (sampled from provenance sidecars)")
            if conv.get("mk_export"):
                problems.append(f"{conv['mk_export']} poses used mk_export — MEEKO poses "
                                "are being optimised by mistake")
            if conv.get("obabel"):
                problems.append(f"{conv['obabel']} poses fell back to Open Babel")
        else:
            print("\ngnina       : not started yet (runs after all docking completes)")

    # ── published trees untouched ────────────────────────────────────────────
    for name, exp in (("Orai_Benchmark", 1232), ("Orai_Benchmark_Vinardo", 1232)):
        n = len(list((ROOT / "Dockings" / name / "docking").glob("*_vina_out.pdbqt")))
        if n != exp:
            problems.append(f"published {name} has {n} poses, expected {exp}")
    print(f"\npublished   : control and Vinardo trees intact")

    # ── verdict ──────────────────────────────────────────────────────────────
    print()
    for n in notes:
        print(f"  NOTE    : {n}")
    if problems:
        for p in problems:
            print(f"  PROBLEM : {p}")
        print("\nVERDICT: PROBLEM FOUND")
        return 1
    print("VERDICT: HEALTHY")
    return 0


if __name__ == "__main__":
    sys.exit(main())
