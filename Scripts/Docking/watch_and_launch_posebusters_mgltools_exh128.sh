#!/usr/bin/env bash
# Wait for the Orai ADFRsuite exh128 docking driver (docking + inline gnina) to finish,
# validate that its gnina layer is sound, then launch PoseBusters on the new arm.
#
# The launch is GATED on validation, not merely on the process exiting. A gnina layer
# that fell back to mk_export or Open Babel would mean the wrong poses were optimised,
# and busting them would produce a plausible-looking but wrong result.
set -uo pipefail

ROOT=/home/manndo/master_dev
SCRATCH=/tmp/claude-1000/-home-manndo-master-dev/291e74c9-2cdf-49c5-8d5e-2c8f77e80cdc/scratchpad
DOCKLOG="$SCRATCH/orai_mgltools_exh128_run.log"
PBLOG="$SCRATCH/orai_mgltools_exh128_posebusters.log"
WLOG="$SCRATCH/posebusters_chain.log"
STATE="$SCRATCH/posebusters_chain.state"
CFG="$ROOT/Scripts/Docking/Posebusters/posebusters_orai_benchmark_mgltools_exh128_config.yaml"
OPTLOG="$ROOT/Dockings/Orai_Benchmark_MGLTools_exh128/mgl_tools/optimization_log.csv"
PY=/home/manndo/anaconda3/envs/vina/bin/python

cd "$ROOT" || exit 1
say() { echo "[$(date '+%F %T')] $*" | tee -a "$WLOG"; }

say "chain armed; waiting for the docking+gnina driver to finish"
echo "WAITING_FOR_DOCKING" > "$STATE"

# ── 1. wait for the driver to exit ───────────────────────────────────────────
while pgrep -f 'python.*run_orai_mgltools_arm\.py' >/dev/null 2>&1; do
    sleep 120
done
say "driver has exited"

# ── 2. did it finish cleanly? ────────────────────────────────────────────────
if ! grep -q 'Done\. Results:' "$DOCKLOG" 2>/dev/null; then
    say "ABORT: driver exited without 'Done. Results:' — it failed or was killed."
    tail -25 "$DOCKLOG" | tee -a "$WLOG"
    echo "ABORTED_DOCKING_INCOMPLETE" > "$STATE"
    exit 1
fi

# ── 3. validate the gnina layer BEFORE busting anything ─────────────────────
say "validating the gnina layer"
"$PY" - <<'PYEOF' | tee -a "$WLOG"
import csv, collections, sys
p = ("/home/manndo/master_dev/Dockings/Orai_Benchmark_MGLTools_exh128/"
     "mgl_tools/optimization_log.csv")
try:
    rows = list(csv.DictReader(open(p)))
except OSError as e:
    print(f"FAIL: cannot read optimization_log.csv: {e}"); sys.exit(1)
conv = collections.Counter(r.get("converter", "?") for r in rows)
st   = collections.Counter(r.get("status", "?") for r in rows)
rk   = collections.Counter(r.get("rank_metric", "?") for r in rows)
print(f"  rows={len(rows)} converter={dict(conv)} status={dict(st)} rank={dict(rk)}")
bad = []
if conv.get("mk_export"):
    bad.append(f"{conv['mk_export']} rows used mk_export — MEEKO poses were optimised")
if conv.get("obabel"):
    bad.append(f"{conv['obabel']} rows fell back to Open Babel — chemistry unreliable")
if set(rk) - {"cnn_affinity"}:
    bad.append(f"unexpected rank_metric {set(rk)}")
if len(rows) < 12000:
    bad.append(f"only {len(rows)} optimiser rows; expected ~12,320")
if bad:
    for b in bad:
        print(f"  FAIL: {b}")
    sys.exit(1)
print("  gnina layer OK")
PYEOF
if [ "${PIPESTATUS[0]}" -ne 0 ]; then
    say "ABORT: gnina layer failed validation; NOT launching PoseBusters."
    echo "ABORTED_GNINA_INVALID" > "$STATE"
    exit 1
fi

# ── 4. don't start on a saturated box ────────────────────────────────────────
for _ in $(seq 1 60); do
    load=$(awk '{print int($1)}' /proc/loadavg)
    [ "$load" -lt 20 ] && break
    say "load $load too high; waiting"
    sleep 60
done

# ── 5. launch ────────────────────────────────────────────────────────────────
say "launching PoseBusters on the new arm"
echo "LAUNCHING_POSEBUSTERS" > "$STATE"
setsid nohup "$PY" -u "$ROOT/Scripts/Docking/Posebusters/run_posebusters.py" \
    -c "$CFG" > "$PBLOG" 2>&1 &
PB=$!
sleep 30
if kill -0 "$PB" 2>/dev/null; then
    say "PoseBusters running, PID $PB -> $PBLOG"
    echo "POSEBUSTERS_RUNNING $PB" > "$STATE"
else
    say "PoseBusters exited within 30 s:"
    tail -25 "$PBLOG" | tee -a "$WLOG"
    echo "POSEBUSTERS_LAUNCH_FAILED" > "$STATE"
    exit 1
fi
