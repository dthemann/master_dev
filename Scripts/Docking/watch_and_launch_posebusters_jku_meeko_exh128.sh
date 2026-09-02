#!/usr/bin/env bash
# Wait for the Orai x JKU ADFRsuite exh128 driver (docking + inline gnina) to finish,
# validate its gnina layer, launch PoseBusters, and then SUPERVISE that run.
#
# The supervision half is not optional. PoseBusters deadlocks during pool drain on this
# box and its own stall watchdog does NOT fire — observed twice (2026-08-29 at
# 8,990/8,991, dead 6h24m; 2026-08-30 at 12,310/12,315, dead 36m). Both logs advertise
# "restart the pool if no pose completes for 1800s"; neither acted. Because
# overwrite:false makes runs resumable, the fix is to kill and relaunch: it picks up
# from the last CSV checkpoint.
set -uo pipefail

ROOT=/home/manndo/master_dev
S=/tmp/claude-1000/-home-manndo-master-dev/291e74c9-2cdf-49c5-8d5e-2c8f77e80cdc/scratchpad
DOCKLOG="$S/orai_jku_meeko_exh128_run.log"
PBLOG="$S/orai_jku_meeko_exh128_posebusters.log"
WLOG="$S/jku_meeko_chain.log"
STATE="$S/jku_meeko_chain.state"
CFG="$ROOT/Scripts/Docking/Posebusters/posebusters_orai_jku_meeko_exh128_config.yaml"
OPT="$ROOT/Dockings/Orai_JKU_Meeko_exh128/mgl_tools/optimization_log.csv"
PY=/home/manndo/anaconda3/envs/vina/bin/python
STALL=1800        # seconds of log silence before we declare the pool dead
MAX_RESTART=4

cd "$ROOT" || exit 1
say(){ echo "[$(date '+%F %T')] $*" | tee -a "$WLOG"; }

say "chain armed; waiting for the JKU docking+gnina driver"
echo "WAITING_FOR_DOCKING" > "$STATE"
while pgrep -f 'python.*run_orai_mgltools_arm\.py.*orai_jku_autodock_meeko' >/dev/null 2>&1; do sleep 30; done
say "driver exited"

if ! grep -q 'Done\. Results:' "$DOCKLOG" 2>/dev/null; then
    say "ABORT: driver did not finish cleanly"; tail -25 "$DOCKLOG" | tee -a "$WLOG"
    echo "ABORTED_DOCKING_INCOMPLETE" > "$STATE"; exit 1
fi

say "validating the gnina layer"
"$PY" - "$OPT" <<'PYEOF' | tee -a "$WLOG"
import csv, collections, sys
rows = list(csv.DictReader(open(sys.argv[1])))
conv = collections.Counter(r.get("converter","?") for r in rows)
st   = collections.Counter(r.get("status","?")    for r in rows)
rk   = collections.Counter(r.get("rank_metric","?") for r in rows)
print(f"  rows={len(rows)} converter={dict(conv)} status={dict(st)} rank={dict(rk)}")
bad=[]
# POLARITY: this is the MEEKO-ligand twin. Meeko poses carry REMARK SMILES, so mk_export is the
# CORRECT and expected converter here — it is what the published JKU_meeko arm used. The failure
# mode to catch is the mirror image: rdkit_template_map would mean ADFRsuite poses were staged
# into this arm by mistake, silently collapsing the pair into an arm compared with itself.
if conv.get("rdkit_template_map"):
    bad.append(f"{conv['rdkit_template_map']} rdkit_template_map rows — ADFRsuite poses in the MEEKO arm")
if conv.get("obabel"):    bad.append(f"{conv['obabel']} obabel rows — chemistry unreliable")
if set(rk)-{"cnn_affinity"}: bad.append(f"unexpected rank_metric {set(rk)}")
if not rows: bad.append("empty optimiser log")
for b in bad: print(f"  FAIL: {b}")
sys.exit(1 if bad else 0)
PYEOF
if [ "${PIPESTATUS[0]}" -ne 0 ]; then
    say "ABORT: gnina layer failed validation; NOT busting."; echo "ABORTED_GNINA_INVALID" > "$STATE"; exit 1
fi

PB_PID=""
# SCOPE TO THIS ARM. A bare `pgrep/pkill -f run_posebusters.py` matches the OTHER JKU
# arm's run too — and the notebook arms both chains concurrently — so an unscoped pkill
# would kill the sibling's PoseBusters, and an unscoped pgrep would report the sibling
# as "this arm is alive". Track our own PID instead.
launch(){ setsid nohup "$PY" -u "$ROOT/Scripts/Docking/Posebusters/run_posebusters.py" -c "$CFG" >> "$PBLOG" 2>&1 &
         PB_PID=$!; }

say "launching PoseBusters"
echo "POSEBUSTERS_RUNNING" > "$STATE"
launch; sleep 30

restarts=0
while true; do
    sleep 60
    if ! kill -0 "$PB_PID" 2>/dev/null; then
        if grep -qE 'PASSED ALL [0-9]+ TESTS' "$PBLOG" 2>/dev/null; then
            say "PoseBusters COMPLETE"; grep -E 'PASSED ALL [0-9]+ TESTS|PoseBusters Results' "$PBLOG" | tail -2 | tee -a "$WLOG"
            echo "COMPLETED" > "$STATE"; exit 0
        fi
        say "process gone without a completion line — treating as failure"
        tail -20 "$PBLOG" | tee -a "$WLOG"; echo "FAILED" > "$STATE"; exit 1
    fi
    silent=$(( $(date +%s) - $(stat -c %Y "$PBLOG") ))
    if [ "$silent" -gt "$STALL" ]; then
        if [ "$restarts" -ge "$MAX_RESTART" ]; then
            say "stalled again after $restarts restarts; giving up"; echo "STALLED_GIVEUP" > "$STATE"; exit 1
        fi
        restarts=$((restarts+1))
        say "STALL: ${silent}s of silence (internal watchdog did not fire). Killing and resuming (restart $restarts/$MAX_RESTART)"
        kill -TERM "$PB_PID" 2>/dev/null; sleep 10; kill -9 "$PB_PID" 2>/dev/null; sleep 5
        launch; sleep 30
    fi
done
