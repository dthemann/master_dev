#!/usr/bin/env bash
# Wait for the running PoseBusters job to finish, then launch the Orai x Benchmark
# ADFRsuite-receptor + ADFRsuite-ligand exh128 + gnina docking arm.
#
# The launch is chained here rather than driven by an external poller so the docking
# starts within ~60 s of the machine freeing up, regardless of anything else.
set -uo pipefail

TARGET_PID=${TARGET_PID:-3281521}
ROOT=/home/manndo/master_dev
SCRATCH=/tmp/claude-1000/-home-manndo-master-dev/291e74c9-2cdf-49c5-8d5e-2c8f77e80cdc/scratchpad
STATE="$SCRATCH/orai_arm_watcher.state"
DOCKLOG="$SCRATCH/orai_mgltools_exh128_run.log"
WATCHLOG="$SCRATCH/orai_arm_watcher.log"
CFG="$ROOT/Scripts/Docking/orai_benchmark_autodock_mgltools_exh128_gnina_config.yaml"
PY=/home/manndo/anaconda3/envs/vina/bin/python

cd "$ROOT" || exit 1
mkdir -p "$SCRATCH"

say() { echo "[$(date '+%F %T')] $*" | tee -a "$WATCHLOG"; }

say "watcher started; waiting on PID $TARGET_PID (PoseBusters)"
echo "WAITING" > "$STATE"

# ── Phase 1: wait for the target to exit ─────────────────────────────────────
# Guard against PID reuse by requiring the cmdline to still look like the job.
while kill -0 "$TARGET_PID" 2>/dev/null; do
    if ! tr '\0' ' ' < "/proc/$TARGET_PID/cmdline" 2>/dev/null | grep -q 'run_posebusters.py'; then
        say "PID $TARGET_PID no longer looks like run_posebusters (PID reuse); proceeding"
        break
    fi
    sleep 60
done
say "PID $TARGET_PID has exited"

# ── Phase 2: let stragglers drain ────────────────────────────────────────────
# run_posebusters forks ~50 workers; they can outlive the parent briefly.
for _ in $(seq 1 30); do
    n=$(pgrep -fc run_posebusters.py 2>/dev/null || echo 0)
    [ "$n" -eq 0 ] && break
    say "waiting for $n PoseBusters worker(s) to drain"
    sleep 30
done

# ── Phase 3: refuse to launch on top of another heavy AutoDock/DiffDock job ──
if pgrep -f 'run_autodock|run_diffdock_arm|run_unidock' >/dev/null 2>&1; then
    say "ABORT: another docking job is running; not launching. Start manually when clear."
    echo "ABORTED_BUSY" > "$STATE"
    exit 1
fi

load=$(awk '{print $1}' /proc/loadavg)
say "load average now $load; launching docking"

# ── Phase 4: launch ──────────────────────────────────────────────────────────
echo "LAUNCHING" > "$STATE"
setsid nohup "$PY" -u "$ROOT/Scripts/Docking/run_orai_mgltools_arm.py" -c "$CFG" \
    > "$DOCKLOG" 2>&1 &
DOCK_PID=$!
sleep 20

if kill -0 "$DOCK_PID" 2>/dev/null; then
    say "docking launched, PID $DOCK_PID -> $DOCKLOG"
    echo "RUNNING $DOCK_PID" > "$STATE"
else
    say "LAUNCH FAILED — process exited within 20 s. Tail of log:"
    tail -20 "$DOCKLOG" | tee -a "$WATCHLOG"
    echo "LAUNCH_FAILED" > "$STATE"
    exit 1
fi

# Surface the guard block so a guard rejection is visible immediately.
grep -E '\[ok\]|REFUSING TO RUN' "$DOCKLOG" 2>/dev/null | tee -a "$WATCHLOG"
say "watcher done; docking is running"
