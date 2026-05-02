"""Structured signal monitor for the EquiBind docking pipeline."""
from __future__ import annotations

from enum import IntEnum
from typing import Optional


class MonitorLevel(IntEnum):
    """Verbosity levels (lower = more verbose)."""
    ALL = 1       # Full trace
    WARNING = 2   # Warnings + critical only
    CRITICAL = 3  # Fatal errors only


class DockingMonitor:
    """Tagged-message monitor; only messages at >= configured level print."""

    _PREFIXES = {
        MonitorLevel.ALL: "    \u00b7",
        MonitorLevel.WARNING: " \u26a0\ufe0f ",
        MonitorLevel.CRITICAL: " \U0001f534",
    }

    def __init__(self, level: MonitorLevel = MonitorLevel.ALL):
        self.level = level
        self.call_count = 0
        self.success_count = 0
        self.fail_count = 0
        self.rejected_count = 0
        self.in_pocket_count = 0
        self.outside_pocket_count = 0

    def _emit(self, lvl: MonitorLevel, msg: str) -> None:
        if lvl >= self.level:
            print(f"{self._PREFIXES.get(lvl, '')} {msg}")

    def info(self, msg: str) -> None:
        self._emit(MonitorLevel.ALL, msg)

    def warning(self, msg: str) -> None:
        self._emit(MonitorLevel.WARNING, msg)

    def critical(self, msg: str) -> None:
        self._emit(MonitorLevel.CRITICAL, msg)

    def equibind_call(self, protein: str, ligand: str, output_dir: str,
                      seed: int, device: str) -> None:
        self.call_count += 1
        self.info(
            f"[CALL #{self.call_count}] EquiBind \u2190 protein={protein}, "
            f"ligand={ligand}, seed={seed}, device={device}, out={output_dir}"
        )

    def equibind_return(self, success: bool, output_sdf: Optional[str],
                        error: str = "") -> None:
        if success:
            self.success_count += 1
            self.info(f"[RECV] EquiBind \u2192 OK  output={output_sdf}")
        else:
            self.fail_count += 1
            self.warning(f"[RECV] EquiBind \u2192 FAIL  error={error[:200]}")

    def pose_accepted_in_pocket(self, pose_id: str, pocket_id: str,
                                distance: float, threshold: float) -> None:
        self.in_pocket_count += 1
        self.info(
            f"[POCKET \u2713] {pose_id} INSIDE {pocket_id}  "
            f"(dist={distance:.1f}\u00c5 \u2264 {threshold:.1f}\u00c5)"
        )

    def pose_rejected_from_pocket(self, pose_id: str, pocket_id: str,
                                  distance: float, threshold: float) -> None:
        self.rejected_count += 1
        self.warning(
            f"[POCKET \u2717] {pose_id} OUTSIDE {pocket_id}  "
            f"(dist={distance:.1f}\u00c5 > {threshold:.1f}\u00c5) \u2192 REJECTED"
        )

    def pose_outside_all_pockets(self, pose_id: str, nearest_pocket: str,
                                 nearest_dist: float, threshold: float) -> None:
        self.outside_pocket_count += 1
        self.warning(
            f"[POCKET \u2717] {pose_id} not in any pocket  "
            f"(nearest={nearest_pocket}, dist={nearest_dist:.1f}\u00c5 > {threshold:.1f}\u00c5)"
        )

    def pose_inside_pocket(self, pose_id: str, pocket_id: str,
                           distance: float, threshold: float,
                           pocket_source: str) -> None:
        self.in_pocket_count += 1
        self.info(
            f"[POCKET \u2713] {pose_id} in {pocket_source} pocket {pocket_id}  "
            f"(dist={distance:.1f}\u00c5 \u2264 {threshold:.1f}\u00c5)"
        )

    def header(self, msg: str) -> None:
        bar = "\u2500" * 70
        print(f"\n{bar}\n  {msg}\n{bar}")

    def section(self, msg: str) -> None:
        self._emit(MonitorLevel.ALL, f"\u2500\u2500 {msg} \u2500\u2500")

    def print_summary(self) -> None:
        bar = "\u2550" * 70
        print(f"\n{bar}")
        print("  SIGNAL MONITOR SUMMARY")
        print(bar)
        print(f"  EquiBind calls:        {self.call_count}")
        print(f"  Successful returns:    {self.success_count}")
        print(f"  Failed returns:        {self.fail_count}")
        print(f"  Poses inside pocket:   {self.in_pocket_count}")
        print(f"  Poses rejected/outside:{self.rejected_count + self.outside_pocket_count}")
        print(f"    \u251c guided rejected:   {self.rejected_count}")
        print(f"    \u2514 unguided outside:  {self.outside_pocket_count}")
        print("\u2550" * 70)


# Singleton monitor — modules import this directly.
monitor = DockingMonitor(level=MonitorLevel.ALL)


def set_level(level: MonitorLevel) -> None:
    monitor.level = level
