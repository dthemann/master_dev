#!/usr/bin/env python3
"""EquiBind 3-phase pose-generation pipeline — entry point."""
from __future__ import annotations

import sys
from pathlib import Path

_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

from equibind_pipeline.config import CFG               # noqa: E402
from equibind_pipeline.gpu_check import validate_gpu   # noqa: E402
from equibind_pipeline.pipeline import run_pipeline    # noqa: E402


def main() -> int:
    CFG.validate()
    CFG.ensure_output_dirs()
    CFG.banner()
    validate_gpu()
    summary = run_pipeline()
    print("\n\u2713 Pipeline complete.")
    print(f"  Poses (success): {summary['totals']['poses_success']}")
    print(f"  Poses (failed):  {summary['totals']['poses_failed']}")
    print(f"  Wall time:       {summary['global_timing']['pipeline_wall_time_s']:.1f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
