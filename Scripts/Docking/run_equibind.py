#!/usr/bin/env python3
"""EquiBind 3-phase pose-generation pipeline — entry point.

NOTE: This script is configured by ``EQ_*`` environment variables only (see
``equibind_pipeline/config.py``); it does NOT read
``equibind_docking_config.yaml`` directly. That YAML is consumed by
``Master_Docking.ipynb``, which maps each key to an ``EQ_*`` var before
launching this script. Editing the YAML therefore has no effect on a bare
``python run_equibind.py`` run — set the ``EQ_*`` vars (or launch via the
notebook) instead.

All settings can be supplied as ``EQ_*`` environment variables (see
``equibind_pipeline/config.py``); the CLI flags below override them for the two
knobs added for the clash study:

  --clamp-mode {on,off,both}    centroid clamping toward the pocket centre.
                                'both' emits __clampON / __clampOFF variants.
  --refine-mode {on,off,both}   docking/force-field re-search of each pose.
                                'both' emits __refRAW / __ref<TOOL> variants.
  --refine-tool {smina,gnina,both}
                                re-search backend. smina = physics; gnina =
                                physics + CNN rescore (writes CNNscore/affinity);
                                'both' emits one refined variant per tool
                                (__refSMINA + __refGNINA).

Examples:
    # Default behaviour (clamp on, no re-search) — unchanged from before.
    python run_equibind.py

    # Measure the effect of clamping AND of an smina re-search in one run
    # (4 output poses per EquiBind pose: clamp{ON,OFF} x ref{RAW,SMINA}).
    python run_equibind.py --clamp-mode both --refine-mode both

    # Disable clamping, always re-search with gnina (CNN-scored).
    python run_equibind.py --clamp-mode off --refine-mode on --refine-tool gnina

    # Emit raw + smina + gnina refined poses for a head-to-head comparison.
    python run_equibind.py --refine-mode both --refine-tool both
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

from equibind_pipeline.config import CFG               # noqa: E402
from equibind_pipeline.gpu_check import validate_gpu   # noqa: E402
from equibind_pipeline.pipeline import run_pipeline     # noqa: E402


def _parse_args(argv=None) -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="EquiBind 3-phase pose-generation pipeline.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    ap.add_argument("--clamp-mode", choices=("on", "off", "both"), default=None,
                    help="Centroid clamping toward the pocket centre "
                         "(default: EQ_CLAMP_MODE or 'on').")
    ap.add_argument("--refine-mode", choices=("on", "off", "both"), default=None,
                    help="Docking/force-field re-search of each pose "
                         "(default: EQ_REFINE_MODE or 'off').")
    ap.add_argument("--refine-tool", choices=("smina", "gnina", "both"), default=None,
                    help="Re-search backend: smina | gnina | both "
                         "(default: EQ_REFINE_TOOL or 'smina').")
    ap.add_argument("--smina-search", choices=("local_only", "minimize"), default=None,
                    help="smina/gnina local-optimisation mode (default: 'local_only').")
    ap.add_argument("--smina-path", default=None,
                    help="Path to the smina executable (default: EQ_SMINA / PATH).")
    ap.add_argument("--gnina-path", default=None,
                    help="Path to the gnina executable (default: EQ_GNINA / PATH).")
    ap.add_argument("--gnina-gpu", dest="gnina_gpu", action="store_const", const=True,
                    default=None,
                    help="Run gnina's CNN on the GPU (default: off / --no_gpu, since "
                         "the GPU is usually busy with EquiBind inference).")
    return ap.parse_args(argv)


def _apply_cli_overrides(args: argparse.Namespace) -> None:
    """CLI flags take precedence over the EQ_* env vars already baked into CFG."""
    if args.clamp_mode is not None:
        CFG.clamp_mode = args.clamp_mode
    if args.refine_mode is not None:
        CFG.refine_mode = args.refine_mode
    if args.refine_tool is not None:
        CFG.refine_tool = args.refine_tool
    if args.smina_search is not None:
        CFG.smina_search = args.smina_search
    if args.smina_path is not None:
        CFG.smina_executable = args.smina_path
    if args.gnina_path is not None:
        CFG.gnina_executable = args.gnina_path
    if args.gnina_gpu is not None:
        CFG.gnina_use_gpu = args.gnina_gpu


def main(argv=None) -> int:
    args = _parse_args(argv)
    _apply_cli_overrides(args)
    CFG.validate()
    CFG.ensure_output_dirs()
    CFG.banner()
    validate_gpu()
    summary = run_pipeline()
    print("\n✓ Pipeline complete.")
    print(f"  Poses (success): {summary['totals']['poses_success']}")
    print(f"  Poses (failed):  {summary['totals']['poses_failed']}")
    print(f"  Wall time:       {summary['global_timing']['pipeline_wall_time_s']:.1f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
