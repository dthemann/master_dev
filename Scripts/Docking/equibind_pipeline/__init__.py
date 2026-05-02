"""EquiBind batch docking pipeline (3-phase CPU/GPU).

Public entry point:
    from equibind_pipeline.pipeline import run_pipeline
    run_pipeline()  # uses defaults from equibind_pipeline.config
"""
from .config import CFG  # re-export for convenience

__all__ = ["CFG"]
