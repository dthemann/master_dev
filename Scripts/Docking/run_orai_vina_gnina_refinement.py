#!/usr/bin/env python3
"""Refine the completed Orai x benchmark AutoDock Vina poses with GNINA."""

from __future__ import annotations

import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from Scripts.Docking import run_orai_vinardo_gnina_refinement as shared_driver


def main() -> int:
    if any(argument.startswith("--source-scoring") for argument in sys.argv[1:]):
        raise RuntimeError(
            "the Vina-specific launcher pins --source-scoring vina; do not override it"
        )
    sys.argv[1:1] = ["--source-scoring", "vina"]
    return shared_driver.main()


if __name__ == "__main__":
    raise SystemExit(main())
