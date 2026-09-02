#!/usr/bin/env python
"""Build the flat MGLTools-prepared ligand PDBQT directory that
``optimize_prepared_ligand_pdbqt_dir`` points at.

The MGLTools whole-protein arm stages its prepared ligands *per complex*
(``<output_dir>/<PDB_CCD>/_staging/ligands/pdbqt/``), but
``run_autodock._resolve_mgl_pose_template`` resolves one flat
``<root>/<ligand_name>.pdbqt``. Without this directory the config's
``optimize_require_template_reconstruction`` preflight fails from a clean
checkout, so this script is the committed way to (re)create it.

MGLTools ligand preparation is deterministic, so the files written here are
byte-identical to the per-complex copies Vina docks; ``--verify`` asserts that
for every complex already staged.

    python Scripts/Docking/build_mgltools_ligand_template_root.py
    python Scripts/Docking/build_mgltools_ligand_template_root.py --verify
"""
from __future__ import annotations

import argparse
import filecmp
import sys
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from Scripts.Docking.run_autodock import get_pdbqt_dir  # noqa: E402
from Scripts.Utilities.prep_docking import run_workflow  # noqa: E402

DEFAULT_CONFIG = (
    REPO_ROOT
    / "Scripts/Docking/autodock_vina_docking_config_full_protein_search_vina_mgltools.yaml"
)
BENCHMARK_DIR = REPO_ROOT / "Data/PoseBuster Benchmark Set"
IDS_FILE = BENCHMARK_DIR / "posebusters_pdb_ccd_ids.txt"
EXPECTED_COUNT = 308


def authoritative_ids() -> list[str]:
    ids = sorted({
        line.strip() for line in IDS_FILE.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    })
    if len(ids) != EXPECTED_COUNT:
        raise SystemExit(f"expected {EXPECTED_COUNT} authoritative IDs, found {len(ids)}")
    return ids


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-c", "--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument(
        "--verify", action="store_true",
        help="also assert byte-identity against every per-complex staged PDBQT",
    )
    args = parser.parse_args()

    cfg = yaml.safe_load(args.config.read_text())
    prepared_root = Path(cfg["optimize_prepared_ligand_pdbqt_dir"])
    if not prepared_root.is_absolute():
        prepared_root = REPO_ROOT / prepared_root
    stage = prepared_root.parent          # <output_dir>/_staging/ligands
    stage.mkdir(parents=True, exist_ok=True)

    ids = authoritative_ids()
    linked = 0
    for pdb_id in ids:
        src = BENCHMARK_DIR / pdb_id / f"{pdb_id}_ligand_start_conf.sdf"
        if not src.is_file():
            raise SystemExit(f"missing authoritative ligand SDF: {src}")
        link = stage / src.name
        if not link.exists():
            link.symlink_to(src.resolve())
            linked += 1
    print(f"staged {len(ids)} ligand SDFs into {stage} ({linked} new symlink(s))")

    # Same call the notebook's Variant-A cell makes per complex, so the flat
    # copies and the per-complex staged copies come from one code path.
    run_workflow(
        input_dir=stage, contains="ligands", output_dir=get_pdbqt_dir(stage),
        process_postfixes=False, convert_ligands_with_meeko=False,
        converter="mgltools", convert_ligands=True, verbose=False,
    )
    written = sorted(prepared_root.glob("*.pdbqt"))
    print(f"prepared {len(written)} ligand PDBQT(s) in {prepared_root}")
    if len(written) != EXPECTED_COUNT:
        raise SystemExit(f"expected {EXPECTED_COUNT} prepared PDBQTs, found {len(written)}")

    if args.verify:
        output_base = REPO_ROOT / cfg["output_dir"]
        checked = mismatched = 0
        for pdb_id in ids:
            staged = (
                output_base / pdb_id / "_staging/ligands/pdbqt"
                / f"{pdb_id}_ligand_start_conf.pdbqt"
            )
            if not staged.is_file():
                continue
            checked += 1
            if not filecmp.cmp(staged, prepared_root / staged.name, shallow=False):
                mismatched += 1
                print(f"  MISMATCH: {staged}")
        print(f"verified {checked} per-complex staged PDBQT(s); {mismatched} mismatch(es)")
        if mismatched:
            return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
