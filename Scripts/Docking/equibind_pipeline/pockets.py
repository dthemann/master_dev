"""PocketInfo dataclass plus fpocket / p2rank result parsers."""
from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

from .config import CFG


@dataclass
class PocketInfo:
    source: str
    pocket_id: int
    unique_id: str
    score: float
    center: Tuple[float, float, float]
    radius: float
    protein_name: str
    extra: dict = field(default_factory=dict)

    def to_dict(self) -> Dict:
        return {
            "source": self.source, "pocket_id": self.pocket_id,
            "unique_id": self.unique_id, "score": self.score,
            "center": list(self.center), "radius": self.radius,
            "protein_name": self.protein_name, "extra": self.extra,
        }


def _parse_fpocket_pocket_centroid(pocket_pdb: Path) -> Tuple[Tuple[float, float, float], float]:
    coords: List[Tuple[float, float, float]] = []
    with open(pocket_pdb, "r") as fh:
        for line in fh:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
                except (ValueError, IndexError):
                    continue
    if not coords:
        return (0.0, 0.0, 0.0), 0.0
    arr = np.asarray(coords, dtype=np.float64)
    centroid = arr.mean(axis=0)
    radius = float(np.linalg.norm(arr - centroid, axis=1).max())
    return (float(centroid[0]), float(centroid[1]), float(centroid[2])), radius


def parse_fpocket_results(fpocket_dir: Path, protein_name: str,
                          n_top: int | None = None) -> List[PocketInfo]:
    if n_top is None:
        n_top = CFG.n_top_pockets
    pockets: List[PocketInfo] = []
    for out_dir in sorted(fpocket_dir.glob(f"{protein_name}*_out")):
        pockets_dir = out_dir / "pockets"
        if not pockets_dir.exists():
            continue
        dir_suffix = out_dir.name[len(protein_name):]
        variant = dir_suffix.replace("_out", "").strip("_") or "raw"

        scores: Dict[int, float] = {}
        info_files = list(out_dir.glob("*_info.txt"))
        if info_files:
            with open(info_files[0]) as fh:
                current_pocket = None
                for line in fh:
                    s = line.strip()
                    if s.startswith("Pocket ") and ":" in s:
                        try:
                            current_pocket = int(s.split()[1])
                        except (ValueError, IndexError):
                            current_pocket = None
                    elif current_pocket is not None and s.startswith("Score"):
                        try:
                            scores[current_pocket] = float(s.split(":")[-1].strip())
                        except ValueError:
                            pass

        for ppdb in sorted(pockets_dir.glob("pocket*_atm.pdb")):
            try:
                pocket_num = int(ppdb.stem.replace("pocket", "").replace("_atm", ""))
            except ValueError:
                continue
            centroid, radius = _parse_fpocket_pocket_centroid(ppdb)
            uid = f"fpocket_{variant}_p{pocket_num:03d}"
            pockets.append(PocketInfo(
                source="fpocket", pocket_id=pocket_num, unique_id=uid,
                score=scores.get(pocket_num, 0.0), center=centroid, radius=radius,
                protein_name=protein_name,
                extra={"pdb_file": str(ppdb), "out_dir": str(out_dir), "variant": variant},
            ))
    pockets.sort(key=lambda p: p.score, reverse=True)
    return pockets[:n_top]


def parse_p2rank_results(p2rank_dir: Path, protein_name: str,
                         n_top: int | None = None) -> List[PocketInfo]:
    if n_top is None:
        n_top = CFG.n_top_pockets
    pockets: List[PocketInfo] = []
    for pred_file in sorted(p2rank_dir.glob(f"{protein_name}*_predictions.csv")):
        base = pred_file.stem.replace("_predictions", "")
        suffix_part = base[len(protein_name):].replace(".pdb", "").strip("_")
        variant = suffix_part or "raw"
        with open(pred_file) as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                # CSV may have whitespace-padded keys
                cleaned = {k.strip(): (v.strip() if isinstance(v, str) else v) for k, v in row.items()}
                try:
                    rank = int(cleaned["rank"])
                    score = float(cleaned["score"])
                    prob = float(cleaned["probability"])
                    cx = float(cleaned["center_x"])
                    cy = float(cleaned["center_y"])
                    cz = float(cleaned["center_z"])
                    name = cleaned.get("name", "")
                    sas = float(cleaned.get("sas_points") or 0)
                except (TypeError, ValueError, KeyError):
                    continue
                radius = max(5.0, np.sqrt(sas) * 0.5)
                uid = f"p2rank_{variant}_p{rank:03d}"
                pockets.append(PocketInfo(
                    source="p2rank", pocket_id=rank, unique_id=uid,
                    score=score, center=(cx, cy, cz), radius=radius,
                    protein_name=protein_name,
                    extra={"probability": prob, "sas_points": sas,
                           "pred_file": str(pred_file), "name": name, "variant": variant},
                ))
    pockets.sort(key=lambda p: p.score, reverse=True)
    return pockets[:n_top]
