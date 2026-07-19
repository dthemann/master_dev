#!/usr/bin/env python3
"""Define the Orai1 transmembrane (TM1 pore) exclusion layer and filter out
docked poses that land inside it.

Biology
-------
Orai1 is a hexameric Ca2+ channel; its pore is lined by the TM1 helix of each of
the 6 subunits (chains A-F here). Going along the pore axis, TM1 runs from the
cytosolic basic region (Arg91, "R91") to the extracellular selectivity filter
(Glu106, "E106").  NOTE: residue 106 in these structures is GLUTAMATE (E106), not
glutamine -- glutamine at that end of TM1 is Q108. The membrane-spanning band
between the R91 ring and the E106 ring is the transmembrane region: it is where
the channel conducts ions, not a druggable surface pocket. Poses that dock there
are therefore not biologically meaningful binding sites and are excluded.

Geometry ("the layer")
----------------------
Per receptor frame we take the six R91 Calpha atoms and the six E106 Calpha atoms
(one per subunit). Their centroids define the two ends of the pore axis:

    p_R91  = mean of the 6 Arg91  CA          (t = 0        , axial_frac = 0)
    p_E106 = mean of the 6 Glu106 CA          (t = L        , axial_frac = 1)
    axis   = unit(p_E106 - p_R91)             (~ the pore / membrane normal)

For any point x, t(x) = (x - p_R91) . axis  is its position along the axis and
radial(x) = |(x - p_R91) - t*axis| its distance from the axis. The exclusion
LAYER is the axial slab  -pad <= t <= L + pad  (a slab perpendicular to the pore
axis, i.e. the membrane-spanning band). r_pore, the mean radial distance of the
12 ring CA atoms, describes the pore/backbone radius and is used only for the
optional "in_pore_lumen" sub-flag; the primary transmembrane test is slab
membership at any radius (a ligand on the lipid-facing TM surface is still in the
membrane).

Because the slab is radially unbounded, drawing it at r_pore would misleadingly
suggest only the conducting channel is excluded. So the layer also records r_outer:
the 99th-percentile radial distance of all protein heavy atoms whose axial position
falls inside the slab -- i.e. how far the hexamer's OUTER SHELL reaches from the axis
in this membrane band (~35 A vs ~8 A for the pore). The visualization (PDB cage /
PyMOL CGO) draws the exclusion boundary out to r_outer, with the pore cylinder kept
as an inner reference and explicit R91 / E106 boundary planes, so the picture matches
the actual test: the whole cross-section of the membrane-spanning band is excluded.

A pose is flagged in_transmembrane when its ligand heavy-atom CENTROID falls in
the slab. frac_atoms_in_layer and majority_in_layer are also reported so a
stricter atom-level rule can be applied without recomputing.

Because the 4 Orai frames (START-Fr0 + MD snapshots Fr300/400/499) are different
conformations in different coordinate systems, the layer is rebuilt per frame from
that frame's own receptor file (the exact `protein_file_used` the poses were
docked into), so poses and layer always share a coordinate system.

Inputs / outputs
----------------
Pose source: the Orai PoseBusters per-pose CSV (`posebusters_filtered_results.csv`),
which lists every generated pose (AutoDock Vina, DiffDock, EquiBind) with its
`pose_file` (per-pose SDF in receptor coords) and `protein_file_used` receptor.

Per dataset out-dir (default posebusters_results/<ds>/transmembrane_filter/):
  tm_exclusion_zones.json     -- per-frame layer parameters (axis, ends, L, r_pore, r_outer)
  <frame>_tm_zone.pdb         -- pseudoatom cage of the band: outer shell + pore walls,
                                 centreline, R91/E106 end-cap rings (load w/ receptor)
  <frame>_tm_zone.pml         -- PyMOL: receptor + translucent outer-shell exclusion
                                 cylinder (r_outer), inner pore reference, R91/E106 planes
  tm_pose_classification.csv  -- EVERY pose + geometry + flags
  tm_excluded_poses.csv       -- poses inside the TM layer (the filtered-out list)
  tm_kept_poses.csv           -- poses that survive the filter
  tm_invalid_poses.csv        -- poses dropped as not placed on the protein
  summary.json                -- counts by method / frame

  fig_tool_tm_loss.png        -- per tool: of all produced poses, how many PB-valid
                                 ones are lost to the TM vs survive (+ tool_tm_loss_summary.csv)
  fig_pbvalid_lost_matrix.png -- Orai frame x ligand matrix of PB-valid poses lost to
                                 the TM (+ pbvalid_lost_matrix.csv)
  fig_rank_survival.png       -- pose rank vs number of PB-valid poses OUTSIDE the TM
                                 (the potentially-valuable poses; + rank_survival.csv)

Toolchains / optimizer variants
-------------------------------
DiffDock and EquiBind each carry raw + optimizer variants (DiffDock raw/smina/gnina;
EquiBind smina/gnina). The per-pose CSVs keep EVERY variant, each tagged with its
`variant` and a display `toolchain` (an asterisk marks a gnina-optimised toolchain:
DiffDock*, EquiBind*). The FIGURES / TABLE / summary report ONE variant per tool,
chosen by --diffdock-variant / --equibind-variant (default gnina), so DiffDock is
not double-counted. EquiBind has no native score, so its poses are ranked by gnina
energy (`gnina_affinity`, most-negative = rank 1) for the rank-survival figure.

Run in the analysis conda env (`vina`):
    conda run -n vina python Scripts/Analysis/orai_transmembrane_exclusion.py
    # gnina toolchains are the default; restrict ligands / pick variants:
    conda run -n vina python Scripts/Analysis/orai_transmembrane_exclusion.py \
        --dataset orai_jku --ligands gsk7975a Synta-66 \
        --diffdock-variant gnina --equibind-variant gnina
"""
from __future__ import annotations

import argparse
import json
import math
import re
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

# Shared, unit-tested statistical helpers (Holm/BH, Wilson CI, G-test of
# independence, Cochran-Armitage trend, ...). Robust path insert so the module
# imports whether the script is run from the repo root or its own directory.
try:
    import os as _os, sys as _sys
    _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
    import stats_utils as su
    _HAVE_SU = True
except Exception:                                   # pragma: no cover
    su = None
    _HAVE_SU = False

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    _HAVE_MPL = True
except Exception:                                   # pragma: no cover
    _HAVE_MPL = False

# Per-tool identity colours (fixed order, shared with orai_pose_cluster_report.py).
TOOL_STYLE = {
    "autodock": ("AutoDock Vina", "#4C72B0", "o"),
    "diffdock": ("DiffDock", "#55A868", "^"),
    "equibind": ("EquiBind", "#C44E52", "s"),
}
# Fate of a produced pose (semantic outcome palette; hatch on the "lost" segment
# so the survive/lost distinction is not carried by colour alone — CVD-safe).
_FATE = [
    ("survive",  "PB-valid · outside TM (kept)",       "#55A868", ""),
    ("lost_tm",  "PB-valid · inside TM (lost)",        "#C44E52", "///"),
    ("unplaced", "PB-valid · off-protein/unevaluable", "#DD8452", ""),
    ("nonvalid", "not PB-valid",                       "#CFCFCF", ""),
]


def _frame_sort_key(protein: str) -> int:
    """MD-time ordinal for an Orai frame label (START-Fr0 < MDSnap-Fr300 < Fr400
    < Fr499). Used to order the frames for the Cochran-Armitage trend test."""
    m = re.search(r"Fr(\d+)", str(protein))
    return int(m.group(1)) if m else 10 ** 9


def _tool_key(method: str) -> str:
    m = str(method).lower()
    for k in ("autodock", "diffdock", "equibind"):
        if m.startswith(k):
            return k
    return m


def _pose_variant(method, optimizer, refine_variant, pose_name) -> str:
    """The optimizer/refinement variant of a pose: 'gnina' | 'smina' | 'raw' (and
    'none' for AutoDock, which has no post-processing). DiffDock optimisation is
    post-hoc (`optimizer` col / optimized_<tool>/ path); EquiBind refinement is
    inline (`refine_variant` col)."""
    base = _tool_key(method)
    if base == "autodock":
        return "none"
    if base == "diffdock":
        v = str(optimizer).strip().lower()
        if v in ("gnina", "smina"):
            return v
        n = str(pose_name)
        if "optimized_gnina" in n or n.endswith("_gnina.sdf"):
            return "gnina"
        if "optimized_smina" in n or n.endswith("_smina.sdf"):
            return "smina"
        return "raw"                          # 'original' / NaN
    if base == "equibind":
        v = str(refine_variant).strip().lower()
        return v if v in ("gnina", "smina") else "raw"
    return "none"


def _toolchain_label(base: str, variant: str) -> str:
    """Display label; an asterisk marks a gnina-optimised toolchain (DiffDock*,
    EquiBind*). Raw = the plain tool name; smina named explicitly."""
    name = TOOL_STYLE.get(base, (str(base).title(),))[0]
    if base == "autodock" or variant in ("none", ""):
        return name
    if variant == "gnina":
        return f"{name}*"
    if variant == "raw":
        return name
    return f"{name} ({variant})"


def _variant_matches(variant: str, selector: str) -> bool:
    """Does a pose variant satisfy a --<tool>-variant selector ('all' = any)."""
    sel = str(selector).lower()
    if sel == "all":
        return True
    if sel in ("raw", "original"):
        return variant == "raw"
    return variant == sel

# Reuse the pipeline's single source of truth for the PoseBusters validity verdict
# (the canonical 20-test allowlist + its bool coercion) so pb_valid here agrees
# EXACTLY with run_posebusters.py / the validity report. The "any boolean-like
# column" heuristic is deliberately NOT used: it sweeps in diagnostic columns such
# as `most_extreme_clash_protein` (always False) and mis-coerces object-dtype
# bool/NaN columns, silently flipping real failures to passes.
_PROJECT_ROOT = Path(__file__).resolve().parents[2]
_PB_DIR = _PROJECT_ROOT / "Scripts" / "Docking" / "Posebusters"
for _p in (str(_PROJECT_ROOT), str(_PB_DIR)):
    if _p not in sys.path:
        sys.path.insert(0, _p)
try:
    from run_posebusters import CANONICAL_TEST_COLUMNS, coerce_test_cols_to_bool
    _HAVE_CANON = True
except Exception:                                   # pragma: no cover
    CANONICAL_TEST_COLUMNS, coerce_test_cols_to_bool = (), None
    _HAVE_CANON = False

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
R91_RESNUM, R91_RESNAME = 91, "ARG"
E106_RESNUM, E106_RESNAME = 106, "GLU"
RING_ATOM = "CA"

# metadata columns carried into the output CSVs when present
_CARRY = [
    "docking_method", "protein", "ligand", "pose_name", "pose_file",
    "autodock_rank", "autodock_affinity", "diffdock_confidence",
    "smina_affinity", "gnina_affinity",
    "refine_variant", "optimizer", "pocket_source", "protein_file_used",
    "protein-ligand_maximum_distance",
]

# Physical-placement gate. A pose only counts as a real binding pose if the ligand
# actually sits on the protein. DiffDock occasionally emits rigid-body-translated
# poses hundreds-to-thousands of Å out of the pocket (internally intact but not
# placed — e.g. the JKU 2abp/boron worst-rank poses). geom_ok (|coord|<=1e4) is too
# loose to catch these, so they would otherwise fall into tm_kept_poses.csv as
# "outside the TM" with absurd axial_frac. We drop them via PoseBusters' own
# protein-ligand_maximum_distance verdict (authoritative), backed by a geometric
# envelope in case that column is absent. Legit poses sit within radial ~58 Å and
# axial_frac ~[-2.2, 2.9]; the corrupt ones run to radial 12000 Å / axial_frac 400.
PLACEMENT_FLAG = "protein-ligand_maximum_distance"
RADIAL_MAX = 80.0                    # Å from the pore axis (receptor envelope + margin)
AXIAL_FRAC_LO, AXIAL_FRAC_HI = -4.0, 5.0   # centroid position along R91→E106 axis

DATASETS = {
    "orai_benchmark": "posebusters_results/orai_benchmark/dock/posebusters_filtered_results.csv",
    "orai_jku": "posebusters_results/orai_jku/dock/posebusters_filtered_results.csv",
}


# ---------------------------------------------------------------------------
# Receptor parsing / zone construction
# ---------------------------------------------------------------------------
def _ring_ca(pdb_path: Path, resnum: int, resname: str, atom: str = RING_ATOM) -> np.ndarray:
    """All `atom` coords for residue (resname, resnum) across every chain."""
    pts = []
    for ln in pdb_path.read_text(errors="ignore").splitlines():
        if not ln.startswith(("ATOM", "HETATM")):
            continue
        if ln[12:16].strip() != atom:
            continue
        if ln[17:20].strip() != resname:
            continue
        try:
            if int(ln[22:26]) != resnum:
                continue
            pts.append([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
        except ValueError:
            continue
    return np.asarray(pts, dtype=float)


def _pca_axis(pdb_path: Path) -> Optional[np.ndarray]:
    """Principal (largest-variance) axis of all CA atoms -- a sanity reference for
    the R91->E106 axis (they should be near-parallel for a channel)."""
    xs = []
    for ln in pdb_path.read_text(errors="ignore").splitlines():
        if ln.startswith(("ATOM", "HETATM")) and ln[12:16].strip() == "CA":
            try:
                xs.append([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
            except ValueError:
                continue
    if len(xs) < 10:
        return None
    a = np.asarray(xs)
    _, vecs = np.linalg.eigh(np.cov((a - a.mean(0)).T))
    return vecs[:, -1]


def _protein_heavy_xyz(pdb_path: Path) -> np.ndarray:
    """All protein heavy-atom coordinates (ATOM records, hydrogens excluded)."""
    pts = []
    for ln in pdb_path.read_text(errors="ignore").splitlines():
        if not ln.startswith("ATOM"):
            continue
        el = ln[76:78].strip() or ln[12:16].strip().lstrip("0123456789")[:1]
        if el.upper() == "H":
            continue
        try:
            pts.append([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
        except ValueError:
            continue
    return np.asarray(pts, dtype=float)


def build_zone(pdb_path: Path, pad: float) -> Dict:
    """Construct the TM exclusion layer for one receptor frame."""
    r91 = _ring_ca(pdb_path, R91_RESNUM, R91_RESNAME)
    e106 = _ring_ca(pdb_path, E106_RESNUM, E106_RESNAME)
    if len(r91) == 0 or len(e106) == 0:
        raise ValueError(f"{pdb_path}: found {len(r91)} R91 and {len(e106)} E106 CA atoms")
    p_lo = r91.mean(0)                      # R91 end  (axial_frac 0)
    p_hi = e106.mean(0)                     # E106 end (axial_frac 1)
    vec = p_hi - p_lo
    length = float(np.linalg.norm(vec))
    axis = vec / length
    # radial pore radius from the 12 ring CA atoms
    ring = np.vstack([r91, e106])
    proj = (ring - p_lo) @ axis
    radial = np.linalg.norm((ring - p_lo) - np.outer(proj, axis), axis=1)
    r_pore = float(radial.mean())
    # OUTER SHELL radius of the transmembrane band. The exclusion test is slab
    # membership at ANY radius (a ligand on the lipid-facing TM surface is still in
    # the membrane), so the pore radius r_pore describes only the conducting channel,
    # NOT the true radial extent of the excluded band. For an honest visualization we
    # measure how far the hexamer reaches out from the axis WITHIN this axial slab:
    # the radial distance of every protein heavy atom whose axial position falls in
    # [t_lo, t_hi]. r_outer uses the 99th percentile so a few flexible loop tips that
    # happen to sit axially in-band (but are not part of the compact TM bundle) do not
    # balloon the drawn boundary; r_outer_max keeps the absolute extent for reference.
    r_outer = r_outer_max = None
    n_heavy_in_slab = 0
    heavy = _protein_heavy_xyz(pdb_path)
    if len(heavy):
        dh = heavy - p_lo
        th = dh @ axis
        rh = np.linalg.norm(dh - np.outer(th, axis), axis=1)
        in_slab = (th >= -pad) & (th <= length + pad)
        if in_slab.any():
            rs = rh[in_slab]
            n_heavy_in_slab = int(in_slab.sum())
            r_outer = round(float(np.percentile(rs, 99)), 3)
            r_outer_max = round(float(rs.max()), 3)
    # sanity: angle between the residue-defined axis and the CA-PCA axis
    pca = _pca_axis(pdb_path)
    axis_pca_deg = None
    if pca is not None:
        c = abs(float(np.clip(np.dot(axis, pca / np.linalg.norm(pca)), -1, 1)))
        axis_pca_deg = round(math.degrees(math.acos(c)), 2)
    return {
        "receptor_file": str(pdb_path),
        "p_R91": p_lo.tolist(), "p_E106": p_hi.tolist(),
        "axis": axis.tolist(), "length": length, "pad": pad,
        "t_lo": -pad, "t_hi": length + pad,
        "r_pore_mean": r_pore, "r_pore_max": float(radial.max()),
        "r_outer": r_outer, "r_outer_max": r_outer_max,
        "n_heavy_in_slab": n_heavy_in_slab,
        "n_R91_subunits": int(len(r91)), "n_E106_subunits": int(len(e106)),
        "axis_vs_pca_deg": axis_pca_deg,
        "residues": {"lower": f"{R91_RESNAME}{R91_RESNUM}", "upper": f"{E106_RESNAME}{E106_RESNUM}"},
    }


def classify_points(zone: Dict, xyz: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Return (t, radial) along/from the pore axis for an (n,3) array of points."""
    p_lo = np.asarray(zone["p_R91"]); axis = np.asarray(zone["axis"])
    d = xyz - p_lo
    t = d @ axis
    radial = np.linalg.norm(d - np.outer(t, axis), axis=1)
    return t, radial


# ---------------------------------------------------------------------------
# Ligand pose parsing (V2000 SDF; heavy atoms only, no RDKit needed)
# ---------------------------------------------------------------------------
def _read_v3000(lines: List[str]) -> List[list]:
    out, in_atom = [], False
    for ln in lines:
        s = ln.strip()
        if s.startswith("M  V30 BEGIN ATOM"):
            in_atom = True; continue
        if s.startswith("M  V30 END ATOM"):
            break
        if in_atom and s.startswith("M  V30"):
            parts = s.split("V30", 1)[1].split()          # idx element x y z aamap ...
            if len(parts) >= 5:
                el = parts[1]
                if el not in ("H", "D"):
                    out.append([float(parts[2]), float(parts[3]), float(parts[4])])
    return out


def _read_heavy_xyz(path: str) -> Optional[np.ndarray]:
    """Heavy-atom coordinates from a single-molecule SDF (V2000 or V3000)."""
    try:
        lines = Path(path).read_text(errors="ignore").splitlines()
        if len(lines) < 4:
            return None
        if "V3000" in lines[3]:
            out = _read_v3000(lines)
        else:
            natoms = int(lines[3][:3])
            out = []
            for ln in lines[4:4 + natoms]:
                el = ln[31:34].strip()
                if el in ("H", "D"):
                    continue
                out.append([float(ln[0:10]), float(ln[10:20]), float(ln[20:30])])
        if not out:
            return None
        return np.asarray(out, dtype=float)
    except (ValueError, IndexError, OSError):
        return None


# worker globals (avoid re-sending zones per task)
_ZONES: Dict[str, Dict] = {}


def _init_worker(zones: Dict[str, Dict]):
    global _ZONES
    _ZONES = zones


def _process_pose(task: Tuple[int, str, str]) -> Dict:
    idx, frame, pose_file = task
    zone = _ZONES.get(frame)
    res = {"_idx": idx, "n_heavy_atoms": 0, "cx": None, "cy": None, "cz": None,
           "axial_t": None, "axial_frac": None, "radial_dist": None,
           "frac_atoms_in_layer": None, "in_transmembrane": None,
           "in_pore_lumen": None, "majority_in_layer": None, "geom_ok": False}
    if zone is None:
        res["geom_error"] = "no_zone_for_frame"
        return res
    xyz = _read_heavy_xyz(pose_file)
    if xyz is None or not np.all(np.isfinite(xyz)) or np.max(np.abs(xyz)) > 1e4:
        res["geom_error"] = "unreadable_or_corrupt"
        return res
    cent = xyz.mean(0)
    t_c, rad_c = classify_points(zone, cent[None, :])
    t_c, rad_c = float(t_c[0]), float(rad_c[0])
    t_atoms, _ = classify_points(zone, xyz)
    in_layer_atoms = (t_atoms >= zone["t_lo"]) & (t_atoms <= zone["t_hi"])
    frac = float(in_layer_atoms.mean())
    in_tm = bool(zone["t_lo"] <= t_c <= zone["t_hi"])
    res.update({
        "n_heavy_atoms": int(len(xyz)),
        "cx": round(float(cent[0]), 3), "cy": round(float(cent[1]), 3),
        "cz": round(float(cent[2]), 3),
        "axial_t": round(t_c, 3),
        "axial_frac": round(t_c / zone["length"], 4),
        "radial_dist": round(rad_c, 3),
        "frac_atoms_in_layer": round(frac, 4),
        "in_transmembrane": in_tm,
        "in_pore_lumen": bool(in_tm and rad_c <= zone["r_pore_mean"]),
        "majority_in_layer": bool(frac >= 0.5),
        "geom_ok": True,
    })
    return res


# ---------------------------------------------------------------------------
# pb_valid derivation
# ---------------------------------------------------------------------------
def _derive_pb_valid(df: pd.DataFrame) -> pd.Series:
    """Canonical PoseBusters validity: passes every canonical test present.
    Uses the exact allowlist + coercion from run_posebusters.py."""
    if not _HAVE_CANON:
        return pd.Series(np.nan, index=df.index)
    checks = [c for c in CANONICAL_TEST_COLUMNS if c in df.columns]
    if not checks:
        return pd.Series(np.nan, index=df.index)
    sub = df[checks].copy()
    coerce_test_cols_to_bool(sub, checks)
    return sub.all(axis=1)


# ---------------------------------------------------------------------------
# Visualization writers
# ---------------------------------------------------------------------------
def _perp_basis(axis: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    tmp = np.array([1.0, 0.0, 0.0]) if abs(axis[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    u = np.cross(axis, tmp); u /= np.linalg.norm(u)
    v = np.cross(axis, u)
    return u, v


def write_zone_pdb(zone: Dict, out: Path, n_ring: int = 48, n_layers: int = 9):
    """Pseudoatom cage tracing the TM exclusion band. The band is an axial slab
    between the R91 and E106 rings, unbounded in radius; here we draw it out to the
    hexamer's OUTER SHELL so the visualization reflects that the exclusion covers the
    whole membrane-spanning cross-section (the lipid-facing surface), not just the
    conducting pore. Elements/chains separate the parts for easy colouring in PyMOL:
      chain O (resn OUT, carbon)   -- outer shell wall  (r_outer)
      chain P (resn POR, oxygen)   -- inner pore wall    (r_pore, reference)
      chain Y (resn AXS, nitrogen) -- pore centreline
      chain B (resn CAP, sulphur)  -- end-cap rings at the R91 and E106 boundary planes
    """
    p_lo = np.asarray(zone["p_R91"]); axis = np.asarray(zone["axis"])
    L = zone["length"]
    r_pore = float(zone["r_pore_mean"])
    r_out = float(zone.get("r_outer") or r_pore)
    u, v = _perp_basis(axis)
    recs = []
    serial = 0

    def _wall(radius: float, el: str, resn: str, ch: str):
        nonlocal serial
        for k in range(n_layers):
            frac = k / (n_layers - 1)
            center = p_lo + frac * L * axis
            for j in range(n_ring):
                th = 2 * math.pi * j / n_ring
                p = center + radius * (math.cos(th) * u + math.sin(th) * v)
                serial += 1
                recs.append((serial, el, resn, ch, serial, p, frac * 100))

    _wall(r_out, "C", "OUT", "O")                 # outer shell wall (the exclusion boundary)
    _wall(r_pore, "O", "POR", "P")                # inner pore wall (reference)
    # pore centreline
    for k in range(n_layers):
        frac = k / (n_layers - 1)
        p = p_lo + frac * L * axis
        serial += 1
        recs.append((serial, "N", "AXS", "Y", serial, p, frac * 100))
    # end-cap rings: concentric rings from the pore out to the shell at the R91 (frac 0)
    # and E106 (frac 1) planes, so the two axial boundaries of the band read as disks
    n_cap_rings = 4
    for frac in (0.0, 1.0):
        center = p_lo + frac * L * axis
        for m in range(n_cap_rings + 1):
            radius = r_pore + (r_out - r_pore) * m / n_cap_rings
            for j in range(n_ring):
                th = 2 * math.pi * j / n_ring
                p = center + radius * (math.cos(th) * u + math.sin(th) * v)
                serial += 1
                recs.append((serial, "S", "CAP", "B", serial, p, frac * 100))
    lines = ["REMARK  Orai1 transmembrane exclusion band (R91 ring <-> E106 ring)",
             f"REMARK  axis {np.round(axis,4).tolist()}  length {L:.2f}  "
             f"r_pore {r_pore:.2f}  r_outer {r_out:.2f}"]
    for serial, el, resn, ch, resi, p, b in recs:
        lines.append(
            f"HETATM{serial % 100000:>5d} {el:>2s}   {resn:>3s} {ch}{resi % 10000:>4d}    "
            f"{p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}  1.00{b:6.2f}          {el:>2s}"
        )
    lines.append("END")
    out.write_text("\n".join(lines) + "\n")


def write_zone_pml(zone: Dict, frame: str, out: Path):
    p_lo = np.asarray(zone["p_R91"]); p_hi = np.asarray(zone["p_E106"])
    axis = np.asarray(zone["axis"])
    r_pore = float(zone["r_pore_mean"])
    r_out = float(zone.get("r_outer") or r_pore)
    receptor = zone["receptor_file"]
    zone_pdb = out.with_suffix(".pdb").name
    # thin end-cap disks (short, fat CGO cylinders) marking the R91 and E106 planes,
    # each spanning the full band cross-section (radius r_outer) so the two axial
    # boundaries of the excluded band are explicit, colour-matched to the residues.
    cap = 0.5
    r91_a = p_lo - 0.5 * cap * axis;  r91_b = p_lo + 0.5 * cap * axis
    e106_a = p_hi - 0.5 * cap * axis; e106_b = p_hi + 0.5 * cap * axis
    # The whole scene is wrapped in a `python ... python end` block so PyMOL runs it as
    # genuine Python: the multi-line CGO list literals with injected coordinates parse
    # reliably (the .pml command interpreter mishandles multi-line list continuation).
    txt = f"""# PyMOL: Orai1 transmembrane exclusion band for {frame}
# usage:  pymol {out.name}
# The excluded band is the axial slab between the R91 and E106 CA rings, unbounded in
# radius. It is drawn out to the hexamer OUTER SHELL (r_outer={r_out:.1f} A) so it is
# clear the exclusion covers the whole membrane-spanning cross-section, not just the
# pore (r_pore={r_pore:.1f} A, drawn as an inner reference).
python
from pymol import cmd
from pymol.cgo import CYLINDER

cmd.load(r"{receptor}", "receptor")
cmd.hide("everything", "receptor")
cmd.show("cartoon", "receptor")
cmd.color("grey80", "receptor")
cmd.set("cartoon_transparency", 0.15, "receptor")
cmd.select("R91", "receptor and resi 91 and resn ARG")
cmd.select("E106", "receptor and resi 106 and resn GLU")
cmd.show("sticks", "R91 or E106")
cmd.color("blue", "R91")
cmd.color("red", "E106")

# OUTER shell exclusion boundary -- spans the whole hexamer cross-section, R91 -> E106
tm_outer = [CYLINDER, {p_lo[0]:.3f}, {p_lo[1]:.3f}, {p_lo[2]:.3f}, {p_hi[0]:.3f}, {p_hi[1]:.3f}, {p_hi[2]:.3f}, {r_out:.3f}, 1.0, 0.55, 0.0, 1.0, 0.55, 0.0]
cmd.load_cgo(tm_outer, "TM_outer_shell")
cmd.set("cgo_transparency", 0.70, "TM_outer_shell")

# INNER pore / conducting channel (reference only)
tm_pore = [CYLINDER, {p_lo[0]:.3f}, {p_lo[1]:.3f}, {p_lo[2]:.3f}, {p_hi[0]:.3f}, {p_hi[1]:.3f}, {p_hi[2]:.3f}, {r_pore:.3f}, 0.60, 0.0, 0.0, 0.60, 0.0, 0.0]
cmd.load_cgo(tm_pore, "TM_pore_channel")
cmd.set("cgo_transparency", 0.35, "TM_pore_channel")

# R91 (blue) and E106 (red) boundary planes -- the axial limits of the excluded band
r91_plane = [CYLINDER, {r91_a[0]:.3f}, {r91_a[1]:.3f}, {r91_a[2]:.3f}, {r91_b[0]:.3f}, {r91_b[1]:.3f}, {r91_b[2]:.3f}, {r_out:.3f}, 0.10, 0.30, 1.0, 0.10, 0.30, 1.0]
cmd.load_cgo(r91_plane, "R91_boundary_plane")
cmd.set("cgo_transparency", 0.80, "R91_boundary_plane")
e106_plane = [CYLINDER, {e106_a[0]:.3f}, {e106_a[1]:.3f}, {e106_a[2]:.3f}, {e106_b[0]:.3f}, {e106_b[1]:.3f}, {e106_b[2]:.3f}, {r_out:.3f}, 1.0, 0.10, 0.10, 1.0, 0.10, 0.10]
cmd.load_cgo(e106_plane, "E106_boundary_plane")
cmd.set("cgo_transparency", 0.80, "E106_boundary_plane")

cmd.zoom("receptor")
cmd.set("two_sided_lighting", 1)
# For a pseudoatom cage instead of solid CGO, load {zone_pdb}:
#   cmd.load(r"{zone_pdb}", "tm_cage")   # chains O=outer P=pore Y=axis B=end-caps
python end
"""
    out.write_text(txt)


# ---------------------------------------------------------------------------
# Pose rank + figures
# ---------------------------------------------------------------------------
def _pose_rank(method, pose_name, pose_file, autodock_rank=None) -> Optional[int]:
    """Native rank of a pose (1 = top-scored). AutoDock = Vina mode / affinity rank,
    DiffDock = confidence rank, EquiBind = unguided sample index. None if unknown."""
    m = _tool_key(method)
    pf = str(pose_file)
    name = Path(pf).name if pf and pf != "nan" else str(pose_name)
    if m == "autodock":
        if autodock_rank is not None and pd.notna(autodock_rank):
            try:
                return int(float(autodock_rank))
            except (TypeError, ValueError):
                pass
        mt = re.search(r"_(?:model|pose)(\d+)", name)
        return int(mt.group(1)) if mt else None
    if m == "diffdock":
        mt = re.search(r"rank(\d+)", name)
        return int(mt.group(1)) if mt else None
    mt = re.search(r"(?:unguided|guided|pose|rank|model)_?(\d+)", name)
    return int(mt.group(1)) if mt else None


def _apply_equibind_gnina_rank(out: pd.DataFrame) -> int:
    """Override EquiBind pose_rank with a gnina-energy rank (most-negative
    gnina_affinity = rank 1) within each (frame, ligand). EquiBind emits no native
    confidence score, so its unguided sample index is not a quality rank; the gnina
    minimisation energy is. Poses with no gnina_affinity get pose_rank = None (not
    gnina-rankable). Returns the number of poses ranked. Non-EquiBind untouched."""
    if "gnina_affinity" not in out.columns:
        return 0
    eb = out["docking_method"].astype(str).str.startswith("equibind")
    if not eb.any():
        return 0
    gaff = pd.to_numeric(out["gnina_affinity"], errors="coerce")
    n = 0
    for _, idx in out[eb].groupby(["protein", "ligand"]).groups.items():
        sub = gaff.loc[idx]
        ranked = sub.dropna().sort_values(kind="stable")     # ascending: most-negative first
        for r, i in enumerate(ranked.index, 1):
            out.at[i, "pose_rank"] = r
            n += 1
        for i in sub[sub.isna()].index:                      # no gnina energy → unrankable
            out.at[i, "pose_rank"] = None
    return n


def make_figures(out: pd.DataFrame, out_dir: Path, name: str,
                 dd_variant: str = "gnina", eb_variant: str = "gnina") -> None:
    """Three PB-valid-focused figures + their data CSVs (see module docstring)."""
    if not _HAVE_MPL:
        print("  matplotlib unavailable — skipping figures")
        return
    df = out.copy()
    df["_base"] = df["docking_method"].map(_tool_key)
    if "variant" not in df.columns:
        df["variant"] = [_pose_variant(m, o, r, pn) for m, o, r, pn in zip(
            df["docking_method"], df.get("optimizer", ""), df.get("refine_variant", ""),
            df.get("pose_name", ""))]
    if "toolchain" not in df.columns:
        df["toolchain"] = [_toolchain_label(b, v) for b, v in zip(df["_base"], df["variant"])]
    # restrict the REPORTING to the selected optimizer variant per tool (default gnina);
    # per-pose CSVs still hold every variant. AutoDock has no variant.
    keep = ((df["_base"] == "autodock")
            | ((df["_base"] == "diffdock") & df["variant"].map(lambda v: _variant_matches(v, dd_variant)))
            | ((df["_base"] == "equibind") & df["variant"].map(lambda v: _variant_matches(v, eb_variant))))
    df = df[keep].copy()
    pv = df["pb_valid"] == True                                    # noqa: E712
    st = df["analysis_status"]
    df["_fate"] = np.where(~pv, "nonvalid",
                  np.where(st == "kept", "survive",
                  np.where(st == "transmembrane", "lost_tm", "unplaced")))
    # ordered toolchains present + a (colour, marker) per toolchain (colour = base tool)
    base_of = dict(zip(df["toolchain"], df["_base"]))
    order = [tc for base in ("autodock", "diffdock", "equibind")
             for tc in sorted(df.loc[df["_base"] == base, "toolchain"].unique())]
    _EXTRA_MK = ["D", "v", "P", "X", "*"]
    tc_style = {}
    for base in ("autodock", "diffdock", "equibind"):
        for i, tc in enumerate([t for t in order if base_of.get(t) == base]):
            mk = TOOL_STYLE[base][2] if i == 0 else _EXTRA_MK[(i - 1) % len(_EXTRA_MK)]
            tc_style[tc] = (TOOL_STYLE[base][1], mk)
    if not order:
        print("  no poses for the selected variant — skipping figures")
        return

    # ── Statistics (shared stats_utils; each block wrapped so a stats failure
    #    degrades to the current test-free figure and never breaks the pipeline).
    #    tm_stats.json next to the figures holds the full numeric results. ────────
    tm_stats: Dict = {
        "dataset": name,
        "unit_of_analysis_caveat": (
            "Counts below are individual POSES — pseudoreplicated: each tool emits "
            "many correlated poses per (frame, ligand) unit. The G-test / "
            "Cochran-Armitage p-values are computed on pose counts and therefore "
            "OVERSTATE significance; the honest independent unit is the "
            "(frame, ligand) pair (Orai x JKU has only ~3 ligands x 4 frames). "
            "Per-(frame,ligand) loss-fraction aggregates and n_units are reported "
            "alongside so the pooled result can be read against the real unit."),
    }
    # independent-unit count (frame x ligand pairs carrying >=1 PB-valid pose)
    try:
        _u = df.loc[pv, ["protein", "ligand"]].drop_duplicates()
        n_units = int(len(_u))
        n_ligands = int(df.loc[pv, "ligand"].nunique())
    except Exception:
        n_units = n_ligands = 0
    # Orai x JKU has only ~3-4 distinct ligands: too few independent chemotypes for
    # formal inference, so its conclusions are labelled exploratory. Benchmark's
    # hundreds of ligands are not. Trigger on few units OR few distinct ligands.
    small_n = (n_units < 12) or (n_ligands < 10)
    tm_stats["n_frame_ligand_units_with_pbvalid"] = n_units
    tm_stats["n_distinct_ligands"] = n_ligands
    tm_stats["small_n_exploratory"] = bool(small_n)

    # ── Figure 1: per-toolchain PB-valid fate (produced → survive vs lost-to-TM) ──
    counts = {tc: df[df["toolchain"] == tc]["_fate"].value_counts() for tc in order}
    totals = {tc: int((df["toolchain"] == tc).sum()) for tc in order}
    ymax = max(totals.values()) if totals else 1

    # toolchain x fate contingency table -> G-test of independence (+ Cramer's V,
    # adjusted residuals); per-tool TM-loss fraction with Wilson CI; per-(frame,
    # ligand) aggregate of the loss fraction (the pseudoreplication-free view).
    _fate_keys = [k for k, _lab, _c, _h in _FATE]     # survive, lost_tm, unplaced, nonvalid
    fig1_ci = {tc: (None, None) for tc in order}       # Wilson CI on lost/PB-valid
    fig1_resid = {tc: None for tc in order}            # adjusted residual, lost_tm cell
    fig1_flag = {tc: False for tc in order}            # |resid| > 2 on lost_tm
    g_title = ""
    perunit_title = ""       # pseudoreplication-free per-(frame,ligand) loss-test caption
    if _HAVE_SU:
        try:
            f1: Dict = {"unit": "poses (pseudoreplicated — see unit_of_analysis_caveat)",
                        "fate_categories": _fate_keys, "toolchains": list(order)}
            # Wilson CI on the per-tool TM-loss fraction (lost_tm / PB-valid)
            per_tool = {}
            for tc in order:
                nvalid = int(((df["toolchain"] == tc) & pv).sum())
                nlost = int(counts[tc].get("lost_tm", 0))
                lo, hi = su.wilson_ci(nlost, nvalid) if nvalid else (float("nan"), float("nan"))
                fig1_ci[tc] = (lo, hi)
                per_tool[tc] = {"pb_valid": nvalid, "lost_tm": nlost,
                                "loss_fraction": (nlost / nvalid) if nvalid else None,
                                "wilson_ci": [lo, hi]}
            f1["per_toolchain_loss_fraction"] = per_tool
            # G-test on the toolchain x fate table (drop all-zero rows/cols so
            # chi2_contingency has valid marginals)
            tbl = np.array([[int(counts[tc].get(k, 0)) for k in _fate_keys]
                            for tc in order], float)
            keep_r = tbl.sum(1) > 0
            keep_c = tbl.sum(0) > 0
            rows_k = [tc for tc, kr in zip(order, keep_r) if kr]
            cols_k = [k for k, kc in zip(_fate_keys, keep_c) if kc]
            sub = tbl[np.ix_(keep_r, keep_c)]
            if sub.shape[0] >= 2 and sub.shape[1] >= 2:
                g = su.gtest_independence(sub)
                resid = np.asarray(g["residuals"])
                f1["gtest"] = {"G": g["G"], "p": g["p"], "df": g["df"],
                               "cramers_v": g["cramers_v"],
                               "n_low_expected": g["n_low_expected"],
                               "rows": rows_k, "cols": cols_k,
                               "adjusted_residuals": resid.tolist()}
                # map the lost_tm-column residual back to each toolchain
                if "lost_tm" in cols_k:
                    jc = cols_k.index("lost_tm")
                    for ir, tc in enumerate(rows_k):
                        z = float(resid[ir, jc])
                        fig1_resid[tc] = z
                        fig1_flag[tc] = abs(z) > 2.0
                    f1["lost_tm_residual_flagged"] = [tc for tc in rows_k if fig1_flag[tc]]
                g_title = (f"G={g['G']:.1f}, {su.fmt_p(g['p'])} {su.p_stars(g['p'])}, "
                           f"Cramer's V={g['cramers_v']:.2f}")
            else:
                f1["gtest"] = {"skipped": "table has <2 non-empty rows/cols"}
            # pseudoreplication-free view: per-(frame,ligand) loss fraction per tool
            agg = {}
            for tc in order:
                s = df[(df["toolchain"] == tc) & pv]
                if s.empty:
                    continue
                fr = (s.assign(_l=(s["analysis_status"] == "transmembrane").astype(float))
                        .groupby(["protein", "ligand"])["_l"].mean())
                agg[tc] = {"n_units": int(fr.size),
                           "mean_per_unit_loss_fraction": float(fr.mean()) if fr.size else None,
                           "per_unit_loss_fraction": {f"{a}|{b}": float(v)
                                                      for (a, b), v in fr.items()}}
            f1["per_frame_ligand_loss_fraction"] = agg
            # pseudoreplication-FREE test: paired across tools on the per-(frame,
            # ligand) loss fraction (one value per unit per tool), listwise-complete
            # → Friedman + Kendall's W + Wilcoxon/Holm; plus a bootstrap CI on each
            # tool's MEAN per-unit loss fraction. Honest companion to the pooled-pose
            # G-test / Wilson CI above (whose unit is the pseudoreplicated pose).
            try:
                unit_mat = pd.DataFrame({tc: pd.Series(agg[tc]["per_unit_loss_fraction"])
                                         for tc in order if tc in agg})
                put: Dict = {"test": "friedman + wilcoxon (paired across tools)",
                             "quantity": "per-(frame,ligand) TM-loss fraction",
                             "unit": "(frame,ligand) pair"}
                boot = {}
                for tc in order:
                    vals = np.array(list(agg.get(tc, {}).get("per_unit_loss_fraction", {}).values()), float)
                    if vals.size:
                        est, blo, bhi = su.bootstrap_ci(vals, statistic=np.mean)
                        boot[tc] = {"mean": float(est), "ci": [float(blo), float(bhi)],
                                    "n_units": int(vals.size)}
                put["per_tool_mean_ci"] = boot
                complete = unit_mat.dropna()
                put["n_units_complete"] = int(len(complete))
                if len(complete) >= 5 and complete.shape[1] >= 3:
                    res = su.paired_continuous(complete, labels=list(complete.columns))
                    put.update(res)
                    om = res["omnibus"]
                    perunit_title = (f"per-(frame,ligand) loss Friedman {su.p_stars(om['p'])} "
                                     f"{su.fmt_p(om['p'])} (W={om['kendall_w']:.2f}, n={om['n']})")
                elif len(complete) >= 5 and complete.shape[1] == 2:
                    a, b = list(complete.columns)
                    rb, p, npair = su.wilcoxon_rankbiserial(complete[a].to_numpy(float),
                                                            complete[b].to_numpy(float))
                    put["pairwise"] = [{"a": a, "b": b, "rank_biserial": float(rb),
                                        "p_raw": float(p), "star": su.p_stars(p), "n": int(npair)}]
                    perunit_title = (f"per-(frame,ligand) loss Wilcoxon {su.p_stars(p)} "
                                     f"{su.fmt_p(p)} (n={npair})")
                else:
                    put["note"] = "n too small — exploratory"
                    perunit_title = "per-(frame,ligand) loss: n too small — exploratory"
                f1["per_frame_ligand_loss_test"] = put
            except Exception as e:                    # pragma: no cover
                f1["per_frame_ligand_loss_test"] = {"error": str(e)}
            if small_n:
                f1["note"] = "exploratory — few independent (frame,ligand) units"
            tm_stats["fig_tool_tm_loss"] = f1
        except Exception as e:                        # pragma: no cover
            print(f"  [warn] fig_tool_tm_loss stats failed: {e}")
            tm_stats["fig_tool_tm_loss"] = {"error": str(e)}

    fig, ax = plt.subplots(figsize=(1.9 * len(order) + 3.5, 6.2))
    x = np.arange(len(order)); bottoms = np.zeros(len(order))
    for key, label, colour, hatch in _FATE:
        vals = np.array([int(counts[tc].get(key, 0)) for tc in order], float)
        ax.bar(x, vals, bottom=bottoms, color=colour, width=0.6, label=label,
               hatch=hatch, edgecolor="white", linewidth=0.7)
        for xi, (v, b) in enumerate(zip(vals, bottoms)):
            if v >= 0.03 * ymax:
                ax.text(xi, b + v / 2, f"{int(v)}", ha="center", va="center",
                        fontsize=8, color="#222")
        bottoms += vals
    for xi, tc in enumerate(order):
        nvalid = int(((df["toolchain"] == tc) & pv).sum())
        nlost = int(counts[tc].get("lost_tm", 0))
        pct = 100 * nlost / nvalid if nvalid else 0
        lo, hi = fig1_ci.get(tc, (None, None))
        ci_txt = (f" [{100 * lo:.0f}–{100 * hi:.0f}]"
                  if lo is not None and lo == lo else "")   # Wilson 95% CI, %
        mark = " ‡" if fig1_flag.get(tc) else ""            # |adj. residual|>2 on lost_tm
        ax.text(xi, totals[tc] + 0.015 * ymax,
                f"produced {totals[tc]}\nPB-valid {nvalid}\n"
                f"lost to TM {nlost} ({pct:.0f}%{ci_txt}){mark}",
                ha="center", va="bottom", fontsize=8, color="#222")
    ax.set_xticks(x); ax.set_xticklabels(order)
    ax.set_ylabel("Number of poses produced")
    _sub = []
    if g_title:
        _sub.append(f"toolchain×fate {g_title}")
    if perunit_title:
        _sub.append(perunit_title + " — pseudoreplication-free unit")
    if any(fig1_flag.values()):
        _sub.append("‡ = |adjusted residual|>2 for the lost-to-TM cell")
    _sub.append("bracket = Wilson 95% CI on % PB-valid lost; counts are poses "
                "(pseudoreplicated)" + (" · exploratory (few units)" if small_n else ""))
    ax.set_title(f"PB-valid poses lost to the transmembrane vs surviving, per toolchain  ({name})",
                 fontsize=11)
    ax.margins(y=0.18)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.07), ncol=2,
              frameon=False, fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    if _sub:
        _half = (len(_sub) + 1) // 2                       # wrap footnote onto 2 lines
        _foot = "\n".join(("   ·   ".join(_sub[:_half]),
                           "   ·   ".join(_sub[_half:]))).strip()
        fig.text(0.5, -0.14, _foot, ha="center", va="top",
                 fontsize=8, color="0.35")
    fig.tight_layout()
    fig.savefig(out_dir / "fig_tool_tm_loss.png", dpi=130, bbox_inches="tight")
    plt.close(fig)

    # ── Figure 2: matrix of PB-valid poses lost to TM, per Orai frame × ligand ──
    lost = df[pv & (st == "transmembrane")]
    mat = lost.groupby(["protein", "ligand"]).size().unstack("ligand", fill_value=0)
    if mat.size == 0:
        mat = pd.DataFrame(0, index=sorted(df["protein"].unique()),
                           columns=sorted(df["ligand"].unique()))
    validm = (df[pv].groupby(["protein", "ligand"]).size().unstack("ligand", fill_value=0)
              .reindex(index=mat.index, columns=mat.columns, fill_value=0))
    MAX_COLS = 24
    col_note = ""
    if mat.shape[1] > MAX_COLS:
        top = mat.sum(0).sort_values(ascending=False).head(MAX_COLS).index
        mat, validm = mat[top], validm[top]
        col_note = f"  (top {MAX_COLS} ligands by loss of {df['ligand'].nunique()})"
    frames, ligs = list(mat.index), list(mat.columns)
    cmap = LinearSegmentedColormap.from_list("loss", ["#f7f7f7", "#C44E52"])
    fig, ax = plt.subplots(figsize=(max(6, 0.62 * len(ligs) + 3),
                                    max(3, 0.6 * len(frames) + 2)))
    im = ax.imshow(mat.values, cmap=cmap, aspect="auto",
                   vmin=0, vmax=max(1, int(mat.values.max())))
    ax.set_xticks(range(len(ligs)))
    ax.set_xticklabels(ligs, rotation=40, ha="right", fontsize=8)
    ax.set_yticks(range(len(frames)))
    ax.set_yticklabels([f.replace("Orai1WT-", "") for f in frames], fontsize=9)
    vmax = mat.values.max() if mat.values.size else 1
    for i in range(len(frames)):
        for j in range(len(ligs)):
            nl, nv = int(mat.values[i, j]), int(validm.values[i, j])
            if nv:
                cell = f"{nl}/{nv}\n{100 * nl / nv:.0f}%"
            else:
                cell = f"{nl}"
            ax.text(j, i, cell, ha="center", va="center", linespacing=1.15,
                    fontsize=7.5, color="white" if nl > 0.6 * vmax else "#333")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03).set_label(
        "PB-valid poses lost to transmembrane", fontsize=9)

    # Cochran-Armitage trend: does the fraction of PB-valid poses lost to the TM
    # rise/fall across the ORDERED MD frames (START-Fr0 → Fr300 → Fr400 → Fr499)?
    # Pooled across toolchains + per toolchain. Scores = the MD-time frame numbers.
    ca_title = ""
    if _HAVE_SU:
        try:
            f2: Dict = {"unit": "poses (pseudoreplicated — see unit_of_analysis_caveat)",
                        "frame_order_scores": "Fr number (MD time)"}
            frames_ord = sorted(df.loc[pv, "protein"].unique(), key=_frame_sort_key)
            scores = [_frame_sort_key(f) for f in frames_ord]

            def _ca(mask):
                succ = [int((mask & pv & (df["protein"] == f) & (st == "transmembrane")).sum())
                        for f in frames_ord]
                tot = [int((mask & pv & (df["protein"] == f)).sum()) for f in frames_ord]
                z, p, sgn = su.cochran_armitage(succ, tot, scores=scores)
                return {"frames": [str(f) for f in frames_ord], "scores": scores,
                        "lost": succ, "pb_valid": tot,
                        "z": None if z != z else float(z),
                        "p": None if p != p else float(p), "sign": int(sgn)}

            true_all = pd.Series(True, index=df.index)
            if len(frames_ord) >= 3:
                pooled = _ca(true_all)
                f2["pooled"] = pooled
                per_tc = {}
                for tc in order:
                    per_tc[tc] = _ca(df["toolchain"] == tc)
                f2["per_toolchain"] = per_tc
                if pooled["p"] is not None:
                    trend = ("↑" if pooled["sign"] > 0 else "↓" if pooled["sign"] < 0 else "flat")
                    ca_title = (f"Cochran–Armitage frame trend (pooled): z={pooled['z']:.2f}, "
                                f"{su.fmt_p(pooled['p'])} {su.p_stars(pooled['p'])} {trend}")
            else:
                f2["skipped"] = f"only {len(frames_ord)} ordered frame(s) — no trend test"
            if small_n:
                f2["note"] = "exploratory — few independent (frame,ligand) units"
            tm_stats["fig_pbvalid_lost_matrix"] = f2
        except Exception as e:                        # pragma: no cover
            print(f"  [warn] fig_pbvalid_lost_matrix stats failed: {e}")
            tm_stats["fig_pbvalid_lost_matrix"] = {"error": str(e)}

    _t2 = (f"PB-valid poses lost to the transmembrane, per Orai frame × ligand{col_note}\n"
           f"cell = lost / total PB-valid (% of PB-valid lost)   ({name})")
    if ca_title:
        _t2 += ("\n" + ca_title
                + ("  · exploratory" if small_n else "")
                + "  · poses (pseudorepl.)")
    ax.set_title(_t2, fontsize=10)
    fig.tight_layout()
    fig.savefig(out_dir / "fig_pbvalid_lost_matrix.png", dpi=130, bbox_inches="tight")
    plt.close(fig)

    # ── Figure 3: rank survival — PB-valid poses outside the TM, by rank ──
    # AutoDock = Vina affinity rank, DiffDock = confidence rank, EquiBind = gnina-
    # energy rank (assigned upstream, since EquiBind has no native score).
    rk = df[pv & df["pose_rank"].notna()].copy()
    fig, ax = plt.subplots(figsize=(8, 5))
    rank_rows = []
    for tc in order:
        sub = rk[rk["toolchain"] == tc]
        if sub.empty:
            continue
        sub = sub.assign(pose_rank=sub["pose_rank"].astype(int))
        total = sub.groupby("pose_rank").size()
        valuable = (sub[sub["analysis_status"] == "kept"].groupby("pose_rank").size()
                    .reindex(total.index, fill_value=0))
        colour, mk = tc_style[tc]
        ax.plot(total.index, valuable.values, "-", marker=mk, color=colour, lw=2, ms=6,
                label=tc)
        ax.plot(total.index, total.values, "--", color=colour, lw=1.2, alpha=0.4)
        for rr in total.index:
            rank_rows.append({"toolchain": tc, "rank": int(rr),
                              "pb_valid": int(total[rr]),
                              "pb_valid_outside_tm": int(valuable[rr])})
    ax.set_xlabel("Pose rank (1 = best-scored: AutoDock affinity / DiffDock confidence / "
                  "EquiBind gnina energy)")
    ax.set_ylabel("Number of PB-valid poses outside the transmembrane")
    ax.set_title(f"Potentially valuable poses by rank\n"
                 f"solid = PB-valid & outside TM   ·   dashed = all PB-valid   ({name})",
                 fontsize=11)
    # cap the x-axis at the AutoDock/DiffDock scoring depth so those tools stay
    # readable; EquiBind's deeper gnina-energy-ranked tail is clipped.
    if not rk.empty:
        ranked = rk[rk["_base"].isin(["autodock", "diffdock"])]["pose_rank"]
        xmax = int((ranked if not ranked.empty else rk["pose_rank"]).astype(int).max())
        ax.set_xlim(0.5, xmax + 0.5)
        ax.set_xticks(range(1, xmax + 1))
    ax.set_ylim(bottom=0)
    ax.grid(True, axis="y", alpha=0.25)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(title="Toolchain", frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(out_dir / "fig_rank_survival.png", dpi=130, bbox_inches="tight")
    plt.close(fig)

    # ── supporting data CSVs ──
    rows = []
    for tc in order:
        sub = df[df["toolchain"] == tc]; c = sub["_fate"].value_counts()
        nvalid = int((sub["pb_valid"] == True).sum())                # noqa: E712
        rows.append({"toolchain": tc, "produced": int(len(sub)), "pb_valid": nvalid,
                     "pbvalid_survive_outside_tm": int(c.get("survive", 0)),
                     "pbvalid_lost_to_tm": int(c.get("lost_tm", 0)),
                     "pbvalid_offprotein_or_uneval": int(c.get("unplaced", 0)),
                     "pct_pbvalid_lost_to_tm": round(100 * c.get("lost_tm", 0) / nvalid, 1)
                     if nvalid else 0.0})
    pd.DataFrame(rows).to_csv(out_dir / "tool_tm_loss_summary.csv", index=False)
    mat.to_csv(out_dir / "pbvalid_lost_matrix.csv")
    pd.DataFrame(rank_rows).to_csv(out_dir / "rank_survival.csv", index=False)

    # fig_rank_survival: descriptive only (no inferential test drawn on the figure,
    # per validation plan §8 — a tool-vs-tool AUC permutation would need many more
    # independent units than Orai x JKU offers). Record the per-rank survival counts
    # so the descriptive claim is recoverable from the sidecar.
    if _HAVE_SU:
        try:
            tm_stats["fig_rank_survival"] = {
                "test": "none — descriptive",
                "note": ("PB-valid poses outside the TM by native pose rank; solid = "
                         "outside TM, dashed = all PB-valid. See rank_survival.csv."),
                "rows": rank_rows}
        except Exception as e:                        # pragma: no cover
            print(f"  [warn] fig_rank_survival stats note failed: {e}")

    # sidecar: full numeric statistical results next to the figures
    if _HAVE_SU:
        try:
            def _san(o):
                if isinstance(o, dict):
                    return {k: _san(v) for k, v in o.items()}
                if isinstance(o, (list, tuple)):
                    return [_san(v) for v in o]
                if isinstance(o, (np.floating,)):
                    return float(o)
                if isinstance(o, (np.integer,)):
                    return int(o)
                if isinstance(o, np.ndarray):
                    return o.tolist()
                return o
            (out_dir / "tm_stats.json").write_text(json.dumps(_san(tm_stats), indent=2))
            print(f"  stats → {out_dir/'tm_stats.json'}")
        except Exception as e:                        # pragma: no cover
            print(f"  [warn] writing tm_stats.json failed: {e}")

    print("  figures → fig_tool_tm_loss.png, fig_pbvalid_lost_matrix.png, fig_rank_survival.png")


# ---------------------------------------------------------------------------
# Per-dataset driver
# ---------------------------------------------------------------------------
def run_dataset(name: str, csv_path: Path, out_dir: Path, pad: float, workers: int,
                ligands: Optional[List[str]] = None,
                dd_variant: str = "gnina", eb_variant: str = "gnina"):
    print(f"\n=== {name} ===")
    print(f"reading {csv_path}")
    df = pd.read_csv(csv_path, low_memory=False)
    # essential columns
    for req in ("protein", "pose_file", "docking_method", "ligand"):
        if req not in df.columns:
            raise SystemExit(f"{csv_path}: missing required column '{req}'")

    # ── optional ligand restriction (case-insensitive substring match) ──
    if ligands:
        pats = [str(p).lower() for p in ligands]
        lig_lower = df["ligand"].astype(str).str.lower()
        mask = lig_lower.apply(lambda s: any(p in s for p in pats))
        matched = sorted(df.loc[mask, "ligand"].unique())
        if not matched:
            avail = sorted(df["ligand"].unique())
            raise SystemExit(f"--ligands {ligands} matched no ligand in {name}. "
                             f"Available ({len(avail)}): {avail}")
        print(f"  --ligands: keeping {int(mask.sum())} of {len(df)} poses across "
              f"{len(matched)} ligand(s): {matched}")
        df = df[mask].copy()

    n0 = len(df)
    df["pb_valid"] = _derive_pb_valid(df)

    # build a zone per frame from that frame's own receptor file
    out_dir.mkdir(parents=True, exist_ok=True)
    zones: Dict[str, Dict] = {}
    for frame, sub in df.groupby("protein"):
        recs = sub["protein_file_used"].dropna().unique() if "protein_file_used" in df.columns else []
        rec = None
        for cand in recs:
            if Path(cand).exists():
                rec = Path(cand); break
        if rec is None:  # fall back to canonical receptor location
            for cand in (Path(f"Data/Receptors/{frame}.pdb"),
                         Path(f"posebusters_results/_{name}_staging/receptors/{frame}.pdb")):
                if cand.exists():
                    rec = cand; break
        if rec is None:
            print(f"  !! no receptor file for frame {frame}; skipping its poses")
            continue
        z = build_zone(rec, pad)
        zones[str(frame)] = z
        _rout = f"{z['r_outer']:.1f}" if z.get("r_outer") is not None else "n/a"
        print(f"  {frame}: R91 z-end {np.round(z['p_R91'],2).tolist()}  "
              f"E106 z-end {np.round(z['p_E106'],2).tolist()}  L={z['length']:.1f}  "
              f"r_pore={z['r_pore_mean']:.1f}  r_outer={_rout}  "
              f"axis|PCA={z['axis_vs_pca_deg']}deg")
        write_zone_pdb(z, out_dir / f"{frame}_tm_zone.pdb")
        write_zone_pml(z, str(frame), out_dir / f"{frame}_tm_zone.pml")
    (out_dir / "tm_exclusion_zones.json").write_text(json.dumps(
        {"pad": pad, "definition": "axial slab between R91 and E106 CA rings; in_transmembrane = ligand centroid in slab",
         "zones": zones}, indent=2))

    # classify every pose (parallel)
    tasks = [(i, str(r.protein), str(r.pose_file)) for i, r in
             enumerate(df.itertuples(index=False))]
    results: List[Optional[Dict]] = [None] * len(tasks)
    print(f"  classifying {len(tasks)} poses on {workers} workers ...")
    with ProcessPoolExecutor(max_workers=workers, initializer=_init_worker,
                             initargs=(zones,)) as ex:
        futs = [ex.submit(_process_pose, t) for t in tasks]
        for fut in as_completed(futs):
            r = fut.result()
            results[r["_idx"]] = r
    geom = pd.DataFrame(results).drop(columns=["_idx"])

    carry = [c for c in _CARRY if c in df.columns]
    if "pb_valid" in df.columns and "pb_valid" not in carry:
        carry.append("pb_valid")
    out = pd.concat([df[carry].reset_index(drop=True), geom.reset_index(drop=True)], axis=1)
    out.insert(0, "dataset", name)

    # optimizer/refinement variant + display toolchain (asterisk marks gnina)
    _op = out["optimizer"] if "optimizer" in out.columns else pd.Series([None] * len(out))
    _rv = out["refine_variant"] if "refine_variant" in out.columns else pd.Series([None] * len(out))
    _pn = out["pose_name"] if "pose_name" in out.columns else pd.Series([None] * len(out))
    _ar = out["autodock_rank"] if "autodock_rank" in out.columns else pd.Series([None] * len(out))
    out["variant"] = [_pose_variant(m, o, r, pn) for m, o, r, pn in
                      zip(out["docking_method"], _op, _rv, _pn)]
    out["toolchain"] = [_toolchain_label(_tool_key(m), v) for m, v in
                        zip(out["docking_method"], out["variant"])]

    # native pose rank (AutoDock Vina affinity rank, DiffDock confidence rank);
    # EquiBind has none, so it is (re)ranked by gnina energy just below.
    out["pose_rank"] = [_pose_rank(m, pn, pf, ar) for m, pn, pf, ar in
                        zip(out["docking_method"], _pn, out["pose_file"], _ar)]
    n_gr = _apply_equibind_gnina_rank(out)
    if n_gr:
        print(f"  ranked {n_gr} EquiBind poses by gnina energy (most-negative = rank 1)")

    # ── physical-placement gate (drop poses not actually on the protein) ──────
    if PLACEMENT_FLAG in df.columns:
        raw = df[PLACEMENT_FLAG].reset_index(drop=True)
        flag_map = {"True": True, "true": True, "False": False, "false": False,
                    True: True, False: False, 1: True, 0: False, 1.0: True, 0.0: False}
        pld = raw if raw.dtype == bool else raw.map(flag_map)
        pld_ok = pld.fillna(False).astype(bool)
    else:
        pld_ok = pd.Series(True, index=out.index)          # no PB flag → geometry only
    geo_sane = (out["radial_dist"] <= RADIAL_MAX) & out["axial_frac"].between(AXIAL_FRAC_LO, AXIAL_FRAC_HI)
    out["placed_on_protein"] = out["geom_ok"].astype(bool) & pld_ok & geo_sane.fillna(False)
    out["analysis_status"] = np.where(
        ~out["geom_ok"].astype(bool), "unevaluable",
        np.where(~out["placed_on_protein"], "invalid_placement",
                 np.where(out["in_transmembrane"] == True, "transmembrane", "kept")))  # noqa: E712

    n_bad = int((out["analysis_status"] == "unevaluable").sum())
    n_invalid = int((out["analysis_status"] == "invalid_placement").sum())
    if n_bad:
        print(f"  {n_bad} poses unevaluable (unreadable/corrupt/no-zone)")
    if n_invalid:
        print(f"  {n_invalid} poses dropped — not placed on protein (translated/exploded off the pocket)")

    class_csv = out_dir / "tm_pose_classification.csv"
    out.to_csv(class_csv, index=False)

    placed = out[out["analysis_status"].isin(["transmembrane", "kept"])].copy()
    excluded = out[out["analysis_status"] == "transmembrane"].copy()
    kept = out[out["analysis_status"] == "kept"].copy()
    invalid = out[out["analysis_status"] == "invalid_placement"].copy()
    _sort = ["docking_method", "protein", "ligand"]
    excluded.sort_values(_sort).to_csv(out_dir / "tm_excluded_poses.csv", index=False)
    kept.sort_values(_sort).to_csv(out_dir / "tm_kept_poses.csv", index=False)
    invalid.sort_values(_sort).to_csv(out_dir / "tm_invalid_poses.csv", index=False)

    # summary — computed over PLACED poses only (kept + transmembrane)
    def _by(col, frame=placed):
        res = {}
        for k, v in frame.groupby(col)["analysis_status"]:
            n = int(v.size); exc = int((v == "transmembrane").sum())
            res[str(k)] = {"n": n, "excluded": exc, "pct": round(100 * exc / max(n, 1), 1)}
        return res
    # reporting subset: the selected optimizer variant per tool (default gnina)
    _base = placed["docking_method"].map(_tool_key)
    sel = ((_base == "autodock")
           | ((_base == "diffdock") & placed["variant"].map(lambda v: _variant_matches(v, dd_variant)))
           | ((_base == "equibind") & placed["variant"].map(lambda v: _variant_matches(v, eb_variant))))
    reported = placed[sel]
    summary = {
        "dataset": name, "n_poses_total": n0,
        "n_evaluated": int(len(placed)),               # placed on protein (kept + TM), ALL variants
        "n_unevaluable": n_bad,
        "n_invalid_placement": n_invalid,              # parsed but not on the protein
        "n_excluded_transmembrane": int(len(excluded)),
        "pct_excluded": round(100 * len(excluded) / max(len(placed), 1), 2),
        "n_in_pore_lumen": int(placed["in_pore_lumen"].sum()),
        "n_excluded_but_pb_valid": int((excluded["pb_valid"] == True).sum())  # noqa: E712
            if "pb_valid" in excluded else None,
        "placement_gate": PLACEMENT_FLAG,
        "reported_diffdock_variant": dd_variant,
        "reported_equibind_variant": eb_variant,
        "by_toolchain": _by("toolchain", reported),    # selected variant, asterisk labels
        "by_method": _by("docking_method"),            # all variants, base tool
        "by_frame": _by("protein"),
        "pad_angstrom": pad,
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2))
    print(f"  -> {len(excluded)}/{len(placed)} placed poses inside the TM layer "
          f"({summary['pct_excluded']}%); {summary['n_in_pore_lumen']} in the pore lumen")
    print(f"  wrote {class_csv}")
    print(f"        {out_dir/'tm_excluded_poses.csv'}  ({len(excluded)})")
    print(f"        {out_dir/'tm_kept_poses.csv'}  ({len(kept)})")
    print(f"        {out_dir/'tm_invalid_poses.csv'}  ({len(invalid)})")
    print(f"  reported toolchains (diffdock={dd_variant}, equibind={eb_variant}):")
    for m, s in summary["by_toolchain"].items():
        print(f"      {m:18s} {s['excluded']:5d}/{s['n']:<5d} ({s['pct']}%) in TM")

    try:
        make_figures(out, out_dir, name, dd_variant, eb_variant)
    except Exception as e:                          # never let plotting kill the run
        print(f"  [warn] figure generation failed: {e}")
    return summary


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dataset", choices=list(DATASETS) + ["all"], default="all")
    ap.add_argument("--pad", type=float, default=0.0,
                    help="Angstrom padding added beyond each residue ring plane (default 0).")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--ligands", nargs="+", default=None, metavar="NAME",
                    help="restrict the analysis to ligands whose name contains any of "
                         "these (case-insensitive substrings); default = all ligands.")
    ap.add_argument("--diffdock-variant", default="gnina",
                    choices=["gnina", "smina", "raw", "original", "all"],
                    help="DiffDock optimizer variant used in the figures/table/summary "
                         "(default gnina → labelled 'DiffDock*'). CSVs keep all variants.")
    ap.add_argument("--equibind-variant", default="gnina",
                    choices=["gnina", "smina", "raw", "all"],
                    help="EquiBind refine variant used in the figures/table/summary "
                         "(default gnina → 'EquiBind*', ranked by gnina energy).")
    ap.add_argument("--out-root", default="posebusters_results",
                    help="root under which <dataset>/transmembrane_filter is written")
    args = ap.parse_args(argv)

    names = list(DATASETS) if args.dataset == "all" else [args.dataset]
    summaries = []
    for name in names:
        csv_path = Path(DATASETS[name])
        if not csv_path.exists():
            print(f"skip {name}: {csv_path} not found")
            continue
        out_dir = Path(args.out_root) / name / "transmembrane_filter"
        summaries.append(run_dataset(name, csv_path, out_dir, args.pad, args.workers,
                                     ligands=args.ligands,
                                     dd_variant=args.diffdock_variant,
                                     eb_variant=args.equibind_variant))
    print("\n==== DONE ====")
    for s in summaries:
        print(f"{s['dataset']}: {s['n_excluded_transmembrane']}/{s['n_evaluated']} "
              f"({s['pct_excluded']}%) excluded as transmembrane")


if __name__ == "__main__":
    main()
