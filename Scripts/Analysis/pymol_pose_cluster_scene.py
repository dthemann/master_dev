#!/usr/bin/env python3
"""Build a PyMOL scene of docked poses + their clusters on an Orai receptor.

For a chosen Orai receptor frame this script gathers every docked pose (from the
clustered ``per_pose.csv`` produced by ``orai_pose_cluster_report.py``) and writes:

  * ``<dataset>__<frame>__poses.pdb`` — one HETATM residue per pose, with the
    tool chain, ligand and cluster encoded in standard PDB fields so the poses
    stay distinguishable even in a bare viewer:

        chain ID   → tool chain      (A = AutoDock Vina, D = DiffDock, E = EquiBind …)
        segID      → ligand          (L01, L02, … ; see the REMARK legend / .pml)
        resName    → ligand CCD code (cosmetic, human-readable label)
        resSeq     → unique pose index within the scene
        tempFactor → cluster id      (per (frame, ligand); ``spectrum b`` colours it)
        occupancy  → PoseBusters validity (1 = valid, 0 = invalid, 0.5 = unknown)

    Bonds are carried as CONECT records taken straight from each SDF bond block,
    so sticks render exactly (no distance guessing).

  * ``<dataset>__<frame>__scene.pml`` — a PyMOL script that loads the matching
    receptor + the poses PDB, hides hydrogens, and builds named selection groups
    for every Tool, Ligand and Cluster, plus labelled pseudo-atoms at each
    cluster centre. Three one-word commands recolour the whole scene:

        by_tool      colour every pose by the tool chain that produced it
        by_ligand    colour every pose by ligand
        by_cluster   give every (ligand, cluster) group its own colour

  * ``<dataset>__<frame>__legend.csv`` — the pose→(tool, ligand, cluster, fields)
    mapping, including any poses skipped (with the reason).

Because clusters are defined per (frame, ligand) and a full benchmark frame holds
~15 000 poses, each scene covers ONE frame and a bounded set of ligands (all of
them for small ligand sets such as the JKU actives; the busiest ``--max-ligands``
otherwise — the cap is always logged, never silent).

Usage
-----
    # JKU actives (3 ligands) — one scene per frame, all ligands
    python Scripts/Analysis/pymol_pose_cluster_scene.py \
        --per-pose-csv posebusters_results/orai_jku/pose_clusters/per_pose.csv \
        --out-dir posebusters_results/orai_jku/pymol_scene

    # Benchmark — one frame, one ligand
    python Scripts/Analysis/pymol_pose_cluster_scene.py \
        --frames Orai1WT-START-Fr0 --ligands 5SAK_ZRY \
        --out-dir posebusters_results/orai_benchmark/pymol_scene

Then in PyMOL:  ``pymol posebusters_results/.../..._scene.pml``
"""
from __future__ import annotations

import argparse
import colorsys
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

# ── tool → (pretty name, colour) — matches Scripts/Analysis/orai_pose_cluster_report.py
TOOL_STYLE = {
    "autodock": ("AutoDock Vina", "#4C72B0"),
    "diffdock": ("DiffDock", "#55A868"),
    "equibind": ("EquiBind", "#C44E52"),
}
_TOOL_FALLBACK = ("#8172B3", "#CCB974", "#64B5CD", "#937860", "#DA8BC3", "#8C8C8C")
# preferred single-letter chain per known tool; others get spare letters
_PREF_CHAIN = {"autodock": "A", "diffdock": "D", "equibind": "E"}
_SPARE_CHAINS = list("BCFGHIJKMNOPQRSTUVWXYZ")

# 20-colour qualitative palette for per-ligand colouring (tab20-ish, hex)
_LIG_PALETTE = [
    "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", "#8C564B",
    "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF", "#AEC7E8", "#FFBB78",
    "#98DF8A", "#FF9896", "#C5B0D5", "#C49C94", "#F7B6D2", "#DBDB8D",
    "#9EDAE5", "#393B79",
]

# PDB fixed-format limits
_MAX_SERIAL = 99999          # ATOM/CONECT serial field is 5 chars
_MAX_RESSEQ = 9999           # resSeq field is 4 chars
# The signed %8.3f x/y/z field holds 8 chars: positive up to "9999.999", but only
# down to "-999.999" on the negative side (the minus sign costs a column). A value
# outside [-999.999, 9999.999] can't be written in fixed-column PDB without shifting
# every field to its right, so such poses are skipped.
_PDB_COORD_MIN = -999.999
_PDB_COORD_MAX = 9999.999


def _coords_fit_pdb(coords: np.ndarray) -> bool:
    return bool(np.all(np.isfinite(coords))
                and coords.min() >= _PDB_COORD_MIN
                and coords.max() <= _PDB_COORD_MAX)


# ════════════════════════════════════════════════════════════════════════
# SDF reading
# ════════════════════════════════════════════════════════════════════════
def read_sdf(path: Path) -> Optional[Tuple[List[Tuple[str, float, float, float]],
                                            List[Tuple[int, int]]]]:
    """Parse the first molecule of a V2000 SDF.

    Returns (atoms, bonds) where atoms = [(element, x, y, z), …] and bonds are
    0-based (i, j) pairs, or None if the file cannot be parsed.
    """
    try:
        lines = path.read_text(errors="ignore").splitlines()
    except OSError:
        return None
    if len(lines) < 4:
        return None
    counts = lines[3]
    if "V2000" not in counts:
        return None                      # V3000 / malformed — skip
    try:
        n_atom = int(counts[0:3])
        n_bond = int(counts[3:6])
    except ValueError:
        return None
    if n_atom <= 0 or len(lines) < 4 + n_atom + n_bond:
        return None
    atoms: List[Tuple[str, float, float, float]] = []
    for i in range(4, 4 + n_atom):
        L = lines[i]
        try:
            x = float(L[0:10]); y = float(L[10:20]); z = float(L[20:30])
        except ValueError:
            return None
        elem = L[31:34].strip() or "C"
        atoms.append((elem, x, y, z))
    bonds: List[Tuple[int, int]] = []
    for i in range(4 + n_atom, 4 + n_atom + n_bond):
        L = lines[i]
        try:
            a = int(L[0:3]) - 1
            b = int(L[3:6]) - 1
        except ValueError:
            continue
        if 0 <= a < n_atom and 0 <= b < n_atom and a != b:
            bonds.append((a, b))
    return atoms, bonds


# ════════════════════════════════════════════════════════════════════════
# PDB writing helpers
# ════════════════════════════════════════════════════════════════════════
def _atom_name_field(elem: str, idx: int) -> str:
    """4-char PDB atom-name field, obeying the 1-char-element leading-space rule.

    A 1-char element normally sits in col 14 with a leading space (" C12"); once the
    element+index would exceed 3 chars (idx>=100) the leading space is dropped so all
    four columns are used ("C100") instead of silently truncating a digit.
    """
    name = f"{elem}{idx}"
    if len(elem) == 1 and len(name) <= 3:
        return f" {name:<3}"[:4]
    return name[:4].ljust(4)


def _pdb_atom_line(serial: int, elem: str, atom_idx: int, resn: str, chain: str,
                   resseq: int, x: float, y: float, z: float,
                   occ: float, bfac: float, segid: str) -> str:
    name = _atom_name_field(elem, atom_idx)
    return (
        "HETATM"
        f"{serial:5d}"
        " "
        f"{name}"                      # cols 13-16 (already 4 chars)
        " "                            # altLoc
        f"{resn[:3]:>3}"
        " "
        f"{chain[:1]:1}"
        f"{resseq:4d}"
        " "                            # iCode
        "   "                          # cols 28-30
        f"{x:8.3f}{y:8.3f}{z:8.3f}"
        f"{occ:6.2f}{bfac:6.2f}"
        "      "                        # cols 67-72
        f"{segid[:4]:<4}"
        f"{elem[:2]:>2}"
    )


def _conect_lines(serials: List[int], bonds: List[Tuple[int, int]]) -> List[str]:
    """CONECT records from 0-based bonds, mapped through per-atom global serials."""
    nbr: Dict[int, List[int]] = {}
    for a, b in bonds:
        nbr.setdefault(serials[a], []).append(serials[b])
        nbr.setdefault(serials[b], []).append(serials[a])
    out = []
    for s in sorted(nbr):
        partners = nbr[s]
        for k in range(0, len(partners), 4):        # max 4 partners per CONECT line
            chunk = partners[k:k + 4]
            out.append("CONECT" + f"{s:5d}" + "".join(f"{p:5d}" for p in chunk))
    return out


# ════════════════════════════════════════════════════════════════════════
# Encoding maps
# ════════════════════════════════════════════════════════════════════════
def build_chain_map(tools: List[str]) -> Dict[str, str]:
    used = set()
    cmap: Dict[str, str] = {}
    for t in tools:                                  # known tools first, stable letters
        if t in _PREF_CHAIN and _PREF_CHAIN[t] not in used:
            cmap[t] = _PREF_CHAIN[t]
            used.add(_PREF_CHAIN[t])
    spare = [c for c in _SPARE_CHAINS if c not in used]
    for t in tools:
        if t in cmap:
            continue
        if spare:
            cmap[t] = spare.pop(0)
        else:                                        # >25 tool chains (never here: 3)
            cmap[t] = "Z"
            print(f"  WARN: ran out of chain letters; tool {t!r} collapsed to chain Z")
        used.add(cmap[t])
    return cmap


def _pb_to_occ(v) -> float:
    """PoseBusters validity -> occupancy: True->1.0, False->0.0, unknown/NaN->0.5.
    Works whether the value is a python/numpy bool, 0/1, or the canonical True/False
    written by load_per_pose."""
    if v is True or v is False:
        return 1.0 if v else 0.0
    try:
        if v == 1:
            return 1.0
        if v == 0:
            return 0.0
    except TypeError:
        pass
    return 0.5


def ligand_ccd(ligand: str) -> str:
    """3-char cosmetic residue code — CCD token if it looks like one, else first alnum."""
    tail = ligand.split("_")[-1]
    if 1 <= len(tail) <= 3 and tail.isalnum():
        return tail.upper()
    alnum = "".join(ch for ch in ligand if ch.isalnum())
    return (alnum[:3] or "LIG").upper()


# ════════════════════════════════════════════════════════════════════════
# Scene assembly for one (frame, ligand-set)
# ════════════════════════════════════════════════════════════════════════
def build_scene(df_frame: pd.DataFrame, frame: str, receptor_pdb: Optional[Path],
                out_dir: Path, dataset: str, max_offset: float,
                max_extent: float) -> Optional[dict]:
    """Write poses.pdb + scene.pml + legend.csv for one receptor frame."""
    ligands = list(dict.fromkeys(df_frame["ligand"].tolist()))   # preserve order
    tools_present = sorted(df_frame["tool"].unique().tolist())
    chain_map = build_chain_map(tools_present)
    lig_tag = {lig: f"L{i + 1:02d}" for i, lig in enumerate(ligands)}
    lig_color = {lig: _LIG_PALETTE[i % len(_LIG_PALETTE)] for i, lig in enumerate(ligands)}

    # receptor centroid — for the max-offset sanity filter
    rec_center = None
    if receptor_pdb and receptor_pdb.exists():
        xs = []
        for L in receptor_pdb.read_text(errors="ignore").splitlines():
            if L.startswith(("ATOM", "HETATM")):
                try:
                    xs.append([float(L[30:38]), float(L[38:46]), float(L[46:54])])
                except ValueError:
                    continue
        if xs:
            rec_center = np.asarray(xs).mean(axis=0)

    atom_lines: List[str] = []
    conect_lines: List[str] = []
    legend_rows: List[dict] = []
    serial = 0
    pose_idx = 0
    # cluster centre accumulation: (ligand, cluster) -> list of pose centroids
    cluster_pts: Dict[Tuple[str, int], List[np.ndarray]] = {}
    # pose index lists per (ligand, cluster) for exact PyMOL selections
    cluster_resi: Dict[Tuple[str, int], List[int]] = {}
    # poses dropped from the scene, per (ligand, cluster) — for on-scene disclosure
    cluster_dropped: Dict[Tuple[str, int], int] = {}
    n_skip_parse = n_skip_offset = n_skip_overflow = n_skip_extent = n_skip_capacity = 0

    def _drop(base, reason, counter_delta):
        cluster_dropped[(base["ligand"], base["cluster"])] = \
            cluster_dropped.get((base["ligand"], base["cluster"]), 0) + 1
        legend_rows.append({**base, "included": False, "reason": reason})

    for row in df_frame.itertuples(index=False):
        pose_file = Path(str(row.pose_file))
        base = {
            "pose_index": None, "tool": row.tool, "ligand": row.ligand,
            "ccd": ligand_ccd(row.ligand), "cluster": int(row.cluster),
            "pb_valid": row.pb_valid, "chain": chain_map[row.tool],
            "segid": lig_tag[row.ligand], "pose_file": str(pose_file),
        }
        parsed = read_sdf(pose_file)
        if parsed is None:
            n_skip_parse += 1
            _drop(base, "sdf_parse_failed", None)
            continue
        atoms, bonds = parsed
        coords = np.asarray([(x, y, z) for _, x, y, z in atoms], dtype=float)
        cen = coords.mean(axis=0)
        if not _coords_fit_pdb(coords):
            n_skip_overflow += 1
            _drop(base, "coord_overflow", None)
            continue
        if rec_center is not None and np.linalg.norm(cen - rec_center) > max_offset:
            n_skip_offset += 1
            _drop(base, f"offset>{max_offset:g}A", None)
            continue
        # geometry-sanity: max atom->centroid radius. Real drug-like poses are
        # <~12 A; broken EquiBind geometry explodes to 50-100+ A and renders as a
        # box-spanning web that wrecks the view framing.
        extent = float(np.linalg.norm(coords - cen, axis=1).max())
        if extent > max_extent:
            n_skip_extent += 1
            _drop(base, f"distorted_geometry(extent={extent:.0f}A)", None)
            continue
        # capacity guards (PDB fixed-format fields)
        if serial + len(atoms) > _MAX_SERIAL or pose_idx + 1 > _MAX_RESSEQ:
            n_skip_capacity += 1
            _drop(base, "pdb_capacity", None)
            continue

        pose_idx += 1
        resseq = pose_idx
        occ = _pb_to_occ(row.pb_valid)
        bfac = float(int(row.cluster))
        pose_serials: List[int] = []
        elem_count: Dict[str, int] = {}
        for elem, x, y, z in atoms:
            serial += 1
            elem_count[elem] = elem_count.get(elem, 0) + 1
            atom_lines.append(_pdb_atom_line(
                serial, elem, elem_count[elem], base["ccd"], base["chain"],
                resseq, x, y, z, occ, bfac, base["segid"]))
            pose_serials.append(serial)
        conect_lines.extend(_conect_lines(pose_serials, bonds))

        key = (row.ligand, int(row.cluster))
        cluster_pts.setdefault(key, []).append(cen)
        cluster_resi.setdefault(key, []).append(resseq)
        legend_rows.append({**base, "pose_index": resseq, "included": True,
                            "reason": "", "n_atoms": len(atoms),
                            "cx": round(float(cen[0]), 3), "cy": round(float(cen[1]), 3),
                            "cz": round(float(cen[2]), 3)})

    n_incl = pose_idx
    if n_incl == 0:
        print(f"  [{frame}] no drawable poses — skipped")
        return None

    out_dir.mkdir(parents=True, exist_ok=True)
    stem = f"{dataset}__{frame}"
    pdb_path = out_dir / f"{stem}__poses.pdb"
    pml_path = out_dir / f"{stem}__scene.pml"
    legend_path = out_dir / f"{stem}__legend.csv"

    # ── poses.pdb ─────────────────────────────────────────────────────────
    hdr = [
        f"REMARK 300 Docked-pose scene for Orai frame {frame} (dataset: {dataset})",
        "REMARK 300 Encoding:  chain=tool  segID=ligand  resName=ligand CCD  "
        "resSeq=pose index",
        "REMARK 300           tempFactor=cluster id   occupancy=PoseBusters valid "
        "(1/0/0.5)",
        "REMARK 300 Tool chains:",
    ]
    for t in tools_present:
        hdr.append(f"REMARK 300   chain {chain_map[t]} = {TOOL_STYLE.get(t, (t.title(),))[0]}")
    hdr.append("REMARK 300 Ligands:")
    for lig in ligands:
        hdr.append(f"REMARK 300   segID {lig_tag[lig]} = {lig} (resName {ligand_ccd(lig)})")
    pdb_path.write_text("\n".join(hdr + atom_lines + conect_lines + ["END", ""]))

    # ── legend.csv ────────────────────────────────────────────────────────
    pd.DataFrame(legend_rows).to_csv(legend_path, index=False)

    # ── scene.pml ─────────────────────────────────────────────────────────
    pml = _render_pml(stem, pdb_path.name, receptor_pdb, frame, dataset,
                      tools_present, chain_map, ligands, lig_tag, lig_color,
                      cluster_pts, cluster_resi, cluster_dropped)
    pml_path.write_text(pml)

    n_clusters = len({k for k in cluster_pts})
    print(f"  [{frame}] {n_incl} poses, {len(ligands)} ligand(s), "
          f"{len(tools_present)} tool(s), {n_clusters} (ligand,cluster) groups "
          f"-> {pml_path.name}")
    skipped = (n_skip_parse + n_skip_offset + n_skip_overflow
               + n_skip_extent + n_skip_capacity)
    if skipped:
        print(f"        skipped {skipped} poses "
              f"(parse={n_skip_parse}, offset>{max_offset:g}A={n_skip_offset}, "
              f"distorted>{max_extent:g}A={n_skip_extent}, "
              f"coord_overflow={n_skip_overflow}, capacity={n_skip_capacity}) "
              f"— see {legend_path.name}")
    return {"frame": frame, "n_poses": n_incl, "pdb": pdb_path, "pml": pml_path,
            "legend": legend_path, "n_ligands": len(ligands),
            "n_clusters": n_clusters, "skipped": skipped}


def _hex_pymol(h: str) -> str:
    """'#55A868' -> '0x55A868' for PyMOL colour arguments."""
    return "0x" + h.lstrip("#")


def _distinct_hex(n: int) -> List[str]:
    """n visually-distinct colours (evenly-spaced hues) as '#RRGGBB'."""
    out = []
    for i in range(max(0, n)):
        r, g, b = colorsys.hsv_to_rgb(i / max(1, n), 0.62, 0.88)
        out.append(f"#{int(r * 255):02X}{int(g * 255):02X}{int(b * 255):02X}")
    return out


def _render_pml(stem: str, pdb_name: str, receptor_pdb: Optional[Path], frame: str,
                dataset: str, tools: List[str], chain_map: Dict[str, str],
                ligands: List[str], lig_tag: Dict[str, str],
                lig_color: Dict[str, str],
                cluster_pts: Dict[Tuple[str, int], List[np.ndarray]],
                cluster_resi: Dict[Tuple[str, int], List[int]],
                cluster_dropped: Dict[Tuple[str, int], int]) -> str:
    def _safe(name: str) -> str:
        return "".join(c if (c.isalnum() or c == "_") else "_" for c in name)

    # unique, readable selection stem per ligand — disambiguate any two ligand
    # names that sanitise to the same string (e.g. differ only in '-' vs '.') with
    # the already-unique segID tag, so selections never silently overwrite.
    lig_stem: Dict[str, str] = {}
    seen: Dict[str, int] = {}
    for lig in ligands:
        s = _safe(lig)
        if s in seen:
            s = f"{s}_{lig_tag[lig]}"
        seen[_safe(lig)] = seen.get(_safe(lig), 0) + 1
        lig_stem[lig] = s

    # ordered (ligand, cluster) groups + a distinct colour each (for by_cluster)
    cluster_order: List[Tuple[str, int]] = []
    for lig in ligands:
        for (_, cl) in sorted([k for k in cluster_resi if k[0] == lig],
                              key=lambda k: k[1]):
            cluster_order.append((lig, cl))
    cluster_col = dict(zip(cluster_order, _distinct_hex(len(cluster_order))))

    lines: List[str] = []
    A = lines.append
    A(f"# PyMOL scene — {dataset} docked poses & clusters on Orai frame {frame}")
    A(f"# Generated by Scripts/Analysis/pymol_pose_cluster_scene.py")
    A("# Open with:  pymol " + f"{stem}__scene.pml")
    A("#")
    A("# Recolour commands (type in the PyMOL prompt):")
    A("#   by_tool     colour poses by docking tool chain")
    A("#   by_ligand   colour poses by ligand")
    A("#   by_cluster  give every (ligand,cluster) group its own colour")
    A("# (tempFactor still holds the per-ligand cluster id, so `spectrum b, rainbow,")
    A("#  poses and segi L01` colours one ligand's clusters if you prefer a gradient.)")
    A("")
    A("reinitialize")
    A("bg_color white")
    A("set ray_shadows, 0")
    A("set orthoscopic, 1")
    A("")
    # receptor
    if receptor_pdb and receptor_pdb.exists():
        A(f"load {receptor_pdb.resolve()}, receptor")
        A("hide everything, receptor")
        A("show cartoon, receptor")
        A("color grey80, receptor")
        A("set cartoon_transparency, 0.55, receptor")
        A("set cartoon_side_chain_helper, 1")
    else:
        A(f"# receptor PDB not found for frame {frame}; load it yourself, e.g.:")
        A(f"#   load Data/Receptors/{frame}.pdb, receptor")
    A("")
    # poses  (pml sits next to the poses PDB -> load by basename)
    A(f"load {pdb_name}, poses")
    A("hide everything, poses")
    A("show sticks, poses")
    A("set stick_radius, 0.14, poses")
    A("hide everything, hydro")              # hydrogens hidden by default (reversible)
    A("set valence, 0")
    A("")
    # ── tool selection group ──────────────────────────────────────────────
    A("# ---- Tools (chain = tool) ----")
    for t in tools:
        sel = f"tool_{_safe(TOOL_STYLE.get(t, (t.title(),))[0])}"
        A(f"select {sel}, poses and chain {chain_map[t]}")
    A("group Tools, tool_*")
    A("")
    # ── ligand selection group ────────────────────────────────────────────
    A("# ---- Ligands (segID = ligand) ----")
    for lig in ligands:
        A(f"select lig_{lig_stem[lig]}, poses and segi {lig_tag[lig]}")
    A("group Ligands, lig_*")
    A("")
    # ── cluster selection group + labelled centres ────────────────────────
    A("# ---- Clusters (per ligand; resSeq lists) ----")
    for (lig, cl) in cluster_order:
        resi = "+".join(str(r) for r in cluster_resi[(lig, cl)])
        A(f"select cl_{lig_stem[lig]}_c{cl}, poses and resi {resi}")
    A("group Clusters, cl_*")
    A("")
    A("# ---- Cluster centres (labelled pseudo-atoms) ----")
    for (lig, cl) in cluster_order:
        pts = cluster_pts[(lig, cl)]
        c = np.mean(pts, axis=0)
        dropped = cluster_dropped.get((lig, cl), 0)
        extra = f" (+{dropped} off-scene)" if dropped else ""
        lbl = f"{ligand_ccd(lig)} c{cl} n={len(pts)}{extra}"
        A(f"pseudoatom ctr_{lig_stem[lig]}_c{cl}, "
          f"pos=[{c[0]:.3f}, {c[1]:.3f}, {c[2]:.3f}], label=\"{lbl}\"")
    A("group ClusterCentres, ctr_*")
    A("hide everything, ctr_*")
    A("show nonbonded, ctr_*")
    A("set label_size, 14")
    A("set label_color, black")
    A("")
    # ── recolour commands ─────────────────────────────────────────────────
    tool_cmds = []
    for t in tools:
        col = _hex_pymol(TOOL_STYLE.get(t, (None, _TOOL_FALLBACK[0]))[1]
                         if t in TOOL_STYLE else
                         _TOOL_FALLBACK[hash(t) % len(_TOOL_FALLBACK)])
        tool_cmds.append(f"color {col}, poses and chain {chain_map[t]}")
    A("alias by_tool, " + "; ".join(tool_cmds))
    lig_cmds = [f"color {_hex_pymol(lig_color[lig])}, poses and segi {lig_tag[lig]}"
                for lig in ligands]
    A("alias by_ligand, " + "; ".join(lig_cmds))
    # explicit per-(ligand,cluster) colouring so distinct spatial clusters never
    # share a colour across ligands (unlike a single global `spectrum b`)
    clu_cmds = [f"color {_hex_pymol(cluster_col[(lig, cl)])}, cl_{lig_stem[lig]}_c{cl}"
                for (lig, cl) in cluster_order]
    A("alias by_cluster, " + "; ".join(clu_cmds))
    A("")
    A("# initial view: colour by tool")
    A("by_tool")
    A("deselect")
    A("orient poses")
    A("zoom poses, 6")
    A("# for the whole channel in context instead, run:  zoom receptor")
    A("")
    return "\n".join(lines) + "\n"


# ════════════════════════════════════════════════════════════════════════
# Loading / filtering
# ════════════════════════════════════════════════════════════════════════
def load_per_pose(csv: Path) -> pd.DataFrame:
    df = pd.read_csv(csv, low_memory=False)
    need = {"frame", "ligand", "tool", "cluster", "pose_file"}
    missing = need - set(df.columns)
    if missing:
        sys.exit(f"ERROR: {csv} is missing required columns: {sorted(missing)}\n"
                 "Expected the per_pose.csv written by orai_pose_cluster_report.py.")
    if "pb_valid" not in df.columns:
        df["pb_valid"] = np.nan
    df["pb_valid"] = _normalize_pb_valid(df["pb_valid"])
    df["tool"] = df["tool"].astype(str).str.lower()
    df["frame"] = df["frame"].astype(str)
    df["ligand"] = df["ligand"].astype(str)
    return df


def _normalize_pb_valid(col: pd.Series) -> pd.Series:
    """Coerce a pb_valid column (bool / 0-1 / 'True'/'False' strings / NaN) to a
    canonical object column of {True, False, np.nan} so occupancy is never silently
    collapsed to 'unknown' just because pandas read the column as strings."""
    truthy = {"true", "1", "1.0", "yes", "t"}
    falsy = {"false", "0", "0.0", "no", "f"}

    def _one(v):
        if v is True or v is False:
            return v
        if isinstance(v, float) and np.isnan(v):
            return np.nan
        s = str(v).strip().lower()
        if s in truthy:
            return True
        if s in falsy:
            return False
        if s in ("", "nan", "none"):
            return np.nan
        return np.nan

    return col.map(_one).astype(object)


def select_ligands(df_frame: pd.DataFrame, requested: Optional[List[str]],
                   max_ligands: int) -> Tuple[List[str], bool]:
    """Return (ligands to draw, capped?). Hard cap at 99 so the L01.. segID tags
    (4-char PDB segID field) stay unique regardless of --max-ligands / --ligands."""
    _HARD = 99
    counts = df_frame.groupby("ligand").size().sort_values(ascending=False)
    all_ligs = counts.index.tolist()
    if requested:
        chosen = [l for l in requested if l in set(all_ligs)]
        if len(chosen) > _HARD:
            print(f"  NOTE: {len(chosen)} ligands requested; capping to the "
                  f"{_HARD} busiest (segID limit).")
            order = {l: i for i, l in enumerate(all_ligs)}
            chosen = sorted(chosen, key=lambda l: order.get(l, 1 << 30))[:_HARD]
            return chosen, True
        return chosen, False
    cap = min(max_ligands, _HARD)
    if len(all_ligs) <= cap:
        return all_ligs, False
    return all_ligs[:cap], True


# ════════════════════════════════════════════════════════════════════════
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--per-pose-csv",
                    default="posebusters_results/orai_benchmark/pose_clusters/per_pose.csv",
                    help="Clustered per_pose.csv from orai_pose_cluster_report.py.")
    ap.add_argument("--receptor-dir", default="Data/Receptors",
                    help="Directory holding <frame>.pdb receptor files.")
    ap.add_argument("--out-dir",
                    default="posebusters_results/orai_benchmark/pymol_scene")
    ap.add_argument("--frames", default="",
                    help="Comma-separated frames to render (default: all present).")
    ap.add_argument("--ligands", default="",
                    help="Comma-separated ligands (default: all, capped by --max-ligands).")
    ap.add_argument("--max-ligands", type=int, default=12,
                    help="Cap ligands per scene when --ligands not given (busiest kept).")
    ap.add_argument("--max-offset", type=float, default=150.0,
                    help="Drop poses whose centroid is >this many A from the receptor "
                         "centroid (corrupt / runaway coordinates).")
    ap.add_argument("--max-extent", type=float, default=30.0,
                    help="Drop poses whose atom->centroid radius exceeds this many A "
                         "(catastrophically distorted geometry, e.g. exploded EquiBind "
                         "poses that render as box-spanning webs). Real poses are <~12 A.")
    ap.add_argument("--dataset-label", default="",
                    help="Label used in filenames/titles (default: inferred from CSV path).")
    ap.add_argument("--render", action="store_true",
                    help="Also ray-trace a <stem>__preview.png of each scene "
                         "(needs PyMOL; best-effort, coloured by tool).")
    args = ap.parse_args(argv)

    csv = Path(args.per_pose_csv)
    if not csv.exists():
        sys.exit(f"ERROR: per-pose CSV not found: {csv}")
    df = load_per_pose(csv)

    dataset = args.dataset_label or _infer_dataset(csv)
    receptor_dir = Path(args.receptor_dir)
    out_dir = Path(args.out_dir)

    frames = ([f.strip() for f in args.frames.split(",") if f.strip()]
              or sorted(df["frame"].unique().tolist()))
    req_ligs = [l.strip() for l in args.ligands.split(",") if l.strip()] or None

    print(f"Loaded {len(df)} pose rows from {csv}")
    print(f"Dataset label: {dataset} | frames: {frames}")

    summaries = []
    for frame in frames:
        df_f = df[df["frame"] == frame].copy()
        if df_f.empty:
            print(f"  [{frame}] no poses in CSV — skipped")
            continue
        ligs, capped = select_ligands(df_f, req_ligs, args.max_ligands)
        if not ligs:
            print(f"  [{frame}] none of the requested ligands present — skipped")
            continue
        if capped:
            print(f"  [{frame}] {df_f['ligand'].nunique()} ligands present; "
                  f"drawing the {len(ligs)} busiest (--max-ligands {args.max_ligands}). "
                  f"Pass --ligands to choose: {', '.join(ligs)}")
        df_f = df_f[df_f["ligand"].isin(ligs)].copy()
        # keep the busiest-first ligand order for stable colour/segID assignment
        order = {l: i for i, l in enumerate(ligs)}
        df_f["_ord"] = df_f["ligand"].map(order)
        df_f = df_f.sort_values(["_ord", "cluster"]).drop(columns="_ord")

        receptor_pdb = receptor_dir / f"{frame}.pdb"
        s = build_scene(df_f, frame, receptor_pdb, out_dir, dataset, args.max_offset,
                        args.max_extent)
        if s:
            summaries.append(s)

    if not summaries:
        print("No scenes written.")
        return 1
    print(f"\nWrote {len(summaries)} scene(s) to {out_dir}/")
    print("Open a scene with, e.g.:")
    print(f"    pymol {summaries[0]['pml']}")
    if args.render:
        print("Rendering preview PNGs ...")
        pngs = render_previews(summaries)
        for p in pngs:
            print(f"    preview: {p}")
    return 0


def _infer_dataset(csv: Path) -> str:
    for part in csv.parts:
        if part.startswith("orai_"):
            return part
    return "orai"


def render_previews(summaries: List[dict], size=(1200, 900)) -> List[Path]:
    """Best-effort headless-PyMOL render of each scene to <stem>__preview.png.

    Returns the PNG paths written. Silently degrades (prints a note) if PyMOL is
    unavailable or rendering fails — the .pml/.pdb are the primary deliverable.
    """
    import os
    try:
        import pymol
        pymol.finish_launching(["pymol", "-qc"])   # quiet, no GUI
        from pymol import cmd
    except Exception as exc:                        # noqa: BLE001
        print(f"  (preview render skipped — PyMOL unavailable: {exc})")
        return []
    # Resolve every path against the ORIGINAL cwd now — cmd.cd() (used below so the
    # .pml can load poses.pdb by basename) mutates the process cwd, which would
    # otherwise corrupt relative paths on later iterations.
    origin = Path(os.getcwd())
    jobs = []
    for s in summaries:
        pml_abs = (origin / s["pml"]).resolve()
        png_abs = pml_abs.with_name(pml_abs.name.replace("__scene.pml", "__preview.png"))
        jobs.append((pml_abs, png_abs))
    pngs: List[Path] = []
    for pml_abs, png in jobs:
        try:
            cmd.reinitialize()
            cmd.cd(str(pml_abs.parent))             # .pml loads poses.pdb by basename
            cmd.do(f"@{pml_abs.name}")
            cmd.set("ray_opaque_background", 1)
            cmd.bg_color("white")
            cmd.show("nb_spheres", "ctr_*")
            cmd.set("nonbonded_size", 0.5)
            cmd.ray(size[0], size[1])
            cmd.png(str(png), dpi=120)
            if png.exists():
                pngs.append(png)
        except Exception as exc:                    # noqa: BLE001
            print(f"  (preview render failed for {pml_abs.name}: {exc})")
    return pngs


if __name__ == "__main__":
    raise SystemExit(main())
