#!/usr/bin/env python3
"""
Generate PyMOL visualization scripts from PandaMap analysis results.
This script reads a PandaMap complex PDB and creates a PyMOL script
that visualizes all protein-ligand interactions.

Uses the EXACT same colors and label format as PandaMap 2D visualization.

Usage:
    python generate_pandamap_pymol.py <complex_pdb> [output_pml]
"""

import sys
import math
from pathlib import Path
from typing import Dict, List, Tuple, Set


# ============================================================================
# PANDAMAP INTERACTION STYLES (EXACT MATCH from pandamap/core.py)
# ============================================================================
PANDAMAP_INTERACTION_STYLES = {
    'hydrogen_bonds': {
        'color': 'green',           # PandaMap: 'green'
        'pymol_color': 'green',
        'linestyle': '-',
        'marker_text': 'H',
        'name': 'Hydrogen Bond'
    },
    'carbon_pi': {
        'color': '#666666',         # PandaMap: '#666666' (gray)
        'pymol_color': '0x666666',
        'linestyle': '--',
        'marker_text': 'C-π',
        'name': 'Carbon-Pi'
    },
    'pi_pi_stacking': {
        'color': '#9370DB',         # PandaMap: '#9370DB' (medium purple)
        'pymol_color': '0x9370DB',
        'linestyle': '--',
        'marker_text': 'π-π',
        'name': 'Pi-Pi'
    },
    'donor_pi': {
        'color': '#FF69B4',         # PandaMap: '#FF69B4' (hot pink)
        'pymol_color': '0xFF69B4',
        'linestyle': '--',
        'marker_text': 'D',
        'name': 'Donor-Pi'
    },
    'amide_pi': {
        'color': '#A52A2A',         # PandaMap: '#A52A2A' (brown)
        'pymol_color': '0xA52A2A',
        'linestyle': '--',
        'marker_text': 'A',
        'name': 'Amide-Pi'
    },
    'hydrophobic': {
        'color': '#808080',         # PandaMap: '#808080' (gray)
        'pymol_color': '0x808080',
        'linestyle': ':',
        'marker_text': 'h',
        'name': 'Hydrophobic'
    },
    'ionic': {
        'color': '#FF4500',         # PandaMap: '#FF4500' (orange-red)
        'pymol_color': '0xFF4500',
        'linestyle': '-',
        'marker_text': 'I',
        'name': 'Ionic'
    },
    'halogen_bonds': {
        'color': '#00CED1',         # PandaMap: '#00CED1' (dark turquoise)
        'pymol_color': '0x00CED1',
        'linestyle': '-',
        'marker_text': 'X',
        'name': 'Halogen Bond'
    },
    'cation_pi': {
        'color': '#FF00FF',         # PandaMap: '#FF00FF' (magenta)
        'pymol_color': 'magenta',
        'linestyle': '--',
        'marker_text': 'C+π',
        'name': 'Cation-Pi'
    },
    'metal_coordination': {
        'color': '#FFD700',         # PandaMap: '#FFD700' (gold)
        'pymol_color': '0xFFD700',
        'linestyle': '-',
        'marker_text': 'M',
        'name': 'Metal Coordination'
    },
    'salt_bridge': {
        'color': '#FF6347',         # PandaMap: '#FF6347' (tomato)
        'pymol_color': '0xFF6347',
        'linestyle': '-',
        'marker_text': 'S',
        'name': 'Salt Bridge'
    },
    'covalent': {
        'color': '#000000',         # PandaMap: '#000000' (black)
        'pymol_color': 'black',
        'linestyle': '-',
        'marker_text': 'COV',
        'name': 'Covalent Bond'
    },
    'alkyl_pi': {
        'color': '#4682B4',         # PandaMap: '#4682B4' (steel blue)
        'pymol_color': '0x4682B4',
        'linestyle': '--',
        'marker_text': 'A-π',
        'name': 'Alkyl-Pi'
    },
    'attractive_charge': {
        'color': '#1E90FF',         # PandaMap: '#1E90FF' (dodger blue)
        'pymol_color': '0x1E90FF',
        'linestyle': '-',
        'marker_text': 'A+',
        'name': 'Attractive Charge'
    },
    'pi_cation': {
        'color': '#FF00FF',         # PandaMap: '#FF00FF' (magenta)
        'pymol_color': 'magenta',
        'linestyle': '--',
        'marker_text': 'π-C+',
        'name': 'Pi-Cation'
    },
    'repulsion': {
        'color': '#DC143C',         # PandaMap: '#DC143C' (crimson)
        'pymol_color': '0xDC143C',
        'linestyle': '-',
        'marker_text': 'R',
        'name': 'Repulsion'
    }
}


# ============================================================================
# RESIDUE NAME MAPPINGS
# ============================================================================
# Map non-standard residue names to standard 3-letter codes
# This handles different force field naming conventions (CHARMM, AMBER, etc.)
RESNAME_MAP = {
    # Histidine variants
    'HSD': 'HIS',  # CHARMM: delta-protonated
    'HSE': 'HIS',  # CHARMM: epsilon-protonated  
    'HSP': 'HIS',  # CHARMM: doubly protonated (charged)
    'HID': 'HIS',  # AMBER: delta-protonated
    'HIE': 'HIS',  # AMBER: epsilon-protonated
    'HIP': 'HIS',  # AMBER: doubly protonated (charged)
    # Cysteine variants
    'CYX': 'CYS',  # Disulfide-bonded cysteine
    'CYM': 'CYS',  # Deprotonated cysteine
    # Aspartate variants
    'ASH': 'ASP',  # Protonated aspartate
    # Glutamate variants
    'GLH': 'GLU',  # Protonated glutamate
    # Lysine variants
    'LYN': 'LYS',  # Neutral lysine
    # Terminal residues
    'NTER': None,  # N-terminal cap
    'CTER': None,  # C-terminal cap
    'ACE': None,   # Acetyl cap
    'NME': None,   # N-methyl cap
}

def standardize_resname(resname: str) -> str:
    """Convert non-standard residue names to standard 3-letter codes."""
    return RESNAME_MAP.get(resname, resname)

def get_resname_variants(resname: str) -> List[str]:
    """Get all possible variants of a residue name for PyMOL selection."""
    variants = [resname]
    # Add reverse mappings
    for variant, standard in RESNAME_MAP.items():
        if standard == resname:
            variants.append(variant)
    return list(set(variants))


# ============================================================================
# AMINO ACID CLASSIFICATIONS (matching PandaMap detection logic)
# ============================================================================
AROMATIC_RESIDUES = {'PHE', 'TYR', 'TRP', 'HIS'}
HBOND_DONORS = {'ARG', 'LYS', 'HIS', 'ASN', 'GLN', 'SER', 'THR', 'TYR', 'TRP'}
HBOND_ACCEPTORS = {'ASP', 'GLU', 'ASN', 'GLN', 'HIS', 'SER', 'THR', 'TYR'}
NEG_CHARGED = {'ASP', 'GLU'}
POS_CHARGED = {'ARG', 'LYS', 'HIS'}
HYDROPHOBIC_RESIDUES = {'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'TRP', 'PRO', 'TYR'}
ALKYL_RESIDUES = {'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PRO'}
AMIDE_RESIDUES = {'ASN', 'GLN'}

# Extended sets including variants
POLAR_AA = {'SER', 'THR', 'ASN', 'GLN', 'HIS', 'TYR', 'CYS'}
CHARGED_AA = {'ASP', 'GLU', 'LYS', 'ARG'}
AROMATIC_AA = {'PHE', 'TYR', 'TRP', 'HIS'}
HYDROPHOBIC_AA = {'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PRO'}

POLAR_AA_ALL = POLAR_AA | {'HSD', 'HSE', 'HSP', 'HID', 'HIE', 'HIP', 'CYX', 'CYM'}
CHARGED_AA_ALL = CHARGED_AA | {'ASH', 'GLH', 'LYN', 'HSP', 'HIP'}
AROMATIC_AA_ALL = AROMATIC_AA | {'HSD', 'HSE', 'HSP', 'HID', 'HIE', 'HIP'}


def infer_interaction_type(resname: str, distance: float) -> str:
    """
    Infer the most likely interaction type based on residue type and distance.
    This attempts to match PandaMap's classification logic.
    """
    std_res = standardize_resname(resname)
    
    # Very close = likely hydrogen bond or ionic
    if distance < 3.5:
        if std_res in NEG_CHARGED or std_res in POS_CHARGED:
            return 'ionic'
        elif std_res in HBOND_DONORS or std_res in HBOND_ACCEPTORS:
            return 'hydrogen_bonds'
    
    # Medium range = could be pi-stacking, hydrophobic, or h-bonds
    if distance < 4.5:
        if std_res in AROMATIC_RESIDUES:
            return 'pi_pi_stacking'
        elif std_res in NEG_CHARGED or std_res in POS_CHARGED:
            return 'salt_bridge'
        elif std_res in HYDROPHOBIC_RESIDUES:
            return 'hydrophobic'
        elif std_res in HBOND_DONORS or std_res in HBOND_ACCEPTORS:
            return 'hydrogen_bonds'
    
    # Longer range = hydrophobic or alkyl-pi
    if std_res in ALKYL_RESIDUES:
        return 'alkyl_pi' if distance < 5.0 else 'hydrophobic'
    elif std_res in AROMATIC_RESIDUES:
        return 'carbon_pi'
    elif std_res in HYDROPHOBIC_RESIDUES:
        return 'hydrophobic'
    
    # Default
    return 'hydrophobic'


def read_complex_pdb(pdb_file: str) -> Tuple[List[Tuple], List[Tuple]]:
    """Read protein atoms and ligand atoms from complex PDB."""
    protein_atoms = []
    ligand_atoms = []
    
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                resname = line[17:20].strip()
                chain = line[21].strip() or 'A'
                resnum = int(line[22:26])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                atom_name = line[12:16].strip()
                element = line[76:78].strip() if len(line) > 76 else atom_name[0]
                protein_atoms.append((chain, resnum, resname, atom_name, element, x, y, z))
            elif line.startswith('HETATM'):
                resname = line[17:20].strip()
                if resname == 'LIG':
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    atom_name = line[12:16].strip()
                    element = line[76:78].strip() if len(line) > 76 else atom_name[0]
                    ligand_atoms.append((x, y, z, atom_name, element))
    
    return protein_atoms, ligand_atoms


def find_nearby_residues(protein_atoms: List[Tuple], ligand_atoms: List[Tuple], 
                          cutoff: float = 5.0) -> Dict[Tuple, float]:
    """Find protein residues within cutoff distance of ligand."""
    nearby = {}
    
    for chain, resnum, resname, atom_name, element, px, py, pz in protein_atoms:
        for lx, ly, lz, _, _ in ligand_atoms:
            dist = math.sqrt((px-lx)**2 + (py-ly)**2 + (pz-lz)**2)
            if dist < cutoff:
                key = (chain, resnum, resname)
                if key not in nearby:
                    nearby[key] = dist
                else:
                    nearby[key] = min(nearby[key], dist)
    
    return nearby


def classify_residues_by_interaction(nearby: Dict[Tuple, float]) -> Dict[str, List[Tuple]]:
    """Classify residues by PandaMap interaction type."""
    interactions_by_type = {itype: [] for itype in PANDAMAP_INTERACTION_STYLES.keys()}
    
    for (chain, resnum, resname), dist in nearby.items():
        std_resname = standardize_resname(resname)
        itype = infer_interaction_type(resname, dist)
        interactions_by_type[itype].append((chain, resnum, resname, std_resname, dist))
    
    # Sort by residue number and remove duplicates
    for itype in interactions_by_type:
        seen = set()
        unique = []
        for item in sorted(interactions_by_type[itype], key=lambda x: x[1]):
            if item[1] not in seen:
                seen.add(item[1])
                unique.append(item)
        interactions_by_type[itype] = unique
    
    return interactions_by_type


def generate_pymol_script(pdb_file: str, interactions_by_type: Dict[str, List[Tuple]], 
                          nearby: Dict[Tuple, float]) -> str:
    """Generate PyMOL script content with PandaMap-synchronized colors."""
    
    pdb_path = Path(pdb_file).resolve()
    complex_name = pdb_path.stem
    
    # Get the main chain
    main_chain = 'A'
    if nearby:
        main_chain = list(nearby.keys())[0][0] or 'A'
    
    # Build script
    script_lines = [
        f'''# ============================================================================
# PyMOL Script - PandaMap Interactions Visualization
# Complex: {complex_name}
# Generated with PandaMap-synchronized colors and labels
# ============================================================================

# Clean slate
delete all
reinitialize

# Load complex
load {pdb_path}, complex

# Global settings
bg_color white
set antialias, 2
set ray_shadows, 0
set depth_cue, 1
set fog_start, 0.4
hide everything

# ============================================================================
# PROTEIN AND LIGAND DISPLAY
# ============================================================================
select protein, complex and polymer.protein
select ligand, complex and resn LIG

# Show binding site region
select binding_region, protein within 12 of ligand
show cartoon, binding_region
color gray90, binding_region
set cartoon_transparency, 0.4, binding_region

# Ligand display
show sticks, ligand
set stick_radius, 0.20, ligand
color atomic, ligand
color gray20, ligand and elem C

# ============================================================================
# INTERACTION-SPECIFIC SELECTIONS AND COLORING
# ============================================================================
'''
    ]
    
    # Generate selections for each interaction type
    interaction_summary = []
    
    for itype, style in PANDAMAP_INTERACTION_STYLES.items():
        residues = interactions_by_type.get(itype, [])
        if not residues:
            continue
        
        resi_list = '+'.join(str(r[1]) for r in residues)
        label_list = ', '.join(f"{r[3]} {r[1]}" for r in residues)  # PandaMap format
        
        selection_name = f"int_{itype}"
        pymol_color = style['pymol_color']
        display_name = style['name']
        
        script_lines.append(f'''
# {display_name} interactions
select {selection_name}, protein and chain {main_chain} and resi {resi_list}
show sticks, {selection_name}
set stick_radius, 0.15, {selection_name}
color {pymol_color}, {selection_name} and elem C
''')
        
        interaction_summary.append(f"{display_name}: {label_list}")
    
    # Add distance objects
    script_lines.append('''
# ============================================================================
# INTERACTION DISTANCE OBJECTS
# ============================================================================
''')
    
    # Hydrogen bonds - green
    if interactions_by_type.get('hydrogen_bonds'):
        resi = '+'.join(str(r[1]) for r in interactions_by_type['hydrogen_bonds'])
        script_lines.append(f'''
# Hydrogen bonds (green, solid)
distance hbonds, ligand and (elem N+O), (protein and chain {main_chain} and resi {resi}) and (elem N+O), 3.5, 2
color green, hbonds
set dash_width, 3, hbonds
set dash_gap, 0.15, hbonds
set dash_radius, 0.08, hbonds
hide labels, hbonds
''')
    
    # Ionic/Salt bridge - orange-red/tomato
    ionic_resi = []
    if interactions_by_type.get('ionic'):
        ionic_resi.extend(str(r[1]) for r in interactions_by_type['ionic'])
    if interactions_by_type.get('salt_bridge'):
        ionic_resi.extend(str(r[1]) for r in interactions_by_type['salt_bridge'])
    if ionic_resi:
        resi = '+'.join(set(ionic_resi))
        script_lines.append(f'''
# Ionic/Salt bridge (orange-red)
distance ionic, ligand and (elem N+O), (protein and chain {main_chain} and resi {resi}) and (elem N+O), 4.0, 2
color 0xFF4500, ionic
set dash_width, 3, ionic
set dash_gap, 0.2, ionic
set dash_radius, 0.08, ionic
hide labels, ionic
''')
    
    # Hydrophobic - gray
    hydro_resi = []
    if interactions_by_type.get('hydrophobic'):
        hydro_resi.extend(str(r[1]) for r in interactions_by_type['hydrophobic'])
    if interactions_by_type.get('alkyl_pi'):
        hydro_resi.extend(str(r[1]) for r in interactions_by_type['alkyl_pi'])
    if hydro_resi:
        resi = '+'.join(set(hydro_resi))
        script_lines.append(f'''
# Hydrophobic contacts (gray, dotted)
distance hydrophobic, ligand and elem C, (protein and chain {main_chain} and resi {resi}) and elem C, 4.5, 0
color 0x808080, hydrophobic
set dash_width, 2, hydrophobic
set dash_gap, 0.4, hydrophobic
set dash_radius, 0.06, hydrophobic
hide labels, hydrophobic
''')
    
    # Pi-stacking - purple
    pi_resi = []
    for pi_type in ['pi_pi_stacking', 'carbon_pi', 'cation_pi', 'pi_cation', 'donor_pi', 'amide_pi']:
        if interactions_by_type.get(pi_type):
            pi_resi.extend(str(r[1]) for r in interactions_by_type[pi_type])
    if pi_resi:
        resi = '+'.join(set(pi_resi))
        script_lines.append(f'''
# Pi-stacking (purple, dashed)
distance pi_stack, ligand, protein and chain {main_chain} and resi {resi}, 5.5, 0
color 0x9370DB, pi_stack
set dash_width, 2.5, pi_stack
set dash_gap, 0.25, pi_stack
set dash_radius, 0.07, pi_stack
hide labels, pi_stack
''')
    
    # Halogen bonds - cyan
    if interactions_by_type.get('halogen_bonds'):
        resi = '+'.join(str(r[1]) for r in interactions_by_type['halogen_bonds'])
        script_lines.append(f'''
# Halogen bonds (cyan)
distance halogen, ligand and (elem CL+BR+I+F), (protein and chain {main_chain} and resi {resi}) and (elem N+O), 3.5, 2
color 0x00CED1, halogen
set dash_width, 3, halogen
set dash_gap, 0.15, halogen
hide labels, halogen
''')
    
    # Labels - PandaMap format: "RESNAME RESNUM"
    all_interacting_resi = set()
    for itype, residues in interactions_by_type.items():
        for r in residues:
            all_interacting_resi.add(r[1])
    
    if all_interacting_resi:
        resi_list = '+'.join(str(r) for r in sorted(all_interacting_resi))
        script_lines.append(f'''
# ============================================================================
# RESIDUE LABELS (PandaMap format: "RESNAME RESNUM")
# ============================================================================
select all_interacting, protein and chain {main_chain} and resi {resi_list}
label all_interacting and name CA, "%s %s" % (resn, resi)
set label_size, 14
set label_font_id, 7
set label_color, black
set label_bg_color, white
set label_bg_transparency, 0.3
set label_position, (1.5, 1.5, 2)
''')
    
    # Camera
    script_lines.append('''
# ============================================================================
# CAMERA AND FINAL SETTINGS
# ============================================================================
zoom ligand, 8
center ligand
orient ligand
turn y, 20
turn x, -10

deselect
''')
    
    # Interaction summary
    script_lines.append('\n# ============================================================================')
    script_lines.append('# INTERACTION SUMMARY (PandaMap-style labels)')
    script_lines.append('# ============================================================================')
    for summary in interaction_summary:
        script_lines.append(f'print("{summary}")')
    
    return '\n'.join(script_lines)


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        print("\nPandaMap Interaction Color Scheme:")
        for itype, style in PANDAMAP_INTERACTION_STYLES.items():
            print(f"  {style['name']:20s}: {style['color']}")
        sys.exit(1)
    
    pdb_file = sys.argv[1]
    output_pml = sys.argv[2] if len(sys.argv) > 2 else None
    
    if not Path(pdb_file).exists():
        print(f"Error: File not found: {pdb_file}")
        sys.exit(1)
    
    # Read and analyze
    print(f"Reading: {pdb_file}")
    protein_atoms, ligand_atoms = read_complex_pdb(pdb_file)
    print(f"  Protein atoms: {len(protein_atoms)}")
    print(f"  Ligand atoms: {len(ligand_atoms)}")
    
    nearby = find_nearby_residues(protein_atoms, ligand_atoms, cutoff=5.0)
    print(f"  Residues within 5Å: {len(nearby)}")
    
    interactions_by_type = classify_residues_by_interaction(nearby)
    
    print(f"\nClassified residues by PandaMap interaction type:")
    for itype, residues in interactions_by_type.items():
        if residues:
            style = PANDAMAP_INTERACTION_STYLES[itype]
            res_list = [f'{r[3]} {r[1]}' for r in residues]
            print(f"  {style['name']}: {', '.join(res_list)}")
    
    # Generate script
    script = generate_pymol_script(pdb_file, interactions_by_type, nearby)
    
    # Output
    if output_pml:
        output_path = Path(output_pml)
    else:
        output_path = Path(pdb_file).with_suffix('.pml')
    
    with open(output_path, 'w') as f:
        f.write(script)
    
    print(f"\nPyMOL script saved: {output_path}")
    print(f"Run in PyMOL with: @{output_path}")


if __name__ == "__main__":
    main()
