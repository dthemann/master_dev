# ============================================================================
# PyMOL Script to Reproduce PandaMap Interactions
# Complex: mgl_Orai1WT-MDSnap-Fr300__gsk7975a-prot-OPT_pose1
# ============================================================================

# Clean slate
delete all
reinitialize

# Load the combined complex file (contains both protein and ligand)
load //wsl.localhost/Ubuntu/home/manndo/MasterProject/pandamap_results/autodock/mgl_Orai1WT-MDSnap-Fr300__gsk7975a-prot-OPT/mgl_Orai1WT-MDSnap-Fr300__gsk7975a-prot-OPT_pose1_complex.pdb, complex

# ============================================================================
# BASIC DISPLAY SETTINGS
# ============================================================================
bg_color white
set antialias, 2
set ray_shadows, 0
hide everything

# ============================================================================
# SELECT PROTEIN AND LIGAND
# ============================================================================
select protein, complex and polymer.protein
select ligand, complex and resn LIG

# ============================================================================
# SHOW BINDING SITE (15Å around ligand)
# ============================================================================
select binding_region, protein within 15 of ligand
show cartoon, binding_region
color gray90, binding_region
set cartoon_transparency, 0.4

# ============================================================================
# SHOW LIGAND AS STICKS
# ============================================================================
show sticks, ligand
set stick_radius, 0.20, ligand
color atomic, ligand
color gray10, ligand and elem C

# ============================================================================
# SELECT INTERACTING RESIDUES (from PandaMap analysis)
# Residues within 5Å of ligand, classified by interaction type
# ============================================================================

# H-bond/Polar residues (GLN489, GLN809)
select hbond_key, protein and chain A and resi 489+809

# Hydrophobic residues (VAL488, ALA492, LEU500, LEU501, ALA503, LEU575, VAL578, PRO582, LEU806)
select hydrophobic_key, protein and chain A and resi 488+492+500+501+503+575+578+582+806

# Charged residues (ASP491, ASP493, ASP495, LYS579, LYS802, LYS807)
select charged_key, protein and chain A and resi 491+493+495+579+802+807

# Aromatic residues (HIS494, TYR496, PHE504)
select aromatic_key, protein and chain A and resi 494+496+504

# ============================================================================
# SHOW AND COLOR SIDE CHAINS
# ============================================================================
show sticks, hbond_key or hydrophobic_key or charged_key or aromatic_key
set stick_radius, 0.15

# Color by interaction type (matching PandaMap color scheme)
color green, hbond_key and elem C        # Green for H-bond donors/acceptors
color cyan, hydrophobic_key and elem C   # Cyan for hydrophobic contacts
color magenta, charged_key and elem C    # Magenta for charged/ionic interactions
color yellow, aromatic_key and elem C    # Yellow for aromatic (pi-stacking)

# ============================================================================
# HYDROGEN BONDS
# Show H-bonds between ligand and polar/charged residues
# ============================================================================
distance hbonds, ligand and (elem N+O), (hbond_key or charged_key) and (elem N+O), 3.5, 2
color green, hbonds
set dash_width, 4, hbonds
set dash_gap, 0.15, hbonds
set dash_radius, 0.12, hbonds
hide labels, hbonds

# ============================================================================
# HYDROPHOBIC CONTACTS
# Show contacts between ligand carbons and hydrophobic residues
# ============================================================================
distance hydrophobic_contacts, ligand and elem C, hydrophobic_key and elem C, 4.5, 0
color cyan, hydrophobic_contacts
set dash_width, 2, hydrophobic_contacts
set dash_gap, 0.3, hydrophobic_contacts
set dash_radius, 0.08, hydrophobic_contacts
hide labels, hydrophobic_contacts

# ============================================================================
# PI-STACKING / AROMATIC INTERACTIONS
# Show interactions with aromatic residues
# ============================================================================
distance pi_interactions, ligand, aromatic_key, 5.5, 0
color orange, pi_interactions
set dash_width, 3, pi_interactions
set dash_gap, 0.2, pi_interactions
set dash_radius, 0.10, pi_interactions
hide labels, pi_interactions

# ============================================================================
# CATION-PI / IONIC INTERACTIONS
# Show interactions with charged residues
# ============================================================================
distance ionic_contacts, ligand, charged_key and (resn LYS and name NZ or resn ARG and name NH* or resn ASP+GLU and name OE*+OD*), 5.0, 0
color deeppurple, ionic_contacts
set dash_width, 3.5, ionic_contacts
set dash_gap, 0.18, ionic_contacts
set dash_radius, 0.10, ionic_contacts
hide labels, ionic_contacts

# ============================================================================
# LABELS (matching PandaMap style)
# ============================================================================
label (hbond_key or hydrophobic_key or charged_key or aromatic_key) and name CA, "%s %s" % (resn, resi)
set label_size, 14
set label_font_id, 7
set label_color, black
set label_bg_color, white
set label_bg_transparency, 0.3
set label_position, (1.5, 1.5, 2)

# ============================================================================
# CAMERA / VIEW SETTINGS
# ============================================================================
zoom ligand, 8
center ligand
orient ligand

# Fine-tune view angle for better visualization
turn y, 25
turn x, -15

# ============================================================================
# LEGEND (create as CGO objects)
# ============================================================================
# Create pseudo-atoms for legend reference
# pseudoatom legend_hbond, pos=[-25, 15, 0], color=green, label="H-bond/Polar"
# pseudoatom legend_hydro, pos=[-25, 13, 0], color=cyan, label="Hydrophobic"
# pseudoatom legend_charged, pos=[-25, 11, 0], color=magenta, label="Charged"
# pseudoatom legend_aromatic, pos=[-25, 9, 0], color=yellow, label="Aromatic"

# ============================================================================
# HIGH QUALITY RENDERING
# ============================================================================
# Uncomment below to render high-quality image:
# set ray_trace_mode, 1
# set ray_opaque_background, on
# ray 2400, 2400
# png pandamap_pymol_interactions.png, dpi=300

# ============================================================================
# CLEANUP - deselect all
# ============================================================================
deselect

print("PandaMap interactions visualization loaded!")
print("Residues shown:")
print("  Green (H-bond/Polar): GLN489, GLN809")
print("  Cyan (Hydrophobic): VAL488, ALA492, LEU500, LEU501, ALA503, LEU575, VAL578, PRO582, LEU806")
print("  Magenta (Charged): ASP491, ASP493, ASP495, LYS579, LYS802, LYS807")
print("  Yellow (Aromatic): HIS494, TYR496, PHE504")
