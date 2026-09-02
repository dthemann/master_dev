# -*- coding: utf-8 -*-
"""Build the FHTW-templated master-thesis defence deck."""
import copy
from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.oxml.ns import qn
from PIL import Image

TPL = "template.pptx"
OUT = "Mann_MasterThesis_Defence_v3.pptx"
MEDIA = "/home/manndo/master_dev/thesis_latex/media/media"

BLUE  = RGBColor(0x00, 0x64, 0x9C)   # FHTW accent5
GREY  = RGBColor(0x72, 0x77, 0x7A)   # FHTW accent1
DARK  = RGBColor(0x1A, 0x1A, 0x1A)
LIGHT = RGBColor(0xEC, 0xF1, 0xF5)

FOOTER = "Whole-Protein Docking with AutoDock Vina, DiffDock and EquiBind  |  Dominik Mann"

L_TITLE, L_ONE, L_TWO = 6, 7, 8      # layout indices in master 0

prs = Presentation(TPL)
LAYOUTS = prs.slide_masters[0].slide_layouts
SW, SH = prs.slide_width, prs.slide_height

# ---------------------------------------------------------------- housekeeping
for i in range(len(prs.slides) - 1, -1, -1):
    rId = prs.slides._sldIdLst[i].rId
    prs.part.drop_rel(rId)
    del prs.slides._sldIdLst[i]


def slidenum_field(tf):
    """Turn a slide-number placeholder into a live <slidenum> field."""
    p = tf.paragraphs[0]._p
    fld = p.makeelement(qn("a:fld"), {"id": "{B7B3D2A1-1C7E-4E4C-9B45-9A2E4C7F0001}",
                                      "type": "slidenum"})
    t = p.makeelement(qn("a:t"), {})
    t.text = "1"
    fld.append(t)
    p.append(fld)


def new_slide(layout_idx, title):
    s = prs.slides.add_slide(LAYOUTS[layout_idx])
    ph = {p.placeholder_format.idx: p for p in s.placeholders}
    if 0 in ph and title is not None:
        tf = ph[0].text_frame
        tf.text = title
        r = tf.paragraphs[0].runs[0]
        r.font.size = Pt(22)
        r.font.bold = True
        r.font.color.rgb = BLUE
    for lp in LAYOUTS[layout_idx].placeholders:
        kind = str(lp.placeholder_format.type)
        if not kind.startswith(("FOOTER", "SLIDE_NUMBER")):
            continue
        el = copy.deepcopy(lp._element)
        s.shapes._spTree.append(el)
        if kind.startswith("FOOTER"):
            new = s.shapes[-1]
            new.text_frame.paragraphs[0].runs[0].text = FOOTER
            new.text_frame.paragraphs[0].runs[0].font.size = Pt(8)
            new.text_frame.paragraphs[0].runs[0].font.color.rgb = GREY
    return s, ph


def drop(shape):
    shape._element.getparent().remove(shape._element)


def textbox(s, x, y, w, h):
    tb = s.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tf = tb.text_frame
    tf.word_wrap = True
    tf.margin_left = tf.margin_right = Emu(0)
    tf.margin_top = tf.margin_bottom = Emu(0)
    return tf


def para(tf, text, size=12, bold=False, colour=DARK, bullet=True,
         space_before=5, space_after=0, italic=False, first=False, indent=0):
    p = tf.paragraphs[0] if first else tf.add_paragraph()
    for e in p._p.findall(qn("a:r")) + p._p.findall(qn("a:br")):
        p._p.remove(e)
    p.space_before = Pt(space_before)
    p.space_after = Pt(space_after)
    p.level = indent
    pPr = p._p.get_or_add_pPr()
    for tag in ("a:buChar", "a:buAutoNum", "a:buNone"):
        for e in pPr.findall(qn(tag)):
            pPr.remove(e)
    if bullet:
        pPr.set("marL", str(int(Inches(0.20 + 0.20 * indent))))
        pPr.set("indent", str(-int(Inches(0.20))))
        bu = pPr.makeelement(qn("a:buChar"), {"char": "–" if indent else "▪"})
        buf = pPr.makeelement(qn("a:buFont"), {"typeface": "Arial"})
        pPr.append(buf)
        pPr.append(bu)
    else:
        pPr.set("marL", "0")
        pPr.set("indent", "0")
        pPr.append(pPr.makeelement(qn("a:buNone"), {}))
    r = p.add_run()
    r.text = text
    r.font.size = Pt(size)
    r.font.bold = bold
    r.font.italic = italic
    r.font.color.rgb = colour
    r.font.name = "Arial"
    return p


def bullets(tf, items, size=12, gap=6):
    for i, it in enumerate(items):
        if isinstance(it, tuple):
            txt, lvl = it
        else:
            txt, lvl = it, 0
        para(tf, txt, size=size - (1 if lvl else 0), first=(i == 0),
             space_before=0 if i == 0 else gap, indent=lvl,
             colour=DARK if lvl == 0 else GREY)


def col_header(tf, text, first=True):
    para(tf, text, size=13, bold=True, colour=BLUE, bullet=False,
         first=first, space_before=0, space_after=4)


def picture(s, path, x, y, max_w, max_h):
    """Place image centred in the (x, y, max_w, max_h) box, aspect preserved."""
    iw, ih = Image.open(path).size
    ar = iw / ih
    w, h = max_w, max_w / ar
    if h > max_h:
        h, w = max_h, max_h * ar
    left = x + (max_w - w) / 2
    top = y + (max_h - h) / 2
    return s.shapes.add_picture(path, Inches(left), Inches(top),
                                Inches(w), Inches(h))


def caveat(s, x, y, w, text, size=9, h=0.44):
    tf = textbox(s, x, y, w, h)
    para(tf, text, size=size, italic=True, colour=GREY, bullet=False,
         first=True, space_before=0)
    return tf


def notes(s, text):
    s.notes_slide.notes_text_frame.text = text.strip()


def chip(s, x, y, w, h, text, fill=LIGHT, line=BLUE, size=9.5, bold=True):
    from pptx.enum.shapes import MSO_SHAPE
    sh = s.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE,
                            Inches(x), Inches(y), Inches(w), Inches(h))
    sh.fill.solid()
    sh.fill.fore_color.rgb = fill
    sh.line.color.rgb = line
    sh.line.width = Pt(1)
    sh.shadow.inherit = False
    tf = sh.text_frame
    tf.word_wrap = True
    tf.margin_left = tf.margin_right = Emu(45720)
    tf.margin_top = tf.margin_bottom = Emu(18000)
    tf.vertical_anchor = MSO_ANCHOR.MIDDLE
    lines = text.split("\n")
    for i, ln in enumerate(lines):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.alignment = PP_ALIGN.CENTER
        p.space_before = Pt(0)
        p.space_after = Pt(0)
        r = p.add_run()
        r.text = ln
        r.font.size = Pt(size if i == 0 else size - 1)
        r.font.bold = bold if i == 0 else False
        r.font.color.rgb = DARK if i == 0 else GREY
        r.font.name = "Arial"
    return sh


def table(s, x, y, w, col_w, rows, hdr_h=0.26, row_h=0.30):
    """Two-row grouped header table with explicit cell formatting."""
    from pptx.enum.text import PP_ALIGN as _A
    n_r, n_c = len(rows), len(col_w)
    shp = s.shapes.add_table(n_r, n_c, Inches(x), Inches(y), Inches(w),
                             Inches(hdr_h * 2 + row_h * (n_r - 2)))
    tbl = shp.table
    tbl.first_row = False
    tbl.horz_banding = False
    for j, cw in enumerate(col_w):
        tbl.columns[j].width = Inches(cw)
    for i in range(n_r):
        tbl.rows[i].height = Inches(hdr_h if i < 2 else row_h)
    return tbl


def cell(tbl, r, c, text, size=9, bold=False, colour=DARK, fill=None,
         align="c"):
    from pptx.enum.text import PP_ALIGN as _A
    cl = tbl.cell(r, c)
    cl.margin_left = cl.margin_right = Emu(45720)
    cl.margin_top = cl.margin_bottom = Emu(9000)
    cl.vertical_anchor = MSO_ANCHOR.MIDDLE
    if fill is None:
        cl.fill.background()
    else:
        cl.fill.solid(); cl.fill.fore_color.rgb = fill
    tf = cl.text_frame
    tf.word_wrap = True
    p = tf.paragraphs[0]
    p.alignment = {"l": _A.LEFT, "c": _A.CENTER, "r": _A.RIGHT}[align]
    r0 = p.add_run(); r0.text = text
    r0.font.size = Pt(size); r0.font.bold = bold
    r0.font.color.rgb = colour; r0.font.name = "Arial"
    return cl


def arrow(s, x, y, w, h):
    from pptx.enum.shapes import MSO_SHAPE
    sh = s.shapes.add_shape(MSO_SHAPE.RIGHT_ARROW, Inches(x), Inches(y),
                            Inches(w), Inches(h))
    sh.fill.solid()
    sh.fill.fore_color.rgb = RGBColor(0xC5, 0xD3, 0xDD)
    sh.line.fill.background()
    sh.shadow.inherit = False
    sh.text_frame.text = ""
    return sh


# ================================================================== 1  TITLE
s, ph = new_slide(L_TITLE, None)
t = ph[0].text_frame
t.word_wrap = True
t.text = "Whole-Protein Docking with AutoDock Vina, DiffDock and EquiBind"
r = t.paragraphs[0].runs[0]
r.font.size = Pt(28); r.font.bold = True; r.font.color.rgb = BLUE
sub = ph[10]
sub.left, sub.top = Emu(252000), Emu(2659502)
sub.width, sub.height = Emu(8640000), Inches(1.05)
tf = sub.text_frame
tf.word_wrap = True
para(tf, "Physical validity and near-nativeness across ranking depth",
     size=14, bullet=False, first=True, space_before=0, colour=DARK)
para(tf, "Dominik Mann   |   MSc Artificial Intelligence Engineering   |   "
         "FH Technikum Wien   |   September 2026",
     size=11, bullet=False, colour=GREY, space_before=6)
notes(s, """
Master's thesis defence deck, 15 slides. Framing for an examiner from computer
science rather than structural biology: docking is treated throughout as a
retrieval-and-ranking problem over 3-D poses, with a physics-based search on one
side and two learned generative or regression models on the other.
Thesis title: Whole-Protein Docking with AutoDock Vina, DiffDock and EquiBind.
""")

# ======================================================= 2  MOTIVATION + RQs
s, ph = new_slide(L_TWO, "Motivation and Research Questions")
tf = ph[15].text_frame; tf.word_wrap = True
col_header(tf, "Why compare these three")
bullets(tf, [
    "Docking predicts where and how a small molecule binds a protein",
    "It is the cheap triage step before any laboratory work",
    "AI docking claims large gains in speed and accuracy",
    "Published comparisons judge a single nominated pose",
    "Search budgets are rarely matched across the competing tools",
    "Physical plausibility of the pose is often not checked at all",
    "No published comparison covers the Orai1 calcium channel",
], size=11.5, gap=6)

tf = ph[16].text_frame; tf.word_wrap = True
col_header(tf, "Research questions")
bullets(tf, [
    "RQ1  Physical validity, near-nativeness and ranking depth",
    "RQ2  What post-docking optimisation actually changes",
    "RQ3  Computational cost per usable result",
    "RQ4  Opportunities and limits of blind whole-protein docking",
], size=11.5, gap=9)
para(tf, "Scope", size=11.5, bold=True, colour=BLUE, bullet=False, space_before=14)
para(tf, "Each question is asked of three configured pipelines, not of "
         "“AI versus physics” as method classes. Every workflow couples a "
         "search or a model to a preparation, a refiner and a ranking rule.",
     size=10, colour=GREY, bullet=False, space_before=3)
notes(s, """
Likely examiner challenge: "isn't this just AI versus physics?" No. Each arm is a
pipeline, and the pipeline cannot be decomposed. The AutoDock arm is Vina search
at exhaustiveness 128 followed by a gnina CNN re-ranking pass, so it already
contains a learned component. The honest claim is a ranking of three configured
workflows on this benchmark.
Second challenge: why these three tools? Vina is the established baseline,
DiffDock is the leading diffusion model and EquiBind the reference one-shot
regression model. Newer co-folding systems were out of scope and are named in
the outlook.
""")

# ============================================================ 3  TOOLS + METRICS
s, ph = new_slide(L_ONE, "Tools, Pipeline and Metrics")
drop(ph[14])
BY, BH = 0.66, 0.72
xs = [0.20, 2.10, 4.00, 5.90, 7.80]
chip(s, xs[0], BY, 1.72, BH, "Ligand + receptor\nwhole protein, no pocket given")
chip(s, xs[1], BY, 1.72, BH, "Pose generation\nVina / DiffDock / EquiBind")
chip(s, xs[2], BY, 1.72, BH, "Local optimisation\nsmina or gnina")
chip(s, xs[3], BY, 1.72, BH, "PoseBusters filter\n22 physics checks")
chip(s, xs[4], BY, 1.72, BH, "Ranking + metrics\nrank-1 and top-k")
for x in (1.94, 3.84, 5.74, 7.64):
    arrow(s, x, BY + 0.28, 0.14, 0.16)

tf = textbox(s, 0.20, 1.62, 4.55, 3.35)
col_header(tf, "Three pipelines, one protocol")
bullets(tf, [
    "AutoDock Vina, stochastic search over an empirical score",
    "DiffDock, diffusion model with a learned confidence rank",
    "EquiBind, one forward pass, produces no score of its own",
    "smina and gnina relax every pose, gnina also re-ranks it",
    "None of the three receives a binding site, the box is the whole protein",
], size=11, gap=6)
para(tf, "* in every figure marks the carried-forward variant: Vina exh128 + gnina, "
         "DiffDock + smina, EquiBind + gnina.",
     size=9, italic=True, colour=GREY, bullet=False, space_before=8)

tf = textbox(s, 5.05, 1.62, 4.75, 3.35)
col_header(tf, "What is measured")
bullets(tf, [
    "Near-native, in-place RMSD ≤ 2 Å, position, orientation and shape",
    "Form, best-fit (Kabsch) RMSD ≤ 1 Å, shape only",
    "PB-valid, passes all 22 reference-free physics and chemistry checks",
    "Headline endpoint, PB-valid and ≤ 2 Å together",
    "Ranking depth, rank-1 against best-of-top-k, precision@1 against recall@k",
    "Cost, one charged basis, CPU-core-seconds / 32 threads + GPU-seconds",
], size=11, gap=6)
notes(s, """
The CS translation to offer out loud: a docking tool is a candidate generator plus
a ranker. RMSD is the distance to the labelled answer, PoseBusters is a hard
feasibility constraint, and rank-1 versus best-of-top-k is exactly precision@1
versus recall@k. Everything downstream is that distinction.
Anticipated question: why refine the AI outputs at all, does that not change what
is being measured? It does, and that is RQ2. Both raw and refined numbers are
reported, and the refined arm is the fair one because Vina minimises its poses
during its own search, so the raw comparison would be the unfair one.
Anticipated question about parity: the EquiBind arm generated its input conformers
with a different generator than the other two arms. This is disclosed in the
appendix and is a genuine residual non-parity.
""")

# ================================================================= 4  DATASETS
s, ph = new_slide(L_ONE, "Two Datasets: Calibrate, Then Apply")
drop(ph[14])
tf = textbox(s, 0.20, 0.66, 4.05, 4.30)
col_header(tf, "Calibration — PoseBusters benchmark")
bullets(tf, [
    "308 high-quality crystal complexes released after 2021",
    "303 analysed, DiffDock produced nothing for 5",
    "Every complex carries a ground-truth ligand pose",
    "Median ligand: 24 heavy atoms, 359 Da, 5 rotatable bonds",
    "Answers RQ1 to RQ3, because accuracy is measurable here",
    ("Post-2021 release reduces temporal overlap with training data "
     "but does not exclude homologous pockets", 1),
], size=11, gap=6)

tf = textbox(s, 4.55, 0.66, 3.15, 4.30)
col_header(tf, "Application — Orai1 channel")
bullets(tf, [
    "Calcium channel, no ligand-bound structure exists",
    "4 molecular-dynamics snapshots × 3 known modulators",
    "Control panel: the 308 benchmark ligands × 4 frames",
    "No RMSD is possible without a reference pose",
    "Endpoint is operational: valid, on the receptor, outside the pore slab",
    ("Effective independent observations: three ligands", 1),
], size=11, gap=6)

picture(s, f"{MEDIA}/image14.png", 7.85, 1.35, 1.95, 1.55)
caveat(s, 7.80, 2.98, 2.00,
       "Orai1 hexamer viewed down the pore. Yellow spheres mark the "
       "transmembrane exclusion slab.", size=8, h=0.75)
notes(s, """
The two datasets do different jobs. The benchmark has labels and therefore
supports accuracy claims. Orai1 has no label and therefore supports only an
operational claim about admissibility and placement.
Anticipated question: is the Orai1 arm underpowered? Yes, and it is reported that
way. Three ligands and four correlated frames of one trajectory. Every cross-panel
contrast is re-tested at ligand level, and almost none survives.
Anticipated question: is cognate re-docking on the benchmark harder or easier than
Orai1? Easier. Each benchmark complex is docked back into its own crystal
receptor, whereas Orai1 uses non-cognate snapshots of a homology model. The
benchmark numbers bound the Orai1 expectation from above.
""")

# ==================================================== 5  RQ2 OPTIMISATION
s, ph = new_slide(L_ONE, "One Optimisation Pass Makes Learned Poses Legal")
drop(ph[14])
tf = textbox(s, 0.20, 0.63, 9.60, 0.68)
bullets(tf, [
    "Raw output fails the physics checks almost everywhere, DiffDock at a 13% "
    "median and EquiBind at 0%",
    "A single local minimisation lifts them to 93% and 80%, gains of +62 and "
    "+58 percentage points",
    "AutoDock Vina has no headroom, it already minimises every pose during "
    "its own search",
], size=10.5, gap=3)
picture(s, f"{MEDIA}/image2.png", 0.55, 1.36, 8.90, 3.30)
caveat(s, 0.20, 4.72, 9.60,
       "Geometry is repaired, placement is not. Accuracy at a fixed rank moves "
       "by at most a few points, so optimisation cannot supply a binding mode "
       "that was never generated.")
notes(s, """
This is the answer to RQ2 and the single largest effect in the thesis, over 70
percentage points of pooled validity between Vina and the raw learned output.
Anticipated and important challenge: is this comparison circular? Partly, and the
thesis says so. The discriminating PoseBusters checks are intermolecular distance
checks, and they inspect coordinates that Vina has already minimised against a
closely related objective. The margin is therefore partly constitutive. The
mitigation is that requiring validity removes at most four of 303 complexes from
the headline endpoint, so it changes the ranking hardly at all.
Second challenge: raw versus rescored hydrogen handling. Both arms are rebuilt
from a bond-order template and PoseBusters adds the hydrogen shell itself, so the
two sides are measured on the same footing.
""")

# ================================ 6  ACCURACY: PLACEMENT AND FORM (merged)
s, ph = new_slide(L_ONE, "Placement and Form Separate the Tools Differently")
drop(ph[14])
for xh, lbl in ((0.20, "Position, orientation and shape   ·   in-place RMSD"),
                (5.10, "Shape only, after superposition   ·   best-fit RMSD")):
    tfh = textbox(s, xh, 0.64, 4.70, 0.24)
    para(tfh, lbl, size=10.5, bold=True, colour=BLUE, bullet=False,
         first=True, space_before=0)
picture(s, f"{MEDIA}/image3.png", 0.20, 0.90, 4.70, 3.34)
picture(s, f"{MEDIA}/image4.png", 5.10, 0.90, 4.70, 3.34)
tf = textbox(s, 0.20, 4.30, 9.60, 0.80)
bullets(tf, [
    "Placement at rank-1: AutoDock 36.6%, DiffDock 34.0%, EquiBind 18.2%, and "
    "the two leaders are not separable there",
    "Form at rank-1: 52.5%, 52.1% and 33.3%, a tie between the leaders at every "
    "depth, so the gap is rigid-body, not conformer quality",
    "Depth lifts both curves, but six in ten rank-1 failures hold no qualifying "
    "pose at any depth",
], size=9.5, gap=2)
para(tf, "Fairness caveat: AutoDock was tuned over five search budgets, the two "
         "learned arms only over a choice of refiner.",
     size=9, italic=True, colour=GREY, bullet=False, space_before=3)
notes(s, """
The two panels are the same 303 complexes and the same poses, scored on two axes.
Left is in-place RMSD, which mixes where the ligand sits with what shape it has.
Right is best-fit RMSD, which superimposes the two molecules first and therefore
reports the internal conformation alone.
Read them together. AutoDock leads on the left from top-5 onwards and ties on the
right everywhere, so the difference between the two leaders is rigid-body, not
conformer quality. DiffDock builds the right shape and seats it less well, which
is consistent with its generative treatment of translation, rotation and torsion.
Be precise about the word placement if it comes up. In-place RMSD is measured in
the crystal's frame with no superposition, so it carries translation, orientation
and conformation together. Best-fit RMSD superimposes first, so it carries
conformation alone. Placement is never measured directly. It is the residual, the
square root of in-place squared minus form squared, and that residual bundles
orientation with position. A ligand sitting in the right pocket but flipped end
for end therefore counts as placement-limited even though it is correctly
localised. The metric that actually tests location is the 4 Angstrom centroid
clustering two slides later.
Anticipated challenge, and a fair one: AutoDock was tuned. It ran a ladder of five
exhaustiveness values and the top rung was carried forward, while the two learned
arms were only given a choice of refiner. The winner's-curse correction is between
0.4 and 0.8 percentage points per family and the largest selection margin is 3.0
points, so it does not overturn the ordering, but the comparison is a tuned
pipeline against two untuned ones.
Anticipated challenge: no held-out split. Correct, variant selection and
evaluation use the same 303 complexes.
Anticipated challenge: is the non-significant rank-1 result equivalence? No. The
minimum detectable difference there is 9.8 points at 80% power and no equivalence
margin was prespecified.
Anticipated challenge: the 1 Angstrom form threshold is not a published
convention. Correct, it was fixed before analysis as half the 2 Angstrom in-place
threshold, and the next slide plus the appendix resolve it across thresholds.
""")

# ============================= 7  ALL CRITERIA AT RANK-1 AND TOP-15 (table)
s, ph = new_slide(L_ONE, "All Criteria at Rank-1 and Top-15")
drop(ph[14])
CRIT = [
    ("Valid near-native, RMSD ≤ 2 Å", "36.6", "65.3", "34.0", "55.1", "18.2", "26.4"),
    ("Strict near-native, RMSD ≤ 1 Å", "22.1", "34.7", "15.5", "30.0", "6.9", "10.2"),
    ("Form, best-fit RMSD ≤ 1 Å", "52.5", "71.9", "52.1", "72.3", "33.3", "51.2"),
    ("Form, best-fit RMSD ≤ 2 Å", "78.9", "93.4", "82.5", "95.7", "61.4", "79.9"),
    ("Crystal pocket reached, ≤ 4 Å", "59", "83", "54", "81", "43", "47"),
    ("Median in-place RMSD (Å)", "1.49", "5.04", "1.47", "1.93", "2.98", "4.07"),
    ("Median best-fit RMSD (Å)", "0.85", "1.45", "0.89", "0.92", "1.12", "1.34"),
    ("Placement-limited poses (%)", "47.9", "77.8", "40.9", "47.9", "59.4", "70.2"),
]
col_w = [3.06] + [1.09] * 6
tbl = table(s, 0.20, 0.70, 9.60, col_w, [None] * (len(CRIT) + 2))
cell(tbl, 0, 0, "Criterion", size=9.5, bold=True, colour=RGBColor(0xFF,0xFF,0xFF),
     fill=BLUE, align="l")
tbl.cell(0, 0).merge(tbl.cell(1, 0))
for k, (name, c0) in enumerate((("AutoDock*", 1), ("DiffDock*", 3),
                                ("EquiBind*", 5))):
    cell(tbl, 0, c0, name, size=9.5, bold=True,
         colour=RGBColor(0xFF, 0xFF, 0xFF), fill=BLUE)
    cell(tbl, 0, c0 + 1, "", fill=BLUE)
    tbl.cell(0, c0).merge(tbl.cell(0, c0 + 1))
    for d, lbl in ((0, "rank-1"), (1, "top-15")):
        cell(tbl, 1, c0 + d, lbl, size=8.5, bold=True, colour=BLUE, fill=LIGHT)
for i, row in enumerate(CRIT):
    r = i + 2
    band = RGBColor(0xF5, 0xF7, 0xF9) if i % 2 else None
    rule = i in (4, 5)
    cell(tbl, r, 0, row[0], size=9, align="l", fill=band,
         bold=False, colour=DARK)
    for j, v in enumerate(row[1:]):
        best = j in (0, 1) and i < 5
        cell(tbl, r, j + 1, v, size=9, fill=band,
             bold=False, colour=DARK)
tf = textbox(s, 0.20, 3.86, 9.60, 0.52)
bullets(tf, [
    "Placement separates the two leaders, form does not, on the same poses "
    "at the same depth",
    "AutoDock's median in-place RMSD degrades from 1.5 to 5.0 Å with depth, "
    "so its ranking is informative only near the top",
], size=10, gap=3)
caveat(s, 0.20, 4.44, 9.60,
       "Selected variants: AutoDock Vina exh128 + gnina, DiffDock + smina, "
       "EquiBind unguided + gnina. The first five rows are the percentage of 303 "
       "complexes holding at least one PoseBusters-valid pose that meets the "
       "criterion within the depth. The last three describe PoseBusters-valid "
       "poses lying within 8 Å of the crystal site. Placement-limited is a "
       "residual of the two RMSD measures and bundles orientation with position, "
       "so it is not a statement about the pocket alone.", h=0.72)
notes(s, """
This is the numeric backup for the previous slide and the one to turn to if the
examiner asks for anything the curves do not show.
What to point out. Rows one and two show that tightening the distance threshold
from 2 to 1 Angstrom widens the AutoDock lead at rank-1, from 2.6 to 6.6 points,
so the leaders are closer at the conventional threshold than at a strict one.
Rows three and four show the form tie surviving both thresholds, and DiffDock is
in fact nominally ahead at 2 Angstrom.
Row five is the coarsest criterion, a cluster centroid within 4 Angstroms of the
crystal pocket, and it is where every tool scores highest. Reaching the pocket is
much easier than nominating the pose.
Rows six to eight are pose-level, not complex-level, and they are the diagnostic
rows. AutoDock's median in-place RMSD moves from 1.49 to 5.04 Angstroms between
rank-1 and top-15 while DiffDock moves from 1.47 to 1.93, so deeper AutoDock ranks
add correctly shaped but badly placed ligands. Its placement-limited share rises
from 48 to 78 percent, which is the same statement in a different currency.
Caveat to volunteer: the last three rows use a near-site cohort, trimmed to poses
within 8 Angstroms of the crystal ligand, so they describe the poses that landed
near the right place and are not comparable with the first five rows.
Second caveat: the in-place and best-fit values each minimise over
symmetry-equivalent atom mappings independently, so the placement decomposition is
approximate. It is also a residual rather than a measurement, and it bundles
orientation with position, so a correctly located but rotated pose lands in the
placement-limited bin.
""")

# ====================================== 8  SITE RECOVERY + CONTACT FIDELITY
s, ph = new_slide(L_ONE, "Finding the Pocket Is the Easy Part")
drop(ph[14])
picture(s, "site_reach.png", 0.20, 0.66, 5.35, 4.20)
tf = textbox(s, 5.75, 0.74, 4.05, 4.15)
col_header(tf, "Site recovery")
bullets(tf, [
    "Retrospectively, the pooled cloud of all three tools contains the "
    "crystal pocket for 95% of complexes",
    "The best cheap blind rule ranks that pocket first in only 57%",
    "A 38-point gap in selection, not in accuracy",
], size=10.5, gap=6)
para(tf, "Contact fidelity", size=13, bold=True, colour=BLUE, bullet=False,
     space_before=13, space_after=4)
bullets_items = [
    "Interaction-fingerprint overlap with the crystal, top-5 union: "
    "DiffDock 0.45, AutoDock 0.42, EquiBind 0.38",
    "The two leaders resemble each other more (0.46) than either "
    "resembles the crystal",
]
for i, it in enumerate(bullets_items):
    para(tf, it, size=10.5, space_before=6 if i else 0)
para(tf, "Different contact hypotheses, so agreement between tools is not "
         "evidence of correctness.",
     size=9.5, italic=True, colour=GREY, bullet=False, space_before=10)
notes(s, """
The retrieval framing lands well here. Pooling every candidate from every tool
almost always contains the right pocket. The unsolved problem is ranking that
pocket first without a reference, which is a re-ranking problem, not a sampling
problem.
Note on the 95% number: it is a cluster-level oracle at a 4 Angstrom centroid
threshold, which is much coarser than a valid pose within 2 Angstrom. It must not
be quoted as 95% pose accuracy, and the slide deliberately labels it as an oracle.
On Jaccard: it counts shared and non-shared contacts, so it is a similarity, not a
recall. At rank-1 the two leaders are not separable on contact F1 either, and
DiffDock's lead exists only for the union of its top five.
""")

# ===================================================== 9  ORAI USABLE YIELD
s, ph = new_slide(L_ONE, "On Orai1, Physical Validity Overstates Usefulness")
drop(ph[14])
tf = textbox(s, 0.20, 0.63, 9.60, 0.68)
bullets(tf, [
    "All 120 AutoDock poses pass the physics checks, yet only 54% are usable",
    "DiffDock is less valid at 93% and more usable at 68%, EquiBind keeps 1 of 120",
    "The failure is placement, not chemistry, the rejected poses sit in the pore",
], size=10.5, gap=3)
picture(s, f"{MEDIA}/image43.png", 0.55, 1.36, 8.90, 3.30)
caveat(s, 0.20, 4.72, 9.60,
       "The reversal is a property of the secondary ranker, not of the search. "
       "Capping AutoDock on its raw Vina order instead of the gnina re-ranking "
       "returns it to 75% and removes the reversal.")
notes(s, """
This is the most quotable Orai1 result and also the most attackable, so present
the caveat yourself before the examiner does.
Three separate caveats. First, the exclusion slab is a modelling hypothesis. It
encodes an outer-pore binding hypothesis and may well reject valid GSK-7975A or
Synta-66 poses, since both act at or near the pore mouth. Second, the two panels
were capped by different selection rules. The experimental ten poses are the
gnina-best of thirty searched modes, the control ten are all that were written.
Matching that moves the AutoDock loss from 28.7 to 17.0 points and its adjusted
p from 5.6e-4 to above 0.2. Third, the between-tool reversal disappears entirely
on the Vina prefix, so it is an artefact of the gnina ranker.
If asked whether the recovered residues validate the published Synta-66 site: no.
A decoy control contacts the same residues at the same rate, so the concordance
is not specific.
""")

# ============================================= 10  ORAI CROSS-TOOL AGREEMENT
s, ph = new_slide(L_ONE, "The Tools Do Not Agree on Where the Ligand Binds")
drop(ph[14])
tf = textbox(s, 0.20, 0.63, 9.60, 0.68)
bullets(tf, [
    "Median distance between two tools' consensus sites is about 30 Å on the "
    "controls and 33 to 44 Å on the modulators",
    "AutoDock and DiffDock land within 5 Å in 14% of control pairs and none of 11",
    "Six-fold channel symmetry does not explain it, folding onto the nearest "
    "copy moves the median only to 27.3 Å",
], size=10.5, gap=3)
picture(s, "orai_agreement_panelA.png", 0.55, 1.36, 8.90, 3.30)
caveat(s, 0.20, 4.72, 9.60,
       "Consensus between independent tools is therefore not usable as a "
       "confidence signal on this target.")
notes(s, """
A negative result, and a useful one. The intuitive safeguard, trusting a pose
because two independent methods agree on it, does not survive contact with a large
elongated membrane protein.
Anticipated question: could this be a coordinate-frame or symmetry artefact? Both
were tested. All three tools dock into the same receptor file per frame, so
within-frame comparisons share one coordinate system with no folding applied, and
the six-fold symmetry fold changes the median by two Angstroms.
Anticipated question: does filtering to valid, outside-slab poses tighten the
clouds? No, not detectably. Filtering improves bootstrap reproducibility but not
spread. Note that no equivalence test was run, so this bounds the filtering effect
rather than proving the scatter is intrinsic.
Statistical honesty point: eleven evaluable experimental pairs. This is
descriptive, not a corrected finding.
""")

# =================================================== 11  RQ3 COST
s, ph = new_slide(L_ONE, "Cost per Usable Pose Spans Two Orders of Magnitude")
drop(ph[14])
picture(s, f"{MEDIA}/image24.png", 0.20, 0.68, 5.20, 4.10)
tf = textbox(s, 5.65, 0.76, 4.15, 4.15)
col_header(tf, "Charged seconds per qualifying pose")
bullets(tf, [
    "EquiBind 2.4 s, DiffDock 49.4 s, AutoDock 137.2 s",
    "One basis for all three, CPU-core-seconds / 32 threads + GPU-seconds",
    "Success rates differ, 80, 169 and 199 complexes out of 303",
    "68.5% of the AutoDock cost is the gnina rescoring pass, not the search",
    "Vina alone needs 21.7 s and no GPU, cheaper than DiffDock",
], size=11, gap=7)
para(tf, "AutoDock buys coverage and admissibility, not economy. EquiBind's "
         "figure is a floor, conditional on succeeding a quarter of the time.",
     size=9.5, italic=True, colour=GREY, bullet=False, space_before=10)
notes(s, """
The accounting basis is the part an examiner will press on, so state it first.
Charging device occupancy rather than wall clock stops work that merely ran with
more concurrency from looking cheaper. The gnina pass ran sixteen ways in
parallel and finished in 0.74 wall-hours, and charging it that way would have
flattered AutoDock to 5.47 total hours instead of 15.00.
Anticipated challenge: cost conditional on success favours whoever fails most.
True, and it is stated. EquiBind's 2.4 s describes only its 80 successful
complexes and it fails on the other 223. The alternative denominators are also
reported and none of them changes the ordering.
Anticipated challenge: the rescorer was invoked once per pose, so each call paid
process start-up on a mostly idle GPU. That is an implementation artefact.
Batching would reduce it by an amount this run does not determine.
""")

# ============================================== 12  ANSWERS RQ1 AND RQ2
s, ph = new_slide(L_TWO, "Answering RQ1 and RQ2")
tf = ph[15].text_frame; tf.word_wrap = True
col_header(tf, "RQ1  Validity, near-nativeness, ranking depth")
bullets(tf, [
    "Rank-1 valid near-native: 36.6%, 34.0% and 18.2%",
    "Best-of-top-15: 65.3%, 55.1% and 26.4%",
    "The leaders separate only from top-5 onwards",
    "Pooled physical validity: 98.9%, 84.5% and 62.0%",
    "Both leaders are decisively ahead of EquiBind everywhere",
    ("Rank-1 is unresolved, not shown to be equivalent", 1),
    ("The ordering is conditional on search effort, AutoDock ran the top "
     "rung of a five-value ladder", 1),
], size=11, gap=5)

tf = ph[16].text_frame; tf.word_wrap = True
col_header(tf, "RQ2  Effect of post-docking optimisation")
bullets(tf, [
    "DiffDock validity rises from 24.4% raw to 84.5%",
    "EquiBind rises from 2.9% to 62.0%",
    "Accuracy at a fixed rank moves by only 0.1 to 4.6 points",
    "gnina re-ranking lifts AutoDock rank-1 from 31.0% to 36.6%",
    "Optimisation repairs geometry and re-orders a pool",
    ("It cannot create a binding mode that was never sampled", 1),
    ("Refiner choice, smina against gnina, was not statistically "
     "separable within either family", 1),
], size=11, gap=5)
notes(s, """
Keep the two answers separable. RQ1 is a ranking of pipelines with an explicit
unresolved case. RQ2 is a mechanism claim, and the mechanism is repair plus
re-ordering, never generation.
The single most likely examiner question on this slide: "so is AutoDock better?"
The defensible answer is that on this benchmark, under blind whole-protein search,
with the variants carried forward, AutoDock has the higher point estimate
everywhere and a statistically separable lead only from top-5 onwards. It is not
a claim about physics-based docking in general, and the AutoDock arm is itself
CNN-rescored, so it is not even a pure physics arm.
""")

# ============================================== 13  ANSWERS RQ3 AND RQ4
s, ph = new_slide(L_TWO, "Answering RQ3 and RQ4")
tf = ph[15].text_frame; tf.word_wrap = True
col_header(tf, "RQ3  Computational cost")
bullets(tf, [
    "2.4 s, 49.4 s and 137.2 s per qualifying pose",
    "Configured AutoDock is dearest on every currency",
    "Learned rescoring is 68.5% of that total",
    "Removing gnina puts Vina below DiffDock on charged time",
    "Vina still has by far the highest CPU cost",
    ("Describes these implementations and this hardware, "
     "not method classes", 1),
], size=11, gap=5)

tf = ph[16].text_frame; tf.word_wrap = True
col_header(tf, "RQ4  Opportunities and limits")
bullets(tf, [
    "The three tools fail in different ways, which favours a staged workflow",
    "Generate, optimise, screen validity and placement, then rank for the "
    "endpoint you care about",
    "Generation and localisation dominate failure, not ranking",
    "Blind whole-protein search is harder than published site-centred "
    "re-docking",
    "Refined DiffDock and AutoDock Vina are complementary",
    ("EquiBind gives too little coverage for routine use in this "
     "configuration", 1),
], size=11, gap=5)
notes(s, """
RQ3 depends entirely on the denominator, so name it before quoting the numbers.
RQ4 is deliberately a workflow answer rather than a winner, because the failure
modes are not the same failure mode. AutoDock fails by mis-ranking a pool that
contains the answer, DiffDock fails by not localising, EquiBind fails by not
generating anything usable.
If pressed on why blind whole-protein search was chosen when Vina normally gets a
box: because DiffDock and EquiBind localise the ligand themselves, so giving Vina
a pocket would have handed it information the others earn. This is stricter than
blind benchmarks that supply a predicted pocket, and a predicted-pocket Vina arm
would be the right intermediate comparison. It was not run.
""")

# ============================================ 14  CONCLUSIONS 1
s, ph = new_slide(L_ONE, "Conclusions: What the Study Establishes")
drop(ph[14])
tf = textbox(s, 0.20, 0.70, 5.85, 4.20)
bullets(tf, [
    "Under blind whole-protein search a tuned classical pipeline still matches "
    "or beats two learned ones",
    "Post-hoc physical optimisation is what makes learned poses admissible, "
    "and it is cheap",
    "Rank-1 and best-of-top-k tell different stories, so report both",
    "Physical validity is an admissibility filter, never evidence of binding",
    "Without a reference structure, validity has to be paired with an explicit "
    "placement rule",
    "Agreement between independent tools is not a confidence signal",
    "The contribution is a separation of sampling, ranking, validity and "
    "localisation failures",
], size=11.5, gap=8)

chip(s, 6.30, 0.90, 3.50, 1.30,
     "Recommended workflow\ngenerate → optimise → screen validity and "
     "placement → rank for the endpoint", size=11)
tf = textbox(s, 6.30, 2.45, 3.50, 2.40)
col_header(tf, "Practical selection")
bullets(tf, [
    "AutoDock Vina and refined DiffDock as complementary hypothesis generators",
    "EquiBind not in this configuration",
    "Retain uncertainty wherever the two disagree",
], size=10.5, gap=6)
notes(s, """
Close on the workflow, not on a winner. The examiner should leave with the four
separated failure modes, since that is the transferable contribution and it
survives the arrival of the next generation of models.
If asked for the single most useful practical finding: one cheap local
minimisation pass converts most learned output from physically impossible to
admissible, at a cost that is negligible next to the generation step itself.
""")

# ============================================ 15  CONCLUSIONS 2
s, ph = new_slide(L_TWO, "Conclusions: Limits and Outlook")
tf = ph[15].text_frame; tf.word_wrap = True
col_header(tf, "What this work does not show")
bullets(tf, [
    "Not AI against physics, three configured pipelines at single settings",
    "Variant selection and evaluation share the same 303 complexes, "
    "with no held-out split",
    "The benchmark DiffDock run is one unseeded draw, no arm was "
    "run across multiple seeds",
    "Post-2021 release does not exclude training-set leakage",
    "Orai1 rests on three ligands and four correlated frames of one trajectory",
    "No endpoint tests binding free energy, protonation or stability",
], size=10.5, gap=6)

tf = ph[16].text_frame; tf.word_wrap = True
col_header(tf, "Where it goes next")
bullets(tf, [
    "Match information and compute against current hybrid and co-folding models",
    "Split by protein, pocket and ligand novelty rather than by release date",
    "Add a predicted-pocket Vina arm as the missing intermediate comparison",
    "Batch the learned rescorer, the GPU was idle on a per-pose invocation",
    "Treat the Orai1 poses as testable hypotheses for mutagenesis, competition "
    "assays and ligand-bound simulation",
    "The experimental follow-up continues at JKU Linz",
], size=10.5, gap=6)
notes(s, """
Lead with the limitations rather than waiting to be asked. The strongest examiner
criticisms are already on this slide, which converts them from attacks into
evidence of judgement.
Two further limits worth having ready but not on the slide. Receptor preparation
stripped waters, ions, metals and cofactors for the AutoDock and DiffDock arms,
and a metal sits within 5 Angstroms of the crystal ligand in 78 of the 303
complexes, so a quarter of the set is scored without a coordinating partner the
native pose depends on. The AutoDock lead at top-15 is +6.7 points on the 225
metal-free complexes against +20.5 on the 78 metal-adjacent ones, so the pooled
contrast averages a narrow gap with a much wider one.
Second, the interaction profiler is blind to chlorine and bromine in its halogen
test, so any halogen-bond result is discounted.
""")

prs.save(OUT)
print("saved", OUT, "slides:", len(prs.slides))
