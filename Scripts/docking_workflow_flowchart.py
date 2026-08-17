"""Horizontal flowchart of the physics-based molecular-docking workflow.

Renders the five sequential stages described in the methods text — protein/ligand
preparation, cavity detection, pose sampling, scoring, and post-docking analysis —
as two rows of cards joined by arrows (a serpentine: three cards across the top
row, two across the bottom, with a wrap-around connector between them). Each card carries a
coloured stage badge, a short title, the source's own italic sub-label (keeping the
"sampling = conformation & position" / "scoring = binding affinity" framing) and a
few concise bullet points with the original citation markers.

Pure matplotlib (no graphviz / external diagram tool), so it runs anywhere the
`vina` env does. Usable two ways:

    # command line
    python Scripts/docking_workflow_flowchart.py -o docking_workflow.png

    # from the notebook (Flow Charts.ipynb)
    from Scripts.docking_workflow_flowchart import build_flowchart
    fig = build_flowchart()          # returns the Matplotlib Figure

Colours are the three Pantone brand colours — 301C blue (#00649c), 376C green
(#8bb31d) and 431C grey (#72777a) — applied per stage at different opacities so
the whole figure stays on-brand while every stage remains distinct. Badge text is
auto-contrasted to its fill and sub-labels are nudged to a legible tone, so every
label clears a comfortable contrast ratio.
"""
from __future__ import annotations

import argparse
import textwrap
from pathlib import Path
from typing import Dict, List, Sequence

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch
from matplotlib.path import Path as MplPath


# --------------------------------------------------------------------------- #
# Content — distilled straight from the methods passage, one dict per stage.   #
# `refs` are the original bracketed citation numbers, rendered as superscripts. #
# --------------------------------------------------------------------------- #
STAGES: List[Dict] = [
    {
        "title": "System preparation",
        "subtitle": "protein & ligand set-up",
        "bullets": [
            ("Strip crystallographic waters, add hydrogens, "
             "assign protonation states", "17"),
            ("Optimise geometry for stable conformations", "15"),
            ("Build 3D ligand, assign partial charges "
             "(Gasteiger–Marsili)", "18"),
            ("Generate multiple conformers for ligand flexibility", "19"),
        ],
    },
    {
        "title": "Cavity detection",
        "subtitle": "binding-site identification",
        "bullets": [
            ("Scan the protein surface for concave regions / buried "
             "volumes using grids or spheres", "6"),
            ("Rank pockets by size, shape and depth", None),
            ("Score volume, hydrophobicity and H-bond potential", "20"),
            ("Blind (whole surface) vs site-directed (known pocket) "
             "docking", "15"),
        ],
    },
    {
        "title": "Pose sampling",
        "subtitle": "sampling — conformation & position",
        "bullets": [
            ("Flexibly place the ligand into the binding site", None),
            ("Explore orientations & conformations "
             "(e.g. Monte Carlo)", None),
            ("Model non-covalent interactions: van der Waals, "
             "electrostatics, H-bonds, water displacement", "21"),
        ],
    },
    {
        "title": "Scoring",
        "subtitle": "scoring — binding affinity",
        "bullets": [
            ("A scoring function estimates the binding affinity of "
             "each pose", None),
            ("Lowest predicted free energy ≈ most stable "
             "complex", "19"),
            ("Scoring differs across packages — software choice "
             "matters (AutoDock, HADDOCK …)", None),
        ],
    },
    {
        "title": "Post-docking analysis",
        "subtitle": "ranking & validation",
        "bullets": [
            ("Rank and validate the surviving poses", None),
            ("Inspect key interactions: H-bonds to catalytic "
             "residues, hydrophobic contacts", None),
            ("Assess complex stability, e.g. MD simulation", "10, 11"),
        ],
    },
]

# --- Brand palette (three Pantone spot colours) ---------------------------- #
BRAND_BLUE = "#00649c"    # Pantone 301C   RGB 0/100/156
BRAND_GREEN = "#8bb31d"   # Pantone 376C   RGB 139/179/29
BRAND_GREY = "#72777a"    # Pantone 431C   RGB 114/120/122

BRAND_GREEN_LT = "#bcd37c"   # a lighter tint of 376C green

# Five stages coloured from the brand palette, each hue reused at a lower opacity
# for its second stage so the pairs stay related but distinct. Each entry is
# (base colour, opacity):
#   1 System preparation    blue        100%
#   2 Cavity detection      blue         80%
#   3 Pose sampling         light green 100%
#   4 Scoring               light green  80%
#   5 Post-docking analysis grey         60%
PALETTE_BRAND = [
    (BRAND_BLUE, 1.00),
    (BRAND_BLUE, 0.80),
    (BRAND_GREEN_LT, 1.00),
    (BRAND_GREEN_LT, 0.80),
    (BRAND_GREY, 0.60),
]
# Single-hue alternative: brand blue stepped through five opacities.
PALETTE_MONO = [(BRAND_BLUE, a) for a in (1.00, 0.82, 0.64, 0.48, 0.34)]

# Chrome / ink tokens (light mode).
INK_PRIMARY = "#0b0b0b"
INK_SECONDARY = "#40403d"
INK_MUTED = "#7c7b75"
SURFACE = "#ffffff"

# Dark-mode ink / surface overrides (used when theme="dark").
DARK = {
    "ink_primary": "#ffffff",
    "ink_secondary": "#d4d3ca",
    "ink_muted": "#9a998f",
    "surface": "#1f1f1d",
    "page": "#141413",
}


def _rgb(hex_color: str):
    return tuple(int(hex_color[i:i + 2], 16) for i in (1, 3, 5))


def _hex(rgb) -> str:
    return "#{:02x}{:02x}{:02x}".format(*(max(0, min(255, round(c))) for c in rgb))


def _blend(fg: str, bg: str, alpha: float) -> str:
    """Composite ``fg`` over ``bg`` at ``alpha`` opacity → an opaque hex."""
    f, b = _rgb(fg), _rgb(bg)
    return _hex(alpha * f[k] + (1 - alpha) * b[k] for k in range(3))


def _contrast(a: str, b: str) -> float:
    la, lb = _luminance(a), _luminance(b)
    hi, lo = max(la, lb), min(la, lb)
    return (hi + 0.05) / (lo + 0.05)


def _readable(fg: str, bg: str, target: float = 3.4) -> str:
    """Nudge ``fg`` toward black (light bg) or white (dark bg) until it clears
    ``target``:1 contrast — keeps a brand hue legible as small text."""
    if _contrast(fg, bg) >= target:
        return fg
    toward = "#ffffff" if _luminance(bg) < 0.2 else "#000000"
    for k in range(1, 21):
        c = _blend(fg, toward, 1 - k * 0.045)
        if _contrast(c, bg) >= target:
            return c
    return toward


def _luminance(hex_color: str) -> float:
    """WCAG relative luminance of an ``#rrggbb`` colour."""
    r, g, b = (int(hex_color[i:i + 2], 16) / 255 for i in (1, 3, 5))

    def lin(c: float) -> float:
        return c / 12.92 if c <= 0.03928 else ((c + 0.055) / 1.055) ** 2.4

    return 0.2126 * lin(r) + 0.7152 * lin(g) + 0.0722 * lin(b)


def _text_on(hex_color: str) -> str:
    """White text where it clears 3:1 on ``hex_color``, else black.

    Keeps white numerals on the deep hues (blue/violet/orange) but flips to
    black on the mid-luminance ones (amber/aqua) where white would be too faint.
    """
    contrast_white = 1.05 / (_luminance(hex_color) + 0.05)
    return "#ffffff" if contrast_white >= 3.0 else "#111111"


def _superscript(refs: str) -> str:
    """Turn a ref string like ``"10, 11"`` into unicode superscripts."""
    table = {"0": "⁰", "1": "¹", "2": "²", "3": "³",
             "4": "⁴", "5": "⁵", "6": "⁶", "7": "⁷",
             "8": "⁸", "9": "⁹", ",": "˒", " ": " "}
    return "".join(table.get(ch, ch) for ch in refs)


def _unit(p, q):
    """Unit vector pointing from ``p`` to ``q`` (zero vector if coincident)."""
    dx, dy = q[0] - p[0], q[1] - p[1]
    d = (dx * dx + dy * dy) ** 0.5
    return (0.0, 0.0) if d == 0 else (dx / d, dy / d)


def _rounded_path(points, radius: float = 0.20) -> "MplPath":
    """Build a :class:`~matplotlib.path.Path` through ``points`` with the
    interior corners rounded off by ``radius`` (quadratic-Bezier fillets).

    ``points`` is an orthogonal polyline; used to route the serpentine
    connector that wraps from the end of the top row to the start of the
    bottom row without cutting diagonally across the cards.
    """
    verts = [points[0]]
    codes = [MplPath.MOVETO]
    for i in range(1, len(points) - 1):
        p_prev, p, p_next = points[i - 1], points[i], points[i + 1]
        vin, vout = _unit(p_prev, p), _unit(p, p_next)
        a = (p[0] - vin[0] * radius, p[1] - vin[1] * radius)
        b = (p[0] + vout[0] * radius, p[1] + vout[1] * radius)
        verts += [a, p, b]
        codes += [MplPath.LINETO, MplPath.CURVE3, MplPath.CURVE3]
    verts.append(points[-1])
    codes.append(MplPath.LINETO)
    return MplPath(verts, codes)


def build_flowchart(
    stages: Sequence[Dict] = STAGES,
    *,
    theme: str = "light",
    palette: str = "categorical",
    show_refs: bool = False,
    title: str = "Physics-based molecular docking workflow",
    subtitle: str = ("Sampling predicts ligand conformation & position; "
                     "scoring estimates binding affinity"),
    width: float = 18.0,
    card_height: float = 5.9,
):
    """Build and return the flowchart :class:`~matplotlib.figure.Figure`.

    All geometry is in inches on an equal-aspect axis so rounded corners stay
    circular and arrows keep their proportions at any figure size.
    """
    n = len(stages)

    # Resolve theme tokens.
    if theme == "dark":
        ink_primary, ink_secondary, ink_muted = (
            DARK["ink_primary"], DARK["ink_secondary"], DARK["ink_muted"])
        surface, page = DARK["surface"], DARK["page"]
    else:
        ink_primary, ink_secondary, ink_muted = INK_PRIMARY, INK_SECONDARY, INK_MUTED
        surface, page = SURFACE, "#f7f7f5"

    # Brand colours applied at their per-stage opacity. Accents are composited
    # over white so an opacity of 1.0 is the true spot colour and lower values
    # read as lighter tints of it; card bodies get a faint wash of the same hue,
    # and the connector arrows use brand grey.
    accents = PALETTE_MONO if palette == "mono" else PALETTE_BRAND
    accent_colours = [_blend(base, "#ffffff", op) for base, op in accents]
    body_tints = [_blend(base, surface, op * 0.12) for base, op in accents]
    arrow = _blend(BRAND_GREY, surface, 0.55)

    # --- layout (inches) --------------------------------------------------- #
    # Two rows (serpentine): the first ceil(n/2) stages fill the top row, the
    # rest the bottom row. Card width is sized to the wider (top) row so every
    # card is identical; the shorter bottom row is left-aligned beneath it and
    # linked by a wrap-around connector routed through the gap between rows.
    margin_x = 0.5
    gap = 0.66                      # horizontal space between cards (holds arrow)
    row_gap = 1.15                  # vertical corridor between the two rows
    title_band = 1.9                # headroom for the figure title + subtitle
    bottom_pad = 0.4

    top_count = (n + 1) // 2
    cols_per_row = top_count
    row_of = lambda i: 0 if i < top_count else 1
    col_of = lambda i: i if i < top_count else i - top_count

    card_w = (width - 2 * margin_x - (cols_per_row - 1) * gap) / cols_per_row
    height = title_band + 2 * card_height + row_gap + bottom_pad

    fig = plt.figure(figsize=(width, height), dpi=100)
    fig.patch.set_facecolor(page)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, width)
    ax.set_ylim(0, height)
    ax.set_aspect("equal")
    ax.axis("off")

    # Per-row vertical extents (top row 0 sits above bottom row 1).
    row_top = [height - title_band, height - title_band - card_height - row_gap]
    row_bottom = [row_top[0] - card_height, row_top[1] - card_height]
    row_ymid = [(row_top[r] + row_bottom[r]) / 2 for r in (0, 1)]

    def card_x(i):
        x_left = margin_x + col_of(i) * (card_w + gap)
        return x_left, x_left + card_w

    # --- figure title ------------------------------------------------------ #
    ax.text(width / 2, height - 0.64, title, ha="center", va="center",
            fontsize=29, fontweight="bold", color=ink_primary)
    ax.text(width / 2, height - 1.28, subtitle, ha="center", va="center",
            fontsize=17.5, style="italic", color=ink_secondary)

    pad = 0.24
    badge_r = 0.34
    mono = FontProperties(family="DejaVu Sans")

    for i, stage in enumerate(stages):
        row = row_of(i)
        x_left, x_right = card_x(i)
        card_top = row_top[row]
        card_bottom = row_bottom[row]
        accent = accent_colours[i]

        # Card: faint hue wash, coloured border, subtle rounded corners.
        card = FancyBboxPatch(
            (x_left, card_bottom), card_w, card_height,
            boxstyle="round,pad=0,rounding_size=0.12",
            facecolor=body_tints[i], edgecolor=accent, linewidth=1.8,
            mutation_aspect=1.0, zorder=2,
        )
        ax.add_patch(card)
        # Thin coloured accent rail down the left edge of the card.
        rail = FancyBboxPatch(
            (x_left, card_bottom), 0.10, card_height,
            boxstyle="round,pad=0,rounding_size=0.05",
            facecolor=accent, edgecolor="none", zorder=3,
        )
        ax.add_patch(rail)

        # Stage badge (numbered circle) at the top-left, with a "STAGE n" kicker.
        bx = x_left + pad + badge_r + 0.04
        by = card_top - pad - badge_r
        ax.add_patch(Circle((bx, by), badge_r, facecolor=accent,
                             edgecolor="none", zorder=4))
        ax.text(bx, by, str(i + 1), ha="center", va="center",
                fontsize=19, fontweight="bold",
                color=_text_on(accent), zorder=5,
                fontproperties=mono)
        ax.text(bx + badge_r + 0.18, by, f"STAGE {i + 1}", ha="left", va="center",
                fontsize=13, fontweight="bold", color=ink_muted,
                fontproperties=mono)

        # Title on its own full-width line under the badge (fits long titles),
        # then the source's italic sub-label.
        title_y = by - badge_r - 0.36
        ax.text(x_left + pad, title_y, stage["title"], ha="left", va="center",
                fontsize=20, fontweight="bold", color=ink_primary)
        sub_y = title_y - 0.48
        ax.text(x_left + pad, sub_y, stage["subtitle"], ha="left", va="center",
                fontsize=14.5, style="italic", color=_readable(accent, surface))

        # Divider under the header block.
        div_y = sub_y - 0.28
        ax.plot([x_left + pad, x_right - pad], [div_y, div_y],
                color=accent, linewidth=1.2, alpha=0.55, zorder=3)

        # Bullets, wrapped to the card's inner width with a hanging indent.
        inner_w = card_w - 2 * pad - 0.10
        wrap_chars = max(18, int(inner_w / 0.116))
        y = div_y - 0.38
        line_h = 0.35
        for text, refs in stage["bullets"]:
            body = text
            if show_refs and refs:
                body = f"{text} {_superscript(refs)}"
            wrapped = textwrap.wrap(body, width=wrap_chars) or [body]
            block = "•  " + "\n   ".join(wrapped)
            ax.text(x_left + pad + 0.06, y, block, ha="left", va="top",
                    fontsize=15.5, color=ink_secondary,
                    linespacing=1.32, zorder=4)
            y -= len(wrapped) * line_h + 0.19

        # Arrow to the next card within the same row.
        if i + 1 < n and row_of(i + 1) == row:
            a0 = x_right + 0.08
            a1 = x_right + gap - 0.08
            ax.add_patch(FancyArrowPatch(
                (a0, row_ymid[row]), (a1, row_ymid[row]),
                arrowstyle="-|>", mutation_scale=24,
                linewidth=2.6, color=arrow, zorder=1,
                shrinkA=0, shrinkB=0,
            ))

    # --- serpentine wrap connector: end of top row -> start of bottom row -- #
    if n > top_count:
        top_left, top_right = card_x(top_count - 1)
        bot_left, bot_right = card_x(top_count)
        top_cx = (top_left + top_right) / 2
        bot_cx = (bot_left + bot_right) / 2
        corridor_y = (row_bottom[0] + row_top[1]) / 2
        pts = [
            (top_cx, row_bottom[0] - 0.02),   # drop out of the last top card
            (top_cx, corridor_y),             # into the inter-row corridor
            (bot_cx, corridor_y),             # sweep left across the corridor
            (bot_cx, row_top[1] + 0.02),      # down into the first bottom card
        ]
        ax.add_patch(FancyArrowPatch(
            path=_rounded_path(pts, radius=0.22),
            arrowstyle="-|>", mutation_scale=26,
            linewidth=2.6, color=arrow, zorder=1,
            shrinkA=0, shrinkB=0,
        ))

    return fig


def main(argv: Sequence[str] | None = None) -> None:
    matplotlib.use("Agg")  # CLI is headless — only ever saving to file
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("-o", "--output", default="docking_workflow_flowchart.png",
                   help="output image path (extension sets the format)")
    p.add_argument("--theme", choices=["light", "dark"], default="light")
    p.add_argument("--palette", choices=["brand", "mono"], default="brand",
                   help="three brand colours across the stages, or brand blue "
                        "stepped through opacities")
    p.add_argument("--refs", action="store_true",
                   help="show the [nn] citation superscripts (off by default)")
    p.add_argument("--dpi", type=int, default=300)
    p.add_argument("--width", type=float, default=18.0,
                   help="figure width in inches")
    p.add_argument("--also-pdf", action="store_true",
                   help="additionally write a vector PDF beside the raster output")
    args = p.parse_args(argv)

    fig = build_flowchart(theme=args.theme, palette=args.palette,
                          show_refs=args.refs, width=args.width)

    out = Path(args.output)
    fig.savefig(out, dpi=args.dpi, bbox_inches="tight",
                facecolor=fig.get_facecolor())
    print(f"wrote {out.resolve()}")
    if args.also_pdf:
        pdf = out.with_suffix(".pdf")
        fig.savefig(pdf, bbox_inches="tight", facecolor=fig.get_facecolor())
        print(f"wrote {pdf.resolve()}")
    plt.close(fig)


if __name__ == "__main__":
    main()
