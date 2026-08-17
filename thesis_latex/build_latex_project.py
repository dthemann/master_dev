from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parent
RAW = ROOT / "body_raw.tex"
FULLTEXT = ROOT.parent / "thesis_revision" / "final_thesis_fulltext.md"
MAIN_OUT = ROOT / "body_main.tex"
APPENDIX_OUT = ROOT / "body_appendix.tex"
BIB_OUT = ROOT / "Literatur.bib"
THESIS_OUT = ROOT / "Thesis.tex"
README_OUT = ROOT / "README.md"


@dataclass(frozen=True)
class TableSpec:
    label: str
    manual_id: str | None
    title: str
    region: str
    landscape: bool = False
    note: str | None = None
    list_title: str | None = None


TABLE_SPECS = [
    TableSpec(
        "tab:literature-comparisons",
        "1",
        "Physics-Based and AI-Based Docking Studies",
        "main",
    ),
    TableSpec(
        "tab:methods-calibration-characteristics",
        "M2",
        "Key Characteristics of the Calibration Benchmark Set",
        "main",
    ),
    TableSpec(
        "tab:methods-experimental-ligands",
        "1",
        "Experimental Orai1 Ligands",
        "main",
    ),
    TableSpec(
        "tab:results-pose-production",
        "R1",
        "Pose Production and Physical-Validity Summary of All Docking Variants",
        "main",
        True,
        r"\%* percentage of all produced poses",
    ),
    TableSpec(
        "tab:results-depth-recovery",
        "R2",
        "PoseBusters Pass-All Recovery by Ranked-Pool Depth",
        "main",
    ),
    TableSpec(
        "tab:results-depth-gain",
        "R4",
        "Decomposition of the Rank-1 to Best-of-Top-15 Pass-All Gain",
        "main",
    ),
    TableSpec(
        "tab:results-near-native-form",
        "R3",
        "Near-Nativeness and Form Recovery at Different Ranking Depths",
        "main",
        True,
    ),
    TableSpec(
        "tab:results-placement-form",
        "3",
        "Placement-versus-Form Distributions by Tool and Ranking Depth",
        "main",
        True,
    ),
    TableSpec(
        "tab:results-cluster-recovery",
        "R9",
        "Crystal-Cluster Reach and Co-Reach by Ranking Depth",
        "main",
    ),
    TableSpec(
        "tab:results-orai-yield",
        "R6",
        "Physical-Validity and Operational-Placement Yield for the Orai1 Experimental Ligands",
        "main",
    ),
    TableSpec(
        "tab:appendix-vina-weights",
        None,
        "AutoDock Vina Scoring-Function Weights Reproduced from Trott and Olson",
        "appendix",
    ),
    TableSpec(
        "tab:appendix-protocols",
        "M1",
        "Docking-Protocol Parameters for the Three Pipelines on the Calibration Benchmark",
        "appendix",
    ),
    TableSpec(
        "tab:appendix-posebusters-checks",
        "M3",
        "PoseBusters Validity Checks",
        "appendix",
    ),
    TableSpec(
        "tab:appendix-refinement-gains",
        "R5",
        "Post-Hoc Geometry-Optimisation Gains by Pose-Rank Band",
        "appendix",
    ),
    TableSpec(
        "tab:appendix-ranking-concordance",
        "R3",
        "Median Per-Complex Kendall Tau-b between Tool Ranking and Oracle Badness",
        "appendix",
    ),
    TableSpec(
        "tab:appendix-oracle-ranks",
        "R7",
        "Median Oracle Rank at Positions 1--15",
        "appendix",
        True,
    ),
    TableSpec(
        "tab:appendix-validity-by-rank",
        "R8",
        "PoseBusters Validity by Rank",
        "appendix",
        True,
    ),
    TableSpec(
        "tab:appendix-ligand-distribution",
        "1",
        "Distribution of Principal Ligand Descriptors in the PoseBusters Subset",
        "appendix",
    ),
    TableSpec(
        "tab:appendix-orai-descriptors",
        "2",
        "Physicochemical Descriptors of the Geometry-Optimised Orai1 Reference Ligands",
        "appendix",
        list_title="Orai1 Ligand Physicochemical Descriptors",
    ),
    TableSpec(
        "tab:appendix-top-k-recovery",
        "A1",
        "Top-k Recovery with Wilson 95\\% Confidence Intervals",
        "appendix",
        True,
    ),
    TableSpec(
        "tab:appendix-refinement-mcnemar",
        "A2",
        "Raw-to-gnina Refinement Effect on Near-Native Recovery",
        "appendix",
        True,
    ),
    TableSpec(
        "tab:appendix-cross-tool-mcnemar",
        "A3",
        "Cross-Tool Comparisons of Near-Native Recovery",
        "appendix",
        True,
    ),
    TableSpec(
        "tab:appendix-accurate-invalid",
        "A4",
        "Accurate-but-Invalid Share by Ranking Depth",
        "appendix",
        True,
    ),
]


REVIEW_QUERIES = [
    (
        "The calibration dataset is derived from the 308 high-resolution",
        "Please document the identities and exclusion reasons for the five source complexes not included in the 303-case shared analysis, and for the ten further cases without comparable timing records.",
    ),
    (
        "Two ligand panels pass through this pipeline.",
        "308 source ligands across four frames would yield 1,232 units. Please list the 30 failed or excluded units, or correct the reported total.",
    ),
    (
        r"\begin{landscape}" + "\n" + r"\scriptsize" + "\n" + r"\setlength{\tabcolsep}{2pt}" + "\n" + r"\begin{longtable}",
        "Table labelled tab:results-pose-production reports 7,337 AutoDock Vina poses, while the Appendix protocol reports 7,407. Reconcile both values against the canonical analysis files.",
    ),
    (
        r"\labelcurrenttable{tab:results-depth-recovery}",
        "Define the exact endpoint for the depth-recovery table and regenerate it from the canonical analysis. Its counts differ from the near-nativeness/form table.",
    ),
    (
        "The following tables give the full per-depth statistics underlying",
        "Appendix top-k recovery reports Vina values of 53.1% at rank 1 and 83.8% at top 15, whereas the main Results report 55.5% and 86.1%. Reconcile the denominator or variant and retain one canonical set.",
    ),
    (
        "For a fixed atom correspondence, rigid superposition gives",
        "Confirm from the analysis code that the same symmetry-equivalent atom mapping was retained for in-place and Kabsch form RMSD. Otherwise describe the exact decomposition as approximate.",
    ),
    (
        r"\labelcurrenttable{tab:appendix-orai-descriptors}",
        "The earlier ligand table and this Appendix table disagree on the identity and descriptors of 2abp-NH2 or 2-APB, Synta-66, and GSK-7975A. Regenerate both from the exact docked structures.",
    ),
]


SUPERSCRIPT_MAP = str.maketrans(
    {
        "⁰": "0",
        "¹": "1",
        "²": "2",
        "³": "3",
        "⁴": "4",
        "⁵": "5",
        "⁶": "6",
        "⁷": "7",
        "⁸": "8",
        "⁹": "9",
        "⁺": "+",
        "⁻": "-",
    }
)
SUBSCRIPT_MAP = str.maketrans(
    {
        "₀": "0",
        "₁": "1",
        "₂": "2",
        "₃": "3",
        "₄": "4",
        "₅": "5",
        "₆": "6",
        "₇": "7",
        "₈": "8",
        "₉": "9",
        "ᵢ": "i",
    }
)


def normalize_unicode(text: str) -> str:
    text = text.replace("W. Lu （陆威）", "W. Lu")
    text = text.replace(
        "√(in-place² − form²)",
        r"\(\sqrt{\text{in-place}^{2}-\text{form}^{2}}\)",
    )
    text = re.sub(
        r"[⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻]+",
        lambda match: r"\textsuperscript{" + match.group(0).translate(SUPERSCRIPT_MAP) + "}",
        text,
    )
    text = re.sub(
        r"[₀₁₂₃₄₅₆₇₈₉ᵢ]+",
        lambda match: r"\textsubscript{" + match.group(0).translate(SUBSCRIPT_MAP) + "}",
        text,
    )
    replacements = {
        "Å": r"\angstrom{}",
        "°": r"\ensuremath{^\circ}",
        "±": r"\ensuremath{\pm}",
        "µ": r"\ensuremath{\mu}",
        "μ": r"\ensuremath{\mu}",
        "·": r"\ensuremath{\cdot}",
        "×": r"\ensuremath{\times}",
        "Δ": r"\ensuremath{\Delta}",
        "α": r"\ensuremath{\alpha}",
        "β": r"\ensuremath{\beta}",
        "δ": r"\ensuremath{\delta}",
        "ε": r"\ensuremath{\varepsilon}",
        "π": r"\ensuremath{\pi}",
        "ρ": r"\ensuremath{\rho}",
        "τ": r"\ensuremath{\tau}",
        "φ": r"\ensuremath{\phi}",
        "χ": r"\ensuremath{\chi}",
        " ": r"\,",
        "→": r"\ensuremath{\rightarrow}",
        "−": r"\ensuremath{-}",
        "≈": r"\ensuremath{\approx}",
        "≤": r"\ensuremath{\leq}",
        "≥": r"\ensuremath{\geq}",
        "√": r"\ensuremath{\sqrt{\vphantom{x}}}",
    }
    for old, new in replacements.items():
        text = text.replace(old, new)
    return text


def convert_research_questions(text: str) -> str:
    pattern = re.compile(
        r"\\begin\{quote\}\n"
        r"How do the quality and stability.*?\n\n"
        r"How does post-docking optimisation.*?\n\n"
        r"How computationally efficient.*?\n\n"
        r"What opportunities and limitations.*?\n"
        r"\\end\{quote\}",
        re.DOTALL,
    )
    replacement = r"""\begin{enumerate}
\item[\textbf{RQ1}] How do the quality and stability of poses predicted by AutoDock Vina, DiffDock and EquiBind differ?
\item[\textbf{RQ2}] How does post-docking optimisation affect the quality and stability of predicted poses?
\item[\textbf{RQ3}] How computationally efficient are physics-based and AI-based docking approaches?
\item[\textbf{RQ4}] What opportunities and limitations do physics-based and AI-based docking tools present?
\end{enumerate}"""
    text, count = pattern.subn(lambda _: replacement, text, count=1)
    if count != 1:
        raise RuntimeError(f"Expected one research-question block, found {count}")
    return text


def figure_label(region: str, manual_id: str) -> str:
    token = manual_id.lower()
    if region == "main" and token == "1":
        return "fig:literature-docking-workflow"
    if region == "main":
        return f"fig:results-{token}"
    if token.startswith("r"):
        return f"fig:appendix-{token}"
    return f"fig:appendix-dataset-{token}"


def convert_figures(
    text: str, region: str, start_number: int
) -> tuple[str, list[tuple[str, str]], int]:
    pattern = re.compile(
        r"(?P<graphic>\\includegraphics\[[^\]]*\]\{[^}]+\})\n\n"
        r"(?P<caption>(?:Figure|Equation)\s+[^\n]+)"
    )
    refs: list[tuple[str, str]] = []
    counter = start_number

    def replacement(match: re.Match[str]) -> str:
        nonlocal counter
        counter += 1
        graphic = match.group("graphic")
        caption = match.group("caption").strip()
        image_match = re.search(r"image(\d+)\.png", graphic)
        if image_match is None:
            raise RuntimeError(f"Could not identify figure image in {graphic}")
        image_number = int(image_match.group(1))
        if caption.startswith("Equation "):
            if image_number != 26:
                raise RuntimeError(f"Unexpected equation screenshot image{image_number}")
            return (
                "% The rasterised duplicate of the Vina scoring equation was omitted; "
                "the native equations below are retained."
            )

        caption_match = re.match(
            r"Figure\s+([A-Z]?\d+[a-z]?)(?:\.\s*|\s*-\s*)(.*)", caption
        )
        if caption_match is None:
            raise RuntimeError(f"Could not parse figure caption: {caption}")
        manual_id = caption_match.group(1)
        caption_text = caption_match.group(2).strip().rstrip(".")
        label = figure_label(region, manual_id)
        refs.append((manual_id, label))
        relative_path = f"media/media/image{image_number}.png"
        return (
            "\\begin{figure}[!htbp]\n"
            "\\centering\n"
            f"\\includegraphics[width=\\linewidth,height=0.78\\textheight,"
            f"keepaspectratio]{{{relative_path}}}\n"
            f"\\caption{{{caption_text}}}\\label{{{label}}}\n"
            "\\end{figure}"
        )

    converted = pattern.sub(replacement, text)
    return converted, refs, counter


def convert_tables(
    text: str,
    specs: list[TableSpec],
) -> tuple[str, list[tuple[str, str]]]:
    pattern = re.compile(
        r"(?P<table>\\begin\{longtable\}.*?\\end\{longtable\})"
        r"(?:\n\n(?P<trailing>Table[^\n]*))?",
        re.DOTALL,
    )
    matches = list(pattern.finditer(text))
    if len(matches) != len(specs):
        raise RuntimeError(
            f"Expected {len(specs)} tables in {specs[0].region}, found {len(matches)}"
        )
    index = 0
    refs: list[tuple[str, str]] = []

    def replacement(match: re.Match[str]) -> str:
        nonlocal index
        spec = specs[index]
        index += 1
        block = match.group("table")
        # The Pandoc minipages are redundant inside fixed-width p columns and
        # become page-wide inside merged multicolumn cells. Replace them with
        # local ragged-right groups before adding the final caption.
        def simplify_minipage(minipage: re.Match[str]) -> str:
            content = minipage.group(1).strip()
            if not content:
                return "{}"
            if r"\\" in content:
                # A line break cannot split the argument of \textbf inside
                # shortstack's internal alignment. Repeat the formatting on
                # each line instead.
                content = re.sub(
                    r"\\textbf\{(.*?)\\\\\s*(.*?)\}(\\strut)?$",
                    lambda found: (
                        rf"\textbf{{{found.group(1)}}}\\"
                        rf"\textbf{{{found.group(2)}}}"
                        + (found.group(3) or "")
                    ),
                    content,
                    flags=re.DOTALL,
                )
                return r"\shortstack[l]{" + content + "}"
            return r"{\raggedright " + content + "}"

        block = block.replace(
            r"\begin{minipage}[b]{\linewidth}\raggedright", "<<<MP_BEGIN>>>"
        )
        block = block.replace(r"\end{minipage}", "<<<MP_END>>>")
        block = re.sub(
            r"<<<MP_BEGIN>>>(.*?)<<<MP_END>>>",
            simplify_minipage,
            block,
            flags=re.DOTALL,
        )
        if "MP_BEGIN" in block or "MP_END" in block or "minipage" in block:
            raise RuntimeError(f"Could not simplify all minipages in {spec.label}")
        block = re.sub(
            r"\\caption\{[^\n]*?\}\\tabularnewline\n?", "", block, count=1
        )
        # Pandoc places the caption inside the repeated longtable header.
        # Keep the caption there, but put the label after the environment so
        # it is evaluated only once when the table spans multiple pages.
        short_caption = (
            f"[{spec.list_title}]"
            if spec.list_title
            else ""
        )
        caption_line = (
            f"\\caption{short_caption}{{{spec.title}}}\\tabularnewline\n"
        )
        block, count = re.subn(
            r"(\\toprule\\noalign\{\}\n)",
            lambda found: caption_line + found.group(1),
            block,
            count=1,
        )
        if count != 1:
            raise RuntimeError(f"Could not insert caption in table {spec.label}")
        caption_start = block.find(caption_line)
        head_end = block.find(r"\endhead", caption_start)
        if caption_start < 0 or head_end < 0:
            raise RuntimeError(f"Could not locate longtable head in {spec.label}")
        first_head_end = block.find(
            r"\endfirsthead",
            caption_start + len(caption_line),
            head_end,
        )
        continued_caption = (
            f"\\caption[]{{{spec.title} (continued)}}\\tabularnewline\n"
        )
        if first_head_end >= 0:
            insertion = first_head_end + len(r"\endfirsthead")
            block = (
                block[:insertion]
                + "\n"
                + continued_caption
                + block[insertion:]
            )
        else:
            column_head = block[caption_start + len(caption_line) : head_end]
            first_and_repeated_head = (
                caption_line
                + column_head
                + "\\endfirsthead\n"
                + continued_caption
                + column_head
                + r"\endhead"
            )
            block = (
                block[:caption_start]
                + first_and_repeated_head
                + block[head_end + len(r"\endhead") :]
            )
        block += f"\n\\labelcurrenttable{{{spec.label}}}"
        if spec.note:
            block += f"\n\\par\\smallskip\\noindent {spec.note}"
        if spec.manual_id:
            refs.append((spec.manual_id, spec.label))
        if spec.landscape:
            block = (
                "\\begin{landscape}\n"
                "\\scriptsize\n"
                "\\setlength{\\tabcolsep}{2pt}\n"
                f"{block}\n"
                "\\end{landscape}"
            )
        return block

    converted = pattern.sub(replacement, text)
    converted = converted.replace(
        "\n\nScoring-function weights of AutoDock Vina, reproduced from "
        "Table 1 of Trott and Olson {[}22{]}.",
        "",
    )
    return converted, refs


def choose_reference(
    manual_id: str,
    candidates: list[tuple[str, str]],
    positions: dict[str, int],
    at: int,
) -> tuple[str, str | None] | None:
    exact = [(token, label) for token, label in candidates if token == manual_id]
    panel = None
    if not exact and re.fullmatch(r"(?:[RMA]?\d+)[A-Z]", manual_id):
        base = manual_id[:-1]
        exact = [(token, label) for token, label in candidates if token == base]
        panel = manual_id[-1]
    if not exact:
        return None
    token, label = min(exact, key=lambda item: abs(positions[item[1]] - at))
    return label, panel


def replace_cross_references(
    text: str,
    figure_refs: list[tuple[str, str]],
    table_refs: list[tuple[str, str]],
) -> str:
    positions = {
        label: text.find(rf"{{{label}}}")
        for _, label in figure_refs + table_refs
    }
    pattern = re.compile(
        r"\b(?P<kind>Figures?|Tables?)\s+"
        r"(?P<ids>[RMA]?\d+[a-zA-Z]?"
        r"(?:\s*(?:,|and|to|--|-)\s*[RMA]?\d+[a-zA-Z]?)*)(?!\w)",
        re.IGNORECASE,
    )

    def replacement(match: re.Match[str]) -> str:
        following = text[match.end() : match.end() + 5].lower()
        if following.startswith(" of "):
            return match.group(0)
        kind = match.group("kind")
        candidates = (
            figure_refs if kind.lower().startswith("figure") else table_refs
        )
        ids = match.group("ids")
        unresolved = False

        def id_replacement(id_match: re.Match[str]) -> str:
            nonlocal unresolved
            manual_id = id_match.group(0)
            chosen = choose_reference(
                manual_id, candidates, positions, match.start()
            )
            if chosen is None:
                unresolved = True
                return manual_id
            label, panel = chosen
            suffix = rf"\,({panel})" if panel else ""
            return rf"\ref{{{label}}}{suffix}"

        converted_ids = re.sub(r"[RMA]?\d+[a-zA-Z]?", id_replacement, ids)
        if unresolved:
            return match.group(0)
        normalized_kind = kind[0].upper() + kind[1:]
        return f"{normalized_kind}~{converted_ids}"

    return pattern.sub(replacement, text)


def convert_citations(text: str) -> str:
    def page_citation(match: re.Match[str]) -> str:
        number = int(match.group(1))
        page = match.group(2)
        return rf"\cite[p.~{page}]{{ref{number:03d}}}"

    text = re.sub(
        r"\{\[\}(\d+), p\. ([^{}]+)\{\]\}", page_citation, text
    )

    def range_citation(match: re.Match[str]) -> str:
        start = int(match.group(1))
        end = int(match.group(2))
        if not (1 <= start <= end <= 72):
            return match.group(0)
        keys = ",".join(f"ref{number:03d}" for number in range(start, end + 1))
        return rf"\cite{{{keys}}}"

    text = re.sub(
        r"\{\[\}(\d+)\{\]\}\s*(?:-|--)\s*\{\[\}(\d+)\{\]\}",
        range_citation,
        text,
    )

    def single_citation(match: re.Match[str]) -> str:
        number = int(match.group(1))
        if not 1 <= number <= 72:
            return match.group(0)
        return rf"\cite{{ref{number:03d}}}"

    text = re.sub(r"\{\[\}(\d+)\{\]\}", single_citation, text)
    return text


def convert_numbered_equations(text: str) -> str:
    pattern = re.compile(r"^\\\((.*?)\\\) \((\d+)\)$", re.MULTILINE)

    def replacement(match: re.Match[str]) -> str:
        equation = match.group(1)
        number = int(match.group(2))
        return (
            "\\begin{equation}\n"
            f"{equation}\n"
            f"\\label{{eq:vina-{number}}}\n"
            "\\end{equation}"
        )

    converted, count = pattern.subn(replacement, text)
    if count != 9:
        raise RuntimeError(f"Expected 9 numbered Vina equations, found {count}")
    converted = converted.replace(
        "Equation 1", r"Equation~\ref{eq:vina-1}"
    )
    return converted


def promote_appendix_headings(text: str) -> str:
    text, count = re.subn(
        r"\\hypertarget\{appendix\}\{%\n"
        r"\\chapter\{Appendix\}\\label\{appendix\}\}",
        r"\\appendix",
        text,
        count=1,
    )
    if count != 1:
        raise RuntimeError("Could not replace the Word Appendix wrapper")

    promoted: list[str] = []
    for line in text.splitlines():
        if re.match(r"\\subsubsection\{", line):
            line = line.replace(r"\subsubsection", r"\subsection", 1)
        elif re.match(r"\\subsection\{", line):
            line = line.replace(r"\subsection", r"\section", 1)
        elif re.match(r"\\section\{", line):
            line = line.replace(r"\section", r"\chapter", 1)
        promoted.append(line)
    return "\n".join(promoted)


def insert_review_queries(text: str) -> str:
    for anchor, query in REVIEW_QUERIES:
        marker = f"% REVIEW QUERY: {query}\n"
        if anchor not in text:
            raise RuntimeError(f"Review-query anchor not found: {anchor[:70]}")
        text = text.replace(anchor, marker + anchor, 1)
    return text


def clean_body(text: str) -> str:
    text = text.replace(
        "/mnt/c/Users/domin/OneDrive/Desktop/Documents/Claude/Projects/"
        "Starting a Business/thesis_latex/media/media/",
        "media/media/",
    )
    text = re.sub(
        r"(?m)^\\?%\* percentage of all produced poses\s*$",
        "",
        text,
    )
    text = text.replace(
        r"Threshold in\angstrom{}ngström",
        r"Threshold (\angstrom{})",
    )
    text = text.replace(
        "Threshold inÅngström",
        r"Threshold (\angstrom{})",
    )
    while re.search(r"\\cite\{[^{}]+\}\\cite\{[^{}]+\}", text):
        text = re.sub(
            r"\\cite\{([^{}]+)\}\\cite\{([^{}]+)\}",
            lambda match: rf"\cite{{{match.group(1)},{match.group(2)}}}",
            text,
        )
    text = text.replace(r"\begin{quote}", "").replace(r"\end{quote}", "")
    text = re.sub(r"\n{3,}", "\n\n", text)
    return normalize_unicode(text).strip() + "\n"


def bib_escape(text: str) -> str:
    text = text.replace("W. Lu （陆威）", "W. Lu")
    text = text.replace("\\", r"\textbackslash{}")
    for old, new in [
        ("&", r"\&"),
        ("%", r"\%"),
        ("_", r"\_"),
        ("#", r"\#"),
        ("{", r"\{"),
        ("}", r"\}"),
        ("$", r"\$"),
        ("~", r"\textasciitilde{}"),
        ("^", r"\textasciicircum{}"),
    ]:
        text = text.replace(old, new)
    text = text.replace("“", "``").replace("”", "''")
    text = text.replace("‘", "`").replace("’", "'")
    text = text.replace("–", "--").replace("—", "---")
    return normalize_unicode(text)


def build_bibliography() -> None:
    source = FULLTEXT.read_text(encoding="utf-8")
    references_text = source.split("# References", 1)[1].split("# Appendix", 1)[0]
    references = re.findall(
        r"^\[(\d+)\]\s+(.+?)(?=\n\n|\Z)", references_text, re.MULTILINE
    )
    if [int(number) for number, _ in references] != list(range(1, 73)):
        raise RuntimeError("Expected a continuous bibliography from [1] to [72]")
    entries = [
        "% Bibliography reconstructed from the 72 formatted references in the revised DOCX.",
        "% The custom misc driver in Thesis.tex prints each verified formatted entry verbatim.",
        "",
    ]
    for number, reference in references:
        key = f"ref{int(number):03d}"
        url_match = re.search(
            r"\s*\[Online\]\.\s*Available:\s*(https?://\S+?)\.?$",
            reference,
        )
        url = None
        if url_match:
            url = url_match.group(1)
            reference = reference[: url_match.start()].rstrip(" ,.")
        entries.extend(
            [
                f"@misc{{{key},",
                f"  sortkey = {{{int(number):03d}}},",
                f"  note = {{{bib_escape(reference)}}}"
                + ("," if url else ""),
                *( [f"  url = {{{url}}}"] if url else [] ),
                "}",
                "",
            ]
        )
    BIB_OUT.write_text("\n".join(entries), encoding="utf-8")


def build_main_tex() -> None:
    all_refs = ",".join(f"ref{number:03d}" for number in range(1, 73))
    thesis = rf"""\PassOptionsToPackage{{hyperfootnotes=false}}{{hyperref}}
\documentclass[Master,english,IEEE]{{twbook}}
\usepackage[utf8]{{inputenc}}
\usepackage[T1]{{fontenc}}
\usepackage{{longtable,booktabs,multirow,calc,etoolbox}}
\usepackage{{pdflscape}}
\usepackage{{footnotehyper}}
\usepackage{{microtype}}

\newcommand{{\angstrom}}{{\ensuremath{{\text{{\AA}}}}}}
\newcounter{{tablebeforeaitools}}
\makeatletter
\newcommand{{\labelcurrenttable}}[1]{{%
  \phantomsection
  \protected@edef\@currentlabel{{\thetable}}%
  \label{{#1}}}}
\makeatother
\providecommand{{\tightlist}}{{%
  \setlength{{\itemsep}}{{0pt}}\setlength{{\parskip}}{{0pt}}}}

\degreecourse{{Artificial Intelligence Engineering}}
\addbibresource{{Literatur.bib}}
\ExecuteBibliographyOptions{{sorting=none}}
\DeclareBibliographyDriver{{misc}}{{%
  \printfield{{note}}%
  \newunit\newblock%
  \printfield{{url}}%
  \finentry}}

\title{{Comparative Analysis of Physics-Based and AI-Based Molecular Docking Tools}}
\author{{Dominik Mann, PhD BSc}}
\studentnumber{{2410585017}}
\supervisor{{FH-Prof. Priv.-Doz. Dipl.-Ing.(FH) Dr.\\Bernhard Knapp}}
\place{{Wien}}

% TODO(author): The source DOCX contains no usable abstract or keyword list.
% Add \outline{{...}}, \keywords{{...}}, \kurzfassung{{...}}, and
% \schlagworte{{...}} here if required before submission.

\begin{{document}}

% Cite all sources in their reviewed numeric order before the first body citation.
\nocite{{{all_refs}}}
\maketitle

\input{{body_main}}

% Transparent disclosure of the assistance used for this revision and conversion.
\aitoolentry{{OpenAI Codex}}{{Consistency review only; not used for drafting any
section of the thesis}}{{Prompt: ``Please conduct a very detailed review of the
thesis. Compare the numbers reported in the tables and figures and make sure
they are in line with the text.'' User requests dated 25--26 July 2026; entire
document}}

\aitoolentry{{Claude Code}}{{Code review and code documentation of the analysis
and docking scripts}}{{Representative prompts: ``Review the code for errors and
bugs''; ``Please add a function description and comments.'' User requests dated
25--26 July 2026; analysis and docking source code}}

\clearpage
\printbibliography
\clearpage
\listoffigures
\clearpage
\listoftables
\clearpage
\setcounter{{tablebeforeaitools}}{{\value{{table}}}}
\listaitools
\setcounter{{table}}{{\value{{tablebeforeaitools}}}}

\clearpage
\input{{body_appendix}}

\end{{document}}
"""
    THESIS_OUT.write_text(thesis, encoding="utf-8")


def build_readme() -> None:
    README_OUT.write_text(
        """# LaTeX thesis project

This project converts the consistency-reviewed Word thesis into the supplied
FH Technikum Wien `twbook` template.

Build in the project directory with:

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error -file-line-error Thesis.tex
```

The main document is `Thesis.tex`; the generated content is split into
`body_main.tex` and `body_appendix.tex`. The 72 reviewed references are in
`Literatur.bib`, and extracted figures are under `media/media/`.

Seven non-printing `% REVIEW QUERY:` comments preserve unresolved source-data
questions from the Word review. The eighth editorial query (mixed figure/table
numbering) is resolved here through automatic LaTeX numbering and references.
The original Pandoc conversion is retained as `body_raw.tex` for auditability.
""",
        encoding="utf-8",
    )


def validate(
    main: str,
    appendix: str,
    main_figures: list[tuple[str, str]],
    appendix_figures: list[tuple[str, str]],
) -> None:
    combined = main + "\n" + appendix
    assertions = {
        "figures": (combined.count(r"\begin{figure}"), 42),
        "tables": (combined.count(r"\begin{longtable}"), 23),
        "footnotes": (combined.count(r"\footnote{"), 4),
        "numbered Vina equations": (combined.count(r"\label{eq:vina-"), 9),
        "main chapters": (main.count(r"\chapter{"), 7),
        "appendix chapters": (appendix.count(r"\chapter{"), 7),
        "figure labels": (len(main_figures) + len(appendix_figures), 42),
        "review queries": (combined.count("% REVIEW QUERY:"), 7),
    }
    failures = [
        f"{name}: got {actual}, expected {expected}"
        for name, (actual, expected) in assertions.items()
        if actual != expected
    ]
    if re.search(r"\{\[\}\d+\{\]\}", combined):
        failures.append("unconverted numeric citation token remains")
    if "Physical-chemical" in combined or "physical-chemistry-based" in combined:
        failures.append("obsolete physics-based terminology remains")
    if r"\chapter{References}" in combined:
        failures.append("manual References chapter remains")
    if r"\chapter{Appendix}" in combined:
        failures.append("generic Appendix chapter remains")
    if "image26.png" in combined:
        failures.append("duplicate scoring-equation screenshot remains")
    if failures:
        raise RuntimeError("Validation failed:\n- " + "\n- ".join(failures))


def main() -> None:
    raw = RAW.read_text(encoding="utf-8").replace("\r\n", "\n")
    motivation = raw.index(r"\hypertarget{motivation}")
    references = raw.index(r"\hypertarget{references}")
    appendix = raw.index(r"\hypertarget{appendix}")
    main = raw[motivation:references]
    appendix_text = raw[appendix:]

    main = convert_research_questions(main)
    appendix_text = promote_appendix_headings(appendix_text)

    main, main_figures, figure_counter = convert_figures(main, "main", 0)
    appendix_text, appendix_figures, figure_counter = convert_figures(
        appendix_text, "appendix", figure_counter
    )
    if figure_counter != 43:
        raise RuntimeError(f"Expected 43 image instances, found {figure_counter}")

    main_specs = [spec for spec in TABLE_SPECS if spec.region == "main"]
    appendix_specs = [spec for spec in TABLE_SPECS if spec.region == "appendix"]
    main, main_tables = convert_tables(main, main_specs)
    appendix_text, appendix_tables = convert_tables(
        appendix_text, appendix_specs
    )

    main = replace_cross_references(main, main_figures, main_tables)
    appendix_text = replace_cross_references(
        appendix_text, appendix_figures, appendix_tables
    )
    main = convert_citations(main)
    appendix_text = convert_citations(appendix_text)
    appendix_text = convert_numbered_equations(appendix_text)

    main = clean_body(main)
    appendix_text = clean_body(appendix_text)
    combined = insert_review_queries(main + "\n% SPLIT APPENDIX\n" + appendix_text)
    main, appendix_text = combined.split("\n% SPLIT APPENDIX\n", 1)
    main = main.strip() + "\n"
    appendix_text = (
        "% REVIEW NOTE: The Word document's mixed manual figure/table numbering "
        "is resolved here with automatic global numbering and cross-references.\n"
        + appendix_text.strip()
        + "\n"
    )

    validate(main, appendix_text, main_figures, appendix_figures)
    MAIN_OUT.write_text(main, encoding="utf-8")
    APPENDIX_OUT.write_text(appendix_text, encoding="utf-8")
    build_bibliography()
    build_main_tex()
    build_readme()
    print(f"Wrote {THESIS_OUT}")
    print(f"Wrote {MAIN_OUT}")
    print(f"Wrote {APPENDIX_OUT}")
    print(f"Wrote {BIB_OUT}")
    print("Validated 7 main chapters, 7 appendix chapters, 42 figures, 23 tables,")
    print("4 anchored footnotes, 9 numbered equations, 72 references, and 7 review queries.")


if __name__ == "__main__":
    main()
