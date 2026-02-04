"""
DOCX -> sectioned LaTeX paper/ scaffold.

Goals:
- Split the manuscript into one .tex per section (Heading 2) + subsections (Heading 3).
- Convert inline figure/table callouts like "Figure 3 | ..." into LaTeX floats that
  reference assets under paper/assets/{figures,tables}/.
- Keep output compilable even when images/tables are missing (prints a boxed placeholder).

Usage (PowerShell):
  python tools/docx_to_paper_latex.py --docx Omega_v1_NatureStyle_EN.docx --out paper
"""

from __future__ import annotations

import argparse
import re
import unicodedata
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Optional, Tuple

from docx import Document  # type: ignore


UNICODE_LATEX_INLINE = {
    "β": r"$\beta$",
    "θ": r"$\theta$",
    "★": r"$\star$",
    "Δ": r"$\Delta$",
    "≤": r"$\le$",
    "≥": r"$\ge$",
    "→": r"$\to$",
    "×": r"$\times$",
}

UNICODE_SEQUENCE_LATEX_INLINE = {
    # Common composed symbol in this manuscript
    "θ★": r"$\theta^\star$",
    # Avoid split math blocks in prose
    "Δw₀": r"$\Delta w_0$",
}

CODE_UNICODE_ASCII = {
    "β": "beta",
    "θ": "theta",
    "★": "star",
    "Δ": "Delta",
    "≤": "<=",
    "≥": ">=",
    "→": "->",
    "×": "x",
    "≈": "~=",
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
    "ⱼ": "j",
    "ₖ": "k",
}

LATEX_SPECIALS = {
    "\\": r"\textbackslash{}",
    "&": r"\&",
    "%": r"\%",
    "$": r"\$",
    "#": r"\#",
    "_": r"\_",
    "{": r"\{",
    "}": r"\}",
    "~": r"\textasciitilde{}",
    "^": r"\textasciicircum{}",
}


def _normalize_ws(s: str) -> str:
    s = s.replace("\u00a0", " ")
    s = re.sub(r"[ \t]+", " ", s)
    return s.strip()


def escape_latex_text(s: str) -> str:
    # Convert inline code spans first so we do NOT inject math inside \texttt{...}.
    def _code_to_texttt(m: re.Match[str]) -> str:
        code = m.group(1)
        for ch, rep in CODE_UNICODE_ASCII.items():
            code = code.replace(ch, rep)
        # Escape specials for LaTeX texttt context.
        code = "".join(LATEX_SPECIALS.get(ch, ch) for ch in code)
        return rf"\texttt{{{code}}}"

    s = re.sub(r"`([^`]+)`", _code_to_texttt, s)

    # Convert composed Unicode sequences before single-character mapping.
    for seq, rep in UNICODE_SEQUENCE_LATEX_INLINE.items():
        s = s.replace(seq, rep)

    # Basic Unicode inline conversions for normal text.
    for ch, rep in UNICODE_LATEX_INLINE.items():
        s = s.replace(ch, rep)

    # Convert markdown-ish emphasis that appears in the DOCX text itself.
    # NOTE: This is conservative and only handles **...**.
    s = re.sub(r"\*\*(.+?)\*\*", r"\\textbf{\1}", s)

    # A few manuscript-specific unicode subscripts/suffixes used in normal prose.
    # Prefer explicit math where possible (e.g., w₀ -> $w_0$).
    s = re.sub(r"([A-Za-z])₀", r"$\1_0$", s)
    s = re.sub(r"([A-Za-z])₁", r"$\1_1$", s)
    s = re.sub(r"([A-Za-z])₂", r"$\1_2$", s)
    s = re.sub(r"([A-Za-z])₃", r"$\1_3$", s)
    s = re.sub(r"([A-Za-z])₄", r"$\1_4$", s)
    s = re.sub(r"([A-Za-z])₅", r"$\1_5$", s)
    s = re.sub(r"([A-Za-z])₆", r"$\1_6$", s)
    s = re.sub(r"([A-Za-z])₇", r"$\1_7$", s)
    s = re.sub(r"([A-Za-z])₈", r"$\1_8$", s)
    s = re.sub(r"([A-Za-z])₉", r"$\1_9$", s)
    s = re.sub(r"([A-Za-z])ᵢ", r"$\1_i$", s)
    s = re.sub(r"([A-Za-z])ⱼ", r"$\1_j$", s)
    s = re.sub(r"([A-Za-z])ₖ", r"$\1_k$", s)
    s = s.replace("≈", r"$\approx$")

    # Escape LaTeX specials (but do not double-escape inside our injected commands).
    #
    # Strategy:
    # - Protect math spans first, then LaTeX commands we injected.
    # - Escape the remaining raw text.
    # - Restore protected spans (iteratively, in case any protection introduced nested tokens).
    protected: List[str] = []

    def _protect(m: re.Match[str]) -> str:
        protected.append(m.group(0))
        return f"@@PROT{len(protected)-1}@@"

    s = re.sub(r"\$[^$]*\$", _protect, s)
    s = re.sub(r"\\[a-zA-Z]+\{[^{}]*\}", _protect, s)
    s = re.sub(r"\\[a-zA-Z]+", _protect, s)

    out_chars: List[str] = []
    for ch in s:
        out_chars.append(LATEX_SPECIALS.get(ch, ch))
    s = "".join(out_chars)

    def _restore(m: re.Match[str]) -> str:
        idx = int(m.group(1))
        return protected[idx]

    # Restore (may need multiple passes if some protected chunks contained protected tokens).
    while re.search(r"@@PROT(\d+)@@", s):
        s = re.sub(r"@@PROT(\d+)@@", _restore, s)
    return s


def slugify(s: str) -> str:
    s = _normalize_ws(s).lower()
    s = unicodedata.normalize("NFKD", s)
    s = "".join(ch for ch in s if not unicodedata.combining(ch))
    s = re.sub(r"[^a-z0-9]+", "_", s).strip("_")
    if not s:
        return "section"
    return s


@dataclass
class Para:
    style: str
    text: str


FIG_RE = re.compile(r"^\s*Figure\s+(\d+)\s*\|\s*(.+?)\s*$", re.IGNORECASE)
TAB_RE = re.compile(r"^\s*Table\s+(\d+)\s*\|\s*(.+?)\s*$", re.IGNORECASE)


def to_latex_figure(fig_no: int, caption: str) -> str:
    caption = escape_latex_text(caption)
    return "\n".join(
        [
            r"\begin{figure}[t]",
            r"  \centering",
            rf"  \maybeincludegraphics[width=\linewidth]{{assets/figures/fig{fig_no}}}",
            rf"  \caption{{{caption}}}",
            rf"  \label{{fig:{fig_no}}}",
            r"\end{figure}",
            "",
        ]
    )


def to_latex_table(tab_no: int, caption: str) -> str:
    caption = escape_latex_text(caption)
    return "\n".join(
        [
            r"\begin{table}[t]",
            r"  \centering",
            rf"  \caption{{{caption}}}",
            rf"  \label{{tab:{tab_no}}}",
            rf"  \maybeincludetable{{assets/tables/table{tab_no}.tex}}",
            r"\end{table}",
            "",
        ]
    )


def iter_doc_paras(doc: Document) -> Iterable[Para]:
    for p in doc.paragraphs:
        text = _normalize_ws(p.text)
        if not text:
            continue
        style = getattr(getattr(p, "style", None), "name", "") or ""
        yield Para(style=style, text=text)


def build_sections(paras: List[Para]) -> Tuple[dict, List[str]]:
    """
    Returns:
      - files: mapping of section_file -> latex content
      - ordered_inputs: list of section files to \input in order
    """
    files: dict[str, List[str]] = {}
    ordered: List[str] = []

    title_lines: List[str] = []
    i = 0
    while i < len(paras) and not paras[i].style.startswith("Heading"):
        title_lines.append(paras[i].text)
        i += 1

    # Abstract block
    abstract_lines: List[str] = []
    if i < len(paras) and paras[i].style.startswith("Heading") and paras[i].text.lower() == "abstract":
        i += 1
        while i < len(paras):
            if paras[i].style.startswith("Heading"):
                break
            abstract_lines.append(paras[i].text)
            i += 1

    frontmatter = []
    if title_lines:
        frontmatter.append(rf"\title{{{escape_latex_text(title_lines[0])}}}")
        if len(title_lines) > 1:
            frontmatter.append(rf"\subtitle{{{escape_latex_text(title_lines[1])}}}")
        author = ""
        date = ""
        for line in title_lines[2:]:
            if "version:" in line.lower() or "date:" in line.lower():
                date = line
            elif not author:
                author = line
        if author:
            frontmatter.append(rf"\author{{{escape_latex_text(author)}}}")
        else:
            frontmatter.append(r"\author{}")
        if date:
            frontmatter.append(rf"\date{{{escape_latex_text(date)}}}")
        else:
            frontmatter.append(r"\date{}")
    else:
        frontmatter.append(r"\title{}")
        frontmatter.append(r"\author{}")
        frontmatter.append(r"\date{}")

    frontmatter.append(r"\maketitle")
    if abstract_lines:
        frontmatter.append(r"\begin{abstract}")
        frontmatter.append(escape_latex_text(" ".join(abstract_lines)))
        frontmatter.append(r"\end{abstract}")
        frontmatter.append("")

    front_path = "sections/00_frontmatter.tex"
    files[front_path] = frontmatter
    ordered.append(front_path)

    current_file: Optional[str] = None
    current_lines: List[str] = []

    def flush():
        nonlocal current_file, current_lines
        if current_file is None:
            return
        files[current_file] = list(current_lines)
        current_file = None
        current_lines = []

    # Skip "Main text" Heading 1 if present.
    while i < len(paras):
        p = paras[i]
        if p.style.startswith("Heading 1") and p.text.lower() == "main text":
            i += 1
            continue
        if p.style.startswith("Heading 1") and "supplementary information" in p.text.lower():
            flush()
            current_file = "sections/90_supplementary.tex"
            ordered.append(current_file)
            current_lines.append(r"\section*{Supplementary Information}")
            current_lines.append(r"\addcontentsline{toc}{section}{Supplementary Information}")
            current_lines.append("")
            i += 1
            continue
        if p.style.startswith("Heading 2"):
            flush()
            slug = slugify(p.text)
            # Ensure stable numbering order.
            num = len([x for x in ordered if x.startswith("sections/0") and x.endswith(".tex")])  # includes frontmatter
            # Frontmatter counts as 0, so start main sections at 1.
            seq = max(1, num)
            current_file = f"sections/{seq:02d}_{slug}.tex"
            ordered.append(current_file)
            current_lines.append(rf"\section{{{escape_latex_text(p.text)}}}")
            current_lines.append("")
            i += 1
            continue
        if p.style.startswith("Heading 3"):
            current_lines.append(rf"\subsection{{{escape_latex_text(p.text)}}}")
            current_lines.append("")
            i += 1
            continue

        # Normal paragraph: detect figure/table caption lines.
        m_fig = FIG_RE.match(p.text)
        if m_fig:
            current_lines.append(to_latex_figure(int(m_fig.group(1)), m_fig.group(2)))
            i += 1
            continue
        m_tab = TAB_RE.match(p.text)
        if m_tab:
            current_lines.append(to_latex_table(int(m_tab.group(1)), m_tab.group(2)))
            i += 1
            continue

        current_lines.append(escape_latex_text(p.text))
        current_lines.append("")
        i += 1

    flush()

    # Ensure references section exists (if DOCX had a heading but no bib).
    if not any("references" in p.lower() for p in ordered):
        ref_path = "sections/80_references.tex"
        files[ref_path] = [
            r"\section{References}",
            r"% TODO: add BibTeX or manual references here.",
            r"\bibliographystyle{naturemag}",
            r"\bibliography{references}",
            "",
        ]
        ordered.append(ref_path)

    return {k: "\n".join(v).rstrip() + "\n" for k, v in files.items()}, ordered


MAIN_TEX = r"""\documentclass[11pt]{article}

% Nature-style single-column submission scaffold (compile with XeLaTeX recommended).
% !TEX program = xelatex

\usepackage[a4paper,margin=1in]{geometry}
\usepackage{graphicx}
\usepackage{booktabs}
\usepackage{etoolbox}
\usepackage{hyperref}
\usepackage{xcolor}
\usepackage{caption}
\usepackage{subcaption}

% A simple \subtitle implementation
\makeatletter
\providecommand{\subtitle}[1]{%
  \apptocmd{\@title}{\par\large #1 \par}{}{}%
}
\makeatother

% Include figure if assets exist; otherwise show placeholder box.
\newcommand{\maybeincludegraphics}[2][]{%
  \IfFileExists{#2.pdf}{\includegraphics[#1]{#2.pdf}}{%
  \IfFileExists{#2.png}{\includegraphics[#1]{#2.png}}{%
  \IfFileExists{#2.jpg}{\includegraphics[#1]{#2.jpg}}{%
  \fbox{Missing file: #2.(pdf|png|jpg)}%
  }}}%
}

\newcommand{\maybeincludetable}[1]{%
  \IfFileExists{#1}{\input{#1}}{\fbox{Missing file: #1}}%
}

\begin{document}

% Auto-generated from Omega_v1_NatureStyle_EN.docx
% Do not edit this file by hand if you plan to re-run the converter.

% INPUTS_START
% INPUTS_END

\end{document}
"""


def write_main_tex(out_dir: Path, ordered_inputs: List[str]) -> None:
    lines = MAIN_TEX.splitlines()
    start = lines.index("% INPUTS_START")
    end = lines.index("% INPUTS_END")
    inputs = [rf"\input{{{p.replace('sections/', 'sections/')}}}" for p in ordered_inputs]
    new_lines = lines[: start + 1] + inputs + lines[end:]
    (out_dir / "main.tex").write_text("\n".join(new_lines) + "\n", encoding="utf-8")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--docx", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    docx_path: Path = args.docx
    out_dir: Path = args.out
    sections_dir = out_dir / "sections"
    sections_dir.mkdir(parents=True, exist_ok=True)

    doc = Document(str(docx_path))
    paras = list(iter_doc_paras(doc))

    files, ordered = build_sections(paras)
    for rel, content in files.items():
        path = out_dir / rel
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")

    write_main_tex(out_dir, ordered)


if __name__ == "__main__":
    main()

