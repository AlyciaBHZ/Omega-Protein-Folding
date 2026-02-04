"""
Import figures/tables from evidence ZIPs into paper/assets/.

This script is intentionally heuristic: it tries to map common naming conventions
to the LaTeX placeholders generated from the DOCX conversion:

  paper/assets/figures/fig1.(png|pdf|jpg)
  paper/assets/tables/table1.tex

Usage (PowerShell):
  python tools/import_zip_assets.py --zip path\\to\\Omega_v1_evidence_repository.zip --paper paper
  python tools/import_zip_assets.py --zip path\\to\\sprint_to_0p9_outputs.zip --paper paper
"""

from __future__ import annotations

import argparse
import csv
import re
import shutil
import zipfile
from pathlib import Path
from typing import Optional


FIG_EXTS = {".png", ".jpg", ".jpeg", ".pdf"}
TABLE_EXTS = {".tex", ".csv", ".tsv"}

FIG_NO_RE = re.compile(r"(?:^|[^a-z0-9])fig(?:ure)?[_\- ]*(\d+)(?:[^0-9]|$)", re.IGNORECASE)
TAB_NO_RE = re.compile(r"(?:^|[^a-z0-9])tab(?:le)?[_\- ]*(\d+)(?:[^0-9]|$)", re.IGNORECASE)


def _safe_copy(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)


def _csv_to_table_tex(src: Path, dst: Path) -> None:
    # A minimal, robust conversion: everything becomes a string cell.
    # Users can refine formatting later.
    delimiter = "\t" if src.suffix.lower() == ".tsv" else ","
    with src.open("r", encoding="utf-8", newline="") as f:
        rows = list(csv.reader(f, delimiter=delimiter))

    if not rows:
        dst.write_text("% Empty table\n", encoding="utf-8")
        return

    ncols = max(len(r) for r in rows)
    colspec = "l" * ncols

    def esc(cell: str) -> str:
        # Keep this very simple; main.tex already loads booktabs.
        return (
            cell.replace("\\", r"\textbackslash{}")
            .replace("&", r"\&")
            .replace("%", r"\%")
            .replace("$", r"\$")
            .replace("#", r"\#")
            .replace("_", r"\_")
            .replace("{", r"\{")
            .replace("}", r"\}")
            .replace("~", r"\textasciitilde{}")
            .replace("^", r"\textasciicircum{}")
        )

    lines = []
    lines.append(r"\begin{tabular}{" + colspec + "}")
    lines.append(r"\toprule")
    header = rows[0] + [""] * (ncols - len(rows[0]))
    lines.append(" & ".join(esc(c) for c in header) + r" \\")
    lines.append(r"\midrule")
    for r in rows[1:]:
        rr = r + [""] * (ncols - len(r))
        lines.append(" & ".join(esc(c) for c in rr) + r" \\")
    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    lines.append("")

    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text("\n".join(lines), encoding="utf-8")


def _extract_zip(zip_path: Path, tmp_dir: Path) -> None:
    if tmp_dir.exists():
        shutil.rmtree(tmp_dir)
    tmp_dir.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(zip_path, "r") as zf:
        zf.extractall(tmp_dir)


def _find_number(pattern: re.Pattern[str], name: str) -> Optional[int]:
    m = pattern.search(name)
    if not m:
        return None
    try:
        return int(m.group(1))
    except ValueError:
        return None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--zip", required=True, type=Path)
    ap.add_argument("--paper", required=True, type=Path)
    args = ap.parse_args()

    zip_path: Path = args.zip
    paper_dir: Path = args.paper

    assets_fig = paper_dir / "assets" / "figures"
    assets_tab = paper_dir / "assets" / "tables"
    tmp_dir = paper_dir / "_import" / zip_path.stem

    _extract_zip(zip_path, tmp_dir)

    # 1) Figures
    for p in tmp_dir.rglob("*"):
        if not p.is_file():
            continue
        ext = p.suffix.lower()
        if ext not in FIG_EXTS:
            continue
        n = _find_number(FIG_NO_RE, p.name)
        if n is None:
            # Keep original name, but normalize to lower-case extension.
            dst = assets_fig / p.name
        else:
            # Prefer png when duplicates exist: do not overwrite a png with a pdf.
            dst = assets_fig / f"fig{n}{ext}"
        if dst.exists() and dst.suffix.lower() == ".png" and ext == ".pdf":
            continue
        _safe_copy(p, dst)

    # 2) Tables
    for p in tmp_dir.rglob("*"):
        if not p.is_file():
            continue
        ext = p.suffix.lower()
        if ext not in TABLE_EXTS:
            continue
        n = _find_number(TAB_NO_RE, p.name)
        if n is None:
            continue
        if ext == ".tex":
            _safe_copy(p, assets_tab / f"table{n}.tex")
        elif ext in {".csv", ".tsv"}:
            _csv_to_table_tex(p, assets_tab / f"table{n}.tex")

    print(f"Imported assets from: {zip_path}")
    print(f"Figures -> {assets_fig}")
    print(f"Tables  -> {assets_tab}")


if __name__ == "__main__":
    main()

