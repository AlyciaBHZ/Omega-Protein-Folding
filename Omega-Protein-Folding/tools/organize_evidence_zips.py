"""
Organize extracted evidence ZIPs into the repo structure.

Inputs (already extracted by user/agent):
  data/raw/Omega_v1_evidence_repository/
  data/raw/sprint_to_0p9_outputs/

Outputs:
  - paper/assets/figures/fig{N}.png  (mapped from known filenames)
  - paper/assets/tables/table{N}.tex (generated from CSVs)
  - scripts/evidence_v1/             (copy python + requirements/specs)
  - data/processed/                  (copy CSVs used by paper)
  - artifacts/reports/               (copy markdown reports)
"""

from __future__ import annotations

import csv
import shutil
from pathlib import Path
from typing import Iterable, List, Tuple


REPO = Path(__file__).resolve().parents[1]
RAW_EVID = REPO / "data" / "raw" / "Omega_v1_evidence_repository"
RAW_SPRINT = REPO / "data" / "raw" / "sprint_to_0p9_outputs"

PAPER = REPO / "paper"
FIG_DIR = PAPER / "assets" / "figures"
TAB_DIR = PAPER / "assets" / "tables"

SCRIPTS_DST = REPO / "scripts" / "evidence_v1"
REPORTS_DST = REPO / "artifacts" / "reports"
DATA_DST = REPO / "data" / "processed"
SPRINT_DST = DATA_DST / "sprint_to_0p9_outputs"


FIG_MAP: List[Tuple[int, str]] = [
    (1, "omega_overview_diagram.png"),
    (2, "imp_operator_usage_bar.png"),
    (3, "ab_hallucination_gap_tmscore_extalpha_q08.png"),
    (4, "omega_projector_L3_1EMA_phratio_vs_tmscore.png"),
    (5, "omega_projector_1EMA_seed0_phason_trajectory.png"),
    (6, "theta_star_AB_phmax_1AKE.png"),
    (7, "theta_star_AB_phmax_1TIM.png"),
    (8, "triple232_oracle_codec_tm_bar.png"),
    # Sprint-to-0.9 figures (prefer evidence pack if present; fall back to sprint zip)
    (9, "sprint_tmscore_distribution.png"),
    (10, "sprint_phason_theta_scatter.png"),
    (11, "sprint_1TIM_phason_strain_trace.png"),
    (12, "sprint_2O5P_hydrophobic_belt_scatter.png"),
]

TABLE_MAP: List[Tuple[int, str]] = [
    # Table 1: impedance-guided benchmarks summary (proxy N=16)
    (1, "imp_adaptive_v2_benchmark_summary.csv"),
    # Table 2: theta-star feasibility A/B
    (2, "theta_star_AB_success_1AKE_1TIM.csv"),
    # Table 3: sprint-to-0.9 A/B
    (3, "sprint_to_0p9_AB_summary.csv"),
]


def _mkdir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def _copy(src: Path, dst: Path) -> None:
    _mkdir(dst.parent)
    shutil.copy2(src, dst)


def _copy_no_overwrite(src: Path, dst: Path) -> Path:
    """
    Windows file systems are case-insensitive, so names like AB_summary vs ab_summary
    can collide. This helper avoids clobbering an existing file: it appends a suffix
    and returns the actual destination used.
    """
    _mkdir(dst.parent)
    if not dst.exists():
        shutil.copy2(src, dst)
        return dst
    stem, suf = dst.stem, dst.suffix
    for k in range(1, 1000):
        alt = dst.with_name(f"{stem}__dup{k}{suf}")
        if not alt.exists():
            shutil.copy2(src, alt)
            return alt
    raise RuntimeError(f"Could not find non-colliding destination for {dst}")


def _read_csv_rows(path: Path) -> List[List[str]]:
    # Try comma, then tab.
    raw = path.read_text(encoding="utf-8")
    for delim in [",", "\t"]:
        rows = list(csv.reader(raw.splitlines(), delimiter=delim))
        if len(rows) >= 1 and max(len(r) for r in rows) >= 2:
            return rows
    # fallback: single column
    return [[x] for x in raw.splitlines() if x.strip()]


def _escape_cell(cell: str) -> str:
    # Minimal LaTeX escaping (tabular context)
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


def csv_to_booktabs_table(csv_path: Path, tex_path: Path, max_rows: int = 25) -> None:
    rows = _read_csv_rows(csv_path)
    if not rows:
        tex_path.write_text("% Empty table\n", encoding="utf-8")
        return

    ncols = max(len(r) for r in rows)
    # Cap rows for readability in the main paper; full CSV remains in data/processed.
    rows = rows[: max_rows + 1]  # header + rows

    def pad(r: List[str]) -> List[str]:
        return r + [""] * (ncols - len(r))

    colspec = "l" * ncols
    out: List[str] = []
    out.append(r"\begin{tabular}{" + colspec + "}")
    out.append(r"\toprule")
    hdr = pad(rows[0])
    out.append(" & ".join(_escape_cell(c) for c in hdr) + r" \\")
    out.append(r"\midrule")
    for r in rows[1:]:
        rr = pad(r)
        out.append(" & ".join(_escape_cell(c) for c in rr) + r" \\")
    out.append(r"\bottomrule")
    out.append(r"\end{tabular}")
    out.append("")

    _mkdir(tex_path.parent)
    tex_path.write_text("\n".join(out), encoding="utf-8")


def first_existing(paths: Iterable[Path]) -> Path:
    for p in paths:
        if p.exists():
            return p
    raise FileNotFoundError("None of the candidate paths exist: " + ", ".join(str(p) for p in paths))


def main() -> None:
    # 0) Ensure dirs
    _mkdir(FIG_DIR)
    _mkdir(TAB_DIR)
    _mkdir(SCRIPTS_DST)
    _mkdir(REPORTS_DST)
    _mkdir(DATA_DST)
    _mkdir(SPRINT_DST)

    # 1) Copy scripts/specs for reproducibility
    for p in (RAW_EVID / "scripts").glob("*"):
        _copy(p, SCRIPTS_DST / p.name)
    for p in (RAW_EVID / "specs").glob("*"):
        _copy(p, SCRIPTS_DST / "specs" / p.name)

    # 2) Copy reports
    if (RAW_EVID / "reports").exists():
        for p in (RAW_EVID / "reports").glob("*.md"):
            _copy(p, REPORTS_DST / p.name)
    if (RAW_SPRINT / "sprint_to_0p9_report.md").exists():
        _copy(RAW_SPRINT / "sprint_to_0p9_report.md", REPORTS_DST / "sprint_to_0p9_report.md")

    # 3) Map figures into LaTeX expected names
    for fig_no, src_name in FIG_MAP:
        src = first_existing(
            [
                RAW_EVID / "figures" / src_name,
                RAW_SPRINT / src_name,
            ]
        )
        _copy(src, FIG_DIR / f"fig{fig_no}{src.suffix.lower()}")

    # 4) Copy CSVs used by paper + generate table tex
    for tab_no, csv_name in TABLE_MAP:
        src = first_existing(
            [
                RAW_EVID / "data" / csv_name,
                RAW_EVID / "data" / csv_name.lower(),
                RAW_EVID / "data" / csv_name.replace("0p9", "to_0p9"),  # tolerate variants
                RAW_SPRINT / csv_name,
                RAW_SPRINT / csv_name.lower(),
            ]
        )
        _copy(src, DATA_DST / csv_name)
        csv_to_booktabs_table(DATA_DST / csv_name, TAB_DIR / f"table{tab_no}.tex")

    # Also copy sprint CSVs for analysis
    for p in RAW_SPRINT.glob("*.csv"):
        # Keep sprint pack files in a dedicated subfolder to prevent collisions.
        _copy_no_overwrite(p, SPRINT_DST / p.name)

    print("Organized evidence into:")
    print(f"- figures: {FIG_DIR}")
    print(f"- tables : {TAB_DIR}")
    print(f"- scripts: {SCRIPTS_DST}")
    print(f"- reports: {REPORTS_DST}")
    print(f"- data   : {DATA_DST}")


if __name__ == "__main__":
    main()

