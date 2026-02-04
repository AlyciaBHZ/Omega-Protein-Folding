#!/usr/bin/env python3
from __future__ import annotations

"""
Parse downloaded QUARK / I-TASSER outputs into committable summaries.

This script is designed to run after you manually download server outputs into:
- data/raw/quark_cache/<target_id>/QAxxxx/
- data/raw/itasser_cache/<target_id>/Sxxxxx/

It writes small summaries under:
- docs/runs/quark_itasser_homology_ablation/

We intentionally keep tarballs and full server outputs untracked.
"""

import argparse
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import pandas as pd


ROOT_MARKER = "docs/runs/quark_itasser_homology_ablation"


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def read_text_maybe(path: Path, max_bytes: int = 2_000_000) -> str:
    if not path.exists() or not path.is_file():
        return ""
    b = path.read_bytes()
    if len(b) > max_bytes:
        b = b[:max_bytes]
    return b.decode("utf-8", errors="ignore")


def iter_files(root: Path) -> Iterable[Path]:
    if not root.exists():
        return []
    for p in root.rglob("*"):
        if p.is_file():
            yield p


@dataclass(frozen=True)
class ItasserSummary:
    job_id: str
    c_scores: Dict[int, float]
    est_tm: float
    est_rmsd: float
    top_templates: List[str]  # raw lines (best-effort)


def parse_itasser_text_blob(txt: str) -> ItasserSummary:
    """
    Best-effort parser that supports multiple I-TASSER result formats.
    We look for:
    - C-score lines
    - Estimated TM-score / RMSD
    - Top template table lines (if present)
    """
    # job id is not always in files; leave empty if absent
    job = ""
    m = re.search(r"\b(S\d{5,})\b", txt)
    if m:
        job = m.group(1)

    # C-score: allow "C-score" or "Cscore"
    c_scores: Dict[int, float] = {}
    for mm in re.finditer(r"Model\s*(\d+)\s*[:\-]?\s*C[\-\s]*score\s*[:=]\s*([\-0-9.]+)", txt, flags=re.I):
        i = int(mm.group(1))
        v = float(mm.group(2))
        c_scores[i] = v

    # Some pages list "C-score ="
    if not c_scores:
        for mm in re.finditer(r"\bC[\-\s]*score\s*[:=]\s*([\-0-9.]+)", txt, flags=re.I):
            # cannot assign model index; store as model1
            try:
                c_scores.setdefault(1, float(mm.group(1)))
            except Exception:
                pass

    est_tm = float("nan")
    est_rmsd = float("nan")
    mm = re.search(r"Estimated\s*TM[\-\s]*score\s*[:=]\s*([0-9.]+)", txt, flags=re.I)
    if mm:
        est_tm = float(mm.group(1))
    mm = re.search(r"Estimated\s*RMSD\s*[:=]\s*([0-9.]+)", txt, flags=re.I)
    if mm:
        est_rmsd = float(mm.group(1))

    # templates: capture a small block following "Top templates" or "threading templates"
    top_templates: List[str] = []
    block = re.search(
        r"(Top\s+\d+\s+(?:threading\s+)?templates.*?)(?:\n\s*\n|\Z)",
        txt,
        flags=re.I | re.S,
    )
    if block:
        lines = [ln.rstrip() for ln in block.group(1).splitlines() if ln.strip()]
        top_templates = lines[:30]

    return ItasserSummary(job_id=job, c_scores=c_scores, est_tm=est_tm, est_rmsd=est_rmsd, top_templates=top_templates)


@dataclass(frozen=True)
class QuarkSummary:
    job_id: str
    ss_line: str
    sa_line: str


def parse_quark_text_blob(txt: str) -> QuarkSummary:
    """
    QUARK pages commonly expose lines labeled:
    - Predicted Secondary Structure
    - Predicted Solvent Accessibility
    We store the first occurrence of each.
    """
    job = ""
    m = re.search(r"\b(QA\d{4,})\b", txt)
    if m:
        job = m.group(1)

    ss_line = ""
    sa_line = ""
    for ln in txt.splitlines():
        if (not ss_line) and re.search(r"Predicted\s+Secondary\s+Structure", ln, flags=re.I):
            ss_line = ln.strip()
            continue
        if (not sa_line) and re.search(r"Predicted\s+Solvent\s+Accessibility", ln, flags=re.I):
            sa_line = ln.strip()
            continue
    return QuarkSummary(job_id=job, ss_line=ss_line, sa_line=sa_line)


def find_itasser_text_files(job_dir: Path) -> List[Path]:
    # Typical names seen in I-TASSER bundles
    pats = [
        "*report*.txt",
        "*summary*.txt",
        "*README*",
        "*.html",
        "*.txt",
    ]
    out: List[Path] = []
    for pat in pats:
        out.extend(sorted(job_dir.glob(pat)))
    # de-dup preserving order
    seen = set()
    out2 = []
    for p in out:
        if p not in seen:
            seen.add(p)
            out2.append(p)
    return out2[:50]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", default=ROOT_MARKER, help="Docs run directory to write summaries into.")
    ap.add_argument("--quark-root", default="data/raw/quark_cache", help="Local untracked QUARK cache root.")
    ap.add_argument("--itasser-root", default="data/raw/itasser_cache", help="Local untracked I-TASSER cache root.")
    args = ap.parse_args()

    root = repo_root()
    run_dir = Path(str(args.run_dir))
    if not run_dir.is_absolute():
        run_dir = (root / run_dir).resolve()
    ensure_dir(run_dir)

    quark_root = Path(str(args.quark_root))
    if not quark_root.is_absolute():
        quark_root = (root / quark_root).resolve()

    it_root = Path(str(args.itasser_root))
    if not it_root.is_absolute():
        it_root = (root / it_root).resolve()

    rows: List[Dict[str, object]] = []

    # QUARK: treat each QAxxxx directory as a job
    if quark_root.exists():
        for job_dir in sorted([p for p in quark_root.rglob("*") if p.is_dir() and re.fullmatch(r"QA\d+", p.name)]):
            txt = ""
            # try to parse an html/text result page if present
            for p in sorted(job_dir.glob("*.html")) + sorted(job_dir.glob("*.txt")):
                txt = read_text_maybe(p)
                if txt:
                    break
            q = parse_quark_text_blob(txt) if txt else QuarkSummary(job_id=job_dir.name, ss_line="", sa_line="")
            rows.append(
                {
                    "method": "QUARK",
                    "job_id": q.job_id or job_dir.name,
                    "target_id": job_dir.parent.name,
                    "ss_line": q.ss_line,
                    "sa_line": q.sa_line,
                }
            )

    # I-TASSER: each Sxxxxx directory as a job
    if it_root.exists():
        for job_dir in sorted([p for p in it_root.rglob("*") if p.is_dir() and re.fullmatch(r"S\d+", p.name)]):
            blob = ""
            for p in find_itasser_text_files(job_dir):
                t = read_text_maybe(p)
                if t:
                    blob = t
                    break
            it = parse_itasser_text_blob(blob) if blob else ItasserSummary(job_id=job_dir.name, c_scores={}, est_tm=float("nan"), est_rmsd=float("nan"), top_templates=[])
            rows.append(
                {
                    "method": "ITASSER",
                    "job_id": it.job_id or job_dir.name,
                    "target_id": job_dir.parent.name,
                    "c_scores": ";".join(f"{k}:{v:.3f}" for k, v in sorted(it.c_scores.items())),
                    "est_tm": it.est_tm,
                    "est_rmsd": it.est_rmsd,
                    "top_templates": " | ".join(it.top_templates[:10]),
                }
            )

    out_csv = run_dir / "server_summaries.csv"
    pd.DataFrame(rows).to_csv(out_csv, index=False)
    print(f"[ok] wrote: {out_csv}")


if __name__ == "__main__":
    main()

