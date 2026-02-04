#!/usr/bin/env python3
from __future__ import annotations

"""
Record TM-align / US-align results into the committable CSV:
  docs/runs/quark_itasser_homology_ablation/tmalign_summary.csv

This supports the plan's step where you paste/compute structural agreement:
- QUARK-A model1 vs I-TASSER-A model1
- QUARK-A vs QUARK-B (homology-fragment dependence)
- I-TASSER-A vs I-TASSER-B (template dependence)

Usage examples:

  python scripts/benchmarks/record_tmalign_results.py \\
    --target-id 1r69 --dataset 4state_reduced \\
    --comparison \"QUARK_A_model1 vs ITASSER_A_model1\" \\
    --tm-score 0.73 --rmsd 2.1 --aligned-len 58 --coverage-q 0.92 --coverage-t 0.89

Or parse a pasted TM-align output:

  python scripts/benchmarks/record_tmalign_results.py \\
    --target-id 1r69 --dataset 4state_reduced \\
    --comparison \"QUARK_A_model1 vs ITASSER_A_model1\" \\
    --from-output-file path/to/tmalign_output.txt
"""

import argparse
import csv
import re
from pathlib import Path
from typing import Dict


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _parse_tmalign_text(txt: str) -> Dict[str, float]:
    """
    Best-effort TM-align/US-align output parsing.
    Returns dict with keys: tm_score, rmsd, aligned_len.
    """
    out: Dict[str, float] = {}

    # RMSD, aligned length: "RMSD=  2.31,   ... aligned length=  58"
    m = re.search(r"RMSD\s*=\s*([0-9.]+).*?aligned\s+length\s*=\s*(\d+)", txt, flags=re.I | re.S)
    if m:
        out["rmsd"] = float(m.group(1))
        out["aligned_len"] = float(m.group(2))

    # TM-score lines vary; capture first TM-score occurrence
    m = re.search(r"TM[\-\s]*score\s*=\s*([0-9.]+)", txt, flags=re.I)
    if m:
        out["tm_score"] = float(m.group(1))

    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--target-id", required=True)
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--comparison", required=True, help="Human-readable comparison label.")
    ap.add_argument("--tm-score", type=float, default=float("nan"))
    ap.add_argument("--rmsd", type=float, default=float("nan"))
    ap.add_argument("--aligned-len", type=float, default=float("nan"))
    ap.add_argument("--coverage-q", type=float, default=float("nan"))
    ap.add_argument("--coverage-t", type=float, default=float("nan"))
    ap.add_argument("--notes", default="")
    ap.add_argument("--from-output-file", default="", help="Parse TM-score/RMSD/aligned_len from a pasted TM-align output file.")
    args = ap.parse_args()

    if str(args.from_output_file).strip():
        p = Path(str(args.from_output_file))
        if not p.is_absolute():
            p = (Path.cwd() / p).resolve()
        txt = p.read_text(encoding="utf-8", errors="ignore")
        d = _parse_tmalign_text(txt)
        if "tm_score" in d:
            args.tm_score = float(d["tm_score"])
        if "rmsd" in d:
            args.rmsd = float(d["rmsd"])
        if "aligned_len" in d:
            args.aligned_len = float(d["aligned_len"])

    root = repo_root()
    csv_path = root / "docs" / "runs" / "quark_itasser_homology_ablation" / "tmalign_summary.csv"
    if not csv_path.exists():
        raise SystemExit(f"Missing CSV template: {csv_path}")

    row = {
        "target_id": str(args.target_id),
        "dataset": str(args.dataset),
        "comparison": str(args.comparison),
        "tm_score": f"{float(args.tm_score):.4f}" if float(args.tm_score) == float(args.tm_score) else "",
        "rmsd": f"{float(args.rmsd):.4f}" if float(args.rmsd) == float(args.rmsd) else "",
        "aligned_len": str(int(float(args.aligned_len))) if float(args.aligned_len) == float(args.aligned_len) else "",
        "coverage_q": f"{float(args.coverage_q):.4f}" if float(args.coverage_q) == float(args.coverage_q) else "",
        "coverage_t": f"{float(args.coverage_t):.4f}" if float(args.coverage_t) == float(args.coverage_t) else "",
        "notes": str(args.notes),
    }

    # append
    with csv_path.open("a", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=[
                "target_id",
                "dataset",
                "comparison",
                "tm_score",
                "rmsd",
                "aligned_len",
                "coverage_q",
                "coverage_t",
                "notes",
            ],
        )
        w.writerow(row)

    print(f"[ok] appended: {csv_path.relative_to(root).as_posix()}")


if __name__ == "__main__":
    main()

