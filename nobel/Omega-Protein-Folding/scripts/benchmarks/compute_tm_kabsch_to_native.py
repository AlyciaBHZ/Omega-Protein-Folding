#!/usr/bin/env python3
from __future__ import annotations

"""
Compute "Kabsch single-pass TM-score" to native for committable models.

IMPORTANT:
- This uses the in-repo implementation in `scripts/pdb_bench/bench_utils.py`:
  one Kabsch alignment + Zhang&Skolnick TM formula.
- This is not TM-align/US-align (no alignment search).

Native structures are downloaded from RCSB as mmCIF into:
  data/raw/pdb_cache/{PDBID}.cif   (gitignored)
"""

import argparse
import csv
from pathlib import Path
from typing import Dict, List

import numpy as np
import requests


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def ensure_native_cif(pdb_id: str) -> Path:
    root = repo_root()
    cache = root / "data" / "raw" / "pdb_cache"
    cache.mkdir(parents=True, exist_ok=True)
    p = cache / f"{pdb_id.upper()}.cif"
    if p.exists():
        return p
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.cif"
    r = requests.get(url, timeout=60)
    r.raise_for_status()
    p.write_bytes(r.content)
    return p


def load_ca(path: Path) -> np.ndarray:
    root = repo_root()
    import sys

    sys.path.insert(0, str((root / "scripts" / "pdb_bench").resolve()))
    from coords_io import load_ca_coords  # type: ignore

    _, ca = load_ca_coords(path)
    return np.asarray(ca, dtype=np.float64)


def tm_kabsch(P: np.ndarray, Q: np.ndarray) -> float:
    root = repo_root()
    import sys

    sys.path.insert(0, str((root / "scripts" / "pdb_bench").resolve()))
    from bench_utils import tm_score  # type: ignore

    return float(tm_score(P, Q))


def upsert_server_summaries(*, target_id: str, updates: Dict[str, float]) -> None:
    """
    Update docs/runs/quark_itasser_homology_ablation/server_summaries.csv
    for rows (target_id, method) with source=tm_kabsch (create if absent).
    """
    root = repo_root()
    p = root / "docs" / "runs" / "quark_itasser_homology_ablation" / "server_summaries.csv"
    rows: List[dict] = []
    if p.exists():
        with p.open("r", encoding="utf-8", newline="") as f:
            rows = list(csv.DictReader(f))

    fieldnames = ["target_id", "method", "tm_score_to_native", "source", "notes"]
    # drop old tm_kabsch rows for this target/method
    out: List[dict] = []
    for r in rows:
        if str(r.get("target_id", "")).lower() == target_id.lower() and str(r.get("method", "")).lower() in updates:
            if str(r.get("source", "")).strip().lower() == "tm_kabsch":
                continue
        out.append(r)

    for method, tm in updates.items():
        out.append(
            {
                "target_id": target_id,
                "method": method,
                "tm_score_to_native": f"{float(tm):.4f}",
                "source": "tm_kabsch",
                "notes": "Kabsch single-pass TM-score (in-repo bench_utils.tm_score) vs RCSB native mmCIF.",
            }
        )

    with p.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, lineterminator="\n")
        w.writeheader()
        for r in out:
            w.writerow({k: r.get(k, "") for k in fieldnames})


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--target-id", required=True, help="e.g. 1r69 or 2cro (PDB id)")
    ap.add_argument("--omega-pdb", required=True)
    ap.add_argument("--quark-pdb", required=True)
    ap.add_argument("--itasser-pdb", required=True)
    args = ap.parse_args()

    root = repo_root()
    tid = str(args.target_id).strip()
    native = ensure_native_cif(tid)

    def rp(s: str) -> Path:
        p = Path(s)
        return p if p.is_absolute() else (root / p).resolve()

    models = {
        "omega": rp(str(args.omega_pdb)),
        "quark": rp(str(args.quark_pdb)),
        "itasser": rp(str(args.itasser_pdb)),
    }

    ca_native = load_ca(native)
    updates: Dict[str, float] = {}
    for method, mp in models.items():
        ca = load_ca(mp)
        if ca.shape != ca_native.shape:
            raise SystemExit(f"Length mismatch {tid} {method}: model {ca.shape[0]} native {ca_native.shape[0]}")
        updates[method] = tm_kabsch(ca, ca_native)
        print(f"[ok] {tid} {method}: tm_kabsch={updates[method]:.4f}")

    upsert_server_summaries(target_id=tid, updates=updates)
    print("[ok] updated server_summaries.csv")


if __name__ == "__main__":
    main()

