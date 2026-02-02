#!/usr/bin/env python3
"""
Generate manuscript tables from the saved CSV outputs.

Usage:
  python scripts/make_tables.py --outdir tables_repro
"""
import argparse
import os
from pathlib import Path

import pandas as pd

def ensure_dir(d):
    os.makedirs(d, exist_ok=True)

def resolve_csv(rel_name: str) -> str:
    """
    Prefer repo-integrated location under data/processed/, but fall back to data/
    for older evidence-pack layouts.
    """
    candidates = [
        Path("data") / "processed" / rel_name,
        Path("data") / rel_name,
    ]
    for p in candidates:
        if p.exists():
            return str(p)
    raise FileNotFoundError(f"Could not find CSV: {rel_name} (tried: {candidates})")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="tables_repro")
    args = ap.parse_args()
    outdir = args.outdir
    ensure_dir(outdir)

    # Table: sprint A/B
    sprint = pd.read_csv(resolve_csv("sprint_to_0p9_AB_summary.csv"))
    sprint.to_csv(os.path.join(outdir, "table_sprint_ab_full.csv"), index=False)

    # Table: theta_star feasibility
    theta = pd.read_csv(resolve_csv("theta_star_AB_success_1AKE_1TIM.csv"))
    theta.to_csv(os.path.join(outdir, "table_theta_star_feasibility.csv"), index=False)

    # Table: triple232 oracle
    try:
        triple = pd.read_csv(resolve_csv("triple232_oracle_codec_tm_table_with_local100.csv"))
    except FileNotFoundError:
        triple = None
        print("WARN: triple232 oracle CSV not found; skipping table_triple232_oracle.csv")
    if triple is not None:
        triple.to_csv(os.path.join(outdir, "table_triple232_oracle.csv"), index=False)

    print(f"Saved tables to {outdir}")

if __name__ == "__main__":
    main()
