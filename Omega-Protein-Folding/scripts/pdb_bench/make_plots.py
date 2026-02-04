#!/usr/bin/env python3
"""
Plot summary figures for the PDB oracle codec benchmark CSV.

Example:
  python scripts/pdb_bench/make_plots.py --csv data/processed/pdb_bench/oracle_codec_bench_n50_seed0_dirscale.csv --tag n50_seed0_dirscale
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True, help="oracle_codec_bench_*.csv")
    ap.add_argument("--tag", default="v1")
    args = ap.parse_args()

    root = Path(__file__).resolve().parents[2]
    csv_path = Path(args.csv)
    if not csv_path.is_absolute():
        csv_path = (root / csv_path).resolve()

    df = pd.read_csv(csv_path)
    out_dir = root / "figures" / "pdb_bench"
    ensure_dir(out_dir)

    # Histogram of TM
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    bins = [i / 50.0 for i in range(0, 51)]  # 0..1 step 0.02
    ax.hist(df["tm_axis12"], bins=bins, alpha=0.6, label="axis12")
    ax.hist(df["tm_pair72"], bins=bins, alpha=0.6, label="pair72")
    ax.hist(df["tm_triple232"], bins=bins, alpha=0.6, label="triple232")
    ax.set_xlabel("TM-score (oracle, Kabsch)")
    ax.set_ylabel("count")
    ax.set_title(f"Oracle codec TM distribution (n={len(df)})")
    ax.legend(frameon=False)
    fig.tight_layout()
    out1 = out_dir / f"oracle_codec_tm_hist_{args.tag}.png"
    fig.savefig(out1, dpi=200)
    plt.close(fig)

    # ECDF
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    for col, label in [
        ("tm_axis12", "axis12"),
        ("tm_pair72", "pair72"),
        ("tm_triple232", "triple232"),
    ]:
        x = df[col].dropna().sort_values().to_numpy()
        y = (np.arange(len(x)) + 1) / (len(x) if len(x) else 1)
        ax.plot(x, y, label=label)
    ax.set_xlabel("TM-score")
    ax.set_ylabel("ECDF")
    ax.set_title("Oracle codec TM ECDF")
    ax.legend(frameon=False)
    fig.tight_layout()
    out2 = out_dir / f"oracle_codec_tm_ecdf_{args.tag}.png"
    fig.savefig(out2, dpi=200)
    plt.close(fig)

    print(f"Wrote: {out1}")
    print(f"Wrote: {out2}")


if __name__ == "__main__":
    main()

