#!/usr/bin/env python3
"""
Plot figures for beam_lift_bench outputs.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))

from bench_utils import ensure_dir  # noqa: E402


def ecdf(x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    x = np.asarray(x, dtype=np.float64)
    x = x[np.isfinite(x)]
    if len(x) == 0:
        return np.array([]), np.array([])
    x = np.sort(x)
    y = (np.arange(len(x)) + 1) / len(x)
    return x, y


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True, help="beam_lift_bench_summary_*.csv")
    ap.add_argument("--tag", default="v1")
    args = ap.parse_args()

    root = Path(__file__).resolve().parents[2]
    csv_path = Path(args.csv)
    if not csv_path.is_absolute():
        csv_path = (root / csv_path).resolve()

    df = pd.read_csv(csv_path)
    out_dir = root / "figures" / "pdb_bench"
    ensure_dir(out_dir)

    # Boxplot ph_rms
    fig, ax = plt.subplots(figsize=(7.8, 4.3))
    ax.boxplot(
        [df["ph_rms_real"].to_numpy(), df["ph_rms_random_mean"].to_numpy()],
        tick_labels=["real", "random_mean"],
        showfliers=False,
    )
    ax.set_ylabel("ph_rms")
    ax.set_title(f"Beam-lift phason proxy (n={len(df)})")
    fig.tight_layout()
    out1 = out_dir / f"beam_lift_phrms_box_{args.tag}.png"
    fig.savefig(out1, dpi=200)
    plt.close(fig)

    # ECDF
    fig, ax = plt.subplots(figsize=(7.8, 4.3))
    for col, label in [("ph_rms_real", "real"), ("ph_rms_random_mean", "random_mean")]:
        x, y = ecdf(df[col].to_numpy())
        ax.plot(x, y, label=label)
    ax.set_xlabel("ph_rms")
    ax.set_ylabel("ECDF")
    ax.set_title("Beam-lift ph_rms ECDF")
    ax.legend(frameon=False)
    fig.tight_layout()
    out2 = out_dir / f"beam_lift_phrms_ecdf_{args.tag}.png"
    fig.savefig(out2, dpi=200)
    plt.close(fig)

    # TM histogram
    fig, ax = plt.subplots(figsize=(7.8, 4.3))
    bins = [i / 50.0 for i in range(0, 51)]
    ax.hist(df["tm_beam"].to_numpy(), bins=bins, alpha=0.8)
    ax.set_xlabel("TM-score (beam-lift, Kabsch)")
    ax.set_ylabel("count")
    ax.set_title("Beam-lift TM distribution")
    fig.tight_layout()
    out3 = out_dir / f"beam_lift_tm_hist_{args.tag}.png"
    fig.savefig(out3, dpi=200)
    plt.close(fig)

    print(f"Wrote: {out1}")
    print(f"Wrote: {out2}")
    print(f"Wrote: {out3}")


if __name__ == "__main__":
    main()

