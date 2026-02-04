#!/usr/bin/env python3
"""
Plot figures for phason_stats outputs.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import sys

# Allow running as a script without installing a package.
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
    ap.add_argument("--summary", required=True, help="phason_stats_summary_*.csv")
    ap.add_argument("--tag", default="n1000_seed0")
    args = ap.parse_args()

    root = Path(__file__).resolve().parents[2]
    summary_path = Path(args.summary)
    if not summary_path.is_absolute():
        summary_path = (root / summary_path).resolve()

    df = pd.read_csv(summary_path)
    out_dir = root / "figures" / "pdb_bench"
    ensure_dir(out_dir)

    for alpha_name in sorted(df["alphabet"].unique().tolist()):
        sub = df[df["alphabet"] == alpha_name]
        if len(sub) == 0:
            continue

        # Violin/box style using matplotlib's boxplot (robust, dependency-free)
        fig, ax = plt.subplots(figsize=(7.8, 4.3))
        data = [
            sub["ph_rms_real"].to_numpy(),
            sub["ph_rms_random_mean"].to_numpy(),
            sub["ph_rms_pert_mean"].to_numpy(),
        ]
        ax.boxplot(data, tick_labels=["real", "random_mean", "perturbed_mean"], showfliers=False)
        ax.set_ylabel("ph_rms (perp-space spread proxy)")
        ax.set_title(f"Phason proxy distribution ({alpha_name}, n={len(sub)})")
        fig.tight_layout()
        out1 = out_dir / f"phason_phrms_box_{alpha_name}_{args.tag}.png"
        fig.savefig(out1, dpi=200)
        plt.close(fig)

        # ECDF plot
        fig, ax = plt.subplots(figsize=(7.8, 4.3))
        for col, label in [
            ("ph_rms_real", "real"),
            ("ph_rms_random_mean", "random_mean"),
            ("ph_rms_pert_mean", "perturbed_mean"),
        ]:
            x, y = ecdf(sub[col].to_numpy())
            ax.plot(x, y, label=label)
        ax.set_xlabel("ph_rms")
        ax.set_ylabel("ECDF")
        ax.set_title(f"Phason proxy ECDF ({alpha_name})")
        ax.legend(frameon=False)
        fig.tight_layout()
        out2 = out_dir / f"phason_phrms_ecdf_{alpha_name}_{args.tag}.png"
        fig.savefig(out2, dpi=200)
        plt.close(fig)

        # N vs ph_rms scatter
        fig, ax = plt.subplots(figsize=(7.8, 4.3))
        ax.scatter(sub["N"], sub["ph_rms_real"], s=10, alpha=0.6, label="real")
        ax.scatter(sub["N"], sub["ph_rms_random_mean"], s=10, alpha=0.4, label="random_mean")
        ax.set_xlabel("N (chain length)")
        ax.set_ylabel("ph_rms")
        ax.set_title(f"Length dependence ({alpha_name})")
        ax.legend(frameon=False)
        fig.tight_layout()
        out3 = out_dir / f"phason_phrms_vs_N_{alpha_name}_{args.tag}.png"
        fig.savefig(out3, dpi=200)
        plt.close(fig)

        print(f"Wrote: {out1}")
        print(f"Wrote: {out2}")
        print(f"Wrote: {out3}")


if __name__ == "__main__":
    main()

