#!/usr/bin/env python3
"""
Plot figures for auric_stats outputs (Fold_m certificate metrics).
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
    ap.add_argument("--summary", required=True, help="auric_stats_summary_*.csv (from phason_stats.py --auric)")
    ap.add_argument("--tag", default="n1000_seed0")
    args = ap.parse_args()

    root = Path(__file__).resolve().parents[2]
    summary_path = Path(args.summary)
    if not summary_path.is_absolute():
        summary_path = (root / summary_path).resolve()

    df = pd.read_csv(summary_path)
    out_dir = root / "figures" / "pdb_bench"
    ensure_dir(out_dir)

    if len(df) == 0:
        raise ValueError("Empty auric summary CSV")

    for alpha_name in sorted(df["alphabet"].unique().tolist()):
        suba = df[df["alphabet"] == alpha_name]
        for ro in sorted(suba["readout"].unique().tolist()):
            sub = suba[suba["readout"] == ro]
            for m in sorted(sub["m"].unique().tolist()):
                subm = sub[sub["m"] == m]
                if len(subm) == 0:
                    continue

                # Boxplot: real vs control means (random/pert)
                fig, ax = plt.subplots(figsize=(7.8, 4.3))
                data = [
                    subm["type_entropy_real"].to_numpy(dtype=np.float64),
                    subm["type_entropy_random_mean"].to_numpy(dtype=np.float64),
                    subm["type_entropy_pert_mean"].to_numpy(dtype=np.float64),
                ]
                ax.boxplot(data, tick_labels=["real", "random_mean", "pert_mean"], showfliers=False)
                ax.set_ylabel("H(type) [bits]")
                ax.set_title(f"Auric type entropy ({alpha_name}, ρ{ro}, m={m}, n={len(subm)})")
                fig.tight_layout()
                out1 = out_dir / f"auric_entropy_box_{alpha_name}_rho{ro}_m{m}_{args.tag}.png"
                fig.savefig(out1, dpi=200)
                plt.close(fig)

                # ECDF of delta vs controls
                fig, ax = plt.subplots(figsize=(7.8, 4.3))
                for col, label in [
                    ("delta_type_entropy_real_vs_random", "δ(real,random)"),
                    ("delta_type_entropy_real_vs_perturbed", "δ(real,perturbed)"),
                ]:
                    x, y = ecdf(subm[col].to_numpy(dtype=np.float64))
                    ax.plot(x, y, label=label)
                ax.set_xlabel("Cliff's δ")
                ax.set_ylabel("ECDF")
                ax.set_title(f"Auric effect size ECDF ({alpha_name}, ρ{ro}, m={m})")
                ax.legend(frameon=False)
                fig.tight_layout()
                out2 = out_dir / f"auric_delta_ecdf_{alpha_name}_rho{ro}_m{m}_{args.tag}.png"
                fig.savefig(out2, dpi=200)
                plt.close(fig)

                print(f"Wrote: {out1}")
                print(f"Wrote: {out2}")

    # Multi-scale curve (median entropy vs m) for each (alphabet, readout)
    for alpha_name in sorted(df["alphabet"].unique().tolist()):
        suba = df[df["alphabet"] == alpha_name]
        for ro in sorted(suba["readout"].unique().tolist()):
            sub = suba[suba["readout"] == ro]
            ms = sorted(sub["m"].unique().tolist())
            if not ms:
                continue
            med_real = []
            med_rand = []
            med_pert = []
            for m in ms:
                subm = sub[sub["m"] == m]
                med_real.append(float(np.nanmedian(subm["type_entropy_real"].to_numpy(dtype=np.float64))))
                med_rand.append(float(np.nanmedian(subm["type_entropy_random_mean"].to_numpy(dtype=np.float64))))
                med_pert.append(float(np.nanmedian(subm["type_entropy_pert_mean"].to_numpy(dtype=np.float64))))
            fig, ax = plt.subplots(figsize=(7.8, 4.3))
            ax.plot(ms, med_real, marker="o", label="real")
            ax.plot(ms, med_rand, marker="o", label="random_mean")
            ax.plot(ms, med_pert, marker="o", label="pert_mean")
            ax.set_xlabel("m (Fold_m window length)")
            ax.set_ylabel("median H(type) [bits]")
            ax.set_title(f"Auric multi-scale entropy ({alpha_name}, ρ{ro})")
            ax.legend(frameon=False)
            fig.tight_layout()
            out3 = out_dir / f"auric_entropy_multiscale_{alpha_name}_rho{ro}_{args.tag}.png"
            fig.savefig(out3, dpi=200)
            plt.close(fig)
            print(f"Wrote: {out3}")


if __name__ == "__main__":
    main()

