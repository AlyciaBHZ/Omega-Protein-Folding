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

    def add_if_present(cols: list[tuple[str, str]], col: str, label: str) -> None:
        if col in df.columns:
            cols.append((col, label))

    for alpha_name in sorted(df["alphabet"].unique().tolist()):
        suba = df[df["alphabet"] == alpha_name]
        for ro in sorted(suba["readout"].unique().tolist()):
            sub = suba[suba["readout"] == ro]
            for m in sorted(sub["m"].unique().tolist()):
                subm = sub[sub["m"] == m]
                if len(subm) == 0:
                    continue

                # Boxplot: real vs null means (random/pert/shuffle/blockshuffle*)
                for metric, ylab, stem in [
                    ("type_entropy", "H(type) [bits]", "entropy"),
                    ("smb_rate_hat", "smb_rate_hat", "smb_rate_hat"),
                ]:
                    cols: list[tuple[str, str]] = []
                    add_if_present(cols, f"{metric}_real", "real")
                    add_if_present(cols, f"{metric}_random_mean", "random_mean")
                    add_if_present(cols, f"{metric}_pert_mean", "pert_mean")
                    add_if_present(cols, f"{metric}_shuffle_mean", "shuffle_mean")
                    # blockshuffle_k{K}_mean columns (if present)
                    for c in sorted([c for c in df.columns if c.startswith(f"{metric}_blockshuffle_k") and c.endswith("_mean")]):
                        cols.append((c, c.replace(f"{metric}_", "")))

                    data = [subm[c].to_numpy(dtype=np.float64) for c, _ in cols]
                    labels = [lab for _, lab in cols]
                    fig, ax = plt.subplots(figsize=(9.6, 4.3))
                    ax.boxplot(data, tick_labels=labels, showfliers=False)
                    ax.set_ylabel(ylab)
                    ax.set_title(f"Auric {metric} ({alpha_name}, rho{ro}, m={m}, n={len(subm)})")
                    fig.tight_layout()
                    out1 = out_dir / f"auric_{stem}_box_{alpha_name}_rho{ro}_m{m}_{args.tag}.png"
                    fig.savefig(out1, dpi=200)
                    plt.close(fig)

                    # ECDF of Cliff's deltas vs nulls (if columns exist)
                    delta_cols: list[tuple[str, str]] = []
                    add_if_present(delta_cols, f"delta_{metric}_real_vs_random", "δ(real,random)")
                    add_if_present(delta_cols, f"delta_{metric}_real_vs_perturbed", "δ(real,perturbed)")
                    add_if_present(delta_cols, f"delta_{metric}_real_vs_shuffle", "δ(real,shuffle)")
                    for c in sorted([c for c in df.columns if c.startswith(f"delta_{metric}_real_vs_blockshuffle_k")]):
                        delta_cols.append((c, c.replace(f"delta_{metric}_real_vs_", "δ(real,") + ")"))
                    fig, ax = plt.subplots(figsize=(9.6, 4.3))
                    for c, lab in delta_cols:
                        x, y = ecdf(subm[c].to_numpy(dtype=np.float64))
                        if len(x) == 0:
                            continue
                        ax.plot(x, y, label=lab)
                    ax.set_xlabel("Cliff's δ")
                    ax.set_ylabel("ECDF")
                    ax.set_title(f"Auric effect size ECDF ({metric}, {alpha_name}, rho{ro}, m={m})")
                    ax.legend(frameon=False, ncols=2)
                    fig.tight_layout()
                    out2 = out_dir / f"auric_{stem}_delta_ecdf_{alpha_name}_rho{ro}_m{m}_{args.tag}.png"
                    fig.savefig(out2, dpi=200)
                    plt.close(fig)

                    print(f"Wrote: {out1}")
                    print(f"Wrote: {out2}")

    # Multi-scale curves (median metric vs m) for each (alphabet, readout)
    for alpha_name in sorted(df["alphabet"].unique().tolist()):
        suba = df[df["alphabet"] == alpha_name]
        for ro in sorted(suba["readout"].unique().tolist()):
            sub = suba[suba["readout"] == ro]
            ms = sorted(sub["m"].unique().tolist())
            if not ms:
                continue
            for metric, ylab, stem in [
                ("type_entropy", "median H(type) [bits]", "entropy"),
                ("smb_rate_hat", "median smb_rate_hat", "smb_rate_hat"),
            ]:
                series: list[tuple[str, str]] = []
                add_if_present(series, f"{metric}_real", "real")
                add_if_present(series, f"{metric}_random_mean", "random_mean")
                add_if_present(series, f"{metric}_pert_mean", "pert_mean")
                add_if_present(series, f"{metric}_shuffle_mean", "shuffle_mean")
                for c in sorted([c for c in df.columns if c.startswith(f"{metric}_blockshuffle_k") and c.endswith("_mean")]):
                    series.append((c, c.replace(f"{metric}_", "")))

                fig, ax = plt.subplots(figsize=(9.6, 4.3))
                for c, lab in series:
                    meds = []
                    for m in ms:
                        subm = sub[sub["m"] == m]
                        meds.append(float(np.nanmedian(subm[c].to_numpy(dtype=np.float64))))
                    ax.plot(ms, meds, marker="o", label=lab)
                ax.set_xlabel("m (Fold_m window length)")
                ax.set_ylabel(ylab)
                ax.set_title(f"Auric multi-scale ({metric}, {alpha_name}, rho{ro})")
                ax.legend(frameon=False, ncols=2)
                fig.tight_layout()
                out3 = out_dir / f"auric_{stem}_multiscale_{alpha_name}_rho{ro}_{args.tag}.png"
                fig.savefig(out3, dpi=200)
                plt.close(fig)
                print(f"Wrote: {out3}")


if __name__ == "__main__":
    main()

