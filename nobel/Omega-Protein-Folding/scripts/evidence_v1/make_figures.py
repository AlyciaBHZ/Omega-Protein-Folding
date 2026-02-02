#!/usr/bin/env python3
"""
Reproduce key figures for the Omega v1 evidence pack.

This script does NOT rerun the Omega folding engine. It regenerates the figures
from the saved CSV outputs included in ../data.

Usage:
  python scripts/make_figures.py --outdir figures_repro
"""
import argparse
import os
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt

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

def fig_sprint_tmscore(df, outpath):
    """
    Plot Sprint-to-0.9 summary.

    Supports two schemas:
      (1) per-run: columns include protein, mode in {A_baseline,B_challenger}, tm
      (2) aggregated: columns include protein, group in {A,B}, tm_best_mean
    """
    fig, ax = plt.subplots(figsize=(8, 4))

    if {"protein", "mode", "tm"}.issubset(df.columns):
        proteins = df["protein"].unique().tolist()
        for i, prot in enumerate(proteins):
            dprot = df[df["protein"] == prot]
            base = dprot[dprot["mode"] == "A_baseline"]["tm"].values
            chal = dprot[dprot["mode"] == "B_challenger"]["tm"].values
            ax.bar(i - 0.15, base.mean() if len(base) > 0 else float("nan"), width=0.25, label="A_baseline" if i == 0 else None)
            ax.bar(i + 0.15, chal.mean() if len(chal) > 0 else float("nan"), width=0.25, label="B_challenger (mean)" if i == 0 else None)
            ax.scatter([i - 0.15] * len(base), base, marker="o")
            ax.scatter([i + 0.15] * len(chal), chal, marker="x")
        ax.set_xticks(range(len(proteins)))
        ax.set_xticklabels(proteins, rotation=0)
        ax.set_ylabel("TM-score")
        ax.set_title("Sprint-to-0.9: TM-score A vs B")
        ax.legend(frameon=False)
        fig.tight_layout()
        fig.savefig(outpath, dpi=200)
        plt.close(fig)
        return

    if {"protein", "group", "tm_best_mean"}.issubset(df.columns):
        proteins = df["protein"].unique().tolist()
        for i, prot in enumerate(proteins):
            dprot = df[df["protein"] == prot]
            base = dprot[dprot["group"] == "A"]["tm_best_mean"].values
            chal = dprot[dprot["group"] == "B"]["tm_best_mean"].values
            ax.bar(i - 0.15, base[0] if len(base) > 0 else float("nan"), width=0.25, label="A_baseline" if i == 0 else None)
            ax.bar(i + 0.15, chal[0] if len(chal) > 0 else float("nan"), width=0.25, label="B_challenger" if i == 0 else None)
        ax.set_xticks(range(len(proteins)))
        ax.set_xticklabels(proteins, rotation=0)
        ax.set_ylabel("TM-score (best mean)")
        ax.set_title("Sprint-to-0.9 (summary): TM-score A vs B")
        ax.legend(frameon=False)
        fig.tight_layout()
        fig.savefig(outpath, dpi=200)
        plt.close(fig)
        return

    raise ValueError(f"Unsupported sprint summary schema. Columns: {list(df.columns)}")

def fig_sprint_phason_scatter(df, outpath):
    # Prefer per-run schema; otherwise fall back to aggregated ph_ratio.
    fig, ax = plt.subplots(figsize=(8, 4))

    if {"protein", "mode", "tm", "ph_max", "theta_star"}.issubset(df.columns):
        proteins = df["protein"].unique().tolist()
        markers = {"A_baseline": "o", "B_challenger": "x"}
        for prot in proteins:
            dprot = df[df["protein"] == prot]
            for mode in ["A_baseline", "B_challenger"]:
                dm = dprot[dprot["mode"] == mode]
                ax.scatter(dm["tm"], dm["ph_max"], marker=markers[mode], label=f"{prot} {mode}")
            theta = dprot["theta_star"].dropna().iloc[0]
            ax.axhline(theta, linewidth=1.0)
        ax.set_xlabel("TM-score")
        ax.set_ylabel("ph_max")
        ax.set_title("Sprint-to-0.9: ph_max vs TM with θ★ lines")
        ax.legend(frameon=False, fontsize=7, ncols=2)
        fig.tight_layout()
        fig.savefig(outpath, dpi=200)
        plt.close(fig)
        return

    if {"protein", "group", "tm_best_mean", "ph_ratio_best_mean"}.issubset(df.columns):
        proteins = df["protein"].unique().tolist()
        for prot in proteins:
            dprot = df[df["protein"] == prot]
            ax.scatter(dprot["tm_best_mean"], dprot["ph_ratio_best_mean"], label=prot)
        ax.set_xlabel("TM-score (best mean)")
        ax.set_ylabel("ph_ratio_best_mean")
        ax.set_title("Sprint-to-0.9 (summary): ph_ratio vs TM")
        ax.legend(frameon=False, fontsize=7)
        fig.tight_layout()
        fig.savefig(outpath, dpi=200)
        plt.close(fig)
        return

    raise ValueError(f"Unsupported sprint phason schema. Columns: {list(df.columns)}")

def fig_triple232_bar(df, outpath):
    # Bar chart comparing codecs per protein
    fig, ax = plt.subplots(figsize=(8,4))
    proteins = df['protein'].tolist()
    x = range(len(proteins))
    ax.bar([i-0.25 for i in x], df['axis12_tm'], width=0.25, label='axis12')
    ax.bar([i for i in x], df['pair72_tm'], width=0.25, label='pair72')
    ax.bar([i+0.25 for i in x], df['triple232_tm'], width=0.25, label='triple232')
    ax.set_xticks(list(x))
    ax.set_xticklabels(proteins)
    ax.set_ylabel("TM-score (oracle)")
    ax.set_title("Oracle expressivity: axis12 vs pair72 vs triple232")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(outpath, dpi=200)
    plt.close(fig)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="figures_repro")
    args = ap.parse_args()
    outdir = args.outdir
    ensure_dir(outdir)

    sprint = pd.read_csv(resolve_csv("sprint_to_0p9_AB_summary.csv"))
    fig_sprint_tmscore(sprint, os.path.join(outdir, "repro_sprint_tmscore.png"))
    fig_sprint_phason_scatter(sprint, os.path.join(outdir, "repro_sprint_phason_theta_scatter.png"))

    try:
        triple = pd.read_csv(resolve_csv("triple232_oracle_codec_tm_table_with_local100.csv"))
    except FileNotFoundError:
        triple = None
        print("WARN: triple232 oracle CSV not found; skipping repro_triple232_oracle_bar.png")
    if triple is not None:
        # Normalize columns
        tdf = triple.rename(columns={
            "axis12_tm": "axis12_tm",
            "pair72_tm": "pair72_tm",
            "triple232_oracle_tm": "triple232_tm",
        })
        if "triple232_tm" not in tdf.columns and "triple232_oracle_tm" in triple.columns:
            tdf["triple232_tm"] = triple["triple232_oracle_tm"]
        fig_triple232_bar(tdf[["protein", "axis12_tm", "pair72_tm", "triple232_tm"]], os.path.join(outdir, "repro_triple232_oracle_bar.png"))

    print(f"Saved figures to {outdir}")

if __name__ == "__main__":
    main()
