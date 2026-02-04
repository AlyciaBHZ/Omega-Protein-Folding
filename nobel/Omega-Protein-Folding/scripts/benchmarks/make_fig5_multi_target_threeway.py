#!/usr/bin/env python3
from __future__ import annotations

"""
Figure 5 (upgraded): multi-target three-way audit (Omega vs QUARK vs I-TASSER).

Targets: 1R69 and 2CRO.

Panels:
- A: Topological accuracy (TM-score vs native) using in-repo tm_kabsch entries
     in docs/runs/quark_itasser_homology_ablation/server_summaries.csv
- B: Physical stability audit (Rosetta constrained FastRelax, nstruct=20):
     grouped boxplots of drift distributions (Cα RMSD to each method's input trace)
     using committable per-struct CSVs in the run folder.
- C: Qualitative overlay (1R69) as an embedded image.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
import pandas as pd


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


RUN_DIR = Path("docs/runs/quark_itasser_homology_ablation")


def _load_tm(targets: list[str]) -> pd.DataFrame:
    root = repo_root()
    p = (root / RUN_DIR / "server_summaries.csv").resolve()
    df = pd.read_csv(p)
    df["target_id"] = df["target_id"].astype(str).str.lower()
    df["method"] = df["method"].astype(str).str.lower()
    df["source"] = df["source"].astype(str).str.lower()
    df = df[df["source"] == "tm_kabsch"].copy()
    df = df[df["target_id"].isin([t.lower() for t in targets])].copy()
    df["tm_score_to_native"] = pd.to_numeric(df["tm_score_to_native"], errors="coerce")
    return df


def _load_per_struct(target: str, method: str, nstruct: int = 20) -> np.ndarray:
    root = repo_root()
    p = (root / RUN_DIR / f"rosetta_relax_{method}_{target}_n{nstruct}_per_struct.csv").resolve()
    df = pd.read_csv(p)
    x = pd.to_numeric(df["rmsd_ca_to_input_trace"], errors="coerce").dropna().to_numpy()
    if x.size == 0:
        raise RuntimeError(f"Empty per-struct drift: {p}")
    return x


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", default="1r69,2cro")
    ap.add_argument("--nstruct", type=int, default=20)
    args = ap.parse_args()

    targets = [t.strip().lower() for t in str(args.targets).split(",") if t.strip()]
    nstruct = int(args.nstruct)
    methods = ["omega", "quark", "itasser"]
    colors = {"omega": "#F28E2B", "quark": "#4E79A7", "itasser": "#59A14F"}

    tm_df = _load_tm(targets)

    fig = plt.figure(figsize=(14.2, 4.8), dpi=180)
    gs = fig.add_gridspec(1, 3, width_ratios=[1.1, 1.6, 1.2], wspace=0.32)

    # Panel A: grouped TM bars (targets on x; methods in legend)
    axA = fig.add_subplot(gs[0, 0])
    x = np.arange(len(targets), dtype=float)
    width = 0.22
    for j, m in enumerate(methods):
        vals = []
        for t in targets:
            row = tm_df[(tm_df["target_id"] == t) & (tm_df["method"] == m)]
            vals.append(float(row["tm_score_to_native"].iloc[0]) if not row.empty else float("nan"))
        axA.bar(x + (j - 1) * width, vals, width=width, color=colors[m], edgecolor="black", linewidth=0.6, label=m)
        for i, v in enumerate(vals):
            if np.isfinite(v):
                axA.text(x[i] + (j - 1) * width, v + 0.02, f"{v:.2f}", ha="center", va="bottom", fontsize=9)
    axA.set_xticks(x, [t.upper() for t in targets])
    axA.set_ylim(0.0, 1.0)
    axA.set_ylabel("TM-score vs native (Kabsch single-pass)")
    axA.set_title("A  Topological accuracy")
    axA.legend(frameon=False, fontsize=9)

    # Panel B: grouped drift boxplots (targets x methods)
    axB = fig.add_subplot(gs[0, 1])
    positions = []
    data = []
    box_colors = []
    xticks = []
    xticklabels = []

    group_gap = 1.2
    method_gap = 0.28
    cur = 0.0
    for t in targets:
        group_center = cur
        for k, m in enumerate(methods):
            pos = cur + (k - 1) * method_gap
            positions.append(pos)
            data.append(_load_per_struct(t, m, nstruct=nstruct))
            box_colors.append(colors[m])
        xticks.append(group_center)
        xticklabels.append(t.upper())
        cur += group_gap

    bp = axB.boxplot(
        data,
        positions=positions,
        widths=0.22,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "black", "linewidth": 1.2},
        whiskerprops={"color": "black", "linewidth": 1.0},
        capprops={"color": "black", "linewidth": 1.0},
    )
    for patch, c in zip(bp["boxes"], box_colors):
        patch.set_facecolor(c)
        patch.set_alpha(0.55)
        patch.set_edgecolor("black")
        patch.set_linewidth(1.0)

    axB.set_xticks(xticks, xticklabels)
    axB.set_ylabel("Cα RMSD drift to input trace (Å)")
    axB.set_title(f"B  Physical stability audit (constrained FastRelax; n={nstruct})")
    axB.grid(True, axis="y", alpha=0.25)

    # Legend for methods (colors)
    handles = [plt.Line2D([0], [0], color=colors[m], lw=8, alpha=0.55) for m in methods]
    axB.legend(handles, methods, frameon=False, fontsize=9, loc="upper left")

    # Panel C: embed overlay image (1R69)
    axC = fig.add_subplot(gs[0, 2])
    overlay = (repo_root() / RUN_DIR / "figures" / "fig5_panelC_overlay_1r69.png").resolve()
    img = mpimg.imread(str(overlay))
    axC.imshow(img)
    axC.axis("off")
    axC.set_title("C  Qualitative overlay (1R69)")

    out_dir = (repo_root() / RUN_DIR / "figures").resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    out_png = out_dir / "fig5_multi_target.png"
    fig.savefig(out_png, bbox_inches="tight")
    plt.close(fig)
    print(f"[ok] wrote: {out_png}")


if __name__ == "__main__":
    main()

