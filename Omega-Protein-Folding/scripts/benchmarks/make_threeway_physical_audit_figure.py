#!/usr/bin/env python3
from __future__ import annotations

"""
Three-way physical audit figure generator.

Produces:
- Panel A (optional): TM-score bars if `server_summaries.csv` contains entries for the target.
- Panel B: Drift distribution boxplots using committable per-struct CSVs (nstruct=20).

Designed for reviewer-facing benchmark figures under:
  docs/runs/quark_itasser_homology_ablation/
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


RUN_DIR = Path("docs/runs/quark_itasser_homology_ablation")


def _load_tm(target_id: str) -> dict[str, float]:
    root = repo_root()
    p = (root / RUN_DIR / "server_summaries.csv").resolve()
    if not p.exists():
        return {}
    df = pd.read_csv(p)
    df = df[df["target_id"].astype(str).str.lower() == target_id.lower()].copy()
    out: dict[str, float] = {}
    for _, r in df.iterrows():
        m = str(r["method"]).strip().lower()
        try:
            if str(r.get("source", "")).strip().lower() != "tm_kabsch":
                continue
            out[m] = float(r["tm_score_to_native"])
        except Exception:
            continue
    return out


def _load_per_struct(target_id: str, nstruct: int) -> dict[str, pd.DataFrame]:
    root = repo_root()
    run_dir = (root / RUN_DIR).resolve()
    out: dict[str, pd.DataFrame] = {}
    for m in ["omega", "quark", "itasser"]:
        p = run_dir / f"rosetta_relax_{m}_{target_id}_n{int(nstruct)}_per_struct.csv"
        if not p.exists():
            continue
        df = pd.read_csv(p)
        df["total_score"] = pd.to_numeric(df["total_score"], errors="coerce")
        df["rmsd_ca_to_input_trace"] = pd.to_numeric(df["rmsd_ca_to_input_trace"], errors="coerce")
        out[m] = df
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--target-id", default="1r69")
    ap.add_argument("--nstruct", type=int, default=20)
    args = ap.parse_args()

    target_id = str(args.target_id)
    nstruct = int(args.nstruct)

    per = _load_per_struct(target_id, nstruct)
    if not all(m in per for m in ["omega", "quark", "itasser"]):
        missing = [m for m in ["omega", "quark", "itasser"] if m not in per]
        raise SystemExit(f"Missing per-struct CSV(s) for {target_id} n={nstruct}: {missing}")

    tm = _load_tm(target_id)
    have_tm = all(m in tm for m in ["omega", "quark", "itasser"])

    methods = ["omega", "quark", "itasser"]
    colors = {"omega": "#F28E2B", "quark": "#4E79A7", "itasser": "#59A14F"}
    cs = [colors[m] for m in methods]

    if have_tm:
        fig = plt.figure(figsize=(10, 4.2), dpi=180)
        gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.25], wspace=0.28)

        axA = fig.add_subplot(gs[0, 0])
        tm_vals = [tm[m] for m in methods]
        axA.bar(methods, tm_vals, color=cs, edgecolor="black", linewidth=0.6)
        axA.set_ylim(0.0, 1.0)
        axA.set_ylabel("TM-score vs native")
        axA.set_title(f"A  Topological accuracy ({target_id.upper()})")
        for i, v in enumerate(tm_vals):
            axA.text(i, v + 0.02, f"{v:.2f}", ha="center", va="bottom", fontsize=9)

        axB = fig.add_subplot(gs[0, 1])
    else:
        fig = plt.figure(figsize=(6.4, 4.2), dpi=180)
        axB = fig.add_subplot(111)

    x = list(range(len(methods)))
    data = [per[m]["rmsd_ca_to_input_trace"].dropna().values for m in methods]
    bp = axB.boxplot(
        data,
        positions=x,
        widths=0.55,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "black", "linewidth": 1.2},
        whiskerprops={"color": "black", "linewidth": 1.0},
        capprops={"color": "black", "linewidth": 1.0},
    )
    axB.set_xticks(x, methods)
    axB.set_ylabel("Cα RMSD drift to input trace (Å)")
    axB.set_title(f"{'B' if have_tm else 'Three-way'}  Physical audit under constrained FastRelax (n={nstruct})")
    axB.grid(True, axis="y", alpha=0.25)

    for patch, c in zip(bp["boxes"], cs):
        patch.set_facecolor(c)
        patch.set_alpha(0.55)
        patch.set_edgecolor("black")
        patch.set_linewidth(1.0)

    best_vals = [float(per[m]["rmsd_ca_to_input_trace"].min()) for m in methods]
    axB.scatter(x, best_vals, s=45, c=cs, edgecolors="black", linewidths=0.7, zorder=3)
    for i, v in enumerate(best_vals):
        axB.text(i, v + 0.02, f"best {v:.2f}Å", ha="center", va="bottom", fontsize=9)

    out_dir = (repo_root() / RUN_DIR / "figures").resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    out_png = out_dir / f"threeway_{target_id.lower()}_panelAB.png"
    fig.savefig(out_png, bbox_inches="tight")
    plt.close(fig)
    print(f"[ok] wrote: {out_png}")


if __name__ == "__main__":
    main()

