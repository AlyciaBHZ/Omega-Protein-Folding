#!/usr/bin/env python3
from __future__ import annotations

"""
Make Figure S5: Physical consistency vs Auric geometry certificate.
"""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def _scatter(ax: plt.Axes, df: pd.DataFrame, *, xcol: str, ycol: str, ylab: str, title: str) -> None:
    native = df[df["label"] == "native"]
    decoy = df[df["label"] != "native"]

    ax.scatter(decoy[xcol], decoy[ycol], s=42, alpha=0.85, c="#1f77b4", edgecolors="none", label="decoys")
    if len(native):
        ax.scatter(
            native[xcol],
            native[ycol],
            s=140,
            marker="*",
            c="#d62728",
            edgecolors="k",
            linewidths=0.6,
            label="native",
            zorder=5,
        )

    ax.set_title(title, fontsize=10)
    ax.set_xlabel("Auric geometry score: $H(\\mathrm{type})$ (rhoA, m=8)")
    ax.set_ylabel(ylab)
    ax.grid(True, alpha=0.25)

    x = df[xcol].to_numpy(dtype=np.float64)
    y = df[ycol].to_numpy(dtype=np.float64)
    mask = np.isfinite(x) & np.isfinite(y)
    r = float(np.corrcoef(x[mask], y[mask])[0, 1]) if int(np.sum(mask)) >= 3 else float("nan")
    ax.text(0.02, 0.98, f"Pearson r={r:.2f}", transform=ax.transAxes, ha="left", va="top", fontsize=9)
    ax.legend(loc="best", frameon=False, fontsize=9)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-csv", default="docs/runs/physics_audit/physics_vs_auric_toy_decoy.csv")
    ap.add_argument("--out-png", default="paper/assets/figures/figS_physics_vs_auric.png")
    args = ap.parse_args()

    root = repo_root()
    in_csv = Path(str(args.in_csv))
    if not in_csv.is_absolute():
        in_csv = (root / in_csv).resolve()
    out_png = Path(str(args.out_png))
    if not out_png.is_absolute():
        out_png = (root / out_png).resolve()
    ensure_dir(out_png.parent)

    df = pd.read_csv(in_csv)
    if len(df) == 0:
        raise ValueError(f"Empty CSV: {in_csv}")

    # Orthogonal-defense view (real decoy audit): x=auric entropy, y=physics_score + boxplot panel.
    if "physics_score" in df.columns and "auric_type_entropy_rhoA" in df.columns:
        df["auric_type_entropy_rhoA"] = pd.to_numeric(df["auric_type_entropy_rhoA"], errors="coerce")
        df["physics_score"] = pd.to_numeric(df["physics_score"], errors="coerce")

        native = df[df["label"] == "native"]
        decoy = df[df["label"] != "native"]

        fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.2), constrained_layout=True)

        ax = axes[0]
        ax.scatter(
            decoy["auric_type_entropy_rhoA"],
            decoy["physics_score"],
            s=16,
            alpha=0.35,
            c="#7f7f7f",
            edgecolors="none",
            label="decoys",
        )
        ax.scatter(
            native["auric_type_entropy_rhoA"],
            native["physics_score"],
            s=180,
            marker="*",
            c="#d62728",
            edgecolors="k",
            linewidths=0.6,
            label="native",
            zorder=5,
        )
        ax.set_title("(S5A) Dual-filter funnel (orthogonal signals)", fontsize=10)
        ax.set_xlabel("Geometric certificate: $H(\\mathrm{type})$ (rhoA, m=8)")
        ax.set_ylabel("Physics score: clashes/N + 1000×bond MSE")
        ax.grid(True, alpha=0.25)

        x = df["auric_type_entropy_rhoA"].to_numpy(dtype=np.float64)
        y = df["physics_score"].to_numpy(dtype=np.float64)
        mask = np.isfinite(x) & np.isfinite(y)
        r = float(np.corrcoef(x[mask], y[mask])[0, 1]) if int(np.sum(mask)) >= 3 else float("nan")
        ax.text(0.02, 0.98, f"Pearson r={r:.2f}", transform=ax.transAxes, ha="left", va="top", fontsize=9)
        ax.legend(loc="best", frameon=False, fontsize=9)

        bx = axes[1]
        bx.set_title("(S5B) Native vs decoy distributions", fontsize=10)
        bx.grid(True, axis="y", alpha=0.25)
        bx.boxplot(
            [
                native["auric_type_entropy_rhoA"].to_numpy(dtype=np.float64),
                decoy["auric_type_entropy_rhoA"].to_numpy(dtype=np.float64),
                native["physics_score"].to_numpy(dtype=np.float64),
                decoy["physics_score"].to_numpy(dtype=np.float64),
            ],
            labels=["Auric (N)", "Auric (D)", "Phys (N)", "Phys (D)"],
            showfliers=False,
        )
        bx.set_ylabel("Value (boxplot; fliers hidden)")

        fig.suptitle("Orthogonality of geometric and physical signals (4state_reduced)", fontsize=12)
        fig.savefig(out_png, dpi=240)
        print(f"[ok] wrote: {out_png}")
        return

    df["bond_mse_1e6"] = df["bond_mse"].astype(float) * 1e6

    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.2), constrained_layout=True)
    _scatter(
        axes[0],
        df,
        xcol="auric_type_entropy_rhoA",
        ycol="bond_mse_1e6",
        ylab="Bond geometry: MSE($d_{CA}-3.8\\AA$) × 1e6",
        title="(S5A) Bond geometry vs Auric",
    )
    _scatter(
        axes[1],
        df,
        xcol="auric_type_entropy_rhoA",
        ycol="clash_per_res",
        ylab="Steric proxy: CA clashes per residue (d<3.5Å, |i-j|>1)",
        title="(S5B) Steric clash proxy vs Auric",
    )
    fig.suptitle("Physical consistency audit (toy native vs decoys)", fontsize=12)
    fig.savefig(out_png, dpi=240)
    print(f"[ok] wrote: {out_png}")


if __name__ == "__main__":
    main()

