#!/usr/bin/env python3
from __future__ import annotations

"""
Figure 5 (draft): Omega vs QUARK vs I-TASSER on 1R69.

Panels:
- A: Topological accuracy (TM-score vs native) from server_summaries.csv
- B: Physical stability audit under Rosetta constrained relax (nstruct=20):
     Drift distribution (Cα RMSD to the *input trace*) across relaxation trajectories.
     If per-struct CSVs are available, render boxplots; otherwise fall back to best-of-20 points.
"""

from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


RUN_DIR = Path("docs/runs/quark_itasser_homology_ablation")


@dataclass(frozen=True)
class Row:
    method: str
    tm: float
    best_energy: float
    best_rmsd: float


def _load_tm(target_id: str) -> dict[str, float]:
    root = repo_root()
    p = (root / RUN_DIR / "server_summaries.csv").resolve()
    df = pd.read_csv(p)
    df = df[df["target_id"].astype(str).str.lower() == target_id.lower()].copy()
    if df.empty:
        raise RuntimeError(f"Missing TM entries for target={target_id} in {p}")
    out: dict[str, float] = {}
    for _, r in df.iterrows():
        out[str(r["method"]).strip().lower()] = float(r["tm_score_to_native"])
    return out


def _load_rosetta_best(target_id: str) -> dict[str, tuple[float, float]]:
    """
    Return method -> (best_total_score, best_rmsd_ca) for nstruct=20.
    """
    root = repo_root()
    # Omega is in the 6-target summary.
    omega_p = (root / RUN_DIR / "rosetta_relax_summary.csv").resolve()
    omega = pd.read_csv(omega_p)
    omega = omega[(omega["method"] == "omega") & (omega["target_id"].astype(str).str.lower() == target_id.lower())].copy()
    if omega.empty:
        raise RuntimeError(f"Missing omega entry for {target_id} in {omega_p}")
    omega_row = omega.iloc[0]
    omega_best = (float(omega_row["best_total_score"]), float(omega_row["best_rmsd_ca"]))

    # QUARK + I-TASSER are in the server-model audit summary for this target.
    srv_p = (root / RUN_DIR / "rosetta_relax_server_models_1r69_n20.csv").resolve()
    srv = pd.read_csv(srv_p)
    out: dict[str, tuple[float, float]] = {
        "omega": omega_best,
    }
    for _, r in srv.iterrows():
        m = str(r["method"]).strip().lower()
        out[m] = (float(r["best_total_score"]), float(r["best_rmsd_ca"]))
    return out


def _load_per_struct(target_id: str) -> dict[str, pd.DataFrame]:
    root = repo_root()
    # Prefer committable copies under the run folder (so Figure 5 is reproducible from a clean clone).
    run_dir = (root / RUN_DIR).resolve()
    committable = {
        "omega": run_dir / f"rosetta_relax_omega_{target_id}_n20_per_struct.csv",
        "quark": run_dir / f"rosetta_relax_quark_{target_id}_n20_per_struct.csv",
        "itasser": run_dir / f"rosetta_relax_itasser_{target_id}_n20_per_struct.csv",
    }
    cache_base = root / "data" / "raw" / "rosetta_cache" / "quark_itasser_homology_ablation"
    cached = {
        "omega": cache_base / "omega" / target_id / "relax_1r69_per_struct.csv",
        "quark": cache_base / "quark" / target_id / "relax_quark_1r69_per_struct.csv",
        "itasser": cache_base / "itasser" / target_id / "relax_itasser_1r69_per_struct.csv",
    }
    paths = {k: (committable[k] if committable[k].exists() else cached[k]) for k in committable}
    out: dict[str, pd.DataFrame] = {}
    for m, p in paths.items():
        if p.exists():
            df = pd.read_csv(p)
            df["total_score"] = pd.to_numeric(df["total_score"], errors="coerce")
            df["rmsd_ca_to_input_trace"] = pd.to_numeric(df["rmsd_ca_to_input_trace"], errors="coerce")
            out[m] = df
    return out


def main() -> None:
    target_id = "1r69"
    tm = _load_tm(target_id)
    best = _load_rosetta_best(target_id)
    per = _load_per_struct(target_id)

    rows = [
        Row("omega", tm["omega"], best["omega"][0], best["omega"][1]),
        Row("quark", tm["quark"], best["quark"][0], best["quark"][1]),
        Row("itasser", tm["itasser"], best["itasser"][0], best["itasser"][1]),
    ]

    # consistent ordering
    methods = [r.method for r in rows]
    tm_vals = [r.tm for r in rows]
    rmsd_vals = [r.best_rmsd for r in rows]
    energy_vals = [r.best_energy for r in rows]

    colors = {"omega": "#F28E2B", "quark": "#4E79A7", "itasser": "#59A14F"}
    cs = [colors[m] for m in methods]

    fig = plt.figure(figsize=(10, 4.2), dpi=180)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.25], wspace=0.28)

    # Panel A: TM bars
    axA = fig.add_subplot(gs[0, 0])
    axA.bar(methods, tm_vals, color=cs, edgecolor="black", linewidth=0.6)
    axA.set_ylim(0.0, 1.0)
    axA.set_ylabel("TM-score vs native")
    axA.set_title("A  Topological accuracy (1R69)")
    for i, v in enumerate(tm_vals):
        axA.text(i, v + 0.02, f"{v:.2f}", ha="center", va="bottom", fontsize=9)

    # Panel B: Physical stability audit (prefer per-struct drift distribution)
    axB = fig.add_subplot(gs[0, 1])
    x = range(len(methods))
    axB.set_xticks(list(x), methods)
    axB.set_ylabel("Cα RMSD drift to input trace (Å)")
    axB.set_title("B  Physical audit under constrained FastRelax (n=20)")
    axB.grid(True, axis="y", alpha=0.25)

    if all(m in per for m in methods):
        data = [per[m]["rmsd_ca_to_input_trace"].dropna().values for m in methods]
        bp = axB.boxplot(
            data,
            positions=list(x),
            widths=0.55,
            patch_artist=True,
            showfliers=False,
            medianprops={"color": "black", "linewidth": 1.2},
            whiskerprops={"color": "black", "linewidth": 1.0},
            capprops={"color": "black", "linewidth": 1.0},
        )
        axB.set_xticks(list(x), methods)
        for patch, c in zip(bp["boxes"], cs):
            patch.set_facecolor(c)
            patch.set_alpha(0.55)
            patch.set_edgecolor("black")
            patch.set_linewidth(1.0)
        # overlay best-of-20 points for drift
        axB.scatter(x, rmsd_vals, s=45, c=cs, edgecolors="black", linewidths=0.7, zorder=3)
        for i, rmsd in enumerate(rmsd_vals):
            axB.text(i, rmsd + 0.02, f"best {rmsd:.2f}Å", ha="center", va="bottom", fontsize=9)
    else:
        axB.scatter(x, rmsd_vals, s=55, c=cs, edgecolors="black", linewidths=0.6)
        for i, rmsd in enumerate(rmsd_vals):
            axB.text(i, rmsd + 0.03, f"{rmsd:.2f}Å", ha="center", va="bottom", fontsize=9)

    out_dir = (repo_root() / RUN_DIR / "figures").resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    out_png = out_dir / "fig5_panelAB_1r69.png"
    fig.savefig(out_png, bbox_inches="tight")
    plt.close(fig)
    print(f"[ok] wrote: {out_png}")


if __name__ == "__main__":
    main()

