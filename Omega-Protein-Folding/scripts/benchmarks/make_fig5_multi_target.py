#!/usr/bin/env python3
from __future__ import annotations

"""
Figure 5 (multi-target): PhasonFold vs QUARK vs I-TASSER physical audit.

Panels:
- A: Kabsch single-pass TM-score bars for 1R69 and 2CRO (model1; from server_summaries.csv).
- B: Drift distributions under coordinate-constrained FastRelax (nstruct=20):
     Cα RMSD to each method's input trace (from committable per-struct CSVs).
- C: Qualitative Cα overlay vs native for 1R69 (pre-rendered panel).

Outputs:
- docs/runs/quark_itasser_homology_ablation/figures/fig5_multi_target.png
- paper/assets/figures/fig5_multi_target.png
"""

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


RUN_DIR = Path("docs/runs/quark_itasser_homology_ablation")


def _read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as f:
        r = csv.DictReader(f)
        return [dict(row) for row in r]


def _load_tm_kabsch(target_ids: list[str]) -> dict[str, dict[str, float]]:
    """
    Return target_id -> method -> tm_score_to_native for source == tm_kabsch.
    """
    p = (repo_root() / RUN_DIR / "server_summaries.csv").resolve()
    rows = _read_csv_rows(p)
    out: dict[str, dict[str, float]] = {}
    for tid in target_ids:
        out[tid.lower()] = {}
    for row in rows:
        tid = str(row.get("target_id", "")).strip().lower()
        if tid not in out:
            continue
        if str(row.get("source", "")).strip().lower() != "tm_kabsch":
            continue
        m = str(row.get("method", "")).strip().lower()
        try:
            out[tid][m] = float(row.get("tm_score_to_native", "nan"))
        except Exception:
            continue
    return out


def _load_drift_per_struct(target_id: str, nstruct: int) -> dict[str, np.ndarray]:
    """
    Return method -> drift array for the given target and nstruct.
    """
    run_dir = (repo_root() / RUN_DIR).resolve()
    out: dict[str, np.ndarray] = {}
    for m in ["omega", "quark", "itasser"]:
        p = run_dir / f"rosetta_relax_{m}_{target_id.lower()}_n{int(nstruct)}_per_struct.csv"
        if not p.exists():
            continue
        rows = _read_csv_rows(p)
        vals: list[float] = []
        for row in rows:
            try:
                v = float(row.get("rmsd_ca_to_input_trace", "nan"))
            except Exception:
                continue
            if np.isfinite(v):
                vals.append(v)
        out[m] = np.asarray(vals, dtype=np.float64)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", default="1r69,2cro")
    ap.add_argument("--nstruct", type=int, default=20)
    args = ap.parse_args()

    target_ids = [t.strip().lower() for t in str(args.targets).split(",") if t.strip()]
    nstruct = int(args.nstruct)

    methods = ["omega", "quark", "itasser"]
    colors = {"omega": "#F28E2B", "quark": "#4E79A7", "itasser": "#59A14F"}
    method_label = {"omega": "PhasonFold", "quark": "QUARK", "itasser": "I-TASSER"}

    # Panel A data
    tm = _load_tm_kabsch(target_ids)
    if not all(t in tm for t in target_ids):
        raise SystemExit("Missing TM entries for some targets.")
    for t in target_ids:
        missing = [m for m in methods if m not in tm[t]]
        if missing:
            raise SystemExit(f"Missing TM(method) for target={t}: {missing}")

    # Panel B data
    drift: dict[str, dict[str, np.ndarray]] = {}
    for t in target_ids:
        per = _load_drift_per_struct(t, nstruct)
        missing = [m for m in methods if m not in per or per[m].size == 0]
        if missing:
            raise SystemExit(f"Missing per-struct drift for target={t} n={nstruct}: {missing}")
        drift[t] = per

    # Panel C image (pre-rendered overlay on 1R69)
    overlay = (repo_root() / RUN_DIR / "figures" / "fig5_panelC_overlay_1r69.png").resolve()
    if not overlay.exists():
        raise SystemExit(f"Missing overlay panel: {overlay}")

    fig = plt.figure(figsize=(14.2, 4.8), dpi=220)
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.45, 1.10], wspace=0.22)

    # --- Panel A: grouped TM bars across targets ---
    axA = fig.add_subplot(gs[0, 0])
    x = np.arange(len(methods), dtype=np.float64)
    width = 0.36 if len(target_ids) == 2 else max(0.18, 0.8 / max(1, len(target_ids)))
    offsets = np.linspace(-(len(target_ids) - 1) * width / 2.0, (len(target_ids) - 1) * width / 2.0, len(target_ids))

    for j, (t, off) in enumerate(zip(target_ids, offsets)):
        vals = [tm[t][m] for m in methods]
        axA.bar(x + off, vals, width=width, label=t.upper(), alpha=0.95, edgecolor="black", linewidth=0.6)

    axA.set_ylim(0.0, 1.0)
    axA.set_xticks(x, [method_label.get(m, m) for m in methods])
    axA.set_ylabel("TM-score vs native")
    axA.set_title("A  Topological accuracy (Kabsch single-pass TM)")
    axA.grid(True, axis="y", alpha=0.2)
    axA.legend(frameon=False, fontsize=9, loc="upper left")

    # --- Panel B: drift distributions across targets ---
    axB = fig.add_subplot(gs[0, 1])
    positions: list[float] = []
    data: list[np.ndarray] = []
    box_colors: list[str] = []
    xtick_pos: list[float] = []
    xtick_lab: list[str] = []

    group_gap = 1.15
    inner_gap = 0.75
    pos = 0.0
    for t in target_ids:
        start = pos
        for m in methods:
            positions.append(pos)
            data.append(drift[t][m])
            box_colors.append(colors[m])
            xtick_pos.append(pos)
            xtick_lab.append(method_label.get(m, m))
            pos += inner_gap
        mid = (start + (pos - inner_gap)) / 2.0
        axB.text(mid, -0.07, t.upper(), ha="center", va="top", transform=axB.get_xaxis_transform(), fontsize=10)
        pos += group_gap

    bp = axB.boxplot(
        data,
        positions=positions,
        widths=0.55,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "black", "linewidth": 1.1},
        whiskerprops={"color": "black", "linewidth": 0.9},
        capprops={"color": "black", "linewidth": 0.9},
    )
    for patch, c in zip(bp["boxes"], box_colors):
        patch.set_facecolor(c)
        patch.set_alpha(0.55)
        patch.set_edgecolor("black")
        patch.set_linewidth(0.9)

    # best-of-nstruct points
    best_vals = [float(np.min(d)) for d in data]
    axB.scatter(positions, best_vals, s=28, c=box_colors, edgecolors="black", linewidths=0.6, zorder=3)

    axB.set_xticks(xtick_pos, xtick_lab, rotation=0)
    axB.set_ylabel("Cα RMSD drift to input trace (Å)")
    axB.set_title(f"B  Physical audit under constrained FastRelax (n={nstruct})")
    axB.grid(True, axis="y", alpha=0.2)

    # --- Panel C: overlay image ---
    axC = fig.add_subplot(gs[0, 2])
    img = plt.imread(str(overlay))
    axC.imshow(img)
    axC.axis("off")
    axC.set_title("C  1R69 Cα overlay (aligned to native)", fontsize=11)

    out_dir = (repo_root() / RUN_DIR / "figures").resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    out_png = out_dir / "fig5_multi_target.png"
    fig.savefig(out_png, bbox_inches="tight")
    plt.close(fig)
    print(f"[ok] wrote: {out_png}")

    # Copy to paper asset folder (so LaTeX build uses a committable path).
    paper_dir = (repo_root() / "paper" / "assets" / "figures").resolve()
    paper_dir.mkdir(parents=True, exist_ok=True)
    paper_png = paper_dir / "fig5_multi_target.png"
    paper_png.write_bytes(out_png.read_bytes())
    print(f"[ok] copied: {paper_png}")


if __name__ == "__main__":
    main()

